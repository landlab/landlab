import numpy as np
import numpy.matlib
import copy as cp
from landlab import Component
from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import _landslide_runout
from landlab.grid.nodestatus import NodeStatus

MAX_HEIGHT_SLOPE = 100 # in m 

class Hylands(Component):
    """Hybrid Landscape evolution: landslide component

    See the publication:
    

    Note: 

    Examples
    ---------
     None Listed
   

    References
    ----------
    **Required Software Citation(s) Specific to this Component**
     None Listed


    **Additional References**

    None Listed

    """

    _name = "Hylands"    
    
    _unit_agnostic = True

    _info = {
        
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },        
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "sediment__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m3/s",
            "mapping": "node",
            "doc": "Sediment flux (volume per unit time of sediment entering each node)",
        },
        "landslide__bed_erosion": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "landslide__bed_erosion",
        },
        "landslide__deposition": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "landslide__bed_erosion",
        }
        
    }

    _cite_as = """@Article{gmd-10-4577-2017,
                  AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
                  TITLE = {HyLands 1.0: a hybrid landscape evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution.},
                  JOURNAL = {Geoscientific Model Development},
                  VOLUME = {13},
                  YEAR = {2020},
                  NUMBER = {9},
                  PAGES = {3863--3886},
                  URL = {https://doi.org/10.5194/gmd-13-3863-2020},
                  DOI = {10.5194/gmd-13-3863-2020}
                  }"""

    def __init__(
        self,
        grid,
        angle_int_frict = 1.0,
        cohesion_eff = 1e4,
        landslides_return_time = 1e5,
        rho_r = 2700,
        grav = 9.81,
        fraction_fines_LS = 0,
        phi = 0,
        max_pixelsize_landslide = 1e9, 
        verbose_landslides = False,        
        landslides_size = None,
        landslides_volume = None,
        landslides_volume_sed = None,
        landslides_volume_bed = None,
    ):
        """Initialize the HyLands model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        angle_int_frict: float, optional
            Materials angle of internal friction in [m/m]
            Default = 1.0
        cohesion_eff : float
            Effective cohesian of material [m L^-1 T^-2].
        landslides_return_time  : float
            Return time for stochastic landslide events to occur
            Default = 1e5
        rho_r : float
            Bulk density rock [m L^-3].
        fraction_fines_LS : float
            Fraction of permanently suspendable fines in bedrock [-].
        phi : float
            Sediment porosity [-].        
        max_pixelsize_landslide : int 
            Maximum size for landslides in number of pixels
            Default = 1e9
        verbose_landslides : bool 
            Print output as number of simulated landslides per timestep
            Default = False       
        landslides_size : list  
            Store the size of landslides through time
            Pass empty list ('[]') to store sizes, pass None to not store any data
            Default = None
        landslides_volume : list 
            Store the volume of landslides through time.
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None
        landslides_volume_sed  : list 
            Store the volume of eroded bedrock through time
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None
        landslides_volume_bed : list  
            Store the volume of eroded sediment landslides through time
            Pass empty list ('[]') to store volumes, pass None to not store any data
            Default = None

        """
        super(Hylands, self).__init__(
            grid)
        
        # Store grid and parameters
        self._angle_int_frict = angle_int_frict
        self._cohesion_eff = cohesion_eff
        self._rho_r = rho_r
        self._grav = grav
        self._fraction_fines_LS = fraction_fines_LS
        self._phi = phi
        self._landslides_return_time = landslides_return_time
        self._max_pixelsize_landslide = max_pixelsize_landslide
        self._verbose_landslides = verbose_landslides
        self._Qs_change_LS =np.zeros(self.grid.number_of_nodes)
        
        # Data structure to store properties of simulated landslides. 
        # If None, data will not be stored
        self.landslides_size=landslides_size
        self.landslides_volume=landslides_volume
        self.landslides_volume_sed=landslides_volume_sed
        self.landslides_volume_bed=landslides_volume_bed    

    def _landslide_erosion(self, dt):
        
        
        
        dx = self.grid.dx
        dx2 = np.square(dx)
        
        # Reset LS Plains
        self.grid.at_node['landslide__bed_erosion'] = np.zeros(self.grid.number_of_nodes)
        self.grid.at_node['landslide__deposition'] = np.zeros(self.grid.number_of_nodes)
               
       
        
        """Calculate bedrock landsliding for a time period 'dt'.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        
        # identify flooded nodes
        flood_status = self.grid.at_node["flood_status_code"]
        flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
      

        # Calculate gradients
        S_slope=self.grid.at_node['topographic__steepest_slope']
        
        H_el=cp.deepcopy(self.grid.at_node['topographic__elevation']\
            -self.grid.at_node['topographic__elevation'][self.grid.at_node['flow__receiver_node']])
        H_el[flooded_nodes]=0    
        
        H_el[H_el>MAX_HEIGHT_SLOPE] = MAX_HEIGHT_SLOPE 
        
        beta_rad = np.arctan(S_slope)
        sc_rad = np.arctan(self._angle_int_frict)
        Hc = (4 * self._cohesion_eff/(self._grav*self._rho_r))\
            *np.divide(np.multiply(np.sin(beta_rad),np.cos(sc_rad)),(1-np.cos(beta_rad-sc_rad)))
        
            
     
        
        # %% Spatial probability 
        p_S = np.divide(H_el,Hc)
        p_S[beta_rad<=sc_rad] = 0
        p_S[p_S>1]=1
        p_L = p_S/self._landslides_return_time
        p =1-np.exp(-p_L*dt)
        slides = np.random.rand(p.size)<p
          
        # Now, find the critical node, which is the receiver of i_slide
        # Critical nodes must be unique (a given node can have more receivers...)
        i_slide=np.unique(self.grid.at_node['flow__receiver_node'][np.where(slides)])
        # Remove boundary nodes
        i_slide = i_slide[(self.grid.node_is_boundary(i_slide) == False)]
       
        
        # %%     
        # output variables
        suspended_Sed = np.float64(0)     
        if self._verbose_landslides:
            print('nbSlides = ' + str( np.sum(i_slide)))
        
        
        storeV_cum=0.0    
        while i_slide.size > 0:   
            
        
            ind = 0
            cP=i_slide[ind] #Critical Node
            cP_el=self.grid.at_node['topographic__elevation'][cP]
            
            # get 8 neighbors and only keep those to active nodes which are upstream
            nb=np.concatenate((self.grid.active_adjacent_nodes_at_node[cP],
                      self.grid.diagonal_adjacent_nodes_at_node[cP]))           
            nb=nb[nb!=-1]            
            nb_up=nb[self.grid.at_node['topographic__elevation'][nb]>cP_el];
            
            X_cP = self.grid.node_x[cP]
            Y_cP = self.grid.node_y[cP]
            
            distToIni_all= np.sqrt(np.add(np.square(X_cP-self.grid.node_x[nb_up]),
                                np.square(Y_cP-self.grid.node_y[nb_up])))
            all_iP_el=self.grid.at_node['topographic__elevation'][nb_up]
            s_slide_all=np.divide((all_iP_el-cP_el),distToIni_all)
            
            nb_up=nb_up[s_slide_all>self._angle_int_frict]
            s_slide_all=s_slide_all[s_slide_all>self._angle_int_frict]
            
            if s_slide_all.size > 0:
                s_slide=max(s_slide_all)
                storeV_bed = np.float64(0.0)
                storeV_sed = np.float64(0.0)
                upstream = 0
                uP=nb_up
                uP = uP[(self.grid.node_is_boundary(uP) == False)] 
                # Fix sliding angle of particular LS
                ang_sl=np.float64((self._angle_int_frict+s_slide)/2)
                
                stall = 0
                
                
                while uP.size > 0 & (upstream<=self._max_pixelsize_landslide and stall<1e4) :
                    #print(uP.size)
                    distToIni= np.sqrt(np.add(np.square(X_cP-self.grid.node_x[uP[0]]),
                                    np.square(Y_cP-self.grid.node_y[uP[0]])))                
                    newEl=cP_el+distToIni*ang_sl
                    stall+=1
                    if newEl<self.grid.at_node['topographic__elevation'][uP[0]]:
                        # Do actual slide
                        upstream=upstream+1;
                        sed_LS_E=np.float64(max(0.0,min(self.grid.at_node["soil__depth"][uP[0]],self.grid.at_node['topographic__elevation'][uP[0]]-newEl)))
                        self.grid.at_node["soil__depth"][uP[0]] -= sed_LS_E
                        self.grid.at_node['topographic__elevation'][uP[0]] = newEl
                        bed_LS_E = max(0, self.grid.at_node['bedrock__elevation'][uP[0]] - (newEl - self.grid.at_node["soil__depth"][uP[0]]))
                        self.grid.at_node['bedrock__elevation'][uP[0]] -= bed_LS_E
                        self.grid.at_node['topographic__elevation'][uP[0]] = newEl
                        
                        vol_sed = sed_LS_E*(1-self._phi)*dx2
                        vol_bed = bed_LS_E*dx2                            
                        storeV_sed=storeV_sed+vol_sed
                        storeV_bed=storeV_bed+vol_bed
                        
                        
                        
                        nb=np.concatenate((self.grid.active_adjacent_nodes_at_node[uP[0]],
                          self.grid.diagonal_adjacent_nodes_at_node[uP[0]]))           
                        nb=nb[nb!=-1]            
                        nb_up=nb[self.grid.at_node['topographic__elevation'][nb]>cP_el]                       
                        uP=[*uP, *nb_up]
                        
                        temp,idx = np.unique(uP,'first')   
                        uP = np.array(uP) 
                        uP = uP[np.sort(idx)]
                        uP = uP[(self.grid.node_is_boundary(uP) == False)] 
                        # if one of the LS pixels also appears in i_slide list, 
                        # remove it there so that no new landslide is initialized 
                        i_slide=i_slide[np.where((i_slide!=uP[0]))]           
                        
                        self.grid.at_node['landslide__bed_erosion'][uP[0]]=sed_LS_E+bed_LS_E                        
                       
                        
                    uP=np.delete(uP, 0, 0)
                
                storeV=storeV_sed+storeV_bed
                storeV_cum += storeV
                if upstream > 0:
                    self._Qs_change_LS[cP] += np.float64((storeV/dt)*(1-self._fraction_fines_LS))
                    suspended_Sed += np.float64((storeV/dt)*self._fraction_fines_LS)                    
                    
                    if self.landslides_size is not None:    
                        self.landslides_size.append(upstream)
                    if self.landslides_volume is not None: 
                        self.landslides_volume.append(storeV)
                    if self.landslides_volume_sed is not None: 
                        self.landslides_volume_sed.append(storeV_sed)
                    if self.landslides_volume_bed is not None: 
                        self.landslides_volume_bed.append(storeV_bed)    
                
            if i_slide.size > 0:
                i_slide=np.delete(i_slide, 0, 0)   
            

        return suspended_Sed
        
   
    
    def _landslide_runout(self,dt):
        
                
        # %% calcualte hillslope _landslide_runout 
        z = self.grid.at_node["topographic__elevation"]
        br = self.grid.at_node["bedrock__elevation"]
        H = self.grid.at_node["soil__depth"]
        stack_rev = np.flip(self.grid.at_node["flow__upstream_node_order"])
        node_status = self.grid.status_at_node
        # Only process core nodes
        stack_rev_sel = stack_rev[(node_status[stack_rev] == NodeStatus.CORE)]        
        receivers = self.grid.at_node["hill_flow__receiver_node"]
        fract = self.grid.at_node["hill_flow__receiver_proportions"]
        slope = cp.deepcopy(self.grid.at_node["hill_topographic__steepest_slope"])    
        
        # keep only steepest slope
        slope = np.max(slope,axis=1)
        slope = np.maximum(slope,np.zeros(slope.shape))
                                      
        Qs_in = self._Qs_change_LS*dt# Qs_in, in m3 per timestep
        
        # L following carretier 
        dx = self.grid.dx
        dx2 = np.square(dx)
        
        L_Hill = numpy.matlib.divide(dx,(1-numpy.matlib.minimum( numpy.matlib.square(numpy.matlib.divide(slope,self._angle_int_frict )),0.999)))        
        Qs_out = np.zeros(z.shape)
        # Node that this results in zero _landslide_runout at watershed divides
        dH_Hill = np.zeros(z.shape)
        H_i_temp = cp.deepcopy(z)
        max_D = np.zeros(z.shape)
        max_dH = np.ones(z.shape)+np.inf
        
        phi = self._phi
        
        
        # Number of receivers
        nb_r = receivers.shape[1]
        
        # Call the cfunc to calculate non-local, non-linear _landslide_runout
        Grid_Status = np.float64(self.grid.status_at_node)        
       
     
        _landslide_runout(                   
            dx,
            phi,
            stack_rev_sel,
            receivers,
            fract,
            Qs_in,
            L_Hill,
            Qs_out,
            dH_Hill,
            H_i_temp,
            max_D,
            max_dH)
            
        Qs_coreNodes = np.sum(Qs_in[self.grid.status_at_node == 0])#/(self.grid.number_of_nodes*dx2)*dt # This value should be close to zero, difference with zero is error
        Qs_leaving = np.sum(Qs_in[self.grid.status_at_node == 1])
        V_leaving = np.sum(Qs_in) #Qs_leaving # in m3 per timestep
        
     
        
        # Change sediment layer
        H[:] += dH_Hill[:]
        z[:] = br[:] + H[:]        
        
        # Reset Qs
        self._Qs_change_LS =np.zeros(self.grid.number_of_nodes)    
        
        self.grid.at_node['landslide__deposition'][:]   =   dH_Hill
        
        
        return dH_Hill,V_leaving,Qs_coreNodes
        


    def run_one_step(self, dt):
        """Advance cubic soil flux component by one time step of size dt.

        Parameters
        ----------
        dt: float (time)
            The imposed timestep.
        """
        
        if self.current_time is None:
            self.current_time=dt
        else:            
            self.current_time+=dt
        
        dt = np.float64(dt)

        #Landslides
        vol_SSY = self._landslide_erosion(dt)     
        dH_Hill,V_leaving,Qs_coreNodes = self._landslide_runout(dt)   
        
        return vol_SSY, V_leaving
            
            


 
