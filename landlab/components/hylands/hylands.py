
import numpy as np
import os
import numpy.matlib
import copy as cp
from landlab import Component
from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import non_local_Depo
# from .nonLocalDeposition_loops import non_local_Depo

from landlab import imshow_grid
from matplotlib import pyplot as plt
from matplotlib.pyplot import pause       

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

    _info = {
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
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "LS__Bed_E": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "LS__Bed_E",
        },
        "LS__Bed_D": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "LS__Bed_E",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
    }

    _cite_as = """    """

    def __init__(
        self,
        grid,
        LS_properties_Size,
        LS_properties_Volume,
        LS_properties_Volume_Reg,
        LS_properties_Volume_Bed,   
        slope_crit = 1.0,
        C_eff = 1e4,
        rho_r = 2700,
        grav = 9.81,
        Ff_LS = 0,
        phi = 0,
        LS_returnTime = 1e5,
        maxPixels = 1e9, 
        if_unstable = "pass",
        output_prefix="terrainbento_output",
        output_dir=os.getcwd(),
        figure_dir=os.getcwd(),
        plot_each = 1e9,
        verbose_LS = False,
        max_H_el = 100,
    ):
        """Initialize the Space model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        slope_crit: float, optional
            Critical slope
            Default = 1.0
        C_eff : float
            Effective cohesian of material [m L^-1 T^-2].
        rho_r : float
            Bulk density rock [m L^-3].
        rho_s : float
            Bulk density sediment [m L^-3].
        F_f_{LS} : float
            Fraction of permanently suspendable fines in bedrock [-].
        phi : float
            Sediment porosity [-].
        if_unstable : string (optional, default is "pass")
            Keyword argument to determine how potential instability due to
            slopes that are too high is handled. Options are "pass", "warn",
            and "raise".

        """
        super(Hylands, self).__init__(
            grid)
        self.LS_properties_Size=LS_properties_Size
        self.LS_properties_Volume=LS_properties_Volume
        self.LS_properties_Volume_Reg=LS_properties_Volume_Reg
        self.LS_properties_Volume_Bed=LS_properties_Volume_Bed
        
        # Store grid and parameters
        self._slope_crit = slope_crit
        self._C_eff = C_eff
        self._rho_r = rho_r
        self._grav = grav
        self._Ff_LS = Ff_LS
        self._phi = phi
        self._LS_returnTime = LS_returnTime
        self._maxPixels = maxPixels
        self._if_unstable = if_unstable
        self._out_file_name = output_prefix
        self._out_file_dir = output_dir 
        self._out_figure_dir = figure_dir
        self._plot_each = plot_each 
        self._verbose_LS = verbose_LS
        self._max_H_el = max_H_el
        self._Qs_change_LS =np.zeros(self.grid.number_of_nodes)
        


       
     

    def LS_Cullman(self, dt, dx2):
        
        # Reset LS Plains
        self.grid.at_node['LS__Bed_E'] = np.zeros(self.grid.number_of_nodes)
        self.grid.at_node['LS__Bed_D'] = np.zeros(self.grid.number_of_nodes)
               
       
        
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
        
        H_el[H_el>self._max_H_el]=self._max_H_el    
        
        beta_rad = np.arctan(S_slope)
        sc_rad = np.arctan(self._slope_crit)
        Hc = (4 * self._C_eff/(self._grav*self._rho_r))\
            *np.divide(np.multiply(np.sin(beta_rad),np.cos(sc_rad)),(1-np.cos(beta_rad-sc_rad)))
        
            
     
        
        # %% Spatial probability v2
        p_S = np.divide(H_el,Hc)
        p_S[beta_rad<=sc_rad] = 0
        p_S[p_S>1]=1
        p_L = p_S/self._LS_returnTime
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
        if self._verbose_LS:
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
            
            nb_up=nb_up[s_slide_all>self._slope_crit]
            s_slide_all=s_slide_all[s_slide_all>self._slope_crit]
            
            if s_slide_all.size > 0:
                s_slide=max(s_slide_all)
                storeV_bed = np.float64(0.0)
                storeV_sed = np.float64(0.0)
                upstream = 0
                uP=nb_up
                uP = uP[(self.grid.node_is_boundary(uP) == False)] 
                # Fix sliding angle of particular LS
                ang_sl=np.float64((self._slope_crit+s_slide)/2)
                
                stall = 0
                
                
                while uP.size > 0 & (upstream<=self._maxPixels and stall<1e4) :
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
                        
                        self.grid.at_node['LS__Bed_E'][uP[0]]=sed_LS_E+bed_LS_E                        
                       
                        
                    uP=np.delete(uP, 0, 0)
                
                storeV=storeV_sed+storeV_bed
                storeV_cum += storeV
                if upstream > 0:
                    self._Qs_change_LS[cP] += np.float64((storeV/dt)*(1-self._Ff_LS))
                    suspended_Sed += np.float64((storeV/dt)*self._Ff_LS)
                    self.LS_properties_Size.append(upstream)
                    self.LS_properties_Volume.append(storeV)
                    self.LS_properties_Volume_Reg.append(storeV_sed)
                    self.LS_properties_Volume_Bed.append(storeV_bed)    
                
            if i_slide.size > 0:
                i_slide=np.delete(i_slide, 0, 0)   
            

        return suspended_Sed
        
   
    
    def deposition(self,dt,dx2):
        
                
        # %% calcualte hillslope deposition 
        # Aliases
        stack = self.grid.at_node["flow__upstream_node_order"]
        receivers = self.grid.at_node["hill_flow__receiver_node"]
        fract = self.grid.at_node["hill_flow__receiver_proportions"]
        slope = cp.deepcopy(self.grid.at_node["hill_topographic__steepest_slope"])
        # el = cp.deepcopy(self.grid.at_node["topographic__elevation"])
        
        z = self.grid.at_node["topographic__elevation"]
        br = self.grid.at_node["bedrock__elevation"]
        H = self.grid.at_node["soil__depth"]
        
        
        # keep only steepest slope
        slope = np.max(slope,axis=1)
        slope = np.maximum(slope,np.zeros(slope.shape))
        
        #core_nodes = DepoGrid_v1.core_nodes
                              
        Qs_in = self._Qs_change_LS*dt# Qs_in, in m3 per timestep
        
        # L following carretier 
        dx = self.grid.dx
        dx2 = np.square(dx)
        
        L_Hill = numpy.matlib.divide(dx,(1-numpy.matlib.minimum( numpy.matlib.square(numpy.matlib.divide(slope,self._slope_crit )),0.999)))        
        Qs_out = np.zeros(z.shape)
        # Node that this results in zero deposition at watershed divides
        dH_Hill = np.zeros(z.shape)
        H_i_temp = cp.deepcopy(z)
        max_D = np.zeros(z.shape)
        max_dH = np.ones(z.shape)+np.inf
        
        phi = self._phi
        
        
        # Number of points
        nb = receivers.shape[0]
        nb_r = receivers.shape[1]
        # Call the cfunc to calculate non-local, non-linear deposition
        Grid_Status = np.float64(self.grid.status_at_node)        
       
     
        non_local_Depo(nb,
                   nb_r,
                   dx,
                   dx2,
                   phi,
                   stack,
                   receivers,
                   fract,
                   Qs_in,
                   L_Hill,
                   Qs_out,
                   dH_Hill,
                   H_i_temp,
                   max_D,
                   max_dH,
                   Grid_Status)
            
        Qs_coreNodes = np.sum(Qs_in[self.grid.status_at_node == 0])#/(self.grid.number_of_nodes*dx2)*dt # This value should be close to zero, difference with zero is error
        Qs_leaving = np.sum(Qs_in[self.grid.status_at_node == 1])
        V_leaving = np.sum(Qs_in) #Qs_leaving # in m3 per timestep
        
     
        
        # Change sediment layer
        H[:] += dH_Hill[:]
        z[:] = br[:] + H[:]        
        
        # Reset Qs
        self._Qs_change_LS =np.zeros(self.grid.number_of_nodes)    
        
        self.grid.at_node['LS__Bed_D'][:]   =   dH_Hill
        
        
        return dH_Hill,V_leaving,Qs_coreNodes
        


    def run_one_step(self, dt, dx2):
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
        vol_SSY = self.LS_Cullman(dt,dx2)     
        dH_Hill,V_leaving,Qs_coreNodes = self.deposition(dt,dx2)   
        
        return vol_SSY, V_leaving, self.LS_properties_Size,self.LS_properties_Volume,self.LS_properties_Volume_Reg,self.LS_properties_Volume_Bed
        
            
            


 
