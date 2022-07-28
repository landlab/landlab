"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab import Component
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import warnings

class PlantGrowth(object):
    """
    Add Intro Stuff here
    """
    _name = "PlantGrowth"

    _unit_agnostic = False


    _cite_as = """
    @article{piercygv,
        author = {Piercy, C.D.; Swannack, T.M.; Carrillo, C.C.; Russ, E.R.; Charbonneau, B. M.]
    }
    """
    #Add all variables to be saved or chosen as outputs here
    _info = {
        "vegetation__total_biomass": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units":"g",
            "mapping":"cell",
            "doc": "Total plant biomass for the plant class at the end of the time step"
        },
    }

    def __init__(
        self,
        grid,
        vegparams={'Corn': 
            {
                'colparams': {}, 
                'dispparams': {}, 
                'growparams': {
                    'growing_season_start': 91,
                    'growing_season_end': 290,
                    'senescence_start': 228,
                    'respiration_coefficient': [0.015,0.015,0.03],
                    'glucose_requirement': [1.444,1.513,1.463],
                    'k_light_extinct':0.02,
                    'light_half_sat':9,
                    'p_max':0.055,
                    'root_to_leaf_coeffs': [0.031,0.951,0],
                    'root_to_stem_coeffs': [-0.107, 1.098, 0.0216],
                    'plant_part_min':[0.01,0.1,0.5]
                }, 
                'mortparams': {
                    's1_days': 365, 
                    's1_name': 'Mortality factor', 
                    's1_pred': [1, 2, 3, 4], 
                    's1_rate': [0, 0.1, 0.9, 1], 
                    's1_weight': [1000, 1, 1, 1000]
                },
                'plant_factors':{
                    'species':'Corn',
                    'growth_form': 1,
                    'monocot_dicot': 'monocot',
                    'angio_gymno': 'angiosperm',
                    'annual_perennial': 'annual',
                    'ptype':'C3'
                }, 
                'sizeparams': {
                    'max_height_stem': 2.5, 
                    'max_mass_stem': 72, 
                    'max_n_stems': 3, 
                    'max_plant_density': 1
                },
                'storparams': {
                    'r_wint_die': 0.25, 
                    'r_wint_stor': 0.25
                    } 
                }
        },
        dt=1,
        plants=pd.DataFrame(columns=['pid','species','cell_index'])
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
        grid: RasterModelGrid
            A Landlab ModelGrid

        vegparams: dict, required,
            a nested dictionary of named vegetation parameters for each
            species or community and process of interest with below sub-dictionaries.
            plant_factors: dict, required,
                dictionary of plant characteristics describing the
                species or community of interest with below keys
                species: string, required,
                    name of species or community used to identify plant
                growth_form: string, required,
                    USDA plant growth habit,
                    graminoid, forb/herb, shrub, tree, nonvascular
                monocot_dicot: string, required,
                    should be monocot or dicot
                angio_gymno: string, required,
                    should be angiosperm or gymnosperm
                annual_perennial: string, required,
                    plant growth duration, annual (1 year) or
                    perennial (multiple years)
            growparams: dict, required,
                dictionary of paramaters required to simulate plant growth
                ptype: string, required,
                    photosythesis type, either 'C3', 'C4', or 'CAM'
                gs_start: int, required,
                    growing season start day of year,
                    must be between 1-365
                gs_end: int, required,
                    growing season end day of year,
                    must be between 1-365
                gs_length: int, required,
                    growing season length in days,
                    must be less than or equal to 365
                senes: int, required,
                    start of senescence period after plant reaches peak biomass,
                    must be between gs_start and gs_end
                res_co: float, required,
                    respiration coefficient
                glu_req: float, required,
                    glucose requirement
                le_k: float, required,
                    light extinction coefficient
                hi: float, required,
                    something
                p_max: float, required,
                    maximum photosyntehtic output

        dt: int, required,
            time step interval

        plants: Pandas DataFrame of plant objects, optional
            with columns
            pid: plant ID
            species: plant species code
            cell_index: index of cell location on grid
            leaf_biomass: plant live leaf biomass in g
            stem_biomass: plant live stem biomass in g
            root_biomass: plant live root biomass in g
        """
        #This will move to GenVeg class
        self._grid=grid
        
        # From chapter 2 of Teh 2006 (pg. 31; equation 2.13)
        (_,self._latitude)=self._grid.xy_of_reference
        self._lat_rad = np.radians(self._latitude)
        
        self.vegparams=vegparams
        self.dt=dt
        self.plants=plants
        
        if self.plants.empty:
            self._estimate_initial_plant_biomass()
 
        #Set constants for PAR formula
        self._wgaus = [0.2778, 0.4444, 0.2778]
        self._xgaus = [0.1127, 0.5, 0.8873]

        #Calculate calculated parameters for speed and error check growth parameter values
        for item in self.vegparams:
            growdict=self.vegparams[item]['growparams']
            #Check validity of growing season beginning and ending days, senescence beginning, and calculate growing season length
            end=growdict['growing_season_end']
            begin=growdict['growing_season_start']
            senes=growdict['senescence_start']
            if (begin > 0) & (begin  < 366) & (end > 0) & (end < 366) & (end > begin) & (senes > begin) & (senes < end):
                length=end-begin+1
            elif (begin < 1) | (begin > 365):
                msg='Growing season beginning must be between 1-365'
                raise ValueError(msg)
            elif (end < 1) | (end > 365):
                msg='Growing season end must be between 1-365'
                raise ValueError(msg)
            elif end < begin:
                msg='Growing season beginning must be before growing season end'
                raise ValueError(msg)
            elif (senes < begin) | (senes > end):
                msg='Start of senescence must be within the growing season'
                raise ValueError(msg)
            else:
                msg='Growing season beginning and end must be integer values between 1-365'
                raise ValueError(msg)
            self.vegparams[item]['growparams']['growing_season_length']=length
            
            #Calculate initial plant part allocation - When does this happen? 
            #during initialization step, first day of growing season, or is this it's own function 
            #that can be called here or at start?

    def plant_df(self):
        """The plants dataframe"""
        return self.plants

    def _grow(self, current_day, dt):
        #Calculate current day-of-year
        jday_td=current_day-np.datetime64(str(current_day.astype('datetime64[Y]'))+'-01-01')
        self._current_jday=jday_td.astype(int)

        #Add photosythesis code here
        for species in self.vegparams:
            factordict=self.vegparams[species]['plant_factors']
            growdict=self.vegparams[species]['growparams']
            _last_veg_root_biomass=self.plants.root_biomass.loc[self.plants['species']==species].to_numpy()
            _last_veg_leaf_biomass=self.plants.leaf_biomass.loc[self.plants['species']==species].to_numpy()
            _last_veg_stem_biomass=self.plants.stem_biomass.loc[self.plants['species']==species].to_numpy()
            _last_veg_biomass=_last_veg_root_biomass+_last_veg_leaf_biomass+_last_veg_stem_biomass
            _gridradiation=pd.DataFrame(data=self._grid["cell"]["radiation__net_flux"],columns=['Radiation'])
            _gridradiation.reset_index(inplace=True)
            _gridradiation=_gridradiation.rename(columns={'index':'cell_index'})
            _rad_join=pd.merge(self.plants,_gridradiation)
            _radiation =_rad_join.Radiation.loc[_rad_join['species']==species].to_numpy()
            ##################################################
            #Calculate Vegetation Structure Metrics Each Day
            ##################################################
            
            #twlvd = total (t) weight(w) leaves(lv), d @ dead and g is living or green
            #totLeafWeight = twlvd + twlvg 
            #total weight of stems where d at end is dead and g is living or green
            #totStemWeight = twstd + twstg 
            #total weight of roots where d at end is dead and g is living or green
            #totRootWeight = twrtd + twrtg 

            ##################################################
            #Growth and Respiration
            ##################################################
                      
            #maintenance respiration
            #if then statement stops respiration at end of growing season
            if self._current_jday == growdict['growing_season_end']:
                respMaint = 0

            if self._current_jday < growdict['growing_season_end'] and self._current_jday >= growdict['growing_season_start']:

                #Calculate current daylength based on current day-of-year and latitude of grid origin
                _declination = np.radians(23.45) * (np.cos(2 * np.pi / 365 * (172 - self._current_jday)))
                _daylength = 24/np.pi*np.arccos(-(np.sin(_declination)*np.sin(self._lat_rad))/(np.cos(_declination)*np.cos(self._lat_rad)))

                #repiration coefficient for lvs, temp dependence from Teh 2006
                #change to read temperature off grid
                kmLVG = growdict['respiration_coefficient'][1] * pow(2,((_radiation - 25)/10))  
                #respiration coefficient for stems, temp depencence from Teh 2006 page 134
                kmSTG = growdict['respiration_coefficient'][2]* pow(2,((_radiation - 25)/10)) 
                #respiration coefficient for roots, temp dependence from Teh 2006 page 134
                kmRTG = growdict['respiration_coefficient'][0] * pow(2,((_radiation - 25)/10)) 
                #maintenance respiration per day from Teh 2006
                rmPrime = (kmLVG * _last_veg_leaf_biomass) + (kmSTG * _last_veg_stem_biomass) + (kmRTG * _last_veg_root_biomass)  
                #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration
                plantAge = _last_veg_biomass 
                #plant age dependence from Teh 2006 page 145    
                respMaint = rmPrime * plantAge  
                
                #glucose requirement for growth save to temp private variable since we are not saving to class
                #from Teh 2006 page 148
                _glu_req = ((_last_veg_leaf_biomass/_last_veg_biomass) * growdict['glucose_requirement'][1]) + \
                           ((_last_veg_stem_biomass/_last_veg_biomass) * growdict['glucose_requirement'][2]) + \
                           ((_last_veg_root_biomass/_last_veg_biomass) * growdict['glucose_requirement'][0]) 

                
                #writes results for daily respiration, plant age, and maintenance respiration
                #I don't know if we need to write these anywhere but we can write to a spearate test grid for now (NOT the main grid)
                #dailyrespiration.iloc[j] = rmPrime  #how to handle

                #dailyplantage[j] = plantAge    #how to handle

                #dailyrespMaint[j] = respMaint   #how to handle
                
                    
                    
            #Enter photosynthesis loop  
                #dailyphoto = []
                #determine how we want to read in light. Ideally this would be dynamic
                dtga=0
                #radiation measured 3x daily, roughly correlates to morning, noon, afternoon
                #change to read solar radiation from grid
                for hr in range(0,3):  
                    #convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on
                    parMicroE = (_radiation) * (868/208.32) #are we leaving as W/m2 or leaving microeinsteins??? FOR FUTURE!!!
                    #from Charisma instructions: tells how much of the light a plant is going to get as PAR in microeinsteins based on how many leaves are on the plant
                    intSolarRad = parMicroE*np.exp(-(growdict['k_light_extinct'])*_last_veg_leaf_biomass)  
                    #amount of light absorbed, per half saturaion constants from Charisma eq. 3. the monod or michaelis/menten function is adequate for describing the photosynthetic response to light
                    intLightpH = intSolarRad/(intSolarRad+growdict['light_half_sat']) 
                    #pMax is the maximum rate of photosynthesis, species specific
                    photosynthesis = (growdict['p_max']) * intLightpH 
                    #calculates gross assimilation of fgross(like APT) via photosynthesis at specific hour calculate growth per day at three times per day, morning, middday, and evenning and this amount is weighted based on how much light is hitting hte plant based on the latitude of your study site
                    #fgross.iloc[hr] = photosynthesis 
                    #weights fgross for specific time of day
                    dtgastep = photosynthesis*self._wgaus[hr]
                    dtga+=dtgastep
                        
                    
                    
                #calculates total biomass gained across plant (twlvg is amount of leaver/green matter): you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
                dtgaCollapsed = dtga*_last_veg_leaf_biomass  
                #total biomass for day length
                assimilatedCH2O = dtgaCollapsed*_daylength 
                #converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
                gphot = assimilatedCH2O*(30/44) 
                
                #I don't think we need to write all of this out except for testing
                ##how do we want to handle something like this?
                #dailyphoto.iloc[j] = gphot  just commenting out for easy purpose in case we need it later

                
                
                #if then statement ends glucose generation at end of growing season 
                if self._current_jday == growdict['growing_season_end']:
                    gphot = 0

                #direct solve method for calculating change in biomass
                #coefficients rename
                aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
                astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
                #Calculate the change in leaf and stem mass per unit change in root mass given the current size of the root biomass
                delta_leaf_unit_root=np.zeros_like(_last_veg_leaf_biomass)
                delta_leaf_unit_root[_last_veg_leaf_biomass!=0]=0.43429448190325176*b1leaf/_last_veg_root_biomass[_last_veg_root_biomass!=0] + \
                            0.37722339402322774*b2leaf*np.log(_last_veg_root_biomass[_last_veg_root_biomass!=0])/_last_veg_root_biomass[_last_veg_root_biomass!=0]
                delta_stem_unit_root=np.zeros_like(_last_veg_stem_biomass)
                delta_stem_unit_root[_last_veg_root_biomass!=0]=0.43429448190325176*b1stem/_last_veg_root_biomass[_last_veg_root_biomass!=0] + \
                            0.37722339402322774*b2stem*np.log(_last_veg_root_biomass[_last_veg_root_biomass!=0])/_last_veg_root_biomass[_last_veg_root_biomass!=0]
                #Calculate the total change in biomass this timestep
                delta_tot=(gphot-respMaint)/_glu_req
                #Calculate the change in root biomass this timestep
                #Solve for delta root: delta_tot=delta_root+delta_leaf_unit_root*delta_root+delta_stem_unit_root*delta_root
                delta_root=delta_tot/(1+delta_leaf_unit_root+delta_stem_unit_root)
                delta_leaf=delta_leaf_unit_root*delta_root
                delta_stem=delta_stem_unit_root*delta_root
                self.plants.root_biomass.loc[self.plants['species']==species]=_last_veg_root_biomass+delta_root
                self.plants.leaf_biomass.loc[self.plants['species']==species]=_last_veg_leaf_biomass+delta_leaf
                self.plants.stem_biomass.loc[self.plants['species']==species]=_last_veg_stem_biomass+delta_stem

        return self.plants

#May be able to eliminate this since using built in solar radiation component
    def _PAR(self, day, lat):
        #required to convert degrees to radians
        degree_to_rad = 0.017453292  
        #rad_to_degree = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))
    
        tmpvec = []
        #Change to use numpy as np 
        declination = (-np.asin((np.sin(23.45*degree_to_rad))*(np.cos(2*np.pi*(day+10)/365))))
        
        #Intermediate variables
        #radians
        sinld = ((np.sin(lat*degree_to_rad))*(np.sin(declination)))
        #radians   
        cosld = np.cos(lat*degree_to_rad)*np.cos(declination)  
        #radians
        aob = (sinld/cosld)  
        
        temp1 = np.asin(aob)
        #calculates daylength based on declination and latitude
        daylength = 12 * (1 + 2 * temp1/np.pi)   
        
        dsinB = 3600 * (daylength * sinld + 24 * cosld * np.sqrt(1 - aob * aob)/np.pi)
        dsinBE = 3600 * (daylength * (sinld + 0.4 * (sinld * sinld + cosld * cosld * 0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * np.sqrt(1 - aob * aob)/np.pi)
        #solar constant
        sc = 1370 * (1 + 0.033 * np.cos(2 * np.pi * day/365))  
        #Daily solar radiation
        dso = sc * dsinB  
        
        for hr in range(0,3):
            #calculates hour in which photosynthesis is applied
            hour1 = 12 + (daylength * 0.5 * (self._xgaus[hr]))  
            print(hour1)
            
            sinb_tmp = sinld + cosld * np.cos(2 * np.pi * (hour1 + 12)/24)
            print(sinb_tmp)
            #calculates sin of solar elevation, max functions prevents values less than 0
            sinB = max(np.Series([0, sinb_tmp]))  
            print(sinB)
            #dso can be replaced with values from FAO chart
            PAR1 = 0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE  
            print(PAR1)
            #convert to correct units
            PAR1 = PAR1 * (868/208.32)  
            print(PAR1)
            #output of function is vector of 3 values that represents time of day
            tmpvec.append(PAR1)  
            
        #returns a vector of light values in MicroEinsteins    
        return tmpvec  

    def _estimate_initial_plant_biomass(self):
        pidval=0
        for cell in range(self._grid.number_of_cells):
            cell_index=cell
            cell_plants=self._grid['cell']['vegetation__plant_species'][cell]
            for plant in cell_plants:
                if plant != 'null':
                    species=plant
                    newrow=[pidval,species,cell_index]
                    self.plants.loc[pidval]=newrow
                    pidval += 1
        
        biomass=np.empty([1,3])
        pidset=np.empty(1)
        #Set coefficients to variable names to pass to solver function
        for item in self.vegparams:
            growdict=self.vegparams[item]['growparams']
            aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
            astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
            coeffs=[aleaf,b1leaf,b2leaf,astem,b1stem,b2stem]
            pid=self.plants.pid.loc[self.plants['species']==item].to_numpy(int)
            pidset=np.concatenate((pidset,pid))
            total_biomass=np.random.rand(pid.shape[0])
            species_biomass=self._init_biomass_allocation(total_biomass, coeffs)
            biomass=np.concatenate((biomass,species_biomass), axis=0)
        biomass=pd.DataFrame(biomass, columns=['leaf_biomass','stem_biomass','root_biomass'])
        biomass['pid']=pidset
        self.plants=self.plants.set_index('pid').join(biomass.set_index('pid'))
        self.plants.reset_index(inplace=True)
    
    def _init_biomass_allocation(self, Tbio,coeffs):
        #Initialize arrays to calculate root, leaf and stem biomass to
        root=[]
        leaf=[]
        stem=[]
        root_all=np.zeros_like(Tbio)
        leaf_all=np.zeros_like(Tbio)
        stem_all=np.zeros_like(Tbio)
        
        #Loop through grid array
        for cell in Tbio:
            #Calculate initial guess on allocation based on total biomass
            zGuess = np.full(3,cell/3)
            #call to system of equations to solve for 
            z=fsolve(self._solverFuncs,zGuess,(coeffs,cell))
            #Transform from log10 values
            zvals=10**z
            root.append(zvals[0])
            leaf.append(zvals[1])
            stem.append(zvals[2])
        
        #Convert to numpy array
        root=np.array(root)
        leaf=np.array(leaf)
        stem=np.array(stem)

        biomass=np.vstack((leaf,stem, root))
        biomass=np.transpose(biomass)
       
        return biomass

    def _solverFuncs(self,x,coeffs,totval):
        r=x[0]
        l=x[1]
        s=x[2]
        T=totval
        F = np.empty([(3)])
        F[0]=10**r+10**l+10**s-T
        F[1]=coeffs[0]+coeffs[1]*r+coeffs[2]*r**2-l
        F[2]=coeffs[3]+coeffs[4]*r+coeffs[5]*r**2-s
        return F
        
    def _mortality(self):
        biomass=self._last_veg_biomass
        type=self._last_veg_type
        for item in self.vegparams:
            mortdict=self.vegparams[item]['mortparams']
            loop=0
            for fact in mortdict['mort_factor']:
                try:
                    pred=self._grid["cell"][fact][:].copy()
                    coeffs=mortdict['coeffs'][loop]
                    surv=1/(1+coeffs[0]*np.exp(-coeffs[1]*pred))
                    surv[np.isnan(surv)]=1.0
                    surv[surv<0]=0
                    rand_dum=np.random.randn(*self._last_veg_biomass.shape)
                    survd=surv^(1/(mortdict['duration'][loop]/self.dt))
                    risk=rand_dum<survd
                    biomass[type==item]=biomass[type==item]*risk[type==item]
                except KeyError:
                    msg=(f'No data available for mortality factor {loop}')
                    raise ValueError(msg)
        return biomass
