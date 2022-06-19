"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab import Component
import numpy as np
from scipy.optimize import fsolve

class PlantGrowth(Component):
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
        "plant__total_biomass": {
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
        current_day=0 
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

        current_day: int, required,
            day of simulation
        """
        
        super().__init__(grid)

        #Check to see if grid contains required vegetation fields
        try:
            self._last_veg_biomass = self._grid["cell"]["vegetation__live_biomass"][:].copy()
        except KeyError:
            msg =("GenVeg requires initial vegetation biomass as an at-cell field.")
            raise ValueError(msg)
        try:
            self._last_veg_n_plant = self._grid["cell"]["vegetation__n_plants"][:].copy()
        except KeyError:
            msg = ("GenVeg requires initial number of plants as an at-cell field.")
            raise ValueError(msg)        
        #We need to decide how plant functional types, species, etc. are handled. I do not want the model
        #to carry all the data forward so I'd prefer it to be embedded in the plant dictionaries. So we would 
        #have a string field on the grid that has the dictionary key name labels. Thoughts?
        try:
            self._last_veg_species=self._grid["cell"]["vegetation__plant_species"][:].copy()
        except KeyError:
            msg =("GenVeg requires initial vegetation species distribution as an at-cell field.")
            raise ValueError(msg)

        #Check to see if grid contains required environmental fields
        try:
            self._air_temp = self._grid['cell']['air__temperature_C'][:].copy()
        except KeyError:
            msg=('GenVeg requires air temperature in Celcius for each time step.')
            raise ValueError(msg)
        
        try:
            self._radiation = self._grid["cell"]["radiation__net_flux"][:].copy()
        except KeyError:
            msg=('GenVeg requires incoming radiation flux for each time step')
            raise ValueError(msg)
        self._last_twg = self._last_veg_biomass/self._last_veg_n_plant
        
        self.vegparams=vegparams
        self.dt=dt
        self.current_day=current_day

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

            #Set coefficients to variable names to pass to solver function
            aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
            astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
            coeffs=[aleaf,b1leaf,b2leaf,astem,b1stem,b2stem]
            self._last_veg_root_biomass, self._last_veg_leaf_biomass, self._last_veg_stem_biomass=self._init_biomass_allocation(item, coeffs)


        
    def run_one_step(self):
        #Add photosythesis code here
        for species in self.vegparams:
            factordict=self.vegparams[species]['plant_factors']
            growdict=self.vegparams[species]['growparams']
            if self.current_day == growdict['gs_start']:
                #allocate biomass - I don't think it how we want to do it in the long term but for now.
                self._last_veg_root_biomass, self._last_veg_leaf_biomass, self._last_veg_stem_biomass=self._init_biomass_allocation(item, coeffs )
    
            ##################################################
            #Calculate Vegetation Structure Metrics Each Day
            ##################################################
            
            totLeafWeight = twlvd + twlvg #twlvd = total (t) weight(w) leaves(lv), d @ dead and g is living or green
            totStemWeight = twstd + twstg #total weight of stems where d at end is dead and g is living or green
            totRootWeight = twrtd + twrtg #total weight of roots where d at end is dead and g is living or green

            ##################################################
            #Growth and Respiration
            ##################################################
                      
            #maintenance respiration
            #if then statement stops respiration at end of growing season
            if self.current_day == growdict['growing_season_end']:
                respMaint = 0

            if self.current_day < growdict['growing_season_end'] and self.current_day >= growdict['growing_season_start']:
                kmLVG = kmLVG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10))  #repiration coefficient for lvs, temp dependence from Teh 2006
    
                kmSTG = kmSTG_prime* pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for stems, temp depencence from Teh 2006 page 134
  
                kmRTG = kmRTG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for roots, temp dependence from Teh 2006 page 134
     
                rmPrime = (kmLVG * twlvg) + (kmSTG * twstg) + (kmRTG * twrtg)  #maintenance respiration per day from Teh 2006
     
                plantAge = twlvg/totLeafWeight #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration

                if np.isnan(plantAge):
                    plantAge = 0
                    
                respMaint = rmPrime * plantAge  #plant age dependence from Teh 2006 page 145
                
                #glucose requirement for growth save to temp private variable since we are not saving to class
                _glu_req = ((self._last_veg_leaf_biomass/self._last_veg_biomass) * growdict['glucose_requirement'][1]) + \
                           ((self._last_veg_stem_biomass/self._last_veg_biomass) * growdict['glucose_requirement'][2]) + \
                           ((self._last_veg_root_biomass/self._last_veg_biomass) * growdict['glucose_requirement'][0]) #from Teh 2006 page 148

                
                #writes results for daily respiration, plant age, and maintenance respiration
                #I don't know if we need to write these anywhere but we can write to a spearate test grid for now (NOT the main grid)
                #dailyrespiration.iloc[j] = rmPrime  #how to handle

                #dailyplantage[j] = plantAge    #how to handle

                #dailyrespMaint[j] = respMaint   #how to handle
                
                    
                    
            #Enter photosynthesis loop  
                #dailyphoto = []
                #determine how we want to read in light. Ideally this would be dynamic
                dtga=0
                for hr in range(0,3):  #radiation measured 3x daily, roughly correlates to morning, noon, afternoon
                    parMicroE = (Light.iloc[j,hr]) * (868/208.32) #convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on
                
                    intSolarRad = parMicroE*np.exp(-(growdict['k_light_extinct'])*self._last_veg_leaf_biomass)  #from Charisma instructions: tells how much of the light a plant is going to get as PAR in microeinsteins based on how many leaves are on the plant

                    intLightpH = intSolarRad/(intSolarRad+Hi) #amoung of light absorbed, per half saturaion constants from Charisma eq. 3. the monod or michaelis/menten function is adequate for describing the photosynthetic response to light

                    photosynthesis = (growdict['p_max']) * intLightpH #pMax is the maximum rate of photosynthesis, species specific

                    #fgross.iloc[hr] = photosynthesis #calculates gross assimilation of fgross(like APT) via photosynthesis at specific hour calculate growth per day at three times per day, morning, middday, and evenning and this amount is weighted based on how much light is hitting hte plant based on the latitude of your study site
                
                    dtgastep = photosynthesis*wgaus[hr] #weights fgross for specific time of day
                    dtga+=dtgastep
                        
                    
                    
                    
                dtgaCollapsed = dtga*self._last_veg_leaf_biomass  #calculates total biomass gained across plant (twlvg is amount of leaver/green matter): you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
                
                assimilatedCH2O = dtgaCollapsed*Light['daylength'][self.current_day] #total biomass for day length

                gphot = assimilatedCH2O*(30/44) #converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
                
                #I don't think we need to write all of this out except for testing
                dailyphoto.iloc[j] = gphot  ##how do we want to handle something like this?

                
                
            #if then statement ends glucose generation at end of growing season 
            if self.current_day == growdict['gs_end']:
                gphot = 0

            #direct solve method for calculating change in biomass
            #coefficients rename
            aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
            astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
            #Calculate the change in leaf and stem mass per unit change in root mass given the current size of the root biomass
            delta_leaf_unit_root=0.43429448190325176*b1leaf/self._last_root_biomass + \
                        0.37722339402322774*b2leaf*np.log(self._last_root_biomass)/self._last_root_biomass
            delta_stem_unit_root=0.43429448190325176*b1stem/self._last_root_biomass + \
                        0.37722339402322774*b2stem*np.log(self._last_root_biomass)/self._last_root_biomass
            #Calculate the total change in biomass this timestep
            delta_tot=(gphot-respMaint)/_glu_req
            #Calculate the change in root biomass this timestep
            #Solve for delta root: delta_tot=delta_root+delta_leaf_unit_root*delta_root+delta_stem_unit_root*delta_root
            delta_root=delta_tot/(1+delta_leaf_unit_root+delta_stem_unit_root)
            delta_leaf=delta_leaf_unit_root*delta_root
            delta_stem=delta_stem_unit_root*delta_root
            self._current_vegetation_biomass=self._last_vegetation_biomass+delta_tot
            self._current_root_biomass=self._last_root_biomass+delta_root
            self._current_leaf_biomass=self._last_leaf_biomass+delta_leaf
            self._current_stem_biomass=self._last_stem_biomass+delta_stem


    def _PAR(self, day, lat):
    
        degree_to_rad = 0.017453292  #required to convert degrees to radians
        #rad_to_degree = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))
    
    
    
        tmpvec = []
        declination = (-math.asin((math.sin(23.45*degree_to_rad))*(math.cos(2*math.pi*(day+10)/365))))
        
        #Intermediate variables
        sinld = ((math.sin(lat*degree_to_rad))*(math.sin(declination)))   #radians
        cosld = math.cos(lat*degree_to_rad)*math.cos(declination)  #radians
        
        aob = (sinld/cosld)  #radians
        
        temp1 = math.asin(aob)
        
        daylength = 12 * (1 + 2 * temp1/math.pi)   #calculates daylength based on declination and latitude
        
        dsinB = 3600 * (daylength * sinld + 24 * cosld * math.sqrt(1 - aob * aob)/math.pi)
        dsinBE = 3600 * (daylength * (sinld + 0.4 * (sinld * sinld + cosld * cosld * 0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * math.sqrt(1 - aob * aob)/math.pi)
        
        sc = 1370 * (1 + 0.033 * math.cos(2 * math.pi * day/365))  #solar constant
        
        dso = sc * dsinB  #Daily solar radiation
        
        for hr in range(0,3):
            hour1 = 12 + (daylength * 0.5 * (xgauss[hr]))  #calculates hour in which photosynthesis is applied
            print(hour1)
            
            sinb_tmp = sinld + cosld * math.cos(2 * math.pi * (hour1 + 12)/24)
            print(sinb_tmp)
            
            sinB = max(pd.Series([0, sinb_tmp]))  #calculates sin of solar elevation, max functions prevents values less than 0
            print(sinB)
            
            PAR1 = 0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE  #dso can be replaced with values from FAO chart
            print(PAR1)
            
            PAR1 = PAR1 * (868/208.32)  #convert to correct units
            print(PAR1)
            
            tmpvec.append(PAR1)  #output of function is vector of 3 values that represents time of day
            
            
        return tmpvec  #returns a vector of light values in MicroEinsteins

    def _init_biomass_allocation(self, species, coeffs):
        #Initialize arrays to calculate root, leaf and stem biomass to
        Tbio=self._last_veg_biomass
        root=np.zeros_like(Tbio)
        leaf=np.zeros_like(Tbio)
        stem=np.zeros_like(Tbio)
        root_all=np.zeros_like(Tbio)
        leaf_all=np.zeros_like(Tbio)
        stem_all=np.zeros_like(Tbio)
        rowind=0

        #Loop through grid array
        for cell in Tbio:
            #Calculate initial guess on allocation based on total biomass
            zGuess = cell/3
            #call to system of equations to solve for 
            z=fsolve(self._solverFuncs,zGuess,(coeffs,cell))
            #Transform from log10 values
            zvals=10**z
            root[cell]=zvals[0]
            leaf[cell]=zvals[1]
            stem[cell]=zvals[2]
            colind+=1
            #Only save values to cells where species matches
            root_all[self._last_veg_species==species]=root[self._last_veg_species==species]
            leaf_all[self._last_veg_species==species]=leaf[self._last_veg_species==species]
            stem_all[self._last_veg_species==species]=stem[self._last_veg_species==species]
        
        return root_all, leaf_all, stem_all

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
