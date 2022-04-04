"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab import Component

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
                'growthparams': {
                    'glucose_requirement': [1.444, 1.513, 1.463], 
                    'growing_season_end': 290, 
                    'growing_season_start': 91, 
                    'k_light_extinct':0.02, 
                    'light_half_sat': 9, 
                    'p_max': 0.055, 
                    'plant_part_allocation': [0.2, 0.3, 0.5], 
                    'plant_part_min': [0.01, 0.1, 0.5], 
                    'ptype': 'C3', 
                    'respiration_coefficient': [0.015, 0.015, 0.03], 
                    'senescence_start': 228
                }, 
                'mortparams': {
                    's1_days': 365, 
                    's1_name': 'Mortality factor', 
                    's1_pred': [1, 2, 3, 4], 
                    's1_rate': [0, 0.1, 0.9, 1], 
                    's1_weight': [1000, 1, 1, 1000]
                },
                'plant_factors': {
                    'annual_perennial': 'C4', 
                    'growth_form': 1, 
                    'species': 'Corn'
                }, 
                'sizeparams': {
                    'max_height_stem': 2.5, 
                    'max_mass_stem': 72, 
                    'max_n_stems': 3, 
                    'max_plant_density': 1, 
                    'total_cs_area_stems': 0.231
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
                    USDA plant growth habit
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
                allocate: list, float <1, required
                    3-element list of proportion of biomass in roots, stems, and leaves
                    Values must be between 0-1 and all sum to 1

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
            msg =("GenVeg requires initial vegetation biomass coverage as an at-cell field.")
            raise ValueError(msg)
        
        try:
            self._last_veg_ftype = self._grid["cell"]["vegetation__plant_functional_type"][:].copy()
        except KeyError:
            msg =("GenVeg requires initial vegetation functional type distribution as an at-cell field.")
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
        
        self.dt=dt
        self.current_day=current_day

        #Calculate calculated parameters for speed and error check growth parameter values
        for item in self.vegparams:
            growdict=self.vegparams[item]['growparams']
            #Check validity of growing season beginning and ending days, senescence beginning, and calculate growing season length
            end=growdict['growing_season_end']
            begin=growdict['growing_season_start']
            senes=growdict['senescence_start']
            if begin > 0 & begin  < 366 & end > 0 & end < 366 & end > begin & senes > begin & senes < end:
                length=end-begin+1
            elif begin < 1 | begin > 365:
                msg='Growing season beginning must be between 1-365'
                raise ValueError(msg)
            elif end < 1 | end > 365:
                msg='Growing season end must be between 1-365'
                raise ValueError(msg)
            elif end < begin:
                msg='Growing season beginning must be before growing season end'
                raise ValueError(msg)
            elif senes < begin | senes > end:
                msg='Start of senescence must be within the growing season'
                raise ValueError(msg)
            else:
                msg='Growing season beginning and end must be integer values between 1-365'
                raise ValueError(msg)
            self.vegparams[item]['growparams']['growing_season_length']=length

            #Check plant part allocation to determine if all parts are between 0 and 1 and all sum to 1
            allo=growdict['plant_part_allocation']
            allosum=0
            for vals in allo:
                if vals <= 1 and vals >=0:
                    allosum += vals
                else:
                    msg='Biomass allocation proportions for plant ID ' + item +' must be between 0 and 1'
                    raise ValueError(msg)
            if allosum != 1.0:
                msg='Biomass allocation for plant ID ' +item +' must sum to 1'
                raise ValueError(msg)
        
    def run_one_step(self):
        #Add photosythesis code here
        if self.current_day == self.gs_start:
            twlvg = 0.25
            twstg = 0.1
            twrtg = 0.15
            
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
        if self.current_day < self.gs_end and self.current_day >= self.gs_start:

            if self.current_day < self.gs_end:
                kmLVG = kmLVG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10))  #repiration coefficient for lvs, temp dependence from Teh 2006
  
                kmSTG = kmSTG_prime* pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for stems, temp depencence from Teh 2006 page 134

                
                kmRTG = kmRTG_prime * pow(2,((Light.iloc[j]['meantemp'] - 25)/10)) #respiration coefficient for roots, temp dependence from Teh 2006 page 134

                
                rmPrime = (kmLVG * twlvg) + (kmSTG * twstg) + (kmRTG * twrtg)  #maintenance respiration per day from Teh 2006

                
                plantAge = twlvg/totLeafWeight #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration

    
                if math.isnan(plantAge):
                    plantAge = 0
                
                respMaint = rmPrime * plantAge  #plant age dependence from Teh 2006 page 145

        
            #if then statement stops respiration at end of growing season
            if self.current_day >= self.gs_end:

                respMaint = 0
            
            #glucose requirement for growth
            self.glu_req = (FracDM_LVG * glucoseReqLVG) + (FracDM_STG * glucoseReqSTG) + (FracDM_RTG * glucoseReqRTG) #from Teh 2006 page 148

            
            #writes results for daily respiration, plant age, and maintenance respiration
            dailyrespiration.iloc[j] = rmPrime  #how to handle

            dailyplantage[j] = plantAge    #how to handle

            dailyrespMaint[j] = respMaint   #how to handle
            
                
                
            #Enter photosynthesis loop  
            
            if self.current_day < self.gs_end:
                #dailyphoto = []
                for hr in range(0,3):  #radiation measured 3x daily, roughly correlates to morning, noon, afternoon
                    parMicroE = (Light.iloc[j,hr]) * (868/208.32) #convert to correct units which is microeinsteins which is the unit measure of light and what this model is based on
                   
                    intSolarRad = parMicroE*math.exp(-(PlantParameters['k'])*twlvg)  #from Charisma instructions: tells how much of the light a plant is going to get as PAR in microeinsteins based on how many leaves are on the plant

                    intLightpH = intSolarRad/(intSolarRad+Hi) #amoung of light absorbed, per half saturaion constants from Charisma eq. 3. the monod or michaelis/menten function is adequate for describing the photosynthetic response to light

                    photosynthesis = (PlantParameters['pMax']) * intLightpH #pMax is the maximum rate of photosynthesis, species specific

                    fgross.iloc[hr] = photosynthesis #calculates gross assimilation of fgross(like APT) via photosynthesis at specific hour calculate growth per day at three times per day, morning, middday, and evenning and this amount is weighted based on how much light is hitting hte plant based on the latitude of your study site
                   
                    dtga.iloc[hr] = fgross[hr]*wgaus[hr] #weights fgross for specific time of day
                    
                
                
                
            dtgaCollapsed = sum(dtga)*twlvg  #calculates total biomass gained across plant (twlvg is amount of leaver/green matter): you feed the model total biomass and then from that we determine how much leaf mass there is and so then basically an average of how much that average leaf will produce multiplied by the number of leaves, this is assuming that all leaves are mature
            
            assimilatedCH2O = dtgaCollapsed*Light['daylength'][self.current_day] #total biomass for day length

            gphot = assimilatedCH2O*(30/44) #converts carbohydrates to glucose where photosynthesis unit is glucose and then we later convert that glucose to biomass in another section
            
            dailyphoto.iloc[j] = gphot  ##how do we want to handle something like this?
            
            
        #if then statement ends glucose generation at end of growing season 
        if self.current_day >= self.gs_end:
            gphot = 0



    def PAR(day, lat):
    
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
