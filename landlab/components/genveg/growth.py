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
        plant_ids,
        growthparams={
            'ftype': {'Corn':['1']},
            'ph_type': {'Corn':['C4']},
            'gs_start': {'Corn':[65]},
            'gs_end': {'Corn':[270]},
            'gs_length': {'Corn':[205]},
            'senes': {'Corn':[209]},
            'res_co': {'Corn':[0.15]},
            'glu_req': {'Corn':[55]},
            'le_k': {'Corn':[30]},
            'hi': {'Corn':[0.1]},
            'p_max': {'Corn':[2.5]}
        },
        sizeparams={},
        dispparams={},
        storparams={},
        colparams={},
        mortparams={},
        dt=1,
        current_day=0 
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
        grid: RasterModelGrid
            A Landlab ModelGrid

        plant_ids: dict, required,
            names of plant communities to be simulated

        growthparams: dict, required,
            dictionary of named growth parameters with below keys.
            A separate dict is stored for each type of plant with 
            the key as the unique plant type ID as a str.

            ftype: dict, int, optional,
                functional type of plant class
            ph_type: dict, string, required,
                photosynthesis type (option are: CAM, C3 or C4)
            gs_start: dict, int, required,
                growing season start day of year
            gs_end: dict, int, required,
                growing season end day of year
            gs_length: dict, int, required,
                growing season length in days
            senes: dict, int, required,
                start of senescence period after plant reaches peak biomass
            res_co: dict, float, required,
                respiration coefficient
            glu_req: dict, float, required,
                glucose requirement
            le_k: dict, float, required,
                light extinction coefficient
            hi: dict, float, required,
                something
            p_max: dict, float, required,
                maximum photosyntehtic output
            allocate: dict, float <1, required
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
        
        self.plant_ids=plant_ids
        self.num_ids=len(self.plant_ids)
        self.growth_params=growthparams
        self.dt=dt
        self.current_day=current_day

        #Verify that all named plant classes have defined growth parameters
        self._check_params(self.growth_params)

        #Calculate calculated parameters for speed and error check growth parameter values
        new_item={}
        for item in self.plant_ids:
            end=self.growth_params['gs_end'][item][0]
            begin=self.growth_params['gs_start'][item][0]
            length=end-begin+1
            new_item.update({item:[length]})
            allo=self.growth_params['allocate'][item]
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
        
        self.growth_params['gs_length']=new_item


    def _check_params(self,check):
        """Check to make sure the number of parameters is consistent
        with the number of vegetation classes in the parameter dictionary
        """
        for items in check:
            num_entries=len(check.get(items).keys())
            if num_entries != self.num_ids:
                msg = ('The number of values in ' + items + ' does not'
                'match the number of plant types')
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
