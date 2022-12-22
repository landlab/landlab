"""
Growth component of GenVeg - this is the main driver of vegetation growth and
is driven by a photosynthesis model. Vegetation growth depends on the availability
of carbohydrate produced by photosynthetically active plant parts. 
"""
from operator import invert
from tabnanny import check
from matplotlib.pyplot import grid
#from landlab.components import Radiation
from landlab.data_record import DataRecord
import numpy as np
import pandas as pd
from scipy.optimize import fsolve
import numpy.lib.recfunctions as nprf
from sympy import symbols, diff, lambdify, log
import warnings

from landlab.components.genveg.species import Species

class PlantGrowth(Species):
    """
    Add Intro Stuff here
    _name = "PlantGrowth"
    _unit_agnostic = False
    _cite_as = 
    @article{piercygv,
        author = {Piercy, C.D.; Swannack, T.M.; Carrillo, C.C.; Russ, E.R.; Charbonneau, B. M.]
    }
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
    """
    def __init__(
        self,
        grid,
        dt,
        current_day,
        start_date,
        plants= np.recarray((0,),
            dtype=[('species','U10'),('pid',int),('cell_index',int),('root_biomass',float),('leaf_biomass',float),('stem_biomass',float)]),
        **kwargs
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
                grid: RasterModelGrid
            A Landlab ModelGrid
        
        dt: int, required,
            time step interval

        plants: Numpy Record array of individual plants, optional
            with columns
            species: plant species code
            pid: plant ID
            root_biomass: plant live root biomass in g
            cell_index: index of cell location on grid
            leaf_biomass: plant live leaf biomass in g
            stem_biomass: plant live stem biomass in g

        **kwargs to send to init
            species_params: dict, optional,
                a nested dictionary of named vegetation parameters for the
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
                duration_params: dict, required,
                    dictionary of parameters defining the growing season,
                    gs_start: int, required,
                        growing season start day of year,
                        must be between 1-365
                    gs_end: int, required,
                        growing season end day of year,
                        must be between 1-365
                    senes: int, required,
                        start of senescence period after plant reaches peak biomass,
                        must be between gs_start and gs_end
                grow_params: dict, required,
                    dictionary of paramaters required to simulate plant growth
                    ptype: string, required,
                        photosythesis type, either 'C3', 'C4', or 'CAM'
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
        """
        #Initialize species object to get correct species parameter list
        species_params=kwargs['species_params']
        super().__init__(species_params)
        self.species_name=self.species_plant_factors['species']
        self._grid=grid
        (_,_latitude)=self._grid.xy_of_reference
        self._lat_rad = np.radians(_latitude)       

        
        self.dt=dt
        self.plants=plants
        self.current_day=current_day
        self.time_ind=1
        
        if len(self.plants)==0:
            self.plants=self._estimate_initial_plant_biomass()
            
        rel_time=(self.current_day-start_date).astype(float)
        #Create empty Datarecord to store plant data
        #Instantiate data record
        self.record_plants = DataRecord(
            self._grid, 
            time = [(self.current_day-start_date).astype(float)], 
            items = {
                "grid_element": np.repeat(['cell'], self.plants['pid'].size).reshape(self.plants['pid'].size,1),
                "element_id": np.reshape(self.plants['cell_index'],(self.plants['pid'].size,1))
            },
            data_vars = {
                'vegetation__species':(['item_id','time'],np.reshape(self.plants['species'],(self.plants['pid'].size,1))),
                'vegetation__root_biomass':(['item_id','time'],np.reshape(self.plants['root_biomass'],(self.plants['pid'].size,1))),
                'vegetation__leaf_biomass':(['item_id','time'],np.reshape(self.plants['leaf_biomass'],(self.plants['pid'].size,1))),
                'vegetation__stem_biomass':(['item_id','time'],np.reshape(self.plants['stem_biomass'],(self.plants['pid'].size,1)))
        },
            attrs={"vegetation__species": 'species name, string',
                "vegetation__root_biomass": 'g',
                "vegetation__leaf_biomass": 'g',
                "vegetation__stem_biomass": 'g',}
        )
        self.plants=nprf.append_fields(self.plants,('item_id'),data=self.record_plants.item_coordinates)
        #Set constants for PAR formula
        self._wgaus = [0.2778, 0.4444, 0.2778]
        self._xgaus = [0.1127, 0.5, 0.8873]
        self.delta_tot=[]

    def inc_growth(self):
        return self.delta_tot

    def species_plants(self):
        return self.plants
    
    def species_grow_params_out(self):
        return self.species_grow_params

    def _grow(self, current_day):
        #Calculate current day-of-year
        jday_td=current_day-np.datetime64(str(current_day.astype('datetime64[Y]'))+'-01-01')
        _current_jday=jday_td.astype(int)
        
        #set up shorthand aliases
        growdict=self.species_grow_params
        durationdict=self.species_duration_params
        _last_biomass=self.plants
        _total_biomass=_last_biomass['leaf_biomass']+_last_biomass['stem_biomass']+_last_biomass['root_biomass']
        _par=self._grid['cell']['radiation__par_tot'][_last_biomass['cell_index']]
        _temperature=self._grid['cell']['air__temperature_C'][_last_biomass['cell_index']]
        
        #define plant process flags
        _photosythesis_period=bool((_current_jday>=durationdict['growing_season_start'])&(_current_jday<=durationdict['growing_season_end']))
        _senescence_period=bool((_current_jday>=durationdict['senescence_start'])&(_current_jday<durationdict['growing_season_end']))
        _dormant_day=bool(_current_jday==durationdict['growing_season_end'])

        ##################################################
        #Growth and Respiration
        ##################################################
                    
        #maintenance respiration
        #if then statement stops respiration at end of growing season

        if _photosythesis_period:
            _declination = np.radians(23.45) * (np.cos(2 * np.pi / 365 * (172 - _current_jday)))
            _daylength = 24/np.pi*np.arccos(-(np.sin(_declination)*np.sin(self._lat_rad))/(np.cos(_declination)*np.cos(self._lat_rad)))

            #repiration coefficient for lvs, temp dependence from Teh 2006
            kmLVG = growdict['respiration_coefficient'][1] * pow(2,((_temperature - 25)/10))  
            #respiration coefficient for stems, temp depencence from Teh 2006 page 134
            kmSTG = growdict['respiration_coefficient'][2]* pow(2,((_temperature - 25)/10)) 
            #respiration coefficient for roots, temp dependence from Teh 2006 page 134
            kmRTG = growdict['respiration_coefficient'][0] * pow(2,((_temperature - 25)/10)) 
            #maintenance respiration per day from Teh 2006
            rmPrime = (kmLVG * _last_biomass['leaf_biomass']) + (kmSTG * _last_biomass['stem_biomass']) + (kmRTG * _last_biomass['root_biomass'])  
            #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration
            #THIS NEEDS TO BE UPDATED
            plantAge = _last_biomass['leaf_biomass']/_last_biomass['leaf_biomass']
            #plant age dependence from Teh 2006 page 145    
            respMaint = rmPrime * plantAge  
            
            #glucose requirement for growth save to temp private variable since we are not saving to class
            #from Teh 2006 page 148
            _glu_req = ((_last_biomass['leaf_biomass']/_total_biomass) * growdict['glucose_requirement'][1]) + \
                       ((_last_biomass['stem_biomass']/_total_biomass) * growdict['glucose_requirement'][2]) + \
                       ((_last_biomass['root_biomass']/_total_biomass) * growdict['glucose_requirement'][0]) 

            #dailyplantage[j] = plantAge    #how to handle

        #Enter photosynthesis loop  
            
            #radiation measured 3x daily, roughly correlates to morning, noon, afternoon
            #rad_est=self._PAR(_declination, _daylength, _current_jday)
            #change to read solar radiation from grid
            if _current_jday == durationdict['growing_season_end']:
                gphot = 0
            else:
                gphot=self.photosynthesize(_par, growdict, _last_biomass, _daylength)
            #if then statement ends glucose generation at end of growing season 


            #direct solve method for calculating change in biomass
            #coefficients rename
            aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
            astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
            #set up sympy equations
            rootsym=symbols('rootsym')              
            dleaf=diff(10**(aleaf+b1leaf*log(rootsym,10)+b2leaf*(log(rootsym,10))**2),rootsym)
            dstem=diff(10**(astem+b1stem*log(rootsym,10)+b2stem*(log(rootsym,10))**2),rootsym)
            #Generate numpy expressions and solve for rate change in leaf and stem biomass per unit mass of root
            fleaf=lambdify(rootsym,dleaf,'numpy')
            fstem=lambdify(rootsym,dstem,'numpy')
            delta_leaf_unit_root=fleaf(_last_biomass['root_biomass'])
            delta_stem_unit_root=fstem(_last_biomass['root_biomass'])
            #Calculate total new biomass and floor at zero
            delta_tot=(gphot-respMaint)/_glu_req
            delta_tot=np.where(delta_tot<0,0,delta_tot)
            #Calculate biomass that goes to root, stem, and leaf
            delta_root=delta_tot/(1+delta_leaf_unit_root+delta_stem_unit_root)
            delta_leaf=delta_leaf_unit_root*delta_root
            delta_stem=delta_stem_unit_root*delta_root
            #save for intermediate output
            self.delta_tot=delta_leaf
            #update biomass in plant array
            self.plants['root_biomass']=_last_biomass['root_biomass']+delta_root
            self.plants['leaf_biomass']=_last_biomass['leaf_biomass']+delta_leaf
            self.plants['stem_biomass']=_last_biomass['stem_biomass']+delta_stem
            
            if _senescence_period:
                self.plants=self.senesce(self.plants)
            if _dormant_day:
                self.plants=self.enter_dormancy(self.plants)

#May be able to eliminate this since using built in solar radiation component
    def _PAR(self, _declination, _daylength, _current_jday):
        lat=self._lat_rad
    
        tmpvec = []
        
        #Intermediate variables
        #radians
        sinld = ((np.sin(lat))*(np.sin(_declination)))
        #radians   
        cosld = np.cos(lat)*np.cos(_declination)  
        #radians
        aob = (sinld/cosld)  
        
        dsinB = 3600 * (_daylength * sinld + 24 * cosld * np.sqrt(1 - aob * aob)/np.pi)
        dsinBE = 3600 * (_daylength * (sinld + 0.4 * (sinld * sinld + cosld * cosld * 0.5)) + 12 * cosld * (2 + 3 * 0.4 * sinld) * np.sqrt(1 - aob * aob)/np.pi)
        #solar constant
        sc = 1370 * (1 + 0.033 * np.cos(2 * np.pi * _current_jday/365))  
        #Daily solar radiation
        dso = sc * dsinB  
        
        for hr in range(0,3):
            #calculates hour in which photosynthesis is applied
            hour1 = 12 + (_daylength * 0.5 * (self._xgaus[hr]))  
            sinb_tmp = sinld + cosld * np.cos(2 * np.pi * (hour1 + 12)/24)

            #calculates sin of solar elevation, max functions prevents values less than 0
            sinB = max([0, sinb_tmp])  

            #dso can be replaced with values from FAO chart
            PAR1 = 0.5 * dso * sinB * (1 + 0.4 * sinB) / dsinBE  

            #convert to correct units
            PAR1 = PAR1 * (868/208.32)  

            #output of function is vector of 3 values that represents time of day
            tmpvec.append(PAR1)  
            
        #returns a vector of light values in MicroEinsteins    
        return tmpvec  
    
    def _estimate_initial_plant_biomass(self):
        #Define datatypes for record array
        dtypes=[('species','U10'),('pid',int),('cell_index',int)]        
        #Create temporary variables to build initial plant list
        pidval=0
        buildlist=[]
        #Loop through grid
        for cell in range(self._grid.number_of_cells):
            cell_index=cell
            cell_plants=self._grid['cell']['vegetation__plant_species'][cell]
            #Loop through list of plants stored on grid
            for plant in cell_plants:
                if plant == self.species_plant_factors['species']:
                    buildlist.append((plant,pidval,cell_index))
                    pidval += 1
        #Create plant record array 
        if len(buildlist)==0:
            plant_array=np.recarray((0,), dtypes)
        else:                   
            plant_array=np.rec.array(buildlist,dtype=dtypes)
        #Create nan arrays and fill undefined variables
        fillnan=np.empty(len(plant_array))
        fillnan[:]=np.nan
        plant_array=nprf.append_fields(plant_array,('leaf_biomass','stem_biomass','root_biomass'),data=[fillnan,fillnan,fillnan])
        
        #Set coefficients to variable names to pass to solver function
        growdict=self.species_grow_params
        aleaf,b1leaf,b2leaf=growdict['root_to_leaf_coeffs']
        astem,b1stem,b2stem=growdict['root_to_stem_coeffs']
        coeffs=[aleaf,b1leaf,b2leaf,astem,b1stem,b2stem]
        init_mass=growdict['init_biomass']
        total_biomass=np.random.rand(plant_array.size)*(init_mass[1]-init_mass[0])+init_mass[0]
        #Call biomass allocation method
        root_bio,leaf_bio,stem_bio=self._init_biomass_allocation(total_biomass, coeffs)
        plant_array['root_biomass']=root_bio
        plant_array['leaf_biomass']=leaf_bio
        plant_array['stem_biomass']=stem_bio

        return plant_array
    
    def _init_biomass_allocation(self, total_biomass, solver_coeffs):
        #Initialize arrays to calculate root, leaf and stem biomass from total
        root=[]
        leaf=[]
        stem=[]
        
        #Loop through grid array
        for total_biomass_in_cell in total_biomass:
            solver_guess = np.full(3,np.log10(total_biomass_in_cell/3))
            
            part_biomass_log10=fsolve(self._solverFuncs,solver_guess,(solver_coeffs,total_biomass_in_cell))
            
            part_biomass=10**part_biomass_log10
            
            root.append(part_biomass[0])
            leaf.append(part_biomass[1])
            stem.append(part_biomass[2])
        
        #Convert to numpy array
        root=np.array(root)
        leaf=np.array(leaf)
        stem=np.array(stem)      
        return root, leaf, stem

    def _solverFuncs(self,solver_guess,solver_coeffs,total_biomass):
        root_part_log10=solver_guess[0]
        leaf_part_log10=solver_guess[1]
        stem_part_log10=solver_guess[2]
        plant_part_biomass_log10 = np.empty([(3)])

        plant_part_biomass_log10[0]=10**root_part_log10+10**leaf_part_log10+10**stem_part_log10-total_biomass
        plant_part_biomass_log10[1]=solver_coeffs[0]+solver_coeffs[1]*root_part_log10+solver_coeffs[2]*root_part_log10**2-leaf_part_log10
        plant_part_biomass_log10[2]=solver_coeffs[3]+solver_coeffs[4]*root_part_log10+solver_coeffs[5]*root_part_log10**2-stem_part_log10
        
        return plant_part_biomass_log10

    #Save plant array output Modify this in future to take user input and add additional parameters that can be saved out
    def save_plant_output(self, rel_time, save_params):
        prev_time=self.record_plants.latest_time
        #new_ids=np.where(np.isin(self.plants['pid'],self.record_plants.dataset['item_id'].values[:,prev_time], invert=True))
        #for i in new_ids:
        #    self.record_plants.add_item(time=[rel_time],
        #        new_item={
        #            'grid_element'=['cell'],
        #            'element_id'=self.plants['cell_index'][i],
        #        },

        #    )
        self.record_plants.add_record(
            time=np.array([rel_time])
        )
        self.record_plants.ffill_grid_element_and_id()

        self.record_plants.dataset['vegetation__species'].values[:,self.time_ind]=self.plants['species']
        self.record_plants.dataset['vegetation__root_biomass'].values[:,self.time_ind]=self.plants['root_biomass']
        self.record_plants.dataset['vegetation__leaf_biomass'].values[:,self.time_ind]=self.plants['leaf_biomass']
        self.record_plants.dataset['vegetation__stem_biomass'].values[:,self.time_ind]=self.plants['stem_biomass']
        #    item_id= np.reshape(self.plants['item_id'],(self.plants['item_id'].size,1)),
        #    new_record={
        #        'vegetation__species':(['item_id','time'],self.plants['species']),
        #        #'vegetation__root_biomass':(['item_id','time'],np.reshape(self.plants['root_biomass'],(self.plants['cell_index'].size,1))),
        #        #'vegetation__leaf_biomass':(['item_id','time'],np.reshape(self.plants['leaf_biomass'],(self.plants['cell_index'].size,1))),
        #        #'vegetation__stem_biomass':(['item_id','time'],np.reshape(self.plants['stem_biomass'],(self.plants['cell_index'].size,1)))
        #    }
        #)
        self.time_ind +=1
        
    
    #Initial recode of simple mortality function        
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
