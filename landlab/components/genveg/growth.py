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
        rel_time,
        _current_jday,
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
        self.time_ind=1
        _in_growing_season,_,_=self.set_event_flags(_current_jday)
        
        if len(self.plants)==0:
            self.plants=self._init_plants_from_grid(_in_growing_season)
            
        #Create empty Datarecord to store plant data
        #Instantiate data record
        self.record_plants = DataRecord(
            self._grid, 
            time = [rel_time], 
            items = {
                "grid_element": np.repeat(['cell'], self.plants['pid'].size).reshape(self.plants['pid'].size,1),
                "element_id": np.reshape(self.plants['cell_index'],(self.plants['pid'].size,1))
            },
            data_vars = {
                'vegetation__species':(['item_id','time'],np.reshape(self.plants['species'],(self.plants['pid'].size,1))),
                'vegetation__root_biomass':(['item_id','time'],np.reshape(self.plants['root_biomass'],(self.plants['pid'].size,1))),
                'vegetation__leaf_biomass':(['item_id','time'],np.reshape(self.plants['leaf_biomass'],(self.plants['pid'].size,1))),
                'vegetation__stem_biomass':(['item_id','time'],np.reshape(self.plants['stem_biomass'],(self.plants['pid'].size,1))),
                'vegetation__storage_biomass':(['item_id','time'],np.reshape(self.plants['storage_biomass'],(self.plants['pid'].size,1))),
                'vegetation__plant_age':(['item_id','time'],np.reshape(self.plants['plant_age'],(self.plants['pid'].size,1)))
            },
            attrs={"vegetation__species": 'species name, string',
                "vegetation__root_biomass": 'g',
                "vegetation__leaf_biomass": 'g',
                "vegetation__stem_biomass": 'g',
                "vegetation__storage_biomass": 'g',
                "vegetation__plant_age": 'days',
            }
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

    def _grow(self, _current_jday):
        #Calculate current day-of-year

        
        #set up shorthand aliases
        growdict=self.species_grow_params
        durationdict=self.species_duration_params
        _last_biomass=self.plants
        _total_biomass=_last_biomass['leaf_biomass']+_last_biomass['stem_biomass']+_last_biomass['root_biomass']+_last_biomass['storage_biomass']
        _par=self._grid['cell']['radiation__par_tot'][_last_biomass['cell_index']]
        _temperature=self._grid['cell']['air__temperature_C'][_last_biomass['cell_index']]
        
        _in_growing_season, _in_senescence_period, _is_dormant_day=self.set_event_flags(_current_jday)


        ##################################################
        #Growth and Respiration
        ##################################################
                    
        #maintenance respiration
        #if then statement stops respiration at end of growing season
        #repiration coefficient temp dependence from Teh 2006
        kmLVG = growdict['respiration_coefficient'][1] * pow(2,((_temperature - 25)/10))  
        kmSTG = growdict['respiration_coefficient'][2]* pow(2,((_temperature - 25)/10)) 
        kmRTG = growdict['respiration_coefficient'][0] * pow(2,((_temperature - 25)/10)) 
        #maintenance respiration per day from Teh 2006
        rmPrime = (kmLVG * _last_biomass['leaf_biomass']) + (kmSTG * _last_biomass['stem_biomass']) + (kmRTG * _last_biomass['root_biomass'])  
        #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration
        #THIS NEEDS TO BE UPDATED
        plantAge = _last_biomass['leaf_biomass']/_last_biomass['leaf_biomass']
        #plant age dependence from Teh 2006 page 145    
        respMaint = rmPrime * plantAge  
        #from Teh 2006 page 148
        _glu_req = ((_last_biomass['leaf_biomass']/_total_biomass) * growdict['glucose_requirement'][1]) + \
                    ((_last_biomass['stem_biomass']/_total_biomass) * growdict['glucose_requirement'][2]) + \
                    ((_last_biomass['root_biomass']/_total_biomass) * growdict['glucose_requirement'][0]) 
        
        
        if _in_growing_season:
            _declination = np.radians(23.45) * (np.cos(2 * np.pi / 365 * (172 - _current_jday)))
            _daylength = 24/np.pi*np.arccos(-(np.sin(_declination)*np.sin(self._lat_rad))/(np.cos(_declination)*np.cos(self._lat_rad)))

            gphot=self.photosynthesize(_par, growdict, _last_biomass, _daylength)

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
            
            if _in_senescence_period:
                self.plants=self.senesce(self.plants)
            if _is_dormant_day:
                self.plants=self.enter_dormancy(self.plants)
        else:
            delta_tot=-respMaint/_glu_req
            self.plants['root_biomass']=(_last_biomass['root_biomass']/_total_biomass)*delta_tot+_last_biomass['root_biomass']
            self.plants['leaf_biomass']=(_last_biomass['leaf_biomass']/_total_biomass)*delta_tot+_last_biomass['leaf_biomass']
            self.plants['stem_biomass']=(_last_biomass['stem_biomass']/_total_biomass)*delta_tot+_last_biomass['stem_biomass']
    
    def _init_plants_from_grid(self,in_growing_season):
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
        plant_array=nprf.append_fields(plant_array,('leaf_biomass','stem_biomass','root_biomass','storage_biomass','plant_age'),data=[fillnan,fillnan,fillnan,fillnan,fillnan])
        plant_array=self.set_initial_biomass(plant_array,in_growing_season)
        return plant_array
        
   
    def set_event_flags(self,_current_jday):
        durationdict=self.species_duration_params
        _in_growing_season=bool((_current_jday>=durationdict['growing_season_start'])&(_current_jday<=durationdict['growing_season_end']))
        _in_senescence_period=bool((_current_jday>=durationdict['senescence_start'])&(_current_jday<durationdict['growing_season_end']))
        _is_dormant_day=bool(_current_jday==durationdict['growing_season_end'])
        return _in_growing_season,_in_senescence_period,_is_dormant_day

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
