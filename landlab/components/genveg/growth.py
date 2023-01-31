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
           species_params={
            'col_params': {}, 
            'disp_params': {},
            'duration_params': {
                'growing_season_start': 91,
                'growing_season_end': 290,
                'senescence_start': 228,
            },
            'grow_params': {
                'respiration_coefficient': [0.015,0.015,0.03],
                'glucose_requirement': [1.444,1.513,1.463],
                'k_light_extinct':0.02,
                'light_half_sat':9,
                'p_max':0.055,
                'root_to_leaf_coeffs': [0.031,0.951,0],
                'root_to_stem_coeffs': [-0.107, 1.098, 0.0216],
                'plant_part_min':[0.01,0.1,0.5]
            }, 
            'mort_params': {
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
                'p_type':'C3'
            }, 
            'size_params': {
                'max_height_stem': 2.5, 
                'max_mass_stem': 72, 
                'max_n_stems': 3, 
                'max_plant_density': 1
            },
            'stor_params': {
                'r_wint_die': 0.25, 
                'r_wint_stor': 0.25
                },
        },
        **kwargs
    ):
        """Instantiate PlantGrowth
        Parameters
        ----------
        grid: RasterModelGrid
            A Landlab ModelGrid
        
        dt: NumPy time delta, required,
            time step interval

        rel_time: int, required,
            number of time steps elapsed

        _current_jday: int, required
            day of the year assuming Jan 1 is 1 

        **kwargs to send to init
            plants: Numpy structured array of individual plants, optional
                with columns
                species: string, plant species names
                pid: int, plant ID
                cell_index: int, index of cell location on grid
                root_biomass: float, title='root', plant live root biomass in g
                leaf_biomass: float, title='stem', plant live leaf biomass in g
                stem_biomass: float, title='stem', plant live stem biomass in g
                storage_biomass: float, title='storage', plant live stem biomass in g
                repro_biomass: float, title='reproductive', plant live reproductive biomass in g
                plant_age: int, plant age in days
        
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
                        graminoid, forb/herb, shrub, tree, vine
                    monocot_dicot: string, required,
                        should be monocot or dicot
                    angio_gymno: string, required,
                        should be angiosperm or gymnosperm
                    annual_perennial: string, required,
                        plant growth duration, annual (1 year) or
                        perennial (multiple years)
                    p_type: string, required,
                        photosythesis type, either 'C3', 'C4', or 'CAM'
                    leaf_retention: string, required,
                        evergreen or deciduous (annuals are deciduous)
                duration_params: dict, required,
                    dictionary of parameters defining the growing season,
                    growing_season_start: int, required,
                        growing season start day of year,
                        must be between 1-365
                    growing_season_end: int, required,
                        growing season end day of year,
                        must be between 1-365
                    senesecence_start: int, required,
                        start of senescence period after plant reaches peak biomass,
                        must be between gs_start and gs_end
                grow_params: dict, required,
                    dictionary of paramaters required to simulate plant growth

                    respiration_coefficient: dict, required,
                        respiration coefficient with keys
                            'root': float
                            'leaf': float
                            'stem': float
                            'storage': float
                            'reproductive': float
                    glucose_requirements: dict, required,
                        glucose requirement
                    le_k: float, required,
                        light extinction coefficient
                    hi: float, required,
                        something
                    p_max: float, required,
                        maximum photosyntehtic output
        """
        #Initialize species object to get correct species parameter list
     

        super().__init__(species_params)
        self.species_name=self.species_plant_factors['species']
        self.populate_biomass_allocation_array()

        self._grid=grid
        (_,_latitude)=self._grid.xy_of_reference
        self._lat_rad = np.radians(_latitude)      
        
        self.dt=dt
        self.time_ind=1

        event_flags=self.set_event_flags(_current_jday)
        _in_growing_season=event_flags.pop('_in_growing_season')
        self.plants=kwargs.get(
            'plant_array',
            self._init_plants_from_grid(_in_growing_season)
        )
        self.call=[]
            
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
                'vegetation__repro_biomass':(['item_id','time'],np.reshape(self.plants['repro_biomass'],(self.plants['pid'].size,1))),
                'vegetation__plant_age':(['item_id','time'],np.reshape(self.plants['plant_age'],(self.plants['pid'].size,1)))
            },
            attrs={"vegetation__species": 'species name, string',
                "vegetation__root_biomass": 'g',
                "vegetation__leaf_biomass": 'g',
                "vegetation__stem_biomass": 'g',
                "vegetation__storage_biomass": 'g',
                "vegetation__repro_biomass":'g',
                "vegetation__plant_age": 'days',
            }
        )
        self.plants['item_id']=self.record_plants.item_coordinates
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
        print(_current_jday)
        print(self.species_plant_factors['species'])
        #set up shorthand aliases
        growdict=self.species_grow_params
        _last_biomass=self.plants
        _new_biomass=_last_biomass
        _total_biomass=self.sum_plant_parts(_last_biomass,parts='total')
        filter=np.where(_total_biomass>0)

        #calculate variables needed to run plant processes
        _par=self._grid['cell']['radiation__par_tot'][_last_biomass['cell_index']]
        _temperature=self._grid['cell']['air__temperature_C'][_last_biomass['cell_index']]
        _declination = np.radians(23.45) * (np.cos(2 * np.pi / 365 * (172 - _current_jday)))
        _daylength = 24/np.pi*np.arccos(-(np.sin(_declination)*np.sin(self._lat_rad))/(np.cos(_declination)*np.cos(self._lat_rad)))

        #Decide what processes happen today
        event_flags=self.set_event_flags(_current_jday)
        processes={
            '_in_growing_season':self.photosynthesize,
            '_in_senescence_period':self.senesce,
            '_in_reproductive_period':self.disperse,
            '_is_emergence_day':self.emerge,
            '_is_dormant_day':self.enter_dormancy
        }
        
        #Run respiration and photosynthesis to determine carb gained
        _glu_req=np.zeros_like(_total_biomass)
        for part in self.all_parts:
            _glu_req[filter]+=(_last_biomass[part][filter]/_total_biomass[filter])*growdict['glucose_requirement'][part]
            delta_biomass_respire=self.respire(_temperature,_last_biomass,_glu_req) 

        if event_flags['_in_growing_season']:
            delta_photo=processes['_in_growing_season'](_par, _last_biomass, _glu_req, _daylength)
            delta_tot=delta_biomass_respire+delta_photo
            #we will need to add a turnover rate estiamte. Not sure where to include it though.
            #delta_tot=delta_tot-self.species_grow_params['biomass_turnover_rate']*self.dt
            _new_biomass=self.allocate_biomass_dynamically(delta_tot)
        else:
            _new_biomass=self.allocate_biomass_proportionately(delta_biomass_respire)
        event_flags.pop('_in_growing_season')

        #Run all other processes that need to occur
        for process in event_flags.items():
            if process[1]:
                _new_biomass=processes[process[0]](_new_biomass)

        self.plants=_new_biomass

    def _init_plants_from_grid(self,in_growing_season):
        ###This method initializes the plants in the PlantGrowth class
        #from the vegetation fields stored on the grid. This method
        #is only called if no initial plant array is parameterized 
        #as part of the PlantGrowth initialization.
        #Required parameters are a boolean inidicating if the plants are
        #in the active growing season.
        ###
        dtypes=[
            ('species','U10'),
            ('pid',int),
            ('cell_index',int),
            (('root','root_biomass'),float),
            (('leaf','leaf_biomass'),float),
            (('stem','stem_biomass'),float),
            (('storage','storage_biomass'), float),
            (('reproductive','repro_biomass'),float),
            ('plant_age',float),
            ('item_id',int)
        ]
        pidval=0
        plantlist=[]
        #Loop through grid cells
        for cell_index in range(self._grid.number_of_cells):
            cell_plants=self._grid['cell']['vegetation__plant_species'][cell_index]
            #Loop through list of plants stored on grid cell
            for plant in cell_plants:
                if plant == self.species_plant_factors['species']:
                    plantlist.append((plant,pidval,cell_index,0.0,0.0,0.0,0.0,0.0,0.0,0))
                    pidval += 1
        plant_array=np.array(plantlist, dtype=dtypes)
        plant_array=self.set_initial_biomass(plant_array,in_growing_season)
        return plant_array
        
    def allocate_biomass_dynamically(self,delta_tot):
        ###This method allocates new biomass according to the size-dependent
        #biomass allocation array calculated upon initiation of the PlantGrowth class.
        #The array is only valid for actively growing plants so this method is only
        #used during the growing season.
        #After initial allocation, storage redistribution, reallocation, and the 
        #minimum size check methods are called to adjust the biomass in each part before
        #saving. 
        #Required parameters are the new net biomass generated for the day
        ###

        #Interpolate values from biomass allocation array
        _new_biomass=self.plants
        delta_leaf_unit_root=np.interp(
            _new_biomass['root_biomass'], 
            self.biomass_allocation_array['prior_root_biomass'],
            self.biomass_allocation_array['delta_leaf_unit_root']
        )
        delta_stem_unit_root=np.interp(
            _new_biomass['root_biomass'], 
            self.biomass_allocation_array['prior_root_biomass'],
            self.biomass_allocation_array['delta_stem_unit_root']
        )
        
        #Calculate allocation           
        delta_root=delta_tot/(1+delta_leaf_unit_root+delta_stem_unit_root)
        delta_leaf=delta_leaf_unit_root*delta_root
        delta_stem=delta_stem_unit_root*delta_root
        
        #update biomass in plant array
        _new_biomass['root_biomass']=self.plants['root_biomass']+delta_root
        _new_biomass['leaf_biomass']=self.plants['leaf_biomass']+delta_leaf
        _new_biomass['stem_biomass']=self.plants['stem_biomass']+delta_stem
        
        #Adjust biomass allocation among storage and growth parts
        _new_biomass=self.redistribute_storage_biomass(_new_biomass)
        _new_biomass=self.adjust_biomass_allocation_towards_ideal(_new_biomass)
       
        #Kill plants that are too small to survive
        _new_biomass=self.kill_small_plants(_new_biomass)
        return _new_biomass

    def redistribute_storage_biomass(self,_new_biomass):
        ###This method redistributes biomass to and from the storage pool of biomass
        #If the growth part biomass is below the minimum and there is carbohydrate
        #avilable in the storage pool, the available biomass is distributed amongst
        #all growth parts. If the growth parts exceed the maximum size, biomass
        #is redistributed to the storage pool.
        #Required parameters are a structured plant array containing the newly updated
        #biomass values for the day.
        ###
        growth_min_biomass=self.species_grow_params['growth_min_biomass']
        growth_max_biomass=self.species_grow_params['growth_max_biomass']
        _growth_parts_biomass=self.sum_plant_parts(_new_biomass,parts='growth')
        filter=np.where(_growth_parts_biomass>0)
        conditions=[
            (_growth_parts_biomass<growth_min_biomass),
            (_growth_parts_biomass>=growth_min_biomass)&(_growth_parts_biomass<=growth_max_biomass),
            (_growth_parts_biomass>growth_max_biomass)
        ]
        adjustments=[
            np.minimum((growth_min_biomass-_growth_parts_biomass),(_new_biomass['storage_biomass']-self.species_grow_params['plant_part_min']['storage'])),
            np.zeros_like(_new_biomass['root_biomass']),
            (growth_max_biomass-_growth_parts_biomass)
        ]
        adjustment=np.select(conditions, adjustments)
        _new_biomass['storage_biomass'][filter]=_new_biomass['storage_biomass'][filter]-adjustment[filter]
        for part in self.growth_parts:
            _new_biomass[part][filter]=_new_biomass[part][filter]+adjustment[filter]*(_new_biomass[part][filter]/_growth_parts_biomass[filter])

        return _new_biomass

    def allocate_biomass_proportionately(self, delta_tot):
        ###This method allocates new net biomass amongst growth parts
        #proportionately based on the relative size of the part. This 
        #method is used outside of the growing season since some plant parts
        #may not be present while the plant is dormant. The storage redistribution
        #method is called after initial biomass allocation to redistribute storage
        #biomass to dormant growth parts as needed.
        #Required parameter is a numpy array of net biomass change to be applied
        # to the plant and it returns the structured array _new_biomass.
        _last_biomass=self.plants
        _new_biomass=_last_biomass
        _total_biomass=self.sum_plant_parts(_last_biomass,parts='growth')
        filter=np.where(_total_biomass!=0)
        for part in self.growth_parts:
            _new_biomass[part][filter]=(
                (_last_biomass[part][filter]/_total_biomass[filter])*
                delta_tot[filter]+_last_biomass[part][filter]
            )
        _new_biomass=self.redistribute_storage_biomass(_new_biomass)
        return _new_biomass

    def populate_biomass_allocation_array(self):
        aleaf,b1leaf,b2leaf=self.species_grow_params['root_to_leaf_coeffs']
        astem,b1stem,b2stem=self.species_grow_params['root_to_stem_coeffs']
        prior_root_biomass=np.arange(
            start=self.species_grow_params['plant_part_min']['root'],
            stop=self.species_grow_params['plant_part_max']['root']+0.1,
            step=0.1
        )
        length_of_array =len(prior_root_biomass)
        place_zeros=np.zeros(length_of_array)
        biomass_allocation_map=np.column_stack((prior_root_biomass,place_zeros,place_zeros,place_zeros,place_zeros,place_zeros))
        biomass_allocation_map=list(map(tuple,biomass_allocation_map))
        self.biomass_allocation_array=np.array(biomass_allocation_map,
            dtype=[('prior_root_biomass',float),('total_biomass',float),('delta_leaf_unit_root',float),('delta_stem_unit_root',float),('leaf_mass_frac',float),('stem_mass_frac',float)])

        #set up sympy equations
        rootsym=symbols('rootsym')              
        dleaf=diff(10**(aleaf+b1leaf*log(rootsym,10)+b2leaf*(log(rootsym,10))**2),rootsym)
        dstem=diff(10**(astem+b1stem*log(rootsym,10)+b2stem*(log(rootsym,10))**2),rootsym)
        #Generate numpy expressions and solve for rate change in leaf and stem biomass per unit mass of root
        fleaf=lambdify(rootsym,dleaf,'numpy')
        fstem=lambdify(rootsym,dstem,'numpy')
        self.biomass_allocation_array['delta_leaf_unit_root']=fleaf(self.biomass_allocation_array['prior_root_biomass'])
        self.biomass_allocation_array['delta_stem_unit_root']=fstem(self.biomass_allocation_array['prior_root_biomass'])
        _leaf_biomasss=10**(aleaf+b1leaf*np.log10(prior_root_biomass)+b2leaf*(np.log10(prior_root_biomass))**2)
        _stem_biomass=10**(astem+b1stem*np.log10(prior_root_biomass)+b2stem*(np.log10(prior_root_biomass))**2)
        self.biomass_allocation_array['total_biomass']=self.biomass_allocation_array['prior_root_biomass']+_leaf_biomasss+_stem_biomass
        self.biomass_allocation_array['leaf_mass_frac']=_leaf_biomasss/self.biomass_allocation_array['total_biomass']
        self.biomass_allocation_array['stem_mass_frac']=_stem_biomass/self.biomass_allocation_array['total_biomass']
    
    def adjust_biomass_allocation_towards_ideal(self, _new_biomass):
        _total_biomass=self.sum_plant_parts(_new_biomass, parts='growth')
        current_leaf_mass_frac=np.divide(_new_biomass['leaf_biomass'], _total_biomass, out=np.zeros_like(_total_biomass), where=~np.isclose(_total_biomass,np.zeros_like(_total_biomass)))
        current_stem_mass_frac=np.divide(_new_biomass['stem_biomass'], _total_biomass, out=np.zeros_like(_total_biomass), where=~np.isclose(_total_biomass,np.zeros_like(_total_biomass)))

        ideal_leaf_mass_frac=np.interp(
            _total_biomass, 
            self.biomass_allocation_array['total_biomass'],
            self.biomass_allocation_array['leaf_mass_frac']
        )
        ideal_stem_mass_frac=np.interp(
            _total_biomass, 
            self.biomass_allocation_array['total_biomass'],
            self.biomass_allocation_array['stem_mass_frac']
        )

        current_diff_leaf=ideal_leaf_mass_frac-current_leaf_mass_frac
        current_diff_stem=ideal_stem_mass_frac-current_stem_mass_frac

        _new_leaf_mass_frac=ideal_leaf_mass_frac-current_diff_leaf*(1-np.exp(-self.species_grow_params['p_max']*self.dt.astype(int)))
        _new_stem_mass_frac=ideal_stem_mass_frac-current_diff_stem*(1-np.exp(-self.species_grow_params['p_max']*self.dt.astype(int)))

        _new_biomass['root_biomass']=(1-_new_leaf_mass_frac-_new_stem_mass_frac)*_total_biomass
        _new_biomass['leaf_biomass']=_new_leaf_mass_frac*_total_biomass
        _new_biomass['stem_biomass']=_new_stem_mass_frac*_total_biomass
        
        return _new_biomass
        
    
    def set_event_flags(self,_current_jday):
        durationdict=self.species_duration_params
        dispersaldict=self.species_dispersal_params
        flags_to_test={
            '_in_growing_season': bool((_current_jday>durationdict['growing_season_start'])&(_current_jday<durationdict['growing_season_end'])),
            '_is_emergence_day': bool(_current_jday==durationdict['growing_season_start']),
            '_in_reproductive_period': bool((_current_jday>dispersaldict['reproduction_start'])&(_current_jday<durationdict['senescence_start'])),
            '_in_senescence_period': bool((_current_jday>=durationdict['senescence_start'])&(_current_jday<durationdict['growing_season_end'])),
            '_is_dormant_day': bool(_current_jday==durationdict['growing_season_end'])
        }
        return flags_to_test

    def kill_small_plants(self, _new_biomass):
        min_size=self.species_grow_params['total_min_biomass']
        total_biomass=self.sum_plant_parts(_new_biomass,parts='total')
        for part in self.all_parts:
            _new_biomass[part][total_biomass<min_size]=0.0
        return _new_biomass

    def sum_plant_parts(self, _new_biomass, parts='total'):
        if parts=='total':
            parts_dict=self.all_parts
        elif parts=='growth':
            parts_dict=self.growth_parts
        _new_tot=np.zeros_like(_new_biomass['root_biomass'])
        for part in parts_dict:
            _new_tot+=_new_biomass[part]
        return _new_tot

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
        self.record_plants.dataset['vegetation__storage_biomass'].values[:,self.time_ind]=self.plants['storage_biomass']
        self.record_plants.dataset['vegetation__repro_biomass'].values[:,self.time_ind]=self.plants['repro_biomass']
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
