"""
Species class definition, composition classes, and factory methods to generate species classes. 
These are used by PlantGrowth to differentiate plant properties and processes for species.
"""
from .habit import *
from .form import *
from .shape import *
from .photosynthesis import *
import numpy as np
from sympy import symbols, diff, lambdify, log

#Define species class that inherits composite class methods
class Species(object):
    def __init__(
        self, 
        species_params
        ):
        self.all_parts=list(species_params['grow_params']['glucose_requirement'].keys())
        self.growth_parts=self.all_parts.copy()
        self.growth_parts.remove('storage')
        self.growth_parts.remove('reproductive')
        self.abg_parts=self.growth_parts.copy()
        self.abg_parts.remove('root')

        self.validate_plant_factors(species_params['plant_factors'])
        self.validate_duration_params(species_params['duration_params'])        
        
        species_params=self.calculate_derived_params(species_params)
        
        self.species_plant_factors=species_params['plant_factors']
        self.species_duration_params=species_params['duration_params']
        self.species_grow_params=species_params['grow_params']
        self.species_dispersal_params=species_params['dispersal_params']  
        self.species_mort_params=species_params['mortality_params']
        self.species_morph_params=species_params['morph_params'] 

        self.populate_biomass_allocation_array()
        
        self.habit=self.select_habit_class(
            self.species_plant_factors['growth_habit'], 
            self.species_plant_factors['duration'],
            self.species_plant_factors['leaf_retention']
            )
        self.form=self.select_form_class(self.species_plant_factors['growth_form'])
        self.shape=self.select_shape_class(self.species_plant_factors['shape'])
        self.photosynthesis=self.select_photosythesis_type(self.species_plant_factors['p_type'])

    def validate_plant_factors(self,plant_factors):
        plant_factor_options={
            'species':[],
            'growth_habit':['forb_herb','graminoid','shrub','tree','vine'],
            'monocot_dicot':['monocot','dicot'],
            'angio_gymno':['angiosperm','gymnosperm'],
            'duration':['annual','perennial'],
            'leaf_retention':['deciduous','evergreen'],
            'growth_form':['bunch','colonizing','multiple_stems','rhizomatous','single_crown','single_stem','stoloniferous','thicket_forming'],
            'shape':['climbing','columnar','conical','decumbent','erect','irregular','oval','prostrate','rounded','semi_erect','vase'],
            'p_type':['C3','C4']
        }

        for key in plant_factors:
            try: 
                opt_list=plant_factor_options[key]
                if opt_list:
                    if plant_factors[key] not in opt_list:
                        msg='Invalid '+str(key)+' option'
                        raise ValueError(msg)
            except:
                msg='Unexpected variable name in species parameter dictionary. Please check input parameter file.'
                raise ValueError(msg)

    def validate_duration_params(self,duration_params):
        if (duration_params['growing_season_start'] < 0) | (duration_params['growing_season_start']  > 366):
            msg='Growing season beginning must be integer values between 1-365'
            raise ValueError(msg)
        elif (duration_params['growing_season_end'] < duration_params['growing_season_start']) | (duration_params['growing_season_end'] >366):
            msg='Growing season end must be between 1-365 and greater than the growing season beginning'
            raise ValueError(msg)
        elif (duration_params['senescence_start'] < duration_params['growing_season_start']) | (duration_params['senescence_start'] > duration_params['growing_season_end']):
            msg='Start of senescence must be within the growing season'
            raise ValueError(msg)

    def calculate_derived_params(self, species_params):
        species_params['morph_params']['max_crown_area']=np.pi/4*species_params['morph_params']['max_shoot_sys_width']**2
        species_params['morph_params']['min_crown_area']=np.pi/4*species_params['morph_params']['min_shoot_sys_width']**2
        species_params['morph_params']['max_root_area']=np.pi/4*species_params['morph_params']['max_root_sys_width']**2
        species_params['morph_params']['min_root_area']=np.pi/4*species_params['morph_params']['min_root_sys_width']**2
        species_params['morph_params']['max_vital_volume']=species_params['morph_params']['max_crown_area']*species_params['morph_params']['max_height']
        species_params['morph_params']['area_per_stem']=species_params['morph_params']['max_crown_area']/species_params['morph_params']['max_n_stems']
        species_params['morph_params']['min_abg_aspect_ratio']=species_params['morph_params']['max_height']/species_params['morph_params']['min_shoot_sys_width']
        species_params['morph_params']['max_abg_aspect_ratio']=species_params['morph_params']['max_height']/species_params['morph_params']['max_shoot_sys_width']

        
        sum_vars=[
            ['max_total_biomass','plant_part_max', self.all_parts],
            ['max_growth_biomass','plant_part_max', self.growth_parts],
            ['max_abg_biomass', 'plant_part_max', self.abg_parts],
            ['min_total_biomass', 'plant_part_min', self.all_parts],
            ['min_growth_biomass','plant_part_min', self.growth_parts],
            ['min_abg_biomass','plant_part_min', self.abg_parts]
        ]
        for sum_var in sum_vars:
            species_params['grow_params'][sum_var[0]]=0
            for part in sum_var[2]:
                species_params['grow_params'][sum_var[0]] += species_params['grow_params'][sum_var[1]][part]
       
        species_params['morph_params']['biomass_packing']=species_params['grow_params']['max_growth_biomass']/species_params['morph_params']['max_vital_volume']

        return species_params

    def select_photosythesis_type(self, p_type):
        photosynthesis_options={
            'C3': C3(),
            'C4': C4(),
            'cam': Cam()
        }
        return photosynthesis_options[p_type]

    def select_habit_class(self, habit_val, duration, retention_val):
        habit={
            'forb_herb':Forbherb(self.species_grow_params, duration),
            'graminoid':Graminoid(self.species_grow_params, duration),
            'shrub':Shrub(self.species_grow_params, duration, retention_val),
            'tree':Tree(self.species_grow_params, duration, retention_val),
            'vine':Vine(self.species_grow_params, duration, retention_val)
        }
        return habit[habit_val]
    
    def select_form_class(self, form_val):
        form={
            'bunch':Bunch(self.species_dispersal_params, self.species_grow_params),
            'colonizing':Colonizing(self.species_dispersal_params, self.species_grow_params),
            'multiple_stems':Multiplestems(self.species_dispersal_params, self.species_grow_params),
            'rhizomatous':Rhizomatous(self.species_dispersal_params, self.species_grow_params),
            'single_crown':Singlecrown(self.species_dispersal_params, self.species_grow_params),
            'single_stem':Singlestem(self.species_dispersal_params, self.species_grow_params),
            'stoloniferous':Stoloniferous(self.species_dispersal_params, self.species_grow_params),
            'thicket_forming':Thicketforming(self.species_dispersal_params, self.species_grow_params)
        }
        return form[form_val]

    def select_shape_class(self,shape_val):
        shape={
            'climbing':Climbing(self.species_morph_params, self.species_grow_params),
            'conical':Conical(self.species_morph_params,self.species_grow_params),
            'decumbent':Decumbent(self.species_morph_params,self.species_grow_params),
            'erect':Erect(self.species_morph_params,self.species_grow_params),
            'irregular':Irregular(self.species_morph_params,self.species_grow_params),
            'oval':Oval(self.species_morph_params,self.species_grow_params),
            'prostrate':Prostrate(self.species_morph_params,self.species_grow_params),
            'rounded':Rounded(self.species_morph_params,self.species_grow_params),
            'semi_erect':Semierect(self.species_morph_params,self.species_grow_params),
            'vase':Vase(self.species_morph_params,self.species_grow_params)
        }
        return shape[shape_val]

    def populate_biomass_allocation_array(self):
        root2leaf=self.species_grow_params['root_to_leaf']
        root2stem=self.species_grow_params['root_to_stem']
        prior_root_biomass=np.arange(
            start=self.species_grow_params['plant_part_min']['root'],
            stop=self.species_grow_params['plant_part_max']['root']+0.1,
            step=0.1
        )
        length_of_array =len(prior_root_biomass)
        place_zeros=np.zeros(length_of_array)
        biomass_allocation_map=np.column_stack((prior_root_biomass,place_zeros,place_zeros,place_zeros,place_zeros,place_zeros,place_zeros))
        biomass_allocation_map=list(map(tuple,biomass_allocation_map))
        self.biomass_allocation_array=np.array(biomass_allocation_map,
            dtype=[('prior_root_biomass',float),('total_biomass',float),('delta_leaf_unit_root',float),('delta_stem_unit_root',float),('leaf_mass_frac',float),('stem_mass_frac',float),('abg_biomass', float)])

        #set up sympy equations
        rootsym=symbols('rootsym')              
        dleaf=diff(10**(root2leaf['a']+root2leaf['b1']*log(rootsym,10)+root2leaf['b2']*(log(rootsym,10))**2),rootsym)
        dstem=diff(10**(root2stem['a']+root2stem['b1']*log(rootsym,10)+root2stem['b2']*(log(rootsym,10))**2),rootsym)
        #Generate numpy expressions and solve for rate change in leaf and stem biomass per unit mass of root
        fleaf=lambdify(rootsym,dleaf,'numpy')
        fstem=lambdify(rootsym,dstem,'numpy')
        self.biomass_allocation_array['delta_leaf_unit_root']=fleaf(self.biomass_allocation_array['prior_root_biomass'])
        self.biomass_allocation_array['delta_stem_unit_root']=fstem(self.biomass_allocation_array['prior_root_biomass'])
        _leaf_biomasss=10**(root2leaf['a']+root2leaf['b1']*np.log10(prior_root_biomass)+root2leaf['b2']*(np.log10(prior_root_biomass))**2)
        _stem_biomass=10**(root2stem['a']+root2stem['b1']*np.log10(prior_root_biomass)+root2stem['b2']*(np.log10(prior_root_biomass))**2)
        self.biomass_allocation_array['total_biomass']=self.biomass_allocation_array['prior_root_biomass']+_leaf_biomasss+_stem_biomass
        self.biomass_allocation_array['leaf_mass_frac']=_leaf_biomasss/self.biomass_allocation_array['total_biomass']
        self.biomass_allocation_array['stem_mass_frac']=_stem_biomass/self.biomass_allocation_array['total_biomass']
        self.biomass_allocation_array['abg_biomass']=_leaf_biomasss+_stem_biomass
    
    def test_output(self):
        return self.habit.duration.emerge_plants

    def branch(self):
        self.form.branch()
        
    
    def sum_plant_parts(self, _new_biomass, parts='total'):
        if parts=='total':
            parts_dict=self.all_parts
        elif parts=='growth':
            parts_dict=self.growth_parts
        elif parts=='aboveground':
            parts_dict=self.abg_parts
        _new_tot=np.zeros_like(_new_biomass['root_biomass'])
        for part in parts_dict:
            _new_tot+=_new_biomass[part]
        return _new_tot
    
    def disperse(self, plants):
        #decide how to parameterize reproductive schedule, make repro event
        #right now we are just taking 20% of available storage and moving to
        filter=np.where(self.sum_plant_parts(plants, parts='growth') >= (self.species_dispersal_params['min_size_dispersal']*self.species_grow_params['max_growth_biomass'])) 
        available_stored_biomass=plants['storage_biomass']-self.species_grow_params['plant_part_min']['storage']
        plants['repro_biomass'][filter]=plants['repro_biomass'][filter]+0.2*(available_stored_biomass[filter])
        plants['storage_biomass'][filter]=plants['storage_biomass'][filter]-0.2*(available_stored_biomass[filter])
        plants=self.form.disperse(plants)
        return plants

    def enter_dormancy(self, plants):
        plants=self.habit.enter_dormancy(plants)
        return plants

    def emerge(self, plants):
        plants=self.habit.duration.emerge(plants)
        return plants
    
    def photosynthesize(self, _par, _last_biomass, _glu_req,_daylength):
        delta_tot=self.photosynthesis.photosynthesize(_par, self.species_grow_params, _last_biomass, _glu_req, _daylength)
        return delta_tot

    def respire(self, _temperature, _last_biomass, _glu_req):
        growdict=self.species_grow_params
        #repiration coefficient temp dependence from Teh 2006
        maint_respire=np.zeros_like(_glu_req)
        #can i create a dictionary with alias?
        for part in self.all_parts:
            maint_respire+=growdict['respiration_coefficient'][part]*_last_biomass[part]

        maint_respire_adj=maint_respire*2**((_temperature - 25)/10)

        delta_biomass_respire=np.zeros_like(_glu_req)
        delta_biomass_respire[_glu_req!=0]=(-maint_respire_adj[_glu_req!=0])/_glu_req[_glu_req!=0]
        return delta_biomass_respire

    def senesce(self, plants):
        plants=self.habit.senesce(plants)
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        growdict=self.species_grow_params
        morphdict=self.species_morph_params
        plants['storage_biomass']=rng.uniform(low=growdict['plant_part_min']['storage'],high=growdict['plant_part_max']['storage'],size=plants.size)
        plants['repro_biomass']=rng.uniform(low=growdict['plant_part_min']['reproductive'], high=growdict['plant_part_max']['reproductive'], size=plants.size)
        crown_area=self.shape.calc_crown_area_from_shoot_width(plants['shoot_sys_width'])
        plants['shoot_sys_height']=self.habit.set_initial_height(morphdict['min_height'], morphdict['max_height'], crown_area.size)
        log_vital_volume=np.log10(crown_area*plants['shoot_sys_height'])
        log_abg_biomass_ideal=(log_vital_volume/np.log10(morphdict['max_vital_volume']))*np.log10(growdict['max_abg_biomass']/1000)
        total_biomass=np.interp(
            ((10**log_abg_biomass_ideal)*1000), 
            self.biomass_allocation_array['abg_biomass'],
            self.biomass_allocation_array['total_biomass']
        )
        plants['root_biomass'],plants['leaf_biomass'],plants['stem_biomass']=self.habit.duration._solve_biomass_allocation(total_biomass)
        plants['root_sys_width'] = self.shape.calc_root_sys_width(plants['shoot_sys_width'], plants['shoot_sys_height'])
        plants['n_stems']=self.form.set_initial_branches(morphdict['max_n_stems'], crown_area.size)
        plants = self.habit.duration.set_initial_biomass(plants, in_growing_season)
        return plants
    
    def update_morphology(self, plants):
        abg_biomass=self.sum_plant_parts(plants, parts='aboveground')
        dims=self.shape.calc_abg_dims_from_biomass(abg_biomass)
        plants['shoot_sys_width']=dims[0]
        plants['shoot_sys_height']=dims[1]
        plants['root_sys_width'] = self.shape.calc_root_sys_width(plants['shoot_sys_width'])
        return plants


