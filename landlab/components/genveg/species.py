"""
Species class definition, composition classes, and factory methods to generate species classes. 
These are used by PlantGrowth to differentiate plant properties and processes for species.
"""
from .habit import *
from .form import *
from .shape import *
from .photosynthesis import *
import numpy as np

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

        self.validate_plant_factors(species_params['plant_factors'])
        self.validate_duration_params(species_params['duration_params'])        
        self.validate_grow_params(species_params['grow_params'])
        
        self.species_plant_factors=species_params['plant_factors']
        self.species_duration_params=species_params['duration_params']
        self.species_grow_params=species_params['grow_params']
        self.species_dispersal_params=species_params['dispersal_params']   

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

    def validate_grow_params(self,grow_params):
        multipart_vars=['respiration_coefficient', 'glucose_requirement','root_to_leaf_coeffs','root_to_stem_coeffs']
        for vars in multipart_vars:
            if len(grow_params[vars])<3:
                msg='Must include respiration coefficients for at least roots, leaves, and stems'
                raise ValueError(msg)

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
            'bunch':Bunch(),
            'colonizing':Colonizing(),
            'multiple_stems':Multiplestems(),
            'rhizomatous':Rhizomatous(),
            'single_crown':Singlecrown(),
            'single_stem':Singlestem(),
            'stoloniferous':Stoloniferous(),
            'thicket_forming':Thicketforming()
        }
        return form[form_val]

    def select_shape_class(self,shape_val):
        shape={
            'climbing':Climbing(),
            'conical':Conical(),
            'decumbent':Decumbent(),
            'erect':Erect(),
            'irregular':Irregular(),
            'oval':Oval(),
            'prostrate':Prostrate(),
            'rounded':Rounded(),
            'semi_erect':Semierect(),
            'vase':Vase()
        }
        return shape[shape_val]

    def test_output(self):
        return self.habit.duration.emerge_plants

    def branch(self):
        self.form.branch()
        
    def disperse(self, plants):
        #decide how to parameterize reproductive schedule, make repro event
        #right now we are just taking 20% of available storage and moving to 
        available_stored_biomass=plants['storage_biomass']-self.species_grow_params['plant_part_min']['storage']
        plants['repro_biomass']=plants['repro_biomass']+0.2*(available_stored_biomass)
        plants['storage_biomass']=plants['storage_biomass']-0.2*(available_stored_biomass)
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
        plants=self.habit.duration.set_initial_biomass(plants,in_growing_season)
        return plants

        
