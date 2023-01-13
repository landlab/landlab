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
    ):
        self.validate_plant_factors(species_params['plant_factors'])
        self.validate_duration_params(species_params['duration_params'])        
        self.validate_grow_params(species_params['grow_params'])
        
        self.species_plant_factors=species_params['plant_factors']
        self.species_duration_params=species_params['duration_params']
        self.species_grow_params=species_params['grow_params']        

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
        
    def disperse(self):
        self.form.dispersal()

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
        kmLVG = growdict['respiration_coefficient'][1] * pow(2,((_temperature - 25)/10))  
        kmSTG = growdict['respiration_coefficient'][2]* pow(2,((_temperature - 25)/10)) 
        kmRTG = growdict['respiration_coefficient'][0] * pow(2,((_temperature - 25)/10)) 
        #maintenance respiration per day from Teh 2006
        rmPrime = (kmLVG * _last_biomass['leaf_biomass']) + (kmSTG * _last_biomass['stem_biomass']) + (kmRTG * _last_biomass['root_biomass'])  
        #calculates respiration adjustment based on aboveground biomass, as plants age needs less respiration
        #THIS NEEDS TO BE UPDATED
        #plantAge = _last_biomass['leaf_biomass']/_last_biomass['leaf_biomass']
        #plant age dependence from Teh 2006 page 145    
        respMaint = rmPrime
        delta_respire=np.zeros_like(_glu_req)
        delta_respire[_glu_req!=0]=(-respMaint[_glu_req!=0])/_glu_req[_glu_req!=0]
        return delta_respire

    def senesce(self, plants):
        plants=self.habit.senesce(plants)
        return plants

    def set_initial_biomass(self, plants, in_growing_season):
        plants=self.habit.duration.set_initial_biomass(plants,in_growing_season)
        return plants
