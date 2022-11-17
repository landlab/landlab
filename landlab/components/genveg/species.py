"""
Species class definition, composition classes, and factory methods to generate species classes. 
These are used by PlantGrowth to differentiate plant properties and processes for species.
"""
from .habit import *
from .form import *
from .shape import *

#Define species class that inherits composite class methods
class Species(object):
    def __init__(
        self, 
        species_params={
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
                },
            },
    ):
        self.validate_plant_factors(species_params['plant_factors'])
        self.validate_grow_params(species_params['growparams'])
        
        self.species_params=species_params
        
        self.habit=self.select_habit_class(
            self.species_params['plant_factors']['growth_habit'], 
            self.species_params['plant_factors']['duration'], 
            self.species_params['plant_factors']['leaf_retention']
            )
        self.form=self.select_form_class(self.species_params['plant_factors']['growth_form'])
        self.shape=self.select_shape_class(self.species_params['plant_factors']['shape'])
        
        self.set_initial_biomass()

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

    def validate_grow_params(self,growparams):
        multipart_vars=['respiration_coefficient', 'glucose_requirement','root_to_leaf_coeffs','root_to_stem_coeffs']
        if (growparams['growing_season_start'] < 0) | (growparams['growing_season_start']  > 366):
            msg='Growing season beginning must be integer values between 1-365'
            raise ValueError(msg)
        elif (growparams['growing_season_end'] < growparams['growing_season_start']) | (growparams['growing_season_end'] >366):
            msg='Growing season end must be between 1-365 and greater than the growing season beginning'
            raise ValueError(msg)
        elif (growparams['senescence_start'] < growparams['growing_season_start']) | (growparams['senescence_start'] > growparams['growing_season_end']):
            msg='Start of senescence must be within the growing season'
            raise ValueError(msg)
        for vars in multipart_vars:
            if len(growparams[vars])<3:
                msg='Must include respiration coefficients for at least roots, leaves, and stems'
                raise ValueError(msg)

    def select_habit_class(self, habit_val, duration, retention):
        habit={
            'forb_herb':Forbherb,
            'graminoid':Graminoid,
            'shrub':Shrub,
            'tree':Tree,
            'vine':Vine
        }
        return habit[habit_val](duration, retention)
    
    def select_form_class(self, form_val):
        form={
            'bunch':Bunch,
            'colonizing':Colonizing,
            'multiple_stems':Multiplestems,
            'rhizomatous':Rhizomatous,
            'single_crown':Singlecrown,
            'single_stem':Singlestem,
            'stoloniferous':Stoloniferous,
            'thicket_forming':Thicketforming
        }
        return form[form_val]()

    def select_shape_class(self,shape_val):
        shape={
            'climbing':Climbing,
            'conical':Conical,
            'decumbent':Decumbent,
            'erect':Erect,
            'irregular':Irregular,
            'oval':Oval,
            'prostrate':Prostrate,
            'rounded':Rounded,
            'semi_erect':Semierect,
            'vase':Vase
        }
        return shape[shape_val]()

    def branch(self):
        self.form.branch()
        
    def disperse(self):
        self.form.dispersal()

    def enter_dormancy(self):
        self.habit.enter_dormancy()

    def emerge(self):
        self.habit.emerge()

    def set_initial_biomass(self):
        min_mass, max_mass=self.habit.initialize_biomass()
        self.species_params['growparams']['init_biomass']=[min_mass,max_mass]
