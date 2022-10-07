"""
Species class definition, composition classes, and factory methods to generate species classes. 
These are used by PlantGrowth to differentiate plant properties and processes for species.
"""

#########################################################
#Define composite classes here that the species will use
#########################################################

#Leaf retention classes and selection method
class Evergreen(object):
    def __init__(self):
        self.keep_green_parts=True

class Deciduous(object):
    def __init__(self):
        self.keep_green_parts=False

#Dispersal classes and selection method
class Clonal(object):
    def __init__(self):
        pass
    def dispersal(self):
        print('Create a new clone')

class Seed(object):
    def __init__(self):
        pass
    def disperse(self):
        print('New plant emerges from seed some distance within parent plant')

class Random(object):
    def __init__(self):
        pass
    def disperse(self):
        print('New plant randomly appears')

#Duration classes and selection method
class Duration(object):
    def __init__(self, **kwargs):
        pass

class Annual(Duration):
    def __init__(self):
        super().__init__()

    def enter_dormancy(self):
        print('I die at dormancy')
    
    def emerge(self,emerge_size=[0.01,0.1,0.5]):
        print('I start as a seedling between 100 and 200% of the minimum size')
        #emerge size should be the min size for annuals
        return self.initialize_biomass(emerge_size)
    
    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5]}):
        #This function provides the range of mass an initial annual plant can have
        #Since initalization of annuals is the same as emergence, use the same function
        init_size_min=sum(grow_params['plant_part_min'])
        init_size_max=2*init_size_min
        return init_size_min, init_size_max

class Perennial(Duration):
    def __init__(self, retention_val):
        self.retention=self.select_retention_class(retention_val)
    
    def select_retention_class(self, retention_val):
        retention={
            'evergreen':Evergreen,
            'deciduous':Deciduous
        }
        return retention[retention_val]

    def enter_dormancy(self):
        print('I move carbs among live parts around during dormancy')
    
    def emerge(self, emerge_min=[0.01,0.1,0.5]):
        print('I will transfer some stored carbohydrate to photosynthetic parts')

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        #This function provides the range of mass an initial perennial plant can have
        init_min_mass=sum(grow_params['plant_part_min'])
        init_max_mass=sum(grow_params['plant_part_max'])
        return [init_min_mass, init_max_mass]


#Growth habit classes and selection method
#Growth habit uses duration properties to assign dormancy and emergence methods
class Habit(object):
    def __init__(self,duration_val, retention_val):
        self.duration=self.select_duration_class(duration_val, retention_val)

    def select_duration_class(self, duration_val, retention_val):
        duration={
            'annual':Annual,
            'perennial':Perennial(retention_val)
        }
        return duration[duration_val]   

class Forbherb(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val, retention_val)
        self.green_parts=('leaf', 'stem')

    def enter_dormancy(self):
        self.duration.enter_dormancy()

    def emerge(self):
        # Use this to move carbohydrate to aboveground biomass from storage organs and roots
        self.duration.emerge()

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        self.duration.initialize_biomass(grow_params)

class Graminoid(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val,retention_val)

        self.green_parts=('leaf','stem')

    def enter_dormancy(self):
        self.duration.enter_dormancy()

    def emerge(self):
        self.duration.emerge()

    def initialize_biomass(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        init_mass_min, init_mass_max=self.duration.initialize_biomass(grow_params)
        return init_mass_min, init_mass_max

class Shrub(Habit):
    def __init__(self, duration_val, retention_val='None'):
        super().__init__(duration_val,retention_val)
        
        self.green_parts=('leaf')

class Tree(Habit):
    def __init__(self, duration_val, retention_val):
        super().__init__(duration_val, retention_val='None')

        self.green_parts=('leaf')

class Vine(Habit):
    pass

#Growth form classes and selection method
class Bunch(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Limited lateral branching due to clumping')

class Colonizing(Random):
    def __init__(self):
        pass
    def branch(self):
        print('No branching annual')

class Multiplestems(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Create two or more main stems at or near soil surface')

class Rhizomatous(Clonal):
    def __init__(self):
        pass
    def branch(self):
        print('Tiller via rhizomes')

class Singlecrown(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Herbaceous plant with one persistent base')

class Singlestem(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Plant develops one stem like a tree or a corn plant')

class Stoloniferous(Clonal):
    def __init__(self):
        pass
    def branch(self):
        print('Branches via stolons')

class Thicketforming(Seed):
    def __init__(self):
        pass
    def branch(self):
        print('Limited lateral branching due to dense thickets')

#Shape and orientation classes and selection method
class Climbing(object):
    pass

class Conical(object):
    pass

class Decumbent(object):
    pass

class Erect(object):
    pass

class Irregular(object):
    pass

class Oval(object):
    pass

class Prostrate(object):
    pass

class Rounded(object):
    pass

class Semierect(object):
    pass

class Vase(object):
    pass


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
