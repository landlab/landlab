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

def retention_select(string_val):
    retention={
        'evergreen':Evergreen,
        'deciduous':Deciduous
    }
    return retention[string_val]()

#Dispersal classes and selection method
class Clonal(object):
    def __init__(self):
        pass
    def dispersal(self):
        print('Create a new clone')

class Seed(object):
    def __init__(self):
        pass
    def dispersal(self):
        print('New plant emerges from seed some distance within parent plant')

class Random(object):
    def __init__(self):
        pass
    def dispersal(self):
        print('New plant randomly appears')

def dispersal_select(string_val):
    dispersal={
        'clonal':Clonal,
        'seed':Seed,
        'random':Random
    }
    return dispersal[string_val]()

#Duration classes and selection method
class Annual(object):
    def __init__(self, retention):
        self.retention=retention_select(retention)
    
    def dormancy(self):
        print('I die at dormancy')
    
    def emergence(self,emerge_size=[0.01,0.1,0.5]):
        print('I start as a seedling between 100 and 200% of the minimum size')
        #emerge size should be the min size for annuals
        emerge_size_min=sum(emerge_size)
        emerge_size_max=2*emerge_size_min
        return emerge_size_min, emerge_size_max
    
    def initial(self, grow_params={'plant_part_min':[0.01,0.1,0.5]}):
        #This function provides the range of mass an initial annual plant can have
        #Since initalization of annuals is the same as emergence, use the same function
        init_mass_min, init_mass_max=self.emergence(grow_params['plant_part_min'])
        return init_mass_min, init_mass_max

class Perennial(object):
    def __init__(self, retention):
        self.retention=retention_select(retention)
    def dormancy(self):
        print('I move carbs among live parts around during dormancy')
    
    def emergence(self, emerge_min=[0.01,0.1,0.5]):
        print('I will transfer some stored carbohydrate to photosynthetic parts')

    def initial(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        #This function provides the range of mass an initial perennial plant can have
        init_min_mass=sum(grow_params['plant_part_min'])
        init_max_mass=sum(grow_params['plant_part_max'])
        return [init_min_mass, init_max_mass]

def duration_select(string_val, retention):
    duration={
        'annual':Annual,
        'perennial':Perennial
    }
    return duration[string_val](retention)

#Growth habit classes and selection method
#Growth habit uses duration properties to assign dormancy and emergence methods
class Forbherb(object):
    def __init__(self, duration, retention='None'):
        self.duration=duration_select(duration, retention)
        self.green_parts=('leaf', 'stem')

    def dormancy(self):
        self.duration.dormancy()

    def emergence(self):
        # Use this to move carbohydrate to aboveground biomass from storage organs and roots
        pass

    def initial_habit(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        self.duration.initial(grow_params)

class Graminoid(object):
    def __init__(self, duration, retention='None'):
        self.duration=duration_select(duration, retention)
        self.green_parts=('leaf','stem')

    def dormancy(self):
        self.duration.dormancy()

    def initial_habit(self, grow_params={'plant_part_min':[0.01,0.1,0.5], 'plant_part_max':[2,2,2]}):
        init_mass_min, init_mass_max=self.duration.initial(grow_params)
        return init_mass_min, init_mass_max

class Shrub(object):
    def __init__(self, duration, retention='None'):
        self.duration=duration_select(duration, retention)
        self.green_parts=('leaf')

class Tree(object):
    def __init__(self, duration, retention='None'):
        self.duration=duration_select(duration, retention)
        self.green_parts=('leaf')

class Vine(object):
    pass

def habit_select(string_val, duration, retention):
    habit={
        'forb_herb':Forbherb,
        'graminoid':Graminoid,
        'shrub':Shrub,
        'tree':Tree,
        'vine':Vine
    }
    return habit[string_val](duration, retention)

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

def form_select(string_val):
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
    return form[string_val]()

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

def shape_select(string_val):
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
    return shape[string_val]()


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
        self.species_params=species_params
        #read vegetation parameters and select 
        self.habit=habit_select(
            self.species_params['plant_factors']['growth_habit'], 
            self.species_params['plant_factors']['duration'], 
            self.species_params['plant_factors']['leaf_retention']
            )
        self.form=form_select(self.species_params['plant_factors']['growth_form'])
        self.shape=shape_select(self.species_params['plant_factors']['shape'])

        growdict=species_params['growparams']
        #Check validity of growing season beginning and ending days, senescence beginning, and calculate growing season length
        end=growdict['growing_season_end']
        begin=growdict['growing_season_start']
        senes=growdict['senescence_start']
        if (begin > 0) & (begin  < 366) & (end > 0) & (end < 366) & (end > begin) & (senes > begin) & (senes < end):
            length=end-begin+1
        elif (begin < 1) | (begin > 365):
            msg='Growing season beginning must be between 1-365'
            raise ValueError(msg)
        elif (end < 1) | (end > 365):
            msg='Growing season end must be between 1-365'
            raise ValueError(msg)
        elif end < begin:
            msg='Growing season beginning must be before growing season end'
            raise ValueError(msg)
        elif (senes < begin) | (senes > end):
            msg='Start of senescence must be within the growing season'
            raise ValueError(msg)
        else:
            msg='Growing season beginning and end must be integer values between 1-365'
            raise ValueError(msg)
        self.species_params['growparams']['growing_season_length']=length
        self.initial_mass()

    def branch(self):
        self.form.branch()
        
    def disperse(self):
        self.form.dispersal()

    def dormancy(self):
        self.habit.dormancy()

    def emerge(self):
        self.habit.emergence()

    def initial_mass(self):
        min_mass, max_mass=self.habit.initial_habit()
        self.species_params['growparams']['init_biomass']=[min_mass,max_mass]
        #return [min_mass, max_mass]


#retention='Deciduous'
#duration='Annual'
#habit='Graminoid'
#form='Rhizomatous'
#shape='Erect'

#ammophila=Species(duration, retention, habit, form, shape)
#ammophila.disperse()
#ammophila.branch()
#ammophila.emerge()