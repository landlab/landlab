from .leaf_retention import *
import numpy as np
from scipy.optimize import fsolve
rng = np.random.default_rng()
#from landlab.components.genveg import PlantGrowth 

#Duration classes and selection method
class Duration(object):
    def __init__(self, species_grow_params, green_parts, retention_val='deciduous'):
        self.min_mass=sum(species_grow_params['plant_part_min'])
        self.new_max_mass=2*self.min_mass
        self.max_mass=sum(species_grow_params['plant_part_max'])
        self.allocation_coeffs=species_grow_params['root_to_leaf_coeffs']+species_grow_params['root_to_stem_coeffs']
        self.green_parts=green_parts
        self.plant_array_dict={
            'root':['root_biomass'],
            'leaf':['leaf_biomass'],
            'stem':['stem_biomass']
        }
        self.retention=self.select_retention_class(retention_val)
    
    def select_retention_class(self, retention_val):
        retention={
            'evergreen':Evergreen(),
            'deciduous':Deciduous()
        }
        return retention[retention_val]

    def set_initial_biomass(self, plants, in_growing_season):
        total_biomass_ideal=rng.uniform(low=self.min_mass,high=self.max_mass,size=plants.size)
        plants['root_biomass'],plants['leaf_biomass'],plants['stem_biomass']=self._solve_biomass_allocation(total_biomass_ideal, self.allocation_coeffs)
        plants['storage_biomass']=np.zeros_like(plants['root_biomass'])
        plants['plant_age']=np.ones_like(plants['root_biomass'])
        if not in_growing_season and not self.retention.keep_green_parts:
            for part in self.green_parts:
                plants[self.plant_array_dict[part]]=np.zeros_like(plants[self.plant_array_dict[part]])
        return plants
    
    def set_new_biomass(self, plants):
        print('I create new plants')
        total_biomass_ideal=rng.uniform(low=self.min_mass,high=self.new_max_mass,size=plants.size)
        plants['root_biomass'],plants['leaf_biomass'],plants['stem_biomass']=self._solve_biomass_allocation(total_biomass_ideal, self.allocation_coeffs)
        plants['storage_biomass']=np.zeros_like(plants['root_biomass'])
        plants['plant_age']=np.zeros_like(plants['root_biomass'])
        return plants
    
    def _solve_biomass_allocation(self, total_biomass, solver_coeffs):
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

class Annual(Duration):
    def __init__(self, species_grow_params):
        green_parts=('root','leaf','stem')
        super().__init__(species_grow_params, green_parts)

    def senesce(self, plants):
        print('I start to lose biomass during senescence periood')
        plants['root_biomass'] = plants['root_biomass'] - (plants['root_biomass'] * 0.02)
        plants['leaf_biomass'] = plants['leaf_biomass'] - (plants['leaf_biomass'] * 0.02)
        plants['stem_biomass'] = plants['stem_biomass'] - (plants['stem_biomass'] * 0.02)
        return plants
    
    def enter_dormancy(self, plants):
        plants['root_biomass'] = np.zeros_like(plants['root_biomass'])
        plants['leaf_biomass'] = np.zeros_like(plants['leaf_biomass'])
        plants['stem_biomass'] = np.zeros_like(plants['stem_biomass']) 
        return plants
    
    def emerge(self, plants):
        print('I emerge from dormancy')
        plants=self.set_new_biomass(plants)
        return plants
    
class Perennial(Duration):
    def __init__(self, species_grow_params, green_parts, retention_val):
        super().__init__(species_grow_params, green_parts, retention_val)
        self.min_part_mass={
            'leaf':species_grow_params['plant_part_min'][1],
            'stem':species_grow_params['plant_part_min'][2]
        }
    
    def senesce(self, plants):
        #copied from annual for testing. This needs to be updated
        plants['root_biomass'] = plants['root_biomass'] - (plants['root_biomass'] * 0.02)
        plants['leaf_biomass'] = plants['leaf_biomass'] - (plants['leaf_biomass'] * 0.02)
        plants['stem_biomass'] = plants['stem_biomass'] - (plants['stem_biomass'] * 0.02)
        return plants

    def enter_dormancy(self, plants):
        print('I kill green parts at end of growing season')
    
    def emerge(self, plants):
        print('I emerge from dormancy')
        #for part in self.green_parts:
        #    plants[self.plant_array_dict[part]]=rng.uniform(low=self.min_part_mass[part],high=self.min_part_mass[part]*2,size=plants.size)
        #    plants['root_biomass']=plants['root_biomass']-plants[self.plant_array_dict[part]]
        plants['leaf_biomass']=rng.uniform(low=self.min_part_mass['leaf'],high=self.min_part_mass['leaf']*2,size=plants.size)
        plants['stem_biomass']=rng.uniform(low=self.min_part_mass['stem'],high=self.min_part_mass['stem']*2,size=plants.size)
        plants['root_biomass']=plants['leaf_biomass']-plants['stem_biomass']
        return plants
