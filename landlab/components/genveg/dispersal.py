import numpy as np
rng = np.random.default_rng()

#Dispersal classes and selection method
class Repro(object):
    def __init__(self, grow_params):
        self.min_size=grow_params['min_growth_biomass']

class Clonal(Repro):
    def __init__(self, disperse_params, grow_params):
        super().__init__(grow_params)
        self.unit_cost=disperse_params['unit_cost_dispersal']
        self.max_dist_dispersal=disperse_params['max_dist_dispersal']
        
    
    def disperse(self, plants):
        #plants['pup_x_loc'] = np.full_like(plants['root'], np.nan)
        #plants['pup_y_loc'] = np.full_like(plants['root'], np.nan)
        #plants['pup_cost'] = np.full_like(plants['root'], np.nan)
        runner_length=np.zeros_like(plants['root'])
        available_carb=plants['reproductive']-2*self.min_size
        runner_length[available_carb>0]=available_carb[available_carb>0]/self.unit_cost
        pup_dist = rng.uniform(low=plants['root_sys_width']/2, high=self.max_dist_dispersal, size=plants.size)
        pup_azimuth = np.deg2rad(rng.uniform(low=0, high=360, size=plants.size))

        filter=np.where(pup_dist < runner_length)
            
        plants['pup_x_loc'][filter] = pup_dist[filter]*np.cos(pup_azimuth[filter])+plants['x_loc'][filter]
        plants['pup_y_loc'][filter] = pup_dist[filter]*np.sin(pup_azimuth)[filter]+plants['y_loc'][filter]
        plants['pup_cost'][filter] = pup_dist[filter]*self.unit_cost
        return plants
    
class Seed(Repro):
    def __init__(self):
        pass
    def disperse(self):
        print('New plant emerges from seed some distance within parent plant')

class Random(Repro):
    def __init__(self):
        pass
    def disperse(self):
        print('New plant randomly appears')