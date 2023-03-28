#Dispersal classes and selection method
class Clonal(object):
    def __init__(self, disperse_params):
        self.specific_length=disperse_params['specific_length']
        self.max_runner_length=disperse_params['max_runner_length']
    
    def disperse(self, plants):
        max_runner_length=(plants['repro']-min_size)/self.specific_length
        pup_dist = rng.uniform(low=plants['root_sys_width']/2, high=self.max_runner_length, size=plants.size)
        pup_azimuth = rng.uniform(low=0, high=360, size=plants.size)
        plants['pup_x_loc']=np.empty_like(plants['repro'])
        plants['pup_y_loc']=np.empty_like(plants['repro'])

        filter=np.where(pup_dist >= max_runner_length)
            
        plants['pup_x_loc'][filter] = pup_dist[filter]*np.cos(pup_azimuth[filter])+plants['x_loc'][filter]
        plants['pup_y_loc'][filter] = pup_dist[filter]*np.sin(pup_azimuth)[filter]+plants['y_loc'][filter]
        plants['pup_runner_length'] = max_runner_length
        return plants
    
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