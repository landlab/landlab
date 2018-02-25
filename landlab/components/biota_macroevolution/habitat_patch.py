"""TODO: Description.
"""

import numpy as np
from uuid import uuid4
from watershed import get_watershed_masks_with_area_threshold


def get_patches_with_drainage_area_threshold(grid, area_threshold, time):
        
    grid.at_node['watershed'] = get_watershed_masks_with_area_threshold(grid,
                area_threshold)
    outlets = np.unique(grid.at_node['watershed'])
    outlets = np.delete(outlets, np.where(outlets == -1))
            
    patches = []
    
    for outlet in outlets:
        watershed_mask = grid.at_node['watershed'] == outlet
        patches.append(HabitatPatch(time, watershed_mask))
        
    return patches


def _get_stream_patch_vectors(self, prior_patches, current_patches, prior_grid,
                              time):
    
    vectors = {'origin': [], 'destinations': [], 'node_was_captured': []}

    ps = np.vstack(prior_patches)
    cs = np.vstack(current_patches)
    
    for cp in current_patches:
        vectors['origin'].append(cp)

        priors_in_cp = np.all([cp == ps, ps], 0)
        indices = np.argwhere(priors_in_cp)[:,0]
        priors_overlapping_cp = current_patches[indices]

        if len(priors_overlapping_cp) == 1:
            pp = priors_overlapping_cp[0]
            currents_in_prior = np.all([pp == cs, cs], 0)
            indices = np.argwhere(currents_in_prior)[:,0]
            overlapping_currents = current_patches[indices]

            if len(overlapping_currents) == 1:
                # cp is actually pp.
                cp = pp

        elif len(priors_overlapping_cp) > 1:
            mean_z = []
            prior_z = prior_grid.at_node['topographic__elevation']
            prior_stream_nodes = prior_grid.at_node['stream']
            for op in priors_overlapping_cp:
                op_nodes = op.at_time[time]
                stream_nodes = np.all([prior_stream_nodes, op_nodes], 0)
                mean_z.append(np.mean(prior_z[stream_nodes]))
                
            cp = priors_overlapping_cp[np.argmin(mean_z)]
        
        vectors['destinations'].append(priors_overlapping_cp)
        
        vectors['node_was_captured'].append(priors_in_cp[indices])
    
    return vectors


class HabitatPatch(object):
    
    def __init__(self, time, nodes):
        """
        """
        
        self.at_time = {time: [nodes]}
        self.identifier = uuid4()
        self.species = []

    def update(self, time, nodes):
                
        self.at_time[time] = nodes
