"""BiotaEvolver StreamHabitatPatch object.
"""

from landlab.components.biota_macroevolution import Zone
import numpy as np
from watershed import get_watershed_masks_with_area_threshold


class StreamHabitatPatch(Zone):

    def __init__(self, time, mask):
        super().__init__(time, mask)

    @staticmethod
    def get_patches_with_area_threshold(grid, area_threshold, time):

        grid.at_node['watershed'] = get_watershed_masks_with_area_threshold(
                grid, area_threshold)
        outlets = np.unique(grid.at_node['watershed'])
        outlets = np.delete(outlets, np.where(outlets == -1))

        patches = []

        for outlet in outlets:
            watershed_mask = grid.at_node['watershed'] == outlet
            patches.append(Zone(time, watershed_mask))

        return patches

    @classmethod
    def _get_habitat_patch_vectors(cls, prior_patches, new_patches, time,
                                   **kwargs):

        vectors = super(cls, StreamHabitatPatch)._get_habitat_patch_vectors(
                prior_patches, new_patches, time, **kwargs)

        return vectors
