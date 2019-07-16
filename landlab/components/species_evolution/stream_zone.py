"""StreamZone SpeciesEvolver object.
"""
import numpy as np

from landlab.utils import get_watershed_masks_with_area_threshold

from .zone import Zone

class StreamZone(Zone):
    """The nodes and attributes of the entities that species populate.

    A portion of the model domain. It is the model entity that recognizes the
    model grid changes relevant to a species.
    """
    _required_params = ['critical_area']

    @classmethod
    def get_zones(cls, grid, object_params):
        """Get zones using a drainage area threshold for streams.

        Requires

        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        object_params: dictionary


        """
        cls._check_params(object_params)

        A_c = object_params['critical_area']

        # Get watershed outlets.

        watershed_masks = get_watershed_masks_with_area_threshold(grid, A_c)
        outlets_with_nulls = np.unique(watershed_masks)
        outlets = np.delete(outlets_with_nulls,
                            np.where(outlets_with_nulls == -1))

        # Create a zone at stream nodes in each watershed.

        zones = [None] * len(outlets)

        stream_mask = grid.at_node['drainage_area'] >= A_c

        for i, outlet in enumerate(outlets):
            watershed_mask = watershed_masks == outlet
            zone_mask = np.all([watershed_mask, stream_mask], 0)
            zones[i] = StreamZone(zone_mask)

        return zones
