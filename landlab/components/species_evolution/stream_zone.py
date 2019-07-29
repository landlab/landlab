#!/usr/bin/env python
"""StreamZoneManager of SpeciesEvolver.
"""
import numpy as np

from landlab.utils import get_watershed_masks_with_area_threshold
from .zone import Zone, ZoneManager


class StreamZoneManager(ZoneManager):
    """The nodes and attributes of the entities that species populate.

    A portion of the model domain. It is the model entity that recognizes the
    model grid changes relevant to a species.
    """

    def __init__(self, grid, min_drainage_area):
        """
        Parameters
        ----------
        grid : ModelGrid
            A Landlab ModelGrid.
        min_drainage_area : float
            The lower limit of the contributing drainage area of streams.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components.species_evolution import StreamZoneManager
        >>> mg = RasterModelGrid((3, 3))

        ZoneManager is instantiated with a model grid and an existing grid
        field.

        >>> zm = StreamZoneManager(mg, 3)
        >>> mg.at_node[zm.field_name].reshape(mg.shape)


        """
        self._grid = grid
        self._A_c = min_drainage_area

    def _create_zones(self):
        """Get zones using a drainage area threshold for streams.
        """
        # Get watershed outlets.

        ws_masks = get_watershed_masks_with_area_threshold(self._grid,
                                                           self._A_c)
        outlets_with_nulls = np.unique(ws_masks)
        outlets = np.delete(outlets_with_nulls,
                            np.where(outlets_with_nulls == -1))

        # Create a zone for each stream network within a watershed.

        zones = [None] * len(outlets)

        stream_mask = self._grid.at_node['drainage_area'] >= self._A_c

        for i, outlet in enumerate(outlets):
            ws_mask = ws_masks == outlet
            zone_mask = np.all([ws_mask, stream_mask], 0)
            zones[i] = Zone(self, zone_mask)

        return zones
