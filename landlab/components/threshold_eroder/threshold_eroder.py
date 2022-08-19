#!/usr/bin/env python3
"""Created on Wed Aug  4 11:00:08 2021.

@author: benjamincampforts
"""


import numpy as np

from landlab import Component, RasterModelGrid

from .cfuncs import _thresholder


class ThresholdEroder(Component):

    """Threshold eroder.

    Threshold eroder that cuts off slopes at a given threshold slope (Sc) and assumes material to dissolve away

    .. math::

        S(S>Sc) = Sc

    To be coupled with FlowDirectorSteepest or PriorityFloodFlowRouter for the calculation of steepest
    slope at each timestep.

    Component written by Benjamin Campforts, 2021

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    slope_crit: float (default=1.)
        Critical slope [L/L]

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ThresholdEroder,PriorityFloodFlowRouter

    Define grid and initial topography:

        - 3x5 grid
        - east and west boundaries are open, north and south are closed
        - Initial topography is plane at base level on the boundaries and
          1m of elevation elsewhere (core)

    >>> mg = RasterModelGrid((5, 5))
    >>> mg.set_closed_boundaries_at_grid_edges(False, False, False, False)
    >>> z = np.array([0., 0., 0., 0., 0.,
    ...               0., 1., 1., 1., 0.,
    ...               0., 1., 10., 1., 0.,
    ...               0., 1., 1., 1., 0.,
    ...               0., 0., 0., 0., 0.])
    >>> _ = mg.add_field("topographic__elevation", z, at="node")

    Instantiate Flow director (steepest slope type) and TL hillslope diffuser

    >>> fdir = PriorityFloodFlowRouter(mg)
    >>> th_ero = ThresholdEroder(
    ...     mg,
    ...     slope_crit=0.6)

    Run the components for ten short timepsteps

    >>> for t in range(2):
    ...     fdir.run_one_step()
    ...     th_ero.run_one_step()


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    """

    _name = "ThresholdEroder"

    _unit_agnostic = True

    _cite_as = """@Article{gmd-13-3863-2020,
                  AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
                  TITLE = {BedrockLandslider 1.0: a hybrid landscape evolution model to simulate the impact of landslides and landslide-derived sediment on landscape evolution.},
                  JOURNAL = {Geoscientific Model Development},
                  VOLUME = {13},
                  YEAR = {2020},
                  NUMBER = {9},
                  PAGES = {3863--3886},
                  URL = {https://doi.org/10.5194/gmd-13-3863-2020},
                  DOI = {10.5194/gmd-13-3863-2020}
                  }"""

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/m",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of soil or weathered bedrock",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of the bedrock surface",
        },
    }

    def __init__(self, grid, slope_crit=1.0):

        """Initialize Threshold Eroder.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        slope_crit: float (default=1.)
            Critical slope [L/L]
        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that ThresholdEroder is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        # Store grid and parameters

        self._slope_crit = slope_crit

        # Create fields:
        # Elevation
        self._elev = self._grid.at_node["topographic__elevation"]
        # Downstream steepest slope at node:
        self._steepest = self._grid.at_node["topographic__steepest_slope"]
        self._r = self._grid.at_node["flow__receiver_node"]
        self._link_to_reciever = grid.at_node["flow__link_to_receiver_node"]
        if isinstance(grid, RasterModelGrid):
            self._link_lengths = grid.length_of_d8
        else:
            self._link_lengths = grid.length_of_link

        self._stack = self.grid.at_node["flow__upstream_node_order"]

        self.initialize_output_fields()

        if "soil__depth" in self._grid.at_node.keys():
            if "bedrock__elevation" not in self._grid.at_node.keys():
                raise Exception(
                    "If soil__depth is provided as a field, also bedrock__elevation mut be provided as a field"
                )
            self._soilFlag = True
            self._soil = self._grid.at_node["soil__depth"]
            self._bed = self._grid.at_node["bedrock__elevation"]
        else:
            self._soilFlag = False

    def erode(self):
        """Erode landscape to threshold and dissolve sediment.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        """
        _thresholder(
            self._stack,
            self._link_to_reciever,
            self._r,
            self._link_lengths,
            self._elev,
            self._slope_crit,
        )

        if self._soilFlag:
            self._bed[:] = np.minimum(self._bed[:], self._elev[:])
            self._soil[:] = self._elev[:] - self._bed[:]

    def run_one_step(self):
        """Advance one timestep.

        Advance threshold erosion component.

        Parameters
        ----------
        """
        self.erode()
