#!/usr/bin/env python3
"""Created on Wed Aug  4 11:00:08 2021.

@author: benjamincampforts
"""
from landlab import Component
from landlab import RasterModelGrid

from .cfuncs import _thresholder


class ThresholdEroder(Component):
    """Threshold eroder.

    Threshold eroder that cuts off slopes at a given threshold slope (Sc) and
    assumes material to dissolve away

    .. math::

        S(S>Sc) = Sc

    To be coupled with :class:`~.flow_director_steepest.FlowDirectorSteepest` or
    :class:`~.priority_flood_flow_router.PriorityFloodFlowRouter` for the
    calculation of steepest slope at each timestep. Note that ThresholdEroder
    run_one_step() cuts off slopes and computes new elevations based on the
    steepest slopes as calculated by the FlowDirectorSteepest or
    PriorityFloodFlowRouter. If slopes over the entire model grid need be set
    to a threshold slope, several iterations of running the flow router and the
    Threshold eroder are required.

    Component written by Benjamin Campforts, 2022

    Parameters
    ----------
    grid : ModelGrid
        Landlab ModelGrid object
    slope_crit: float, optional
        Critical slope [L/L]

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import ThresholdEroder, PriorityFloodFlowRouter

    Define grid and initial topography:

    - 3x5 grid
    - east and west boundaries are open, north and south are closed
    - Initial topography is plane at base level on the boundaries and
      1m of elevation elsewhere (core)

    >>> mg = RasterModelGrid((5, 5))
    >>> mg.set_closed_boundaries_at_grid_edges(False, False, False, False)
    >>> z = np.array(
    ...     [
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...         [0.0, 1.0, 1.0, 1.0, 0.0],
    ...         [0.0, 1.0, 10.0, 1.0, 0.0],
    ...         [0.0, 1.0, 1.0, 1.0, 0.0],
    ...         [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     ]
    ... )
    >>> _ = mg.add_field("topographic__elevation", z, at="node")

    Instantiate Flow director (steepest slope type) and TL hillslope diffuser

    >>> fdir = PriorityFloodFlowRouter(mg)
    >>> th_ero = ThresholdEroder(mg, slope_crit=0.6)

    Run the components for ten short timepsteps

    >>> for t in range(2):
    ...     fdir.run_one_step()
    ...     th_ero.run_one_step()
    ...


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    """

    _name = "ThresholdEroder"

    _unit_agnostic = True

    _cite_as = """
    @Article{gmd-13-3863-2020,
      AUTHOR = {Campforts B., Shobe C.M., Steer P., Vanmaercke M., Lague D., Braun J.},
      TITLE = {BedrockLandslider 1.0: a hybrid landscape evolution model to
               simulate the impact of landslides and landslide-derived sediment on
               landscape evolution.},
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
        slope_crit: float, optional
            Critical slope [L/L]
        """
        super().__init__(grid)

        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that ThresholdEroder is compatible "
                "with route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        # Store grid and parameters
        self._slope_crit = slope_crit

        # Link lengths depending on raster type:
        if isinstance(grid, RasterModelGrid):
            self._link_lengths = grid.length_of_d8
        else:
            self._link_lengths = grid.length_of_link

        # Create fields
        self.initialize_output_fields()

        if (
            "soil__depth" in self._grid.at_node
            and "bedrock__elevation" not in self._grid.at_node
        ):
            raise ValueError(
                "If soil__depth is provided as a field, "
                "bedrock__elevation must also be provided as a field"
            )

    def erode(self):
        """Erode landscape to threshold and dissolve sediment."""
        _thresholder(
            self.grid.at_node["flow__upstream_node_order"],
            self.grid.at_node["flow__link_to_receiver_node"],
            self._grid.at_node["flow__receiver_node"],
            self._link_lengths,
            self._grid.at_node["topographic__elevation"],
            self._slope_crit,
        )

        if "soil__depth" in self._grid.at_node:
            self._grid.at_node["bedrock__elevation"].clip(
                None,
                self._grid.at_node["topographic__elevation"],
                out=self._grid.at_node["bedrock__elevation"],
            )
            self._grid.at_node["soil__depth"][:] = (
                self._grid.at_node["topographic__elevation"]
                - self._grid.at_node["bedrock__elevation"]
            )

    def run_one_step(self):
        """Advance threshold erosion component one timestep."""
        self.erode()
