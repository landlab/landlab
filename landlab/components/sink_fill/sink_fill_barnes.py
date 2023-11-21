#!/usr/env/python

"""fill_sinks_barnes.py.

Fill sinks in a landscape to the brim, following the Barnes et al.
(2014) algorithms.
"""


import numpy as np

from landlab.components import LakeMapperBarnes
from landlab.utils.return_array import return_array_at_node


class SinkFillerBarnes(LakeMapperBarnes):
    """Uses the Barnes et al (2014) algorithms to replace pits in a topography
    with flats, or optionally with very shallow gradient surfaces to allow
    continued draining.

    This component is NOT intended for use iteratively as a model runs;
    rather, it is to fill in an initial topography. If you want to repeatedly
    fill pits as a landscape develops, you are after the LakeMapperBarnes
    component. If you want flow paths on your filled landscape, manually run a
    FlowDirector and FlowAccumulator for yourself.

    The locations and depths etc. of the fills will be tracked, and properties
    are provided to access this information.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnes, R., Lehman, C., Mulla, D. (2014). Priority-flood: An optimal
    depression-filling and watershed-labeling algorithm for digital elevation
    models. Computers and Geosciences  62(C), 117 - 127.
    https://dx.doi.org/10.1016/j.cageo.2013.04.024

    **Additional References**

    None Listed

    """

    _name = "SinkFillerBarnes"

    _unit_agnostic = True

    _cite_as = """
    @article{BARNES2014117,
        title = {Priority-flood: An optimal depression-filling and
                 watershed-labeling algorithm for digital elevation models},
        journal = "Computers & Geosciences",
        volume = "62",
        pages = "117 - 127",
        year = "2014",
        issn = "0098-3004",
        doi = "https://doi.org/10.1016/j.cageo.2013.04.024",
        url = "http://www.sciencedirect.com/science/article/pii/S0098300413001337",
        author = "Richard Barnes and Clarence Lehman and David Mulla",
        keywords = "Pit filling, Terrain analysis, Hydrology, Drainage network, Modeling, GIS"
        }"""

    _info = {
        "sediment_fill__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Depth of sediment added at eachnode",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        method="D8",
        fill_flat=False,
        ignore_overfill=False,
    ):
        """Initialise the component.

        Parameters
        ----------
        grid : ModelGrid
            A grid.
        surface : field name at node or array of length node
            The surface to fill.
        method : {'Steepest', 'D8'}
            Whether or not to recognise diagonals as valid flow paths, if a raster.
            Otherwise, no effect.
        fill_flat : bool
            If True, pits will be filled to perfectly horizontal. If False, the new
            surface will be slightly inclined (at machine precision) to give
            steepest descent flow paths to the outlet, once they are calculated.
        ignore_overfill : bool
            If True, suppresses the Error that would normally be raised during
            creation of a gentle incline on a fill surface (i.e., if not
            fill_flat). Typically this would happen on a synthetic DEM where more
            than one outlet is possible at the same elevation. If True, the
            was_there_overfill property can still be used to see if this has
            occurred.

        """
        if "flow__receiver_node" in grid.at_node and grid.at_node[
            "flow__receiver_node"
        ].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that SinkFillerBarnes is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        # Most of the functionality of this component is directly inherited
        # from SinkFillerBarnes, so
        super().__init__(
            grid,
            surface=surface,
            method=method,
            fill_flat=fill_flat,
            fill_surface=surface,
            redirect_flow_steepest_descent=False,
            reaccumulate_flow=False,
            ignore_overfill=ignore_overfill,
            track_lakes=True,
        )
        # note we will always track the fills, since we're only doing this
        # once... Likewise, no need for flow routing; this is not going to
        # get used dynamically.
        self._supplied_surface = return_array_at_node(grid, surface).copy()
        # create the only new output field:
        self._sed_fill_depth = self._grid.add_zeros(
            "node", "sediment_fill__depth", clobber=True
        )

    def run_one_step(self):
        """Fills the surface to remove all pits.

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import SinkFillerBarnes, FlowAccumulator
        >>> mg = RasterModelGrid((5, 6))
        >>> for edge in ("left", "top", "bottom"):
        ...     mg.status_at_node[mg.nodes_at_edge(edge)] = mg.BC_NODE_IS_CLOSED
        ...
        >>> mg.at_node["topographic__elevation"] = [
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ...     [0.0, 2.1, 1.1, 0.6, 1.6, 0.0],
        ...     [0.0, 2.0, 1.0, 0.5, 1.5, 0.0],
        ...     [0.0, 2.2, 1.2, 0.7, 1.7, 0.0],
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ... ]
        >>> z = mg.at_node["topographic__elevation"]
        >>> z_init = z.copy()
        >>> sfb = SinkFillerBarnes(mg, method="Steepest")  # , surface=z

        TODO: once return_array_at_node is fixed, this example should also
        take surface... GIVE IT surface=z  !!

        >>> sfb.run_one_step()
        >>> z_out = [
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ...     [0.0, 2.1, 1.5, 1.5, 1.6, 0.0],
        ...     [0.0, 2.0, 1.5, 1.5, 1.5, 0.0],
        ...     [0.0, 2.2, 1.5, 1.5, 1.7, 0.0],
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ... ]
        >>> np.allclose(z.reshape(mg.shape), z_out)
        True
        >>> fill_map = [
        ...     [-1, -1, -1, -1, -1, -1],
        ...     [-1, -1, 16, 16, -1, -1],
        ...     [-1, -1, 16, 16, -1, -1],
        ...     [-1, -1, 16, 16, -1, -1],
        ...     [-1, -1, -1, -1, -1, -1],
        ... ]
        >>> np.all(sfb.fill_map.reshape(mg.shape) == fill_map)
        True
        >>> np.all(sfb.fill_at_node == (sfb.fill_map > -1))
        True
        >>> sfb.was_there_overfill  # everything fine with slope adding
        False

        >>> fr = FlowAccumulator(mg, flow_director="D4")  # routing now works
        >>> fr.run_one_step()
        >>> np.all(mg.at_node["flow__sink_flag"][mg.core_nodes] == 0)
        True
        >>> drainage_area = [
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ...     [0.0, 1.0, 2.0, 3.0, 1.0, 1.0],
        ...     [0.0, 1.0, 4.0, 9.0, 10.0, 10.0],
        ...     [0.0, 1.0, 2.0, 1.0, 1.0, 1.0],
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ... ]
        >>> np.allclose(mg.at_node["drainage_area"].reshape(mg.shape), drainage_area)
        True

        Test two pits and a flat fill:

        >>> from collections import deque
        >>> z[:] = mg.node_x.max() - mg.node_x
        >>> z[[10, 23]] = 1.1  # raise "guard" exit nodes
        >>> z[7] = 2.0  # is a lake on its own
        >>> z[9] = 0.5
        >>> z[15] = 0.3
        >>> z[14] = 0.6  # [9, 14, 15] is a lake
        >>> z[22] = 0.9  # a non-contiguous lake node also draining to 16
        >>> z_init = z.copy()
        >>> sfb = SinkFillerBarnes(mg, method="Steepest", fill_flat=True)
        >>> sfb.run_one_step()
        >>> sfb.fill_dict == {8: deque([7]), 16: deque([15, 9, 14, 22])}
        True
        >>> sfb.number_of_fills
        2
        >>> sfb.fill_outlets == [16, 8]
        True
        >>> np.allclose(sfb.fill_areas, np.array([4.0, 1.0]))  # same order
        True

        Unlike the LakeMapperBarnes equivalents, fill_depths and fill_volumes
        are always available through this component:

        >>> fill_depths = [
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ...     [0.0, 1.0, 0.0, 0.5, 0.0, 0.0],
        ...     [0.0, 0.0, 0.4, 0.7, 0.0, 0.0],
        ...     [0.0, 0.0, 0.0, 0.0, 0.1, 0.0],
        ...     [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ... ]
        >>> np.allclose(sfb.fill_depths.reshape(mg.shape), fill_depths)
        True
        >>> np.allclose(sfb.fill_volumes, np.array([1.7, 1.0]))
        True

        Note that with a flat fill, we can't drain the surface. The surface
        is completely flat, so each and every core node within the fill
        becomes a sink.

        >>> fr.run_one_step()
        >>> where_is_filled = np.where(sfb.fill_map > -1, 1, 0)
        >>> np.all(
        ...     mg.at_node["flow__sink_flag"][mg.core_nodes]
        ...     == where_is_filled[mg.core_nodes]
        ... )
        True

        (Note that the fill_map does not think that the perimeter nodes are
        sinks, since they haven't changed elevation. In contrast, the
        FlowAccumulator *does* think they are, because these nodes are where
        flow terminates.)
        """
        if "flow__receiver_node" in self._grid.at_node and self._grid.at_node[
            "flow__receiver_node"
        ].size != self._grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has "
                "not verified that SinkFillerBarnes is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        super().run_one_step()

        self._sed_fill_depth[:] = self._surface - self._supplied_surface

    @property
    def fill_dict(self):
        """Return a dictionary where the keys are the outlet nodes of each
        filled area, and the values are deques of nodes within each.

        Items are not returned in ID order.
        """
        return super().lake_dict

    @property
    def fill_outlets(self):
        """Returns the outlet for each filled area, not necessarily in ID
        order."""
        return super().lake_outlets

    @property
    def number_of_fills(self):
        """Return the number of individual filled areas."""
        return super().number_of_lakes

    @property
    def fill_map(self):
        """Return an array of ints, where each filled node is labelled with its
        outlet node ID.

        Nodes not in a filled area are labelled with
        BAD_INDEX_VALUE (default -1).
        """
        return super().lake_map

    @property
    def fill_at_node(self):
        """Return a boolean array, True if the node is filled, False
        otherwise."""
        return super().lake_at_node

    @property
    def fill_depths(self):
        """Return the change in surface elevation at each node this step."""
        return self._sed_fill_depth

    @property
    def fill_areas(self):
        """A nlakes-long array of the area of each fill.

        The order is the same as that of the keys in fill_dict, and of
        fill_outlets.
        """
        return super().lake_areas

    @property
    def fill_volumes(self):
        """A nlakes-long array of the volume of each fill.

        The order is the same as that of the keys in fill_dict, and of
        fill_outlets.
        """
        fill_vols = np.empty(self.number_of_fills, dtype=float)
        col_vols = self._grid.cell_area_at_node * self._sed_fill_depth
        for i, fillnodes in enumerate(self.fill_dict.values()):
            fill_vols[i] = col_vols[fillnodes].sum()
        return fill_vols
