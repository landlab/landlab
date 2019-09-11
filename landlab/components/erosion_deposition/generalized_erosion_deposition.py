import numpy as np

from landlab import Component, RasterModelGrid
from landlab.utils.return_array import return_array_at_node

DEFAULT_MINIMUM_TIME_STEP = 0.001  # default minimum time step duration


class _GeneralizedErosionDeposition(Component):
    """ Base class for erosion-deposition type components.

    More documenation here.
    """

    _name = "ErosionDeposition"

    _input_var_names = set((
        "flow__receiver_node",
        "flow__upstream_node_order",
        "topographic__steepest_slope",
        "surface_water__discharge",
    ))

    _output_var_names = set(("topographic__elevation", "sediment__flux"))

    _var_units = {
        "flow__receiver_node": "-",
        "flow__upstream_node_order": "-",
        "topographic__steepest_slope": "-",
        "surface_water__discharge": "m**2/s",
        "topographic__elevation": "m",
        "sediment__flux": "m3/s",
    }

    _var_mapping = {
        "flow__receiver_node": "node",
        "flow__upstream_node_order": "node",
        "topographic__steepest_slope": "node",
        "surface_water__discharge": "node",
        "topographic__elevation": "node",
        "sediment__flux": "node",
    }

    _var_doc = {
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "topographic__steepest_slope": "Topographic slope at each node",
        "surface_water__discharge": "Water discharge at each node",
        "topographic__elevation": "Land surface topographic elevation",
        "sediment__flux": "TODO",
    }

    def __init__(
        self,
        grid,
        m_sp,
        n_sp,
        phi,
        F_f,
        v_s,
        discharge_field="surface_water__discharge",
        dt_min=DEFAULT_MINIMUM_TIME_STEP,
    ):
        """Initialize the ErosionDeposition model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        m_sp : float
            Discharge exponent (units vary)
        n_sp : float
            Slope exponent (units vary)
        phi : float
            Sediment porosity [-].
        F_f : float
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load. Defaults to zero.
        v_s : float
            Effective settling velocity for chosen grain size metric [L/T].
        discharge_field : float, field name, or array
            Discharge [L^2/T].
        dt_min : float, optional
            Only applies when adaptive solver is used. Minimum timestep that
            adaptive solver will use when subdividing unstable timesteps.
            Default values is 0.001. [T].
        """
        super(_GeneralizedErosionDeposition, self).__init__(grid)

        self._flow_receivers = grid.at_node["flow__receiver_node"]
        self._stack = grid.at_node["flow__upstream_node_order"]
        self._topographic__elevation = grid.at_node["topographic__elevation"]
        self._slope = grid.at_node["topographic__steepest_slope"]
        self._link_to_reciever = grid.at_node["flow__link_to_receiver_node"]
        self._cell_area_at_node = grid.cell_area_at_node

        if isinstance(grid, RasterModelGrid):
            self._link_lengths = grid.length_of_d8
        else:
            self._link_lengths = grid.length_of_link

        if "sediment__flux" in grid.at_node:
            self._qs = grid.at_node["sediment__flux"]
        else:
            self._qs = grid.add_zeros("sediment__flux", at="node", dtype=float)

        self._q = return_array_at_node(grid, discharge_field)

        # Create arrays for sediment influx at each node, discharge to the
        # power "m", and deposition rate
        self._qs_in = np.zeros(grid.number_of_nodes)
        self._Q_to_the_m = np.zeros(grid.number_of_nodes)
        self._S_to_the_n = np.zeros(grid.number_of_nodes)
        self._depo_rate = np.zeros(grid.number_of_nodes)

        # store other constants
        self._m_sp = float(m_sp)
        self._n_sp = float(n_sp)
        self._phi = float(phi)
        self._v_s = float(v_s)
        self._dt_min = dt_min
        self._F_f = float(F_f)

        if phi >= 1.0:
            raise ValueError("Porosity must be < 1.0")

        if F_f > 1.0:
            raise ValueError("Fraction of fines must be <= 1.0")

        if phi < 0.0:
            raise ValueError("Porosity must be > 0.0")

        if F_f < 0.0:
            raise ValueError("Fraction of fines must be > 0.0")

        self._verify_output_fields()

    def _update_flow_link_slopes(self):
        """Updates gradient between each core node and its receiver.

        Used to update slope values between sub-time-steps, when we do not
        re-run flow routing.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import FlowAccumulator
        >>> from landlab.components.erosion_deposition.generalized_erosion_deposition import _GeneralizedErosionDeposition
        >>> rg = RasterModelGrid((3, 4))
        >>> z = rg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = rg.x_of_node + rg.y_of_node
        >>> fa = FlowAccumulator(rg, flow_director='FlowDirectorD8')
        >>> fa.run_one_step()
        >>> rg.at_node['topographic__steepest_slope'][5:7]
        array([ 1.41421356,  1.41421356])
        >>> sp = _GeneralizedErosionDeposition(rg, phi=0.1, v_s=0.001,
        ...                                    m_sp=0.5, n_sp=1.0, F_f=0)
        >>> z *= 0.1
        >>> sp._update_flow_link_slopes()
        >>> rg.at_node['topographic__steepest_slope'][5:7]
        array([ 0.14142136,  0.14142136])
        """
        self._slope[:] = (
            self._topographic__elevation
            - self._topographic__elevation[self._flow_receivers]
        ) / self._link_lengths[self._link_to_reciever]

    def _calc_hydrology(self):
        self._Q_to_the_m[:] = np.power(self._q, self._m_sp)
