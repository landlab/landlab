import numpy as np
from six import string_types

from landlab import Component
from landlab import RasterModelGrid
from landlab.utils.return_array import return_array_at_node
from .cfuncs import calculate_qs_in
DEFAULT_MINIMUM_TIME_STEP = 0.001 # default minimum time step duration

class _GeneralizedErosionDeposition(Component):
    """

    """

    _name= 'ErosionDeposition'

    _input_var_names = (
        'flow__receiver_node',
        'flow__upstream_node_order',
        'topographic__steepest_slope',
        'drainage_area',
    )

    _output_var_names = (
        'topographic__elevation'
    )

    _var_units = {
        'flow__receiver_node': '-',
        'flow__upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'drainage_area': 'm**2',
        'topographic__elevation': 'm',
    }

    _var_mapping = {
        'flow__receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'drainage_area': 'node',
        'topographic__elevation': 'node',
    }

    _var_doc = {
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'topographic__steepest_slope':
            'Topographic slope at each node',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'topographic__elevation':
            'Land surface topographic elevation',
    }

    def __init__(self, grid, m_sp=None, n_sp=None,
                 phi=None, F_f=None, v_s=None,
                 dt_min=DEFAULT_MINIMUM_TIME_STEP):
        """Initialize the ErosionDeposition model.

        Parameters
        ----------
        grid
        m_sp
        n_sp
        phi
        F_f
        v_s
        dt_min : float, optional
            Default values is 0.001.
        """
        super(_GeneralizedErosionDeposition, self).__init__(grid)

        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.topographic__elevation = grid.at_node['topographic__elevation']
        self.slope = grid.at_node['topographic__steepest_slope']
        self.link_to_reciever = grid.at_node['flow__link_to_receiver_node']

        if isinstance(grid, RasterModelGrid):
            self.link_lengths = grid.length_of_d8
        else:
            self.link_lengths = grid.length_of_link

        try:
            self.qs = grid.at_node['sediment__flux']
        except KeyError:
            self.qs = grid.add_zeros(
                'sediment__flux', at='node', dtype=float)
        try:
            self.q = grid.at_node['surface_water__discharge']
        except KeyError:
            self.q = grid.add_zeros(
                'surface_water__discharge', at='node', dtype=float)

        # Create arrays for sediment influx at each node, discharge to the
        # power "m", and deposition rate
        self.qs_in = np.zeros(grid.number_of_nodes)
        self.Q_to_the_m = np.zeros(grid.number_of_nodes)
        self.S_to_the_n = np.zeros(grid.number_of_nodes)
        self.depo_rate = np.zeros(grid.number_of_nodes)

        # store other constants
        self.m_sp = float(m_sp)
        self.n_sp = float(n_sp)
        self.phi = float(phi)
        self.v_s = float(v_s)
        self.dt_min = dt_min
        self.F_f = float(F_f)

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
        self.slope[:] = ((self.topographic__elevation -
                          self.topographic__elevation[self.flow_receivers]) /
                         self.link_lengths[self.link_to_reciever])



    #three choices for erosion methods:
    def simple_stream_power(self):
        """Use non-threshold stream power.

        simple_stream_power uses no entrainment or erosion thresholds,
        and uses either q=A^m or q=Q^m depending on discharge method. If
        discharge method is None, default is q=A^m.
        """
        if self.discharge_method is None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        else:
            if self.discharge_method == 'area_field':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.drainage_area, self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if isinstance(self.discharge_field, string_types):
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
            else:
                raise TypeError('Unrecognized discharge_method')

        #TODO: FIGURE OUT WHY WE BOTHER TO CALC AND STORE BOTH ES/ER AND
        #**_EROSION_TERM?


    def threshold_stream_power(self):
        """Use stream power with entrainment/erosion thresholds.

        threshold_stream_power works the same way as simple SP but includes
        user-defined thresholds for sediment entrainment and bedrock erosion.
        """
        self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.discharge_method is None:
            self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
        else:
            if self.discharge_method == 'area_field':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.drainage_area, self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if isinstance(self.discharge_field, string_types):
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')


    def stochastic_hydrology(self):
        """Allows custom area and discharge fields, no default behavior.

        stochastic_hydrology forces the user to supply either an array or
        field name for either drainage area or discharge, and will not
        default to q=A^m.
        """
        self.Q_to_the_m = np.zeros(len(self.grid.at_node['drainage_area']))
        if self.discharge_method is None:
            raise TypeError('Supply a discharge method to use stoc. hydro!')
        else:
            if self.discharge_method == 'drainage_area':
                if self.area_field is not None:
                    if isinstance(self.area_field, string_types):
                        self.drainage_area = self._grid.at_node[self.area_field]
                    elif len(self.area_field) == self.grid.number_of_nodes:
                        self.drainage_area = np.array(self.area_field)
                    else:
                        raise TypeError('Supplied type of area_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
                self.Q_to_the_m[:] = np.power(self.grid.at_node['drainage_area'], self.m_sp)
            elif self.discharge_method == 'discharge_field':
                if self.discharge_field is not None:
                    if isinstance(self.discharge_field, string_types):
                        self.q[:] = self._grid.at_node[self.discharge_field]
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    elif len(self.discharge_field) == self.grid.number_of_nodes:
                        self.q[:] = np.array(self.discharge_field)
                        self.Q_to_the_m[:] = np.power(self.q, self.m_sp)
                    else:
                        raise TypeError('Supplied type of discharge_field ' +
                                'was not recognised, or array was ' +
                                'not nnodes long!')
            else:
                raise ValueError('Specify discharge method for stoch hydro!')
