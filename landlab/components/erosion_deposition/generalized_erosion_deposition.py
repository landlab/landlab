import numpy as np
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

        """
        super(_GeneralizedErosionDeposition, self).__init__(grid)

        self.flow_receivers = grid.at_node['flow__receiver_node']
        self.stack = grid.at_node['flow__upstream_node_order']
        self.elev = grid.at_node['topographic__elevation']
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
