import numpy as np


def _default_axis_names(n_dims):
    """Returns a tuple of the default axis names."""
    _DEFAULT_NAMES = ('z', 'y', 'x')
    return _DEFAULT_NAMES[- n_dims:]


def _default_axis_units(n_dims):
    """Returns a tuple of the default axis units."""
    return ('-', ) * n_dims


class BaseGrid(object):
    def __init__(self, axis_name=None, axis_units=None, node_status=None):
        self._axis_name = axis_name or _default_axis_names(self.ndim)
        self._axis_units = axis_units or _default_axis_units(self.ndim)

        self._num_nodes = 0
        self._num_cells = 0
        self._num_links = 0
        self._num_faces = 0

        self._node_id_at_cells = np.empty(self._num_nodes)

        StatusGrid.__init__(self, node_status)
        super(BaseGrid, self).__init__(node_status)

    @property
    def ndim(self):
        return 2

    @property
    def number_of_cells(self):
        """Number of cells.
        """
        return self._num_cells

    @property
    def number_of_nodes(self):
        """Number of nodes.
        """
        return self._num_nodes

    @property
    def node_id_at_cells(self):
        return self._node_id_at_cells

    def _reset_list_of_active_links(self):
        """
        Creates or resets a list of active links. We do this by sweeping
        through the given lists of from and to nodes, and checking the status
        of these as given in the node_status list. A link is active if both its
        nodes are active interior points, or if one is an active interior and
        the other is an active boundary.
        """
        status_at_link_start = self.node_status[self._node_id_at_link_start]
        status_at_link_end = self.node_status[self._node_id_at_link_end]

        (self._active_link_ids, ) = np.where(
            link_is_active(status_at_link_start, status_at_link_end))

        self._number_of_active_links = len(self._active_link_ids)
        self._number_of_active_faces = len(self._active_link_ids)
        self._node_id_at_active_link_start = self._node_id_at_link_start[self._active_link_ids]
        self._node_id_at_active_link_end = self._node_id_at_link_end[self._active_link_ids]

        # Set up active inlink and outlink matrices
        self._setup_active_inlink_and_outlink_matrices()


