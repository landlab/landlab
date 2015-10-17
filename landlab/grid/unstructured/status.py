import numpy as np

# Define the boundary-type codes

#: Indicates a node is *core*.
CORE_NODE = 0

#: Indicates a boundary node is has a fixed values.
FIXED_VALUE_BOUNDARY = 1

#: Indicates a boundary node is has a fixed gradient.
FIXED_GRADIENT_BOUNDARY = 2

#: Indicates a boundary node is wrap-around.
TRACKS_CELL_BOUNDARY = 3

#: Indicates a boundary node is closed
CLOSED_BOUNDARY = 4

BOUNDARY_STATUS_FLAGS_LIST = [
    FIXED_VALUE_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    TRACKS_CELL_BOUNDARY,
    CLOSED_BOUNDARY,
]
BOUNDARY_STATUS_FLAGS = set(BOUNDARY_STATUS_FLAGS_LIST)


class StatusGrid(object):

    def __init__(self, node_status):
        self._node_status = np.array(node_status)

    @property
    def node_status(self):
        """Status of nodes

        Return an array of the status of a grid's nodes. The node status can
        be one of the following:
        - `CORE_NODE`
        - `FIXED_VALUE_BOUNDARY`
        - `FIXED_GRADIENT_BOUNDARY `
        - `TRACKS_CELL_BOUNDARY`
        - `CLOSED_BOUNDARY `
        """
        return self._node_status

    @node_status.setter
    def node_status(self, node_status):
        self._node_status = node_status

    def active_nodes(self):
        """
        Node IDs of all active (core & open boundary) nodes.
        core_nodes will return just core nodes.
        """
        (active_node_ids, ) = np.where(self.node_status != CLOSED_BOUNDARY)
        return active_node_ids

    def core_nodes(self):
        """Node IDs of all core nodes.
        """
        (core_node_ids, ) = np.where(self.node_status == CORE_NODE)
        return core_node_ids

    def boundary_nodes(self):
        """Node IDs of all boundary nodes.
        """
        (boundary_node_ids, ) = np.where(self.node_status != CORE_NODE)
        return boundary_node_ids

    def open_boundary_nodes(self):
        """Node id of all open boundary nodes.
        """
        (open_boundary_node_ids, ) = np.where(
            (self.node_status != CLOSED_BOUNDARY) &
            (self.node_status != CORE_NODE))
        return open_boundary_node_ids

    def closed_boundary_nodes(self):
        """Node id of all closed boundary nodes.
        """
        (closed_boundary_node_ids, ) = np.where(
            self.node_status == CLOSED_BOUNDARY)
        return closed_boundary_node_ids

    def fixed_gradient_boundary_nodes(self):
        """Node id of all fixed gradient boundary nodes
        """
        (fixed_gradient_boundary_node_ids, ) = np.where(
            self.node_status == FIXED_GRADIENT_BOUNDARY)
        return fixed_gradient_boundary_node_ids

    def fixed_value_boundary_nodes(self):
        """Node id of all fixed value boundary nodes
        """
        (fixed_value_boundary_node_ids, ) = np.where(
            self.node_status == FIXED_VALUE_BOUNDARY)
        return fixed_value_boundary_node_ids
