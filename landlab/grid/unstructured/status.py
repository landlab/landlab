import numpy as np

from ..nodestatus import NodeStatus


class StatusGrid(object):
    def __init__(self, node_status):
        self._node_status = np.array(node_status)

    @property
    def node_status(self):
        """Status of nodes.

        Return an array of the status of a grid's nodes. The node status can
        be one of the following:
        - ``NodeStatus.CORE``
        - ``NodeStatus.FIXED_VALUE``
        - ``NodeStatus.FIXED_GRADIENT``
        - ``NodeStatus.LOOPED``
        - ``NodeStatus.CLOSED``
        """
        return self._node_status

    @node_status.setter
    def node_status(self, node_status):
        self._node_status = node_status

    def active_nodes(self):
        """Node IDs of all active (core & open boundary) nodes.

        core_nodes will return just core nodes.
        """
        (active_node_ids,) = np.where(self.node_status != NodeStatus.CLOSED)
        return active_node_ids

    def core_nodes(self):
        """Node IDs of all core nodes."""
        (core_node_ids,) = np.where(self.node_status == NodeStatus.CORE)
        return core_node_ids

    def boundary_nodes(self):
        """Node IDs of all boundary nodes."""
        (boundary_node_ids,) = np.where(self.node_status != NodeStatus.CORE)
        return boundary_node_ids

    def open_boundary_nodes(self):
        """Node id of all open boundary nodes."""
        (open_boundary_node_ids,) = np.where(
            (self.node_status != NodeStatus.CLOSED)
            & (self.node_status != NodeStatus.CORE)
        )
        return open_boundary_node_ids

    def closed_boundary_nodes(self):
        """Node id of all closed boundary nodes."""
        (closed_boundary_node_ids,) = np.where(self.node_status == NodeStatus.CLOSED)
        return closed_boundary_node_ids

    def fixed_gradient_boundary_nodes(self):
        """Node id of all fixed gradient boundary nodes."""
        (fixed_gradient_boundary_node_ids,) = np.where(
            self.node_status == NodeStatus.FIXED_GRADIENT
        )
        return fixed_gradient_boundary_node_ids

    def fixed_value_boundary_nodes(self):
        """Node id of all fixed value boundary nodes."""
        (fixed_value_boundary_node_ids,) = np.where(
            self.node_status == NodeStatus.FIXED_VALUE
        )
        return fixed_value_boundary_node_ids
