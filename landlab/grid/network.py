#! /usr/bin/env python

from ..graph import Graph
from ..field import GraphFields


class NetworkModelGrid(Graph, GraphFields):

    """A ModelGrid of just nodes and links.
    
    Parameters
    ----------
    yx_of_node : tuple of ndarray
        Node y and x coordinates.
    links : array of tuple of int
        Nodes at link tail and head.
    
    Examples
    --------
    >>> from landlab.grid.network import NetworkModelGrid
    >>> y_of_node = (0, 1, 2, 2)
    >>> x_of_node = (0, 0, -1, 1)
    >>> nodes_at_link = ((1, 0), (2, 1), (3, 1))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.x_of_node
    array([ 0.,  0., -1.,  1.])
    >>> grid.y_of_node
    array([ 0.,  1.,  2.,  2.])
    >>> grid.nodes_at_link
    array([[0, 1],
           [2, 1],
           [1, 3]])
    """

    def __init__(self, yx_of_node, links):
        Graph.__init__(self, yx_of_node, links=links)
        GraphFields.__init__(self,
                             {'node': self.number_of_nodes,
                              'link': self.number_of_links},
                             default_group='node')
