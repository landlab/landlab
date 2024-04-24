"""Define a graph of nodes-links-patches and its dual.

This class should not be used directly. Instead, it should be used as a
base class when defining other types of graphs.
"""

import inspect
from functools import cached_property

import numpy as np

from ..core.utils import as_id_array
from .graph import Graph
from .graph_convention import ConventionConverter
from .sort.sort import reverse_one_to_one


class DualGraphMeta(type):
    def __init__(cls, name, bases, dct):
        type.__init__(cls, name, bases, dct)

        converter = ConventionConverter("cfc")
        # for name, prop in inspect.getmembers(cls, inspect.isdatadescriptor):
        for name, prop in inspect.getmembers(
            cls, lambda o: inspect.isdatadescriptor(o) or inspect.ismethoddescriptor(o)
        ):
            new_name = converter.conform(name, "nlp")
            if hasattr(cls, new_name):
                continue

            fdoc = inspect.getdoc(prop)
            if fdoc:
                fdoc = inspect.cleandoc(
                    """{}

                    See Also
                    --------
                    :attr:`~{}`
                    """.format(
                        converter.conform(fdoc.splitlines()[0], "nlp"), name
                    )
                )

            setattr(
                cls,
                new_name,
                property(lambda x, name=name: getattr(x._dual, name), None, None, fdoc),
            )


class DualGraph(metaclass=DualGraphMeta):
    @property
    def dual(self):
        return self._dual

    @property
    def node_at_cell(self):
        return self.ds["node_at_cell"].values

    @property
    def nodes_at_face(self):
        return self.ds["nodes_at_face"].values

    @cached_property
    def cell_at_node(self):
        return reverse_one_to_one(self.node_at_cell, minlength=self.number_of_nodes)

    @cached_property
    def link_at_face(self):
        return self._create_link_at_face()

    def _create_link_at_face(self):
        link_at_nodes = {}
        for link, pair in enumerate(self.nodes_at_link):
            # pair.sort()
            link_at_nodes[tuple(np.sort(pair))] = link

        link_at_face = np.full((self.number_of_faces,), -1, dtype=int)
        # for face, pair in enumerate(self._nodes_at_face):
        for face, pair in enumerate(self.nodes_at_face):
            # pair.sort()
            link_at_face[face] = link_at_nodes[tuple(np.sort(pair))]
        self._link_at_face = link_at_face

        return self._link_at_face

    @cached_property
    def face_at_link(self):
        return reverse_one_to_one(self.link_at_face, minlength=self.number_of_links)

    def sort(self):
        from .sort.ext.remap_element import remap_graph_element

        sorted_nodes, sorted_links, sorted_patches = Graph.sort(self)
        sorted_corners, sorted_faces, sorted_cells = self.dual.sort()

        with self.thawed():
            self.node_at_cell[:] = self.node_at_cell[sorted_cells]
            self.nodes_at_face[:] = self.nodes_at_face[sorted_faces]

            remap_graph_element(
                as_id_array(self.node_at_cell), as_id_array(np.argsort(sorted_nodes))
            )
            remap_graph_element(
                as_id_array(self.nodes_at_face).reshape((-1,)),
                as_id_array(np.argsort(sorted_nodes)),
            )
