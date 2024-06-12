import numpy as np

from landlab.graph.structured_quad.structured_quad import (
    StructuredQuadLayoutCython,
    StructuredQuadLayoutPython,
)


class TimeQuad:
    param_names = ["implementation", "number-of-rows"]
    params = [[], []]

    implementation = {
        "cython": StructuredQuadLayoutCython,
        "python": StructuredQuadLayoutPython,
    }

    def time_layout(self, implementation, size):
        getattr(self.implementation[implementation], self.method)((size, size))


class TimeQuadMethod:
    params = [["cython", "python"], list(2 ** np.arange(7, 11))]


class TimeQuadLinksAtPatch(TimeQuadMethod, TimeQuad):
    method = "links_at_patch"


class TimeQuadNodesAtLink(TimeQuadMethod, TimeQuad):
    method = "nodes_at_link"


class TimeQuadHorizontalLinks(TimeQuadMethod, TimeQuad):
    method = "horizontal_links"


class TimeQuadVerticalLinks(TimeQuadMethod, TimeQuad):
    method = "vertical_links"


class TimeQuadPerimeterNodes(TimeQuadMethod, TimeQuad):
    method = "perimeter_nodes"


class TimeQuadLinksAtNode(TimeQuadMethod, TimeQuad):
    method = "links_at_node"


class TimeQuadPatchesAtLink(TimeQuadMethod, TimeQuad):
    method = "patches_at_link"


class TimeQuadLinkDirsAtNode(TimeQuadMethod, TimeQuad):
    method = "link_dirs_at_node"


class TimeQuadPatchesAtNode(TimeQuadMethod, TimeQuad):
    method = "patches_at_node"
