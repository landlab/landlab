import numpy as np


def get_area_of_patch(graph, out=None):
    from .ext.of_patch import calc_area_at_patch

    if out is None:
        out = np.empty(graph.number_of_patches, dtype=float)

    calc_area_at_patch(
        graph.nodes_at_patch,
        np.ascontiguousarray(graph.x_of_node),
        np.ascontiguousarray(graph.y_of_node),
        out,
    )

    return out


def get_centroid_of_patch(graph, out=None):
    from .ext.of_patch import calc_centroid_at_patch

    if out is None:
        out = np.empty((graph.number_of_patches, 2), dtype=float)

    # calc_centroid_at_patch(
    #     graph.nodes_at_patch,
    #     np.ascontiguousarray(graph.x_of_node),
    #     np.ascontiguousarray(graph.y_of_node),
    #     out,
    # )

    links_at_patch = graph.ds["links_at_patch"].values
    calc_centroid_at_patch(
        links_at_patch,
        # graph.links_at_patch,
        np.ascontiguousarray(graph.xy_of_link[:, 0]),
        np.ascontiguousarray(graph.xy_of_link[:, 1]),
        # graph.ds["nodes_at_patch"].values,
        # graph.nodes_at_patch,
        # np.ascontiguousarray(graph.x_of_node),
        # np.ascontiguousarray(graph.y_of_node),
        out,
    )

    return out
