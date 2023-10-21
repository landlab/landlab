#! /usr/bin/env python
"""
Source Tracking Algorithm
+++++++++++++++++++++++++
.. autosummary::

    ~landlab.utils.source_tracking_algorithm.convert_arc_flow_directions_to_landlab_node_ids
    ~landlab.utils.source_tracking_algorithm.track_source
    ~landlab.utils.source_tracking_algorithm.find_unique_upstream_hsd_ids_and_fractions

Authors: Sai Nudurupati & Erkan Istanbulluoglu

Ref 1: 'The Landlab LandslideProbability Component User Manual' @
https://github.com/RondaStrauch/pub_strauch_etal_esurf/blob/master/LandslideComponentUsersManual.pdf

+----------+-------------------------------------------------------------------+
| Notation | Definition                                                        |
+==========+===================================================================+
| MD       | Modeling Domain - Raster grid that is being analyzed/worked upon. |
+----------+-------------------------------------------------------------------+
+ HSD      | Hydrologic Source Domain - Grid that is at least as coarse as MD. |
|          | For more info, refer Ref 1                                        |
+----------+-------------------------------------------------------------------+

"""
import copy
from collections import Counter

import numpy as np


def convert_arc_flow_directions_to_landlab_node_ids(grid, flow_dir_arc):
    """Convert Arc flow_directions to RasterModelGrid node ids.

    This function receives flow directions (D8) from ESRI ArcGIS and converts
    them to Landlab's RasterModelGrid node id. ESRI ArcGIS D8 flow directions
    are either of the eight valid output directions relating to the eight
    adjacent cells into which flow could travel. The valid output directions
    are powers of 2 starting from 2^0 (1) in the Eastern neighbor going
    clockwise to 2^7 (128) at Northeastern neighbor. For more information
    refer 'https://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/
    how-flow-direction-works.htm'

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    flow_dir_arc: ndarray of int, shape (n_nodes, )
        flow directions derived from ESRII ArcGIS.

    Returns
    -------
    receiver_nodes: ndarray of int, shape (n_nodes, )
        downstream node at each node. Note that this array gives the
        receiver nodes only for the core nodes. For non-core
        nodes, a zero is used.
    """
    r_arc_raw = np.log2(flow_dir_arc)
    r_arc_raw = r_arc_raw.astype("int")
    neigh_ = grid.adjacent_nodes_at_node
    diag_ = grid.diagonals_at_node
    neigh_ = np.fliplr(neigh_)
    diag_ = np.fliplr(diag_)
    a_n = np.hsplit(neigh_, 4)
    a_d = np.hsplit(diag_, 4)
    neighbors = np.hstack(
        (a_n[-1], a_d[0], a_n[0], a_d[1], a_n[1], a_d[2], a_n[2], a_d[3])
    )
    # Now neighbors has node ids of neighboring nodes in cw order starting at
    # right, hence the order of neighbors = [r, br, b, bl, l, tl, t, tr]
    receiver_nodes = np.zeros(grid.number_of_nodes, dtype=int)
    receiver_nodes[grid.core_nodes] = np.choose(
        r_arc_raw[grid.core_nodes], np.transpose(neighbors[grid.core_nodes])
    )
    return receiver_nodes


# %%
# Source Routing Algorithm
# Note 1: This algorithm works on core nodes only because core nodes
# have neighbors that are real values and not -1s.
# Note 2: Nodes in the following comments in this section refer to core nodes.
def track_source(grid, hsd_ids, flow_directions=None):
    """Track all contributing upstream core nodes for each core node.

    This algorithm traverses the grid based on information of flow directions
    at nodes and at every node identifies all the nodes upstream of a given
    node. The algorithm creates a dictionary with an entry for each node;
    a node's entry in the dictionary will contain a list with the node_ids
    of all upstream nodes. Thus this method permits identification of the
    source area contributing to each and every node in the model grid. This
    function is different from a standard flow accumulation routine in that
    it not only calculates the amount of flow at each node, but records the
    IDs of all upstream nodes. However, similar to a standard
    flow accumulation routine, it produces an at_node array of the amount
    of flow passing through the node. It also differs from a standard
    flow accumulation routing in that it permits the mapping of flow inputs
    from a coarser grid to to a finer model grid.

    In its present implementation, the algorithm has not been optimized
    for efficient time use. Its methods are brute force and it should be
    expected to be time intensive. It is not recommended to be run frequently
    in a modeling exercise. Due to its intensive nature, this algorithm may
    fail with large watersheds (a present, the development team has not
    derived a maximum stable watershed size).

    This function was initially developed to find contributing area of a
    30 m grid (MD), where the quantitative data that we were interested in was
    available in significantly coarser resolution (called Hydrologic Source
    Domain (HSD)). Therefore, we started working with re-sampled HSD,
    that is at the same resolution as MD, and represents exactly the same
    landscape. Alternatively, one can use the node ids of MD
    (grid.nodes.flatten()) as input for hsd_ids.

    For more information, refer Ref 1.

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    hsd_ids: ndarray of int, shape (n_nodes, )
        array that maps the nodes of the grid to, possibly coarser,
        Hydrologic Source Domain (HSD) grid ids.
    flow_directions: ndarray of int, shape (n_nodes, ), optional.
        downstream node at each node. Alternatively, this data can be
        provided as a nodal field 'flow__receiver_node' on the grid.

    Returns
    -------
    (hsd_upstr, flow_accum): (dictionary, ndarray of shape (n_nodes))
        'hsd_upstr' maps each grid node to corresponding
        contributing upstream hsd_ids. hsd_upstr.keys() will return
        node_ids of the grid. hsd_upstr.values() will return lists of
        all upstream contributing hsd_ids, including repitions of hsd_ids,
        at corresponding node_ids.
        'flow_accum' is an array of the number of upstream contributing
        nodes at each node.
    """
    if flow_directions is None:
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that the source tracking utility is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        r = grid.at_node["flow__receiver_node"]
    else:
        r = flow_directions
    z = grid.at_node["topographic__elevation"]
    core_nodes = grid.core_nodes
    core_elev = z[core_nodes]
    # Sort all nodes in the descending order of elevation
    sor_z = core_nodes[np.argsort(core_elev, kind="stable")[::-1]]
    # Create a list to record all nodes that have been visited
    # To store nodes that have already been counted
    alr_counted = []
    flow_accum = np.zeros(grid.number_of_nodes, dtype=int)
    hsd_upstr = {}
    # Loop through all nodes
    for i in sor_z:
        # Check 1: Check if this node has been visited earlier. If yes,
        # then skip to next node
        if i in alr_counted:
            continue
        # Check 2: If the visited node is a sink
        if r[i] == i:
            hsd_upstr.update({i: [hsd_ids[i]]})
            flow_accum[i] += 1.0
            alr_counted.append(i)
            continue
        # Check 3: Now, if the node is not a sink and hasn't been visited, it
        # belongs to a stream segment. Hence, all the nodes in the stream will
        # have to betraversed.
        # stream_buffer is a list that will hold the upstream contributing
        # node information for that particular segment until reaching outlet.
        stream_buffer = []
        j = i
        switch_i = True
        a = 0.0
        # Loop below will traverse the segment of the stream until an outlet
        # is reached.
        while True:
            # Following if loop is to execute the contents once the first node
            # in the segment is visited.
            if not switch_i:
                j = r[j]
                if j not in core_nodes:
                    break
            # If this node is being visited for the first time,
            # this 'if statement' will executed.
            if flow_accum[j] == 0.0:
                a += 1.0
                alr_counted.append(j)
                stream_buffer.append(hsd_ids[j])
            # Update number of upstream nodes.
            flow_accum[j] += a
            # If the node is being visited for the first time, the dictionary
            # 'hsd_upstr' will be updated.
            if j in hsd_upstr:
                hsd_upstr[j] += copy.copy(stream_buffer)
            # If the node has been already visited, then the upstream segment
            # that was not accounted for in the main stem, would be added to
            # all downstream nodes, one by one, until the outlet is reached.
            else:
                hsd_upstr.update({j: copy.copy(stream_buffer)})
            # If the outlet is reached, the 'while' loop will be exited.
            if r[j] == j:
                break
            # This will be executed only for the first node of the
            # stream segment.
            if switch_i:
                switch_i = False
    return (hsd_upstr, flow_accum)


# %%
# Algorithm to calculate coefficients of each upstream HSD ID
def find_unique_upstream_hsd_ids_and_fractions(hsd_upstr):
    """Finds unique entries in hsd_upstr.values()

    This function operates on hsd_upstr.values(), that are lists of hsd_ids.
    Two new Python dictionaries, 'unique_ids' and 'fractions' are created.

    unique_ids.keys() = hsd_upstr.keys()
    unique_ids.values()[i] = list of unique entries in hsd_upstr.values()[i]

    fractions.keys() = hsd_upstr.keys()
    fractions.values()[i] = (number of entries of each unique_id.values()[i]/
    length of hsd_upstr.values()[i]) for each unique_id.values()[i] in the
    same order.

    Note that 'hsd_upstr' is the output of track_source(). You can use
    an alternative input. In that case, please refer to the documentation
    of track_source() or refer source_tracking_algorithm_user_manual for
    more information.

    Parameters
    ----------
    hsd_upstr: dictionary
        'hsd_upstr' maps each MD grid node to corresponding
        contributing upstream HSD ids.

    Returns
    -------
    (unique_ids, fractions): (dictionary, dictionary)
        Tuple of data. 'unique_ids' maps each MD node with all upstream HSD
        ids without repitition. 'fractions' maps each MD node with the
        fractions of contributions of the corresponding upstream HSD ids in
        the same order as uniques_ids[node_id].
    """
    unique_ids = {}  # Holds unique upstream HSD ids
    C = {}  # Holds corresponding total numbers
    fractions = {}  # Holds corresponding fractions of contribution
    for ke in hsd_upstr.keys():
        cnt = Counter()
        for num in hsd_upstr[ke]:
            cnt[num] += 1
        unique_ids.update({ke: cnt.keys()})
        buf = []
        for k in cnt.keys():
            buf.append(cnt[k])
        C.update({ke: buf})
        e = [s / float(sum(buf)) for s in buf]
        fractions.update({ke: e})
    return (unique_ids, fractions)
