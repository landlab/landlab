#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions to read shapefiles and create a NetworkModelGrid."""
import pathlib

import numpy as np
import shapefile as ps
from shapefile import ShapefileException

from landlab.graph.graph import NetworkGraph
from landlab.grid.network import NetworkModelGrid


def _read_shapefile(file, dbf):
    if isinstance(file, pathlib.PurePath):
        file = str(file)
    try:
        sf = ps.Reader(file)
    except ShapefileException:
        try:
            sf = ps.Reader(shp=file, dbf=dbf)
        except ShapefileException:
            raise ShapefileException(("Bad file path provided to read_shapefile."))
    return sf


def read_shapefile(
    file,
    dbf=None,
    store_polyline_vertices=True,
    points_shapefile=None,
    points_dbf=None,
    link_fields=None,
    node_fields=None,
    link_field_conversion=None,
    node_field_conversion=None,
    link_field_dtype=None,
    node_field_dtype=None,
    threshold=0.0,
):
    """Read shapefile and create a NetworkModelGrid.

    There are a number of assumptions that are requied about the shapefile.
        * The shape file must be a polyline shapefile.
        * All polylines must be their own object (e.g. no multi-part
          polylines).
        * Polyline endpoints match perfectly.

    You might notice that there is no ``write_shapefile`` function. If this is
    something you need for your work, please make a GitHub issue to start this
    process.

    Parameters
    ----------
    file: str or file-like
        File path or file-like of a valid polyline shapefile
    dbf: file-like, optional
        If file is file-like, the dbf must also be passed.
    store_polyline_vertices: bool, optional
        If True (default), store the vertices of the polylines in
        the at_link fields ``x_of_polyline`` and ``y_of_polyline``.
    points_shapefile: str or file-like
        File path or file-like of a valid point shapefile.
    points_dbf: file-like, optional
        If file is file-like, the dbf must also be passed.
    link_fields: list, optional
        List of polyline shapefile attributes to import as landlab at-link
        fields. Default is to import all.
    node_fields: list, optional
        List of point shapefile attributes to import as landlab at-node
        fields. Default is to import all.
    link_field_conversion: dict, optional
        Dictionary mapping polyline shapefile field names to desired at link
        field names. Default is no remapping.
    node_field_conversion: dict, optional
        Dictionary mapping node shapefile field names to desired at node field
        names. Default is no remapping.
    link_field_dtype: dict, optional
        Dictionary mapping node shapefile field names to desired dtype. Default
        is no change to dtype.
    node_field_dtype: dict, optional
        Dictionary mapping node shapefile field names to desired dtype. Default
        is no change to dtype.
    threshold: float, optional
        Maximum distance between a point in the point shapefile and a polyline
        junction in the polyline shapefile. Units are the same as in the
        shapefiles. Default is zero (requiring perfect overlap).

    Returns
    -------
    grid : NetworkModelGrid instance
        The network model grid will have nodes at the endpoints of the
        polylines, and links that connect these nodes. Any fields
        associated with the shapefile will be added as at-link fields. If a
        point shapefile is provided those values will be added as at-node
        fields.

    Examples
    --------
    First, we make a simple shapefile

    >>> from io import BytesIO
    >>> import shapefile
    >>> shp = BytesIO()
    >>> shx = BytesIO()
    >>> dbf = BytesIO()
    >>> w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    >>> w.shapeType = 3
    >>> w.field("spam", "N")
    >>> w.line([[[5,5],[10,10]]])
    >>> w.record(37)
    >>> w.line([[[5,0],[5,5]]])
    >>> w.record(100)
    >>> w.line([[[5,5],[0,10]]])
    >>> w.record(239)
    >>> w.close()

    Now create a NetworkModelGrid with read_shapefile:

    >>> from landlab.io import read_shapefile
    >>> grid = read_shapefile(shp, dbf=dbf)
    >>> grid.nodes
    array([0, 1, 2, 3])
    >>> grid.x_of_node
    array([  5.,   5.,   0.,  10.])
    >>> grid.y_of_node
    array([  0.,   5.,  10.,  10.])
    >>> grid.nodes_at_link
    array([[0, 1],
           [2, 1],
           [1, 3]])
    >>> assert "spam" in grid.at_link
    >>> grid.at_link["spam"]
    array([100, 239, 37])

    Next lets also include a points file. First create both shapefiles.

    >>> shp = BytesIO()
    >>> shx = BytesIO()
    >>> dbf = BytesIO()
    >>> w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)
    >>> w.shapeType = 3
    >>> w.field("spam", "N")
    >>> w.line([[[5,5],[10,10]]])
    >>> w.record(37)
    >>> w.line([[[5,0],[5,5]]])
    >>> w.record(100)
    >>> w.line([[[5,5],[0,10]]])
    >>> w.record(239)
    >>> w.close()

    >>> p_shp = BytesIO()
    >>> p_shx = BytesIO()
    >>> p_dbf = BytesIO()
    >>> p_w = shapefile.Writer(shp=p_shp, shx=p_shx, dbf=p_dbf)
    >>> p_w.shapeType = 1
    >>> p_w.field("eggs", "N")
    >>> p_w.point(5, 0)
    >>> p_w.record(2)
    >>> p_w.point(5, 5)
    >>> p_w.record(4)
    >>> p_w.point(0, 10)
    >>> p_w.record(8)
    >>> p_w.point(10, 10)
    >>> p_w.record(6)
    >>> p_w.close()

    Now read in both files together.

    >>> grid = read_shapefile(shp,dbf=dbf,points_shapefile=p_shp,points_dbf=p_dbf)
    >>> grid.nodes
    array([0, 1, 2, 3])
    >>> grid.x_of_node
    array([  5.,   5.,   0.,  10.])
    >>> grid.y_of_node
    array([  0.,   5.,  10.,  10.])
    >>> grid.nodes_at_link
    array([[0, 1],
           [2, 1],
           [1, 3]])
    >>> assert "spam" in grid.at_link
    >>> grid.at_link["spam"]
    array([100, 239, 37])
    >>> assert "eggs" in grid.at_node
    >>> grid.at_node["eggs"]
    array([2, 4, 8, 6])
    """
    sf = _read_shapefile(file, dbf)

    link_field_conversion = link_field_conversion or dict()
    node_field_conversion = node_field_conversion or dict()
    link_field_dtype = link_field_dtype or dict()
    node_field_dtype = node_field_dtype or dict()

    if sf.shapeType != 3:
        raise ValueError(
            (
                "landlab.io.shapefile read requires a polyline "
                "type shapefile. The provided shapefile does "
                "not meet these requirements."
            )
        )

    if points_shapefile:
        psf = _read_shapefile(points_shapefile, points_dbf)
        if psf.shapeType != 1:
            raise ValueError(
                (
                    "landlab.io.shapefile read requires a point "
                    "type shapefile. The provided shapefile does "
                    "not meet these requirements."
                )
            )

    # get record information, the first element is ('DeletionFlag', 'C', 1, 0)
    # which we will ignore.
    records = sf.fields[1:]

    # initialize data structures for node (x,y) tuples,
    # link (head_node_id, tail_node_id) tuples, and a dictionary of at-link
    # fields.
    # besides the at-link fields on the shapefile, we'll also store an array of
    # x and y of the full polyline segment'.

    node_xy = []
    links = []
    fields = {rec[0]: [] for rec in records}

    # store which link fields to retain
    link_fields_to_retain = link_fields or list(fields.keys())
    if store_polyline_vertices:
        link_fields_to_retain.append("x_of_polyline")
        link_fields_to_retain.append("y_of_polyline")
    if store_polyline_vertices:
        fields["x_of_polyline"] = []
        fields["y_of_polyline"] = []

    record_order = [rec[0] for rec in records]

    # itterate through shapes and records
    shapeRecs = sf.shapeRecords()
    for sr in shapeRecs:

        # if not a multi-part polyline:
        if len(sr.shape.parts) == 1:

            # get all the points on the polyline and deconstruct into x and y
            points = sr.shape.points
            x, y = zip(*points)

            # construct the (x,y) tuples of the head and tail nodes of each
            # polyline. Note here, that head and tail just refer to starting and
            # ending, they will be re-oriented if necessary by landlab.

            head_node_xy = (x[0], y[0])
            tail_node_xy = (x[-1], y[-1])

            # we should expect that the head node and tail node of later links will
            # already be part of the model grid. So we check, and add the nodes,
            # if they don't already exist.

            if head_node_xy not in node_xy:
                node_xy.append(head_node_xy)

            if tail_node_xy not in node_xy:
                node_xy.append(tail_node_xy)

            # get the index of the head and tail node index.
            head_node__node_id = node_xy.index(head_node_xy)
            tail_node__node_id = node_xy.index(tail_node_xy)

            # append the head and tail node ids to the link array
            links.append((head_node__node_id, tail_node__node_id))

            for i in range(len(sr.record)):
                field_name = record_order[i]
                fields[field_name].append(sr.record[i])

            if store_polyline_vertices:
                fields["x_of_polyline"].append(x)
                fields["y_of_polyline"].append(y)

        else:
            raise ValueError(
                (
                    "landlab.io.shapefile currently does not support ",
                    "reading multipart polyline shapefiles.",
                )
            )

    # Create a Network Model Grid.
    x_of_node, y_of_node = zip(*node_xy)

    # We want to ensure that we maintain sorting, so start by creating an
    # unsorted network graph and sorting.
    # The sorting is important to ensure that the fields are assigned to
    # the correct links.
    graph = NetworkGraph((y_of_node, x_of_node), links=links, sort=False)
    sorted_nodes, sorted_links, sorted_patches = graph.sort()

    # use the sorting information to make a new network model grid.
    grid = NetworkModelGrid(
        (np.asarray(y_of_node)[sorted_nodes], np.asarray(x_of_node)[sorted_nodes]),
        np.vstack((graph.node_at_link_head, graph.node_at_link_tail)).T,
    )

    # add values to fields.
    for field_name in fields:
        if field_name in link_fields_to_retain:
            mapped_field_name = link_field_conversion.get(field_name, field_name)
            mapped_dtype = link_field_dtype.get(field_name, None)
            grid.at_link[mapped_field_name] = np.asarray(
                fields[field_name], dtype=mapped_dtype
            )[sorted_links]

    # if a points shapefile is added, bring in and use.
    if points_shapefile:
        # get ready to store fields.
        psf_records = psf.fields[1:]
        psf_record_order = [rec[0] for rec in psf_records]
        psf_fields = {rec[0]: [] for rec in psf_records}

        # store which node fields to retain
        node_fields_to_retain = node_fields or list(psf_fields.keys())

        # we don't need to store node xy, just need to store which index each
        # node maps to on the new grid.
        psf_node_mapping = -1 * np.ones(grid.x_of_node.shape, dtype=int)

        # loop through each node
        psf_shapeRecs = psf.shapeRecords()
        for node_idx, sr in enumerate(psf_shapeRecs):
            # find the closest
            point_x = sr.shape.points[0][0]
            point_y = sr.shape.points[0][1]
            x_diff = grid.x_of_node - point_x
            y_diff = grid.y_of_node - point_y

            dist = np.sqrt(x_diff ** 2 + y_diff ** 2)

            # check that the distance is small.
            if np.min(dist) > threshold:
                msg = (
                    "landlab.io.shapefile: a point in the points shapefile "
                    "is {dist} away from the closet polyline junction in the ".format(
                        dist=np.min(dist)
                    ),
                    "polyline shapefile. This is larger than the threshold"
                    "value of {thresh}. This may mean that the threshold".format(
                        thresh=threshold
                    ),
                    "or that something is wrong with the points file.",
                )
                raise ValueError(msg)

            ind = np.nonzero(dist == np.min(dist))[0]
            # verify that there is only one closest.

            if psf_node_mapping[ind[0]] >= 0:
                msg = (
                    "landlab.io.shapefile requires that the points file "
                    "have a 1-1 mapping to the polylines file. More than one "
                    "at-node point provided maps to the node with Landlab ID "
                    "{ind}, (x,y). This point has coordinates of ({x}, {y})".format(
                        ind=ind[0], x=point_x, y=point_y
                    )
                )
                raise ValueError(msg)

            psf_node_mapping[ind[0]] = node_idx

            for rec_idx in range(len(sr.record)):
                field_name = psf_record_order[rec_idx]
                psf_fields[field_name].append(sr.record[rec_idx])

        if np.any(psf_node_mapping < 0):
            msg = (
                "landlab.io.shapefile requires that the points file "
                "contain the same number of points as polyline junctions. "
                "The points file contains fewer points than polyline junctions."
            )
            raise ValueError(msg)

        # add values to nodes.
        for field_name in psf_fields:
            if field_name in node_fields_to_retain:
                mapped_field_name = node_field_conversion.get(field_name, field_name)
                mapped_dtype = node_field_dtype.get(field_name, None)

                grid.at_node[mapped_field_name] = np.asarray(
                    psf_fields[field_name], dtype=mapped_dtype
                )[psf_node_mapping]

    return grid
