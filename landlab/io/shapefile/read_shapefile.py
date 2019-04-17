#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions to read shapefiles and create a NetworkModelGrid.
"""
import shapefile as ps
from shapefile import ShapefileException

from landlab.grid.network import NetworkModelGrid


def read_shapefile(file, dbf=None, store_polyline_vertices=True):
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
    file : str or file-like
        File path or file-like of a valid shapefile
    dbf : file-like, optional
        If file is file-like, the dbf must also be passed.
    store_polyline_vertices: bool, optional
        If True (default), store the vertices of the polylines in
        the at_link fields ``x_of_polyline`` and ``y_of_polyline``.

    Returns
    -------
    grid : NetworkModelGrid instance
        The network model grid will have nodes at the endpoints of the
        polylines, and links that connect these nodes. Any fields
        associated with the shapefile will be added as at-link fields.

    Examples
    --------
    First, we make a simple shapefile

    >>> from six import BytesIO
    >>> import shapefile
    >>> shp = BytesIO()
    >>> shx = BytesIO()
    >>> dbf = BytesIO()

    >>> w = shapefile.Writer(shp=shp, shx=shx, dbf=dbf)

    >>> w.shapeType = 3
    >>> w.field("spam", "N")
    >>> w.line([[[5,0],[5,5]]])
    >>> w.record(100)
    >>> w.line([[[5,5],[0,10]]])
    >>> w.record(239)
    >>> w.line([[[5,5],[10,10]]])
    >>> w.record(37)
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
           [1, 2],
           [1, 3]])
    >>> assert "spam" in grid.at_link
    >>> grid.at_link["spam"]
    array([100, 239, 37])
    """
    try:
        sf = ps.Reader(file)
    except ShapefileException:
        try:
            sf = ps.Reader(shp=file, dbf=dbf)
        except ShapefileException:
            raise ShapefileException(("Bad file path provided to read_shapefile."))

    if sf.shapeType != 3:
        raise ValueError(
            (
                "landlab.io.shapefile read requires a polyline "
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

    # Create a Network Model Grid
    x_of_node, y_of_node = zip(*node_xy)
    grid = NetworkModelGrid((y_of_node, x_of_node), links)
    for field_name in fields:
        grid.at_link[field_name] = fields[field_name]

    return grid
