
import numpy as np

from six.moves.urllib.request import urlopen
from six import BytesIO, StringIO

from osgeo import gdal
from osgeo import ogr
from osgeo import osr

from landlab.io import read_esri_ascii
from landlab.io.esri_ascii import add_halo

_URL = "http://opentopo.sdsc.edu/otr/getdem"


def read_wgs84_from_from_opentopography(demtype,
                                        west,
                                        south,
                                        east,
                                        north,
                                        **kwargs):
    """
    Read data from OpenTopography Restful Server in WGS84 lat/lon coordinates.

    Parameters
    ----------
    demtype: string
        SRTM GL3 (90m) is "SRTMGL3", SRTM GL1 (30m) is "SRTMGL1",
        SRTM GL1 (Ellipsoidal) is "SRTMGL1_E", ALOS World 3D 30m is "AW3D30",
        ALOS World 3D (Ellipsoidal) is "AW3D30_E".
    west: float
        The western longitude of the bounding box.
    south: float
        The southern latitude of the bounding box.
    east: float
        The eastern longitude of the bounding box.
    north: float
        The northern latitude of the bounding box.

    Optionally, the ``grid``, ``halo``, and ``name`` used by
    **read_esri_ascii** may be provided as keyword arguments.

    Returns
    -------
    grid, z

    Example
    -------
    >>> from landlab.io import read_wgs84_from_from_opentopography
    >>> mg, z = read_wgs84_from_from_opentopography()
    >>> z[0]
    >>> mg.shape()

    """
    # input checking.
    if demtype not in ["SRTMGL3",
                       "SRTMGL1",
                       "SRTMGL1_E",
                       "AW3D30",
                       "AW3D30_E"]:
        raise ValueError("demtype not valid.")

    if west > east:
        raise ValueError("coordinate for west is greater than for east.")

    if south > north:
        raise ValueError("coordinate for south is greater than north")

    url = (_URL + "?" +
           "demtype=" + demtype + "&" +
           "west=" + str(west) + "&" +
           "south=" + str(south) + "&"
           "east=" + str(east) + "&"
           "north=" + str(north) + "&"
           "outputFormat=AAIGrid")

    f = urlopen(url)

    b = BytesIO(f.read())
    file_like = StringIO(b.getvalue().decode("UTF-8"))

    mg, z = read_esri_ascii(file_like, **kwargs)
    # deal with no-data values correctly.

    return mg, z


def read_utm_from_from_opentopography(demtype,
                                      west,
                                      south,
                                      east,
                                      north,
                                      resolution,
                                      zone,
                                      **kwargs):
    """
    Read data from OpenTopography Restful Server in UTM coordinates.


    Parameters
    ----------
    demtype: string
        SRTM GL3 (90m) is "SRTMGL3", SRTM GL1 (30m) is "SRTMGL1",
        SRTM GL1 (Ellipsoidal) is "SRTMGL1_E", ALOS World 3D 30m is "AW3D30",
        ALOS World 3D (Ellipsoidal) is "AW3D30_E".
    west: float
        The western longitude of the bounding box.
    south: float
        The southern latitude of the bounding box.
    east: float
        The eastern longitude of the bounding box.
    north: float
        The northern latitude of the bounding box.
    resolution: float
        The resolution at which to reproject from WGS84 lat/lon into UTM m.
    zone: int
        UTM zone number.
    utm_north: bool, optional
        Flag to indicate if the UTM zone is north (default) or south.

    Optionally, the ``grid``, ``halo``, and ``name`` used by
    **read_esri_ascii** may be provided as keyword arguments.

    Returns
    -------
    grid, z

    Example
    -------
    >>> from landlab.io import read_wgs84_from_from_opentopography
    >>> mg, z = read_wgs84_from_from_opentopography()
    >>> z[0]
    >>> mg.shape()
    """
    # get the correct UTM zone.
    if utm_north is True:
        epsg = 32600 + int(zone)
    else:
        epsg = 32700 + int(zone)

    # lat-lon bounding box
    # based on min/max of gdal.Translate
    source = osr.SpatialReference()
    source.ImportFromEPSG(epsg)

    target = osr.SpatialReference()
    target.ImportFromEPSG(4326)

    transform = osr.CoordinateTransformation(source, target)

    wkt = ("POLYGON ((" +
           str(east) + " " + str(north) + ", " +
           str(east) + " " + str(south) + ", " +
           str(west) + " " + str(south) + ", " +
           str(west) + " " + str(north) + ", " +
           str(east) + " " + str(north) + "))")

    box = ogr.CreateGeometryFromWkt(wkt)
    box.Transform(transform)
    ring = box.GetGeometryRef(0)
    points = ring.GetPoints(0)
    x, y = zip(*points)


    # get wgs84 values from opentopo
    grid, z_wgs84 = read_wgs84_from_from_opentopography(demtype,
                                                        east=np.max(x),
                                                        south=np.min(y),
                                                        west=np.min(x),
                                                        north=np.max(y))

    # reproject to UTM with gdal.Warp (specifying resolution)
    # create in_ds and out_ds

    # or get this directly from the esri ascii objecdt.
    ds = gdal_array.OpenArray(z_wgs84.reshape(grid.shape))
    # indicate the correct bounds.

    gdal.Warp("",
              in_ds,
              format="VRT",
              srcSRS='EPSG:4326', dstSRS='EPSG:'+str(epsg),
              xRes=res, yRes=res,
              outputBounds=(minx, miny, maxx, maxy),
              outputType=gdal.GDT_CFloat64)

    # get "data"
    # data =
    # add halo.
    if "halo" in kwargs:
        halo = kwargs.pop("halo")
        if halo > 0:
            data, shape = add_halo()

    if "grid" in kwargs:
        grid = kwargs.pop("grid")
        if (grid.number_of_node_rows != shape[0]) or (
            grid.number_of_node_columns != shape[1]
        ):
            raise MismatchGridDataSizeError(
                shape[0] * shape[1],
                grid.number_of_node_rows * grid.number_of_node_columns,
            )

    if grid is None:
        grid = RasterModelGrid(shape, spacing=spacing, origin=origin)
    if "name" in kwargs:
        name = kwargs.pop("name")
        grid.add_field("node", name, data)

    return (grid, data)
