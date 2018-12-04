
from six.moves.urllib.request import urlopen
from six import BytesIO, StringIO

from landlab.io import read_esri_ascii

_URL = 'http://opentopo.sdsc.edu/otr/getdem'


def get_wgs84_from_from_opentopography(demtype, west, south, east, north, **kwargs):
    """


    Parameters
    ----------
    demtype: string
        SRTM GL3 (90m) is 'SRTMGL3', SRTM GL1 (30m) is 'SRTMGL1',
        SRTM GL1 (Ellipsoidal) is 'SRTMGL1_E', ALOS World 3D 30m is 'AW3D30',
        ALOS World 3D (Ellipsoidal) is 'AW3D30_E'.

    west: float
        The bounding box of data to pull in the order
        (west, south, east, north)
        Where each coordinate is in WGS84 Longitude.
    south: float
    east: float
    north: float

    Returns
    -------
    grid, z

    Example
    -------

    """
    # input checking.
    if demtype not in ['SRTMGL3', 'SRTMGL1', 'SRTMGL1_E', 'AW3D30', 'AW3D30_E']:
        raise ValueError('')

    if west > east:
        raise ValueError('')

    if south > north:
        raise ValueError('')

    url = (_URL + '?' +
           'demtype=' + demtype + '&' +
           'west=' + str(west) + '&' +
           'south=' + str(south) + '&'
           'east=' + str(east) + '&'
           'north=' + str(north) + '&'
           'outputFormat=AAIGrid')

    f = urlopen(url)

    b = BytesIO(f.read())
    file_like = StringIO(b.getvalue().decode('UTF-8'))

    mg, z = read_esri_ascii(file_like, **kwargs)
    # deal with no-data values correctly.

    return mg, z


def interpolate_to_utm_grid(grid, field, x, y, z):
    """
    """
    pass
