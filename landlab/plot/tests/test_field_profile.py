from numpy import append, array
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.plot import FieldProfiler


def test_field_profile():
    mg = RasterModelGrid((6, 6), 10)

    z = array([
        -9999., -9999., -9999., -9999., -9999., -9999.,
        -9999.,    28.,     0.,    25.,    28., -9999.,
        -9999.,    30.,     3.,    14.,    30., -9999.,
        -9999.,    32.,    11.,    25.,    32., -9999.,
        -9999.,    34.,    32.,    34.,    36., -9999.,
        -9999., -9999., -9999., -9999., -9999., -9999.])

    field = 'topographic__elevation'

    mg.at_node[field] = z
    mg.set_watershed_boundary_condition_outlet_id(2, z, nodata_value=-9999.)

    # values of row 2
    cn = mg.core_nodes
    values = mg.at_node[field][cn][mg.y_of_node[cn] == 2 * mg.dy]

    # profile of row 2 using coordinates
    fp = FieldProfiler(mg, field, [(10.7, 19.3), (38.0, 18.4)])
    assert_equal(fp.field_value, values)

    # profile of row 2 using nodes
    fp = FieldProfiler(mg, field, [13, 16])
    assert_equal(fp.field_value, values)

    # multisegment profile
    new_node = 19
    values = append(mg.at_node[field][new_node], values)
    fp = FieldProfiler(mg, field, [new_node, 13, 16])
    assert_equal(fp.field_value, values)

    # profile of row 2 using nodes
    values = mg.at_node[field][cn][mg.y_of_node[cn] == 2 * mg.dy]
    values = append(values, values[-1])
    fp = FieldProfiler(mg, field, [13, 16], sample_spacing = 9)
    assert_equal(fp.field_value, values)
