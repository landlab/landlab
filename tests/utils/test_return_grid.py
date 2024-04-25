import numpy as np
import pytest

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.utils.return_array import return_array_at_link
from landlab.utils.return_array import return_array_at_node


def test_no_field():
    mg = RasterModelGrid((10, 10))
    with pytest.raises(FieldError):
        return_array_at_node(mg, "spam")
    with pytest.raises(FieldError):
        return_array_at_link(mg, "spam")


def test_return_array():
    mg = RasterModelGrid((10, 10))
    node_vals = np.arange(mg.number_of_nodes)
    out = return_array_at_node(mg, node_vals)

    np.testing.assert_array_equal(np.arange(mg.number_of_nodes), out)

    link_vals = np.arange(mg.number_of_links)
    out = return_array_at_link(mg, link_vals)

    np.testing.assert_array_equal(np.arange(mg.number_of_links), out)
