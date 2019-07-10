import pytest

import numpy as np
from landlab import FieldError, RasterModelGrid
from landlab.utils.return_array import (return_array_at_node,
                                        return_array_at_link)


def test_no_field():
    mg = RasterModelGrid((10, 10))
    with pytest.raises(FieldError):
        return_array_at_node(mg, "spam")
    with pytest.raises(FieldError):
        return_array_at_link(mg, "spam")

def test_return_array():
    mg = RasterModelGrid((10, 10))
    node_vals = np.arange(mg.number_of_nodes)
    return_array_at_node(mg, node_vals)
    link_vals = np.arange(mg.number_of_links)
    return_array_at_link(mg, link_vals)
