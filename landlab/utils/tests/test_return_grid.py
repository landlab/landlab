import pytest

from landlab import FieldError, RasterModelGrid
from landlab.utils.return_array import return_array_at_node


def test_no_field():
    mg = RasterModelGrid((10, 10))
    with pytest.raises(FieldError):
        return_array_at_node(mg, "spam")
