import pytest

from landlab.field.errors import GroupSizeError


def test_error_as_str():
    with pytest.raises(GroupSizeError, match="number of node elements has changed"):
        raise GroupSizeError("node", 4, 5)
