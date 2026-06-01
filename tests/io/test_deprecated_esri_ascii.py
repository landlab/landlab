import re

import pytest

from landlab.io._deprecated_esri_ascii import DataSizeError
from landlab.io._deprecated_esri_ascii import KeyTypeError
from landlab.io._deprecated_esri_ascii import KeyValueError
from landlab.io._deprecated_esri_ascii import MismatchGridDataSizeError
from landlab.io._deprecated_esri_ascii import MismatchGridXYLowerLeft
from landlab.io._deprecated_esri_ascii import MismatchGridXYSpacing


@pytest.mark.parametrize(
    "err, args, match",
    (
        (DataSizeError, (3, 5), "3 != 5"),
        (KeyTypeError, ("foo", "int"), "unable to convert 'foo' to int"),
        (KeyValueError, ("foo", "bar"), "'foo': bar"),
        (MismatchGridDataSizeError, (5, 4), "(data size) 5 != 4 (grid size)"),
        (
            MismatchGridXYLowerLeft,
            ((0, 1), (1, 2)),
            "(data lower-left) (0, 1) != (1, 2) (grid lower-left)",
        ),
        (MismatchGridXYSpacing, (1, 5), "(data dx) 1 != 5 (grid dx)"),
    ),
)
def test_error_as_string(err, args, match):
    with pytest.raises(err, match=re.escape(match)):
        raise err(*args)
