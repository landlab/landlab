import re

import pytest

from landlab.grid.warnings import DeprecatedSignature


def test_deprecated_signature():
    with pytest.raises(
        DeprecatedSignature, match="You are using a deprecated calling signature"
    ):
        raise DeprecatedSignature("foobar")


def test_deprecated_signature_with_new():
    with pytest.raises(DeprecatedSignature, match=re.escape("foobar('foo', bar=2)")):
        raise DeprecatedSignature("foobar", new=(("foo",), {"bar": 2}))
