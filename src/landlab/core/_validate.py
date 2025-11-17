import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import DTypeLike
from numpy.typing import NDArray

from landlab.core.errors import ValidationError


def require_between(
    value: ArrayLike,
    a_min: float | None = None,
    a_max: float | None = None,
    *,
    inclusive_min: bool = True,
    inclusive_max: bool = True,
) -> ArrayLike:
    arr = np.asarray(value)

    if a_min is not None:
        cmp = np.less if inclusive_min else np.less_equal
        op = ">=" if inclusive_min else ">"
        if np.any(cmp(arr, a_min)):
            raise ValidationError(f"value must be {op} {a_min}")

    if a_max is not None:
        cmp = np.greater if inclusive_max else np.greater_equal
        op = "<=" if inclusive_max else "<"
        if np.any(cmp(arr, a_max)):
            raise ValidationError(f"value must be {op} {a_max}")

    return value


def require_positive(value: ArrayLike) -> ArrayLike:
    return require_between(value, a_min=0.0, a_max=None, inclusive_min=False)


def require_nonnegative(value: ArrayLike) -> ArrayLike:
    return require_between(value, a_min=0.0, a_max=None, inclusive_min=True)


def require_negative(value: ArrayLike) -> ArrayLike:
    return require_between(value, a_min=None, a_max=0.0, inclusive_max=False)


def require_nonpositive(value: ArrayLike) -> ArrayLike:
    return require_between(value, a_min=None, a_max=0.0, inclusive_max=True)


def validate_array(
    array: NDArray,
    *,
    dtype: DTypeLike | None = None,
    shape: tuple[int, ...] | None = None,
    writable: bool | None = None,
    contiguous: bool | None = None,
):
    """Validate an array to satisfy requirements.

    Parameters
    ----------
    array : ndarray
        The array to be validated.
    dtype : data-type, optional
        The required data type.
    shape : tuple of int, optional
        The required shape.
    writable : bool, optional
        Require the array to be writable.
    contiguous : bool, optional
        Require the array to be c-contiguous.

    Returns
    -------
    array : ndarray
        The validated array.

    Raises
    ------
    ValidationError
        If the array is invalid.

    Examples
    --------
    >>> import numpy as np
    >>> validate_array(np.array([1, 2, 3, 4]), shape=(2, 2))
    Traceback (most recent call last):
    ...
    landlab.core.errors.ValidationError: incorrect shape: expected (2, 2), got (4,)
    """
    if shape is not None and array.shape != shape:
        raise ValidationError(f"incorrect shape: expected {shape}, got {array.shape}")

    if dtype is not None and array.dtype != np.dtype(dtype):
        raise ValidationError(f"incorrect type: expected {dtype}, got {array.dtype}")

    if writable and not array.flags.writeable:
        raise ValidationError("array is not writable")

    if contiguous and not array.flags.c_contiguous:
        raise ValidationError("array is not contiguous")

    return array
