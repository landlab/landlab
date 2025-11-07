import numpy as np
from numpy.typing import DTypeLike
from numpy.typing import NDArray

from landlab.core.errors import ValidationError


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
        raise ValidationError(
            f"array has incorrect shape: array should be {shape}, got {array.shape}"
        )

    if dtype is not None and array.dtype != np.dtype(dtype):
        raise ValidationError(
            f"array has incorrect type: array should be {dtype}, got {array.dtype}"
        )

    if writable and not array.flags.writeable:
        raise ValidationError("array must be writable.")

    if contiguous and not array.flags.c_contiguous:
        raise ValidationError("array must be contiguous.")

    return array
