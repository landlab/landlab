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
    """
    if shape is not None and array.shape != shape:
        raise ValidationError(
            f"array has incorrect shape: array should be {shape}, got {array.shape}"
        )

    if dtype is not None and array.dtype != np.dtype(dtype):
        raise ValidationError(
            "array has incorrect type: array should by {dtype}, got {array.dtype}"
        )

    if writable is not None and not array.flags.writeable:
        raise ValidationError("array must be writable.")

    if contiguous is not None and not array.flags.c_contiguous:
        raise ValidationError("array must be contiguous.")

    return array
