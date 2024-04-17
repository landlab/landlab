from enum import IntEnum
from enum import unique


@unique
class FloodStatus(IntEnum):
    """Codes for depression status"""

    _UNFLOODED = 0
    _PIT = 1
    _CURRENT_LAKE = 2
    _FLOODED = 3
