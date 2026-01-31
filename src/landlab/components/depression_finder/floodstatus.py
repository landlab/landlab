from enum import IntEnum
from enum import unique


@unique
class FloodStatus(IntEnum):
    """Codes for depression status"""

    UNFLOODED = 0
    PIT = 1
    CURRENT_LAKE = 2
    FLOODED = 3
