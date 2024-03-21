#! /usr/bin/env python
from enum import IntEnum
from enum import unique


@unique
class NodeStatus(IntEnum):
    """Define the boundary-type codes"""

    #: Indicate a node is *core*.
    CORE = 0
    #: Indicate a boundary node is has a fixed values.
    FIXED_VALUE = 1
    #: Indicate a boundary node is has a fixed gradient.
    FIXED_GRADIENT = 2
    #: Indicate a boundary node is wrap-around.
    LOOPED = 3
    #: Indicate a boundary node is closed
    CLOSED = 4
