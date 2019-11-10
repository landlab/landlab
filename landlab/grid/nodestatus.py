#! /usr/bin/env python
from enum import IntEnum, unique


@unique
class NodeStatus(IntEnum):
    """Define the boundary-type codes"""
    CORE = 0
    FIXED_VALUE = 1
    FIXED_GRADIENT = 2
    LOOPED = 3
    CLOSED = 4


# Define the boundary-type codes

#: Indicates a node is *core*.
CORE_NODE = 0

#: Indicates a boundary node is has a fixed values.
FIXED_VALUE_BOUNDARY = 1

#: Indicates a boundary node is has a fixed gradient.
FIXED_GRADIENT_BOUNDARY = 2

#: Indicates a boundary node is wrap-around.
LOOPED_BOUNDARY = 3

#: Indicates a boundary node is closed
CLOSED_BOUNDARY = 4

BOUNDARY_STATUS_FLAGS_LIST = [
    FIXED_VALUE_BOUNDARY,
    FIXED_GRADIENT_BOUNDARY,
    LOOPED_BOUNDARY,
    CLOSED_BOUNDARY,
]
BOUNDARY_STATUS_FLAGS = set(BOUNDARY_STATUS_FLAGS_LIST)
