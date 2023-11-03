from enum import IntFlag


class LinkOrientation(IntFlag):
    """Define the link orientations."""

    E = 1
    ENE = 2
    NNE = 4
    N = 8
    NNW = 16
    ESE = 32
