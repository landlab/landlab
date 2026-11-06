#!/usr/bin/env python3
"""
Created on Fri Oct 23 20:00:53 2026

@author: matheuswanderleydealmeida
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class FireEvent:
    """A single fire event.

    Attributes
    ----------
    year : float
        Model time at which the fire occurred.
    xy_of_ignition : tuple of float
        Coordinates of the ignition point, in grid units.
    burned_nodes : tuple of int
        IDs of the nodes that burned.
    fire_size : float
        Burned area, in grid units squared.
    severity_factor : float
        Fraction of vegetation removed at each burned node.
    aridity : float
        Aridity index at the time of the fire.
    vegetation_removed_fraction : float
        Mean fraction of vegetation removed across the burned nodes that
        had vegetation before the fire.

    Examples
    --------
    >>> event = FireEvent(
    ...     year=1.0,
    ...     xy_of_ignition=(2.0, 3.0),
    ...     burned_nodes=(12, 13),
    ...     fire_size=2.0,
    ...     severity_factor=0.64,
    ...     aridity=0.5,
    ...     vegetation_removed_fraction=0.64,
    ... )
    >>> event.burned_nodes
    (12, 13)
    """

    year: float
    xy_of_ignition: tuple[float, float]
    burned_nodes: tuple[int, ...]
    fire_size: float
    severity_factor: float
    aridity: float
    vegetation_removed_fraction: float


class FireEventLog:
    """Collect :class:`FireEvent` records in memory.

    An instance is callable, so it can be passed directly as the ``log``
    argument of :class:`WildfireGenerator`. Any other callable that accepts
    a :class:`FireEvent` (for example, one that writes to a CSV file) can be
    used instead.

    Examples
    --------
    >>> log = FireEventLog()
    >>> len(log)
    0
    >>> log(FireEvent(1.0, (2.0, 3.0), (12,), 1.0, 0.64, 0.5, 0.64))
    >>> len(log)
    1
    >>> log.records[0].xy_of_ignition
    (2.0, 3.0)
    """

    def __init__(self):
        self._records = []

    def __call__(self, event: FireEvent):
        self._records.append(event)

    def __len__(self):
        return len(self._records)

    def __iter__(self):
        return iter(self._records)

    @property
    def records(self):
        """All recorded fire events, oldest first."""
        return self._records
