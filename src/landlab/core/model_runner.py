from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from requireit import require_less_than
from requireit import require_positive


@dataclass(frozen=True, slots=True)
class Clock:
    """Define the time domain and default time step for a model run.

    Parameters
    ----------
    start : float, optional
        Initial model time.
    stop : float, optional
        Final model time. It must be greater than ``start``.
    step : float, optional
        Positive, finite default time-step duration.

    Examples
    --------
    >>> clock = Clock(start=2.0, stop=8.0, step=0.5)
    >>> clock.duration
    6.0
    """

    start: float = 0.0
    stop: float = np.inf
    step: float = 1.0

    def __post_init__(self) -> None:
        require_less_than(self.start, self.stop, name="start")
        require_positive(self.step, name="step")
        if np.isinf(self.step):
            raise ValueError("step must be finite")

    @property
    def duration(self) -> float:
        return self.stop - self.start
