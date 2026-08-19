import math
from collections.abc import Callable
from collections.abc import Iterator
from numbers import Integral

from requireit import require_between
from requireit import require_instance
from requireit import require_positive


def iter_time_steps(duration: float, *, dt: float | None = None) -> Iterator[float]:
    """Yield fixed-size time steps that evenly span a requested duration.

    Split *duration* into equally-sized substeps, so that no substeps
    are longer that *dt*.

    Parameters
    ----------
    duration : float
        Total amount of time to advance.
    dt : float, optional
        Maximum time-step size. If not given, use *duration* as a
        single time step.

    Yields
    ------
    float
        The next time-step size.

    Raises
    ------
    ValueError
        If *duration* is negative, or if *dt* is not finite and
        positive.

    Examples
    --------
    >>> from landlab.core.component_utils import iter_time_steps

    >>> list(iter_time_steps(10.0, dt=2.5))
    [2.5, 2.5, 2.5, 2.5]

    A duration that doesn't divide evenly is split into equal substeps,
    each no longer than *dt*, rather than leaving a short final step.

    >>> list(iter_time_steps(10.0, dt=3.0))
    [2.5, 2.5, 2.5, 2.5]

    If *dt* isn't given, *duration* is used as a single time step.

    >>> list(iter_time_steps(5.0))
    [5.0]
    """
    duration = require_between(
        duration, 0.0, math.inf, inclusive_max=False, name="duration"
    )

    if duration == 0.0:
        return

    dt = duration if dt is None else dt

    dt = require_between(
        dt, 0.0, math.inf, inclusive_min=False, inclusive_max=False, name="dt"
    )

    n_steps = math.ceil(duration / dt)

    step = duration / n_steps
    for _ in range(n_steps):
        yield step


def iter_adaptive_time_steps(
    duration: float,
    *,
    calc_dt: Callable[[], float | None],
    max_steps: int | None = None,
    rtol: float = 1e-12,
) -> Iterator[float]:
    """Yield adaptive time steps that advance up to a requested duration.

    Repeatedly call *calc_dt* to obtain the next stable time-step size,
    capping each step so that the total does not exceed *duration*.
    Note that iteration may stop within the tolerance specified by *rtol*.

    Parameters
    ----------
    duration : float
        Total amount of time to advance.
    calc_dt : callable
        Called with no arguments before each substep to obtain the
        current stable time-step size. A return value of ``None`` signals
        that iteration should be stopped before *duration* is reached.
        A return value of `inf` advances to *duration*.
    max_steps : int, optional
        Maximum number of substeps to yield before raising a
        ``RuntimeError``.
    rtol : float, optional
        Stop once the remaining time is no greater than
        ``rtol * duration``. Consequently, the yielded time steps may
        sum to slightly less than *duration*.

    Yields
    ------
    float
        The next time-step size.

    Raises
    ------
    ValueError
        If *duration* or *rtol* are out of range.
    RuntimeError
        If *max_steps* is exceeded, or if a returned step is either
        invalid or is too small, relative to the elapsed time, to make
        further progress.

    Examples
    --------
    >>> from landlab.core.component_utils import iter_adaptive_time_steps

    >>> steps = iter([2.0, 2.0, 2.0, 1.0])
    >>> list(iter_adaptive_time_steps(7.0, calc_dt=lambda: next(steps)))
    [2.0, 2.0, 2.0, 1.0]

    Return ``None`` from *calc_dt* to stop before *duration* is reached.

    >>> steps = iter([2.0, 2.0, None])
    >>> list(iter_adaptive_time_steps(10.0, calc_dt=lambda: next(steps)))
    [2.0, 2.0]

    >>> list(iter_adaptive_time_steps(10.0, calc_dt=lambda: 3.0))
    [3.0, 3.0, 3.0, 1.0]
    """
    duration = require_between(
        duration, 0.0, math.inf, inclusive_max=False, name="duration"
    )

    rtol = require_between(rtol, 0.0, 1.0, inclusive_max=False, name="rtol")
    if max_steps is not None:
        require_instance(max_steps, Integral, name="max_steps")
        require_positive(max_steps, name="max_steps")

    if duration == 0.0:
        return

    tol = rtol * duration
    elapsed = 0.0
    n_steps = 0

    while duration - elapsed > tol:
        if max_steps is not None and n_steps >= max_steps:
            raise RuntimeError(f"unable to reach {duration!r} in {max_steps} substeps")

        step = calc_dt()
        if step is None:
            break
        if step <= 0.0 or math.isnan(step):
            raise RuntimeError("step must be positive or None")

        this_dt = min(step, duration - elapsed)

        new_elapsed = elapsed + this_dt
        if new_elapsed == elapsed:
            raise RuntimeError(
                "time step is too small relative to the elapsed time"
                " to make further progress due to floating-point precision"
            )

        yield this_dt

        elapsed = new_elapsed
        n_steps += 1
