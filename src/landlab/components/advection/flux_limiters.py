#!/usr/bin/env python3
"""Flux limiter functions for advection solver.

There are many flux-limiter functions. For second-order TVD
schemes, there is an envelope of acceptable values of ``phi(r)``::

    r <= phi(r) <= 2r, (0 <= r <= 1)
    phi(1) = 1
    1 <= phi(r) <= r, (1 <= r <= 2)
    1 <= phi(r) <= 2, (r > 2)
"""

import numpy as np


def flux_lim_vanleer(r):
    """Apply Van Leer flux-limiter function."""
    return (r + np.abs(r)) / (1.0 + np.abs(r))
