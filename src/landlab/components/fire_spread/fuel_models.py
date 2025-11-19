from __future__ import annotations

from typing import TypedDict


class FuelModelParams(TypedDict):
    name: str
    w0: float  # fuel loading (tons/acre)
    delta: float  # fuel bed depth (m)
    sigma: float  # surface-area-to-volume ratio (1/ft → converted later)
    h: float  # heat of combustion (BTU/lb)
    tau: float  # time constant (s)


ANDERSON_13: dict[int, FuelModelParams] = {
    1: {
        "name": "Short grass (1 ft)",
        "w0": 0.034,
        "delta": 0.30,
        "sigma": 3500,
        "h": 8000,
        "tau": 900,
    },
    2: {
        "name": "Timber grass & understory",
        "w0": 0.20,
        "delta": 1.0,
        "sigma": 3000,
        "h": 8000,
        "tau": 1800,
    },
    3: {
        "name": "Tall grass (2.5 ft)",
        "w0": 0.138,
        "delta": 0.76,
        "sigma": 3500,
        "h": 8000,
        "tau": 1800,
    },
    4: {
        "name": "Chaparral (6 ft)",
        "w0": 1.0,
        "delta": 2.0,
        "sigma": 2000,
        "h": 8000,
        "tau": 3600,
    },
    5: {
        "name": "Brush (2 ft)",
        "w0": 0.30,
        "delta": 0.6,
        "sigma": 2000,
        "h": 8000,
        "tau": 1800,
    },
    6: {
        "name": "Dormant brush, hardwood slash",
        "w0": 0.60,
        "delta": 0.8,
        "sigma": 2500,
        "h": 8000,
        "tau": 3600,
    },
    7: {
        "name": "Southern rough",
        "w0": 0.40,
        "delta": 0.8,
        "sigma": 2500,
        "h": 8000,
        "tau": 3600,
    },
    8: {
        "name": "Closed timber litter",
        "w0": 0.20,
        "delta": 0.3,
        "sigma": 2000,
        "h": 8000,
        "tau": 1800,
    },
    9: {
        "name": "Hardwood litter",
        "w0": 0.25,
        "delta": 0.3,
        "sigma": 2500,
        "h": 8000,
        "tau": 1800,
    },
    10: {
        "name": "Timber (litter + understory)",
        "w0": 0.60,
        "delta": 1.0,
        "sigma": 2000,
        "h": 8000,
        "tau": 3600,
    },
    11: {
        "name": "Light logging slash",
        "w0": 0.60,
        "delta": 1.0,
        "sigma": 1500,
        "h": 8000,
        "tau": 3600,
    },
    12: {
        "name": "Medium logging slash",
        "w0": 1.2,
        "delta": 1.5,
        "sigma": 1500,
        "h": 8000,
        "tau": 3600,
    },
    13: {
        "name": "Heavy logging slash",
        "w0": 2.0,
        "delta": 2.0,
        "sigma": 1500,
        "h": 8000,
        "tau": 3600,
    },
}
