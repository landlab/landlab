# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 12:06:29 2019

Temporary file for use while writing tests


@author: pfeif
"""
import numpy as np


def _calculate_alluvium_depth(stored_volume,width_of_upstream_links, length_of_upstream_links, width_of_downstream_link, length_of_downstream_link, porosity):
    """Calculate alluvium depth based on adjacent link inactive parcel volumes.

    Parameters
    ----------
    stored_volume : float
        Total volume of inactive parcels in this link.
    width_of_upstream_links : float
        Channel widths of upstream links.
    length_of_upstream_link : float
        Channel lengths of upstream links.
    width_of_downstream_link : float
        Channel widths of downstream links.
    length_of_downstream_link : float
        Channel lengths of downstream links.
    porosity: float
        Channel bed sediment porosity.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_alluvium_depth
    >>> import pytest
    >>> _calculate_alluvium_depth(100,np.array([0.5,1]),np.array([10,10]) 1, 10, 0.2)
    10.0
    >>>_calculate_alluvium_depth(x, x, x, x, x)
    0.0001
    >>> with pytest.raises(ValueError):
    ...     _calculate_alluvium_depth(x, x, x, x, x)

    """
    import numpy as np
    alluvium__depth = (
        2
        * stored_volume
        / (
            np.sum(width_of_upstream_links * length_of_upstream_links)
            + width_of_downstream_link * length_of_downstream_link
        )
        /(1-porosity)
    )
    # NOTE: Jon, porosity was left out in earlier version of the LL component,
    # but it seems it should be in here. Check me: is the eqn correct?

    if alluvium__depth < 0.0:
        raise ValueError("NST Alluvium Depth Negative")

    return alluvium__depth

alluvium__depth = _calculate_alluvium_depth(24,np.array([0.1,3]),np.array([10,10]), 1, 1, 0.5)

# %%

def _calculate_reference_shields_stress(
                        thing
                        ):
    """Calculate reference shields stress (taursg) using the sand content of
    the bed surface, as per Wilcock and Crowe (2003).

    Parameters
    ----------
    fluid_density : float
        Density of fluid (generally, water).
    R: float
        Specific weight..?
    mean_active_grain_size: float
        Mean grain size of the 'active' sediment parcels.
    frac_sand: float
        Fraction of the bed surface grain size composed of sand sized parcels.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter _calculate_reference_shields_stress
    >>> import pytest
    >>> _calculate_reference_shields_stress(x, x, x, x, x)
    1.0
    >>>_calculate_reference_shields_stress(x, x, x, x, x)
    0.0001
    >>> with pytest.raises(ValueError):
    ...     _calculate_reference_shields_stress(x, x, x, x, x)

    """

    taursg = (
        fluid_density
        * R
        * self.g
        * mean_active_grain_size
        * (0.021 + 0.015 * np.exp(-20.0 * frac_sand))
    )

    if taursg < 0.0:
        raise ValueError("NST reference Shields stress is negative")
    if taursg > 0.05:
        raise ValueError("NST reference Shields stress is unreasonably high")
    return taursg




def _calculate_parcel_volume_post_abrasion(
                        starting_volume,
                        travel_distance,
                        abrasion_rate
                        ):
    """Calculate parcel volumes after abrasion, according to classic
    Sternberg exponential abrasion.

    Parameters
    ----------
    starting_volume : float
        Starting volume of each parcel.
    travel_distance: float
        Travel distance for each parcel during this timestep, in ___.
    abrasion_rate: float
        Mean grain size of the 'active' sediment parcels.

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_volume_post_abrasion
    >>> import pytest
    >>> _calculate_parcel_volume_post_abrasion(x, x, x, x, x)
    1.0
    >>>_calculate_parcel_volume_post_abrasion(x, x, x, x, x)
    0.0001
    >>> with pytest.raises(ValueError):
    ...     _calculate_parcel_volume_post_abrasion(x, x, x, x, x)

    """

    volume = starting_volume * np.exp(travel_distance * (-abrasion_rate))

    if volume < 0.0:
        raise ValueError("NST parcel volume is negative")

    return volume


def _calculate_parcel_grain_diameter_post_abrasion(
                        starting_diameter,
                        pre_abrasion_volume,
                        post_abrasion_volume
                        ):
    """Calculate parcel grain diameters after abrasion, according to classic
    Sternberg exponential abrasion.

    Parameters
    ----------
    starting_diameter : float
        Starting volume of each parcel.
    pre_abrasion_volume: float
        xxx
    post_abrasion_volume: float
        xxx

    Examples
    --------
    >>> from landlab.components.network_sediment_transporter.network_sediment_transporter import _calculate_parcel_grain_diameter_post_abrasion
    >>> import pytest
    >>> _calculate_parcel_volume_post_abrasion(x, x, x, x, x)
    1.0
    >>>_calculate_parcel_volume_post_abrasion(x, x, x, x, x)
    0.0001
    >>> with pytest.raises(ValueError):
    ...     _calculate_parcel_volume_post_abrasion(x, x, x, x, x)

    """

    grain_size = (starting_diameter
            * (post_abrasion_volume
            / pre_abrasion_volume
            ) ** (1 / 3)
            )

    return grain_size
