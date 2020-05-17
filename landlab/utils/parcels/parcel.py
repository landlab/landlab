import numpy as np


def make_sediment(
    grid,
    time,
    time_arrival_in_link,
    abrasion_rate,
    density,
    active_layer,
    location_in_link,
    D,
    volume,
    number_of_parcels,
):
    """One line description.

        More info goes here.

        Note that this is unit agnostic, but that it is designed to work with
        the :py:class:`~landlab.components.NetworkSedimentTransporter` which
        requires mks units.


        When input is a dictionary it has the form:

        {"distribution": "name of distribution",
        other keyword argument for that distribution in numpy random.}

        Parameters
        ----------
        grid : model grid
        time : float
            Time at which sediment is created.
        number_of_parcels

        abrasion_rate : float or dict (see above)
            must be greater than zero
        density : float
            must be greater than zero.



        location_in_link
        D
        volume

        time_arrival_in_link : float or dict
            must be greater than zero

        active_layer :
            Double check that a user specifying this is meaningful.
            Or is it overwritten by NST run one step.

        Examples
        --------
        >>> from landlab.utils.parcels import make_sediment

        """

    if isinstance(abrasion_rate, dict):
        distribution = abrasion_rate.pop("distribution")
        function = np.random.__dict__[distribution]
        values = function(size=number_of_parcels, **abrasion_rate)
    else:
        values = abrasion_rate * np.ones(number_of_parcels)

    # given attributes, create the dictionary that data record wants for
    # add items.

    new_parcels = {
        "grid_element": newpar_grid_elements,
        "element_id": newpar_element_id,
    }

    new_variables = {
        "starting_link": (["item_id"], new_starting_link),
        "abrasion_rate": (["item_id"], new_abrasion_rate),
        "density": (["item_id"], new_density),
        "lithology": (["item_id"], new_lithology),
        "time_arrival_in_link": (["item_id", "time"], new_time_arrival_in_link),
        "active_layer": (["item_id", "time"], new_active_layer),
        "location_in_link": (["item_id", "time"], new_location_in_link),
        "D": (["item_id", "time"], new_D),
        "volume": (["item_id", "time"], new_volume),
    }

    items = {"time": [time], "new_item": new_parcels, "new_item_spec": new_variables}
    return items
