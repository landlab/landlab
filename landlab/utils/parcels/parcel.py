import numpy as np

from landlab.components.network_sediment_transporter.network_sediment_transporter import (
    _INACTIVE,
)


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

    .. code-block:: python

      {"distribution": "name of distribution",
       #other keyword argument for that distribution in numpy random.
       }

    Right now we are only making the core required inputs for
    the NST. But eventually want to be able to expand to other
    arbitratry inputs (e.g., lithology, osl characteristics,
    etc)


    Parameters
    ----------
    grid : model grid
    time : float
        Time at which sediment is created.
    number_of_parcels

    which_links :
        information about which links parcels are placed on (default is all?)

    abrasion_rate : float or dict (see above)
        must be greater than zero
    density : float
        must be greater than zero.

    location_in_link
    D
    volume

    # talk with AMP and JC about what the right way to specify the order of
    # being added.

    # check with AMP regarding which of these are actually required.
    # is time arrival in link just the current model time?

    # need to think about time arrival at link carefully because time arrival
    # at link is what controls when a parcel is brought into the active layer
    # or not. When active, you move downstream,
    # So you can be deep in the link, and close to the end of the link, and not
    # move for a long time.

    # Have eric help with creating additional attributes, he will have good
    # ideas based on his work with layers.

    Examples
    --------
    >>> from landlab.utils.parcels import make_sediment

    Make one example that uses all default values.

    Make one example that uses all scalars

    Make one example that uses numpy.random

    In unit tests ensure volume is always correct, even with
    wierd distributions.
    Check that thing that are floats are always floats
    Check (same with ints)


    """
    # Part 1:
    # For each required attribute, do something like this to create the
    # required values.

    if isinstance(abrasion_rate, dict):
        distribution = abrasion_rate.pop("distribution")
        function = np.random.__dict__[distribution]
        values = function(size=number_of_parcels, **abrasion_rate)
    else:
        values = abrasion_rate * np.ones(number_of_parcels)

    # Part 2: Some attributes values are pre-constrained.
    # time_arrival_in_link is the time (though need to be careful with FILO)

    # active_layer can be set to _INACTIVE (it will be overwritten in run one step).

    # Part 3: Deal with any additional attributes that are not required:

    # Part 4:
    # given attributes, create the dictionary that data record wants for
    # add items.
    # The following code is pasted directly from Allison's example script.
    # so all the variables/names will need to be created.

    # create new parcel grid elements and element IDs.

    # elements is an array of "link"

    # based on which_links determine which link IDs parcels are placed on.

    # Note that the syntax of how these three variables (new_parcels,
    # new_variables) is very touchy b/c the DataRecord is touchy.
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
