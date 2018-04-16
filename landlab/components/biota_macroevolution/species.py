"""Species BiotaEvolver object.
"""

from landlab.components.biota_macroevolution import BiotaEvolverObject, Zone
import numpy as np
import pandas as pd


class Species(BiotaEvolverObject):
    """A BiotaEvolver species.

    Species contains

    A universally unique identifier (UUID) is assigned to the species at
    initialization. The id is passed to child species.
    """

    subtype = 'base'

    def __init__(self, initial_time, initial_zones,
                 parent_species=-1):
        """Initialize a species.

        Parameters
        ----------
        initial_time : float
            Initial time of the species.
        initial_zones : Zone list
            A list of BiotaEvolver Zone objects of the species at the initial
            time.
        parent_species : BiotaEvolver Species
            The parent species object. An id of -1 indicates no parent species.
        """
        BiotaEvolverObject.__init__(self)

        self.record = pd.DataFrame(columns=['zones'])

        # Set parameters.
        self.parent_species = parent_species

        # Set initial zone(s).
        if isinstance(initial_zones, list):
            z = initial_zones
        else:
            z = [initial_zones]
        self.record.loc[initial_time, 'zones'] = z

    def __str__(self):
        return '<{} at {}>'.format(self.__class__.__name__, hex(id(self)))

    def run_macroevolution_processes(self, time, zone_paths):
        """ Run disperal, speciation, and extinction processes.

        Extinction is not explicitly implemented in this method. The base class
        of species leaves extinction to the disappearance of the range of a
        species.

        Parameters
        ----------
        time : float

        zone_paths : Pandas DataFrame

        Returns
        -------
        surviving_species : BiotaEvolver Species list
            The species that exist after the macroevolution processes run. This
            may include self and/or child species of self, or None if no
            species will persist in `time`.
        """
        surviving_species = []

        # Disperse and speciate.

        for v in zone_paths.itertuples():
            if v.path_type in [Zone.ONE_TO_ONE, Zone.MANY_TO_ONE]:

                self.record.loc[time, 'zones'] = v.destinations
                surviving_species.append(self)

            elif v.path_type in [Zone.ONE_TO_MANY, Zone.MANY_TO_MANY]:

                for d in v.destinations:
                    child_species = Species(time, d, parent_species=self)
                    surviving_species.append(child_species)

        surviving_species = np.array(list(set(surviving_species)))

        return list(surviving_species)

    @property
    def clade(self):
        return self._identifier[0]
