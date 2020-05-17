from landlab import Component
from landlab.data_record import DataRecord
from landlab.utils.parcels import make_sediment


class SyntheticPulseParcelCreator(Component):

    _name = "SyntheticPulseParcelCreator"

    _unit_agnostic = True

    _info = {}

    def __init__(self, grid, parcels, abrasion_rate):
        """Create one or more synthetic pulses that add parcels to a network

        More description here.
        And more description here.

        Note that this is unit agnostic, but that it is designed to work with
        the :py:class:`~landlab.components.NetworkSedimentTransporter` which
        requires mks units.

        Parameters
        ----------
        grid : NetworkModelGrid
        parcels : DataRecord


        information about when parcels come in:
            - float (specific time)
            - array (set of times)
            - some stochastic method (look at the Precipitation Distribution
            component?)

        Add all the parcel attributes that are necessary for make_sediment here.

        abrasion_rate

        Examples
        --------
        >>> from landlab.components import SyntheticPulseParcelCreator

        """
        # initial input type checking.
        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)

        if not isinstance(parcels, DataRecord):
            msg = ""
            raise ValueError(msg)

        # run super. this will check for required inputs specified by _info
        # will also provide the powers of a COmponent.
        super().__init__(grid)

        # set initial time.
        self.current_time = 0.0

        # save all the information about pulse parcel attributes.
        self.abrasion_rate = abrasion_rate

        # set and save the up ability to deterimine if a pulse is occuring or
        # not. this might include some if-else depending on how its specified.

        # set up information that controls the number of parcels necessary.

    # The following is an example of how @property and a setter might be used to
    # support different options for input.
    @property
    def abrasion_rate(self):
        return self._abrasion_rate

    @abrasion_rate.setter
    def abrasion_rate(self, val):
        """Abrasion rate

        Describe the characterstics

        Parameters
        ----------
        val : float or dict
            If dict then must be ... describe.
            Must be greater than zero.

        """
        # check that the value passed is correctly.
        # either a float that is greater than zero, or a dictionary in which #the first element is something in numpy.random.

        if isinstance(val, dict):
            if val["distribution"] not in np.random.__dict__:
                raise ValueError
        else:
            if not isinstance(val, float):
                raise ValueError
            if val < 0:
                raise ValueError

        self._abrasion_rate = val

    # You can make any functions you want that make init and run one step easier.

    def run_one_step(self, dt):
        """Advance the SyntheticPulseParcelCreator forward by ``dt``.

        Parameters
        ----------
        dt : float
            Timestep duration
        """

        # figure out how many parcels we want. (presently with dummy value.)
        number_of_parcels = 100

        # figure out whether a pulse is occuring.(presently with dummy value.)
        pulse_this_timestep = True

        # increment time (sometimes it makes sense to do this now, sometimes at
        # the end of the run one step. )
        self.current_time += dt

        if pulse_this_timestep:

            items = make_sediment(
                self.grid,
                self.current_time,
                abrasion_rate=self.abrasion_rate,
                D=self.D,
                density=self.density,
                volume=self.volume,
                number_of_parcels=number_of_parcels,
            )

            self.parcels.add_items(items)
