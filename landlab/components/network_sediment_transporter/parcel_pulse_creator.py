from landlab import Component
from landlab.data_record import DataRecord
from landlab.utils.parcels import make_sediment

#comment from jeff, try 2
class SyntheticPulseParcelCreator(Component):

    _name = "SyntheticPulseParcelCreator"

    _unit_agnostic = True

    _info = {}

    def __init__(self, grid, parcels):
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


        Examples
        --------
        >>> from landlab.components import SyntheticPulseParcelCreator

        """
        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)

        # run super. this will check for required inputs specified by _info
        super().__init__(grid)
        self.current_time = 0.0

        # save all the attributes we need .

    def run_one_step(self, dt):
        """Advance the SyntheticPulseParcelCreator forward by ``dt``.

        Parameters
        ----------
        dt : float
            Timestep duration
        """

        # pulse_this_timestep = (find out if there is a pulse, true/false)

        self.current_time += dt

        if pulse_this_timestep:

            items = make_sediment(self.grid, self.current_time)

            parcels.add_items(items)
