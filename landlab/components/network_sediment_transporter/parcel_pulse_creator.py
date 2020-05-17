from landlab import Component
from landlab.data_record import DataRecord


class SyntheticPulseParcelCreator(Component):
    # metadata goes here.

    def __init__(self, grid, parcels):
        """Create one or more synthetic pulses that add parcels to a network

        More description here.

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

    def run_one_step(self, dt):
        """Advance the SyntheticPulseParcelCreator forward by ``dt``.

        Parameters
        ----------
        dt : float
            Timestep duration
        """
        pass
