import numpy as np

from landlab.data_record import DataRecord
from landlab.grid.network import NetworkModelGrid

class SedimentPulserBase:

    
    """
    This utility is the base class of the SedimentPulserAtLinks and 
    SedimentPulserEachPulse utilities.
    
    The SedimentPulserAtLinks and SedimentPulserEachPulse utilities run the
    landlab DataRecord "add_item" method on a DataRecord configured for the 
    NetworkSedimentTransporter component. 
    

    Parameters
    ----------
    grid : ModelGrid
        landlab *ModelGrid* to place sediment parcels on.
    parcels: landlab DataRecord 
        Tracks parcel location and variables
    D50: float, optional
        median grain size [m]
    D_sd: float, optional
        standard deviation of grain sizes [m]
    rho_sediment : float, optional
        Sediment grain density [kg / m^3].
    parcel_volume : float, optional
        parcel volume used for all parcels that do not have a specified volume
    abrasion_rate: float, optional
        rate that grain size decreases with distance along channel [mm/km?]
    
    
    

    Examples
    --------
    >>> from landlab import NetworkModelGrid
    >>> from landlab.components.network_sediment_transporter.sediment_pulser_base import SedimentPulserBase

    >>> y_of_node = (0, 100, 200, 200, 300, 400, 400, 125)
    >>> x_of_node = (0, 0, 100, -50, -100, 50, -150, -100)
    >>> nodes_at_link = ((1, 0), (2, 1), (1, 7), (3, 1), (3, 4), (4, 5), (4, 6))
    >>> grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    >>> grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    >>> grid.at_link["channel_slope"] = np.full(grid.number_of_links, .01)  # m / m
    >>> grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m
    >>> make_pulse_base = SedimentPulserBase(grid)
    >>> make_pulse.parcels
    None
    
    SedimentPulserBase does not have any methods for adding a pulse
    
    >>> a_pulse = make_pulse_base()
    NotImplementedError: the base component has no call method
    

    """
    def __init__(
        self,
        grid,
        parcels = None,
        D50 = 0.05,
        D_sd = 0.03,
        rho_sediment = 2650.0,
        parcel_volume = 0.5,
        abrasion_rate = 0.0
        ):
        
        self._grid = grid
        self._parcels = parcels
        self._D50 = D50
        self._D_sd = D_sd
        self._rho_sediment = rho_sediment
        self._parcel_volume = parcel_volume
        self._abrasion_rate = abrasion_rate
 
        if not isinstance(grid, NetworkModelGrid):
            msg = "NetworkSedimentTransporter: grid must be NetworkModelGrid"
            raise ValueError(msg)
            
    def __call__(self):
        """__call__ is not implemented for this component."""
        raise NotImplementedError("the base component has no call method")

    def calc_lognormal_distribution_parameters(self, mu_x, sigma_x):
        '''        
        determine mean and standard deviation of the underlying normal distribution
        of a sample that is lognormally distributed following Maidment, 1990, 
        Chapter 18, eq. 18.2.6 
    
        Parameters
        ----------
        mu_x : float
            mean grain size.
        sigma_x : float
            standard deviation of grain sizes.
    
        Returns
        -------
        mu_y : float
            mean of natural log of grain size
        sigma_y : float
            standard deviation of natural log of grain sizes.
    
        '''
        sigma_y = (np.log(((sigma_x**2)/(mu_x**2))+1))**(1/2)
        mu_y = np.log(mu_x)-(sigma_y**2)/2        
    
        
        return mu_y, sigma_y