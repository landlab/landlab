from landlab import NetworkModelGrid
from landlab.utils.parcels import SedimentPulser
import numpy as np
import scipy.stats

def test_pulser_defaults_equal():
    y_of_node = (0, 100, 200)
    x_of_node = (0, 0, 100)
    nodes_at_link = ((1, 0), (2, 1), (1, 2))
    #create network model grid
    grid = NetworkModelGrid((y_of_node, x_of_node), nodes_at_link)
    grid.at_link["channel_width"] = np.full(grid.number_of_links, 1.0)  # m
    grid.at_link["channel_slope"] = np.full(grid.number_of_links, .01)  # m / m
    grid.at_link["reach_length"] = np.full(grid.number_of_links, 100.0)  # m
    
    
    # instantiate pulse
    make_pulse = SedimentPulser(grid)
    
    # obtain values for all defaults of initializer
    default_rho_sediment = make_pulse._rho_sediment
    default_std_dev = make_pulse._std_dev
    default_time_to_pulse = make_pulse._time_to_pulse
    #assign d50 of lognormal sediment distribution 
    d50 = 1.5
    
    # "run model" for single time step (will pulse once)
    for time in range(1):
        # call pulser
        parcels = make_pulse(d50, time)
    
    # obtain values for all defaults of call
    # can we access these somehow from 'make_pulse' so they aren't hard-coded?
    default_n_parcels = 100 
    default_link_to_pulse = 0
    default_abrasion_rate = 0.0
    default_time_arrival = 0.0 #for this case we only have 1 time step 
    
    # check values-----------
    np.testing.assert_equal(
        parcels.dataset.element_id.values,
        np.expand_dims(np.full(n_parcels, default_link_to_pulse),1)
    )
    
    np.testing.assert_equal(
        parcels.dataset.starting_link.values,
        np.expand_dims(np.full(n_parcels, default_link_to_pulse),1)
    )
    
    # test 
    np.testing.assert_equal(
        parcels.dataset.abrasion_rate.values,
        np.expand_dims(np.full(n_parcels, default_abrasion_rate),1)
    )
    
    # test starting location
    np.testing.assert_equal(
        parcels.dataset.time_arrival_in_link.values,
        np.expand_dims(np.full(n_parcels, default_time_arrival),1)
    )
    
    # test D distribution
    ks_test_stat, pval = scipy.stats.kstest(
        np.ravel(parcels.dataset.D),
        "lognorm",
        scipy.stats.lognorm.fit(np.ravel(parcels.dataset.D))
    )
    
    if pval < 0.05:
        raise ValueError('D distribution not lognormal')
        # print('D distribution not lognormal')
        
    # test location in link
    
    # test volumes








    
