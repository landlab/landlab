"""
Calculate recurrence interval (return period) given a daily flow timeseries
and a desired return period.

.. codeauthor:: Jon Czuba

Created on May 18, 2020

"""

import numpy as np
from scipy import interpolate


def calc_recurrence_interval_flow(flow_cms, return_period_years):
    """
    Calculate recurrence interval (return period) given a daily flow 
    timeseries and a desired return period.

    Parameters
    ----------
    flow_cms : float
        Timeseries of flow discharge, cubic meters per second (cms).
    return_period_years : float
        Return period (recurrence interval) at which to calculate flow, years

    Returns
    -------
    recurrence_interval_flow : float
        Flow at the specified return period (recurrence interval), cms

    """
    
    #verify return period is at least 1 day
    if not 1/365.25 <= return_period_years:
        msg = "Return period must be specified in years and greater than 1 day"
        raise ValueError(msg)
    
    #check length of flow_cms relative to return period
    multiples_of_return_period = np.size(flow_cms)/365/return_period_years
    if multiples_of_return_period < 10:
        msg = "Length of flow record is <10 times the desired return period"
        raise ValueError(msg)


    #sort flow values in descending order
    flow_cms_descend = np.flip(np.sort(flow_cms))

    #only pull non-nan values
    #Qs=Qs(~isnan(Qs));

    number_of_flows = np.size(flow_cms_descend)

    rank = np.arange(1,number_of_flows+1)

    percent_exceedence = ( rank / (number_of_flows + 1) ) * 100

    #recurrence_interval_daily = 1 / (percent_exceedence/100)
    #recurrence_interval_years = recurrence_interval_daily / 365.25

    #plot flow duration curve
    #plt.loglog(percent_exceedence, flow_gage_descend)

    #useful info to understand variables
    # 2-yr flow
    # days in 2 years = 365.25*2 = 730.5
    # fraction exceedance = 1/(365.25*2) = 0.0013689253935660506
    # percent exceedance = 1/(365.25*2)*100 = 0.13689253935660506

    f = interpolate.interp1d(percent_exceedence, np.log10(flow_cms_descend))
    xnew = 1/(365.25*return_period_years)*100
    recurrence_interval_flow = 10 ** f(xnew)   # use interpolation function returned by `interp1d`
    
    return recurrence_interval_flow