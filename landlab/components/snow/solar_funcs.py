"""
This file defines a set of functions for computing both shortwave
and longwave radiation.  Most of the functions are from Appendix
E of Dingman (2002).  The equation for optical air mass is from
Kasten and Young (1989).

Dingman, S.L. (2002) Physical Hydrology, Appendix E: Radiation
  on Sloping Surfaces, book, 2nd edition.
Kasten and Young (1989) Optical air mass equation.
Whitman, A.M. (2003) A Simple Expression for the Equation of Time,
  an online document at:
  http://www.sunspot.noao.edu/sunspot/pr/answerbook/expl-5.html
Marks and Dozier (1992) WRR paper
Boehner (2009) Ch.8 in Geomorphometry book
See also papers by:  
   Glen Liston and David Tarboton's (w/ C. Luce)
"""
#-----------------------------------------------------------------------
#
#  Copyright (c) 2005-2023, Scott D. Peckham
#
#  Sep 2023.  Added Day_Length() function to match Day_Length_Slope().
#             Added Solar_Elevation_Angle() for completeness.
#             Added Saturation_Vapor_Pressure() & Vapor_Pressure();
#               these were in met_base.py but not here.
#  Aug 2023.  Fixed bug in Optical_Air_Mass() function.
#             Updated value in Solar_Constant() function.
#  Jul 2021.  Added year keyword to True_Solar_Noon(), which is
#             then passed to Equation_Of_Time().
#             Note that met_base.py calls True_Solar_Noon().
#             New version of Earth_Perihelion that uses Python
#             dictionary and covers 1981 to 2060.
#             Small changes to Julian_Day() and Day_Angle().
#             Day_Angle() does not yet handle leap years.
#  May 2020.  Fixed small bug in Day_Angle().
#             Updated Current_Year() to use datetime.
#  Jul 2010.  Cleaned up, removed old GUI routines.
#             Replaced T_air & RH as args to radiation routines w/ W_p.
#  Jan 2009.  Converted from IDL.
#  Mar 2007.  Updates.
#  Feb 2007.  Added longwave radiation tools.
#  Aug 2005.  Updates.
#  Jul 2005.  Created; shortwave radiation only.
#
#------------------------------------------------------------------------
#
#  Note:  Earth_Perihelion() uses a table of values that only
#         spans the years 1981 to 2060.  For 2001 to 2100,
#         as well as other centuries, see:
#         http://astropixels.com/ephemeris/perap2001.html
#
#  Note:  NB!  Not yet ready to create grids of lats & 
#         lons for DEM with fixed-length pixels (e.g. UTM).
#
#  Note:  Functions that should be double-checked include:
#           Vernal_Equinox, Earth_Perihelion and
#           ET_Radiation_Flux_Slope
#------------------------------------------------------------------------
#
#----------------
#   Functions
#----------------
#   Current_Year()
#   Solar_Constant
#   Day_Angle
#   Eccentricity_Correction
#   Declination
#   Earth_Angular_Velocity    # (See Earth_Rotation_Rate)
#   Zenith_Angle
#   Solar_Elevation_Angle
#   Sunrise_Offset
#   Sunset_Offset
#   Day_Length
#   ET_Radiation_Flux
#------------------------------
#   Saturation_Vapor_Pressure
#   Vapor_Pressure
#   Dew_Point                  
#   Precipitable_Water_Content
#   Optical_Air_Mass
#   Dust_Attenuation
#   Atmospheric_Transmissivity
#   Direct_Radiation_Flux
#------------------------------
#   Scattering_Attenuation
#   Diffuse_Radiation_Flux
#------------------------------
#   Global_Radiation_Flux
#   BS_Radiation_Flux
#------------------------------
#   Longitude_Offset
#   Equivalent_Latitude
#   Noon_Offset_Slope
#   Sunrise_Offset_Slope
#   Sunset_Offset_Slope
#   Day_Length_Slope
#   ET_Radiation_Flux_Slope
#   Clear_Sky_Radiation
#------------------------------
#   Julian_Day
#   Days_Per_Year
#   Earth_Rotation_Rate
#   Earth_Tilt_Angle
#   Earth_Orbit_Eccentricity
#   Vernal_Equinox
#   Earth_Perihelion
#   Equation_of_Time
#   True_Solar_Noon
#------------------------------
#   Latitude_Grid
#   Longitude_Grid
#
#------------------------------------------------------------------------

import numpy as np
import os
from datetime import datetime

#------------------------------------------------------------------------
def Current_Year():

    # currentSecond = datetime.now().second
    # currentMinute = datetime.now().minute
    # currentHour   = datetime.now().hour
    # currentDay    = datetime.now().day
    # currentMonth  = datetime.now().month
    # currentYear   = datetime.now().year

    return np.int16( datetime.now().year )

#   Current_Year()    
#------------------------------------------------------------------------
def Solar_Constant():

    #---------------------------------------------------------
    # See: Wikipedia: Solar Constant (2023-08-29)
    # It varies from 1361 at solar min to 1362 at solar max
    # as measured by satellite.
    # The value 1367 is widely used, see:
    # https://www.sciencedirect.com/science/article/pii/S0360544210007565
    # Dingman may have used: 1367.
    #---------------------------------------------------------
    return np.float64(1361.5)   # [Watts / m^2]
    
#   Solar_Constant()
#------------------------------------------------------------------------
def Day_Angle( Julian_day, DEGREES=False ):

    #---------------------------------------------------------
    # Notes:  The Julian day does not need to be an integer;
    #         decimal values can be used for more precision.
    #---------------------------------------------------------

    #-------------------------------------    
    # Use this if Julian Day starts at 1
    #-------------------------------------
    ## angle = (2 * np.pi) * (Julian_day - np.float64(1)) / np.float64(365)

    #-------------------------------------    
    # Use this if Julian Day starts at 0
    #-------------------------------------
    # Don't use Days_Per_Year() here.
    #-----------------------------------------------------------
    # We should be using 366 vs. 365 for leap years, but would
    # then need to pass year to every Day_Angle() call.
    #-----------------------------------------------------------
    angle = (2 * np.pi) * Julian_day / np.float64(365)
            
    if (DEGREES):    
        angle = angle * (np.float64(180) / np.pi)
    
    return angle
    
#   Day_Angle()
#------------------------------------------------------------------------
def Eccentricity_Correction( day_angle ):

    #----------------------------------------
    # Note: Range is about (0.966, 1.035).
    # See Jupyter notebook for Meteorology.
    #----------------------------------------
    E0 = np.float64(1.000110) + \
         (np.float64(0.034221) * np.cos(day_angle)) + \
         (np.float64(0.001280) * np.sin(day_angle)) + \
         (np.float64(0.000719) * np.cos(np.float64(2) * day_angle)) + \
         (np.float64(0.000077) * np.sin(np.float64(2) * day_angle))
    
    return E0   # [unitless, ratio of two lengths]
    
#   Eccentricity_Correction()
#------------------------------------------------------------------------
def Declination( day_angle, DEGREES=False, DMS=False ):

    ########################################################
    # NB! Make sure that DEGREES and DMS default to False.
    ########################################################

    #-----------------------------------------------------------
    # Note:  The declination reaches its lowest value of -23.5
    #        degrees on the Winter Solstice (Dec. 21/22) and
    #        reaches its highest value of 23.5 degrees on the
    #        Summer Solstice (June 21/22).  It is zero for
    #        both the Vernal Equinox (Mar. 20/21) and the
    #        Autumnal Equinox (Sept. 22/23).  The value of
    #        23.4397 degrees is the fixed tilt angle of the
    #        Earth's axis from from the plane of the ecliptic.
    #-----------------------------------------------------------  
    delta = np.float64(0.006918) - \
            (np.float64(0.399912) * np.cos(day_angle)) + \
            (np.float64(0.070257) * np.sin(day_angle)) - \
            (np.float64(0.006758) * np.cos(np.float64(2) * day_angle)) + \
            (np.float64(0.000907) * np.sin(np.float64(2) * day_angle)) - \
            (np.float64(0.002697) * np.cos(np.float64(3) * day_angle)) + \
            (np.float64(0.001480) * np.sin(np.float64(3) * day_angle))
    
    #------------------------------------
    # Convert from radians to degrees ?
    #------------------------------------
    if (DEGREES):    
        delta = delta * (np.float64(180) / np.pi)
    
    #----------------------------------------
    # Convert from radians to "decimal DMS"
    #----------------------------------------
    if (DMS):    
        delta = delta * (np.float64(180) / np.pi)  # [decimal degrees]
        deg = np.int16(delta)
        min = np.int16((delta - deg) * np.float64(60))
        sec = np.int16(((delta - deg) * np.float64(60) - min) * np.float64(60))
        delta = deg + (min / np.float64(100)) + (sec / np.float64(10000))   # [decimal DMS, DD.MMSS]
    
    return delta
    
#   Declination()
#------------------------------------------------------------------------
def Earth_Angular_Velocity():

    #---------------------------------------------------
    # Notes:  Compare to Earth_Rotation_Rate function.
    #---------------------------------------------------
    deg_per_hour = np.float64(360) / np.float64(24)   #(equals 15)
    rad_per_hour = deg_per_hour * (np.pi / np.float64(180))
    
    return rad_per_hour
    
#   Earth_Angular_Velocity()
#------------------------------------------------------------------------
def Zenith_Angle( lat_deg, declination, th ):

    #----------------------------------------------------------
    # Notes: lat_deg has units of DEGREES and declination
    #        must have units of RADIANS.

    #        th is number of hours before (-) or or after (+)
    #        the true solar noon.

    #        Note that zenith_angle < 0 at night.
    #        Sunrise and sunset occur when zenith angle is
    #        equal to pi/2, so cos(Z)=0 and we can then
    #        solve for time offsets as th.
    #----------------------------------------------------------
    omega = Earth_Angular_Velocity()    # [radians / hour]
    lat_rad = lat_deg * (np.pi / np.float64(180))
    term1 = np.sin(lat_rad) * np.sin(declination)
    term2 = np.cos(lat_rad) * np.cos(declination) * np.cos(omega * th)
    
    return np.arccos(term1 + term2)   # [radians]
    
#   Zenith_Angle()
#------------------------------------------------------------------------
def Solar_Elevation_Angle( lat_deg, declination, th, DEGREES=False ):

    #---------------------------------------------------------- 
    # Note: This is the complement of the solar zenith angle.
    #       At sunrise and sunset, this angle is zero.
    #       See notes for that function for more info.
    #----------------------------------------------------------
    theta = Zenith_Angle( lat_deg, declination, th )
    angle = (np.pi/2 - theta)   # [radians]
    if (DEGREES):
        angle *= (180 / np.pi)  # [degrees]
    return angle
    
#   Solar_Elevation_Angle()
#------------------------------------------------------------------------
def Sunrise_Offset( lat_deg, declination ):

    #----------------------------------------------------------
    # Notes: lat and declination must have units of RADIANS.
    #        time has units of hours before true solar noon.

    #        If (abs(lat_deg) gt 66.5) we are above Arctic
    #        circle or below Antarctic circle.  This can also
    #        happen if lat_deg is an "equivalent latitude"
    #        for a slope.  In this case, the absolute value
    #        of the argument to ACOS can exceed one and there
    #        is either no sunrise or no sunset at the given
    #        location.
    #----------------------------------------------------------
    omega = Earth_Angular_Velocity()
    lat_rad = lat_deg * (np.pi / np.float64(180))
    
    #--------------------------------------------------
    # See Notes for TF_Tan function in utils_TF.py ??
    #--------------------------------------------------
    arg  = -np.float64(1) * np.tan(lat_rad) * np.tan(declination)
    arg  = np.minimum( np.maximum(-1, arg), 1 )
    time = (-np.float64(1) * np.arccos(arg) / omega)
    
    return time
    
#   Sunrise_Offset()
#------------------------------------------------------------------------
def Sunset_Offset( lat_deg, declination ):

    #----------------------------------------------------------
    # Notes: lat and declination must have units of RADIANS.
    #        time has units of hours before true solar noon.

    #        If (abs(lat_deg) gt 66.5) we are above Arctic
    #        circle or below Antarctic circle.  This can also
    #        happen if lat_deg is an "equivalent latitude"
    #        for a slope.  In this case, the absolute value
    #        of the argument to ACOS can exceed one and there
    #        is either no sunrise or no sunset at the given
    #        location.
    #----------------------------------------------------------
    omega   = Earth_Angular_Velocity()
    lat_rad = lat_deg * (np.pi / np.float64(180))
    
    #--------------------------------------------------
    # See Notes for TF_Tan function in utils_TF.py ??
    #--------------------------------------------------
    arg  = -np.float64(1) * np.tan(lat_rad) * np.tan(declination)
    arg  = np.minimum( np.maximum(-1, arg), 1 )
    time = (np.arccos(arg) / omega)
    
    return time
    
#   Sunset_Offset()
#------------------------------------------------------------------------
def Day_Length( lat_deg, Julian_day):

    day_angle   = Day_Angle( Julian_day )
    declination = Declination( day_angle, DEGREES=False )

    t_sr = Sunrise_Offset(lat_deg, declination)
    t_ss = Sunset_Offset(lat_deg,  declination)
    
    return (t_ss - t_sr)  # [hours]
    
#   Day_Length()
#------------------------------------------------------------------------
def ET_Radiation_Flux( lat_deg, Julian_day, th ):

    #------------------------------------------------------------
    # Notes:  This is the instantaneous extraterrestrial
    #         radiation flux on a horizontal plane at a time
    #         th hours before (-) or after (+) true solar noon.
    #------------------------------------------------------------
    I_sc  = Solar_Constant()            # [W / m^2]
    omega = Earth_Angular_Velocity()    # [radians / hour]
    #---------------------------------
    Gamma = Day_Angle(Julian_day)              # [radians]
    delta = Declination(Gamma)                 # [radians]
    E0    = Eccentricity_Correction(Gamma)     # [unitless]
    lat_rad = lat_deg * (np.pi / np.float64(180))
    #------------------------------------------------------------
    term1 = np.cos(delta) * np.cos(lat_rad) * np.cos(omega * th)
    term2 = np.sin(delta) * np.sin(lat_rad)
    K_ET  = I_sc * E0 * (term1 + term2)
 
    TEST = False
    if (TEST):
        print('Julian day =', Julian_day )
        print('min(lat_deg), max(lat_deg) =', lat_deg.min(), lat_deg.max() )
        print('I_sc =', I_sc)
        print('min(omega), max(omega) =', omega.min(), omega.max() )
        print('min(Gamma), max(Gamma) =', Gamma.min(), Gamma.max() )
        print('min(delta), max(delta) =', delta.min(), delta.max() )
        print('min(E0), max(E0)       =', E0.min(), E0.max() )
        print('min(K_ET), max(K_ET)   =', K_ET.min(), K_ET.max() )
        print()
           
    #-------------------------------------------------------------
    # NB! During local nightime hours, K_ET < 0.
    #     At equator (lat_deg=0), this occurs at th<-6 and th>6.
    #     When negative, return zero as shown.
    #-------------------------------------------------------------
    np.maximum( K_ET, 0.0, K_ET )    # in-place
    return K_ET    # [Watts / m^2]

        
#     K_ET_min = K_ET.min()
#     K_ET_max = K_ET.max()
#     if (K_ET_min < 0):
#         print('------------------------------------------')
#         print('ERROR in ET_Radiation_Flux():')
#         print('Incoming radiation flux should be > 0.')
#         print('min, max =', K_ET_min, K_ET_max)
#         print('------------------------------------------')
#         print()
#     return K_ET    # [Watts / m^2]
    
#   ET_Radiation_Flux()
#------------------------------------------------------------------------
def Saturation_Vapor_Pressure(T, method='BRUTSAERT', MBAR=False):

        if (method == 'BRUTSAERT'):  
            #------------------------------
            # Use Brutsaert (1975) method
            #------------------------------
            term1 = (np.float64(17.3) * T) / (T + np.float64(237.3))
            e_sat = np.float64(0.611) * np.exp(term1)  # [kPa]
        else:    
            #-------------------------------
            # Use Satterlund (1979) method     #### DOUBLE CHECK THIS (7/26/13)
            #-------------------------------
            term1 = np.float64(2353) / (T + np.float64(273.15))
            e_sat = np.float64(10) ** (np.float64(11.4) - term1)   # [Pa]
            e_sat = (e_sat / np.float64(1000))  # [kPa]

        #-----------------------------------
        # Convert units from kPa to mbars?
        #-----------------------------------
        if (MBAR):    
            e_sat = (e_sat * np.float64(10))   # [mbar]
        return e_sat
        
#   Saturation_Vapor_Pressure()
#------------------------------------------------------------------------
def Vapor_Pressure(T, rel_humidity, MBAR=False):

        #-------------------------------------------------         
        # RH is in [0,1], so e gets same units as e_sat.
        # So we never need to convert units of e.
        #-------------------------------------------------
        e_sat = Saturation_Vapor_Pressure(T, MBAR=MBAR)   
        e = (rel_humidity * e_sat)
        return e

#   Vapor_Pressure()
#------------------------------------------------------------------------
def Dew_Point( T, rel_humidity ):

    #---------------------------------------------------------
    # Notes:  Temps are in degrees C, and vapor pressure
    #         units are kPa.  Relative humidity is unitless.
    #---------------------------------------------------------    
    vp  = Vapor_Pressure(T, rel_humidity, MBAR=False)  # [kPa]
    top = np.log(vp) + np.float64(0.4926)
    bot = np.float64(0.0708) - np.float64(0.00421) * np.log(vp)
    Td  = (top / bot)
    
    return Td
    
#   Dew_Point()
#------------------------------------------------------------------------
def Precipitable_Water_Content(T, rel_humidity):

    #---------------------------------------
    # Note: Wp > 0 always, even if Td < 0.
    #---------------------------------------
    Td = Dew_Point(T, rel_humidity)   # [degrees C]
    Wp = np.float64(1.12) * np.exp(np.float64(0.0614) * Td)    # [centimeters]
    
    return Wp  # [centimeters]
    
#   Precipitable_Water_Content()
#------------------------------------------------------------------------
def Optical_Air_Mass( lat_deg, declination, th ):

    #------------------------------------------------------------
    # Notes: This is a dimensionless number that gives the
    #        relative path length (greater than 1) that
    #        radiation must travel through the atmosphere as
    #        the result of not entering at a right angle.

    #        th is number of hours before (-) or or after (+)
    #        the true solar noon.

    #        Dingman gives only a table (Figure E-4, p. 605)
    #        of daily average values as a function of lat
    #        and declination.

    #        The approximation formula used here is widely
    #        used and is from Kasten and Young (1989).
    #------------------------------------------------------------
    # Note: This version was used up to 2023-08-29, but is
    #       wrong because Z was not converted to degrees in
    #       term1.  In term1, note that 96.07995 = (90 + c),
    #       in degrees, but Z is in radians.  Note also that
    #       sin(gamma_rad) = cos(Z_rad).  
    #------------------------------------------------------------    
#     Z = Zenith_Angle(lat_deg, declination, th)  # [radians]
#     term1 = (np.float64(96.07995) - Z) ** (-np.float64(1.6364))
#     denom = np.cos(Z) + (np.float64(0.50572) * term1)
#     M_opt = (np.float64(1) / denom)
    
    #----------------------------------------------------------
    # Kasten (1965) gave an approximation formula for the
    # "relative optical air mass" that was widely used.
    # (The relative optical air mass is dimensionless.)
    # This formula was a function of the "solar elevation
    # angle", denoted by gamma, which is the complement of
    # the solar zenith angle.  The formula had 3 fitting
    # parameters: a, b and c.  Kasten and Young (1989) gave
    # an improved approximation formula (equation 3) in
    # which only the 3 fitting parameters were different.
    # This monotonic, rapidly decreasing function has:
    #    f[0] = 37.9196, and f[90] = 0.999712.   
    #----------------------------------------------------------
    a = 0.50572
    b = 6.07995  # [degrees]
    c = 1.6364
#     KASTEN_1965 = False
#     if (KASTEN_1965):
#         a = 0.1500
#         b = 3.885  # [degrees]
#         c = 1.253

    Z = Zenith_Angle(lat_deg, declination, th)  # [radians]
    Z_deg = Z * (180 / np.pi)  # [degrees]
    gamma = (90.0 - Z_deg)     # [degrees]
    #------------------------------------------------------
    # NB!  Zenith angle = pi/2 (or 90 degrees) at sunrise
    #      and sunset and is > 90 degrees at night.
    #      Therefore, gamma < 0 at night, but formula
    #      doesn't allow gamma < 0.
    #      To avoid error, set gamma = 0 at night.
    #      Can use gamma.size if 0D or 2D ndarray.
    #      Check if gamma is an ndarray. (2023-08-29)
    #------------------------------------------------------
    if (isinstance(gamma, np.ndarray)):
        np.maximum( gamma, 0.0, gamma)  # in-place
    else:
        gamma = max(gamma, 0)
 
    term1 = np.sin( gamma * (np.pi / 180) )
    term2 = a / (gamma + b)**c
    M_opt = (np.float64(1) / (term1 + term2))

    return M_opt   # [unitless]
    
#   Optical_Air_Mass()
#------------------------------------------------------------------------
def Dust_Attenuation():

    #----------------------------------------------------------
    # Notes:  Typical clear-sky values are between 0 and 0.2.
    #         Bolsenga (1964) cites values of:
    #             0.00 to 0.05,   remote stations
    #             0.03 to 0.10,   moderate-sized cities
    #             0.10 to 0.13,   larger metro. areas
    #         See Dingman, p. 604-605.
    #----------------------------------------------------------
    return np.float64(0.08)
    
#   Dust_Attenuation()
#------------------------------------------------------------------------
def Atmospheric_Transmissivity( lat_deg, Julian_day, W_p, 
                                th, gamma_dust=None ):

    #------------------------------------------------------------
    # Notes:  W_p is precipitable water content in centimeters,
    #         which depends on air temp and relative humidity.
    #         W_p = Precipitable_Water_Content(T, rel_humidity).
    #
    #         lat_deg is latitude in decimal degrees.
    #         Atmospheric Trans. = tau = unitless and in [0,1].
    #------------------------------------------------------------
    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()
        
    Gamma  = Day_Angle(Julian_day)   # [radians]
    delta  = Declination(Gamma)      # [radians]
    #----------------------------------------------------
    # Note: 0 <= M_opt <= 37.9196, and a_sa & b_sa < 0.
    #----------------------------------------------------    
    a_sa   = -np.float64(0.1240) - (np.float64(0.0207) * W_p)
    b_sa   = -np.float64(0.0682) - (np.float64(0.0248) * W_p)
    M_opt  = Optical_Air_Mass(lat_deg, delta, th) # [unitless]
    tau_sa = np.exp(a_sa + (b_sa * M_opt))
    tau    = (tau_sa - gamma_dust)          # [unitless]
    
    return np.minimum( np.maximum(tau, 0), 1)
    
#   Atmospheric_Transmissivity()
#------------------------------------------------------------------------
def Direct_Radiation_Flux( lat_deg, Julian_day, W_p,
                           th, gamma_dust=None ):

    #------------------------------------------------------------
    # Notes:  W_p is precipitable water content in centimeters,
    #         which depends on air temp and relative humidity.
    #         W_p = Precipitable_Water_Content(T, rel_humidity).
    #
    #         lat_deg is latitude in decimal degrees.
    #         Atmospheric Trans. = tau = unitless and in [0,1].
    #------------------------------------------------------------
    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()
        
    tau   = Atmospheric_Transmissivity(lat_deg, Julian_day,
                                       W_p, th, gamma_dust)
    K_ET  = ET_Radiation_Flux(lat_deg, Julian_day, th)
    K_dir = (tau * K_ET)
    
    return K_dir
    
#   Direct_Radiation_Flux()
#------------------------------------------------------------------------
def Scattering_Attenuation( lat_deg, Julian_day, W_p,
                            th, gamma_dust=None ):

    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()

    Gamma = Day_Angle(Julian_day)   # [radians]
    delta = Declination(Gamma)      # [radians]
    #----------------------------------------------------
    a_s   = -np.float64(0.0363) - (np.float64(0.0084) * W_p)
    b_s   = -np.float64(0.0572) - (np.float64(0.0173) * W_p)
    M_opt = Optical_Air_Mass(lat_deg, delta, th)  # [unitless]
    tau_s = np.exp(a_s + (b_s * M_opt))
    gam_s = (1 - tau_s) + gamma_dust
    
    return gam_s
    
#   Scattering_Attenuation()
#------------------------------------------------------------------------
def Diffuse_Radiation_Flux( lat_deg, Julian_day, W_p,
                            th, gamma_dust=None ):

    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()

    gam_s = Scattering_Attenuation(lat_deg, Julian_day, W_p,
                                   th, gamma_dust)  # [unitless]
    K_ET  = ET_Radiation_Flux(lat_deg, Julian_day, th)
    K_dif = np.float64(0.5) * gam_s * K_ET

    return K_dif   # [Watts / meter^2]
    
#   Diffuse_Radiation_Flux()
#------------------------------------------------------------------------
def Global_Radiation_Flux( lat_deg, Julian_day, W_p,
                           th, gamma_dust=None ):

    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()
        
    K_dir = Direct_Radiation_Flux(lat_deg, Julian_day, W_p,
                                  th, gamma_dust)
    
    K_dif = Diffuse_Radiation_Flux(lat_deg, Julian_day, W_p,
                                   th, gamma_dust)
    
    K_global = (K_dir + K_dif)
    
    return K_global   # [Watts / m^2]
    
#   Global_Radiation_Flux()
#------------------------------------------------------------------------
def BS_Radiation_Flux( lat_deg, Julian_day, W_p,
                       albedo, th, gamma_dust=None):

    #----------------------------------------------------------
    # Notes:  Compute the backscattered radiation flux.
    #         A table of typical albedos is given by Dingman,
    #         Table D-2 on page 584.
    #----------------------------------------------------------
    if (gamma_dust is None):    
        gamma_dust = Dust_Attenuation()
        
    gam_s = Scattering_Attenuation(lat_deg, Julian_day, W_p,
                                   th, gamma_dust)
    
    K_global = Global_Radiation_Flux(lat_deg, Julian_day, W_p,
                                     th, gamma_dust)
    
    #----------------
    # For debugging
    #----------------
    #print 'min(gam_s),    max(gam_s)    =', gam_s.min(),    gam_s.max()
    #print 'min(albedo),   max(albedo)   =', albedo.min(),   albedo.max()
    #print 'min(K_global), max(K_global) =', K_global.min(), K_global.max()
    
    K_bs = np.float64(0.5) * gam_s * albedo * K_global
    
    return K_bs
    
#   BS_Radiation_Flux()
#------------------------------------------------------------------------
def Longitude_Offset( lat_deg, alpha, beta ):

    #-------------------------------------------------------------
    # Notes:  beta  = "slope angle" satisfies slope = tan(beta).
    #         alpha = "aspect_angle" or azimuth is measured
    #                 clockwise from north.
    #         Both angles have units of radians.
    #         Returned value, dlon, has units of radians.

    #         If (alpha eq 0) or (beta eq 0) then (term1 eq 0)
    #         and this offset will be 0.
    #-------------------------------------------------------------
    lat_rad = lat_deg * (np.pi / np.float64(180))
    term1   = np.sin(beta) * np.sin(alpha)
    term2   = np.cos(beta) * np.cos(lat_rad)
    term3   = np.sin(beta) * np.sin(lat_rad) * np.cos(alpha)
    dlon    = np.arctan(term1 / (term2 - term3))
    
    return dlon  # [radians]
    
#   Longitude_Offset()
#------------------------------------------------------------------------
def Equivalent_Latitude( lat_deg, alpha, beta, DEGREES=False ):

    #-------------------------------------------------------------
    # Notes:  beta  = "slope angle" satisfies slope = tan(beta).
    #         alpha = "aspect_angle" or azimuth is measured
    #                 clockwise from north.
    #         Both angles have units of radians.

    #         Note that if (beta eq 0), then lat_eq (in deg) is
    #         always equal to lat_deg.  Also, beta will always
    #         be in the range [0, pi/2].
    #-------------------------------------------------------------
    lat_rad = lat_deg * (np.pi / np.float64(180))
    term1   = np.sin(beta) * np.cos(alpha) * np.cos(lat_rad)
    term2   = np.cos(beta) * np.sin(lat_rad)
    
    eq_lat   = np.arcsin(term1 + term2)
    
    #--------------------------------------
    # Convert to degrees?  Sunrise_Offset
    # function requires lat in degrees.
    #--------------------------------------
    if (DEGREES):    
        eq_lat = eq_lat * (np.float64(180) / np.pi)
        return eq_lat    # [degrees]
    else:    
        return eq_lat    # [radians]

#   Equivalent_Latitude()
#------------------------------------------------------------------------
def Noon_Offset_Slope( lat_deg, alpha, beta ):

    dlon  = Longitude_Offset(lat_deg, alpha, beta)  # [radians]
    omega  = Earth_Angular_Velocity()
    t_noon = -np.float64(1) * dlon / omega
    
    return t_noon
    
#   Noon_Offset_Slope()
#------------------------------------------------------------------------
def Sunrise_Offset_Slope( lat_deg, Julian_day, alpha, beta ):

    #-----------------------------------------------------------
    # Notes:  beta  = "slope angle" satisfies slope = tan(beta).
    #         alpha = "aspect_angle" or azimuth is measured
    #                 clockwise from north.
    #         Both angles have units of radians.
    #-----------------------------------------------------------
    Gamma = Day_Angle(Julian_day)         # [radians]
    delta = Declination(Gamma)            # [radians]
    eq_lat_deg = Equivalent_Latitude(lat_deg, alpha, beta, DEGREES=True)
    #----------------------------------------------------
    t_noon = Noon_Offset_Slope(lat_deg, alpha, beta)
    t_sr   = Sunrise_Offset(eq_lat_deg, delta)
    t_sr   = t_sr + t_noon   # (plus sign is correct)
    
    #-----------------------------------------------
    # This is what Dingman does in his spreadsheet
    #-----------------------------------------------
    t_sr = np.maximum( t_sr, Sunrise_Offset(lat_deg, delta) )
    
    return t_sr   # [hours]
    
#   Sunrise_Offset_Slope()
#------------------------------------------------------------------------
def Sunset_Offset_Slope( lat_deg, Julian_day, alpha, beta ):

    #-----------------------------------------------------------
    # Notes:  beta  = "slope angle" satisfies slope = tan(beta).
    #         alpha = "aspect_angle" or azimuth is measured
    #                 clockwise from north.
    #         Both angles have units of radians.
    #-----------------------------------------------------------
    Gamma = Day_Angle(Julian_day)         # [radians]
    delta = Declination(Gamma)            # [radians]
    eq_lat_deg = Equivalent_Latitude(lat_deg, alpha, beta, DEGREES=True)
    #----------------------------------------------------
    t_noon = Noon_Offset_Slope(lat_deg, alpha, beta)
    t_ss   = Sunset_Offset(eq_lat_deg, delta)
    t_ss   = t_ss + t_noon
    
    #-----------------------------------------------
    # This is what Dingman does in his spreadsheet
    #-----------------------------------------------
    t_ss = np.minimum( t_ss, Sunset_Offset(lat_deg, delta) )
    
    return t_ss   # [hours]
    
#  Sunset_Offset_Slope()
#------------------------------------------------------------------------
def Day_Length_Slope( lat_deg, Julian_day, alpha, beta ):

    t_sr = Sunrise_Offset_Slope(lat_deg, Julian_day, alpha, beta)
    t_ss = Sunset_Offset_Slope(lat_deg, Julian_day, alpha, beta)
    
    return (t_ss - t_sr)  # [hours]
    
#   Day_Length_Slope()
#------------------------------------------------------------------------
def ET_Radiation_Flux_Slope( lat_deg, Julian_day, th, alpha, beta ):

    #-------------------------------------------------------------
    # Notes:  This is the instantaneous extraterrestrial
    #         radiation flux on a sloping plane.
    #-------------------------------------------------------------
    # Notes:  beta  = "slope angle" satisfies slope = tan(beta).
    #         alpha = "aspect_angle" or azimuth is measured
    #                 clockwise from north.
    #         Both angles have units of radians.
    #-------------------------------------------------------------
    I_sc = Solar_Constant()             # [W / m^2]
    omega = Earth_Angular_Velocity()    # [radians / hour]
    #---------------------------------
    Gamma = Day_Angle(Julian_day)                # [radians]
    delta = Declination(Gamma)                   # [radians]
    E0 = Eccentricity_Correction(Gamma)          # [unitless]
    #--------------------------------------------------------------
    lat_eq = Equivalent_Latitude(lat_deg, alpha, beta) #[radians]
    dlon = Longitude_Offset(lat_deg, alpha, beta)    #[radians]
    #--------------------------------------------------------------
    term1 = np.cos(delta) * np.cos(lat_eq)
    term2 = np.cos((omega * th) + dlon)
    term3 = np.sin(lat_eq) * np.sin(delta)
    K_ET = I_sc * E0 * ((term1 * term2) + term3)
    
    #----------------
    # For debugging
    #-------------------------------------------
    # NaNs should only occur because alpha and
    # beta grids have them along the edges.
    #-------------------------------------------
    #print 'min(I_sc),  max(I_sc)  =', I_sc.min(),  I_sc.max()
    #print 'min(E0),    max(E0)    =', E0.min(),    E0.max()
    #print 'min(term1), max(term1) =', term1.min(), term1.max()
    #print 'min(term2), max(term2) =', term2.min(), term2.max()
    #print 'min(term3), max(term3) =', term3.min(), term3.max()
    #print ' '
    
    #----------------------------
    # This shouldn't be needed.
    #----------------------------
    K_ET = np.maximum( K_ET, 0 )
    
    return K_ET    # [Watts / m^2]
    
#   ET_Radiation_Flux_Slope()
#------------------------------------------------------------------------
def Clear_Sky_Radiation( lat_deg, Julian_day, W_p,
                         TSN_offset, alpha, beta,
                         albedo, gamma_dust ):

    #--------------------------------------------------------
    # Notes:  I think K_cs is the same as the Qnet required
    #         for the energy-balance routines in TopoFlow.
    #         Both have units of Watts/m^2.
    #--------------------------------------------------------
    #** if (n_elements(gamma_dust) eq 0) then $
    #**     gamma_dust = Dust_Attenuation()
    #--------------------------------------------------------
    tau   = Atmospheric_Transmissivity(lat_deg, Julian_day,
                                       W_p, TSN_offset,
                                       gamma_dust)
    K_ET  = ET_Radiation_Flux_Slope(lat_deg, Julian_day, TSN_offset,
                                    alpha, beta)
    K_dif = Diffuse_Radiation_Flux(lat_deg, Julian_day, W_p,
                                   TSN_offset, gamma_dust)
    K_bs  = BS_Radiation_Flux(lat_deg, Julian_day, W_p,
                              albedo, TSN_offset, gamma_dust)
    
    K_cs  = (tau * K_ET) + K_dif + K_bs
    
    #----------------
    # For debugging
    #----------------
    TEST = False
    if (TEST):
        print( 'min(alpha), max(alpha) =', alpha.min(), alpha.max() )
        print( 'min(beta),  max(beta)  =', beta.min(),  beta.max() )
        print( 'min(tau),   max(tau)   =', tau.min(),   tau.max() )
        print( 'min(K_ET),  max(K_ET)  =', K_ET.min(),  K_ET.max() )
        print( 'min(K_dif), max(K_dif) =', K_dif.min(), K_dif.max() )
        print( 'min(K_bs),  max(K_bs)  =', K_bs.min(),  K_bs.max() )
        print( 'min(K_cs),  max(K_cs)  =', K_cs.min(),  K_cs.max() )
        print()
    
    #-----------------------------------------------
    # Set K_cs to zero between (local) dusk & dawn
    # NB!  These next two variables are GRIDS.
    #-----------------------------------------------
    T_sr = Sunrise_Offset_Slope(lat_deg, Julian_day, alpha, beta)
    T_ss = Sunset_Offset_Slope(lat_deg,  Julian_day, alpha, beta)

    #------------------------------------------------
    # Use of LE & GE also takes care of case where
    # Tsr = T_ss = 0, when abs(eq_lat_deg) gt 66.5.
    #------------------------------------------------
    # Without WHERE call, cProfile time was reduced
    # from 0.604 to
    #------------------------------------------------    
    dark = np.logical_or( (TSN_offset <= T_sr), (TSN_offset >= T_ss) )  
    K_cs[ dark ] = np.float64(0)
            
    #------------------------------------------------
    # Use of LE & GE also takes care of case where
    # Tsr = T_ss = 0, when abs(eq_lat_deg) gt 66.5.
    #------------------------------------------------
#     dark   = np.where( np.logical_or((TSN_offset <= T_sr),
#                                      (TSN_offset >= T_ss)) )
#     n_dark = np.size( dark[0] )
#     if (n_dark != 0):    
#         K_cs[ dark ] = np.float64(0)
    
    return K_cs   # [Watts / m^2]
    
#   Clear_Sky_Radiation()
#------------------------------------------------------------------------
def Julian_Day( month_num, day_num, hour_num=None, year=None ):

    #-----------------------------------------------------------
    # NB!  month_num is an integer between 1 and 12 inclusive.
    #      day_num is day of the month.
    #      hour_num is in [0,24], in UTC time zone.
    #      This function can handle leap years, given year.
    #         Julian_Day(1,1,0)    = 0.0  # Jan 1 starts at midnight.
    #         Julian_Day(1,1,24)   = 1.0
    #         Julian_Day(2,1,0)    = 31.0
    #         Julian_Day(12,31,0)  = 364.0
    #         Julian_Day(12,31,24) = 365.0
    #         Julian_Day(12,31,24, year=2024) = 366.0
    #      We are using Julian_Day to compute Day_Angle.   
    #-----------------------------------------------------------
    # NB!  For our purposes we need the Julian day of the year,
    #      not the Julian day from the start of Julian period.
    #      See:  https://en.wikipedia.org/wiki/Julian_day
    #      https://calendars.fandom.com/wiki/Julian_day_number
    #------------------------------------------------------------
    if (year is None) or ((year % 4) != 0):
        month_days = np.array([0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])        
    else:
        # Leap year     
        month_days = np.array([0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) 

    #-------------------------------------------------------   
    # This method gives JD in [0,365]
    # Note that JD < 1 between midnight of 12-31 (start of
    # the day 01-01) to midnight of 01-01 (start of the
    # day 01-02).  That is, JD=1 when one, full 24-hour
    # day has elapsed in the given year, and JD=365 when
    # a full 24-day has elapsed on the last day of the
    # year (or 366 for leap years).
    #-------------------------------------------------------
    JD = np.sum(month_days[:month_num]) + np.maximum(day_num - 1, 0)
    
    #---------------------------------------------------
    # (2020-05-10) Online source, Parkin (2017), says:
    #        JD(Jan. 1) = 1, JD(Feb. 1) = 32.
    #---------------------------------------------------
    # JD = np.sum(month_days[:month_num]) + np.maximum(day_num,1)
 
    if (hour_num is not None):
        JD = JD + (hour_num / np.float64(24))

        #----------------------------------------------------            
        # Wikipedia says to subtract 12; day starts at noon
        #----------------------------------------------------
        # JD = JD + (hour_num - 12) / np.float64(24)   # (2021-07-29)  

    return JD
    
#   Julian_Day()
#------------------------------------------------------------------------
def Days_Per_Year( SIDEREAL=False ):

    #----------------------------------------------------------
    # Notes: A day is typically defined as a mean solar day
    #        and contains exactly 24 hours.  A year is then
    #        the number of these days required for the Earth
    #        to trace out one complete orbit, and is given by
    #        365.24219 days (365.2422).  A year is also known
    #        as a tropical year.

    #        A sidereal day is the length of time that it
    #        takes for the Earth to spin through 360 degrees
    #        on its axis.

    #        1 mean solar day = 86400 seconds = 24 hours

    #        1 sidereal day   = 86164.09 secs = 23.934470 hrs.

    #        This function returns the number of solar days
    #        in one year (complete orbit) by default, but can
    #        also return the number of sidereal days in one
    #        year by setting the SIDEREAL keyword.
    #
    #        An approximate value of 365.2425 for the number
    #        solar days per year is used in accordance with
    #        the Gregorian calendar.
    #----------------------------------------------------------
    if (SIDEREAL):    
        #-------------------------------------
        # Return the number of sidereal days
        #-------------------------------------
        n_days = np.float64(366.2425)
    else:    
        #----------------------------------
        # Return the number of solar days
        # (Gregorian calendar convention)
        #----------------------------------
        n_days = np.float64(365.2425)
    
    return n_days   # [days]
    
#   Days_Per_Year()
#------------------------------------------------------------------------
def Earth_Rotation_Rate( PER_HOUR=False ):

    #------------------------------------------------------------
    # Notes:  Compare to Earth_Angular_Velocity function.
    #         The default is to return the rotation rate
    #         as the total number of radians rotated in one
    #         Earth orbit (tropical year).  About 2294.8863.
    #------------------------------------------------------------
    #         1 mean solar day = 86400 seconds = 24 hours
    #         1 sidereal day   = 86164.09 secs = 23.934470 hrs.
    #------------------------------------------------------------
    if (PER_HOUR):    
        Omega = np.float64(2) * np.pi / np.float64(24)   # [radians / hour]
    else:    
        DPY = Days_Per_Year(SIDEREAL=True)
        Omega = np.float64(2) * np.pi * DPY   # [radians / year]
    
    return Omega
    
#   Earth_Rotation_Rate()
#------------------------------------------------------------------------
def Earth_Tilt_Angle( DEGREES=False ):

    #------------------------------------------------------
    # Note:  The Earth's tilt angle is slowly decreasing.
    #        It is also known as the "obliquity".
    #------------------------------------------------------
    angle = np.float64(23.4397)  # [degrees]
    
    if not(DEGREES):    
        angle = angle * (np.pi / np.float64(180))  #[radians]
    
    return angle  # [degrees or radians]
    
#   Earth_Tilt_Angle()
#------------------------------------------------------------------------
def Earth_Orbit_Eccentricity():

    #----------------------------------------------------
    # Notes:  Return the eccentricity of Earth's orbit,
    #         which measures the difference between its
    #         elliptical orbit and a circular orbit.
    #         It is computed as:  e = (b-a)/a, where
    #         a and b are the semi-major and semi-minor
    #         axes of the elliptical orbit.
    #----------------------------------------------------
    return np.float64(0.016713)
    
#   Earth_Orbit_Eccentricity()
#------------------------------------------------------------------------
def Vernal_Equinox( year ):

    #---------------------------------------------------------------
    # Notes: This function assumes that vernal equinoxes from one
    #        year to the next are separated by exactly 365.2425
    #        days, or one tropical year, in accordance with the
    #        Gregorian calendar.  The difference between using this
    #        value and a "true" tropical year of 365.2422 days will
    #        be about 2.88 hours in 400 years.

    #        The time of vernal equinox for the year 2000 A.D. is
    #        March 20, 7:36 GMT [NASA Reference Publication 1349,
    #        Oct. 1994].

    #        We assume here that the vernal equinox for year 2000
    #        will be on March 20, 7:30, or 79.3125 days from 2000
    #        January 1, hour 0.  Vernal equinoxes for other years
    #        are returned as the total number of days since
    #        2000 January 1, hour 0.

    #        Note that:   79.3125 = 31 + 29 + 19 + 7.5/24.
    #---------------------------------------------------------------
    DPY = Days_Per_Year()
    VE_2000 = np.float64(79.3125)
    VE = VE_2000 + DPY * (year - np.float64(2000))
    
    return VE
    
#   Vernal_Equinox()
#------------------------------------------------------------------------
def Earth_Perihelion( year=None ):

    #-------------------------------------------------------------
    # NOTES:  Perihelion refers to the point along the orbit of
    #         a planet around the sun where it is closest to the
    #         sun.  For Earth, this typically occurs between the
    #         dates of January 2 and January 5.  This function
    #         returns the time when this event occurs as a
    #         Julian date in the given year.
    #-------------------------------------------------------------
    # Notes:  New version that goes through 2050 (2021-07-28).
    # See: http://www.astropixels.com/ephemeris/perap2001.html
    #-------------------------------------------------------------
    if (year is None):    
        year = Current_Year()
        
    if (year < 1981) or (year > 2060):
        print('WARNING: Earth_Perihelion is not available')
        print('         for the year:', year)
        print('         Will use current year instead.') 
        year = Current_Year()
                   
    #--------------------------------------------    
    # Store data in a dictionary (1981 to 2060)
    #--------------------------------------------
    Tp_dict = {
    1981:(2,2),  1982:(4,11), 1983:(2,15), 1984:(3,22), 1985:(3,20),
    1986:(2,5),  1987:(4,23), 1988:(3,0),  1989:(1,22), 1990:(4,17),
    1991:(3,3),  1992:(3,15), 1993:(4,3),  1994:(2,6),  1995:(4,11),
    1996:(4,7),  1997:(2,0),  1998:(4,21), 1999:(3,13), 2000:(3,5),
    2001:(4,9),  2002:(2,14), 2003:(4,5),  2004:(4,18), 2005:(2,1),
    2006:(4,15), 2007:(3,20), 2008:(3,0),  2009:(4,15), 2010:(3,0),
    2011:(3,19), 2012:(5,0),  2013:(2,5),  2014:(4,12), 2015:(4,7), 
    2016:(2,23), 2017:(4,14), 2018:(3,6),  2019:(3,5),  2020:(5,8),
    2021:(2,14), 2022:(4,7),  2023:(4,16), 2024:(3,1),  2025:(4,13),
    2026:(3,17), 2027:(3,3),  2028:(5,12), 2029:(2,18), 2030:(3,10),
    2031:(4,21), 2032:(3,5),  2033:(4,12), 2034:(4,5),  2035:(3,1),
    2036:(5,14), 2037:(3,4),  2038:(3,5),  2039:(5,7),  2040:(3,12),
    2041:(3,22), 2042:(4,9),  2043:(2,22), 2044:(5,13), 2045:(3,15),
    2046:(3,1),  2047:(5,12), 2048:(3,18), 2049:(3,10), 2050:(4,20),
    2051:(3,6),  2052:(5,9),  2053:(3,22), 2054:(2,18), 2055:(5,12),
    2056:(4,4),  2057:(3,3),  2058:(5,4),  2059:(3,11), 2060:(4,23) }

    #----------------------------------------
    # Get day and hour from table for given
    # year and convert to a Julian day.
    #----------------------------------------
    (Tp_day, Tp_hour) = Tp_dict[ year ]
    Tp_Julian_Day     = Julian_Day(1, Tp_day, Tp_hour)
    return Tp_Julian_Day
    
#   Earth_Perihelion()
#------------------------------------------------------------------------
# def Earth_Perihelion_old( year=None ):

    #-------------------------------------------------------------
    # NOTES:  Perihelion refers to the point along the orbit of
    #         a planet around the sun where it is closest to the
    #         sun.  For Earth, this typically occurs between the
    #         dates of January 2 and January 5.  This function
    #         returns the time when this event occurs as a
    #         Julian date in the given year.
    #-------------------------------------------------------------
#     if (year is None):    
#         year = Current_Year()
#         
# #     if (year < 1992) or (year > 2020):    
# #         year = Current_Year()
#     
#     #--------------------------------------------------------------
#     # Use published values from a table for the years 1992-2020.
#     #-------------------------------------------------------------
#     Tp_years = np.arange(29, dtype='int16') + 1992
#     
#     Tp_days  = np.array([3, 4, 2, 4, 4, 2, 4, 3, 3, 4, 2, 4, 4, 2, 4, 3, \
#                       3, 4, 3, 3, 5, 2, 4, 4, 2, 4, 3, 3, 5])
#     
#     Tp_hours = np.array([15, 3, 6, 11, 7, 0, 21, 13, 5, 9, 14, 5, 18, 1, 15, \
#                       20, 0, 15, 0, 19, 0, 5, 12, 7, 23, 14, 6, 5, 8])
#     
#     #----------------------------------------
#     # Get day and hour from table for given
#     # year and convert to a Julian day.
#     #----------------------------------------
#     w       = np.where( Tp_years == year )
#     Tp_vals = Julian_Day(1, Tp_days, Tp_hours)
#     Tp_JD   = Tp_vals[w]
#     
#     return Tp_JD
#     
# #   Earth_Perihelion_old()
#------------------------------------------------------------------------
def Equation_Of_Time( Julian_day, year=None,
                      DEGREES=False, DMS=False ):

    ############################################################
    # NB!  Should DEGREES and DMS both be False by default ??
    ############################################################

    #-----------------------------------------------------
    # Notes: The so-called "equation of time" gives the
    #        time difference between true solar noon and
    #        local clock noon, without accounting for any
    #        arbitrary time zone adjustments.  The latter
    #        are determined by humans and would introduce
    #        a whole-number offset from local clock noon
    #        in hours.

    #        The equation of time is closely related to
    #        the figure-8-shaped "analemma".

    #        Note that TE equals zero at 4 different
    #        times during the year.

    #        To test this against tables, try this:
    #            IDL>  JD = dindgen(365) + 1d
    #            IDL>  TE = Equation_of_Time(JD)
    #            IDL>  plot, JD, TE
    #            IDL>  Gamma = Day_Angle(JD)
    #            IDL>  delta = Declination(Gamma)
    #----------------------------------------------------- 
    if (year is None):   
        year = Current_Year()
        
#     if (year < 1981) or (year > 2060):    
#         year = Current_Year()
    
    #--------------------------------
    # Eccentricity of Earth's orbit
    # Computed as: e = (b-a)/a
    #--------------------------------
    ## print '### Computing e...'
    e = Earth_Orbit_Eccentricity()
    
    #----------------------------------
    # Earth's tilt angle or obliquity
    # (which is slowly decreasing)
    #----------------------------------
    ## print '### Computing eps...'
    eps = Earth_Tilt_Angle()
    
    #---------------------------------------
    # Number of mean solar days that occur
    # in one complete Earth orbit
    #---------------------------------------
    ## print '### Computing days_per_year...'
    days_per_year = Days_Per_Year()
    
    #------------------------------------
    # Get Julian date of the perihelion
    #------------------------------------
    ## print '### Computing Tp_JD...'
    Tp_JD = Earth_Perihelion(year)
    
    #----------------------------------------------------
    # Compute the mean anomaly, or the angular distance
    # from perihelion that would be travelled by a
    # uniformly moving (mean) sun.  It is zero at the
    # perihelion.
    #----------------------------------------------------
    twopi = np.float64(2) * np.pi
    M = (twopi / days_per_year) * (Julian_day - Tp_JD)    # [radians]
    M = (M + twopi) % twopi
    
    #------------------------------------
    # Get "longitude of the perihelion"
    #----------------------------------------------------
    # This is the angle between the semi-major axis
    # (the line of apsides) and a line between the Sun
    # and the Earth at the time of the Vernal Equinox.
    # Note that celestial longitudes are measured from
    # the Vernal Equinox (analogous to prime meridian).
    #----------------------------------------------------
    # omega is roughly equal to 4.9358 radians,
    # or -77.20 degrees or 282.8 degrees.
    #--------------------------------------------
    ## print '### Computing VE_JD...'
    ## year0 = np.int16(2000)
    ## VE_JD = Vernal_Equinox(year0)
    VE_JD = Vernal_Equinox( year )
    PT = (np.float64(365) + Tp_JD) - VE_JD    # [days, about 287]
    omega = twopi * (PT / days_per_year)   # [radians]
    
    #--------------------------------------
    # Compute "mean longitude of the sun"
    #--------------------------------------
    ## print '### Computing L...'
    L = (M + omega)   # [radians]
    
    #-----------------------------
    # Compute "equation of time"
    #-----------------------------
    TE = (-2.0 * e * np.sin(M)) + (np.sin(2 * L) * (eps / 2) ** 2.0)
    
    #------------------------------------
    # Convert from radians to degrees ?
    #------------------------------------
    if (DEGREES):    
        TE = TE * (np.float64(180) / np.pi)
        return TE    # [degrees]
    
    #----------------------------------------
    # Convert from radians to "decimal DMS"
    #----------------------------------------
    if (DMS):    
        TE = TE * (np.float64(180) / np.pi)  # [decimal degrees]
        deg = np.int16(TE)
        min = np.int16((TE - deg) * np.float64(60))
        sec = np.int16(((TE - deg) * np.float64(60) - min) * np.float64(60))
        TE = deg + (min / np.float64(100)) + (sec / np.float64(10000))
        return TE   # [decimal DMS, DD.MMSS]
    
    #----------------------------------------
    # Earth's rotation rate (angular speed)
    #----------------------------------------
    spin_rate = Earth_Rotation_Rate(PER_HOUR=True)    # [radians / hour]
    
    #--------------------------------
    # Convert from radians to hours
    #--------------------------------
    TE = (TE / spin_rate)   # [hours]
    return TE   # [hours]
    
#   Equation_of_Time()
#------------------------------------------------------------------------
def True_Solar_Noon( Julian_day, longitude, GMT_offset=None,
                     DST_offset=None, year=None):

    #------------------------------------------------------------
    # Notes: We need to know the local clock time when True
    #        Solar Noon occurs, since some of our equations
    #        depend on the time offset in hours from True Solar
    #        Noon. Note that TE may be negative or positive.

    #        The GMT_offset at the location of interest
    #        should be entered as an integer between 0 and 12,
    #        with negative values for locations west of the
    #        prime meridian.  LC is longitude correction and
    #        should be negative for longitudes east of the
    #        time zone's central meridian and positive other-
    #        wise, since solar noon will occur earlier for
    #        locations further to the east.  Note that some
    #        countries, like Iceland, may lie entirely outside
    #        of the time zone strip (i.e. the 15-degree wide
    #        strip of longitudes) that they set their clocks by.
    #-------------------------------------------------------------
    # NB!    Should we add or subtract TE below ??  **************
    #------------------------------------------------------------
    # NB!    The effect of Daylight Savings Time can be obtained
    #        by choosing an adjacent time zone number, or by
    #        using the optional DST_offset argument.  Be aware
    #        that different countries use different conventions.
    #------------------------------------------------------------
    # Notes: Added year keyword on 2021-07-28.
    #------------------------------------------------------------    
    if (year is None):
        year = Current_Year()

    time_zone_center_lon = GMT_offset * np.float64(15)  # [degrees]
    lon_diff = (time_zone_center_lon - longitude)    # [degrees]
    LC = (lon_diff / np.float64(15))      # [hours]
    ## print '### Computing TE...'
    TE = Equation_Of_Time(Julian_day, year=year)  # [hours]
    T_noon = np.float64(12) + LC + TE     # [hours; 24-hour military time]
    
    #-----------------------------
    # Add or subtract 1 hour for
    # Daylight Savings Time ?
    #-----------------------------
    if (DST_offset is not None):
        T_noon = T_noon + DST_offset
    
    return T_noon
    
#   True_Solar_Noon()
#------------------------------------------------------------------------
def Latitude_Grid( info ):

    #-----------------------------
    # Create a grid of latitudes
    #-----------------------------
    if (info.pixel_geom == 0):    
        #----------------------------------------
        # Geographic coords, fixed-angle pixels
        # Compute lats for pixel centers.
        #----------------------------------------
        dy = (info.yres / np.float64(3600))  #[arcsecs -> degrees]
        lats = (np.arange(info.nrows, dtype='float64') * dy) + info.y_south_edge + (dy/2)
        lats = np.flipud(lats)  ## Need FLIPUD vs. ROT90 to reverse 1D arrays.
        ones = np.ones([info.ncols], dtype='float64')
        lat_deg = np.outer( lats, ones )
        ## lat_deg = (transpose(matrixmultiply(transpose(ones), transpose(lats))))
        
        #** print,'min(lat_deg), max(lat_deg) = ', min(lat_deg), max(lat_deg)
        
    else:    
        #----------------------------------
        # UTM coords, fixed-length pixels
        #--------------------------------------------------
        # Must convert UTM coords to lats, which is a           ;***********
        # complicated procedure. Do this with RiverTools?
        #--------------------------------------------------
        lat_deg = np.float64(0)
        msg = np.array(['SORRY: ', ' ',
                     'The DEM for this data set uses UTM coordinates ',
                     'and TopoFlow cannot yet convert UTM coordinates ',
                     'to Geographic coordinates (lon and lat). ', ' ',
                     'A latitude value of 0.0 will be returned.', ' '])
        for line in msg:
            print(line)
        ## GUI_Error_Message(msg)

    #-----------------------------------------
    # Convert to 4-byte floats vs. doubles ?
    #-----------------------------------------
    ## lat_deg = np.float32( lat_deg )
    
    return lat_deg
    
#   Latitude_Grid()
#------------------------------------------------------------------------
def Longitude_Grid( info ):

    #------------------------------
    # Create a grid of longitudes
    #------------------------------
    if (info.pixel_geom == 0):    
        #----------------------------------------
        # Geographic coords, fixed-angle pixels
        # Compute lons for pixel centers.
        #----------------------------------------
        dx = (info.xres / np.float64(3600))  #[arcsecs -> degrees]
        lons = (np.arange(info.ncols, dtype='float64') * dx) + info.x_west_edge + (dx/2)
        ones = np.ones([info.nrows], dtype='float64')
        lon_deg = np.outer( ones, lons )
        ## lon_deg = (transpose(matrixmultiply(transpose(lons), transpose(ones))))
        
        #** print,'min(lon_deg), max(lon_deg) = ', min(lon_deg), max(lon_deg)
        
    else:    
        #----------------------------------
        # UTM coords, fixed-length pixels
        #--------------------------------------------------
        # Must convert UTM coords to lons, which is a           ;***********
        # complicated procedure. Do this with RiverTools?
        #--------------------------------------------------
        lon_deg = np.float64(0)
        msg = np.array(['SORRY: ', ' ',
                     'The DEM for this data set uses UTM coordinates ',
                     'and TopoFlow cannot yet convert UTM coordinates ',
                     'to Geographic coordinates (lon and lat). ', ' ',
                     'A longitude value of 0.0 will be returned.', ' '])
        for line in msg:
            print(line)
        ## GUI_Error_Message(msg)

    #-----------------------------------------
    # Convert to 4-byte floats vs. doubles ?
    #-----------------------------------------
    ## lon_deg = np.float32( lon_deg )
    
    return lon_deg
    
#   Longitude_Grid()
#------------------------------------------------------------------------


