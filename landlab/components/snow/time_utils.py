
#------------------------------------------------------------------------
#  Copyright (c) 2021, Scott D. Peckham
#
# Jul. 2021.  Changes to support MINT netCDF and datetimes.
#             Added get_time_dtype(), get_time_letter(),
#             get_duration() (from tf_utils), get_current_datetime().
# June 2020.  Created for TopoFlow calibration notebooks.
#             Started from code in: balto_gui.py and calibrate.py.
#
#------------------------------------------------------------------------

import datetime
from dateutil.relativedelta import relativedelta   ## 2022-02-18
import numpy as np
import sys

#------------------------------------------------------------------------
#  get_time_dtype()      # 2021-07-15
#  get_time_letter()     # 2021-07-22
#  get_duration()
#  get_current_datetime()
#  get_end_datetime()    # 2021-10-13
#
#  get_datetime_str()    # 2021-07-29
#  standardize_datetime_str()
#  get_datetime_obj_from_str()
#  get_datetime_obj_from_one_str()
#  pad_with_zeros()
#  split_datetime_str()
#  split_date_str()
#  split_time_str()
#  get_time_since_from_datetime()
#  get_datetime_from_time_since()
#  get_month_difference()
#
#  convert_times_from_hhmm_to_minutes()
#  convert_times_from_minutes_to_hhmm()
#  convert_times_from_datetime_to_minutes()
#  convert_times_from_minutes_to_datetime()
#
#  get_duration()
#
#  save_averaged_time_series()      ##### NOT FINISHED
#
#  COMMENTED OUT FOR NOW.
#  get_actual_time_units()
#  get_time_delta_str()
#  get_start_datetime_obj()
#  get_end_datetime_obj()
#  get_dt_from_datetime_str()
#
#--------------------------------------------------------------------
def get_time_dtype( time_units ):

    #----------------------------------------------------------
    # Note: Added on 2021-07-15 so that time vector in netCDF
    #       files can be returned as an array of datetimes.
    #       Callers: ncgs_files.py, ncts_files.py,
    #       ncps_files.py, nccs_files.py.
    #--------------------------------------------------------------
    # Numpy now supports a datetime64 "data type" that can be
    # used with np.arange() to create an array of datetimes.
    # The resolution is determined by a bracketed letter from:
    # Y=year, M=month, D=day, h=hour, m=minute, s=second
    # https://numpy.org/doc/stable/reference/arrays.datetime.html
    #--------------------------------------------------------------
    umap = {'seconds':'datetime64[s]', 'minutes':'datetime64[m]',
            'hours':'datetime64[h]',   'days':'datetime64[D]',
            'months':'datetime64[M]',  'years':'datetime64[Y]'}
    return umap[ time_units ]

#   get_time_dtype()
#--------------------------------------------------------------------
def get_time_letter( time_units ):
 
    umap = {'seconds':'s', 'minutes':'m', 'hours':'h',
            'days':'D', 'months':'M',  'years':'Y'}
    return umap[ time_units ]
 
#   get_time_letter()      
#--------------------------------------------------------------------
def get_duration(start_date=None, start_time=None,
                 end_date=None, end_time=None,
                 dur_units=None, REPORT=False):
                 ### time_step_secs=None):
                 
    #------------------------------------------------    
    # Note:  Compute time span between 2 datetimes.
    #------------------------------------------------
    # Next block is used for testing.
    #------------------------------------------------    
    if (start_date is None): start_date = '2014-01-01'
    if (start_time is None): start_time = '00:00:00'
    if (end_date is None):   end_date   = '2015-01-01'
    if (end_time is None):   end_time   = '00:00:00'
    if (dur_units is None):  dur_units  = 'days'
    #------------------------------------------------          
#     if (end_date is None):   end_date   = '2014-12-31'
#     if (end_time is None):   end_time   = '23:30:00'
#     if (time_step_secs is None):  time_step_secs  = '1800'


    date1 = start_date.split('-')
    y1 = int(date1[0])
    m1 = int(date1[1])  # NOTE:  int('08') = 8
    d1 = int(date1[2])
    #------------------------------  
    time1 = start_time.split(':')
    h1  = int(time1[0])
    mm1 = int(time1[1])
    s1  = int(time1[2])
    #------------------------------
    date2 = end_date.split('-')
    y2 = int(date2[0])
    m2 = int(date2[1])
    d2 = int(date2[2])
    #------------------------------   
    time2 = end_time.split(':')
    h2  = int(time2[0])
    mm2 = int(time2[1])
    s2  = int(time2[2])
    #-----------------------------------------------------------   
    start_obj     = datetime.datetime(y1, m1, d1, h1, mm1, s1)
    end_obj       = datetime.datetime(y2, m2, d2, h2, mm2, s2)
    duration_obj  = (end_obj - start_obj)
    duration_secs = duration_obj.total_seconds()
    
    #-------------------------------------------------------
    # If end_date is really the beginning of the last
    # interval, need to add time_step to get full duration
    #-------------------------------------------------------
    ## duration_secs += int(time_step_secs)

    #-----------------------------------------    
    # Convert duration to dur_units provided
    #-----------------------------------------
    if (dur_units == 'seconds'):
        duration = duration_secs
    elif (dur_units == 'minutes'):
        duration = (duration_secs / 60.0)
    elif (dur_units == 'hours'):
        duration = (duration_secs / 3600.0)
    elif (dur_units == 'days'):
        duration = (duration_secs / 86400.0)
    elif (dur_units == 'months'):
        #----------------------------------------------------
        # Different months have different number of seconds
        #----------------------------------------------------
        duration = get_month_difference(start_obj, end_obj)
    elif (dur_units == 'years'):
        duration = (duration_secs / 31536000.0)
    else:
        print('Unknown duration units = ' + dur_units + '.')
        print('Returning duration in hours.')
        duration = (duration_secs / 3600.0)
        
    if (REPORT):
        print( 'duration =', duration, '[' + dur_units + ']' )

    return duration
    
    #-----------------------------------------      
    # Alternate approach, where dur_units is
    # determined and then returned
    #-----------------------------------------   
#     if (duration_secs < 60):
#         duration  = duration_secs
#         dur_units = 'seconds'
#     elif (duration_secs < 3600):
#         duration  = divmod( duration_secs, 60 )[0]
#         dur_units = 'minutes'
#     elif (duration_secs < 86400):
#         duration  = divmod( duration_secs, 3600 )[0]
#         dur_units = 'hours'
#     elif (duration_secs <  31536000):          
#         duration = divmod( duration_secs, 86400 )[0]
#         dur_units = 'days'
#     else:
#         duration = divmod( duration_secs, 86400 )[0]
#         dur_units = 'days'
#               
#     return (duration, dur_units)
     
#   get_duration()   
#--------------------------------------------------------------------
def get_current_datetime( start_datetime, time, time_units ):
    
    start_datetime_obj = get_datetime_obj_from_one_str( start_datetime )
    #-----------------------------------------------------------
    # Modified next function to use dateutil.relativedelta vs.
    # datetime.timedelta to get support for months and years.
    # This was needed to write datetimes to netCDF when
    # time_units was equal to 'months'. (2022-02-18)
    #-----------------------------------------------------------
    datetime = get_datetime_from_time_since(start_datetime_obj,
                            time, units=time_units)
                            
    #-----------------------------------------------------
    # The return type is: <class 'datetime.datetime'>
    # Wrap result with "str()" before writing to netCDF.
    #-----------------------------------------------------
    return datetime

#   get_current_datetime()
#--------------------------------------------------------------------
def get_end_datetime( start_datetime, n_steps, time_units ):

    start_datetime_obj = get_datetime_obj_from_one_str( start_datetime )
    end_datetime_obj   = get_datetime_from_time_since(start_datetime_obj,
                             n_steps, units=time_units)
    return end_datetime_obj

#   get_end_datetime()
#---------------------------------------------------------------------
def get_datetime_str(y, m1, d, h, m2, s):

    m1_str = pad_with_zeros(m1, 2)
    d_str  = pad_with_zeros(d, 2)
    date_str = str(y) + '-' + m1_str + '-' + d_str
    #-----------------------------------------------
    h_str  = pad_with_zeros(h, 2)
    m2_str = pad_with_zeros(m2, 2)
    s_str  = pad_with_zeros(s, 2)
    time_str = h_str + ':' + m2_str + ':' + s_str
    #-----------------------------------------------    
    datetime_str = date_str + ' ' + time_str
    return datetime_str

#   get_datetime_str()
#---------------------------------------------------------------------
def standardize_datetime_str( datetime_str ):

	#---------------------------------------------------
	# Note: Handle any delim between date & time, such
	#       as '', ' ' or 'T' (from Ankush).
	#---------------------------------------------------  
    date_str     = datetime_str[0:10]   # (2015-10-01)
    time_str     = datetime_str[-8:]    # (00:00:00)
    datetime_str = date_str + ' ' + time_str
    return datetime_str
    
#   standardize_datetime_str()
#--------------------------------------------------------------------
def get_datetime_obj_from_str( date_str, time_str='00:00:00' ):

    #---------------------------------------------------
    # date_str = 'YYYY-MM-DD', time_str = 'HH:MM:SS'
    #---------------------------------------------------
    ## e.g. d1 = str(self.datetime_end_date.value)
    ## e.g. t1 = self.datetime_end_time.value

    (y, m1, d) = split_date_str(date_str)
    (h, m2, s) = split_time_str(time_str)
    if( y <= 0 ):
        # msg  = 'Year cannot be < 1 in start date.\n'
        # msg += 'Changed year from ' + str(y) + ' to 1.'
        # self.datetime_notes.value = msg
        print('Year cannot be < 1 in start date.')
        print('Changed year from ' + str(y) + ' to 1.')
        print()
        y = 1
    datetime_obj = datetime.datetime(y, m1, d, h, m2, s) 
    return datetime_obj
    
#   get_datetime_obj_from_str()
#--------------------------------------------------------------------                       
def get_datetime_obj_from_one_str( datetime_str ):

    (date, time) = split_datetime_str( datetime_str )
    (y, m1,  d)  = split_date_str( date )
    (h, m2, s)   = split_time_str( time )
    datetime_obj = datetime.datetime(y, m1, d, h, m2, s)
    return datetime_obj

#   get_datetime_obj_from_one_str()
#------------------------------------------------------------------------
def pad_with_zeros(num, target_len):
  
    num_string = str( int(num) )  # int removes decimal part
    n = len( num_string )
    m = (target_len - n)
    num_string = ('0'*m) + num_string
    return num_string

#   pad_with_zeros()
#--------------------------------------------------------------------  
def split_datetime_str(datetime_obj, datetime_sep=' ',
                       ALL=False):

    #-----------------------------------------------     
    # Note: Still works if datetime_obj is string.
    #-----------------------------------------------
    datetime_str = str(datetime_obj)
    parts = datetime_str.split( datetime_sep )
    ## print('## datetime_str =', datetime_str )
    ## print('## parts =', str(parts) )
    
    date_str = parts[0]
    time_str = parts[1]
    if not(ALL):
        return (date_str, time_str)
    else:
        (y,m1,d) = split_date_str( date_str )
        (h,m2,s) = split_time_str( time_str )
        return (y,m1,d,h,m2,s)

#   split_datetime_str()
#--------------------------------------------------------------------  
def split_date_str(date_str, date_sep='-'):

    date_parts = date_str.split( date_sep )
    year  = int(date_parts[0])
    month = int(date_parts[1])   # NOTE:  int('08') = 8
    day   = int(date_parts[2])
    return (year, month, day)
     
#   split_date_str()
#--------------------------------------------------------------------  
def split_time_str(time_str, time_sep=':'):

    time_parts = time_str.split( time_sep )
    hour   = int(time_parts[0])
    minute = int(time_parts[1])
    second = int(time_parts[2])
    return (hour, minute, second)

#   split_time_str()
#--------------------------------------------------------------------                       
def get_time_since_from_datetime(origin_datetime_obj,
                                 datetime_obj, units='days'):

    #-------------------------------------------------
    # Compute time duration between datetime objects
    #-------------------------------------------------
    origin_obj    = origin_datetime_obj
    duration_obj  = (datetime_obj - origin_obj)
    duration_secs = duration_obj.total_seconds()
    #---------------------------------------------------
    # There is not a fixed number of seconds per month
    # Also 52 (weeks/year) * 7 (days/week) = 364.
    #---------------------------------------------------
    secs_per_unit_map = {
    'years':31536000.0, 'weeks':604800.0, 'days':86400.0,
    'hours':3600.0, 'minutes':60.0, 'seconds':1 }          
    secs_per_unit = secs_per_unit_map[ units ]      
    duration = (duration_secs / secs_per_unit )
    time_since = duration  # (in units provided)

    return time_since
        
#   get_time_since_from_datetime()
#--------------------------------------------------------------------                       
def get_datetime_from_time_since(origin_datetime_obj,
                           time_since, units='days'):

    # For testing
#         print('## type(times_since) =', type(time_since) )
#         print('## time_since =', time_since )
#         print('## int(time_since) =', int(time_since) )

    #---------------------------------------------------   
    # Note: Use: dateutil.relativedelta.relativedelta:
    #       https://dateutil.readthedocs.io/en/stable/
    #-----------------------------------------------------------
    # Note: Unlike datetime.timedelta, dateutil.relativedelta
    #       accepts numpy types (e.g. np.int16, np.float32),
    #       and supports weeks, months, & years.  A month is
    #       not a fixed length of time, unlike most other
    #       date/time units.
    #-----------------------------------------------------------
    delta = None
    time_since2 = float(time_since)  # (may not be needed)
    #------------------------------------------------------
    if (units == 'years'):
        delta = relativedelta( years=time_since2 )  
    if (units == 'months'):
        delta = relativedelta( months=time_since2 )  
    if (units == 'weeks'):
        delta = relativedelta( weeks=time_since2 )    
    if (units == 'days'):
        delta = relativedelta( days=time_since2 )
    if (units == 'hours'):
        delta = relativedelta( hours=time_since2 )
    if (units == 'minutes'):
        delta = relativedelta( minutes=time_since2 )
    if (units == 'seconds'):
        delta = relativedelta( seconds=time_since2 )
    #-----------------------------------------------------------
    # This should also work.
    # delta = eval("relativedelta(" + units + "=timesince2""))
    #-----------------------------------------------------------
    if (delta is None):
        msg = '### ERROR: Units: ' + units + ' not supported.'
        return origin_obj
        
    #---------------------------------------------------   
    # Note: datetime.timedelta() can take integer or
    #       float arguments, and the arguments can be
    #       very large numbers.  However, it does not
    #       accept any numpy types, whether float or
    #       int (e.g. np.int16, np.float32).
    #  https://docs.python.org/3/library/datetime.html
    #---------------------------------------------------
    # Note: This routine doesn't work for 'months' &
    #       'years', but see get_current_datetime()
    #---------------------------------------------------
#     delta = None
#     time_since2 = float(time_since)  ## No numpy types
#     #------------------------------------------------------    
#     if (units == 'days'):
#         delta = datetime.timedelta( days=time_since2 )
#     if (units == 'hours'):
#         delta = datetime.timedelta( hours=time_since2 )
#     if (units == 'minutes'):
#         delta = datetime.timedelta( minutes=time_since2 )
#     if (units == 'seconds'):
#         delta = datetime.timedelta( seconds=time_since2 )
#     #------------------------------------------------------
#     if (delta is None):
#         msg = '### ERROR: Units: ' + units + ' not supported.'
#         return origin_obj

    # For testing
    ## print('#### delta =', delta)
    
    #---------------------------------------------        
    # Create new datetime object from time_since
    # This is used by: get_current_datetime().
    #---------------------------------------------
    # The return type is:
    # <class 'datetime.datetime'>
    #---------------------------------------------    
    origin_obj = origin_datetime_obj
    new_dt_obj = (origin_obj + delta)
    return new_dt_obj

#   get_datetime_from_time_since()
#--------------------------------------------------------------------  
def get_month_difference(start_datetime_obj, end_datetime_obj ):

    #-------------------------------------------
    # Example 0: 2017-09 to 2017-09
    # months = (2017-2017)*12 = 0
    # months = (months - 9) = (0-9) = -0
    # months = (months + 9) = 0   (as index)
    #-------------------------------------------           
    # Example 1: 2017-09 to 2018-02
    # 9:10, 10:11, 11:12, 12:1, 1:2 = 5 (if same days)
    # months = (2018-2017)*12  = 12
    # months = (months - 9) = 3
    # months = (months + 2) = 3 + 2 = 5
    #-------------------------------------------
    start_year = start_datetime_obj.year
    end_year   = end_datetime_obj.year
    months     = (end_year - start_year) * 12
    #-------------------------------------------
    start_month = start_datetime_obj.month
    end_month   = end_datetime_obj.month
    months = months - start_month
    months = months + end_month
    ## months = months + 1  # (no: get 1 if dates same)
    ## print('month difference =', months)
    return float(months)   # see get_duration() re: float()

#   get_month_difference()
#--------------------------------------------------------------------
#--------------------------------------------------------------------
def convert_times_from_hhmm_to_minutes( times_hhmm ):

    #----------------------------------------------------      
    # Notes:  For simple time strings with only 'hhmm',
    #         return the time in minutes.
    #----------------------------------------------------    
    n_times   = len( times_hhmm )
    times_min = np.zeros( n_times )
    for k in range(n_times):
        hhmm = times_hhmm[k]
        hour = np.int16( hhmm[:2] )
        min  = np.int16( hhmm[2:] )
        times_min[k] = (hour * 60) + min

    return times_min

#   convert_times_from_hhmm_to_minutes()
#---------------------------------------------------------------------
def convert_times_from_minutes_to_hhmm( times_min ):

    n_times   = len( times_min )
    times_hhmm = np.zeros( n_times )
    for k in range(n_times):
        hour = int( times_min[k] / 60 )
        min  = int( times_min[k] % 60 )
        #-----------------------------------        
        hh   = str(hour)
        if (len(hh) == 1):  hh = ('0' + hh)
        #-----------------------------------
        mm   = str(min)
        if (len(mm) == 1):  mm = ('0' + mm)
        #-----------------------------------
        times_hhmm[k] = int(hh + mm)

    return times_hhmm

#   convert_times_from_minutes_to_hhmm()
#---------------------------------------------------------------------
def convert_times_from_datetime_to_minutes( times_datetime,
                  origin_datetime_obj=None ):
    
    if (origin_datetime_obj is None):
        print('ERROR: origin_datetime_obj is required.')
        return 

    n_times   = len( times_datetime )
    times_min = np.zeros( n_times )
    for k in range(n_times):
        datetime_str = times_datetime[k]
        #---------------------------------------------------
        # Note: Handle any delim between date & time, such
        #       as '', ' ' or 'T' (from Ankush).
        #--------------------------------------------------- 
        datetime_str = standardize_datetime_str( datetime_str )        
        datetime_obj = get_datetime_obj_from_one_str( datetime_str )
        times_min[k] = get_time_since_from_datetime(origin_datetime_obj,
                                datetime_obj, units='minutes')

    return times_min

#   convert_times_from_datetime_to_minutes()
#---------------------------------------------------------------------
def convert_times_from_minutes_to_datetime( times_min,
                  origin_datetime_obj=None ):
    
    if (origin_datetime_obj is None):
        print('ERROR: origin_datetime_obj is required.')
        return 

    n_times   = len( times_min )
    times_datetime = np.zeros( n_times, dtype='U30' )   # string array
    for k in range(n_times):
        time_since = times_min[k]
        times_datetime[k] = get_datetime_from_time_since(origin_datetime_obj,
                                         time_since, units='minutes')

    return times_datetime

#   convert_times_from_minutes_to_datetime()
#-------------------------------------------------------------------
#---------------------------------------------------------------------
# def save_averaged_time_series( values, n_days=1, PLOT=False,
#                   start_datetime='2015-10-01 00:00:00'):
# 
#   THIS IS NOT FINISHED YET.
#
# 	#------------------------------------
# 	# Convert hourly time series values
# 	# to n-day average values
# 	#------------------------------------
# 	## start_date = '2015-10-01'
# 	# n_days  = 12
# 
# 	#------------------------	
# 	# Hours to average over
# 	#------------------------
# 	n_hours = 24 * n_days
# 	n_avg   = np.int32(values.size / n_hours)
# 	v_avg   = np.zeros(n_avg, dtype='float32')
# 	t_avg   = np.zeros(n_avg, dtype='float32')    ###########
# 	i1 = 0
# 	i2 = n_hours
# 
# 	hours = 0
# 
# 	for k in range(n_avg):
# 		v_avg[k] = (values[i1:i2].sum() / n_hours)
# 		t_avg[k] = ''
# 		hours += n_hours
# 		i1 += n_hours
# 		i2 += n_hours
# 
#     #-----------------------------------------
#     # Write time series to multi-column text
#     #-----------------------------------------
#     # from topoflow.utils import text_ts_files
#     
#     #-------------------------------------
#     # Option to plot the new time series
#     #-------------------------------------
#     if (PLOT):
#         from topoflow.utils import visualize as tfvis
# 	    tfvis.plot_data(t_avg, v_avg)
# 
# #   save_averaged_time_series()
#---------------------------------------------------------------------
#     def get_actual_time_units(self):
# 
# #         secs_per_unit_list = [1, 60.0, 3600.0, 86400, 31536000.0, -1]
# #         next_unit_factor = [60.0, 60.0, 24.0, 365.0, -1, -1]
# 
#         units_list = ['second', 'minute', 'hour',
#                        'day', 'year', 'None']   # ascending, skip month
# 
#         for units in units_list:
#              if (self.time_units_str.startswith(units)):
#                  break
#         if (units != None):
#             units += 's'   # (make units plural now; not before)
#         else:
#             print('ERROR: No match found for units.')
#             return
#         self.time_units = units
# 
#     #   get_actual_time_units()
#     #--------------------------------------------------------------------
#     def get_time_delta_str(self):
# 
#         ## print('### self.time_var.size =', self.time_var.size )
#         ## print('###')
#         
#         #-----------------------------------
#         # Check size of the time_var array
#         #-----------------------------------
#         if (self.time_var.size == 1):
#             dt = 0
#             self.time_delta = '0000-00-00 00:00:00'
#             # print('At top of get_time_delta_str():')
#             # print('self.time_var.size =', self.time_var.size )
#             # print('self.time_delta =', self.time_delta )
#             return
#         if (self.time_var.size > 1):  
#             dt  = (self.time_var[1] - self.time_var[0])
#             print('dt1 =', dt)
#         if (self.time_var.size > 3):
#             dt2 = (self.time_var[2] - self.time_var[1])  ###
#             dt3 = (self.time_var[3] - self.time_var[2])  ###
#             print('dt2 =', dt2)  # check if evenly spaced
#             print('dt3 =', dt3)
#                 
#         #---------------------------------------------------        
#         # Note: Actual time units were stripped from units
#         #       string and saved as self.time_units.
#         #       A full units attribute string may be:
#         #        'hour since 0000-00-00 00:00:00'
#         #---------------------------------------------------
#         units_list = ['seconds', 'minutes', 'hours',
#                       'days', 'years', 'None']  # ascending, skip month
#         secs_per_unit_list = [1, 60.0, 3600.0, 86400, 31536000.0, -1]
#         next_unit_factor   = [60.0, 60.0, 24.0, 365.0, -1, -1]
#         units       = self.time_units
#         units_index = units_list.index( units )
#         #----------------------------------------
#         if (units == 'years'):
#             s = self.pad_with_zeros(dt,4)
#         else:
#             if (len(str(dt)) <= 2):
#                 s = self.pad_with_zeros(dt,2)
#             else:
#                 #-------------------------------
#                 # Must convert units to get dt
#                 # down to 1 or 2 digits.
#                 #-------------------------------
#                 old_dt    = dt
#                 old_units = units
#                 k = units_index
#                 n = len( str(int(dt)) )
#                 while (n > 2) and (units != 'None'):
#                     k     = k + 1
#                     dt    = (dt / next_unit_factor[k-1])
#                     units = units_list[k]
#                     n     = len( str(int(dt)) )
#                 if (units == 'None'):
#                     print('#####################################')
#                     print('ERROR in get_time_delta_str():')
#                     print('      dt has too many digits.')
#                     print('#####################################')
#                     return
#                 else:
#                     # Note that any remainder has been dropped.
#                     s = self.pad_with_zeros(dt,2)
#                     print('Old dt and units =', old_dt, old_units)
#                     print('New dt and units =', dt, units)
#                     print('Remainder not retained yet.')
#         #----------------------------------------------
#         if (units == 'years'):
#             td = (s + '-00-00 00:00:00')
# #         if (units == 'months'):
# #             td= ('0000-' + s + '-00 00:00:00')
#         if (units == 'days'):
#             td = ('0000-00-' + s + ' 00:00:00')
#         if (units == 'hours'):
#             td = ('0000-00-00 ' + s + ':00:00')
#         if (units == 'minutes'):
#             td = ('0000-00-00 00:' + s + ':00')
#         if (units == 'seconds'):
#             td = ('0000-00-00 00:00:' + s)
#         #------------------------------------------------
#         self.time_delta = td
#         # print('At bottom of get_time_delta_str():')
#         # print('self.time_delta =', td)
#         # print()
# 
#     #   get_time_delta_str()
#     #--------------------------------------------------------------------
#     def get_start_datetime_obj(self):
# 
#         #---------------------------------------
#         # d1.value is a datetime "date object"
#         # t1.value is a time string: 00:00:00
#         #---------------------------------------
#         d1 = self.datetime_start_date
#         t1 = self.datetime_start_time
#         if (d1.value is None):
#             return None
#         date_str = str(d1.value)
#         time_str = t1.value   # (already string)
#         ## print('In get_start_datetime_obj():')
#         ## print('date_str =', date_str)
#         ## print('time_str =', time_str)
#         
#         datetime_obj = self.get_datetime_obj_from_str(date_str, time_str)
#         return datetime_obj
#     
#     #   get_start_datetime_obj()
#     #--------------------------------------------------------------------
#     def get_end_datetime_obj(self):
# 
#         #---------------------------------------
#         # d1.value is a datetime "date object"
#         # t1.value is a time string: 00:00:00
#         #---------------------------------------
#         d1 = self.datetime_end_date
#         t1 = self.datetime_end_time
#         if (d1.value is None):
#             return None
#         date_str = str(d1.value)
#         time_str = t1.value   # (already string)
#         ## print('In get_end_datetime_obj():')
#         ## print('date_str =', date_str)
#         ## print('time_str =', time_str)
#         
#         datetime_obj = self.get_datetime_obj_from_str(date_str, time_str)
#         return datetime_obj
# 
#     #   get_end_datetime_obj()
#---------------------------------------------------------------------

