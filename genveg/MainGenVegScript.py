"""
GenVeg is a generalized vegetation growth model that simulates growth due to
photosynthesis as well as processes that control vegetation distribution such
as mortality, senescence, dispersal, and colonization. The model utilizes a 
mixture of existing and new landlab components to setup the model, simulate
relevant vegetation community dynamics with the local environment, and analyze
results of the simulation. 
"""

"""
Model control setup
"""
##Import components
import pathlib
from landlab import RasterModelGrid
from landlab.data_record import DataRecord
from landlab.components.genveg import VegParams
import numpy as np
import pandas as pd
from datetime import date
import pathlib

##USER INPUTS
##Set up model simulation parameters
#Time control
sim_start_date = date(2010,1,1)
tot_sim_len_y = 30
veg_time_step_d = 1
env_time_step_d = 1

#Filename containing vegetation parameters
VegInputFile = 'GenVeg_params_inputs_1col.xlsx'

#Filename of environmental driver inputs
#temperature and latitude

## Set up model grid
#Define grid size - can import these parameters from ESRI ASCII
#For testing, assume 100x100 cell grid with 2-m grid cells
pg = RasterModelGrid((100, 100), 2.0)
pg_size = pg.number_of_cell_columns*pg.number_of_cell_rows

#Define data stored in grid
#Initialize biomass to zero
_ = pg.add_zeros('vegetation__live_biomass',at='cell', units='g')
#Initialize plant functional type as grass or bare
_ = pg.add_field('vegetation__plant_functional_type',np.random.choice([0,3],pg_size), at='cell')

#Define environmental drivers
_ = pg.add_field('temperature__mean_air_temp', np.random.uniform(low=0, high=55,size=pg.number_of_nodes), at='node', units='C' )

## Set up time loop list using real dates
#Define model end and timesteps
sim_end = date(sim_start_date.year+tot_sim_len_y,sim_start_date.month,sim_start_date.day)
dt = min([veg_time_step_d,env_time_step_d])
day_freq = str(dt)+'D'

#Construct model run date list, list of steps, and zip to tuple of step number and dates
date_counter = pd.date_range(start=sim_start_date,end=sim_end,freq=day_freq).strftime('%Y-%m-%d').tolist()
step_counter=range(len(date_counter))
model_counter=tuple(zip(step_counter,date_counter))

#Read in vegetation parameter data sets
fdir=pathlib.Path(__file__).parent
inp_file=pathlib.Path(fdir,VegInputFile)
print(inp_file)

data = VegParams(fpath=inp_file,processes=['plantsize','dispersal','mortality','colonization','storage'])  # Uses GenVeg function to read formatted Excel File

#Read in environmental time series data
#Here we will read in the vegetation driver data sets
pass

#Start model loop here using vegetation time step
#Will need to change structure for 
for i in model_counter:
    # # Update objects
    # Env driver check

    # Calculate Day of Year (DOY)
    Julian = np.int(np.floor((current_time - np.floor(current_time)) * 365.))
    # Generate seasonal storms
    # for Dry season
    if Julian < data['doy__start_of_monsoon'] or Julian > data[
            'doy__end_of_monsoon']:
        PD_D.update()
        P[i] = PD_D.get_storm_depth()
        Tr[i] = PD_D.get_precipitation_event_duration()
        Tb[i] = PD_D.get_interstorm_event_duration()
    # Wet Season—Jul to Sep—NA Monsoon
    else:
        PD_W.update()
        P[i] = PD_W.get_storm_depth()
        Tr[i] = PD_W.get_precipitation_event_duration()
        Tb[i] = PD_W.get_interstorm_event_duration()

    # Spatially distribute PET and its 30-day-mean (analogous to degree day)
    grid['cell']['surface__potential_evapotranspiration_rate'] = (
        (np.choose(grid['cell']['vegetation__plant_functional_type'],
                   PET_[Julian])) * Rad_Factor[Julian])
    grid['cell']['surface__potential_evapotranspiration_30day_mean'] = (
        (np.choose(grid['cell']['vegetation__plant_functional_type'],
                   EP30[Julian])) * Rad_Factor[Julian])

    # Assign spatial rainfall data
    grid['cell']['rainfall__daily_depth'] = P[i] * np.ones(
        grid.number_of_cells)
    
    #Update Tb and Tr in SM and VEG
    VEG.Tb, VEG.Tr = Tb[i], Tr[i]
    SM.Tb, SM.Tr = Tb[i], Tr[i]

    # Update soil moisture component
    current_time = SM.update()

    # Decide whether its growing season or not
    if Julian != 364:
        if EP30[Julian + 1, 0] > EP30[Julian, 0]:
            PET_threshold = 1
            # 1 corresponds to ETThresholdup (begin growing season)
        else:
            PET_threshold = 0
            # 0 corresponds to ETThresholddown (end growing season)

    # Update vegetation component
    #VEG.update(PETthreshold_switch=PET_threshold, Tb=Tb[i], Tr=Tr[i])
    VEG.update()

    # Update yearly cumulative water stress data
    WS += (grid['cell']['vegetation__water_stress']) * Tb[i] / 24.

    # Record time (optional)
    Time[i] = current_time

    # Cellular Automata
    if (current_time - time_check) >= 1.:
        if yrs % 5 == 0:
            print('Elapsed time = {time} years'.format(time=yrs))
        VegType[yrs] = grid['cell']['vegetation__plant_functional_type']
        grid['cell']['vegetation__cumulative_water_stress'] = WS / Tg
        vegca.update()
        SM.initialize()
        VEG.initialize()
        time_check = current_time
        WS = 0
        yrs += 1
VegType[yrs] = grid['cell']['vegetation__plant_functional_type']


# Calculate approximate number of storms per year - hold over from CatGRASS
fraction_wet = (data['doy__end_of_monsoon'] -
                data['doy__start_of_monsoon']) / 365.
fraction_dry = 1 - fraction_wet
no_of_storms_wet = (8760 * (fraction_wet) /
                    (data['mean_interstorm_wet'] + data['mean_storm_wet']))
no_of_storms_dry = (8760 * (fraction_dry) /
                    (data['mean_interstorm_dry'] + data['mean_storm_dry']))
n = int(n_years * (no_of_storms_wet + no_of_storms_dry))