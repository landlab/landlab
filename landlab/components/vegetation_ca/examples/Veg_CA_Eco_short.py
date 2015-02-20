
# Authors: Sai Nudurupati & Erkan Istanbulluoglu, 12Dec14

# Import required libraries

import os
import numpy as np
import matplotlib.pyplot as plt
from landlab.io import read_esri_ascii
from landlab import RasterModelGrid as rmg
from landlab.components.uniform_precip.generate_uniform_precip \
                        import PrecipitationDistribution
from landlab.components.radiation.radiation_field import Radiation
from landlab.components.PET.potential_evapotranspiration_field \
                        import PotentialEvapotranspiration
from landlab.components.soil_moisture.soil_moisture_multi_pft_new  \
                        import SoilMoisture
from landlab.components.single_vegetation.vegetation_multi_pft_new \
                        import Vegetation
from landlab.components.vegetation_ca.CA_Veg_new  import VegCA
import matplotlib as mpl
import time

# GRASS = 0; SHRUB = 1; TREE = 2; BARE = 3;
# SHRUBSEEDLING = 4; TREESEEDLING = 5

## Function that converts text file to a dictionary
def txt_data_dict( InputFile ):
    f = open( InputFile )
    data1 = {}
    for line in f:
       if line.strip() != '' and line[0] != '#':
           m, n = line.split(':')
           line = f.next()
           e = line[:].strip()
           if e[0].isdigit():
               if e.find('.') != -1:
                   data1[m.strip()] = float(line[:].strip())
               else:
                   data1[m.strip()] = int(line[:].strip())
           else:
               data1[m.strip()] = line[:].strip()
    f.close()
    return data1.copy()


## Point to the input DEM
_DEFAULT_INPUT_FILE_1 = os.path.join(os.path.dirname(__file__),
                                 'DEM_10m.asc')

InputFile = 'Inputs_Vegetation_CA.txt'
data = txt_data_dict( InputFile ) # Create dictionary that holds the inputs

## Importing Grid and Elevations from DEM
(grid1,elevation) = read_esri_ascii(_DEFAULT_INPUT_FILE_1)
grid1['node']['Elevation'] = elevation
grid = rmg(5,4,5)
grid['node']['Elevation'] = 1700. * np.ones(grid.number_of_nodes)
grid1['cell']['VegetationType'] = np.random.randint(0,6,grid1.number_of_cells)
grid['cell']['VegetationType'] = np.arange(0,6)

# Create radiation, soil moisture and Vegetation objects
PD_D = PrecipitationDistribution(mean_storm = data['mean_storm_dry'],  \
                    mean_interstorm = data['mean_interstorm_dry'],
                    mean_storm_depth = data['mean_storm_depth_dry'])
PD_W = PrecipitationDistribution(mean_storm = data['mean_storm_wet'],  \
                    mean_interstorm = data['mean_interstorm_wet'],
                    mean_storm_depth = data['mean_storm_depth_wet'])
Rad = Radiation( grid )
PET_Tree = PotentialEvapotranspiration( grid1, method = data['PET_method'], \
                    MeanTmaxF = data['MeanTmaxF_tree'],
                    DeltaD = data['DeltaD'] )
PET_Shrub = PotentialEvapotranspiration( grid1, method = data['PET_method'], \
                    MeanTmaxF = data['MeanTmaxF_shrub'],
                    DeltaD = data['DeltaD'] )
PET_Grass = PotentialEvapotranspiration( grid1, method = data['PET_method'], \
                    MeanTmaxF = data['MeanTmaxF_grass'],
                    DeltaD = data['DeltaD'] )

SM = SoilMoisture( grid, data )   # Soil Moisture object
VEG = Vegetation( grid, data )    # Vegetation object
vegca = VegCA( grid1, data )      # Cellular automaton object

##########
n = data['n_short']   # Defining number of storms the model will be run
##########

## Create arrays to store modeled data
P = np.empty(n)    # Record precipitation
Tb = np.empty(n)    # Record inter storm duration
Tr = np.empty(n)    # Record storm duration
Time = np.empty(n) # To record time elapsed from the start of simulation

CumWaterStress = np.empty([n/50, grid1.number_of_cells]) # Cum Water Stress
CumWS = np.empty([n/50, grid.number_of_cells]) # Cum Water Stress of different plant types
VegType = np.empty([n/50, grid1.number_of_cells],dtype = int)
PETn_ = np.empty([n, grid.number_of_cells])
EP30_ = np.empty([n, grid.number_of_cells])
AET_ = np.empty([n, grid.number_of_cells])
SM_ = np.empty([n, grid.number_of_cells])
SMini_ = np.empty([n, grid.number_of_cells])
SMavg_ = np.empty([n, grid.number_of_cells])
LAI = np.empty([n, grid.number_of_cells])
WaterStress = np.empty([n, grid.number_of_cells])
ETmax = np.empty([n, grid.number_of_cells])
LAI_fr = np.empty([n, grid.number_of_cells])
MAP = np.zeros(n)
Julian_ = np.empty(n)

PET_ = np.zeros([365, grid.number_of_cells])
Rad_Factor = np.empty([365,grid.number_of_cells])
EP30 = np.empty([365, grid.number_of_cells]) # 30 day average PET to determine season
PET_threshold = 0  # Initializing PET_threshold to ETThresholddown

## Initializing inputs for Soil Moisture object
grid['cell']['LiveLeafAreaIndex'] = 1.6 * np.ones( grid.number_of_cells )
SM._SO = 0.59 * np.ones(grid.number_of_cells) # Initializing Soil Moisture
Tg = 365    # Growing season in days


## Calculate current time in years such that Julian time can be calculated
current_time = 0                   # Start from first day of June

for i in range(0,365):
    Rad.update( float(i)/365.25)
    PET_Tree.update( float(i)/365.25 )
    PET_Shrub.update( float(i)/365.25 )
    PET_Grass.update( float(i)/365.25 )
    PET_[i] = [PET_Grass._PET_value, PET_Shrub._PET_value, PET_Tree._PET_value,
                0., PET_Shrub._PET_value, PET_Tree._PET_value]
    Rad_Factor[i] = grid['cell']['RadiationFactor']
    if i < 30:
        if i == 0:
            EP30[0] = PET_[0]
        else:
            EP30[i] = np.mean(PET_[:i],axis = 0)
    else:
        EP30[i] = np.mean(PET_[i-30:i], axis = 0)


time_check = 0.               # Buffer to store current_time at previous storm
i_check = 0                   #
yrs = 0

Start_time = time.clock()     # Recording time taken for simulation
WS = 0.

## Run Time Loop
for i in range(0, n):
    # Update objects
    Julian = np.floor( (current_time - np.floor( current_time)) * 365.)
    if Julian < 182 or Julian > 273:  # Dry Season
        PD_D.update()
        P[i] = PD_D.storm_depth
        Tr[i] = PD_D.storm_duration
        Tb[i] = PD_D.interstorm_duration
    else:                             # Wet Season - Jul to Sep - NA Monsoon
        PD_W.update()
        P[i] = PD_W.storm_depth
        Tr[i] = PD_W.storm_duration
        Tb[i] = PD_W.interstorm_duration

    grid['cell']['PotentialEvapotranspiration'] = PET_[Julian]
    grid['cell']['PotentialEvapotranspiration30'] = EP30[Julian]
    current_time = SM.update( current_time, P = P[i], Tr = Tr[i], Tb = Tb[i] )
    PETn_[i] = SM._PET

    if Julian != 364:
        if EP30[Julian + 1,0] > EP30[Julian,0]:
            PET_threshold = 1  # 1 corresponds to ETThresholdup
        else:
            PET_threshold = 0  # 0 corresponds to ETThresholddown

    VEG.update(PotentialEvapotranspirationThreshold = PET_threshold)

    WS += (grid['cell']['WaterStress'])*Tb[i]/24.
    PETn_[i] = grid['cell']['PotentialEvapotranspiration']
    AET_[i] = grid['cell']['ActualEvapotranspiration']
    SM_[i] = grid['cell']['SaturationFraction']
    SMini_[i] = SM._Sini
    SMavg_[i] = (SMini_[i]+SM_[i])/2.
    Time[i] = current_time
    WaterStress[i] = grid['cell']['WaterStress']
    MAP[yrs] += P[i]
    Julian_[i] = Julian
    LAI[i] = grid['cell']['LiveLeafAreaIndex']
    ETmax[i] = SM._ETmax
    LAI_fr[i] = SM._fr

    # Cellular Automata
    if (current_time - time_check) >= 1.:
        VegType[yrs] = grid1['cell']['VegetationType']
        WS_ = np.choose(VegType[yrs],WS)
        CumWaterStress[yrs] = WS_/Tg
        grid1['cell']['CumulativeWaterStress'] = CumWaterStress[yrs]
        vegca.update()
        CumWS[yrs] = WS/Tg   # Record of Cumulative Water Stress of different plant types
        time_check = current_time
        WS = 0
        i_check = i
        yrs += 1


Final_time = time.clock()

VegType[yrs] = grid1['cell']['VegetationType']
Time_Consumed = (Final_time - Start_time)/60.    # in minutes

## Evaluating Cumulative WaterStress of each plant functional type
CWS = CumWS.round(decimals=2)
# Separate arrays according to their PFTs
CWS_G = CWS[0:yrs,0]
CWS_S = CWS[0:yrs,1]
CWS_T = CWS[0:yrs,2]
CWS_SS = CWS[0:yrs,4]
CWS_TS = CWS[0:yrs,5]
# Sort them in ascending order
CWS_G.sort()
CWS_S.sort()
CWS_T.sort()
CWS_SS.sort()
CWS_TS.sort()
# Weibull probability of exceedance (Ref. Erkan's class notes - Hw1 Physical Hydrology)
Q_G = [(1. - (m/(yrs+1.))) for m in range(1,len(CWS_G)+1)]
Q_S = [(1. - (m/(yrs+1.))) for m in range(1,len(CWS_S)+1)]
Q_T = [(1. - (m/(yrs+1.))) for m in range(1,len(CWS_T)+1)]
Q_SS = [(1. - (m/(yrs+1.))) for m in range(1,len(CWS_SS)+1)]
Q_TS = [(1. - (m/(yrs+1.))) for m in range(1,len(CWS_TS)+1)]


## Plotting
cmap = mpl.colors.ListedColormap(   \
                    [ 'green', 'red', 'black', 'white', 'red', 'black' ] )
bounds = [-0.5,0.5,1.5,2.5,3.5,4.5,5.5]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

# Grass
pic = 0
plt.figure(pic)
plt.subplot(5,1,1)
plt.bar(Time,P,width = 0.01)
plt.xlim([0, Time[n-1]])
plt.ylabel('Rainfall (mm)')
plt.title('Grass')
plt.subplot(5,1,2)
plt.plot(Time, PETn_[:,0])
plt.ylabel('PET (mm)')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,3)
plt.plot(Time, AET_[:,0])
plt.ylabel('AET (mm)')
plt.xlim([0, Time[n-1]])
plt.ylim([0, max(AET_[:,0])])
plt.subplot(5,1,4)
plt.plot(Time, SMavg_[:,0])
plt.ylabel('S - avg')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,5)
plt.plot(np.arange(1,yrs+1),CumWS[0:yrs,0])
plt.ylabel('Cum WS')
plt.xlim([0, yrs])
plt.xlabel('Time in Years')

pic += 1
plt.figure(pic)
plt.subplot(5,1,1)
plt.bar(Time,P,width = 0.01)
plt.xlim([0, Time[n-1]])
plt.ylabel('Rainfall (mm)')
plt.title('Shrub')
plt.subplot(5,1,2)
plt.plot(Time, PETn_[:,1])
plt.ylabel('PET (mm)')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,3)
plt.plot(Time, AET_[:,1])
plt.ylabel('AET (mm)')
plt.xlim([0, Time[n-1]])
plt.ylim([0, max(AET_[:,1])])
plt.subplot(5,1,4)
plt.plot(Time, SMavg_[:,1])
plt.ylabel('S - avg')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,5)
plt.plot(np.arange(1,yrs+1),CumWS[0:yrs,1])
plt.ylabel('Cum WS')
plt.xlim([0, yrs])
plt.xlabel('Time in Years')

pic += 1
plt.figure(pic)
plt.subplot(5,1,1)
plt.bar(Time,P,width = 0.01)
plt.xlim([0, Time[n-1]])
plt.ylabel('Rainfall (mm)')
plt.title('Tree')
plt.subplot(5,1,2)
plt.plot(Time, PETn_[:,2])
plt.ylabel('PET (mm)')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,3)
plt.plot(Time, AET_[:,2])
plt.ylabel('AET (mm)')
plt.xlim([0, Time[n-1]])
plt.ylim([0, max(AET_[:,2])])
plt.subplot(5,1,4)
plt.plot(Time, SMavg_[:,2])
plt.ylabel('S - avg')
plt.xlim([0, Time[n-1]])
plt.subplot(5,1,5)
plt.plot(np.arange(1,yrs+1),CumWS[0:yrs,2])
plt.ylabel('Cum WS')
plt.xlim([0, yrs])
plt.xlabel('Time in Years')

pic += 1
plt.figure(pic)
plt.subplot(3,1,1)
plt.plot(Time, LAI[:,0])
plt.title('Live Leaf Area Index')
plt.ylabel('Grass')
plt.subplot(3,1,2)
plt.plot(Time, LAI[:,1])
plt.ylabel('Shrub')
plt.subplot(3,1,3)
plt.plot(Time, LAI[:,2])
plt.ylabel('Tree')
plt.xlabel('Time in Years')

pic += 1
plt.figure(pic)
plt.plot(CWS_G, Q_G, '-g', label = 'Grass')
plt.hold(True)
plt.plot(CWS_S, Q_S, '-r', label = 'Shrub')
plt.hold(True)
plt.plot(CWS_T, Q_T, '-k', label = 'Tree')
plt.hold(True)
plt.plot(CWS_SS, Q_SS, '--r', label = 'Shrub Seedling')
plt.hold(True)
plt.plot(CWS_TS, Q_TS, '--k', label = 'Tree Seedling')
plt.xlabel('Cumulative Water Stress')
plt.ylabel('Probability of exceedance of Water Stress')
plt.xlim([0,1])
plt.legend(loc=0)

pic += 1
plt.figure(pic)
plt.plot(PETn_[:,2], ETmax[:,2], '.')
plt.hold(True)
plt.plot(range(0,int(PETn_[:,2].max()+1)),range(0,int(PETn_[:,2].max()+1)),'-k')
plt.xlim([0.,PETn_.max()])
plt.ylim([0.,ETmax.max()])
plt.xlabel('Tmax in mm/day')
plt.ylabel('ETmax in mm/day')
plt.title('Tree - Tmax vs ETmax')

pic += 1
plt.figure(pic)
plt.plot(Time, PETn_[:,2], '-k', label = 'Tmax')
plt.hold(True)
plt.plot(Time, ETmax[:,2], '-r', label = 'ETmax')
plt.legend(loc=0)
plt.xlabel('Time in years')
plt.ylabel(' Transpiration mm/day ')
plt.show()
