######### soil_moisture_dynamics_driver.py ###########
##
## Soil Moisture Dynamics Model:
## 
##    This Model simulates 'soil moisture dynamics' using following 
## eco-hydrologic components: 
##   
##    1. SM_spatial_DEM_PET.py - Evaluates soil moisture decay every storm
##    2. ET_PriestlyTaylor_Generated.py - Estimates Priestly Taylor potential
##                                                evapotranspiration
##    3. generate_uniform_precip.py - Generate storms
##    4. radiation.py - Estimate influence of aspect on potential -
##                                               evapotranspiration
##    
##    This Model also demonstrates use of cool landlab utilities that import
## a DEM (read_esri_ascii.py) and plot variables on raster grid (imshow.py). 
##
## Sai Nudurupati and E.I. 11/17/2013
#############################################

# Import all required components and tools

from numpy import *
from landlab.components.soil_moisture.SM_spatial_DEM_PET import *
from landlab.components.PET.ET_PriestlyTaylor_Generated import *
from landlab.components.radiation.radiation import *
from landlab.io import read_esri_ascii
import time         # We will need this to check performance of the code

# Point to the input DEM
_DEFAULT_INPUT_FILE_2 = os.path.join(os.path.dirname(__file__),
                                 'saihugoRAW.asc')

""" Create an Soil Moisture Object """
SM = SoilMoisture()

""" Initialize Soil Moisture """
SM.initialize()

""" Set number of storms for the simulation """
SM.set_strms( 5 )
SM._length = SM._iterate_storm

""" Create an instance of Storm Class """
PD = PrecipitationDistribution()

""" Initialize Storm - Create first storm """
PD.initialize()

""" Create an instance of Potential Evapotranspiration (pet) Class """
PET = ET()

""" Initialize potential evaporation class """
PET.initialize()
pet = 0     # local variable initialization

""" Importing Grid and Elevations from DEM """
(RMG,elev) = read_esri_ascii(_DEFAULT_INPUT_FILE_2)

""" Using a dummy variable to capture shape of a multi-dimensional 
    array spread over (all) cells on the raster grid """
dummy = RMG.cell_vector_to_raster(RMG.zeros(centering = 'cell'))

""" Initializing OPTIONAL arrays to store values of variables and plotting """

# following arrays are constant in space but variable in time
""" Store precipitation depth """
P_ = zeros( SM._length, dtype = float )
""" Store effective precipitation depth """ 
Peff_ = zeros( SM._length, dtype = float )
""" Store local time after every storm """
TIME = zeros( SM._length, dtype = float )

# following arrays store data that are variable both in space & time
""" Store soil moisture """
S_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )
""" Store actual evapotranspiration """
ETA_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )
""" Store drainage """
D_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )
""" Store runoff """
Ro_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )
""" Store water stress """
WS_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )
""" Store radiation factor that records influence of aspect """
RADF_ = zeros( (SM._length,dummy.shape[0],dummy.shape[1]), dtype = float )

""" Initializing 'Initial Soil Moisture' """
SO = 0.5 * ones((dummy.shape))   

""" Create an instance of Radiation Class """
RAD = Radiation()

""" Initialize Radiation class """
RAD.initialize()

# Transfer Elevation data
RAD._Z = RMG.empty()
RAD._Z = elev  

#######
""" Evaluated time taken for the simulation """
Time_start = time.clock()    # Real time for code performance measure
print 'Clock starts ticking now!'

""" Loop through the storms """
for i in range( 0, SM._iterate_storm ):
    if i != 0:
       """ Create a new storm """ 
       PD.update()
    """ Obtain storm data """
    P = PD.storm_depth
    Tb = PD.interstorm_duration
    Tr = PD.storm_duration
    """ Assume fully vegetated field for the time being """
    Fveg = 1.0 #self._soil_Fveg
    
    """ Obtain 'pet' for the storm at half the time between storms.
        Assuming 'pet' will remain constant between the storms"""
    PET.update((TIME[i-1]+(1.0/2.0)*(Tb+Tr)/(24*365.4)))
    pet = max(PET._ETp, 6) # Adjustment of PET to be atleast 6 mm/day
    
    """ Evaluate influence of aspect at the same time when 'pet' is
        obtained.Assuming influence of aspect remains the same 
        between storms """
    RAD.set_current_time( (TIME[i-1]+(1.0/2.0)*(Tb+Tr)/(24*365.4)) )
    RAD.update(RMG)
    RADF_[i] = RMG.cell_vector_to_raster(RAD._Si)

    """ Distribute 'pet' spatially over grid considering aspect influence """
    SM.set_PET((pet*RADF_[i]))    
    

    """ Update Soil Moisture """
    SM.update( P, Tb, Tr, Fveg, SO)
    
    """ Record outputs """
    P_[i] = P
    Peff_[i] = SM.get_Peff()  
    S_[i] = SM.get_S()
    ETA_[i] = SM.get_ET()     
    Ro_[i] = SM.get_Ro()
    D_[i] = SM._D             
    TIME[i] = SM.get_time()      
    WS_[i] = SM.get_WaterStress()
    SO = S_[i]

#######
Time_end = time.clock()                               # Code Performance time
Time_elapsed = Time_end - Time_start                  # Code Performance time
print 'Time Elapsed = ', Time_elapsed, 'secs'        # Code Performance time

""" Plot output - urrently plots soil moisture state at the end of all storms"""

SM.plott( P_, S_[i], ETA_[i], D_[i], Ro_[i], TIME, RMG )