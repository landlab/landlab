from matplotlib.pyplot import figure, show
from numpy import arange, sin, pi
from scipy import interpolate
import matplotlib.pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
from landlab import RasterModelGrid, Component


class SeaLevelTimeSeries(Component):

    _name = 'Sea Level Changer'

    _time_units = 'y'

    _input_var_names = ()

    _output_var_names = (
        'sea_level__elevation',
    )

    _var_units = {
        'sea_level__elevation': 'm',
    }

    _var_mapping = {
        'sea_level__elevation': 'grid',
    }

    _var_doc = {
        'sea_level__elevation': 'Sea level elevation',
    }

    def __init__(self, grid, filepath, kind='linear', start=0., **kwds):
        """Generate sea level values.

        Parameters
        ----------
        grid: ModelGrid
            A landlab grid.
        filepath: str
            Name of csv-formatted sea-level file.
        kind: str, optional
            Kind of interpolation as a string (one of 'linear',
            'nearest', 'zero', 'slinear', 'quadratic', 'cubic').
            Default is 'linear'.
        """
        super(SeaLevelTimeSeries, self).__init__(grid, **kwds)

        data = np.loadtxt(filepath, delimiter=',')
        self._sea_level = interp1d(data[:, 0], data[:, 1], kind=kind,
                                   copy=True, assume_sorted=True,
                                   bounds_error=True)

        self._time = start

    @property
    def time(self):
        return self._time

    def run_one_step(self, dt):
        self._time += dt
        self.grid.at_grid['sea_level__elevation'] = self._sea_level(self.time)


class SinusoidalSeaLevel(SeaLevelTimeSeries):

    def __init__(self, grid, wave_length=1., amplitude=1., phase=0.,
                 start=0., **kwds):
        """Generate sea level values.

        Parameters
        ----------
        grid: ModelGrid
            A landlab grid.
        """
        wave_length /= 2. * np.pi
        super(SeaLevelTimeSeries, self).__init__(grid, **kwds)

        self._sea_level = lambda time: np.sin((time - phase) / wave_length) * amplitude

        self._time = start


def sea_level_type(dictionary):
    from landlab.components.submarine_diffusion.sea_level import sea_level_file
    sl_type = dictionary['sl_type']
    if sl_type == 'sinusoid':
        return sea_level_function(dictionary)
    else:
        sl_file_name = dictionary['sl_file_name']
        return sea_level_file(sl_file_name, dictionary)
       

def sea_level_function(dictionary):
    """
    t is an array of x values (ex. arange(0,1000, pi/4))
    a is the amplitude of the sin function
    p is th phase shift
    t is the title of the graph (String)
    xt is the title of the x axis (string)
    xy is the title of the y axis (STRING)
    Function starts at 0 or P.
    Fs is the period of the function. (10,000)
    """
    p = dictionary['sea_level_phase'] 
    a = dictionary['sea_level_amplitude'] 
    Fs = dictionary ['sea_level_period']
    start_time = dictionary['start_time']
    run_duration = dictionary['run_duration']
    dt = dictionary['dt']
    t = arange(start_time, start_time + run_duration, dt)
    sl_array = a * sin((2*pi * (p +t))/ Fs)
    """
    fig = plt.figure()
    plt.plot(t, sl_array)
    fig.suptitle(ti)
    plt.xlabel(xt)
    plt.ylabel(xy)
    plt.show()
    """
    return t, sl_array
    
def sea_level_file(filename, dictionary):
    """
    reading in the file above
    x is an array of x values (ex. x = arange(0, 10))
    y is an array of y values (ex. y = np.exp(x/2.0))
    start time (default should be 0)
    dt
    run duration
    
    Note: the array of x values can be pretermined to a set of values. Goes backwards so the start will be at -12500 years
    There will be a sea level array that stores these values
    """
    start_time = dictionary['start_time']
    run_duration = dictionary['run_duration']
    dt = dictionary['dt']
    
    xes = []
    ys = []
    with open(filename) as f:
         for line in f:
             x, y = line.split()
             xes.append(x)
             ys.append(y)
    x = []
    for item in xes:
        x.append(float(item))
    y = []
    for item in ys:
        y.append(float(item))    
    
    f = interpolate.interp1d(x, y, kind= 'cubic')
    times = arange(start_time, start_time + run_duration, dt)
    sl_array = f(times)
    plt.plot(x, y, 'o', times, sl_array, '-')
    plt.show()
    return times, sl_array
