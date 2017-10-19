from netCDF4 import Dataset
import numpy

def calculate_subsidence (dictionary) :
    
    if len (x) < number_of_node_rows:
        while len (x) < number_of_node_rows:
            x = x.append(x[len(x)- 1] + dx)
            y = y.append(y[len(y) - 1])
            return x, y
        print ('file is too short, last value was repeated')
        
    f = interpolate.interp1d(x, y, kind= 'cubic')
    distance = arange(0, dx, (number_of_node_rows -1) * dx)
    subsidence_array = f(xnew)
    
    plt.plot(x, y, 'o', distance, subsidence_array, '-')
    plt.show()
    # plots the interpolated arrays
    return distance, subsidence_array