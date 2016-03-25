""" generate_detachment_ltd_transport.py

 This component simulates detachment-limited sediment
 transport conditions based Whipple and Tucker, (2002).

Written by Jordan Adams

E = K*f(qs)*(A**m)*(S**n)

"""

from landlab import Component, ModelParameterDictionary, CLOSED_BOUNDARY
import pylab
import numpy as np
from matplotlib import pyplot as plt
from math import atan, degrees
import os


_DEFAULT_INPUT_FILE = os.path.join(os.path.dirname(__file__),
                                  'detach_ltd_input.txt')

class DetachmentLtdErosion(Component):
    '''  Landlab component that simulates detachment-limited river erosion.

This component calculates changes in elevation in response to vertical incision. For now
this does not assume an erosion threshold (8/29/14).


    '''
    _name = 'DetachmentLtdErosion'

    _input_var_names = set([
        'elevation',
        'slope_at_nodes',
        'water_discharge_at_nodes',
    ])
    _output_var_names = set([
        'elevation',
    ])

    _var_units = {
        'elevation': 'm',
        'water_discharge_at_nodes': 'm3/s',
        'slope_at_nodes': 'm/m',
    }

    def __init__(self, grid, **kwds):
        self.m = kwds.pop('m', 0.5)
        self.n = kwds.pop('n', 1.0)
        self.K = kwds.pop('K', 0.000005)
        self.f_qs = kwds.pop('f_qs', 1.0)
        self.U = kwds.pop('U', 0.0)
        self.rg = grid
        self.dzdt = self.rg.zeros(centering='node')
        self.summed_dzdt = self.rg.zeros(centering='node')



    def change_elev(self, z, S, q):
        self.z = z
        self.slope = S
        self.Q = q
        self.Q_to_m = self.Q**self.m
        self.S_to_n = self.slope**self.n
        self.interior_nodes = (
            np.where(self.rg.status_at_node != CLOSED_BOUNDARY)[0])

        self.dzdt = self.U - (self.K*self.f_qs*self.Q_to_m*self.S_to_n)


        self.z += self.dzdt
        #self.summed_dzdt += self.dzdt
        return self.z, self.dzdt

    def plot_elev_changes(self):
        plt.figure('Elevation Changes')

        # This call selects the color map and interval which it will be blocked by
        cmap=plt.get_cmap('RdYlGn', 10)


        # This converts the new slope vector (in degrees) to a raster
        fx = self.rg.node_vector_to_raster(self.summed_dzdt)

        # All slopes beneath the minimum are shown in white.
        cmap.set_over('white')

        # Creating the grid image to show.
        im2 = pylab.imshow(fx, cmap=cmap, extent=[0, self.rg.number_of_node_columns *self.rg.dx,0, self.rg.number_of_node_rows * self.rg.dx])

        # Setting minimum and maximum values for the colorbar
        #pylab.clim(vmax=-0.0000000001)

        # Creating the colorbar instance
        cb = pylab.colorbar(im2)

        # Set titles for the figure and the colorbar.
        cb.set_label('dzdt', fontsize=12)
        pylab.title('dzdt')
        plt.show()

