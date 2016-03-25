#! /usr/env/python
"""

Python implementation of drainag area calculation.

Note that this does not properly handle pits.  In the case of pit, two cells
can point to each other.  If this happens, the code would enter an infinte loop.
Have a temporary fix for this, but this needs improving.

Last updated NG 6/2013

"""
from __future__ import print_function

#from landlab.model_grid import RasterModelGrid
from numpy import *


class CalcDrainageArea(object):
    """
    This class finds the drainage area at all active nodes on a mesh.
    Note that the drainage area is initially set to zero.
    Drainage area is not calculated at boundary points, so these points
    should have a drainage area value of zero.
    """

    def __init__(self, num_cells):
        """
        This sets the num_cells parameter.
        This class assumes that the number of cells does not change
        after a class item has been instantiated.
        """

        self.num_cells = num_cells
        self.initialize()

        # print 'CalcDrainageArea.__init__'

    def initialize(self):
        """
        This sets up the drainage area vector.
        It is initialized to zero for all nodes.
        Is this method really needed?  Can we just put this in __init__?
        """

        self.drainarea = zeros(self.num_cells)

        # print 'CalcDrainageArea.__initialize__'

    def calc_DA(self, mg, fd):
        """
        This calculates the drainage area at each cell by moving downstream
        from each cell and adding the area of that cell to the drainage area of
        the downstream cells.  This algorithm includes the local area of Cell A
        in Cell A's drainage area.  In otherwords, the cells that only drain
        themselves have a drainage area that is the cell area.
        This method assumes that the flow direction vector was updated before
        the method was called.

        Method inputs: the model grid and the flow direction vector
        Method returns: the drainage area vector
        """

        # loop through each cell
        for i in range(0, self.num_cells):
            # only update drainage area in interior cells
            if mg.is_interior(i):
                # add local cell area to its drainage area
                self.drainarea[i] = self.drainarea[i] + mg.cellarea
                prev_cell = i
                next_cell = fd[i]
                # move down the flow direction, adding the cell area to the
                # drainage area of all of the downstream cells
                # stop if two cells have flow directions pointing to each other
                # this can happen if there is a pit, which would trigger an
                # infinite loop.  This s not a good fix, but a temporary fix.
                while (mg.is_interior(next_cell) and fd[next_cell] != prev_cell):
                    self.drainarea[next_cell] = self.drainarea[
                        next_cell] + mg.cellarea
                    prev_cell = next_cell
                    next_cell = fd[next_cell]

        return self.drainarea
