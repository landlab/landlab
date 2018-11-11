# -*- coding: utf-8 -*-
"""
hex_grid_types.py

Example showing the four types of hex grid.

Created on Sun Nov 16 09:25:04 2014

@author: gtucker
"""

from landlab import HexModelGrid
from numpy import arange
from pylab import figure, show, title

# Case 1: Make and display a hex-shaped grid with horizontal rows of nodes

# Create the grid
hg1 = HexModelGrid(5, 3, 1.0, orientation='horizontal', shape='hex')

# Make some data
d = hg1.add_zeros('node', 'mydata')
d[:] = arange(hg1.number_of_nodes)

# Display the grid
figure(1)
hg1.hexplot(d)
title('hexagon shape, horizontal orientation')
show()


# Case 2: Make and display a hex-shaped grid with vertical columns of nodes

# Create the grid
hg2 = HexModelGrid(3, 5, 1.0, orientation='vertical', shape='hex')

# Make some data
d = hg2.add_zeros('node', 'mydata')
d[:] = arange(hg2.number_of_nodes)

# Display the grid
figure(2)
hg2.hexplot(d)
title('hexagon shape, vertical orientation')
show()


# Case 3: Make and display a rectangular-shaped grid with horizontal rows
# of nodes

# Create the grid
hg3 = HexModelGrid(5, 5, 1.0, orientation='horizontal', shape='rect')

# Make some data
d = hg3.add_zeros('node', 'mydata')
d[:] = arange(hg3.number_of_nodes)

# Display the grid
figure(3)
hg3.hexplot(d)
title('rectangular shape, horizontal orientation')
show()


# Case 4: Make and display a rectangular-shaped grid with vertical rows of
# nodes

# Create the grid
hg4 = HexModelGrid(5, 5, 1.0, orientation='vertical', shape='rect')

# Make some data
d = hg4.add_zeros('node', 'mydata')
d[:] = arange(hg4.number_of_nodes)

# Display the grid
figure(4)
hg4.hexplot(d)
title('rectangular shape, vertical orientation')
show()
