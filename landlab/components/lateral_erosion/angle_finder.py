#plotter for runs with saved output

import numpy as np
from pylab import * 
from landlab import RasterModelGrid
import matplotlib.pyplot as plt
import time
import math
#reload(flow_accumulation)
#reload(raster)

    
def angle_finder(grid, dn, cn, rn):
	xcoord=grid.node_axis_coordinates(axis=0)
	ycoord=grid.node_axis_coordinates(axis=1)
	#print "xcoord", xcoord
	
	sl1=(ycoord[cn]-ycoord[dn])/(xcoord[cn]-xcoord[dn])
	sl2=(ycoord[rn]-ycoord[cn])/(xcoord[rn]-xcoord[cn])
	#print "slope1", sl1
	#print "slope2", sl2
	angle1=math.degrees(math.atan(sl1))
	angle2=math.degrees(math.atan(sl2))
	#print "angle1", angle1
	#print "angle2", angle2
	angle_diff=angle2-angle1
	#if (angle_diff==45. or angle_diff==135.):
	#	if angle1==0.0:
	#		lat_node=
	
	#print "anglediff", angle_diff
	#print delta
	#ang_lat=(abs(angle_diff)-315.)*sign(angle_diff)
	#print "angle_diff=abs(angle2-angle1)", abs(angle_diff)-135.
	#print sign(angle_diff)
	#print "ang lat", ang_lat
	#lat_x=math.degrees(mat.sine(ang_lat))
	#lat_y=math.degrees(mat.sine(ang_lat))
	angle_diff=abs(angle2-angle1)
	#print "angle diff", angle_diff
	return angle_diff
if __name__ == '__main__':
    main()