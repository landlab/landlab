#! /usr/env/python
"""Set boundary conditions for watersheds from a DEM.

Classes for Landlab that deal with setting boundary conditions if running
on a DEM.  This takes advantage of methods in the grid, but combines them
to make it easier for a user.
"""
from landlab import FIXED_VALUE_BOUNDARY
import numpy as np
import pylab


class WatershedBoundaryConditions(object):

    """Set the boundary conditions of a watershed extracted from a DEM."""

    def set_bc_node_coords(self, grid, node_data, nodata_value, outlet_row,
                           outlet_column):
        """Set the boundary conditions for a watershed.

        Assumes that outlet is already known.

        This must be passed the grid, node_data and nodata_value,
        and the values of the outlet_row and outlet_column.
        """
        # set no data nodes to inactive boundaries
        grid.set_nodata_nodes_to_inactive(node_data, nodata_value)

        # find the id of the outlet node
        outlet_node = grid.grid_coords_to_node_id(outlet_row, outlet_column)
        # set the boundary condition (fixed value) at the outlet_node
        grid.status_at_node[outlet_node] = FIXED_VALUE_BOUNDARY

    def set_bc_node_id(self, grid, node_data, nodata_value, outlet_node):
        """
        Sets the boundary conditions for a watershed.
        Assumes that outlet is already known.

        This must be passed the grid, node_data and nodata_value,
        and the id of the outlet node.
        """

        #set no data nodes to inactive boundaries
        grid.set_nodata_nodes_to_inactive(node_data, nodata_value)

        #set the boundary condition (fixed value) at the outlet_node
        grid.status_at_node[outlet_node] = FIXED_VALUE_BOUNDARY

    def set_bc_find_outlet(self, grid, node_data, nodata_value):
        """
        Finds the node adjacent to a boundary node with the smallest value.
        This node is set as the outlet.

        This must be passed the grid, node_data and nodata_value.
        """
        #for this to be a watershed, need to make sure that there is a ring
        #of no data values around the outside of the watershed, barring the
        #outlet location.  So enforce that all outer nodes
        #are inactive boundaries now, then set the outlet location later.
        #By enforcing the ring of closed values first, then fixing the outlet
        #later, it should be OK if the outlet is on the outer ring.
        grid.set_inactive_boundaries(True, True, True, True)

        #set no data nodes to inactive boundaries
        #this may be redundant, but must do in case there are no data
        #values that are not on the outer boundary
        grid.set_nodata_nodes_to_inactive(node_data, nodata_value)

        #This method works well if the watershed topography is already
        #established.  If it's not, then this is an ineffiient method, but
        #seems likely that one would only call this if the watershed
        #topography was already established.

        #ng this could maybe be generalized (?) so you don't need this if/else
        #first case is if nodata_value is minimum value
        if (min(node_data) == nodata_value):
            #min value is nodata_value, so need to find values
            #that are not no_data

            #locs is a list that contains locations where
            #node data is greater than the nodata value
            locs=list(np.where(node_data>nodata_value)[0])

            #now find minimum of the data values
            min_val=np.min(node_data[locs])

            #now find where minimum values are
            min_locs=list(np.where(node_data==min_val)[0])

            #check all the locations with the minimum value to see if one
            #is adjacent to a boundary location.  If so, that will be the
            #watershed outlet.  If none of these points qualify, then
            #increase the minimum value and check again.  Keep checking
            #until a point next to the boundary is found.
            #
            #NG I think the only way this would become an infinite loop
            #is if there are no interior nodes.
            not_found=True
            while not_found:
                #now check these locations to see if any are next to
                #a boundary node
                local_not_found=True
                i=0
                while (i<len(min_locs) and local_not_found):
                    if grid.has_boundary_neighbor(min_locs[i]):
                        local_not_found = False
                        #outlet_loc contains the index of the outlet location
                        #in the node_data array
                        outlet_loc=min_locs[i]
                    else:
                        i += 1

                #checked all of the min vals, (so done with inner while)
                #and none of the min values were outlet candidates
                if local_not_found:
                    #need to find the next largest minimum value
                    #first find the locations of all values greater
                    #than the old minimum
                    #not done with outer while
                    locs=list(np.where(node_data>min_val)[0])
                    #now find new minimum of these values
                    min_val=np.min(node_data[locs])
                    min_locs=list(np.where(node_data==min_val)[0])
                else:
                    #if locally found, it is also globally found
                    #so done with outer while
                    not_found = False
        else:
            #entering this else b/c nodata_value is not the minimum value
            #can I generalize this?
            #will the no data value ever be anything other than the minimum
            #value?
            #this code is exactly the same as above, just have to check that
            #every min value is not the no data value.
            min_val=np.min(node_data)
            if (min_val == nodata_value):
                #found place where no data value is
                #find locations that have values greater than
                #the no data value, then find the min val of these locations
                locs=list(np.where(node_data>nodata_value)[0])
                min_val = np.min(node_data[locs])

            #now find where minimum values are
            min_locs=list(np.where(node_data==min_val)[0])
            #check all the locations with the minimum value to see if one
            #is adjacent to a boundary location.  If so, that will be the
            #watershed outlet.  If none of these points qualify, then
            #increase the minimum value and check again.  Keep checking
            #until a point next to the boundary is found.
            #
            #NG I think the only way this would become an infinite loop
            #is if there are no interior nodes.
            not_found=True
            while not_found:
                #now check these locations to see if any are next to
                #a boundary node

                local_not_found=True
                i=0
                while (i<len(min_locs) and local_not_found):
                    if grid.has_boundary_neighbor(min_locs[i]):
                        local_not_found = False
                        #outlet_loc contains the index of the outlet location
                        #in the node_data array
                        outlet_loc=min_locs[i]
                    else:
                        i += 1

                #checked all of the min vals, (so done with inner while)
                #and none of the min values were outlet candidates
                if local_not_found:
                    #need to find the next largest minimum value
                    #first find the locations of all values greater
                    #than the old minimum
                    #not done with outer while
                    locs=list(np.where(node_data>min_val)[0])
                    #now find new minimum of these values
                    min_val=np.min(node_data[locs])
                    if (min_val == nodata_value):
                        #found place where no data value is
                        #find locations that have values greater than
                        #the no data value, then find the min val of these locations
                        locs=list(np.where(node_data>nodata_value)[0])
                        min_val = np.min(node_data[locs])
                    #now find where minimum values are
                    min_locs=list(np.where(node_data==min_val)[0])
                else:
                    #if locally found, it is also globally found
                    #so done with outer while
                    not_found = False

        #set outlet boundary condition
        grid.status_at_node[outlet_loc] = FIXED_VALUE_BOUNDARY
        #x=grid.node_x[outlet_loc]
        #y=node_y[outlet_loc]
        #print "outlet_loc ", outlet_loc," x ",x," y ",y
        return outlet_loc
