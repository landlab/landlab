#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Create a RockBlock object with different properties. 

"""

import numpy as np
from landlab.layers import EventLayers


class RockBlock(object):
    """Create a RockBlock object.
    
    A RockBlock is a three dimentional representation of material operated on 
    by landlab components. Material can be removed through erosion or added to
    through deposition. 
    
    It is constructed by specifying a series of thicknesses and a series of 
    rock type IDs. Thicknesses and IDs are both specified in order of  closest 
    to the surface to furthest from the surface. Thicknesses can either be a 
    single value (cooresponding to a layer of uniform thickness) or a number-of
    -nodes length array (cooresponding to a non-uniform layer).
    
    Additionally, an attribute dictionary specifies the attributes of each 
    rock type. This dictionary is expected to have the form of:
     
    .. code-block:: python

        attrs = {'K_sp': {1: 0.001,
                          2: 0.0001},
                 'D': {1: 0.01, 
                       2: 0.001}}
    
    Where ``'K_sp'`` and ``'D'`` are attributes to track, and ``1`` and ``2`` 
    are rock type IDs. The rock type IDs can be any type that is valid as a 
    python dictionary key. 
    
    Methods
    -------
    add_rock_type
    add_rock_attribute
    update_rock_attribute
    add_layer
    get_surface_values
    run_one_step
        
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.layers import RockBlock
    >>> mg = RasterModelGrid(3, 3)
    >>> thickness = [1, 2, 4, 1]
    >>> ids = [1, 2, 1, 2]
    >>> attrs = {'K_sp': {1: 0.001,
    ...                   2: 0.0001},
    ...          'D': {1: 0.01, 
    ...                2: 0.001}}
    >>> rb = RockBlock(mg, thickness, ids, attrs)
    >>> rb.get_surface_values('K_sp')
    array(0.001, 0.001)
    """

    _name = 'RockBlock'

    _cite_as = """ """

    def __init__(self, grid,  thicknesses, ids, attrs):
        """Create a RockBlock.

        More here
        
        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        thicknesses, top first first... 
        ids
        attr
        
        Examples
        --------
        >>> 
        
        """
        # save reference to the grid and the last time steps's elevation. 
        self._grid = grid
        self.last_elevation = self._grid['node']['topographic__elevation'][:].copy()

        # save inital information about thicknesses, layers, attributes, and ids. 
        self._init_thicknesses = np.asarray(thicknesses)
        self._ids = np.asarray(ids)
        self._attrs = attrs
        self._number_of_init_layers = self._init_thicknesses.shape[0]
        self.attributes = list(attrs.keys()) 
        self.ids = np.unique(self._ids)
        
        # assert that thicknesses and ids are correct and consistent shapes
        if self._init_thicknesses.shape != self._ids.shape:
            raise ValueError(('Thicknesses and IDs provided to RockBlock are ',
                              'inconsistent with each other.'))
        if self._init_thicknesses.shape[1] != self._grid.number_of_nodes:
            raise ValueError(('Thicknesses provided to RockBlock are ',
                              'inconsistent with the ModelGrid.'))
            
        # assert that attrs are pointing to fields (or create them)
        for at in self.attributes:
            if at not in grid.at_node:
                self._grid.add_empty('node', at) 
                
        # add a field for the rock type id
        self._grid.add_empty('node', 'rock_type__id')

        # verify that and that all IDs have attributes. 
        self._check_attribute_dictionary()
        
        # create a EventLayers instance
        self._layers = EventLayers(grid.number_of_nodes)
        
        # From bottom to top, add layers to the RockBlock with attributes.
        for i in range(self._number_of_init_layers-1, -1, -1):
            self.add_layer(self._init_thicknesses[i, :], self._ids[i, :])
    
        # update values for attributes at the topographic surface. 
        self._update_surface_values()
        
        # check that rock exists at the surface everywhere. 
        self._check_thickness()
    
    def _check_attribute_dictionary(self):
        """
        Check compatibility of RockBlock and attribute dictionary. 
        
        Parameters
        ----------
        at_dict : dict
            Attribute dictionary
            
        """
        for at in self.attributes:
            for i in self.ids:
                if i not in self._attrs[at]:
                    msg = ('A rock type with ID value ' + str(i) + 'was '
                           'specified in instantiation of RockBlock. No value '
                           'for this ID was provided in attribute ' + at + '.')
                    raise ValueError(msg)
                    
    def _check_thickness(self):
        """Check thickness."""
        # verify that rockblock has some thickness in each stack, otherwise 
        # raise a Runtime error
        if np.any(self._layers.thickness == 0):
            msg = ('RockBlock has a thickness of zero at at least one node.')
            raise RuntimeError(msg)
                        
    def _update_surface_values(self):
        """Update surface values"""
        # Update surface values for each attribute.
        self._grid['node']['rock_type__id'][:] = self.get_surface_values('rock_type__id')
        for at in self.attributes:            
            self._grid['node'][at][:] = self.get_surface_values(at)
        
    def add_layer(self, thickness, ids=None):
        """
        Add a new layer to RockBlock.
        
        Parameters
        ----------
        thickness : float or array
            Float or number-of-nodes array of thicknesses. Positive values add
            to RockBlock while negative values erode RockBlock. 
        ids : single value or itterable, optional if only erosion occurs
            Rock type ID for new deposits. Can be single value or an number-
            of-nodes array. 
            
        Examples
        --------
        
        """
        new_layer_attributes = {'rock_type__id': ids}
        for at in self.attributes:
            try:
                layer_value = list(map(self._attrs[at].get, ids))
            except TypeError: # if not iterable
                layer_value = self._attrs[at].get(ids)
                
            new_layer_attributes[at] = layer_value
        self._layers.add(thickness, **new_layer_attributes)
    
    def get_surface_values(self, at):
        """
        Get RockBlock attribute values at the topographic surface. 
        
        Parameters
        ----------
        at : str
            Name of attribute
        
        Examples
        --------
        
        """
        return self._layers.surface_values(at)
    
    def add_rock_attribute(self, attrs):
        """
        Add new attribute to RockBlock
        
        Parameters
        ----------
        attrs : dict
            Rock attribute dictionary for the new attribute(s).
        
        """
    
    def add_rock_type(self, attrs):
        """
        Add rock type to RockBlock.
        
        Parameters
        ----------
        attrs : dict
            Rock attribute dictionary for the new rock type(s).
            
        Examples
        --------
        
        """
        # Check that the new rock type has all existing attributes
        
        # do this with sets. 
        for at in self.attributes:
            if at not in attrs:
                msg = 'The new rock type is missing attribute ' + str(at) + '.'
                raise ValueError(msg)
        # And no new attributes
        for at in attrs:
            if at not in self.attributes:
                msg = ('The new rock type has an attribute (' + str(at) + ') '
                       'that no other rock type has. This is not permitted.')
                raise ValueError(msg)
        
        new_ids = []
        for at in attrs:
            att_dict = attrs[at]
            rids = att_dict.keys()
            for rid in rids:
                if rid in self._ids:
                    msg = ('Rock type ID ' + str(rid) + ' for attribute ' 
                           '' + str(at) + ' has already been added. This is '
                           'not allowed')
                    raise ValueError(msg)
                else:
                    new_ids.append(rid)
                    self._attrs[at][rid] = att_dict[rid]
        self._ids.extend(list(set(new_ids)))
    
    def update_rock_attribute(self, at, rock_id, value):
        """
        Update rock type attribute.
        
        Paramters
        ---------
        
        
        Examples
        
        """
        if at not in self.attributes:
            msg = ('RockBlock cannot update the value of ' + str(at) + 'as '
                   'this attribute does not exist.')
            raise ValueError(msg)
        
        if rock_id not in self._attrs[at]:
            msg = ('RockBlock cannot update the value of rock type '
                   '' + str(rock_id) + 'for attribute ' + str(at) + ' as this '
                   'rock type is not yet defined.')
            raise ValueError(msg)
        
        self._attrs[at][rock_id] = value
        
    def run_one_step(self, dz_advection=0, layer_id=None):
        """
        Update RockBlock.
        
        Parameters
        ----------
        dz_advenction : float or array, optional
            Change in rock elevation due to advection by some external process.
        layer_id : 
            
        Examples
        --------
        
        """
        # calculate amount of erosion
        elevation_change = (self._grid['node']['topographic__elevation'] - 
                            (self.last_elevation + dz_advection))
        
        
        if layer_id not in self._ids:
            if np.any(elevation_change>0):
                msg = 'You are depositing without a valid rock type'
                raise RuntimeError(msg)
            
        self.add_layer(elevation_change)
        
        # update surface values
        self._update_surface_values()
        
        # check that rock exists at the surface everywhere. 
        self._check_thickness()
        
        # update the last elevation. 
        self.last_elevation = self._grid['node']['topographic__elevation'][:].copy()
