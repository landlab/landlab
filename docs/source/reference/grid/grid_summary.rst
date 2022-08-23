Grid Cheatsheet
===============

Getting Information about a Grid
--------------------------------

The following attributes, properties, and methods provide data about the grid,
its geometry, and the connectivity among the various elements. Each grid
element has an ID number, which is also its position in an array that
contains information about that type of element. For example, the *x*
coordinate of node 5 would be found at `grid.node_x[5]`.

The naming of grid-element arrays is *attribute*`_at_`*element*, where
*attribute* is the name of the data in question, and *element* is the element
to which the attribute applies. For example, the property `node_at_cell`
contains the ID of the node associated with each cell. For example,
`node_at_cell[3]` contains the *node ID* of the node associated with cell 3.
The *attribute* is singular if there is only one value per element; for
example, there is only one node associated with each cell. It is plural when
there are multiple values per element; for example, the `faces_at_cell` array
contains multiple faces for each cell. Exceptions to these general rules are
functions that return indices of a subset of all elements of a particular type.
For example, you can obtain an array with IDs of only the core nodes using
`core_nodes`, while `active_links` provides an array of IDs of active links
(only). Finally, attributes that represent a measurement of something, such as
the length of a link or the surface area of a cell, are described using `_of_`,
as in the example `area_of_cell`.

    
.. contents::
  :local:


Grid Elements
~~~~~~~~~~~~~

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, cats in grids |dictsort %}
  
  .. tab:: {{ grid }}
    
    {% for cat, label in [('info-node', 'Nodes'), ('info-link', 'Links'), ('info-patch', 'Patches')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in cats[cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %} 
  
Dual-Grid Elements
~~~~~~~~~~~~~~~~~~

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, cats in grids |dictsort %}
  
  .. tab:: {{ grid }}
    
    {% for cat, label in [('info-corner', 'Corners'), ('info-face', 'Faces'), ('info-cell', 'Cells')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in cats[cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %}  


Grid
~~~~

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, cats in grids |dictsort %}
  
  .. tab:: {{ grid }}
    
    {% for cat, label in [('boundary-condition', 'Boundary Conditions'), ('subset', 'Subsetting')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in cats[cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %}  


Fields
~~~~~~

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, cats in grids |dictsort %}
  
  .. tab:: {{ grid }}
          
    {% for cat, label in [('field-add', 'New'), ('field-io', 'Access'), ('map', 'Mappers'), ('gradient', 'Gradients'), ('surface', 'Analysis')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in cats[cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %}  


Deprecated and Uncategorized
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. jinja:: llcats
  
  .. currentmodule:: landlab
    
  {% for grid, cats in grids |dictsort %}
  
  .. tab:: {{ grid }}
    
    {% for cat, label in [('uncategorized', 'Uncategorized'), ('deprecated', 'Deprecated')] %}
    
      .. tab:: {{label}}
      
        .. autosummary::
          :nosignatures:
        
          {% for func in cats[cat] %}
            ~{{func}}      
          {% endfor %}
    {% endfor %}
  {% endfor %}
  

Data Fields
~~~~~~~~~~~


:class:`~.ModelGrid` inherits from the :class:`~.GraphFields` class. This
provides `~.ModelGrid`, and its subclasses, with the ability to, optionally,
store data values that are associated with the different types grid elements
(nodes, cells, etc.). In particular, as part of ``ModelGrid.__init__()``,
data field *groups* are added to the `ModelGrid` that provide containers to
put data fields into. There is one group for each of the eight grid elements
(node, cell, link, face, core_node, core_cell, active_link, and active_face).

To access these groups, use the same methods as accessing groups with
`~.GraphFields`. ``ModelGrid.__init__()`` adds the following attributes to
itself that provide access to the values groups:


.. jinja:: llcats

  .. currentmodule:: landlab
    
  {% for grid, cats in grids.items() %}
  
  .. tab:: {{ grid }}
    
      .. tab:: Access
        
        Each of these attributes returns a ``dict``-like object whose keys are value
        names as strings and values are numpy arrays that gives quantities at
        grid elements.
          
        .. autosummary::
          :nosignatures:
          
          ~landlab.{{grid}}.at_node
          ~landlab.{{grid}}.at_cell
          ~landlab.{{grid}}.at_link
          ~landlab.{{grid}}.at_face
          ~landlab.{{grid}}.at_patch
          ~landlab.{{grid}}.at_corner
            
      .. tab:: New
      
        :class:`~.ModelGrid` inherits several useful methods for creating new data
        fields and adding new data fields to a ModelGrid instance. Methods to add or
        create a new data array follow the ``numpy`` syntax for creating arrays. The
        folowing methods create and, optionally, initialize new arrays. These arrays
        are of the correct size but a new field will not be added to the field:
            
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["field-add"] %}
            ~{{func}}      
          {% endfor %}
          
      .. tab:: Element Mapping
        
        These methods allow mapping of values defined on one grid element type onto a
        second, e.g., mapping upwind node values onto links, or mean link values onto
        nodes.
        
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["map"] %}
            ~{{func}}      
          {% endfor %}
      
      .. tab:: Element Mapping
            
        Landlab is designed to easily calculate gradients in quantities across the
        grid, and to construct fluxes and flux divergences from them. Because these
        calculations tend to be a little more involved than property lookups, the
        methods tend to start with `calc_`.
        
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["map"] %}
            ~{{func}}      
          {% endfor %}
      
      .. tab:: Modifying Fields
                
        :class:`~.ModelGrid` inherits several useful methods for creating new data
        fields and adding new data fields to a ModelGrid instance. Methods to add or
        create a new data array follow the ``numpy`` syntax for creating arrays. The
        folowing methods create and, optionally, initialize new arrays. These arrays
        are of the correct size but a new field will not be added to the field:
        
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["field-io"] %}
            ~{{func}}      
          {% endfor %}
      
      .. tab:: Field
          
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["info-field"] %}
            ~{{func}}      
          {% endfor %}
      
      .. tab:: Boundary
              
        .. autosummary::
          :nosignatures:
        
          {% for func in cats["boundary-condition"] %}
            ~{{func}}      
          {% endfor %}
    
  {% endfor %}