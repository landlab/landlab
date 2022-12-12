Fields
======

:class:`~.ModelGrid` inherits several useful methods for creating new data
fields and adding new data fields to a :class:`~.ModelGrid` instance. Methods to add or
create a new data array follow the ``numpy`` syntax for creating arrays. The
following methods create and, optionally, initialize new arrays. The size of the
new array is determined by the *at* keyword, which indicates on which element
the array is defined. Methods with the prefix ``add_`` will add the newly created
field to the grid, otherwise a newly-created array is **not** added to the grid.

.. tab:: Create

    .. autosummary::
        :nosignatures:

        ~landlab.grid.base.ModelGrid.add_empty
        ~landlab.grid.base.ModelGrid.add_field
        ~landlab.grid.base.ModelGrid.add_ones
        ~landlab.grid.base.ModelGrid.add_zeros
        ~landlab.grid.base.ModelGrid.delete_field
        ~landlab.grid.base.ModelGrid.empty
        ~landlab.grid.base.ModelGrid.ones
        ~landlab.grid.base.ModelGrid.zeros

.. tab:: Access

    .. autosummary::
        :nosignatures:

        ~landlab.grid.base.ModelGrid.at_cell
        ~landlab.grid.base.ModelGrid.at_corner
        ~landlab.grid.base.ModelGrid.at_face
        ~landlab.grid.base.ModelGrid.at_link
        ~landlab.grid.base.ModelGrid.at_node
        ~landlab.grid.base.ModelGrid.at_patch
        ~landlab.grid.base.ModelGrid.field_values
        ~landlab.grid.base.ModelGrid.return_array_or_field_values

.. tab:: Info

    .. autosummary::
        :nosignatures:

        ~landlab.grid.base.ModelGrid.size
        ~landlab.grid.base.ModelGrid.keys
        ~landlab.grid.base.ModelGrid.has_group
        ~landlab.grid.base.ModelGrid.has_field
        ~landlab.grid.base.ModelGrid.field_units
        ~landlab.grid.base.ModelGrid.field_values
        ~landlab.grid.base.ModelGrid.groups
