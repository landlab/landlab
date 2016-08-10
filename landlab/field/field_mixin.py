#! /usr/bin/env python
from .grouped import ModelDataFields, GroupSizeError


_GROUPS = ('node', 'cell', 'link', 'face', 'core_node', 'core_cell',
           'active_link', 'active_face', )


class ModelDataFieldsMixIn(ModelDataFields):

    """Mix-in that provides un-sized fields.

    Inherit from this class to provide fields that do not have to be given
    a size on instantiation. The size of a particular group will only be set
    after the first time a group is accessed (a field is added to a group or
    an array is created for a group). The size of that group will be the
    size as returned by the `number_of_elements` method.

    This mix-in assumes it is being added to a class that implements a method
    that, given an element name as a string, returns the number of elements
    for that group. A `RasterModelGrid` is an excellent example of this.

    >>> from landlab import RasterModelGrid
    >>> grid = RasterModelGrid((4, 5))
    >>> grid.number_of_elements('node')
    20

    Examples
    --------

    A `RasterModelGrid` is an example of a class that inherits from
    `ModelDataFieldsMixIn`.

    >>> from landlab import RasterModelGrid, CLOSED_BOUNDARY
    >>> grid = RasterModelGrid((4, 5))

    After we initially create the grid, it has a number of active links.
    However, since we haven't assigned or created a field based on the number,
    the field is un-sized.

    >>> grid.number_of_elements('active_link')
    17
    >>> grid.at_active_link.size is None
    True

    If we change the status of some of the nodes, the number of active links
    changes. If we now create an array based on the current number of active
    links, that number *is now fixed*.

    >>> grid.status_at_node[:5] = CLOSED_BOUNDARY
    >>> grid.number_of_elements('active_link')
    14
    >>> grid.ones(at='active_link', dtype=int)
    array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])
    >>> grid.at_active_link.size
    14

    Because the number of active links is fixed in our collection of fields,
    if we change the number of active links in our grid, an exception is
    raised if we try to create a new active link array.

    >>> grid.status_at_node[-5:] = CLOSED_BOUNDARY
    >>> grid.number_of_elements('active_link')
    11
    >>> grid.ones(at='active_link', dtype=int)
    ...     # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    GroupSizeError: number of active_link elements has changed. (was = 14, now=11)
    """

    def __init__(self, **kwds):
        super(ModelDataFieldsMixIn, self).__init__(**kwds)
        for group in _GROUPS:
            ModelDataFields.new_field_location(self, group)

    def new_field_location(self, group, size=None):
        raise AttributeError(
            "'ModelDataFieldsMixIn' object has no attribute "
            "'new_field_location'")

    def empty(self, *args, **kwds):
        """Array, filled with unititialized values, for a given element.

        Returns a numpy array of uninitialized values that is the same length
        as the number of nodes in the grid. Use the *centering* keyword to
        return an array for other elements of the grid. *centering* is a
        string that is one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> len(grid.empty())
        20
        """
        if len(args) == 0:
            group = kwds.pop('at', kwds.pop('centering', 'node'))
        else:
            group = args[0]

        n_elements = self.number_of_elements(group)

        if self[group].size is None:
            self[group].size = n_elements
        elif self[group].size != n_elements:
            raise GroupSizeError(group, self[group].size, n_elements)

        return ModelDataFields.empty(self, group, **kwds)

    def ones(self, *args, **kwds):
        """Array, filled with ones, for a given element.

        Returns a numpy array of ones that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.zeros(dtype=int) # doctest: +NORMALIZE_WHITESPACE
        array([0, 0, 0, 0, 0,
               0, 0, 0, 0, 0,
               0, 0, 0, 0, 0,
               0, 0, 0, 0, 0])
        >>> grid.zeros('cell') # doctest: +NORMALIZE_WHITESPACE
        array([ 0., 0., 0.,
                0., 0., 0.])
        """
        if len(args) == 0:
            group = kwds.pop('at', kwds.pop('centering', 'node'))
        else:
            group = args[0]

        n_elements = self.number_of_elements(group)

        if self[group].size is None:
            self[group].size = n_elements
        elif self[group].size != n_elements:
            raise GroupSizeError(group, self[group].size, n_elements)

        return ModelDataFields.ones(self, group, **kwds)

    def zeros(self, *args, **kwds):
        """Array, filled with zeros, for a given element.

        Returns a numpy array of zeros that is the same length as the number
        of nodes in the grid. Use the *centering* keyword to return an
        array for other elements of the grid. *centering* is a string that is
        one of *node*, *cell*, *link*, or *face*.

        All other keywords are the same as for the numpy zeros function.

        Parameters
        ----------
        centering : str, optional
            Grid element on which the values are defined.

        Returns
        -------
        ndarray
            A newly-allocated array.

        Examples
        --------
        >>> from landlab import RasterModelGrid
        >>> grid = RasterModelGrid((4, 5))
        >>> grid.zeros()
        array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
                0.,  0.,  0.,  0.,  0.,  0.,  0.])
        """
        if len(args) == 0:
            group = kwds.pop('at', kwds.pop('centering', 'node'))
        else:
            group = args[0]

        n_elements = self.number_of_elements(group)

        if self[group].size is None:
            self[group].size = n_elements
        elif self[group].size != n_elements:
            raise GroupSizeError(group, self[group].size, n_elements)

        return ModelDataFields.zeros(self, group, **kwds)
