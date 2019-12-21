#! /usr/bin/env python
from .grouped import GroupSizeError, ModelDataFields


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
    """

    def __init__(self, **kwds):
        super(ModelDataFieldsMixIn, self).__init__(**kwds)

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

        LLCATS: FIELDADD
        """
        if len(args) == 0:
            group = kwds.pop("at", kwds.pop("centering", "node"))
        else:
            group = args[0]

        if group == "grid":
            raise ValueError(
                "empty is not supported for at='grid', if you "
                "want to create a field at the grid, use\n"
                "grid.at_grid['value_name']=value\n"
                "instead."
            )

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

        LLCATS: FIELDADD
        """
        if len(args) == 0:
            group = kwds.pop("at", kwds.pop("centering", "node"))
        else:
            group = args[0]

        if group == "grid":
            raise ValueError(
                "ones is not supported for at='grid', if you "
                "want to create a field at the grid, use\n"
                "grid.at_grid['value_name']=value\n"
                "instead.\nAlternatively, if you want ones"
                "of the shape stored at_grid, use np.array(1)."
            )

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

        LLCATS: FIELDADD
        """
        if len(args) == 0:
            group = kwds.pop("at", kwds.pop("centering", "node"))
        else:
            group = args[0]

        if group == "grid":
            raise ValueError(
                "zeros is not supported for at='grid', if you "
                "want to create a field at the grid, use\n"
                "grid.at_grid['value_name']=value\n"
                "instead.\nAlternatively, if you want zeros"
                "of the shape stored at_grid, use np.array(0)."
            )

        n_elements = self.number_of_elements(group)

        if self[group].size is None:
            self[group].size = n_elements
        elif self[group].size != n_elements:
            raise GroupSizeError(group, self[group].size, n_elements)

        return ModelDataFields.zeros(self, group, **kwds)
