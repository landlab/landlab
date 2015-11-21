#! /usr/bin/env python
from .grouped import ModelDataFields


_GROUPS = ('node', 'cell', 'link', 'face', 'core_node', 'core_cell',
           'active_link', 'active_face', )


class ModelDataFieldsMixIn(ModelDataFields):
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

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
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

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
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

        if self[group].size is None:
            self[group].size = self.number_of_elements(group)
        return ModelDataFields.zeros(self, group, **kwds)
