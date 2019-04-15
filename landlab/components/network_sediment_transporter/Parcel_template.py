import numpy as np

from landlab.data_record import DataRecord


class Parcel(DataRecord):

    # Your collection will inherit from ItemCollection so:
    #
    # - it will be instantiated like ItemCollection, which takes as input the
    # grid, a data dictionary, the grid element your items will live on, and
    # the element id where the items are (initially). See ItemCollection
    # documentation for more.
    #
    # - you can add as many fields as is appropriate for your collection in
    # the data dictionary
    #
    # - it will inherit the methods that are defined in ItemCollection
    # (_check_sizes, etc.) AND all methods and attributes related to
    # Pandas DataFrames (at, groupby, etc.).

    """
    Documentation and tests of MyCollection here
    """

    def __init__(self, grid, param_1=[], param_2=[], param_3=[], param_4=[]):
        """
        Documentation on init here
        """

        # Store the grid:
        self._grid = grid

        # New parameters could be defined here:
        param_5 = param_2 + param_3

        # All fields should have the same length (one value for each item):
        self._nb_of_items = len(param_1)

        # Create the data dictionary:
        my_data = {
            "data_name_1": param_1,
            "data_name_2": param_2,
            "data_name_3": param_3,
            "data_name_5": param_5,
            "data_name_6": np.zeros(self._nb_of_items),
            "data_name_7": np.full(self._nb_of_items, np.NaN),
        }

        # In this example, the positions (element ID) of the items at the
        # time of initialization are stored in the array param_4:

        my_element_id = param_4

        # If your items leave on nodes and you know the initial (x,y)
        # coordinates of your items, you can use the grid method
        # 'find_nearest_node' to find the node ID corresponding to the closest
        # node to each item coordinates (only for RasterModelGrid, only for
        # nodes):

        # my_element_id = self._grid.find_nearest_node(
        #         (clast_x[:], clast_y[:]))

        # Build ItemCollection containing your data:
        DataRecord.__init__(
            self,
            self._grid,
            data=my_data,
            grid_element="node",
            element_id=my_element_id,
        )

        # In this example, all elements live on nodes.
        # If not all your elements live on nodes, you can pass an array (or a
        # parameter in __init__) instead of 'node' in grid_element above.
        # E.g. in __init__(..., param_8=['node', 'link, 'link', 'node',...])
        # and then in the lines above ItemCollection.__init__(self,...
        #                                               grid_element=param_8,
        #                                                ...)

        # Using python 3:

    #        super().__init__(
    #                self,
    #                self._grid,
    #                data=my_data,
    #                grid_element='node',
    #                element_id=my_element_id)

    def method_1(self, other_input):
        """
        Documentation on method_1 here
        """
        # Do things
