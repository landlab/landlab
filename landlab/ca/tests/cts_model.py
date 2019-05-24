#!/usr/env/python
import time

from matplotlib.pyplot import axis
from numpy import random

from landlab.ca.celllab_cts import CAPlotter, Transition
from landlab.io.native_landlab import save_grid

_DEBUG = False


class CTSModel(object):
    """
    Implement a generic CellLab-CTS model.

    This is the base class from which models should inherit.
    """

    def __init__(
        self,
        grid_size=(5, 5),
        report_interval=5.0,
        grid_orientation="vertical",
        node_layout="rect",
        show_plots=False,
        cts_type="oriented_hex",
        run_duration=1.0,
        output_interval=1.0e99,
        plot_every_transition=False,
        **kwds
    ):

        self.initialize(
            grid_size,
            report_interval,
            grid_orientation,
            node_layout,
            show_plots,
            cts_type,
            run_duration,
            output_interval,
            plot_every_transition,
            **kwds
        )

    def initialize(
        self,
        grid_size=(5, 5),
        report_interval=5.0,
        grid_orientation="vertical",
        node_layout="rect",
        show_plots=False,
        cts_type="oriented_hex",
        run_duration=1.0,
        output_interval=1.0e99,
        plot_every_transition=False,
        **kwds
    ):

        # Remember the clock time, and calculate when we next want to report
        # progress.
        self.current_real_time = time.time()
        self.next_report = self.current_real_time + report_interval
        self.report_interval = report_interval

        # Interval for output
        self.output_interval = output_interval

        # Duration for run
        self.run_duration = run_duration

        # Create a grid
        self.create_grid_and_node_state_field(
            grid_size[0], grid_size[1], grid_orientation, node_layout, cts_type
        )

        # Create the node-state dictionary
        ns_dict = self.node_state_dictionary()

        # Initialize values of the node-state grid
        nsg = self.initialize_node_state_grid()

        # Create the transition list
        xn_list = self.transition_list()

        # Create the CA object
        if cts_type == "raster":
            from landlab.ca.raster_cts import RasterCTS

            self.ca = RasterCTS(self.grid, ns_dict, xn_list, nsg)
        elif cts_type == "oriented_raster":
            from landlab.ca.oriented_raster_cts import OrientedRasterCTS

            self.ca = OrientedRasterCTS(self.grid, ns_dict, xn_list, nsg)
        elif cts_type == "hex":
            from landlab.ca.hex_cts import HexCTS

            self.ca = HexCTS(self.grid, ns_dict, xn_list, nsg)
        else:
            from landlab.ca.oriented_hex_cts import OrientedHexCTS

            self.ca = OrientedHexCTS(self.grid, ns_dict, xn_list, nsg)

        # Initialize graphics
        self._show_plots = show_plots
        if show_plots:
            self.initialize_plotting(**kwds)

    def create_grid_and_node_state_field(
        self, num_rows, num_cols, grid_orientation, node_layout, cts_type
    ):
        """Create the grid and the field containing node states."""

        if cts_type == "raster" or cts_type == "oriented_raster":
            from landlab import RasterModelGrid

            self.grid = RasterModelGrid(shape=(num_rows, num_cols), xy_spacing=1.0)
        else:
            from landlab import HexModelGrid

            self.grid = HexModelGrid(
                num_rows,
                num_cols,
                xy_spacing=1.0,
                orientation=grid_orientation,
                node_layout=node_layout,
            )

        self.grid.add_zeros("node", "node_state", dtype=int)

    def node_state_dictionary(self):
        """Create and return a dictionary of all possible node (cell) states.

        This method creates a default set of states (just two); it is a
        template meant to be overridden.
        """
        ns_dict = {0: "on", 1: "off"}
        return ns_dict

    def transition_list(self):
        """Create and return a list of transition objects.

        This method creates a default set of transitions (just two); it is a
        template meant to be overridden.
        """
        xn_list = []
        xn_list.append(Transition((0, 1, 0), (1, 0, 0), 1.0))
        xn_list.append(Transition((1, 0, 0), (0, 1, 0), 1.0))
        return xn_list

    def write_output(self, grid, outfilename, iteration):
        """Write output to file (currently netCDF)."""
        filename = outfilename + str(iteration).zfill(4) + ".nc"
        save_grid(grid, filename)

    def initialize_node_state_grid(self):
        """Initialize values in the node-state grid.

        This method should be overridden. The default is random "on" and "off".
        """
        num_states = 2
        for i in range(self.grid.number_of_nodes):
            self.grid.at_node["node_state"][i] = random.randint(num_states)
        return self.grid.at_node["node_state"]

    def initialize_plotting(self, **kwds):
        """Create and configure CAPlotter object."""
        self.ca_plotter = CAPlotter(self.ca, **kwds)
        self.ca_plotter.update_plot()
        axis("off")

    def run_for(self, dt):

        self.ca.run(self.ca.current_time + dt, self.ca.node_state)


if __name__ == "__main__":
    ctsm = CTSModel(show_plots=True)
    ctsm.run_for(1.0)
    ctsm.ca_plotter.update_plot()
