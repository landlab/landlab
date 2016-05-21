#! /usr/bin/env python

"""Create mp4 animation of Landlab output.

This component allows creation of mp4 animations of output from Landlab.
It does so by stitching together output from the conventional Landlab
static plotting routines from plot/imshow.py.

It is compatible with all Landlab grids, though cannot handle an evolving grid
as yet.

Initialize the video object vid at the start of your code, then simply call
vid.add_frame(grid, data) each time you want to add a frame. At the end of
the model run, call vid.produce_video().

CAUTION: This component may prove *very* memory-intensive. It is recommended
that the total number of frames included in the output multiplied by the
number of pixels (nodes) in the image not exceed XXXXXXXXX.

Due to some issues with codecs in matplotlib, at the moment on .gif output
movies are recommended. If this irritates you, you can modify your own
PYTHONPATH to allow .mp4 compilation (try a google search for the warning
raised by this method for some hints). These (known) issues are apparently
likely to resolve themselves in a future release of matplotlib.
"""
import six
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from landlab.plot import imshow


class VideoPlotter(object):

    """
    Create animations of landlab output.

    Create Landlab movies.

    Parameters
    ----------
    grid : RasterModelGrid
        A RasterModelGrid.
    data_centering : {'node', 'cell'}, optional
        Where data are centered.
    start : float, optional
        Model time at which filming starts.
    stop : float, optional
        Model time at which filming stops.
    step : float, optional
        Model time frequency at which frames are made.
    """

    def __init__(self, grid, data_centering='node', start=None, stop=None,
                 step=None):
        """Create Landlab movies.

        Parameters
        ----------
        grid : RasterModelGrid
            A RasterModelGrid.
        data_centering : {'node', 'cell'}, optional
            Where data are centered.
        start : float, optional
            Model time at which filming starts.
        stop : float, optional
            Model time at which filming stops.
        step : float, optional
            Model time frequency at which frames are made.
        """
        self.initialize(grid, data_centering, start, stop, step)

    def initialize(self, grid, data_centering, start, stop, step):
        """Set up the plotter.

        A copy of the grid is required.

        *data_centering* controls the type of data the video will be plotting.

        It can be set to: `"node"` (default) or `"cell"`.

        *start*, *stop,* and *step* control when a frame is added. They are
        absolute times in the model run. All are optional.

        Parameters
        ----------
        grid : RasterModelGrid
            A RasterModelGrid.
        data_centering : {'node', 'cell'}, optional
            Where data are centered.
        start : float
            Model time at which filming starts.
        stop : float
            Model time at which filming stops.
        step : float
            Model time frequency at which frames are made.
        """
        options_for_data_centering = ['node', 'cell']

        if data_centering not in options_for_data_centering:
            raise ValueError('data_centering not valid')

        self.grid = grid
        #self.image_list = []
        self.data_list = []

        # this controls the intervals at which to plot
        self.last_remainder = float('inf')
        self.last_t = float('-inf')
        if start is None:
            start = float('-inf')
        if stop is None:
            stop = float('inf')
        self.step_control_tuple = (start, stop, step)

        # initialize the plots for the vid
        if data_centering == 'node':
            self.centering = 'n'
            self.plotfunc = imshow.imshow_grid_at_node
        elif data_centering == 'cell':
            self.centering = 'c'
            self.plotfunc = imshow.imshow_cell_grid

        self.randomized_name = ("my_animation_" +
                                str(int(np.random.random() * 10000)))
        self.fig = plt.figure(self.randomized_name)  # randomized name

    def add_frame(self, grid, data, elapsed_t, **kwds):
        """Add a frame to the video.

        data can be either the data to plot (nnodes, or appropriately lengthed
        numpy array), or a string for grid field access.

        kwds can be any of the usual plotting keywords, e.g., cmap.
        """
        if type(data) == str:
            if self.centering == 'n':
                data_in = grid.at_node[data]
            elif self.centering == 'c':
                data_in = grid.at_cell[data]
        else:
            data_in = data

        self.kwds = kwds

        if self.last_t < elapsed_t:
            try:
                normalized_elapsed_t = elapsed_t - self.start_t
            except AttributeError:
                self.start_t = elapsed_t
                normalized_elapsed_t = 0.
        else:  # time has apparently gone "backwards"; reset the module
            #...note a *forward* jump in time wouldn't register
            self.clear_module()
            self.start_t = elapsed_t
            normalized_elapsed_t = 0.

        # we're between start & stop
        if self.step_control_tuple[0] <= elapsed_t < self.step_control_tuple[1]:
            if not self.step_control_tuple[2]:  # no step provided
                six.print_('Adding frame to video at elapsed time %f' %
                           elapsed_t)
                self.data_list.append(data_in.copy())
            else:
                excess_fraction = (normalized_elapsed_t %
                                   self.step_control_tuple[2])
                # Problems with rounding errors make this double check
                # necessary
                if excess_fraction < self.last_remainder or np.allclose(excess_fraction, self.step_control_tuple[2]):
                    six.print_(
                        'Adding frame to video at elapsed time %f' % elapsed_t)
                    self.data_list.append(data_in.copy())
                self.last_remainder = excess_fraction
        self.last_t = elapsed_t

    def produce_video(self, interval=200, repeat_delay=2000,
                      filename='video_output.gif', override_min_max=None):
        """Finalize and save the video of the data.

        Parameters
        ----------
        interval : int, optional
            Interval between frames in milliseconds.
        repeat_delay : int, optional
            Repeat delay before restart in milliseconds.
        filename : str, optional
            Name of the file to save in the present working directory. At
            present, only *.gifs* will implement reliably without
            tweaking Python's PATHs.
        override_min_max : tuple of float
            Minimum and maximum for the scale on the plot as (*min*, *max*).
        """
        six.print_("Assembling video output, may take a while...")
        plt.figure(self.randomized_name)
        # find the limits for the plot:
        if not override_min_max:
            self.min_limit = np.amin(self.data_list[0])
            self.max_limit = np.amax(self.data_list[0])

            if len(self.data_list) <= 1:
                raise ValueError('Animation must have at least one frame.')

            # assumes there is more than one frame in the loop
            for i in self.data_list[1:]:
                self.min_limit = min((self.min_limit, np.amin(i)))
                self.max_limit = max((self.max_limit, np.amax(i)))
        else:
            self.min_limit = override_min_max[0]
            self.max_limit = override_min_max[1]

        self.fig.colorbar(
            self.plotfunc(self.grid, self.data_list[0],
                          limits=(self.min_limit, self.max_limit),
                          allow_colorbar=False, **self.kwds))

        ani = animation.FuncAnimation(self.fig, _make_image,
                                      frames=self._yield_image,
                                      interval=interval, blit=True,
                                      repeat_delay=repeat_delay)
        ani.save(filename, fps=1000. / interval)
        plt.close()

    def _yield_image(self):
        """
        Helper function designed to generate image_list items for plotting,
        rather than storing them all.
        """

        for i in self.data_list:
            # yield self.grid.node_vector_to_raster(i)
            yield (i, self.plotfunc, (self.min_limit, self.max_limit),
                   self.grid, self.kwds)

    def clear_module(self):
        """Clear internally held data.

        Wipe all internally held data that would cause trouble if module
        were to be rerun without being reinstantiated.
        """
        self.data_list = []


def _make_image(yielded_tuple):
    yielded_raster_data = yielded_tuple[0]
    plotfunc = yielded_tuple[1]
    limits_in = yielded_tuple[2]
    grid = yielded_tuple[3]
    kwds = yielded_tuple[4]
    im = plotfunc(grid, yielded_raster_data, limits=limits_in,
                  allow_colorbar=False, **kwds)
    return im
