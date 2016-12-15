#! /usr/env/python
"""

Example of a simple diffusion model that uses the LinearDiffuser.

Created July 2013 GT
Last updated August 2013 GT

"""

from landlab.components.diffusion import LinearDiffuser
from landlab.grid import create_and_initialize_grid
from landlab import ModelParameterDictionary
import pylab
import numpy


def display_model(grid, elevation):

    # Convert z to a 2D array
    zr = grid.node_vector_to_raster(elevation, flip_vertically=True)

    # Create a shaded image
    pylab.close()  # clear any pre-existing plot
    im = pylab.imshow(zr, cmap=pylab.cm.RdBu)  # display a colored image

    # add contour lines with labels
    # add a color bar on the side
    if numpy.amax(zr)>numpy.amin(zr):
        cset = pylab.contour(zr)
        pylab.clabel(cset, inline=True, fmt='%1.1f', fontsize=10)
        pylab.colorbar(im)

    # add a title
    pylab.title('Simulated topography with uplift and diffusion')

    # Display the plot
    pylab.show()


def main():

    # INITIALIZE

    # Name of parameter input file
    input_file_name = 'test_inputs_for_diffusion_model.txt'

    # Open input file and read run-control parameters
    mpd = ModelParameterDictionary(input_file_name)
    run_duration = mpd.get('RUN_DURATION', ptype=float)

    # Create and initialize a grid
    mg = create_and_initialize_grid(mpd)

    # Create and initialize a diffusion component
    dc = LinearDiffuser(mg)
    dc.initialize(mpd)

    # RUN

    # Run the diffusion component until it's time for the next output
    dc.run_until(run_duration)

    # FINALIZE

    # Display results to screen
    mg.imshow('node', 'landscape_surface__elevation')
    import pylab
    pylab.show()

    from landlab.io.netcdf import write_netcdf
    write_netcdf('diffusion_example.nc', mg)


if __name__ == "__main__":
    main()
