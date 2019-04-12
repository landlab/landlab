

def plot_network_average(grid, parcels, time, color=attribute_name, line_width=xxx, **kwargs):
  this would plot lines from the saved shapefile attrite on grid and make
  color based on the given attribute.
  kwargs would have to do with color map, etc

  have multiple options to save or return plot instance.


def plot_network_parcels(grid, parcels, time, color=xxx, shape=xxx, size=xxs):
  this would plot lines from saved shapefile,
  then calculate distance down each grid segment,
  then place each parcel along the squiggly shapefile line correctly with the
  correct color, size, etc.

  have multiple options to save or return plot instance.


def other_plots():!
    # something that works nicely with bokeh or holloviews and has a time sliderbar.
