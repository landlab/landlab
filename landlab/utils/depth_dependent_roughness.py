"""
Creates a field of Manning's n value, where each value is dependent on the
water depth at a given node (e.g. Jain et al., 2005, Mugler et al., 2011 and
Rengers et al., 2016). This "effectively creates separate values of Manning's n
for hillslopes and channels" (Rengers et al., 2016).

This component creates -or- overwrites a field on the grid called
'mannings_n' where each node is assigned a Manning's n value based on the
minimum Manning's n value for the landscape, the local water depths, an index
(or threshold) water depth above which all Manning's n values are considered
constant, and a vegetation drag coefficent (for more on vegetation drag and
the impact on surface roughness, see Wu et al., 1999 in the Journal of
Hydraulic Engineering.)

This can be used iteratively inside a driver loop, to update Manning's n values
as water depths change in another component (e.g. OverlandFlow)

.. codeauthor:: Jordan Adams

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> grid = RasterModelGrid((5, 5))
>>> grid.at_node['surface_water__depth'] = np.array(
...     [5., 5., 5., 5., 5.,
...     4., 4., 4., 4., 4.,
...     3., 3., 3., 3., 3.,
...     2., 2., 2., 2., 2.,
...     1., 1., 1., 1., 1.])
>>> depth_dependent_mannings_n(grid, index_flow_depth=2.)
>>> grid.at_node['mannings_n']
array([ 0.06      ,  0.06      ,  0.06      ,  0.06      ,  0.06      ,
        0.06      ,  0.06      ,  0.06      ,  0.06      ,  0.06      ,
        0.06      ,  0.06      ,  0.06      ,  0.06      ,  0.06      ,
        0.06      ,  0.06      ,  0.06      ,  0.06      ,  0.06      ,
        0.07559526,  0.07559526,  0.07559526,  0.07559526,  0.07559526])
"""

import numpy as np

from landlab import FieldError


def depth_dependent_mannings_n(
    grid,
    water_depths="surface_water__depth",
    min_mannings_n=0.06,
    index_flow_depth=0.003,
    veg_drag_exponent=(-1.0 / 3.0),
):
    """
    Method to create or overwrite a Manning's n field
    (grid.at_node['mannings_n']) with spatially variable n values based on
    water depths.

    Parameters
    ----------
    grid : A Landlab RasterModelGrid instance
        A Landlab grid - only works with RasterModelGrid instances as of
        1/31/17.
    water_depths : array or Landlab field of floats
        Array of values, with length of number of nodes, water depths
        at all grid node locations. (m)
    min_mannings_n : float
        This is the minimum Manning's n coefficient for a given landscape,
        following Chow, 1959. (s m^(-1./3.))
    index_flow_depth : float
        The flow depth above which it is assumed that Manning's n is
        constant. (m)
    veg_drag_exponent :
        An exponent related to vegetation drag effects, which increases
        effective Manning's n at low flow conditions.
    """

    # Looks for a field called 'mannings_n' attached to the grid instance. If
    # one is found, a FieldError is thrown but ignored. This method
    # REWRITES over the existing Manning's n fields after the calculation.
    try:
        grid.add_zeros("mannings_n", at="node")

    except FieldError:
        pass

    # Identifies locations where water depth is lower than the value supplied
    # through keyword index_flow_depth.
    (locs_less,) = np.where(grid.at_node["surface_water__depth"] <= index_flow_depth)

    # Identifies locations where water depth is greater than the value
    # supplied through keyword index_flow_depth.
    (locs_more,) = np.where(grid.at_node["surface_water__depth"] > index_flow_depth)

    # At all locations lower than the index flow depth (assumed to be shallow
    # flow on hillslopes), a new Manning's n value is calculated to that
    # incorporates effects of vegetation drag. These Manning's n values will
    # be greater than the supplied Manning's n (keyword: min_mannings_n)
    grid.at_node["mannings_n"][locs_less] = (
        min_mannings_n
        * (grid.at_node["surface_water__depth"][locs_less] / index_flow_depth)
        ** veg_drag_exponent
    )

    # Resets the field with the new Manning's n values.
    grid.at_node["mannings_n"][locs_more] = min_mannings_n
