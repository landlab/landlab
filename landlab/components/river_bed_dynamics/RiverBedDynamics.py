"""Landlab component that calculates temporal and spatial changes in river
bed elevation and grain size distribution. Also, this component predicts
fractional or total sediment transport based on the specified bed load
transport model. Hydraulic properties are obtained from an external flow
component. We recommend coupling it with the OverlandFlow component from
Adams et al, 2017.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

Examples
--------

Let's import all the required libraries

>>> import numpy as np
>>> from landlab import RasterModelGrid, imshow_grid
>>> from landlab.components import RiverBedDynamics
>>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

Create a grid on which to calculate sediment transport

>>> grid = RasterModelGrid((5, 5))

The grid will need some data to run the RiverBedDynamics component.
To check the names of the fields that provide input use the *input_var_names*
class property.

>>> RiverBedDynamics.input_var_names
('surface_water__depth', 'surface_water__velocity', 'topographic__elevation')

Create fields of data for each of these input variables. When running a
complete simulation some of these variables will be created by the flow model.
Notice that surface water depth and velocity are required at links. However,
specifying these variables at nodes is easier and then we can map the fields
onto links. By doing so, we don't have to deal with links numbering. When this
component is coupled to OverlandFlow there is no need to map fields because it
is done automatically within the component.

We start by creating the topography data

>>> grid.at_node["topographic__elevation"] = [
...     [1.07, 1.06, 1.00, 1.06, 1.07],
...     [1.08, 1.07, 1.03, 1.07, 1.08],
...     [1.09, 1.08, 1.07, 1.08, 1.09],
...     [1.09, 1.09, 1.08, 1.09, 1.09],
...     [1.09, 1.09, 1.09, 1.09, 1.09],
... ]

We set the boundary conditions

>>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

And check the node status

>>> grid.status_at_node.reshape(grid.shape)
array([[4, 4, 1, 4, 4],
       [4, 0, 0, 0, 4],
       [4, 0, 0, 0, 4],
       [4, 0, 0, 0, 4],
       [4, 4, 4, 4, 4]], dtype=uint8)

Which tell us that there is one outlet located on the 3rd node

The topography data can be display using

>>> imshow_grid(grid, "topographic__elevation")

Now we add some water into the watershed. In this case is specified in nodes

>>> grid.at_node["surface_water__depth"] = [
...     [0.102, 0.102, 0.102, 0.102, 0.102],
...     [0.102, 0.102, 0.102, 0.102, 0.102],
...     [0.102, 0.102, 0.102, 0.102, 0.102],
...     [0.102, 0.102, 0.102, 0.102, 0.102],
...     [0.102, 0.102, 0.102, 0.102, 0.102],
... ]

Now, we give the water a velocity.

>>> grid.at_node["surface_water__velocity"] = [
...     [0.25, 0.25, 0.25, 0.25, 0.25],
...     [0.25, 0.25, 0.25, 0.25, 0.25],
...     [0.25, 0.25, 0.25, 0.25, 0.25],
...     [0.25, 0.25, 0.25, 0.25, 0.25],
...     [0.25, 0.25, 0.25, 0.25, 0.25],
... ]

Note that in this example, we are attempting to specify a vector at a node
using a single value. This is done intentionally to emphasize the process.
The component will interpret this as the vector's magnitude, and, given its
location in the grid, it will manifest different components. When using
OverlandFlow, there is no need to specify a velocity because it is a
byproduct of the component.

For the purpose of this illustration, we will make an assumption that the
conditions remain identical to the previous time step.

By default, when creating our grid we used a spacing of 1 m in the x and y
directions. Therefore, the discharge is 0.0255 m3/s. Discharge is always in
units of m3/s.

Now we map nodes into links when it is required

>>> grid.at_link["surface_water__depth"] = map_mean_of_link_nodes_to_link(
...     grid, "surface_water__depth"
... )
>>> grid.at_link["surface_water__velocity"] = map_mean_of_link_nodes_to_link(
...     grid, "surface_water__velocity"
... )

We will assume, for the sake of demonstration, that we have two sectors with
different bed surface grain size (GSD). We can tell the component the location
of these two different GSD within the watershed (labeled as 0 and 1). This will
be included during the instantiation

>>> gsd_loc = [
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
...     [0, 1.0, 1.0, 1.0, 0],
... ]

We assign a GSD to each location. The structure of this array is:
First column contains the grain sizes in milimiters
Second column is location 0 in 'bed_grainSizeDistribution__location'
Third column is location 1 in 'bed_grainSizeDistribution__location', and so on

>>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]

Instantiate the `RiverBedDynamics` component to work on the grid, and run it.

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Parker1990",
...     bed_surface__gsd_location_node=gsd_loc,
...     output_vector=True,
...     track_stratigraphy=True,
... )
>>> rbd.run_one_step()

After running RiverBedDynamics, we can check if the different GSD locations
were correctly included

>>> gsd_loc_rbd = rbd._bed_surface__gsd_location_node.reshape(grid.shape)
>>> gsd_loc_rbd
array([[ 0.,  1.,  1.,  1.,  0.],
       [ 0.,  1.,  1.,  1.,  0.],
       [ 0.,  1.,  1.,  1.,  0.],
       [ 0.,  1.,  1.,  1.,  0.],
       [ 0.,  1.,  1.,  1.,  0.]])

Let's check at the calculated net bedload. Hereinafter, most of the results will show
between 3 to 6 decimals. This is just to avoid roundoff errors problems

>>> qb = rbd._sed_transport__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001598, -0.004793,  0.001598,  0.      ],
       [ 0.      ,  0.      ,  0.001597,  0.      ,  0.      ],
       [ 0.      ,  0.      , -0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Let's check at the calculated shear_stress

>>> shearStress = rbd._surface_water__shear_stress_link
>>> np.round(shearStress, decimals=3)
array([  0.   ,  60.017,  60.017,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,  40.011,  40.011,   0.   ,  10.003,
        10.003,  40.011,  10.003,  10.003,   0.   ,  10.003,  10.003,
         0.   ,   0.   ,  10.003,  10.003,  10.003,   0.   ,   0.   ,
        10.003,  10.003,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ])

Considering the link upstream the watershed exit, link Id 15, we can obtain the
bed load transport rate

>>> qb_l15 = rbd._sed_transport__bedload_rate_link[15]
>>> formatted_number = "{:.4e}".format(qb_l15)
>>> print(formatted_number)
-1.5977e-03

Therefore, the bed load transport rate according to Parker 1990 surface-based
equation is 1.598 * 10^-3 m2/s. Negative means that is going in the negative
Y direction

The GSD at this place is:

>>> qb_gsd_l15 = rbd._sed_transport__bedload_gsd_link[15]
>>> np.round(qb_gsd_l15, decimals=3)
array([ 0.475,  0.525])

Which in cummulative percentage is equivalent to

==== =======
D mm % Finer
==== =======
32   100.000
16   52.498
8    0.000
==== =======


Grain sizes are always given in mm.
We can also check the bed load grain size distribution in all links

>>> qb_gsd_l = rbd._sed_transport__bedload_gsd_link
>>> np.round(qb_gsd_l, decimals=3)
array([[ 0.   ,  0.   ],
       [ 0.485,  0.515],
       [ 0.485,  0.515],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.475,  0.525],
       [ 0.475,  0.525],
       [ 0.   ,  0.   ],
       [ 0.541,  0.459],
       [ 0.282,  0.718],
       [ 0.475,  0.525],
       [ 0.282,  0.718],
       [ 0.541,  0.459],
       [ 0.   ,  0.   ],
       [ 0.282,  0.718],
       [ 0.282,  0.718],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.282,  0.718],
       [ 0.282,  0.718],
       [ 0.282,  0.718],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.282,  0.718],
       [ 0.282,  0.718],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ],
       [ 0.   ,  0.   ]])

Zeros indicate that there is no sediment transport of that grain size
at that location.

After the flow acted on the bed and sediment transport occured we can check
the new topographic elevation field

>>> z = grid.at_node["topographic__elevation"].reshape(grid.shape)
>>> np.round(z, decimals=3)
array([[ 1.07 ,  1.06 ,  1.007,  1.06 ,  1.07 ],
       [ 1.08 ,  1.068,  1.037,  1.068,  1.08 ],
       [ 1.09 ,  1.08 ,  1.068,  1.08 ,  1.09 ],
       [ 1.09 ,  1.09 ,  1.08 ,  1.09 ,  1.09 ],
       [ 1.09 ,  1.09 ,  1.09 ,  1.09 ,  1.09 ]])

Let's take a look at bed load transport rate when we use the different bedload equations
For the defaul MPM model we get:

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="MPM",
...     bed_surface__gsd_location_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transport__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-5.4212e-04

For Fernandez Luque and Van Beek:

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="FLvB",
...     bed_surface__gsd_location_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transport__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-1.3681e-04

For Wilcock and Crowe 2003:

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WilcockAndCrowe",
...     bed_surface__gsd_location_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transport__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-8.3608e-06

The previous example, covers a relatively complete case. For demonstration purposes
let's see some other options that show how to use or change the default setting.
If the grain size distribution is not specified, what value will river bed dynamics use?

>>> rbd = RiverBedDynamics(grid)
>>> rbd._gsd
array([[ 32, 100],
       [ 16,  25],
       [  8,   0]])

The sand content can be calculated from a grain size distribution

>>> gsd = [[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [2, 10], [1, 0]]
>>> rbd = RiverBedDynamics(grid, gsd=gsd, bedload_equation="MPM")
>>> sandContent = rbd._bed_surface__sand_fraction_node[20]
>>> print(float("{:.4f}".format(round(sandContent, 4))))
0.1

But it is different if we use Parker 1990, because it removes sand content

>>> rbd = RiverBedDynamics(grid, gsd=gsd, bedload_equation="Parker1990")
>>> rbd._bed_surface__sand_fraction_node[20]
0.0

What happens if we give it a set of wrong optional fields. The following
fields will have only two elements, which is different than the number of
nodes and links

>>> qb_imposed = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, sed_transport__bedload_rate_imposed_link=qb_imposed)
>>> rbd._sed_transport__bedload_rate_imposed_link
array([ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,
        0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.])

>>> gsd_loc = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surface__gsd_location_node=gsd_loc)
>>> rbd._bed_surface__gsd_location_node.reshape(grid.shape)
array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.]])

>>> gsd_fix = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surface__gsd_fixed_node=gsd_fix)
>>> rbd._bed_surface__gsd_fixed_node.reshape(grid.shape)
array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.]])

>>> elev_fix = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surface__elevation_fixed_node=elev_fix)
>>> rbd._bed_surface__elevation_fixed_node.reshape(grid.shape)
array([[ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.]])

>>> vel_n_1 = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, surface_water__velocity_prev_time_link=vel_n_1)
>>> rbd._surface_water__velocity_prev_time_link
array([ 0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,  0.25,
        0.25,  0.25,  0.25,  0.25])

For sed_transport__bedload_gsd_imposed_link, for simplicity, we only show
links 0, 5, and 10

>>> qb_gsd_imposed = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, sed_transport__bedload_gsd_imposed_link=qb_gsd_imposed)
>>> rbd._sed_transport__bedload_gsd_imposed_link[[0, 5, 10], :]
array([[ 0.,  0.],
       [ 0.,  0.],
       [ 0.,  0.]])

In summary, in all these cases the wrong given value is override by default values.
But, if the size of the array is correct the specified condition is used.
For simplicity, we only show links 0, 5, and 10

>>> qb_imposed = np.full(grid.number_of_links, 1)
>>> rbd = RiverBedDynamics(grid, sed_transport__bedload_rate_imposed_link=qb_imposed)
>>> rbd._sed_transport__bedload_rate_imposed_link[[0, 5, 10]]
array([1, 1, 1])

>>> gsd_loc = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surface__gsd_location_node=gsd_loc)
>>> rbd._bed_surface__gsd_location_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> gsd_fix = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surface__gsd_fixed_node=gsd_fix)
>>> rbd._bed_surface__gsd_fixed_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> elev_fix = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surface__elevation_fixed_node=elev_fix)
>>> rbd._bed_surface__elevation_fixed_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> vel_n_1 = np.full(grid.number_of_links, 1)
>>> rbd = RiverBedDynamics(grid, surface_water__velocity_prev_time_link=vel_n_1)
>>> rbd._surface_water__velocity_prev_time_link[[0, 5, 10]]
array([1, 1, 1])

>>> qb_gsd_imposed = np.ones((grid.number_of_links, 2))  # 2 comes from gsd.shape[0]-1
>>> rbd = RiverBedDynamics(grid, sed_transport__bedload_gsd_imposed_link=qb_gsd_imposed)
>>> rbd._sed_transport__bedload_gsd_imposed_link[[0, 5, 10]]
array([[ 1.,  1.],
       [ 1.,  1.],
       [ 1.,  1.]])

Using the hydraulics radius is also possible. Let's compare the shear stress with and
without that option

>>> rbd = RiverBedDynamics(grid)
>>> rbd.run_one_step()
>>> shearStress = rbd._shear_stress[15]
>>> np.round(shearStress, decimals=3)
-15.532

>>> rbd = RiverBedDynamics(grid, use_hydraulics_radius_in_shear_stress=True)
>>> rbd.run_one_step()
>>> shearStress = rbd._shear_stress[15]
>>> print(float("{:.3f}".format(round(shearStress, 3))))
-12.892

So, there is an important difference between the two ways of calculating it.

For a complete list of all the fields created during the execution of river bed dynamics
use:

>>> fields = RiverBedDynamics.get_available_fields()
"""
import copy
import os
import shutil
import time

import numpy as np
import pandas as pd
import scipy.constants
from scipy.interpolate import interp1d

from landlab import Component


class RiverBedDynamics(Component):

    """Predicts the evolution of a river bed.

    Landlab component that predicts the evolution of a river bed
    considering changes in elevation and grain size distribution in response to
    bed load transport according to the Exner equation and the transfer
    function of Toro-Ecobar et al., (1996).

    To estimate temporal and spatial changes in river bed properties, this
    component predicts the bedload transport rate and fractional transport at
    each link using unsteady shear stress. Time-varying hydraulic
    variables are obtained from a surface flow, for instance,
    :class:`~OverlandFlow`.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Not required but recommended

    Adams, J., Gasparini, N., Hobley, D., Tucker, G., Hutton, E., Nudurupati,
    S., Istanbulluoglu, E. (2017). The Landlab v1. 0 OverlandFlow component:
    a Python tool for computing shallow-water flow across watersheds.
    Geoscientific Model Development  10(4), 1645.
    https://dx.doi.org/10.5194/gmd-10-1645-2017

    **Additional References**

    G. Parker (1990) Surface-based bedload transport relation for gravel
    rivers, Journal of Hydraulic Research, 28:4, 417-436,
    DOI: 10.1080/00221689009499058

    Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model for
    mixed-size sediment. Journal of hydraulic engineering, 129(2), 120-128.
    DOI: 10.1061/(ASCE)0733-9429(2003)129:2(120)

    Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport,
    Proceedings, 2nd Congress, International Association of Hydraulic Research,
    Stockholm: 39-64.

    Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of
    bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
    https://doi.org/10.1080/00221687609499677

    Mueller, E. R., J. Pitlick, and J. M. Nelson (2005), Variation in the
    reference Shields stress for bed load transport in gravelbed streams and
    rivers, Water Resour. Res., 41, W04006, doi:10.1029/2004WR003692

    Carlos M. Toro-Escobar, Chris Paola & Gary Parker (1996) Transfer function
    for the deposition of poorly sorted gravel in response to streambed
    aggradation, Journal of Hydraulic Research, 34:1, 35-53,
    DOI: 10.1080/00221689609498763
    """

    _name = "RiverBedDynamics"

    _unit_agnostic = False

    _info = {
        "surface_water__depth": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "link",
            "doc": "Depth of water on the surface",
        },
        "surface_water__velocity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "Speed of water flow above the surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        gsd=None,
        rho=1000,
        rho_s=2650,
        bedload_equation="MPM",
        variable_critical_shear_stress=False,
        use_hydraulics_radius_in_shear_stress=False,
        lambda_p=0.35,
        outlet_boundary_condition="zeroGradient",
        evolve_bed=True,
        update_bed_surface_GSD=True,
        track_stratigraphy=False,
        number_cycles_to_process_stratigraphy=10,
        new_surface_layer_thickness=1,
        output_vector=False,
        dt=1,
        alpha=1.0,
        bed_surface__elevation_fixed_node=None,
        sed_transport__bedload_rate_imposed_link=None,
        sed_transport__bedload_gsd_imposed_link=None,
        bed_surface__gsd_fixed_node=None,
        bed_surface__gsd_location_node=None,
        surface_water__velocity_prev_time_link=None,
        current_t=0.0,
    ):
        """Calculates the evolution of a river bed based on bed load transport
        and fractional rates on links. An external flow hydraulics solver, such
        as OverlandFlow, is required to predict flow variables in time and space.
        The shear stress used in sediment transport equations takes into account
        the time and spatial variability of water depth and flow velocity.

        This component adjusts topographic elevation and grain size
        distribution at each node within a Landlab grid.

        Parameters
        ----------
        grid : RasterModelGrid
            area Landlab grid.
        gsd : ndarray of float
            Grain size distribution. Must contain as many GSDs as there are
            different indexes in GSD Location.
        rho : float, optional
            Fluid density. Defaults to the density of water at 1000 kg/m^3.
        rho_s : float, optional
            Sediment density. Defaults to sediment density of 2650 kg/m^3.
        bedload_equation : String, optional
            Selects the bedload transport equation. Options include:

            * ``'MPM'`` for Meyer Peter and Muller,
            * ``'FLvB'`` for Fernandez Luque & van Beek (1976),
            * ``'Parker1990'`` for Parker 1990,
            * ``'WilcockAndCrowe'`` for Wilcock and Crowe 2003.
        variable_critical_shear_stress: bool, optional
            If ``True``, the critical shear stress in ``'MPM'`` or ``'FLvB'`` will be
            obtained using Mueller et al. (2005) equation.
        use_hydraulics_radius_in_shear_stress: bool, optional
            If ``True``, shear stress calculations will use the hydraulic radius.
            Defaults to False, which uses water depth.
        lambda_p : float, optional
            Sediment porosity. Default value is 0.35.
        outlet_boundary_condition : str, optional
            Sets the boundary condition at the watershed outlet. Options are:
            ``'zeroGradient'`` (maps the outlet to the upstream node, default) or
            ``'fixedValue'`` (does not change the value at the outlet during the
            run).
        evolve_bed : bool, optional
            If ``True``, the bed evolves according to the local bed load transport
            rate and GSD. If ``False``, the bed load transport rate and GSD are
            calculated and output.
        update_bed_surface_GSD : bool, optional
            If ``True``, the bed GSD evolves according to the interaction between
            the current surface GSD, and the bed load GSD and rate.
        track_stratigraphy : bool, optional
            If ``True``, the component stores the GSD of each layer at every node.
            This is computationally demanding as it needs to read/write data
            at every time step. Recommended only when cycles of erosion and
            deposition are frequent and very important.
        number_cycles_to_process_stratigraphy : int, optional
            If ``track_stratigraphy`` is ``True``, data will be read and stored every
            ``'number_cycles_to_process_stratigraphy'`` time steps. Must be larger
            or equal to 1. All data are written to the hard drive. It does not use
            *Landlab* layers.
        dt: float, optional
            Time step in seconds. When this component is coupled to a flow model,
            it is dynamically updated.
        alpha : float, optional
            An upwinding coefficient for a central difference scheme when
            updating the bed GSD - default value is 1.0 - a value of 0.5
            generates a central differences scheme.
        bed_surface__elevation_fixed_node: ndarray of int, optional
            Sets a node as a fixed elevation.
            Units: - , mapping: node
        sed_transport__bedload_rate_imposed_link: ndarray of float, optional
            Sets the sediment transport rate per unit width where sediment supply is imposed.
            Units: m^2/s, mapping: link
        bed_surface__gsd_fixed_node: ndarray of int, optional
            Sets the locations (nodes) that do not change its GSD in time
            Units: - , mapping: node
        bed_surface__gsd_location_node: ndarray of int, optional
            Sets the location at each node in which the GSD applies
            Units: - , mapping: node
        sed_transport__bedload_gsd_imposed_link: ndarray of float, optional
            Sets the sediment transport GSD where sediment supply is imposed. It's
            size should be columns: number of links rows: Number of grain sizes - 1
            Units: - , mapping: link
        surface_water__velocity_prev_time_link: ndarray of float, optional
            Speed of water flow above the surface in the previous time step
            Units: m/s, mapping: link
        current_t : float, optional
            Current simulation time or elapsed time. It does not update automatically
            Units: s
        """
        super().__init__(grid)

        self._g = scipy.constants.g  # Acceleration due to gravity (m/s**2).
        self._rho = rho
        self._rho_s = rho_s
        self._R = (rho_s - rho) / rho
        self._shear_stress_star_rsgo = 0.0386  # Reference dimensionless shear
        # stress for the median size
        self._beta = 0.0951  # Coefficient for the hiding function

        if gsd is None:
            self._gsd = np.array([[32, 100], [16, 25], [8, 0]])
        else:
            self._gsd = np.array(gsd)

        # self._bedload_equation = bedload_equation

        if bedload_equation == "Parker1990":
            self._bedload_equation = self.bedload_equation_Parker_1990
        elif bedload_equation == "MPM":
            self._bedload_equation = self.bedload_equation_MeyerPeter_Muller
        elif bedload_equation == "FLvB":
            self._bedload_equation = self.bedload_equation_FernandezLuque_VanBeek
        elif bedload_equation == "WilcockAndCrowe":
            self._bedload_equation = self.bedload_equation_Wilcock_Crowe_2003

        self._variable_critical_shear_stress = variable_critical_shear_stress
        self._use_hydraulics_radius_in_shear_stress = (
            use_hydraulics_radius_in_shear_stress
        )
        self._lambda_p = lambda_p
        self._alpha = alpha
        self._outlet_boundary_condition = outlet_boundary_condition

        self._grid._dt = dt
        self._current_t = current_t

        # True if it is the first iteration
        self._first_iteration = True

        # Initialize the field such that exist before the initial bed properties
        # function is called

        if bed_surface__gsd_location_node is None:
            self._bed_surface__gsd_location_node = self._grid.zeros(at="node")
        else:
            bed_surface__gsd_location_node = np.array(
                bed_surface__gsd_location_node
            ).flatten()
            if bed_surface__gsd_location_node.size > 0:
                if bed_surface__gsd_location_node.size == self._grid.number_of_nodes:
                    self._bed_surface__gsd_location_node = np.array(
                        bed_surface__gsd_location_node
                    ).flatten()
                else:
                    self._bed_surface__gsd_location_node = np.zeros(
                        self._grid.number_of_nodes
                    )

        # Initialize the bed surface grain size properties. At the beginning,
        # and only for the first timestep, the Grain Size Distribution (GSD)
        # from the input file is used. For subsequent timesteps, the
        # previously calculated GSD is used.
        self.define_initial_bed_properties()

        # Finds the Id of the links that are next to a border cell
        self.links_at_border_cells()

        # This flag changes to True if there are links with fixed bedload GSD
        self._fixed_link_calculate = False

        # This flag is used to activate or deactivate the bed evolution part
        # of the component.
        self._evolve_bed = evolve_bed

        # This flag is used to activate or deactivate the bed GSD updating part
        # of the component.
        self._update_bed_surface_GSD = update_bed_surface_GSD

        # This flag is used to activate or deactivate option to store the GSD
        # of individual layers in each node.
        self._track_stratigraphy = track_stratigraphy

        # This value is used to merge all deposited layers in a new subsurface
        # layer
        self._new_surface_layer_thickness = new_surface_layer_thickness

        # When new_surface_layer_thickness is reached this flag is used to
        # record and read the data
        self._compute_stratigraphy = False

        # When new_surface_layer_thickness is reached this flag is used to
        # update subsurface data
        self._update_subsurface = False
        self._update_subsurface_eroded = False
        self._subsurface_changed = False

        # If true exports shear stress and velocity as vectors on nodes at each
        # time step
        self._output_vector = output_vector

        # Sets initial values to enter into the write and read stratigraphy
        # routine
        self._stratigraphy_cycle = 0
        self._number_cycles_to_process_stratigraphy = (
            number_cycles_to_process_stratigraphy
        )

        # Creating optional grid fields at time zero - If these fields were not
        # defined before instantiation it will create it and fill all values
        # with zeros. Velocity at previous time simply copies the input velocity.
        # An error will be raised if the defined flied size is not correct

        if sed_transport__bedload_rate_imposed_link is None:
            self._sed_transport__bedload_rate_imposed_link = np.zeros(
                self._grid.number_of_links
            )
        else:
            sed_transport__bedload_rate_imposed_link = np.array(
                sed_transport__bedload_rate_imposed_link
            ).flatten()
            if sed_transport__bedload_rate_imposed_link.size > 0:
                if (
                    sed_transport__bedload_rate_imposed_link.size
                    == self._grid.number_of_links
                ):
                    self._sed_transport__bedload_rate_imposed_link = (
                        sed_transport__bedload_rate_imposed_link
                    ).flatten()
                else:
                    self._sed_transport__bedload_rate_imposed_link = np.zeros(
                        self._grid.number_of_links
                    )

        if sed_transport__bedload_gsd_imposed_link is None:
            self._sed_transport__bedload_gsd_imposed_link = np.zeros(
                (self._grid.number_of_links, self._gsd.shape[0] - 1)
            )
        else:
            sed_transport__bedload_gsd_imposed_link = np.array(
                sed_transport__bedload_gsd_imposed_link
            )

            if sed_transport__bedload_gsd_imposed_link.shape[0] > 0:
                if (
                    sed_transport__bedload_gsd_imposed_link.shape[0]
                    == self._grid.number_of_links
                ) and (
                    sed_transport__bedload_gsd_imposed_link.shape[1]
                    == self._gsd.shape[0] - 1
                ):
                    self._sed_transport__bedload_gsd_imposed_link = np.array(
                        sed_transport__bedload_gsd_imposed_link
                    )
                else:
                    self._sed_transport__bedload_gsd_imposed_link = np.zeros(
                        (self._grid.number_of_links, self._gsd.shape[0] - 1)
                    )

        if bed_surface__gsd_fixed_node is None:
            self._bed_surface__gsd_fixed_node = self._grid.zeros(at="node")
        else:
            bed_surface__gsd_fixed_node = np.array(
                bed_surface__gsd_fixed_node
            ).flatten()
            if bed_surface__gsd_fixed_node.size > 0:
                if bed_surface__gsd_fixed_node.size == self._grid.number_of_nodes:
                    self._bed_surface__gsd_fixed_node = np.array(
                        bed_surface__gsd_fixed_node
                    ).flatten()
                else:
                    self._bed_surface__gsd_fixed_node = np.zeros(
                        self._grid.number_of_nodes
                    )

        if bed_surface__elevation_fixed_node is None:
            self._bed_surface__elevation_fixed_node = np.zeros(
                self._grid.number_of_nodes
            )
        else:
            bed_surface__elevation_fixed_node = np.array(
                bed_surface__elevation_fixed_node
            ).flatten()
            if bed_surface__elevation_fixed_node.size > 0:
                if bed_surface__elevation_fixed_node.size == self._grid.number_of_nodes:
                    self._bed_surface__elevation_fixed_node = (
                        bed_surface__elevation_fixed_node
                    ).flatten()
                else:
                    self._bed_surface__elevation_fixed_node = np.zeros(
                        self._grid.number_of_nodes
                    )

        if surface_water__velocity_prev_time_link is None:
            self._surface_water__velocity_prev_time_link = copy.deepcopy(
                self._grid["link"]["surface_water__velocity"]
            )
        else:
            surface_water__velocity_prev_time_link = np.array(
                surface_water__velocity_prev_time_link
            ).flatten()
            if surface_water__velocity_prev_time_link.size > 0:
                if (
                    surface_water__velocity_prev_time_link.size
                    == self._grid.number_of_links
                ):
                    self._surface_water__velocity_prev_time_link = (
                        surface_water__velocity_prev_time_link
                    ).flatten()
                else:
                    self._surface_water__velocity_prev_time_link = copy.deepcopy(
                        self._grid["link"]["surface_water__velocity"]
                    )

        # Initializating fields at time zero
        # Volumetric bed load transport rate per unit width
        self._sed_transport__bedload_rate_link = self._grid.zeros(at="link")
        # Net sediment transport rate at a node
        self._sed_transport__net_bedload_node = self._grid.zeros(at="node")
        # bed load grain size distribution
        self._sed_transport__bedload_gsd_link = self._grid.zeros(at="link")
        # Shear stress applied by the surface water on the bed surface - Pa
        self._surface_water__shear_stress_link = self._grid.zeros(at="link")

        # Define faces normal vector
        self._normal = -self._grid.link_dirs_at_node

        # Makes a copy of the original bed surface elevation and maps into links
        self._topographic__elevation_original_node = self._grid["node"][
            "topographic__elevation"
        ].copy()
        self._topographic__elevation_original_link = (
            self._grid.map_mean_of_link_nodes_to_link(
                self._topographic__elevation_original_node
            )
        )

        self._topographic__elevation_subsurface_link = (
            self._grid.map_mean_of_link_nodes_to_link(
                self._topographic__elevation_original_node
            )
        )

        # Initialize the field
        self._bed_surface__surface_thickness_new_layer_link = np.zeros_like(
            self._topographic__elevation_original_node
        )

        # Defines some variables for storing stratigraphy
        self._x = 0.5 * (
            self._grid.node_x[self._grid.node_at_link_head]
            + self._grid.node_x[self._grid.node_at_link_tail]
        )
        self._y = 0.5 * (
            self._grid.node_y[self._grid.node_at_link_head]
            + self._grid.node_y[self._grid.node_at_link_tail]
        )

        # folders location
        self._cwd = os.getcwd()
        self._stratigraphy_temp_files_path = "stratigraphyTempFiles"
        self._stratigraphy_raw_data_path = "stratigraphyRawData"

    def define_initial_bed_properties(self):
        """This method performs the initial setup of the bed grain size
        distribution properties. It reads the input data and populates the
        necessary variables.

        This configuration is only performed during the first time step. Subsequent
        time steps will utilize the bed information, which has already
        been calculated or updated and formatted appropriately.
        """

        grain_size_D = self._gsd[:, 0]  # Grain sizes
        grain_size_frequency = self._gsd[:, 1:]  # Grain sizes frequency cumulative

        # Number of locations with different GSD during the first time step
        number_gsd_locations = self._gsd.shape[1] - 1

        # Calculates the sand fraction in each location
        if np.min(grain_size_D) < 2:  # Grain smaller than 2 mm are sand
            sand_fraction_0 = self.calculate_sand_fraction()
            # Flag is set as True to adjust the GSD and calculate a sand-free
            # GSD when using Parker 1990 eq
            adjust_gsd_flag = True
        else:
            sand_fraction_0 = np.zeros(self._gsd.shape[1] - 1)
            # Flag is set as False, there is no need to adjust the GSD
            adjust_gsd_flag = False

        # If Parker Eq is used we remove sand content
        if adjust_gsd_flag:
            # We add 2mm into GSD and update the GSD
            id2mm = np.argmin(grain_size_D >= 2)  # Location where 2 mm will be placed

            if self._bedload_equation == self.bedload_equation_Parker_1990:
                # Adds 2 mm and removes sand and smaller grains
                grain_size_D = np.concatenate([grain_size_D[0:id2mm], [2]])
                grain_size_frequency = np.concatenate(
                    [grain_size_frequency[0:id2mm, :], [sand_fraction_0]]
                )
                gravel_fraction = (
                    grain_size_frequency[0, :] - grain_size_frequency[-1, :]
                ) / 100
                # Calculates cumulative grain sizes frequency
                grain_size_frequency0 = np.flip(
                    -np.diff(grain_size_frequency / 100, axis=0) / gravel_fraction,
                    axis=0,
                )
                grain_size_frequency = np.abs(
                    np.flip(np.cumsum(grain_size_frequency0, axis=0), axis=0) * 100
                )
                grain_size_frequency = np.vstack(
                    (grain_size_frequency, np.zeros([1, grain_size_frequency.shape[1]]))
                )
                sand_fraction_0 = sand_fraction_0 * 0
            else:
                # Adds 2 mm but does not remove sand and smaller grains
                grain_size_D = np.concatenate(
                    [grain_size_D[0:id2mm], [2], grain_size_D[id2mm:]]
                )
                grain_size_frequency = np.concatenate(
                    [
                        grain_size_frequency[0:id2mm, :],
                        [sand_fraction_0],
                        grain_size_frequency[id2mm:, :],
                    ]
                )

        # Copies the grain sizes after removing sand if Parker Eq. is used.
        # All grains will be copies if other equation is used.
        self._grain_size_D_original = grain_size_D

        # Grain sizes frequency - Now is not cumulative anymore and has the
        # same dimensions as the equivalent grain_size_D or grain_size_D_equivalent
        grain_size_frequency = np.abs(-np.diff(grain_size_frequency / 100, axis=0))
        # Equivalent grain sizes
        grain_size_D_equivalent = (grain_size_D[0:-1] * grain_size_D[1:]) ** 0.5
        sand_fraction = np.zeros_like(
            self._bed_surface__gsd_location_node
        )  # Sand fraction at each node

        # Bed grain sizes frequency in each node
        grain_size_D_equiv_freq = np.zeros(
            [self._grid.number_of_nodes, grain_size_D.shape[0] - 1]
        )
        for i in range(number_gsd_locations):
            (id_gsdLocation,) = np.where(self._bed_surface__gsd_location_node == i)
            sand_fraction[id_gsdLocation] = sand_fraction_0[i] / 100
            grain_size_D_equiv_freq[id_gsdLocation, :] = grain_size_frequency[:, i]

        # Calculating median_size_D50 at each node - Grain sizes in Psi scale
        grain_size_Psi_scale_D = np.flip(np.log2(grain_size_D), axis=0)
        # Frequency of each grain size in each node
        grain_size_D_equiv_freq_cum = np.hstack(
            (
                np.zeros(
                    [
                        self._bed_surface__gsd_location_node.shape[0],
                        1,
                    ]
                ),
                np.cumsum(np.flip(grain_size_D_equiv_freq, axis=1), axis=1),
            )
        )

        # Finds the index of the grain size smaller than 50% in each node
        i0 = np.argmin(grain_size_D_equiv_freq_cum <= 0.5, axis=1) - 1
        # Finds the index of the grain size larger than 50% in each node
        i1 = np.argmax(grain_size_D_equiv_freq_cum > 0.5, axis=1)
        nodes_list = np.arange(0, self._grid.number_of_nodes)

        grain_size_Psi_scale_D50 = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equiv_freq_cum[nodes_list, i1]
                - grain_size_D_equiv_freq_cum[nodes_list, i0]
            )
        ) * (0.5 - grain_size_D_equiv_freq_cum[nodes_list, i0])
        median_size_D50 = (
            2**grain_size_Psi_scale_D50
        )  # Median grain size in each node

        # Calculating the geometric mean and standard deviation
        # Equivalent grain sizes in Psi scale
        grain_size_D_equivalent_Psi = np.log2(grain_size_D_equivalent)
        grain_size_D_equivalent_Psi_mean = np.sum(
            grain_size_D_equiv_freq * grain_size_D_equivalent_Psi, axis=1
        )
        grain_size_geom_mean = 2**grain_size_D_equivalent_Psi_mean
        # Equivalent grain sizes in Psi scale at each node
        grain_size_D_equivalent_Psi = np.tile(
            grain_size_D_equivalent_Psi, (self._grid.number_of_nodes, 1)
        )
        grain_size_D_equivalent_Psi_mean = np.reshape(
            grain_size_D_equivalent_Psi_mean, [self._grid.number_of_nodes, 1]
        )
        grain_size_geo_std = 2 ** np.sqrt(
            np.sum(
                ((grain_size_D_equivalent_Psi - grain_size_D_equivalent_Psi_mean) ** 2)
                * grain_size_D_equiv_freq,
                axis=1,
            )
        )

        # We save all data into the grid so it can be shared with different
        # functions and components
        # Bed grain sizes frequency in each node
        self._bed_surface__gsd_node = grain_size_D_equiv_freq
        self._bed_surface__gsd_original_node = grain_size_D_equiv_freq.copy()
        self._bed_subsurface__gsd_node = copy.deepcopy(grain_size_D_equiv_freq)
        self._bed_surface__median_size_node = median_size_D50
        self._bed_surface__geom_mean_size_node = grain_size_geom_mean
        self._bed_surface__geo_std_size_node = grain_size_geo_std
        self._bed_surface__sand_fraction_node = sand_fraction

        # GSD properties is mapped from nodes onto links
        self._bed_surface__gsd_link = self.map_mean_of_nodes_to_link(
            grain_size_D_equiv_freq
        )
        self._bed_surface__median_size_link = self.map_mean_of_nodes_to_link(
            median_size_D50
        )
        self._bed_surface__geom_mean_size_link = self.map_mean_of_nodes_to_link(
            grain_size_geom_mean
        )
        self._bed_surface__geo_std_size_link = self.map_mean_of_nodes_to_link(
            grain_size_geo_std
        )
        self._bed_surface__sand_fraction_link = self.map_mean_of_nodes_to_link(
            sand_fraction
        )

        self._bed_surface__gsd_original_link = copy.deepcopy(
            self._bed_surface__gsd_link
        )
        self._bed_subsurface__gsd_link = copy.deepcopy(self._bed_surface__gsd_link)

        self.calculate_act_layer_thick()

        # Sets the initial active layer thickness as zero meter deep.
        self._bed_surface__act_layer_thick_prev_time_node = np.zeros(
            self._grid.number_of_nodes
        )
        self._bed_surface__act_layer_thick_prev_time_link = np.zeros(
            self._grid.number_of_links
        )

        # Every equation will work with equivalent grain sizes.
        # We've dropped the 'Eq' in this definition.
        self._grain_size = grain_size_D_equivalent

    def links_at_border_cells(self):
        self._links_at_border_cells = np.sort(
            np.hstack(
                (
                    self._grid.links_at_node[:, 2][self._grid.nodes_at_right_edge],
                    self._grid.links_at_node[:, 3][self._grid.nodes_at_top_edge],
                    self._grid.links_at_node[:, 0][self._grid.nodes_at_left_edge],
                    self._grid.links_at_node[:, 1][self._grid.nodes_at_bottom_edge],
                )
            )
        )

    def run_one_step(self):
        """The component can be divided into two parts. In the first part, all bed
        load transport and GSD calculations, including shear stress estimates,
        are conducted. In the second part, bed GSD and bed elevation can evolve.

        **First part**

        Calculates shear stress and bed load transport rates across the grid.

        For one time step, this generates the shear stress across a given
        grid by accounting for the local water depth, bed surface slope, and
        water velocity at links. Then, based on this shear stress and the local
        bed surface grain size distribution, the bed load transport rate is
        calculated at each link and mapped onto each node. Bed load grain size
        distributions are calculated when using Parker's 1990 or Wilcock and
        Crowe's 2003 equations. Meyer-Peter and Muller and Fernandez Luque and
        van Beek models will only calculate the total bed load transport.

        Outputs the following bed surface properties over time at every link
        in the input grid: geometric mean size, GSD, median size, sand fraction,
        standard deviation size.

        Also outputs the bed load GSD, bed load rate, and shear stress over
        time at every link in the input grid. The net bed load is output
        over time at every node.

        **Second Part**

        Erode grid topography. Starts at self._evolve_bed

        For one time step, this erodes the grid topography according to
        Exner equation::

            (1-λp) ∂Z/∂t = - (∂qbx/∂x + ∂qby/∂y)

        Simplifying, we get::

            ∂Z/∂t = - (1 / (1-λp)) * (∂Qb/∂A)
            Z_t+1 = -(Δt * ΔQb)/(1-λp) + Z_t

        The grid field ``"topographic__elevation"`` is altered each time step.
        """

        if self._first_iteration:
            # Identify nodes with fixed elevations and surface GSD
            self.fixed_nodes_info()
            # Identify links with fixed bedload GSD
            self.fixed_links_info()
            # Identify the node upstream of the outlet
            self.outlet_nodes_info()
            self.Parker_1990_interp_points()

        # Unsteady shear stress is calculated every time step
        self.shear_stress()

        self._bedload_equation()

        # Calculates the net bedload transport from links into nodes
        self.calculate_net_bedload()

        # Maps bed load vectors from links (m2/s) onto nodes (m3/s) preserving
        # vector components.
        if self._output_vector:
            (
                self._shear_stress_vector,
                self._shear_stress_magnitude,
            ) = RiverBedDynamics.vector_mapper(self.grid, self._shear_stress)
            self._u_vector, self._u_magnitude = RiverBedDynamics.vector_mapper(
                self.grid, self._shear_stress
            )
            self.bedload_mapper()

        if self._evolve_bed:
            # This routine is entered only after the second time step to avoid
            # re-calculation.
            # When the bed has undergone changes, new bed properties are calculated.
            # This updated information will then be used when calculating bedload rates.

            # Update the bed elevation
            self.update_bed_elevation()

            # Update the bed grain size distribution
            if (self._bedload_equation == self.bedload_equation_Parker_1990) or (
                self._bedload_equation == self.bedload_equation_Wilcock_Crowe_2003
            ):
                if self._update_bed_surface_GSD is True:
                    self.update_bed_surface_gsd()
                    self.update_bed_surface_properties()
                self.map_gsd_from_link_to_node()

        self._first_iteration = False
        self._stratigraphy_cycle += 1

    def fixed_nodes_info(self):
        """Search and identify nodes defined as having a fixed elevation"""

        # Gives the ID of the outlet node -  First we make sure that we use only integers
        self._bed_surface__elevation_fixed_node = np.round(
            self._bed_surface__elevation_fixed_node
        ).astype(int)
        self._bed_surface__gsd_fixed_node = np.round(
            self._bed_surface__gsd_fixed_node
        ).astype(int)
        (self._fixed_nodes_id,) = np.where(self._bed_surface__elevation_fixed_node == 1)
        (self._fixed_surface_gsd_nodes_id,) = np.where(
            self._bed_surface__gsd_fixed_node == 1
        )

        # All connecting links to these nodes will be set as fixed too
        self._fixed_links = np.unique(
            self._grid.links_at_node[self._fixed_surface_gsd_nodes_id]
        )
        self._fixed_links = self._fixed_links[np.where(self._fixed_links >= 0)]

    def fixed_links_info(self):
        """Search and identify links defined as having a fixed bed load GSD
        This function is used only when the bedload gsd in imposed on a link.
        For testing purposes, the a case where this function is called is included
        below. It is based on the main example.

        Examples
        --------

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
        >>> grid.at_node["surface_water__velocity"] = np.full(
        ...     grid.number_of_nodes, 0.25
        ... )
        >>> grid.at_link["surface_water__depth"] = map_mean_of_link_nodes_to_link(
        ...     grid, "surface_water__depth"
        ... )
        >>> grid.at_link["surface_water__velocity"] = map_mean_of_link_nodes_to_link(
        ...     grid, "surface_water__velocity"
        ... )
        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.07, 1.06, 1.00, 1.06, 1.07],
        ...     [1.08, 1.07, 1.03, 1.07, 1.08],
        ...     [1.09, 1.08, 1.07, 1.08, 1.09],
        ...     [1.09, 1.09, 1.08, 1.09, 1.09],
        ...     [1.09, 1.09, 1.09, 1.09, 1.09],
        ... ]
        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> gsd_loc = [
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ... ]
        >>> gsd = [[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [2, 10], [1, 0]]
        >>> qb_imposed_gsd = np.zeros(
        ...     (grid.number_of_links, np.array(gsd).shape[0] - 1)
        ... )
        >>> qb_imposed_gsd[29, :] = np.array([0.15, 0.15, 0.2, 0.2, 0.15, 0.15])
        >>> rbd = RiverBedDynamics(
        ...     grid,
        ...     gsd=gsd,
        ...     bedload_equation="Parker1990",
        ...     bed_surface__gsd_location_node=gsd_loc,
        ...     output_vector=True,
        ...     track_stratigraphy=True,
        ...     sed_transport__bedload_gsd_imposed_link=qb_imposed_gsd,
        ... )
        >>> rbd.run_one_step()

        Let's check which link is has an imposed bedload gsd

        >>> rbd._fixed_surface_gsd_link_id[0]
        29
        """
        # Gives the ID of the links
        if np.max(self._sed_transport__bedload_gsd_imposed_link > 0.0):
            (a0, a1) = np.where(self._sed_transport__bedload_gsd_imposed_link > 0.0)
            self._fixed_surface_gsd_link_id = a0
            self._fixed_surface_gsd_link_fractions = a1

            # We use this flag to avoid extra calculations when there is no fixed
            # bed load GSD
            if a0.shape[0] > 0:
                self._fixed_link_calculate = True

    def calculate_sand_fraction(self):
        grain_size_Psi_scale_D = np.flip(np.log2(self._gsd[:, 0]), axis=0)
        grain_size_D_equiv_freq_cum = np.flip(self._gsd[:, 1:], axis=0)
        sand_fraction = np.zeros(grain_size_D_equiv_freq_cum.shape[1])
        i = np.max(np.where(grain_size_Psi_scale_D <= 1))
        for j in np.arange(0, grain_size_D_equiv_freq_cum.shape[1]):
            sand_fraction[j] = (
                (1 - grain_size_Psi_scale_D[i])
                / (grain_size_Psi_scale_D[i + 1] - grain_size_Psi_scale_D[i])
            ) * (
                grain_size_D_equiv_freq_cum[i + 1, j]
                - grain_size_D_equiv_freq_cum[i, j]
            ) + grain_size_D_equiv_freq_cum[
                i, j
            ]

        # self._bed_surface__sand_fraction_node = sand_fraction
        return sand_fraction

    def map_mean_of_nodes_to_link(self, grain_size_variable):
        return 0.5 * (
            grain_size_variable[self.grid.nodes_at_link[:, 0]]
            + grain_size_variable[self.grid.nodes_at_link[:, 1]]
        )

    def calculate_act_layer_thick(self):
        """First for nodes"""

        grain_size_D = self._grain_size_D_original  # Grain sizes
        # Grain sizes in Psi scale
        grain_size_Psi_scale_D = np.flip(np.log2(grain_size_D), axis=0)

        grain_size_D_equiv_freq = self._bed_surface__gsd_node
        frequency_grain_size_list = np.arange(self._grid.number_of_nodes)

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equiv_freq,
                np.zeros([grain_size_D_equiv_freq.shape[0], 1]),
            )
        )

        grain_size_D_equiv_freq_cum = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equiv_freq_cum <= 0.9, axis=1) - 1
        i1 = np.argmax(grain_size_D_equiv_freq_cum > 0.9, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equiv_freq_cum[frequency_grain_size_list, i1]
                - grain_size_D_equiv_freq_cum[frequency_grain_size_list, i0]
            )
        ) * (0.9 - grain_size_D_equiv_freq_cum[frequency_grain_size_list, i0])
        surface_D90 = 2**grain_size_Psi_scale_DX
        self._bed_surface__act_layer_thick_node = 2 * surface_D90 / 1000

        # Now for links

        grain_size_D_equiv_freq = self._bed_surface__gsd_link
        frequency_grain_size_list = np.arange(self._grid.number_of_links)

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equiv_freq,
                np.zeros([grain_size_D_equiv_freq.shape[0], 1]),
            )
        )

        grain_size_D_equiv_freq_cum = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equiv_freq_cum <= 0.9, axis=1) - 1
        i1 = np.argmax(grain_size_D_equiv_freq_cum > 0.9, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equiv_freq_cum[frequency_grain_size_list, i1]
                - grain_size_D_equiv_freq_cum[frequency_grain_size_list, i0]
            )
        ) * (0.9 - grain_size_D_equiv_freq_cum[frequency_grain_size_list, i0])
        surface_D90 = 2**grain_size_Psi_scale_DX  # 90th percentile
        self._bed_surface__act_layer_thick_link = 2 * surface_D90 / 1000

    def outlet_nodes_info(self):
        """Search and identify the node upstream the outlet to apply boundary
        conditions

        This function is automatically called but for verification purposes we
        will test it. We will explore the horizontal link id of the outlet
        in different cases

        Examples
        --------

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))
        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
        >>> grid.at_node["surface_water__velocity"] = np.full(
        ...     grid.number_of_nodes, 0.25
        ... )
        >>> grid.at_link["surface_water__depth"] = map_mean_of_link_nodes_to_link(
        ...     grid, "surface_water__depth"
        ... )
        >>> grid.at_link["surface_water__velocity"] = map_mean_of_link_nodes_to_link(
        ...     grid, "surface_water__velocity"
        ... )

        In this topography the outlet is at the left edge

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.07, 1.08, 1.09, 1.09, 1.09],
        ...     [1.06, 1.07, 1.08, 1.09, 1.09],
        ...     [1.00, 1.03, 1.07, 1.08, 1.09],
        ...     [1.06, 1.07, 1.08, 1.09, 1.09],
        ...     [1.07, 1.08, 1.09, 1.09, 1.09],
        ... ]
        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> rbd = RiverBedDynamics(grid)
        >>> rbd.run_one_step()
        >>> rbd._outlet_links_horizontal
        array([[18],
               [19]])

        In this topography the outlet is at the top edge

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.09, 1.09, 1.09, 1.09, 1.09],
        ...     [1.09, 1.09, 1.08, 1.09, 1.09],
        ...     [1.09, 1.08, 1.07, 1.08, 1.09],
        ...     [1.08, 1.07, 1.03, 1.07, 1.08],
        ...     [1.07, 1.06, 1.00, 1.06, 1.07],
        ... ]
        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> rbd = RiverBedDynamics(grid)
        >>> rbd.run_one_step()
        >>> rbd._outlet_links_horizontal
        array([[28, 37],
               [29, 38]])

        In this topography the outlet is at the right edge

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.09, 1.09, 1.09, 1.08, 1.07],
        ...     [1.09, 1.09, 1.08, 1.07, 1.06],
        ...     [1.09, 1.08, 1.07, 1.03, 1.00],
        ...     [1.09, 1.09, 1.08, 1.07, 1.06],
        ...     [1.09, 1.09, 1.09, 1.08, 1.07],
        ... ]
        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> rbd = RiverBedDynamics(grid)
        >>> rbd.run_one_step()
        >>> rbd._outlet_links_horizontal
        array([[20],
               [21]])
        """

        # Gives the ID of the outlet node
        (self._out_id,) = np.where(self._grid.status_at_node == 1)

        # The outlet can be located on the bottom, top, right, or left edge of
        # the Raster. Depending on the edge, the upstream node is then identify
        # upstream1 is upstream of the upstream node

        if np.isin(self._out_id, self._grid.nodes_at_right_edge).any():
            self._upstream_out_id = self._out_id - 1
            up_upstream_out_id = self._upstream_out_id - 1
            up_up_upstream_out_id = up_upstream_out_id - 1
            out_direction = 0

        elif np.isin(self._out_id, self._grid.nodes_at_top_edge).any():
            self._upstream_out_id = self._out_id - self._grid._shape[1]
            up_upstream_out_id = self._upstream_out_id - self._grid._shape[1]
            up_up_upstream_out_id = up_upstream_out_id - self._grid._shape[1]
            out_direction = 1

        elif np.isin(self._out_id, self._grid.nodes_at_left_edge).any():
            self._upstream_out_id = self._out_id + 1
            up_upstream_out_id = self._upstream_out_id + 1
            up_up_upstream_out_id = up_upstream_out_id + 1
            out_direction = 2

        elif np.isin(self._out_id, self._grid.nodes_at_bottom_edge).any():
            self._upstream_out_id = self._out_id + self._grid._shape[1]
            up_upstream_out_id = self._upstream_out_id + self._grid._shape[1]
            up_up_upstream_out_id = up_upstream_out_id + self._grid._shape[1]
            out_direction = 3

        outlet_nodes = np.array((self._out_id, self._upstream_out_id))
        outlet_nodes = np.sort(outlet_nodes.flatten())
        outlet_links = np.unique(self._grid.links_at_node[outlet_nodes])
        outlet_links = outlet_links[np.where(outlet_links >= 0)]

        # 4 comes from the number of levels that are affected by the OverlandFlow boundary
        # conditions. 4 rows or columns of links.
        outlet_nodes_4 = np.array((up_upstream_out_id, up_up_upstream_out_id))
        outlet_nodes_4 = np.sort(outlet_nodes_4.flatten())
        outlet_links_4 = np.unique(self._grid.links_at_node[outlet_nodes_4])
        outlet_links_4 = outlet_links_4[np.where(outlet_links_4 >= 0)]
        outlet_links_4 = outlet_links_4[~np.in1d(outlet_links_4, outlet_links)]

        outlet_links_horizontal = outlet_links[
            np.in1d(outlet_links, self._grid.horizontal_links)
        ]
        outlet_links_vertical = outlet_links[
            np.in1d(outlet_links, self._grid.vertical_links)
        ]
        outlet_links_4_horizontal = outlet_links_4[
            np.in1d(outlet_links_4, self._grid.horizontal_links)
        ]
        outlet_links_4_vertical = outlet_links_4[
            np.in1d(outlet_links_4, self._grid.vertical_links)
        ]

        number_outlets = self._out_id.shape[0]

        # 3 comes from the number of levels that we are adjusting.
        if out_direction == 0:  # Outlet is at the right
            self._outlet_links_horizontal = np.reshape(
                outlet_links_horizontal, [2, number_outlets]
            )
            self._outlet_links_vertical = np.reshape(
                outlet_links_vertical, [number_outlets + 1, 2]
            )
            self._outlet_links_4_horizontal = np.reshape(
                outlet_links_4_horizontal, [2, number_outlets]
            )
            self._outlet_links_4_vertical = np.reshape(
                outlet_links_4_vertical, [number_outlets + 1, 2]
            )

        elif out_direction == 1:  # Outlet is at the top
            self._outlet_links_horizontal = np.reshape(
                outlet_links_horizontal, [2, number_outlets + 1]
            ).T
            self._outlet_links_vertical = np.reshape(
                outlet_links_vertical, [number_outlets, 2]
            ).T
            self._outlet_links_4_horizontal = np.reshape(
                outlet_links_4_horizontal, [2, number_outlets + 1]
            ).T
            self._outlet_links_4_vertical = np.reshape(
                outlet_links_4_vertical, [number_outlets, 2]
            ).T

        elif out_direction == 2:  # Outlet is at the  left edge
            self._outlet_links_horizontal = np.flip(
                np.reshape(outlet_links_horizontal, [2, number_outlets]), axis=1
            )
            self._outlet_links_vertical = np.flip(
                np.reshape(outlet_links_vertical, [number_outlets + 1, 2]), axis=1
            )
            self._outlet_links_4_horizontal = np.flip(
                np.reshape(outlet_links_4_horizontal, [2, number_outlets]), axis=1
            )
            self._outlet_links_4_vertical = np.flip(
                np.reshape(outlet_links_4_vertical, [number_outlets + 1, 2]), axis=1
            )

        elif out_direction == 3:  # Outlet is at the bottom edge
            self._outlet_links_horizontal = np.flip(
                np.reshape(outlet_links_horizontal, [2, number_outlets + 1]).T, axis=1
            )
            self._outlet_links_vertical = np.flip(
                np.reshape(outlet_links_vertical, [number_outlets, 2]).T, axis=1
            )
            self._outlet_links_4_horizontal = np.flip(
                np.reshape(outlet_links_4_horizontal, [2, number_outlets + 1]).T, axis=1
            )
            self._outlet_links_4_vertical = np.flip(
                np.reshape(outlet_links_4_vertical, [number_outlets, 2]).T, axis=1
            )

    def map_gsd_from_link_to_node(self, location="bed_surface"):
        """Map the bed surface grain size distribution from links to nodes.

        Given that all our calculations are conducted in links we implemented
        this function to display results in a raster or in nodes.
        default type is bedload, alternative type='bedload'
        """
        gsd_nodes = np.zeros([self.grid.number_of_nodes, self._gsd.shape[0] - 1])

        if location == "bed_surface":
            gsd_l = self._bed_surface__gsd_link

        else:
            gsd_l = self._sed_transport__bedload_gsd_link

        # split into grain sizes
        for i in range(self._gsd.shape[1] - 1):
            for j in range(self.grid.number_of_nodes):
                gsd_nodes_0 = gsd_l[self.grid.links_at_node[j]][:, i]
                mask = gsd_nodes_0 > 0.0
                gsd_nodes_0 = gsd_nodes_0[mask]
                if gsd_nodes_0.size > 0:
                    gsd_nodes[j, i] = np.mean(gsd_nodes_0)

        if location == "bed_surface":
            self._bed_surface__gsd_node = gsd_nodes
            # Revert any changes to the fixed GSD nodes and fixed elevations nodes
            fixed_nodes = np.unique(
                np.hstack((self._fixed_surface_gsd_nodes_id, self._fixed_nodes_id))
            )
            self._bed_surface__gsd_node[
                fixed_nodes
            ] = self._bed_surface__gsd_original_node[fixed_nodes]
        else:
            self._sed_transport__bedload_gsd_node = gsd_nodes

    def update_bed_surface_properties(self):
        """Calculates the updated GSD properties"""
        number_links = self._grid.number_of_links
        grain_size_D = self._grain_size_D_original  # Grain sizes
        grain_size_D_equivalent = self._grain_size  # Equivalent grain sizes

        # Equivalent grain sizes frequency
        grain_size_D_equiv_freq = self._bed_surface__gsd_link
        # Grain sizes frequency
        grain_size_frequency = np.hstack(
            (
                grain_size_D_equiv_freq,
                np.zeros([grain_size_D_equiv_freq.shape[0], 1]),
            )
        )

        # Calculates median_size_D50 at each node based on the updated GSD
        # Grain sizes in Psi scale
        grain_size_Psi_scale_D = np.flip(np.log2(grain_size_D), axis=0)
        grain_size_D_equiv_freq_cum = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each link

        # Finds the index of the grain size smaller than 50% in each link
        i0 = np.argmin(grain_size_D_equiv_freq_cum <= 0.5, axis=1) - 1
        # Finds the index of the grain size larger than 50% in each link
        i1 = np.argmax(grain_size_D_equiv_freq_cum > 0.5, axis=1)
        linkList = np.arange(number_links)

        grain_size_Psi_scale_D50 = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equiv_freq_cum[linkList, i1]
                - grain_size_D_equiv_freq_cum[linkList, i0]
            )
        ) * (0.5 - grain_size_D_equiv_freq_cum[linkList, i0])
        median_size_D50 = (
            2**grain_size_Psi_scale_D50
        )  # Median grain size in each link

        # Calculates the geometric mean and standard deviation
        # Equivalent grain sizes in Psi scale
        grain_size_D_equivalent_Psi = np.log2(grain_size_D_equivalent)
        grain_size_D_equivalent_Psi_mean = np.sum(
            grain_size_D_equiv_freq * grain_size_D_equivalent_Psi, axis=1
        )
        # Geometric mean size in each link
        grain_size_geom_mean = 2**grain_size_D_equivalent_Psi_mean

        # Equivalent grain sizes in Psi scale at each link
        grain_size_D_equivalent_Psi = np.tile(
            grain_size_D_equivalent_Psi, (number_links, 1)
        )
        grain_size_D_equivalent_Psi_mean = np.reshape(
            grain_size_D_equivalent_Psi_mean, [number_links, 1]
        )
        # Standard deviation at each node
        grain_size_geo_std = 2 ** np.sqrt(
            np.sum(
                ((grain_size_D_equivalent_Psi - grain_size_D_equivalent_Psi_mean) ** 2)
                * grain_size_D_equiv_freq,
                axis=1,
            )
        )

        self._bed_surface__median_size_link = median_size_D50
        self._bed_surface__geom_mean_size_link = grain_size_geom_mean
        self._bed_surface__geo_std_size_link = grain_size_geo_std

        self._bed_surface__act_layer_thick_prev_time_link = copy.deepcopy(
            self._bed_surface__act_layer_thick_link
        )
        self.calculate_act_layer_thick()

    def shear_stress(self):
        """Unsteady shear stress calculated at links

        Shear stress is calculated as::

            τ = rho * g * h * sf

        where sf is the unsteady friction slope and is calculated as::

            sf = S0 - dh/ds - U/g du/ds - 1/g du/dt

        Alternatively, τ can be calculated as::

            τ = rho * g * rh * sf

        but need to be manually selected (for consistency in the code)

        The term ds indicates a certain direction, X and Y in this case.
        All derivatives are first or second order approximations in each
        direction. Most of grid parameters are read again because they may have
        changed in time in OverlandFlow or when the bed evolved.
        """

        # Reads the current topographic elevation
        z = self._grid.at_node["topographic__elevation"]

        # First term in the friction slope - gradient of bed elevation
        # S0 = - dz/ds
        self._dz_ds = -self._grid.calc_grad_at_link(z)

        # Second term in the friction slope - gradient of water depth
        h = self._grid["node"]["surface_water__depth"]
        dh_ds = self._grid.calc_grad_at_link(h)

        h_links = self._grid.at_link["surface_water__depth"]

        # Third term in the friction slope - gradient of flow velocity
        # Velocity at current time step
        self._u = self._grid["link"]["surface_water__velocity"]

        # Velocity gradients are calculated in each direction
        du_ds = np.zeros_like(self._dz_ds)

        # In the X direction only horizontal gradients are valid. Here we use
        # the map_mean_of_horizontal_links_to_node method. However, this tool
        # is valid only for the X direction
        u_nodes_horizontal = self._grid.map_mean_of_horizontal_links_to_node(self._u)
        du_ds[self._grid.horizontal_links] = self._grid.calc_grad_at_link(
            u_nodes_horizontal
        )[self._grid.horizontal_links]

        # In the Y direction an extra step is required.
        u_nodes_vertical = self._grid.map_mean_of_vertical_links_to_node(self._u)
        u_nodes_vertical = np.flip(
            np.flip(
                np.reshape(
                    u_nodes_vertical, (self._grid._shape[0], self._grid._shape[1])
                )
            ),
            axis=1,
        )

        du_ds_vertical = np.zeros(
            (u_nodes_vertical.shape[0] - 1, u_nodes_vertical.shape[1])
        )

        for i in np.arange(u_nodes_vertical.shape[0] - 2, -1, -1):
            du_ds_vertical[i, :] = (
                u_nodes_vertical[i, :] - u_nodes_vertical[i + 1, :]
            ) / self._grid.dy

        # Three operations are in the right side. Here we recover landLab
        # indexing order
        du_ds[self._grid.vertical_links] = np.flip(du_ds_vertical.T, axis=1).flatten(
            order="F"
        )

        # Fourth term in the friction slope - rate of change of velocity in time
        u_at_prev_time = self._surface_water__velocity_prev_time_link
        du_dt = (self._u - u_at_prev_time) / self._grid._dt

        # Friction slope calculation at links including unsteady effects
        sf = self._dz_ds - dh_ds - (self._u / self._g) * du_ds - 1 / self._g * du_dt

        # And finally, the shear stress at links including unsteady effects
        if self._use_hydraulics_radius_in_shear_stress is True:
            # uses hydraulics ratio, so shear stress is calculated as
            # τ = rho * g * rh * sf

            # Different grid sizes (dx~=dy) are possible
            area = np.zeros_like(h_links)
            area[self._grid.horizontal_links] = (
                h_links[self._grid.horizontal_links] * self._grid.dx
            )
            area[self._grid.vertical_links] = (
                h_links[self._grid.vertical_links] * self._grid.dy
            )

            perimeter = np.zeros_like(h_links)
            perimeter[self._grid.horizontal_links] = (
                self._grid.dx + 2 * h_links[self._grid.horizontal_links]
            )
            perimeter[self._grid.vertical_links] = (
                self._grid.dy + 2 * h_links[self._grid.vertical_links]
            )
            rh = area / perimeter  # rh = wetted area / wetted perimeter
            self._shear_stress = self._rho * self._g * rh * sf
        else:
            # Equation is τ = rho * g * h * sf """
            self._shear_stress = self._rho * self._g * h_links * sf

        # Apply a boundary condition of zero flux at link next to borders
        self._shear_stress[self._links_at_border_cells] = 0

        # Direction of flux will be recovered at the end of the bedload
        # transport routine
        self._shear_stress_total = np.abs(self._shear_stress)

        # Now we write the shear stress field into the grid so other components
        # or postprocess functions can access it
        self._surface_water__shear_stress_link = self._shear_stress_total

    def Parker_1990_interp_points(self):
        """The following variables: po, oo, so are the interpolation points for
        the omega and sigma functions in Parker 1990"""

        po = [
            [0.6684, 0.7639, 0.8601, 0.9096, 0.9615, 1.000],
            [1.055, 1.108, 1.197, 1.302, 1.407, 1.529],
            [1.641, 1.702, 1.832, 1.937, 2.044, 2.261],
            [2.499, 2.732, 2.993, 3.477, 4.075, 4.469],
            [5.016, 6.158, 7.821, 10.06, 14.38, 19.97],
            [25.79, 38.57, 68.74, 91.95, 231.2, 2320],
        ]
        self._po = np.array(po).flatten()

        oo = [
            [1.011, 1.011, 1.01, 1.008, 1.004, 0.9997],
            [0.9903, 0.9789, 0.9567, 0.9273, 0.8964, 0.8604],
            [0.8287, 0.8123, 0.7796, 0.7554, 0.7326, 0.6928],
            [0.6585, 0.6345, 0.615, 0.5877, 0.564, 0.5523],
            [0.5395, 0.5209, 0.5045, 0.4917, 0.479, 0.4712],
            [0.4668, 0.462, 0.4578, 0.4564, 0.4541, 0.4527],
        ]
        self._oo = np.array(oo).flatten()

        so = [
            [0.8157, 0.8157, 0.8182, 0.8233, 0.8333, 0.8439],
            [0.8621, 0.8825, 0.9214, 0.9723, 1.025, 1.083],
            [1.13, 1.153, 1.196, 1.225, 1.25, 1.287],
            [1.313, 1.333, 1.352, 1.38, 1.403, 1.414],
            [1.426, 1.444, 1.458, 1.469, 1.48, 1.486],
            [1.49, 1.493, 1.497, 1.498, 1.499, 1.5],
        ]
        self._so = np.array(so).flatten()

    def bedload_equation_Parker_1990(self):
        """Surface-based bedload transport equation of Parker 1990

        G. Parker (1990) Surface-based bedload transport relation for gravel rivers,
        Journal of Hydraulic Research, 28:4, 417-436, DOI: 10.1080/00221689009499058
        """

        # Variables definition - All these variables are updated each time step.
        # Therefore, we read them again to ensure they are updated

        grain_size_frequency = self._bed_surface__gsd_link
        grain_size_geom_mean = self._bed_surface__geom_mean_size_link
        grain_size_geo_std = self._bed_surface__geo_std_size_link

        shear_stress_star_sg = self._shear_stress_total / (
            self._rho * self._R * self._g * (grain_size_geom_mean / 1000)
        )
        self._phi_sgo = shear_stress_star_sg / self._shear_stress_star_rsgo

        self.strain_functions()
        self._omega = 1 + (np.log2(grain_size_geo_std) / self._sigma0) * (
            self._omega0 - 1
        )

        grain_size_Di_Dsg = np.tile(
            self._grain_size, (self._grid.number_of_links, 1)
        ) / (np.reshape(grain_size_geom_mean, [self._grid.number_of_links, 1]))
        phi_sgo = np.reshape(self._phi_sgo, [self._phi_sgo.shape[0], 1])
        omega = np.reshape(self._omega, [self._omega.shape[0], 1])
        phi_i = omega * phi_sgo * (grain_size_Di_Dsg) ** -self._beta

        qb_G = np.zeros_like(phi_i)

        # There are three intervals where qb_G is evaluated
        (id0, id1) = np.where(phi_i > 1.59)
        if id0.shape[0] > 0:
            qb_G[id0, id1] = 5474 * (1 - 0.853 / phi_i[id0, id1]) ** 4.5

        (id0, id1) = np.where((phi_i >= 1) & (phi_i <= 1.59))
        if id0.shape[0] > 0:
            qb_G[id0, id1] = np.exp(
                14.2 * (phi_i[id0, id1] - 1) - 9.28 * (phi_i[id0, id1] - 1) ** 2
            )

        (id0, id1) = np.where(phi_i < 1)
        if id0.shape[0] > 0:
            qb_G[id0, id1] = phi_i[id0, id1] ** 14.2

        qb_Gf = grain_size_frequency * qb_G
        qb_Gf_sum = np.sum(qb_Gf, axis=1)

        # Total bedload transport rate
        qb_W_star_s = 0.00218 * qb_Gf_sum
        self._sed_transport__bedload_rate_link = (
            ((np.sqrt(self._shear_stress_total / self._rho)) ** 3 * qb_W_star_s)
            / (self._R * self._g)
        ) * np.sign(self._shear_stress)

        # Now that bed load rate has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate is True:
            link_id_fixed_links = np.unique(self._fixed_surface_gsd_link_id)

            self._sed_transport__bedload_rate_link[
                link_id_fixed_links
            ] = self._sed_transport__bedload_rate_imposed_link[link_id_fixed_links]

        # During the calculation of fractional sediment transport, there might
        # be instances where certain fractions will be zero. This situation
        # could trigger a division by zero warning. However, we are aware that
        # in these specific cases, there's no transport occurring. Thus, such a
        # warning is not indicative of an actual problem in the calculation.
        # To avoid unnecessary warnings, we've deactivated them.

        np.seterr(divide="ignore", invalid="ignore")

        # Frational bedload transport rate
        p = np.zeros_like(qb_Gf)
        qb_Gf_sum = np.transpose(np.tile(qb_Gf_sum, (grain_size_frequency.shape[1], 1)))
        p = qb_Gf / qb_Gf_sum
        id0 = np.isnan(p)
        p[id0] = 0
        self._sed_transport__bedload_gsd_link = p

        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there is no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate:
            linkId_0 = self._fixed_surface_gsd_link_id
            linkId_1 = self._fixed_surface_gsd_link_fractions

            self._sed_transport__bedload_gsd_link[
                linkId_0, linkId_1
            ] = self._sed_transport__bedload_gsd_imposed_link[linkId_0, linkId_1]

    def strain_functions(self):
        omega0 = np.zeros_like(self._phi_sgo)
        sigma0 = np.zeros_like(self._phi_sgo)

        # There are three intervals where we can interpolate Omega and Sigma
        (link_id_phi_sgo,) = np.where(self._phi_sgo <= 0.7639)
        if link_id_phi_sgo.shape[0] > 0:
            omega0[link_id_phi_sgo] = self._oo[0]
            sigma0[link_id_phi_sgo] = self._so[0]

        (link_id_phi_sgo,) = np.where(self._phi_sgo > 231.2)
        if link_id_phi_sgo.shape[0] > 0:
            omega0[link_id_phi_sgo] = np.interp(self._phi_sgo, self._po, self._oo)[
                link_id_phi_sgo
            ]
            sigma0[link_id_phi_sgo] = np.interp(self._phi_sgo, self._po, self._oo)[
                link_id_phi_sgo
            ]

        (link_id_phi_sgo,) = np.where(
            (self._phi_sgo > 0.7639) & (self._phi_sgo < 231.2)
        )
        if link_id_phi_sgo.shape[0] > 0:
            foo = interp1d(self._po, self._oo, kind="cubic")
            fso = interp1d(self._po, self._so, kind="cubic")
            omega0[link_id_phi_sgo] = foo(self._phi_sgo[link_id_phi_sgo])
            sigma0[link_id_phi_sgo] = fso(self._phi_sgo[link_id_phi_sgo])

        self._omega0 = omega0
        self._sigma0 = sigma0

    def bedload_equation_Wilcock_Crowe_2003(self):
        """Surface-based bedload transport equation of Wilcock and Crowe 2003

        Wilcock, P. R., & Crowe, J. C. (2003). Surface-based transport model
        for mixed-size sediment. Journal of hydraulic engineering, 129(2),
        120-128.
        """

        # Variables definition - All these variables are updated each time step.
        # Therefore, we read them again to ensure they are updated
        grain_size_frequency = self._bed_surface__gsd_link
        grain_size_geom_mean = self._bed_surface__geom_mean_size_link
        sand_fraction = self._bed_surface__sand_fraction_link

        shear_stress_star_sg = self._shear_stress_total / (
            self._rho * self._R * self._g * (grain_size_geom_mean / 1000)
        )
        shear_stress_star_rsg0 = 0.021 + 0.015 * np.exp(-20 * sand_fraction)
        phi_sg0 = shear_stress_star_sg / shear_stress_star_rsg0

        phi_i = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))

        qb_G = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))
        p = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))

        for i in np.arange(0, self._grain_size.shape[0]):
            b = 0.67 / (1 + np.exp(1.5 - self._grain_size[i] / grain_size_geom_mean))
            phi_i[:, i] = phi_sg0 * (self._grain_size[i] / (grain_size_geom_mean)) ** (
                -b
            )

        # There are two intervals where qb_G is evaluated
        (id0, id1) = np.where(phi_i >= 1.35)
        if id0.shape[0] > 0:
            qb_G[id0, id1] = 14 * (1 - 0.894 / (phi_i[id0, id1] ** 0.5)) ** 4.5

        (id0, id1) = np.where(phi_i < 1.35)
        if id0.shape[0] > 0:
            qb_G[id0, id1] = 0.002 * phi_i[id0, id1] ** 7.5

        qb_W_star_i = grain_size_frequency * qb_G
        qb_W_star = np.sum(qb_W_star_i, axis=1)

        # Total bedload transport rate is calculated at each link
        # Notice that qb_W_star includes fi (Eq 2 in the paper)
        self._sed_transport__bedload_rate_link = (
            ((np.sqrt(self._shear_stress_total / self._rho)) ** 3 * qb_W_star)
            / (self._R * self._g)
        ) * np.sign(self._shear_stress)

        # During the calculation of fractional sediment transport, there might
        # be instances where certain fractions will be zero. This situation
        # could trigger a division by zero warning. However, we are aware that
        # in these specific cases, there's no transport occurring. Thus, such a
        # warning is not indicative of an actual problem in the calculation.
        # To avoid unnecessary warnings, we've deactivated them.

        np.seterr(divide="ignore", invalid="ignore")
        p = qb_W_star_i / np.reshape(qb_W_star, [qb_W_star.shape[0], 1])
        id0 = np.isnan(p)
        p[id0] = 0
        self._sed_transport__bedload_gsd_link = p

        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there is no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate is True:
            linkId_0 = self._fixed_surface_gsd_link_id
            linkId_1 = self._fixed_surface_gsd_link_fractions

            self._sed_transport__bedload_gsd_link[
                linkId_0, linkId_1
            ] = self._sed_transport__bedload_gsd_imposed_link[linkId_0, linkId_1]

    def bedload_equation_MeyerPeter_Muller(self):
        """Surface-based bedload transport equation of Meyer-Peter and Müller

        Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport,
        Proceedings, 2nd Congress, International Association of Hydraulic
        Research, Stockholm: 39-64.
        """
        median_size_D50 = self._bed_surface__median_size_link
        shear_stress_star = self._shear_stress_total / (
            self._rho * self._R * self._g * (median_size_D50 / 1000)
        )
        self._critical_shear_stress_star = 0.047
        if self._variable_critical_shear_stress is True:
            self._critical_shear_stress_star = self.critical_shear_stress_star()
        qb_star = np.where(
            shear_stress_star - self._critical_shear_stress_star > 0,
            8 * np.abs(shear_stress_star - self._critical_shear_stress_star) ** (3 / 2),
            0,
        )

        self._sed_transport__bedload_rate_link = (
            qb_star
            * (
                np.sqrt(self._R * self._g * (median_size_D50 / 1000))
                * (median_size_D50 / 1000)
            )
        ) * np.sign(self._shear_stress)

    def bedload_equation_FernandezLuque_VanBeek(self):
        """Surface-based bedload transport equation of Fernandez Luque and
        van Beek

        Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of
        bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
        """
        median_size_D50 = self._bed_surface__median_size_link
        shear_stress_star = self._shear_stress_total / (
            self._rho * self._R * self._g * (median_size_D50 / 1000)
        )
        self._critical_shear_stress_star = 0.045
        if self._variable_critical_shear_stress is True:
            self._critical_shear_stress_star = self.critical_shear_stress_star()
        qb_star = np.where(
            shear_stress_star - self._critical_shear_stress_star > 0,
            5.7
            * np.abs(shear_stress_star - self._critical_shear_stress_star) ** (3 / 2),
            0,
        )

        # Total bedload transport rate is calculated at each link
        self._sed_transport__bedload_rate_link = (
            qb_star
            * (
                np.sqrt(self._R * self._g * (median_size_D50 / 1000))
                * (median_size_D50 / 1000)
            )
        ) * np.sign(self._shear_stress)

    def critical_shear_stress_star(self):
        """Uses the current bed slope at links for calculating a slope-dependent
        critical shear stress. This assumes that the flow direction always goes
        in the direction of descending slope. However, there might be local
        cases where this assumption does not hold true
        """
        bed_slope = np.abs(self._dz_ds)
        # Mueller et al. (2005) equation is used in steep slopes
        return np.where(
            bed_slope > 0.03, 2.18 * bed_slope + 0.021, self._critical_shear_stress_star
        )

    def calculate_net_bedload(self):
        """Calculates the net volumetric bedload coming from all links (m2/s)
        onto nodes (m3/s).

        This method takes the volumetric bedload entering and exiting through a
        face and determines the net volumetric bedload on a given node.
        """

        # Reads and modify the field self._sed_transport__bedload_rate_link to account
        # for links where sediment supply is imposed.
        (self._id_link_upstream_sediment_supply,) = np.nonzero(
            self._sed_transport__bedload_rate_imposed_link
        )
        self._sed_transport__bedload_rate_link[
            self._id_link_upstream_sediment_supply
        ] = self._sed_transport__bedload_rate_imposed_link[
            self._id_link_upstream_sediment_supply
        ]

        qb_x = (
            np.sum(
                self._sed_transport__bedload_rate_link[
                    self._grid.links_at_node[:, [0, 2]]
                ]
                * self._normal[:, [0, 2]],
                axis=1,
            )
            * self._grid.dy
        )
        qb_y = (
            np.sum(
                self._sed_transport__bedload_rate_link[
                    self._grid.links_at_node[:, [1, 3]]
                ]
                * self._normal[:, [1, 3]],
                axis=1,
            )
            * self._grid.dx
        )

        net_qb_nodes = qb_x + qb_y

        ## At the boundary, there is no exiting link, so we assume a zero flux
        # exiting. This assumption is overridden in the Exner equation, where a
        # zero gradient boundary condition is employed.

        net_qb_nodes[self._grid.boundary_nodes] = 0

        self._sed_transport__net_bedload_node = net_qb_nodes

    def update_bed_elevation(self):
        """Applies the Exner equation and boundary conditions to predict
        the change in bed surface elevation.
        """

        # Bed elevation is updated using Exner equation
        area = self._grid.dx * self._grid.dy
        d_qb = self._sed_transport__net_bedload_node
        z0 = copy.deepcopy(self._grid["node"]["topographic__elevation"])
        dz = -self._grid._dt / ((1 - self._lambda_p) * area) * d_qb

        self._grid["node"]["topographic__elevation"] += dz

        # The outlet node has been modifed, but not according to the specified
        # boundary condition. Here, we return the outlet to the previous state
        self._grid["node"]["topographic__elevation"][self._out_id] = z0[self._out_id]

        # Fixed nodes may have been modifed, but not according to the specified
        # condition. Here, we return the nodes to the orignal state
        self._grid["node"]["topographic__elevation"][
            self._fixed_nodes_id
        ] = self._topographic__elevation_original_node[self._fixed_nodes_id]

        # Now we can apply the boundary conditions
        if self._outlet_boundary_condition == "zeroGradient":
            dz_outlet = dz[self._upstream_out_id]
        elif self._outlet_boundary_condition == "fixedValue":
            dz_outlet = 0

        self._grid["node"]["topographic__elevation"][self._out_id] = (
            z0[self._out_id] + dz_outlet
        )

        # Now we map data into links to update Bed GSD
        self._grid["link"][
            "topographic__elevation"
        ] = self._grid.map_mean_of_link_nodes_to_link(
            self._grid["node"]["topographic__elevation"]
        )

        # We correct the water depth given the change in bed elevation, preserving WSE
        # Keeps the discharge, so mass conservation is preserved
        dz_current_time = z0 - self._grid["node"]["topographic__elevation"]
        (self._id_eroded_nodes,) = np.where(dz_current_time < 0)
        self._id_eroded_nodes = self._id_eroded_nodes.astype(int)

        if self._id_eroded_nodes.shape[0] > 0:
            self._grid["node"]["surface_water__depth"][self._id_eroded_nodes] += np.abs(
                dz_current_time[self._id_eroded_nodes]
            )

        # Here we register how deep the deposited/eroded layer is
        if self._track_stratigraphy:
            self._bed_surface__surface_thickness_new_layer_link = (
                self._grid["link"]["topographic__elevation"]
                - self._topographic__elevation_subsurface_link
            )

            # Checks if deposited material needs to be updated
            (self._id_deep_links,) = np.where(
                self._bed_surface__surface_thickness_new_layer_link
                > self._new_surface_layer_thickness
            )
            self._id_deep_links = self._id_deep_links[
                np.in1d(self._id_deep_links, self._grid.active_links)
            ]
            self._id_deep_links = self._id_deep_links[
                ~np.in1d(self._id_deep_links, self._fixed_links)
            ]

            if self._id_deep_links.shape[0] > 0:
                self._compute_stratigraphy = True
                self._update_subsurface = True

            # Checks if eroded material needs to be updated
            (self._id_eroded_links,) = np.where(
                self._bed_surface__surface_thickness_new_layer_link
                < -self._new_surface_layer_thickness
            )
            self._id_eroded_links = self._id_eroded_links[
                np.in1d(self._id_eroded_links, self._grid.active_links)
            ]
            self._id_eroded_links = self._id_eroded_links[
                ~np.in1d(self._id_eroded_links, self._fixed_links)
            ]

            if self._id_eroded_links.shape[0] > 0:
                self._compute_stratigraphy = True
                self._update_subsurface_eroded = True

    def update_bed_surface_gsd(self):
        """Uses the fractional Exner equation to update the bed GSD"""
        # Here we create a number of variables that will be used in the
        # following definitions to make the code a little bit cleaner
        # No deep copies are done to avoid unnecesary repetition.

        number_links = self._grid.number_of_links
        num_col = self._grid.number_of_node_columns
        gsd_F = self._bed_surface__gsd_link
        gsd_Fs = self._bed_subsurface__gsd_link
        pl = self._sed_transport__bedload_gsd_link
        qbT = (
            self._sed_transport__bedload_rate_link
        )  # total bed load transport rate at each link
        active_layer_La = np.reshape(
            self._bed_surface__act_layer_thick_link, [number_links, 1]
        )
        active_layer_La_old = np.reshape(
            self._bed_surface__act_layer_thick_prev_time_link,
            [number_links, 1],
        )

        lps = self._lambda_p
        dx = self._grid.dx
        dy = self._grid.dy
        alpha = self._alpha
        dt = self._grid._dt

        dv = 2 * num_col - 1

        qbT = np.reshape(qbT, [number_links, 1])
        qbTdev = np.zeros([number_links, 1])
        qbTdevNeg = np.zeros([number_links, 1])
        qjj1dev = np.zeros_like(pl)
        qjj1devNeg = np.zeros_like(pl)

        # Horizontal Links
        hlL = np.arange(
            0, number_links - num_col + 2, 2 * num_col - 1
        )  # Horizontal Links at left edge
        hlR = np.arange(
            num_col - 2, number_links, 2 * num_col - 1
        )  # Horizontal Links at right edge
        hl_Id = np.in1d(
            self._grid.horizontal_links, np.hstack((hlL, hlR))
        )  # Links at the top row
        hl = self._grid.horizontal_links[~hl_Id]

        # Vertical Links
        vl = self._grid.vertical_links[num_col:-num_col]  # Links within middle region
        vlB = self._grid.vertical_links[0:num_col]  # Links at the bottom row
        vlT = self._grid.vertical_links[-num_col:None]  # Links at the top row

        # First we assume that everywhere flow direction is in the positive direction
        qbTdev[hl] = (
            alpha * (qbT[hl] - qbT[hl - 1]) / dy
            + (1 - alpha) * (qbT[hl + 1] - qbT[hl]) / dy
        )
        qbTdev[hlL] = (qbT[hlL + 1] - qbT[hlL]) / dy
        qbTdev[hlR] = (qbT[hlR] - qbT[hlR - 1]) / dy

        qjj1dev[hl, :] = (
            alpha * (qbT[hl] * pl[hl, :] - qbT[hl - 1] * pl[hl - 1, :]) / dy
            + (1 - alpha) * (qbT[hl + 1] * pl[hl + 1] - qbT[hl] * pl[hl]) / dy
        )
        qjj1dev[hlL, :] = (qbT[hlL + 1] * pl[hlL + 1, :] - qbT[hlL] * pl[hlL, :]) / dy
        qjj1dev[hlR, :] = (qbT[hlR] * pl[hlR, :] - qbT[hlR - 1] * pl[hlR - 1, :]) / dy

        qbTdev[vl] = (
            alpha * (qbT[vl] - qbT[vl - dv]) / dx
            + (1 - alpha) * (qbT[vl + dv] - qbT[vl]) / dx
        )
        qbTdev[vlB] = (qbT[vlB + dv] - qbT[vlB]) / dx
        qbTdev[vlT] = (qbT[vlT] - qbT[vlT - dv]) / dx

        qjj1dev[vl, :] = (
            alpha * (qbT[vl] * pl[vl, :] - qbT[vl - dv] * pl[vl - dv, :]) / dx
            + (1 - alpha) * (qbT[vl + dv] * pl[vl + dv] - qbT[vl] * pl[vl]) / dx
        )
        qjj1dev[vlB, :] = (qbT[vlB + dv] * pl[vlB + dv, :] - qbT[vlB] * pl[vlB, :]) / dx
        qjj1dev[vlT, :] = (qbT[vlT] * pl[vlT, :] - qbT[vlT - dv] * pl[vlT - dv, :]) / dx

        # Now we correct for flow at locations where it is flowing towards the
        # negative direction
        (hl_neg_id,) = np.where(qbT[hl][:, 0] < 0)
        (hlL_neg_id,) = np.where(qbT[hlL][:, 0] < 0)
        (hlR_neg_id,) = np.where(qbT[hlR][:, 0] < 0)

        qbTdevNeg[hl] = (
            alpha * (qbT[hl] - qbT[hl + 1]) / dy
            + (1 - alpha) * (qbT[hl - 1] - qbT[hl]) / dy
        )
        qbTdevNeg[hlL] = (qbT[hlL] - qbT[hlL - 1]) / dy
        qbTdevNeg[hlR] = (qbT[hlR - 1] - qbT[hlR]) / dy

        qbTdev[hl[hl_neg_id]] = -qbTdevNeg[hl[hl_neg_id]]
        qbTdev[hlL[hlL_neg_id]] = -qbTdevNeg[hlL[hlL_neg_id]]
        qbTdev[hlR[hlR_neg_id]] = -qbTdevNeg[hlR[hlR_neg_id]]

        qjj1devNeg[hl, :] = (
            alpha * (qbT[hl] * pl[hl, :] - qbT[hl + 1] * pl[hl + 1, :]) / dy
            + (1 - alpha) * (qbT[hl - 1] * pl[hl - 1] - qbT[hl] * pl[hl]) / dy
        )
        qjj1devNeg[hlL, :] = (
            qbT[hlL] * pl[hlL, :] - qbT[hlL - 1] * pl[hlL - 1, :]
        ) / dy
        qjj1devNeg[hlR, :] = (
            qbT[hlR - 1] * pl[hlR - 1, :] - qbT[hlR] * pl[hlR, :]
        ) / dy

        qjj1dev[hl[hl_neg_id]] = -qjj1devNeg[hl[hl_neg_id]]
        qjj1dev[hlL[hlL_neg_id]] = -qjj1devNeg[hlL[hlL_neg_id]]
        qjj1dev[hlR[hlR_neg_id]] = -qjj1devNeg[hlR[hlR_neg_id]]

        (vl_neg_id,) = np.where(qbT[vl][:, 0] < 0)
        (vlB_neg_id,) = np.where(qbT[vlB][:, 0] < 0)
        (vlT_neg_id,) = np.where(qbT[vlT][:, 0] < 0)

        qbTdevNeg[vl] = (
            alpha * (qbT[vl] - qbT[vl + dv]) / dx
            + (1 - alpha) * (qbT[vl - dv] - qbT[vl]) / dx
        )
        qbTdevNeg[vlB] = (qbT[vlB] - qbT[vlB + dv]) / dx
        qbTdevNeg[vlT] = (qbT[vlT - dv] - qbT[vlT]) / dx

        qbTdev[vl[vl_neg_id]] = -qbTdevNeg[vl[vl_neg_id]]
        qbTdev[vlB[vlB_neg_id]] = -qbTdevNeg[vlB[vlB_neg_id]]
        qbTdev[vlT[vlT_neg_id]] = -qbTdevNeg[vlT[vlT_neg_id]]

        qjj1devNeg[vl, :] = (
            alpha * (qbT[vl] * pl[vl, :] - qbT[vl + dv] * pl[vl + dv, :]) / dx
            + (1 - alpha) * (qbT[vl - dv] * pl[vl - dv] - qbT[vl] * pl[vl]) / dx
        )
        qjj1devNeg[vlB, :] = (
            qbT[vlB] * pl[vlB, :] - qbT[vlB + dv] * pl[vlB + dv, :]
        ) / dx
        qjj1devNeg[vlT, :] = (
            qbT[vlT - dv] * pl[vlT - dv, :] - qbT[vlT] * pl[vlT, :]
        ) / dx

        qjj1dev[vl[vl_neg_id]] = -qjj1devNeg[vl[vl_neg_id]]
        qjj1dev[vlB[vlB_neg_id]] = -qjj1devNeg[vlB[vlB_neg_id]]
        qjj1dev[vlT[vlT_neg_id]] = -qjj1devNeg[vlT[vlT_neg_id]]

        # Correction done

        gsd_FIexc = copy.deepcopy(gsd_Fs)
        (id0,) = np.where(qbTdev[:, 0] <= 0)
        gsd_FIexc[id0, :] = 0.7 * gsd_F[id0, :] + 0.3 * pl[id0, :]

        qjj2dev = gsd_FIexc * np.reshape(qbTdev, [number_links, 1])

        gsd_Fnew = gsd_F + dt * (-qjj1dev + qjj2dev) / (1 - lps) / active_layer_La
        if (self._track_stratigraphy is False) and (self._first_iteration is False):
            gsd_Fnew = gsd_Fnew + (gsd_FIexc - gsd_F) / active_layer_La * (
                active_layer_La - active_layer_La_old
            )
        elif (self._track_stratigraphy is True) and (self._first_iteration is False):
            if self._subsurface_changed is True:
                gsd_Fnew = gsd_Fnew + (gsd_FIexc - gsd_F) / active_layer_La * (
                    active_layer_La - active_layer_La_old
                )
                self._subsurface_changed = False

        (id0, id1) = np.where(gsd_Fnew <= 0)
        gsd_Fnew[id0, id1] = 0
        gsd_Fnew = gsd_Fnew / np.reshape(np.sum(gsd_Fnew, axis=1), [number_links, 1])
        gsd_Fnew = np.nan_to_num(gsd_Fnew)

        # Given the way in which OverlandFLow calculates flow near the outlets
        # we need to correct changes to GSD that are caused by the sudden drop in
        # water depth. This affects always the 3 rows or columns of links near the
        # outlet. Here we identify those 3 rows or columns

        # Corrects horizontal links

        # In very small grids we may not have enough links, so we limit the
        # correction to those that have more enough data

        if self._outlet_links_4_horizontal.shape[1] > 1:
            for i in np.arange(0, self._outlet_links_horizontal.shape[0]):
                ds = dx
                m = (
                    gsd_Fnew[self._outlet_links_4_horizontal[i][1]]
                    - gsd_Fnew[self._outlet_links_4_horizontal[i][0]]
                ) / ds
                b = gsd_Fnew[self._outlet_links_4_horizontal[i][0]]
                gsd_Fnew[self._outlet_links_horizontal[i][0]] = m * 2 * ds + b
                gsd_Fnew[self._outlet_links_horizontal[i][1]] = m * 3 * ds + b

        # Corrects vertical links

        # In very small grids we may not have enough links, so we limit the
        # correction to those that have more enough data
        if self._outlet_links_4_vertical.shape[1] > 1:
            for i in np.arange(0, self._outlet_links_vertical.shape[0]):
                ds = dy
                m = (
                    gsd_Fnew[self._outlet_links_4_vertical[i][1]]
                    - gsd_Fnew[self._outlet_links_4_vertical[i][0]]
                ) / ds
                b = gsd_Fnew[self._outlet_links_4_vertical[i][0]]
                gsd_Fnew[self._outlet_links_vertical[i][0]] = m * 2 * ds + b
                gsd_Fnew[self._outlet_links_vertical[i][1]] = m * 3 * ds + b

        # Now we update the bed surface GSD
        self._bed_surface__gsd_link = copy.deepcopy(gsd_Fnew)

        # If the bed is eroded below the original elevation it restores this
        # initial GSD. First looks if there is any such link
        (id_eroded_links,) = np.where(
            self._grid["link"]["topographic__elevation"]
            < self._topographic__elevation_original_link
        )

        if id_eroded_links.shape[0] > 0:
            self._bed_surface__gsd_link[
                id_eroded_links
            ] = self._bed_surface__gsd_original_link[id_eroded_links]

            # Now, an eroded node cannot return to the original GSD if starts depositing again.
            # It will only use the original GSD if erodes deeper than the maximum that has
            # been eroded
            self._topographic__elevation_original_link[id_eroded_links] = copy.deepcopy(
                self._grid["link"]["topographic__elevation"][id_eroded_links]
            )

        # Revert any changes to the fixed GSD nodes
        self._bed_surface__gsd_link[
            self._fixed_links
        ] = self._bed_surface__gsd_original_link[self._fixed_links]

        if self._track_stratigraphy is True:
            if (
                (self._first_iteration is True)
                or (
                    self._stratigraphy_cycle
                    >= self._number_cycles_to_process_stratigraphy
                )
                or (self._compute_stratigraphy is True)
            ):
                self.stratigraphy()

    def stratigraphy(self):
        """This function controls how the stratigraphy is beign updated

        Example where stratigraphy with deposition and erosion is tested below.
        This is a very slow test

        Examples
        --------

        >>> from shutil import rmtree
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics

        >>> grid = RasterModelGrid((8, 3), xy_spacing=100)

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.12, 1.00, 1.12, 1.12, 1.01, 1.12, 1.12, 1.01],
        ...     [1.12, 1.12, 1.01, 1.12, 1.12, 1.01, 1.12, 1.12],
        ...     [1.01, 1.12, 1.12, 1.01, 1.12, 1.12, 1.12, 1.12],
        ... ]

        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.40)
        >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.40)
        >>> grid.at_link["surface_water__velocity"] = np.full(
        ...     grid.number_of_links, 0.40
        ... )
        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> gsd = [[8, 100], [4, 90], [2, 0]]

        >>> fixed_nodes_id = np.array((1, 4))
        >>> fixed_nodes = np.zeros(grid.number_of_nodes)
        >>> fixed_nodes[fixed_nodes_id] = 1

        >>> fixed_bed_gsd_nodes = np.zeros(grid.number_of_nodes)
        >>> fixed_bed_gsd_nodes[fixed_nodes_id] = 1

        >>> in_qb = -0.002
        >>> in_l = np.array((28, 33))
        >>> qb = np.full(grid.number_of_links, 0.0)
        >>> qb[in_l] = in_qb

        >>> rbd = RiverBedDynamics(
        ...     grid,
        ...     gsd=gsd,
        ...     bedload_equation="Parker1990",
        ...     outlet_boundary_condition="fixedValue",
        ...     bed_surface__elevation_fixed_node=fixed_nodes,
        ...     sed_transport__bedload_rate_imposed_link=qb,
        ...     bed_surface__gsd_fixed_node=fixed_bed_gsd_nodes,
        ...     track_stratigraphy=True,
        ...     new_surface_layer_thickness=0.02,
        ...     number_cycles_to_process_stratigraphy=2,
        ... )

        >>> for t in range(1300):
        ...     rbd._current_t = t
        ...     rbd.run_one_step()
        ...
        >>> z = grid.at_node["topographic__elevation"][16]
        >>> np.round(z, decimals=3)
        1.05

        >>> rbd._sed_transport__bedload_rate_imposed_link[in_l] = 0.001
        >>> for t in range(1300, 2605 + 1300):
        ...     rbd._current_t = t
        ...     rbd.run_one_step()
        ...
        >>> z = grid.at_node["topographic__elevation"][16]
        >>> np.round(z, decimals=3)
        1.01

        We will delete any temp folder created during docstring tests

        >>> rmtree("stratigraphyRawData")
        >>> rmtree("stratigraphyTempFiles")
        """

        # Here we create a number of variables that will be used in the
        # following definitions to make the code a little bit cleaner

        z = self._grid["link"]["topographic__elevation"]
        grain_size_D = self._grain_size_D_original
        gsd_F = self._bed_surface__gsd_link
        gsd_Fs = self._bed_subsurface__gsd_link
        active_layer_La = self._bed_surface__act_layer_thick_link
        dzl = self._bed_surface__surface_thickness_new_layer_link

        if self._first_iteration:
            # Creates a folder to store results and a file for each active link and node
            if os.path.exists(self._stratigraphy_temp_files_path):
                shutil.rmtree(self._stratigraphy_temp_files_path)

            os.mkdir(self._stratigraphy_temp_files_path)

            if os.path.exists(self._stratigraphy_raw_data_path):
                shutil.rmtree(self._stratigraphy_raw_data_path)

            os.mkdir(self._stratigraphy_raw_data_path)

            # Now goes to the raw data folder to store the properties at time
            # zero. The format is unfriendly but after the simulation is done a
            # function will take care of making it more friendly.

            os.chdir(self._stratigraphy_raw_data_path)

            for i in self._grid.active_links:
                filename = "link_" + str(i) + ".txt"
                data = np.hstack((self._current_t, z[i], gsd_F[i, :], gsd_Fs[i, :]))
                data = np.reshape(data, [1, data.shape[0]])
                with open(filename, "ab") as f:
                    np.savetxt(f, data, "%.3f")

        # Here we store data
        os.chdir(self._cwd)
        os.chdir(self._stratigraphy_temp_files_path)
        for i in self._grid.active_links:
            filename = "link_" + str(i) + ".txt"
            data = np.hstack(
                (
                    self._current_t,
                    z[i],
                    dzl[i],
                    active_layer_La[i],
                    gsd_F[i, :],
                    gsd_Fs[i, :],
                )
            )
            data = np.reshape(data, [1, data.shape[0]])

            # Try up to 10 times
            for _attempt in range(10):
                try:
                    with open(filename, "ab") as f:
                        np.savetxt(f, data, "%.3f")
                    break  # Successfully written, exit the loop
                except PermissionError:
                    time.sleep(1)  # Wait for 1 second before retrying

        # Here we update stratigraphy in case of deposition
        # only for new layers
        os.chdir(self._cwd)
        os.chdir(self._stratigraphy_temp_files_path)
        if self._update_subsurface:
            for i in self._id_deep_links:
                link_data = np.loadtxt("link_" + str(i) + ".txt")
                mean_subsurface_gsd = np.mean(
                    link_data[:, 4 : 4 + grain_size_D.shape[0] - 1], axis=0
                )  # 4 is the number of elements before gsd_F GSD
                self._bed_subsurface__gsd_link[i, :] = mean_subsurface_gsd / np.sum(
                    mean_subsurface_gsd
                )
                os.remove("link_" + str(i) + ".txt")

            # Now updates the raw data for the links with deposition
            os.chdir(self._cwd)
            os.chdir(self._stratigraphy_raw_data_path)

            for i in self._id_deep_links:
                filename = "link_" + str(i) + ".txt"
                data = np.hstack(
                    (
                        self._current_t,
                        z[i],
                        gsd_F[i, :],
                        self._bed_subsurface__gsd_link[i, :],
                    )
                )
                data = np.reshape(data, [1, data.shape[0]])
                with open(filename, "ab") as f:
                    np.savetxt(f, data, "%.3f")

            self._bed_surface__surface_thickness_new_layer_link[self._id_deep_links] = (
                self._bed_surface__surface_thickness_new_layer_link[self._id_deep_links]
                - self._new_surface_layer_thickness
            )
            self._topographic__elevation_subsurface_link[self._id_deep_links] = (
                self._topographic__elevation_subsurface_link[self._id_deep_links]
                + self._new_surface_layer_thickness
            )

            self._subsurface_changed = True

        # Here we update the stratigraphy in case of erosion - only for new layers
        os.chdir(self._cwd)
        if self._update_subsurface_eroded is True:
            os.chdir(self._stratigraphy_temp_files_path)
            for i in self._id_eroded_links:
                os.remove("link_" + str(i) + ".txt")

            os.chdir(self._cwd)
            os.chdir(self._stratigraphy_raw_data_path)

            for i in self._id_eroded_links:
                link_data = np.loadtxt("link_" + str(i) + ".txt")
                if (link_data.shape[0] > 1) and (link_data.ndim > 1):
                    self._bed_subsurface__gsd_link[i, :] = link_data[
                        -2, (2 + self._gsd.shape[0] - 1) : None
                    ]
                    data = np.hstack(
                        (
                            self._current_t,
                            z[i],
                            gsd_F[i, :],
                            self._bed_subsurface__gsd_link[i, :],
                        )
                    )
                    data = np.reshape(data, [1, data.size])
                with open("link_" + str(i) + ".txt", "ab") as f:
                    np.savetxt(f, data, "%.3f")
                self._bed_surface__gsd_link[i, :] = self._bed_subsurface__gsd_link[i, :]

            self._bed_surface__surface_thickness_new_layer_link[
                self._id_eroded_links
            ] = (
                self._bed_surface__surface_thickness_new_layer_link[
                    self._id_eroded_links
                ]
                + self._new_surface_layer_thickness
            )
            self._topographic__elevation_subsurface_link[self._id_eroded_links] = (
                self._topographic__elevation_subsurface_link[self._id_eroded_links]
                - self._new_surface_layer_thickness
            )

            self._subsurface_changed = True

        os.chdir(self._cwd)
        self._stratigraphy_cycle = 0
        self._compute_stratigraphy = False
        self._update_subsurface = False
        self._update_subsurface_eroded = False

    def bedload_mapper(self):
        """Maps the magnitude of the volumetric bedload per unit width values
        from links (m2/s) onto nodes (m3/s).

        This method takes the bedload transport rates from links and calculates
        the magnitude at a given node.
        """

        qb_x_r = (
            self._sed_transport__bedload_rate_link[self._grid.links_at_node[:, 0]]
            * self._grid.dy
        )
        (link_id,) = np.where(self._grid.links_at_node[:, 0] < 0)
        qb_x_r[link_id] = 0

        qb_x_l = (
            self._sed_transport__bedload_rate_link[self._grid.links_at_node[:, 2]]
            * self._grid.dy
        )
        (link_id,) = np.where(self._grid.links_at_node[:, 2] < 0)
        qb_x_l[link_id] = 0

        qb_y_u = (
            self._sed_transport__bedload_rate_link[self._grid.links_at_node[:, 1]]
            * self._grid.dx
        )
        (link_id,) = np.where(self._grid.links_at_node[:, 1] < 0)
        qb_y_u[link_id] = 0

        qb_y_l = (
            self._sed_transport__bedload_rate_link[self._grid.links_at_node[:, 3]]
            * self._grid.dx
        )
        (link_id,) = np.where(self._grid.links_at_node[:, 3] < 0)
        qb_y_l[link_id] = 0

        qb_x = 0.5 * (qb_x_r + qb_x_l)
        qb_y = 0.5 * (qb_y_u + qb_y_l)

        self._bedloadRate_vector = np.transpose(np.vstack((qb_x, qb_y)))
        self._bedloadRate_magnitude = np.sqrt(qb_x**2 + qb_y**2)

    def calculate_DX(self, fX, mapped_in="link"):
        """Calculate the bed surface and bed load grain size corresponding to
        any fraction. For example, 50%, which is the median_size_D50.
        In that case use fX = 0.5

        This method takes the user specified fraction, from 0 to 1, and outputs
        the corresponding grain size in nodes or links. By default the link
        option is used. Use ``mapped_in="node"`` to calculate at nodes

        Examples
        --------

        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.07, 1.06, 1.00, 1.06, 1.07],
        ...     [1.08, 1.07, 1.03, 1.07, 1.08],
        ...     [1.09, 1.08, 1.07, 1.08, 1.09],
        ...     [1.09, 1.09, 1.08, 1.09, 1.09],
        ...     [1.09, 1.09, 1.09, 1.09, 1.09],
        ... ]

        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
        >>> grid.at_node["surface_water__velocity"] = np.full(
        ...     grid.number_of_nodes, 0.25
        ... )
        >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
        >>> grid.at_link["surface_water__velocity"] = np.full(
        ...     grid.number_of_links, 0.25
        ... )

        >>> gsd_loc = [
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ... ]

        >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]

        >>> rbd = RiverBedDynamics(
        ...     grid,
        ...     gsd=gsd,
        ...     bedload_equation="Parker1990",
        ...     bed_surface__gsd_location_node=gsd_loc,
        ... )
        >>> rbd.run_one_step()

        For simplicity, only a few nodes and links will be selected for displaying.

        >>> nodes = np.arange(2, 25, 5)
        >>> links = grid.active_links[np.in1d(grid.active_links, grid.vertical_links)]

        >>> (surface_D90, D90_bedload) = rbd.calculate_DX(0.9)
        >>> np.round(surface_D90[links], 3)
        array([[ 27.849],
               [ 27.858],
               [ 27.866],
               [ 27.858],
               [ 27.858],
               [ 27.858],
               [ 27.858]])

        >>> np.round(D90_bedload[links], 3)
        array([[  0.   ],
               [ 25.032],
               [ 27.655],
               [ 25.032],
               [ 25.032],
               [ 25.032],
               [ 25.032]])

        >>> (surface_D90, D90_bedload) = rbd.calculate_DX(0.9, mapped_in="node")
        >>> np.round(surface_D90[nodes], 3)
        array([[ 28.084],
               [ 27.858],
               [ 27.86 ],
               [ 27.858],
               [ 28.086]])

        >>> np.round(D90_bedload[nodes], 3)
        array([[ 27.737],
               [ 27.655],
               [ 25.945],
               [ 25.032],
               [  0.   ]])
        """
        nodes_list = np.arange(self._grid.number_of_nodes)
        grain_size_D = self._grain_size_D_original  # Grain sizes
        grain_size_Psi_scale_D = np.flip(
            np.log2(grain_size_D), axis=0
        )  # Grain sizes in Psi scale

        if mapped_in == "link":
            frequency_grain_size_list = np.arange(self._grid.number_of_links)
        else:
            frequency_grain_size_list = nodes_list

        def compute_grain_size_DX(grain_size_equiv_freq):
            grain_size_frequency = np.hstack(
                (grain_size_equiv_freq, np.zeros([grain_size_equiv_freq.shape[0], 1]))
            )
            grain_size_equiv_freq_cum = np.cumsum(
                np.flip(grain_size_frequency, axis=1), axis=1
            )
            i0 = np.argmin(grain_size_equiv_freq_cum <= fX, axis=1) - 1
            i1 = np.argmax(grain_size_equiv_freq_cum > fX, axis=1)

            grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
                (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
                / (
                    grain_size_equiv_freq_cum[frequency_grain_size_list, i1]
                    - grain_size_equiv_freq_cum[frequency_grain_size_list, i0]
                )
            ) * (fX - grain_size_equiv_freq_cum[frequency_grain_size_list, i0])
            return 2**grain_size_Psi_scale_DX

        grain_size_DX_surface = compute_grain_size_DX(
            self._bed_surface__gsd_link
            if mapped_in == "link"
            else self._bed_surface__gsd_node
        )
        if mapped_in != "link":
            self.map_gsd_from_link_to_node(location="bedload")
        grain_size_DX_bedload = compute_grain_size_DX(
            self._sed_transport__bedload_gsd_link
            if mapped_in == "link"
            else self._sed_transport__bedload_gsd_node
        )

        grain_size_DX_surface = grain_size_DX_surface.reshape(
            grain_size_DX_surface.size, 1
        )
        grain_size_DX_bedload = grain_size_DX_bedload.reshape(
            grain_size_DX_bedload.size, 1
        )
        return grain_size_DX_surface, grain_size_DX_bedload

    def format_gsd(self, bedload_gsd):
        """Give a more friendly format for the bed surface or bed load GSD.
        Reads a bed load GSD, from links or nodes, and returns the GSD in
        cumulative percetage

        Examples
        --------

        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.07, 1.06, 1.00, 1.06, 1.07],
        ...     [1.08, 1.07, 1.03, 1.07, 1.08],
        ...     [1.09, 1.08, 1.07, 1.08, 1.09],
        ...     [1.09, 1.09, 1.08, 1.09, 1.09],
        ...     [1.09, 1.09, 1.09, 1.09, 1.09],
        ... ]

        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )

        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
        >>> grid.at_node["surface_water__velocity"] = np.full(
        ...     grid.number_of_nodes, 0.25
        ... )
        >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
        >>> grid.at_link["surface_water__velocity"] = np.full(
        ...     grid.number_of_links, 0.25
        ... )

        >>> gsd_loc = [
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ...     [0, 1.0, 1.0, 1.0, 0],
        ... ]

        >>> gsd = [[32, 100, 100], [16, 25, 50], [8, 0, 0]]

        >>> timeStep = 1  # time step in seconds

        >>> rbd = RiverBedDynamics(
        ...     grid,
        ...     gsd=gsd,
        ...     dt=timeStep,
        ...     bedload_equation="Parker1990",
        ...     bed_surface__gsd_location_node=gsd_loc,
        ... )
        >>> rbd.run_one_step()
        >>> rbd.format_gsd(rbd._sed_transport__bedload_gsd_link)
                   8      16     32
        Link_0   0.0   0.000    0.0
        Link_1   0.0  51.521  100.0
        Link_2   0.0  51.521  100.0
        Link_3   0.0   0.000    0.0
        Link_4   0.0   0.000    0.0
        Link_5   0.0   0.000    0.0
        Link_6   0.0   0.000    0.0
        Link_7   0.0   0.000    0.0
        Link_8   0.0   0.000    0.0
        Link_9   0.0   0.000    0.0
        Link_10  0.0  52.498  100.0
        Link_11  0.0  52.498  100.0
        Link_12  0.0   0.000    0.0
        Link_13  0.0  45.945  100.0
        Link_14  0.0  71.774  100.0
        Link_15  0.0  52.498  100.0
        Link_16  0.0  71.774  100.0
        Link_17  0.0  45.945  100.0
        Link_18  0.0   0.000    0.0
        Link_19  0.0  71.774  100.0
        Link_20  0.0  71.774  100.0
        Link_21  0.0   0.000    0.0
        Link_22  0.0   0.000    0.0
        Link_23  0.0  71.774  100.0
        Link_24  0.0  71.774  100.0
        Link_25  0.0  71.774  100.0
        Link_26  0.0   0.000    0.0
        Link_27  0.0   0.000    0.0
        Link_28  0.0  71.774  100.0
        Link_29  0.0  71.774  100.0
        Link_30  0.0   0.000    0.0
        Link_31  0.0   0.000    0.0
        Link_32  0.0   0.000    0.0
        Link_33  0.0   0.000    0.0
        Link_34  0.0   0.000    0.0
        Link_35  0.0   0.000    0.0
        Link_36  0.0   0.000    0.0
        Link_37  0.0   0.000    0.0
        Link_38  0.0   0.000    0.0
        Link_39  0.0   0.000    0.0

        It also works for nodes.

        >>> rbd.map_gsd_from_link_to_node(location="bedload")
        >>> rbd.format_gsd(rbd._sed_transport__bedload_gsd_node)
                   8      16     32
        Node_0   0.0   0.000    0.0
        Node_1   0.0  51.521  100.0
        Node_2   0.0  51.521  100.0
        Node_3   0.0  51.521  100.0
        Node_4   0.0   0.000    0.0
        Node_5   0.0  45.945  100.0
        Node_6   0.0  62.136  100.0
        Node_7   0.0  52.498  100.0
        Node_8   0.0  62.136  100.0
        Node_9   0.0  45.945  100.0
        Node_10  0.0  45.945  100.0
        Node_11  0.0  71.774  100.0
        Node_12  0.0  66.955  100.0
        Node_13  0.0  71.774  100.0
        Node_14  0.0  45.945  100.0
        Node_15  0.0   0.000    0.0
        Node_16  0.0  71.774  100.0
        Node_17  0.0  71.774  100.0
        Node_18  0.0  71.774  100.0
        Node_19  0.0   0.000    0.0
        Node_20  0.0   0.000    0.0
        Node_21  0.0   0.000    0.0
        Node_22  0.0   0.000    0.0
        Node_23  0.0   0.000    0.0
        Node_24  0.0   0.000    0.0
        """

        if bedload_gsd.shape[0] == self._grid.number_of_links:
            indexText = "Link_"
        else:
            indexText = "Node_"

        bedload_gsd = np.hstack((bedload_gsd, np.zeros([bedload_gsd.shape[0], 1])))
        bedload_gsd = np.cumsum(np.fliplr(bedload_gsd), axis=1) * 100
        bedload_gsd = np.around(bedload_gsd, 3)

        grain_size_D_ascending = np.sort(self._grain_size_D_original)
        columns = ["".join(item) for item in grain_size_D_ascending.astype(str)]

        df = pd.DataFrame(
            bedload_gsd,
            columns=columns,
            index=[indexText + str(i) for i in range(bedload_gsd.shape[0])],
        )
        return df

    @staticmethod
    def get_available_fields():
        """Return a list of available fields and their units.

        Examples
        --------

        >>> from pprint import pprint
        >>> from landlab.components import RiverBedDynamics

        >>> pprint(RiverBedDynamics.get_available_fields())
        [('rbd._bed_subsurface__gsd_link', '[mm,%]'),
         ('rbd._bed_subsurface__gsd_node', '[mm,%]'),
         ('rbd._bed_surface__act_layer_thick_link', '[m]'),
         ('rbd._bed_surface__act_layer_thick_node', '[m]'),
         ('rbd._bed_surface__act_layer_thick_prev_time_link', '[m]'),
         ('rbd._bed_surface__act_layer_thick_prev_time_node', '[m]'),
         ('rbd._bed_surface__elevation_fixed_node', '[m]'),
         ('rbd._bed_surface__geo_std_size_link', '[mm]'),
         ('rbd._bed_surface__geo_std_size_node', '[mm]'),
         ('rbd._bed_surface__geom_mean_size_link', '[mm]'),
         ('rbd._bed_surface__geom_mean_size_node', '[mm]'),
         ('rbd._bed_surface__gsd_fixed_node', '[mm,%]'),
         ('rbd._bed_surface__gsd_link', '[mm,%]'),
         ('rbd._bed_surface__gsd_node', '[mm,%]'),
         ('rbd._bed_surface__gsd_original_link', '[mm,%]'),
         ('rbd._bed_surface__gsd_original_node', '[mm,%]'),
         ('rbd._bed_surface__median_size_link', '[mm]'),
         ('rbd._bed_surface__median_size_node', '[mm]'),
         ('rbd._bed_surface__sand_fraction_link', '[-]'),
         ('rbd._bed_surface__sand_fraction_node', '[-]'),
         ('rbd._bed_surface__surface_thickness_new_layer_link', '[m]'),
         ('rbd._sed_transport__bedload_gsd_imposed_link', '[mm,%]'),
         ('rbd._sed_transport__bedload_gsd_link', '[mm,%]'),
         ('rbd._sed_transport__bedload_rate_imposed_link', '[m^2/s]'),
         ('rbd._sed_transport__bedload_rate_link', '[m^2/s]'),
         ('rbd._sed_transport__net_bedload_node', '[m^2/s]'),
         ('rbd._surface_water__shear_stress_link', '[Pa]'),
         ('rbd._surface_water__velocity_prev_time_link', '[m/s]'),
         ('rbd._topographic__elevation_original_link', '[m]'),
         ('rbd._topographic__elevation_original_node', '[m]'),
         ('rbd._topographic__elevation_subsurface_link', '[m]')]
        """

        # Define available fields and their units
        fields = sorted(
            [
                ("rbd._bed_subsurface__gsd_link", "[mm,%]"),
                ("rbd._bed_subsurface__gsd_node", "[mm,%]"),
                ("rbd._bed_surface__act_layer_thick_link", "[m]"),
                ("rbd._bed_surface__act_layer_thick_node", "[m]"),
                ("rbd._bed_surface__act_layer_thick_prev_time_link", "[m]"),
                ("rbd._bed_surface__act_layer_thick_prev_time_node", "[m]"),
                ("rbd._bed_surface__elevation_fixed_node", "[m]"),
                ("rbd._bed_surface__geom_mean_size_link", "[mm]"),
                ("rbd._bed_surface__geom_mean_size_node", "[mm]"),
                ("rbd._bed_surface__geo_std_size_link", "[mm]"),
                ("rbd._bed_surface__geo_std_size_node", "[mm]"),
                ("rbd._bed_surface__gsd_fixed_node", "[mm,%]"),
                ("rbd._bed_surface__gsd_link", "[mm,%]"),
                ("rbd._bed_surface__gsd_node", "[mm,%]"),
                ("rbd._bed_surface__gsd_original_link", "[mm,%]"),
                ("rbd._bed_surface__gsd_original_node", "[mm,%]"),
                ("rbd._bed_surface__median_size_link", "[mm]"),
                ("rbd._bed_surface__median_size_node", "[mm]"),
                ("rbd._bed_surface__sand_fraction_link", "[-]"),
                ("rbd._bed_surface__sand_fraction_node", "[-]"),
                ("rbd._bed_surface__surface_thickness_new_layer_link", "[m]"),
                ("rbd._sed_transport__bedload_gsd_imposed_link", "[mm,%]"),
                ("rbd._sed_transport__bedload_gsd_link", "[mm,%]"),
                ("rbd._sed_transport__bedload_rate_link", "[m^2/s]"),
                ("rbd._sed_transport__net_bedload_node", "[m^2/s]"),
                ("rbd._sed_transport__bedload_rate_imposed_link", "[m^2/s]"),
                ("rbd._surface_water__shear_stress_link", "[Pa]"),
                ("rbd._surface_water__velocity_prev_time_link", "[m/s]"),
                ("rbd._topographic__elevation_original_link", "[m]"),
                ("rbd._topographic__elevation_original_node", "[m]"),
                ("rbd._topographic__elevation_subsurface_link", "[m]"),
            ]
        )

        return fields

    @staticmethod
    def vector_mapper(grid, vector):
        """Map vector, in this case shear stress or velocity, values from
        links onto nodes preserving the components.

        This method takes the vectors values on links and determines the
        vectors components in nodes.

        Examples
        --------

        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.07, 1.06, 1.00, 1.06, 1.07],
        ...     [1.08, 1.07, 1.03, 1.07, 1.08],
        ...     [1.09, 1.08, 1.07, 1.08, 1.09],
        ...     [1.09, 1.09, 1.08, 1.09, 1.09],
        ...     [1.09, 1.09, 1.09, 1.09, 1.09],
        ... ]

        >>> grid.set_watershed_boundary_condition(
        ...     grid.at_node["topographic__elevation"]
        ... )
        >>> grid.at_node["surface_water__depth"] = np.full(grid.number_of_nodes, 0.102)
        >>> grid.at_node["surface_water__velocity"] = np.full(
        ...     grid.number_of_nodes, 0.25
        ... )
        >>> grid.at_link["surface_water__depth"] = np.full(grid.number_of_links, 0.102)
        >>> grid.at_link["surface_water__velocity"] = np.full(
        ...     grid.number_of_links, 0.25
        ... )

        >>> rbd = RiverBedDynamics(grid)
        >>> rbd.run_one_step()

        >>> (velocityVector, velocityMagnitude) = RiverBedDynamics.vector_mapper(
        ...     grid, rbd._u
        ... )
        >>> velocityVector_x = velocityVector[:, 0]
        >>> velocityVector_x.reshape(grid.shape)
        array([[ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125],
               [ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125],
               [ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125],
               [ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125],
               [ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125]])

        >>> velocityVector_y = velocityVector[:, 1]
        >>> velocityVector_y.reshape(grid.shape)
        array([[ 0.125,  0.125,  0.125,  0.125,  0.125],
               [ 0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ],
               [ 0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ],
               [ 0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ],
               [ 0.125,  0.125,  0.125,  0.125,  0.125]])
        """
        vector_x_r = vector[grid.links_at_node[:, 0]]
        (link_id,) = np.where(grid.links_at_node[:, 0] < 0)
        vector_x_r[link_id] = 0

        vector_x_l = vector[grid.links_at_node[:, 2]]
        (link_id,) = np.where(grid.links_at_node[:, 2] < 0)
        vector_x_l[link_id] = 0

        vector_y_u = vector[grid.links_at_node[:, 1]]
        (link_id,) = np.where(grid.links_at_node[:, 1] < 0)
        vector_y_u[link_id] = 0

        vector_y_l = vector[grid.links_at_node[:, 3]]
        (link_id,) = np.where(grid.links_at_node[:, 3] < 0)
        vector_y_l[link_id] = 0

        vector_x = 0.5 * (vector_x_r + vector_x_l)
        vector_y = 0.5 * (vector_y_u + vector_y_l)

        return np.transpose(np.vstack((vector_x, vector_y))), np.sqrt(
            vector_x**2 + vector_y**2
        )
