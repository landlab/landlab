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
>>> import copy
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

Let's save a copy of this topography, we will use it later.

>>> z0 = copy.deepcopy(grid.at_node["topographic__elevation"])

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
...     bed_surf__gsd_loc_node=gsd_loc,
...     track_stratigraphy=True,
... )
>>> rbd.run_one_step()

After running RiverBedDynamics, we can check if the different GSD locations
were correctly included

>>> gsd_loc_rbd = rbd._bed_surf__gsd_loc_node.reshape(grid.shape)
>>> gsd_loc_rbd
array([[0, 1, 1, 1, 0],
       [0, 1, 1, 1, 0],
       [0, 1, 1, 1, 0],
       [0, 1, 1, 1, 0],
       [0, 1, 1, 1, 0]])

Let's check at the calculated shear_stress

>>> shearStress = rbd._surface_water__shear_stress_link
>>> np.round(shearStress, decimals=3)
array([  0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,  40.011,  40.011,   0.   ,   0.   ,
        10.003,  40.011,  10.003,   0.   ,   0.   ,  10.003,  10.003,
         0.   ,   0.   ,  10.003,  10.003,  10.003,   0.   ,   0.   ,
        10.003,  10.003,   0.   ,   0.   ,   0.   ,   0.   ,   0.   ,
         0.   ,   0.   ,   0.   ,   0.   ,   0.   ])

Notice that links at borders have zero shear stress. Let's check at the calculated
net bedload. Hereinafter, most of the results will show between 3 to 6 decimals.
This is just to avoid roundoff errors problems

>>> qb = rbd._sed_transp__net_bedload_node.reshape(grid.shape)
>>> np.around(qb, decimals=6)
array([[ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.001598, -0.004793,  0.001598,  0.      ],
       [ 0.      ,  0.      ,  0.001597,  0.      ,  0.      ],
       [ 0.      ,  0.      , -0.      ,  0.      ,  0.      ],
       [ 0.      ,  0.      ,  0.      ,  0.      ,  0.      ]])

Considering the link upstream the watershed exit, link Id 15, we can obtain the
bed load transport rate

>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> formatted_number = "{:.4e}".format(qb_l15)
>>> print(formatted_number)
-1.5977e-03

Therefore, the bed load transport rate according to Parker 1990 surface-based
equation is 1.598 * 10^-3 m2/s. Negative means that is going in the negative
Y direction (towards bottom or towards south)

The GSD at this place is:

>>> qb_gsd_l15 = rbd._sed_transp__bedload_gsd_link[15]
>>> np.round(qb_gsd_l15, decimals=3)
array([0.475, 0.525])

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

>>> qb_gsd_l = rbd._sed_transp__bedload_gsd_link
>>> np.round(qb_gsd_l, decimals=3)
array([[0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.475, 0.525],
       [0.475, 0.525],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.282, 0.718],
       [0.475, 0.525],
       [0.282, 0.718],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.282, 0.718],
       [0.282, 0.718],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.282, 0.718],
       [0.282, 0.718],
       [0.282, 0.718],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.282, 0.718],
       [0.282, 0.718],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ],
       [0.   , 0.   ]])

Zeros indicate that there is no sediment transport of that grain size
at that location.

After the flow acted on the bed and sediment transport occured we can check
the new topographic elevation field

>>> z = grid.at_node["topographic__elevation"].reshape(grid.shape)
>>> np.round(z, decimals=3)
array([[1.07 , 1.06 , 1.007, 1.06 , 1.07 ],
       [1.08 , 1.068, 1.037, 1.068, 1.08 ],
       [1.09 , 1.08 , 1.068, 1.08 , 1.09 ],
       [1.09 , 1.09 , 1.08 , 1.09 , 1.09 ],
       [1.09 , 1.09 , 1.09 , 1.09 , 1.09 ]])

Let's take a look at bed load transport rate when we use the different bedload equations.
First, let's recover the original topography.

>>> grid.at_node["topographic__elevation"] = z0.copy()

For the defaul MPM model we get:

>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="MPM",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-2.2970e-03

For Fernandez Luque and Van Beek:

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="FLvB",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-1.6825e-03

For Wong and Parker:

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WongAndParker",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-1.1326e-03

For Huang:

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="Huang",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-1.1880e-03

For Wilcock and Crowe 2003:

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(
...     grid,
...     gsd=gsd,
...     bedload_equation="WilcockAndCrowe",
...     bed_surf__gsd_loc_node=gsd_loc,
... )
>>> rbd.run_one_step()
>>> qb_l15 = rbd._sed_transp__bedload_rate_link[15]
>>> print("{:.4e}".format(qb_l15))
-5.3081e-04

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
>>> sandContent = rbd._bed_surf__sand_fract_node[20]
>>> print(float("{:.4f}".format(round(sandContent, 4))))
0.1

But it is different if we use Parker 1990, because it removes sand content

>>> rbd = RiverBedDynamics(grid, gsd=gsd, bedload_equation="Parker1990")
>>> rbd._bed_surf__sand_fract_node[20]
0.0

What happens if we give it a set of wrong optional fields. The following
fields will have only two elements, which is different than the number of
nodes and links

>>> qb_imposed = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, sed_transp__bedload_rate_fix_link=qb_imposed)
>>> rbd._sed_transp__bedload_rate_fix_link
array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
       0., 0., 0., 0., 0., 0.])

>>> gsd_loc = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surf__gsd_loc_node=gsd_loc)
>>> rbd._bed_surf__gsd_loc_node.reshape(grid.shape)
array([[0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0]])

>>> gsd_fix = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surf__gsd_fix_node=gsd_fix)
>>> rbd._bed_surf__gsd_fix_node.reshape(grid.shape)
array([[0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0]])

>>> elev_fix = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, bed_surf__elev_fix_node=elev_fix)
>>> rbd._bed_surf__elev_fix_node.reshape(grid.shape)
array([[0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0],
       [0, 0, 0, 0, 0]])

>>> vel_n_1 = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, surface_water__velocity_prev_time_link=vel_n_1)
>>> rbd._surface_water__velocity_prev_time_link
array([0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
       0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
       0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
       0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25])

For sed_transp__bedload_gsd_fix_link, for simplicity, we only show
links 0, 5, and 10

>>> qb_gsd_imposed = np.array([1, 2])
>>> rbd = RiverBedDynamics(grid, sed_transp__bedload_gsd_fix_link=qb_gsd_imposed)
>>> rbd._sed_transp__bedload_gsd_fix_link[[0, 5, 10], :]
array([[0., 0.],
       [0., 0.],
       [0., 0.]])

In summary, in all these cases the wrong given value is override by default values.
But, if the size of the array is correct the specified condition is used.
For simplicity, we only show links 0, 5, and 10

>>> qb_imposed = np.full(grid.number_of_links, 1)
>>> rbd = RiverBedDynamics(grid, sed_transp__bedload_rate_fix_link=qb_imposed)
>>> rbd._sed_transp__bedload_rate_fix_link[[0, 5, 10]]
array([1., 1., 1.])

>>> gsd_loc = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surf__gsd_loc_node=gsd_loc)
>>> rbd._bed_surf__gsd_loc_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> gsd_fix = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surf__gsd_fix_node=gsd_fix)
>>> rbd._bed_surf__gsd_fix_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> elev_fix = np.full(grid.number_of_nodes, 1)
>>> rbd = RiverBedDynamics(grid, bed_surf__elev_fix_node=elev_fix)
>>> rbd._bed_surf__elev_fix_node.reshape(grid.shape)
array([[1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1],
       [1, 1, 1, 1, 1]])

>>> vel_n_1 = np.full(grid.number_of_links, 1)
>>> rbd = RiverBedDynamics(grid, surface_water__velocity_prev_time_link=vel_n_1)
>>> rbd._surface_water__velocity_prev_time_link[[0, 5, 10]]
array([1., 1., 1.])

>>> qb_gsd_imposed = np.ones((grid.number_of_links, 2))  # 2 comes from gsd.shape[0]-1
>>> rbd = RiverBedDynamics(grid, sed_transp__bedload_gsd_fix_link=qb_gsd_imposed)
>>> rbd._sed_transp__bedload_gsd_fix_link[[0, 5, 10]]
array([[1., 1.],
       [1., 1.],
       [1., 1.]])

Using the hydraulics radius is also possible. Let's compare the shear stress with and
without that option. First, without including the hydraulics radius.

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(grid)
>>> rbd.run_one_step()
>>> print(np.around(rbd._shear_stress[15], decimals=3))
-40.011

Now, we will consider the hydraulics radius

>>> grid.at_node["topographic__elevation"] = z0.copy()
>>> rbd = RiverBedDynamics(grid, use_hydraulics_radius_in_shear_stress=True)
>>> rbd.run_one_step()
>>> print(np.around(rbd._shear_stress[15], decimals=3))
-33.232

So, there is an important difference between the two ways of calculating it.

"""

import numpy as np
import scipy.constants

from landlab import Component

from . import _bedload_eq_MPM_style as MPM_style
from . import _bedload_eq_Parker_1990 as Parker1990
from . import _bedload_eq_Wilcock_Crowe_2003 as WilcockAndCrowe2003
from . import _initialize_fields as initialize
from . import _initialize_gsd as initialize_gsd
from . import _nodes_and_links_info as info
from . import _stratigraphy as stratigraphy
from . import _utilities as utilities


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
        dt=1,
        alpha=1.0,
        bed_surf__elev_fix_node=None,
        bed_surf__gsd_fix_node=None,
        sed_transp__bedload_rate_fix_link=None,
        sed_transp__bedload_gsd_fix_link=None,
        bed_surf__gsd_loc_node=None,
        surface_water__velocity_prev_time_link=None,
        current_t=0.0,
        track_stratigraphy=False,
        num_cycles_to_process_strat=10,
        bed_surf_new_layer_thick=1,
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
            * ``'WongAndParker'`` for Wong and Parker 2006.
            * ``'Huang'`` for Huang 2010.
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
        dt: float, optional
            Time step in seconds. When this component is coupled to a flow model,
            it is dynamically updated.
        alpha : float, optional
            An upwinding coefficient for a central difference scheme when
            updating the bed GSD - default value is 1.0 - a value of 0.5
            generates a central differences scheme.
        bed_surf__elev_fix_node: ndarray of int, optional
            Sets a node as a fixed elevation, this means that it's elevation does
            not change during the simulation. Use 1 to denote a fix node and 0 for
            nodes that can change elevation
            Units: - , mapping: node
        sed_transp__bedload_rate_fix_link: ndarray of float, optional
            Sets the sediment transport rate per unit width as fixed at a given link.
            This means that its value will not change unless is manually redefined.
            When defining it, use 0 for non-fixed and other float numbers to specify the
            bed load rate. Specifing actually a zero bed load transport is not supported,
            but a very small value could be used instead.
            Units: m^2/s, mapping: link
        bed_surf__gsd_fix_node: ndarray of int, optional
            Sets a node as fixed gsd, this means that it's gsd does not change during
            the simulation.
            Units: - , mapping: node
        bed_surf__gsd_loc_node: ndarray of int, optional
            Sets the location at each node in which the GSD applies
            Units: - , mapping: node
        sed_transp__bedload_gsd_fix_link: ndarray of float, optional
            Sets the sediment transport GSD where sediment supply is imposed. It's
            size should be columns: number of links; rows: Number of grain sizes - 1
            Units: - , mapping: link
        surface_water__velocity_prev_time_link: ndarray of float, optional
            Speed of water flow above the surface in the previous time step
            Units: m/s, mapping: link
        current_t : float, optional
            Current simulation time or elapsed time. It does not update automatically
            Units: s
        track_stratigraphy : bool, optional
            If ``True``, the component stores the GSD of each layer at every node.
            This is computationally demanding as it needs to read/write data
            at every time step. Recommended only when cycles of erosion and
            deposition are frequent and very important.
        num_cycles_to_process_strat : int, optional
            If ``track_stratigraphy`` is ``True``, data will be read and stored every
            ``'num_cycles_to_process_strat'`` time steps. Must be larger
            or equal to 1. All data are written to the hard drive. It does not use
            *Landlab* layers.
        bed_surf_new_layer_thick : float, optional
            When this thickness is reached in a deposition or erosion process, meaning
            that the new bed at a node is "bed_surf_new_layer_thick" higher or
            deeper a new layer is created, in which the time changes in gsd are
            consolidated and becomes a single value.
        """
        super().__init__(grid)

        self._g = scipy.constants.g  # Acceleration due to gravity (m/s**2).
        self._rho = rho
        self._rho_s = rho_s
        self._R = (rho_s - rho) / rho

        self._normal = -self._grid.link_dirs_at_node  # Define faces normal vector

        if gsd is None:
            self._gsd = np.array([[32, 100], [16, 25], [8, 0]])
        else:
            self._gsd = np.array(gsd)

        self._bedload_equation = bedload_equation

        self._variable_critical_shear_stress = variable_critical_shear_stress
        self._use_hydraulics_radius_in_shear_stress = (
            use_hydraulics_radius_in_shear_stress
        )
        self._lambda_p = lambda_p
        self._alpha = alpha
        self._outlet_boundary_condition = outlet_boundary_condition

        self._grid._dt = dt
        self._current_t = current_t

        # Initialize required fields and identify nodes and links with
        # fixed surface elevation and/or gsd
        self._bed_surf__gsd_loc_node = initialize.field_at_node(
            grid, bed_surf__gsd_loc_node
        )

        self._bed_surf__elev_fix_node = initialize.field_at_node(
            grid, bed_surf__elev_fix_node
        ).astype(int)
        self._bed_surf__elev_fix_node_id = np.where(
            (self._bed_surf__elev_fix_node) == 1
        )

        self._bed_surf__elev_fix_link_id = np.unique(
            self._grid.links_at_node[self._bed_surf__elev_fix_node_id]
        )

        self._bed_surf__gsd_fix_node = initialize.field_at_node(
            grid, bed_surf__gsd_fix_node
        ).astype(int)
        self._bed_surf__gsd_fix_node_id = np.where(self._bed_surf__gsd_fix_node == 1)

        # Now self._bed_surf__gsd_fix_node changes from nodes to links
        bed_surf__gsd_fix_link = np.unique(
            grid.links_at_node[self._bed_surf__gsd_fix_node_id]
        )
        self._bed_surf__gsd_fix_link = bed_surf__gsd_fix_link[
            np.where(bed_surf__gsd_fix_link >= 0)
        ]

        self._surface_water__velocity_prev_time_link = initialize.velocity_at_link(
            grid, surface_water__velocity_prev_time_link
        )
        self._sed_transp__bedload_rate_link = self._grid.zeros(
            at="link"
        )  # Volumetric bed load transport rate per unit width
        self._sed_transp__net_bedload_node = self._grid.zeros(at="node")
        self._sed_transp__bedload_gsd_link = self._grid.zeros(at="link")
        self._surface_water__shear_stress_link = self._grid.zeros(at="link")

        # Initialize the bed surface grain size properties using inputs
        self.define_initial_bed_properties()

        self._sed_transp__bedload_rate_fix_link = initialize.field_at_link(
            grid, sed_transp__bedload_rate_fix_link
        )
        self._sed_transp__bedload_rate_fix_link_id = np.where(
            self._sed_transp__bedload_rate_fix_link != 0
        )[0]
        self._sed_transp__bedload_gsd_fix_link = initialize.gsd_at_link(
            self._grid, sed_transp__bedload_gsd_fix_link, self._gsd
        )
        self._sed_transp__bedload_gsd_fix_link_id = info.fixed_links(
            self._sed_transp__bedload_gsd_fix_link
        )

        # Identify the node upstream of the outlet
        # Used for boundary conditions when updating bed elevation and bed gsd
        (
            self._out_id,
            self._upstream_out_id,
            self._outlet_links,
            self._closed_nodes,
            self._boundary_links,
        ) = info.outlet_nodes(self._grid)

        # This flag is used to activate or deactivate the bed GSD updating part
        # of the component.
        self._update_bed_surf_GSD = False

        # Activates option to store the GSD of individual layers in each node.
        self._track_stratigraphy = track_stratigraphy

        # Threshold to deposit layers in a new subsurface layer
        self._bed_surf_new_layer_thick = bed_surf_new_layer_thick

        # When bed_surf_new_layer_thick is reached this flag is used to
        self._update_stratigraphy = False  # record and read the data
        self._update_subsurface_deposited = (
            False  # update subsurface data when there is deposition
        )
        self._update_subsurface_eroded = (
            False  # update subsurface data when there is erosion
        )

        # Sets initial values to enter into the write and read stratigraphy routine
        self._stratigraphy_cycle = 0
        self._num_cycles_to_process_strat = num_cycles_to_process_strat

        # Makes a copy of the original bed surface elevation and maps into links
        self._grid["link"]["topographic__elevation"] = (
            self._grid.map_mean_of_link_nodes_to_link(
                self._grid["node"]["topographic__elevation"]
            )
        )
        self._topogr__elev_orig_node = self._grid["node"][
            "topographic__elevation"
        ].copy()
        self._topogr__elev_orig_link = self._grid.map_mean_of_link_nodes_to_link(
            self._topogr__elev_orig_node
        )
        self._topogr__elev_subsurf_link = self._topogr__elev_orig_link.copy()
        self._bed_surf__thick_new_layer_link = np.zeros_like(
            self._topogr__elev_orig_link
        )

        # Check that parker 1990 or Wilcock and Crowe were selected when
        # tracking stratigraphy
        stratigraphy.checks_correct_equation_to_track_stratigraphy(self)

        # Creates links dictionary to store changes in stratigraphy
        stratigraphy.create_links_dictionary(self)

    def define_initial_bed_properties(self):
        """This method performs the initial setup of the bed grain size distribution properties.
        It reads the input data and populates the necessary variables.

        This configuration is only performed during the first time step. Subsequent time steps
        will utilize the bed information, which has already been calculated or updated and
        formatted appropriately.
        """
        # Adds the 2mm fraction
        self._gsd = initialize_gsd.adds_2mm_to_gsd(self._gsd)

        # Removes sand fractions in case Parker 1999 is selected
        self._gsd = initialize_gsd.remove_sand_from_gsd(
            self._gsd, self._bedload_equation
        )

        # Maps inputs to nodes - self._gs contains the equivalent grain sizes D_eq
        (
            sand_fraction,
            gs_D_equiv_freq,
            self._gs,
        ) = initialize_gsd.map_initial_bed_properties_to_nodes(
            self._gsd, self._bed_surf__gsd_loc_node
        )

        median_size_D50 = self.calculate_DX(
            gs_D_equiv_freq, 0.5
        )  # Median grain size in each node
        (gs_geom_mean, gs_geo_std) = self.calculate_gsd_geo_mean_and_geo_std(
            gs_D_equiv_freq
        )

        # Bed grain sizes frequency in each node
        self._bed_surf__gsd_node = gs_D_equiv_freq
        self._bed_surf__gsd_orig_node = gs_D_equiv_freq.copy()
        self._bed_subsurf__gsd_node = gs_D_equiv_freq.copy()
        self._bed_surf__median_size_node = median_size_D50
        self._bed_surf__geom_mean_size_node = gs_geom_mean
        self._bed_surf__geo_std_size_node = gs_geo_std
        self._bed_surf__sand_fract_node = sand_fraction

        def map_mean_of_nodes_to_link(var, r_node, l_node):
            return 0.5 * (var[r_node] + var[l_node])

        # GSD properties is mapped from nodes onto links
        r_node = self.grid.nodes_at_link[:, 0]
        l_node = self.grid.nodes_at_link[:, 1]

        self._bed_surf__gsd_link = map_mean_of_nodes_to_link(
            gs_D_equiv_freq, r_node, l_node
        )
        self._bed_surf__median_size_link = map_mean_of_nodes_to_link(
            median_size_D50, r_node, l_node
        )
        self._bed_surf__geom_mean_size_link = map_mean_of_nodes_to_link(
            gs_geom_mean, r_node, l_node
        )
        self._bed_surf__geo_std_size_link = map_mean_of_nodes_to_link(
            gs_geo_std, r_node, l_node
        )
        self._bed_surf__sand_fract_link = map_mean_of_nodes_to_link(
            sand_fraction, r_node, l_node
        )

        self._bed_surf__gsd_orig_link = self._bed_surf__gsd_link.copy()
        self._bed_subsurf__gsd_link = self._bed_surf__gsd_link.copy()

        self._bed_surf__act_layer_thick_link = (
            2 * self.calculate_DX(self._bed_surf__gsd_link, 0.9) / 1000
        )
        self._bed_surf__act_layer_thick_prev_time_link = (
            self._bed_surf__act_layer_thick_link.copy()
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

        Changes grid topography. Starts at self.update_bed_elevation()

        For one time step, this erodes the grid topography according to
        Exner equation::

            (1-λp) ∂Z/∂t = - (∂qbx/∂x + ∂qby/∂y)

        Simplifying, we get::

            ∂Z/∂t = - (1 / (1-λp)) * (∂Qb/∂A)
            Z_t+1 = -(Δt * ΔQb)/(1-λp) + Z_t

        The grid field ``"topographic__elevation"`` is altered each time step.
        """
        self.shear_stress()  # Shear stress calculation
        self.bedload_equation()  # Bedload calculation
        self.calculate_net_bedload()  # Calculates bedload transport from links into nodes
        self.update_bed_elevation()  # Changes bed elevation
        stratigraphy.checks_erosion_or_deposition(self)
        stratigraphy.evolve(self)
        self.update_bed_surf_gsd()  # Changes bed surface grain size distribution
        self.update_bed_surf_properties()  # Updates gsd properties

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
        self._shear_stress[self._boundary_links] = 0

        # Direction of flux will be recovered at the end of the bedload
        # transport routine
        self._surface_water__shear_stress_link = np.abs(self._shear_stress)

    def bedload_equation(self):
        if self._bedload_equation == "Parker1990":
            (
                self._sed_transp__bedload_rate_link,
                self._sed_transp__bedload_gsd_link,
            ) = Parker1990.bedload_equation(self)
        elif self._bedload_equation == "WilcockAndCrowe":
            (
                self._sed_transp__bedload_rate_link,
                self._sed_transp__bedload_gsd_link,
            ) = WilcockAndCrowe2003.bedload_equation(self)
        else:
            self._sed_transp__bedload_rate_link = MPM_style.bedload_equation(self)

    def calculate_net_bedload(self):
        """Calculates the net volumetric bedload coming from all links (m2/s)
        onto nodes (m3/s).

        This method takes the volumetric bedload entering and exiting through a
        face and determines the net volumetric bedload on a given node.
        """

        qb_x = (
            np.sum(
                self._sed_transp__bedload_rate_link[self._grid.links_at_node[:, [0, 2]]]
                * self._normal[:, [0, 2]],
                axis=1,
            )
            * self._grid.dy
        )
        qb_y = (
            np.sum(
                self._sed_transp__bedload_rate_link[self._grid.links_at_node[:, [1, 3]]]
                * self._normal[:, [1, 3]],
                axis=1,
            )
            * self._grid.dx
        )

        self._sed_transp__net_bedload_node = qb_x + qb_y

        ## At the boundary, there is no exiting link, so we assume a zero flux
        # exiting. This assumption is overridden in the Exner equation, where a
        # zero gradient boundary condition can be used.
        self._sed_transp__net_bedload_node[self._grid.boundary_nodes] = 0

    def update_bed_elevation(self):
        """Applies the Exner equation and boundary conditions to predict
        the change in bed surface elevation.

        ∂Z/∂t = - (1 / (1-λp)) * (∂Qb/∂A)
        Z_t+1 = -(Δt * ΔQb/ΔA)/(1-λp) + Z_t

        """
        # Variable definitions
        area = self._grid.dx * self._grid.dy
        d_qb = self._sed_transp__net_bedload_node  # Represents ΔQb
        z = self._grid["node"]["topographic__elevation"]
        z0 = self._grid["node"]["topographic__elevation"].copy()

        # Applies Exner Equation
        dz = -self._grid._dt / ((1 - self._lambda_p) * area) * d_qb
        z += dz

        # Boundary conditions
        # Return outlet, closed, and fixed nodes to the previous condition
        z[self._out_id] = z0[self._out_id]
        z[self._bed_surf__elev_fix_node_id] = z0[self._bed_surf__elev_fix_node_id]
        z[self._closed_nodes] = z0[self._closed_nodes]

        # Now we can apply the boundary conditions
        if self._outlet_boundary_condition == "zeroGradient":
            dz_outlet = dz[self._upstream_out_id]
        elif self._outlet_boundary_condition == "fixedValue":
            dz_outlet = 0

        z[self._out_id] = z0[self._out_id] + dz_outlet

        # Now we map data into links to update Bed GSD
        self._grid["link"]["topographic__elevation"] = (
            self._grid.map_mean_of_link_nodes_to_link(z)
        )
        self._grid["node"]["topographic__elevation"] = z

    def update_bed_surf_gsd(self):
        """Uses the fractional Exner equation to update the bed GSD"""
        # Here we create a number of variables that will be used in the
        # following definitions to make the code a little bit cleaner

        if self._track_stratigraphy:
            n_links = self._grid.number_of_links
            n_cols = self._grid.number_of_node_columns
            gsd_F = self._bed_surf__gsd_link
            gsd_Fs = self._bed_subsurf__gsd_link
            pl = self._sed_transp__bedload_gsd_link
            qbT = self._sed_transp__bedload_rate_link
            la = np.reshape(self._bed_surf__act_layer_thick_link, [n_links, 1])
            la0 = np.reshape(
                self._bed_surf__act_layer_thick_prev_time_link,
                [n_links, 1],
            )
            z = self._grid["link"]["topographic__elevation"]
            z0 = self._topogr__elev_orig_link

            lps = self._lambda_p
            dx = self._grid.dx
            dy = self._grid.dy
            alpha = self._alpha
            dt = self._grid._dt

            dv = 2 * n_cols - 1

            qbT = np.reshape(qbT, [n_links, 1])
            qbTdev = np.zeros([n_links, 1])
            qjj1dev = np.zeros_like(pl)

            # Horizontal Links
            # Border links to apply 1st order approximation
            hlL = self._grid.links_at_node[self._grid.nodes_at_left_edge][:, 0]
            hlR = self._grid.links_at_node[self._grid.nodes_at_right_edge][:, 2]
            # rest of horizontal links to apply 1st or second order approximation
            hl_Id = np.isin(self._grid.horizontal_links, np.hstack((hlL, hlR)))
            hl = self._grid.horizontal_links[~hl_Id]

            # Determine flow direction at each horizontal link location
            hl_pos = hl[np.where(qbT[hl][:, 0] >= 0)]
            hl_neg = hl[np.where(qbT[hl][:, 0] < 0)]
            hlL_pos = hlL[np.where(qbT[hlL][:, 0] >= 0)]
            hlL_neg = hlL[np.where(qbT[hlL][:, 0] < 0)]
            hlR_pos = hlR[np.where(qbT[hlR][:, 0] >= 0)]
            hlR_neg = hlR[np.where(qbT[hlR][:, 0] < 0)]

            qbTdev[hl_pos] = (
                alpha * (qbT[hl_pos] - qbT[hl_pos - 1]) / dy
                + (1 - alpha) * (qbT[hl_pos + 1] - qbT[hl_pos]) / dy
            )
            qbTdev[hl_neg] = (
                alpha * (qbT[hl_neg] - qbT[hl_neg + 1]) / dy
                + (1 - alpha) * (qbT[hl_neg - 1] - qbT[hl_neg]) / dy
            )
            qbTdev[hlL_pos] = (qbT[hlL_pos + 1] - qbT[hlL_pos]) / dy
            qbTdev[hlL_neg] = (qbT[hlL_neg] - qbT[hlL_neg - 1]) / dy
            qbTdev[hlR_pos] = (qbT[hlR_pos] - qbT[hlR_pos - 1]) / dy
            qbTdev[hlR_neg] = (qbT[hlR_neg - 1] - qbT[hlR_neg]) / dy

            qjj1dev[hl_pos, :] = (
                alpha
                * (qbT[hl_pos] * pl[hl_pos, :] - qbT[hl_pos - 1] * pl[hl_pos - 1, :])
                / dy
                + (1 - alpha)
                * (qbT[hl_pos + 1] * pl[hl_pos + 1] - qbT[hl_pos] * pl[hl_pos])
                / dy
            )
            qjj1dev[hl_neg, :] = (
                alpha
                * (qbT[hl_neg] * pl[hl_neg, :] - qbT[hl_neg + 1] * pl[hl_neg + 1, :])
                / dy
                + (1 - alpha)
                * (qbT[hl_neg - 1] * pl[hl_neg - 1] - qbT[hl_neg] * pl[hl_neg])
                / dy
            )
            qjj1dev[hlL_pos, :] = (
                qbT[hlL_pos + 1] * pl[hlL_pos + 1, :] - qbT[hlL_pos] * pl[hlL_pos, :]
            ) / dy
            qjj1dev[hlL_neg, :] = (
                qbT[hlL_neg] * pl[hlL_neg, :] - qbT[hlL_neg - 1] * pl[hlL_neg - 1, :]
            ) / dy
            qjj1dev[hlR_pos, :] = (
                qbT[hlR_pos] * pl[hlR_pos, :] - qbT[hlR_pos - 1] * pl[hlR_pos - 1, :]
            ) / dy
            qjj1dev[hlR_neg, :] = (
                qbT[hlR_neg - 1] * pl[hlR_neg - 1, :] - qbT[hlR_neg] * pl[hlR_neg, :]
            ) / dy

            # Vertical Links
            # Border links to apply 1st order approximation
            vlB = self._grid.links_at_node[self._grid.nodes_at_bottom_edge][:, 1]
            vlT = self._grid.links_at_node[self._grid.nodes_at_top_edge][:, 3]
            # rest of vertical links to apply 1st or second order approximation
            vl_Id = np.isin(self._grid.vertical_links, np.hstack((vlB, vlT)))
            vl = self._grid.vertical_links[~vl_Id]

            vl_pos = vl[np.where(qbT[vl][:, 0] >= 0)]
            vl_neg = vl[np.where(qbT[vl][:, 0] < 0)]
            vlB_pos = vlB[np.where(qbT[vlB][:, 0] >= 0)]
            vlB_neg = vlB[np.where(qbT[vlB][:, 0] < 0)]
            vlT_pos = vlT[np.where(qbT[vlT][:, 0] >= 0)]
            vlT_neg = vlT[np.where(qbT[vlT][:, 0] < 0)]

            qbTdev[vl_pos] = (
                alpha * (qbT[vl_pos] - qbT[vl_pos - dv]) / dx
                + (1 - alpha) * (qbT[vl_pos + dv] - qbT[vl_pos]) / dx
            )
            qbTdev[vl_neg] = (
                alpha * (qbT[vl_neg] - qbT[vl_neg + dv]) / dx
                + (1 - alpha) * (qbT[vl_neg - dv] - qbT[vl_neg]) / dx
            )
            qbTdev[vlB_pos] = (qbT[vlB_pos + dv] - qbT[vlB_pos]) / dx
            qbTdev[vlB_neg] = (qbT[vlB_neg] - qbT[vlB_neg + dv]) / dx
            qbTdev[vlT_pos] = (qbT[vlT_pos] - qbT[vlT_pos - dv]) / dx
            qbTdev[vlT_neg] = (qbT[vlT_neg - dv] - qbT[vlT_neg]) / dx

            # Very small rounding errors, for example -1e-20 is considered negative
            # but in fact is zero. That affects gsd_FIexc
            qbTdev = np.round(qbTdev, 8)

            gsd_FIexc = gsd_Fs.copy()
            (id0,) = np.where(qbTdev[:, 0] <= 0)
            gsd_FIexc[id0, :] = 0.7 * gsd_F[id0, :] + 0.3 * pl[id0, :]

            qjj1dev[vl_pos, :] = (
                alpha
                * (qbT[vl_pos] * pl[vl_pos, :] - qbT[vl_pos - dv] * pl[vl_pos - dv, :])
                / dx
                + (1 - alpha)
                * (qbT[vl_pos + dv] * pl[vl_pos + dv] - qbT[vl_pos] * pl[vl_pos])
                / dx
            )
            qjj1dev[vl_neg, :] = (
                alpha
                * (qbT[vl_neg] * pl[vl_neg, :] - qbT[vl_neg + dv] * pl[vl_neg + dv, :])
                / dx
                + (1 - alpha)
                * (qbT[vl_neg - dv] * pl[vl_neg - dv] - qbT[vl_neg] * pl[vl_neg])
                / dx
            )
            qjj1dev[vlB_pos, :] = (
                qbT[vlB_pos + dv] * pl[vlB_pos + dv, :] - qbT[vlB_pos] * pl[vlB_pos, :]
            ) / dx
            qjj1dev[vlB_neg, :] = (
                qbT[vlB_neg] * pl[vlB_neg, :] - qbT[vlB_neg + dv] * pl[vlB_neg + dv, :]
            ) / dx
            qjj1dev[vlT_pos, :] = (
                qbT[vlT_pos] * pl[vlT_pos, :] - qbT[vlT_pos - dv] * pl[vlT_pos - dv, :]
            ) / dx
            qjj1dev[vlT_neg, :] = (
                qbT[vlT_neg - dv] * pl[vlT_neg - dv, :] - qbT[vlT_neg] * pl[vlT_neg, :]
            ) / dx

            qjj2dev = gsd_FIexc * np.reshape(qbTdev, [n_links, 1])
            gsd_Fnew = gsd_F + dt * (-qjj1dev + qjj2dev) / (1 - lps) / la

            # skips the very first time step
            if self._current_t > 0:
                gsd_Fnew += (gsd_FIexc - gsd_F) / la * (la - la0)

            # Renormalize gsd
            gsd_Fnew[np.where(gsd_Fnew <= 0)] = 0  # In case a negative value appears
            gsd_Fnew = gsd_Fnew / np.reshape(np.sum(gsd_Fnew, axis=1), [n_links, 1])
            gsd_Fnew = np.nan_to_num(gsd_Fnew)

            # Boundary conditions
            # Links around the outlet and upstream the outlet are reverse to their
            # original state. This is specially important if OverlandFlow is used,
            # because the sudden drop in water depth.

            gsd_Fnew[self._outlet_links] = gsd_F[self._outlet_links]

            # Revert any changes to fixed GSD and closed links
            gsd_Fnew[self._bed_surf__gsd_fix_link] = gsd_F[self._bed_surf__gsd_fix_link]
            gsd_Fnew[self._boundary_links] = gsd_F[self._boundary_links]

            # Restores the gsd in links that previously had deposition and now are
            # eroded below the original elevation
            id_eroded_links = np.where(z < z0)[0]

            if id_eroded_links.size > 0:
                gsd_Fnew[id_eroded_links] = gsd_F[id_eroded_links]
                # Updates the "original" bed elevation to account for repetitive
                # deposition/erosion
                z0[id_eroded_links] = z[id_eroded_links]

            # Now we update the bed surface GSD and original bed elevation
            self._topogr__elev_orig_link = z0.copy()
            self._bed_surf__gsd_link = gsd_Fnew.copy()
            self._bed_surf__gsd_node = utilities.map_gsd_from_link_to_node(self)

    def update_bed_surf_properties(self):
        """Calculates the updated GSD properties"""

        self._bed_surf__median_size_link = self.calculate_DX(
            self._bed_surf__gsd_link, 0.5
        )  # Median grain size in each node
        (
            self._bed_surf__geom_mean_size_link,
            self._bed_surf__geo_std_size_link,
        ) = self.calculate_gsd_geo_mean_and_geo_std(self._bed_surf__gsd_link)
        self._bed_surf__act_layer_thick_prev_time_link = (
            self._bed_surf__act_layer_thick_link.copy()
        )
        self._bed_surf__act_layer_thick_link = (
            2 * self.calculate_DX(self._bed_surf__gsd_link, 0.9) / 1000
        )

        self._bed_surf__median_size_node = self.calculate_DX(
            self._bed_surf__gsd_node, 0.5
        )
        (
            self._bed_surf__geom_mean_size_node,
            self._bed_surf__geo_std_size_node,
        ) = self.calculate_gsd_geo_mean_and_geo_std(self._bed_surf__gsd_node)

    def calculate_DX(self, gs_D_equiv_freq, fX):
        """Calculate the grain size corresponding to any fraction.
        For example, 50%, which is the median_size_D50. In that case use fX = 0.5

        This method takes the user specified fraction, from 0 to 1, and outputs
        the corresponding grain size in nodes or links, which is considered
        by taking the size of the passed argument"""

        gs_D = self._gsd[:, 0]
        if gs_D_equiv_freq.shape[0] == self.grid.number_of_links:
            freq_gs_list = np.arange(self.grid.number_of_links)
        else:
            freq_gs_list = np.arange(self.grid.number_of_nodes)

        gs_Psi_scale_D = np.flip(np.log2(gs_D), axis=0)  # Psi scale

        gs_freq = np.hstack((gs_D_equiv_freq, np.zeros([gs_D_equiv_freq.shape[0], 1])))
        gs_D_equiv_freq_cum = np.cumsum(np.flip(gs_freq, axis=1), axis=1)
        i0 = np.argmin(gs_D_equiv_freq_cum <= fX, axis=1) - 1
        i1 = np.argmax(gs_D_equiv_freq_cum > fX, axis=1)

        gs_Psi_scale_DX = gs_Psi_scale_D[i0] + (
            (gs_Psi_scale_D[i1] - gs_Psi_scale_D[i0])
            / (
                gs_D_equiv_freq_cum[freq_gs_list, i1]
                - gs_D_equiv_freq_cum[freq_gs_list, i0]
            )
        ) * (fX - gs_D_equiv_freq_cum[freq_gs_list, i0])

        gs_DX = 2**gs_Psi_scale_DX

        return gs_DX

    def calculate_gsd_geo_mean_and_geo_std(self, gs_D_equiv_freq):
        """Calculates the geometric mean and standard deviation in links or nodes
        depending on the input"""

        gs_D_eq_Psi = np.log2(self._gs)
        gs_D_eq_Psi_mean = np.sum(gs_D_equiv_freq * gs_D_eq_Psi, axis=1)
        gs_geom_mean = 2**gs_D_eq_Psi_mean
        # Equivalent grain sizes in Psi scale

        if gs_D_equiv_freq.shape[0] == self.grid.number_of_links:
            n = self.grid.number_of_links
        else:
            n = self.grid.number_of_nodes

        gs_D_eq_Psi = np.tile(gs_D_eq_Psi, (n, 1))
        gs_D_eq_Psi_mean = np.reshape(gs_D_eq_Psi_mean, [n, 1])
        gs_geo_std = 2 ** np.sqrt(
            np.sum(
                ((gs_D_eq_Psi - gs_D_eq_Psi_mean) ** 2) * gs_D_equiv_freq,
                axis=1,
            )
        )

        return gs_geom_mean, gs_geo_std

    def stratigraphy_write_evolution(self):
        """Writes the stratigraphy time evolution into a csv file

        The example below is a shorter version of the one used in _stratigraphy.py

        Examples
        --------

        As per usual, we define import the required libraries and create a grid
        and configure the mandatory fields.

        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> import os

        >>> grid = RasterModelGrid((8, 3), xy_spacing=100)

        >>> grid.at_node["topographic__elevation"] = [
        ...     [1.12, 1.00, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.01, 1.12],
        ...     [1.12, 1.12, 1.12],
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

        >>> fixed_nodes = np.zeros(grid.number_of_nodes)
        >>> fixed_nodes[[1, 4]] = 1

        >>> fixed_bed_gsd_nodes = np.zeros(grid.number_of_nodes)
        >>> fixed_bed_gsd_nodes[[1, 4]] = 1

        >>> qb = np.full(grid.number_of_links, 0.0)
        >>> qb[[28, 33]] = -0.002

        >>> rbd = RiverBedDynamics(
        ...     grid,
        ...     gsd=gsd,
        ...     bedload_equation="Parker1990",
        ...     outlet_boundary_condition="fixedValue",
        ...     bed_surf__elev_fix_node=fixed_nodes,
        ...     bed_surf__gsd_fix_node=fixed_bed_gsd_nodes,
        ...     sed_transp__bedload_rate_fix_link=qb,
        ...     track_stratigraphy=True,
        ...     bed_surf_new_layer_thick=0.02,
        ...     num_cycles_to_process_strat=2,
        ... )

        We will run the model for 1299 s. This is exactly the time required for the
        first link to reach a deposition of 2 cm (Notice bed_surf_new_layer_thick=0.02).
        The evolution can is written to a file called Stratigraphy_evolution.csv

        >>> for t in range(1300):
        ...     rbd._current_t = t
        ...     rbd.run_one_step()
        ...

        >>> rbd.stratigraphy_write_evolution()

        In this case we will delete the file to keep Landlab clean

        >>> os.remove("Stratigraphy_evolution.csv")

        """
        stratigraphy.write_evolution(self)
