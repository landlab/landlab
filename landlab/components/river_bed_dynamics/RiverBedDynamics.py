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

>>> import numpy as np
>>> import copy
>>> from landlab import RasterModelGrid
>>> from landlab.components import RiverBedDynamics
>>> from landlab import imshow_grid
>>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

Create a grid on which to calculate sediment transport

>>> grid = RasterModelGrid((5, 5))

The grid will need some data to run the RiverBedDynamics component.
To check the names of the fields that provide input use the *input_var_names*
class property.

>>> RiverBedDynamics.input_var_names
('bed_surface__grain_size_distribution_location',
 'surface_water__depth',
 'surface_water__velocity',
 'surface_water__velocity_previous_time',
 'topographic__elevation')

Create fields of data for each of these input variables. When running a
complete simulation some of these variables will be created by the flow model.
Notice that surface water depth and velocity are required at links. However,
specifying these variables at nodes is easier and then we can map the fields
onto links. By doing so, we don't have to deal with links numbering. When this
component is coupled to OverlandFlow there is no need to map fields because it
is done automatically within the component.

We create the bed surface grain size (GSD) distribution location. We assume
that there are two different GSD within the watershed (labeled as 0 and 1)

>>> grid.at_node['bed_surface__grain_size_distribution_location'] = np.array([
... 0, 1., 1., 1., 0,
... 0, 1., 1., 1., 0,
... 0, 1., 1., 1., 0,
... 0, 1., 1., 1., 0,
... 0, 1., 1., 1., 0,])

Now we create the topography data

>>> grid.at_node['topographic__elevation'] = np.array([
... 1.07, 1.06, 1.00, 1.06, 1.07,
... 1.08, 1.07, 1.03, 1.07, 1.08,
... 1.09, 1.08, 1.07, 1.08, 1.09,
... 1.09, 1.09, 1.08, 1.09, 1.09,
... 1.09, 1.09, 1.09, 1.09, 1.09,])

We copy this elevation to define the original bed surface elevation

>>> grid.at_node['topographic__elevation_original'] = \
    copy.deepcopy(grid.at_node['topographic__elevation'])

We set the boundary conditions

>>> grid.set_watershed_boundary_condition(grid.at_node['topographic__elevation'])

And check the node status

>>> grid.status_at_node
array([4, 4, 1, 4, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 4, 4, 4, 4,
       4, 4], dtype=uint8)

Which tell us that there is one outlet located on the 3rd node

the topography data can be display using
>>> imshow_grid(grid,'topographic__elevation')

Now we add some water into the watershed. In this case is specified in nodes

>>> grid.at_node['surface_water__depth'] = np.array([
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,
... 0.102, 0.102, 0.102, 0.102, 0.102,])

There are other most efficient ways to fill the 'surface_water__depth', but for
demonstration purposes we show the extended version. A more efficient way to
set the previous field could be:
grid.at_node['surface_water__depth'] = np.full(grid.number_of_nodes,0.102)

Now, we give the water a velocity.

>>> grid.at_node['surface_water__velocity'] = np.array([
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,
... 0.25, 0.25, 0.25, 0.25, 0.25,])

Note that in this example, we are attempting to specify a vector at a node
using a single value. This is done intentionally to emphasize the process.
The component will interpret this as the vector's magnitude, and, given its
location in the grid, it will manifest different components. When using
OverlandFlow, there is no need to specify a velocity because it is a
byproduct of the component.

For the purpose of this illustration, we will make an assumption that the
conditions remain identical to the previous time step.

>>> grid.at_node['surface_water__velocity_previous_time'] = \
    grid.at_node['surface_water__velocity']

By default, when creating our grid we used a spacing of 1 m in the x and y
directions. Therefore, the discharge is 0.0255 m3/s. Discharge is always in
units of m3/s.

Now we map nodes into links when it is required

>>> grid['link']['surface_water__depth'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__depth')
>>> grid['link']['surface_water__velocity'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__velocity')
>>> grid['link']['surface_water__velocity_previous_time'] = \
    map_mean_of_link_nodes_to_link(grid,'surface_water__velocity_previous_time')

We assign a GSD to each location. The structure of this array is:
First column contains the grain sizes in milimiters
Second column is location 0 in 'bed_grainSizeDistribution__location'
Third column is location 1 in 'bed_grainSizeDistribution__location', and so on

>>> gsd = np.array([[32, 100, 100], [16, 25, 50], [8, 0, 0]])

Provide a time step, usually an output of OverlandFlow, but it can be
overridden with a custom value.

>>> timeStep = 1 # time step in seconds

Instantiate the `RiverBedDynamics` component to work on the grid, and run it.

>>> RBD = RiverBedDynamics(grid, gsd = gsd, dt = timeStep, bedload_equation = 'Parker1990')
>>> RBD.run_one_step()

After instantiating RiverBedDynamis, new fields have been added to the grid.
Use the *output_var_names* property to see the names of the fields
that have been changed.

>>> RBD.output_var_names
('bed_surface__geometric_mean_size',
 'bed_surface__grain_size_distribution',
 'bed_surface__median_size',
 'bed_surface__sand_fraction',
 'bed_surface__standard_deviation_size',
 'sediment_transport__bedload_grain_size_distribution',
 'sediment_transport__bedload_rate',
 'sediment_transport__net_bedload',
 'surface_water__shear_stress')

The `sediment_transport__bedload_rate` field is defined at nodes.

>>> RBD.var_loc('sediment_transport__net_bedload')
'node'

>>> grid.at_node['sediment_transport__net_bedload'] # doctest: +NORMALIZE_WHITESPACE
array([  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         1.59750963e-03,  -4.79298187e-03,   1.59750963e-03,
         0.00000000e+00,   0.00000000e+00,   1.50988732e-07,
         1.59720766e-03,   1.50988732e-07,   0.00000000e+00,
         0.00000000e+00,   3.01977464e-07,  -1.50988732e-07,
         3.01977464e-07,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
         0.00000000e+00])

The 'surface_water__shear_stress' field is defined at links.

>>> RBD.var_loc('surface_water__shear_stress')
'link'

>>> grid.at_link['surface_water__shear_stress'] # doctest: +NORMALIZE_WHITESPACE
array([  0.      ,  60.016698,  60.016698,   0.      ,   0.      ,
         0.      ,   0.      ,   0.      ,   0.      ,   0.      ,
        40.011132,  40.011132,   0.      ,  10.002783,  10.002783,
        40.011132,  10.002783,  10.002783,   0.      ,  10.002783,
        10.002783,   0.      ,   0.      ,  10.002783,  10.002783,
        10.002783,   0.      ,   0.      ,  10.002783,  10.002783,
         0.      ,   0.      ,   0.      ,   0.      ,   0.      ,
         0.      ,   0.      ,   0.      ,   0.      ,   0.      ])

Considering the link upstream the watershed exit, link Id 15, we can obtain the
bed load transport rate

>>> grid.at_link['sediment_transport__bedload_rate'][15]
-0.0015976606223666004

Therefore, the bed load transport rate according to Parker 1991 surface-based
equation is 1.598 * 10^-3 m2/s. Negative means that is going in the negative
Y direction

The GSD at this place is:

>>> grid.at_link['sediment_transport__bedload_grain_size_distribution'][15]
array([ 0.47501858,  0.52498142])

Which in cummulative percentage is equivalent to
D mm    % Finer
32      100.000
16      52.498
8       0.000

Grain sizes are always given in mm.
We can also check the bed load grain size distribution in all links

>>> grid.at_link['sediment_transport__bedload_grain_size_distribution']
array([[ 0.        ,  0.        ],
       [ 0.48479122,  0.51520878],
       [ 0.48479122,  0.51520878],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.47501858,  0.52498142],
       [ 0.47501858,  0.52498142],
       [ 0.        ,  0.        ],
       [ 0.54055384,  0.45944616],
       [ 0.28225526,  0.71774474],
       [ 0.47501858,  0.52498142],
       [ 0.28225526,  0.71774474],
       [ 0.54055384,  0.45944616],
       [ 0.        ,  0.        ],
       [ 0.28225526,  0.71774474],
       [ 0.28225526,  0.71774474],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.28225526,  0.71774474],
       [ 0.28225526,  0.71774474],
       [ 0.28225526,  0.71774474],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.28225526,  0.71774474],
       [ 0.28225526,  0.71774474],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ],
       [ 0.        ,  0.        ]])

Zeros indicate that there is no sediment transport of that grain size
at that location.

After the flow acted on the bed and sediment transport occured we can check
the new topographic elevation field

>>> grid.at_node['topographic__elevation']  # doctest: +NORMALIZE_WHITESPACE
array([ 1.07      ,  1.06      ,  1.00737382,  1.06      ,  1.07      ,
        1.08      ,  1.06754229,  1.03737382,  1.06754229,  1.08      ,
        1.09      ,  1.07999977,  1.06754276,  1.07999977,  1.09      ,
        1.09      ,  1.08999954,  1.08000023,  1.08999954,  1.09      ,
        1.09      ,  1.09      ,  1.09      ,  1.09      ,  1.09      ])
"""
import copy
import os
import shutil
import time

import numpy as np
import pandas as pd
import scipy.constants
from scipy.interpolate import interp1d

from landlab import Component, FieldError


class RiverBedDynamics(Component):

    """Landlab component that predicts the evolution of a river bed
    considering changes in elevation and grain size distribution in response to
    bed load transport according to the Exner equation and the transfer
    function of Toro-Ecobar et al., (1996).

    To estimate temporal and spatial changes in river bed properties, this
    component predicts the bedload transport rate and fractional transport at
    each link using unsteady shear stress. Time-varying hydraulic
    variables are obtained from a surface flow, for instance, OverlandFlow.

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
        "bed_surface__elevation_fixed": {
            "dtype": int,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Sets a node as a fixed elevation",
        },
        "bed_surface__geometric_mean_size": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm",
            "mapping": "link",
            "doc": "Bed surface geometric mean grain size",
        },
        "bed_surface__grain_size_distribution": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "link",
            "doc": "Bed surface grain size fraction",
        },
        "bed_surface__grain_size_distribution_fixed": {
            "dtype": int,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Determine locations that do not change its GSD in time",
        },
        "bed_surface__grain_size_distribution_location": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Accounts for the initial grain size distribution that can \
                vary in space. Each value corresponds to a given GSD location.",
        },
        "bed_surface__median_size": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm",
            "mapping": "node",
            "doc": "Bed surface median grain size",
        },
        "bed_surface__sand_fraction": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Bed surface sand content",
        },
        "bed_surface__standard_deviation_size": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm",
            "mapping": "node",
            "doc": "Bed surface standard deviation grain size",
        },
        "sediment_transport__bedload_grain_size_distribution": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "link",
            "doc": "Bed load grain size distribution",
        },
        "sediment_transport__bedload_grain_size_distribution_imposed": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "-",
            "mapping": "link",
            "doc": "Sets the sediment transport GSD where \
                  sediment supply is imposed",
        },
        "sediment_transport__bedload_rate": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m^2/s",
            "mapping": "link",
            "doc": "Volumetric bed load transport rate per unit width",
        },
        "sediment_transport__net_bedload": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m^3/s",
            "mapping": "node",
            "doc": "Net sediment transport rate at a node",
        },
        "sediment_transport__sediment_supply_imposed": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m^2/s",
            "mapping": "link",
            "doc": "Sets the sediment transport rate per unit width where \
                  sediment supply is imposed",
        },
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
            "doc": "Average velocity of the surface water",
        },
        "surface_water__velocity_previous_time": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "Average velocity of the surface water in the previous time step",
        },
        "surface_water__shear_stress": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "Pa",
            "mapping": "link",
            "doc": "Shear stress applied by the surface water on the bed \
                surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        gsd=np.array([[32, 100], [16, 25], [8, 0]]),  # Initializes the gsd array.
        rho=1000,  # Sets the fluid density (kg/m**3).
        rho_s=2650,  # Sets the sediment density (kg/m**3).
        bedload_equation="MPM",  # Selects the bedload equation.
        variable_critical_shear_stress=False,  # If set to True, allows the
        # critical shear stress to vary with the bed slope.
        use_hydraulics_radius_in_shear_stress=False,  # If set to True, uses
        # the hydraulic radius, not water depth, to
        # calculate shear stresses.
        lambda_p=0.35,  # Sets the sediment porosity.
        outlet_boundary_condition="zeroGradient",  # Sets the boundary condi-
        # tion at the watershed outlet. Available
        # options are: 'zeroGradient' and 'fixedValue'.
        evolve_bed=True,  # If set to True, allows changes in bed surface
        # elevation and GSD in response to bed load transport.
        update_bed_surface_GSD=True,  # If set to True, updates the bed surface
        # GSD considering feedbacks from bed load rate and
        # bed load GSD.
        track_stratigraphy=False,  # If set to True, records the GSD of each
        # sediment layer at each node on the hard drive.
        # This function does not use Landlab layers.
        number_cycles_to_process_stratigraphy=10,  # Accesses the read/write
        # stratigraphy routine after a set number of time
        # steps.
        new_surface_layer_thickness=1,  # Sets the threshold thickness for a
        # deposited surface layer to become a new subsurface
        # layer.
        output_vector=False,  # If set to True, enables exporting shear stress
        # and velocity vectors.
        dt=1,  # Sets the time step (s). When coupled to
        # OverlandFlow, this value is adjusted dynamically.
        alpha=1.0,  # Sets the upwinding coefficient for the central
        # difference scheme used when updating the bed GSD.
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
        grid :  RasterModelGrid
            A Landlab grid.
        gsd :   Array float, mandatory
            Grain size distribution. Must contain as many GSDs as there are
            different indexes in GSDLocation.
        rho :   float, optional
            Fluid density. Defaults to the density of water at 1000 kg/m^3.
        rho_s : float, optional
            Sediment density. Defaults to sediment density of 2650 kg/m^3.
        bedload_equation :    String, optional
            Selects the bedload transport equation. Options include:
            'MPM' for Meyer Peter and Muller,
            'FLvB' for Fernandez Luque & van Beek (1976),
            'Parker1990' for Parker 1990,
            'WilcockAndCrowe' for Wilcock and Crowe 2003.
            Default is 'MPM'.
        variable_critical_shear_stress: Boolean, optional
            If True, the critical shear stress in 'MPM' or 'FLvB' will be
            obtained using Mueller et al. (2005) equation.
        use_hydraulics_radius_in_shear_stress: Boolean, optional
            If True, shear stress calculations will use the hydraulic radius.
            Defaults to False, which uses water depth.
        lambda_p : float, optional
            Sediment porosity. Default value is 0.35.
        outlet_boundary_condition : str, optional
            Sets the boundary condition at the watershed outlet. Options are:
            'zeroGradient' (maps the outlet to the upstream node, default) or
            'fixedValue' (does not change the value at the outlet during the
            run).
        evolve_bed : Boolean, optional
            If True, the bed evolves according to the local bed load transport
            rate and GSD. If False, the bed load transport rate and GSD are
            calculated and output.
        update_bed_surface_GSD : Boolean, optional
            If True, the bed GSD evolves according to the interaction between
            the current surface GSD, and the bed load GSD and rate.
        track_stratigraphy : Boolean, optional
            If True, the component stores the GSD of each layer at every node.
            This is computationally demanding as it needs to read/write data
            at every time step. Recommended only when cycles of erosion and
            deposition are frequent and very important.
        number_cycles_to_process_stratigraphy : int, optional
            If track_stratigraphy is True, data will be read and stored every
            'number_cycles_to_process_stratigraphy' time steps. Must be larger
            or equal to 1. All data is written to the hard drive. It does not use
            Landlab layers.
        dt: float, optional
            Time step in seconds. When this component is coupled to a flow model,
            it is dynamically updated.
        alpha : float, optional
            An upwinding coefficient for a central difference scheme when
            updating the bed GSD - default value is 1.0 - a value of 0.5
            generates a central differences scheme.

        """
        super().__init__(grid)

        self._g = scipy.constants.g  # Acceleration due to gravity (m/s**2).
        self._rho = rho
        self._rho_s = rho_s
        self._R = (rho_s - rho) / rho
        self._shear_stress_star_rsgo = 0.0386  # Reference dimensionless shear
        # stress for the median size
        self._beta = 0.0951  # Coefficient for the hiding function
        self._gsd = gsd
        self._bedload_equation = bedload_equation
        self._variable_critical_shear_stress = variable_critical_shear_stress
        self._use_hydraulics_radius_in_shear_stress = (
            use_hydraulics_radius_in_shear_stress
        )
        self._lambda_p = lambda_p
        self._alpha = alpha
        self._outlet_boundary_condition = outlet_boundary_condition

        self._gsd_location = self._grid.at_node[
            "bed_surface__grain_size_distribution_location"
        ]

        self._grid._dt = dt

        # True if it is the first iteration
        self._first_iteration = True

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

        # The following variables: po, oo, so are the interpolation points for
        # the omega and sigma functions in Parker 1990
        self._po = np.array(
            [
                0.6684,
                0.7639,
                0.8601,
                0.9096,
                0.9615,
                1,
                1.055,
                1.108,
                1.197,
                1.302,
                1.407,
                1.529,
                1.641,
                1.702,
                1.832,
                1.937,
                2.044,
                2.261,
                2.499,
                2.732,
                2.993,
                3.477,
                4.075,
                4.469,
                5.016,
                6.158,
                7.821,
                10.06,
                14.38,
                19.97,
                25.79,
                38.57,
                68.74,
                91.95,
                231.2,
                2320,
            ]
        )
        self._oo = np.array(
            [
                1.011,
                1.011,
                1.01,
                1.008,
                1.004,
                0.9997,
                0.9903,
                0.9789,
                0.9567,
                0.9273,
                0.8964,
                0.8604,
                0.8287,
                0.8123,
                0.7796,
                0.7554,
                0.7326,
                0.6928,
                0.6585,
                0.6345,
                0.615,
                0.5877,
                0.564,
                0.5523,
                0.5395,
                0.5209,
                0.5045,
                0.4917,
                0.479,
                0.4712,
                0.4668,
                0.462,
                0.4578,
                0.4564,
                0.4541,
                0.4527,
            ]
        )
        self._so = np.array(
            [
                0.8157,
                0.8157,
                0.8182,
                0.8233,
                0.8333,
                0.8439,
                0.8621,
                0.8825,
                0.9214,
                0.9723,
                1.025,
                1.083,
                1.13,
                1.153,
                1.196,
                1.225,
                1.25,
                1.287,
                1.313,
                1.333,
                1.352,
                1.38,
                1.403,
                1.414,
                1.426,
                1.444,
                1.458,
                1.469,
                1.48,
                1.486,
                1.49,
                1.493,
                1.497,
                1.498,
                1.499,
                1.5,
            ]
        )

        # Creating optional grid fields at time zero- If the grid field was not
        # defined before instantiation it will create it and fill all values
        # with zeros, otherwise, it will print a message and use the given values
        try:
            self._grid["link"]["sediment_transport__bedload_rate"] = grid.add_zeros(
                "sediment_transport__bedload_rate",
                at="link",
                units=self._info["sediment_transport__bedload_rate"]["units"],
            )
        except FieldError:
            print("'sediment_transport__bedload_rate' at links - Initialized")

        try:
            self._grid["node"]["bed_surface__elevation_fixed"] = grid.add_zeros(
                "bed_surface__elevation_fixed",
                at="node",
                units=self._info["bed_surface__elevation_fixed"]["units"],
            )
        except FieldError:
            print("'bed_surface__elevation_fixed' at nodes - Initialized")

        try:
            self._grid["link"][
                "sediment_transport__sediment_supply_imposed"
            ] = grid.add_zeros(
                "sediment_transport__sediment_supply_imposed",
                at="link",
                units=self._info["sediment_transport__sediment_supply_imposed"][
                    "units"
                ],
            )
        except FieldError:
            print(
                "'sediment_transport__sediment_supply_imposed' at links - Initialized"
            )

        try:
            self._grid["node"][
                "bed_surface__grain_size_distribution_fixed"
            ] = grid.add_zeros(
                "bed_surface__grain_size_distribution_fixed",
                at="node",
                units=self._info["bed_surface__grain_size_distribution_fixed"]["units"],
            )
        except FieldError:
            print("'bed_surface__grain_size_distribution_fixed' at nodes - Initialized")

        try:
            self._grid["link"][
                "sediment_transport__bedload_grain_size_distribution_imposed"
            ] = grid.add_zeros(
                "sediment_transport__bedload_grain_size_distribution_imposed",
                at="link",
                units=self._info[
                    "sediment_transport__bedload_grain_size_distribution_imposed"
                ]["units"],
            )
        except FieldError:
            print(
                "'sediment_transport__bedload_grain_size_distribution_imposed' at links - Initialized"
            )

        # Define faces normal vector
        self._normal = -(self._grid.link_dirs_at_node)

        # Makes a copy of the original bed surface elevation and maps into links
        self._grid["node"]["topographic__elevation_original"] = copy.deepcopy(
            self._grid["node"]["topographic__elevation"]
        )
        self._grid["link"][
            "topographic__elevation_original"
        ] = self._grid.map_mean_of_link_nodes_to_link(
            self._grid["node"]["topographic__elevation_original"]
        )
        self._grid["link"][
            "topographic__elevation_subsurface"
        ] = self._grid.map_mean_of_link_nodes_to_link(
            self._grid["node"]["topographic__elevation_original"]
        )
        self._grid["link"]["bed_surface__surface_thickness_new_layer"] = (
            self._grid["link"]["topographic__elevation_original"] * 0
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
        """This method performs the initial setup of the bed grain size dis-
        tribution properties. It reads the input data and populates the nece-
        ssary variables.

        This configuration is only performed during the first time step. Subse-
        quent time steps will utilize the bed information, which has already
        been calculated or updated and formatted appropriately."""

        grain_size_D = self._gsd[:, 0]  # Grain sizes
        grain_size_frequency = self._gsd[:, 1:]  # Grain sizes frequency cumulative

        # Number of locations with different GSD during the first time step
        number_gsd_locations = self._gsd.shape[1] - 1

        # Calculates the sand fraction in each location
        if np.min(grain_size_D) < 2:  # Grain smaller than 2 mm are sand
            sand_fraction_0 = self.calculate_sand_fraction()
            # Flag is set as True to adjust the GSD and calculate a sand-free
            # GSD when using Parker 1990 eq
            self._adjust_gsd_flag = True
        else:
            sand_fraction_0 = np.zeros(self._gsd.shape[1] - 1)
            # Flag is set as False, there is no need to adjust the GSD
            self._adjust_gsd_flag = False

        # If Parker Eq is used we remove sand content
        if self._adjust_gsd_flag is True:
            # We add 2mm into GSD and update the GSD
            id2mm = np.argmin(grain_size_D >= 2)  # Location where 2 mm will be placed

            if self._bedload_equation == "Parker1990":
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
        grain_size_D_equivalent = (
            grain_size_D[0:-1] * grain_size_D[1:]
        ) ** 0.5  # Equivalent grain sizes
        sand_fraction = np.zeros_like(self._gsd_location)  # Sand fraction at each node

        # Bed grain sizes frequency in each node
        grain_size_D_equivalent_frequency = np.zeros(
            [self._grid.number_of_nodes, grain_size_D.shape[0] - 1]
        )
        for i in np.arange(0, number_gsd_locations):
            (id_gsdLocation,) = np.where(self._gsd_location == i)
            sand_fraction[id_gsdLocation] = sand_fraction_0[i] / 100
            grain_size_D_equivalent_frequency[id_gsdLocation, :] = grain_size_frequency[
                :, i
            ]

        # Calculating D50 at each node - Grain sizes in Psi scale
        grain_size_Psi_scale_D = np.flip(np.log2(grain_size_D), axis=0)
        # Frequency of each grain size in each node
        grain_size_D_equivalent_frequency_cumulative = np.hstack(
            (
                np.zeros([self._gsd_location.shape[0], 1]),
                np.cumsum(np.flip(grain_size_D_equivalent_frequency, axis=1), axis=1),
            )
        )

        # Finds the index of the grain size smaller than 50% in each node
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= 0.5, axis=1) - 1
        # Finds the index of the grain size larger than 50% in each node
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > 0.5, axis=1)
        nodes_list = np.arange(0, self._grid.number_of_nodes)

        grain_size_Psi_scale_D50 = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[nodes_list, i1]
                - grain_size_D_equivalent_frequency_cumulative[nodes_list, i0]
            )
        ) * (0.5 - grain_size_D_equivalent_frequency_cumulative[nodes_list, i0])
        D50 = 2**grain_size_Psi_scale_D50  # Median grain size in each node

        # Calculating the geometric mean and standard deviation
        Psi_grain_size_D_equivalent = np.log2(
            grain_size_D_equivalent
        )  # Equivalent grain sizes in Psi scale
        Psi_grain_size_D_equivalent_mean = np.sum(
            grain_size_D_equivalent_frequency * Psi_grain_size_D_equivalent, axis=1
        )
        grain_size_geometric_mean = 2**Psi_grain_size_D_equivalent_mean
        # Equivalent grain sizes in Psi scale at each node
        Psi_grain_size_D_equivalent = np.tile(
            Psi_grain_size_D_equivalent, (self._grid.number_of_nodes, 1)
        )
        Psi_grain_size_D_equivalent_mean = np.reshape(
            Psi_grain_size_D_equivalent_mean, [self._grid.number_of_nodes, 1]
        )
        grain_size_geometric_standard_deviation = 2 ** np.sqrt(
            np.sum(
                ((Psi_grain_size_D_equivalent - Psi_grain_size_D_equivalent_mean) ** 2)
                * grain_size_D_equivalent_frequency,
                axis=1,
            )
        )

        # We save all data into the grid so it can be shared with different
        # functions and components
        # Bed grain sizes frequency in each node
        self._grid["node"][
            "bed_surface__grain_size_distribution"
        ] = grain_size_D_equivalent_frequency
        self._grid["node"][
            "bed_surface__grain_size_distribution_Original"
        ] = copy.deepcopy(grain_size_D_equivalent_frequency)
        self._grid["node"]["bed_subsurface__grain_size_distribution"] = copy.deepcopy(
            grain_size_D_equivalent_frequency
        )
        self._grid["node"]["bed_surface__median_size"] = D50
        self._grid["node"][
            "bed_surface__geometric_mean_size"
        ] = grain_size_geometric_mean
        self._grid["node"][
            "bed_surface__geometric_standard_deviation_size"
        ] = grain_size_geometric_standard_deviation

        # GSD properties is mapped from nodes onto links
        self._grid["link"][
            "bed_surface__grain_size_distribution"
        ] = self.map_mean_of_nodes_to_link(grain_size_D_equivalent_frequency)

        self._grid["link"]["bed_surface__median_size"] = self.map_mean_of_nodes_to_link(
            D50
        )

        self._grid["link"][
            "bed_surface__geometric_mean_size"
        ] = self.map_mean_of_nodes_to_link(grain_size_geometric_mean)

        self._grid["link"][
            "bed_surface__geometric_standard_deviation_size"
        ] = self.map_mean_of_nodes_to_link(grain_size_geometric_standard_deviation)

        self._grid["link"][
            "bed_surface__sand_fraction"
        ] = self.map_mean_of_nodes_to_link(sand_fraction)

        self._grid["link"][
            "bed_surface__grain_size_distribution_Original"
        ] = copy.deepcopy(self._grid["link"]["bed_surface__grain_size_distribution"])
        self._grid["link"]["bed_subsurface__grain_size_distribution"] = copy.deepcopy(
            self._grid["link"]["bed_surface__grain_size_distribution"]
        )

        self.calculate_active_layer_thickness()

        # Sets the initial active layer thickness as zero meter deep.
        self._grid["node"][
            "bed_surface__active_layer_thickness_previous_time"
        ] = copy.deepcopy(self._grid["node"]["bed_surface__active_layer_thickness"] * 0)
        self._grid["link"][
            "bed_surface__active_layer_thickness_previous_time"
        ] = copy.deepcopy(self._grid["link"]["bed_surface__active_layer_thickness"] * 0)

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

        --- First part ---
        Calculates shear stress and bed load transport rates across the grid.

        For one time step, this generates the shear stress across a given
        grid by accounting for the local water depth, bed surface slope, and
        water velocity at links. Then, based on this shear stress and the local
        bed surface grain size distribution, the bed load transport rate is
        calculated at each link and mapped onto each node. Bed load grain size
        distributions are calculated when using Parker's 1991 or Wilcock and
        Crowe's 2003 equations. Meyer-Peter and Muller and Fernandez Luque and
        van Beek models will only calculate the total bed load transport.

        Outputs the following bed surface properties over time at every link
        in the input grid: geometric mean size, GSD, median size, sand fraction,
        standard deviation size.

        Also outputs the bed load GSD, bed load rate, and shear stress over
        time at every link in the input grid. The net bed load is output
        over time at every node."""

        if self._first_iteration is True:
            # Identify nodes with fixed elevations and surface GSD
            self.fixed_nodes_info()
            # Identify links with fixed bedload GSD
            self.fixed_links_info()
            # Identify the node upstream of the outlet
            self.outlet_nodes_info()

        # Unsteady shear stress is calculated every time step
        self.shear_stress()

        # Selects one bedload transport model to conduct all calculations
        if self._bedload_equation == "Parker1990":
            self.bedload_equation_Parker_1990()
        elif self._bedload_equation == "MPM":
            self.bedload_equation_MeyerPeter_Muller()
        elif self._bedload_equation == "FLvB":
            self.bedload_equation_FernandezLuque_VanBeek()
        elif self._bedload_equation == "WilcockAndCrowe":
            self.bedload_equation_Wilcock_Crowe_2003()

        # Calculates the net bedload transport from links into nodes
        self.calculate_net_bedload()

        # Maps bed load vectors from links (m2/s) onto nodes (m3/s) preserving
        # vector components.
        if self._output_vector is True:
            (
                self._shear_stress_vector,
                self._shear_stress_magnitude,
            ) = self.vector_mapper(self._shear_stress)
            self._u_vector, self._u_magnitude = self.vector_mapper(self._u)
            self._qb_vector, self._qb_magnitude = self.bedload_mapper()

        """Second Part ---

        Erode grid topography.

        For one time step, this erodes the grid topography according to
        Exner equation.

        (1-λp) ∂Z/∂t = - (∂qbx/∂x + ∂qby/∂y)
        Simplifying, we get:
        ∂Z/∂t = - (1 / (1-λp)) * (∂Qb/∂A)
        Z_t+1 = -(Δt * ΔQb)/(1-λp) + Z_t

        The grid field 'topographic__elevation' is altered each time step.
        """
        if self._evolve_bed is True:
            # This routine is entered only after the second time step to avoid
            # re-calculation.
            # When the bed has undergone changes, new bed properties are calculated.
            # This updated information will then be used when calculating bedload rates.

            # Update the bed elevation
            self.update_bed_elevation()

            # Update the bed grain size distribution
            if (self._bedload_equation == "Parker1990") or (
                self._bedload_equation == "WilcockAndCrowe"
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
        self._grid["node"]["bed_surface__elevation_fixed"] = np.round(
            self._grid["node"]["bed_surface__elevation_fixed"]
        ).astype(int)
        self._grid["node"]["bed_surface__grain_size_distribution_fixed"] = np.round(
            self._grid["node"]["bed_surface__grain_size_distribution_fixed"]
        ).astype(int)
        (self._fixed_nodes_id,) = np.where(
            self._grid["node"]["bed_surface__elevation_fixed"] == 1
        )
        (self._fixed_surface_gsd_nodes_id,) = np.where(
            self._grid["node"]["bed_surface__grain_size_distribution_fixed"] == 1
        )

        # All connecting links to these nodes will be set as fixed too
        self._fixed_links = np.unique(
            self._grid.links_at_node[self._fixed_surface_gsd_nodes_id]
        )
        self._fixed_links = self._fixed_links[np.where(self._fixed_links >= 0)]

    def fixed_links_info(self):
        """Search and identify links defined as having a fixed bed load GSD"""

        # Gives the ID of the links
        if np.max(
            self._grid["link"][
                "sediment_transport__bedload_grain_size_distribution_imposed"
            ]
            > 0
        ):
            (a0, a1) = np.where(
                self._grid["link"][
                    "sediment_transport__bedload_grain_size_distribution_imposed"
                ]
                > 0
            )
            self._fixed_surface_gsd_nodes_id_row = a0
            self._fixed_surface_gsd_nodes_id_col = a1

            # We use this flag to avoid extra calculations when there is no fixed
            # bed load GSD
            if a0.shape[0] > 0:
                self._fixed_link_calculate = True

    def calculate_sand_fraction(self):
        grain_size_Psi_scale_D = np.flip(np.log2(self._gsd[:, 0]), axis=0)
        grain_size_D_equivalent_frequency_cumulative = np.flip(self._gsd[:, 1:], axis=0)
        sand_fraction = np.zeros(grain_size_D_equivalent_frequency_cumulative.shape[1])
        i = np.max(np.where(grain_size_Psi_scale_D <= 1))
        for j in np.arange(0, grain_size_D_equivalent_frequency_cumulative.shape[1]):
            sand_fraction[j] = (
                (1 - grain_size_Psi_scale_D[i])
                / (grain_size_Psi_scale_D[i + 1] - grain_size_Psi_scale_D[i])
            ) * (
                grain_size_D_equivalent_frequency_cumulative[i + 1, j]
                - grain_size_D_equivalent_frequency_cumulative[i, j]
            ) + grain_size_D_equivalent_frequency_cumulative[
                i, j
            ]

        return sand_fraction

    def map_mean_of_nodes_to_link(self, grain_size_variable):
        return 0.5 * (
            grain_size_variable[self.grid.nodes_at_link[:, 0]]
            + grain_size_variable[self.grid.nodes_at_link[:, 1]]
        )

    def calculate_active_layer_thickness(self):
        """First for nodes"""

        grain_size_D = self._grain_size_D_original  # Grain sizes
        grain_size_Psi_scale_D = np.flip(
            np.log2(grain_size_D), axis=0
        )  # Grain sizes in Psi scale

        grain_size_D_equivalent_frequency = self._grid["node"][
            "bed_surface__grain_size_distribution"
        ]
        frequency_grain_size_list = np.arange(0, self._grid.number_of_nodes)

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equivalent_frequency,
                np.zeros([grain_size_D_equivalent_frequency.shape[0], 1]),
            )
        )  # Grain sizes freq

        grain_size_D_equivalent_frequency_cumulative = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= 0.9, axis=1) - 1
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > 0.9, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i1
                ]
                - grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i0
                ]
            )
        ) * (
            0.9
            - grain_size_D_equivalent_frequency_cumulative[
                frequency_grain_size_list, i0
            ]
        )
        D90_surface = 2**grain_size_Psi_scale_DX
        self._grid["node"]["bed_surface__active_layer_thickness"] = (
            2 * D90_surface / 1000
        )

        """ Now for links"""

        grain_size_D_equivalent_frequency = self._grid["link"][
            "bed_surface__grain_size_distribution"
        ]
        frequency_grain_size_list = np.arange(0, self._grid.number_of_links)

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equivalent_frequency,
                np.zeros([grain_size_D_equivalent_frequency.shape[0], 1]),
            )
        )  # Grain sizes freq

        grain_size_D_equivalent_frequency_cumulative = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= 0.9, axis=1) - 1
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > 0.9, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i1
                ]
                - grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i0
                ]
            )
        ) * (
            0.9
            - grain_size_D_equivalent_frequency_cumulative[
                frequency_grain_size_list, i0
            ]
        )
        D90_surface = 2**grain_size_Psi_scale_DX  # 90th percentile
        self._grid["link"]["bed_surface__active_layer_thickness"] = (
            2 * D90_surface / 1000
        )

    def outlet_nodes_info(self):
        """Search and identify the node upstream the outlet to apply boundary
        conditions"""

        # Gives the ID of the outlet node
        (self._out_id,) = np.where(self._grid.status_at_node == 1)

        # nCols = self._grid.number_of_node_columns
        # nRows = self._grid.number_of_node_rows
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

    def map_gsd_from_link_to_node(self):
        """Maps the bed surface grain size distribution from links to nodes.

        Given that the all our calculations are conducted in links we implemented
        this function to display results in a raster or in nodes."""

        grain_size_D_equivalent_frequency = self._grid["link"][
            "bed_surface__grain_size_distribution"
        ]

        grain_size_D_equivalent_frequency_nodes = 0.25 * (
            grain_size_D_equivalent_frequency[self._grid.links_at_node[:, 0], :]
            + grain_size_D_equivalent_frequency[self._grid.links_at_node[:, 1], :]
            + grain_size_D_equivalent_frequency[self._grid.links_at_node[:, 2], :]
            + grain_size_D_equivalent_frequency[self._grid.links_at_node[:, 3], :]
        )

        grain_size_D_equivalent_frequency_nodes = (
            grain_size_D_equivalent_frequency_nodes
            / np.reshape(
                np.sum(grain_size_D_equivalent_frequency_nodes, axis=1),
                [self._grid.number_of_nodes, 1],
            )
        )

        # Revert any changes to the fixed GSD nodes and fixed elevations nodes
        fixed_nodes = np.unique(
            np.hstack((self._fixed_surface_gsd_nodes_id, self._fixed_nodes_id))
        )

        grain_size_D_equivalent_frequency_nodes[fixed_nodes] = self._grid["node"][
            "bed_surface__grain_size_distribution_Original"
        ][fixed_nodes]
        self._grid["node"]["bed_surface__grain_size_distribution"] = copy.deepcopy(
            grain_size_D_equivalent_frequency_nodes
        )

    def update_bed_surface_properties(self):
        """Calculates the updated GSD properties"""
        number_links = self._grid.number_of_links
        grain_size_D = self._grain_size_D_original  # Grain sizes
        grain_size_D_equivalent = self._grain_size  # Equivalent grain sizes

        # Equivalent grain sizes frequency
        grain_size_D_equivalent_frequency = self._grid["link"][
            "bed_surface__grain_size_distribution"
        ]
        # Grain sizes frequency
        grain_size_frequency = np.hstack(
            (
                grain_size_D_equivalent_frequency,
                np.zeros([grain_size_D_equivalent_frequency.shape[0], 1]),
            )
        )

        """ Calculates D50 at each node based on the updated GSD """
        grain_size_Psi_scale_D = np.flip(
            np.log2(grain_size_D), axis=0
        )  # Grain sizes in Psi scale
        grain_size_D_equivalent_frequency_cumulative = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each link

        # Finds the index of the grain size smaller than 50% in each link
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= 0.5, axis=1) - 1
        # Finds the index of the grain size larger than 50% in each link
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > 0.5, axis=1)
        linkList = np.arange(0, number_links)

        grain_size_Psi_scale_D50 = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[linkList, i1]
                - grain_size_D_equivalent_frequency_cumulative[linkList, i0]
            )
        ) * (0.5 - grain_size_D_equivalent_frequency_cumulative[linkList, i0])
        D50 = 2**grain_size_Psi_scale_D50  # Median grain size in each link

        # Calculates the geometric mean and standard deviation
        # Equivalent grain sizes in Psi scale
        Psi_grain_size_D_equivalent = np.log2(grain_size_D_equivalent)
        Psi_grain_size_D_equivalent_mean = np.sum(
            grain_size_D_equivalent_frequency * Psi_grain_size_D_equivalent, axis=1
        )
        grain_size_geometric_mean = (
            2**Psi_grain_size_D_equivalent_mean
        )  # Geometric mean size in each link

        # Equivalent grain sizes in Psi scale at each link
        Psi_grain_size_D_equivalent = np.tile(
            Psi_grain_size_D_equivalent, (number_links, 1)
        )
        Psi_grain_size_D_equivalent_mean = np.reshape(
            Psi_grain_size_D_equivalent_mean, [number_links, 1]
        )
        # Standard deviation at each node
        grain_size_geometric_standard_deviation = 2 ** np.sqrt(
            np.sum(
                ((Psi_grain_size_D_equivalent - Psi_grain_size_D_equivalent_mean) ** 2)
                * grain_size_D_equivalent_frequency,
                axis=1,
            )
        )

        self._grid["link"]["bed_surface__median_size"] = D50
        self._grid["link"][
            "bed_surface__geometric_mean_size"
        ] = grain_size_geometric_mean
        self._grid["link"][
            "bed_surface__geometric_standard_deviation_size"
        ] = grain_size_geometric_standard_deviation

        self._grid["link"][
            "bed_surface__active_layer_thickness_previous_time"
        ] = copy.deepcopy(self._grid["link"]["bed_surface__active_layer_thickness"])
        self.calculate_active_layer_thickness()

    def shear_stress(self):
        """Unsteady shear stress calculated at links according to
        τ = rho * g * h * Sf

        where Sf is the unsteady friction slope and is calculated as
        Sf = S0 - dh/ds - U/g du/ds - 1/g du/dt

        Alternatively, tau can be calculated as τ = rho * g * Rh * Sf
        but need to be manually selected (for consistency in the code)

        The term ds indicates a certain direction, X and Y in this case.
        All derivatives are first or second order approximations in each
        direction. Most of grid parameters are read again because they may have
        changed in time in OverlandFlow or when the bed evolved."""

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
                (
                    np.reshape(
                        u_nodes_vertical, (self._grid._shape[0], self._grid._shape[1])
                    )
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
        u_at_previous_time = self._grid["link"]["surface_water__velocity_previous_time"]
        du_dt = (self._u - u_at_previous_time) / self._grid._dt

        # Friction slope calculation at links including unsteady effects
        Sf = self._dz_ds - dh_ds - (self._u / self._g) * du_ds - 1 / self._g * du_dt

        # And finally, the shear stress at links including unsteady effects
        if self._use_hydraulics_radius_in_shear_stress is True:
            # uses hydraulics ratio, so shear stress is calculated as
            # τ = Rho * g * Rh * Sf

            # Different grid sizes (dx~=dy) are possible
            A = np.zeros_like(h_links)
            A[self._grid.horizontal_links] = (
                h_links[self._grid.horizontal_links] * self._grid.dx
            )
            A[self._grid.vertical_links] = (
                h_links[self._grid.vertical_links] * self._grid.dy
            )

            P = np.zeros_like(h_links)
            P[self._grid.horizontal_links] = (
                self._grid.dx + 2 * h_links[self._grid.horizontal_links]
            )
            P[self._grid.vertical_links] = (
                self._grid.dy + 2 * h_links[self._grid.vertical_links]
            )
            Rh = (
                A / P
            )  # Rh = wetted area / wetted perimeter#Rh = wetted area / wetted perimeter
            self._shear_stress = self._rho * self._g * Rh * Sf
        else:
            # Equation is τ = Rho * g * h * Sf """
            self._shear_stress = self._rho * self._g * h_links * Sf

        # Apply a boundary condition of zero flux at link next to borders
        self._shear_stress[self._links_at_border_cells] = 0

        # Direction of flux will be recovered at the end of the bedload
        # transport routine
        self._shear_stress_total = np.abs(self._shear_stress)

        # Now we write the shear stress field into the grid so other components
        # or postprocess functions can access it
        self._grid["link"]["surface_water__shear_stress"] = self._shear_stress_total

    def bedload_equation_Parker_1990(self):
        """Surface-based bedload transport equation of Parker 1990

        G. Parker (1990) Surface-based bedload transport relation for gravel rivers,
        Journal of Hydraulic Research, 28:4, 417-436, DOI: 10.1080/00221689009499058
        """

        # Variables definition - All these variables are updated each time step.
        # Therefore, we read them again to ensure they are updated

        grain_size_frequency = self._grid["link"][
            "bed_surface__grain_size_distribution"
        ]
        grain_size_geometric_mean = self._grid["link"][
            "bed_surface__geometric_mean_size"
        ]
        grain_size_geometric_standard_deviation = self._grid["link"][
            "bed_surface__geometric_standard_deviation_size"
        ]

        shear_stress_star_sg = self._shear_stress_total / (
            self._rho * self._R * self._g * (grain_size_geometric_mean / 1000)
        )
        self._phi_sgo = shear_stress_star_sg / self._shear_stress_star_rsgo

        self.strainFunctions()
        self._omega = 1 + (
            np.log2(grain_size_geometric_standard_deviation) / self._sigma0
        ) * (self._omega0 - 1)

        Di_Dsg = np.tile(self._grain_size, (self._grid.number_of_links, 1)) / (
            np.reshape(grain_size_geometric_mean, [self._grid.number_of_links, 1])
        )
        phi_sgo = np.reshape(self._phi_sgo, [self._phi_sgo.shape[0], 1])
        omega = np.reshape(self._omega, [self._omega.shape[0], 1])
        phi_i = omega * phi_sgo * (Di_Dsg) ** -self._beta

        G = np.zeros_like(phi_i)

        # There are three intervals where G is evaluated
        (id0, id1) = np.where(phi_i > 1.59)
        if id0.shape[0] > 0:
            G[id0, id1] = 5474 * (1 - 0.853 / phi_i[id0, id1]) ** 4.5

        (id0, id1) = np.where((phi_i >= 1) & (phi_i <= 1.59))
        if id0.shape[0] > 0:
            G[id0, id1] = np.exp(
                14.2 * (phi_i[id0, id1] - 1) - 9.28 * (phi_i[id0, id1] - 1) ** 2
            )

        (id0, id1) = np.where(phi_i < 1)
        if id0.shape[0] > 0:
            G[id0, id1] = phi_i[id0, id1] ** 14.2

        Gf = grain_size_frequency * G
        Gf_sum = np.sum(Gf, axis=1)

        # Total bedload transport rate
        W_star_s = 0.00218 * Gf_sum
        self._grid["link"]["sediment_transport__bedload_rate"] = (
            ((np.sqrt(self._shear_stress_total / self._rho)) ** 3 * W_star_s)
            / (self._R * self._g)
        ) * np.sign(self._shear_stress)

        # Now that bed load rate has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate is True:
            link_id_fixed_links = np.unique(self._fixed_surface_gsd_nodes_id_row)

            self._grid["link"]["sediment_transport__bedload_rate"][
                link_id_fixed_links
            ] = self._grid["link"]["sediment_transport__sediment_supply_imposed"][
                link_id_fixed_links
            ]

        # During the calculation of fractional sediment transport, there might
        # be instances where certain fractions will be zero. This situation
        # could trigger a division by zero warning. However, we are aware that
        # in these specific cases, there's no transport occurring. Thus, such a
        # warning is not indicative of an actual problem in the calculation.
        # To avoid unnecessary warnings, we've deactivated them.

        np.seterr(divide="ignore", invalid="ignore")

        # Frational bedload transport rate
        p = np.zeros_like(Gf)
        Gf_sum = np.transpose(np.tile(Gf_sum, (grain_size_frequency.shape[1], 1)))
        p = Gf / Gf_sum
        id0 = np.isnan(p)
        p[id0] = 0
        self._grid["link"]["sediment_transport__bedload_grain_size_distribution"] = p

        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate is True:
            linkId_0 = self._fixed_surface_gsd_nodes_id_row
            linkId_1 = self._fixed_surface_gsd_nodes_id_col

            self._grid["link"]["sediment_transport__bedload_grain_size_distribution"][
                linkId_0, linkId_1
            ] = self._grid["link"][
                "sediment_transport__bedload_grain_size_distribution_imposed"
            ][
                linkId_0, linkId_1
            ]

    def strainFunctions(self):
        omega0 = np.zeros_like(self._phi_sgo)
        sigma0 = np.zeros_like(self._phi_sgo)

        # There are three intervals where we can interpolate Omega and Sigma
        (I,) = np.where(self._phi_sgo <= 0.7639)
        if I.shape[0] > 0:
            omega0[I] = self._oo[0]
            sigma0[I] = self._so[0]

        (I,) = np.where(self._phi_sgo > 231.2)
        if I.shape[0] > 0:
            omega0[I] = np.interp(self._phi_sgo, self._po, self._oo)[I]
            sigma0[I] = np.interp(self._phi_sgo, self._po, self._oo)[I]

        (I,) = np.where((self._phi_sgo > 0.7639) & (self._phi_sgo < 231.2))
        if I.shape[0] > 0:
            foo = interp1d(self._po, self._oo, kind="cubic")
            fso = interp1d(self._po, self._so, kind="cubic")
            omega0[I] = foo(self._phi_sgo[I])
            sigma0[I] = fso(self._phi_sgo[I])

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
        grain_size_frequency = self._grid["link"][
            "bed_surface__grain_size_distribution"
        ]
        grain_size_geometric_mean = self._grid["link"][
            "bed_surface__geometric_mean_size"
        ]
        sand_fraction = self._grid["link"]["bed_surface__sand_fraction"]

        shear_stress_star_sg = self._shear_stress_total / (
            self._rho * self._R * self._g * (grain_size_geometric_mean / 1000)
        )
        shear_Stress_star_rsg0 = 0.021 + 0.015 * np.exp(-20 * sand_fraction)
        phi_sg0 = shear_stress_star_sg / shear_Stress_star_rsg0

        phi_i = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))

        G = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))
        p = np.zeros((phi_sg0.shape[0], self._grain_size.shape[0]))

        for i in np.arange(0, self._grain_size.shape[0]):
            b = 0.67 / (
                1 + np.exp(1.5 - self._grain_size[i] / grain_size_geometric_mean)
            )
            phi_i[:, i] = phi_sg0 * (
                self._grain_size[i] / (grain_size_geometric_mean)
            ) ** (-b)

        # There are two intervals where G is evaluated
        (id0, id1) = np.where(phi_i >= 1.35)
        if id0.shape[0] > 0:
            G[id0, id1] = 14 * (1 - 0.894 / (phi_i[id0, id1] ** 0.5)) ** 4.5

        (id0, id1) = np.where(phi_i < 1.35)
        if id0.shape[0] > 0:
            G[id0, id1] = 0.002 * phi_i[id0, id1] ** 7.5

        W_star_i = grain_size_frequency * G
        W_star = np.sum(W_star_i, axis=1)

        # Total bedload transport rate is calculated at each link
        # Notice that W_star includes fi (Eq 2 in the paper)
        self._grid["link"]["sediment_transport__bedload_rate"] = (
            ((np.sqrt(self._shear_stress_total / self._rho)) ** 3 * W_star)
            / (self._R * self._g)
        ) * np.sign(self._shear_stress)

        # During the calculation of fractional sediment transport, there might
        # be instances where certain fractions will be zero. This situation
        # could trigger a division by zero warning. However, we are aware that
        # in these specific cases, there's no transport occurring. Thus, such a
        # warning is not indicative of an actual problem in the calculation.
        # To avoid unnecessary warnings, we've deactivated them.

        np.seterr(divide="ignore", invalid="ignore")
        p = W_star_i / np.reshape(W_star, [W_star.shape[0], 1])
        id0 = np.isnan(p)
        p[id0] = 0
        self._grid["link"]["sediment_transport__bedload_grain_size_distribution"] = p

        # Now that bedload GSD has been calculated in all links we replace those that
        # were imposed. If there there are no imposed bedload GSD this will do nothing.
        if self._fixed_link_calculate is True:
            linkId_0 = self._fixed_surface_gsd_nodes_id_row
            linkId_1 = self._fixed_surface_gsd_nodes_id_col

            self._grid["link"]["sediment_transport__bedload_grain_size_distribution"][
                linkId_0, linkId_1
            ] = self._grid["link"][
                "sediment_transport__bedload_grain_size_distribution_imposed"
            ][
                linkId_0, linkId_1
            ]

    def bedload_equation_MeyerPeter_Muller(self):
        """Surface-based bedload transport equation of Meyer-Peter and Müller

        Meyer-Peter, E. and Müller, R., 1948, Formulas for Bed-Load Transport,
        Proceedings, 2nd Congress, International Association of Hydraulic
        Research, Stockholm: 39-64.
        """
        D50 = self._grid["link"]["bed_surface__median_size"]
        shear_stress_star = self._shear_stress_total / (
            self._rho * self._R * self._g * (D50 / 1000)
        )
        self._critical_shear_stress_star = 0.047
        if self._variable_critical_shear_stress is True:
            self._critical_shear_stress_star = self.critical_shear_stress_star()
        qb_star = np.where(
            shear_stress_star - self._critical_shear_stress_star > 0,
            8 * np.abs(shear_stress_star - self._critical_shear_stress_star) ** (3 / 2),
            0,
        )

        self._grid["link"]["sediment_transport__bedload_rate"] = (
            qb_star * (np.sqrt(self._R * self._g * (D50 / 1000)) * (D50 / 1000))
        ) * np.sign(self._shear_stress)

    def bedload_equation_FernandezLuque_VanBeek(self):
        """Surface-based bedload transport equation of Fernandez Luque and
        van Beek

        Fernandez Luque, R. and R. van Beek, 1976, Erosion and transport of
        bedload sediment, Journal of Hydraulic Research, 14(2): 127-144.
        """
        D50 = self._grid["link"]["bed_surface__median_size"]
        shear_stress_star = self._shear_stress_total / (
            self._rho * self._R * self._g * (D50 / 1000)
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
        self._grid["link"]["sediment_transport__bedload_rate"] = (
            qb_star * (np.sqrt(self._R * self._g * (D50 / 1000)) * (D50 / 1000))
        ) * np.sign(self._shear_stress)

    def critical_shear_stress_star(self):
        """Uses the current bed slope at links for calculating a slope-dependent
        critical shear stress. This assumes that the flow direction always goes
        in the direction of descending slope. However, there might be local
        cases where this assumption does not hold true"""
        bed_slope = np.abs(self._dz_ds)
        # Mueller et al. (2005) equation is used in steep slopes
        return np.where(
            bed_slope > 0.03, 2.18 * bed_slope + 0.021, self._critical_shear_stress_star
        )

    def calculate_net_bedload(self):
        """Calculates the net volumetric bedload coming from all links (m2/s)
        onto nodes (m3/s).

        This method takes the volumetric bedload entering and exiting through a
        face and determines the net volumetric bedload on a given node."""

        # Reads and modify the field sediment_transport__bedload_rate to account
        # for links where sediment supply is imposed.
        (self._Id_link_upstream_sediment_supply,) = np.nonzero(
            self._grid["link"]["sediment_transport__sediment_supply_imposed"]
        )
        self._grid["link"]["sediment_transport__bedload_rate"][
            self._Id_link_upstream_sediment_supply
        ] = self._grid["link"]["sediment_transport__sediment_supply_imposed"][
            self._Id_link_upstream_sediment_supply
        ]

        qb_x = (
            np.sum(
                self._grid["link"]["sediment_transport__bedload_rate"][
                    self._grid.links_at_node[:, [0, 2]]
                ]
                * self._normal[:, [0, 2]],
                axis=1,
            )
            * self._grid.dy
        )
        qb_y = (
            np.sum(
                self._grid["link"]["sediment_transport__bedload_rate"][
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

        self._grid["node"]["sediment_transport__net_bedload"] = net_qb_nodes

    def update_bed_elevation(self):
        """Applies the Exner equation and boundary conditions to predict
        the change in bed surface elevation"""

        # Bed elevation is updated using Exner equation
        A = self._grid.dx * self._grid.dy
        d_qb = self._grid["node"]["sediment_transport__net_bedload"]
        z0 = copy.deepcopy(self._grid["node"]["topographic__elevation"])
        dz = -self._grid._dt / ((1 - self._lambda_p) * A) * d_qb

        self._grid["node"]["topographic__elevation"] += dz

        # The outlet node has been modifed, but not according to the specified
        # boundary condition. Here, we return the outlet to the previous state
        self._grid["node"]["topographic__elevation"][self._out_id] = z0[self._out_id]

        # Fixed nodes may have been modifed, but not according to the specified
        # condition. Here, we return the nodes to the orignal state
        self._grid["node"]["topographic__elevation"][self._fixed_nodes_id] = self._grid[
            "node"
        ]["topographic__elevation_original"][self._fixed_nodes_id]

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

        if self._id_eroded_nodes.shape[0] > 0:
            self._grid["node"]["surface_water__depth"][self._id_eroded_nodes] += np.abs(
                dz_current_time[self._id_eroded_nodes]
            )

        if (
            self._track_stratigraphy is True
        ):  # Here we register how deep the deposited/eroded layer is
            self._grid["link"]["bed_surface__surface_thickness_new_layer"] = (
                self._grid["link"]["topographic__elevation"]
                - self._grid["link"]["topographic__elevation_subsurface"]
            )

            # Checks if deposited material needs to be updated
            (self._id_deep_links,) = np.where(
                self._grid["link"]["bed_surface__surface_thickness_new_layer"]
                > self._new_surface_layer_thickness
            )
            self._id_deep_links = self._id_deep_links[
                np.in1d(self._id_deep_links, self._grid.active_links)
            ]
            self._id_deep_links = self._id_deep_links[
                ~np.in1d(self._id_deep_links, self._fixed_links)
            ]

            if self._id_deep_links.shape[0] > 0:
                if (
                    np.min(
                        self._grid["link"]["bed_surface__active_layer_thickness"][
                            self._id_deep_links
                        ]
                    )
                    > self._new_surface_layer_thickness
                ):
                    print(
                        "Warning - New surface layer thickness is too thin compared to the active layer thickness"
                    )
                self._compute_stratigraphy = True
                self._update_subsurface = True

            # Checks if eroded material needs to be updated
            (self._id_eroded_links,) = np.where(
                self._grid["link"]["bed_surface__surface_thickness_new_layer"]
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
        nCols = self._grid.number_of_node_columns
        F = self._grid["link"]["bed_surface__grain_size_distribution"]
        Fs = self._grid["link"]["bed_subsurface__grain_size_distribution"]
        pl = self._grid["link"]["sediment_transport__bedload_grain_size_distribution"]
        qbT = self._grid["link"][
            "sediment_transport__bedload_rate"
        ]  # total bed load transport rate at each link
        La = np.reshape(
            self._grid["link"]["bed_surface__active_layer_thickness"], [number_links, 1]
        )
        Laold = np.reshape(
            self._grid["link"]["bed_surface__active_layer_thickness_previous_time"],
            [number_links, 1],
        )

        lps = self._lambda_p
        dx = self._grid.dx
        dy = self._grid.dy
        alpha = self._alpha
        dt = self._grid._dt

        dv = 2 * nCols - 1

        qbT = np.reshape(qbT, [number_links, 1])
        qbTdev = np.zeros([number_links, 1])
        qbTdevNeg = np.zeros([number_links, 1])
        qjj1dev = np.zeros_like(pl)
        qjj1devNeg = np.zeros_like(pl)

        # Horizontal Links
        hlL = np.arange(
            0, number_links - nCols + 2, 2 * nCols - 1
        )  # Horizontal Links at left edge
        hlR = np.arange(
            nCols - 2, number_links, 2 * nCols - 1
        )  # Horizontal Links at right edge
        hl_Id = np.in1d(
            self._grid.horizontal_links, np.hstack((hlL, hlR))
        )  # Links at the top row
        hl = self._grid.horizontal_links[~hl_Id]

        # Vertical Links
        vl = self._grid.vertical_links[nCols:-nCols]  # Links within middle region
        vlB = self._grid.vertical_links[0:nCols]  # Links at the bottom row
        vlT = self._grid.vertical_links[-nCols:None]  # Links at the top row

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

        # Now we correct for flow at locations where it is flowing towards the negative direction
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

        FIexc = copy.deepcopy(Fs)
        (id0,) = np.where(qbTdev[:, 0] <= 0)
        FIexc[id0, :] = 0.7 * F[id0, :] + 0.3 * pl[id0, :]

        qjj2dev = FIexc * np.reshape(qbTdev, [number_links, 1])

        Fnew = F + dt * (-qjj1dev + qjj2dev) / (1 - lps) / La
        if (self._track_stratigraphy is False) and (self._first_iteration is False):
            Fnew = Fnew + (FIexc - F) / La * (La - Laold)
        elif (self._track_stratigraphy is True) and (self._first_iteration is False):
            if self._subsurface_changed is True:
                Fnew = Fnew + (FIexc - F) / La * (La - Laold)
                self._subsurface_changed = False

        (id0, id1) = np.where(Fnew <= 0)
        Fnew[id0, id1] = 0
        Fnew = Fnew / np.reshape(np.sum(Fnew, axis=1), [number_links, 1])
        Fnew = np.nan_to_num(Fnew)

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
                    Fnew[self._outlet_links_4_horizontal[i][1]]
                    - Fnew[self._outlet_links_4_horizontal[i][0]]
                ) / ds
                b = Fnew[self._outlet_links_4_horizontal[i][0]]
                Fnew[self._outlet_links_horizontal[i][0]] = m * 2 * ds + b
                Fnew[self._outlet_links_horizontal[i][1]] = m * 3 * ds + b

        # Corrects vertical links

        # In very small grids we may not have enough links, so we limit the
        # correction to those that have more enough data
        if self._outlet_links_4_vertical.shape[1] > 1:
            for i in np.arange(0, self._outlet_links_vertical.shape[0]):
                ds = dy
                m = (
                    Fnew[self._outlet_links_4_vertical[i][1]]
                    - Fnew[self._outlet_links_4_vertical[i][0]]
                ) / ds
                b = Fnew[self._outlet_links_4_vertical[i][0]]
                Fnew[self._outlet_links_vertical[i][0]] = m * 2 * ds + b
                Fnew[self._outlet_links_vertical[i][1]] = m * 3 * ds + b

        # Now we update the bed surface GSD
        self._grid["link"]["bed_surface__grain_size_distribution"] = copy.deepcopy(Fnew)

        # If the bed is eroded below the original elevation it restores this
        # initial GSD. First looks if there is any such link
        (id_eroded_links,) = np.where(
            self._grid["link"]["topographic__elevation"]
            < self._grid["link"]["topographic__elevation_original"]
        )

        if id_eroded_links.shape[0] > 0:
            self._grid["link"]["bed_surface__grain_size_distribution"][
                id_eroded_links
            ] = self._grid["link"]["bed_surface__grain_size_distribution_Original"][
                id_eroded_links
            ]

            # Now, an eroded node cannot return to the original GSD if starts depositing again.
            # It will only use the original GSD if erodes deeper than the maximum that has been eroded
            self._grid["link"]["topographic__elevation_original"][
                id_eroded_links
            ] = copy.deepcopy(
                self._grid["link"]["topographic__elevation"][id_eroded_links]
            )

        # Revert any changes to the fixed GSD nodes
        self._grid["link"]["bed_surface__grain_size_distribution"][
            self._fixed_links
        ] = self._grid["link"]["bed_surface__grain_size_distribution_Original"][
            self._fixed_links
        ]

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
        # Here we create a number of variables that will be used in the
        # following definitions to make the code a little bit cleaner

        z = self._grid["link"]["topographic__elevation"]
        grain_size_D = self._grain_size_D_original
        F = self._grid["link"]["bed_surface__grain_size_distribution"]
        Fs = self._grid["link"]["bed_subsurface__grain_size_distribution"]
        La = self._grid["link"]["bed_surface__active_layer_thickness"]
        dzl = self._grid["link"]["bed_surface__surface_thickness_new_layer"]

        if self._first_iteration is True:
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
                data = np.hstack((self.t, z[i], F[i, :], Fs[i, :]))
                data = np.reshape(data, [1, data.shape[0]])
                with open(filename, "ab") as f:
                    np.savetxt(f, data, "%.3f")

        # Here we store data
        os.chdir(self._cwd)
        os.chdir(self._stratigraphy_temp_files_path)
        for i in self._grid.active_links:
            filename = "link_" + str(i) + ".txt"
            data = np.hstack((self.t, z[i], dzl[i], La[i], F[i, :], Fs[i, :]))
            data = np.reshape(data, [1, data.shape[0]])
            # with open(filename, 'ab') as f:
            #    np.savetxt(f, data,'%.3f')
            # Try up to 3 times
            for attempt in range(3):
                try:
                    with open(filename, "ab") as f:
                        np.savetxt(f, data, "%.3f")
                    break  # Successfully written, exit the loop
                except PermissionError:
                    print(
                        f"PermissionError when trying to write to {filename}. Retry in 1 second."
                    )
                    time.sleep(1)  # Wait for 1 second before retrying

        # Here we update stratigraphy in case of deposition - only for new layers
        os.chdir(self._cwd)
        os.chdir(self._stratigraphy_temp_files_path)
        if self._update_subsurface is True:
            # print('Updating stratigraphy at time ',self.t, 's. New layer(s) created at links ',self._id_deep_links)
            for i in self._id_deep_links:
                link_data = np.loadtxt("link_" + str(i) + ".txt")
                mean_subsurface_gsd = np.mean(
                    link_data[:, 4 : 4 + grain_size_D.shape[0] - 1], axis=0
                )  # 4 is the number of elements before F GSD
                self._grid["link"]["bed_subsurface__grain_size_distribution"][
                    i, :
                ] = mean_subsurface_gsd / np.sum(mean_subsurface_gsd)
                os.remove("link_" + str(i) + ".txt")

            # Now updates the raw data for the links with deposition
            os.chdir(self._cwd)
            os.chdir(self._stratigraphy_raw_data_path)

            for i in self._id_deep_links:
                filename = "link_" + str(i) + ".txt"
                data = np.hstack(
                    (
                        self.t,
                        z[i],
                        F[i, :],
                        self._grid["link"]["bed_subsurface__grain_size_distribution"][
                            i, :
                        ],
                    )
                )
                data = np.reshape(data, [1, data.shape[0]])
                with open(filename, "ab") as f:
                    np.savetxt(f, data, "%.3f")

            self._grid["link"]["bed_surface__surface_thickness_new_layer"][
                self._id_deep_links
            ] = (
                self._grid["link"]["bed_surface__surface_thickness_new_layer"][
                    self._id_deep_links
                ]
                - self._new_surface_layer_thickness
            )
            self._grid["link"]["topographic__elevation_subsurface"][
                self._id_deep_links
            ] = (
                self._grid["link"]["topographic__elevation_subsurface"][
                    self._id_deep_links
                ]
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
                if link_data.shape[0] > 1:
                    self._grid["link"]["bed_subsurface__grain_size_distribution"][
                        i, :
                    ] = link_data[
                        link_data.shape[0] - 2, 9:None
                    ]  # 9 is the number of elements before Fs GSD
                    os.remove("link_" + str(i) + ".txt")
                    link_data[-2, 0:2] = np.array((self.t, z[i]))  # Removes last layer
                    with open("link_" + str(i) + ".txt", "ab") as f:
                        np.savetxt(f, link_data[0:-1, :], "%.3f")
                    self._grid["link"]["bed_surface__grain_size_distribution"][
                        i, :
                    ] = self._grid["link"]["bed_subsurface__grain_size_distribution"][
                        i, :
                    ]

            self._grid["link"]["bed_surface__surface_thickness_new_layer"][
                self._id_eroded_links
            ] = (
                self._grid["link"]["bed_surface__surface_thickness_new_layer"][
                    self._id_eroded_links
                ]
                + self._new_surface_layer_thickness
            )
            self._grid["link"]["topographic__elevation_subsurface"][
                self._id_eroded_links
            ] = (
                self._grid["link"]["topographic__elevation_subsurface"][
                    self._id_eroded_links
                ]
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
        the magnitude at a given node."""

        qb_x_r = (
            self._grid["link"]["sediment_transport__bedload_rate"][
                self._grid.links_at_node[:, 0]
            ]
            * self._grid.dy
        )
        (I,) = np.where(self._grid.links_at_node[:, 0] < 0)
        qb_x_r[I] = 0

        qb_x_l = (
            self._grid["link"]["sediment_transport__bedload_rate"][
                self._grid.links_at_node[:, 2]
            ]
            * self._grid.dy
        )
        (I,) = np.where(self._grid.links_at_node[:, 2] < 0)
        qb_x_l[I] = 0

        qb_y_u = (
            self._grid["link"]["sediment_transport__bedload_rate"][
                self._grid.links_at_node[:, 1]
            ]
            * self._grid.dx
        )
        (I,) = np.where(self._grid.links_at_node[:, 1] < 0)
        qb_y_u[I] = 0

        qb_y_l = (
            self._grid["link"]["sediment_transport__bedload_rate"][
                self._grid.links_at_node[:, 3]
            ]
            * self._grid.dx
        )
        (I,) = np.where(self._grid.links_at_node[:, 3] < 0)
        qb_y_l[I] = 0

        qb_x = 0.5 * (qb_x_r + qb_x_l)
        qb_y = 0.5 * (qb_y_u + qb_y_l)

        self._bedloadRate_vector = np.transpose(np.vstack((qb_x, qb_y)))
        self._bedloadRate_magnitude = np.sqrt(qb_x**2 + qb_y**2)

    def vector_mapper(self, vector):
        """Maps a vector, in this case shear stress or velocity, values from
        links onto nodes preserving the components.

        This method takes the vectors values on links and determines the
        vectors components in nodes.

        Examples
        --------
        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> import copy
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab import imshow_grid
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link
        >>> from matplotlib import pyplot as plt

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node['bed_surface__grain_size_distribution_location'] = np.array([
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,])

        >>> grid.at_node['topographic__elevation'] = np.array([
        ... 1.07, 1.06, 1.00, 1.06, 1.07,
        ... 1.08, 1.07, 1.03, 1.07, 1.08,
        ... 1.09, 1.08, 1.07, 1.08, 1.09,
        ... 1.09, 1.09, 1.08, 1.09, 1.09,
        ... 1.09, 1.09, 1.09, 1.09, 1.09,])

        >>> grid.at_node['topographic__elevation_original'] = \
            copy.deepcopy(grid.at_node['topographic__elevation'])

        >>> grid.set_watershed_boundary_condition(grid.at_node['topographic__elevation'])

        >>> grid.at_node['surface_water__depth'] = np.full(grid.number_of_nodes,0.102)

        >>> grid.at_node['surface_water__velocity'] = np.array([
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,])

        >>> grid.at_node['surface_water__velocity_previous_time'] = \
            grid.at_node['surface_water__velocity']

        >>> grid['link']['surface_water__depth'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__depth')
        >>> grid['link']['surface_water__velocity'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity')
        >>> grid['link']['surface_water__velocity_previous_time'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity_previous_time')

        >>> gsd = np.array([[32, 100, 100], [16, 25, 50], [8, 0, 0]])

        >>> timeStep = 1 # time step in seconds

        >>> RBD = RiverBedDynamics(grid, gsd = gsd, dt = timeStep, bedload_equation = 'Parker1990')
        >>> RBD.run_one_step()

        If we want to plot the velocity vector on top of the surface water
        depth we can do:

        >>> (velocityVector, velocityMagnitude) = RBD.vector_mapper(RBD._u)
        >>> velocityVector_x = velocityVector[:,0]
        >>> velocityVector_x
        array([ 0.125,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,
        0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,
        0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,
        0.125])

        >>> velocityVector_y = velocityVector[:,1]
        >>> velocityVector_y
        array([ 0.125,  0.125,  0.125,  0.125,  0.125,  0.25 ,  0.25 ,  0.25 ,
        0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.25 ,
        0.25 ,  0.25 ,  0.25 ,  0.25 ,  0.125,  0.125,  0.125,  0.125,
        0.125])

        # A figure can be visualized using:
        # >>> imshow_grid(grid, 'surface_water__depth',cmap='Blues',vmin=0,vmax=0.5)
        # >>> plt.quiver(grid.x_of_node,grid.y_of_node,velocityVector_x,velocityVector_y,scale = 10)
        # >>> plt.show(block=False)

        """
        vector_x_r = vector[self._grid.links_at_node[:, 0]]
        (I,) = np.where(self._grid.links_at_node[:, 0] < 0)
        vector_x_r[I] = 0

        vector_x_l = vector[self._grid.links_at_node[:, 2]]
        (I,) = np.where(self._grid.links_at_node[:, 2] < 0)
        vector_x_l[I] = 0

        vector_y_u = vector[self._grid.links_at_node[:, 1]]
        (I,) = np.where(self._grid.links_at_node[:, 1] < 0)
        vector_y_u[I] = 0

        vector_y_l = vector[self._grid.links_at_node[:, 3]]
        (I,) = np.where(self._grid.links_at_node[:, 3] < 0)
        vector_y_l[I] = 0

        vector_x = 0.5 * (vector_x_r + vector_x_l)
        vector_y = 0.5 * (vector_y_u + vector_y_l)

        return np.transpose(np.vstack((vector_x, vector_y))), np.sqrt(
            vector_x**2 + vector_y**2
        )

    def calculate_DX(self, fX, mapped_in="link"):
        """ Calculates the bed surface and bed load grain size corresponding to
        any fraction. For example, 50%, which is the D50 or median grain size.
        In that case use fX = 0.5

        This method takes the user specified fraction, from 0 to 1, and outputs
        the corresponding grain size in nodes or links. By default the link
        option is used. Use mapped_in = 'node' to calculate at nodes

        Examples
        --------

        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> import copy
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab import imshow_grid
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node['bed_surface__grain_size_distribution_location'] = np.array([
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,])

        >>> grid.at_node['topographic__elevation'] = np.array([
        ... 1.07, 1.06, 1.00, 1.06, 1.07,
        ... 1.08, 1.07, 1.03, 1.07, 1.08,
        ... 1.09, 1.08, 1.07, 1.08, 1.09,
        ... 1.09, 1.09, 1.08, 1.09, 1.09,
        ... 1.09, 1.09, 1.09, 1.09, 1.09,])

        >>> grid.at_node['topographic__elevation_original'] = \
            copy.deepcopy(grid.at_node['topographic__elevation'])

        >>> grid.set_watershed_boundary_condition(grid.at_node['topographic__elevation'])

        >>> grid.at_node['surface_water__depth'] = np.full(grid.number_of_nodes,0.102)

        >>> grid.at_node['surface_water__velocity'] = np.array([
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,])

        >>> grid.at_node['surface_water__velocity_previous_time'] = \
            grid.at_node['surface_water__velocity']

        >>> grid['link']['surface_water__depth'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__depth')
        >>> grid['link']['surface_water__velocity'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity')
        >>> grid['link']['surface_water__velocity_previous_time'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity_previous_time')

        >>> gsd = np.array([[32, 100, 100], [16, 25, 50], [8, 0, 0]])

        >>> timeStep = 1 # time step in seconds

        >>> RBD = RiverBedDynamics(grid, gsd = gsd, dt = timeStep, bedload_equation = 'Parker1990')
        >>> RBD.run_one_step()

        >>> RBD.calculate_DX(0.9)
        (array([ 28.72809865,  27.85762504,  27.85762504,  28.72809865,
         29.17511936,  27.85761803,  27.84919993,  27.85761803,
         29.17511936,  28.64080227,  27.85762504,  27.85762504,
         28.64080227,  29.17511991,  27.85761803,  27.86610585,
         27.85761803,  29.17511991,  28.64080227,  27.85761803,
         27.85761803,  28.64080227,  29.17511991,  27.85761803,
         27.85761803,  27.85761803,  29.17511991,  28.64080227,
         27.85761803,  27.85761803,  28.64080227,  29.17511963,
         27.85761803,  27.85762504,  27.85761803,  29.17511963,
         28.64080227,  27.85761803,  27.85761803,  28.64080227]),
         array([  0.        ,  27.73672675,  27.73672675,   0.        ,
          0.        ,   0.        ,   0.        ,   0.        ,
          0.        ,   0.        ,  27.65525847,  27.65525847,
          0.        ,  28.14885923,  25.03215792,  27.65525847,
         25.03215792,  28.14885923,   0.        ,  25.03215792,
         25.03215792,   0.        ,   0.        ,  25.03215792,
         25.03215792,  25.03215792,   0.        ,   0.        ,
         25.03215792,  25.03215792,   0.        ,   0.        ,
          0.        ,   0.        ,   0.        ,   0.        ,
          0.        ,   0.        ,   0.        ,   0.        ]))

        """
        nodes_list = np.arange(0, self._grid.number_of_nodes)
        number_links = self._grid.number_of_links
        links_list = np.arange(0, number_links)

        grain_size_D = self._grain_size_D_original  # Grain sizes
        grain_size_Psi_scale_D = np.flip(
            np.log2(grain_size_D), axis=0
        )  # Grain sizes in Psi scale

        # First for the bed surface

        # Equivalent grain sizes frequency
        if mapped_in == "link":
            grain_size_D_equivalent_frequency = self._grid["link"][
                "bed_surface__grain_size_distribution"
            ]
            frequency_grain_size_list = links_list
        else:
            grain_size_D_equivalent_frequency = self._grid["node"][
                "bed_surface__grain_size_distribution"
            ]
            frequency_grain_size_list = nodes_list

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equivalent_frequency,
                np.zeros([grain_size_D_equivalent_frequency.shape[0], 1]),
            )
        )  # Grain sizes freq

        grain_size_D_equivalent_frequency_cumulative = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= fX, axis=1) - 1
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > fX, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i1
                ]
                - grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i0
                ]
            )
        ) * (
            fX
            - grain_size_D_equivalent_frequency_cumulative[
                frequency_grain_size_list, i0
            ]
        )
        DX_surface = 2**grain_size_Psi_scale_DX

        # Now for the bed load

        # Equivalent grain sizes frequency

        if mapped_in == "link":
            grain_size_D_equivalent_frequency = self._grid["link"][
                "sediment_transport__bedload_grain_size_distribution"
            ]
        else:
            grain_size_D_equivalent_frequency = self._grid["node"][
                "sediment_transport__bedload_grain_size_distribution"
            ]

        grain_size_frequency = np.hstack(
            (
                grain_size_D_equivalent_frequency,
                np.zeros([grain_size_D_equivalent_frequency.shape[0], 1]),
            )
        )  # Grain sizes freq

        grain_size_D_equivalent_frequency_cumulative = np.cumsum(
            np.flip(grain_size_frequency, axis=1), axis=1
        )  # Cumulative GSD in each node
        # Finds the index of the grain size smaller (i0) and larger (i1) than
        # the fraction requested (fX) in each node
        i0 = np.argmin(grain_size_D_equivalent_frequency_cumulative <= fX, axis=1) - 1
        i1 = np.argmax(grain_size_D_equivalent_frequency_cumulative > fX, axis=1)

        grain_size_Psi_scale_DX = grain_size_Psi_scale_D[i0] + (
            (grain_size_Psi_scale_D[i1] - grain_size_Psi_scale_D[i0])
            / (
                grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i1
                ]
                - grain_size_D_equivalent_frequency_cumulative[
                    frequency_grain_size_list, i0
                ]
            )
        ) * (
            fX
            - grain_size_D_equivalent_frequency_cumulative[
                frequency_grain_size_list, i0
            ]
        )
        DX_bedload = 2**grain_size_Psi_scale_DX

        return DX_surface, DX_bedload

    def formats_gsd_outputs(self, bedload_gsd):
        """ Gives a more friendly format for the bed surface or bed load GSD.
        Reads a bed load GSD, from links or nodes, and returns the GSD in
        cumulative percetage

        Examples
        --------

        Let us copy part of the example from the beginning

        >>> import numpy as np
        >>> import copy
        >>> from landlab import RasterModelGrid
        >>> from landlab.components import RiverBedDynamics
        >>> from landlab import imshow_grid
        >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link

        >>> grid = RasterModelGrid((5, 5))

        >>> grid.at_node['bed_surface__grain_size_distribution_location'] = np.array([
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,
        ... 0, 1., 1., 1., 0,])

        >>> grid.at_node['topographic__elevation'] = np.array([
        ... 1.07, 1.06, 1.00, 1.06, 1.07,
        ... 1.08, 1.07, 1.03, 1.07, 1.08,
        ... 1.09, 1.08, 1.07, 1.08, 1.09,
        ... 1.09, 1.09, 1.08, 1.09, 1.09,
        ... 1.09, 1.09, 1.09, 1.09, 1.09,])

        >>> grid.at_node['topographic__elevation_original'] = \
            copy.deepcopy(grid.at_node['topographic__elevation'])

        >>> grid.set_watershed_boundary_condition(grid.at_node['topographic__elevation'])

        >>> grid.at_node['surface_water__depth'] = np.full(grid.number_of_nodes,0.102)

        >>> grid.at_node['surface_water__velocity'] = np.array([
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,
        ... 0.25, 0.25, 0.25, 0.25, 0.25,])

        >>> grid.at_node['surface_water__velocity_previous_time'] = \
            grid.at_node['surface_water__velocity']

        >>> grid['link']['surface_water__depth'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__depth')
        >>> grid['link']['surface_water__velocity'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity')
        >>> grid['link']['surface_water__velocity_previous_time'] = \
            map_mean_of_link_nodes_to_link(grid,'surface_water__velocity_previous_time')

        >>> gsd = np.array([[32, 100, 100], [16, 25, 50], [8, 0, 0]])

        >>> timeStep = 1 # time step in seconds

        >>> RBD = RiverBedDynamics(grid, gsd = gsd, dt = timeStep, bedload_equation = 'Parker1990')
        >>> RBD.run_one_step()
        >>> RBD.formats_gsd_outputs(grid['link']['sediment_transport__bedload_grain_size_distribution'])
                   8         16     32
        Link_0   0.0   0.000000    0.0
        Link_1   0.0  51.520878  100.0
        Link_2   0.0  51.520878  100.0
        Link_3   0.0   0.000000    0.0
        Link_4   0.0   0.000000    0.0
        Link_5   0.0   0.000000    0.0
        Link_6   0.0   0.000000    0.0
        Link_7   0.0   0.000000    0.0
        Link_8   0.0   0.000000    0.0
        Link_9   0.0   0.000000    0.0
        Link_10  0.0  52.498142  100.0
        Link_11  0.0  52.498142  100.0
        Link_12  0.0   0.000000    0.0
        Link_13  0.0  45.944616  100.0
        Link_14  0.0  71.774474  100.0
        Link_15  0.0  52.498142  100.0
        Link_16  0.0  71.774474  100.0
        Link_17  0.0  45.944616  100.0
        Link_18  0.0   0.000000    0.0
        Link_19  0.0  71.774474  100.0
        Link_20  0.0  71.774474  100.0
        Link_21  0.0   0.000000    0.0
        Link_22  0.0   0.000000    0.0
        Link_23  0.0  71.774474  100.0
        Link_24  0.0  71.774474  100.0
        Link_25  0.0  71.774474  100.0
        Link_26  0.0   0.000000    0.0
        Link_27  0.0   0.000000    0.0
        Link_28  0.0  71.774474  100.0
        Link_29  0.0  71.774474  100.0
        Link_30  0.0   0.000000    0.0
        Link_31  0.0   0.000000    0.0
        Link_32  0.0   0.000000    0.0
        Link_33  0.0   0.000000    0.0
        Link_34  0.0   0.000000    0.0
        Link_35  0.0   0.000000    0.0
        Link_36  0.0   0.000000    0.0
        Link_37  0.0   0.000000    0.0
        Link_38  0.0   0.000000    0.0
        Link_39  0.0   0.000000    0.0

        It also works for nodes.

        """

        if bedload_gsd.shape[0] == self._grid.number_of_links:
            indexText = "Link_"
        else:
            indexText = "Node_"

        bedload_gsd = np.hstack((bedload_gsd, np.zeros([bedload_gsd.shape[0], 1])))
        bedload_gsd = np.cumsum(np.fliplr(bedload_gsd), axis=1) * 100

        D_ascending = np.sort(self._grain_size_D_original)
        columns = ["".join(item) for item in D_ascending.astype(str)]

        df = pd.DataFrame(
            bedload_gsd,
            columns=columns,
            index=[indexText + str(i) for i in range(bedload_gsd.shape[0])],
        )
        return df
