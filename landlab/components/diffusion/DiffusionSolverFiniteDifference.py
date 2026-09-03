#! /usr/env/python
"""Simplified Landlab diffusion component 

Created July 2023. Updates / edits made by Isamar Cortes 
"""


import numpy as np

from landlab import Component, FieldError, LinkStatus, NodeStatus, RasterModelGrid

_ALPHA = 0.15  # time-step stability factor
# ^0.25 not restrictive enough at meter scales w S~1 (possible cases)


class DiffusionSolverFiniteDifference(Component):
    """This component implements linear diffusion of a Landlab raster grid.

    


    Component assumes grid does not deform. If the boundary conditions on the
    grid change after component instantiation, be sure to also call
    :func:`updated_boundary_conditions` to ensure these are reflected in the
    component (especially if fixed_links are present).

    The method keyword allows control of the way the solver works. Options
    other than 'simple' will make the component run slower, but give second
    order solutions that incorporate information from more distal nodes (i.e.,
    the diagonals). Note that the option 'resolve_on_patches' can result in
    somewhat counterintuitive behaviour - this option tells the component to
    treat the diffusivity as a field **with directionality to it** (i.e., that
    the diffusivites are defined on links). Thus if all links have the same
    diffusivity value, with this flag active "effective" diffusivities
    at the nodes will be *higher* than this value (by a factor of root 2) as
    the diffusivity at each patch will be the mean vector sum of that at the
    bounding links.

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> mg = RasterModelGrid((9, 9))
    >>> z = mg.add_zeros("topographic__elevation", at="node")
    >>> z.reshape((9, 9))[4, 4] = 1.
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> ld = LinearDiffuser(mg, linear_diffusivity=1.)
    >>> for i in range(1):
    ...     ld.run_one_step(1.)
    >>> np.isclose(z[mg.core_nodes].sum(), 1.)
    True
    >>> mg2 = RasterModelGrid((5, 30))
    >>> z2 = mg2.add_zeros("topographic__elevation", at="node")
    >>> z2.reshape((5, 30))[2, 8] = 1.
    >>> z2.reshape((5, 30))[2, 22] = 1.
    >>> mg2.set_closed_boundaries_at_grid_edges(True, True, True, True)
    >>> kd = mg2.node_x/mg2.node_x.mean()
    >>> ld2 = LinearDiffuser(mg2, linear_diffusivity=kd)
    >>> for i in range(10):
    ...     ld2.run_one_step(0.1)
    >>> z2[mg2.core_nodes].sum() == 2.
    True
    >>> z2.reshape((5, 30))[2, 8] > z2.reshape((5, 30))[2, 22]
    True

    An example using links:

    >>> mg1 = RasterModelGrid((10, 10), xy_spacing=100.)
    >>> mg2 = RasterModelGrid((10, 10), xy_spacing=100.)
    >>> z1 = mg1.add_zeros("topographic__elevation", at="node")
    >>> z2 = mg2.add_zeros("topographic__elevation", at="node")
    >>> dt = 1.
    >>> nt = 10
    >>> kappa_links = mg2.add_ones("surface_water__discharge", at="link")
    >>> kappa_links *= 10000.
    >>> dfn1 = LinearDiffuser(mg1, linear_diffusivity=10000.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='surface_water__discharge')
    >>> for i in range(nt):
    ...     z1[mg1.core_nodes] += 1.
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn1.run_one_step(dt)
    ...     dfn2.run_one_step(dt)
    >>> np.allclose(z1, z2)
    True
    >>> z2.fill(0.)
    >>> dfn2 = LinearDiffuser(mg2, linear_diffusivity='surface_water__discharge',
    ...                       method='resolve_on_patches')
    >>> for i in range(nt):
    ...     z2[mg2.core_nodes] += 1.
    ...     dfn2.run_one_step(dt)
    >>> np.all(z2[mg2.core_nodes] < z1[mg2.core_nodes])
    True

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Culling, W. (1963). Soil Creep and the Development of Hillside Slopes.
    The Journal of Geology  71(2), 127-161. https://dx.doi.org/10.1086/626891

    """

    _name = "DiffusionSolverFiniteDifference"

    _unit_agnostic = True

    _info = {
        "diffused__quantity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/s",
            "mapping": "node",
            "doc": "Volume flux per unit width along links",
        },
        "diffusion__flux": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2/s",
            "mapping": "link",
            "doc": "Volume flux per unit width along links",
        },
       
    }

    def __init__(self, grid, linear_diffusivity=0.01):
        """
        Parameters
        ----------
        grid : RasterGrid
            A grid.
        """
        super().__init__(grid)
        
        self.diffusion_coef = linear_diffusivity

        self.grid.add_zeros("diffusion__flux", at="link")
        
    
    def run_one_step(self, dt):
        qs = self.grid.at_link["diffusion__flux"]
        value_to_diffuse = self.grid.at_node["diffused__quantity"]
        
        time_remaining = dt
        while time_remaining > 0:
            local_dt = self.time_step
            
            if time_remaining - local_dt < 0:
                local_dt = time_remaining
            
            g = self.grid.calc_grad_at_link(value_to_diffuse)
            qs[self.grid.active_links] = -self.diffusion_coef * g[self.grid.active_links]
            dqdx = self.grid.calc_flux_div_at_node(qs)
            value_to_diffuse[self.grid.core_nodes] = value_to_diffuse[self.grid.core_nodes] + (dsdt[self.grid.core_nodes] * local_dt)

            time_remaining -= local_dt
 
    @property
    def time_step(self):
        """Returns internal time-step size (as a property)."""
        return 0.2 * self.grid.dx ** 2 / self.diffusion_coef
