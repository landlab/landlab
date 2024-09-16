"""GroundwaterDupuitPercolator Component.

@author: G Tucker, D Litwin, K Barnhart
"""

from warnings import warn

import numpy as np

from landlab import Component
from landlab.grid.mappers import map_mean_of_link_nodes_to_link
from landlab.grid.mappers import map_value_at_max_node_to_link
from landlab.utils import return_array_at_link
from landlab.utils import return_array_at_node


# regularization functions used to deal with numerical demons of seepage
def _regularize_G(u, reg_factor):
    """Smooths transition of step function with an exponential.

    0<=u<=1.
    """
    return np.exp(-(1 - u) / reg_factor)


def _regularize_R(u):
    """ramp function on u."""
    return u * np.greater_equal(u, 0)


def _update_thickness(dt, h0, b, f, dqdx, n, r):
    """analytical solution for the linearized governing equation."""
    out = b * (
        1
        - r
        * np.log(
            1
            + np.exp((1 - (h0 + ((f - dqdx) * dt) / n) / b) / r)
            * (1 - np.exp(-(b - h0) / (b * r)))
        )
    )
    out[f <= dqdx] = (h0 + (1 / n * (f - dqdx)) * dt)[f <= dqdx]
    return out


def get_link_hydraulic_conductivity(grid, K):
    """Returns array of hydraulic conductivity on links, allowing for aquifers
    with laterally anisotropic hydraulic conductivity.

    Parameters
    ----------
    K: (2x2) array of floats (m/s)
        The hydraulic conductivity tensor:
        [[Kxx, Kxy],[Kyx,Kyy]]
    """

    u = grid.unit_vector_at_link
    K_link = np.zeros(len(u))
    for i in range(len(u)):
        K_link[i] = np.dot(np.dot(u[i, :], K), u[i, :])
    return K_link


class GroundwaterDupuitPercolator(Component):
    r"""
    Simulate groundwater flow in a shallow unconfined aquifer.

    The GroundwaterDupuitPercolator solves the Boussinesq equation for
    flow in an unconfined aquifer over an impermeable aquifer base and
    calculates groundwater return flow to the surface. This method uses the
    Dupuit-Forcheimer approximation. This means that the model assumes the
    aquifer is laterally extensive in comparison to its thickness, such that
    the vertical component of flow is negligible. It also assumes that the
    capillary fringe is small, such that the water table can be modeled as a
    free surface. Please consider the applicability of these assumptions when
    using this model. For more details, see component documentation
    :ref:`here <dupuit_theory>`.

    Examples
    --------
    Import the grid class and component

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import GroundwaterDupuitPercolator

    Initialize the grid and component

    >>> grid = RasterModelGrid((10, 10), xy_spacing=10.0)
    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> abe = grid.add_zeros("aquifer_base__elevation", at="node")
    >>> elev[:] = 5.0
    >>> gdp = GroundwaterDupuitPercolator(grid)

    Run component forward.

    >>> dt = 1e4
    >>> for i in range(100):
    ...     gdp.run_one_step(dt)
    ...

    When the model generates groundwater return flow, the surface water flux
    out of the domain can be calculated only after a FlowAccumulator is run.
    Below is a more complex example that demonstrates this case.

    >>> from landlab.components import FlowAccumulator

    Set boundary conditions and initialize grid

    >>> grid = RasterModelGrid((5, 41), xy_spacing=10.0)
    >>> grid.set_closed_boundaries_at_grid_edges(True, True, False, True)

    Make a sloping, 3 m thick aquifer, initially fully saturated

    >>> elev = grid.add_zeros("topographic__elevation", at="node")
    >>> elev[:] = grid.x_of_node / 100 + 3
    >>> base = grid.add_zeros("aquifer_base__elevation", at="node")
    >>> base[:] = grid.x_of_node / 100
    >>> wt = grid.add_zeros("water_table__elevation", at="node")
    >>> wt[:] = grid.x_of_node / 100 + 3

    Initialize components

    >>> gdp = GroundwaterDupuitPercolator(grid, recharge_rate=1e-7)
    >>> fa = FlowAccumulator(grid, runoff_rate="surface_water__specific_discharge")

    Advance timestep. Default units are meters and seconds, though the component
    is unit agnostic.

    >>> dt = 1e3
    >>> for i in range(1000):
    ...     gdp.run_one_step(dt)
    ...

    Calculate surface water flux out of domain

    >>> fa.run_one_step()
    >>> np.testing.assert_almost_equal(gdp.calc_sw_flux_out(), 0.0005077)


    Notes
    -----
    Below is a summary of the theory and numerical implementation of
    the ``GroundwaterDupuitPercolator``. A complete description can be found
    :ref:`here <dupuit_theory>`.

    Groundwater discharge per unit length, :math:`q`, is calculated as:

    .. math::
        q = -K_{sat} h \big( \nabla z \big) \cos^2 (\alpha)

    where :math:`K_{sat}` is the saturated hydraulic conductivity, :math:`h` is
    the aquifer thickness, and :math:`\alpha` is the slope angle of the aquifer base.

    Surface water discharge per unit area, :math:`q_s`, is calculated as:

    .. math::
        q_s = \mathcal{G}_r \bigg( \frac{h}{d} \bigg) \mathcal{R} \big(-\nabla \cdot q + f \big)

    where :math:`\mathcal{G}_r` is a smoothed step function, :math:`\mathcal{R}`
    is the ramp function, :math:`d` is the regolith thickness, and :math:`f` is
    the recharge rate.

    The evolution of aquifer thickness is then given by:

    .. math::
        n \frac{\partial h}{\partial t} = f - q_s - \nabla \cdot q

    where :math:`n` is the drainable porosity.

    A semi-analytical approach is used to update aquifer thickness :math:`h`, in which
    the differential equation above is linearized by assuming :math:`\nabla \cdot q` varies
    minimally over the duration of a timestep. In this approach, :math:`q` is
    reported at the beginning of the timestep, while :math:`q_s` is reported at the
    end of the timestep.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Litwin, D. G., Tucker, G.E., Barnhart, K. R., Harman, C. J. (2020).
    GroundwaterDupuitPercolator: A Landlab component for groundwater flow.
    Journal of Open Source Software, 5(46), 1935, https://doi.org/10.21105/joss.01935.

    **Additional References**

    Marçais, J., de Dreuzy, J. R. & Erhel, J. Dynamic coupling of subsurface
    and seepage flows solved within a regularized partition formulation.
    Advances in Water Resources 109, 94–105 (2017).

    Childs, E. C. Drainage of Groundwater Resting on a Sloping Bed. Water
    Resources Research 7, 1256–1263 (1971).
    """

    _name = "GroundwaterDupuitPercolator"

    _cite_as = """@article{litwin2020groundwater,
      doi = {10.21105/joss.01935},
      url = {https://doi.org/10.21105/joss.01935},
      year = {2020},
      publisher = {The Open Journal},
      volume = {5},
      number = {46},
      pages = {1935},
      author = {David Litwin and Gregory Tucker and Katherine Barnhart and Ciaran Harman},
      title = {GroundwaterDupuitPercolator: A Landlab component for groundwater flow},
      journal = {Journal of Open Source Software}
    }"""

    _unit_agnostic = True

    _info = {
        "aquifer__thickness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "thickness of saturated zone",
        },
        "aquifer_base__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of impervious layer",
        },
        "aquifer_base__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of the aquifer base in the link direction",
        },
        "average_surface_water__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "average surface water specific discharge over variable timesteps",
        },
        "groundwater__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m2/s",
            "mapping": "link",
            "doc": "discharge per width in link dir",
        },
        "groundwater__velocity": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "link",
            "doc": "velocity of groundwater in link direction",
        },
        "hydraulic__gradient": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/m",
            "mapping": "link",
            "doc": "gradient of water table in link direction",
        },
        "surface_water__specific_discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m/s",
            "mapping": "node",
            "doc": "rate of seepage to surface",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water_table__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "elevation of water table",
        },
    }

    def __init__(
        self,
        grid,
        hydraulic_conductivity=0.001,
        porosity=0.2,
        recharge_rate=1.0e-8,
        regularization_f=1e-2,
        courant_coefficient=0.5,
        vn_coefficient=0.8,
        callback_fun=lambda *args, **kwargs: None,
        **callback_kwds,
    ):
        r"""
        Parameters
        ----------
        grid: ModelGrid
            Landlab ModelGrid object
        hydraulic_conductivity: float, field name, array of float or function.
            the aquifer saturated hydraulic conductivity, m/s.
            If function is given, it should take a landlab ModelGrid and return
            an array of floats at link. This may be used if the lateral hydraulic
            conductivity is not vertically homogenous and the effective hydraulic
            conductivity needs to be modified based upon on the position of the
            water table. See component tests for example.
            Default = 0.001 m/s
        porosity: float, field name or array of float
            the drainable porosity of the aquifer [-]
            Default = 0.2
        recharge_rate: float, field name, or array of float
            Rate of recharge, m/s
            Default = 1.0e-8 m/s
        regularization_f: float
            factor controlling the smoothness of the transition between
            surface and subsurface flow
            Default = 0.01
        courant_coefficient: float (-)
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over groundwater flow
            velocity. This parameter is only used with
            ``run_with_adaptive_time_step_solver`` and must be greater than
            zero.
            Default = 0.5
        vn_coefficient: float (-)
            The multiplying factor C for the condition :math:`dt >= C*dx^2/(4D)`,
            where :math:`D = Kh/n` is the diffusivity of the Boussinesq
            equation. This arises from a von Neumann stability analysis of
            the Boussinesq equation when the hydraulic gradient is small.
            This parameter is only used with ``run_with_adaptive_time_step_solver``
            and must be greater than zero.
            Default = 0.8
        callback_fun: function(grid, recharge_rate, substep_dt, \*\*kwargs)
            Optional function that will be executed at the end of each sub-timestep
            in the run_with_adaptive_time_step_solver method. Intended purpose
            is to write output not otherwise visible outside of the method call.
            The function should have three required arguments:
            grid: the ModelGrid instance used by GroundwaterDupuitPercolator
            recharge_rate: an array at node that is the specified recharge rate
            substep_dt: the length of the current substep determined internally
            by run_with_adaptive_time_step_solver to meet stability criteria.
        callback_kwds: any additional keyword arguments for the provided callback_fun.
        """
        super().__init__(grid)

        # Shorthand
        self._cores = grid.core_nodes

        # Create fields:

        self._elev = self._grid.at_node["topographic__elevation"]
        self._base = self._grid.at_node["aquifer_base__elevation"]

        self.initialize_output_fields()

        self._wtable = self._grid.at_node["water_table__elevation"]
        self._wtable[grid.closed_boundary_nodes] = 0

        self._thickness = self._grid.at_node["aquifer__thickness"]
        self._thickness[:] = self._wtable - self._base
        self._thickness[grid.closed_boundary_nodes] = 0

        self._hydr_grad = self._grid.at_link["hydraulic__gradient"]
        self._base_grad = self._grid.at_link["aquifer_base__gradient"]

        self._q = self._grid.at_link["groundwater__specific_discharge"]
        self._qs = self._grid.at_node["surface_water__specific_discharge"]
        self._qsavg = self.grid.at_node["average_surface_water__specific_discharge"]

        self._vel = self._grid.at_link["groundwater__velocity"]

        # Convert parameters to fields if needed, and store a reference
        self.K = hydraulic_conductivity
        self.recharge = recharge_rate
        self.n = porosity
        self._r = regularization_f

        # save courant_coefficient (and test)
        self.courant_coefficient = courant_coefficient
        self.vn_coefficient = vn_coefficient

        # set callback function
        self._callback_kwds = callback_kwds
        self.callback_fun = callback_fun

    @property
    def callback_fun(self):
        r"""callback function for adaptive timestep solver

        Parameters
        ----------
        callback_fun: function(grid, recharge_rate, substep_dt, \*\*callback_kwds)
            Optional function that will be executed at the end of each sub-timestep
            in the run_with_adaptive_time_step_solver method. Intended purpose
            is to write output not otherwise visible outside of the method call.
            The function should have three required arguments:
            grid: the ModelGrid instance used by GroundwaterDupuitPercolator
            recharge_rate: an array at node that is the specified recharge rate
            substep_dt: the length of the current substep determined internally
            by run_with_adaptive_time_step_solver to meet stability criteria.
        """
        return self._callback_fun

    @callback_fun.setter
    def callback_fun(self, new_val):
        try:  # New style callback function.
            new_val(self._grid, self.recharge, 0.0, **self._callback_kwds)
            self._callback_fun = new_val
        except TypeError as exc:  # Nonfunctional callback function.
            raise ValueError(
                f"{str(exc)}: Please supply a callback function with the form "
                "function(grid, recharge_rate, substep_dt, **kwargs)"
            ) from exc

    @property
    def courant_coefficient(self):
        """Courant coefficient for adaptive time step.

        Parameters
        ----------
        courant_coefficient: float (-)
            The muliplying factor on the condition that the timestep is
            smaller than the minimum link length over groundwater flow
            velocity. This parameter is only used with
            ``run_with_adaptive_time_step_solver`` and must be greater than
            zero.
        """
        return self._courant_coefficient

    @courant_coefficient.setter
    def courant_coefficient(self, new_val):
        if new_val <= 0:
            raise ValueError("courant_coefficient must be > 0.")
        self._courant_coefficient = new_val

    @property
    def vn_coefficient(self):
        """Coefficient for the diffusive timestep condition in
        the adaptive timestep solver.

        Parameters
        ----------
        vn_coefficient: float (-)
            The multiplying factor C for the condition dt >= C*dx^2/(4D),
            where D = Kh/n is the diffusivity of the Boussinesq
            equation. This arises from a von Neumann stability analysis of
            the Boussinesq equation when the hydraulic gradient is small.
            This parameter is only used with ``run_with_adaptive_time_step_solver``
            and must be greater than zero.
        """
        return self._vn_coefficient

    @vn_coefficient.setter
    def vn_coefficient(self, new_val):
        """set coefficient for the diffusive timestep condition in
        the adaptive timestep solver."""
        if new_val <= 0:
            raise ValueError("vn_coefficient must be > 0.")
        self._vn_coefficient = new_val

    @property
    def K(self):
        """hydraulic conductivity at link (m/s)"""
        if self._kfunc:
            self._K = return_array_at_link(self._grid, self._func(self._grid))
        return self._K

    @K.setter
    def K(self, new_val):
        """set hydraulic conductivity at link (m/s)"""
        if callable(new_val):
            self._kfunc = True

            if (
                not isinstance(new_val(self._grid), np.ndarray)
                and len(new_val(self._grid)) == self._grid.number_of_links
            ):
                raise TypeError(
                    "If a function is provided it must take a ModelGrid and "
                    "return an array of length number_of_links."
                )
            else:
                self._func = new_val
                self._K = return_array_at_link(self._grid, self._func(self._grid))
        else:
            self._kfunc = False
            self._K = return_array_at_link(self._grid, new_val)

    @property
    def recharge(self):
        """recharge rate (m/s)"""
        return self._recharge

    @recharge.setter
    def recharge(self, new_val):
        """set recharge rate (m/s)"""
        self._recharge = return_array_at_node(self._grid, new_val)

    @property
    def n(self):
        """drainable porosity of the aquifer (-)"""
        return self._n

    @n.setter
    def n(self, new_val):
        """set aquifer drainable porosity (-)"""
        self._n = return_array_at_node(self._grid, new_val)
        self._n_link = map_mean_of_link_nodes_to_link(self._grid, self._n)

    @property
    def number_of_substeps(self):
        """
        The number of substeps used by the run_with_adaptive_time_step_solver
        method in the latest method call.
        """
        if self._num_substeps:
            return self._num_substeps
        else:
            warn("The method run_with_adaptive_time_step_solver has not been used")

        return self._num_substeps

    def calc_recharge_flux_in(self):
        """Calculate flux into the domain from recharge.

        Includes recharge that may immediately become saturation excess
        overland flow. (m3/s)
        """
        return np.sum(
            self._grid.cell_area_at_node[self._cores] * self._recharge[self._cores]
        )

    def calc_gw_flux_out(self):
        """Groundwater flux through open boundaries may be positive (out of the
        domain) or negative (into the domain).

        This function determines the correct sign for specific discharge
        based upon this convention, and sums the flux across the
        boundary faces. (m3/s)
        """
        # get links at open boundary nodes
        open_nodes = self._grid.open_boundary_nodes
        links_at_open = self._grid.links_at_node[open_nodes]
        link_dirs_at_open = self._grid.active_link_dirs_at_node[open_nodes]

        # find active links at open boundary nodes
        active_links_at_open = links_at_open[link_dirs_at_open != 0]
        active_link_dirs_at_open = link_dirs_at_open[link_dirs_at_open != 0]

        # get groundwater specific discharge at these locations
        q_at_open_links = self._grid.at_link["groundwater__specific_discharge"][
            active_links_at_open
        ]

        # get cell widths at these locations
        faces = self._grid.face_at_link[active_links_at_open]
        face_widths = self._grid.length_of_face[faces]

        # get volume flux out at these locations
        gw_volume_flux_rate_out = (
            q_at_open_links * active_link_dirs_at_open * face_widths
        )

        return np.sum(gw_volume_flux_rate_out)

    def calc_sw_flux_out(self):
        """Surface water flux out of the domain through seepage and saturation
        excess.

        Note that model does not allow for reinfiltration.  (m3/s)
        """
        return np.sum(
            self._grid.at_node["surface_water__discharge"][
                self._grid.open_boundary_nodes
            ]
        )

    def calc_gw_flux_at_node(self):
        """Calculate the sum of the groundwater flux leaving a node.

        (m2/s)
        """
        gw = (
            self._grid.at_link["groundwater__specific_discharge"][
                self._grid.links_at_node
            ]
            * self._grid.link_dirs_at_node
        )
        gw_out = -np.sum(gw * (gw < 0), axis=1)
        return gw_out

    def calc_total_storage(self):
        """calculate the current water storage in the aquifer (m3)"""
        return np.sum(
            self._n[self._cores]
            * self._grid.cell_area_at_node[self._cores]
            * self._grid.at_node["aquifer__thickness"][self._cores]
        )

    def run_one_step(self, dt):
        """Advance component by one time step of size dt.

        Parameters
        ----------
        dt: float
            The imposed timestep.
        """

        # check water table above surface
        if (self._wtable > self._elev).any():
            warn(
                "water table above elevation surface. "
                "Setting water table elevation here to "
                "elevation surface"
            )
            self._wtable[self._wtable > self._elev] = self._elev[
                self._wtable > self._elev
            ]
            self._thickness[self._cores] = (self._wtable - self._base)[self._cores]

        # Calculate base gradient
        self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._base
        )[self._grid.active_links]
        cosa = np.cos(np.arctan(self._base_grad))

        # Calculate hydraulic gradient
        self._hydr_grad[self._grid.active_links] = (
            self._grid.calc_grad_at_link(self._wtable) * cosa
        )[self._grid.active_links]

        # Calculate groundwater velocity
        self._vel[:] = -self._K * self._hydr_grad

        # Aquifer thickness at links (upwind)
        hlink = (
            map_value_at_max_node_to_link(
                self._grid, "water_table__elevation", "aquifer__thickness"
            )
            * cosa
        )

        # Calculate specific discharge
        self._q[:] = hlink * self._vel

        # Groundwater flux divergence
        dqdx = self._grid.calc_flux_div_at_node(self._q)

        # Regolith thickness
        reg_thickness = self._elev - self._base

        # update thickness from analytical
        self._thickness[self._cores] = _update_thickness(
            dt, self._thickness, reg_thickness, self._recharge, dqdx, self._n, self._r
        )[self._cores]
        self._thickness[self._thickness < 0] = 0.0

        # Recalculate water surface height
        self._wtable[:] = self._base + self._thickness

        # Calculate surface discharge at nodes
        self._qs[:] = _regularize_G(
            self._thickness / reg_thickness, self._r
        ) * _regularize_R(self._recharge - dqdx)

    def run_with_adaptive_time_step_solver(self, dt):
        """
        Advance component by one time step of size dt, subdividing the timestep
        into substeps as necessary to meet stability conditions.
        Note this method returns the fluxes at the last substep, but also
        returns a new field, average_surface_water__specific_discharge, that is
        averaged over all subtimesteps. To return state during substeps,
        provide a callback_fun.

        Parameters
        ----------
        dt: float
            The imposed timestep.
        """

        # check water table above surface
        if (self._wtable > self._elev).any():
            warn(
                "water table above elevation surface. "
                "Setting water table elevation here to "
                "elevation surface"
            )
            self._wtable[self._wtable > self._elev] = self._elev[
                self._wtable > self._elev
            ]
            self._thickness[self._cores] = (self._wtable - self._base)[self._cores]

        # Calculate base gradient
        self._base_grad[self._grid.active_links] = self._grid.calc_grad_at_link(
            self._base
        )[self._grid.active_links]
        cosa = np.cos(np.arctan(self._base_grad))

        # Initialize reg_thickness
        reg_thickness = self._elev - self._base

        # Initialize for average surface discharge
        qs_cumulative = np.zeros_like(self._elev)

        # Initialize variable timestep
        remaining_time = dt
        self._num_substeps = 0

        while remaining_time > 0.0:
            # Calculate hydraulic gradient
            self._hydr_grad[self._grid.active_links] = (
                self._grid.calc_grad_at_link(self._wtable) * cosa
            )[self._grid.active_links]

            # Calculate groundwater velocity
            self._vel[:] = -self._K * self._hydr_grad

            # Aquifer thickness at links (upwind)
            hlink = (
                map_value_at_max_node_to_link(
                    self._grid, "water_table__elevation", "aquifer__thickness"
                )
                * cosa
            )

            # Calculate specific discharge
            self._q[:] = hlink * self._vel

            # Groundwater flux divergence
            dqdx = self._grid.calc_flux_div_at_node(self._q)

            # calculate criteria for timestep
            dt_vn = self._vn_coefficient * np.min(
                np.divide(
                    (self._n_link * self._grid.length_of_link**2),
                    (4 * self._K * hlink),
                    where=hlink > 0,
                    out=np.ones_like(self._q) * 1e15,
                )
            )

            dt_courant = self._courant_coefficient * np.min(
                np.divide(
                    self._grid.length_of_link,
                    abs(self._vel / self._n_link),
                    where=abs(self._vel) > 0,
                    out=np.ones_like(self._q) * 1e15,
                )
            )
            substep_dt = min([dt_courant, dt_vn, remaining_time])
            # 0 = courant limited, 1 = vn limited, 2 = not limited
            # print(np.argmin(np.array([self._dt_courant, self._dt_vn, remaining_time])))

            # update thickness from analytical
            self._thickness[self._cores] = _update_thickness(
                substep_dt,
                self._thickness,
                reg_thickness,
                self._recharge,
                dqdx,
                self._n,
                self._r,
            )[self._cores]
            self._thickness[self._thickness < 0] = 0.0

            # Recalculate water surface height
            self._wtable[:] = self._base + self._thickness

            # Calculate surface discharge at nodes
            self._qs[:] = _regularize_G(
                self._thickness / reg_thickness, self._r
            ) * _regularize_R(self._recharge - dqdx)

            # add cumulative sw discharge in substeps
            qs_cumulative += self._qs * substep_dt

            # calculate the time remaining and advance count of substeps
            remaining_time -= substep_dt
            self._num_substeps += 1

            # run callback function if supplied
            self._callback_fun(
                self._grid, self._recharge, substep_dt, **self._callback_kwds
            )

        self._qsavg[:] = qs_cumulative / dt
