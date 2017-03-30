# -*- coding: utf-8 -*-
"""
This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np
from six.moves import range
from landlab import RasterModelGrid, Component, FieldError, INACTIVE_LINK, \
    CLOSED_BOUNDARY, CORE_NODE
import inspect
from landlab.utils.decorators import use_file_name_or_kwds


class DischargeDiffuser(Component):
    """Diffuse sediment proportional to an implicit water discharge value.

    This class implements Voller, Hobley, and Paola's scheme for sediment
    diffusion, where the diffusivity of the sediment is proportional to the
    local discharge of water. The method works by solving for a potential
    field describing the water discharge at all nodes on the grid, which
    enforces both mass conservation and flow downhill along topographic
    gradients. This routine is designed to construct sediment fans.

    Note that both the water and sediment discharges are calculated together
    within the component.

    The algorithm uses a pair of rules that look like:

        q_water = K * S**c

        q_sed = a * q_water * (S**n - S_crit**n)**m

    where S_crit is a critical slope threshold, and c, a, m, and n are constant
    parameters which reflect the chosen transport law. In this context, note
    that S, q_water, and q_sed are all vector terms.

    As of 28/3/17, this component needs the following enhancements:
        * m, n powers on the slope terms.
        * demonstration of effectiveness on non-flat topos
        * Optional water depth to fill (probably a second component)
        * Forcing to follow a specific flow law
        * Forcing to follow a specific transport law
        * BC control (as of now, all boundaries are open and fixed grad (?))

    The primary method of this class is :func:`run_one_step`.

    Construction::

        DischargeDiffuser(grid, ...)

    Notes
    -----
    This is a "research grade" component, and is subject to dramatic change
    with little warning. No guarantees are made regarding its accuracy or
    utility. It is not recommended for user use yet!

    Parameters
    ----------
    grid : ModelGrid
        A grid.


    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> import numpy as np
    >>> mg = HexModelGrid(4, 6, dx=2., shape='rect', orientation='vertical')
    >>> z = mg.add_zeros('node', 'topographic__elevation')
    >>> Q_in = mg.add_ones('node', 'water__unit_flux_in')
    >>> z += mg.node_y.copy()
    >>> potfr = PotentialityFlowRouter(mg)
    >>> potfr.run_one_step()
    >>> Q_at_core_nodes = np.array(
    ...     [ 17.02012846,  16.88791903,  13.65746194,  14.85578934,
    ...       11.41908145,  11.43630865,   8.95902559,  10.04348075,
    ...        6.28696459,   6.44316089,   4.62478522,   5.29145188])
    >>> np.allclose(mg.at_node['surface_water__discharge'][mg.core_nodes],
    ...             Q_at_core_nodes)
    True
    """
    _name = 'DischargeDiffuser'

    _input_var_names = ('topographic__elevation',
                        'water__discharge_in',
                        'sediment__discharge_in',
                        )

    _output_var_names = ('topographic__elevation',
                         'surface_water__discharge',
                         'flow__potential',
                         )

    _var_units = {'topographic__elevation': 'm',
                  'water__unit_flux_in': 'm/s',
                  'surface_water__discharge': 'm**3/s',
                  'flow__potential': 'm**3/s',
                  'surface_water__depth': 'm',
                  }

    _var_mapping = {'topographic__elevation': 'node',
                    'water__unit_flux_in': 'node',
                    'surface_water__discharge': 'node',
                    'flow__potential': 'node',
                    'surface_water__depth': 'node',
                    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'water__unit_flux_in': (
            'External volume water per area per time input to each node ' +
            '(e.g., rainfall rate)'),
        'surface_water__discharge': (
            'Magnitude of volumetric water flux out of each node'),
        'flow__potential': (
            'Value of the hypothetical field "K", used to force water flux ' +
            'to flow downhill'),
        'surface_water__depth': (
            'If Manning or Chezy specified, the depth of flow in the cell, ' +
            'calculated assuming flow occurs over the whole surface'),
                  }

    _min_slope_thresh = 1.e-24
    # if your flow isn't connecting up, this probably needs to be reduced

    def __init__(self, grid, threshold, flow_law_slope_exponent=1.,
                 sediment_flux_equation='powerlaw', a=1., m=1., n=1.,
                 Cf=25., rho_sed=2700., rho_fl=1000., g=9.81,
                 stab_thresh=1.e-10, **kwds):
        """Initialize flow router.

        Parameters
        ----------
        grid : RasterModelGrid
            The grid must be a raster.
        threshold : float
            The critical threshold below which no transport occurs. If
            sediment_flux_equation == 'powerlaw', this is a slope value.
            If sediment_flux_equation == 'MPM', this is a critical shear
            stress.
        flow_law_slope_exponent : float
            The power, c, on slope in the relation q_w = K * S**c. This is
            1. in Voller's original derivation, which equates to a time-
            averaged sheetflow solution. If <1., this will result in a
            channel closure condition, and can lead to preferential flow on
            grid diagonals.
        sediment_flux_equation : {'powerlaw', 'MPM'}
            The equation form to use for sediment transport. If 'powerlaw',
            the component uses q_sed = a * q_water * (S**n - S_crit**n)**m,
            where a, m, and n are user-settable. If 'MPM', the coefficients
            in that equation are automatically adjusted to mimic the Meyer-
            Peter Muller equation. In this case, 'threshold', above, is in
            units of shear stress, and a check is performed to ensure that
            'flow_law_slope_exponent' == 0.5, the physically reasonable case.
        a : float
            The prefactor for the sed equation, if
            'sediment_flux_equation' == 'powerlaw'.
        m : float
            The power on each slope term for the sed equation, if
            'sediment_flux_equation' == 'powerlaw'.
        n : float
            The power on the thresholded slope value for the sed equation, if
            'sediment_flux_equation' == 'powerlaw'.
        """
        if RasterModelGrid in inspect.getmro(grid.__class__):
            assert grid.number_of_node_rows >= 3
            assert grid.number_of_node_columns >= 3
            self._raster = True
        else:
            self._raster = False

        assert self._raster is True  # ...for now

        self._grid = grid
        self._slope = threshold

        if np.isclose(flow_law_slope_exponent, 1.):
            self._simple_flow = True
        else:
            self._simple_flow = False
            self._c = flow_law_slope_exponent

        if (np.isclose(a, 1.) and np.isclose(m, 1.) and np.isclose(n, 1.)):
            self._simple_powers = True
        else:
            self._simple_powers = False
        if sediment_flux_equation == 'powerlaw':
            self._powerlaw = True
            self._m = m
            self._n = n
            self._a = a
        else:
            assert sediment_flux_equation == 'MPM'
            self._powerlaw = False  # assume MPM for now
            assert self._simple_powers
            self._m = 1.5
            self._n = 2./3.
            self._a = 8. * Cf**0.5 / ((rho_sed/rho_fl) - 1.)
            tautoslope_prefactor = 1./(rho_fl * (g**2. * Cf)**(1./3.))
            # this prefactor needs also to include K**-(2./3.)
            self._slope *= tautoslope_prefactor

        # hacky fix because water__discharge is defined on both links and nodes
        for out_field in self._output_var_names:
            if self._var_mapping[out_field] == 'node':
                try:
                    self.grid.add_zeros(self._var_mapping[out_field],
                                        out_field, dtype=float)
                except FieldError:
                    pass
            else:
                pass
            try:
                self.grid.add_zeros('node', 'surface_water__discharge',
                                    dtype=float)
            except FieldError:
                pass

        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns

        # stab_thresh needs to not vary according to nnodes, so:
        self._stab_thresh = stab_thresh / ni / nj

        self._K = grid.zeros('node', dtype=float)
        self._Knew = grid.zeros('node', dtype=float)
        self._prevK = grid.zeros('node', dtype=float)
        self._znew = grid.zeros('node', dtype=float)
        # discharge across north, south, west, and east face of control volume
        self._Qn = np.zeros((ni, nj), dtype='float')
        self._Qs = np.zeros((ni, nj), dtype='float')
        self._Qw = np.zeros((ni, nj), dtype='float')
        self._Qe = np.zeros((ni, nj), dtype='float')

        # coefficenst used in solition of flow conductivity K
        self._app = np.zeros((ni, nj), dtype='float')
        self._apz = np.zeros((ni, nj), dtype='float')
        self._aww = np.zeros((ni, nj), dtype='float')
        self._awp = np.zeros((ni, nj), dtype='float')
        self._awz = np.zeros((ni, nj), dtype='float')
        self._aee = np.zeros((ni, nj), dtype='float')
        self._aep = np.zeros((ni, nj), dtype='float')
        self._aez = np.zeros((ni, nj), dtype='float')
        self._ass = np.zeros((ni, nj), dtype='float')
        self._asp = np.zeros((ni, nj), dtype='float')
        self._asz = np.zeros((ni, nj), dtype='float')
        self._ann = np.zeros((ni, nj), dtype='float')
        self._anp = np.zeros((ni, nj), dtype='float')
        self._anz = np.zeros((ni, nj), dtype='float')

        self._slx = np.empty((ni, nj), dtype=float)
        self._sly = np.empty((ni, nj), dtype=float)
        self._Qsed_w = np.empty((ni, nj), dtype=float)
        self._Qsed_e = np.empty((ni, nj), dtype=float)
        self._Qsed_n = np.empty((ni, nj), dtype=float)
        self._Qsed_s = np.empty((ni, nj), dtype=float)

        self._dxbyy = grid.dx/grid.dy
        self._dybyx = grid.dy/grid.dx

        # define this slice set for clarity:
        self._all = slice(0, None, 1)
        self._allpad = slice(1, -1, 1)
        self._low = slice(0, -1, 1)
        self._lowpad = slice(0, -2, 1)
        self._high = slice(1, None, 1)
        self._highpad = slice(2, None, 1)
        self._east = (self._all, self._high)
        self._eastpad = (self._allpad, self._highpad)
        self._west = (self._all, self._low)
        self._westpad = (self._allpad, self._lowpad)
        self._south = (self._low, self._all)
        self._southpad = (self._lowpad, self._allpad)
        self._north = (self._high, self._all)
        self._northpad = (self._highpad, self._allpad)
        self._centpad = (self._allpad, self._allpad)
        self._eastedge = (self._all, slice(-1, None, 1))
        self._westedge = (self._all, slice(0, 1, 1))
        self._southedge = (slice(0, 1, 1), self._all)
        self._northedge = (slice(-1, None, 1), self._all)
        self._SWpad = (self._lowpad, self._lowpad)
        self._NEpad = (self._highpad, self._highpad)
        self._SEpad = (self._lowpad, self._highpad)
        self._NWpad = (self._highpad, self._lowpad)

    def run_one_step(self, dt, **kwds):
        """
        """
        grid = self.grid
        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns
        z = grid.at_node['topographic__elevation']
        Qsp = grid.at_node['water__discharge_in'].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns))
        Qsource = grid.at_node['sediment__discharge_in'].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns))

        # make the fields (ni, nj) to enable efficient slicing
        eta = z.reshape((ni, nj))
        K = self._K.reshape((ni, nj))
        Knew = self._Knew.reshape((ni, nj))

        # begin the component stability loop:
        current_internal_time = 0.
        breakcheck = False
        while not breakcheck:
            # pad eta
            pad_eta = np.pad(eta, ((1, 1), (1, 1)), 'edge')

            # must do water part 1st for stab analysis to work OK:
            # do the water routing on links
            # These calculations are based on assuming that the flow is a sheet
            # flow that can be characterized with a potential equation. If this
            # flow is isotropic (an assumption we should revisit) with this
            # model the flow discharge in the x-direction (for example) can be
            # calculated as a constant (K the 'flow conductivity') times the
            # component of the sediment slope in that direction. It helps to
            # define a 'slope velocity' u, with components ustar=-deta/dx and
            # vstar=-deta/dx which allows us to write down the following
            # advection like gov. eq. for the flow  ----div(Ku)+Q=0---where Q
            # represents external flow inputs

            # Since we can readily determine u from the current sediment
            # topography, We solve this equation for K using an upwind scheme

            # Build upwinded coefficients. Vals only 0 if if flow is in upwind
            # direction
            # note cols/rows which don't get updated will always remain as 0,
            # which is right provided we want no flow BCs
# replace w more explicit treatment to check
            # # N
            # eta_diff = (-pad_eta[self._centpad] +
            #             pad_eta[self._northpad]) * self._dxbyy
            # self._ann[:] = eta_diff.clip(0.)
            # self._anp[:] = (-eta_diff).clip(0.)
            # # S
            # eta_diff = (-pad_eta[self._centpad] +
            #             pad_eta[self._southpad]) * self._dxbyy
            # self._ass[:] = eta_diff.clip(0.)
            # self._asp[:] = (-eta_diff).clip(0.)
            # # E
            # eta_diff = (-pad_eta[self._centpad] +
            #             pad_eta[self._eastpad]) * self._dybyx
            # self._aee[:] = eta_diff.clip(0.)
            # self._aep[:] = (-eta_diff).clip(0.)
            # # W
            # eta_diff = (-pad_eta[self._centpad] +
            #             pad_eta[self._westpad]) * self._dybyx
            # self._aww[:] = eta_diff.clip(0.)
            # self._awp[:] = (-eta_diff).clip(0.)
            # 
            # self._app[:] = self._awp + self._aep + self._asp + self._anp
            # 
            # # zero elevation treatment
            # # at a zero elevation we use a simple averaging approach
            # # this rationale is questionable - a propagation across flats may
            # # be preferable
            # self._app += self._min_slope_thresh
            # # this is VV's treatment for flats; now replaced with a simple-
            # # minded addition of small vals to app and axx ->
            # # flats = np.abs(self._app) < self._flat_thresh
            # # apz[flats] = 4
            # # for NSEW in (awz, aez, asz, anz):
            # #     NSEW[flats] = 1
            # # NOTE when we do not have a zero elevation condition the
            # # coefficients a*z are the upwind coefficents
            # 
            # if self._simple_flow:  # c == 1
            #     apz = self._app
            #     awz = self._aww
            #     aez = self._aee
            #     asz = self._ass
            #     anz = self._ann
            #     awpz = self._awp
            #     aepz = self._aep
            #     aspz = self._asp
            #     anpz = self._anp
            # else:  # Manning/Chezy style routing
            #     # second correction here accounts for fact it should be
            #     # slope, not d_eta, inside the term
            #     awz = self._aww**self._c * self.grid.dy**(1.-self._c)
            #     aez = self._aee**self._c * self.grid.dy**(1.-self._c)
            #     asz = self._ass**self._c * self.grid.dx**(1.-self._c)
            #     anz = self._ann**self._c * self.grid.dx**(1.-self._c)
            #     awpz = self._awp**self._c * self.grid.dy**(1.-self._c)
            #     aepz = self._aep**self._c * self.grid.dy**(1.-self._c)
            #     aspz = self._asp**self._c * self.grid.dx**(1.-self._c)
            #     anpz = self._anp**self._c * self.grid.dx**(1.-self._c)
            #     apz = awpz + aepz + aspz + anpz + self._min_slope_thresh

# replace here:
            # these are true gradients
            # N
            eta_diff = (-pad_eta[self._centpad] +
                        pad_eta[self._northpad]) / grid.dy
            self._ann[:] = eta_diff.clip(0.)
            self._anp[:] = (-eta_diff).clip(0.)
            # S
            eta_diff = (-pad_eta[self._centpad] +
                        pad_eta[self._southpad]) / grid.dy
            self._ass[:] = eta_diff.clip(0.)
            self._asp[:] = (-eta_diff).clip(0.)
            # E
            eta_diff = (-pad_eta[self._centpad] +
                        pad_eta[self._eastpad]) / grid.dx
            self._aee[:] = eta_diff.clip(0.)
            self._aep[:] = (-eta_diff).clip(0.)
            # W
            eta_diff = (-pad_eta[self._centpad] +
                        pad_eta[self._westpad]) / grid.dx
            self._aww[:] = eta_diff.clip(0.)
            self._awp[:] = (-eta_diff).clip(0.)

            ##self._app[:] = self._awp + self._aep + self._asp + self._anp

            # zero elevation treatment
            # at a zero elevation we use a simple averaging approach
            # this rationale is questionable - a propagation across flats may
            # be preferable
            ##self._app += self._min_slope_thresh
            # this is VV's treatment for flats; now replaced with a simple-
            # minded addition of small vals to app and axx ->
            # flats = np.abs(self._app) < self._flat_thresh
            # apz[flats] = 4
            # for NSEW in (awz, aez, asz, anz):
            #     NSEW[flats] = 1
            # NOTE when we do not have a zero elevation condition the
            # coefficients a*z are the upwind coefficents

            if self._simple_flow:  # c == 1
                # these are now discharge prefactors
                # self._app[:] = (self._awp*grid.dy + self._aep*grid.dy +
                #                 self._asp*grid.dx + self._anp*grid.dx)
                # self._app += self._min_slope_thresh
                # _app is no longer well defined
                awz = self._aww * grid.dy
                aez = self._aee * grid.dy
                asz = self._ass * grid.dx
                anz = self._ann * grid.dx
                awpz = self._awp * grid.dy
                aepz = self._aep * grid.dy
                aspz = self._asp * grid.dx
                anpz = self._anp * grid.dx
                apz = awpz + aepz + aspz + anpz + self._min_slope_thresh
            else:  # Manning/Chezy style routing
                # second correction here accounts for fact it should be
                # slope, not d_eta, inside the term
                awz = self._aww**self._c * self.grid.dy
                aez = self._aee**self._c * self.grid.dy
                asz = self._ass**self._c * self.grid.dx
                anz = self._ann**self._c * self.grid.dx
                awpz = self._awp**self._c * self.grid.dy
                aepz = self._aep**self._c * self.grid.dy
                aspz = self._asp**self._c * self.grid.dx
                anpz = self._anp**self._c * self.grid.dx
                apz = awpz + aepz + aspz + anpz + self._min_slope_thresh

            # Solve upwind equations for nodal K
            # this involves iteration to a stable solution
            # calc the new K based on incoming discharges
            mismatch = 1.
            # if VV's 4/1/1/1/1 flat scheme is used, this thresh is too low
            # current approach with tiny addition to app and axx works better
            while mismatch > self._stab_thresh:
                # in here, without VV's method, we can use axx and app instead
                # of axz and apz
                Knew.fill(self._min_slope_thresh)
                Knew[self._east] += awz[self._east] * K[self._west]
                Knew[self._westedge] += awz[self._westedge] * K[
                    self._westedge]
                Knew[self._west] += aez[self._west] * K[self._east]
                Knew[self._eastedge] += aez[self._eastedge] * K[
                    self._eastedge]
                Knew[self._north] += asz[self._north] * K[self._south]
                Knew[self._southedge] += asz[self._southedge] * K[
                    self._southedge]
                Knew[self._south] += anz[self._south] * K[self._north]
                Knew[self._northedge] += anz[self._northedge] * K[
                    self._northedge]
                Knew += Qsp
                Knew /= apz
                mismatch = np.sum(np.square(Knew - K))
                K[:] = Knew

            Kpad = np.pad(K, ((1, 1), (1, 1)), 'edge')
            self._Qw[:] = awz * Kpad[self._westpad]
            self._Qw -= awpz * K
            self._Qe[:] = aez * Kpad[self._eastpad]
            self._Qe -= aepz * K
            self._Qs[:] = asz * Kpad[self._southpad]
            self._Qs -= aspz * K
            self._Qn[:] = anz * Kpad[self._northpad]
            self._Qn -= anpz * K

            self._K = K  # ...to make it accessible

            # STABILITY ANALYSIS:
            # not immediately clear if messing with powers in equations will
            # alter this... This version derived from linear equations
            # assume rectangular grid (can't be warped)
            # establish the max flux across a face:
            EWmax = max(np.fabs(self._Qw).max(),
                        np.fabs(self._Qe).max())/self.grid.dy
            NSmax = max(np.fabs(self._Qn).max(),
                        np.fabs(self._Qs).max())/self.grid.dx
            maxwaterflux = max(EWmax, NSmax)
            if np.isclose(maxwaterflux, 0.):
                dt_stab = dt
                breakcheck = True
            else:
                dt_stab = self.grid.dx * self.grid.dy / (
                    2. * self._a * maxwaterflux)
                dt_stab = min(dt_stab, dt)
                if current_internal_time + dt_stab >= dt:
                    breakcheck = True
                    dt_stab = dt - current_internal_time

            # modify S if we're using MPM:
            if self._powerlaw is True:
                thresh = self._slope
            else:
                thresh = self._slope * K**(-2./3.)

            # do the sediment diffusion
            for dir in ('W', 'E', 'S', 'N'):
                self._grad_on_link(pad_eta, dir, self._slx, self._sly)
                Cslope = np.sqrt(self._slx**2 + self._sly**2)
                self._link_sed_flux_from_slope(Cslope, thresh, dir,
                                               a=self._a, m=self._m, n=self._n)

            try:
                Qin = Qsource.reshape((ni, nj))
            except AttributeError:
                Qin = float(Qsource)  # if both fail, we're in trouble
            eta[:] += dt_stab*(
                self._Qsed_e + self._Qsed_n + self._Qsed_w + self._Qsed_s +
                Qin)/grid.dx/grid.dy

            # increment the component time:
            current_internal_time += dt_stab

    @property
    def discharges_at_links(self):
        """Return the discharges at links.

        Note that if diagonal routing, this will return number_of_d8_links.
        Otherwise, it will be number_of_links.
        """
        return self._discharges_at_link

    def _grad_on_link(self, padded_eta, direction, slx, sly):
        """
        Updates slx and sly with link gradient values according to `direction`.

        eta = elevations in grid form
        direction = {'E', 'N', 'S', 'W'}

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 4), (0.5, 1.))
        >>> z = mg.add_zeros('node', 'topographic__elevation')
        >>> z[:] = np.array([[1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6]])
        >>> zpad = np.pad(z.reshape((3, 4)), ((1, 1), (1, 1)), 'edge')
        >>> dd = DischargeDiffuser(mg, 0.25)
        >>> dd._grad_on_link(zpad, 'W', dd._slx, dd._sly)
        """
        core = (slice(1, -1, 1), slice(1, -1, 1))
        if direction == 'W':
            slx[:] = (
                padded_eta[self._westpad] -
                padded_eta[self._centpad])/self.grid.dx
            sly[:] = padded_eta[self._SWpad]
            sly -= padded_eta[self._NWpad]
            sly += padded_eta[self._southpad]
            sly -= padded_eta[self._northpad]
            sly *= 0.25
            sly /= self.grid.dy

        elif direction == 'E':
            slx[:] = (
                padded_eta[self._eastpad] -
                padded_eta[self._centpad])/self.grid.dx
            sly[:] = padded_eta[self._SEpad]
            sly -= padded_eta[self._NEpad]
            sly += padded_eta[self._southpad]
            sly -= padded_eta[self._northpad]
            sly *= 0.25
            sly /= self.grid.dy

        elif direction == 'S':
            sly[:] = (
                padded_eta[self._southpad] -
                padded_eta[self._centpad])/self.grid.dy
            slx[:] = padded_eta[self._SWpad]
            slx -= padded_eta[self._SEpad]
            slx += padded_eta[self._westpad]
            slx -= padded_eta[self._eastpad]
            slx *= 0.25
            slx /= self.grid.dx

        elif direction == 'N':
            sly[:] = (
                padded_eta[self._northpad] -
                padded_eta[self._centpad])/self.grid.dy
            slx[:] = padded_eta[self._NWpad]
            slx -= padded_eta[self._NEpad]
            slx += padded_eta[self._westpad]
            slx -= padded_eta[self._eastpad]
            slx *= 0.25
            slx /= self.grid.dx

        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")

    def _link_sed_flux_from_slope(self, S_val, S_thresh,
                                  direction, a=1, m=1., n=1.):
        """
        Update the sed flux array for a given link dir, assuming a critical S.
        """
        if direction == 'W':
            dir_sed_flux = self._Qsed_w
            dir_water_flux = self._Qw
            thisslice = (self._east)
            deadedge = (self._westedge)
        elif direction == 'E':
            dir_sed_flux = self._Qsed_e
            dir_water_flux = self._Qe
            thisslice = (self._west)
            deadedge = (self._eastedge)
        elif direction == 'N':
            dir_sed_flux = self._Qsed_n
            dir_water_flux = self._Qn
            thisslice = (self._south)
            deadedge = (self._northedge)
        elif direction == 'S':
            dir_sed_flux = self._Qsed_s
            dir_water_flux = self._Qs
            thisslice = (self._north)
            deadedge = (self._southedge)
        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")
        if (np.isclose(m, 1.) and np.isclose(n, 1.)):
            slope_diff = (S_val - S_thresh).clip(0.)
        else:
            # VV's analysis in unpub MS shows e.g. MPM solution is close to
            # the linear response in most cases, so here we just crudely
            # modify the slopes with the requisite powers, without propagating
            # the nonlinearity throughout the whole method
            slope_diff = ((S_val**n - S_thresh**n).clip(0.))**m
        if np.isclose(a, 1.):
            dir_sed_flux[thisslice] = (dir_water_flux[thisslice] *
                                       slope_diff[thisslice])
        else:
            dir_sed_flux[thisslice] = a * (dir_water_flux[thisslice] *
                                           slope_diff[thisslice])
        dir_sed_flux[deadedge] = 0.


if __name__ == '__main__':
    import numpy as np
    from landlab import RasterModelGrid, imshow_grid_at_node
    S_crit = 0.25
    mg = RasterModelGrid((40, 40), (0.25, 0.25))
    z = mg.add_zeros('node', 'topographic__elevation')
    # source_nodes = np.sqrt(mg.node_x**2 + mg.node_y**2) < 3.
    # source_count = np.sum(source_nodes.astype(float))
    source_nodes = 0
    source_count = 1.
    Qw_in = mg.add_zeros('node', 'water__discharge_in')
    Qs_in = mg.add_zeros('node', 'sediment__discharge_in')
    # z -= mg.node_x/10000.
    # z -= mg.node_y/10000.
    Qw_in[source_nodes] = 0.5*np.pi/source_count
    Qs_in[source_nodes] = (1. - S_crit)*0.5*np.pi/source_count
    # Qw_in[0] = 0.5*np.pi  # 820 puts it about in the middle
    # Qs_in[0] = (1. - S_crit)*0.5*np.pi
    dd = DischargeDiffuser(mg, S_crit, flow_law_slope_exponent=1.,
                           sediment_flux_equation='powerlaw')
    for i in range(501):  # 501
        dd.run_one_step(0.08)  # 0.08
    imshow_grid_at_node(mg, 'topographic__elevation')
