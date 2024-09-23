"""This is an implementation of Vaughan Voller's experimental boundary method
reduced complexity flow router. Credit: Voller, Hobley, Paola.

Created on Fri Feb 20 09:32:27 2015

@author: danhobley (SiccarPoint), after volle001@umn.edu
"""

import numpy as np

from landlab import Component
from landlab import RasterModelGrid


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

    The algorithm uses a rule that looks like:

        q_sed = q_water * (S - S_crit)

    where S_crit is a critical slope threshold. [MODIFY THIS]

    It is VITAL you initialize this component AFTER setting boundary
    conditions.

    The primary method of this class is :func:`run_one_step`.

    Notes
    -----
    This is a "research grade" component, and is subject to dramatic change
    with little warning. No guarantees are made regarding its accuracy or
    utility. It is not recommended for user use yet!

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    None Listed

    """

    _name = "DischargeDiffuser"

    _unit_agnostic = True

    _info = {
        "flow__potential": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": (
                "Value of the hypothetical field 'K', used to force water "
                "flux to flow downhill"
            ),
        },
        "sediment__discharge_in": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Sediment discharge into a node.",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__discharge_in": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Incoming water discharge at node.",
        },
    }

    _min_slope_thresh = 1.0e-24
    # if your flow isn't connecting up, this probably needs to be reduced

    def __init__(self, grid, slope=0.25, flat_thresh=1.0e-4):
        """
        Parameters
        ----------
        grid : ModelGrid
            A grid.
        """
        super().__init__(grid)

        if isinstance(grid, RasterModelGrid):
            assert grid.number_of_node_rows >= 3
            assert grid.number_of_node_columns >= 3
            self._raster = True
        else:
            self._raster = False

        assert self._raster is True  # ...for now

        self._slope = slope
        self._flat_thresh = flat_thresh

        # hacky fix because water__discharge is defined on both links and nodes
        self.initialize_output_fields()

        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns

        self._K = grid.zeros("node", dtype=float)
        self._Knew = grid.zeros("node", dtype=float)
        self._prevK = grid.zeros("node", dtype=float)
        self._znew = grid.zeros("node", dtype=float)
        # discharge across north, south, west, and east face of control volume
        self._Qn = np.zeros((ni, nj), dtype="float")
        self._Qs = np.zeros((ni, nj), dtype="float")
        self._Qw = np.zeros((ni, nj), dtype="float")
        self._Qe = np.zeros((ni, nj), dtype="float")

        # coefficenst used in solition of flow conductivity K
        self._app = np.zeros((ni, nj), dtype="float")
        self._apz = np.zeros((ni, nj), dtype="float")
        self._aww = np.zeros((ni, nj), dtype="float")
        self._awp = np.zeros((ni, nj), dtype="float")
        self._awz = np.zeros((ni, nj), dtype="float")
        self._aee = np.zeros((ni, nj), dtype="float")
        self._aep = np.zeros((ni, nj), dtype="float")
        self._aez = np.zeros((ni, nj), dtype="float")
        self._ass = np.zeros((ni, nj), dtype="float")
        self._asp = np.zeros((ni, nj), dtype="float")
        self._asz = np.zeros((ni, nj), dtype="float")
        self._ann = np.zeros((ni, nj), dtype="float")
        self._anp = np.zeros((ni, nj), dtype="float")
        self._anz = np.zeros((ni, nj), dtype="float")

        self._slx = np.empty((ni, nj), dtype=float)
        self._sly = np.empty((ni, nj), dtype=float)
        self._Qsed_w = np.empty((ni, nj), dtype=float)
        self._Qsed_e = np.empty((ni, nj), dtype=float)
        self._Qsed_n = np.empty((ni, nj), dtype=float)
        self._Qsed_s = np.empty((ni, nj), dtype=float)

    def run_one_step(self, dt):
        """Run forward a duration of time, dt.

        Parameters
        ----------
        dt: float
        """
        grid = self._grid
        ni = grid.number_of_node_rows
        nj = grid.number_of_node_columns
        z = grid.at_node["topographic__elevation"]
        Qsp = grid.at_node["water__discharge_in"].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns)
        )
        Qsource = grid.at_node["sediment__discharge_in"].reshape(
            (grid.number_of_node_rows, grid.number_of_node_columns)
        )

        # #####STABILITY ANALYSIS GOES HERE
        dt_stab = dt

        # elevation at current and new time
        # Note a horizonal surface is the initial condition
        eta = z.reshape((ni, nj))
        K = self._K.reshape((ni, nj))
        Knew = self._Knew.reshape((ni, nj))
        # etan = self._znew.reshape((grid.number_of_node_rows,
        #                            grid.number_of_node_columns))

        # pad eta
        pad_eta = np.pad(eta, ((1, 1), (1, 1)), "edge")
        # do the sediment diffusion
        for dir in ("W", "E", "S", "N"):
            self._grad_on_link(pad_eta, dir)
            Cslope = np.sqrt(self._slx**2 + self._sly**2)
            self._link_sed_flux_from_slope(Cslope, self._slope, dir)

        try:
            Qin = Qsource.reshape((ni, nj))
        except AttributeError:
            Qin = float(Qsource)  # if both fail, we're in trouble
        eta[:] += (
            dt_stab
            * (self._Qsed_e + self._Qsed_n + self._Qsed_w + self._Qsed_s + Qin)
            / grid.dx
            / grid.dy
        )

        # do the water routing on links
        # These calculations are based on assuming that the flow is a sheet
        # flow that can be characterized with a potential equation. If this
        # flow is isotropic (an assumption we should revisit) with this model
        # the flow discharge in the x-direction (for example) can be calculated
        # as a constant (K the 'flow conductivity') times the component of the
        # sediment slope in that direction. It helps to define a 'slope
        # velocity' u, with components ustar=-deta/dx and vstar=-deta/dx which
        # allows us to write down the following advection like gov. eq. for
        # the flow  ----div(Ku)+Q=0---where Q represents external flow inputs

        # Since we can readily determine u from the current sediment topography
        # We solve this equation for K using an upwind scheme

        # Build upwinded coefficients. Vals only 0 if if flow is in upwind dir
        # note cols/rows which don't get updated will always remain as 0,
        # which is right provided we want no flow BCs
        eta_diff = -eta[:-1, :] + eta[1:, :]
        self._ann[:-1, :] = eta_diff.clip(0.0)
        self._anp[:-1, :] = (-eta_diff).clip(0.0)
        eta_diff = -eta[1:, :] + eta[:-1, :]
        self._ass[1:, :] = eta_diff.clip(0.0)
        self._asp[1:, :] = (-eta_diff).clip(0.0)
        eta_diff = -eta[:, :-1] + eta[:, 1:]
        self._aee[:, :-1] = eta_diff.clip(0.0)
        self._aep[:, :-1] = (-eta_diff).clip(0.0)
        eta_diff = -eta[:, 1:] + eta[:, :-1]
        self._aww[:, 1:] = eta_diff.clip(0.0)
        self._awp[:, 1:] = (-eta_diff).clip(0.0)

        self._app[:] = self._awp + self._aep + self._asp + self._anp

        apz = self._app.copy()
        awz = self._aww.copy()
        aez = self._aee.copy()
        asz = self._ass.copy()
        anz = self._ann.copy()
        # zero elevation treatment
        # at a zero elevation we use a simple averaging approach
        # this rationale is questionable - a propagation across flats may be
        # preferable
        flats = np.abs(self._app) < self._flat_thresh
        apz[flats] = 4
        for NSEW in (awz, aez, asz, anz):
            NSEW[flats] = 1
        # NOTE when we do not have a zero elevation condition the
        # coefficients a*z are the upwind coefficents

        # Solve upwind equations for nodal K
        # this involves iteration to a stable solution
        # #####IMPLEMENT IT
        # calc the new K based on incoming discharges
        for _ in range(1):
            Knew[:, 1:] += awz[:, 1:] + K[:, :-1]
            Knew[:, 0] += awz[:, 0] + K[:, 0]
            Knew[:, :-1] += aez[:, :-1] + K[:, 1:]
            Knew[:, -1] += awz[:, -1] + K[:, -1]
            Knew[1:, :] += asz[1:, :] + K[:-1, :]
            Knew[0, :] += asz[0, :] + K[0, :]
            Knew[:-1, :] += anz[:-1, :] + K[1:, :]
            Knew[-1, :] += asz[-1, :] + K[-1, :]
            Knew += Qsp
            Knew /= apz
            K[:] = Knew

        Kpad = np.pad(K, ((1, 1), (1, 1)), "edge")
        self._Qw += self._aww * Kpad[1:-1, :-2]
        self._Qw -= self._awp * K
        self._Qe += self._aee * Kpad[1:-1, 2:]
        self._Qe -= self._aep * K
        self._Qs += self._ass * Kpad[:-2, 1:-1]
        self._Qs -= self._asp * K
        self._Qn += self._ann * Kpad[2:, 1:-1]
        self._Qn -= self._anp * K

    @property
    def discharges_at_links(self):
        """Return the discharges at links.

        Note that if diagonal routing, this will return
        number_of_d8_links. Otherwise, it will be number_of_links.
        """
        return self._discharges_at_link

    def _grad_on_link(self, padded_eta, direction):
        """Updates slx and sly with link gradient values according to
        `direction`.

        eta = elevations in grid form
        direction = {'E', 'N', 'S', 'W'}

        Examples
        --------
        >>> import numpy as np
        >>> from landlab import RasterModelGrid
        >>> mg = RasterModelGrid((3, 4), xy_spacing=(1.0, 0.5))
        >>> z = mg.add_zeros("topographic__elevation", at="node")
        >>> z = mg.add_zeros("water__discharge_in", at="node")
        >>> z = mg.add_zeros("sediment__discharge_in", at="node")
        >>> z[:] = np.array([[1, 2, 3, 4, 2, 3, 4, 5, 3, 4, 5, 6]])
        >>> zpad = np.pad(z.reshape((3, 4)), ((1, 1), (1, 1)), "edge")
        >>> dd = DischargeDiffuser(mg, 0.25)
        >>> dd._grad_on_link(zpad, "W")
        """
        core = (slice(1, -1, 1), slice(1, -1, 1))
        if direction == "W":
            self._slx[:] = (padded_eta[1:-1, :-2] - padded_eta[core]) / self._grid.dx
            self._sly[:] = padded_eta[:-2, :-2]
            self._sly -= padded_eta[2:, :-2]
            self._sly += padded_eta[:-2, 1:-1]
            self._sly -= padded_eta[2:, 1:-1]
            self._sly *= 0.25
            self._sly /= self._grid.dy

        elif direction == "E":
            self._slx[:] = (padded_eta[1:-1, 2:] - padded_eta[core]) / self._grid.dx
            self._sly[:] = padded_eta[:-2, 2:]
            self._sly -= padded_eta[2:, 2:]
            self._sly += padded_eta[:-2, 1:-1]
            self._sly -= padded_eta[2:, 1:-1]
            self._sly *= 0.25
            self._sly /= self._grid.dy

        elif direction == "S":
            self._sly[:] = (padded_eta[:-2, 1:-1] - padded_eta[core]) / self._grid.dy
            self._slx[:] = padded_eta[:-2, :-2]
            self._slx -= padded_eta[:-2, 2:]
            self._slx += padded_eta[1:-1, :-2]
            self._slx -= padded_eta[1:-1, 2:]
            self._slx *= 0.25
            self._slx /= self._grid.dx

        elif direction == "N":
            self._sly[:] = (padded_eta[2:, 1:-1] - padded_eta[core]) / self._grid.dy
            self._slx[:] = padded_eta[2:, :-2]
            self._slx -= padded_eta[2:, 2:]
            self._slx += padded_eta[1:-1, :-2]
            self._slx -= padded_eta[1:-1, 2:]
            self._slx *= 0.25
            self._slx /= self._grid.dx

        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")

    def _link_sed_flux_from_slope(self, S_val, S_thresh, direction):
        """Update the sed flux array for a given link dir, assuming a critical
        S."""
        if direction == "W":
            dir_sed_flux = self._Qsed_w
            dir_water_flux = self._Qw
            thisslice = (slice(0, -1, 1), slice(1, -1, 1))
            deadedge = (slice(0, -1, 1), slice(0, 1, 1))
        elif direction == "E":
            dir_sed_flux = self._Qsed_e
            dir_water_flux = self._Qe
            thisslice = (slice(0, -1, 1), slice(1, -2, 1))
            deadedge = (slice(0, -1, 1), slice(-2, -1, 1))
        elif direction == "N":
            dir_sed_flux = self._Qsed_n
            dir_water_flux = self._Qn
            thisslice = (slice(0, -2, 1), slice(0, -1, 1))
            deadedge = (slice(-2, -1, 1), slice(0, -1, 1))
        elif direction == "S":
            dir_sed_flux = self._Qsed_s
            dir_water_flux = self._Qs
            thisslice = (slice(1, -1, 1), slice(0, -1, 1))
            deadedge = (slice(0, 1, 1), slice(0, -1, 1))
        else:
            raise NameError("direction must be {'E', 'N', 'S', 'W'}")
        slope_diff = (S_val - S_thresh).clip(0.0)
        dir_sed_flux[thisslice] = dir_water_flux[thisslice] * slope_diff[thisslice]
        dir_sed_flux[deadedge] = 0.0

    def diffuse_sediment(self, Qw_in, Qsed_in):
        """ """
        pass


if __name__ == "__main__":
    from landlab import imshow_grid_at_node

    S_crit = 0.25
    mg = RasterModelGrid((20, 20), 0.5)
    mg.add_zeros("topographic__elevation", at="node")
    Qw_in = mg.add_zeros("water__discharge_in", at="node")
    Qs_in = mg.add_zeros("sediment__discharge_in", at="node")
    Qw_in[0] = 0.5 * np.pi
    Qs_in[0] = (1.0 - S_crit) * 0.5 * np.pi
    dd = DischargeDiffuser(mg, S_crit)
    for _ in range(5):  # 501
        dd.run_one_step(0.01)  # 0.08
    imshow_grid_at_node(mg, "topographic__elevation")
