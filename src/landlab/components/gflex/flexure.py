"""Landlab wrapper for the gFlex lithospheric flexure model.

Original wrapper by Daniel Hobley (SiccarPoint), 2015.
Redesigned 2026: A. Wickert with Claude Code (Anthropic).
"""

import numpy as np
import scipy.constants

from landlab import Component
from landlab import FieldError
from landlab import RasterModelGrid

try:
    import gflex
except ImportError:
    NO_GFLEX = True
else:
    NO_GFLEX = False


class gFlex(Component):
    """Landlab wrapper for the gFlex lithospheric flexure model.

    Computes the steady-state flexural isostatic deflection of the
    lithosphere given a surface load stress field.  The deflection is
    written to ``lithosphere__vertical_displacement`` on the grid.
    Applying that displacement to topography (or any other field) is the
    caller's responsibility.

    Install gFlex with::

        pip install gflex

    For full documentation see https://github.com/awickert/gFlex.

    Note that gFlex component **demands lengths in meters**, including
    the grid spacing.

    Examples
    --------

    NB: these tests are not actually run as our automated testing becomes
    confused if gFlex is not installed on the testing machine!

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import gFlex
    >>> mg = RasterModelGrid((10, 10), xy_spacing=25000.0)
    >>> stress = mg.add_zeros("surface_load__stress", at="node", dtype=float)
    >>> stress.view().reshape(mg.shape)[3:7, 3:7] += 1.0e6
    >>> gf = gFlex(mg, bc_east="free", bc_north="periodic", bc_south="periodic")  # doctest: +SKIP
    >>> gf.run_one_step()  # doctest: +SKIP
    >>> mg.at_node["lithosphere__vertical_displacement"]  # doctest: +SKIP

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Wickert, A. (2016). Open-source modular solutions for flexural isostasy:
    gFlex v1.0. Geoscientific Model Development  9(3), 997-1017.
    https://dx.doi.org/10.5194/gmd-9-997-2016

    **Additional References**

    None Listed

    """

    _name = "gFlex"

    _unit_agnostic = True

    _cite_as = """
    @article{wickert2016open,
      author = {Wickert, A. D.},
      title = {{Open-source modular solutions for flexural isostasy: gFlex v1.0}},
      issn = {1991-959X},
      doi = {10.5194/gmd-9-997-2016},
      pages = {997--1017},
      number = {3},
      volume = {9},
      journal = {Geoscientific Model Development},
      year = {2016}
    }
    """

    _info = {
        "lithosphere__elastic_thickness": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Elastic thickness of the lithosphere. When present, read "
                "on every run_one_step() call; update via BMI set_value() "
                "to vary T_e between coupling steps."
            ),
        },
        "lithosphere__vertical_displacement": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Vertical displacement of the lithosphere surface due to "
                "flexural isostasy. Downward displacement is negative."
            ),
        },
        "surface_load__stress": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "Pa",
            "mapping": "node",
            "doc": "Magnitude of stress exerted by surface load",
        },
    }

    def __init__(
        self,
        grid,
        Youngs_modulus=65e9,
        Poissons_ratio=0.25,
        rho_mantle=3300.0,
        rho_fill=0.0,
        elastic_thickness=35000.0,
        method="fd",
        quiet=True,
        bc_west="no_outside_loads",
        bc_east="no_outside_loads",
        bc_north="no_outside_loads",
        bc_south="no_outside_loads",
        g=scipy.constants.g,
    ):
        """Constructor for the gFlex Landlab component.

        Parameters
        ----------
        Youngs_modulus : float
            Young's modulus for the lithosphere [Pa].
        Poissons_ratio : float
            Poisson's ratio for the lithosphere.
        rho_mantle : float
            Density of the mantle [kg m⁻³].
        rho_fill : float
            Density of the infilling material (air, water, sediment…)
            [kg m⁻³].
        elastic_thickness : float, str, or array-like
            Elastic thickness of the lithosphere [m]. May be a scalar
            float, the name of an existing node field on the grid, or a
            NumPy array of shape ``grid.shape``. To vary T_e between
            timesteps, create a ``lithosphere__elastic_thickness`` node
            field and update it before each ``run_one_step()`` call.
        method : str
            Solution method: ``'fd'`` (finite difference, supports
            variable T_e and all BC types), ``'fft'``, ``'sas'``, or
            ``'sas_ng'``.
        quiet : bool
            Suppress gFlex log output.
        bc_west, bc_east, bc_north, bc_south : str or dict
            Boundary condition for each grid edge. Any string accepted by
            ``gflex.VALID_BC_STRINGS_2D`` is valid. Canonical names:
            ``'no_outside_loads'`` (alias ``'infinite'``, default),
            ``'zero_displacement_zero_slope'`` (alias ``'clamped'``),
            ``'zero_displacement_zero_moment'`` (alias ``'pinned'``),
            ``'zero_moment_zero_shear'`` (alias ``'free'``),
            ``'zero_slope_zero_shear'`` (alias ``'mirror'``),
            ``'periodic'``. Periodic boundaries must be paired
            (west–east or north–south). A dict of prescribed values
            (e.g. ``{"displacement": arr, "slope": arr}``) may be
            passed for inhomogeneous (nested-domain) BCs.
        g : float
            Acceleration due to gravity [m s⁻²].
        """
        super().__init__(grid)

        if not isinstance(grid, RasterModelGrid):
            raise TypeError(f"gFlex requires a RasterModelGrid; got {type(grid)}")

        if NO_GFLEX:
            raise ImportError(
                "gFlex not installed! For installation instructions see "
                "https://github.com/awickert/gFlex"
            )

        self._flex = gflex.F2D()
        flex = self._flex

        flex.dx = grid.dx
        flex.dy = grid.dy
        flex.method = method.lower()
        flex.quiet = quiet
        flex.E = float(Youngs_modulus)
        flex.nu = float(Poissons_ratio)
        flex.rho_m = float(rho_mantle)
        flex.rho_fill = float(rho_fill)
        flex.g = float(g)

        for name, val in (
            ("bc_west", bc_west),
            ("bc_east", bc_east),
            ("bc_north", bc_north),
            ("bc_south", bc_south),
        ):
            if not isinstance(val, dict) and val not in gflex.VALID_BC_STRINGS_2D:
                raise ValueError(
                    f"{name}={val!r} is not a valid boundary condition. "
                    f"Choose from: {sorted(gflex.VALID_BC_STRINGS_2D)}"
                )

        if (bc_west == "periodic") != (bc_east == "periodic"):
            raise ValueError("bc_west and bc_east must both be 'periodic', or neither.")
        if (bc_north == "periodic") != (bc_south == "periodic"):
            raise ValueError(
                "bc_north and bc_south must both be 'periodic', or neither."
            )

        flex.bc_west = bc_west
        flex.bc_east = bc_east
        flex.bc_north = bc_north
        flex.bc_south = bc_south

        Te_in = elastic_thickness
        try:
            flex.T_e = float(Te_in)
        except (TypeError, ValueError):
            try:
                flex.T_e = grid.at_node[Te_in].view().reshape(grid.shape)
            except TypeError:
                flex.T_e = np.asarray(Te_in, dtype=float).reshape(grid.shape)

        flex.qs = grid.at_node["surface_load__stress"].view().reshape(grid.shape)
        flex.initialize()

        if isinstance(flex.T_e, float):
            flex.cache_factorization = True

        self._grid.add_zeros(
            "lithosphere__vertical_displacement",
            at="node",
            dtype=float,
            clobber=True,
        )

    def run_one_step(self):
        """Compute the flexural isostatic deflection for the current load.

        Reads ``surface_load__stress`` (and ``lithosphere__elastic_thickness``
        if present) from the grid, runs the gFlex solver, and writes the
        result to ``lithosphere__vertical_displacement``.
        """
        try:
            self._flex.T_e = (
                self._grid.at_node["lithosphere__elastic_thickness"]
                .view()
                .reshape(self._grid.shape)
            )
        except FieldError:
            pass

        self._flex.qs = (
            self._grid.at_node["surface_load__stress"]
            .view()
            .reshape(self._grid.shape)
        )
        self._flex.run()
        self._grid.at_node["lithosphere__vertical_displacement"][
            :
        ] = self._flex.w.view().ravel()
