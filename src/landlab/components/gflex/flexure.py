"""This is a Landlab wrapper for A Wickert's gFlex flexure model (Wickert et
al., submitted to Geoscientific Model Development). The most up-to-date version
of his code can be found at github.com/awickert/gFlex.

This Landlab wrapper will use a snapshot of that code, which YOU need to
install on your own machine.
A stable snapshot of gFlex is hosted on PyPI, which is the recommended version
to install.
If you have pip (the Python package install tool), simply run
'pip install gFlex' from a command prompt.
Alternatively, you can download and unpack the code (from github, or with PyPI,
pypi.python.org/pypi/gFlex/), then run 'python setup.py install'.

Created on Thu Feb 19 18:47:11 2015

@author: daniel.hobley (SiccarPoint @Github)

...following AW's run_in_script_2D.py.
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
    """This is a Landlab wrapper for A Wickert's gFlex flexure model (Wickert
    et al., 2016, Geoscientific Model Development). The most up-to-date version
    of his code can be found at github.com/awickert/gFlex.

    This Landlab wrapper will use a snapshot of that code, which YOU need to
    install on your own machine.
    A stable snapshot of gFlex is hosted on PyPI, which is the recommended
    version to install.
    If you have pip (the Python package install tool), simply run
    'pip install gFlex' from a command prompt.
    Alternatively, you can download and unpack the code (from github, or with
    PyPI, pypi.python.org/pypi/gFlex/), then run 'python setup.py install'.

    Note that gFlex maintains its own internal version if the grid, but this
    should not affect performance.

    This component will modify the topographic__elevation field only if one
    already exists. Note that the gFlex component **demands lengths in
    meters**, including the grid dimensions.
    The component also recognises the gFlex specific parameters 'Method',
    'PlateSolutionType', 'Solver', and 'Quiet'. See the gFlex software
    documentation for more details.

    Examples
    --------

    NB: these tests are not actually run as our automated testing becomes
    confused if gFlex is not installed on the testing machine!

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import gFlex
    >>> mg = RasterModelGrid((10, 10), xy_spacing=25000.0)
    >>> z = mg.add_zeros("topographic__elevation", at="node", dtype=float)
    >>> stress = mg.add_zeros("surface_load__stress", at="node", dtype=float)
    >>> stress.view().reshape(mg.shape)[3:7, 3:7] += 1.0e6
    >>> gf = gFlex(
    ...     mg, BC_E="0Moment0Shear", BC_N="Periodic", BC_S="Periodic"
    ... )  # doctest: +SKIP
    >>> gf.run_one_step()  # doctest: +SKIP

    N-S profile across flexed plate:

    >>> z.reshape(mg.shape)[:, 5]  # doctest: +SKIP
    array([-4.54872677, -4.6484927 , -4.82638669, -5.03001546, -5.15351385,
           -5.15351385, -5.03001546, -4.82638669, -4.6484927 , -4.54872677])

    W-E profile, noting the free BC to the east side:

    >>> z.reshape(mg.shape)[5, :]  # doctest: +SKIP
    array([-0.43536739, -1.19197738, -2.164915  , -3.2388464 , -4.2607558 ,
           -5.15351385, -5.89373366, -6.50676947, -7.07880156, -7.63302576])

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
        "lithosphere_surface__elevation_increment": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "The change in elevation of the top of the lithosphere (the "
                "land surface) in one timestep"
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
        "topographic__elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        Youngs_modulus=6.5e11,
        Poissons_ratio=0.25,
        rho_mantle=3300.0,
        rho_fill=0.0,
        elastic_thickness=35000.0,
        Method="FD",
        Solver="direct",
        PlateSolutionType="vWC1994",
        quiet=True,
        BC_W="0Displacement0Slope",
        BC_E="0Displacement0Slope",
        BC_N="0Displacement0Slope",
        BC_S="0Displacement0Slope",
        g=scipy.constants.g,
    ):
        """Constructor for Wickert's gFlex in Landlab.

        Parameters
        ----------
        Youngs_modulus : float
            Young's modulus for the lithosphere.
        Poissons_ratio : float
            Poisson's ratio for the lithosphere.
        rho_mantle : float (kg*m**-3)
            The density of the mantle.
        rho_fill : float (kg*m**-3)
            The density of the infilling material (air, water...)
        elastic_thickness : float (m)
            The elastic thickness of the lithosphere.
        BC_W, BC_E, BC_N, BC_S : {'0Displacement0Slope', '0Moment0Shear',
                                  'Periodic'}
            The boundary condition status of each grid edge, following gFlex's
            definitions. Periodic boundaries must be paired (obviously).
        g : float (m*s**-2)
            The acceleration due to gravity.
        """
        super().__init__(grid)

        assert isinstance(grid, RasterModelGrid)

        if NO_GFLEX:
            raise ImportError(
                "gFlex not installed! For installation instructions see "
                + "gFlex on GitHub: https://github.com/awickert/gFlex"
            )
        BC_options = (
            "0Displacement0Slope",
            "0Moment0Shear",
            "0Slope0Shear",
            "Periodic",
        )

        # instantiate the module:
        self._flex = gflex.F2D()
        flex = self._flex

        # set up the grid variables:

        flex.dx = grid.dx
        flex.dy = grid.dy

        # we assume these properties are fixed in this relatively
        # straightforward implementation, but they can still be set if you
        # want:
        flex.Method = Method
        flex.PlateSolutionType = PlateSolutionType
        flex.Solver = Solver
        flex.Quiet = quiet

        flex.E = float(Youngs_modulus)
        flex.nu = float(Poissons_ratio)
        flex.rho_m = float(rho_mantle)
        flex.rho_fill = float(rho_fill)
        flex.g = float(g)
        flex.BC_W = BC_W
        flex.BC_E = BC_E
        flex.BC_S = BC_S
        flex.BC_N = BC_N
        for i in (flex.BC_E, flex.BC_W, flex.BC_N, flex.BC_S):
            assert i in BC_options

        if BC_W == "Periodic":
            assert BC_E == "Periodic"
        if BC_E == "Periodic":
            assert BC_W == "Periodic"
        if BC_N == "Periodic":
            assert BC_S == "Periodic"
        if BC_S == "Periodic":
            assert BC_N == "Periodic"

        Te_in = elastic_thickness
        try:
            flex.Te = float(Te_in)
        except ValueError:
            try:
                flex.Te = grid.at_node[Te_in].view().reshape(grid.shape)
            except TypeError:
                flex.Te = Te_in.view().reshape(grid.shape)
            self._input_var_names.add(Te_in)
            self._output_var_names.add(Te_in)

        # set up the link between surface load stresses in the gFlex component
        # and the LL grid field:
        flex.qs = grid.at_node["surface_load__stress"].view().reshape(grid.shape)

        # create a holder for the "pre-flexure" state of the grid, to allow
        # updating of elevs:
        self._pre_flex = np.zeros(grid.number_of_nodes, dtype=float)

        # create the primary output field:
        self._grid.add_zeros(
            "lithosphere_surface__elevation_increment",
            at="node",
            dtype=float,
            clobber=True,
        )

    def flex_lithosphere(self):
        """Executes (& finalizes, from the perspective of gFlex) the core
        method of gFlex.

        Note that flexure of the lithosphere proceeds to steady state in
        a single timestep.
        """
        self._flex.qs = (
            self._grid.at_node["surface_load__stress"].view().reshape(self._grid.shape)
        )
        self._flex.initialize()
        self._flex.run()
        self._flex.finalize()

        self._grid.at_node["lithosphere_surface__elevation_increment"][
            :
        ] = self._flex.w.view().ravel()

        try:
            self._grid.at_node["topographic__elevation"]
            # a topo exists...
        except FieldError:
            pass
        else:
            topo_diff = (
                self._grid.at_node["lithosphere_surface__elevation_increment"]
                - self._pre_flex
            )
            self._grid.at_node["topographic__elevation"] += topo_diff
            self._pre_flex += topo_diff

    def run_one_step(self):
        """Flex the lithosphere to find its steady state form.

        The standardized run method for this component.

        Parameters
        ----------
        None
        """
        self._flex_lithosphere()
