import numpy as np

from landlab.components.erosion_deposition.erosion_deposition import ErosionDeposition
from landlab.components.erosion_deposition.generalized_erosion_deposition import (
    DEFAULT_MINIMUM_TIME_STEP,
)
from landlab.utils import return_array_at_node

TIME_STEP_FACTOR = 0.5


class SharedStreamPower(ErosionDeposition):
    """Shared Stream Power Model in the style of Hergarten (2021).

    Implements the Shared Stream Power Model in the style of Hergarten (2021).
    Designed to simultaneously model
    river incision and sediment transport in order to seamlessly
    transition between detachment limited to transport limited erosion.
    Mathematically equivalent to the linear decline model from Davy and Lague (2009),
    which is used by the base class, ErosionDeposition. In addition, this component is
    designed to work with varying runoff rates, and can update the discharge
    and other parameters effected by discharge with each timestep.

    Here is the equation for erosion without a threshold::

        E = k_bedrock * A**m_sp * S**n_sp - k_bedrock / k_transport * Qs / A

    where ``Q`` is water discharge, ``Qs`` is sediment flux, ``S`` is slope, ``m_sp``
    and ``n_sp`` are scaling exponents, coefficient ``k_bedrock`` is the erodibility of
    the bedrock and coefficient ``k_transport`` is the ability to transport sediment.

    The first term, ``k_bedrock * A**m_sp * S**n_sp``, is the incision term, and the
    second term, ``k_bedrock / k_transport * Qs / A``, is the transport term. Note that
    ``k_bedrock / k_transport`` determines the relative amount of incision and
    sediment transport. ``k_bedrock`` modifies the incision term.

    The equivalent equation used by ErosionDeposition from Davy & Lague (2009) is::

        E = K * q**m_sp * S**n_sp - v_s * Qs / q

    where ``K`` is sediment erodibility, ``v_s`` is the settling velocity for sediment,
    and ``q`` is the water discharge.

    To translate the shared stream power input for ErosionDeposition, we use the
    equations::

        q = Ar
        K = k_bedrock / r**m_sp
        v_s = k_bedrock / k_transport

    It is important to note that the second two equations were derived only
    for calibrating the model, and do not necessarily correlate to landscape evolution
    in nature.

    To write the final equation we define the incision term as omega::

        omega = k_bedrock * A**m_sp * S**n_sp

    and incorporate ``sp_crit``, the critical stream power needed to erode bedrock,
    giving::

        E = omega * (1 - exp(omega / sp_crit) ) - k_bedrock / k_transport * Qs / A


    Written by A. Thompson.

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Hergarten, S. (2021). The influence of sediment transport on stationary
    and mobile knickpoints in river profiles. Journal of Geophysical Research:
    Earth Surface, 126, e2021JF006218. https://doi.org/10.1029/2021JF006218

    **Additional References**

    Barnhart, K., Glade, R., Shobe, C., Tucker, G. (2019). Terrainbento 1.0: a
    Python package for multi-model analysis in long-term drainage basin
    evolution. Geoscientific Model Development  12(4), 1267--1297.
    https://dx.doi.org/10.5194/gmd-12-1267-2019

    Davy, P., Lague, D. (2009). Fluvial erosion/transport equation of landscape
    evolution models revisited Journal of Geophysical Research  114(F3),
    F03007. https://dx.doi.org/10.1029/2008jf001146

    Examples
    ---------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.components import DepressionFinderAndRouter
    >>> from landlab.components import ErosionDeposition
    >>> from landlab.components import FastscapeEroder
    >>> np.random.seed(seed=5000)

    Define grid and initial topography:

    * 5x5 grid with baselevel in the lower left corner
    * All other boundary nodes closed
    * Initial topography is plane tilted up to the upper right + noise

    >>> grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    >>> grid.at_node["topographic__elevation"] = (
    ...     grid.y_of_node / 10
    ...     + grid.x_of_node / 10
    ...     + np.random.rand(grid.number_of_nodes) / 10
    ... )
    >>> grid.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=True,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> grid.set_watershed_boundary_condition_outlet_id(
    ...     0, grid.at_node["topographic__elevation"], -9999.0
    ... )
    >>> fsc_dt = 100.0
    >>> ed_dt = 1.0

    Check initial topography

    >>> grid.at_node["topographic__elevation"].reshape(grid.shape)
    array([[0.02290479, 1.03606698, 2.0727653 , 3.01126678, 4.06077707],
           [1.08157495, 2.09812694, 3.00637448, 4.07999597, 5.00969486],
           [2.04008677, 3.06621577, 4.09655859, 5.04809001, 6.02641123],
           [3.05874171, 4.00585786, 5.0595697 , 6.04425233, 7.05334077],
           [4.05922478, 5.0409473 , 6.07035008, 7.0038935 , 8.01034357]])

    Instantiate Fastscape eroder, flow router, and depression finder

    >>> fr = FlowAccumulator(grid, flow_director="D8")
    >>> df = DepressionFinderAndRouter(grid)
    >>> fsc = FastscapeEroder(grid, K_sp=0.001, m_sp=0.5, n_sp=1)

    Burn in an initial drainage network using the Fastscape eroder:

    >>> for _ in range(100):
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     flooded = np.where(df.flood_status == 3)[0]
    ...     fsc.run_one_step(dt=fsc_dt)
    ...     grid.at_node["topographic__elevation"][0] -= 0.001  # uplift
    ...

    Instantiate the SharedStreamPower component:

    >>> ssp = SharedStreamPower(
    ...     grid,
    ...     k_bedrock=0.00001,
    ...     k_transport=0.001,
    ...     m_sp=0.5,
    ...     n_sp=1.0,
    ...     sp_crit=0,
    ... )

    Now run the E/D component for 2000 short timesteps:

    >>> for _ in range(2000):  # E/D component loop
    ...     fr.run_one_step()
    ...     df.map_depressions()
    ...     ssp.run_one_step(dt=ed_dt)
    ...     grid.at_node["topographic__elevation"][0] -= 2e-4 * ed_dt
    ...

    Now we test to see if topography is right:

    >>> np.around(grid.at_node["topographic__elevation"], decimals=3).reshape(
    ...     grid.shape
    ... )
    array([[-0.477,  1.036,  2.073,  3.011,  4.061],
           [ 1.082, -0.08 , -0.065, -0.054,  5.01 ],
           [ 2.04 , -0.065, -0.065, -0.053,  6.026],
           [ 3.059, -0.054, -0.053, -0.035,  7.053],
           [ 4.059,  5.041,  6.07 ,  7.004,  8.01 ]])
    """

    _name = "SharedStreamPower"

    _unit_agnostic = True

    def __init__(
        self,
        grid,
        k_bedrock=0.001,
        k_transport=0.001,
        runoff_rate=1.0,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit=0.0,
        F_f=0.0,
        discharge_field="surface_water__discharge",
        solver="basic",
    ):
        """Initialize the Shared Stream Power model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        k_bedrock : str, or array_like, optional
            Erodibility for bedrock (units vary).
        k_transport : str, or array_like, optional
            Ability to transport sediment (units vary).
        runoff_rate : float, optional
            Runoff rate. Scales Q = Ar. [m/yr]
        m_sp : float, optional
            Discharge exponent (units vary).
        n_sp : float, optional
            Slope exponent (units vary).
        sp_crit : str or array_like
            Critical stream power to erode substrate [E/(TL^2)]
        F_f : float, optional
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load.
        discharge_field : str or array_like, optional
            Discharge [L^2/T]. The default is to use the grid field
            ``"surface_water__discharge"``, which is simply drainage area
            multiplied by the default rainfall rate (1 m/yr). To use custom
            spatially/temporally varying rainfall, use 'water__unit_flux_in'
            to specify water input to the FlowAccumulator.
        solver : {"basic", "adaptive"}, optional
            Solver to use. Options at present include:

            1. ``"basic"`` (default): explicit forward-time extrapolation.
               Simple but will become unstable if time step is too large.
            2. ``"adaptive"``: adaptive time-step solver that estimates a
               stable step size based on the shortest time to "flattening"
               among all upstream-downstream node pairs.
        """
        self._discharge_field = discharge_field
        self._runoff_rate = runoff_rate
        self._k_bedrock = k_bedrock
        self._k_transport = k_transport

        # convert shared stream power inputs to erosion deposition inputs
        v_s = self.k_bedrock * self.runoff_rate / self.k_transport
        K_s = self.k_bedrock / self.runoff_rate**m_sp

        # instantiate ErosionDeposition
        super().__init__(
            grid,
            K=K_s,
            v_s=v_s,
            m_sp=m_sp,
            n_sp=n_sp,
            sp_crit=sp_crit,
            F_f=F_f,
            discharge_field=discharge_field,
            solver=solver,
            dt_min=DEFAULT_MINIMUM_TIME_STEP,
        )

    @property
    def k_bedrock(self):
        """Erodibility for bedrock (units vary)."""
        if isinstance(self._k_bedrock, str):
            return self.grid.at_node[self._k_bedrock]
        else:
            return self._k_bedrock

    @property
    def k_transport(self):
        """Ability to transport sediment (units vary)."""
        if isinstance(self._k_transport, str):
            return self.grid.at_node[self._k_transport]
        else:
            return self._k_transport

    @property
    def runoff_rate(self):
        """Runoff rate. Scales Q = Ar. [m/yr]"""
        if isinstance(self._runoff_rate, str):
            return self.grid.at_node[self._runoff_rate]
        else:
            return self._runoff_rate

    def update_runoff(self, new_runoff=1.0):
        """Update runoff variables.

        Updates ``runoff_rate``, ``K``, ``v_s``, and ``"water__unit_flux_in"`` for a new
        runoff rate. Works only if discharge field is set to ``"water__unit_flux_in"``.

        Parameters
        ----------
        new_runoff : str or array_like
            New runoff rate.
        """
        if self._discharge_field != "water__unit_flux_in":
            ValueError(
                "The SharedStreamPower's update_runoff method can only be used"
                "when discharge field is set to water__unit_flux_in (got"
                f" {self._discharge_field})"
            )

        self._runoff_rate = new_runoff

        self._K = self.k_bedrock / self.runoff_rate**self.m_sp
        self._v_s = self.k_bedrock * self.runoff_rate / self.k_transport
        np.multiply(
            self.runoff_rate,
            self._grid.at_node["drainage_area"],
            out=self._grid.at_node["water__unit_flux_in"],
        )
        self._q = return_array_at_node(self._grid, self._discharge_field)
