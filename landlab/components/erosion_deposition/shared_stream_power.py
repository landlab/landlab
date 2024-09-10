from landlab.components.erosion_deposition.erosion_deposition import ErosionDeposition
from landlab.components.erosion_deposition.generalized_erosion_deposition import (
    DEFAULT_MINIMUM_TIME_STEP,
)
from landlab.utils import return_array_at_node

TIME_STEP_FACTOR = 0.5


class SharedStreamPower(ErosionDeposition):
    r"""
    Implements the Shared Stream Power Model in the style of Hergarten (2021).
    Designed to simultaneously model
    river incision and sediment transport in order to seamlessly
    transition between detachment limited to transport limited erosion.
    Mathematically equivalent to the linear decline model from Davy and Lague (2009),
    which is used by the base class, ErosionDeposition. In addition, this component is
    designed to work with varying runoff rates, and can update the discharge
    and other parameters effected by discharge with each timestep.

    Here is the equation for erosion without a threshold:

        E = K_d * A^m_sp * S^n_sp - K_d/K_t * Qs / A

    where Q is water discharge, Qs is sediment flux, S is slope, m_sp and n_sp are scaling
    exponents, coefficient K_d is the erodibility of the bedrock and coefficient K_t is the
    ability to transport sediment.

    The first term, K_d * A ** m_sp * S ** n_sp, is the incision term, and the second term,
    K_d/K_t * Qs / A, is the transport term. Note that K_d/K_t determines the relative amount of
    incision and sediment transport. K_d modifies the incision term.

    The equivalent equation used by ErosionDeposition from Davy & Lague (2009) is:

        E = K * q^m_sp * S^n_sp - v_s * Qs / q

    Where K is sediment erodibility, v_s is the settling velocity for sediment, and q is the
    water discharge.

    To translate the shared stream power input for ErosionDeposition, we use the equations:

        q = Ar
        K = K_d / r^m_sp
        v_s = K_d/K_t

    It is important to note that the second two equations were derived only
    for calibrating the model,and do not necessarily correlate to landscape evolution
    in nature.

    To write the final equation we define the incision term as omega:

        omega = K_d * A ** m_{sp} * S ** {n_sp}

    And incorporate sp_crit, the critical stream power needed to erode bedrock, giving:

        E = omega(1 - exp(omega / sp_crit) ) - K_d/K_t * Qs / A


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

    """

    _name = "SharedStreamPower"

    _unit_agnostic = True

    def __init__(
        self,
        grid,
        K_d=0.001,
        K_t=0.001,
        r=1.0,
        m_sp=0.5,
        n_sp=1.0,
        sp_crit=0.0,
        F_f=0.0,
        discharge_field="surface_water__discharge",
        solver="basic",
        **kwds,
    ):
        """Initialize the Shared Stream Power model.

        Parameters
        ----------
        grid : ModelGrid
            Landlab ModelGrid object
        K_d : float, field name, or array
            Erodibility for bedrock (units vary).
        K_t : float, field name, or array
            Ability to transport sediment (units vary).
        r : float
            Runoff rate. Scales Q = Ar. Default is 1.0 [m/yr]
        m_sp : float
            Discharge exponent (units vary). Default is 0.5.
        n_sp : float
            Slope exponent (units vary). Default is 1.0.
        sp_crit : float, field name, or array
            Critical stream power to erode substrate [E/(TL^2)]
        F_f : float
            Fraction of eroded material that turns into "fines" that do not
            contribute to (coarse) sediment load. Defaults to zero.
        discharge_field : float, field name, or array
            Discharge [L^2/T]. The default is to use the grid field
            'surface_water__discharge', which is simply drainage area
            multiplied by the default rainfall rate (1 m/yr). To use custom
            spatially/temporally varying rainfall, use 'water__unit_flux_in'
            to specify water input to the FlowAccumulator.
        solver : string
            Solver to use. Options at present include:
                (1) 'basic' (default): explicit forward-time extrapolation.
                    Simple but will become unstable if time step is too large.
                (2) 'adaptive': adaptive time-step solver that estimates a
                    stable step size based on the shortest time to "flattening"
                    among all upstream-downstream node pairs.

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
            -5x5 grid with baselevel in the lower left corner
            -all other boundary nodes closed
            -Initial topography is plane tilted up to the upper right + noise

        >>> nr = 5
        >>> nc = 5
        >>> dx = 10
        >>> mg = RasterModelGrid((nr, nc), xy_spacing=10.0)
        >>> _ = mg.add_zeros("node", "topographic__elevation")
        >>> mg.at_node["topographic__elevation"] += (
        ...     mg.node_y / 10 + mg.node_x / 10 + np.random.rand(len(mg.node_y)) / 10
        ... )
        >>> mg.set_closed_boundaries_at_grid_edges(
        ...     bottom_is_closed=True,
        ...     left_is_closed=True,
        ...     right_is_closed=True,
        ...     top_is_closed=True,
        ... )
        >>> mg.set_watershed_boundary_condition_outlet_id(
        ...     0, mg.at_node["topographic__elevation"], -9999.0
        ... )
        >>> fsc_dt = 100.0
        >>> ed_dt = 1.0

        Check initial topography

        >>> mg.at_node["topographic__elevation"].reshape(mg.shape)
        array([[0.02290479, 1.03606698, 2.0727653 , 3.01126678, 4.06077707],
               [1.08157495, 2.09812694, 3.00637448, 4.07999597, 5.00969486],
               [2.04008677, 3.06621577, 4.09655859, 5.04809001, 6.02641123],
               [3.05874171, 4.00585786, 5.0595697 , 6.04425233, 7.05334077],
               [4.05922478, 5.0409473 , 6.07035008, 7.0038935 , 8.01034357]])

        Instantiate Fastscape eroder, flow router, and depression finder

        >>> fr = FlowAccumulator(mg, flow_director="D8")
        >>> df = DepressionFinderAndRouter(mg)
        >>> fsc = FastscapeEroder(mg, K_sp=0.001, m_sp=0.5, n_sp=1)

        Burn in an initial drainage network using the Fastscape eroder:

        >>> for x in range(100):
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     flooded = np.where(df.flood_status == 3)[0]
        ...     fsc.run_one_step(dt=fsc_dt)
        ...     mg.at_node["topographic__elevation"][0] -= 0.001  # uplift
        ...

        Instantiate the SharedStreamPower component:

        >>> ssp = SharedStreamPower(
        ...     mg, K_d=0.00001, K_t=0.001, m_sp=0.5, n_sp=1.0, sp_crit=0
        ... )

        Now run the E/D component for 2000 short timesteps:

        >>> for x in range(2000):  # E/D component loop
        ...     fr.run_one_step()
        ...     df.map_depressions()
        ...     ssp.run_one_step(dt=ed_dt)
        ...     mg.at_node["topographic__elevation"][0] -= 2e-4 * ed_dt
        ...

        Now we test to see if topography is right:

        >>> np.around(mg.at_node["topographic__elevation"], decimals=3).reshape(
        ...     mg.shape
        ... )
        array([[-0.477,  1.036,  2.073,  3.011,  4.061],
               [ 1.082, -0.08 , -0.065, -0.054,  5.01 ],
               [ 2.04 , -0.065, -0.065, -0.053,  6.026],
               [ 3.059, -0.054, -0.053, -0.035,  7.053],
               [ 4.059,  5.041,  6.07 ,  7.004,  8.01 ]])
        """

        self.discharge_field = discharge_field
        self.r = r
        self.K_d = K_d
        self.K_t = K_t
        self.m_sp = m_sp

        # convert shared stream power inputs to erosion deposition inputs
        vs_ = self.K_d * self.r / self.K_t
        K_s = self.K_d / self.r**self.m_sp

        # instantiate ErosionDeposition
        super().__init__(
            grid,
            K=K_s,
            v_s=vs_,
            m_sp=m_sp,
            n_sp=n_sp,
            sp_crit=sp_crit,
            F_f=F_f,
            discharge_field=self.discharge_field,
            solver=solver,
            dt_min=DEFAULT_MINIMUM_TIME_STEP,
            **kwds,
        )

    def update_runoff(self, new_r=1.0):
        """
        Updates r, K, v_s, and "water__unit_flux_in" for a new runoff rate. Works only if
        discharge field is set to "water__unit_flux_in."

        Parameters
        ----------
        new_r : float, field name, or array
            Updated runoff rate. Default is 1.0.

        """
        if not self.discharge_field == "water__unit_flux_in":
            ValueError(
                "The SharedStreamPower's update_runoff method can only be used"
                "when discharge field is set to water__unit_flux_in"
            )

        self.r = new_r
        self.K = self.K_d / self.r**self.m_sp
        self.v_s = self.K_d * self.r / self.K_t
        self._grid.at_node["water__unit_flux_in"] = (
            self.r * self._grid.at_node["drainage_area"]
        )
        self._q = return_array_at_node(self._grid, self.discharge_field)
