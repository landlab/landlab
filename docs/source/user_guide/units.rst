.. _unit_agnostic:

How Landlab Is/Is Not Unit Agnostic
-----------------------------------

This page describes how Landlab handles units. Our approach is intended to
balance usability, effective communication of component expectations to users,
and a low barrier for developers.

Component inputs and outputs with units may be fields (arrays stored on the
grid), or arguments/keyword arguments to the component instantiation and run
functions.

All components require that a user is consistent between the space/time/mass
units used within one and across multiple components. For example, if you were
to couple a ``StreamPowerEroder`` and a ``LinearDiffuser`` your x, y, and z
coordinates, the field `topographic__elevation` and all input parameters would
need to share a common length unit. Further, input parameters with units that
include time and your time step would need to share a common time unit.

Components that require that you are consistent with units but do not care
whether you use feet or meters for your length unit are called "unit agnostic".
You can find out if a component is unit agnostic by querying the attribute:
``Component.unit_agnostic`` which will return ``True`` or ``False``.

Unit agnostic components will provide specific units in the field and parameter
metadata. However, no computation within a unit agnostic component assumes a
specific unit. If you provide a consistent set of inputs, you can use whichever
units you prefer. Note that this may require specifying ALL keyword arguments
(e.g., if gravitational acceleration as 9.81 m/s^2 is a default value, you must
provide 32.2 ft/s^2 as an input if you want to use feet as your length unit).

In contrast to unit-agnostic components, non-unit-agnostic components REQUIRE
that a specific set of units be used. This typically occurs because built into
the source code of the computation are assumptions or conversions about the
units of inputs.

When it doubt, the best approach is to open a GitHub issue.

Below is a list of non-unit agnostic components current as of February 2020:

* :py:class:`~landlab.components.GroundwaterDupuitPercolator`
* :py:class:`~landlab.components.KinwaveImplicitOverlandFlow`
* :py:class:`~landlab.components.KinwaveOverlandFlowModel`
* :py:class:`~landlab.components.LandslideProbability`
* :py:class:`~landlab.components.LateralEroder`
* :py:class:`~landlab.components.OverlandFlow`
* :py:class:`~landlab.components.OverlandFlowBates`
* :py:class:`~landlab.components.PotentialEvapotranspiration`
* :py:class:`~landlab.components.PotentialityFlowRouter`
* :py:class:`~landlab.components.Radiation`
* :py:class:`~landlab.components.SedDepEroder`
* :py:class:`~landlab.components.SoilMoisture`
* :py:class:`~landlab.components.SoilInfiltrationGreenAmpt`
* :py:class:`~landlab.components.SpatialPrecipitationDistribution`
* :py:class:`~landlab.components.VegCA`
* :py:class:`~landlab.components.Vegetation`
