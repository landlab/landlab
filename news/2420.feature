Modernised the `gFlex` component for the gFlex v2.x API with several
breaking changes: output field renamed to `lithosphere__vertical_displacement`
(total deflection; the CSDMS standard name used by the gFlex BMI); load
input field renamed to `load__normal_component_of_stress` (material-agnostic
normal stress; callers convert ice/water/sediment thickness to Pa before
calling `run_one_step`); optional `lithosphere__elastic_thickness` node field
added for BMI-compatible T_e updates between coupling steps; BC parameters
renamed `bc_west/east/north/south`; default BC changed to `no_outside_loads`
(infinite-plate assumption); `topographic__elevation` interaction removed
(applying deflection to topography is now the caller's responsibility);
`Solver` and `PlateSolutionType` parameters removed (no longer supported by
gFlex v2); LU factorisation cached across calls for scalar T_e; BC validation
delegates to `gflex.VALID_BC_STRINGS_2D`, automatically accepting new aliases.
Requires `gflex>=2.0.0`.