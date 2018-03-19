# Change Log
All notable changes to landlab will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

This file was auto-generated using `scripts/make_changelog.py`.


## [HEAD] 2018-03-19

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]

## [v1.3.0] 2018-03-14

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* Clean up Sphinx documentation

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.2.0] 2017-10-19

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.1.0] 2017-06-26

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.3] 2017-03-04

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.2] 2016-11-24

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.1] 2016-08-25

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/patch-flowaccum-reclimit [#375]
* Merge remote-tracking branch 'refs/remotes/origin/master' into release
* Merge remote-tracking branch 'refs/remotes/origin/master' into release

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.8] 2016-07-07

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.7] 2016-07-07

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.6] 2016-06-16

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.5] 2016-06-14

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.4] 2016-06-14

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/pot-fr-modernise [#344]
* landlab/mcflugen/fix-landlab-test-function [#345]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.3] 2016-06-06

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.2] 2016-06-06

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.12] 2016-08-02

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/patch-flowaccum-reclimit [#375]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.11] 2016-07-19

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* Merge remote-tracking branch 'refs/remotes/origin/master' into release

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.10] 2016-07-14

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* Merge remote-tracking branch 'refs/remotes/origin/master' into release

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta.1] 2016-05-18

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-beta] 2016-05-13

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0-alpha] 2016-05-13

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v1.0.0] 2016-08-19

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/patch-flowaccum-reclimit [#375]
* Merge remote-tracking branch 'refs/remotes/origin/master' into release
* Merge remote-tracking branch 'refs/remotes/origin/master' into release
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.5.0] 2016-03-12

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.2.3] 2016-03-10

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.2.2] 2016-02-09

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.2.1] 2016-02-08

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.2.0] 2015-12-20

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.41] 2015-12-13

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.40] 2015-12-13

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.39] 2015-10-28

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/SiccarPoint/D4routing [#176]
* landlab/sed_fill [#173]
* landlab/sed_fill [#168]
* landlab/mcflugen/add-gradients-module [#169]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.38] 2015-10-21

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.37] 2015-10-20

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.36] 2015-10-20

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.35] 2015-10-20

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.34] 2015-10-19

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.33] 2015-10-11

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.32] 2015-10-11

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.31] 2015-10-11

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.30] 2015-10-11

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.29] 2015-09-13

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.28] 2015-09-09

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]
* landlab/mcflugen/update-mappers [#142]
* landlab/fix-issue-128 [#129]
* landlab/mcflugen/add-netcdf-tests [#138]
* landlab/mcflugen/add-better-testing-docs [#137]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.27] 2015-07-29

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]
* landlab/mcflugen/update-mappers [#142]
* landlab/fix-issue-128 [#129]
* landlab/mcflugen/add-netcdf-tests [#138]
* landlab/mcflugen/add-better-testing-docs [#137]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.26] 2015-07-20

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]
* landlab/mcflugen/update-mappers [#142]
* landlab/fix-issue-128 [#129]
* landlab/mcflugen/add-netcdf-tests [#138]
* landlab/mcflugen/add-better-testing-docs [#137]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.25] 2015-07-16

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]
* landlab/mcflugen/update-mappers [#142]
* landlab/fix-issue-128 [#129]
* landlab/mcflugen/add-netcdf-tests [#138]
* landlab/mcflugen/add-better-testing-docs [#137]
* resolved conflicts in lake mapper
* ierge branch 'master' of https://github.com/landlab/landlab

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.24] 2015-07-12

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit
* landlab/mcflugen/fix-component-imports [#214]
* landlab/SiccarPoint/component-tests [#209]
* landlab/SiccarPoint/component-tests [#204]
* landlab/mcflugen/add-value-name-decorator [#199]
* landlab/SiccarPoint/component-introspection [#200]
* landlab/mcflugen/fix-voronoi-cell-areas [#202]
* landlab/mcflugen/fix-testing-script [#198]
* landlab/mcflugen/add-neighbor-at-node [#195]
* landlab/mcflugen/fix-appveyor-builds [#185]
* landlab/mcflugen/clean-up-imshow [#178]
* landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182]
* landlab/SiccarPoint/var-doc [#187]
* landlab/gtucker/cleanup [#186]
* landlab/mcflugen/add-gradients-module [#169]
* landlab/SiccarPoint/delete-fields [#167]
* landlab/voronoi_stream_power [#158]
* landlab/mcflugen/make-node-status-private [#152]
* landlab/mcflugen/code-clean-up [#148]
* landlab/mcflugen/add-appveyor-slack-notifications [#149]
* Merge remote-tracking branch 'origin/link_status'
* landlab/link_status [#143]
* landlab/mcflugen/fix-failing-docs-build [#145]
* landlab/mcflugen/update-mappers [#142]
* landlab/fix-issue-128 [#129]
* landlab/mcflugen/add-netcdf-tests [#138]
* landlab/mcflugen/add-better-testing-docs [#137]
* resolved conflicts in lake mapper
* ierge branch 'master' of https://github.com/landlab/landlab
* landlab/mcflugen/add-python-3-support-with-six [#121]
* landlab/mcflugen/fix-setuptools-error [#126]

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

## [v0.1.23] tag

### Removed
* Removed inlink and outlink matrices.
* Removed deprecated raster_steepest_descent module.
* Removed corner_node_at_cell
* Removed old and unused _route_flow_old from lake_mapper
* Removed unused code from flow_direction_DN

### Fixed
* Fixed incorrect signatures of decorated methods in docs.
* Fixed Travis build errors with Python version conflicts.
* Fixed values not being cached (#614)
* Fixed component documentation not building (issue #575)
* Fixed netcdf4 import error

### Changed
* landlab/barnhark/updates_to_channel_profile [#637]
* landlab/barnhark/typo_in_imshow [#636]
* landlab/barnhark/add_component_docs [#634]
* landlab/gt-ca-uplift [#581]
* landlab/barnhark/make_stream_profiler [#605]
* landlab/mcflugen/remove-old-flux-div [#619]
* Simplified continuous integration and versioning.
* landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions [#612]
* landlab/barnhark/fix_stream_power_type_check [#610]
* Clean up API for diagonals.
* landlab/gt-taylor-fix [#606]
* landlab/mcflugen/fix-travis-ioerror [#607]
* landlab/barnhark/depth_dependent_boundary_conditions [#601]
* landlab/mcflugen/tidy-green-ampt [#591]
* landlab/barnhark/improving_cubic_flux [#582]
* Clean up Sphinx documentation
* landlab/margauxmouchene/test2 [#546]
* landlab/gt-fastscape-q [#574]
* amanaster2/master [#572]
* landlab/barnhark/kwargs_depth_dependent_diffuser [#553]
* landlab/gt-lattice-uplifter [#539]
* landlab/gt-add-phi-to-space-adaptive [#551]
* landlab/barnhark/cubic_nl_diffuser_kwargs [#550]
* landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548]
* landlab/gt-cmap-in-hexplot [#544]
* landlab/SiccarPoint/uniform_precip [#517]
* landlab/mcflugen/fix-greenampt-issue-530 [#535]
* landlab/mcflugen/add-logging-function [#504]
* landlab/gt-try-dyn-ts-space [#529]
* landlab/barnhark/get_set_state_methods_for_grid [#525]
* landlab/fixing_small_bug_in_erosion_deposition [#528]
* landlab/barnhark/eroder_depo_with_n_less_than_one [#523]
* landlab/barnhark/cubic_timestepper [#519]
* landlab/barnhark/addressing_brent_method_index_error [#510]
* landlab/gt-edit-erodep [#516]
* cmshobe/cmshobe/make-erosion-deposition-component [#511]
* landlab/barnhark/lake_mapper_faster [#512]
* nathanlyons/master [#505]
* cmshobe/cmshobe/minor_fixes_to_space [#509]
* cmshobe/cmshobe/change-hybrid-to-SPACE [#506]
* cmshobe/cmshobe/fix-hybrid-q-mechanics [#502]
* RondaStrauch/master [#480]
* landlab/barnhark/use_newton_fastscape [#492]
* landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
* landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497]
* landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
* landlab/barnhark/cython_hybrid_alluviaum [#494]
* cmshobe/cmshobe/fix_hybrid_q_options [#488]
* landlab/barnhark/smallchangestohybrid [#487]
* landlab/gt-add-stretched-expo [#485]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#481]
* landlab/mcflugen/add-graph-class [#477]
* landlab/barnhark/accumulator_efficiency [#476]
* landlab/barnhark/making_flow_accumulator_faster [#474]
* landlab/barnhark/fixing_kinwave_flow_issue [#471]
* cmshobe/cmshobe_fixes_to_hybrid_alluv [#469]
* landlab/gt-implicit-kinwave [#464]
* cmshobe/cmshobe/make_hybrid_alluv_initis [#467]
* Glader011235/master [#465]
* landlab/nicgaspar/diffusion_not_depositing [#463]
* landlab/kbarnhart/make_raster_netcdf [#462]
* cmshobe/cmshobe_hybrid_alluvium_model [#461]
* cmshobe/cmshobe_hybrid_alluvium_model [#460]
* Merge remote-tracking branch 'origin/master'
* Merge remote-tracking branch 'origin/SiccarPoint/pot-fr'
* landlab/kbarnhart/consistent_parameter_names [#459]
* landlab/gt-stream-power-K [#457]
* landlab/gt-fix-fastscape-variable-k [#456]
* landlab/gt-create-depth-dep-cubic-diffuser [#452]
* landlab/mcflugen/add-py36-builds [#453]
* landlab/kbarnhart/stream_power_error [#450]
* landlab/gt-fix-issue-448 [#449]
* landlab/mcflugen/fix-issue-428 [#447]
* landlab/jadams15/depth_slope_product [#445]
* landlab/SiccarPoint/fix_429 [#430]
* landlab/SiccarPoint/add-docs [#442]
* landlab/gt-fix-issue-431 [#433]
* landlab/gt-add-Q-stream-power-smooth-thresh [#443]
* landlab/SiccarPoint/auto-build-docs [#437]
* landlab/jadams15/spatially_variable_roughness [#438]
* landlab/kbarnhart/make_nd_fields [#434]
* landlab/kbarnhart/improvements_to_set_watershed_boundary [#426]
* landlab/gt-float64-fastscape [#427]
* landlab/gt-more-cts-cython [#378]
* landlab/gt-smooth-threshold-stream-power [#418]
* landlab/gt-tweak-cubic-diff [#416]
* landlab/gt-fix-init_typo [#415]
* landlab/jk-move-old-rst [#412]
* landlab/gt-merge-rg-cubic [#414]
* cmshobe/cmshobe-drainage-density [#398]
* fix minor conflict in raster.py
* landlab/SiccarPoint/grid_docs [#329]
* landlab/SiccarPoint/diagonal_link_lengths [#328]
* landlab/mcflugen/remove-deprecations [#327]
* landlab/SiccarPoint/imshow_grid_returns_im [#326]
* landlab/SiccarPoint/last-minute-deprecation [#324]
* landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323]
* landlab/SiccarPoint/patch-methods [#322]
* landlab/SiccarPoint/tweak-plotter [#321]
* landlab/saisiddu/Version_1_final [#320]
* landlab/SiccarPoint/modernise-field-names [#319]
* landlab/SiccarPoint/modernise-components [#314]
* landlab/SiccarPoint/most-egregious-diagonals [#315]
* landlab/gt-calc-of-to-at [#316]
* landlab/saisiddu/Version_1_final [#317]
* landlab/jadams15/uniform_precip_changes [#310]
* landlab/saisiddu/Version_1 [#311]
* landlab/mcflugen/moved-slope-methods [#313]
* landlab/SiccarPoint/Horn-slope [#309]
* landlab/mcflugen/remove-craters [#312]
* landlab/mcflugen/fix-docs-not-building [#304]
* landlab/SiccarPoint/grid_trawl [#307]
* landlab/nicgaspar/watershed_boundary_condition [#306]
* landlab/SiccarPoint/slopes [#305]
* landlab/gt-fix-diffuser-bug [#294]
* landlab/gt-update-gradients [#303]
* landlab/doc-component-reorg [#296]
* landlab/gt-fix-ca-tectonics [#297]
* landlab/gt-flux-divergence [#295]
* landlab/jk_cleanup_grid_docs [#289]
* landlab/SiccarPoint/fastscape-threshold [#290]
* landlab/SiccarPoint/component-modernisation [#288]
* landlab/gt_fix_faces_at_cell [#282]
* landlab/sed-flux-dep [#277]
* landlab/SiccarPoint/chi [#273]
* landlab/SiccarPoint/plotter_modernisation [#274]
* landlab/jk_rearrange_index [#275]
* landlab/SiccarPoint/steepness-index [#271]
* landlab/mcflugen/fix-issue-268 [#269]
* landlab/mcflugen/add-py35-support [#270]
* landlab/SiccarPoint/fix-issue-250 [#251]
* landlab/SiccarPoint/stream_power_standard [#256]
* landlab/mcflugen/fix-travis-not-running-all-tests [#265]
* landlab/SiccarPoint/dynamic-docstring-import [#258]
* landlab/SiccarPoint/enhance-mappers [#263]
* landlab/SiccarPoint/enhance-mappers [#262]
* Merged fix for deployment from AppVeyor to PyPI.
* landlab/SiccarPoint/enhance-mappers [#255]
* landlab/jk_copy_init_docstring [#248]
* landlab/gtucker/node_link_connectivity [#253]
* landlab/mcflugen/add-bmi-bridge [#246]
* landlab/gt-handle-flooded-nodes-in-stream-power [#247]
* landlab/jk_include_init_docstring [#244]
* landlab/mcflugen/fix-issue-242 [#243]
* Changed to deploy on release branch.
* landlab/SiccarPoint/fix-issue-237 [#239]
* landlab/mcflugen/fix-flexure-init [#231]
* landlab/jadams15/fix_node_links [#234]
* merge commit

### Added
* Added additional tests for SoilInfiltrationGreenAmpt.
* Added citation tracker for components.
* Added nodes_at_link attribute to ModelGrid.
* Added event layers to track stratigraphy

