=====================
Changelog for landlab
=====================


2.3.1 (unreleased)
------------------

- Nothing changed yet.


2.3.0 (2021-03-19)
------------------

- Fixed a bug when adding a missing at_grid field when testing components (#1286)

- Cleaned up landlab metadata files (#1285)

- Removed versioneer, we'll use zest.releaser from now on the manage versions (#1285)

- Added release and prerelease github actions (#1275)

- Fixed a bug in the FlowAccumulator to update pit present logic to also include node
  flood status (#1277)

- Fixed an error in the streampower notebook (#1260)

- Added building and testing of landlab with Python 3.9 (#1274)

- Added additional references for 2020/21 (#1273)

- Added simple submarine diffuser component (#1269)

- Changed to use GitHub Actions for CI (#1270)

- Added tutorial for river input to LEMs (#1258)

- Added a tidal flow component (#1225)

- Added a reference to the papers and presentations list (#1255)

- Fixed bug in Flow router/depression finder which incorrectly directed flow (#1248)

- Added ExponentialWeathererIntegrated component (#1254)

- Added docs for the flow accumulator (#1251)

- Added Wavefront OBJ output (#1241)

- Added two 2020 publications (#1243)

- Small changes to JOSS paper (#1237)

- Added recharge to the GroundwaterDupuitPercolator callback (#1223)

- NetworkSedimentTtransporter JOSS paper fixes (#1235)

- Fixed documentation bugs (#1233)

1.5.1 (2018-06-22)
------------------

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]

1.5.0 (2018-06-18)
------------------

.. _fixed-1:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]

.. _changed-1:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]

1.4.0 (2018-05-03)
------------------

.. _fixed-2:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]

.. _added-1:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]

.. _changed-2:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]

1.3.1 (2018-03-24)
------------------

.. _fixed-3:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]

.. _added-2:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]

.. _changed-3:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  Set versioneer to ignore ``v`` prefix in tags [Eric Hutton]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]

1.3.0 (2018-03-14)
------------------

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-4:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]

.. _added-3:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-4:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]

1.2.0 (2017-10-19)
------------------

.. _removed-1:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-5:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]

.. _added-4:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-5:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]

1.1.0 (2017-06-26)
------------------

.. _removed-2:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-6:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-5:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-6:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]

1.0.3 (2017-03-04)
------------------

.. _removed-3:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-7:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-6:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-7:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.2 (2016-11-24)
------------------

.. _removed-4:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-8:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-7:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-8:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.1 (2016-08-25)
------------------

.. _removed-5:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-9:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-8:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-9:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/patch-flowaccum-reclimit [#375] [Dan Hobley]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]

1.0.0-beta.8 (2016-07-07)
-------------------------

.. _removed-6:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-10:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-9:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-10:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.7 (2016-07-07)
-------------------------

.. _removed-7:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-11:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-10:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-11:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.6 (2016-06-16)
-------------------------

.. _removed-8:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-12:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-11:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-12:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.5 (2016-06-14)
-------------------------

.. _removed-9:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-13:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-12:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-13:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.4 (2016-06-14)
-------------------------

.. _removed-10:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-14:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-13:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-14:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/pot-fr-modernise [#344] [Jordan Adams]
-  landlab/mcflugen/fix-landlab-test-function [#345] [Dan Hobley]

1.0.0-beta.3 (2016-06-06)
-------------------------

.. _removed-11:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-15:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-14:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-15:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.2 (2016-06-06)
-------------------------

.. _removed-12:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-16:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-15:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-16:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0-beta.12 (2016-08-02)
--------------------------

.. _removed-13:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-17:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-16:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-17:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/patch-flowaccum-reclimit [#375] [Dan Hobley]

1.0.0-beta.11 (2016-07-19)
--------------------------

.. _removed-14:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-18:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-17:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-18:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]

1.0.0-beta.10 (2016-07-14)
--------------------------

.. _removed-15:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-19:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-18:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-19:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]

1.0.0-beta.1 (2016-05-18)
-------------------------

.. _removed-16:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-20:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-19:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-20:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]

1.0.0-beta (2016-05-13)
-----------------------

.. _removed-17:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-21:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-20:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-21:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]

1.0.0-alpha (2016-05-13)
------------------------

.. _removed-18:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-22:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-21:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-22:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]

1.0.0 (2016-08-19)
------------------

.. _removed-19:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-23:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-22:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-23:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/patch-flowaccum-reclimit [#375] [Dan Hobley]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]
-  Merge remote-tracking branch ‘refs/remotes/origin/master’ into
   release [saisiddu]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]

0.5.0 (2016-03-12)
------------------

.. _removed-20:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-24:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-23:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-24:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.2.3 (2016-03-10)
------------------

.. _removed-21:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-25:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-24:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-25:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.2.2 (2016-02-09)
------------------

.. _removed-22:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-26:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-25:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-26:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.2.1 (2016-02-08)
------------------

.. _removed-23:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-27:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-26:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-27:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.2.0 (2015-12-20)
------------------

.. _removed-24:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-28:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-27:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-28:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.1.41 (2015-12-13)
-------------------

.. _removed-25:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-29:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-28:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-29:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]

0.1.40 (2015-12-13)
-------------------

.. _removed-26:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-30:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-29:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-30:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]

0.1.39 (2015-10-28)
-------------------

.. _removed-27:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-31:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-30:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-31:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/SiccarPoint/D4routing [#176] [Jordan Adams]
-  landlab/sed_fill [#173] [Eric Hutton]
-  landlab/sed_fill [#168] [Jordan Adams]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]

0.1.38 (2015-10-21)
-------------------

.. _removed-28:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-32:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-31:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-32:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]

0.1.37 (2015-10-20)
-------------------

.. _removed-29:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-33:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-32:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-33:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]

0.1.36 (2015-10-20)
-------------------

.. _removed-30:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-34:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-33:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-34:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]

0.1.35 (2015-10-20)
-------------------

.. _removed-31:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-35:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-34:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-35:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]

0.1.34 (2015-10-19)
-------------------

.. _removed-32:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-36:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-35:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-36:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]

0.1.33 (2015-10-11)
-------------------

.. _removed-33:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-37:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-36:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-37:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]

0.1.32 (2015-10-11)
-------------------

.. _removed-34:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-38:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-37:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-38:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]

0.1.31 (2015-10-11)
-------------------

.. _removed-35:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-39:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-38:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-39:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]

0.1.30 (2015-10-11)
-------------------

.. _removed-36:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-40:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-39:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-40:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]

0.1.29 (2015-09-13)
-------------------

.. _removed-37:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-41:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-40:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-41:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]

0.1.28 (2015-09-09)
-------------------

.. _removed-38:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-42:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-41:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-42:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]
-  landlab/mcflugen/update-mappers [#142] [Eric Hutton]
-  landlab/fix-issue-128 [#129] [Eric Hutton]
-  landlab/mcflugen/add-netcdf-tests [#138] [Eric Hutton]
-  landlab/mcflugen/add-better-testing-docs [#137] [Eric Hutton]

0.1.27 (2015-07-29)
-------------------

.. _removed-39:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-43:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-42:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-43:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]
-  landlab/mcflugen/update-mappers [#142] [Eric Hutton]
-  landlab/fix-issue-128 [#129] [Eric Hutton]
-  landlab/mcflugen/add-netcdf-tests [#138] [Eric Hutton]
-  landlab/mcflugen/add-better-testing-docs [#137] [Eric Hutton]

0.1.26 (2015-07-20)
-------------------

.. _removed-40:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-44:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-43:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-44:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]
-  landlab/mcflugen/update-mappers [#142] [Eric Hutton]
-  landlab/fix-issue-128 [#129] [Eric Hutton]
-  landlab/mcflugen/add-netcdf-tests [#138] [Eric Hutton]
-  landlab/mcflugen/add-better-testing-docs [#137] [Eric Hutton]

0.1.25 (2015-07-16)
-------------------

.. _removed-41:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-45:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-44:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-45:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]
-  landlab/mcflugen/update-mappers [#142] [Eric Hutton]
-  landlab/fix-issue-128 [#129] [Eric Hutton]
-  landlab/mcflugen/add-netcdf-tests [#138] [Eric Hutton]
-  landlab/mcflugen/add-better-testing-docs [#137] [Eric Hutton]
-  resolved conflicts in lake mapper [gregtucker]
-  ierge branch ‘master’ of https://github.com/landlab/landlab
   [gregtucker]

0.1.24 (2015-07-12)
-------------------

.. _removed-42:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-46:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-45:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-46:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
-  landlab/mcflugen/fix-component-imports [#214] [Dan Hobley]
-  landlab/SiccarPoint/component-tests [#209] [Eric Hutton]
-  landlab/SiccarPoint/component-tests [#204] [Eric Hutton]
-  landlab/mcflugen/add-value-name-decorator [#199] [Dan Hobley]
-  landlab/SiccarPoint/component-introspection [#200] [Jenny Knuth]
-  landlab/mcflugen/fix-voronoi-cell-areas [#202] [Dan Hobley]
-  landlab/mcflugen/fix-testing-script [#198] [Eric Hutton]
-  landlab/mcflugen/add-neighbor-at-node [#195] [saisiddu]
-  landlab/mcflugen/fix-appveyor-builds [#185] [Eric Hutton]
-  landlab/mcflugen/clean-up-imshow [#178] [saisiddu]
-  landlab/mcflugen/fix-reset_lists_of_nodes_cells [#182] [Greg Tucker]
-  landlab/SiccarPoint/var-doc [#187] [Eric Hutton]
-  landlab/gtucker/cleanup [#186] [Eric Hutton]
-  landlab/mcflugen/add-gradients-module [#169] [Eric Hutton]
-  landlab/SiccarPoint/delete-fields [#167] [Dan Hobley]
-  landlab/voronoi_stream_power [#158] [Dan Hobley]
-  landlab/mcflugen/make-node-status-private [#152] [Eric Hutton]
-  landlab/mcflugen/code-clean-up [#148] [Eric Hutton]
-  landlab/mcflugen/add-appveyor-slack-notifications [#149] [Eric
   Hutton]
-  Merge remote-tracking branch ‘origin/link_status’ [mcflugen]
-  landlab/link_status [#143] [Eric Hutton]
-  landlab/mcflugen/fix-failing-docs-build [#145] [Eric Hutton]
-  landlab/mcflugen/update-mappers [#142] [Eric Hutton]
-  landlab/fix-issue-128 [#129] [Eric Hutton]
-  landlab/mcflugen/add-netcdf-tests [#138] [Eric Hutton]
-  landlab/mcflugen/add-better-testing-docs [#137] [Eric Hutton]
-  resolved conflicts in lake mapper [gregtucker]
-  ierge branch ‘master’ of https://github.com/landlab/landlab
   [gregtucker]
-  landlab/mcflugen/add-python-3-support-with-six [#121] [Eric Hutton]
-  landlab/mcflugen/fix-setuptools-error [#126] [Eric Hutton]

v0.1.23 (2015-05-15)
--------------------

.. _removed-43:

Removed
```````

-  Removed inlink and outlink matrices. [Eric Hutton]
-  Removed deprecated raster_steepest_descent module. [Eric Hutton]
-  Removed corner_node_at_cell [Eric Hutton]
-  Removed old and unused \_route_flow_old from lake_mapper [Eric
   Hutton]
-  Removed unused code from flow_direction_DN [Eric Hutton]

.. _fixed-47:

Fixed
`````

-  Fixed bug in Flexure1D when using “flexure” method [Eric Hutton]
-  Fixed unit test failures related to masked arrays (#710) [Eric
   Hutton]
-  Fixed failed Travis builds being reported as passing [Eric Hutton]
-  Fixed doctest for graph.adjacent_nodes_at_node [Eric Hutton]
-  Fixed names of packages deployed to Anaconda Cloud [Eric Hutton]
-  Fixed incorrect signatures of decorated methods in docs. [Eric
   Hutton]
-  Fixed Travis build errors with Python version conflicts. [Eric
   Hutton]
-  Fixed values not being cached (#614) [Eric Hutton]
-  Fixed component documentation not building (issue #575) [Eric Hutton]
-  Fixed netcdf4 import error [Eric Hutton]

.. _added-46:

Added
`````

-  Added CONTRIBUTING.md document [Eric Hutton]
-  Added script to create a nicely formatted changelog [Eric Hutton]
-  Added 1D Flexure component [Eric Hutton]
-  Added cite_as function to get landlab component citations [Eric
   Hutton]
-  Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
   [Eric Hutton]
-  Added additional tests for SoilInfiltrationGreenAmpt. [Eric Hutton]
-  Added citation tracker for components. [Eric Hutton]
-  Added nodes_at_link attribute to ModelGrid. [Eric Hutton]
-  Added event layers to track stratigraphy [Eric Hutton]

.. _changed-47:

Changed
```````

-  amanaster2/master [#733] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/major_cleanup_to_space_and_erodepo_init [#709] [Katy
   Barnhart]
-  landlab/gt/fix-doctest-issue-726 [#728] [Greg Tucker]
-  landlab/gt/ca-top-hit-bug [#720] [Greg Tucker]
-  landlab/barnhark/space_cell_area [#719] [Greg Tucker]
-  landlab/barnhark/use_field_name_array_or_float [#683] [Katy Barnhart]
-  landlab/barnhark/give_hex_models_watershed_methods [#685] [Katy
   Barnhart]
-  landlab/SiccarPoint/fix-issue-702 [#706] [Katy Barnhart]
-  Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility [#658]
   [Katy Barnhart]
-  landlab/barnhark/revert_channel_profiler [#695] [Katy Barnhart]
-  landlab/barnhark/space_rounding [#698] [Katy Barnhart]
-  landlab/barnhark/add_docs_to_normal_fault [#677] [Katy Barnhart]
-  landlab/barnhark/space_type_updates [#669] [Katy Barnhart]
-  landlab/barnhark/minor_changes_to_normal_fault [#663] [Katy Barnhart]
-  landlab/gt-debug-ca-propswap [#661] [Greg Tucker]
-  landlab/barnhark/space_hex [#655] [Katy Barnhart]
-  landlab/barnhark/add_kwargs [#645] [Katy Barnhart]
-  landlab/barnhark/normal_fault_kwargs [#649] [Katy Barnhart]
-  landlab/barnhark/normal_fault [#640] [Katy Barnhart]
-  landlab/barnhark/exponential_weatherer_docstring [#643] [Katy
   Barnhart]
-  landlab/nathanlyons/watershed [#545] [Nathan Lyons]
-  landlab/barnhark/updates_to_channel_profile [#637] [Katy Barnhart]
-  landlab/barnhark/typo_in_imshow [#636] [Katy Barnhart]
-  landlab/barnhark/add_component_docs [#634] [Katy Barnhart]
-  landlab/gt-ca-uplift [#581] [Greg Tucker]
-  landlab/barnhark/make_stream_profiler [#605] [Katy Barnhart]
-  landlab/mcflugen/remove-old-flux-div [#619] [Dan Hobley]
-  Simplified continuous integration and versioning. [Eric Hutton]
-  landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
   [#612] [Katy Barnhart]
-  landlab/barnhark/fix_stream_power_type_check [#610] [Katy Barnhart]
-  Clean up API for diagonals. [Eric Hutton]
-  landlab/gt-taylor-fix [#606] [Katy Barnhart]
-  landlab/mcflugen/fix-travis-ioerror [#607] [Nathan Lyons]
-  landlab/barnhark/depth_dependent_boundary_conditions [#601] [Katy
   Barnhart]
-  landlab/mcflugen/tidy-green-ampt [#591] [Jordan Adams]
-  landlab/barnhark/improving_cubic_flux [#582] [Katy Barnhart]
-  Clean up Sphinx documentation [Eric Hutton]
-  landlab/margauxmouchene/test2 [#546] [margauxmouchene]
-  landlab/gt-fastscape-q [#574] [Greg Tucker]
-  amanaster2/master [#572] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/kwargs_depth_dependent_diffuser [#553] [Katy
   Barnhart]
-  landlab/gt-lattice-uplifter [#539] [Greg Tucker]
-  landlab/gt-add-phi-to-space-adaptive [#551] [Greg Tucker]
-  landlab/barnhark/cubic_nl_diffuser_kwargs [#550] [Katy Barnhart]
-  landlab/barnhark/no_kwargs_in_dd_cubic_diffuser [#548] [Katy
   Barnhart]
-  landlab/gt-cmap-in-hexplot [#544] [Greg Tucker]
-  landlab/SiccarPoint/uniform_precip [#517] [Dan Hobley]
-  landlab/mcflugen/fix-greenampt-issue-530 [#535] [Katy Barnhart]
-  landlab/mcflugen/add-logging-function [#504] [Eric Hutton]
-  landlab/gt-try-dyn-ts-space [#529] [Katy Barnhart]
-  landlab/barnhark/get_set_state_methods_for_grid [#525] [Greg Tucker]
-  landlab/fixing_small_bug_in_erosion_deposition [#528] [Greg Tucker]
-  landlab/barnhark/eroder_depo_with_n_less_than_one [#523] [Greg
   Tucker]
-  landlab/barnhark/cubic_timestepper [#519] [Greg Tucker]
-  landlab/barnhark/addressing_brent_method_index_error [#510] [Katy
   Barnhart]
-  landlab/gt-edit-erodep [#516] [Katy Barnhart]
-  cmshobe/cmshobe/make-erosion-deposition-component [#511] [Greg
   Tucker]
-  landlab/barnhark/lake_mapper_faster [#512] [Greg Tucker]
-  nathanlyons/master [#505] [Nicole M Gasparini]
-  cmshobe/cmshobe/minor_fixes_to_space [#509] [Katy Barnhart]
-  cmshobe/cmshobe/change-hybrid-to-SPACE [#506] [Katy Barnhart]
-  cmshobe/cmshobe/fix-hybrid-q-mechanics [#502] [Katy Barnhart]
-  RondaStrauch/master [#480] [Sai Siddhartha Nudurupati]
-  landlab/barnhark/use_newton_fastscape [#492] [Katy Barnhart]
-  landlab/barnhark/improve_streampower_smooth_thresh_stability [#499]
   [Greg Tucker]
-  landlab/barnhark/dynamic_timestep_cubic_flux_diffuser [#497] [Greg
   Tucker]
-  landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient [#490]
   [Katy Barnhart]
-  landlab/barnhark/cython_hybrid_alluviaum [#494] [Greg Tucker]
-  cmshobe/cmshobe/fix_hybrid_q_options [#488] [Katy Barnhart]
-  landlab/barnhark/smallchangestohybrid [#487] [Greg Tucker]
-  landlab/gt-add-stretched-expo [#485] [Katy Barnhart]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#481] [Katy Barnhart]
-  landlab/mcflugen/add-graph-class [#477] [Greg Tucker]
-  landlab/barnhark/accumulator_efficiency [#476] [Greg Tucker]
-  landlab/barnhark/making_flow_accumulator_faster [#474] [Greg Tucker]
-  landlab/barnhark/fixing_kinwave_flow_issue [#471] [Greg Tucker]
-  cmshobe/cmshobe_fixes_to_hybrid_alluv [#469] [Greg Tucker]
-  landlab/gt-implicit-kinwave [#464] [Greg Tucker]
-  cmshobe/cmshobe/make_hybrid_alluv_initis [#467] [Katy Barnhart]
-  Glader011235/master [#465] [Katy Barnhart]
-  landlab/nicgaspar/diffusion_not_depositing [#463] [Jordan Adams]
-  landlab/kbarnhart/make_raster_netcdf [#462] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#461] [Katy Barnhart]
-  cmshobe/cmshobe_hybrid_alluvium_model [#460] [Greg Tucker]
-  Merge remote-tracking branch ‘origin/master’ [SiccarPoint]
-  Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
   [SiccarPoint]
-  landlab/kbarnhart/consistent_parameter_names [#459] [Katy Barnhart]
-  landlab/gt-stream-power-K [#457] [Greg Tucker]
-  landlab/gt-fix-fastscape-variable-k [#456] [Katy Barnhart]
-  landlab/gt-create-depth-dep-cubic-diffuser [#452] [Katy Barnhart]
-  landlab/mcflugen/add-py36-builds [#453] [Eric Hutton]
-  landlab/kbarnhart/stream_power_error [#450] [Greg Tucker]
-  landlab/gt-fix-issue-448 [#449] [Dan Hobley]
-  landlab/mcflugen/fix-issue-428 [#447] [Jordan Adams]
-  landlab/jadams15/depth_slope_product [#445] [Jordan Adams]
-  landlab/SiccarPoint/fix_429 [#430] [Katy Barnhart]
-  landlab/SiccarPoint/add-docs [#442] [Katy Barnhart]
-  landlab/gt-fix-issue-431 [#433] [Dan Hobley]
-  landlab/gt-add-Q-stream-power-smooth-thresh [#443] [Katy Barnhart]
-  landlab/SiccarPoint/auto-build-docs [#437] [Dan Hobley]
-  landlab/jadams15/spatially_variable_roughness [#438] [Jordan Adams]
-  landlab/kbarnhart/make_nd_fields [#434] [Greg Tucker]
-  landlab/kbarnhart/improvements_to_set_watershed_boundary [#426] [Katy
   Barnhart]
-  landlab/gt-float64-fastscape [#427] [Greg Tucker]
-  landlab/gt-more-cts-cython [#378] [Greg Tucker]
-  landlab/gt-smooth-threshold-stream-power [#418] [Greg Tucker]
-  landlab/gt-tweak-cubic-diff [#416] [Greg Tucker]
-  landlab/gt-fix-init_typo [#415] [Greg Tucker]
-  landlab/jk-move-old-rst [#412] [Greg Tucker]
-  landlab/gt-merge-rg-cubic [#414] [Greg Tucker]
-  cmshobe/cmshobe-drainage-density [#398] [Katy Barnhart]
-  fix minor conflict in raster.py [Greg Tucker]
-  landlab/SiccarPoint/grid_docs [#329] [Dan Hobley]
-  landlab/SiccarPoint/diagonal_link_lengths [#328] [Eric Hutton]
-  landlab/mcflugen/remove-deprecations [#327] [Eric Hutton]
-  landlab/SiccarPoint/imshow_grid_returns_im [#326] [Dan Hobley]
-  landlab/SiccarPoint/last-minute-deprecation [#324] [Eric Hutton]
-  landlab/SiccarPoint/BAD-INDEX-is-minus1 [#323] [Eric Hutton]
-  landlab/SiccarPoint/patch-methods [#322] [Eric Hutton]
-  landlab/SiccarPoint/tweak-plotter [#321] [Eric Hutton]
-  landlab/saisiddu/Version_1_final [#320] [Eric Hutton]
-  landlab/SiccarPoint/modernise-field-names [#319] [Dan Hobley]
-  landlab/SiccarPoint/modernise-components [#314] [Eric Hutton]
-  landlab/SiccarPoint/most-egregious-diagonals [#315] [Dan Hobley]
-  landlab/gt-calc-of-to-at [#316] [Greg Tucker]
-  landlab/saisiddu/Version_1_final [#317] [Eric Hutton]
-  landlab/jadams15/uniform_precip_changes [#310] [Dan Hobley]
-  landlab/saisiddu/Version_1 [#311] [Dan Hobley]
-  landlab/mcflugen/moved-slope-methods [#313] [Dan Hobley]
-  landlab/SiccarPoint/Horn-slope [#309] [Eric Hutton]
-  landlab/mcflugen/remove-craters [#312] [Eric Hutton]
-  landlab/mcflugen/fix-docs-not-building [#304] [Dan Hobley]
-  landlab/SiccarPoint/grid_trawl [#307] [Eric Hutton]
-  landlab/nicgaspar/watershed_boundary_condition [#306] [Jordan Adams]
-  landlab/SiccarPoint/slopes [#305] [Dan Hobley]
-  landlab/gt-fix-diffuser-bug [#294] [Dan Hobley]
-  landlab/gt-update-gradients [#303] [Greg Tucker]
-  landlab/doc-component-reorg [#296] [Greg Tucker]
-  landlab/gt-fix-ca-tectonics [#297] [Greg Tucker]
-  landlab/gt-flux-divergence [#295] [Greg Tucker]
-  landlab/jk_cleanup_grid_docs [#289] [Greg Tucker]
-  landlab/SiccarPoint/fastscape-threshold [#290] [Jordan Adams]
-  landlab/SiccarPoint/component-modernisation [#288] [Greg Tucker]
-  landlab/gt_fix_faces_at_cell [#282] [Greg Tucker]
-  landlab/sed-flux-dep [#277] [Dan Hobley]
-  landlab/SiccarPoint/chi [#273] [Greg Tucker]
-  landlab/SiccarPoint/plotter_modernisation [#274] [Greg Tucker]
-  landlab/jk_rearrange_index [#275] [Greg Tucker]
-  landlab/SiccarPoint/steepness-index [#271] [nicgaspar]
-  landlab/mcflugen/fix-issue-268 [#269] [Dan Hobley]
-  landlab/mcflugen/add-py35-support [#270] [saisiddu]
-  landlab/SiccarPoint/fix-issue-250 [#251] [Dan Hobley]
-  landlab/SiccarPoint/stream_power_standard [#256] [Eric Hutton]
-  landlab/mcflugen/fix-travis-not-running-all-tests [#265] [Eric
   Hutton]
-  landlab/SiccarPoint/dynamic-docstring-import [#258] [Greg Tucker]
-  landlab/SiccarPoint/enhance-mappers [#263] [Dan Hobley]
-  landlab/SiccarPoint/enhance-mappers [#262] [Dan Hobley]
-  Merged fix for deployment from AppVeyor to PyPI. [mcflugen]
-  landlab/SiccarPoint/enhance-mappers [#255] [Greg Tucker]
-  landlab/jk_copy_init_docstring [#248] [Jenny Knuth]
-  landlab/gtucker/node_link_connectivity [#253] [Dan Hobley]
-  landlab/mcflugen/add-bmi-bridge [#246] [Greg Tucker]
-  landlab/gt-handle-flooded-nodes-in-stream-power [#247] [Dan Hobley]
-  landlab/jk_include_init_docstring [#244] [Jenny Knuth]
-  landlab/mcflugen/fix-issue-242 [#243] [Eric Hutton]
-  Changed to deploy on release branch. [mcflugen]
-  landlab/SiccarPoint/fix-issue-237 [#239] [Dan Hobley]
-  landlab/mcflugen/fix-flexure-init [#231] [Jordan Adams]
-  landlab/jadams15/fix_node_links [#234] [Eric Hutton]
-  merge commit [Jenny Knuth]
