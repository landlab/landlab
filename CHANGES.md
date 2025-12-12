# Release Notes

```{towncrier-draft-entries} Not yet released
```

<!-- towncrier-draft-entries:: Not yet released -->

<!-- towncrier release notes start -->

## 2.10.1 (2025-08-10)


### ✨ New Components

- Adds a new 2D flow solver to Landlab. RiverFlowDynamics implements a semi-implicit,
  semi-Lagrangian finite-volume approximation of the depth-averaged shallow water
  equations, originally proposed by Casulli and Cheng in 1992, along with
  subsequent related work. [#1979](https://github.com/landlab/landlab/issues/1979)
- Added ConcentrationTrackerForSpace companion component to the SpaceLargeScaleEroder. The component calculates the movement of sediment property concentrations by fluvial transport across the Landlab grid. [#2020](https://github.com/landlab/landlab/issues/2020)


### 🍰 New Features

- Updated `KinwaveImplicitOverlandFlow` to accept arrays of floats in
  addition to scalars for both the `runoff_rate` and `roughness` keywords. [#2098](https://github.com/landlab/landlab/issues/2098)
- Added `neighbor_to_arrow` function that converts neighbor node indices into
  arrow symbols. [#2178](https://github.com/landlab/landlab/issues/2178)


### 🛠️ Bug Fixes

- Fixed an issue that caused a `TypeError` to be raised when printing numpy arrays using
  older versions of numpy. [#2018](https://github.com/landlab/landlab/issues/2018)
- Fix dimensional error in calculation of Courant number. [#2163](https://github.com/landlab/landlab/issues/2163)


### 📖 Documentation Enhancements

- Update notebooks that demonstrate Mesa + Landlab. [#2032](https://github.com/landlab/landlab/issues/2032)
- Fixed an issue with the citation section of the docs that caused the
  tabs to not render correctly. [#2049](https://github.com/landlab/landlab/issues/2049)
- Moved *Lithology* text from the top-level reference page to the *layers* page,
  which describes different types of layering that *Landlab* supports. [#2050](https://github.com/landlab/landlab/issues/2050)
- Updated the links to the landlab documentation to point to landlab.csdms.io. [#2117](https://github.com/landlab/landlab/issues/2117)
- Fixed a number of broken links in the documentation. [#2118](https://github.com/landlab/landlab/issues/2118)
- Fixed some formatting in the developer installation instructions. [#2155](https://github.com/landlab/landlab/issues/2155)
- Update list of publications in "usedby" [#2170](https://github.com/landlab/landlab/issues/2170)


### 🔩 Other Changes and Additions

- Minor improvements to plot_network_and_parcels to allow straightforward compatability with imshowgrid. [#1986](https://github.com/landlab/landlab/issues/1986)
- Deprecated passing the `at` keyword as the first argument to
  functions like, for example, :meth:`~.GraphFields.add_ones`. [#1999](https://github.com/landlab/landlab/issues/1999)
- Added support for Python 3.13 and dropped support for 3.10. [#2017](https://github.com/landlab/landlab/issues/2017)
- Removed unnecessary toctree directives from the docs. [#2048](https://github.com/landlab/landlab/issues/2048)
- Fixed a *jinja* error when building the docs using newer versions of
  *sphinx*. [#2051](https://github.com/landlab/landlab/issues/2051)
- Removed the current module from *sphinx* autosummary tables to be compatible
  with newer version of *sphinx*. [#2052](https://github.com/landlab/landlab/issues/2052)
- Changed the *mixed-line-ending* *pre-commit* hook to ensure all files
  use a line feed as the end-of-line character. [#2083](https://github.com/landlab/landlab/issues/2083)
- Removed the benchmarks from our unit tests, as they'll be moved to a separate
  repostory. This reduces the amount of time our continuous integration
  tests take to run. [#2087](https://github.com/landlab/landlab/issues/2087)
- Added a new GitHub Actions job that checks the documentation for broken links. [#2121](https://github.com/landlab/landlab/issues/2121)
- Simplified the github actions job that builds the docs by using
  the setup-python action rather than setup-miniconda. [#2122](https://github.com/landlab/landlab/issues/2122)
- Improved the speed of the continuous integration tests by switching to using
  the *setup-python* action and running the *richdem* tests as their own job
  using a *conda* environment. [#2129](https://github.com/landlab/landlab/issues/2129)
- Added ubuntu arm64 runners to our continuous integration builds. This drastically
  speeds up building for arm64 as compared to using emulation via QEMU. [#2132](https://github.com/landlab/landlab/issues/2132)
- Changed from *coveralls* to *codecov* for tracking test coverage. [#2133](https://github.com/landlab/landlab/issues/2133)
- Added a mass balance check in the unit tests for deAlmeida OverlandFlow component. [#2166](https://github.com/landlab/landlab/issues/2166)
- Added support for parallel builds. The build process now utilizes
  multiple CPU cores where available, significantly reducing build time. [#2200](https://github.com/landlab/landlab/issues/2200)

## 2.9.1 (2024-10-10)


### ✨ New Components

- Added shared stream power component that models mixed bedrock-alluvial erosion in channels. [#1976](https://github.com/landlab/landlab/issues/1976)


### 🍰 New Features

- Added new, faster, aggregator functions for the {class}`~.DataRecord`. These new functions,
  {func}`~.aggregate_items_at_link_count`, {func}`~.aggregate_items_at_link_sum`, and
  {func}`~.aggregate_items_at_link_mean` are written in *cython* and are several orders
  of magnitude faster than the previous method of using *merge*/*group_by*. [#1648](https://github.com/landlab/landlab/issues/1648)
- Improved perfomance when creating graphs by using openmp for
  parallelism. [#1940](https://github.com/landlab/landlab/issues/1940)
- Added {func}`~.esri_ascii.dump` and {func}`~.esri_ascii.load` functions for reading and writing
  grids to/from ESRI ASCII format. [#1948](https://github.com/landlab/landlab/issues/1948)
- Added transport dependent bed material abrasion to the {class}`~.NetworkSedimentTransporter`, as in Pfeiffer et al. (2022). [#1952](https://github.com/landlab/landlab/issues/1952)


### 🛠️ Bug Fixes

- Fixed a bug in {class}`~.NetworkSedimentTransporter` where porosity was
  corrected for twice in calculated topographic change at nodes. [#1921](https://github.com/landlab/landlab/issues/1921)
- Update {class}`~.ConcentrationTrackerForDiffusion` run methods to to solve two situations that cause the mass balance to fail. [#1967](https://github.com/landlab/landlab/issues/1967)
- Fixed an issue when using :func:`~.jaggedarray.unravel` when using
  an `int32` array for offsets. [#2015](https://github.com/landlab/landlab/issues/2015)


### 📖 Documentation Enhancements

- Significantly reduced the amount of time needed to build the documentation. [#1987](https://github.com/landlab/landlab/issues/1987)
- Autogenerate Landlab's API reference docs using *sphinx-apidoc*. [#1989](https://github.com/landlab/landlab/issues/1989)
- Cleaned up the documentation by removing obsolete documentation and moving
  posters and presentations to figshare. [#1993](https://github.com/landlab/landlab/issues/1993)
- Added a notebook gallery that shows all of the *Landlab* tutorial and teaching
  notebooks on the documentation site. [#1994](https://github.com/landlab/landlab/issues/1994)
- Converted documentation files from ReStructuredText to MyST. [#2003](https://github.com/landlab/landlab/issues/2003)


### 🔩 Other Changes and Additions

- Speed up the Flexure component through parallelism using openmp. [#1660](https://github.com/landlab/landlab/issues/1660)
- Provide effective sediment erosion and deposition rates for {class}`~.SpaceLargeScaleEroder`. [#1931](https://github.com/landlab/landlab/issues/1931)
- Require ``numpy<2`` temporarily while we upgrade to support ``numpy>=2``. [#1945](https://github.com/landlab/landlab/issues/1945)
- Updated *Landlab* to be compatible with *numpy* v2. [#1956](https://github.com/landlab/landlab/issues/1956)
- Enhanced the continuous integration workflows to automatically
  cancel any ongoing tests when new changes are pushed to the same branch. [#1971](https://github.com/landlab/landlab/issues/1971)
- Changed the Landlab package to a src-layout. [#1985](https://github.com/landlab/landlab/issues/1985)
- Decreased the run time of the continuous integraion tests by about a factor of two. [#1995](https://github.com/landlab/landlab/issues/1995)
- Moved the nbmake markers stored in the notebooks' metadata into their own
  namespace for newer versions of nbmake. [#1996](https://github.com/landlab/landlab/issues/1996)
- Added a configurable parameter to the `landlab authors` command that
  allows for a user-specified start string for the credits list in the
  authors file. [#2005](https://github.com/landlab/landlab/issues/2005)
- Changed the changelog and towncrier to use myst files. [#2006](https://github.com/landlab/landlab/issues/2006)
- Modified *Landlab*'s test workflow to run the tests using wheels built
  using cibuildwheel. [#2012](https://github.com/landlab/landlab/issues/2012)
- Combined the `release` and prerelease `workflows` with the `test`
  workflow. [#2013](https://github.com/landlab/landlab/issues/2013)

## 2.8.0 (2024-05-12)

### ✨ New Components

- Added new component {class}`~.ConcentrationTrackerForDiffusion`
  for tracking hillslope sediment properties. ([#1662](https://github.com/landlab/landlab/issues/1662))
- Added new component {class}`~.MassWastingRunout`
  for predicting the hazard extent, sediment transport and topographic change associated with the runout of a landslide. ([#1830](https://github.com/landlab/landlab/issues/1830))

### 🍰 New Features

- Fixed the Radiation component by computing fields with ASCE standard formulas, added
  increased test coverage for both field computations and standard unit testing.
  Min, max, and avg daily temp are also three optional, newly added
  user-defined arguments for the component. ([#1755](https://github.com/landlab/landlab/issues/1755))
- Added a new grid type, {class}`~.IcosphereGlobalGrid` (plus underlying graph
  machinery, etc.). ([#1808](https://github.com/landlab/landlab/issues/1808))
- Added a new function, calc_net_face_flux_at_cell, that computes the
  net flux of a quantity into each of a RasterModelGrid's cells. This
  function uses openmp to parallelize its calculations. ([#1900](https://github.com/landlab/landlab/issues/1900))
- Added a new *vtk* writer, `landlab.io.legacy_vtk.dump` that is
  able to write *Landlab* grids that have three spatial coordinates.
  This function is also able to write both the main grid (*nodes* and
  *patches*) as well as the dual grid (*corners* and *cells*). ([#1932](https://github.com/landlab/landlab/issues/1932))

### 🛠️ Bug Fixes

- Fixed a bug when ordering links at patches with patches composed of varying
  numbers of links. ([#1807](https://github.com/landlab/landlab/issues/1807))
- Fixed a bug where SpaceLargeScaleEroder deviates from analytical solution for mixed bedrock-alluvial river in a portion of the parameter space. ([#1901](https://github.com/landlab/landlab/issues/1901))
- Fixed a bug that caused a `ModuleNotFoundError` for *pkg_config* on
  Python 3.12. ([#1927](https://github.com/landlab/landlab/issues/1927))

### 📖 Documentation Enhancements

- Update list of publications in USEDBY.rst. ([#1928](https://github.com/landlab/landlab/issues/1928))

### 🔩 Other Changes and Additions

- Removed the broken `map_link_vector_to_nodes` method from
  `ModelGrid`. As a replacement, use
  {func}`~.map_link_vector_components_to_node_raster` for raster grids, and
  {func}`~.map_link_vector_components_to_node_hex` for hex grids. ([#1786](https://github.com/landlab/landlab/issues/1786))
- Fixed the path to the requirements file needed by *readthedocs*. ([#1797](https://github.com/landlab/landlab/issues/1797))
- Fixed an issue that caused with the CI to fail when building *multidict* on
  Mac and Python 3.12. ([#1850](https://github.com/landlab/landlab/issues/1850))
- Fixed warnings caused by using xarray.Dataset.dims rather than
  xarray.Dataset.sizes. ([#1910](https://github.com/landlab/landlab/issues/1910))
- Added a list of Landlab's extensions to setup.py that must be
  maintained manually. This replaces the old, and somewhat buggy,
  method of conducting a recursive glob for pyx files. ([#1915](https://github.com/landlab/landlab/issues/1915))
- Added a new linter, *cython-lint*, that checks for lint in cython
  files. ([#1924](https://github.com/landlab/landlab/issues/1924))
- Changed the *numpy* printing options from the legacy 1.13 format
  to the latest, default, version. ([#1929](https://github.com/landlab/landlab/issues/1929))
- Removed duplicate shapefile modules. ([#1933](https://github.com/landlab/landlab/issues/1933))

## 2.7.0 (2023-11-04)

### ✨ New Components

- Added new component {class}`~.GravelBedrockEroder` to model rock-cutting gravel rivers. ([#1505](https://github.com/landlab/landlab/issues/1505))
- Added new component {class}`~.AdvectionSolverTVD` for advection using
  a Total Variation Diminishing method. ([#1582](https://github.com/landlab/landlab/issues/1582))

### 📚 New Tutorial Notebooks

- Added a tutorial notebook for the {class}`~.AdvectionSolverTVD` component. ([#1582](https://github.com/landlab/landlab/issues/1582))

### 🍰 New Features

- Added two new mapping functions to assist numerical advection schemes:
  {func}`~.map_node_to_link_linear_upwind` and {func}`~.map_node_to_link_lax_wendroff`. ([#1570](https://github.com/landlab/landlab/issues/1570))
- Added {attr}`.RasterModelGrid.orientation_of_link` and {attr}`.HexModelGrid.orientation_of_link`
  attributes to get orientation codes for links. ([#1573](https://github.com/landlab/landlab/issues/1573))
- Added {attr}`.RasterModelGrid.parallel_links_at_link` and {attr}`.HexModelGrid.parallel_links_at_link`
  attributes. ([#1576](https://github.com/landlab/landlab/issues/1576))
- AdvectionSolverTVD can now handle advection of multiple fields ([#1632](https://github.com/landlab/landlab/issues/1632))
- Refactor ListricKinematicExtender to use AdvectionSolverTVD ([#1635](https://github.com/landlab/landlab/issues/1635))
- Add output function for legacy VTK files ([#1643](https://github.com/landlab/landlab/issues/1643))
- Added an `rng` keyword to the {class}`~.NetworkSedimentTransporter` utilities
  that allows a user to control the random number generator used. ([#1722](https://github.com/landlab/landlab/issues/1722))
- Added an `alpha` keyword to {func}`~.plot.imshow_grid` that allows a user to set
  the transparency value for image plots. ([#1735](https://github.com/landlab/landlab/issues/1735))
- Added the ability for {class}`~.OverlandFlow` to accept an array
  for the `rainfall_intensity` keyword. ([#1743](https://github.com/landlab/landlab/issues/1743))

### 🛠️ Bug Fixes

- Fixed a bug that prevented the {class}`~.DrainageDensity` component from
  working on hex grids. ([#1266](https://github.com/landlab/landlab/issues/1266))
- Fixed a boundary condition issue on D8 flow accumulation in the {class}`~.PriorityFloodFlowRouter`. ([#1542](https://github.com/landlab/landlab/issues/1542))
- Fixed broken link to header image in `notebooks/tutorials/syllabus.ipynb`. ([#1556](https://github.com/landlab/landlab/issues/1556))
- Update obsolete function name in raster_gradients.calc_slope_at_node ([#1606](https://github.com/landlab/landlab/issues/1606))
- Fixed a bug in {class}`~.SpaceLargeScaleEroder` where it would overwrite parts
  of the *sediment\_\_influx* field with zeros. ([#1638](https://github.com/landlab/landlab/issues/1638))
- Fixed a bug where the `colorbar_label` keyword of {func}`~.imshow.imshow_grid`
  was being ignored for non-raster grids. ([#1654](https://github.com/landlab/landlab/issues/1654))
- Fixed  errors introduced with *argsort* from *numpy* v1.25. These were the result of
  vectorized versions of the quicksort algorithm used on some architectures. ([#1670](https://github.com/landlab/landlab/issues/1670))
- Fixed an issue with the agent based modeling tutorial notebooks that
  caused a "too many values to unpack" error with *mesa* v2. ([#1674](https://github.com/landlab/landlab/issues/1674))
- Fixed an issue with {class}`~.PriorityFloodFlowRouter` where flooded nodes were not updated properly.
  This is fixed by setting the `flood_status_code` to 3 (i.e. `_FLOODED`) ([#1683](https://github.com/landlab/landlab/issues/1683))
- Fixed a bug that caused an incorrect Python version to be used in *Landlab*'s
  continuous integration tests. ([#1754](https://github.com/landlab/landlab/issues/1754))

### 📖 Documentation Enhancements

- Added links in docs and README to open Landlab tutorials on EarthscapeHub. ([#1556](https://github.com/landlab/landlab/issues/1556))
- Removed out-dated installation instructions from the documentation. ([#1592](https://github.com/landlab/landlab/issues/1592))
- Add a tutorial notebook on bringing Landlab raster NetCDF output into Paraview for visualization and animation. ([#1646](https://github.com/landlab/landlab/issues/1646))
- Fixed an error that caused the documentation build to fail with an error
  saying that the documentation was not using `furo.css` as the stylesheet. ([#1696](https://github.com/landlab/landlab/issues/1696))
- Add tutorial on bringing Landlab .obj output into Blender ([#1698](https://github.com/landlab/landlab/issues/1698))
- Updated the installation instructions to include options to fetch dependencies
  from, and only from, *conda-forge*. ([#1704](https://github.com/landlab/landlab/issues/1704))
- Reformatted all *doctests* and *reStructuredText* *code-blocks* to conform
  to [black](https://github.com/psf/black), giving the code across all of
  our documentation a consistent format. To keep things formatted correctly,
  added [blackdoc](https://github.com/keewis/blackdoc) to our linters ([#1785](https://github.com/landlab/landlab/issues/1785))

### 🔩 Other Changes and Additions

- Removed the `on_diagonals` method from the {class}`~.LinearDiffuser` component. ([#1236](https://github.com/landlab/landlab/issues/1236))
- Moved unversioned requirements into *requirements.in* files and pinned
  requirements into *requirements.txt* files. ([#1546](https://github.com/landlab/landlab/issues/1546))
- Set up [dependabot](https://docs.github.com/en/code-security/dependabot/dependabot-version-updates/about-dependabot-version-updates)
  to track and update dependencies. ([#1546](https://github.com/landlab/landlab/issues/1546))
- Added pre-commit hooks for delinting the notebooks and removed newly-found
  lint. ([#1585](https://github.com/landlab/landlab/issues/1585))
- Changed the target branch for *dependabot* pull requests to *dependencies*
  and added a GitHub action that automatically keeps the *dependencies* branch
  up-to-date with *master*. ([#1602](https://github.com/landlab/landlab/issues/1602))
- Added two new references to list of publications. ([#1603](https://github.com/landlab/landlab/issues/1603))
- Added better error reporting and input validation for the LinearDiffser
  component. ([#1607](https://github.com/landlab/landlab/issues/1607))
- Added Cython 3.x (beta) to the build-system for compiling extension modules. ([#1639](https://github.com/landlab/landlab/issues/1639))
- Fixed an issue with a missing package, *hypothesis*, not being installed when
  the notebook tests were run through *nox*. ([#1644](https://github.com/landlab/landlab/issues/1644))
- Added getters for several {class}`~.BedrockLandslider` input parameters. ([#1651](https://github.com/landlab/landlab/issues/1651))
- Added getter for several {class}`~SpaceLargeScaleEroder` input parameters. ([#1652](https://github.com/landlab/landlab/issues/1652))
- Modified the *TaylorDiffuser* components, {class}`~.DepthDependentTaylorDiffuser` and
  {class}`~.TaylorNonLinearDiffuser` to use the shortest link instead of `dx` in calculatting
  time steps. ([#1694](https://github.com/landlab/landlab/issues/1694))
- Changed the continuous integraion to use *micromamba* rather than *miniconda*. ([#1703](https://github.com/landlab/landlab/issues/1703))
- Updated *Landlab* for *matplotlib* 3.7.2. Removed calls to newly deprecated
  `get_cmap` and fixed some notebook errors. ([#1714](https://github.com/landlab/landlab/issues/1714))
- Removed unused requirements for building the documentation. ([#1720](https://github.com/landlab/landlab/issues/1720))
- Fixed a flaky test with the {class}`~.lateral_erosion.lateral_erosion.LateralEroder` where it would occasionally
  fail to reach the steady state solution. ([#1722](https://github.com/landlab/landlab/issues/1722))
- Fixed a flaky test with the `sediment_pulser_at_links.ipynb` notebook where it
  would occasionally hang. ([#1722](https://github.com/landlab/landlab/issues/1722))
- Fixed incorrect doctests for `parallel_links_at_link` and
  `orientation_of_link`. ([#1738](https://github.com/landlab/landlab/issues/1738))
- Fixed an issue with *Landlab*'s environment file that caused an error when
  trying to run the tutorial notebooks through *Binder*. ([#1758](https://github.com/landlab/landlab/issues/1758))
- Updated the *readthedocs* configuration file to exclude the
  [now invalid](https://blog.readthedocs.com/drop-support-system-packages)
  `system_packages` option. ([#1762](https://github.com/landlab/landlab/issues/1762))
- Updated the *isort* configuration to identify *landlab* as a first-party
  package to prevent it from moving *landlab* imports into the third-party
  section. ([#1763](https://github.com/landlab/landlab/issues/1763))
- Updated *dependabot* to only manage *Landlab* direct dependencies and changed
  our CI to ensure we are running with those pinned dependencies. ([#1773](https://github.com/landlab/landlab/issues/1773))
- Added support for Python 3.12 and dropped Python 3.9. ([#1782](https://github.com/landlab/landlab/issues/1782))
- Removed the unused and broken *cython* functions `reorient_links` and
  `get_angle_of_links` from the `remap_element` *cython* module. ([#1788](https://github.com/landlab/landlab/issues/1788))
- Fixed flaky tests of the {class}`~.SedimentPulserAtLinks` and
  {class}`~.SedimentPulserEachParcel` components by testing them using a random seed. ([#1794](https://github.com/landlab/landlab/issues/1794))
- Added a tool that builds a list of *Landlab* contributors and updates the
  `AUTHORS.rst` and `.mailmap` files. ([#1795](https://github.com/landlab/landlab/issues/1795))

## 2.6.0 (2023-02-16)

### ✨ New Components

- Added two {class}`SedimentPulser <.SedimentPulserBase>` components ({class}`~.SedimentPulserAtLinks`,
  {class}`~.SedimentPulserEachParcel`) that allow the user to efficiently add sediment
  parcels to the {class}`~.DataRecord` while using the
  {class}`~.NetworkSedimentTransporter` component. ([#1208](https://github.com/landlab/landlab/issues/1208))

- Added a set of {class}`BedParcelInitializer <.BedParcelInitializerBase>` components
  ({class}`~.BedParcelInitializerDischarge`, {class}`~.BedParcelInitializerDepth`,
  {class}`~.BedParcelInitializerArea`, {class}`~.BedParcelInitializerUserD50`) that
  allow the user to efficiently create initial river bed sediment conditions for use
  in the {class}`~.NetworkSedimentTransporter` component. ([#1208](https://github.com/landlab/landlab/issues/1208))

- Added a new component, {class}`~.GravelRiverTransporter`, that models
  gravel transport and abrasion in a gridded network of river segments. ([#1439](https://github.com/landlab/landlab/issues/1439))

- Added a new component, {class}`~.AreaSlopeTransporter`.

  The {class}`~.AreaSlopeTransporter` is a generic transport-limited landscape evolution component that models the time rate of change of elevation at a set of grid nodes, each of which has a defined contributing drainage area 𝐴 (field drainage_area) and a local steepest-descent slope gradient, 𝑆, defined from the node itself toward one of its neighboring nodes. The drainage area and slope can be computed with a drainage-routing component such as {class}`~.FlowAccumulator` or {class}`~.PriorityFloodFlowRouter`. The component is designed to function as an integral part of a transport-limited landscape evolution model in the spirit of the Willgoose et al. "SIBERIA" model. ([#1502](https://github.com/landlab/landlab/issues/1502))

### 📚 New Tutorial Notebooks

- Added tutorial notebooks for the new
  {class}`BedParcelInitializer <.BedParcelInitializerBase>` and
  {class}`SedimentPulser <.SedimentPulserBase>` components, all associated with the
  {class}`~.NetworkSedimentTransporter`. ([#1208](https://github.com/landlab/landlab/issues/1208))
- Added a tutorial notebook that demonstrates use of the new {class}`~.GravelRiverTransporter` component. ([#1439](https://github.com/landlab/landlab/issues/1439))

### 🍰 New Features

- Updated the `NetworkSedimentTransporter` component to allow the user to
  specify a minimum acceptable channel slope threshold. ([#1208](https://github.com/landlab/landlab/issues/1208))
- Added the `calculate_window_statistic` utility that calculates local grid node statistics within a moving window. ([#1263](https://github.com/landlab/landlab/issues/1263))
- Added the `at` keyword to the `imshow_grid` functions so that they now
  use the same pattern as many other *landlab* functions. ([#1424](https://github.com/landlab/landlab/issues/1424))
- The `plot_graph` function now can take lists of graph elements rather than only comma-separated strings. ([#1425](https://github.com/landlab/landlab/issues/1425))
- Added a new keyword, `axes` to `plot_graph` to allow plotting within an
  existing axes. ([#1425](https://github.com/landlab/landlab/issues/1425))
- Enhanced the `plot_graph` function: allow the `with_id` keyword to
  accept a list of elements that should have included IDs, fill in patches and
  cells. ([#1425](https://github.com/landlab/landlab/issues/1425))
- Added an `imshow` method to all *landlab* model grids that is a wrapper for
  the `imshow_grid` function. ([#1430](https://github.com/landlab/landlab/issues/1430))
- Updated the `BedrockLandslider` component so that a user can now specify a
  threshold slope to determine the transport length within the deposition part
  of the component. ([#1431](https://github.com/landlab/landlab/issues/1431))
- Added the `ThresholdEroder` component that erodes material to a user-suplied maximum slope. ([#1440](https://github.com/landlab/landlab/issues/1440))
- Added a new class of grid, *FramedVoronoiGrid* which is an elaborated version of the VoronoiDelaunayGrid. The user input parameters to automatically calculate the positions of the nodes. The boundary nodes are automatically fixed, in a not random way. The core nodes are first positioned in a rectangular pattern, and then moved by a random distance in such a way that a minimal distance between nodes is respected. This minimal distance is convenient when we have to run diffusion or river incision processes on the grid, which can become unstable for two small distances between nodes (depending on the timestep of the run). ([#1450](https://github.com/landlab/landlab/issues/1450))
- Enhance possibilities for .pyx compilation through setup.py update. Now include the tests and compile using Python 3. Compilation instructions (C++, multithreading openmp, macros) can be added at top of the .pyx and .pxd files. See use case with files linked to the future FlowRouter component (including tests). ([#1467](https://github.com/landlab/landlab/issues/1467))
- Enhance Exponential weatherer so that it takes spatially explicit input values for soil production maximum rate and soil production decay depth. ([#1529](https://github.com/landlab/landlab/issues/1529))

### 🛠️ Bug Fixes

- Fixed an issue in the NetworkSedimentTransporter tutorial notebooks related to
  deprecated xarray dataset syntax in the calc_aggregate_value method of `DataRecord` ([#1208](https://github.com/landlab/landlab/issues/1208))
- Fixed a bug in notebooks that use *bmi-topography* where an incorrect API key was being used. ([#1410](https://github.com/landlab/landlab/issues/1410))
- Fixed a bug in `plot_graph` where patch and cell polygons were not drawn. ([#1428](https://github.com/landlab/landlab/issues/1428))
- Fixed a bug where `plot_graph` would incorrectly include the last
  node/corner with patches/cells that had fewer links/faces than the maximum of
  the graph. ([#1428](https://github.com/landlab/landlab/issues/1428))
- Fixed an issue related to flow re-routing on an irregular Voronoi-Delaunay grid. ([#1442](https://github.com/landlab/landlab/issues/1442))
- Fixed the ABM tutorial notebooks that were using an older syntax for the
  *Mesa* *remove_agent* method. ([#1444](https://github.com/landlab/landlab/issues/1444))
- Fixed a bug where the *tests* folder was also being installed in
  *site-packages*. ([#1445](https://github.com/landlab/landlab/issues/1445))
- Fixed a bug where the `SpaceLargeScaleEroder` was only able to accept a scalar value for the erodibility coefficient.
  Now it is able to accept either a scalar or an array. ([#1477](https://github.com/landlab/landlab/issues/1477))
- Fixed a bug in *imshowhs_grid* where, when a no-data drape was provided, the plot was
  inverted in the north-south directions. ([#1484](https://github.com/landlab/landlab/issues/1484))
- Fixed a bug in *imshowhs_grid* where the hillshade base layer was not plotting data from rows and columns adjacent to boundary nodes. ([#1484](https://github.com/landlab/landlab/issues/1484))
- Fixed a bug in the *HyLandsTutorial* notebook where the *BedrockLandslider*'s
  *topographic\_\_elevation* field was not being updated correctly. ([#1490](https://github.com/landlab/landlab/issues/1490))
- Fixed a bug in *imshowhs_grid* that caused the axis tick marks to be slightly in the wrong position. ([#1492](https://github.com/landlab/landlab/issues/1492))
- Fixed a bug in *imshowhs_grid* where boundary nodes were not indicated even if requested. ([#1492](https://github.com/landlab/landlab/issues/1492))
- Fixed an issue when plotting the colorbar in the `plot_drainage` function
  using *matplotlib* 3.6. ([#1493](https://github.com/landlab/landlab/issues/1493))
- Fixed usages of `plt.gca` that used keywords to create new axes objects. With
  *matplotlib* 3.6, the way to do this is with `plt.axes`. ([#1494](https://github.com/landlab/landlab/issues/1494))
- Fixed a possible memory leak caused by using the *lru_cache* decorator of
  class methods. ([#1514](https://github.com/landlab/landlab/issues/1514))
- Fixed a bug that, when using randomly positioned nodes, sometimes resulted in a voronoi
  diagram that contained cells without any vertices. ([#1516](https://github.com/landlab/landlab/issues/1516))

### 📖 Documentation Enhancements

- Combined multiple "Tectonics" sections on the component documentation page. ([#1415](https://github.com/landlab/landlab/issues/1415))
- Fixed the broken links to the openearthscape JupyterHub. ([#1419](https://github.com/landlab/landlab/issues/1419))
- Cleaned up the indexing of field names used and provided by all landlab
  components. ([#1476](https://github.com/landlab/landlab/issues/1476))
- Cleaned up the categorization of all the landlab grid methods. ([#1476](https://github.com/landlab/landlab/issues/1476))
- Updated the installation instructions for the tutorial notebooks to better
  describe how to install the tutorial dependencies. ([#1526](https://github.com/landlab/landlab/issues/1526))
- Added additional publications to the list in the documentation. ([#1538](https://github.com/landlab/landlab/issues/1538))

### 🔩 Other Changes and Additions

- Added a pull request template that contains a checklist of items for
  contributors to complete. ([#1340](https://github.com/landlab/landlab/issues/1340))
- Added a citation file, using the Citation File Format, that describes how to cite the *landlab* code base. ([#1342](https://github.com/landlab/landlab/issues/1342))
- Added a short script that can be used to download a set of *landlab* notebooks
  that are compatible with a specified version of *landlab*. ([#1408](https://github.com/landlab/landlab/issues/1408))
- Moved static project metadata into pyproject.toml. ([#1409](https://github.com/landlab/landlab/issues/1409))
- Fixed an issue where notebooks that download DEMs from OpenTopography were
  failing with an error about a missing API key. ([#1410](https://github.com/landlab/landlab/issues/1410))
- Fixed some failing read_shapefile tests related to a new version of pyshp by requiring pyshp != v2.3.0. ([#1418](https://github.com/landlab/landlab/issues/1418))
- Fixed some typos in the doctest for the `StreamPowerEroder`. ([#1426](https://github.com/landlab/landlab/issues/1426))
- Added  *water_surface\_\_elevation* as a field in the
  `LinearDiffusionOverlandFlowRouter`. ([#1433](https://github.com/landlab/landlab/issues/1433))
- Fixed doctests that were failing because "0"s were being printed as "-0"s. ([#1435](https://github.com/landlab/landlab/issues/1435))
- Added a GitHub Actions workflow to the continuous integration that checks to
  see if a pull request contains a news fragment. ([#1446](https://github.com/landlab/landlab/issues/1446))
- Update tutorial template notebook remove obsolute "%" magic and edit description of link to tutorials page. ([#1457](https://github.com/landlab/landlab/issues/1457))
- Added unit tests for the cython function, *adjust_flow_receivers*, used by the *FlowDirectorSteepest* component. ([#1459](https://github.com/landlab/landlab/issues/1459))
- Added a *nox* file to help with routine project maintenance tasks like, for
  example, running the tests, and checking for coding style. ([#1469](https://github.com/landlab/landlab/issues/1469))
- Added two new *nox* sessions: *requirements* and *nuke*. *requirements*
  recreates the various requirements files while *nuke* does an extra bit of
  cleaning. ([#1474](https://github.com/landlab/landlab/issues/1474))
- Fixed an issue that prevented the docs from building due to a compatibility
  issue with *sphinxcontrib.towncrier* and *towncrier* v22.8. ([#1480](https://github.com/landlab/landlab/issues/1480))
- Changed `FramedVoronoiGrid` to accept a single seed for the *seed* keyword. ([#1495](https://github.com/landlab/landlab/issues/1495))
- Modified to skip the doctests for `ExampleData` and `write_esri_ascii` that created
  files in the user's working directory. These doctests are now repeated as unit tests
  that clean up after themselves. ([#1496](https://github.com/landlab/landlab/issues/1496))
- Improved the error message that's reported when a user attempts to add a field
  to a grid that already contains a field with that name. ([#1500](https://github.com/landlab/landlab/issues/1500))
- Allow the cumulative_subsidence_depth field in ListricKinematicExtender to clobber a pre-existing field, which is needed if the caller has read in a pre-existing saved grid. ([#1510](https://github.com/landlab/landlab/issues/1510))
- Fixed a broken *pre-commit* hook that caused an error when checking for lint
  with *flake8*. ([#1512](https://github.com/landlab/landlab/issues/1512))
- Added *flake8-comprehension* to the *flake8* *pre-commit* hook to identify
  comprehension-related lint. ([#1512](https://github.com/landlab/landlab/issues/1512))
- Added additional linters via pre-commit hooks and removed the newly discovered
  lint. ([#1514](https://github.com/landlab/landlab/issues/1514))
- In the Tutorials doc, updated the URL to download the `notebook.py` script from GitHub. ([#1520](https://github.com/landlab/landlab/issues/1520))
- Updated code to work with *numpy* v1.24 and *scipy* v1.10. ([#1521](https://github.com/landlab/landlab/issues/1521))
- Removed the *richdem* package as a mandatory requirement for *landlab*. Users
  needing to use *richdem* (i.e. the `PriorityFloodFlowRouter`) must now install it
  separately. ([#1523](https://github.com/landlab/landlab/issues/1523))
- Updated the pre-commit hooks (most notably flake8 and its plugins) and removed
  newly-found lint. ([#1524](https://github.com/landlab/landlab/issues/1524))
- Updated *Landlab*'s CI to use Python 3.11 and to drop testing with Python 3.8. ([#1527](https://github.com/landlab/landlab/issues/1527))
- Updated `landlab.__version__` to match that of the latest release. ([#1531](https://github.com/landlab/landlab/issues/1531))
- Removed obsolete files from the top-level directory of the repository. ([#1534](https://github.com/landlab/landlab/issues/1534))
- Updated the ci workflows to use a newer version of cibuildwheel when building
  wheels for releases and pre-releases. ([#1536](https://github.com/landlab/landlab/issues/1536))
- Fixed a test failure in the `PriorityFlood_realDEMs.ipynb` notebook by
  constraining bmi-topography to versions other than 0.8.1. ([#1539](https://github.com/landlab/landlab/issues/1539))
- Changed the ci testing of the notebooks to use nbmake. ([#1541](https://github.com/landlab/landlab/issues/1541))
- Increased the stacklevel for warnings from 1 (the default) to 2 to provide
  more information to the user. ([#1545](https://github.com/landlab/landlab/issues/1545))

## 2.5.0 (2022-04-15)

### ✨ New Components

- `CarbonateProducer` Grow carbonate strata using growth function of Bosscher and Schlager (1992). ([#1284](https://github.com/landlab/landlab/issues/1284))
- `DimensionlessDischarge`, that calculates the dimensionless discharge value, debris flow threshold value, and boolean for predicted debris flow for stream segments. ([#1377](https://github.com/landlab/landlab/issues/1377))
- `LinearDiffusionOverlandFlowRouter`: overland flow using the linearized diffusion-wave approximation. ([#1383](https://github.com/landlab/landlab/issues/1383))

### 📚 New Tutorial Notebooks

- Added a notebook that shows how to use USGS NHDPlus HR datasets with the
  `NetworkSedimentTransporter` component ([#1345](https://github.com/landlab/landlab/issues/1345))
- Addded a new notebook that demonstrates ways to create a `NetworkModelGrid` from a DEM fetched from *OpenTopography* using the *Topography* utility. ([#1400](https://github.com/landlab/landlab/issues/1400))

### 🍰 New Features

- Added the ability for a user to add layers at grid elements other than cells (i.e.
  nodes, links, etc.).  Previously, the *at_layer* variables could only be at cell elements. ([#1292](https://github.com/landlab/landlab/issues/1292))
- Added the ability to define the units of a field when creating a grid from a file
  through the `create_grid` function. ([#1358](https://github.com/landlab/landlab/issues/1358))
- Added the `network_grid_from_raster` function that creates a `NetworkModelGrid`
  from a `RasterModelGrid`. This function extracts channel segments from the
  source grid to become links of the newly-created grid. ([#1360](https://github.com/landlab/landlab/issues/1360))
- Added *sediment\_\_influx* and *sediment\_\_outflux* fields to the `ErosionDeposition`,
  `LateralEroder`, `SpaceLargeScaleEroder`, and `Space` components. ([#1370](https://github.com/landlab/landlab/issues/1370))
- Added `ticks_km`, `cbar_ticks_color` keywords to the `imshowhs_grid` function for more control of colorbar ticks. ([#1397](https://github.com/landlab/landlab/issues/1397))
- Added control on location of the ylabels of colorbars in the `imshowhs_grid` function using the `y_label_offSet_var_1` and `y_label_offSet_var_2` keywords. ([#1397](https://github.com/landlab/landlab/issues/1397))
- Added a new utility, *plot_layers*, that plots sediment layers along with sea level and bedrock. ([#1398](https://github.com/landlab/landlab/issues/1398))

### 🛠️ Bug Fixes

- Clip active layer thickness to zero in the NetworkSedimentTransporter component. This
  eliminates an `invalid value encountered in power` warning. ([#1356](https://github.com/landlab/landlab/issues/1356))
- Allow *landlab* to be installed without the *richdem* package in the case that
  *richdem* is not available for a particular platform or Python version. ([#1379](https://github.com/landlab/landlab/issues/1379))
- Resolved instabilities related to the use of very small `H*` values when using the `Space_Large_Scale_Eroder`. ([#1397](https://github.com/landlab/landlab/issues/1397))
- Fixed a broken reference in the `PriorityFloodFlowDirector` where the gradient of the hillslopes are being updated. ([#1397](https://github.com/landlab/landlab/issues/1397))
- Fixed a bug that incorrectly diagnosed if the *richdem* engine was installed and working correctly. ([#1399](https://github.com/landlab/landlab/issues/1399))

### 📖 Documentation Enhancements

- Added missing documentation files for `BedrockLandslider` and `SpaceLargeScaleEroder`. ([#1373](https://github.com/landlab/landlab/issues/1373))
- Set up *\[towncrier\](https://towncrier.readthedocs.io/)*
  to update and manage the *landlab* changelog. New fragments are placed in the
  `news/` folder. ([#1396](https://github.com/landlab/landlab/issues/1396))

### 🔩 Other Changes and Additions

- Added an OpenTopography API key to notebooks that use *bmi-topography* to fetch
  data from OpenTopography. ([#1384](https://github.com/landlab/landlab/issues/1384))
- Updated the coding style to conform to new version of black. This was, primarily,
  hugging the `**` operator. ([#1385](https://github.com/landlab/landlab/issues/1385))
- The notebooks are tested only with Python 3.9. ([#1399](https://github.com/landlab/landlab/issues/1399))
- Added Python 3.10 to continuous integration tests and dropped Python 3.7. ([#1399](https://github.com/landlab/landlab/issues/1399))
- Speed up our continuous integration tests by about 2x by running them in parallel using *pytest-xdist*. ([#1399](https://github.com/landlab/landlab/issues/1399))
- Turn off *hypothesis* deadline setting globally when running continuous
  integration tests. ([#1401](https://github.com/landlab/landlab/issues/1401))
- Updated the documentation to build with newer versions of *Sphinx*. ([#1404](https://github.com/landlab/landlab/issues/1404))
- Added several new *landlab*-using references. ([#1407](https://github.com/landlab/landlab/issues/1407))

## 2.4.1 (2021-12-02)

### 📚 New Tutorial Notebooks

- Added two ABM tutorial notebooks ([#1364](https://github.com/landlab/landlab/issues/1364))

### 🔩 Other Changes and Additions

- fixed a bug that causes release workflows to not be triggered ([#1371](https://github.com/landlab/landlab/issues/1371))
- Fixed the building of source distributions for prerelease and release
  workflows ([#1372](https://github.com/landlab/landlab/issues/1372))

## 2.4.0 (2021-11-29)

### 🔩 Other Changes and Additions

- Changed GitHub actions to use cibuildwheel for building wheels ([#1368](https://github.com/landlab/landlab/issues/1368))

## 2.4.0b0 (2021-11-28)

### ✨ New Components

- ListricKinematicExtender: Simulate Extensional Tectonic Motion on a Listric Fault Plane ([#1283](https://github.com/landlab/landlab/issues/1283))
- PriorityFloodFlowRouter and SpaceLargeScaleEroder ([#1352](https://github.com/landlab/landlab/issues/1352))
- Added BedrockLandslider component ([#1362](https://github.com/landlab/landlab/issues/1362))

### 📚 New Tutorial Notebooks

- Added tutorial notebook for depth dependent taylor diffuser ([#1306](https://github.com/landlab/landlab/issues/1306))

- Added tutorial notebook for chi finder ([#1307](https://github.com/landlab/landlab/issues/1307))

- Added tutorial notebook for kinwave impl ([#1308](https://github.com/landlab/landlab/issues/1308))

- Added tutorial notebook for taylor diffuser ([#1309](https://github.com/landlab/landlab/issues/1309))

- Added notebook tutorials for two components (both written by Jordan Adams):
  \* `DepthSlopeProductErosion`,
  \* `DetachmentLtdErosion`

  Added a tutorial showing how to "D4 pit fill" a DEM, and a version of the simple `hugo_site.asc` DEM that has been pit-filled. ([#1313](https://github.com/landlab/landlab/issues/1313))

- Added tutorial notebook for Space component ([#1314](https://github.com/landlab/landlab/issues/1314))

- Added tutorial notebook for erosiondeposition ttl ([#1315](https://github.com/landlab/landlab/issues/1315))

- Added tutorial notebook for erodep ([#1317](https://github.com/landlab/landlab/issues/1317))

- Added tutorial notebook for StreamPowerSmoothThresholdEroder ([#1331](https://github.com/landlab/landlab/issues/1331))

### 🍰 New Features

- Infer data types of fields when reading from shape files ([#1357](https://github.com/landlab/landlab/issues/1357))

### 🛠️ Bug Fixes

- Fixed ability to pass a masked array to imshow_grid_at_node ([#1297](https://github.com/landlab/landlab/issues/1297))
- Fixed xarray 'axis' keyword error in map function ([#1300](https://github.com/landlab/landlab/issues/1300))
- Fixed a missing absolute value in Courant condition in dupuit_percolator ([#1311](https://github.com/landlab/landlab/issues/1311))
- Fixed pits and division by zero in lateral_erosion component ([#1353](https://github.com/landlab/landlab/issues/1353))

### 📖 Documentation Enhancements

- Updated installation instructions ([#1287](https://github.com/landlab/landlab/issues/1287))
- Minor updates to documentation ([#1290](https://github.com/landlab/landlab/issues/1290))
- Run the link checker on docs ([#1336](https://github.com/landlab/landlab/issues/1336))
- Fixed documentation errors in green ampt component ([#1343](https://github.com/landlab/landlab/issues/1343))
- Added new references to landlab ([#1344](https://github.com/landlab/landlab/issues/1344))
- Added a link to launch landlab notebooks on the CSDMS JupyterHub ([#1347](https://github.com/landlab/landlab/issues/1347))

### 🔩 Other Changes and Additions

- Fixed warnings related to unnecessary use of numpy number types ([#1323](https://github.com/landlab/landlab/issues/1323))
- Changed continuous integration to always check the docs build ([#1336](https://github.com/landlab/landlab/issues/1336))
- Added a pre-commit configuration file ([#1338](https://github.com/landlab/landlab/issues/1338))
- Drop the "file:" prefix when referencing pip requirements files ([#1339](https://github.com/landlab/landlab/issues/1339))
- Removed usages of np.int from Cython code ([#1354](https://github.com/landlab/landlab/issues/1354))
- Check that notebooks are both clean and blackened as part of continuous integration ([#1355](https://github.com/landlab/landlab/issues/1355))

## 2.3.0 (2021-03-19)

### ✨ New Components

- Added a tidal flow component ([#1225](https://github.com/landlab/landlab/issues/1225))
- Added ExponentialWeathererIntegrated component ([#1254](https://github.com/landlab/landlab/issues/1254))
- Added simple submarine diffuser component ([#1269](https://github.com/landlab/landlab/issues/1269))

### 📚 New Tutorial Notebooks

- Added tutorial for river input to LEMs ([#1258](https://github.com/landlab/landlab/issues/1258))

### 🍰 New Features

- Added recharge to the GroundwaterDupuitPercolator callback ([#1223](https://github.com/landlab/landlab/issues/1223))
- Added Wavefront OBJ output ([#1241](https://github.com/landlab/landlab/issues/1241))

### 🛠️ Bug Fixes

- Fixed bug in Flow router/depression finder which incorrectly directed flow ([#1248](https://github.com/landlab/landlab/issues/1248))
- Fixed an error in the streampower notebook ([#1260](https://github.com/landlab/landlab/issues/1260))
- Fixed a bug in the FlowAccumulator to update pit present logic to also include node flood status ([#1277](https://github.com/landlab/landlab/issues/1277))
- Fixed a bug when adding a missing at_grid field when testing components ([#1286](https://github.com/landlab/landlab/issues/1286))

### 📖 Documentation Enhancements

- Fixed documentation bugs ([#1233](https://github.com/landlab/landlab/issues/1233))
- Added two 2020 publications ([#1243](https://github.com/landlab/landlab/issues/1243))
- Added docs for the flow accumulator ([#1251](https://github.com/landlab/landlab/issues/1251))
- Added a reference to the papers and presentations list ([#1255](https://github.com/landlab/landlab/issues/1255))
- Added additional references for 2020 and 2021 ([#1273](https://github.com/landlab/landlab/issues/1273))

### 🔩 Other Changes and Additions

- NetworkSedimentTtransporter JOSS paper fixes ([#1235](https://github.com/landlab/landlab/issues/1235))
- Small changes to JOSS paper ([#1237](https://github.com/landlab/landlab/issues/1237))
- Changed to use GitHub Actions for CI ([#1270](https://github.com/landlab/landlab/issues/1270))
- Added building and testing of landlab with Python 3.9 ([#1274](https://github.com/landlab/landlab/issues/1274))
- Added release and prerelease github actions ([#1275](https://github.com/landlab/landlab/issues/1275))
- Cleaned up landlab metadata files; Removed versioneer, we'll use zest.releaser from now on the manage versions ([#1285](https://github.com/landlab/landlab/issues/1285))

## 1.5.1 (2018-06-22)

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]

## 1.5.0 (2018-06-18)

(fixed-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]

(changed-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]

## 1.4.0 (2018-05-03)

(fixed-2-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]

(added-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]

(changed-2-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]

## 1.3.1 (2018-03-24)

(fixed-3-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]

(added-2-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]

(changed-3-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- Set versioneer to ignore `v` prefix in tags \[Eric Hutton\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]

## 1.3.0 (2018-03-14)

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-4-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]

(added-3-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-4-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]

## 1.2.0 (2017-10-19)

(removed-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-5-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]

(added-4-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-5-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]

## 1.1.0 (2017-06-26)

(removed-2-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-6-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-5-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-6-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]

## 1.0.3 (2017-03-04)

(removed-3-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-7-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-6-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-7-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.2 (2016-11-24)

(removed-4-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-8-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-7-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-8-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.1 (2016-08-25)

(removed-5-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-9-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-8-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-9-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/patch-flowaccum-reclimit \[#375\] \[Dan Hobley\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]

## 1.0.0-beta.8 (2016-07-07)

(removed-6-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-10-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-9-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-10-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.7 (2016-07-07)

(removed-7-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-11-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-10-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-11-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.6 (2016-06-16)

(removed-8-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-12-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-11-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-12-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.5 (2016-06-14)

(removed-9-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-13-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-12-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-13-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.4 (2016-06-14)

(removed-10-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-14-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-13-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-14-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/pot-fr-modernise \[#344\] \[Jordan Adams\]
- landlab/mcflugen/fix-landlab-test-function \[#345\] \[Dan Hobley\]

## 1.0.0-beta.3 (2016-06-06)

(removed-11-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-15-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-14-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-15-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.2 (2016-06-06)

(removed-12-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-16-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-15-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-16-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0-beta.12 (2016-08-02)

(removed-13-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-17-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-16-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-17-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/patch-flowaccum-reclimit \[#375\] \[Dan Hobley\]

## 1.0.0-beta.11 (2016-07-19)

(removed-14-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-18-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-17-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-18-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]

## 1.0.0-beta.10 (2016-07-14)

(removed-15-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-19-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-18-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-19-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]

## 1.0.0-beta.1 (2016-05-18)

(removed-16-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-20-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-19-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-20-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]

## 1.0.0-beta (2016-05-13)

(removed-17-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-21-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-20-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-21-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]

## 1.0.0-alpha (2016-05-13)

(removed-18-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-22-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-21-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-22-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]

## 1.0.0 (2016-08-19)

(removed-19-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-23-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-22-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-23-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-flowaccum-reclimit \[#375\] \[Dan Hobley\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]
- Merge remote-tracking branch ‘refs/remotes/origin/master’ into
  release \[saisiddu\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]

## 0.5.0 (2016-03-12)

(removed-20-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-24-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-23-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-24-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.2.3 (2016-03-10)

(removed-21-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-25-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-24-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-25-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.2.2 (2016-02-09)

(removed-22-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-26-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-25-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-26-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.2.1 (2016-02-08)

(removed-23-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-27-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-26-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-27-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.2.0 (2015-12-20)

(removed-24-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-28-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-27-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-28-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.1.41 (2015-12-13)

(removed-25-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-29-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-28-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-29-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]

## 0.1.40 (2015-12-13)

(removed-26-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-30-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-29-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-30-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]

## 0.1.39 (2015-10-28)

(removed-27-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-31-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-30-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-31-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/SiccarPoint/D4routing \[#176\] \[Jordan Adams\]
- landlab/sed_fill \[#173\] \[Eric Hutton\]
- landlab/sed_fill \[#168\] \[Jordan Adams\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]

## 0.1.38 (2015-10-21)

(removed-28-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-32-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-31-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-32-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]

## 0.1.37 (2015-10-20)

(removed-29-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-33-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-32-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-33-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]

## 0.1.36 (2015-10-20)

(removed-30-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-34-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-33-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-34-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]

## 0.1.35 (2015-10-20)

(removed-31-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-35-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-34-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-35-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]

## 0.1.34 (2015-10-19)

(removed-32-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-36-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-35-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-36-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]

## 0.1.33 (2015-10-11)

(removed-33-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-37-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-36-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-37-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]

## 0.1.32 (2015-10-11)

(removed-34-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-38-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-37-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-38-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]

## 0.1.31 (2015-10-11)

(removed-35-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-39-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-38-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-39-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]

## 0.1.30 (2015-10-11)

(removed-36-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-40-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-39-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-40-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]

## 0.1.29 (2015-09-13)

(removed-37-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-41-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-40-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-41-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]

## 0.1.28 (2015-09-09)

(removed-38-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-42-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-41-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-42-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]
- landlab/mcflugen/update-mappers \[#142\] \[Eric Hutton\]
- landlab/fix-issue-128 \[#129\] \[Eric Hutton\]
- landlab/mcflugen/add-netcdf-tests \[#138\] \[Eric Hutton\]
- landlab/mcflugen/add-better-testing-docs \[#137\] \[Eric Hutton\]

## 0.1.27 (2015-07-29)

(removed-39-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-43-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-42-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-43-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]
- landlab/mcflugen/update-mappers \[#142\] \[Eric Hutton\]
- landlab/fix-issue-128 \[#129\] \[Eric Hutton\]
- landlab/mcflugen/add-netcdf-tests \[#138\] \[Eric Hutton\]
- landlab/mcflugen/add-better-testing-docs \[#137\] \[Eric Hutton\]

## 0.1.26 (2015-07-20)

(removed-40-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-44-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-43-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-44-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]
- landlab/mcflugen/update-mappers \[#142\] \[Eric Hutton\]
- landlab/fix-issue-128 \[#129\] \[Eric Hutton\]
- landlab/mcflugen/add-netcdf-tests \[#138\] \[Eric Hutton\]
- landlab/mcflugen/add-better-testing-docs \[#137\] \[Eric Hutton\]

## 0.1.25 (2015-07-16)

(removed-41-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-45-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-44-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-45-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]
- landlab/mcflugen/update-mappers \[#142\] \[Eric Hutton\]
- landlab/fix-issue-128 \[#129\] \[Eric Hutton\]
- landlab/mcflugen/add-netcdf-tests \[#138\] \[Eric Hutton\]
- landlab/mcflugen/add-better-testing-docs \[#137\] \[Eric Hutton\]
- resolved conflicts in lake mapper \[gregtucker\]
- ierge branch ‘master’ of <https://github.com/landlab/landlab>
  \[gregtucker\]

## 0.1.24 (2015-07-12)

(removed-42-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-46-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-45-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-46-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
- landlab/mcflugen/fix-component-imports \[#214\] \[Dan Hobley\]
- landlab/SiccarPoint/component-tests \[#209\] \[Eric Hutton\]
- landlab/SiccarPoint/component-tests \[#204\] \[Eric Hutton\]
- landlab/mcflugen/add-value-name-decorator \[#199\] \[Dan Hobley\]
- landlab/SiccarPoint/component-introspection \[#200\] \[Jenny Knuth\]
- landlab/mcflugen/fix-voronoi-cell-areas \[#202\] \[Dan Hobley\]
- landlab/mcflugen/fix-testing-script \[#198\] \[Eric Hutton\]
- landlab/mcflugen/add-neighbor-at-node \[#195\] \[saisiddu\]
- landlab/mcflugen/fix-appveyor-builds \[#185\] \[Eric Hutton\]
- landlab/mcflugen/clean-up-imshow \[#178\] \[saisiddu\]
- landlab/mcflugen/fix-reset_lists_of_nodes_cells \[#182\] \[Greg Tucker\]
- landlab/SiccarPoint/var-doc \[#187\] \[Eric Hutton\]
- landlab/gtucker/cleanup \[#186\] \[Eric Hutton\]
- landlab/mcflugen/add-gradients-module \[#169\] \[Eric Hutton\]
- landlab/SiccarPoint/delete-fields \[#167\] \[Dan Hobley\]
- landlab/voronoi_stream_power \[#158\] \[Dan Hobley\]
- landlab/mcflugen/make-node-status-private \[#152\] \[Eric Hutton\]
- landlab/mcflugen/code-clean-up \[#148\] \[Eric Hutton\]
- landlab/mcflugen/add-appveyor-slack-notifications \[#149\] \[Eric
  Hutton\]
- Merge remote-tracking branch ‘origin/link_status’ \[mcflugen\]
- landlab/link_status \[#143\] \[Eric Hutton\]
- landlab/mcflugen/fix-failing-docs-build \[#145\] \[Eric Hutton\]
- landlab/mcflugen/update-mappers \[#142\] \[Eric Hutton\]
- landlab/fix-issue-128 \[#129\] \[Eric Hutton\]
- landlab/mcflugen/add-netcdf-tests \[#138\] \[Eric Hutton\]
- landlab/mcflugen/add-better-testing-docs \[#137\] \[Eric Hutton\]
- resolved conflicts in lake mapper \[gregtucker\]
- ierge branch ‘master’ of <https://github.com/landlab/landlab>
  \[gregtucker\]
- landlab/mcflugen/add-python-3-support-with-six \[#121\] \[Eric Hutton\]
- landlab/mcflugen/fix-setuptools-error \[#126\] \[Eric Hutton\]

## v0.1.23 (2015-05-15)

(removed-43-1)=

### Removed

- Removed inlink and outlink matrices. \[Eric Hutton\]
- Removed deprecated raster_steepest_descent module. \[Eric Hutton\]
- Removed corner_node_at_cell \[Eric Hutton\]
- Removed old and unused \_route_flow_old from lake_mapper \[Eric
  Hutton\]
- Removed unused code from flow_direction_DN \[Eric Hutton\]

(fixed-47-1)=

### Fixed

- Fixed bug in Flexure1D when using “flexure” method \[Eric Hutton\]
- Fixed unit test failures related to masked arrays (#710) \[Eric
  Hutton\]
- Fixed failed Travis builds being reported as passing \[Eric Hutton\]
- Fixed doctest for graph.adjacent_nodes_at_node \[Eric Hutton\]
- Fixed names of packages deployed to Anaconda Cloud \[Eric Hutton\]
- Fixed incorrect signatures of decorated methods in docs. \[Eric
  Hutton\]
- Fixed Travis build errors with Python version conflicts. \[Eric
  Hutton\]
- Fixed values not being cached (#614) \[Eric Hutton\]
- Fixed component documentation not building (issue #575) \[Eric Hutton\]
- Fixed netcdf4 import error \[Eric Hutton\]

(added-46-1)=

### Added

- Added CONTRIBUTING.md document \[Eric Hutton\]
- Added script to create a nicely formatted changelog \[Eric Hutton\]
- Added 1D Flexure component \[Eric Hutton\]
- Added cite_as function to get landlab component citations \[Eric
  Hutton\]
- Added adjacent_nodes_at_node, adjacent_corners_at_corner to Graph.
  \[Eric Hutton\]
- Added additional tests for SoilInfiltrationGreenAmpt. \[Eric Hutton\]
- Added citation tracker for components. \[Eric Hutton\]
- Added nodes_at_link attribute to ModelGrid. \[Eric Hutton\]
- Added event layers to track stratigraphy \[Eric Hutton\]

(changed-47-1)=

### Changed

- amanaster2/master \[#733\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/major_cleanup_to_space_and_erodepo_init \[#709\] \[Katy
  Barnhart\]
- landlab/gt/fix-doctest-issue-726 \[#728\] \[Greg Tucker\]
- landlab/gt/ca-top-hit-bug \[#720\] \[Greg Tucker\]
- landlab/barnhark/space_cell_area \[#719\] \[Greg Tucker\]
- landlab/barnhark/use_field_name_array_or_float \[#683\] \[Katy Barnhart\]
- landlab/barnhark/give_hex_models_watershed_methods \[#685\] \[Katy
  Barnhart\]
- landlab/SiccarPoint/fix-issue-702 \[#706\] \[Katy Barnhart\]
- Giuseppecipolla95/Giuseppecipolla95/make_stream_length_utility \[#658\]
  \[Katy Barnhart\]
- landlab/barnhark/revert_channel_profiler \[#695\] \[Katy Barnhart\]
- landlab/barnhark/space_rounding \[#698\] \[Katy Barnhart\]
- landlab/barnhark/add_docs_to_normal_fault \[#677\] \[Katy Barnhart\]
- landlab/barnhark/space_type_updates \[#669\] \[Katy Barnhart\]
- landlab/barnhark/minor_changes_to_normal_fault \[#663\] \[Katy Barnhart\]
- landlab/gt-debug-ca-propswap \[#661\] \[Greg Tucker\]
- landlab/barnhark/space_hex \[#655\] \[Katy Barnhart\]
- landlab/barnhark/add_kwargs \[#645\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault_kwargs \[#649\] \[Katy Barnhart\]
- landlab/barnhark/normal_fault \[#640\] \[Katy Barnhart\]
- landlab/barnhark/exponential_weatherer_docstring \[#643\] \[Katy
  Barnhart\]
- landlab/nathanlyons/watershed \[#545\] \[Nathan Lyons\]
- landlab/barnhark/updates_to_channel_profile \[#637\] \[Katy Barnhart\]
- landlab/barnhark/typo_in_imshow \[#636\] \[Katy Barnhart\]
- landlab/barnhark/add_component_docs \[#634\] \[Katy Barnhart\]
- landlab/gt-ca-uplift \[#581\] \[Greg Tucker\]
- landlab/barnhark/make_stream_profiler \[#605\] \[Katy Barnhart\]
- landlab/mcflugen/remove-old-flux-div \[#619\] \[Dan Hobley\]
- Simplified continuous integration and versioning. \[Eric Hutton\]
- landlab/barnhark/improving_flow_accumulator_lake_mapper_interactions
  \[#612\] \[Katy Barnhart\]
- landlab/barnhark/fix_stream_power_type_check \[#610\] \[Katy Barnhart\]
- Clean up API for diagonals. \[Eric Hutton\]
- landlab/gt-taylor-fix \[#606\] \[Katy Barnhart\]
- landlab/mcflugen/fix-travis-ioerror \[#607\] \[Nathan Lyons\]
- landlab/barnhark/depth_dependent_boundary_conditions \[#601\] \[Katy
  Barnhart\]
- landlab/mcflugen/tidy-green-ampt \[#591\] \[Jordan Adams\]
- landlab/barnhark/improving_cubic_flux \[#582\] \[Katy Barnhart\]
- Clean up Sphinx documentation \[Eric Hutton\]
- landlab/margauxmouchene/test2 \[#546\] \[margauxmouchene\]
- landlab/gt-fastscape-q \[#574\] \[Greg Tucker\]
- amanaster2/master \[#572\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/kwargs_depth_dependent_diffuser \[#553\] \[Katy
  Barnhart\]
- landlab/gt-lattice-uplifter \[#539\] \[Greg Tucker\]
- landlab/gt-add-phi-to-space-adaptive \[#551\] \[Greg Tucker\]
- landlab/barnhark/cubic_nl_diffuser_kwargs \[#550\] \[Katy Barnhart\]
- landlab/barnhark/no_kwargs_in_dd_cubic_diffuser \[#548\] \[Katy
  Barnhart\]
- landlab/gt-cmap-in-hexplot \[#544\] \[Greg Tucker\]
- landlab/SiccarPoint/uniform_precip \[#517\] \[Dan Hobley\]
- landlab/mcflugen/fix-greenampt-issue-530 \[#535\] \[Katy Barnhart\]
- landlab/mcflugen/add-logging-function \[#504\] \[Eric Hutton\]
- landlab/gt-try-dyn-ts-space \[#529\] \[Katy Barnhart\]
- landlab/barnhark/get_set_state_methods_for_grid \[#525\] \[Greg Tucker\]
- landlab/fixing_small_bug_in_erosion_deposition \[#528\] \[Greg Tucker\]
- landlab/barnhark/eroder_depo_with_n_less_than_one \[#523\] \[Greg
  Tucker\]
- landlab/barnhark/cubic_timestepper \[#519\] \[Greg Tucker\]
- landlab/barnhark/addressing_brent_method_index_error \[#510\] \[Katy
  Barnhart\]
- landlab/gt-edit-erodep \[#516\] \[Katy Barnhart\]
- cmshobe/cmshobe/make-erosion-deposition-component \[#511\] \[Greg
  Tucker\]
- landlab/barnhark/lake_mapper_faster \[#512\] \[Greg Tucker\]
- nathanlyons/master \[#505\] \[Nicole M Gasparini\]
- cmshobe/cmshobe/minor_fixes_to_space \[#509\] \[Katy Barnhart\]
- cmshobe/cmshobe/change-hybrid-to-SPACE \[#506\] \[Katy Barnhart\]
- cmshobe/cmshobe/fix-hybrid-q-mechanics \[#502\] \[Katy Barnhart\]
- RondaStrauch/master \[#480\] \[Sai Siddhartha Nudurupati\]
- landlab/barnhark/use_newton_fastscape \[#492\] \[Katy Barnhart\]
- landlab/barnhark/improve_streampower_smooth_thresh_stability \[#499\]
  \[Greg Tucker\]
- landlab/barnhark/dynamic_timestep_cubic_flux_diffuser \[#497\] \[Greg
  Tucker\]
- landlab/barnhark/switching_mfd_and_dinf_from_slope_to_gradient \[#490\]
  \[Katy Barnhart\]
- landlab/barnhark/cython_hybrid_alluviaum \[#494\] \[Greg Tucker\]
- cmshobe/cmshobe/fix_hybrid_q_options \[#488\] \[Katy Barnhart\]
- landlab/barnhark/smallchangestohybrid \[#487\] \[Greg Tucker\]
- landlab/gt-add-stretched-expo \[#485\] \[Katy Barnhart\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#481\] \[Katy Barnhart\]
- landlab/mcflugen/add-graph-class \[#477\] \[Greg Tucker\]
- landlab/barnhark/accumulator_efficiency \[#476\] \[Greg Tucker\]
- landlab/barnhark/making_flow_accumulator_faster \[#474\] \[Greg Tucker\]
- landlab/barnhark/fixing_kinwave_flow_issue \[#471\] \[Greg Tucker\]
- cmshobe/cmshobe_fixes_to_hybrid_alluv \[#469\] \[Greg Tucker\]
- landlab/gt-implicit-kinwave \[#464\] \[Greg Tucker\]
- cmshobe/cmshobe/make_hybrid_alluv_initis \[#467\] \[Katy Barnhart\]
- Glader011235/master \[#465\] \[Katy Barnhart\]
- landlab/nicgaspar/diffusion_not_depositing \[#463\] \[Jordan Adams\]
- landlab/kbarnhart/make_raster_netcdf \[#462\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#461\] \[Katy Barnhart\]
- cmshobe/cmshobe_hybrid_alluvium_model \[#460\] \[Greg Tucker\]
- Merge remote-tracking branch ‘origin/master’ \[SiccarPoint\]
- Merge remote-tracking branch ‘origin/SiccarPoint/pot-fr’
  \[SiccarPoint\]
- landlab/kbarnhart/consistent_parameter_names \[#459\] \[Katy Barnhart\]
- landlab/gt-stream-power-K \[#457\] \[Greg Tucker\]
- landlab/gt-fix-fastscape-variable-k \[#456\] \[Katy Barnhart\]
- landlab/gt-create-depth-dep-cubic-diffuser \[#452\] \[Katy Barnhart\]
- landlab/mcflugen/add-py36-builds \[#453\] \[Eric Hutton\]
- landlab/kbarnhart/stream_power_error \[#450\] \[Greg Tucker\]
- landlab/gt-fix-issue-448 \[#449\] \[Dan Hobley\]
- landlab/mcflugen/fix-issue-428 \[#447\] \[Jordan Adams\]
- landlab/jadams15/depth_slope_product \[#445\] \[Jordan Adams\]
- landlab/SiccarPoint/fix_429 \[#430\] \[Katy Barnhart\]
- landlab/SiccarPoint/add-docs \[#442\] \[Katy Barnhart\]
- landlab/gt-fix-issue-431 \[#433\] \[Dan Hobley\]
- landlab/gt-add-Q-stream-power-smooth-thresh \[#443\] \[Katy Barnhart\]
- landlab/SiccarPoint/auto-build-docs \[#437\] \[Dan Hobley\]
- landlab/jadams15/spatially_variable_roughness \[#438\] \[Jordan Adams\]
- landlab/kbarnhart/make_nd_fields \[#434\] \[Greg Tucker\]
- landlab/kbarnhart/improvements_to_set_watershed_boundary \[#426\] \[Katy
  Barnhart\]
- landlab/gt-float64-fastscape \[#427\] \[Greg Tucker\]
- landlab/gt-more-cts-cython \[#378\] \[Greg Tucker\]
- landlab/gt-smooth-threshold-stream-power \[#418\] \[Greg Tucker\]
- landlab/gt-tweak-cubic-diff \[#416\] \[Greg Tucker\]
- landlab/gt-fix-init_typo \[#415\] \[Greg Tucker\]
- landlab/jk-move-old-rst \[#412\] \[Greg Tucker\]
- landlab/gt-merge-rg-cubic \[#414\] \[Greg Tucker\]
- cmshobe/cmshobe-drainage-density \[#398\] \[Katy Barnhart\]
- fix minor conflict in raster.py \[Greg Tucker\]
- landlab/SiccarPoint/grid_docs \[#329\] \[Dan Hobley\]
- landlab/SiccarPoint/diagonal_link_lengths \[#328\] \[Eric Hutton\]
- landlab/mcflugen/remove-deprecations \[#327\] \[Eric Hutton\]
- landlab/SiccarPoint/imshow_grid_returns_im \[#326\] \[Dan Hobley\]
- landlab/SiccarPoint/last-minute-deprecation \[#324\] \[Eric Hutton\]
- landlab/SiccarPoint/BAD-INDEX-is-minus1 \[#323\] \[Eric Hutton\]
- landlab/SiccarPoint/patch-methods \[#322\] \[Eric Hutton\]
- landlab/SiccarPoint/tweak-plotter \[#321\] \[Eric Hutton\]
- landlab/saisiddu/Version_1_final \[#320\] \[Eric Hutton\]
- landlab/SiccarPoint/modernise-field-names \[#319\] \[Dan Hobley\]
- landlab/SiccarPoint/modernise-components \[#314\] \[Eric Hutton\]
- landlab/SiccarPoint/most-egregious-diagonals \[#315\] \[Dan Hobley\]
- landlab/gt-calc-of-to-at \[#316\] \[Greg Tucker\]
- landlab/saisiddu/Version_1_final \[#317\] \[Eric Hutton\]
- landlab/jadams15/uniform_precip_changes \[#310\] \[Dan Hobley\]
- landlab/saisiddu/Version_1 \[#311\] \[Dan Hobley\]
- landlab/mcflugen/moved-slope-methods \[#313\] \[Dan Hobley\]
- landlab/SiccarPoint/Horn-slope \[#309\] \[Eric Hutton\]
- landlab/mcflugen/remove-craters \[#312\] \[Eric Hutton\]
- landlab/mcflugen/fix-docs-not-building \[#304\] \[Dan Hobley\]
- landlab/SiccarPoint/grid_trawl \[#307\] \[Eric Hutton\]
- landlab/nicgaspar/watershed_boundary_condition \[#306\] \[Jordan Adams\]
- landlab/SiccarPoint/slopes \[#305\] \[Dan Hobley\]
- landlab/gt-fix-diffuser-bug \[#294\] \[Dan Hobley\]
- landlab/gt-update-gradients \[#303\] \[Greg Tucker\]
- landlab/doc-component-reorg \[#296\] \[Greg Tucker\]
- landlab/gt-fix-ca-tectonics \[#297\] \[Greg Tucker\]
- landlab/gt-flux-divergence \[#295\] \[Greg Tucker\]
- landlab/jk_cleanup_grid_docs \[#289\] \[Greg Tucker\]
- landlab/SiccarPoint/fastscape-threshold \[#290\] \[Jordan Adams\]
- landlab/SiccarPoint/component-modernisation \[#288\] \[Greg Tucker\]
- landlab/gt_fix_faces_at_cell \[#282\] \[Greg Tucker\]
- landlab/sed-flux-dep \[#277\] \[Dan Hobley\]
- landlab/SiccarPoint/chi \[#273\] \[Greg Tucker\]
- landlab/SiccarPoint/plotter_modernisation \[#274\] \[Greg Tucker\]
- landlab/jk_rearrange_index \[#275\] \[Greg Tucker\]
- landlab/SiccarPoint/steepness-index \[#271\] \[nicgaspar\]
- landlab/mcflugen/fix-issue-268 \[#269\] \[Dan Hobley\]
- landlab/mcflugen/add-py35-support \[#270\] \[saisiddu\]
- landlab/SiccarPoint/fix-issue-250 \[#251\] \[Dan Hobley\]
- landlab/SiccarPoint/stream_power_standard \[#256\] \[Eric Hutton\]
- landlab/mcflugen/fix-travis-not-running-all-tests \[#265\] \[Eric
  Hutton\]
- landlab/SiccarPoint/dynamic-docstring-import \[#258\] \[Greg Tucker\]
- landlab/SiccarPoint/enhance-mappers \[#263\] \[Dan Hobley\]
- landlab/SiccarPoint/enhance-mappers \[#262\] \[Dan Hobley\]
- Merged fix for deployment from AppVeyor to PyPI. \[mcflugen\]
- landlab/SiccarPoint/enhance-mappers \[#255\] \[Greg Tucker\]
- landlab/jk_copy_init_docstring \[#248\] \[Jenny Knuth\]
- landlab/gtucker/node_link_connectivity \[#253\] \[Dan Hobley\]
- landlab/mcflugen/add-bmi-bridge \[#246\] \[Greg Tucker\]
- landlab/gt-handle-flooded-nodes-in-stream-power \[#247\] \[Dan Hobley\]
- landlab/jk_include_init_docstring \[#244\] \[Jenny Knuth\]
- landlab/mcflugen/fix-issue-242 \[#243\] \[Eric Hutton\]
- Changed to deploy on release branch. \[mcflugen\]
- landlab/SiccarPoint/fix-issue-237 \[#239\] \[Dan Hobley\]
- landlab/mcflugen/fix-flexure-init \[#231\] \[Jordan Adams\]
- landlab/jadams15/fix_node_links \[#234\] \[Eric Hutton\]
- merge commit \[Jenny Knuth\]
