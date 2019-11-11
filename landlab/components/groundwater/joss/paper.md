---
title: 'GroundwaterDupuitPercolator: A Landlab component for groundwater flow'
tags:
  - Python
  - Landlab
authors:
  - name: David G. Litwin
    orcid: 0000-0002-8097-4029
    affiliation:
  - name: Gregory E. Tucker
    orcid: 0000-0003-0364-5800
    affiliation:
  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation:
  - name: Ciaran J. Harman
    orcid: 0000-0002-3185-002X
    affiliation:
affiliations:
  - name: Johns Hopkins University, Department of Environmental Health and Engineering
  - index:
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index:
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index:
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index:
  - name: University of Colorado at Boulder, Institute for Arctic and Alpine Research
    index:
  - name: Johns Hopkins University, Department of Earth and Planetary Science
  - index:   
date:
bibliography: papers.bib
---

# Summary
A large portion of the water that enters a catchment as precipitation percolates through soil and rock before exiting to water bodies or returning to the atmosphere as evapotranspiration. In many places, the discharge of water stored in the subsurface is a primary source of streamflow, and thus controls the ways in which catchments respond to stochastic variations in precipitation and climate [@beck_global_2013, @jasechko_pronounced_2014]. Previous studies have shown the importance of groundwater for a diverse range of processes, from transpiration to solute export [@maxwell_connections_2016, @van_verseveld_role_2009], and across diverse timescales, from rainfall-runoff response to landscape evolution [@huang_modelling_2006, @sklash_role_1979]. Of particular importance to landscape evolution, groundwater can be an important control on the occurrence of overland flow, as the interaction of the water table with the ground surface controls the spatial extent of saturation and groundwater return flow [@dunne_partial_1970].
Variably-saturated groundwater flow is often assumed to be governed by the Richards equation, which describes how water content and/or total energy potential evolve in an idealized porous medium due to fluxes of water driven by gradients in total potential. Numerical solutions to the Richards equation are computationally expensive [e.g. @kirkland_algorithms_1992], often limiting their applications. For computational efficiency, we use the widely applied Dupuit-Forcheimer approximation, which simplifies the Richards equation when aquifers are laterally extensive in comparison to their thickness, and the capillary fringe above the water table is relatively thin [e.g. @brutsaert_hydrology:_2005, @troch_hillslope-storage_2003]. In this case, the water table is modeled as a free surface with groundwater flow driven by gradients in water table elevation. When the model assumptions are valid, this greatly reduces the model complexity while still producing water table elevations and discharges comparable to Richards equation solutions [@hilberts_hillslope-storage_2004].
Implementations of the Dupuit-Forcheimer model often encounter numerical instabilities as the water table intersects the surface and groundwater return flow occurs through a seepage face, due to the presence of a discontinuity in the energy gradient from inside the hillslope (where it is determined by the water table) to the seepage face (where it is determined by the topography). Marçais et al. (2017) introduced a regularization that smooths the transition between subsurface flow and surface flow. Although introduced for numerical stability and not based on physical principles, this smoothing may reproduce the effect of subgrid heterogeneity where saturation within a grid cell is unlikely to be homogenous or well-reproduced by a binary condition (saturated vs unsaturated).
GroundwaterDupuitPercolator:
The GroundwaterDupuitPercolator solves the governing groundwater flow equations with an explicit, forward in time finite volume method, using the @marcais_dynamic_2017 regularization at seepage faces. While the explicit method limits the maximum timestep that can be used without jeopardizing stability, the model includes an adaptive timestep solver that subdivides the user-provided timestep in order to satisfy a Courant-type stability criterion.
The GroundwaterDupuitPercolator can be implemented on both regular (e.g. rectangular and hexagonal) and irregular grids determined by the user. Recharge, hydraulic conductivity, and porosity may be specified as single values, uniform over the model domain, or as vectors on the nodes (recharge, porosity) or links (hydraulic conductivity) of the grid. Link hydraulic conductivity can also be specified from a two-dimensional hydraulic conductivity tensor using an included function. For mass balance calculations, the model includes methods to determine the total groundwater storage in the model domain, the total recharge flux in, and total groundwater and surface water fluxes leaving through the boundaries.
The GroundwaterDupuitPercolator is implemented in Landlab, a Python-based open source Earth surface modeling toolkit [@Hobley2017Creative]. Landlab has a modular framework, which allows for easy coupling of different process components to meet the needs of the modeler. For example, the surface water seepage from the GroundwaterDupuitPercolator can be passed to the FlowAccumulator module to route overland flow and calculate discharge at nodes. A diverse array of components are available, yielding many possibilities for model coupling that have not yet been explored. Given the importance of groundwater for many Earth surface processes, this component is an important contribution to the Landlab environment.  


# Acknowledgements


# References
