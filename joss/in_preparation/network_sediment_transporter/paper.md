---
title: 'NetworkSedimentTransporter: A Landlab component for bed material transport through river networks'
tags:
  - Python
  - Landlab
authors:
  - name: Allison M. Pfeiffer
    orcid: 0000-0002-3974-132X
    affiliation: 1

  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 2, 3, 4

  - name: Jonathan A. Czuba
    orcid: 0000-0002-9485-2604
    affiliation: 5

  - name: Eric W. H. Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 6

affiliations:
  - name: Western Washington University, Geology Department
    index: 1
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 2
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 3
  - name: "Present affiliation: U.S. Geological Survey, Landslide Hazards Program, 1711 Illinois St., Golden, CO 80401"
    index: 4
  - name: Virginia Tech, Department of Biological Systems Engineering and The Global Change Center
    index: 5
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 6
date: 27 Apr 2020
bibliography: papers.bib
---

# Summary

Coarse sediment (sand, gravel, and cobbles) moves downstream through river networks. The transport rate of any particular sediment grain on the river bed surface is a function of both the hydraulics of that reach of river and the size distribution of the other grains in the reach. As sediment moves through a river system, grains may be deposited or eroded, burying and exposing other grains, and in the process changing the elevation and slope of each segment of river. This process of river channel evolution through the process of sediment transport is referred to as morphodynamics [@ParkerEbook]. Computational morphodynamic models allow for the prediction of sediment pulse transport, such as that which occurs after dam removal [@Cuietal2006a; @Cuietal2006b; @Cui2007a] or landsliding events [@An2017; @Benda&Dunne1997], as well as the prediction of changes in river channel bed surface grain size [@Fergusonetal2015].

Most computational morphodynamic models take an Eulerian approach, which tracks changes in bed elevation through time as a function of the spatial gradient in sediment flux [e.g., @ParkerEbook]. These models directly compute bed elevation change and sediment flux throughout the domain. One of the major drawbacks with Eulerian morphodynamic models is the difficulty in being able to 'tag' individual sediment particles to answer questions about how an individual sediment particle/input may move, when it might arrive, and what affect it will have on river morphology when it arrives downstream. To overcome this drawback and to more easily extend morphodynamic models to entire river networks, recent work has focused on developing river-network based Lagrangian sediment transport models, which track the locations of individual sediment units on a river network.  

A more comprehensive overview of river-network based sediment transport models is described by @Czubaetal2017. Of most relevance to the work described herein, @Czuba2018 introduced a network-based, Lagrangian bed material morphodynamic model that tracks the motion of individual units (referred to as “parcels”) of sediment through a river network. The model presented by @Czuba2018 has been applied to post-wildfire debris-flow sediment movement through a river network in Utah [@Murphyetal2019]. Czuba's approach improves on the existing morphodynamic models by: (1) accounting for the full river network, rather than a single longitudinal profile, (2) allowing the user to ‘tag’ particular sediment inputs and track their fate through time. Despite its advances, this existing network sediment transport model, however, has two notable drawbacks: 1) it written in a proprietary scripting language (MATLAB), and 2) it is not explicitly designed to be interoperable with other Earth-surface models, such as streamflow or landslide models.

Here, we present software that overcomes these two drawbacks, translating and expanding upon the network sediment transport model of @Czuba2018 in Landlab, a modular, Python-based package for the modeling of Earth-surface dynamics. Landlab is an open-source Python package for modeling Earth-surface processes [@Hobley2017Creative; @Barnhart2020Short]. It was designed as a modular framework, hosting a variety of process components such as flow routing, hillslope diffusion, and stream power erosion that function on a common set of landscape model grids. The ``NetworkSedimentTransporter`` is the newest of these components. We first describe computational infrastructure built in order to create the ``NetworkSedimentTransporter`` and then describe the new component itself.

The creation of the ``NetworkSedimentTransporter`` required the addition of two new data structures in the Landlab framework. First, the ``NetworkModelGrid``, which represents the model domain as connected nodes and links. Second, the ``DataRecord``, which stores a generic set of items in time and on the model grid. It is used here to store all attributes associated with the sediment parcels that move into, through, and out of the model domain.

In the ``NetworkSedimentTransporter``, sediment is represented as "parcels"-a quantity of sediment grains with common attributes such as grain diameter, lithology, and density. Each parcel is transported, buried, and eroded as a coherent unit. The river network is represented as a series of links and nodes on a ``NetworkModelGrid``. Each time the ``NetworkSedimentTransporter`` is run forward in time, the set of parcels that are in active transport is identified based on the flow conditions and bed surface grain size in each link, transport distances are calculated for all active parcels based on the @WilcockCrowe2003 equations, and parcels move through links on the network by updating their locations based on their transport distances [@Czuba2018]. As a result of parcel redistribution, the elevation of nodes and slope of the links evolves [@Czubaetal2017; @Czuba2018].

Our implementation is not a direct translation of the model implemented in MATLAB and described in @Czuba2018. Here we add three new elements to the model: sediment density that varies across parcels, downstream bed-material abrasion, and enhanced capabilities for specifying the active layer thickness.

The use of the ``DataRecord`` attributes to store density and the abrasion-rate coefficient permits different values for each sediment parcel. The density influences which parcels are mobile and how far they move each timestep. Variable (rather than constant) density permits better representing study sites with lithologic variation. Similarly, different rock types may abrade at different rates. The abrasion-rate is calculated as the loss of particle mass (or volume, because density is constant within each parcel) during transport downstream as:

$W_x = W_0 \exp \left(\alpha x \right)$

Where $x$ is the downstream transport distance, $\alpha$ is the abrasion rate (for mass loss), and $W_x$ and $W_0$ are the resulting and original sediment parcel masses, respectively. The model tracks parcel volumes (not masses) so the actual implementation replaces $W_x$ and $W_0$ with volumes (e.g., $W_0=V_0\rho_s$, where $V_0$ is the original sediment parcel volume and $\rho_s$ is the rock density of the sediment in the parcel); however, the form of the equation for mass or volume is equivalent for a parcel with a constant sediment density (i.e., the $\rho_s$ on both sides of the equation cancel out). Furthermore, once a volume reduction of each parcel is computed, the model also updates the associated reduction in parcel sediment grain size as:

$D_x = D_0 \left(\frac{V_x}{V_0}\right)^{1/3}$

Where $D_x$ and $D_0$ are the resulting and original sediment parcel diameters, respectively.

Our final modification to @Czuba2018 is enhancing the methods used for calculating variable active layer thickness. Many sediment transport models [e.g., @Cui2007b; @Czuba2018] represent the mobile portion of the grains on the riverbed at any given time as an "active layer" of constant thickness. All grains in this layer are transported, whereas all grains below this layer are immobile. Within ``NetworkSedimentTransporter`` the user has the option to specify active layer thickness as a constant value or a multiple of the mean grain size in each link. Alternatively, we incorporated the formulation of @Wongetal2007 to calculate an active layer thickness for each link in the network at each timestep as a function of Shields stress and median grain diameter.

The ``NetworkSedimentTransporter`` component of Landlab is capable of routing mixed grain size sediment through river networks to answer questions about how sediment pulses move through river networks and when, where, and how they affect downstream reaches. The accessibility of this code within the Landlab framework will make it easier for future users to modify and contribute to its continual evolution.

Source code for ``NetworkSedimentTransporter`` is available as part of the [Landlab Python package](https://github.com/landlab/landlab) and can be found in
the [``NetworkSedimentTransporter`` component](https://github.com/landlab/landlab/tree/release/landlab/components/network_sediment_transporter). The first release version of Landlab that includes the ``NetworkSedimentTransporter`` component is tagged as v2.1.0.

The Landlab project maintains a separate repository containing tutorials that introduce core concepts and the use of individual components. In addition to the
source code, a set of Jupyter Notebooks introducing the use of NetworkSedimentTransporter
are now part of the Landlab tutorials repository:
- [Part  1: Introduction with a synthetic network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter.ipynb)
- [Part  2: Using a shapefile-based river network](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_sediment_transporter_shapefile_network.ipynb)
- [Part  3: Plotting options](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/network_sediment_transporter/network_plotting_examples.ipynb)

# Acknowledgements

Barnhart supported by an NSF EAR Postdoctoral Fellowship (NSF Award Number 1725774). Czuba was partially supported by NSF-EAR (1848672), Virginia Agricultural Experiment Station, and USDA Hatch program (1017457). Pfeiffer was supported by the NCED II Synthesis Postdoctoral program and NSF-PREEVENTS (NSF Award Number 1663859 to PI Istanbulluoglu). Landlab is supported by the National Science Foundation (NSF Award Numbers 1147454, 1148305, 1450409, 1450338, and 1450412) and by the Community Surface Dynamics Modeling System (NSF Award Numbers 1226297 and 1831623). The authors thank Associate Editor Kristen Thyng, along with Zoltán Sylvester and Evan Goldstein for their thorough review of this contribution in the midst of a pandemic.

# References
