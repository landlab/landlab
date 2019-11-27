---
title: 'NetworkSedimentTransporter: A Landlab submodule for bed material transport through river networks'
tags:
  - Python
  - Landlab
authors:
  - name: Allison M. Pfeiffer
    orcid: 0000-0002-3974-132X
    affiliation: 1

  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 2, 3

  - name: Jon Czuba
    orcid: 0000-0002-9485-2604
    affiliation: 4

  - name: Eric Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 3, 5

affiliations:
  - name: Western Washington University, Geology Department
    index: 1
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 2
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 3
  - name: Virginia Tech
    index: 4
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 5
date: 0 November 2019
bibliography: papers.bib
---

# Summary

Coarse sediment moves downstream through river networks. The transport rate of any particular sediment grain on the river bed surface is a function of both the hydraulics of that reach of river and the size distribution of the other grains in the reach. As sediment moves through a river system, grains may be deposited or eroded, burying and exposing other grains, and in the process changing the elevation and slope of each segment of river. This process of river channel evolution through the process of sediment transport is referred to as morphodynamics [@cite]. Computational morphodynamic models allow for the prediction of sediment pulse transport, such as that which occurs after dam removal [@citeCui] or landsliding events[@An2017], as well as the prediction of changes in river channel bed surface grain size.

Over the past two decades, computational morphodynamic models have ...
[**JON, could you write two sentences here about the history of morphodynamic models?** ]
Cui et al. [2006] ; Ferguson et al. [2015]
Enter network models:CASCADE...
Czuba (2018) introduced a network-based, Lagrangian bed material morphodynamic model that tracks the motion of individual units (referred to as “parcels”) of sediment through the river network. This approach improves on the existing morphodynamic models by: (1) accounting for the full river network, rather than a single longitudinal profile, (2) allowing the user to ‘tag’ particular sediment inputs and track their fate through time. These existing network sediment transport models, however, have two notable drawbacks: 1) they are written in a proprietary scripting language, and 2) they are not explicitly designed to be interoperable with other earth surface models, such as streamflow or landslide models.

Here, we present software that overcomes these two drawbacks, translating network sediment transport model of @Czuba2018 into Landlab, a modular, Python-based package for the modeling of Earth surface dynamics. Landlab is an Open Source Python package for modeling earth surface processes [@Hobley2017Creative]. It was designed as a modular framework, hosting a variety of process components such as flow routing, hillslope diffusion, and stream power erosion that function on a common set of landscape model grids. The original implementation of this model was in MatLab. ``NetworkSedimentTransporter`` is the newest of these components. We first describe computational infrastructure created in order to create the ``NetworkSedimentTransporter`` and then describe the new componet itself.

The creation of the ``NetworkSedimentTransporter`` required the creation of two new data structures in the Landlab framework. First, the ``NetworkModelGrid``, which represents the model domain as connected nodes and links. Second, the ``DataRecord`` stores a generic set of items in time and on the model grid. It is used here to represent the sediment that move into, through, and out of the modeled domain.

In the ``NetworkSedimentTransporter``, sediment is represented as "parcels"-a quantity of sediment grains with common attributes such as grain diameter, lithology, and density. The parcel transports, is buried, and is eroded as a coherent unit. The river network is represented a series of links and nodes on a ``NetworkModelGrid``. Each time the ``NetworkSedimentTransporter`` is run forward in time, the set of parcels that are in active transport is identified based on hydrologic conditions, parcels transport distances are determined, and the parcels move through links on the network. As a result of parcel redistribution, the elevation of nodes and slope of the links evolves.

The ``NetworkSedimentTransporter`` adds new functionality to the original implementation  by @Czuba2018. Specifically, it permits for variable sediment parcel density and bed material abrasion. The latter is calculating the loss of particle mass during transport downstream as:

$W_x = W_0 e^{\left(\alpha x \right)}$

Where $x$ is the downstream transport distance, $\alpha$ is the abrasion rate, and $W_x$ and $W_0$ are the resulting and original sediment parcel masses, respectively.

The final element of new functionality is a method for calculating the [active layer thickness]. As in many sediment transport models, @Czuba2018 represents the mobile portion of the grains on the riverbed at any given time as an "active layer" of constant thickness. All grains in this layer are transported, while all grains below this layer remain unmoved. Here, we incorporate the formulation of Wong et al. (2007) to calculate an active layer thickness for each link in the network at each timestep as a function of Shields stress and median grain diameter.
[Anything else we added?]

Source code for ``NetworkSedimentTransporter`` is available as part of the [Landlab python package](https://github.com/landlab/landlab) and can be found in
the [``NetworkSedimentTransporter`` submodule](https://github.com/landlab/landlab/tree/release/landlab/components/network_sediment_transporter). The first release version of Landlab that includes the ``NetworkSedimentTransporter`` submodule is tagged as v2.??.

The Landlab project maintains a separate repository containing tutorials that
introduce core concepts and the use of individual submodules. In addition to the
source code, a [Jupyter Notebook introducing the use of NetworkSedimentTransporter](https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/xxxxxxxxxx.ipynb)
is now part of the Landlab tutorials repository.

# Acknowledgements

The authors thank ____________ for _________

Barnhart supported by an NSF EAR Postdoctoral Fellowship (NSF Award Number 1725774).

# References
