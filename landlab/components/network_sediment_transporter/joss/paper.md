---
title: 'NetworkSedimentTransporter: A Landlab submodule for coarse bed material transport through river networks'
tags:
  - Python
  - Landlab
authors:
  - name: Allison M. Pfeiffer
    orcid: --
    affiliation: 1

  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 2, 3

  - name: Jon Czuba
    orcid: ---
    affiliation: 4

affiliations:
  - name: Western Washington University, Geology Department
    index: 1
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 2
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 3
  - name: Virginia Tech
    index: 4

date: 0 November 2019
bibliography: papers.bib
---

# Summary

Coarse sediment transports downstream through river networks. The transport rate of any particular sediment grain on the river bed surface is a function of both the hydraulics of that reach of river and the size distribution of the other grains in the reach. As sediment transports, grains may be deposited or eroded, burying and exposing other grains, and in the process changing the elevation and slope of each segment of river. This process of river channel evolution through the process of sediment transport is referred to as morphodynamics. Computational morphodynamic models allow for the prediction of sediment pulse transport, such as that which occurs after dam removal or catastrophic landsliding events, as well as the prediction of changes in river channel bed surface grain size.

Over the past two decades, computational morphodynamic models have ...
[**JON, could you write this paragraph?** ]
Cui et al. [2006] developed and tested a model for sediment transport following dam removal, separately tracking the grain size of the active layer, subsurface, and bedload.

Ferguson et al. [2015] developed (**did they really develop, or did that happen in an earlier paper?**)
a similar model workflow to model the propagation of placer mine tailings through 500 km of the Fraser River in the absence of knowledge of initial channel bed grain size.
Enter network models:
CASCADE toolbox..
Czuba (2018) introduced a network-based, Lagrangian bed material morphodynamic model that tracks the motion of individuals “parcels” of sediment through the river network. This approach improves on the existing morphodynamic models by: (1) accounting for the full river network, rather than a single longitudinal profile, (2) allowing the user to ‘tag’ particular sediment inputs and track their fate through time.


These existing models, however, have two notable drawbacks: 1) they are written in a proprietary scripting language, and 2) they are not explicitly designed to be interoperable with other earth surface models, such as streamflow or landslide models. Here, we present software that overcomes these two drawbacks, translating the Czuba network sediment transport model into Landlab, a modular, Python-based package for the modeling of earth surface dynamics.  

The model presented here is nearly identical to the Czuba model, with a small number of minor added functionalities:
Abrasion
Active layer depth that varies with shear stress



Paragraph about Landlab:
Open Source Python package
Earth surface dynamics
[@Hobley2017Creative].
modular framework built of process components, primarily used for landscape evolution
New: Network model grid.



This contribution describes ``NetworkSedimentTransporter``, the Landlab submodule designed to model the motion of sediment "parcels" through the river network. Sediment 'parcels' are represented as items within the DataRecord, an xarray dataset-based landlab utility [cite]. Each parcel represents a package of sediment grains of common attributes such as grain diameter, lithology, and density. The parcel transports, is buried, and is eroded as a coherent unit. The river network is represented by the NetworkModelGrid, a landlab [something???], in which segments of river are represented as links, which are joined as nodes.   As parcels transport through links on the network, the elevation of nodes and slope of the links evolves according to addition and removal (deposition and erosion) of parcels from the surrounding links.



*********[All I did here was replace  words...]
Source code for ``NetworkSedimentTransporter`` is available as part of the [Landlab python package](https://github.com/landlab/landlab) and can be found in
the [``NetworkSedimentTransporter`` submodule](https://github.com/landlab/landlab/tree/release/landlab/components/network_sediment_transporter).
The ``NetworkSedimentTransporter`` submodule is documented using Docstrings, and the documentation can be found on the Landlab [ReadTheDocs site] (https://landlab.readthedocs.io/en/release/landlab.components.lithology.html).
Unit and docstring tests provide XXX% coverage of this submodule. [Pull Request # XXX](https://github.com/landlab/landlab/pull/XXX) brought the ``NetworkSedimentTransporter`` submodule into the core Landlab source code. The first release version of Landlab that includes the ``NetworkSedimentTransporter`` submodule is tagged as v2.??.


The Landlab project maintains a separate repository containing tutorials that
introduce core concepts and the use of individual submodules. In addition to the
source code, a [Jupyter Notebook introducing the use of NetworkSedimentTransporter](https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/xxxxxxxxxx.ipynb)
is now part of the Landlab tutorials repository. This tutorial was brought into
the repository with [Pull Request #
XX](https://github.com/landlab/tutorials/pull/19).

# Acknowledgements

The authors thank ____________ for _________

# References
