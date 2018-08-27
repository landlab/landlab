---
title: 'Lithology: 'A Landlab submodule for spatially variable rock properties'
tags:
  - Python
  - Landlab
authors:
  - name: Katherine R. Barnhart
    orcid: 0000-0001-5682-455X
    affiliation: 1, 2
  - name: Eric Hutton
    orcid: 0000-0002-5864-6459
    affiliation: 3, 4
  - name: Nicole M. Gasparini
    orcid: 0000-0002-0803-3697
    affiliation: 5
  - name: Gregory E. Tucker
    orcid: 0000-0003-0364-5800
    affiliation: 1, 2
affiliations:
  - name: University of Colorado at Boulder, Department of Geological Sciences
    index: 1
  - name: University of Colorado at Boulder, Cooperative Institute for Research in Environmental Sciences
    index: 2
  - name: University of Colorado at Boulder, Community Surface Dynamics Modeling System Integration Facility
    index: 3
  - name: University of Colorado at Boulder, Institute for Arctic and Alpine Research
    index: 4
  - name: Tulane University, Department of Earth and Environmental Sciences
    index: 5
date: 16 August 2018
bibliography: papers.bib
---

# Summary

The surface of the Earth reflects the competing advection of rock by tectonic processes and the erosion of rock by wind, water, and ice. The inherent strength and fracture spacing of rock influence erosion rates by changing the rate and physical process by which rock is detached and turned into mobile material. Variations in the rock properties over space and with depth reflect the legacy of sedimentary deposition and tectonic deformation. Long-term landscape evolution modeling experiments that include the impact of spatially and temporally variable rock characteristics can be used to identify the impact of rock strength patterns on other geologic observables such as topography, erosion rates, and detrital mineral records (e.g. [@Forte2016Complexites]). Identifying these relationships allows for better interpretations of the geologic record.

Landlab is an Open Source Python package that provides a framework for the development of 2D numerical models, typically in Earth surface dynamics ([@Hobley2017Creative]). Landlab was designed as a  modular framework in which different process components can be mixed and matched to construct a model based on a user's needs. Much of the prior work on spatially variable lithology in landscape evolution was done using a modified version of the channel-hillslope integrated landscape development (CHILD, [@tucker2001channel]) model (e.g.[@Forte2016]). This modified version of CHILD is not available in an open repository. As Landlab was designed to superceed CHILD, there is a need for a Landlab submodule that can treat spatial variations in rock materials.  

This contribution describes ``lithology``, a Landlab submodule designed to support the representation of 3D variations in rock material properties within the Landlab framework. It includes two classes: ``Lithology`` is a generic representation of spatially varying rock material and ``LithoLayers`` is a derived class that treats parallel layers of material variations. In both classes, each rock type may have multiple attributes, be removed through erosion, or added to through deposition. Two options for the underlying data-structure are supported: event layers, in which the data structure stores each time-step as an event, or material layers, in which material must be present in each layer. This second option is more memory efficient but does not record the transient dynamics of erosion and deposition.

Source code for ``Lithology`` and ``Litholayers`` is available as part of the [Landlab python package](https://github.com/landlab/landlab). The ``lithology`` submodule is documented using Docstrings and the documentation can be found on the [Landlab ReadTheDocs site](https://landlab.readthedocs.io/en/latest/) (note: this URL will be updated to link to the submodule documentation instead of the package documentation). Unit and docstring tests provide 100% coverage of this submodule. [Pull Request # 674](https://github.com/landlab/landlab/pull/674) brought the ``lithology`` submodule into the core Landlab source code (note: as of initial JOSS submission this PR is still active as we anticipate making revisions before merging). A tagged version of Landlab that includes the ``lithology`` submodule is archived in Zenodo with the linked DOI:[@ZenodoLithologySourceCode] (note: doi not active for initial JOSS submission).

The Landlab project maintains a separate repository containing tutorials introducing core concepts and the use of individual submodules. In addition to the source code a Jupyter Notebook introducing the use of Lithology and Litholayers is now part of the Landlab tutorials repository. This tutorial was brought into the repository with [Pull Request # 19](https://github.com/landlab/tutorials/pull/19)(note: as of initial JOSS submission this PR is still active as we anticipate making revisions before merging). A tagged version of Landlab that includes this functionality has been archived with the following DOI: [@ZenodoLithologyNotebook] (note: doi not active for initial JOSS submission).

# Acknowledgements
The authors thank Adam Forte, Matt Rossi, and Brian Yanites for helpful discussions during the development of this code and the accompanying Jupyter notebooks.

# References
