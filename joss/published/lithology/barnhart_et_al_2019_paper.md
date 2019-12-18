---
title: 'Lithology: A Landlab submodule for spatially variable rock properties'
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
    affiliation: 1, 2, 3
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

The surface of the Earth reflects the competing advection of rock by tectonic
processes and the erosion of rock by wind, water, and ice. Rock properties
influence erosion rates by changing the processes responsible for erosion and
the rate at which rock is weathered, detached, and turned into mobile sediment.
Variations in the rock properties over space and with depth reflect the legacy
of sedimentary deposition and tectonic deformation. Long-term landscape
evolution modeling experiments that include the impact of spatially and
temporally variable rock characteristics can be used to identify the impact of
rock strength patterns on other geologic observables such as topography, erosion
rates, and detrital mineral records [e.g., @Forte2016Complexites;
@Perne2017Steady]. Identifying these relationships allows for better
interpretations of the geologic record.

Landlab is an Open Source Python package that provides a framework for the
development of 2D numerical models, typically in Earth surface dynamics
[@Hobley2017Creative]. Landlab was designed as a  modular framework in which
different process components can be mixed and matched to construct a model based
on a user's needs. Prior work on spatially variable lithology in landscape
evolution has been done using a modified version of the channel-hillslope
integrated landscape development [CHILD, @tucker2001channel] model [e.g.,
@Forte2016Complexites] and the FastScape V5 model [@braun2013very;
@Perne2017Steady]. To provide these capabilities within the Landlab framework,
there is a need for a Landlab submodule that can treat spatial variations in
rock materials.  

This contribution describes ``Lithology``, a Landlab submodule designed to
support the representation of 3D variations in rock material properties within
the Landlab framework. It includes two classes: ``Lithology`` is a generic
representation of spatially varying rock material, and ``LithoLayers`` is a
derived class that treats parallel layers of material variations. In both
classes, each rock type may have multiple attributes. Rock layers may be removed
through erosion, or added to through deposition. Two options for the underlying
datastructure are supported: "event layers", in which the data structure stores
each time-step as an event even if there is no material in the layer, or
"material layers", in which entries in the datastructure represent contiguous
material of the same property, but not necessarily the same age. This second
option is more memory efficient but does not record the transient dynamics of
erosion and deposition.

Source code for ``Lithology`` and ``Litholayers`` is available as part of the
[Landlab python package](https://github.com/landlab/landlab) and can be found in
the [``Lithology``
submodule](https://github.com/landlab/landlab/tree/release/landlab/components/lithology).
The ``Lithology`` submodule is documented using Docstrings, and the
documentation can be found on the Landlab ReadTheDocs site. One page exists for
the [Lithology
component](https://landlab.readthedocs.io/en/release/landlab.components.lithology.html)
and a second for the [LithoLayers
component](https://landlab.readthedocs.io/en/release/landlab.components.litholayers.html).
Unit and docstring tests provide 100% coverage of this submodule. [Pull Request #
674](https://github.com/landlab/landlab/pull/674) brought the ``Lithology``
submodule into the core Landlab source code. The first release version of
Landlab that includes the ``Lithology`` submodule is tagged as v1.5.4. The
concept DOI for Landlab is archived in Zenodo with the linked DOI: [@ZenodoLithologySourceCode]
and the archive for this manuscript points to the Zenodo archive of v1.5.4.

The Landlab project maintains a separate repository containing tutorials that
introduce core concepts and the use of individual submodules. In addition to the
source code, a [Jupyter Notebook introducing the use of Lithology and
Litholayers](https://nbviewer.jupyter.org/github/landlab/tutorials/blob/release/lithology/lithology_and_litholayers.ipynb)
is now part of the Landlab tutorials repository. This tutorial was brought into
the repository with [Pull Request #
19](https://github.com/landlab/tutorials/pull/19). The first release version of
the Landlab tutorials that includes this notebooks is tagged as v1.5.4 and is
archived in Zenodo with the linked DOI: [@ZenodoLithologyNotebook].

# Acknowledgements

The authors thank Adam Forte, Matt Rossi, and Brian Yanites for helpful
discussions during the development of this code and the accompanying Jupyter
notebooks.

# References
