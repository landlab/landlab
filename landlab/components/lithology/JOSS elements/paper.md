---
title: Lithology and Litholayers: 'Landlab components for spatially variable rock properties'
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
date: 15 August 2018
bibliography: papers.bib
---

# Summary

The surface of the Earth reflects the competing advection of rock by tectonic processes and the erosion of rock by wind, water, and ice. The inherent strength and fracture spacing of rock influence erosion rates by changing the rate and physical process by which rock is detached and turned into mobile material. Variations in the rock properties over space and with depth reflect the legacy of sedimentary deposition and tectonic deformation.

The Lithology and LithoLayers components in the Landlab toolkit were designed to permit Landlab users to implement variations in rock properties in Landlab models. Landlab is an Open Source Python package that provides a framework for the development of 2D numerical models, typically in Earth surface dynamics [@HobleyLandlab]. Landlab was designed as a  modular framework in which different process components can be mixed and matched to construct a model based on a user's needs.

Lithology is a three dimensional representation of layered. Layers may have spatially variable thickness and multiple attributes. Material can be removed through erosion or added to through deposition.

Material or Event Layers.

LithoLayers is a

Source code for Lithology and Litholayers is available as part of the [Landlab python package](). A tagged version of Landlab that includes this functionality is archived in Zenodo with the linked DOI:[@ZenodoSourceCode]. In addition to the source code a Jupyter Notebook introducing the use of Lithology and Litholayers is now part of the Landlab tutorials repository. A tagged version of Landlab that includes this functionality has been archived with the following DOI: [@ZenodoNotebook].

# Acknowledgements
The authors thank Adam Forte, Matt Rossi, and Brian Yanites for helpful discussions during the development of this code and the accompanying Jupyter notebooks.

# References
