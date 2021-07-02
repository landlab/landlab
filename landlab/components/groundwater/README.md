## Welcome to the README for the GroundwaterDupuitPercolator component submodule.

[![status](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd/status.svg)](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd)

The GroundwaterDupuitPercolator is a component in Landlab for simulating shallow
subsurface flow. A [paper describing it](https://joss.theoj.org/papers/6936ca6851c622de48b2c5f6cf45a7bd)
was published in February 2020 in the Journal of Open Source Software. Here we
summarize installation, documentation, tutorials, tests, and getting help with
this component.

As this component lives within the larger Landlab package ecosystem, most of the
information below provides links into the [main Landlab documentation](https://landlab.readthedocs.io/).

### Installation
To use this component, you will need to install Landlab. Two options for
installation are available:
[a pre-packaged binary](https://landlab.readthedocs.io/en/master/install/index.html)
distributed through PyPI or conda-forge and a
[source code installation](https://landlab.readthedocs.io/en/master/development/install/index.html#developer-install).

The dependencies of the Landlab package are described [here](https://landlab.readthedocs.io/en/master/development/practices/dependencies.html).  

### Documentation
The documentation specific to this component is housed within the Landlab
documentation. There are two pages in the documentation that are most relevant
to this component:
- [The component API](https://landlab.readthedocs.io/en/master/reference/components/groundwater.html).
- [A page](https://landlab.readthedocs.io/en/master/reference/components/dupuit_theory.html#dupuit-theory)
describing the theory and numerical implementation of this component.

If you are new to Landlab and components, we recommend that you also look at the
[User Guide](https://landlab.readthedocs.io/en/master/user_guide/index.html),
in particular, the page on the [model grid](https://landlab.readthedocs.io/en/master/user_guide/grid.html), and [components](https://landlab.readthedocs.io/en/master/user_guide/components.html).

### Tutorials
There is a [Jupyter notebook in the Landlab Tutorials repository](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/groundwater/groundwater_flow.ipynb)
that describes the use of the `GroundwaterDupuitPercolator`.
The link takes you to a binder instance of this notebook. Its filepath within
the repository is `notebooks/tutorials/groundwater/groundwater_flow.ipynb`

A directory of all Landlab notebooks can be found (as a binder instance) [here](https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb)

### Tests of this Component
Along with the rest of the Landlab package, this component uses
[`pytest`](https://docs.pytest.org/en/latest/)
to  discover and run its tests. General information about running the Landlab
tests can be found [here](https://landlab.readthedocs.io/en/master/development/install/test.html#testing).

If you want to run the tests locally, you will need to use a
[source code installation](https://landlab.readthedocs.io/en/master/development/install/index.html#developer-install).

### Getting Help
If you have any questions, comments, issues, or bugs related to this submodule,
please [open an Issue](https://github.com/landlab/landlab/issues/new) so we can
respond.
