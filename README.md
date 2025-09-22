![[DOI][doi-link]][doi-badge]
![[Documentation][rtd-link]][rtd-badge]
![[Coverage][coveralls-link]][coveralls-badge]
![[Testing][test-link]][test-badge]
![[pre-commit.ci status][lint-link]][lint-badge]


[coveralls-badge]: https://coveralls.io/repos/landlab/landlab/badge.png
[coveralls-link]: https://coveralls.io/r/landlab/landlab
[doi-badge]: https://zenodo.org/badge/DOI/10.5281/zenodo.3776837.svg
[doi-link]: https://doi.org/10.5281/zenodo.3776837
[lint-badge]: https://results.pre-commit.ci/badge/github/landlab/landlab/master.svg
[lint-link]: https://results.pre-commit.ci/latest/github/landlab/landlab/master
[rtd-badge]:https://readthedocs.org/projects/landlab/badge/?version=latest
[rtd-link]: https://landlab.csdms.io
[test-badge]: https://github.com/landlab/landlab/actions/workflows/test.yml/badge.svg
[test-link]: https://github.com/landlab/landlab/actions/workflows/test.yml

# Landlab

## What does Landlab do?

<!-- start-intro -->

Landlab is an open-source Python-language package for numerical modeling of
Earth surface dynamics. It contains

- A gridding engine which represents the model domain. Regular and irregular
  grids are supported.
- A library of process components, each of which represents a physical process
  (e.g., generation of rain, erosion by flowing water). These components have
  a common interface and can be combined based on a user's needs.
- Utilities that support general numerical methods, file input/output, and
  visualization.

In addition Landlab contains a set of Jupyter notebook tutorials providing
an introduction to core concepts and examples of use.

Landlab was designed for disciplines that quantify Earth surface dynamics such
as geomorphology, hydrology, glaciology, and stratigraphy. It can also be used
in related fields. Scientists who use this type of model often build
their own unique model from the ground up, re-coding the basic building blocks
of their landscape model rather than taking advantage of codes that have
already been written. Landlab saves practitioners from the need for this kind
of re-invention by providing standardized components that they can re-use.

Watch the webinar [Landlab Toolkit Overview](https://csdms.colorado.edu/wiki/Presenters-0407)
at CSDMS to learn more.

<!-- end-intro -->

______________________________________________________________________

[Read the documentation on ReadTheDocs!](https://landlab.csdms.io/)

______________________________________________________________________

## Installation

To install the latest release of *landlab* using *pip*, simply run the following
in your terminal of choice:

```bash
$ pip install landlab
```

For a full description of how to install *Landlab*, including using *mamba*/*conda*,
please see the documentation for our [installation instructions].

## Source code

If you would like to modify or contribute code to *Landlab* or use the very latest
development version, please see the documentation that describes how to
[install landlab from source].

## Are there any examples of using Landlab I can look at?

The Landlab package contains a directory, `landlab/notebooks`, with
Jupyter Notebooks describing core concepts and giving examples of using components.
The file `landlab/notebooks/welcome.ipynb` provides a table of contents to
the notebooks and is the recommended starting place.
Additionally, there are a set of notebooks curated to teach physical processes
located in the directory `landlab/notebooks/teaching`.

### Run on Binder

To launch an instance of
Binder and [explore the notebooks click here].

To launch a Binder instance that goes straight to the [teaching notebooks click here].

### Run on EarthscapeHub

The Landlab notebooks can also be run on [EarthscapeHub].
Visit this link to learn how to sign up for a free account.
Try the example notebooks on the [explore][explore-hub-example] or
[jupyter][jupyter-hub-example] Hub instance.
Or, use the teaching notebooks on the [explore][explore-hub-teaching] or
[jupyter][jupyter-hub-teaching] Hub instance.
Be sure to run all notebooks with the *CSDMS* kernel.

[explore-hub-example]: https://explore.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fwelcome.ipynb&branch=master&depth=1
[explore-hub-teaching]: https://explore.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fteaching%2Fwelcome_teaching.ipynb&branch=master&depth=1
[jupyter-hub-example]: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fwelcome.ipynb&branch=master&depth=1
[jupyter-hub-teaching]: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fteaching%2Fwelcome_teaching.ipynb&branch=master&depth=1


## License

*landlab* is licensed under the MIT License.

## Citing Landlab

If you use any portion of Landlab, please see the documentation for our
[citation guidelines].

## Contact

<!-- start-contact -->

The recommended way to contact the Landlab team is with a
[GitHub Issue](https://github.com/landlab/landlab/issues).

- **Bug reports**: Please make an Issue describing the bug so we can address it, or work
  with you to address it. Please try to provide a [minimal, reproducible example](https://stackoverflow.com/help/minimal-reproducible-example).
- **Documentation**: If something in our documentation is not clear to you, please make an
  issue describing the what isn't clear. Someone will tag
  the most appropriate member of the core Landlab team. We will work to clarify
  your question and revise the documentation so that it is clear for the next user.

Keep in touch with the latest *landlab* news by following us on [Twitter](https://twitter.com/landlabtoolkit).

During workshops and clinics, we sometimes use the
[Landlab Slack channel](https://landlab.slack.com).

<!-- end-contact -->

[citation guidelines]: https://landlab.csdms.io/about/citing.html
[earthscapehub]: https://csdms.colorado.edu/wiki/JupyterHub
[explore the notebooks click here]: https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/welcome.ipynb
[install landlab from source]: https://landlab.csdms.io/install/
[installation instructions]: https://landlab.csdms.io/installation.html
[teaching notebooks click here]: https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/teaching/welcome_teaching.ipynb
