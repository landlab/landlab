(development)=

# Developer Guide

If you're intending to make changes to the Landlab code base, or want to
develop your own components, this set of pages will provide information you
need.

## Quick Start Tutorial
New contributors are encouraged to watch the tutorial video for an overview of
the general workflow for how to contribute to the Landlab GitHub repository. This tutorial video is part of the CSDMS webinar series:
[*"From issue to pull request: how to contribute to CSDMS' open-source community code
repositories"*](https://youtu.be/YuFhdzVZUFI).

Please click the title above to watch the video,
and a summary of the steps covered in this tutorial can be found
[here](https://csdms.colorado.edu/wiki/Presenters-0697).

## Key Development Recommendations

The most important things to remember are to:

- Start by creating a fork of Landlab's repository;
- Make changes to the source code on a development branch, not the default
  `master` branch; and
- Keep your fork's `master` and development branches up to date with changes
  in the main Landlab repository.

## Supported Python Versions

*Landlab* supports the [current Python version](https://devguide.python.org/versions/)
plus the two previous versions.

If you need to introduce a new dependency, that dependency must be compatible
with the oldest supported Python version and be available on Linux, Mac, and Windows.

```{toctree}
:hidden: true
:maxdepth: 2

contribution/index
practices/index
package_organization
```

## Quick Links For Package Maintenance

There are a few pages that are particularly important for the maintenance of
the package. These are:

- {ref}`directory organization <organization>`
- {ref}`dependency organization <dependencies>`, and
- {ref}`release workflow <dev-releases>`

If package maintainers change any of these, the prior pages likely need
updating.
