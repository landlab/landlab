(style-enforcement)=

# Style Guidelines and Enforcement

## General Coding Style

- Please stick to the coding style described by
  [PEP8](https://www.python.org/dev/peps/pep-0008/). PEP8 is one of
  the standard worldwide stylistic conventions for coding in Python.
- Class and function docstrings should follow the [numpydoc
  conventions](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard).
- Further, Landlab-specific advice for developing your own components
  can be found in the {ref}`component development guide <dev-contributing>`.

Many modern text editors, e.g., [Atom](https://atom.io), provide interactive tools to
examine your code as you edit it and highlight where it diverges from
standards.

## Style enforcement

All tools used for development are specified in our `environment-dev.yml` file.
If you followed the
{ref}`developer installation instructions <install>` you have
everything you need in the `landlab_dev` conda environment.

Currently we check for all the [flake8
violations](https://pycodestyle.readthedocs.io/en/latest/intro.html#error-codes)
and [pycodestyle
violations](http://flake8.pycqa.org/en/latest/user/error-codes.html)
except for (as defined in our `setup.cfg`)

- E203: whitespace before
- E501: line too long (n > 88 characters)
- W503: line break before binary operator

To format files to meet these standards, we recommend using
[isort](https://pypi.org/project/isort/) +
[black](https://github.com/psf/black).

You can run

```bash
$ make pretty
```

from the main landlab directory in a terminal to run both
[isort](https://pypi.org/project/isort/) and
[black](https://github.com/psf/black) on your code.

To check if your files meet the standards that are enforced under
continuous integration, we use
[flake8](http://flake8.pycqa.org/en/latest/). You can run

```bash
make lint
```

from the top level directory or

```bash
flake8 <file-to-check>
```
