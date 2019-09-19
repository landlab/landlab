[Landlab](http://landlab.github.io) | 
[[ About | About ]] |
[[ Examples | Examples ]] |
[[ User Guide | User-Guide ]] |
[ Reference Manual](http://landlab.readthedocs.org/en/latest/#developer-documentation) |
[[ Tutorials| Tutorials ]] |
[[ FAQs | FAQs ]] 

# Style enforcement 
Currently we check for all the [*flake8* violations](https://pycodestyle.readthedocs.io/en/latest/intro.html#error-codes) and [*pycodestyle* violations](http://flake8.pycqa.org/en/latest/user/error-codes.html) except for (as defined in our `setup.cfg`)
* E203: whitespace before ‘:’
* E501: line too long (n > 88 characters)
* W503: line break before binary operator

To format files to meet these standards, we recommend using [isort](https://pypi.org/project/isort/) + [black](https://github.com/ambv/black).

You can run 
```
make pretty
```

from the main landlab directory to run both [isort](https://pypi.org/project/isort/) and [black](https://github.com/ambv/black) on your code.

To check if your files meet the standards that are enforced under continuous integration, install [*flake8*](http://flake8.pycqa.org/en/latest/) and run 
```
make lint
```
from the top level directory or
```
flake8 <file-to-check>
```

# General Coding Style

* Please stick to the coding style described by [PEP8](http://www.python.org/dev/peps/pep-0008/). PEP8 is one of the standard worldwide stylistic conventions for coding in Python.

* Class and function docstrings should follow the [numpydoc conventions](https://numpydoc.readthedocs.io/en/latest/format.html#docstring-standard).
  
* Further, Landlab-specific advice for developing your own components can be found in the [component development guide](https://github.com/landlab/landlab/wiki/Develop-your-own-component).

If you want to check how well you are doing, please look at our [Landscape page](https://landscape.io). Landscape will grade the health of the landlab code with every push to GitHub. Many modern text editors, e.g., [Atom](https://atom.io), provide interactive tools to examine your code as you edit it and highlight where it diverges from standards.