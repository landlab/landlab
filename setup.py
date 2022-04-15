#! /usr/bin/env python

import os
import re

import pkg_resources
from setuptools import Extension, find_packages, setup


numpy_incl = pkg_resources.resource_filename("numpy", "core/include")


def read(filename):
    with open(filename, "r", encoding="utf-8") as fp:
        return fp.read()


def read_requirements(filename):
    """Read requirements from a file, adding dependency specifications where needed.

    Dependency specifications could be included in the requirements file but,
    since the requirements file is also used by conda and conda doesn't
    understand dependency specifications, we add them here.
    """
    lines = [line.strip() for line in read(filename).splitlines()]
    requirements = []
    for requirement in lines:
        if requirement.startswith(("richdem", "bmi-topography")):
            requirement += "; sys_platform != 'win32' or python_version < '3.10'"
        requirements.append(requirement)

    return requirements


long_description = "\n\n".join(
    [
        read("README.rst"),
        read("AUTHORS.rst"),
        read("CHANGES.rst"),
    ]
)


def find_extensions(path="."):
    extensions = []
    for root, _dirs, files in os.walk(os.path.normpath(path)):
        extensions += [
            os.path.join(root, fname) for fname in files if fname.endswith(".pyx")
        ]
    return [
        Extension(re.sub(re.escape(os.path.sep), ".", ext[: -len(".pyx")]), [ext])
        for ext in extensions
    ]


install_requires = read_requirements("requirements.txt")
extras_require = {
    "tests": read_requirements("requirements-testing.txt"),
    "docs": read_requirements("requirements-docs.txt"),
    "dev": read_requirements("requirements-dev.txt"),
    "notebooks": read_requirements("requirements-notebooks.txt"),
}

setup(
    name="landlab",
    version="2.5.1.dev0",
    author="Eric Hutton",
    author_email="eric.hutton@colorado.edu",
    url="https://github.com/landlab",
    description="Plugin-based component modeling tool.",
    long_description=long_description,
    python_requires=">=3.7",
    install_requires=install_requires,
    extras_require=extras_require,
    include_package_data=True,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    keywords=[
        "bmi",
        "component modeling",
        "earth science",
        "gridding engine",
        "model coupling",
        "numerical modeling",
    ],
    packages=find_packages(),
    package_data={
        "": [
            "tests/*txt",
            "data/*asc",
            "data/*nc",
            "data/*shp",
            "test/*shx",
            "data/*dbf",
            "preciptest.in",
            "test_*/*nc",
            "test_*/*asc",
        ]
    },
    entry_points={"console_scripts": ["landlab=landlab.cmd.landlab:main"]},
    include_dirs=[numpy_incl],
    ext_modules=find_extensions("landlab"),
)
