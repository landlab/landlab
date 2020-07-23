#! /usr/bin/env python

import os
import re
from distutils.extension import Extension

import pkg_resources
from setuptools import Extension, find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

import versioneer

numpy_incl = pkg_resources.resource_filename("numpy", "core/include")


def find_extensions(path="."):
    extensions = []
    for root, dirs, files in os.walk(os.path.normpath(path)):
        extensions += [
            os.path.join(root, fname) for fname in files if fname.endswith(".pyx")
        ]
    return [
        Extension(re.sub(re.escape(os.path.sep), ".", ext[: -len(".pyx")]), [ext])
        for ext in extensions
    ]


def register(**kwds):
    import httplib, urllib

    data = urllib.urlencode(kwds)
    header = {
        "Content-type": "application/x-www-form-urlencoded",
        "Accept": "text/plain",
    }
    conn = httplib.HTTPConnection("csdms.colorado.edu")
    conn.request("POST", "/register/", data, header)


def register_landlab():
    try:
        from sys import argv
        import platform

        data = {
            "name": "landlab",
            "version": __version__,
            "platform": platform.platform(),
            "desc": ";".join(argv),
        }
        register(**data)
    except Exception:
        pass


class install_and_register(install):
    def run(self):
        install.run(self)
        register_landlab()


class develop_and_register(develop):
    def run(self):
        develop.run(self)
        register_landlab()


setup(
    name="landlab",
    version=versioneer.get_version(),
    author="Eric Hutton",
    author_email="eric.hutton@colorado.edu",
    url="https://github.com/landlab",
    description="Plugin-based component modeling tool.",
    long_description=open("README.rst").read(),
    python_requires=">=3.6",
    setup_requires=["cython", "numpy"],
    install_requires=open("requirements.txt", "r").read().splitlines(),
    include_package_data=True,
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Cython",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering :: Physics",
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
    cmdclass=versioneer.get_cmdclass(
        {"install": install_and_register, "develop": develop_and_register}
    ),
    entry_points={"console_scripts": ["landlab=landlab.cmd.landlab:main"]},
    include_dirs=[numpy_incl],
    ext_modules=find_extensions("landlab"),
)
