#! /usr/bin/env python

import os
import re

import pkg_resources
from setuptools import Extension, setup


numpy_incl = pkg_resources.resource_filename("numpy", "core/include")


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


setup(include_dirs=[numpy_incl], ext_modules=find_extensions("landlab"))
