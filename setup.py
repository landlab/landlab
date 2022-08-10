""" setup.py is called by pip install -e . to compile the new cython files
into c or cpp platform-dependent executables.
It can also be called this way:
python setup.py build_ext --inplace.

To force all cython files to compile or prevent (as a function of their
modification), edit the force parameter to True or False in setup_cython.yml.
"""
from setuptools import Extension, setup
from Cython.Build import cythonize
import yaml
import os
import re
import numpy

with open("setup_cython.yml") as f:
    setup_cython = yaml.load(f, Loader=yaml.loader.SafeLoader)
roots = []
filenames = []

"""
To get the paths and filenames of the extensions to update,
we walk through the files included in walk_paths. If the directory names
or filenames contains one of the exclude_words_in_path or exclude_words_in_filename
they are not included in the list of extensions.
Otherwise, and if the filenames end by .pyx:
- either force is set to True, and the file is added to the extensions
- of force set to False, the file is added either if
-      - there is no c or cpp file
-      - there is a c or cpp file, but the timestamp is older than the one
-        of the .pyx file
"""
setup_select = setup_cython["cython_compile"]["file_selection"]
for walk_path in setup_select["walk_paths"]:
    for root, _dirs, fnames in os.walk(os.path.normpath(walk_path)):
        exclude = False
        if setup_select["exclude_words_in_path"] is not None:
            for word in setup_select["exclude_words_in_path"]:
                if word in root:
                    exclude = True
                    break
            if exclude:
                continue

        for fname in fnames:
            if not fname.endswith(".pyx"):
                continue
            if setup_select["exclude_words_in_filename"] is not None:
                for word in setup_select["exclude_words_in_filename"]:
                    if word in fname:
                        exclude = True
                        break
                if exclude:
                    continue

            pyx_path = os.path.join(root, fname)
            if setup_select["force"]:
                roots += [root]
                filenames += [fname]
                continue
            c_path = ""
            for ext in [".c", ".cpp"]:
                if os.path.exists(pyx_path[: -len(".pyx")] + ext):
                    c_path = pyx_path[: -len(".pyx")] + ext
                    break
            if c_path == "" or os.path.getmtime(c_path) < os.path.getmtime(pyx_path):
                roots += [root]
                filenames += [fname]

"""
We construct a list of Extension objects from the .pyx files retained
above. If the .pyx file is not associated to a setup_xxx.yml conf file,
no compilation params are given to the extension.
Otherwise, the params of the setup_xxx.yml file are given, with
the potential conversion of the define_macros params into a tuple
(expected by the class Extension).
A list of possible params to write in the setup_xxx.yml file is
yielded here:
https://setuptools.pypa.io/en/latest/userguide/ext_modules.html
sources (the .pyx file and possibly other stuff written in c/cpp
must be included in this params file if you need to specify
other param for compilation.
"""
module_list = []
for i in range(len(filenames)):
    pyx_path = os.path.join(roots[i], filenames[i])
    params = {"sources": [pyx_path]}
    yaml_path = os.path.join(roots[i], "setup_" + filenames[i][: -len(".pyx")] + ".yml")
    ext_name = re.sub(re.escape(os.path.sep), ".", pyx_path[: -len(".pyx")]).lstrip(".").lstrip("/")

    if os.path.exists(yaml_path):
        with open(yaml_path) as f:
            v = yaml.load(f, Loader=yaml.loader.SafeLoader)
        if "cython" in v.keys() and "extension" in v["cython"].keys():
            w = v["cython"]["extension"]
            # tuple required for macros but not implementable in yaml
            if "define_macros" in w.keys():
                z = w["define_macros"]
                if z is not None:
                    for j in range(len(z)):
                        z[j] = (z[j][0], z[j][1])
        params = w

    ext = Extension(name=ext_name, **params)
    module_list.append(ext)

"""
Now we launch the setup. Within this setup we cythonize all modules with the general params
given in the setup_cython.yml conf file
The possible params to set in the conf file are here:
https://cython.readthedocs.io/en/latest/src/userguide/source_files_and_compilation.html
Note that our script implements compiler directives but not Cython compiler options
"""
params = setup_cython["cython_compile"]["cythonize_arguments"]
setup(
    ext_modules=cythonize(module_list=module_list, **params),
    include_dirs=[numpy.get_include()],
    zip_safe=False,
)
