#! /bin/bash
# Build landlab wheels using Anaconda.
#
# This script must be run in the top-level landlab folder (the one with
# requirements.txt). To build a wheel and upload it to PyPI,
#
#     $ dist_to_pypi.sh PYTHON_VERSION [...]
#
# For example,
#
#     $ dist_to_pypi.sh 2.6 2.7
#
# will build wheels for Python 2.6 and 2.7 and upload them to PyPI.

if [ "$#" == "0" ]; then
  echo "Usage: dist_to_pypi.sh PYTHON_VERSION [...]"
  exit 1
fi

PYTHON_VERSIONS=$@
for v in $@; do
  shift
done

for v in $PYTHON_VERSIONS; do
  echo "Building a wheel for py$v using $(which python)"

  conda create -n dist --file=requirements.txt python=$v
  source activate dist

  if [ "$v" == "2.6" ]; then
    conda install scipy=0.14
  fi

  pip install wheel

  python setup.py develop
  python setup.py bdist_wheel

  source deactivate
  conda env remove -n dist
done

conda create -n dist python=2.7
source activate dist
conda install numpy

python setup.py sdist

pip install delocate twine
delocate-wheel dist/*whl
twine upload dist/*

source deactivate
conda env remove -n dist
