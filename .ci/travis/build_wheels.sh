#!/bin/bash
set -e -x

# Install a system package required by our library
# yum install -y atlas-devel

export PATH="/io/anaconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no
conda create -n test_env python
source activate test_env
conda info -a && conda list

pip install numpy
pip install wheel /io/ -w wheelhouse

# Bundle external shared libraries into the wheels
# for whl in wheelhouse/*.whl; do
#     auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
# done

# Install packages and test
# for PYBIN in /opt/python/*/bin/; do
#     "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#     (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
# done
