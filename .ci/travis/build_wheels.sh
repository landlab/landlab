#!/bin/bash
set -e -x

# Install a system package required by our library
# yum install -y atlas-devel

# for PYBIN in /opt/python/*3*/bin/; do
#   ${PYBIN}/pip install numpy
#   ${PYBIN}/pip wheel /io/ -w wheelhouse
# done

/opt/python/cp36-cp36m/bin/pip install numpy
/opt/python/cp36-cp36m/bin/pip wheel /io/ -w /io/wheelhouse

/opt/python/cp37-cp37m/bin/pip install numpy
/opt/python/cp37-cp37m/bin/pip wheel /io/ -w /io/wheelhouse

ls /io/wheelhouse

# Bundle external shared libraries into the wheels
# for whl in wheelhouse/*.whl; do
#     auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
# done

# Install packages and test
# for PYBIN in /opt/python/*/bin/; do
#     "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
#     (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
# done
