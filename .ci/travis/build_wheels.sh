#!/bin/bash
set -e -x

PYBIN=" \
  /opt/python/cp36-cp36m/bin \
  /opt/python/cp37-cp37m/bin \
"

for bindir in $PYBIN; do
  ${bindir}/pip install numpy
  ${bindir}/pip wheel /io/ -w /io/wheelhouse
done

# Bundle external shared libraries into the wheels
for whl in /io/wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

ls /io/wheelhouse
