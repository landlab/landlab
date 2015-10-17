#! /bin/bash

PYTHON=${PYTHON:-python}

run_test()
{
  mkdir -p _test
  cd _test

  INSTALLDIR=$($PYTHON -c "import os; import landlab; print(os.path.dirname(landlab.__file__))")

  coverage run --source=$INSTALLDIR ../scripts/test-installed-landlab.py && (coverage report && cp .coverage ..)
}

run_test
