if [[ "$TRAVIS_TAG" == v* ]]; then
  pip install twine wheel
  if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
    python setup.py sdist
  else
    python setup.py bdist_wheel
  fi
  twine upload -u mcflugen -p$PYPI_PASS dist/*
fi
