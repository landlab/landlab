if [[ "$TRAVIS_TAG" == v* ]]; then
  pip install twine wheel
  python setup.py bdist_wheel
  twine upload -u mcflugen -p$PYPI_PASS dist/*
fi
