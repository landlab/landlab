if [[ "$BUILD_STR" == "" ]]; then
  echo "Building release version"
  export CHANNEL="main"
else
  echo "Building dev version"
  export CHANNEL="dev"
fi

echo "Building package for Python $PYTHON_VERSION, numpy $NUMPY_VERSION"
conda config --set anaconda_upload no

echo "Uploading to $CHANNEL"
anaconda -t $ANACONDA_TOKEN upload --force --user landlab --channel $CHANNEL $HOME/miniconda/conda-bld/**/landlab*bz2

echo "Done."
