if [[ "$TRAVIS_TAG" == v* ]]; then
  export CHANNEL="main"
else
  export CHANNEL="dev"
fi

echo "Uploading to $CHANNEL"
anaconda -t $ANACONDA_TOKEN upload --force --user landlab --channel $CHANNEL $HOME/miniconda/conda-bld/**/landlab*bz2

echo "Done."
