## macOS and OpenMP

The *macOS* wheels for *Landlab* are built **without OpenMP support**. This
keeps the wheels easy to install and avoids requiring an *OpenMP* runtime at
install time.

If you would like *OpenMP* support on *macOS*, you must install an *OpenMP*
runtime and build *Landlab* from source.

````{tab} conda
```bash
conda create -n omp llvm-openmp
conda activate omp

export OPENMP_PREFIX="$CONDA_PREFIX"

pip install landlab
```
````

````{tab} homebrew
```bash
brew install libomp

export OPENMP_PREFIX="$(brew --prefix libomp)"

pip install landlab
```
````
