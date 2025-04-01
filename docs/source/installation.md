(basic-install)=

# Install

:::{important}
The following commands will install *landlab* into your current environment. Although
not necessary, we **highly recommend** you install landlab into its own
{ref}`virtual environment <virtual-environments>`.
:::

In order to use *landlab* you will first need Python. While not
necessary, we recommend using the
[Anaconda Python distribution](https://www.anaconda.com/download)
as it provides a large number of third-party packages useful for
scientific computing.

To install *Landlab*, simply run the following in your terminal of choice:

````{tab} mamba
```bash
mamba install landlab -c nodefaults -c conda-forge --override-channels
```
````

````{tab} conda
```bash
conda install landlab -c nodefaults -c conda-forge --override-channels
```
````

````{tab} pip
```bash
pip install landlab
```
````

If you would like the very latest development version of *landlab* or want to modify
or contribute code to the *landlab* project, you will need to do a
{ref}`developer installation <install>` of *landlab* from source.
