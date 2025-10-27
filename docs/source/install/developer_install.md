(install)=

# Developer Install

:::{important}
The following commands will install *landlab* into your current environment. Although
not necessary, we **highly recommend** you install landlab into its own
{ref}`virtual environment <virtual-environments>`.
:::

If you will be modifying code or contributing new code to *landlab*, you will first
need to get *landlab*'s source code and then install *landlab* from that code.

## Source Install

*Landlab* is actively being developed on GitHub, where the code is freely available.
If you would like to modify or contribute code, you can either clone our
repository

````{tab} ssh
```bash
git clone git@github.com:landlab/landlab.git
```
````

````{tab} https
```bash
git clone https://github.com/landlab/landlab.git
```
````

or download a [zip file](https://github.com/landlab/landlab/archive/refs/heads/master.zip):

```bash
curl -OL https://github.com/landlab/landlab/archive/refs/heads/master.zip
```

Once you have a copy of the source code, you can install it into your current
Python environment by first installing *Landlab* dependencies and then building
and installing *Landlab*.

### Install dependencies

*Landlab*'s dependencies are listed in *requirements.in*.

````{tab} mamba
```bash
cd landlab
mamba install --file=requirements.in -c nodefaults -c conda-forge --override-channels
```
````

````{tab} conda
```bash
cd landlab
conda install --file=requirements.in -c nodefaults -c conda-forge --override-channels
```
````

````{tab} pip
```bash
cd landlab
pip install -r requirements.in
```
````

### Build and install *Landlab*

*Landlab*'s build process includes compiling Python extensions, which requires
you to have a C++ compiler installed. *Linux* will usually already have one,
on *Mac* you can use *XCode*, and on *Windows* you will need to install *MSVC*.
For help on installing *MSVC*, you may want to refer to the *conda-forge* page
on [compiling code on Windows](https://conda-forge.org/docs/maintainer/knowledge_base/#particularities-on-windows)
or the [Python wiki page for Windows compilers](https://wiki.python.org/moin/WindowsCompilers).

If you are using *conda*/*mamba*, set up your compilers to build libraries
compatible with other installed packages,

````{tab} mamba
```bash
mamba install compilers -c nodefaults -c conda-forge --override-channels
```
````

````{tab} conda
```bash
conda install compilers -c nodefaults -c conda-forge --override-channels
```
````

With compilers set up and dependencies installed, build and install *Landlab*,

```bash
pip install -e .
```

## Developer Tools

Once you start developing with *Landlab*, we recommend that you use [nox]  to
automate common tasks such as, for example, running the tests, building the docs, and
finding lint.

```bash
pip install nox
```

The following list shows how to use [nox] for some of the more common tasks:

- Run the tests:

  ```bash
  nox -s test
  ```

- Run the tests on the notebooks:

  ```bash
  nox -s test-notebooks
  ```

- Build the docs:

  ```bash
  nox -s docs-build
  ```

- Run the linters:

  ```bash
  nox -s lint
  ```

- To get a complete list of the available targets:

  ```bash
  nox -l
  ```

[nox]: https://nox.thea.codes/en/stable/
