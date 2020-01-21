.. _organization:

Package Organization
--------------------

Below is a tree view description of the key files and directories of the
Landlab package. Note that not all files are shown and that top level files
have been reorganized out of alphabetical order for clarity. Directories are
shown down to two levels.

::

      landlab
      |
      | ## Top Level Files ##
      |
      | # Package Description
      |
      ├── README.rst
      ├── CHANGELOG.md
      ├── CONTRIBUTING.md
      ├── LICENSE.txt
      ├── MANIFEST.in
      |
      | # Dependency Requirements
      |
      ├── requirements.txt
      ├── requirements-dev.txt
      ├── requirements-notebooks.txt
      ├── requirements-testing.txt
      |
      | # Compilation
      |
      ├── setup.cfg
      ├── setup.py
      ├── Makefile
      |
      | # Environment Specifications
      |
      ├── environment-dev.yml # Developer conda environment
      ├── environment.yml # User/Notebook conda environment
      ├── readthedocs.yml # Documentation buidling conda environment
      |
      | # Continuous Integration
      |
      ├── appveyor.yml # Specification for Appveyor Continuous Integration
      ├── conftest.py
      |
      |
      | ## Directories ##
      |
      ├── docs # The documentation
      │   ├── scipy-sphinx-theme
      │   └── source
      |
      ├── joss # Files associated with JOSS publications
      |
      ├── landlab # The core Landlab package
      │   ├── bmi
      │   ├── ca
      │   ├── cmd
      │   ├── components
      │   ├── core
      │   ├── data_record
      │   ├── field
      │   ├── framework
      │   ├── graph
      │   ├── grid
      │   ├── io
      │   ├── layers
      │   ├── plot
      │   ├── testing
      │   ├── utils
      │   └── values
      |
      ├── notebooks # Jupyter notebooks
      │   ├── teaching
      │   └── tutorials
      |
      ├── scripts
      |
      └── tests # The unit tests, structured to mirror landlab/landab
      |   ├── ca
      |   ├── components
      |   ├── core
      |   ├── data_record
      |   ├── field
      |   ├── graph
      |   ├── grid
      |   ├── io
      |   ├── layers
      |   ├── plot
      |   ├── utils
      |   └── values
      ├── END
