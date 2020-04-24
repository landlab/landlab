.. _troubleshooting:

Troubleshooting Installation/Testing Problems
=============================================

First, make sure you:
---------------------

-  Follow the installation instructions that match your use of Landlab:

   -  :ref:`Install Landlab using Anaconda <anaconda_install>`
      if you are just going to **use** Landlab
   -  :ref:`Install from source code <developer_install>`
      if you are going to **modify** Landlab's source code

-  Follow the instructions relating to the system you are using
   (Mac/Linux or Windows)
-  Install the correct dependencies if you are doing the Developer
   install: they are listed in the `environment_dev.yml` file and you should
   point to this file when :ref:`installing <developer_install>`)
-  Make Anaconda your default Python (when installing Anaconda, or
   :ref:`edit paths afterwards <correcting_install_paths>`).

Try the following commands

On Mac:

.. code-block:: bash

  $ which python
  $ which ipython

On Windows:

.. code-block:: bash

  $ where python
  $ where ipython

python and ipython should show the same path and clearly refer to
Anaconda, looking something like ``/anaconda/bin/python`` \* Update:
Python, Landlab...


MOST COMMON ISSUES:
-------------------

1. Landlab install fails
~~~~~~~~~~~~~~~~~~~~~~~~

Specifications / dependencies conflicts:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Symptoms: when launching the install an error message shows something
   like

``Solving environment: failed``

``UnsatisfiableError: The following specifications were found to be in conflict:...``

-  Possible fix: create and work in a `conda
   environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_

-  Related issues:
   `#625 <https://github.com/landlab/landlab/issues/625>`_,
   `#631 <https://github.com/landlab/landlab/issues/631>`_,
   `#635 <https://github.com/landlab/landlab/issues/635>`_

Issues with PATH configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Symptoms: command not found, conda not found, landlab not found,…
-  Fix: :ref:`Correcting your install
   paths <correcting_install_paths>`
-  Related issues:
   `#537 <https://github.com/landlab/landlab/issues/537>`_ (Windows),
   `#538 <https://github.com/landlab/landlab/issues/538>`_ (Mac),
   `#589 <https://github.com/landlab/landlab/issues/589>`_ (Windows)

Developer version not installing on Windows with Python 3.6.5
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  Symptoms: Installation fails and requires Visual C++ tools to be
   installed
-  Fix: see Issues
   `#768 <https://github.com/landlab/landlab/issues/768>`_\ and
   `#773 <https://github.com/landlab/landlab/issues/773>`_

Landlab for Linux 32-bit (Issue `#617 <https://github.com/landlab/landlab/issues/617>`_):
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sorry, Landlab is not available (yet) for Linux 32-bit!

2. Landlab runs fail
~~~~~~~~~~~~~~~~~~~~

To test your install, use a command prompt, navigate to your landlab
directory and type:

.. code-block:: bash

  $ conda install pytest
  $ pytest

This tests whatever version of landlab you have installed. If this
produces errors for which you cannot find a fix in this page, please
contact the development team by `creating a GitHub
issue <https://github.com/landlab/landlab/issues>`_. Provide the
results of the above test and all the information you can (operating
system, type of install you have done, etc.)

AttributeError
^^^^^^^^^^^^^^

-  Symptoms: Test fails with:

.. code-block:: python

  AttributeError: 'LandlabTester' object has no attribute 'check_fpu_mode'

-  This is a known error, probably due to numpy version 1.14., to check
   your version:

.. code-block:: python

  import numpy
  numpy.__version__

Although tests currently fail, you should be able to run Landlab
normally, you can ignore this error. \* Similar issues:
`#616 <https://github.com/landlab/landlab/issues/616>`_,
`#625 <https://github.com/landlab/landlab/issues/625>`_,
`#630 <https://github.com/landlab/landlab/issues/630>`_

3. Landlab + Spyder problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spyder does not launch
^^^^^^^^^^^^^^^^^^^^^^

-  Symptoms: you cannot launch Spyder or it crashes whenever you try to
   open a file
-  Fix: restart, upgrade, reset… see the `Spyder Troubleshooting
   page <https://github.com/spyder-ide/spyder/wiki/Troubleshooting-Guide-and-FAQ>`_.
   Some security/antimalware software have been known to block Spyder
   from launching: try temporarily disabling it or changing your
   settings. On Windows, some launching problems may be due to conflicts
   with other programs that installed the Qt libraries in your computer
   (you can identify them searching for files with the following
   pattern: ``Qt*.dll``): you can either uninstall these programs or
   ``pip install pyqt5``, which should fix the Spyder launching problem.

Import problem
^^^^^^^^^^^^^^

-  Symptoms: command not found, conda not found, landlab not found,…
-  Fix: likely a path issue, so try :ref:`Correcting your install
   paths <correcting_install_paths>`
-  Related Issues:
   `#538 <https://github.com/landlab/landlab/issues/538>`_

4. Landlab + Jupyter Notebook problems
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  Related issues:
   `#540 <https://github.com/landlab/landlab/issues/540>`_

5. Tutorial fails
~~~~~~~~~~~~~~~~~

Import error
^^^^^^^^^^^^

-  Symptom: ``Cannot import name 'ModuleName'`` This is because a module
   used by the tutorial is not included in the latest release: do the
   :ref:`developer install <developer_install>`
   or wait for the next `release <https://github.com/landlab/landlab/releases>`_.

Other known pb:
~~~~~~~~~~~~~~~

-  **Cygwin** may create some problems: see issue
   `#349 <https://github.com/landlab/landlab/issues/349>`_
-  **GRASS GIS** and Landlab might not like each other too much:
   :ref:`Install GRASS GIS after installing Landlab <grass>`


If this does not solve your issue:
----------------------------------

-  Try a clean install: uninstall and reinstall conda, etc.
-  Use a `conda
   environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_
-  Contact the Landlab development team: `create an issue on
   GitHub <https://github.com/landlab/landlab/issues>`_. Please
   provide all the information you can: the system you are operating on,
   the install you have done, the command that produced the error, what
   the error message is… We will get back to you quickly!
