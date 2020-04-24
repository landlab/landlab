.. _tutorials:

Tutorials
=========

The Landlab Tutorials provide examples of Landlab core concepts and component
introductions. Tutorials exist as interactive Jupyter notebooks that contain
alternating cells of computer code and text that explain the code. In addition
to Landlab Tutorials that exemplify Landlab, notebooks intended to teach and
learn surface dynamics are the
:ref:`Landlab Teaching Tutorials <teaching_tutorials>`.

Launch notebooks online
-----------------------

Landlab Notebooks can be accessed online with the following link:
`Binder <https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/welcome.ipynb>`_.
Here the notebooks are provided within a binder online environment that
includes Landlab.

The welcome page on Binder provides onward links to most of our tutorials.
If you're a newbie you might want to skip directly to a recommended syllabus
for learning Landlab
`here <https://mybinder.org/v2/gh/landlab/landlab/release?filepath=notebooks/tutorials/syllabus.ipynb>`_.

Launch notebooks locally
------------------------

How you installed Landlab determines the steps to launch notebooks that use a
copy of Landlab on your computer.

User installations
``````````````````

If you installed using a prepackaged binary, the most common method, the
notebooks can be run in a conda environment.
:ref:`These instructions <conda_environment>` describe how to create a conda
environment for the notebooks.

Once the conda environment has been created, it must be activated and then the
notebooks can be launched:

.. code-block:: bash

   $ conda activate landlab_notebooks
   $ jupyter notebook notebooks/welcome.ipynb

Developer installations
```````````````````````

If you installed using the
:ref:`developer install instructions <developer_install>`, the tutorials exist
on your computer within the ``notebooks`` folder of your local Landlab clone.

If you installed Landlab within an environment, the environment must be
activated and then the notebooks can be launched. The following code assumes
that the environment name is `landlab_dev`:

.. code-block:: bash

   $ conda activate landlab_dev
   $ jupyter notebook notebooks/welcome.ipynb

Contributing new or modified tutorials
--------------------------------------

If you write a Landlab Tutorial or Gist, please contribute it via a pull request
to the master branch of the Landlab repository. See this
:ref:`page <ongoing_development>` about contributing to Landlab, and
:ref:`reach out for help <contact>` when needed.

Landlab clinics and workshops
-----------------------------

For more examples and tutorials, see also our :ref:`Clinics & workshops
page <clinics_workshops>`.
