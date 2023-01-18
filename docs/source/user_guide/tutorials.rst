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

To run the tutorials locally, you will first need to install *landlab*
on your computer. If you have not already done so, please see the *landlab*
:ref:`installation guide<install>`.
Because several of the notebooks depend on packages that cannot be
installed using *pip*, we recommend you use *conda* to install the
tutorial notebook requirements.

Get the notebooks
`````````````````

If you have the *landlab* source code, you already have the notebooks (they are
in the *notebooks/* folder).

If you don't have a copy of the source code, you can run the
`notebooks.py`_
script to fetch the set of notebooks that matches your version of *landlab*.
This can be done either by running,

.. code-block:: bash

   $ curl -L https://raw.githubusercontent.com/landlab/landlab/master/notebooks.py | python -

or by downloading the script and running the following from the terminal,

.. code-block:: bash

   $ python -m notebooks

.. _notebooks.py: https://github.com/landlab/landlab/blob/master/notebooks.py

Install dependencies
````````````````````

The dependencies required to run the notebooks are listed in the file, *requirements-notebooks.txt*
and can be installed with *conda*,

.. code-block:: bash

   $ conda install --file=requirements-notebooks.txt


Run the tutorials
`````````````````

It's now time to run the tutorials. Opening the welcome page is a good place to start
but you can also open individual notebooks as well.

.. code-block:: bash

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
