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
`Binder <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/welcome.ipynb>`_.
Here the notebooks are provided within a binder online environment that
includes Landlab.

The welcome page on Binder provides onward links to most of our tutorials.
If you're a newbie you might want to skip directly to a recommended syllabus
for learning Landlab
`here <https://mybinder.org/v2/gh/landlab/landlab/master?filepath=notebooks/tutorials/syllabus.ipynb>`_.

.. _tutorials_EarthscapeHub:

Launch notebooks on EarthscapeHub
---------------------------------

Landlab notebooks can also be run on `EarthscapeHub`_.

If you have a *lab* login for class,
start the Landlab welcome notebook `here`__.
Likewise,
if you have a *jupyter* login associated with your CSDMS membership,
start the welcome notebook `here`__.

To directly access the recommended syllabus for learning Landlab
on *lab*, start `here`__, and for *jupyter*, start `here`__.

Note that software is grouped into kernels on EarthscapeHub.
Landlab is installed in the kernel titled "CSDMS";
select this kernel to run the tutorial notebooks.

.. _EarthscapeHub: https://csdms.colorado.edu/wiki/JupyterHub
.. __: https://lab.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fwelcome.ipynb%3Fautodecode&branch=master
.. __: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Fwelcome.ipynb%3Fautodecode&branch=master
.. __: https://lab.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Ftutorials%2Fsyllabus.ipynb%3Fautodecode&branch=master
.. __: https://jupyter.openearthscape.org/hub/user-redirect/git-pull?repo=https%3A%2F%2Fgithub.com%2Flandlab%2Flandlab&urlpath=lab%2Ftree%2Flandlab%2Fnotebooks%2Ftutorials%2Fsyllabus.ipynb%3Fautodecode&branch=master

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

If you don't have a copy of the source code, you can run the `notebooks.py`_
script to fetch the set of notebooks that matches your version of *landlab*.
This can be done either by running,

.. code-block:: bash

   curl -L https://raw.githubusercontent.com/landlab/landlab/master/notebooks.py | python -

or by downloading the script and running the following from the terminal,

.. code-block:: bash

   python -m notebooks

The download script will create a folder called *landlab-<VERSION>*, where *<VERSION>*
is the version of the notebooks you have requested (e.g. *master* or *2.5.0*). If
you would like to get a specific version of the notebooks, which may not match your
version of *Landlab*, you can specify that as a command line argument. For example,

.. code-block:: bash

   curl -L https://raw.githubusercontent.com/landlab/landlab/master/notebooks.py | python - 2.4.1

.. _notebooks.py: https://github.com/landlab/landlab/blob/master/notebooks.py

Install dependencies
````````````````````

Within the downloaded folder is a file named *requirements-notebooks.txt* that
contains a list of requirements needed to run the notebook tutorials.

.. important::

  The following will install the requirements into your current environment. Although
  not necessary, we **highly recommend** you install these into a separate
  :ref:`virtual environment <virtual_environments>`.

Use *mamba* or *conda* to install the requirements.

.. tab:: mamba

  .. code-block:: bash

     mamba install --file=requirements-notebooks.txt

.. tab:: conda

  .. code-block:: bash

     conda install --file=requirements-notebooks.txt


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
