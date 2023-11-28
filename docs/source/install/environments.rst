.. _virtual_environments:

====================
Virtual Environments
====================

A virtual environment is a self-contained directory tree that contains a Python installation for a particular
version of Python along with additional packages. It solves the problem of one application's
package requirements conflicting with another's.

Two popular tools used for creating virtual environments are the built-in *venv* module and *conda*.
For virtual environments created using *conda*, you can use either *conda* or *pip* to install additional
packages, while *venv*-created environments should stick with *pip*.

.. tab:: conda

    .. code-block:: bash

        conda create -n landlab
        conda activate landlab

.. tab:: venv

    .. code-block:: bash

        python -m venv .venv
        source .venv/bin/activate

Note that you will need to activate this environment every time you want to use it in a new shell.

Helpful links on managing virtual environments:

* `conda environments <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_.
* `venv environments <https://docs.python.org/3/tutorial/venv.html>`_.
