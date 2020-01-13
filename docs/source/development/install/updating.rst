.. _updating:

====================
Keep Landlab Updated
====================


To take advantage of new features and new library additions, we recommend you
**update Landlab** fairly frequently. From time to time you should fetch
commits to the "trunk" version of Landlab that you don't have in your forked
working copy. This is also important to keep your branch up to date with
improvements/bug elimination happening as part of core Landlab development.

Before you begin, we recommend you review this `GitHub page describing syncing
a fork <https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/syncing-a-fork>`_.

Through the command line, you should be able to update your Landlab repository
from master using the following command

.. code-block:: bash

   $ git fetch upstream

A common reason why this might not work is if you haven't set the remotes
correctly. This
`GitHub page <https://help.github.com/en/github/using-git/managing-remote-repositories>`_ on
managing remotes might be helpful.

If we have added a dependency, you may see an import error. Similarly, if we
have added a dependency, you may see some tests break. A likely fix is to
ensure you have all the requirements specified in the Landlab development
environment described by the `environment-dev.yml` file

.. code-block:: bash

   $ conda env update landlab_dev  --file environment-dev.yml
   $ conda activate landlab_dev
   $ pip install -e .

The first of these installs any new requirements for development. The second
line activates the environment (called "landlab_dev").  The third will verify
that all requirements are met and recompile Landlab.

If the Cython code has changed since you last updated, you will probably see
errors associated with files named `cfuncs.pyx` (or similar). These changes
require that you recompile your code. This can be done by executing the
following from within the top level landlab folder (which contains the
`setup.py` file)

.. code-block:: bash

   $ pip install -e .

If none of these options work, please make an issue.

Additional issues
-----------------
If you have been developing code and making commits to the ``master`` branch on
your fork, this may not work. It is for this reason that we highly recommend
you develop on a development branch.

Still having problems? This probably means that some time early in our
development cycle you installed Landlab with one of our old procedures. The clue
will be that you still have a (very out of date!) copy of the Landlab code
base somewhere on your machine. Another possibility is that you've previously
tried a
:ref:`developer install <developer_install>`.
This procedure will also work in this case.

Try this:

In a terminal, navigate to the top level directory of
that old code, the one that contains the file *setup.py*.
This is likely to be *your_home_dir*/landlab, if you installed with git
and left all the defaults as is.
Then

.. code-block:: bash

    $ pip uninstall landlab #just to be on the safe side, may get errors again
    $ python setup.py develop -u

This should remove the install, **if** you installed as a developer.

Still getting error messages? This means we're going to have to excise the
old Landlab install "by hand". You're looking to remove any reference to
Landlab that lives inside *your_python_install*/lib/python2.7/site-packages.
**Do this only after you've exhausted other possibilities, above**, as
packages like pip will get annoyed with you if you start manually deleting
their files if they installed them in the first place. To minimize the risk,
once again make sure you have just run

.. code-block:: bash

    $ pip uninstall landlab

Then find your Python directory with

.. code-block:: bash

    $ which python

Find that folder, ignoring everything after and including the subfolder
*bin*. Instead, go to *your_install*/lib/python2.7/site-packages. In here,
you should find one (or more) folders referring to Landlab, e.g.,
*landlab* or *landlab.egg-link*, or some other reference to
*landlab.egg*. Delete these. Leave everything else as it is!

If you are running an Anaconda distribution, we now recommend you replace your
pip install with a conda install. Simply do this

.. code-block:: bash

     $ pip uninstall landlab
     $ conda install landlab -c landlab

You will then be able to update Landlab along with the rest of your conda
packages

.. code-block:: bash

     $ conda update --all -c landlab

If you prefer to remain with pip, try another pip install

.. code-block:: bash

     $ pip install landlab

This should now take.

*Still* having problems? There are probably multiple
versions of Python on your machine interfering with each other. Solve
that problem first (see
:ref:`here <correcting_install_paths>` for
some help), then return to trying to install Landlab.

**Without using command line:**

You can't do the equivalent in the GitHub app, but you can do it through the
website—though it is a bit more involved. Navigate to the
`master Landlab fork <https://github.com/landlab/landlab>`_. Hit the big green
"new pull request" button. Hit the "compare across forks" hyperlink. You now
want to set the *base fork* dropdown menu to your local Landlab
(`UserName/landlab`), and leave the *head fork* as it is (`landlab/landlab`).
Leave both the *base* and *compare* boxes set to `master`, as they should
already be. The site will then have a think, then report back to you whether it
can perform an automatic merge (in green or yellow, just below the dropdown
boxes), and also what the differences between these two versions are (scroll
down).

Take a quick look at the commits list to reassure yourself you got the "base"
and "head" versions the right way round (DH struggles with this every time…).
Presuming you got the green "able to merge" dialogue, now just fill in the pull
request box with a title (e.g., "updating local Landlab fork from master") and
hit "Create pull request." This will take you to a page showing the open
request. Confirm it and actually perform the code update in your fork by
scrolling down and hitting the green "Merge pull request" button. You'll then
need to sync your local version with the fork that just updated on the Github
servers—easiest through the app, using the "Sync" button in the top right (make
sure you're in the local master branch).

If you didn't get the green "able to merge" button, it means some of your local
changes conflict with what's in the main Landlab repo. Don't panic! This is
fine. Just follow the instructions to press ahead. Conflicts can be resolved
manually later in the process—conflicting sections of the code will be tagged
with a distinctive >>>>>>> symbology, and you can find and modify them by hand.
See :ref:`here<what_do_if_cant_merge_pr>`
for more documentation on this.
