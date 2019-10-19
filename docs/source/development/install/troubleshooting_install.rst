.. _dev_troubleshooting:

====================================================
Troubleshoot Installation and Development Challenges
====================================================

I updated my working version and now it is broken. What do I do?
----------------------------------------------------------------

One possibility is that the landlab requirements changed between when
you originally installed landlab and when you updated landlab. To
address this, re-run the following lines and then test the installation.

.. code-block:: bash

   $ conda install --yes --file=requirements.txt
   $ python setup.py develop

.. _what_do_if_cant_merge_pr:

What do I do if my pull request cannot be automatically merged?
---------------------------------------------------------------

Get the latest upstream/master and go to the master branch. Remember,
*do not develop here*. Always develop in a feature branch. Merge the
lastest upstream master with your master:

.. code-block:: bash

   $ git fetch upstream
   $ git checkout master
   $ git merge upstream/master

Go to the branch on which you are developing and merge the lastest
upstream master with your branch:

.. code-block:: bash

   $ git checkout <branch_name>
   $ git merge upstream/master

Fix the conflicts. Do this by hand or with a merge editor. This is where
you decide how to integrate the conflicting changes. Since only you know
what and why you made the changes you did, this can only be done by you:

.. code-block:: bash

   $ git merge tool

After everything has been fixed, commit the changes and push the changes
to the repository. The pull request will automatically be updated:

.. code-block:: bash

   $ git commit
   $ git push

Most of these steps have equivalents in the Github app. Use the
"changes" pane to identify where conflicts exist in your version, then
resolve them one by one. When you're done, commit then sync the
un-conflicted version's changes as if they were any other.

I'm seeing errors about Cython when I try to run my code/import Landlab. It used to be fine.
--------------------------------------------------------------------------------------------

Very occasionally, local code updates or rebasing can break the compiled
code that lives in your local developer's install. *Provided you used to
have a fully working Landlab install*, you can fix this by just calling
again from the main Landlab local folder

.. code-block:: bash

   $ python setup.py develop

as described above in the main text. If this is happening when you call
this install function rather than when you try to actually run some
code/import Landlab, see immediately below.

I see errors about Cython when I try to *install* Landlab
---------------------------------------------------------

If you see errors referring to Cython when you try to run
``python setup.py develop``, it indicates you have a problem with your
local compilers. This can happen both the first time you ever try this,
or also subsequently, apparently at random. On a Mac, check first that
you have the free Apple Xcode app (get it from the app store). If you do
have it already, typically this means Xcode has updated itself (this can
happen automatically without your knowledge!) and needs you to
re-authorize its permissions. Open the Xcode app manually, follow the
instructions it will give you, then try the install for Landlab again.
On a PC? Try updating Anaconda.

I'm still confused
------------------

If you are having problems when installing, testing or running Landlab,
please visit our :ref:`Troubleshooting page <troubleshooting>`.

The Landlab development team will be happy to hear from you. Find contact
information :ref:`here <contact>`.
When reporting your problem (in either place) we recommend that you provide a
minimal, complete, and verifiable example which will help the development team
and involved users reproduce your problem and determine a solution.
`This page from Stack Overflow <https://stackoverflow.com/help/minimal-reproducible-example>`_ provides
some background on how to make a minimal, complete, and verifiable example.
