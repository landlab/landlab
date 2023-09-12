.. _landlab_develop_with_git:

================
Develop with Git
================

*This information is intended to help users who have downloaded the
Landlab source code from Github, and intend to develop new components for
Landlab, or otherwise make changes to the code base. Here we illustrate
a few additional code management techniques using git and Github, which
we hope will help you develop for Landlab in a more efficient, faster,
more organised way.*

If you've installed Landlab from the source code, you've already
interacted with the git code management and version control systems, and
with Github, which is built onto the basic git framework. However, this
same framework offers you powerful ways to organise and track the change
you make to the code after install, and to plug any changes you make
back into the main code base for Landlab.

Your **fork** of Landlab is your starting point for all modifications to
the code. However, you can create *branches* from your master fork in
which you will actually modify code. The great advantage of this system
is that you can have several different branches running in parallel,
each aiming to complete a particular objective. The master version of
the code is not changed until you *merge* your completed branch back
into the master. If you decide that changes you're making are actually
unhelpful or redundant, you can easily roll them back - either by
individually reverting each committed change, or by deleting the whole
branch.

`This <http://rogerdudler.github.io/git-guide/>`_ is a nice, simple
guide to the git command line tools. All this functionality is also
available through the Github app.

*Note: Github provides a nice GUI interface for using git called "Github
Desktop". This provides all the functionality described below in a form
that avoids the command line, if that sounds like something you'd
prefer.*

Branching with git
------------------

Making a new branch
```````````````````

Before making any changes to your code, you should create a new branch.

Update your mirror with any upstream changes you don't have:

.. code-block:: bash

   $ git fetch upstream

If using git at the command line, make the new branch like this:

.. code-block:: bash

   $ git branch name-of-branch upstream/master
   $ git checkout name-of-branch

From the Github app, you can do the same thing with the "branching"
button in the top line. A nice feature of the app is it shows you
graphically how your active branch currently differs from the master
branch - when you branched it, and when changes were committed to both
branches.

You will probably want to choose a descriptive name for your new branch
so that you and others will remember what it is you are intending to do
with your branch (for example, bugfix-for-that-major-problem, or
add-that-cool-feature).

You can create branches from branches!

If changes appear in the master and you want to also have them in your
active branch, you can update your active branch from the master with

.. code-block:: bash

   $ git pull

or with the "Update from master" button in the app.

Pushing changes from your local machine to your fork on Github
``````````````````````````````````````````````````````````````

If you want to keep a copy of the files you have modified or created on
your branches on your public GitHub page for Landlab (you probably do as
this will serve as a file backup) you need to tell git to push changes
to your github repo. This is done with the following command:

.. code-block:: bash

   $ git push --set-upstream origin name-of-branch

In the app, the same functionality is achieved by first "publishing"
your branch (creating it on your page within the github.com central
repository), then using the "sync" button to send the files to the
central github.com servers.

On your Landlab GitHub page you will now be able to toggle between your
various branches to see the code you have committed. The app also lets
you see the structure of your branches on your local machine.

Committing changes and merging branches back in
```````````````````````````````````````````````

Changes you make to your code are "saved" in git when you commit them to
your branch. Save your files, then at the command line

.. code-block:: bash

   $ git commit -m "Text describing the changes"

Again, the app provides the same functionality, but with the added bonus
that it shows you what the changes you've made since your last commit
actually are. Click on the "X Uncommitted Changes" button at the top
centre, and see which files have changed and what's happened. Pick the
files for which you want to store the changes as part of this commit,
type text describing the change in the boxes, then hit "Commit to
my-branch-name". If you go to History, you can see a record of all past
changes in the branch. You can then use git to travel "back in time" and
review what the code was like at any time in the past!

Note that files that you do not explicitly ask git to track (either by
clicking the checkbox next to the file in the GUI or by using
``git add``) are not tracked — and thus not sent to the github.com
central servers when you push changes.

Once you're happy with your branch, and the code is fully functional
again, it's time to merge it back into the master. This procedure
generally works best if you first pull any changes from the main Landlab
master branch (not just the master branch on your fork) into your active
branch, and resolve any conflicts there (so you don't mess up the
master).
Once you've done that pull, in git at the command line, make the master
your active branch again then *merge* the branch:

.. code-block:: bash

   $ git checkout master
   $ git merge my-branch-name

In the app, create the merge by making a "pull request" using the button
in the top right. The process is fairly self explanatory, and provides a
preview of whether there will be any conflicts. Once you've created the
merge, click through the hyperlink and merge it into the master on the
website using the prominent green button. On your local machine, sync
your master branch to pick up the changes locally.

Pulling changes from your fork to the Landlab master fork
---------------------------------------------------------

Once you've completed whatever modifications you were working on with
Landlab, we'd like to incorporate your changes back into the main code
of Landlab so everyone can benefit from your enhancements. This is done
by creating a *pull request* from your fork into the Landlab master
fork. This is basically the inverse process you use to update your fork
from the master fork
(but in this case, one of us will review your changes before it gets
merged in).

Perform this procedure through the Github website. Go to the github page
for your fork, and click the green "New pull request" button at the top.
The next page shows you which branch on which fork (the "head") will be
merged into which other branch and fork (the "base"). These details
should all be correct as shown. There may well be conflicts reported on
this page. If there are, consider
updating your fork from the master fork
before finalising the request. Once you're ready to go, click the next
"Create pull request" green button. You'll be redirected to a discussion
page for your request, and it will be visible to all of the admins for
the main Landlab fork - one of whom will review your changes and
actually make the merge (*please don't do this yourself!*).

If you're confused by this process, just create the request, and one of
us will see it and come to help you. You can create comments on your
request from the website at any time after you've made it.


Troubleshooting
---------------

What do I do if my pull request cannot be automatically merged?
```````````````````````````````````````````````````````````````

Get the latest upstream/master and go to the `master` branch. Remember,
*do not develop here*.  Always develop in a feature branch. Merge the lastest
upstream master with your master

.. code-block:: bash

  $ git fetch upstream
  $ git checkout master
  $ git merge upstream/master

Go to the branch on which you are developing and merge the lastest upstream
master with your branch

.. code-block:: bash

  $ git checkout <branch_name>
  $ git merge upstream/master

Fix the conflicts. Do this by hand or with a merge editor. This is where you
decide how to integrate the conflicting changes. Since only you know what and
why you made the changes you did, this can only be done by you

.. code-block:: bash


  $ git mergetool

After everything has been fixed, commit the changes and push the changes to
the repository.  The pull request will automatically be updated

.. code-block:: bash

  $ git commit
  $ git push
