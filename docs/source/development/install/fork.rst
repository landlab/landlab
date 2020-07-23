.. _dev_install_fork:

=========================
Create of Fork of Landlab
=========================

Getting the code from Github
----------------------------

Landlab development takes place in your own *fork* of the main Landlab
repository. A fork is a *mirror* of the repository and is hosted on your
personal GitHub account. You will use this fork for developing new
Landlab features. Your changes will migrate to the core repository (for
review and merging) by requesting that the main repository "pull" in
your changes. This is known as a pull request and is facilitated through
the GitHub website.

You will only need to do this once for each project to which you want to
contribute. Github has some great documentation on `how to create a
fork <https://help.github.com/en/github/getting-started-with-github/fork-a-repo>`_. We outline
below the basic steps as applied to Landlab.

Create a Github account
```````````````````````

1. You can create a GitHub account by going to the `GitHub
   website <https://github.com>`_.
2. If you haven't already, `Install
   Git <https://help.github.com/en/github/getting-started-with-github/set-up-git>`_. Note that if
   you are using a mac with OS lower than 10.7, we have found it
   difficult to setup git, and you may want to upgrade your OS.
   Alternatively, you can do this in a more user-friendly fashion by
   installing the `Github graphical user
   interface <https://desktop.github.com>`_ for your OS, which comes
   packaged with git for command line by default, and also lets you
   perform many git functions (updating your code, pushing changes back
   to github…) without having to deal with the command line interface.
3. Configure your account to allow write access. If you choose to
   install the git GUI, then it will set-up an SSH key for you. If you
   install on the command line, you might need some help with this, see
   `Generating SSH Keys on
   GitHub <https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh>`_.
4. If you want to ensure the ssh key is set correctly, go to your home
   page on github and hit the account settings (wrench and screwdriver
   button in upper right corner). On this page hit the SSH keys tab on
   the left. This should show that you have a key for whatever computer
   you are currently working on. Note that you may have more than one
   key if you have installed git on more than one computer with the same
   user account.

Creating your own fork of Landlab
`````````````````````````````````

The following steps will create a fork of the Landlab repository under
your GitHub account.

1. Sign in to your GitHub account.
2. Go to the `Landlab home page <https://github.com/landlab/landlab>`_
   on GitHub.
3. Click the fork button in the upper-right corner of the page.

Once completed, you will be redirected to the home page for your own
copy of Landlab (github.com/UserName/landlab).

Cloning your fork to your computer
``````````````````````````````````

You can clone the fork (that lives on the Github website) locally to
your computer either using the Github app (the GUI), or directly from
the command line using git. If you've never used git before, the app is
probably the way to go.

Using the GUI/app
~~~~~~~~~~~~~~~~~

1. Sign in to github on the GUI.
2. Hit on your account on the left side of the GUI.
3. This will show you your fork of Landlab. Hit the clone option next to
   the fork. This will download the Landlab package to your local
   computer. If you are on a windows machine, this will put Landlab in
   the Documents/GitHub folder. If you are on a mac, you are given the
   option of where to put the downloaded Landlab package.

Alternatively, you can also clone your fork in the GUI directly through
the website. Navigate to the page for your fork on the web
(UserName/landlab) and hit the "Clone in Desktop" button
(monitor-with-downward-arrow) on the top/right hand side of the page.
Make sure you have the GUI installed and set up on your machine before
you try this for the most pain-free results. Also ensure you clone your
fork, not the main development version of Landlab (i.e., the version at
github.com/your-user-name/landlab, not the version at
github.com/landlab/landlab).

Using the command line
~~~~~~~~~~~~~~~~~~~~~~

Use the following commands from a terminal of your choice.

.. code-block:: bash

   $ git clone git@github.com:your-user-name/landlab.git
   $ cd landlab
   $ git remote add upstream git://github.com/landlab/landlab.git
