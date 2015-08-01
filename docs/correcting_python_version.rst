.. _correcting_python_version:

==================================
Dealing with Python on your system
==================================

Known issues, and their symptoms
================================

Based on our experience of guiding users through the install process so far, we have
found that the vast majority of technical issues arise due to conflicts between various
different copies of and versions of Python that tend to accumulate on laptops used
frequently for scientific computing. Very many scientific computing applications are
built onto a Python framework (notably various flavours of GIS software), and in some
cases this can cause your machine to become confused over which version it should be using
at a given time.

Other causes of this problem may occur if you have ever messed about with your `.profile`,
`.bash_profile`, or `.bashrc` files (Mac/Linux) or your Environment Variables (PC).
In such situations, you're likely to need :ref:`The Hard Way <the_hard_way>`.

These kind of issues will also almost always happen if you're not running bash
(the most commonly used UNIX shell). Find out which you're running with
``echo $SHELL`` from the prompt. Anything other than `/bin/bash` being returned
here will give you a bad time with your Python install (both Anaconda and
Canopy assume you're running bash). We're going to assume that if you've modified
your shell, you probably know enough to fix this problem yourself. However, a clue
to get you started is that you'll need to modify the unix path in the profile file
for your shell (...which is named differently according to which shell it is -
get Googling!) so that it also points to your Anaconda binaries. Note the syntax
for the shell commands to do this will probably also be different in each shell!


SYMPTOMS:
>>>>>>>>>

#. Typing ``which python`` and ``which ipython`` gives versions of Python that appear to
   point to different locations, or one or other of thses commands gives an error.
#. The install procedure appeared to go fine, including the test import in the Python
   shell described at the end of the instructions, but now you can't import Landlab from
   the Python prompt provided inside the Canopy or Anaconda app.
#. As above, but you can't import Landlab from inside an iPython shell or from an iPython
   Notebook.
#. Landlab appears to import OK (though there may be warnings or errors at the time of
   import), but when you try to use it you see lots of errors and warnings indicating it
   can't find Python modules it needs (e.g., Numpy, SciPy).

You may have some subset of the above symptoms occurring together, but one is normally
enough to indicate Python conflict is the problem.


Resolving the problem:
>>>>>>>>>>>>>>>>>>>>>>

There are two ways to fix this problem; an easy way, and a hard way. Try the easy one
first!

.. _the_easy_way:

The Easy Way
------------

Both Canopy and Anaconda perform tests to see if you have other versions of Python
already on your system when you install. This is why they ask you if you want to set
Canopy/Anaconda as the default Python package at install time. It's possible you clicked
`No` when you should have clicked `yes`, or for some other reason the software failed to
set the first time.

If you are trying to install Canopy, you're in luck; you can also do this after install.
Go to `Preferences`. Under `General`, you will see text telling you Canopy is not the
default environment, with a button to push saying `Set as default`. Press it.

Unfortunately, in Anaconda, you have to set this up more yourself. See :ref:`The Hard Way
<the_hard_way>`.

.. _the_hard_way:

The Hard Way
------------

If the above didn't work, you're going to have to doctor your system environment variables
yourself. What you do is system-dependent:

**MAC/LINUX**. (This assumes you're running a bash shell.)
Open a terminal. It should automatically open in your home directory, but
just to be on the safe side, type::

    cd ~

Now get a list of all the files in this folder, including the hidden ones::

    ls -a

Have a thorough look at this list. You're looking for files called `.profile` and
`.bash_profile`. Which you have will depend on your system configuration. If you have
*both*, this is probably the root cause of your problem; `.bash_profile` overrides
`.profile` if both are present. In this case, open `.profile` and copy everything you
find inside across to `.bash_profile`. If you've already tried The Easy Way, profile
will probably contain several (repeated?) entries clearly referring to Canopy or "EPD"
having tried to reset itself as default. You can safely delete any repeated text in this
file. (Remember, # is the "comment out" symbol in the Bash terminal). Make sure you
don't actually delete anything else, and we recommend you don't actually remove the
.profile file. Use your favorite Unix editor to make the changes - almost all machines
have `nano`, which is super basic, but works. Once you're done, save (Ctrl-X, say yes to
save in `nano`), then quit and restart any apps that are using Python - including the
terminal app itself. Test if this solved your problem.

If you only have one of `.bash_profile` or `.profile`, open it in an editor and have a
look. You're looking for any obvious references to Python on any of the uncommented-out
lines. You'll probably see several lines that look something like::

    PATH="Users/yourusername/Library/...:${PATH}"
    export PATH

These lines are modifying your environment property PATH, by appending the location of
whatever thing it is the line refers to to the start of PATH, which is a big long
composite string containing lots of addresses. The key here is to understand that *things
nearer the front of the PATH list get precedence*. So if you see lots of `export PATH`
commands in this file that all refer to Python, your system is seeing *only the last
entry*. So find any references to the distribution you want, and duplicate them **at the
end of this file**. Alternatively, you can simply put the comment out symbol # in front
of the lines that are referring to the path you see when you type ``which python``.

In the rare case that you don't see any reference to Anaconda at all (with Canopy, The
Easy Way will have already embedded the necessary lines in this file), you will need to
add the lines yourself. Add the following to the end of the file::

    export PATH=~/anaconda/bin:$PATH

(Assuming you installed Anaconda in the default location).

You can check the path that your machine is actually seeing by typing::

    echo $PATH

from the terminal.


**PC**. On a PC, the same principle of modifying your environment variables applies, but
you access them differently. Go to the Control Panel, then System. On Windows 8, you then
want Advanced System Settings, though this will be similar on older OSes. Go to
Advanced, then to the `Environment Variables...` button. Under User Variables, see if
there is an entry called PATH. If there is, we will modify it. If there isn't, we will
create one. It is **VERY IMPORTANT** that you do not modify any existing text,
*especially* under `System Variables` below.

As is the situation for Mac, above, the system reads these PATH strings from left to
right, and stops once it has found what it is looking for. It also reads User before
System variables. Hence, we want to add new
strings to the left hand (start) of the existing text, if there is any.

First, scan the existing string(s) (including under System) to see if there is any
reference to the Python distribution you are trying to set as default already there.
e.g., my User PATH (running Anaconda cleanly) currently reads::

    C:\Users\Dan\AppData\Local\Continuum\Anaconda;C:\Users\Dan\AppData\Local\Continuum\Anaconda\Scripts

If you find a reference or references like this to the version you're currently trying to
run, copy the text, and add it (repeated) at the start of the User string. Copy this
syntax - semicolons separate paths.

If you can't find any reference to your chosen version (Canopy/Anaconda), you'll need to
add the PATH yourself. For Anaconda, assuming you installed it in the default directory,
add the above string. For Canopy, use the "Set Canopy as default" option ("the easy way"),
which really should work.
See `this page <http://docs.enthought.com/canopy/quick-start/install_windows.html>`_ for
more information on the PATHs used by Canopy if you're still struggling.

Note that modifying the User Variables will only affect the current user account. Add the
text - carefully!! - to the System Variables if you want the changes for all users.

Other issues
============

1. Other install issues often mean that some component of your Python distribution is out
of date. A very common culprit is setuptools, which - extremely frustratingly - isn't
updated by a *conda update -all* call for Anaconda. Other packages can also cause this
kind of problem if out of date. An example of a setuptools related error we've seen
recently ends with::

    error: unknown file type '.pyx' (from 'landlab/components/flexure/cfuncs.pyx')

...combined with warnings referencing a problem with PEP 440.

To our knowledge, this issue only arises for developer installs.

Resolve the issue by updating your distribution. For Anaconda, from a terminal just run::

    > conda update --all
    > conda update setuptools
