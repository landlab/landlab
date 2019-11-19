.. _correcting_install_paths:

========================
Correcting Install Paths
========================

Based on our experience of guiding users through the install process so far, we
have found that the vast majority of technical issues arise due to incorrect
paths to your Python distribution. Basically, your computer does not know where
to look for the correct file or program and you need to tell it where to look
using your system's `$PATH` variable.

Symptoms
--------

1. Typing ``conda`` after an Anaconda install gives ``command not found``.
2. Typing ``which python`` and ``which ipython``
   point to different locations, the files can't be found, or one or other of
   these commands gives an error.
3. The install procedure appeared to go fine, including the test import in the
   Python shell
   :ref:`described at the end of the instructions <test_landlab_install>`,
   but now you can't import Landlab from the Python prompt provided inside Anaconda.
4. As above, but you can't import Landlab from inside an iPython shell or from
   an iPython Notebook.
5. Landlab appears to import OK (though there may be warnings or errors at the
   time of import), but when you try to use it you see lots of errors and
   warnings indicating it can't find Python modules it needs (e.g., NumPy,
   SciPy).

You may have some subset of the above symptoms occurring together, but one is
normally enough to indicate a Python conflict is the problem.

The Solution
------------

Path problems can be solved by editing your environmental variables. How you do
this is system dependent. A rough outline for how to do this on MAC/LINUX or
PC follows.

.. _the_hard_way:

MAC/LINUX
---------

To set your path variables in a MAC/LINUX system, you will need to manually
edit the ``.profile`` or ``.rc`` files for your terminal shell.

The most commonly used shell is ``bash``. Both Anaconda and Canopy assume
you're running the Bash shell and will put the path to your install there.

To find out which shell you're running, open a terminal window and type
``echo $SHELL`` from the prompt. Anything other than ``/bin/bash`` being
returned will require a little more manual set up as described below.

Now that you know which shell you are using, find out where your computer is
looking for your Python files. From the prompt, type ``echo $PATH``.

You will see your PATH, a big long composite string containing lots of
addresses that your computer searches—in order—to find the requested file. If
you ask for ``python`` or ``conda``, for example, it will return the first
``python`` or ``conda`` it finds in its path.

The key here is to understand that *things nearer the front of the PATH get
precedence*.

**Your goal** is to simply have your *new* distribution appear near the
beginning of the ``PATH`` that results when you type ``echo $PATH``.
And…that's not so hard!

Here is how to edit your files to give you the correct path:

Go to your home directory by typing ``cd ~`` from the prompt.

Now get a list of all the files in this folder, including the hidden ones::

    ls -la

Have a thorough look at this list.

If you are using Bash, you're looking for files called ``.profile`` and
``.bash_profile``. Which you have will depend on your system configuration. If
you have *both*, this is probably the root cause of your problem;
``.bash_profile`` overrides ``.profile`` if both are present. In this case,
open ``.profile`` and copy everything you find inside across to
``.bash_profile``.

For example, your profile likely contains several (possibly repeated) entries
clearly referring to Canopy or "EPD" or Anaconda. You can safely comment out
any old installs or repeated text in this file. (Remember, # is the "comment
out" symbol in the Bash terminal). Make sure you don't actually delete anything
else, and we recommend you don't actually remove the ``.profile`` file.

Use your favorite Unix editor to make the changes—almost all machines have
``nano``, which is super basic, but works. Once you're done, save
(Ctrl-X, say yes to save in ``nano``), then quit and restart any apps that are
using Python—including the terminal app itself.
:ref:`Test <test_landlab_install>` if this solved your problem.

If you only have one of ``.bash_profile`` or ``.profile``, open it in an editor
and have a look. You're looking for any obvious references to Python or your
distribution (Anaconda, Canopy) on any of the uncommented-out lines. You'll
probably see several lines that look something like::

    PATH="Users/<your user name>/Library/...:${PATH}"
    export PATH

These lines are modifying your environment property PATH, by prepending the new
path, ``Users/<your user name> /Library/...``, to the old path, ``${PATH}``.

Lines closer to the end of the file can prepend folders closer to the beginning
of your PATH.

In the rare case that you don't see any reference to Anaconda at all, you will
need to add the lines yourself. Add, for example, the following to the file::

    export PATH=~/anaconda/bin:$PATH

As of this writing (11/2016), the default path for an Anaconda installation on
OSX was in the User directory, which means the path command should look
something like this::

    export PATH="/Users/<your user name here>/anaconda/bin:$PATH"

If you are using a shell other than Bash, you can likely find the correct line
of instruction in ``.bash_profile`` and then simply copy it and place it in
your shell's ``.profile`` and ``.rc`` files so that they point to your Anaconda
binaries. However, your shell's ``.profile`` and ``.rc`` files will have
different names, so refer to
`this wikipedia page <https://en.wikipedia.org/wiki/Unix_shell>`_
for some guidance on the names of configuation files for popular shells.

Note the syntax for the shell commands to do this will probably also be
different in each shell! But since both Anaconda and Canopy assume you are
using Bash, another clue will be to start with a copy of the (correct) path
they placed in your `.profile` or `.bash_profile` file and copy that line to
both the ``.profile`` and ``.rc`` files of the shell you are using.

**The skinny:** Find any references to the distribution you want and make sure
they are inserted near the beginning of the path for the shell you are using.
Alternatively, you can simply put the comment-out symbol ``#`` in front of the
lines that are referring to the incorrect path you see when you type
``which python``.

Type ::

     > which python
     > which ipython

In both cases the path should be the same and reference your distribution.


PC
``

On a PC, the same principle of modifying your environment variables applies,
but you access them differently. Go to the Control Panel, then System. On
**Windows 8**, you then want Advanced System Settings, though this will be
similar on older OSes. Go to Advanced, then to the `Environment Variables...`
button. Under User Variables, see if there is an entry called PATH. If there
is, we will modify it. If there isn't, we will create one. It is
**VERY IMPORTANT** that you do not modify any existing text, *especially*
under `System Variables` below.

As is the situation for Mac, above, the system reads these PATH strings from
left to right, and stops once it has found what it is looking for. It also
reads User before System variables. Hence, we want to add new strings to the
left hand (start) of the existing text, if there is any.

First, scan the existing string(s) (including under System) to see if there is
any reference to the Python distribution you are trying to set as default
already there. e.g., my User PATH (running Anaconda cleanly) currently reads::

    C:\Users\Dan\AppData\Local\Continuum\Anaconda;C:\Users\Dan\AppData\Local\Continuum\Anaconda\Scripts

If you find a reference or references like this to the version you're currently
trying to run, copy the text, and add it (repeated) at the start of the User
string. Copy this syntax—semicolons separate paths.

If you can't find any reference to your chosen version (Canopy/Anaconda),
you'll need to add the PATH yourself. For Anaconda, assuming you installed it
in the default directory, add the above string. For Canopy, use the "Set
Canopy as default" option ("the easy way"), which really should work. See
`this page <http://docs.enthought.com/canopy/quick-start/install_windows.html>`_
for more information on the PATHs used by Canopy if you're still struggling.

If you are on **Windows 10**, you need to make sure you see these paths.

If you installed for a single user::

    C:\Users\your_user_name\Anaconda3
    C:\Users\your_user_name\Anaconda3\Scripts

If you installed for all users::

    C:\ProgramData\Anaconda3
    C:\ProgramData\Anaconda3\Scripts

Note, if you aren't sure how you installed, just search for 'Anaconda3' on the
main drive to find where it was installed.

Note that modifying the User Variables will only affect the current user
account. Add the text—carefully!!—to the System Variables if you want the
changes for all users.

Type ::

     > where python
     > where ipython

In both cases the path should be the same and reference your distribution.

Other issues
------------

Other install issues often mean that some component of your Python distribution
is out of date. A very common culprit is ``setuptools``, which—extremely
frustratingly—isn't updated by a ``conda update --all`` call for Anaconda.
Other packages can also cause this kind of problem if out of date. An example
of a ``setuptools`` related error we've seen recently ends with::

    error: unknown file type '.pyx' (from 'landlab/components/flexure/cfuncs.pyx')

...combined with warnings referencing a problem with PEP 440.

To our knowledge, this issue only arises for developer installs.

Resolve the issue by updating your distribution. For Anaconda, from a terminal just run::

    > conda update --all
    > conda update setuptools

Finally, if you are still having problems, you can use the nuclear option and
start again from scratch.

For example, your Anaconda distribution is contained in one folder. You can
move this folder to the trash and install a fresh version following the
directions on the `Anaconda <https://www.anaconda.com/distribution/>`_ site.

Update ``conda`` and ``pip``, uninstall Landlab, and then install a fresh copy.
