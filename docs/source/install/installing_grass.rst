.. _grass:

=====================
Landlab and GRASS GIS
=====================

Landlab has some known issues when operating alongside GRASS GIS, the
open-source GIS software package. If you're familiar with GRASS, this
will probably come as no surprise to you, since GRASS has known issues a
lot of the time. We believe these difficulties only apply to Apple Macs,
but it's always possible these exist on other platforms too.

The specific problem relates to GRASS needing to bind to your machine's
32bit system Python when installed on a Mac, as opposed to the Conda
python distribution you will install during Landlab setup. If you
install GRASS before Anaconda, you should not have any problems, but if
you install Anaconda first, you will not have a good day. However, there
is a (relatively) simple work-around for this if you do end up in this
situation, and can't face removing Anaconda entirely while you install
GRASS.

Note that on modern Macs, there is also a second critical issue which
will prevent you from installing (and indeed launching) GRASS. GRASS is
fundamentally incompatible with an inbuilt security widget that is
embedded in both OSX El Capitan and Sierra. It is nonetheless possible
to disable this security, however (see below).

Installing GRASS after Anaconda
-------------------------------

(Before trying this, test that you definitely cannot install GRASS.
There's always a chance that you already had the troublesome libraries
and won't need to do this.)

The key issue is that GRASS GIS has a python library called ``wxpython``
as one of its dependencies, but that ``wxpython`` only exists as a 32bit
release for Mac system architecture. Conversely, Anaconda is only easily
accessible as a 64bit release. Happily, this is only an issue while
installing GRASS, so if you're smart you can set the system back to its
original Python during install, then switch it back to Anaconda once you
have GRASS up and running.

This is done by altering the PATH variable, which Anaconda will have
silently modified during installation. This can be found in the (hidden)
file in your system profile called - most likely - ``.bash_profile``.
Simply open the file in a text editor, comment out (#) the line or lines
that add Anaconda to the path (mine reads simply
``export PATH="/Users/daniel/anaconda/bin:$PATH"``), then save again.
This may be enough to allow install, but in order to make sure GRASS can
continue to open once you restore the path, also add a line reading
``export GRASS_PYTHON=/System/Library/Frameworks/Python.framework/Versions/Current/bin/pythonw``
BEFORE the commented out line setting the PATH.

You can find out a bit more about messing with your system PATH
:ref:`here <correcting_install_paths>`.

Once you've done this, you can follow the normal GRASS instructions,
installing both its framework packages and GRASS itself. Fire GRASS up
once and check it starts. If it starts up now, you can safely reinstate
the line you commented out in the profile file, and GRASS and Anaconda
will run happily side by side.

Except…

GRASS won't run on OSX El Capitan and later
-------------------------------------------

This isn't really our problem, but we suspect it will nonetheless catch
out many of our users. Apparently, Apple have added an inbuilt layer of
security ("CSR") to their newest operating systems, which treats GRASS
as hostile software and stops it running. A solution which you may be
able to find online is to simply disable this security. This is done by
rebooting into recovery mode, opening a terminal, and running
``csrutil disable``, then restarting as normal. *Allegedly*, doing this
creates no additional security risks beyond what have been present on
older OSX releases for a long time, but PROCEED AT YOUR OWN RISK. You'll
have to leave it turned off to get GRASS to run. (Reinstate it with
``csrutil enable``, as you'd expect.)

More details on this piece of black magic can be found here:
http://grassmac.wikidot.com
