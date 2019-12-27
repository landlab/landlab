.. _compile_in_windows:

==============================================================
Installing compilers for building Landlab on a Windows machine
==============================================================

On a Mac, all the necessary C++ compilers needed by Python to build new C++
modules from source code are provided simply by installing the Xcode app from
the app store. However, on a PC, things are a bit more involved. This document
will guide you through the process of getting the necessary compilers, and
getting Python talking to them.

This process is only necessary if you want to build Landlab directly from its
source code on your PC. Some parts of Landlab use compiled C++ code to
accelerate things, via Cython. If you are using a binary install, as described
:ref:`here <anaconda_install_landlab>`, you don't need to do this.

Getting the compilers
---------------------

Anaconda provides some of the necessary compilers. Install it first.

The easiest way to discover if you have other compiler issues is to attempt an
install! Follow the
:ref:`developer installation instructions <developer_install>`
until you attempt to perform the developer installation ::

.. code-block:: bash

   $ python setup.py develop

If you have compiler issues, you will see the message::

  ...building 'landlab.components.flexure.cfuncs' extension
  error: MS VisualC++ 9.0 is required (Unable to find vcvarsall.bat).
  Get it from https://www.microsoft.com/en-us/download/details.aspx?id=44266

Simply go to `that specified Microsoft website <https://www.microsoft.com/en-us/download/details.aspx?id=44266>`_ and from
there download the software needed using the obvious download button.

Once you have it, the above install command should then work fine.
