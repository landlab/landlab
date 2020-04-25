.. _compile_in_windows:

==============================================================
Installing compilers for building Landlab on a Windows machine
==============================================================

On a Mac, all the necessary C++ compilers needed by Python to build new C++
modules from source code are provided simply by installing the Xcode app from
the app store. On a PC you will get these compilers from Build Tools for Visual
Studio.

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
until you attempt to perform the developer installation

.. code-block:: bash

 $ python setup.py develop

If you have compiler issues, you will see an error message. On April 24, 2019
on a clean install on Windows 10 this error message said: ::

  building 'landlab.ca.cfuncs' extension
  error: Microsoft Visual C++ 14.0 is required. Get it with "Build Tools for
  Visual Studio": https://visualstudio.microsoft.com/downloads/

Go to the `specified website <https://visualstudio.microsoft.com/downloads/>`_
and download the Microsoft Visual C++ compiler. In April 2019 this was done by
first downloading "Build Tools for Visual Studio". Then when installing Build
Tools there were options for what exactly should be installed. One of these was
the Visual Studio C++ compilers. Once installed, a restart was required and
landlab compiled as expected. 
