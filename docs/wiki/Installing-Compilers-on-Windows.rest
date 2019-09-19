`Landlab <http://landlab.github.io>`_ |
[[ About | About ]] |
[[ Examples | Examples ]] |
[[ User Guide | User-Guide ]] |
`Reference Manual <http://landlab.readthedocs.org/en/latest/#developer-documentation>`_ |
[[ Tutorials| Tutorials ]] |
[[ FAQs | FAQs ]]

.. _compile_in_windows:

==============================================================
Installing compilers for building Landlab on a Windows machine
==============================================================

On a Mac, all the necessary C++ compilers needed by Python to build new C++ modules from
source code are provided simply by installing the Xcode app from the app store. However,
on a PC, things are a bit more involved. This document will guide you through the process
of getting the necessary compilers, and getting Python talking to them.

This process is only necessary if you want to build Landlab directly from its source code on your PC. Some parts of Landlab use compiled C++ code to accelerate things, via Cython. If you are using a binary install, as described [[ here | Installing-Landlab-with-Anaconda#now-to-install-landlab ]], you don't need to do this.
    
Getting the compilers
=====================

Anaconda provides some of the necessary compilers. Install it first.

The easiest way to discover if you have other compiler issues is to attempt an install! Follow the [[ main instructions | Installing-Landlab-from-source-code-("developer-install") ]] until you attempt to perform the main install::

  > python setup.py develop

If you have compiler issues, you will see the message::

  ...building 'landlab.components.flexure.cfuncs' extension
  error: MS VisualC++ 9.0 is required (Unable to find vcvarsall.bat).
  Get it from http://aka.ms/vcpython27

Simply go to `that specified Microsoft website <http://aka.ms/vcpython27>`_ and from
there download the software needed using the obvious download button.

Once you have it, the above install command should then work fine.