.. _compile_in_windows:

==============================================================
Installing compilers for building Landlab on a Windows machine
==============================================================

On a mac, all the necessary C++ compilers needed by Python to build new C++ modules from
source code are provided simply by installed the Xcode app from the app store. However,
on a PC, things are a bit more involved. This document will guide you through the process
of getting the necessary compilers, and getting Python talking to them.

.. note::

    This process is only necessary if you want to build Landlab directly from its source
    code on your PC. This is necessary because some parts of Landlab use compiled C++
    code to accelerate things, using Cython. If you are using a binary install, as
    described :ref:`here <install>`, you don't need to do this.
    

