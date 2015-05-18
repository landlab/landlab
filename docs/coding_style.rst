Getting started with Landlab
=======================

Landlab is coded in Python and exploits and includes as dependencies a number of widely used scientific Python packages - in particular, numpy and scipy. The decision to code in Python was explicitly made to lower the bar to entry for Landlab; to increase the flexibility and reusability of the code base; and to increase development speed both for the core development team and for future users. The choice of Python also means that developers using Landlab can take advantage of that language’s affinity for rapid development timescales.

Other advantages of this choice include high portability between platforms, open source language, numerous existing scientific libraries that prevent developers having to “reinvent the wheel”, and support for selective optimization of time-critical parts of the code base in Cython (see below).

Dependencies
------------

To run Landlab, you will need Python and several of its packages on your machine. These will also need to be at least certain versions to run. These dependencies are:
Python 2.7
Numpy 1.8 or greater
Scipy 0.12 or greater

We recommend installing a Python distribution like the Enthought Python Distribution or Anaconda; these also have the advantage of providing you with an interactive development environment (IDE) in which to view and edit the code. However, as long as you meet the dependency requirements above, several other options will work. Find out more about installing Python here.



