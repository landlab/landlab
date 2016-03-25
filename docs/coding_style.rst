Why Python?
===========

Landlab is coded in Python and exploits and includes as dependencies a number of widely used scientific Python packages - in particular, numpy and scipy. The decision to code in Python was explicitly made to lower the bar to entry for Landlab; to increase the flexibility and reusability of the code base; and to increase development speed both for the core development team and for future users. The choice of Python also means that developers using Landlab can take advantage of that language’s affinity for rapid development timescales.

Other advantages of this choice include high portability between platforms, open source language, numerous existing scientific libraries that prevent developers having to “reinvent the wheel”, and support for selective optimization of time-critical parts of the code base in Cython (see :ref:`cython`).

Dependencies
============

To run Landlab, you will need Python and several of its packages on your machine. These will also need to be at least certain versions to run. These dependencies are:

* Python 2.7
* Numpy 1.8 or greater
* Scipy 0.12 or greater

We recommend installing a Python distribution like the Enthought Python Distribution or Anaconda; these also have the advantage of providing you with an interactive development environment (IDE) in which to view and edit the code. However, as long as you meet the dependency requirements above, several other options will work. Find out more about installing Python here: :ref:`install`.


Getting to know Python
======================

We recommend you approach Landlab with a basic working knowledge of the Python coding language. A good, concise, complete beginner’s guide that will get you to enough knowledge to be able to get started can be found here. We like `the Software Carpentry intro <http://software-carpentry.org/v4/python/>`_  and `these Python notebooks <http://nbviewer.ipython.org/github/jrjohansson/scientific-python-lectures/tree/master/>`_ as more comprehensive introductions to Python.

If you’re already familiar with Matlab, you will probably feel fairly at home with Python fairly quickly. However, there are some critical differences. Important things to remember are:
Python’s indexing is inclusive at the start and exclusive at the end (in contrast to Matlab). e.g.,

>>> numpy.arange(0,100)
    
will give an array of 100 numbers, starting at 0 and ending at 99.
    
Python doesn’t use parentheses or brackets to delimit code blocks (functions, loops, if statements, etc). Instead it uses a colon to declare the start of a code block, then spaced indenting (normally 4 spaces) to demark which lines of code belong in the block. e.g.,

>>> def myfunction(input_param):
...     if type(input_param) == str:
...         print “The input to the function said: “, input_param
...     else:
...         print “The input parameter wasn’t a string.”
...         print “It was “, input param
    
Lines don’t need to end with the semicolon to suppress output; Python won’t print output unless you explicitly call print.

Finally, but importantly, Python doesn't use the hat (.^) as its raise-to-the-power symbol. Instead, it uses a double star (**). Simple, but typically very frustrating for a day or two during transition! There's also the numpy method, np.square, which if you're using arrays typically outperforms the ** operator.

We have a very short tutorial on Python and numpy from the point of view of Landab (and the key differences with Matlab) `here <http://nbviewer.ipython.org/github/landlab/drivers/blob/master/notebooks/Python_intro.ipynb>`_.

Landlab is also written in an object-oriented coding style. Many of the elements of Landlab that you will interact with - grids, components, utilities - are Python objects, which means they contain both data that describe the object and functions that you can call that operate on the object. Think of the object as a container in which is stored everything relevant to that part of the code, so it can be accessed easily. You can read a bit more about Python objects `here <http://learnpythonthehardway.org/book/ex40.html>`_ for general information. We give more detail on what this means in terms of running our code :ref:`later in this guide <drive_a_model>`.


Numpy and Scipy, and Efficient Coding Style
===========================================

Numpy and scipy are the workhorse scientific computing packages of Python. They provide fast, efficient, and surprisingly comprehensive data structures and numerical methods that we (and you) can exploit to make coding in Landlab faster and easier.

In particular, Landlab makes extensive use of the numpy array data structure. Almost all data input and output from Landlab is in this form (see :ref:`Landlab Fields <fields>` for more information). These arrays allow operations to happen much faster on the data than would be possible in a pure Python data structure like a list or dictionary. (This is possible because arrays partially suppress some of Python’s inbuilt type checking and memory management, and impose a more ordered structure on the way raw data is held in your computer’s memory).

However, in order to exploit the speed gains that numpy can give you, you’ll need to adopt a coding style quite different to what would be natural in, say, C++ (or likely, Matlab). A typical bottleneck in Python code occurs when looping over data, and numpy arrays typically let you avoid doing this. So if you find yourself about to write something like this:

>>> for i in range(len(myarray)):
        myoutputarray[i] = myoutputarray[i] + myarray[i]

Don’t! Try to develop a coding style where each line operates on the whole array at once, like:

>>> myoutputarray += myarray

In particular, it can be very tempting to use loops to apply a condition over a whole array. Try not to do this! e.g.,

>>> for i in myarray:
        if i < 0:
            i=0
            

This will be really slow. Do this instead:

>>> myarray[myarray<0] = 0

You can read a lot more about writing efficient code using numpy on a large number of websites. For example, `UW's astronomy department has a great online intro <http://www.astro.washington.edu/users/vanderplas/Astr599/notebooks/11_EfficientNumpy>`_.
We also strongly recommend the book “High Performance Python” by Gorelick and Ozsvald, published by O’Reilly, if you’re looking for a more comprehensive treatment.

.. _cython:

Cython
======

If you go poking around in the Landlab source code code, you will discover that not all of Landlab is written in pure Python. Some of it is written in Cython. Cython is a very closely related programming language to Python, and indeed, all code written in pure Python is automatically also Cython code. Cython is probably best thought of as a cross between C++ and Python, which aims to combine the flexibility of Python with the brute power and granular control over your code that C++ provides. e.g., if there are sections of code where looping through an array is unavoidable, Cython provides a way of significantly accelerating the speed of this code. Cython code largely looks like straightforward Python, but may have type declarations or other C++-like features within it.

From the user’s perspective, the most important thing to note is that Cython is a compiled language. (This isn’t true of Python, which is an interpreted - i.e., compiled at run time - language.) We provide the pre-compiled executables you will need to run Landlab when you install, and this should be sufficient for the vast majority of users.

However, note that if as a developer you intend to modify any of the sections of code that we provide to you, you will probably need to recompile that code on your machine before the changes take effect. See the :ref:`development guide <dev_install>` for lots more information on this.


