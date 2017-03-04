A Brief Introduction to Using Python, Numpy, Scipy, and other Landlab Dependencies
==================================================================================

Getting to know Python
----------------------

We recommend you approach Landlab with a basic working knowledge of the Python coding language. A good, concise, complete beginner’s guide that will get you to enough knowledge to be able to get started can be found here. We like XXX and XXX as more comprehensive introductions to Python.

If you’re already familiar with Matlab, you will probably feel fairly at home with Python fairly quickly. However, there are some critical differences. Important things to remember are:

Python’s indexing is inclusive at the start and exclusive at the end (in contrast to Matlab). e.g.::

    >>> numpy.arange(0,100)

will give an array of 100 numbers, starting at 0 and ending at 99.

Python doesn’t use parentheses or brackets to delimit code blocks (functions, loops, if statements, etc). Instead it uses a colon to declare the start of a code block, then spaced indenting (normally 4 spaces) to demark which lines of code belong in the block. e.g.::

    def myfunction(input_param):
        if type(input_param) == str:
            print “The input to the function said: “, input_param
        else:
            print “The input parameter wasn’t a string.”
            print “It was “, input param
    
Lines don’t need to end with the semicolon to suppress output; Python won’t print output unless you explicitly call print.

Landlab is also written in an object-oriented coding style. Many of the elements of Landlab that you will interact with - grids, components, utilities - are Python objects, which means they contain both data that describe the object and functions that you can call that operate on the object. Think of the object as a container in which is stored everything relevant to that part of the code, so it can be accessed easily. You can read a bit more about Python objects here for general information. We give more detail on what this means in terms of running our code later in this guide (hotlink to parts about designing a driver).


Numpy and Scipy, and Efficient Coding Style
-------------------------------------------

Numpy and scipy are the workhorse scientific computing packages of Python. They provide fast, efficient, and surprisingly comprehensive data structures and numerical methods that we (and you) can exploit to make coding in Landlab faster and easier.

In particular, Landlab makes extensive use of the numpy array data structure. Almost all data input and output from Landlab is in this form (see “Landlab Fields” for more information). These arrays allow operations to happen much faster on the data than would be possible in a pure Python data structure like a list or dictionary. (This is possible because arrays partially suppress some of Python’s inbuilt type checking and memory management, and impose a more ordered structure on the way raw data is held in your computer’s memory).

However, in order to exploit the speed gains that numpy can give you, you’ll need to adopt a coding style quite different to what would be natural in, say, C++ (or likely, Matlab). A typical bottleneck in Python code occurs when looping over data, and numpy arrays typically let you avoid doing this. So if you find yourself about to write something like this::

    >>> for i in range(len(myarray)):
            myoutputarray[i] = myoutputarray[i] + myarray[i]

Don’t! Try to develop a coding style where each line operates on the whole array at once, like::

    >>> myoutputarray += myarray

In particular, it can be very tempting to use loops to apply a condition over a whole array. Try not to do this! e.g.,::

  >>> for i in myarray:
          if i < 0:
              i=0

This will be really slow. Do this instead::

  >>> myarray[myarray < 0] = 0

You can read a lot more about writing efficient code using numpy on a large number of websites. Try starting here: http://www.astro.washington.edu/users/vanderplas/Astr599/notebooks/11_EfficientNumpy
We also strongly recommend the book “High Performance Python” by Gorelick and Ozsvald, published by O’Reilly, if you’re looking for a more comprehensive treatment.

Cython
------

If you go poking around in the Landlab source code code, you will discover that not all of Landlab is written in pure Python. Some of it is written in Cython. Cython is a very closely related programming language to Python, and indeed, all code written in pure Python is automatically also Cython code. Cython is probably best thought of as a cross between C++ and Python, which aims to combine the flexibility of Python with the brute power and granular control over your code that C++ provides. e.g., if there are sections of code where looping through an array is unavoidable, Cython provides a way of significantly accelerating the speed of this code. Cython code largely looks like straightforward Python, but may have type declarations or other C++-like features within it.

From the user’s perspective, the most important thing to note is that Cython is a compiled language. (This isn’t true of Python, which is an interpreted - i.e., compiled at run time - language.) We provide the pre-compiled executables you will need to run Landlab when you install, and this should be sufficient for the vast majority of users.

However, note that if you intend to modify any of the sections of code that we provide to you, you will probably need to recompile that code on your machine before the changes take effect. See XXX for lots more information on this.
