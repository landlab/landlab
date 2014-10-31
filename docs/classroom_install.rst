.. _classroom_install:

===============================
The Landlab "Classroom Install"
===============================


For all academic or research use, we strongly recommend that users follow the main 
:ref:`install instructions <install>` for Landlab, using Github or Git. However, this
process is not ideal for a teaching environment, where getting a stable, working version
of landlab onto a large number of machines as quickly and easily as possible will be a
priority. This is the reason for the *classroom install* detailed here.

.. note:
    
    Although faster to perform once, we emphasise that this method only provides a single
    "snapshot" of Landlab. You (or your students) will be unable to make changes to the
    Landlab code, and will not be able to update Landlab to future (more stable!) 
    versions.
    
    If you do need to update Landlab after having done this once (e.g., running a future
    class on computers that have already had Landlab installed on them in a previous 
    year), delete the existing Landlab files from the hard drive, then repeat the whole of
    this procedure.


Step 1: Check or get your Python distribution
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

You will still need a Python distribution fulfilling our system requirements (Python 2.7,
Numpy 1.8 or greater, Scipy 0.12 or greater). As before for the full install, we recommend
a scientific Python distribution like Canopy or Anaconda. See the :ref:`main page 
<install>` for full instructions.

If you have to use a departmental lab for this install, you can check the versioning
of existing packages on that machine from a command prompt::

    python --version            #must be >=2.7
    python                      #launch a python shell
    >>>import scipy, numpy 
    >>>scipy.__version__        #must be >=0.12
    >>>numpy.__version__        #must be >=1.8


Step 2: Download the code
>>>>>>>>>>>>>>>>>>>>>>>>>

Go to `the Github Landlab page <https://github.com/landlab/landlab>`_. There's no need
to register for an account. On the right hand side of the page, click the **Download
ZIP** button. This should start an automatic download of the zipped version of the
Landlab code files.

Unzip the file into an appropriate file location on your system.


Step 3: Setup Landlab in Python
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Open a command prompt and navigate to the directory into which you unzipped the downloaded
folder. Check it contains `setup.py`. Now simply type::

    python setup.py install

Check the install worked by opening Python at the prompt and testing the Landlab import::

    python
    >>>import landlab       #no errors indicates a successful import
    
