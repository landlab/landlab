.. note::

    The following instructions assume you have a working version of `Git
    <http://git-scm.com/>`_ installed on your system. Git is a
    distributed version control system (DVCS) and source code management
    system. For an introduction to Git and DVCS, see the official
    `git documentation <http://git-scm.com/documentation>`_. Installing the
    Github graphical user interface (see below) will give you the necessary
    git tools.


.. _source-install:

The Landlab code lives in the `Github <https://github.com>`_ online code repository. For install, 
you have two choices. Firstly, you can sign up to Github as a user of that website, 
download their third party GUI, and use that to get a copy of the code. 
Alternatively, you can manage the code acquisition directly through the command line 
on your machine using the Git text interface. Using Github takes slightly longer, 
but their graphical interface is arguably more straightforward - especially for updating
Landlab once you have it installed.

.. _gui-install:

With GitHub GUI
>>>>>>>>>>>>>>>

#. Go to the `Github webpage <https://github.com>`_. The homepage will prompt you to sign
   up for an account. Do so! Choose the free plan. Note down your user name and password.
#. Install the `GitHub app 
   <https://help.github.com/articles/set-up-git>`_. Follow the directions for
   installing the native app for your operating system.
     * `Mac <https://mac.github.com>`_
     * `Windows <https://windows.github.com>`_
     * Linux: Follow the command-line :ref:`installation instructions
       <command-line-install>`.
#. Open the app. You need to provide it with your user name and password to allow it to
   interact smoothly with the website. You should be prompted to do so when it boots up
   for the first time. If not, go to Preferences and enter your sign-in details. Click 
   through the remainder of the options, skipping the "add repositories" step.
#. With your browser, go to the `landlab page
   <https://github.com/landlab/landlab>`_ on GitHub and click the "Clone in
   Desktop" button (midway down the right hand side of the page). This will automatically
   cause your machine to switch back to the Github app and begin the download process. 
   Pick a location to store the Landlab files on your hard drive, and click through.
   Download will begin.
#. Now, leave the Github app and open a command prompt (PC) or Terminal (Mac/Unix). 
   Navigate to the root directory of your Landlab download (reminder: change directory
   in a prompt/terminal using the command ``cd``, then the name of the subfolder; 
   ``cd ..`` takes you up one folder level). This root directory will contain a file
   called `setup.py` (check with ``dir`` (PC) or ``ls`` (Mac/Linux)).
   From this directory, type at the prompt::

        > python setup.py develop

.. note::
    
    This command tells your install of Python that `landlab` is a Python module that 
    you have now installed on your system, and where to look for the files it needs
    to run. Using the keyword `develop` warns Python that the code you have saved 
    on your disc might change from time to time. This
    means that should you so desire, you can make changes to the code, add 
    functionality, add your own modules, or otherwise tinker with the .py files you
    will find in the directories that Github has placed on your system. Importantly,
    however, it also allows to you quickly and easily use Github to download more
    up-to-date versions of Landlab - which may contain bug fixes, etc. For more on
    updating your installation of Landlab, click :ref:`here <updating_landlab>`.
        
    
#. Finally, test everything worked. From the same command line, type::
    
       > python
    
   An interactive Python window will open in the command line; the prompt will look like
   ``>>>``. From here, enter::
    
        >>> import landlab
    
   If you are returned to the >>> prompt after a few moments, everything is fine. If you
   see an error message, you might have some problems with your install. See the 
   :ref:`install FAQ page <install_FAQ>` for a list of known install issues, and their 
   solutions. 
   
   Leave the Python shell by typing::
   
        >>> exit()
      

.. _command-line-install:

With Git
>>>>>>>>

.. note::

    This assumes that you already have Git on your machine. To check, open a command 
    prompt and type ``git``. If you have it, you will see usage instructions. If you
    don't, you will see an error message.

#. Using the command prompt, clone landlab from the master repository. This is 
   hosted on `github.com <http://www.github.com>`_. The files will be added inside 
   whichever directory you are in when you enter this command.::

    > git clone https://github.com/landlab/landlab.git

#. Navigate From the root directory of your landlab clone (the folder that contains
   `setup.py`). From your likely current location this will probably just be 
   ``cd landlab``. From here, enter::

    > python setup.py develop

#. Finally, test everything worked. From the same command line, type::
    
     > python
    
   An interactive Python window will open in the command line; the prompt will look like
   ``>>>``. From here, enter::
    
      >>> import landlab
    
   If you are returned to the >>> prompt after a few moments, everything is fine. If you
   see an error message, you might have some problems with your install. See the 
   :ref:`install FAQ page <install_FAQ>` for a list of known install issues, and their 
   solutions. 
   
   Leave the Python shell by typing::
   
      >>> exit()

You can find more details about installing Landlab as a developer :ref:`here 
<dev_guide>`.
