[Landlab](http://landlab.github.io) |
[[About | About]] |
[[Examples | Examples]] |
[[User Guide | User-Guide]] |
[Reference Manual](http://landlab.readthedocs.org/en/latest/#developer-documentation) |
[[Tutorials| Tutorials ]] |
[[FAQs |FAQs]]

<p>
    This should work for <a href="https://www.continuum.io/downloads">Anaconda</a> users with Windows 7+, Mac OS 10.6+, or Ubuntu
    Linux (only the latest version has been tested).
  </p>

  <p>
    Once you have <a href="https://github.com/landlab/landlab/wiki/Installing-Python">a full Python distribution on your machine</a>, it is vital to
    check that it has been successfully set as the default copy of Python on
    your system. Open a command prompt (Terminal on a Mac, or Command Prompt
    on a PC) and type the lines below (note the <code>></code> indicates that you are on a
    command line):
  </p>

  MAC:
<p>
   <code>
     > which python
   </code>
  </p>
  <p>
   <code>
     > which ipython
   </code>
  </p>
  <p>

 WINDOWS:
<p>
   <code>
     > where python
   </code>
  </p>
  <p>
   <code>
     > where ipython
   </code>
  </p>
  <p>
    In each case, both commands should return the same path, and it should clearly refer to Anaconda (or
    Canopy). Details will depend on your operating system but it could look something like this:
  </p>

  <p>
   <code>
     /anaconda/bin/python
   </code>
  </p>
  <p>
    If you don’t see reference to your newly installed
    distribution (i.e., <code>/anaconda</code>), <a href="https://github.com/landlab/landlab/wiki/Correcting-Install-Paths">click here</a> to resolve the problem.
  </p>

</p>
   Make sure you have the latest version installed (close anaconda before doing this):
 
</p>

   <code>
    conda update --all
   </code>
  <p>
  <p>
    Once the path to both <code>python</code> and
    <code>ipython</code> point to your new distribution, open the
    Python editor in Anaconda called Spyder.
  </p>

  <p>
    On the Spyder toolbar, go to <code>Tools → Open</code> command prompt to open the command
    line.
  <p>

  </p>
    Alternatively you can open a standard terminal window, such as an
    xterm (X11.app) or terminal window (Terminal.app) on a Mac, or a command
    prompt on a Windows machine. If you do use a standard terminal and run into
    problems, make sure you have <a href="https://github.com/landlab/landlab/wiki/Correcting-Install-Paths">resolved your path issues.</a>
  </p>

<h3>Now to install Landlab!</h3>
  <p>You can either install with the conda or the pip package managers. Conda is recommended, as it reduces the chances of versioning conflicts. Try to remember which you choose to avoid confusion when updating later! (If you installed landlab prior to May 19th 2016, you will have used pip).
  </p>

  <p>
    Type either (for conda install):
    <p>
        <code>
          > conda install landlab -c landlab -c conda-forge
        </code>
    </p>
  </p>

  <p>
    ...or (for pip install, not recommended for entry level users):
    <p>
        <code> 
            > pip install landlab 
        </code>
    </p>
  </p>

  <h3>Test Landlab install</h3>
<p> Once Landlab has been successfully installed, on the
    Python shell line, check to make sure it is up-to-date (note that those are
    double underscores around version; also note that you may need to close and
    reopen Anaconda before typing the below commands):
  </p>

  <p>
    <code>
      > import landlab
    </code>
  </p>
  <p>
    <code>
      > landlab.__version__
    </code>
  </p>

   <p>
     The version number should be greater than 1. You can check the version number of the most recent release <a href= "https://github.com/landlab/landlab/releases">here</a>.
   </p>




<h3> Install/Test problems</h3>

  <p>
     If you are having problems when installing, testing or running Landlab, please visit our <a href= "https://github.com/landlab/landlab/wiki/Troubleshooting">Troubleshooting page</a>.
   </p>