.. _developer_install:

===================================
Developer Installation Instructions
===================================

If you intend to modify the Landlab code base, including design of new
components, you will want to install Landlab directly from the source
code. This way, when you modify the code, the changes you have made will
be reflected immediately in your Landlab library. This section provides
information on how to do this.

Installation Steps and Resources
--------------------------------

.. toctree::
   :maxdepth: 2

   fork
   install
   test
   updating
   troubleshooting_install
   installing_windows_compiler

**Note:** For dev work, we actively recommend Anaconda over the
Enthought Python Distribution, especially on Windows machines. This is
because it ships with a working compiler already associated with Python,
whereas the EPD does not. On a Mac, this is less important as the Xcode
app (available through iTunes) gives you the necessary compilers
instead—install it now if you don't have it! If you choose to use the
EPD on a Windows machine, however, you'll need to install separately
either Visual Basic or MinGW and successfully associate them with your
Python install. See :ref:`this page<compile_in_windows>` on Windows Compilers.

Make a GitHub Issue to contact the development team if you're really struggling.
But unless you're really invested in Canopy and the EPD, uninstalling it and
replacing with Anaconda is probably the more stress-free way to go.*

*Either way, you'll need a working C++ compiler running alongside Python
to be able to perform a full developer install. You'll see errors
referring to* :ref:`Cython <cython>` *if you
don't have working compiler when calling* ``python setup.py develop``
*(see* :ref:`the developer install instructions <dev_install_install>` *).*

Working with your local version of Landlab
------------------------------------------

Obviously, feel free to just dive into modifying the code, but your life
in the future will be a bit easier if you follow some basic
recommendations for good work flow with git forks and branches. Even if
you have a working knowledge of using git in a collaborative project, we
highly recommend that you review :ref:`this section<landlab_develop_with_git>`
of the documentation to get a sense of how to track
modifications to your version of Landlab in a way that makes it easy to
(a) get updates to Landlab made by the development team and other
contributors and (b) contribute improvements and new features you
develop back to the community.

For information about our in-house code
formatting conventions and standards, see :ref:`here<style_enforcement>`.
