.. _cite_as:

Landlab Citation Registry
-------------------------

A relatively new interface also automates the process of extracting citations
for landlab.

.. code-block:: python

  import landlab

  # the code of your model built with landlab goes here, then add a call to

  landlab.registry.format_citations() # will produce a Bibtex-formatted
                                      # citations for all Landlab components
                                      # that you currently have instantiated.
