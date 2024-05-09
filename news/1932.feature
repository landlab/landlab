
Added a new *vtk* writer, ``landlab.io.legacy_vtk.dump`` that is
able to write *Landlab* grids that have three spatial coordinates.
This function is also able to write both the main grid (*nodes* and
*patches*) as well as the dual grid (*corners* and *cells*).
