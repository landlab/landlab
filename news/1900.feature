
Added a new function, calc_net_face_flux_at_cell, that computes the
net flux of a quantity into each of a RasterModelGrid's cells. This
function uses openmp to parallelize its calculations.
