from multiprocessing import Pool

import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport fabs
from libc.stdlib cimport abs


DTYPE = np.double
ctypedef np.double_t DTYPE_t


@cython.boundscheck(False)
def subside_parallel_row(np.ndarray[DTYPE_t, ndim=1] w,
                         np.ndarray[DTYPE_t, ndim=1] load,
                         np.ndarray[DTYPE_t, ndim=1] r,
                         DTYPE_t alpha,
                         DTYPE_t gamma_mantle):
  cdef int ncols = w.shape[0]
  cdef double inv_c = 1. / (2. * np.pi * gamma_mantle * alpha ** 2.)
  cdef double c
  cdef int i
  cdef int j

  for i in range(ncols):
    if fabs(load[i]) > 1e-6:
      c = load[i] * inv_c
      for j in range(ncols):
        w[j] += - c * r[abs(j - i)]


def subside_grid(np.ndarray[DTYPE_t, ndim=2] w,
                 np.ndarray[DTYPE_t, ndim=2] load,
                 np.ndarray[DTYPE_t, ndim=2] r,
                 DTYPE_t alpha, DTYPE_t gamma_mantle):
  cdef int nrows = w.shape[0]
  cdef int i
  cdef int j

  for i in range(nrows):
    for j in range(nrows):
      subside_parallel_row(w[j], load[i], r[abs(j - i)], alpha, gamma_mantle)


def subside_grid_strip(np.ndarray[DTYPE_t, ndim=2] load,
                       np.ndarray[DTYPE_t, ndim=2] r,
                       DTYPE_t alpha, DTYPE_t gamma_mantle, strip_range):
  (start, stop) = strip_range

  cdef np.ndarray w = np.zeros((stop - start, load.shape[1]), dtype=DTYPE)
  cdef i

  for i in range(load.shape[0]):
    for j in range(start, stop):
      subside_parallel_row(w[j - start], load[i], r[abs(j - i)], alpha,
                           gamma_mantle)

  return w, strip_range


def tile_grid_into_strips(grid, n_strips):
    rows_per_strip = grid.shape[0] // n_strips

    starts = np.arange(0, grid.shape[0], rows_per_strip)
    stops = starts + rows_per_strip
    stops[-1] = grid.shape[0]

    return zip(starts, stops)


def _subside_grid_strip_helper(args):
  return subside_grid_strip(*args)


def subside_grid_in_parallel(np.ndarray[DTYPE_t, ndim=2] w,
                             np.ndarray[DTYPE_t, ndim=2] load,
                             np.ndarray[DTYPE_t, ndim=2] r,
                             DTYPE_t alpha, DTYPE_t gamma_mantle, n_procs):
    if n_procs == 1:
        return subside_grid(w, load, r, alpha, gamma_mantle)

    strips = tile_grid_into_strips(w, n_procs)

    args = [(load, r, alpha, gamma_mantle, strip) for strip in strips]

    pool = Pool(processes=n_procs)

    results = pool.map(_subside_grid_strip_helper, args)
    for dz, strip in results:
        start, stop = strip
        w[start:stop] += dz
