from multiprocessing import Pool

import numpy as np
cimport numpy as np
cimport cython


_RHO_MANTLE = 3300.
_GRAVITY = 9.81


DTYPE = np.double
ctypedef np.double_t DTYPE_t

@cython.boundscheck(False)
def subside_parallel_row(np.ndarray[DTYPE_t, ndim=1] w,
                         np.ndarray[DTYPE_t, ndim=1] load,
                         np.ndarray[DTYPE_t, ndim=1] r,
                         alpha):
  cdef int ncols = w.size
  cdef double inv_c = 1. / (2. * np.pi * _RHO_MANTLE * _GRAVITY * alpha ** 2.)
  cdef double c
  cdef int i
  cdef int j

  for i in xrange(ncols):
    if abs(load[i]) > 1e-6:
      c = load[i] * inv_c
      for j in xrange(ncols):
        w[j] += - c * r[abs(j - i)]


def subside_grid(np.ndarray[DTYPE_t, ndim=2] w,
                 np.ndarray[DTYPE_t, ndim=2] load,
                 np.ndarray[DTYPE_t, ndim=2] r,
                 alpha):
  cdef int nrows = w.shape[0]
  cdef i, j

  for i in xrange(nrows):
    for j in xrange(nrows):
      subside_parallel_row(w[j], load[i], r[abs(j - i)], alpha)


def subside_grid_strip(np.ndarray[DTYPE_t, ndim=2] load,
                       np.ndarray[DTYPE_t, ndim=2] r,
                       alpha, strip_range):
  (start, stop) = strip_range

  cdef np.ndarray w = np.zeros((stop - start, load.shape[1]), dtype=DTYPE)
  cdef i

  for i in xrange(load.shape[0]):
    for j in xrange(start, stop):
      subside_parallel_row(w[j - start], load[i], r[abs(j - i)], alpha)

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
                             alpha, n_procs):
    if n_procs == 1:
        return subside_grid(w, load, r, alpha)

    strips = tile_grid_into_strips(w, n_procs)

    args = [(load, r, alpha, strip) for strip in strips]

    pool = Pool(processes=n_procs)

    results = pool.map(_subside_grid_strip_helper, args)
    for dz, strip in results:
        start, stop = strip
        w[start:stop] += dz
