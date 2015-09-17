#! /usr/env/python
"""
just a little script for testing various bits and pieces of model_grid
"""
from __future__ import print_function
from landlab import RasterModelGrid


def main():
    nr = 3
    nc = 4
    ncells = nr * nc
    mg = RasterModelGrid(nr, nc, 1.0)

    for i in xrange(0, ncells):
        if i > 0:
            fid = mg.get_face_connecting_cell_pair(i, i - 1)
            print('Cells', i, 'and', i - 1, 'are connected by face', fid)
        if i < (ncells - 1):
            fid = mg.get_face_connecting_cell_pair(i, i + 1)
            print('Cells', i, 'and', i + 1, 'are connected by face', fid)
        if i > (nc - 1):
            fid = mg.get_face_connecting_cell_pair(i, i - nc)
            print('Cells', i, 'and', i - nc, 'are connected by face', fid)
        if i < (ncells - nc):
            fid = mg.get_face_connecting_cell_pair(i, i + nc)
            print('Cells', i, 'and', i + nc, 'are connected by face', fid)

    my_bc = mg.create_boundary_condition()
    print('I created a BC that looks like:')
    print(my_bc)
    for i in arange(0, size(my_bc.boundary_code)):
        print((i, my_bc.boundary_code[i], my_bc.boundary_gradient[i],
               my_bc.tracks_cell[i]))


if __name__ == '__main__':
    main()
