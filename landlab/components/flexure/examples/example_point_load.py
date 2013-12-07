#! /usr/bin/env python
"""


"""

from landlab.components.flexure import FlexureComponent
from landlab import RasterModelGrid


def add_load_to_middle_of_grid(grid, load):
    shape = grid.shape

    load_array = grid.field_values(
        'node', 'lithosphere__overlying_pressure').view()
    load_array.shape = shape
    load_array[shape[0] / 2, shape[1] / 2] = load


def main():
    (n_rows, n_cols) = (100, 100)
    (dy, dx) = (10e3, 10e3)

    grid = RasterModelGrid(n_rows, n_cols, dx)

    flex = FlexureComponent(grid, method='flexure')

    add_load_to_middle_of_grid(grid, 1e9)

    flex.update()

    grid.imshow('node', 'lithosphere__elevation', symmetric_cbar=True,
                show=True) 


if __name__ == '__main__':
    main()
