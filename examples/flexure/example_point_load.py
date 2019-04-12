#! /usr/bin/env python
"""


"""

from landlab import RasterModelGrid
from landlab.components.flexure import Flexure


def add_load_to_middle_of_grid(grid, load):
    shape = grid.shape

    load_array = grid.field_values(
        "node", "lithosphere__overlying_pressure_increment"
    ).view()
    load_array.shape = shape
    load_array[shape[0] / 2, shape[1] / 2] = load


def main():
    (n_rows, n_cols) = (100, 100)
    spacing = (10e3, 10e3)

    grid = RasterModelGrid(n_rows, n_cols, spacing[1])

    flex = Flexure(grid, method="flexure")

    add_load_to_middle_of_grid(grid, 1e7)

    flex.update()

    grid.imshow(
        "node",
        "lithosphere_surface__elevation_increment",
        symmetric_cbar=True,
        show=True,
    )


if __name__ == "__main__":
    main()
