#! /usr/bin/env python
"""Calculate slope aspects on a `RasterModelGrid`."""


def _one_line_slopes(input_array, grid, vals):
    node = input_array[0]
    diagonals = input_array[5:]
    neighbors = input_array[1:5]

    if not grid.status_at_node[node] == 0:
        raise IndexError("One or more of the provided nodes was closed!")

    try:
        slope_we = (
            (vals[diagonals[1]] + 2.0 * vals[neighbors[2]] + vals[diagonals[2]])
            - (vals[diagonals[0]] + 2.0 * vals[neighbors[0]] + vals[diagonals[3]])
        ) / (8.0 * grid.dx)
        slope_sn = (
            (vals[diagonals[2]] + 2.0 * vals[neighbors[3]] + vals[diagonals[3]])
            - (vals[diagonals[1]] + 2.0 * vals[neighbors[:, 1]] + vals[diagonals[0]])
        ) / (8.0 * grid.dy)
        return slope_we, slope_sn
    except IndexError:
        C = vals[node]
        weighting_verticals = 4.0
        weighting_horizontals = 4.0
        try:
            vertical_grad = (vals[neighbors[3]] - vals[neighbors[1]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                vertical_grad = (C - vals[neighbors[1]]) / grid.dy
            except IndexError:
                try:
                    vertical_grad = (vals[neighbors[3]] - C) / grid.dy
                except IndexError:
                    vertical_grad = 0.0
                    weighting_verticals -= 2.0
        try:
            horizontal_grad = (vals[neighbors[2]] - vals[neighbors[0]]) / (
                2.0 * grid.dx
            )
        except IndexError:
            try:
                horizontal_grad = (C - vals[neighbors[0]]) / grid.dx
            except IndexError:
                try:
                    horizontal_grad = (vals[neighbors[2]] - C) / grid.dx
                except IndexError:
                    horizontal_grad = 0.0
                    weighting_horizontals -= 2.0
        try:
            left_grad = (vals[diagonals[2]] - vals[diagonals[1]]) / (2.0 * grid.dx)
        except IndexError:
            try:
                C = vals[neighbors[2]]
            except IndexError:
                left_grad = 0.0
                weighting_verticals -= 1.0
            else:
                try:
                    left_grad = (C - vals[diagonals[1]]) / grid.dx
                except IndexError:
                    left_grad = (vals[diagonals[2]] - C) / grid.dx
        try:
            right_grad = (vals[diagonals[3]] - vals[diagonals[0]]) / (2.0 * grid.dx)
        except IndexError:
            try:
                C = vals[neighbors[0]]
            except IndexError:
                right_grad = 0.0
                weighting_verticals -= 1.0
            else:
                try:
                    right_grad = (C - vals[diagonals[0]]) / grid.dx
                except IndexError:
                    right_grad = (vals[diagonals[3]] - C) / grid.dx
        try:
            top_grad = (vals[diagonals[1]] - vals[diagonals[0]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                C = vals[neighbors[1]]
            except IndexError:
                top_grad = 0.0
                weighting_horizontals -= 1.0
            else:
                try:
                    top_grad = (C - vals[diagonals[0]]) / grid.dy
                except IndexError:
                    top_grad = (vals[diagonals[1]] - C) / grid.dy
        try:
            bottom_grad = (vals[diagonals[2]] - vals[diagonals[3]]) / (2.0 * grid.dy)
        except IndexError:
            try:
                C = vals[neighbors[3]]
            except IndexError:
                bottom_grad = 0.0
                weighting_horizontals -= 1.0
            else:
                try:
                    bottom_grad = (C - vals[diagonals[3]]) / grid.dy
                except IndexError:
                    bottom_grad = (vals[diagonals[2]] - C) / grid.dy

        slope_we = (
            top_grad + 2.0 * horizontal_grad + bottom_grad
        ) / weighting_horizontals
        slope_sn = (left_grad + 2.0 * vertical_grad + right_grad) / weighting_verticals

        return slope_we, slope_sn
