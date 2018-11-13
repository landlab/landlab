from landlab import RasterModelGrid


def bench_gradient_across_faces():
    rmg = RasterModelGrid(1000, 1000)
    node_values = rmg.zeros()
    rmg.calc_grad_across_cell_faces(node_values)


def bench_gradient_across_corners():
    rmg = RasterModelGrid(1000, 1000)
    node_values = rmg.zeros()
    rmg.calc_grad_across_cell_corners(node_values)


def bench_max_gradients():
    rmg = RasterModelGrid(1000, 1000)
    node_values = rmg.zeros()
    (grads, nodes) = rmg.calculate_max_gradient_across_adjacent_cells(
        node_values, method="d8", return_node=True
    )
