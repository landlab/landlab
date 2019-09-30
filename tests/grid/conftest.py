DIAGONAL_PROPERTIES = (
    "diagonals_at_node",
    "diagonal_dirs_at_node",
    "diagonal_adjacent_nodes_at_node",
    "d8_adjacent_nodes_at_node",
    "nodes_at_diagonal",
    "nodes_at_d8",
    "d8s_at_node",
    "d8_dirs_at_node",
    # 'd8_status_at_node',
    "length_of_diagonal",
    "length_of_d8",
    "status_at_diagonal",
    "diagonal_status_at_node",
    "active_diagonals",
    "active_diagonal_dirs_at_node",
    "status_at_d8",
    "active_d8",
    "active_d8_dirs_at_node",
)


def pytest_generate_tests(metafunc):
    if "diagonal_property" in metafunc.fixturenames:
        metafunc.parametrize("diagonal_property", DIAGONAL_PROPERTIES)
    elif "random_xy" in metafunc.fixturenames:
        from numpy.random import random_sample

        metafunc.parametrize(
            "random_xy",
            (
                tuple(-1e3 * random_sample(2)),
                tuple(1e3 * random_sample(2)),
                tuple(1e3 * (random_sample(2) - 0.5)),
            ),
        )
