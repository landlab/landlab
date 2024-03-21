import numpy as np

from landlab import RasterModelGrid
from landlab.components import Flexure


def parameterized(names, params):
    def decorator(func):
        func.param_names = names
        func.params = params
        return func

    return decorator


class TimeFlexure:
    """Change number of loads while maintaining grid size."""

    params = ([], [])
    param_names = ["grid_size", "n_loads"]

    def setup(self, size, n_loads):
        if n_loads > size:
            raise NotImplementedError()

        grid = RasterModelGrid((size, size), xy_spacing=10.0)
        grid.add_zeros("lithosphere__overlying_pressure_increment", at="node")

        self.loads = grid.zeros(at="node")
        self.loads[0 :: size // n_loads] = 1e9

        self.out = np.zeros((size, size))

        self.flex = Flexure(grid, method="flexure")

    def time_flexure(self, size, n_loads):
        self.flex.subside_loads(self.loads, out=self.out)


class FlexureOneLoad(TimeFlexure):
    params = (list(2 ** np.arange(5, 11)), [1])


class FlexureManyLoads(TimeFlexure):
    params = ([2**8], list(2 ** np.arange(0, 11)))
