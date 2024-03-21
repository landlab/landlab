import numpy as np

from landlab import RasterModelGrid
from landlab.grid.divergence import calc_flux_div_at_node as calc_flux_div_at_node_slow
from landlab.grid.gradients import (
    calc_diff_at_link as calc_diff_at_link_slow,
    calc_grad_at_link as calc_grad_at_link_slow,
)
from landlab.grid.raster_divergence import calc_flux_div_at_node
from landlab.grid.raster_gradients import calc_diff_at_link, calc_grad_at_link


class TimeCalc:
    params = [[], []]
    param_names = ["module", "function"]

    func = {
        "gradients": {
            "calc_diff_at_link": calc_diff_at_link_slow,
            "calc_grad_at_link": calc_grad_at_link_slow,
            "calc_flux_div_at_node": calc_flux_div_at_node_slow,
        },
        "raster": {
            "calc_diff_at_link": calc_diff_at_link,
            "calc_flux_div_at_node": calc_flux_div_at_node,
            "calc_grad_at_link": calc_grad_at_link,
        },
    }

    def time_calculation(self, mod_name, func_name):
        self.func[mod_name][func_name](self.grid, self.value_at_node, out=self.out)


class TimeCalcAtLink(TimeCalc):
    params = [["gradients", "raster"], ["calc_diff_at_link", "calc_grad_at_link"]]

    def setup(self, mod_name, func_name):
        self.grid = RasterModelGrid((400, 5000), (1.0, 2.0))
        self.value_at_node = np.random.uniform(size=self.grid.number_of_links)
        self.out = self.grid.empty(at="link")


class TimeCalcAtNode(TimeCalc):
    params = [["gradients", "raster"], ["calc_flux_div_at_node"]]

    def setup(self, mod_name, func_name):
        self.grid = RasterModelGrid((400, 5000), (1.0, 2.0))
        self.value_at_node = np.random.uniform(size=self.grid.number_of_links)
        self.out = self.grid.empty(at="node")
