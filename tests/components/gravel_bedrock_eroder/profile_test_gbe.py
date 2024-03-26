import cProfile

import numpy as np

from landlab import HexModelGrid
from landlab.components import FlowAccumulator, GravelBedrockEroder


class BigantrSimulator:
    def __init__(
        self,
        nrows,
        ncols,
        dx,
        rand_elev,
        init_sed,
        ngrains,
        abrasion_coefs,
        gravel_fracs,
        pluck_coef,
        uplift_rate,
        update_interval=50000.0,
    ):
        self.grid = HexModelGrid((nrows, ncols), spacing=dx, node_layout="rect")

        self.grid.status_at_node[
            self.grid.perimeter_nodes
        ] = self.grid.BC_NODE_IS_CLOSED
        self.grid.status_at_node[
            self.grid.y_of_node == 0.0
        ] = self.grid.BC_NODE_IS_FIXED_VALUE

        self.elev = self.grid.add_zeros("topographic__elevation", at="node")
        self.rock = self.grid.add_zeros("bedrock__elevation", at="node")
        self.sed = self.grid.add_zeros("soil__depth", at="node")

        np.random.seed(1)
        self.rock[:] += rand_elev * np.random.rand(self.grid.number_of_nodes) - init_sed
        self.sed[:] = init_sed
        self.elev[:] = self.rock + self.sed

        self.flow_router = FlowAccumulator(self.grid)
        self.eroder = GravelBedrockEroder(
            grid=self.grid,
            plucking_coefficient=pluck_coef,
            number_of_sediment_classes=ngrains,
            coarse_fractions_from_plucking=gravel_fracs,
            abrasion_coefficients=abrasion_coefs,
        )

        self.uprate = uplift_rate
        self.current_time = 0.0

        self.update_interval = update_interval
        self.next_update = update_interval

    def update(self, dt):
        """Advance model by one step of size dt."""
        self.rock[self.grid.core_nodes] += self.uprate * dt
        self.elev[:] = self.rock + self.sed
        self.flow_router.run_one_step()
        self.eroder.run_one_step(dt)

    def update_until(self, to_time, dt):
        """Update model from current time to to_time in steps of dt."""
        remaining_time = to_time - self.current_time
        while remaining_time > 0.0:
            this_dt = min(dt, remaining_time)
            self.update(this_dt)
            remaining_time -= this_dt
        self.current_time = to_time


def main():
    model = BigantrSimulator(
        nrows=34,
        ncols=18,
        dx=1000.0,
        rand_elev=1.0,
        init_sed=10.0,
        ngrains=1,
        abrasion_coefs=0.0001,
        gravel_fracs=0.1,
        pluck_coef=0.01,
        uplift_rate=0.0001,
    )

    model.update_until(1.0e5, 100.0)


if __name__ == "__main__":
    cProfile.run("main()", sort="tottime")
