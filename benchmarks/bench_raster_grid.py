from landlab import RasterModelGrid


class TimeRasterModelGrid:
    param_names = ["number-of-rows", "number-of-columns"]
    params = [
        [128, 256, 512, 1024],
        [128, 256, 512, 1024],
    ]

    def time_creation(self, n_rows, n_cols):
        RasterModelGrid((n_rows, n_cols))
