import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

import landlab


def test_imshow_grid():
    rmg = landlab.RasterModelGrid((4, 5))

    pp = PdfPages("test.pdf")

    values = np.arange(rmg.number_of_nodes)
    landlab.plot.imshow_grid(rmg, values, values_at="node", limits=(0, 20))
    pp.savefig()

    plt.clf()
    rmg.status_at_node[7] = landlab.CLOSED_BOUNDARY
    values = np.arange(rmg.number_of_cells)
    landlab.plot.imshow_grid(rmg, values, values_at="cell", symmetric_cbar=True)
    pp.savefig()
    pp.close()
