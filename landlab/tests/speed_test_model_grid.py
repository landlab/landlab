"""
simple script to run speed tests of various functions in model grid
"""

from __future__ import print_function
import landlab as ll
import time


def main():
    mg = ll.RasterModelGrid(20, 30, 1.0)

    nt = 1000

    s = mg.zeros(centering='node')
    g = mg.zeros(centering='active_link')
    divg = mg.zeros(centering='node')

    start_time = time.time()

    for i in range(nt):

        g = mg.calc_grad_of_active_link(s, g)

    time1 = time.time()

    time2 = time.time()

    for i in range(nt):

        divg = mg.calculate_flux_divergence_at_nodes(g, divg)

    time3 = time.time()

    for i in range(nt):

        divg = mg.calculate_flux_divergence_at_nodes_slow(g, divg)

    time4 = time.time()

    for i in range(nt):

        divg = mg.calculate_flux_divergence_at_active_cells(g)

    time5 = time.time()

    for i in range(nt):

        divg = mg.calculate_flux_divergence_at_active_cells_slow(g)

    time6 = time.time()

    print('Elapsed time with fast gradient algo: ' + str(time1 - start_time))
    print('Elapsed time with slow gradient algo: ' + str(time2 - time1))
    print('Elapsed time with fast node-divergence algo: ' + str(time3 - time2))
    print('Elapsed time with slow node-divergence algo: ' + str(time4 - time3))
    print('Elapsed time with fast activecell-divergence algo: ' + str(time5 - time4))
    print('Elapsed time with slow activecell-divergence algo: ' + str(time6 - time5))


if __name__ == '__main__':
    main()
