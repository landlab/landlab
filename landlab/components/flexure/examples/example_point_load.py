#! /usr/bin/env python

import pylab

from landlab.components.flexure import Flexure


def main():

    shape = (100, 100)
    spacing = (10e3, 10e3)

    flex = Flexure(shape, spacing, (0., 0.))

    load = flex['lithosphere__overlying_pressure']
    load[shape[0]/2, shape[1]/2] = 1e9

    dz = flex['lithosphere__elevation']

    flex.update()

    dz.shape = shape

    pylab.imshow(dz)
    pylab.colorbar()
    #pylab.plot(dz[shape[0]/2,:])
    pylab.show()


if __name__ == '__main__':
    main()
