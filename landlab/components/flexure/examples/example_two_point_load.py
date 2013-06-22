#! /usr/bin/env python

import pylab

from landlab.components.flexure import Flexure

_SHAPE = (100, 100)
_SPACING = (10e3, 10e3)
_LOAD_LOCS = [
    (_SHAPE[0]/2, _SHAPE[1]/2),
    (_SHAPE[0]/4, _SHAPE[1]/4),
]


def main():

    flex = Flexure(_SHAPE, _SPACING, (0., 0.))

    load = flex['lithosphere__overlying_pressure']
    for loc in _LOAD_LOCS:
        load[loc] = 1e9

    dz = flex['lithosphere__elevation']

    flex.update()

    dz.shape = _SHAPE

    #pylab.imshow(dz)
    #pylab.colorbar()
    pylab.plot(dz[_SHAPE[0] * 3 / 8, :])
    pylab.show()


if __name__ == '__main__':
    main()
