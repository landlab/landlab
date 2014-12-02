#! /usr/bin/env python

import sys, os

sys.path.pop(0)

from optparse import OptionParser

parser = OptionParser('usage: %prog [options] -- [nosetests options]')
parser.add_option('-v', '--verbose', action='count', dest='verbose',
                  default=1, help='increase verbosity [%default]')
parser.add_option('--doctests', action='store_true', dest='doctests',
                  default=False, help='Run doctests in module [%default]')
parser.add_option('--coverage', action='store_true', dest='coverage',
                  default=False, help='report coverage of landlab [%default]')
parser.add_option('-m', '--mode', action='store', dest='mode', default='fast',
                  help='"fast", "full", or something that can be passed to '
                       'nosetests -A [%default]')

(options, args) = parser.parse_args()


import landlab


result = landlab.test(label=options.mode, verbose=options.verbose,
                      doctests=options.doctests, coverage=options.coverage,
                      extra_argv=args)

if result.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)
