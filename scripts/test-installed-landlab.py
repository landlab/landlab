#! /usr/bin/env python
from __future__ import print_function

import sys, os

sys.path.pop(0)

from optparse import OptionParser

parser = OptionParser('usage: %prog [options] -- [nosetests options]')
parser.add_option('-v', '--verbose', action='count', dest='verbose',
                  default=1, help='increase verbosity [%default]')
parser.add_option('--no-doctests', action='store_false', dest='doctests',
                  default=True,
                  help='Do not run doctests in module [%default]')
parser.add_option('--coverage', action='store_true', dest='coverage',
                  default=False, help='report coverage of landlab [%default]')
parser.add_option('-m', '--mode', action='store', dest='mode', default='fast',
                  help='"fast", "full", or something that can be passed to '
                       'nosetests -A [%default]')

(options, args) = parser.parse_args()

try:
    import landlab
except ImportError:
    print('Unable to import landlab. You may not have landlab installed.')
    print('Here is your sys.path')
    print(os.linesep.join(sys.path))
    raise


result = landlab.test(label=options.mode, verbose=options.verbose,
                      doctests=options.doctests, coverage=options.coverage,
                      extra_argv=args, raise_warnings='release')

if result.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)
