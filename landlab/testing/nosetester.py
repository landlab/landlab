from __future__ import print_function

import os
import sys
import doctest
from numpy.testing import Tester
from nose.plugins import doctests
from nose.plugins.base import Plugin


def show_system_info():
    import landlab
    import nose

    print('landlab version %s' % landlab.__version__)
    landlab_dir = os.path.dirname(landlab.__file__)
    print('landlab is installed in %s' % landlab_dir)

    print('Python version %s' % sys.version.replace('\n', ''))
    print('nose version %d.%d.%d' % nose.__versioninfo__)


class LandlabDoctest(doctests.Doctest):
    name = 'landlabdoctest'
    score = 1000
    doctest_ignore = ('setup.py', )

    def options(self, parser, env=os.environ):
        Plugin.options(self, parser, env)
        self.doctest_tests = True
        self._doctest_result_var = None

    def configure(self, options, config):
        options.doctestOptions = ["+ELLIPSIS", "+NORMALIZE_WHITESPACE"]
        super(LandlabDoctest, self).configure(options, config)


class LandlabTester(Tester):
    excludes = ['examples']

    def __init__(self, package=None, raise_warnings='develop'):
        package_name = None
        if package is None:
            f = sys._getframe(1)
            package_path = f.f_locals.get('__file__', None)
            if package_path is None:
                raise AssertionError
            package_path = os.path.dirname(package_path)
            package_name = f.f_locals.get('__name__', None)
        elif isinstance(package, type(os)):
            package_path = os.path.dirname(package.__file__)
            package_name = getattr(package, '__name__', None)
        else:
            package_path = str(package)

        self.package_path = os.path.abspath(package_path)

        # Find the package name under test; this name is used to limit coverage
        # reporting (if enabled).
        if package_name is None:
            #package_name = get_package_name(package_path)
            package_name = 'landlab'
        self.package_name = package_name

        # Set to "release" in constructor in maintenance branches.
        self.raise_warnings = raise_warnings

    def _get_custom_doctester(self):
        return LandlabDoctest()

    def test(self, **kwds):
        kwds.setdefault('verbose', 2)
        kwds.setdefault('doctests', 1)
        kwds.setdefault('coverage', 1)
        kwds.setdefault('extra_argv', ['-x'])
        show_system_info()
        return super(LandlabTester, self).test(**kwds)
