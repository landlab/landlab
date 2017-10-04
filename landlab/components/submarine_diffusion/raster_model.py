#! /usr/bin/env python
import argparse

import yaml

from landlab import RasterModelGrid, CLOSED_BOUNDARY
from landlab.core import load_params
from landlab.plot import imshow_grid
from landlab.io.netcdf import write_raster_netcdf
from landlab.bmi.bmi_bridge import TimeStepper


class RasterModel(object):

    def __init__(self, grid=None, clock=None):
        self._clock = TimeStepper(**clock)
        self._grid = RasterModelGrid.from_dict(grid)

        self._components = ()

    @property
    def grid(self):
        return self._grid

    @property
    def clock(self):
        return self._clock

    def run_one_step(self, dt=None, output=None):
        """Run each component for one time step."""
        dt = dt or self.clock.step
        self.clock.advance(step=dt)

        for component in self._components:
            component.run_one_step(dt)

        if output:
            write_raster_netcdf(output, self.grid, append=True)

    def run(self, output=None):
        """Run the model until complete."""
        if output:
            write_raster_netcdf(output, self.grid, append=False)

        try:
            while 1:
                self.run_one_step(output=output)
        except StopIteration:
            pass

    @classmethod
    def argparser(cls, *args, **kwds):
        parser = argparse.ArgumentParser(*args, **kwds)

        parser.add_argument('file', nargs='?', help='model configuration file')
        parser.add_argument('--output', help='output file')
        parser.add_argument('--with-citations', action='store_true',
                            help='Print citations for components used')
        parser.add_argument('--verbose', action='store_true', help='be verbose')
        parser.add_argument('--dry-run', action='store_true',
                            help='do not actually run the model')
        parser.add_argument('--set', action='append', default=[],
                            help='set model parameters')
        parser.add_argument('--plot', choices=cls.LONG_NAME.keys(), default=None,
                            help='value to plot')

        return parser

    @classmethod
    def main(cls):
        parser = cls.argparser()

        args = parser.parse_args()

        # params = cls.DEFAULT_PARAMS
        # if args.file:
        #     params_from_file = load_params(args.file)
        #     for group in params.keys():
        #         params[group].update(params_from_file.get(group, {}))

        # params_from_cl = load_params_from_strings(args.set)
        # for group in params.keys():
        #     params[group].update(params_from_cl.get(group, {}))

        params = load_model_params(param_file=args.file,
                                   defaults=cls.DEFAULT_PARAMS,
                                   dotted_params=args.set)

        if args.verbose:
            print(yaml.dump(params, default_flow_style=False))

        model = cls(**params)

        if args.with_citations:
            from landlab.core.model_component import registry
            print(registry.format_citations())

        if not args.dry_run:
            model.run(output=args.output)

            if args.plot:
                imshow_grid(model.grid,
                            model.grid.at_node[cls.LONG_NAME[args.plot]],
                            at='node', show=True, cmap='Blues')


def load_params_from_strings(values):
    params = dict()
    for param in values:
        dotted_name, value = param.split('=')
        params.update(dots_to_dict(dotted_name, yaml.load(value)))

    return params


def dots_to_dict(name, value):
    base = {}
    level = base
    names = name.split('.')
    for k in names[:-1]:
        level[k] = dict()
        level = level[k]
    level[names[-1]] = value
    return base


def dict_to_dots(d):
    dots = []
    for names in walk_dict(d):
        dots.append('.'.join(names[:-1]) + '=' + str(names[-1]))
    return dots


def walk_dict(indict, prev=None):
    prev = prev[:] if prev else []

    if isinstance(indict, dict):
        for key, value in indict.items():
            if isinstance(value, dict):
                for d in walk_dict(value, [key] + prev):
                    yield d
            elif isinstance(value, list) or isinstance(value, tuple):
                yield prev + [key, value]
                # for v in value:
                #     for d in walk_dict(v, [key] + prev):
                #         yield d
            else:
                yield prev + [key, value]
    else:
        yield indict


def load_model_params(param_file=None, defaults=None, dotted_params=()):
    params = defaults or {}

    if param_file:
        params_from_file = load_params(param_file)
        dotted_params = dict_to_dots(params_from_file) + dotted_params
        # for group in params.keys():
        #     params[group].update(params_from_file.get(group, {}))

    params_from_cl = load_params_from_strings(dotted_params)
    for group in params.keys():
        params[group].update(params_from_cl.get(group, {}))

    return params
