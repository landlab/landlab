from __future__ import print_function, absolute_import

import os
import sys


def get_all_components():
    from landlab.components import COMPONENTS
    from landlab.core.model_component import Component

    components = []
    for cls in COMPONENTS:
        if issubclass(cls, Component):
            components.append(cls)

    return components


def get_all_components_by_name():
    return dict([(cls.__name__, cls) for cls in get_all_components()])


def get_components(*args):
    """Get components by name.
    
    Parameters
    ----------
    names : list of str, optional
        Component names.

    Returns
    -------
    list of class
        Components with any of the given names.
    """
    if len(args) == 0 or len(args[0]) == 0:
        components = get_all_components()
    else:
        components_by_name = get_all_components_by_name()
        components = []
        for name in args[0]:
            try:
                components.append(components_by_name[name])
            except KeyError:
                print('{name}: not a component'.format(name=name),
                      file=sys.stderr)

    return components


def get_users_of(var):
    """Get components that use a variable."""
    users = []
    for cls in get_all_components():
        try:
            if var in cls.input_var_names:
                users.append(cls.__name__)
        except (AttributeError, TypeError):
            print('Warning: {name}: unable to get input vars'.format(
                name=cls.__name__), file=sys.stderr)

    return users

def get_providers_of(var):
    """Get components that provide a variable."""
    providers = []
    for cls in get_all_components():
        try:
            if var in cls.output_var_names:
                providers.append(cls.__name__)
        except (AttributeError, TypeError):
            print('Warning: {name}: unable to get output vars'.format(
                name=cls.__name__), file=sys.stderr)

    return providers


def used_by(classes):
    """Get variables used by components."""
    used = []
    for cls in classes:
        try:
            used += cls.input_var_names
        except TypeError:
            pass

    return used


def provided_by(classes):
    """Get variables provided by components."""
    provided = []
    for cls in classes:
        try:
            provided += cls.output_var_names
        except TypeError:
            pass

    return provided


def _test_input_var_names(cls):
    errors = []
    try:
        names = cls.input_var_names
    except AttributeError:
        errors.append('no input_var_names attribute')
    else:
        if not isinstance(names, tuple):
            errors.append('input_var_names is not a tuple')
    return errors


def _test_output_var_names(cls):
    errors = []
    try:
        names = cls.output_var_names
    except AttributeError:
        errors.append('no output_var_names attribute')
    else:
        if not isinstance(names, tuple):
            errors.append('output_var_names is not a tuple')
    return errors


def validate_component(cls):
    from landlab.core.model_component import Component

    errors = []
    if not issubclass(cls, Component):
        errors.append('not a subclass of Component')
    else:
        errors += _test_input_var_names(cls)
        errors += _test_output_var_names(cls)

    return errors


def validate(args):
    failures = 0
    classes = get_components(args.name)
    for cls in classes:
        errors = validate_component(cls)
        if errors:
            failures += 1
            for error in errors:
                print('Error: {name}: {error}'.format(
                    name=cls.__name__, error=error))

    return failures


def used(args):
    for name in used_by(get_components(args.name)):
        print('{name}'.format(name=name))


def provided(args):
    for name in provided_by(get_components(args.name)):
        print('{name}'.format(name=name))


def provides(args):
    for name in get_providers_of(args.var):
        print('{name}'.format(name=name))


def uses(args):
    for name in get_users_of(args.var):
        print('{name}'.format(name=name))


def list(args):
    for cls in get_all_components():
        print(cls.__name__)


def main():
    import argparse
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='subcommands',
                                       description='valid subcommands',
                                       help='additional help')
    parser_list = subparsers.add_parser('list',
                                        help='list landlab components')
    parser_list.set_defaults(func=list)

    parser_used_by = subparsers.add_parser(
        'used_by', help='variables used by a component')
    parser_used_by.add_argument('name', nargs='*', help='component name')
    parser_used_by.set_defaults(func=used)

    parser_provided_by = subparsers.add_parser(
        'provided_by', help='variables provided by a component')
    parser_provided_by.add_argument('name', nargs='*', help='component name')
    parser_provided_by.set_defaults(func=provided)

    parser_uses = subparsers.add_parser(
        'uses', help='landlab components that use a variable')
    parser_uses.add_argument('var', help='variable name')
    parser_uses.set_defaults(func=uses)

    parser_provides = subparsers.add_parser(
        'provides', help='landlab components that provide a variable')
    parser_provides.add_argument('var', help='variable name')
    parser_provides.set_defaults(func=provides)

    parser_validate = subparsers.add_parser(
        'validate', help='validate a landlab component')
    parser_validate.add_argument('name', nargs='*', help='component name')
    parser_validate.set_defaults(func=validate)

    args = parser.parse_args()
    rtn = args.func(args)
    if rtn is not None and rtn != 0:
        parser.exit(status=1)
