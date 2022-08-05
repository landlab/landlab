import os
import sys
from functools import partial

import rich_click as click
from jinja2 import Environment

from landlab import (
    ModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
    HexModelGrid,
    RadialModelGrid,
)
from landlab.core.utils import (
    # CATEGORIES,
    get_category_from_class,
    # get_categories_from_class,
    get_funcs_by_category,
)


GRIDS = [
    ModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
    HexModelGrid,
    RadialModelGrid,
]


click.rich_click.ERRORS_SUGGESTION = (
    "Try running the '--help' flag for more information."
)
click.rich_click.ERRORS_EPILOGUE = (
    "To find out more, visit https://github.com/landlab/landlab"
)
click.rich_click.STYLE_ERRORS_SUGGESTION = "yellow italic"
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = False
click.rich_click.SHOW_METAVARS_COLUMN = True
click.rich_click.USE_MARKDOWN = True

out = partial(click.secho, bold=True, file=sys.stderr)
err = partial(click.secho, fg="red", file=sys.stderr)


@click.group()  # chain=True)
@click.version_option()
@click.option(
    "--cd",
    default=".",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="chage to directory, then execute",
)
@click.option(
    "-s",
    "--silent",
    is_flag=True,
    help="Suppress status status messages, including the progress bar.",
)
@click.option(
    "-v", "--verbose", is_flag=True, help="Also emit status messages to stderr."
)
def landlab(cd, silent, verbose) -> None:
    os.chdir(cd)


@landlab.command()
@click.argument("category", type=str)
@click.option("--title", help="section title", default=None)
@click.option(
    "-t", "--template", help="template file", type=click.File(mode="r"), default="-"
)
def render(category, title, template):
    title = title or category

    # env = Environment(loader=FileSystemLoader("templates"))
    # template = env.get_template(template)

    _template = Environment().from_string(template.read())

    grids = {
        ".".join([grid.__module__, grid.__name__]): get_category_from_class(
            grid, category
        )
        for grid in GRIDS
    }

    print(_template.render(title=title, category=category, grids=grids))


@landlab.command()
@click.argument("grid", type=str)
@click.option("--title", help="section title", default=None)
@click.option(
    "-t", "--template", help="template file", type=click.File(mode="r"), default="-"
)
def render_grid(grid, title, template):
    title = title or grid

    # env = Environment(loader=FileSystemLoader("templates"))
    # template = env.get_template(template)

    _template = Environment().from_string(template.read())

    # categories = {category: get_funcs_by_category(grid) for category in CATEGORIES}
    categories = get_funcs_by_category(RasterModelGrid)

    print(_template.render(title=title, grid=grid, categories=categories))


@landlab.command()
@click.argument("category", nargs=-1)
@click.option("--with-prefix", is_flag=True)
def filter(category, with_prefix):
    for _category in category:
        for grid in GRIDS:
            for func in get_category_from_class(grid, _category):
                prefix = f"{grid.__module__}.{grid.__name__}." if with_prefix else ""
                print(f"{prefix}{func}")


@landlab.command()
def list():
    for cls in get_all_components():
        print(cls.__name__)


@landlab.command()
@click.argument("component", type=str, nargs=-1)
def used_by(component):
    for name in _used_by(get_components(component)):
        print(name)


@landlab.command()
@click.argument("component", type=str, nargs=-1)
def provided_by(component):
    for name in _provided_by(get_components(component)):
        print(name)


@landlab.command()
@click.argument("var", type=str)
def uses(var):
    for name in get_users_of(var):
        print(name)


@landlab.command()
@click.argument("var", type=str)
def provides(var):
    for name in get_providers_of(var):
        print(name)


@landlab.command()
@click.argument("component", type=str, nargs=-1)
def validate(component):
    failures = 0
    classes = get_components(component)
    for cls in classes:
        out(cls.__name__)
        errors = _validate_component(cls)
        if errors:
            failures += 1
            for error in errors:
                err(f"Error: {cls.__name__}: {error}")
    if failures:
        click.Abort()
    else:
        out("💥 All good! 💥")


def get_all_components():
    from landlab.components import COMPONENTS
    from landlab.core.model_component import Component

    components = []
    for cls in COMPONENTS:
        if issubclass(cls, Component):
            components.append(cls)

    return components


def get_all_components_by_name():
    return {cls.__name__: cls for cls in get_all_components()}


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
                print(f"{name}: not a component", file=sys.stderr)

    return components


def get_users_of(var):
    """Get components that use a variable."""
    users = []
    for cls in get_all_components():
        try:
            if var in cls.input_var_names:
                users.append(cls.__name__)
        except (AttributeError, TypeError):
            print(
                f"Warning: {cls.__name__}: unable to get input vars",
                file=sys.stderr,
            )

    return users


def get_providers_of(var):
    """Get components that provide a variable."""
    providers = []
    for cls in get_all_components():
        try:
            if var in cls.output_var_names:
                providers.append(cls.__name__)
        except (AttributeError, TypeError):
            print(
                f"Warning: {cls.__name__}: unable to get output vars",
                file=sys.stderr,
            )

    return providers


def _used_by(classes):
    """Get variables used by components."""
    used = []
    for cls in classes:
        try:
            used += cls.input_var_names
        except TypeError:
            pass

    return used


def _provided_by(classes):
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
        errors.append("no input_var_names attribute")
    else:
        if not isinstance(names, tuple):
            errors.append("input_var_names is not a tuple")
    return errors


def _test_output_var_names(cls):
    errors = []
    try:
        names = cls.output_var_names
    except AttributeError:
        errors.append("no output_var_names attribute")
    else:
        if not isinstance(names, tuple):
            errors.append("output_var_names is not a tuple")
    return errors


def _validate_component(cls):
    from landlab.core.model_component import Component

    errors = []
    if not issubclass(cls, Component):
        errors.append("not a subclass of Component")
    else:
        errors += _test_input_var_names(cls)
        errors += _test_output_var_names(cls)

    return errors


def _validate(args):
    failures = 0
    classes = get_components(args.name)
    for cls in classes:
        errors = _validate_component(cls)
        if errors:
            failures += 1
            for error in errors:
                print(f"Error: {cls.__name__}: {error}")

    return failures
