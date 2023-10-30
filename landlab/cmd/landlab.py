import contextlib
import inspect
import os
import sys
import textwrap
from collections import defaultdict
from functools import partial

import numpy as np
import rich_click as click

from landlab import (
    FramedVoronoiGrid,
    HexModelGrid,
    ModelGrid,
    RadialModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
)

GRIDS = [
    ModelGrid,
    RasterModelGrid,
    VoronoiDelaunayGrid,
    HexModelGrid,
    RadialModelGrid,
    FramedVoronoiGrid,
]

CATEGORIES = {
    "boundary-condition",
    "connectivity",
    "deprecated",
    "field-add",
    "field-io",
    "gradient",
    "info-cell",
    "info-corner",
    "info-face",
    "info-field",
    "info-grid",
    "info-link",
    "info-node",
    "info-patch",
    "map",
    "quantity",
    "subset",
    "surface",
    "uncategorized",
}


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


@landlab.command(name="list")
def _list():
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


@landlab.group(chain=True)
@click.pass_context
def index(ctx):
    pass


@index.command()
@click.pass_context
def grids(ctx):
    verbose = ctx.parent.parent.params["verbose"]
    silent = ctx.parent.parent.params["silent"]

    index = {"grids": {}}
    for cls in GRIDS:
        index["grids"][cls.__name__] = _categorize_class(cls)
        index["grids"][cls.__name__]["field-io"] += [
            f"{cls.__module__}.{cls.__name__}.at_node",
            f"{cls.__module__}.{cls.__name__}.at_link",
            f"{cls.__module__}.{cls.__name__}.at_patch",
            f"{cls.__module__}.{cls.__name__}.at_corner",
            f"{cls.__module__}.{cls.__name__}.at_face",
            f"{cls.__module__}.{cls.__name__}.at_cell",
        ]

    print("")
    print("# Generated using `landlab index grids`")
    print("[grids]")
    for grid, cats in sorted(index["grids"].items()):
        print("")
        print(f"[grids.{grid}]")
        for cat, funcs in sorted(cats.items()):
            print(f"{cat} = [")
            print(
                textwrap.indent(
                    os.linesep.join([repr(f) + "," for f in sorted(funcs)]), "  "
                )
            )
            print("]")

    if verbose and not silent:
        summary = defaultdict(int)
        for cats in index["grids"].values():
            for cat, funcs in cats.items():
                summary[cat] += len(funcs)

        out("[summary]")
        out(f"grids = [{', '.join(sorted(index['grids']))}]")
        out(f"entries = {sum(summary.values())}")
        out("")
        out("[summary.categories]")
        for cat in sorted(summary):
            out(f"{cat} = {summary[cat]}")


@index.command()
@click.pass_context
def components(ctx):
    verbose = ctx.parent.parent.params["verbose"]
    silent = ctx.parent.parent.params["silent"]

    from sphinx.util.docstrings import prepare_docstring

    index = {"components": {}}
    for cls in get_all_components():
        if verbose and not silent:
            out(f"indexing: {cls.__name__}")
        index["components"][cls.__name__] = {
            "name": f"{cls.__module__}.{cls.__name__}",
            "unit_agnostic": cls._unit_agnostic,
            "info": cls._info,
            "summary": prepare_docstring(cls.__doc__)[0],
        }

    print("")
    print("# Generated using `landlab index components`")
    print("[components]")
    for component, info in sorted(index["components"].items()):
        print("")
        print(f"[components.{component}]")
        print(f"name = {info['name']!r}")
        print(f"unit_agnostic = {'true' if info['unit_agnostic'] else 'false'}")
        print(f"summary = {info['summary']!r}")

        for name, values in sorted(info["info"].items()):
            print("")
            print(f"[components.{component}.info.{name}]")
            print(f"doc = {values['doc']!r}")
            print(f"dtype = {str(np.dtype(values['dtype']))!r}")
            print(f"intent = {values['intent']!r}")
            print(f"mapping = {values['mapping']!r}")
            print(f"optional = {'true' if values['optional'] else 'false'}")
            print(f"units = {values['units']!r}")

    if not silent:
        out("[summary]")
        out(f"count = {len(index['components'])}")


@index.command()
@click.pass_context
def fields(ctx):
    verbose = ctx.parent.parent.params["verbose"]
    silent = ctx.parent.parent.params["silent"]

    fields = defaultdict(lambda: defaultdict(list))
    for cls in get_all_components():
        if verbose and not silent:
            out(f"checking {cls.__name__}... {len(cls._info)} fields")

        for name, desc in cls._info.items():
            fields[name]["desc"].append(desc["doc"])
            if desc["intent"].startswith("in"):
                fields[name]["used_by"].append(f"{cls.__module__}.{cls.__name__}")
            if desc["intent"].endswith("out"):
                fields[name]["provided_by"].append(f"{cls.__module__}.{cls.__name__}")

    print("")
    print("# Generated using `landlab index fields`")
    print("[fields]")
    for field, info in sorted(fields.items()):
        print("")
        print(f"[fields.{field}]")
        print(f"desc = {info['desc'][0]!r}")
        if info["used_by"]:
            # used_by = [repr(f) for f in info["used_by"]]
            # print(f"used_by = [{', '.join(used_by)}]")
            print("used_by = [")
            for component in sorted(info["used_by"]):
                print(f"  {component!r},")
            print("]")
        else:
            print("used_by = []")
        if info["provided_by"]:
            print("provided_by = [")
            for component in sorted(info["provided_by"]):
                print(f"  {component!r},")
            print("]")

        else:
            print("provided_by = []")

    if not silent:
        out("[summary]")
        out(f"count = {len(fields)}")


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
        with contextlib.suppress(TypeError):
            used += cls.input_var_names

    return used


def _provided_by(classes):
    """Get variables provided by components."""
    provided = []
    for cls in classes:
        with contextlib.suppress(TypeError):
            provided += cls.output_var_names

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


def _categorize_class(cls):
    funcs = {cat: [] for cat in CATEGORIES}

    for name, func in inspect.getmembers(cls):
        if not name.startswith("_"):
            full_name = ".".join([cls.__module__, cls.__name__, name])
            for cat in _extract_landlab_category(inspect.getdoc(func)):
                funcs[cat].append(full_name)
    return funcs


def _extract_landlab_category(s: str):
    from sphinx.util.docstrings import separate_metadata

    return [
        cat.strip() or "uncategorized"
        for cat in separate_metadata(s)[1].get("landlab", "").split(",")
    ]
