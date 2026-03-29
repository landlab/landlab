import argparse
import contextlib
import inspect
import itertools
import os
import pathlib
import re
import sys
import textwrap
from collections import Counter
from collections import defaultdict
from collections.abc import Iterable
from functools import partial

import numpy as np

from landlab import FramedVoronoiGrid
from landlab import HexModelGrid
from landlab import ModelGrid
from landlab import RadialModelGrid
from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid
from landlab.cmd.authors import AuthorList
from landlab.cmd.authors import AuthorsConfig
from landlab.cmd.authors import AuthorsSubprocessError
from landlab.cmd.authors import GitLog

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


out = partial(print, file=sys.stderr)
err = partial(print, file=sys.stderr)


def _abort() -> None:
    raise SystemExit(1)


def _is_verbose(args) -> bool:
    return getattr(args, "verbose", False)


def _is_silent(args) -> bool:
    return getattr(args, "silent", False)


def _authors_config(args) -> AuthorsConfig:
    params = {}
    for name in ("authors_file", "credits_file"):
        value = getattr(args, name, None)
        if value:
            params[name] = value
    return AuthorsConfig(**params)


def _list(args=None):
    for cls in get_all_components():
        print(cls.__name__)


def used_by(args):
    for name in _used_by(get_components(args.component)):
        print(name)


def provided_by(args):
    for name in _provided_by(get_components(args.component)):
        print(name)


def uses(args):
    for name in get_users_of(args.var):
        print(name)


def provides(args):
    for name in get_providers_of(args.var):
        print(name)


def validate(args):
    failures = 0
    classes = get_components(args.component)
    for cls in classes:
        out(cls.__name__)
        errors = _validate_component(cls)
        if errors:
            failures += 1
            for error in errors:
                err(f"Error: {cls.__name__}: {error}")
    if failures:
        _abort()
    else:
        out("💥 All good! 💥")


def authors_create(args):
    """Create a database of contributors."""
    verbose = _is_verbose(args)
    silent = _is_silent(args)
    config = _authors_config(args)
    credits_file = pathlib.Path(config["credits_file"])

    if verbose and not silent:
        config_str = textwrap.indent(str(config), prefix="  ")
        out("using the following configuration:")
        out(f"{config_str}")

    git_log = GitLog("%an, %ae")
    try:
        names_and_emails = git_log()
    except AuthorsSubprocessError as error:
        err(error)
        raise SystemExit(1) from error
    else:
        if verbose and not silent:
            out(f"{git_log}")

    if not silent and args.update_existing:
        if not credits_file.is_file():
            err(f"nothing to update ({credits_file})")
        else:
            out(f"updating existing author credits ({credits_file})")

    authors = (
        AuthorList.from_toml(credits_file)
        if args.update_existing and credits_file.is_file()
        else AuthorList()
    )

    authors.update(AuthorList.from_csv(names_and_emails))
    lines = [author.to_toml() for author in sorted(authors, key=lambda item: item.name)]

    print((2 * os.linesep).join(lines))


def authors_build(args):
    """Build an authors file."""
    verbose = _is_verbose(args)
    silent = _is_silent(args)
    config = _authors_config(args)

    if verbose and not silent:
        config_str = textwrap.indent(str(config), prefix="  ")
        out("using the following configuration:")
        out(f"{config_str}")

    exclude = args.exclude
    authors_file = pathlib.Path(args.authors_file)
    author_format = args.author_format
    credits_file = pathlib.Path(config["credits_file"])
    start_string = args.start_string

    git_log = GitLog("%aN")
    try:
        commit_authors = git_log()
    except AuthorsSubprocessError as error:
        err(error)
        raise SystemExit(1) from error
    else:
        if verbose and not silent:
            out(f"{git_log}")

    intro = (
        _read_until(authors_file, until=start_string) if authors_file.is_file() else ""
    )
    if len(intro) == 0:
        err(f"empty or missing authors file ({authors_file})")

    authors = (
        AuthorList.from_toml(credits_file) if credits_file.is_file() else AuthorList()
    )

    if len(authors) == 0:
        err(f"missing or empty credits file ({credits_file})")

    commits = Counter()
    for author in commit_authors.splitlines():
        canonical_name = authors.find_author(author.strip()).name
        commits[canonical_name] += 1

    lines = [intro]
    for author in sorted(authors, key=lambda a: commits[a.name], reverse=True):
        github = _guess_github_user(author)
        if github is None:
            github = "landlab"
        author.github = github
        if not exclude_matches_any(author.names, exclude):
            lines.append(
                author_format.format(
                    name=author.name, github=author.github, email=author.email
                )
            )

    print(os.linesep.join(lines))


def _read_until(path_to_file, until=None):
    """Read a file until a line starting with a given string.

    Parameters
    ----------
    path_to_file : str or path-like
        The file to read.
    until : str, optional
        Read lines until reaching a line that starts with ``until``.
        If not provided, read the entire file.

    Returns
    -------
    str
        The contents of the file up to, and including, the search
        string.
    """
    with open(path_to_file) as fp:
        if until is None:
            return fp.read()

        lines = []
        for line in fp.readlines():
            if line.startswith(until):
                lines.append(line)
                break
            lines.append(line)
    return "".join(lines)


def _guess_github_user(author):
    """Guess an author's github username."""
    github = None

    try:
        github = author.github
    except AttributeError:
        pass

    try:
        github = author._extras["github"]
    except KeyError:
        pass

    if github is None:
        for email in sorted(author.emails):
            if email.endswith("github.com"):
                github, _ = email.split("@")
                if "+" in github:
                    _, github = github.split("+")
                break

    if github is None:
        for name in sorted(author.names):
            if name.isalnum():
                github = name
                break

    return github


def authors_mailmap(args):
    """Create a mailmap file from an author list."""
    verbose = _is_verbose(args)
    silent = _is_silent(args)
    config = _authors_config(args)
    credits_file = config["credits_file"]

    if verbose and not silent:
        config_str = textwrap.indent(str(config), prefix="  ")
        out("using the following configuration:")
        out(f"{config_str}")
        out(f"reading author list: {credits_file}")

    print(
        textwrap.dedent(
            """
            # Prevent git from showing duplicate names with commands like "git shortlog"
            # See the manpage of git-shortlog for details.
            # The syntax is:
            #
            #   Name that should be used <email that should be used> Bad name <bad email>
            #
            # You can skip Bad name if it is the same as the one that should be used,
            # and is unique.
            #
            # This file is up-to-date if the command,
            #
            #   git log --format="%aN <%aE>" | sort -u
            #
            # gives no duplicates.
            """
        ).lstrip()
    )
    authors = AuthorList.from_toml(credits_file)
    for author in authors:
        good_name, good_email = author.name, author.email
        for bad_name, bad_email in itertools.product(
            sorted(author.names), sorted(author.emails)
        ):
            print(f"{good_name} <{good_email}> {bad_name} <{bad_email}>")


def authors_list(args):
    verbose = _is_verbose(args)
    silent = _is_silent(args)
    exclude = args.exclude

    git_log = GitLog("%aN")
    try:
        commit_authors = git_log()
    except AuthorsSubprocessError as error:
        err(error)
        raise SystemExit(1) from error
    else:
        if verbose and not silent:
            out(f"{git_log}")

    if verbose and not silent:
        out(f"reading authors from {args.file}")
    authors = AuthorList.from_toml(args.file)

    commits = Counter()
    for author in commit_authors.splitlines():
        canonical_name = authors.find_author(author.strip()).name
        commits[canonical_name] += 1

    for name, n_commits in sorted(commits.items(), key=lambda x: x[1], reverse=True):
        author = authors.find_author(name)
        if not exclude_matches_any(author.names, exclude):
            print(f"{author.name} <{author.email}> ({n_commits})")


def exclude_matches_any(names: Iterable[str], exclude: str):
    exclude_re = re.compile(exclude)
    for name in names:
        if exclude_re.search(name):
            return True
    return False


def index_grids(args):
    verbose = _is_verbose(args)
    silent = _is_silent(args)

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
        summary = Counter()
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


def index_components(args):
    verbose = _is_verbose(args)
    silent = _is_silent(args)

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


def index_fields(args):
    verbose = _is_verbose(args)
    silent = _is_silent(args)

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


def get_components(args):
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
    if len(args) == 0:
        components = get_all_components()
    else:
        components_by_name = get_all_components_by_name()
        components = []
        for name in args:
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="landlab",
        epilog="To find out more, visit https://github.com/landlab/landlab",
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s",
    )
    parser.add_argument(
        "--cd",
        default=".",
        type=pathlib.Path,
        help="change to directory, then execute",
    )
    parser.add_argument(
        "-s",
        "--silent",
        action="store_true",
        help="Suppress status status messages, including the progress bar.",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="Also emit status messages to stderr.",
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    p = subparsers.add_parser("list")
    p.set_defaults(func=_list)

    p = subparsers.add_parser("used_by")
    p.add_argument("component", nargs="*", type=str)
    p.set_defaults(func=used_by)

    p = subparsers.add_parser("provided_by")
    p.add_argument("component", nargs="*", type=str)
    p.set_defaults(func=provided_by)

    p = subparsers.add_parser("uses")
    p.add_argument("var", type=str)
    p.set_defaults(func=uses)

    p = subparsers.add_parser("provides")
    p.add_argument("var", type=str)
    p.set_defaults(func=provides)

    p = subparsers.add_parser("validate")
    p.add_argument("component", nargs="*", type=str)
    p.set_defaults(func=validate)

    authors = subparsers.add_parser(
        "authors", description="Commands for working with lists of authors."
    )
    authors.add_argument("--authors-file", type=str, help="existing authors file")
    authors.add_argument(
        "--credits-file",
        default=".credits.toml",
        type=str,
        help="The file that contains a list of authors",
    )
    authors_subparsers = authors.add_subparsers(dest="authors_command", required=True)

    p = authors_subparsers.add_parser("create")
    p.add_argument(
        "--update-existing", dest="update_existing", action="store_true", default=True
    )
    p.add_argument("--no-update-existing", dest="update_existing", action="store_false")
    p.set_defaults(func=authors_create)

    p = authors_subparsers.add_parser("build")
    p.add_argument("--exclude", required=True)
    p.add_argument("--author-format", required=True)
    p.add_argument("--start-string", required=True)
    p.set_defaults(func=authors_build)

    p = authors_subparsers.add_parser("mailmap")
    p.set_defaults(func=authors_mailmap)

    p = authors_subparsers.add_parser("list")
    p.add_argument(
        "--file",
        default="authors.toml",
        type=str,
        help="existing authors file",
    )
    p.add_argument("--exclude", required=True)
    p.set_defaults(func=authors_list)

    index = subparsers.add_parser("index")
    index_subparsers = index.add_subparsers(dest="index_command", required=True)

    p = index_subparsers.add_parser("grids")
    p.set_defaults(func=index_grids)

    p = index_subparsers.add_parser("components")
    p.set_defaults(func=index_components)

    p = index_subparsers.add_parser("fields")
    p.set_defaults(func=index_fields)

    return parser


def main(args=None) -> int:
    parser = build_parser()
    ns = parser.parse_args(args)

    os.chdir(ns.cd)
    ns.func(ns)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
