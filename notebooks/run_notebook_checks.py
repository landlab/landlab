#! /usr/bin/env python
import pathlib
from dataclasses import dataclass
from functools import partial

import click


out = partial(click.secho, bold=True, err=True)
err = partial(click.secho, fg="red", err=True)


def collect_notebooks(src):
    p = pathlib.Path(src)
    if p.is_dir():
        return set([_p.absolute() for _p in iter_notebooks_in_dir(p, src)])
    elif is_a_notebook(p):
        return set([p.absolute()])
    else:
        raise ValueError("{0}: not a directory or a notebook".format(src))


def is_a_notebook(path):
    return path.is_file() and path.suffix == ".ipynb"


def iter_notebooks_in_dir(path, root):
    for s in path.iterdir():
        p = pathlib.Path(s)

        if p.is_dir() and p.name not in (".git", ".ipynb_checkpoints"):
            yield from iter_notebooks_in_dir(p, root)
        elif is_a_notebook(p):
            # elif p.is_file() and p.suffix == ".ipynb":
            yield p


def _notebook_cell_is_clean(cell):
    """Check if a single notebook cell is clean."""
    return cell["cell_type"] != "code" or (
        not cell["outputs"] and not cell["execution_count"]
    )


def _notebook_check_is_clean(path_to_notebook):
    """Check that a notebook is clean.

    A notebook is considered clean if all of its code cells have no output and
    have not execution count.

    Parameters
    ----------
    path_to_notebook : str or Path
        Path to the notebook to check.

    Returns
    -------
    bool
        ``True`` if the notebook is clean, otherwise ``False``.
    """
    import nbformat

    nb = nbformat.read(path_to_notebook, nbformat.current_nbformat)
    for cell in nb.cells:
        if not _notebook_cell_is_clean(cell):
            return False
    return True


@dataclass
class Report:
    clean_count: int = 0
    dirty_count: int = 0

    def __str__(self):
        report = []
        s = "s" if self.dirty_count != 1 else ""
        report.append(click.style(f"{self.dirty_count} messy notebook{s}", bold=True))
        report.append(click.style(f"{self.clean_count} clean"))

        return ", ".join(report)


@click.command()
@click.option(
    "-v", "--verbose", is_flag=True, help="Also emit status messages to stderr."
)
@click.option(
    "-q",
    "--quiet",
    is_flag=True,
    help=(
        "Don't emit non-error messages to stderr. Errors are still emitted, "
        "silence those with 2>/dev/null."
    ),
)
@click.argument(
    "path",
    type=click.Path(
        exists=True, file_okay=True, dir_okay=True, readable=True, path_type=None
    ),
)
@click.pass_context
def notebook_is_clean(ctx, verbose, quiet, path):
    error_msg = "Oh no! 💥 💔 💥"

    clean = []
    unclean = []
    for notebook in collect_notebooks(path):
        if _notebook_check_is_clean(notebook):
            clean.append(str(notebook))
        else:
            unclean.append(str(notebook))

    report = Report(clean_count=len(clean), dirty_count=len(unclean))

    if verbose:
        for notebook in sorted(clean + unclean):
            out(notebook) if notebook in clean else err(notebook)

    if not quiet:
        out(error_msg if report.dirty_count else "All done! ✨ 🍰 ✨")
        out(str(report))

    for notebook in unclean:
        print(notebook)

    ctx.exit(len(unclean))


if __name__ == "__main__":
    notebook_is_clean()
