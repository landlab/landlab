import pathlib
import subprocess

import pytest
import yaml

from run_notebook_checks import _notebook_check_is_clean


_exclude_file = pathlib.Path(__file__).absolute().parent / "exclude.yml"
with open(_exclude_file, "r") as fp:
    _EXCLUDE = dict([(item["file"], item["reason"]) for item in yaml.safe_load(fp)])


def _notebook_check(notebook):
    """Check a notebook for errors.
    Parameters
    ----------
    notebook : Notebook node
        Path the to notebook to execute.
    Returns
    -------
    errors : list of str
        A list of the errors encountered in the notebook cells.
    """
    errors = [
        output
        for cell in notebook.cells
        if "outputs" in cell
        for output in cell["outputs"]
        if output.output_type == "error"
    ]

    return errors


def _notebook_run(path_to_notebook):
    """Execute a notebook via nbconvert and collect output.

    Parameters
    ----------
    path_to_notebook : str or Path
        Path to the notebook to execute.

    Returns
    -------
    nb : NotebookNode
        The parsed notebook object
    """
    import uuid

    import nbformat

    unique_name = pathlib.Path(
        "{prefix}{suffix}".format(prefix=str(uuid.uuid4()), suffix=".ipynb")
    )

    try:
        subprocess.check_call(
            [
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--ExecutePreprocessor.kernel_name=python",
                "--ExecutePreprocessor.timeout=-1",
                "--output",
                str(unique_name),
                "--output-dir=.",
                str(path_to_notebook),
            ]
        )
    except subprocess.CalledProcessError:
        raise
    else:
        nb = nbformat.read(unique_name, nbformat.current_nbformat)
    finally:
        try:
            unique_name.unlink()
        except FileNotFoundError:
            pass

    return nb


@pytest.mark.notebook
def test_notebook_is_clean(notebook):
    assert _notebook_check_is_clean(notebook)


@pytest.mark.notebook
def test_notebook(tmpdir, notebook):
    try:
        pytest.skip(_EXCLUDE[pathlib.Path(notebook).name])
    except KeyError:
        pass

    with tmpdir.as_cwd():
        nb = _notebook_run(notebook)
        assert _notebook_check(nb) == []
