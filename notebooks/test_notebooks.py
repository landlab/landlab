import os
import pathlib
import subprocess
import tempfile

import pytest
import yaml


_exclude_file = pathlib.Path(__file__).absolute().parent / "exclude.yml"
with open(_exclude_file, "r") as fp:
    _EXCLUDE = dict([(item["file"], item["reason"]) for item in yaml.safe_load(fp)])


def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
    import nbformat

    _, notebook = os.path.split(path)
    base, ext = os.path.splitext(notebook)

    with tempfile.NamedTemporaryFile("w", suffix=".ipynb") as fp:
        args = [
            "jupyter",
            "nbconvert",
            "--to",
            "notebook",
            "--execute",
            "--ExecutePreprocessor.kernel_name=python",
            "--ExecutePreprocessor.timeout=-1",
            "--output",
            fp.name,
            "--output-dir=.",
            path,
        ]
        subprocess.check_call(args)

        nb = nbformat.read(fp.name, nbformat.current_nbformat, encoding="UTF-8")

    errors = [
        output
        for cell in nb.cells
        if "outputs" in cell
        for output in cell["outputs"]
        if output.output_type == "error"
    ]

    return nb, errors


@pytest.mark.notebook
def test_notebook(tmpdir, notebook):
    try:
        pytest.skip(_EXCLUDE[pathlib.Path(notebook).name])
    except KeyError:
        pass

    with tmpdir.as_cwd():
        nb, errors = _notebook_run(notebook)
        assert not errors
