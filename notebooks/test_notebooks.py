import os
import pathlib
import subprocess
import tempfile

import nbformat
import pytest


_EXCLUDE = [
    "animate-landlab-output.ipynb",
    "cellular_automaton_vegetation_flat_domain.ipynb",
    "cellular_automaton_vegetation_DEM.ipynb",
    "stream_power_channels_class_notebook.ipynb",
]


def _notebook_run(path):
    """Execute a notebook via nbconvert and collect output.
       :returns (parsed nb object, execution errors)
    """
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
            "--ExecutePreprocessor.timeout=None",
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
    if pathlib.Path(notebook).name in _EXCLUDE:
        pytest.skip("notebook marked as failing")

    with tmpdir.as_cwd():
        nb, errors = _notebook_run(notebook)
        assert not errors
