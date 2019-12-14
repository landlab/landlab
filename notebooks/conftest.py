import os

import pytest


_TEST_DIR = os.path.abspath(os.path.dirname(__file__))
_EXCLUDE = [
    "animate-landlab-output.ipynb",
    "cellular_automaton_vegetation_flat_domain.ipynb",
    "cellular_automaton_vegetation_DEM.ipynb",
]


def _all_notebooks(path="."):
    notebooks = []
    for root, dirs, files in os.walk(path):
        if ".ipynb_checkpoints" in root:
            continue
        for file in files:
            if file.endswith(".ipynb") and (file not in _EXCLUDE):
                notebooks.append(os.path.join(root, file))
    return notebooks


def pytest_generate_tests(metafunc):
    if "notebook" in metafunc.fixturenames:
        metafunc.parametrize("notebook", _all_notebooks(_TEST_DIR))


def pytest_addoption(parser):
    parser.addoption(
        "--run-notebook", action="store_true", default=False, help="run notebook tests"
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "notebook: mark test as a notebook to run")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-notebook"):
        skip_notebook = pytest.mark.skip(reason="need --run-notebook option to run")
        for item in items:
            if "notebook" in item.keywords:
                item.add_marker(skip_notebook)
