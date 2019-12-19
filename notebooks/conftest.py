import os
import pathlib

import pytest


_TEST_DIR = pathlib.Path(__file__).absolute().parent


def collect_notebooks(src):
    p = pathlib.Path(src)
    if p.is_dir():
        return set([_p.absolute() for _p in iter_notebooks_in_dir(p, src)])
    else:
        raise ValueError("{0}: not a directory".format(src))


def iter_notebooks_in_dir(path, root):
    for s in path.iterdir():
        p = pathlib.Path(s)
        normalized_path = "/" + p.resolve().relative_to(root).as_posix()

        if p.is_dir():
            normalized_path += "/"
        if p.is_dir() and p.name not in (".git", ".ipynb_checkpoints"):
            yield from iter_notebooks_in_dir(p, root)
        elif p.is_file() and p.suffix == ".ipynb":
            yield p


def pytest_generate_tests(metafunc):
    if "notebook" in metafunc.fixturenames:
        metafunc.parametrize(
            "notebook", sorted([str(name) for name in collect_notebooks(_TEST_DIR)])
        )


def pytest_configure(config):
    config.addinivalue_line("markers", "notebook: mark test as a notebook to run")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--run-notebook"):
        skip_notebook = pytest.mark.skip(reason="need --run-notebook option to run")
        for item in items:
            if "notebook" in item.keywords:
                item.add_marker(skip_notebook)
