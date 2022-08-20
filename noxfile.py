import os
import pathlib
import shutil

import nox

PROJECT = "landlab"


@nox.session
def tests(session: nox.Session) -> None:
    """Run the tests."""
    session.install("pytest")
    session.install(".[dev,testing]")
    session.run("pytest", "--cov=landlab", "-vvv")
    session.run("coverage", "report", "--ignore-errors", "--show-missing")
    # "--fail-under=100",


@nox.session
def notebooks(session: nox.Session) -> None:
    """Run the notebooks."""
    session.install(".[dev,notebooks,testing]")
    session.run("pytest", "notebooks", "--run-notebook", "-n", "auto", "-vvv")


@nox.session
def cli(session: nox.Session) -> None:
    """Test the command line interface."""
    session.install(".")
    session.run("landlab", "--help")
    session.run("landlab", "--version")
    session.run("landlab", "index", "--help")
    session.run("landlab", "list", "--help")
    session.run("landlab", "provided-by", "--help")
    session.run("landlab", "provides", "--help")
    session.run("landlab", "used-by", "--help")
    session.run("landlab", "uses", "--help")
    session.run("landlab", "validate", "--help")


@nox.session
def lint(session: nox.Session) -> None:
    """Look for lint."""
    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files")

    # towncrier(session)


@nox.session
def towncrier(session: nox.Session) -> None:
    """Check that there is a news fragment."""
    session.install("towncrier")
    session.run("towncrier", "check", "--compare-with", "origin/master")


@nox.session
def docs(session: nox.Session) -> None:
    """Build the docs."""
    session.install(".[docs]")

    session.chdir("docs")
    if os.path.exists("build"):
        shutil.rmtree("build")
    session.run("sphinx-build", "-b", "html", "-W", "./source", "build/html")


@nox.session
def requirements(session: nox.Session) -> None:
    session.install("tomli")

    with open("requirements.txt", "w") as fp:
        session.run("python", "requirements.py", stdout=fp)

    for extra in ["dev", "docs", "notebooks", "testing"]:
        with open(f"requirements-{extra}.txt", "w") as fp:
            session.run("python", "requirements.py", extra, stdout=fp)


@nox.session
def build(session: nox.Session) -> None:
    """Build sdist and wheel dists."""
    session.install("pip")
    session.install("build")
    session.run("python", "--version")
    session.run("pip", "--version")
    session.run("python", "-m", "build", "--outdir", "./wheelhouse")


@nox.session
def release(session):
    """Tag, build and publish a new release to PyPI."""
    session.install("zest.releaser[recommended]")
    session.install("zestreleaser.towncrier")
    session.run("fullrelease")


@nox.session
def publish_testpypi(session):
    """Publish wheelhouse/* to TestPyPI."""
    session.run("twine", "check", "wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "wheelhouse/*.tar.gz",
    )


@nox.session
def publish_pypi(session):
    """Publish wheelhouse/* to PyPI."""
    session.run("twine", "check", "wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "wheelhouse/*.tar.gz",
    )


@nox.session(python=False)
def clean(session):
    """Remove all .venv's, build files and caches in the directory."""
    root = pathlib.Path(__file__).parent
    folders = (
        [root] if not session.posargs else [pathlib.Path(f) for f in session.posargs]
    )

    for folder in folders:
        with session.chdir(folder):
            _clean()


@nox.session(python=False)
def nuke(session):
    root = pathlib.Path(__file__).parent

    folders = (
        [root] if not session.posargs else [pathlib.Path(f) for f in session.posargs]
    )

    for folder in folders:
        with session.chdir(folder):
            _clean() and _nuke()


def _clean():
    patterns_to_clean = ["*.py[co]", "*.so", "__pycache__"]

    shutil.rmtree("build", ignore_errors=True)
    shutil.rmtree("wheelhouse", ignore_errors=True)
    shutil.rmtree(f"{PROJECT}.egg-info", ignore_errors=True)
    shutil.rmtree(".pytest_cache", ignore_errors=True)
    shutil.rmtree(".venv", ignore_errors=True)

    for pattern in patterns_to_clean:
        _clean_pattern(pattern)

    return True


def _nuke():
    patterns_to_nuke = [".ipynb_checkpoints"]

    for pattern in patterns_to_nuke:
        _clean_pattern(pattern)

    return True


def _clean_pattern(pattern):
    nox_dir = pathlib.Path(".nox")

    for p in pathlib.Path(".").rglob(pattern):
        if nox_dir in p.parents:
            continue
        if p.is_dir():
            # print(f"$ rm -r {str(p)}")
            p.rmdir()
        else:
            # print(f"$ rm {str(p)}")
            p.unlink()
