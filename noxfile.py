import os
import pathlib
import shutil
import tempfile

import nox

PROJECT = "landlab"
ROOT = pathlib.Path(__file__).parent


@nox.session(venv_backend="mamba")
def test(session: nox.Session) -> None:
    """Run the tests."""
    session.conda_install("--file", "requirements.txt")
    session.conda_install("--file", "requirements-testing.txt")
    session.conda_install("richdem")
    session.install("-e", ".", "--no-deps")

    args = [
        "-n",
        "auto",
        "--cov",
        PROJECT,
        "-vvv",
        # "--dist", "worksteal",  # this is not available quite yet
    ] + session.posargs

    if "CI" in os.environ:
        args.append(f"--cov-report=xml:{ROOT.absolute()!s}/coverage.xml")
    session.run("pytest", *args)

    if "CI" not in os.environ:
        session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-notebooks", venv_backend="mamba")
def test_notebooks(session: nox.Session) -> None:
    """Run the notebooks."""
    args = [
        "pytest",
        "notebooks",
        "--nbmake",
        "--nbmake-kernel=python3",
        "--nbmake-timeout=3000",
        "-n",
        "auto",
        "-vvv",
    ] + session.posargs

    session.install("git+https://github.com/mcflugen/nbmake.git@mcflugen/add-markers")
    session.conda_install("richdem")
    session.conda_install(
        "pytest",
        "pytest-xdist",
        "--file",
        "requirements-notebooks.txt",
        "--file",
        "requirements.txt",
    )
    session.install("-e", ".", "--no-deps")

    session.run(*args)


@nox.session(name="test-cli")
def test_cli(session: nox.Session) -> None:
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


@nox.session(name="build-index")
def build_index(session: nox.Session) -> None:
    session.install(".[docs]")

    with open(ROOT / "docs" / "index.toml", "w") as fp:
        print("# This file was automatically generated with:", file=fp, flush=True)
        print("#     nox -s build-index", file=fp, flush=True)
        session.run(
            "landlab", "--silent", "index", "grids", "fields", "components", stdout=fp
        )


@nox.session(name="build-docs", venv_backend="mamba")
def build_docs(session: nox.Session) -> None:
    """Build the docs."""
    session.conda_install("richdem")
    session.install(".[docs]")

    clean_docs(session)
    session.run(
        "sphinx-build",
        "-b",
        "html",
        "-W",
        "--keep-going",
        str(ROOT / "docs/source"),
        str(ROOT / "build/html"),
    )


@nox.session(name="build-requirements")
def build_requirements(session: nox.Session) -> None:
    """Create requirements files from pyproject.toml."""
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
    session.run("python", "-m", "build", "--outdir", "./build/wheelhouse")


@nox.session
def release(session):
    """Tag, build and publish a new release to PyPI."""
    session.install("zest.releaser[recommended]")
    session.install("zestreleaser.towncrier")
    session.run("fullrelease")


@nox.session(name="publish-testpypi")
def publish_testpypi(session):
    """Publish wheelhouse/* to TestPyPI."""
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "--repository-url",
        "https://test.pypi.org/legacy/",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(name="publish-pypi")
def publish_pypi(session):
    """Publish wheelhouse/* to PyPI."""
    session.run("twine", "check", "build/wheelhouse/*")
    session.run(
        "twine",
        "upload",
        "--skip-existing",
        "build/wheelhouse/*.tar.gz",
    )


@nox.session(python=False)
def clean(session):
    """Remove all .venv's, build files and caches in the directory."""
    for folder in _args_to_folders(session.posargs):
        with session.chdir(folder):
            shutil.rmtree("build", ignore_errors=True)
            shutil.rmtree("build/wheelhouse", ignore_errors=True)
            shutil.rmtree(f"{PROJECT}.egg-info", ignore_errors=True)
            shutil.rmtree(".pytest_cache", ignore_errors=True)
            shutil.rmtree(".venv", ignore_errors=True)

            for pattern in ["*.py[co]", "__pycache__"]:
                _clean_rglob(pattern)


@nox.session(python=False, name="clean-checkpoints")
def clean_checkpoints(session):
    """Remove jupyter notebook checkpoint files."""
    for folder in _args_to_folders(session.posargs):
        with session.chdir(folder):
            _clean_rglob("*-checkpoint.ipynb")
            _clean_rglob(".ipynb_checkpoints")


@nox.session(python=False, name="clean-docs")
def clean_docs(session: nox.Session) -> None:
    """Clean up the docs folder."""
    session.chdir(ROOT / "build")
    if os.path.exists("html"):
        shutil.rmtree("html")


@nox.session(python=False, name="clean-ext")
def clean_ext(session: nox.Session) -> None:
    """Clean shared libraries for extension modules."""
    for folder in _args_to_folders(session.posargs):
        with session.chdir(folder):
            _clean_rglob("*.so")


@nox.session(python=False)
def nuke(session):
    """Run all clean sessions."""
    clean_checkpoints(session)
    clean_docs(session)
    clean(session)
    clean_ext(session)


@nox.session(name="list-wheels")
def list_wheels(session):
    print(os.linesep.join(_get_wheels(session)))


@nox.session(name="list-ci-matrix")
def list_ci_matrix(session):
    def _os_from_wheel(name):
        if "linux" in name:
            return "linux"
        elif "macos" in name:
            return "macos"
        elif "win" in name:
            return "windows"

    for wheel in _get_wheels(session):
        print(f"- cibw-only: {wheel}")
        print(f"  os: {_os_from_wheel(wheel)}")


def _get_wheels(session):
    platforms = session.posargs or ["linux", "macos", "windows"]
    session.install("cibuildwheel")

    wheels = []
    for platform in platforms:
        with tempfile.TemporaryFile("w+") as fp:
            session.run(
                "cibuildwheel",
                "--print-build-identifiers",
                "--platform",
                platform,
                stdout=fp,
            )
            fp.seek(0)
            wheels += fp.read().splitlines()
    return wheels


def _args_to_folders(args):
    return [ROOT] if not args else [pathlib.Path(f) for f in args]


def _clean_rglob(pattern):
    nox_dir = pathlib.Path(".nox")

    for p in pathlib.Path(".").rglob(pattern):
        if nox_dir in p.parents:
            continue
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()
