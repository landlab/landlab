import json
import os
import pathlib
import shutil
import sys

import nox
from packaging.requirements import Requirement

PROJECT = "landlab"
ROOT = pathlib.Path(__file__).parent
PYTHON_VERSION = "3.12"
PATH = {
    "build": ROOT / "build",
    "docs": ROOT / "docs",
    "nox": pathlib.Path(".nox"),
    "requirements": ROOT / "requirements",
    "root": ROOT,
}


@nox.session(python=PYTHON_VERSION, venv_backend="conda")
def test(session: nox.Session) -> None:
    """Run the tests."""
    os.environ["WITH_OPENMP"] = "1"

    session.log(f"CC = {os.environ.get('CC', 'NOT FOUND')}")

    if sys.platform.startswith("darwin") and session.python == "3.12":
        session.log("installing multidict from conda-forge.")
        session.conda_install("multidict")

    session.install(
        "-r",
        PATH["requirements"] / "required.txt",
        "-r",
        PATH["requirements"] / "testing.txt",
    )

    session.conda_install("richdem", channel=["nodefaults", "conda-forge"])
    session.install("-e", ".", "--no-deps")

    check_package_versions(session, files=["required.txt", "testing.txt"])

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


@nox.session(name="test-notebooks", python=PYTHON_VERSION, venv_backend="conda")
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

    os.environ["WITH_OPENMP"] = "1"

    if sys.platform.startswith("darwin") and session.python == "3.12":
        session.log("installing multidict from conda-forge")
        session.conda_install("multidict")

    session.install(
        "-r",
        PATH["requirements"] / "required.txt",
        "-r",
        PATH["requirements"] / "testing.txt",
        "-r",
        PATH["requirements"] / "notebooks.txt",
    )
    session.conda_install("richdem", channel=["nodefaults", "conda-forge"])
    session.install("git+https://github.com/mcflugen/nbmake.git@mcflugen/add-markers")
    session.install("-e", ".", "--no-deps")

    check_package_versions(
        session, files=["required.txt", "testing.txt", "notebooks.txt"]
    )

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
    skip_hooks = [] if "--no-skip" in session.posargs else ["check-manifest", "pyroma"]

    session.install("pre-commit")
    session.run("pre-commit", "run", "--all-files", env={"SKIP": ",".join(skip_hooks)})


@nox.session
def towncrier(session: nox.Session) -> None:
    """Check that there is a news fragment."""
    session.install("towncrier")
    session.run("towncrier", "check", "--compare-with", "origin/master")


@nox.session(name="build-index")
def build_index(session: nox.Session) -> None:
    index_file = ROOT / "docs" / "index.toml"
    header = """
# This file was automatically generated with:
#     nox -s build-index
    """.strip()

    session.install("sphinx")
    session.install(".")

    with open(index_file, "w") as fp:
        print(header, file=fp, flush=True)
        session.run(
            "landlab", "--silent", "index", "components", "fields", "grids", stdout=fp
        )
    session.log(f"generated index at {index_file!s}")


@nox.session(name="build-docs")
def build_docs(session: nox.Session) -> None:
    """Build the docs."""

    session.install(
        "-r",
        PATH["requirements"] / "docs.txt",
        "-r",
        PATH["requirements"] / "required.txt",
    )
    session.install("-e", ".", "--no-deps")

    check_package_versions(session, files=["required.txt", "docs.txt"])

    PATH["build"].mkdir(exist_ok=True)
    session.run(
        "sphinx-build",
        "-b",
        "html",
        "-W",
        "--keep-going",
        PATH["docs"] / "source",
        PATH["build"] / "html",
    )
    session.log(f"generated docs at {PATH['build'] / 'html'!s}")


@nox.session(name="check-versions")
def check_package_versions(session, files=("required.txt",)):
    output_lines = session.run("pip", "list", "--format=json", silent=True).splitlines()

    installed_version = {
        p["name"].lower(): p["version"] for p in json.loads(output_lines[0])
    }

    for file_ in files:
        required_version = {}
        with (PATH["requirements"] / file_).open() as fp:
            for line in fp.readlines():
                requirement = Requirement(line)
                required_version[requirement.name.lower()] = requirement.specifier

        mismatch = set()
        for name, version in required_version.items():
            if name not in installed_version or not version.contains(
                installed_version[name]
            ):
                mismatch.add(name)

        session.log(f"Checking installed package versions for {file_}")
        for name in sorted(required_version):
            print(f"[{name}]")
            print(f"requested = {str(required_version[name])!r}")
            if name in installed_version:
                print(f"installed = {installed_version[name]!r}")
            else:
                print("installed = false")

        if mismatch:
            session.warn(
                f"There were package version mismatches for packages required in {file_}"
            )


@nox.session
def locks(session: nox.Session) -> None:
    """Create lock files."""
    folders = session.posargs or [".", "docs", "notebooks"]

    session.install("pip-tools")

    def upgrade_requirements(src, dst="requirements.txt"):
        with open(dst, "wb") as fp:
            session.run("pip-compile", "--upgrade", src, stdout=fp)

    for folder in folders:
        with session.chdir(ROOT / folder):
            upgrade_requirements("requirements.in", dst="requirements.txt")

    for folder in folders:
        session.log(f"updated {ROOT / folder / 'requirements.txt'!s}")

    # session.install("conda-lock[pip_support]")
    # session.run("conda-lock", "lock", "--mamba", "--kind=lock")


@nox.session(name="sync-requirements", python=PYTHON_VERSION, venv_backend="conda")
def sync_requirements(session: nox.Session) -> None:
    """Sync requirements.in with pyproject.toml."""
    with open("requirements.in", "w") as fp:
        session.run(
            "python",
            "-c",
            """
import os, tomllib
with open("pyproject.toml", "rb") as fp:
    print(os.linesep.join(sorted(tomllib.load(fp)["project"]["dependencies"])))
""",
            stdout=fp,
        )


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
    if (PATH["build"] / "html").is_dir():
        with session.chdir(PATH["build"]):
            shutil.rmtree("html")

    if PATH["build"].is_dir():
        session.chdir(PATH["build"])
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
        wheels += session.run(
            "cibuildwheel",
            "--print-build-identifiers",
            "--platform",
            platform,
            silent=True,
        ).splitlines()
    return wheels


def _args_to_folders(args):
    return [ROOT] if not args else [pathlib.Path(f) for f in args]


def _clean_rglob(pattern):
    for p in pathlib.Path(".").rglob(pattern):
        if PATH["nox"] in p.parents:
            continue
        if p.is_dir():
            p.rmdir()
        else:
            p.unlink()


@nox.session
def credits(session):
    """Update the various authors files."""
    from landlab.cmd.authors import AuthorsConfig

    config = AuthorsConfig()

    with open(".mailmap", "wb") as fp:
        session.run(
            "landlab", "--silent", "authors", "mailmap", stdout=fp, external=True
        )

    contents = session.run(
        "landlab",
        "--silent",
        "authors",
        "create",
        "--update-existing",
        external=True,
        silent=True,
    )
    with open(config["credits_file"], "w") as fp:
        print(contents, file=fp, end="")

    contents = session.run(
        "landlab", "--silent", "authors", "build", silent=True, external=True
    )
    with open(config["authors_file"], "w") as fp:
        print(contents, file=fp, end="")
