import glob
import json
import os
import pathlib
import shutil
import tempfile
from collections import defaultdict

import nox

PROJECT = "landlab"
ROOT = pathlib.Path(__file__).parent


@nox.session(python="3.11", venv_backend="mamba")
def test(session: nox.Session) -> None:
    """Run the tests."""
    os.environ["WITH_OPENMP"] = "1"

    # session.conda_install("c-compiler", "cxx-compiler")
    session.log(f"CC = {os.environ.get('CC', 'NOT FOUND')}")
    session.conda_install("--file", "requirements.in")
    session.conda_install("--file", "requirements-testing.in")
    session.conda_install("richdem")
    session.install("-e", ".", "--no-deps")

    args = [
        "-n",
        "auto",
        "--cov",
        f"src/{PROJECT}",
        "-vvv",
        # "--dist", "worksteal",  # this is not available quite yet
    ] + session.posargs

    if "CI" in os.environ:
        args.append(f"--cov-report=xml:{ROOT.absolute()!s}/coverage.xml")
    session.run("pytest", *args)

    if "CI" not in os.environ:
        session.run("coverage", "report", "--ignore-errors", "--show-missing")


@nox.session(name="test-notebooks", python="3.11", venv_backend="mamba")
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

    # session.conda_install("c-compiler", "cxx-compiler")
    session.conda_install("richdem", channel=["nodefaults", "conda-forge"])
    session.conda_install(
        "pytest",
        "pytest-xdist",
        "--file",
        "notebooks/requirements.in",
        "--file",
        "requirements-testing.in",
        "--file",
        "requirements.in",
        channel=["nodefaults", "conda-forge"],
    )
    session.install("git+https://github.com/mcflugen/nbmake.git@mcflugen/add-markers")
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


@nox.session(name="build-gallery-index")
def build_notebook_index(session: nox.Session) -> None:
    docs_dir = ROOT / "docs" / "source"

    for gallery in ["tutorials", "teaching"]:
        index = collect_notebooks(docs_dir / gallery)

        for subdir, entries in sorted(index.items()):
            title = pathlib.Path(subdir).stem.replace("_", " ").title()

            path_to_notebooks = pathlib.Path(gallery) / subdir

            lines = [
                f"{title}",
                f"{'-'*len(title)}",
                "",
                ".. nbgallery::",
                "    :glob:",
                "",
            ] + [f"    /{path_to_notebooks / v!s}" for v in entries]

            generated_dir = docs_dir / "generated" / path_to_notebooks
            generated_file = generated_dir / "_index.rst"

            generated_dir.mkdir(parents=True, exist_ok=True)
            with open(generated_file, "w") as fp:
                print(os.linesep.join(lines), file=fp)
            session.log(generated_file)


def collect_notebooks(path_to_notebooks):
    paths = pathlib.Path(path_to_notebooks)

    index = defaultdict(list)

    for p in (p for p in paths.iterdir() if p.is_dir()):
        subdir = p.relative_to(path_to_notebooks)

        if glob.glob(str(p / "*.ipynb")) + glob.glob(str(p / "*.md")):
            index[subdir] += ["*"]
        if glob.glob(str(p / "**/*.ipynb")) + glob.glob(str(p / "**/*.md")):
            index[subdir] += ["**"]

    return index


@nox.session(name="summarize-linkcheck")
def summarize_linkcheck_output(session: nox.Session) -> None:
    path_to_file = ROOT / "build" / "linkcheck" / "output.json"

    broken = _summarize_broken_links(path_to_file)

    if broken:
        print(broken)
    else:
        session.log("Nothing to summarize, there were no broken links")


def _summarize_broken_links(path_to_file) -> str:
    with open(path_to_file) as fp:
        entries = [json.loads(line) for line in fp.readlines()]

    broken = defaultdict(list)
    for entry in entries:
        if entry["status"] == "broken":
            broken[entry["filename"]] += [entry["uri"]]

    tables = []
    for filename in sorted(broken):
        tables.append(
            os.linesep.join(
                [
                    f"[{filename!r}]",
                    "broken = [",
                ]
                + [f"  {uri!r}," for uri in sorted(broken[filename])]
                + ["]"]
            )
        )
    return (2 * os.linesep).join(tables)


@nox.session(name="check-links", python="3.11", venv_backend="mamba")
def check_links(session: nox.Session) -> None:
    """Check for broken links in the docs."""
    builder = "linkcheck"

    build_dir = _build_docs(session, builders=[builder], success_codes=(0, 1))

    log_file = build_dir / builder / "output.json"

    broken = summarize_linkcheck_output(log_file)

    session.log(f"{log_file!s}")
    if broken:
        print(broken)
        session.error("Broken links were found")


@nox.session(name="build-docs", python="3.11", venv_backend="mamba")
def build_docs(session: nox.Session) -> None:
    """Build the docs."""
    builder = "html"

    build_dir = _build_docs(session, builders=[builder])

    session.log(f"generated docs at {build_dir / {builder} !s}")


def _build_docs(session, builders=("html",), success_codes=(0,)):
    build_dir = ROOT / "build"
    docs_dir = ROOT / "docs"

    build_notebook_index(session)

    session.conda_install(
        "pandoc",
        "--file",
        "requirements.in",
        channel=["nodefaults", "conda-forge"],
    )
    # session.conda_install(f"--file={docs_dir / 'requirements.in'!s}")
    # session.conda_install("--file", "requirements.in")
    session.install("-r", docs_dir / "requirements.in")
    session.install("-e", ".")

    build_dir.mkdir(exist_ok=True)
    for builder in builders:
        session.run(
            "sphinx-build",
            "-b",
            builder,
            "-W",
            "--keep-going",
            "--jobs",
            "auto",
            docs_dir / "source",
            build_dir / builder,
            success_codes=success_codes,
        )
        session.log(f"generated docs at {build_dir / builder !s}")

    return build_dir


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


@nox.session(name="sync-requirements", python="3.11", venv_backend="conda")
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
            shutil.rmtree(f"src/{PROJECT}.egg-info", ignore_errors=True)
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
    build_dir = ROOT / "build"

    if (build_dir / "html").is_dir():
        with session.chdir(build_dir):
            shutil.rmtree("html")

    if (ROOT / "build").is_dir():
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
