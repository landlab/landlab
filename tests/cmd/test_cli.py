import tomllib

import pytest

from landlab.cmd.landlab import main


def test_version():
    with pytest.raises(SystemExit):
        main(["--version"])


def test_help():
    with pytest.raises(SystemExit):
        main(["--help"])


@pytest.mark.parametrize(
    "subcommand",
    (
        "list",
        "used_by",
        "provided_by",
        "uses",
        "provides",
        "validate",
        "authors",
        "index",
        ("authors", "create"),
        ("authors", "build"),
        ("authors", "mailmap"),
        ("authors", "list"),
    ),
)
def test_subcommand_help(subcommand):
    if isinstance(subcommand, str):
        subcommand = (subcommand,)
    with pytest.raises(SystemExit):
        main(subcommand + ("--help",))


def test_list(capsys):
    main(["list"])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())


def test_used_by(capsys):
    main(["used_by"])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())


def test_provided_by(capsys):
    main(["provided_by"])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())


def test_provides(capsys):
    main(["provides", "topographic__elevation"])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())


def test_validate(capsys):
    main(["validate"])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())


@pytest.mark.parametrize("arg", ("components", "fields", "grids"))
def test_index(capsys, arg):
    main(["index", arg])
    index = tomllib.loads(capsys.readouterr().out)
    assert isinstance(index, dict)


@pytest.mark.parametrize("arg", ("mailmap",))
def test_authors(capsys, arg):
    main(["authors", arg])
    lines = capsys.readouterr().out.strip().splitlines()
    assert all(name for name in lines if name.strip())
