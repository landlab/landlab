import itertools
import os
import subprocess
from collections import ChainMap

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


def load_config():
    config = load_first_of(["authors.toml", "pyproject.toml"], "rb")
    return config


def load_first_of(files, mode="rb"):
    files = [files] if isinstance(files, str) else files

    config = {}
    for name in files:
        try:
            with open(name, mode=mode) as fp:
                config = tomllib.load(fp)["tool"]["rollcall"]
        except FileNotFoundError:
            pass
        except KeyError:
            pass
        else:
            break
    try:
        config.pop("author")
    except KeyError:
        pass
    return config


class AuthorsConfig:
    _files = ["authors.toml", "pyproject.toml"]
    _defaults = {
        "authors_file": "AUTHORS.rst",
        "ignore": (),
        "author_format": "{name}",
    }

    def __init__(self, **kwds):

        self._config = ChainMap(
            {k: v for k, v in kwds.items() if k in self._defaults},
            load_first_of(self._files, "rb"),
            self._defaults,
        )

    def __getitem__(self, key):
        return self._config[key]

    def __str__(self):
        lines = ["[tool.rollcall]"]
        for key, value in sorted(self._config.items(), key=lambda item: item[0]):
            lines.append(f"{key} = {value!r}")
        return os.linesep.join(lines)


class GitLog:
    def __init__(self, format):
        self._format = format
        self._args = ["git", "log", f"--format={self._format}"]

    def __call__(self):
        process = subprocess.run(
            self._args,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            stdin=subprocess.PIPE,
        )
        return process

    def __str__(self):
        return " ".join(self._args)

    def __repr__(self):
        return f"GitLog({self._format!r})"


class Author:
    def __init__(self, name, email, aliases=None, alternate_emails=None):
        self._name = name
        self._email = email
        self._aliases = set() if not aliases else set(aliases)
        self._alternate_emails = (
            set() if not alternate_emails else set(alternate_emails)
        )
        self._aliases.discard(self.name)
        self._alternate_emails.discard(self.email)

    @classmethod
    def from_dict(cls, attrs):
        return cls(
            attrs["name"],
            attrs["email"],
            aliases=attrs.get("aliases", None),
            alternate_emails=attrs.get("alternate_emails", None),
        )

    @property
    def name(self):
        return self._name

    @property
    def names(self):
        return (self.name,) + self.aliases

    @property
    def email(self):
        return self._email

    @property
    def emails(self):
        return (self.email,) + self.alternate_emails

    @property
    def aliases(self):
        return tuple(self._aliases)

    @property
    def alternate_emails(self):
        return tuple(self._alternate_emails)

    def add_alias(self, alias):
        if alias != self.name:
            self._aliases.add(alias)

    def to_toml(self):
        lines = [
            "[[tool.rollcall.author]]",
            f"name = {self.name!r}",
            f"email = {self.email!r}",
        ]
        lines += (
            ["aliases = ["]
            + [f"  {alias!r}," for alias in sorted(self.aliases)]
            + ["]"]
        )
        lines += (
            ["alternate_emails = ["]
            + [f"  {email!r}," for email in sorted(self.alternate_emails)]
            + ["]"]
        )

        return os.linesep.join(lines)

    def update(self, other):
        if other.name != self.name:
            self._aliases.add(other.name)
        if other.email != self.email:
            self._alternate_emails.add(other.email)
        self._aliases |= set(other.aliases)
        self._alternate_emails |= set(other.alternate_emails)

    def __repr__(self):
        aliases = None if not self.aliases else self.aliases
        alternate_emails = None if not self.alternate_emails else self.alternate_emails
        return f"Author({self.name!r}, {self.email!r}, aliases={aliases!r}, alternate_emails={alternate_emails!r})"


class AuthorList:
    def __init__(self, authors=None):
        authors = [] if authors is None else authors
        self._name = {}
        self._email = {}

        for author in authors:
            for name in author.names:
                self._name[name] = author
            for email in author.emails:
                self._email[email] = author

    def __iter__(self):
        names = set([author.name for author in self._name.values()])
        for name in sorted(names):
            yield self._name[name]

    def __len__(self):
        names = set([author.name for author in self._name.values()])
        return len(names)

    def find_author(self, name_or_email):
        if name_or_email in self._name:
            author = self._name[name_or_email]
        elif name_or_email in self._email:
            author = self._email[name_or_email]
        else:
            raise KeyError(f"unknown author: {name_or_email}")
        return author

    @classmethod
    def from_toml(cls, toml_file):
        with open(toml_file, "rb") as fp:
            authors = [
                Author.from_dict(author)
                for author in tomllib.load(fp)["tool"]["rollcall"]["author"]
            ]

        return cls(authors=authors)

    @classmethod
    def from_csv(cls, csv):
        authors = cls()
        for line in csv.splitlines():
            name, email = (item.strip() for item in line.rsplit(",", maxsplit=1))
            authors.add(name, email)
        return authors

    def update(self, other):
        for author in other:
            for name, email in itertools.product(author.names, author.emails):
                self.add(name, email)

    def add(self, name, email):
        have_name = name in self._name
        have_email = email in self._email

        new_author = Author(name, email)

        if not (have_name or have_email):
            author = new_author
            self._name[name] = new_author
            self._email[email] = new_author
        elif have_name and not have_email:
            author = self._name[name]
            author.update(new_author)
            self._email[email] = author
        elif have_email and not have_name:
            author = self._email[email]
            author.update(new_author)
            self._name[name] = author
        elif self._name[name].name != self._email[email].name:
            other = self._email[email]
            author = self._name[name]
            author.update(other)
            for name in other.names:
                self._name[name] = author
            for email in other.emails:
                self._email[email] = author
        else:
            author = self._name[name]
