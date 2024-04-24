"""Get the landlab tutorial and teaching notebooks.

Run this script to fetch the set of *landlab* notebooks compatible with
an installed version of landlab.

Usage
-----

Get notebooks for a currently installed *landlab*,

    $ python -m notebooks

To get notebooks for a particular version of *landlab*, provide
a version number as an argument. For example,

    $ python -m notebooks 2.5.0
"""

import argparse
import os
import pathlib
import sys
import tarfile
from urllib.error import HTTPError
from urllib.parse import urljoin
from urllib.request import urlopen

from packaging.version import Version


def main(version=None):
    if version in ("master", "latest", "dev"):
        version = None

    if not version or (version := Version(version)).is_devrelease:
        tag = "master"
    else:
        tag = "v" + version.base_version

    notebooks = NotebookFetcher(tag)

    out(f"fetching notebooks for landlab {notebooks.version}")
    out(f"{notebooks.url}")

    try:
        stream = notebooks.open()
    except NotebookError as error:
        err(str(error))
        return -1

    with tarfile.open(fileobj=stream, mode="r|gz") as tfile:
        base = NotebookExtractor(tfile).extract()

    out(f"notebooks have been extracted into {base}")
    out("To run the notebooks first install the required dependencies:")
    out("")
    out(f"    $ conda install --file={base}/requirements-notebooks.txt")
    out("")
    out("and then open the welcome notebook:")
    out("")
    out(f"    $ jupyter notebook {base}/notebooks/welcome.ipynb")
    print(base)

    return 0


class NotebookError(Exception):
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return self._msg


class NotebookFetcher:
    URL = "https://github.com/landlab/landlab/archive/refs"

    def __init__(self, version):
        self._version = version

    @property
    def version(self):
        return self._version

    @property
    def url(self):
        return urljoin(NotebookFetcher.URL, f"{self.version}.tar.gz")

    def open(self):
        try:
            stream = urlopen(self.url)
        except HTTPError as error:
            if error.code == 404:
                msg = f"unable to find notebooks for requested landlab version ({self.version})"
            else:
                msg = f"unable to fetch notebooks ({error.reason})"
            raise NotebookError(msg) from error
        else:
            return stream


class NotebookExtractor:
    def __init__(self, tfile):
        self._tfile = tfile
        self._names = []

    def extract(self):
        self._tfile.extractall(members=self._notebooks())
        return self.base

    def _notebooks(self):
        for tarinfo in self._tfile:
            parts = pathlib.Path(tarinfo.name).parts
            if (len(parts) > 1 and parts[1] == "notebooks") or (
                parts[-1] == "requirements-notebooks.txt"
            ):
                self._names.append(tarinfo.name)
                yield tarinfo

    @property
    def names(self):
        return sorted(self._names)

    @property
    def base(self):
        return pathlib.Path(os.path.commonprefix(self.names))


def out(text):
    print("\033[1m" + text + "\033[0m", file=sys.stderr)


def err(text):
    print("\033[91m" + text + "\033[0m", file=sys.stderr)


if __name__ == "__main__":
    try:
        import landlab
    except ModuleNotFoundError:
        version = ""
    else:
        version = landlab.__version__

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "version",
        metavar="VERSION",
        nargs="?",
        default=version,
        help="a landlab version (e.g. 2.5.0)",
    )
    args = parser.parse_args()

    sys.exit(main(args.version))
