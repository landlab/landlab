#! /usr/bin/env python
"""Fetch landlab tutorial and teaching notebooks for a given version of landlab.
"""
import argparse
import json
import os
import pathlib
import sys
import tarfile
import urllib
from packaging.version import Version
from urllib.error import HTTPError
from urllib.parse import urljoin
from urllib.request import urlopen


def main(version=None):
    notebooks = NotebookFetcher(version)

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
    print(base)

    return 0


class NotebookError(Exception):
    def __init__(self, msg):
        self._msg = msg

    def __str__(self):
        return self._msg


class NotebookFetcher:
    """Fetch the landlab notebooks.

    Parameters
    ----------
    version : str, optional
        A landlab version to get the notebooks for. If not provided,
        fetch notebooks for the currently installed version of landlab.
        If landlab's not installed, get the latest development notebooks.
    """

    URL = "https://github.com/landlab/landlab/archive/refs"
    API_URL = "https://api.github.com/repos/landlab/landlab/tags"

    def __init__(self, version=None):
        if not version:
            version = NotebookFetcher.landlab_version()

        if not version:
            self._version = "master"
        else:
            version = NotebookFetcher.get_closest_version(version)
            self._version = "v" + Version(version).base_version

    @staticmethod
    def landlab_version():
        """The currently installed version of landlab.

        Returns
        -------
        version : str
            The version of landlab that is installed or ``None`` if not installed.
        """
        try:
            from landlab import __version__
        except ModuleNotFoundError:
            return None
        else:
            return __version__

    @property
    def version(self):
        return self._version

    @property
    def url(self):
        return urljoin(NotebookFetcher.URL, f"{self.version}.tar.gz")

    def open(self):
        if Version(self._version) < Version("2"):
            raise NotebookError("notebooks are not available for landlab < 2")

        try:
            stream = urlopen(self.url)
        except HTTPError as error:
            if error.code == 404:
                msg = f"unable to find notebooks for requested landlab version ({self.version})"
            else:
                msg = f"unable to fetch notebooks ({error.reason})"
            raise NotebookError(msg)
        else:
            return stream

    @staticmethod
    def get_available_versions():
        response = urllib.request.urlopen(NotebookFetcher.API_URL)
        if response.status == 200:
            tags = json.loads(response.read())
        else:
            tags = []

        versions = []
        for tag in tags:
            if tag["name"].startswith("v"):
                versions.append(Version(tag["name"]))

        return sorted(versions)

    @staticmethod
    def get_closest_version(ver):
        import bisect

        versions = sorted(NotebookFetcher.get_available_versions())

        if Version(ver) <= versions[0]:
            closest = versions[0]
        else:
            closest = versions[bisect.bisect_right(versions, Version(ver)) - 1]
        return str(closest)


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
            if len(parts) > 1 and parts[1] == "notebooks":
                self._names.append(tarinfo.name)
                yield tarinfo
            elif parts[-1] == "requirements-notebooks.txt":
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
    parser = argparse.ArgumentParser(description="fetch landlab notebooks")
    parser.add_argument("version", nargs="?", default="")
    args = parser.parse_args()

    sys.exit(main(args.version))
