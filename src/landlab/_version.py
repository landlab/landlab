from importlib.metadata import PackageNotFoundError
from importlib.metadata import version

try:
    __version__ = version("landlab")
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
