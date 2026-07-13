from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("phosx")
except PackageNotFoundError:
    # package is not installed (e.g. running from a source tree without an install)
    __version__ = "0.0.0"
