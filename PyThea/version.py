# Try to import _version or else go to default
# The modification is made so during the cli call of PyThea the version
# the version imports correctly. I will explore this in the future.
# This is a modified script from sunpy
try:
    try:
        from PyThea._version import version
    except ImportError:
        from ._version import version
except Exception:
    import warnings
    warnings.warn(
        f'could not determine {__name__.split(".")[0]} package version; '
        f'this indicates a broken installation')
    del warnings

    version = '0.0.0'


# We use Version to define major, minor, micro, but ignore any suffixes.
def split_version(version):
    pieces = [0, 0, 0]

    try:
        from packaging.version import Version

        v = Version(version)
        pieces = [v.major, v.minor, v.micro]

    except Exception:
        pass

    return pieces


major, minor, bugfix = split_version(version)

del split_version  # clean up namespace.

release = 'dev' not in version
