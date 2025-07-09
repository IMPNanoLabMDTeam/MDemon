import inspect

from .. import _READERS
from ..lib import util


def get_reader_for(filename, format=None):
    """ """
    if inspect.isclass(format):
        return format

    # Only guess if format is not provided
    if format is None:
        format = util.guess_format(filename)
    format = format.upper()
    return _READERS[format]
