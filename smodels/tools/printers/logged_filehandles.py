"""
.. module:: logged_filehandles
   :synopsis: monkey patching open to track which files
   haven been opened

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import builtins
import inspect

opened_by = {}

_original_open = builtins.open

def logged_open(file, *args, **kwargs):
    """ an open method that logs the filename """
    frame = inspect.currentframe().f_back

    if frame is not None:
        caller = (
            f"{frame.f_globals.get('__name__', '?')}."
            f"{frame.f_code.co_name}:"
            f"{frame.f_lineno}"
        )
    else:
        caller = "<unknown>"
    print ( f"[logged_open] file:{file} args:{args} kwargs:{kwargs} caller {caller}" )

    opened_by[str(file)] = caller

    return _original_open(file, *args, **kwargs)

def redirect_open():
    """ redirect the 'open' method to logged_open """
    builtins.open = logged_open

def restore_original ():
    """ restore the original 'open' method """
    builtins.open = _original_open
