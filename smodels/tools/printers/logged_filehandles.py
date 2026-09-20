"""
.. module:: logged_filehandles
   :synopsis: monkey patching open to track which files
   haven been opened

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import builtins
import inspect
from smodels.base.exceptions import SModelSBaseError as SModelSError

opened_by = {}

_original_open = builtins.open

def logged_open(file, *args, **kwargs):
    """ an open method that logs the filename """
    if "r" in args[0]: ## open for reading, we dont care
        return _original_open(file, *args, **kwargs)
    frame = inspect.currentframe().f_back

    if frame is not None:
        caller = f"{frame.f_globals.get('__name__', '?')}"
    else:
        caller = "<unknown>"
    sfile = str(file)
    if sfile in opened_by: 
        if caller != opened_by[sfile]:
            line = f"{sfile} is opened by {opened_by[sfile]} as well as {caller}"
            raise SModelSError ( line )
    opened_by[sfile] = caller
    return _original_open(file, *args, **kwargs)

def redirect_open():
    """ redirect the 'open' method to logged_open """
    builtins.open = logged_open

def restore_original ():
    """ restore the original 'open' method """
    builtins.open = _original_open
