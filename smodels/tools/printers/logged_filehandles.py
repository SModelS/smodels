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
    self = frame.f_locals.get("self")
    function = frame.f_code.co_name
    if self is not None:
        caller = f"{type(self).__name__}.{function}"
    else:
        caller = function

    sfile = str(file)
    # print ( f"[logged_filehandles] {sfile}: {caller}" )
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
