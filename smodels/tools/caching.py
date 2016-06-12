"""
.. module:: tools.caching
   :synopsis: The memoize technique, for caching.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from functools import wraps                                                           
from smodels.tools.physicsUnits import pb, GeV, fb

def _toString ( arg ):
    try:
        return "%.2f " % arg.asUnit(fb)
    except Exception,e:
        pass
    try:
        return "%.3f " % arg.asNumber(GeV)
    except Exception,e:
        pass
    try:
        return "%.2f " % arg.asUnit(1/fb)
    except Exception,e:
        pass
    if type(arg) == float:
        return "%.2f " % arg
    if type(arg) == int:
        return "%d " % arg
    if type(arg) == str:
        return "%s " % arg
    if type(arg) in [ list, tuple ]:
        argstring=""
        for newarg in arg:
            argstring += _toString ( newarg )
        return argstring
    return "%s " % ( str(arg) )
    
    
_cache = {}

def clearCache():
    """
    Clear the cache
    """
    _cache={}
            
def _memoize(func):
    """
    Serves as a wrapper to cache the results of func, since this is a
    computationally expensive function.
    
    """
    @wraps(func)
    def _wrap(*args):
        """
        Wrapper for the function to be memoized
        """ 
        argstring = _toString ( args )
        if argstring not in _cache:
            _cache[argstring] = func(*args)
        return _cache[argstring]
    return _wrap
