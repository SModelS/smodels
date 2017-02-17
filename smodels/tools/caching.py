"""
.. module:: caching
   :synopsis: The memoize technique, for caching.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from functools import wraps                                                           
from smodels.tools.physicsUnits import pb, GeV, fb, IncompatibleUnitsError

def _toString ( arg ):
    try:
        return "%.2f" % arg.asNumber(fb)
    except (AttributeError,IncompatibleUnitsError) as e:
        pass
    try:
        return "%.3f" % arg.asNumber(GeV)
    except (AttributeError,IncompatibleUnitsError) as e:
        pass
    try:
        return "%.2f" % arg.asNumber(1/fb)
    except (AttributeError,IncompatibleUnitsError) as e:
        pass
    if type(arg) == float:
        return "%.2f" % arg
    if type(arg) == int:
        return "%d" % arg
    if type(arg) == str:
        return "%s" % arg
    if type(arg) in [ list, tuple ]:
        argstring=""
        for newarg in arg:
            argstring += _toString ( newarg ) + " "
        argstring=argstring[:-1]
        return argstring
    return "%s" % ( str(arg) )
    
class Cache:
    """ a class for storing results from interpolation """
    _cache = {}
    _cache_order=[]
    n_stored = 1000 ## number of interpolations we keep per result

    @staticmethod
    def size():
        return len(Cache._cache_order)

    @staticmethod
    def _clear_garbage ():
        """
        every once in a while we clear the garbage, i.e. 
        discard half of the cached interpolations
        """
        if len(Cache._cache_order)<Cache.n_stored:
            return
        ns2 = int ( Cache.n_stored / 2 )
        for i in range ( ns2 ):
            Cache._cache.pop ( Cache._cache_order[i] ) ## remove 
        Cache._cache_order=Cache._cache_order[ ns2 : ]

    @staticmethod
    def reset ():
        """ completely reset the cache """
        Cache._cache = {}
        Cache._cache_order = []

    @staticmethod
    def add ( key, value ):
        if key in Cache._cache_order:
            return value
        Cache._cache_order.append ( key )
        Cache._cache[key] = value
        Cache._clear_garbage()
        return value

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
        return Cache.add ( argstring, func(*args) )
    return _wrap
