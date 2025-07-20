"""
.. module:: caching
   :synopsis: The memoize technique, for caching.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from functools import wraps
from unum import IncompatibleUnitsError
from smodels.base.physicsUnits import GeV, fb

import numpy as np
from collections.abc import Iterable
from functools import lru_cache,wraps

def roundObj(obj, digits : int):

    if isinstance(obj,Iterable):
        new_obj = tuple([np.round(x,digits) for x in obj])
    else:
        new_obj = np.round(obj,digits)
    return new_obj

def roundCache(argname = None, argpos : int = 0,
               digits : int = 5, maxsize: int = 128):
    """
    Returns the cached function called with the argument defined
    by argname and in the position argpos rounded to the
    amount of desired digits.
    """

    def roundCacheDec(function):
        func_cache = lru_cache(maxsize=maxsize)(function)

        @wraps(func_cache)
        def wrapper(*args, **kwargs):
            rounded_kwargs = kwargs
            rounded_args = list(args)
            if argname in rounded_kwargs:
                rounded_kwargs[argname] = roundObj(rounded_kwargs[argname],digits)
            elif argpos < len(rounded_args):
                rounded_args[argpos] = roundObj(rounded_args[argpos],digits)
            rounded_args = tuple(rounded_args)
            #print ( f"@@CA roundCache {rounded_args} {kwargs} {func_cache(*rounded_args, **rounded_kwargs)} {function(*rounded_args, **rounded_kwargs)}" )
            #print ( f"@@CA guess1: function(0.,None,True,False,None) {function(0.,None,True,False,None)}" )
            #print ( f"@@CA guess2: function(0.,None,True,False,0.) {function(0.,None,True,False,0.)}" )
            return func_cache(*rounded_args, **rounded_kwargs)

        return wrapper

    return roundCacheDec

def _toString(arg):
    try:
        return "%.2f" % arg.asNumber(fb)
    except (AttributeError, IncompatibleUnitsError):
        pass
    try:
        return "%.3f" % arg.asNumber(GeV)
    except (AttributeError, IncompatibleUnitsError):
        pass
    try:
        return "%.2f" % arg.asNumber(1/fb)
    except (AttributeError, IncompatibleUnitsError):
        pass
    if type(arg) == float:
        return "%.2f" % arg
    if type(arg) == int:
        return "%d" % arg
    if type(arg) == str:
        return "%s" % arg
    if type(arg) in [list, tuple]:
        argstring = ""
        for newarg in arg:
            argstring += _toString(newarg) + " "
        argstring = argstring[:-1]
        return argstring
    return "%s" % (str(arg))


class Cache:
    """ a class for storing results from interpolation """
    _cache = {}
    _cache_order = []
    n_stored = 1000  # number of interpolations we keep per result

    @staticmethod
    def size():
        return len(Cache._cache_order)

    @staticmethod
    def _clear_garbage():
        """
        every once in a while we clear the garbage, i.e.
        discard half of the cached interpolations
        """
        if len(Cache._cache_order) < Cache.n_stored:
            return
        ns2 = int(Cache.n_stored / 2)
        for i in range(ns2):
            Cache._cache.pop(Cache._cache_order[i])  # remove
        Cache._cache_order = Cache._cache_order[ns2:]

    @staticmethod
    def reset():
        """ completely reset the cache """
        Cache._cache = {}
        Cache._cache_order = []

    @staticmethod
    def add(key, value):
        if key in Cache._cache_order:
            return value
        Cache._cache_order.append(key)
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
        argstring = _toString(args)
        return Cache.add(argstring, func(*args))
    return _wrap
