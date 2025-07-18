"""
.. module:: caching
   :synopsis: Extend the functools caching facility to deal with functions with float arguments
              which need to be rounded when caching.

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
            
                
            return func_cache(*rounded_args, **rounded_kwargs)

        return wrapper
    
    return roundCacheDec
