#!/usr/bin/env python3

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
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType
from functools import lru_cache, wraps

def roundObj(obj, digits : int):
    """ round <obj> to <digits> digits """
    if isinstance(obj,Iterable):
        new_obj = tuple([np.round(x,digits) for x in obj])
    else:
        new_obj = float(np.round(obj,digits))
    return new_obj

def roundCache(argname = None, argpos : int = 0, digits : int = 8,
        maxsize: int = 128, verbose : bool = False,
        turnoff : bool = False ):
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
            if verbose:
                print ( f"[cache] mu={args[1]} kwargs {kwargs} returns {func_cache(*rounded_args, **rounded_kwargs)}" )
                print ( f"[cache]              orig {function(*rounded_args, **rounded_kwargs)}" )
            if turnoff:
                return function(*rounded_args, **rounded_kwargs )
            return func_cache(*rounded_args, **rounded_kwargs)

        return wrapper

    return roundCacheDec

if __name__ == "__main__":
    from typing import Union, Text
    class TestCase:
        def __init__ ( self ):
            self.bla = 0
        @roundCache(argname='mu',argpos=1,digits=5)
        def myfunc ( self, mu : float = 0., evaluationType : NllEvalType = observed,
                         asimov : Union[None,float] = None ):
            print ( "calling myfunc" )
            ret = mu
            if evaluationType == apriori:
                ret += 10.
            if evaluationType == aposteriori:
                ret += 20.
            if asimov != None:
                ret += 1000.*(asimov+1000.)
            if asimov == None:
                ret = None
            return ret


    test = TestCase()
    print ( test.myfunc ( 1. ) )
    print ( test.myfunc ( 1. ) )
    print ( test.myfunc ( mu=1. ) )
    print ( test.myfunc ( mu=1. ) )
    print ( test.myfunc ( 1.,evaluationType=apriori, asimov=0. ) )
    print ( test.myfunc ( 1.,evaluationType=apriori, asimov=0. ) )
    print ( test.myfunc ( 1.,apriori, asimov=0. ) )
    print ( test.myfunc ( 1.,apriori, asimov=0. ) )
