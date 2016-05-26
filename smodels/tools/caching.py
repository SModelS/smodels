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
    #    print "argstring=",argstring
        return argstring
    return "%s " % ( str(arg) )
            
def _memoize(func):
    """
    Serves as a wrapper to cache the results of func, since this is a
    computationally expensive function.
    
    """
    cache = {}
    @wraps(func)
    def _wrap(*args):
        """
        Wrapper for the function to be memoized
        """ 
        argstring = _toString ( args )
        if argstring not in cache:
            cache[argstring] = func(*args)
        # print "[_wrap] len of cache is",len(cache)
        # print "[_wrap] args=",args
        ## print "[_wrap] argstring=",argstring
        # import unum
        # unum.Unum.VALUE_FORMAT = "%0.2E"
        # print "[_wrap] ret=",cache[argstring] ## .asNumber(pb)
        return cache[argstring]
    return _wrap
