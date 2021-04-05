"""
.. module:: unitTestHelpers
   :synopsis: helper functions for the unit tests

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os, sys
import unum
import numpy as np
import redirector
from smodels.tools.runSModelS import run
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.tools.smodelsLogging import logger, setLogLevel, getLogLevel

def equalObjs(obj1,obj2,allowedDiff,ignore=[], where=None, fname=None,
               fname2=None ):
    """
    Compare two objects.
    The numerical values are compared up to the precision defined by allowedDiff.

    :param obj1: First python object to be compared
    :param obj2: Second python object to be compared
    :param allowedDiff: Allowed % difference between two numerical values
    :param ignore: List of keys to be ignored
    :param where: keep track of where we are, for easier debugging.
    :param fname: the filename of obj1
    :param fname2: the filename of obj2
    :return: True/False
    """
    if type(fname)==str:
        fname = fname.replace( os.getcwd(), "." )
    if type(obj1) in [ float, int ] and type ( obj2) in [ float, int ]:
        obj1,obj2=float(obj1),float(obj2)

    if type(obj1) != type(obj2):
        logger.warning("Data types differ: (%s,%s) <-> (%s,%s) in %s:%s" %(obj1,type(obj1),obj2,type(obj2),where,fname))
        return False

    if isinstance(obj1,unum.Unum):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff.asNumber() < allowedDiff
    elif isinstance(obj1,float):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        if diff > allowedDiff:
            logger.error ( "values %s and %s differ by %s in %s:%s" % ( obj1, obj2, diff, where, fname) )
        return diff < allowedDiff
    elif isinstance(obj1,str):
        if obj1 != obj2:
            logger.error ( "strings ``%s'' and ``%s'' differ in %s:%s" % ( obj1, obj2, where, fname ) )
        return obj1 == obj2
    elif isinstance(obj1,dict):
        for key in obj1:
            if key in ignore: continue
            if not key in obj2:
                if where == None:
                    where = "unspecified"
                if fname2 == None:
                    fname2 = "unspecified"
                logger.warning("Key ``%s'' missing in %s:%s" % (key, where, fname2 ) )
                return False
            if not equalObjs(obj1[key],obj2[key],allowedDiff, ignore=ignore, where=key, fname = fname, fname2 = fname2 ):
                return False
    elif isinstance(obj1,list):
        if len(obj1) != len(obj2):
            logger.warning('Lists differ in length:\n   %i (this run)\n and\n   %i (default)' %\
                                (len(obj1),len(obj2)))
            return False
        for ival,val in enumerate(obj1):
            if not equalObjs(val,obj2[ival],allowedDiff, fname = fname,
                             fname2 = fname2 ):
                #logger.warning('Lists differ:\n   %s (this run)\n and\n   %s (default)' %\
                #                (str(val),str(obj2[ival])))
                return False
    else:
        return obj1 == obj2

    return True

def importModule ( filename ):
    """ import a module, but giving the filename """
    if sys.version_info[0]==2:
        import imp
        ## python2, use imp
        with open( filename, 'rb') as fp: ## imports file with dots in name
            output_module = imp.load_module("output",fp, filename, \
                    ('.py', 'rb', imp.PY_SOURCE) )
        return output_module.smodelsOutput
    ### python3, use importlib
    import importlib
    spec = importlib.util.spec_from_file_location( "output", filename )
    output_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(output_module)
    return output_module.smodelsOutput

def runMain( filename, timeout = 0, suppressStdout=True, development=False,
             inifile = "testParameters.ini", overridedatabase = None ):
    """ run SModelS proper
    :param filename: slha file
    :param timeout: timeout for the operation, given in seconds
    :param suppressStdout: if True, then redirect stdout and stderr to /dev/null
    :param development: development run (FIXME what does that entail?)
    :param inifile: the config file to be used
    :param overridedatabase: if not None, then use the provided database,
           else use databaseLoader.database
    :returns: printer output
    """
    to = None
    oldlevel = getLogLevel()
    level = 'debug'
    if suppressStdout:
        level = 'error'
        to = os.devnull
    database = None
    if overridedatabase != None:
        database = overridedatabase
    else:
        from databaseLoader import database ## to make sure the db exists
    with redirector.stdout_redirected ( to = to ):
        out = join( iDir(), "test/unitTestOutput" )
        setLogLevel ( level )
        run(filename, parameterFile=join ( iDir(), "test/%s" % inifile ),
             outputDir= out, db= database, timeout = timeout,
             development = development)
        setLogLevel ( oldlevel )
        sfile = join(iDir(),"test/unitTestOutput/%s.py" % basename(filename))
        return sfile


def compareSummaries(outA,outB,allowedDiff):


    fA = np.genfromtxt(outA,dtype=None,encoding='utf-8',
                  skip_header=3,names=True)

    fB = np.genfromtxt(outB,dtype=None,encoding='utf-8',
                  skip_header=3,names=True)

    if sorted(fA['filename']) != sorted(fB['filename']):
        logger.error("Filenames differ:\n %s\n and\n %s" %(sorted(fA['filename']),sorted(fB['filename'])))
        return False

    for fname in fA['filename']:
        ptA = fA[fA['filename'] == fname][0]
        ptB = fB[fB['filename'] == fname][0]
        for col in fA.dtype.names:
            if ptA[col] == ptB[col]:
                continue
            elif isinstance(ptA[col],(float,int)):
                diff = 2.*abs(ptA[col]-ptB[col])/abs(ptA[col]+ptB[col])
                if diff > allowedDiff:
                    logger.error("values for %s differ by %s in %s" % ( col, diff, fname) )
                    return False
            else:
                logger.error("values for %s differ in %s" % ( col, fname))
                return False
    return True
