#!/usr/bin/env python

"""
.. module:: testRunSModelS
   :synopsis: Tests runSModelS

.. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>

"""

import sys,shutil,os
sys.path.insert(0,"../")
import unittest
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from runSModelS import main
import logging as logger
import unum


def equalObjs(obj1,obj2,allowedDiff,ignore=[]):
    """
    Compare two objects.
    The numerical values are compared up to the precision defined by allowedDiff.
    
    :param obj1: First python object to be compared 
    :param obj2: Second python object to be compared
    :param allowedDiff: Allowed % difference between two numerical values 
    :param ignore: List of keys to be ignored
    :return: True/False    
    """      
    
    if type(obj1) != type(obj2):
        logger.info("Data types differ (%s,%s)" %(type(obj1),type(obj2)))
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
        return diff < allowedDiff
    elif isinstance(obj1,str):
        return obj1 == obj2
    elif isinstance(obj1,dict):    
        if len(obj1) != len(obj2):
            error = "Dictionaries have distinct lengths (%i,%i)" %(len(obj1),len(obj2))
            return error
        for key in obj1:
            if key in ignore: continue
            if not key in obj2:
                error = "Key %s missing" %key
                return error
            if not equalObjs(obj1[key],obj2[key],allowedDiff):
                error = 'Objects differ:\n   %s\n and\n   %s' %(str(obj1[key]),str(obj2[key]))
                return error
    elif isinstance(obj1,list):
        for ival,val in enumerate(sorted(obj1)):
            if not equalObjs(val,sorted(obj2)[ival],allowedDiff):
                error = 'Objects differ:\n   %s \n and\n   %s' %(str(val),str(sorted(obj2)[ival]))
                return error
    else:
        return obj1 == obj2
            
    return True


class RunSModelSTest(unittest.TestCase):
    def runMain(self, filename ):
        suppressStdout = True
        if suppressStdout:
            a=sys.stdout
            sys.stdout = open ( "stdout.log", "w" )
        out = join ( iDir(), "test/unitTestOutput" )
        main(filename, parameterFile=join ( iDir(), "test/testParameters.ini" ),
             outputDir= out, verbosity = 'error' )
        sfile = join ( iDir(), 
                "test/unitTestOutput/%s.py" % basename ( filename ) )
        return sfile

    def testGoodFile(self):        
        filename = join ( iDir(), "inputFiles/slha/gluino_squarks.slha" )
        outputfile = self.runMain(filename )
        shutil.copyfile(outputfile,'./output.py')
        from gluino_squarks_default import smodelsOutputDefault
        from output import smodelsOutput
        ignoreFields = ['input file','smodels version']
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.01,ignore=ignoreFields)
        self.assertEqual(equals,True)
        os.remove('./output.py')
        os.remove('./output.pyc')

    def testBadFile(self):
   
        filename = join ( iDir(), "inputFiles/slha/I_dont_exist.slha" )
        outputfile = self.runMain (filename )
        shutil.copyfile(outputfile,'./bad_output.py')
        from bad_default import smodelsOutputDefault
        from bad_output import smodelsOutput
        ignoreFields = ['input file','smodels version']
        equals = equalObjs(smodelsOutput,smodelsOutputDefault,allowedDiff=0.,ignore=ignoreFields)
        self.assertEqual(equals,True)
        os.remove('./bad_output.py')
        os.remove('./bad_output.pyc')

if __name__ == "__main__":
    unittest.main()
