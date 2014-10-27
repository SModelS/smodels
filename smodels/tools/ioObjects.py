#!/usr/bin/env python

"""
    .. module:: ioObjects
    :synopsis: Definitions of input/output parameters which are read from parameter.in
    
    .. moduleauthor:: Ursula Laa <Ursula.Laa@assoc.oeaw.ac.at>    
    .. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
"""

import os
from smodels.theory import lheReader
from smodels.theory.printer import Printer

class ResultList(Printer):
    """
    Class that collects experimental constraints and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresultonly = None, describeTopo = None):
        self.outputarray = outputarray
        self.bestresultonly = bestresultonly
        self.describeTopo = describeTopo

    def addResult(self, res, maxcond):
        mCond = res.getmaxCondition()
        if mCond == 'N/A': return
        if mCond > maxcond: return
        self.outputarray.append(res)
        return

    def getR(self, res):
        #calculate R value
        return res.value[0].value / res.analysis.getUpperLimitFor(res.mass)

    def sort(self):
        #reverse sort outputarray by R value
        self.outputarray = sorted(self.outputarray, key=self.getR, reverse=True)
        return

    def getBestResult(self):
        #return best result
        self.sort()
        return self.outputarray[0]

    def isEmpty(self):
        #check if outputarray is empty
        return len(self.outputarray)==0

    def formatData(self):
        """
        to access printout format
        """
        return self.formatResultsData()

class OutputStatus(Printer):
    """
    Object that holds all status information and has a predefined printout 
    """
    def __init__(self, slhastatus, warnings, databaseVersion):
        self.slhastatus = slhastatus
        self.warnings = warnings
        self.databaseVersion = databaseVersion
        self.statusStrings = {-1: "#could not run the decomposition",
                              -3: "#no cross sections above sigmacut found",
                              -4: "#database not found",
                              -2: "#bad input slha, did not run decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}
        self.status = self.initiateStatus()

    def initiateStatus(self):
        """
        evaluate initial status from slhastatus
        """
        if not self.databaseVersion: return -4
        if self.slhastatus == -1 or self.slhastatus == -3: return -2
        return 0

    def updateStatus(self, status):
        self.status = status
        return

    def addWarning(self, warning):
        self.warnings += warning
        return

    def formatData(self):
        """
        to access printout format
        """
        return self.formatStatusData()


class LheStatus(Printer):
    """
    Object to check if input lhe file is ok
    """

    def __init__(self, filename, massgap=None, maxcond=None):
        self.filename = filename
        #additional information for printer output
        self.massgap = massgap
        self.maxcond = maxcond
        self.status = self.evaluateStatus()

    def formatData(self):
        return self.formatLHEData()

    def evaluateStatus(self):
        if not os.path.exists(self.filename):
            #set status flag to -3, as in slha checks for missing input file
            return -3, "Inputfile %s not found" %self.filename
        lhe = lheReader.LheReader(self.filename)
        if not lhe.metainfo["nevents"]:
            return -1, "No events found in inputfile %s" %self.filename
        return 1, "Input file ok"
