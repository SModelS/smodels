#!/usr/bin/env python

"""
.. module:: CrossSection
        :synopsis: A class that encapsulates the result of the computation of the reference cross section.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from Tools.PhysicsUnits import addunit, rmvunit
import copy
import logging
logger = logging.getLogger(__name__)
import LHEReader

class XSectionInfo:
    """A simple class to store information about the cross-section (center of mass, order and label)"""
    def __init__ (self):
        self.sqrts = None
        self.order = None
        self.label = None

    def __eq__(self,other):     
        if type(other) != type(XSectionInfo()): return False
        if other.sqrts != self.sqrts: return False
        if other.order != self.order: return False
        if other.label != self.label: return False
        return True
 
    def __ne__(self,other):     
        if type(other) != type(XSectionInfo()): return True
        if other.sqrts != self.sqrts: return True
        if other.order != self.order: return True
        return False


class XSection:
    """A simple class to store the information about a single cross-section (value, paritcle ids, center of mass, order and label)
         order = 0 (LO), 1 (NLO) or 2 (NLL)."""
    def __init__ (self):
        self.info = XSectionInfo()
        self.value = None
        self.pid = (None,None)

    def __mul__(self,other):
        newXsec = copy.deepcopy(self)
        if type(other) == type(1.):
            newXsec.value = newXsec.value*other
        else:
            print "[XSection.mul]: Xsections can only be multiplied by floats"
            return False
        return newXsec
        
    def __add__(self,other):
        if type(other) == type(XSection()):
            if self.info == other.info:
                res = copy.deepcopy(self)
                res.value += other.value
                return res
        print "[XSection.add]: Trying to add",type(other),"to a XSection objetc"
        return False

    def __eq__(self,other):
        if type(other) != type(XSection()): return False
        if other.info != self.info: return False
        if other.value != self.value: return False
        if other.pid != self.pid: return False
        return True

    def __ne__(self,other):
        if type(other) != type(XSection()): return True
        if other.info != self.info: return True
        if other.value != self.value: return True
        if other.pid != self.pid: return True
        return False    

    def __str__ (self):
        """cross-section information in string format"""
        st = 'label: '+self.info.label+', value:'+str(self.value)
        return st
    
    def zeroXSec(self):
        """
        Replaces the cross-section value by zero
        """
        self.values = addunit(0., 'fb')

         
class XSectionList:
    """A simple class to store a list of cross-sections to be used"""        
    
    def __init__ (self,infoList=None):
        """Creates a list of XSection objects from the input string with None cross-section values.
        If infoList is defined, create entries with zero cross-sections according to infoList"""
        self.XSections=[]
        
        if infoList:
            for info in infoList:
                newentry = XSection()
                newentry.value = addunit(0.,'fb')
                newentry.pid = (None,None)
                newentry.info = copy.deepcopy(info)
                self.XSections.append(newentry)
                

    def __mul__(self,other):
        newList = copy.deepcopy(self)
        for ixsec,xsec in enumerate(newList.XSections): newList.XSections[ixsec] = xsec*other
        return newList
    
    def __iter__(self):
        return iter(self.XSections)

    def __str__(self):
        return str([str(xsec) for xsec in self])


    def getXsecsFor(self,item):
        """ Returns a list of XSection objects for item (label, pid, sqrts) """
        xsecList = XSectionList()
        for xsec in self:
            if type(item) == type(xsec.info.label) and item == xsec.info.label:
                xsecList.XSections.append(xsec)
            elif type(item) == type(xsec.info.sqrts) and item == xsec.info.sqrts:    
                xsecList.XSections.append(xsec)
            elif type(item) == type(xsec.pid) and item == xsec.pid:    
                xsecList.XSections.append(xsec)
            elif type(item) == type(1) and (item in xsec.pid):
                xsecList.XSections.append(xsec)

        return xsecList
    
    def zeroXSecs(self):
        """
        Replaces the cross-section values in the list by zero
        """
        for xsec in self:  xsec.value = addunit(0.,'fb')

    def delete(self,XSec):
        """Deletes the cross-section entry from the list"""
        for ixsec,xsec in enumerate(self.XSections):
            if xsec == XSec: self.XSections.pop(ixsec)

    def getInfo(self):
        """Gets the basic info about the cross-sections appearing in the list (order,value and label). Returns a list of XSectionInfo objects """
        allInfo = []
        for xsec in self:
            info = xsec.info
            if not info in allInfo: allInfo.append(info)
        return allInfo

    def getLabels(self):
        """Gets all the labels appearing in the list."""
        allLabels = []
        allInfo = self.getInfo()
        for info in allInfo: allLabels.append(info.label)
        return list(set(allLabels)) 

    def getPIDpairs(self):
        """Gets all the particle ID pairs appearing in the list."""
        allPidPairs = []
        for xsec in self: allPidPairs.append(xsec.pid)
        return list(set(allPidPairs)) 

    def getPIDs(self):
        """Gets all the particle IDs appearing in the list."""
        allPids = []
        for xsec in self: allPids.extend(xsec.pid)
        return list(set(allPids)) 


    def getMaxXsec(self):
        """Gets the maximum cross-section value appearing in the list."""
        maxxsec= addunit(0.,'fb')
        for xsec in self:
            if xsec.value > maxxsec: maxxsec = xsec.value
        return maxxsec


    def getDictionary(self,groupBy="pids"):
        """ Converts the list of XSection objects to a nested dictionary. First level keys are the particles IDs (if groupBy=pids) or labels
        (if groupBy=labels) and values are the cross-section labels or particle IDs and the cross-section value.
        If groupBy=pids and a single pid is present, return a simple dictionary with the cross-sections for the pid"""

        XsecDic = {}

        if groupBy == "pids":
            allPids = self.getPIDpairs()
            for pid in allPids:
                XsecDic[pid] = {}
                Xsecs = self.getXsecsFor(pid)
                for xsec in Xsecs: XsecDic[pid][xsec.info.label] = xsec.value
            if len(allPids) == 1: XsecDic = XsecDic[allPids[0]]  #Return standard weight dictionary                

        elif groupBy == "labels":
            allLabels = self.getLabels()
            for label in allLabels:
                XsecDic[label] = {}
                Xsecs = self.getXsecsFor(label)
                for xsec in Xsecs:
                    XsecDic[label][xsec.pid] = xsec.value
                    

        return XsecDic
    
    def combineWith(self,newList):
        """ Adds a new list of cross-sections. If the new cross-sections already appear (have same order and sqrts), add its
        value to the original value, otherwise append it to the list.
        The particle IDs are ignored when adding cross-sections. Hence they are set to (None,None) if any cross-sections are combined"""
        
        for newXsec in newList:
            if not newXsec.info in self.getInfo():
                self.XSections.append(newXsec)                
            else:    
                for oldXsec in self:
                    if newXsec.info == oldXsec.info:
                        oldXsec.value = copy.deepcopy(oldXsec.value + newXsec.value)
                        oldXsec.pid = (None,None)
        

def getXsecFromSLHAFile(slhafile,UseXSecs=None):
    """ obtain cross-sections from input SLHA file. 
        :param slhafile: SLHA input file with cross-sections
        :param UseXsecs: if defined enables the user to select cross-sections to use./
        Must be a XSecInfoList object
        :returns: a XSectionList object     
    """


    XsecsInFile = XSectionList()    #To store information about all cross-sections in the SLHA file
    slha = open(slhafile, 'r')
    lines = slha.readlines()
    xsecblock = False
    for l in lines:
        if l.startswith("#") or len(l)<2: continue
        if 'XSECTION' in l:
            xsecblock = True
            sqrtS =    eval(l.split()[1])/1000.#Values in the SLHA file are in GeV
            pids = (eval(l.split()[5]),eval(l.split()[6]))
            continue
        if not xsecblock: continue    #ignores other entries
        cs_order = eval(l.split()[1])
        cs = addunit(eval(l.split()[6]),'fb')
        wlabel = str(int(sqrtS))+' TeV'
        if cs_order == 0: wlabel += ' (LO)'
        elif cs_order == 1: wlabel += ' (NLO)'
        elif cs_order == 2: wlabel += ' (NLL)'
        else:
            print '[SLHADecomposer] unknown QCD order in XSECTION line', l
            return False
        xsec = XSection()
        xsec.info.sqrts = addunit(sqrtS,'TeV')
        xsec.info.order = cs_order
        xsec.info.label = wlabel
        xsec.value = cs
        xsec.pid = pids
#Do not add xsecs which do not match the user required ones:        
        if UseXSecs and not xsec.info in UseXSecs: continue 
        else: XsecsInFile.XSections.append(xsec)

    slha.close()

    return XsecsInFile


def getXsecFromLHEFile(lhefile):
    """ obtain cross-sections from input LHE file. 
        :param lhefile: LHE input file with unweighted MC events        
        :returns: a XSectionList object     
    """                       

    XsecsInFile = XSectionList()    #To store information about all cross-sections in the LHE file
    reader = LHEReader.LHEReader(lhefile)
    if not rmvunit(reader.metainfo["totalxsec"],'fb'):
        logger.error("Cross-section information not found in LHE file.")
        return False
    elif not reader.metainfo["nevents"]:
        logger.error("Total number of events information not found in LHE file.")
        return False
    elif not rmvunit(reader.metainfo["sqrts"],'TeV'):
        logger.error("Center-of-mass energy information not found in LHE file.")
        return False    

#Common cross-section info
    totxsec = reader.metainfo["totalxsec"]
    nevts = reader.metainfo["nevents"]
    sqrtS = reader.metainfo["sqrts"]
    event_cs = totxsec/float(nevts)
    
#Get all mom pids:
    allpids = []    
    for Event in reader: allpids.append(tuple(Event.getMom()))
    allpids = set(allpids)    
    for pids in allpids:
        xsec = XSection()
        xsec.info.sqrts = sqrtS
        if reader.metainfo.has_key("cs_order"): xsec.info.order = reader.metainfo["cs_order"]
        else: xsec.info.order = 0  #Assume LO xsecs, if not defined in the reader
        wlabel = str(int(rmvunit(sqrtS,'TeV')))+' TeV'
        if xsec.info.order == 0: wlabel += ' (LO)'
        elif xsec.info.order == 1: wlabel += ' (NLO)'
        elif xsec.info.order == 2: wlabel += ' (NLL)'
        xsec.info.label = wlabel
        xsec.value = event_cs
        xsec.pid = pids
        XsecsInFile.XSections.append(xsec)

    reader.close()

    return XsecsInFile
