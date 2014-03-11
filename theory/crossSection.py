"""
.. module:: theory.crossSection
   :synopsis: A class that encapsulates the result of the computation of the
   reference cross section.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from tools.PhysicsUnits import addunit, rmvunit
import logging
import LHEReader

logger = logging.getLogger(__name__)


class XSectionInfo(object):
    """
    A simple class to store information about the cross-section (center of
    mass, order and label).
    
    """
    def __init__ (self):
        self.sqrts = None
        self.order = None
        self.label = None


    def __eq__(self, other):     
        if type(other) != type(XSectionInfo()):
            return False
        if other.sqrts != self.sqrts:
            return False
        if other.order != self.order:
            return False
        if other.label != self.label:
            return False
        return True
 
 
    def __ne__(self, other):
        if type(other) != type(XSectionInfo()):
            return True
        if other.sqrts != self.sqrts:
            return True
        if other.order != self.order:
            return True
        return False
    
    
    def copy(self):
        """
        Generates an independent copy of itself. Faster than deepcopy.
        
        """
        newinfo = XSectionInfo()
        newinfo.sqrts = self.sqrts
        newinfo.order = self.order
        newinfo.label = self.label[:]
        return newinfo


class XSection(object):
    """
    A simple class to store the information about a single cross-section
    (value, paritcle ids, center of mass, order and label).
    
    order = 0 (LO), 1 (NLO) or 2 (NLL).
    
    """
    def __init__ (self):
        self.info = XSectionInfo()
        self.value = None
        self.pid = (None, None)


    def __mul__(self, other):
        newXsec = self.copy()
        if type(other) == type(1.):
            newXsec.value = newXsec.value*other
        else:
            logger.error("Xsections can only be multiplied by floats")
            return False
        return newXsec
    
    
    def __rmul__(self, other):
        return self*other
        
        
    def __add__(self, other):
        if type(other) == type(XSection()):
            if self.info == other.info:
                res = self.copy()
                res.value += other.value
                return res
        logger.error("Trying to add", type(other), "to a XSection objetc")
        return False


    def __eq__(self, other):
        if type(other) != type(XSection()):
            return False
        if other.info != self.info:
            return False
        if other.value != self.value:
            return False
        if other.pid != self.pid:
            return False
        return True


    def __ne__(self, other):
        if type(other) != type(XSection()):
            return True
        if other.info != self.info:
            return True
        if other.value != self.value:
            return True
        if other.pid != self.pid:
            return True
        return False    


    def __str__ (self):
        """
        Cross-section information in string format.
        
        """
        st = self.info.label+':'+str(self.value)
        return st
    
    
    def copy(self):
        """
        Generates an independent copy of itself. Faster than deepcopy.
        
        """
        newXsec = XSection()
        newXsec.info = self.info.copy()
        newXsec.value = self.value
        newXsec.pid = tuple(list(self.pid)[:])
        return newXsec
    
    
    def zeroXSec(self):
        """
        Replaces the cross-section value by zero.
        
        """
        self.value = addunit(0., 'fb')

         
class XSectionList(object):
    """
    A simple class to store a list of cross-sections to be used.
    
    """
    def __init__ (self, infoList=None):
        """
        Creates a list of XSection objects from the input string with None
        cross-section values. If infoList is defined, create entries with zero
        cross-sections according to infoList.
        
        """
        self.xSections=[]
        
        if infoList:
            for info in infoList:
                newentry = XSection()
                newentry.value = addunit(0., 'fb')
                newentry.pid = (None, None)
                newentry.info = info.copy()
                self.add(newentry)
                

    def __mul__(self, other):
        newList = self.copy()
        for ixsec, xsec in enumerate(newList):
            newList[ixsec] = xsec*other
        return newList
    
    
    def __rmul__(self, other):
        return self*other
    
    
    def __iter__(self):
        return iter(self.xSections)
    
    
    def __getitem__(self, index):
        return self.xSections[index]
    
    
    def __setitem__(self, index, xsec):
        if type(xsec) != type(XSection()):
            logger.error("Input object must be a XSection() object")
            return False
        else:
            self.xSections[index] = xsec
    
    
    def __len__(self):
        return len(self.xSections)        


    def __str__(self):
        return str([str(xsec) for xsec in self])
    
    
    def copy(self):
        """
        Generates an independent copy of itself. Faster than deepcopy.
        
        """
        newList = XSectionList()
        for xsec in self.xSections:
            newList.xSections.append(xsec.copy())
        return newList    
    
    
    def add(self, newxsec):
        """
        Append a XSection object to the list.
        
        """
        if type(newxsec) != type(XSection()):
            logger.error("Input object must be a XSection() object")
            return False
        else:
            self.xSections.append(newxsec.copy())
            
        
    def addValue(self, newxsec):
        """
        Adds a XSection object to the list. If the XSection object already
        exists, add to its values, otherwise append the object.
        
        """
        if type(newxsec) != type(XSection()):
            logger.error("Input object must be a XSection() object")
            return False
        else:
            exists = False
            for iXSec, xSec in enumerate(self.xSections):
                if xSec.info == newxsec.info \
                        and sorted(xSec.pid) == sorted(newxsec.pid):
                    self.xSections[iXSec].value = xSec.value + newxsec.value
                    break
            if not exists:
                self.add(newxsec)    


    def getXsecsFor(self, item):
        """
        Returns a list of XSection objects for item (label, pid, sqrts).
        
        """
        xsecList = XSectionList()
        for xsec in self:
            if type(item) == type(xsec.info.label) and item == xsec.info.label:
                xsecList.add(xsec)
            elif type(item) == type(xsec.info.sqrts) \
                    and item == xsec.info.sqrts:
                xsecList.add(xsec)
            elif type(item) == type(xsec.pid) and item == xsec.pid:
                xsecList.add(xsec)
            elif type(item) == type(1) and (item in xsec.pid):
                xsecList.add(xsec)
        return xsecList
    
    
    def zeroXSecs(self):
        """
        Replaces the cross-section values in the list by zero.
        
        """
        for xsec in self:
            xsec.value = addunit(0., 'fb')


    def delete(self, xSec):
        """
        Deletes the cross-section entry from the list.
        
        """
        for ixsec, xsec in enumerate(self):
            if xsec == xSec:
                self.xSections.pop(ixsec)


    def getInfo(self):
        """
        Gets the basic info about the cross-sections appearing in the list
        (order, value and label). Returns a list of XSectionInfo objects.
        
        """
        allInfo = []
        for xsec in self:
            info = xsec.info
            if not info in allInfo:
                allInfo.append(info)
        return allInfo


    def getLabels(self):
        """
        Gets all the labels appearing in the list.
        
        """
        allLabels = []
        allInfo = self.getInfo()
        for info in allInfo:
            allLabels.append(info.label)
        return list(set(allLabels)) 


    def getPIDpairs(self):
        """
        Gets all the particle ID pairs appearing in the list.
        
        """
        allPidPairs = []
        for xsec in self:
            allPidPairs.append(xsec.pid)
        return list(set(allPidPairs)) 


    def getPIDs(self):
        """
        Gets all the particle IDs appearing in the list.
        
        """
        allPids = []
        for xsec in self:
            allPids.extend(xsec.pid)
        return list(set(allPids)) 


    def getMaxXsec(self):
        """
        Gets the maximum cross-section value appearing in the list.
        
        """
        maxxsec = addunit(0., 'fb')
        for xsec in self:
            if xsec.value > maxxsec:
                maxxsec = xsec.value
        return maxxsec


    def getDictionary(self, groupBy="pids"):
        """
        Converts the list of XSection objects to a nested dictionary. First
        level keys are the particles IDs (if groupBy=pids) or labels
        (if groupBy=labels) and values are the cross-section labels or particle
        IDs and the cross-section value. If groupBy=pids and a single pid is
        present, return a simple dictionary with the cross-sections for the
        pid.
        
        """
        xSecDictionary = {}

        if groupBy == "pids":
            allPids = self.getPIDpairs()
            for pid in allPids:
                xSecDictionary[pid] = {}
                xSecs = self.getXsecsFor(pid)
                for xsec in xSecs:
                    xSecDictionary[pid][xsec.info.label] = xsec.value
            if len(allPids) == 1:
                # Return standard weight dictionary
                xSecDictionary = xSecDictionary[allPids[0]]                   

        elif groupBy == "labels":
            allLabels = self.getLabels()
            for label in allLabels:
                xSecDictionary[label] = {}
                xSecs = self.getXsecsFor(label)
                for xsec in xSecs:
                    xSecDictionary[label][xsec.pid] = xsec.value

        return xSecDictionary
    
    
    def combineWith(self, newXsecs):
        """
        Adds a new list of cross-sections. If the new cross-sections already
        appear (have same order and sqrts), add its value to the original
        value, otherwise append it to the list. The particle IDs are ignored
        when adding cross-sections. Hence they are set to (None, None) if any
        cross-sections are combined.
        
        """        
        newList = newXsecs
        if type(newXsecs) == type(XSection()):
            newList = [newXsecs]
        for newXsec in newList:
            if not newXsec.info in self.getInfo():
                self.add(newXsec)
            else:
                for oldXsec in self:
                    if newXsec.info == oldXsec.info:
                        oldXsec.value = oldXsec.value + newXsec.value
                        oldXsec.pid = (None, None)
        

def getXsecFromSLHAFile(slhafile, useXSecs=None):
    """
    Obtain cross-sections from input SLHA file. 
    
    :param slhafile: SLHA input file with cross-sections
    :param UseXsecs: if defined enables the user to select cross-sections to
    use. Must be a XSecInfoList object
    :returns: a XSectionList object    
     
    """
    # To store information about all cross-sections in the SLHA file
    xSecsInFile = XSectionList()    
    slha = open(slhafile, 'r')
    lines = slha.readlines()
    xsecblock = False
    for l in lines:
        if l.startswith("#") or len(l)<2:
            continue
        if 'XSECTION' in l:
            xsecblock = True
            # Values in the SLHA file are in GeV
            sqrtS =    eval(l.split()[1])/1000.
            pids = (eval(l.split()[5]), eval(l.split()[6]))
            continue
        if not xsecblock:
            continue    #ignores other entries
        csOrder = eval(l.split()[1])
        cs = addunit(eval(l.split()[6]), 'fb')
        wlabel = str(int(sqrtS))+' TeV'
        if csOrder == 0:
            wlabel += ' (LO)'
        elif csOrder == 1:
            wlabel += ' (NLO)'
        elif csOrder == 2:
            wlabel += ' (NLL)'
        else:
            logger.error("Unknown QCD order in XSECTION line", l)
            return False
        xsec = XSection()
        xsec.info.sqrts = addunit(sqrtS, 'TeV')
        xsec.info.order = csOrder
        xsec.info.label = wlabel
        xsec.value = cs
        xsec.pid = pids
        # Do not add xsecs which do not match the user required ones:
        if useXSecs and not xsec.info in useXSecs:
            continue 
        else:
            xSecsInFile.add(xsec)

    slha.close()

    return xSecsInFile


def getXsecFromLHEFile(lhefile, addEvents=True):
    """
    Obtain cross-sections from input LHE file.
    
    :param lhefile: LHE input file with unweighted MC events
    :param addEvents: if True, add cross-sections with the same mothers,
    otherwise return the event weight for each pair of mothers
    :returns: a XSectionList object
    
    """
    # To store information about all cross-sections in the LHE file
    xSecsInFile = XSectionList()    
    reader = LHEReader.LHEReader(lhefile)
    if not rmvunit(reader.metainfo["totalxsec"],'fb'):
        logger.error("Cross-section information not found in LHE file.")
        return False
    elif not reader.metainfo["nevents"]:
        logger.error("Total number of events information not found in LHE " +
                     "file.")
        return False
    elif not rmvunit(reader.metainfo["sqrts"],'TeV'):
        logger.error("Center-of-mass energy information not found in LHE " +
                     "file.")
        return False    

    # Common cross-section info
    totxsec = reader.metainfo["totalxsec"]
    nevts = reader.metainfo["nevents"]
    sqrtS = reader.metainfo["sqrts"]
    eventCs = totxsec/float(nevts)
    
    # Get all mom pids
    allpids = []    
    for event in reader:
        allpids.append(tuple(sorted(event.getMom())))
    pids = set(allpids)        
    # Generate zero cross-sections for all independent pids    
    for pid in pids:
        xsec = XSection()
        xsec.info.sqrts = sqrtS
        if reader.metainfo.has_key("cs_order"):
            xsec.info.order = reader.metainfo["cs_order"]
        else:
            # Assume LO xsecs, if not defined in the reader
            xsec.info.order = 0  
        wlabel = str(int(rmvunit(sqrtS,'TeV')))+' TeV'
        if xsec.info.order == 0:
            wlabel += ' (LO)'
        elif xsec.info.order == 1:
            wlabel += ' (NLO)'
        elif xsec.info.order == 2:
            wlabel += ' (NLL)'
        xsec.info.label = wlabel
        xsec.value = addunit(0.,'fb')        
        xsec.pid = pid
        # If addEvents = False, set cross-section value to event weight
        if not addEvents:
            xsec.value = eventCs 
        xSecsInFile.add(xsec)      
    # If addEvents = True, sum up event weights with same mothers
    if addEvents:
        for pid in allpids:
            for ixsec, xsec in enumerate(xSecsInFile.xSections):
                if xsec.pid == pid:
                    xSecsInFile.xSections[ixsec].value += eventCs
    
    reader.close()

    return xSecsInFile
