"""
.. module:: crossSection
   :synopsis: Encapsulates the result of the computation of the reference
              cross section.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import TeV, pb, GeV
from smodels.theory import lheReader
import pyslha
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger

## orders in perturbation theory
LO,NLO,NLL,NNLL = range(4)

def orderToString ( order, short=False, raise_error=False ):
    """ return the string that describes the perturbation order
    :param short: if true, return a short version of string
    :param raise_error: if true, raise exception if order is not know
    """
    if order == LO:
        return "LO"
    if order == NLO:
        return "NLO"
    if order == NLL:
        if short:
            return "NLL"
        return "NLO+NLL"
    if order == NNLL:
        if short:
            return "NNLL"
        return "NLO+NLL+NNLL"
    if raise_error:
        line = "Unknown QCD order %d" % order
        raise SModelSError ( line )
    return "?"

def stringToOrder ( strng ):
    """ from a string describing the order return the perturbation order """
    order={ "LO": LO, "NLO": NLO, "NLL": NLL, "NNLL": NNLL,
            "NLO+NLL": NLL, "NLO+NLL+NNLL": NNLL }
    if strng in order:
        return order[strng]
    return -1

class XSectionInfo(object):
    """
    An instance of this class represents information regarding a cross section.

    This class is used to store information of a cross section (center of
    mass, order and label).

    """
    def __init__(self, sqrts = None, order = None, label = None ):
        """
        Constructor.
        :param: sqrts  center of mass energy, with unit
        :param: order perturbation order of xsec computation
        :param: label, a string that describes the xsec computation
        """
        self.sqrts = sqrts
        self.order = order
        self.label = label


    def __eq__(self, other):
        if type(other) != type(self):
            return False
        if other.sqrts != self.sqrts:
            return False
        if other.order != self.order:
            return False
        return True

    def __hash__( self ):
        ## for now lets just make it a string
        order = self.order
        if order == None:
            order=9
        return int(self.sqrts.asNumber(GeV)) + order

    def __str__(self):
        if not self.order:
            return str(self.sqrts)
        return "%s (%s)" % ( self.sqrts, self.order )

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
        Generate an independent copy of self.

        Faster than deepcopy.

        """
        newinfo = XSectionInfo()
        newinfo.sqrts = self.sqrts
        newinfo.order = self.order
        newinfo.label = self.label[:]
        return newinfo


class XSection(object):
    """
    An instance of this class represents a cross section.

    This class is used to store the information of a single cross section
    (value, particle ids, center of mass, order and label).

    order = 0 (LO), 1 (NLO), 2 (NLL), or 3 (NNLL).

    """
    def __init__(self):
        """
        Initializes the object to store a cross section value.
        All initial info is set to None.
        """
        self.info = XSectionInfo()
        self.value = None
        self._pid = (None, None)

    @property
    def pid(self):
        return self._pid

    @pid.setter
    def pid(self,pn):
        if None in pn:
            self._pid = pn
            return
        self._pid = tuple ( sorted (pn) )

    def __mul__(self, other):
        """
        Multiplies the value of the cross section by the factor other (should
        be a float or int).
        """

        newXsec = self.copy()
        if isinstance(other,float) or isinstance(other,int):
            newXsec.value = newXsec.value * float(other)
        else:
            print (other,type(other))
            logger.error("Xsections can only be multiplied by floats")
            raise SModelSError()
        return newXsec


    def __rmul__(self, other):
        """
        Right multiplication (see left multiplication).
        """
        return self * other


    def __add__(self, other):
        """
        Returns a copy of self with the value of other added to its value.
        """

        if type(other) == type(XSection()):
            if self.info == other.info:
                res = self.copy()
                res.value += other.value
                return res
        logger.error("Trying to add " + type(other) + " to a XSection object")
        raise SModelSError()


    def __eq__(self, other):
        """
        Compare two XSection objects. Returns True if .info and type and value and
        pid are equal.
        """

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
        """
        Compare two XSection objects. Returns True if .info or type or value or
        pid is not equal.
        """

        if type(other) != type(XSection()):
            return True
        if other.info != self.info:
            return True
        if other.value != self.value:
            return True
        if other.pid != self.pid:
            return True
        return False


    def __str__(self):
        """
        Generate cross section information in string format.
        """
        st = self.info.label + ':' + str(self.value)+ " " + str(self.pid)
        return st

    def niceStr(self):
        """
        Generates a more human readable string. The string format is:
        Sqrts: self.info.sqrts,  Weight: self.value
        """
        st = 'Sqrts: '+str(self.info.sqrts) + ', Weight:' + str(self.value)
        return st

    def copy(self):
        """
        Generates an independent copy of self.

        Faster than deepcopy.

        """
        newXsec = XSection()
        newXsec.info = self.info.copy()
        newXsec.value = self.value
        newXsec.pid = tuple(list(self.pid)[:])
        return newXsec


    def _zeroXSec(self):
        """
        Replace the cross section value by zero.

        """
        self.value = 0. * pb


class XSectionList(object):
    """
    An instance of this class represents a list of cross sections.

    This class is used to store a list of cross sections.
    The list is sorted by cross section, highest cross section first.

    """
    def __init__(self, infoList=None):
        """
        If infoList is defined, create entries with zero cross sections
        according to infoList. infoList must be a list of XSectionInfo objects.

        """
        self.xSections = []

        if infoList:
            for info in infoList:
                newentry = XSection()
                newentry.value = 0. * pb
                newentry.pid = (None, None)
                newentry.info = info.copy()
                self.add(newentry)

    def __mul__(self, other):

        newList = self.copy()
        for ixsec, xsec in enumerate(newList):
            newList[ixsec] = xsec * other
        return newList


    def __rmul__(self, other):

        return self*other

    def __add__(self,other):

        newList = self.copy()
        if type(other) != type(self):
            logger.warning("Trying to add a XSectionList and a "+str(type(other)))
            return self

        newList += other

        return  newList

    def __iadd__(self, newXsecs):
        """
        Add a new list of cross sections.

        If the new cross sections already appear (have same order and sqrts),
        add its value to the original value, otherwise append it to the list.
        The particle IDs are ignored when adding cross sections. Hence, they
        are set to (None, None) if any cross sections are combined.

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
                        if newXsec.pid != oldXsec.pid:
                            oldXsec.pid = (None, None)

        return self

    def __iter__(self):
        return iter(self.xSections)


    def __getitem__(self, index):
        if len(self)<=index:
            txt = "Index in XSectionList out of bounds: idx(%d)>=length(%d). " \
                   "(Are there cross sections given in the input?)" % \
                   ( index, len(self) )
            logger.error ( txt )
            raise SModelSError ( txt )
        return self.xSections[index]


    def __setitem__(self, index, xsec):
        if type(xsec) != type(XSection()):
            logger.error("Input object must be a XSection() object")
            raise SModelSError()
        else:
            self.xSections[index] = xsec


    def __len__(self):
        return len(self.xSections)


    def __str__(self):
        return str([str(xsec) for xsec in self])

    def niceStr(self):
        st = ""
        for xsec in self:
            st += xsec.niceStr()+'\n'
        return st


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
            raise SModelSError()
        else:
            self.xSections.append(newxsec.copy())

    def _addValue(self, newxsec):
        """
        Add a XSection object to the list.

        If the XSection object already exists, add to its values, otherwise
        append the object.

        """
        if type(newxsec) != type(XSection()):
            logger.error("Input object must be a XSection() object")
            raise SModelSError()
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
        Return a list of XSection objects for item (label, pid, sqrts).

        """
        xsecList = XSectionList()
        for xsec in self:
            if type(item) == type(xsec.info.label) and item == xsec.info.label:
                xsecList.add(xsec)
            elif type(item) == type(xsec.info.sqrts) \
                    and abs((item - xsec.info.sqrts).asNumber(TeV)) < 1e-2:
                xsecList.add(xsec)
            elif type(item) == type(xsec.pid) and item == xsec.pid:
                xsecList.add(xsec)
            elif type(item) == type(1) and (item in xsec.pid):
                xsecList.add(xsec)
        return xsecList


    def _zeroXSecs(self):
        """
        Replace the cross section values in the list by zero.

        """
        for xsec in self:
            xsec.value = 0. * pb


    def delete(self, xSec):
        """
        Delete the cross section entry from the list.

        """
        for ixsec, xsec in enumerate(self):
            if xsec == xSec:
                self.xSections.pop(ixsec)


    def getInfo(self):
        """
        Get basic info about the cross sections appearing in the list (order,
        value and label).

        :returns: list of XSectionInfo objects

        """
        allInfo = []
        for xsec in self:
            info = xsec.info
            if not info in allInfo:
                allInfo.append(info)
        return allInfo


    def _getLabels(self):
        """
        Get all labels appearing in the list.

        """
        allLabels = []
        allInfo = self.getInfo()
        for info in allInfo:
            allLabels.append(info.label)
        return list(set(allLabels))


    def getPIDpairs(self):
        """
        Get all particle ID pairs appearing in the list.

        """
        allPidPairs = []
        for xsec in self:
            allPidPairs.append(xsec.pid)
        return list(set(allPidPairs))


    def getPIDs(self):
        """
        Get all particle IDs appearing in the list.

        """
        allPids = []
        for xsec in self:
            allPids.extend(xsec.pid)
        return sorted(list(set(allPids)))


    def getMaxXsec(self):
        """
        Get the maximum cross section value appearing in the list.

        """
        maxxsec = 0. * pb
        for xsec in self:
            if xsec.value > maxxsec:
                maxxsec = xsec.value
        return maxxsec

    def getMinXsec(self):
        """
        Get minimum cross section value appearing in the list.

        """
        if len(self) > 0:
            minxsec = self.xSections[0].value
        else: return False
        for xsec in self:
            if xsec.value < minxsec:
                minxsec = xsec.value
        return minxsec


    def getDictionary(self, groupBy="pids"):
        """
        Convert the list of XSection objects to a nested dictionary.

        First level keys are the particles IDs (if groupBy == pids) or labels
        (if groupBy == labels) and values are the cross section labels or
        particle IDs and the cross section value.

        """
        xSecDictionary = {}

        if groupBy == "pids":
            allPids = self.getPIDpairs()
            for pid in allPids:
                xSecDictionary[pid] = {}
                xSecs = self.getXsecsFor(pid)
                for xsec in xSecs:
                    xSecDictionary[pid][xsec.info.label] = xsec.value

        elif groupBy == "labels":
            allLabels = self._getLabels()
            for label in allLabels:
                xSecDictionary[label] = {}
                xSecs = self.getXsecsFor(label)
                for xsec in xSecs:
                    xSecDictionary[label][xsec.pid] = xsec.value

        return xSecDictionary

    def removeLowerOrder(self):
        """
        Keep only the highest order cross section for each process in the list.

        Remove order information and set default labels.

        """

        newList = XSectionList()
        for pids in self.getPIDpairs():
            xsecs = self.getXsecsFor(pids)
            for i, ixsec in enumerate(xsecs):
                newxsec = ixsec.copy()
                removeXsec = False
                isqrts = ixsec.info.sqrts
                iorder = ixsec.info.order
                # Check if the xsec appear with the same sqrts but at a higher
                # order
                for j, jxsec in enumerate(xsecs):
                    if i == j:
                        continue
                    jsqrts = jxsec.info.sqrts
                    jorder = jxsec.info.order
                    if jsqrts == isqrts and jorder > iorder:
                        removeXsec = True
                        break
                if not removeXsec:
                    # Erase cross section labels and information
                    newxsec.info.label = str(newxsec.info.sqrts)
                    newxsec.info.order = None
                    newList.add(newxsec)

        if len(self) != len(newList):
            logger.debug("Ignoring %i lower order cross sections",
                           (len(self) - len(newList)))
        self.xSections = newList.xSections

    def removeDuplicates(self):
        """
        If there are two entries for the same process, center of mass energy and order,
        keep only one (the one with highest value).
        """

        newList = XSectionList()
        for pids in self.getPIDpairs():
            xsecs = self.getXsecsFor(pids)
            #Make sure xsecs are sorted by sqrts,order and value:
            xsecs.xSections = sorted(xsecs.xSections,
                key=lambda xsec: (xsec.info.sqrts,xsec.info.order,xsec.value.asNumber(pb )))
            for i, ixsec in enumerate(xsecs):
                keepXsec = True
                isqrts = ixsec.info.sqrts
                iorder = ixsec.info.order
                # Check if there is a duplicate entry (with higher value):
                for j, jxsec in enumerate(xsecs):
                    if i >= j:
                        continue
                    jsqrts = jxsec.info.sqrts
                    jorder = jxsec.info.order
                    if jsqrts == isqrts and jorder == iorder:
                        keepXsec = False
                        break
                if keepXsec:
                    newList.add(ixsec.copy())

        if len(self) != len(newList):
            logger.warning("Removing %i duplicate cross sections",
                           (len(self) - len(newList)))
            self.xSections = newList.xSections

    def order(self):
        """
        Order the cross section in the list by their PDG pairs
        """

        self.xSections = sorted(self.xSections, key=lambda xsec: xsec.pid)

    def __lt__(self, other):
        return self.xSections[0].pid < other.xSections[0].pid

    def sort ( self ):
        """ sort the xsecs by the values """
        self.xSections = sorted(self.xSections,
                key=lambda xsec: xsec.value.asNumber(pb ), reverse=True )


def getXsecFromSLHAFile(slhafile, useXSecs=None, xsecUnit = pb):
    """
    Obtain cross sections for pair production of R-odd particles from input SLHA file.
    The default unit for cross section is pb.

    :parameter slhafile: SLHA input file with cross sections
    :parameter useXSecs: if defined enables the user to select cross sections to
                     use. Must be a XSecInfoList object
    :parameter xsecUnit: cross section unit in the input file (must be a Unum unit)
    :returns: a XSectionList object

    """
    # Store information about all cross sections in the SLHA file
    xSecsInFile = XSectionList()
    f=pyslha.readSLHAFile(slhafile)
    for production in f.xsections:
        process = f.xsections.get ( production )
        for pxsec in process.xsecs:
            csOrder = pxsec.qcd_order
            wlabel = str( int ( pxsec.sqrts / 1000) ) + ' TeV'
            wlabel += f" ({orderToString(csOrder,True,True)})"
            xsec = XSection()
            xsec.info.sqrts = pxsec.sqrts/1000. * TeV
            xsec.info.order = csOrder
            xsec.info.label = wlabel
            xsec.value = pxsec.value * pb
            xsec.pid = production[2:]
            # Do not add xsecs which do not match the user required ones:
            if (useXSecs and not xsec.info in useXSecs):
                continue
            else: xSecsInFile.add(xsec)

    #Make sure duplicates are removed.
    xSecsInFile.removeDuplicates()

    return xSecsInFile


def getXsecFromLHEFile(lhefile, addEvents=True):
    """
    Obtain cross sections from input LHE file.

    :parameter lhefile: LHE input file with unweighted MC events
    :parameter addEvents: if True, add cross sections with the same mothers,
                      otherwise return the event weight for each pair of mothers
    :returns: a XSectionList object

    """
    # Store information about all cross sections in the LHE file
    xSecsInFile = XSectionList()
    reader = lheReader.LheReader(lhefile)
    if not type ( reader.metainfo["totalxsec"] ) == type ( pb) :
        logger.error("cross section information not found in LHE file.")
        raise SModelSError()
    elif not reader.metainfo["nevents"]:
        logger.error("Total number of events information not found in LHE " +
                     "file.")
        raise SModelSError()
    elif not type ( reader.metainfo["sqrts"] ) == type ( TeV ):
        logger.error("Center-of-mass energy information not found in LHE " +
                     "file.")
        raise SModelSError()

    # Common cross section info
    totxsec = reader.metainfo["totalxsec"]
    nevts = reader.metainfo["nevents"]
    sqrtS = reader.metainfo["sqrts"]
    eventCs = totxsec / float(nevts)


    # Get all mom pids
    allpids = []
    for event in reader:
        allpids.append(tuple(sorted(event.getMom())))
    pids = set(allpids)
    # Generate zero cross sections for all independent pids
    for pid in pids:
        xsec = XSection()
        xsec.info.sqrts = sqrtS
        if 'cs_order' in reader.metainfo:
            xsec.info.order = reader.metainfo["cs_order"]
        else:
            # Assume LO xsecs, if not defined in the reader
            xsec.info.order = 0
        wlabel = str( sqrtS / TeV ) + ' TeV'
        wlabel += f' ({orderToString(xsec.info.order,True,False)})'
        xsec.info.label = wlabel
        xsec.value = 0. * pb
        xsec.pid = pid
        # If addEvents = False, set cross section value to event weight
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

    #Make sure duplicates are removed.
    xSecsInFile.removeDuplicates()

    return xSecsInFile
