"""
.. module:: unitTestHelpers
   :synopsis: helper functions for the unit tests

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
sys.path.append('../')
import unum
import re
import pyslha
from xml.etree import ElementTree
import numpy as np
import redirector
from smodels.tools.runSModelS import run
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.base.smodelsLogging import logger, setLogLevel, getLogLevel
from smodels.base.physicsUnits import fb
from smodels.decomposition.theorySMS import TheorySMS
from smodels.experiment.expSMS import ExpSMS

def canonNameToVertNumb(topoDict,cName):
    """
    Convert canonical name to old vertnumb notation (in string format)

    :param topoDict: TopologyDict object
    :param cName: Canonical name for the desired topology

    :return: Vertnumb in string format
    """

    if cName not in topoDict:
        return None
    sms = topoDict[cName][0]
    evenParticles = sms.treeToBrackets()[0]
    vertnumb = str([len(v) for v in evenParticles[0]])
    vertnumb += str([len(v) for v in evenParticles[1]])
    vertnumb = vertnumb.replace(' ','')

    return vertnumb

def theorySMSFromString(stringEl,model,prodXSec = 1.0*fb,
                              maxWeight = 1.0*fb, intermediateState=None,
                              finalState=None):

    """
    Facility to construct a TheorySMS object from a string (only needed for unit tests)
    """

    # Hack to create a theory element from a string:
    expSMS = ExpSMS.from_string(stringEl, model=model,finalState=finalState,
                                intermediateState=intermediateState)
    sms = TheorySMS()
    for nodeIndex,node in zip(expSMS.nodeIndices,expSMS.nodes):
        sms.add_node(node,nodeIndex)
    sms.add_edges_from(expSMS.edgeIndices)
    sms.prodXSec = prodXSec
    sms.maxWeight = maxWeight
    sms.setGlobalProperties()

    return sms

def sortExptRes ( exptRes ):
    """ the experimental results may be in different orders.
        sort by AnalysisId+datasetid+TxNames """
    exptRes.sort ( key = lambda x: x["AnalysisID"]+str(x["DataSetID"])+str(x["TxNames"] ) )
    return exptRes

def sortSModelSOutput ( smodelsOutput ):
    smodelsOutput["ExptRes"] = sortExptRes ( smodelsOutput["ExptRes"] )
    return smodelsOutput

def flattenElement(elStr):

    oldStr = elStr[elStr.find('[')+1:elStr.rfind(']')]
    newStr = oldStr.replace('[','').replace(']','')
    newStr = ','.join(sorted([x for x in newStr.split(',') if x.strip()]))
    newElStr = elStr.replace(oldStr,newStr)
    return newElStr


def equalObjs(obj1, obj2, allowedRelDiff, ignore=[], where=None, fname=None,
              fname2=None, checkBothOrders=True):
    """
    Compare two objects.
    The numerical values are compared up to the precision defined by allowedRelDiff.

    :param obj1: First python object to be compared
    :param obj2: Second python object to be compared
    :param allowedRelDiff: Allowed % difference between two numerical values
    :param ignore: List of keys to be ignored
    :param where: keep track of where we are, for easier debugging.
    :param fname: the filename of obj1
    :param fname2: the filename of obj2
    :param checkBothOrders: If True, check if obj1 == obj2 and obj2 == obj1.
    :param version3: If True, tries to take into account differences of output
                     between version 2 and version 3
    :param allowedAbsDiff: If the relative difference is larger than allowedRelDiff,
                           check if the absolute difference is within allowedAbsDiff

    :return: True/False
    """

    if type(fname) == str:
        fname = fname.replace(os.getcwd(), ".")
    if type(obj1) in [float, int] and type(obj2) in [float, int]:
        obj1, obj2 = float(obj1), float(obj2)

    if type(obj1) != type(obj2):
        if (where == 'Mass (GeV)' or where == 'Width (GeV)'):
            if obj1 is None or obj2 is None:
                return True
        logger.warning("Data types differ: (%s,%s) <-> (%s,%s) in ''%s'':%s" % (obj1, type(obj1), obj2, type(obj2), where, fname ))
        return False

    if isinstance(obj1, unum.Unum):
        if obj1 == obj2:
            return True
        abs_diff = abs(obj1-obj2)
        rel_diff = 2.*abs_diff/abs(obj1+obj2)
        # For numbers with units, do not check for absolute difference
        ret = (rel_diff.asNumber() < allowedRelDiff)
        if not ret:
            logger.error( f"values {obj1} and {obj2} differ by {rel_diff}" )
        return ret
    elif isinstance(obj1, float):
        if obj1 == obj2:
            return True
        abs_diff = abs(obj1-obj2)
        rel_diff = 2.*abs_diff/abs(obj1+obj2)
        ret = (rel_diff < allowedRelDiff)
        if not ret:
            logger.error("values %s and %s differ by %s in ''%s'': %s != %s" % (obj1, obj2, rel_diff, where, fname, fname2))
        return ret
    elif isinstance(obj1, str):
        obj1 = obj1.replace(" ","")  # Remove blanks
        obj2 = obj2.replace(" ","")  # Remove blanks
        if obj1 != obj2:
            logger.error("strings ``%s'' and ``%s'' differ in %s:%s" % (obj1, obj2, where, fname))
        return obj1 == obj2
    elif isinstance(obj1, dict):
        for key in obj1:
            if key in ignore:
                continue
            if key not in obj2:
                if where is None:
                    where = "unspecified"
                if fname2 is None:
                    fname2 = "unspecified"
                deffile = f" (default file {fname2})"
                if fname2 == "unspecified":
                    deffile = ""
                logger.warning("Key ``%s'' missing in %s:%s%s" % (key, where, fname, deffile ))
                return False
            if not equalObjs(obj1[key], obj2[key], allowedRelDiff, ignore=ignore,
                             where=key, fname=fname, fname2=fname2):
                logger.warning("Objects differ for %s" % (key))
                return False
    elif isinstance(obj1, list):
        if len(obj1) != len(obj2):
            logger.warning('Lists for %s differ in length:\n   %i (this run)\n and\n   %i (default)' %
                           (where,len(obj1), len(obj2)))
            return False
        for ival, val in enumerate(obj1):
            if not equalObjs(val, obj2[ival], allowedRelDiff, fname=fname, ignore=ignore,
                             fname2=fname2):
                # logger.warning('Lists differ:\n   %s (this run)\n and\n   %s (default)' %\
                #                (str(val),str(obj2[ival])))
                return False
    else:
        return obj1 == obj2

    # Now check for the opposite order of the objects
    if checkBothOrders:
        if not equalObjs(obj2, obj1, allowedRelDiff, ignore, where,
                         fname2, fname, checkBothOrders=False):
            logger.error("Objects %s and %s differ in %s: %s != %s" % (obj1, obj2, where, fname, fname2))
            return False
    return True


def sortXML(xmltree):
    for el in xmltree:
        sortXML(el)
    xmltree[:] = sorted(xmltree, key=lambda el: [el.tag, ElementTree.tostring(el)])


def compareXML(xmldefault, xmlnew, allowedRelDiff, ignore=[]):

    if len(xmldefault) != len(xmlnew):
        logger.warning( f"lengths of document {len(xmldefault)} != {len(xmlnew)}" )
        return False
    for i, el in enumerate(xmldefault):
        newel = xmlnew[i]
        if len(el) != len(newel):
            logger.warning( f"lengths of elements {el.tag} and {newel.tag} differ ({len(el)} != {len(newel)})" )
            return False
        if len(el) == 0:
            if el.tag in ignore:
                continue
            if el.tag != newel.tag:
                logger.warning( f"tags {el.tag} and {newel.tag} differ" )
                return False

            if el.text == newel.text:
                continue

            if type(el.text) == str and "[" not in el.text:
                try:
                    el.text = eval(el.text)
                    newel.text = eval(newel.text)
                except (TypeError, NameError, SyntaxError):
                    el.text = el.text.replace(" ","")
                    newel.text = newel.text.replace(" ","")

            if isinstance(el.text, float) and isinstance(newel.text, float):
                diff = 2.*abs(el.text-newel.text)/abs(el.text+newel.text)
                if diff > allowedRelDiff:
                    logger.warning( f"values {el.text} and {newel.text} differ" )
                    return False
            elif newel.text != el.text:
                logger.warning( f"texts {el.text} and {newel.text} differ" )
                return False
        else:
            if not compareXML(el, newel, allowedRelDiff, ignore):
                return False

    return True


def compareSLHA(slhadefault, slhanew, allowedRelDiff):

    newData = pyslha.read(slhanew, ignorenomass=True, ignorenobr=True, ignoreblocks=["SMODELS_SETTINGS"])
    defaultData = pyslha.read(slhadefault, ignorenomass=True, ignorenobr=True, ignoreblocks=["SMODELS_SETTINGS"])
    defaultBlocks = sorted([defaultData.blocks[b].name for b in defaultData.blocks])
    newBlocks = sorted([newData.blocks[b].name for b in newData.blocks])
    if defaultBlocks != newBlocks:
        logger.error(f'Block structure differs! {defaultBlocks} != {newBlocks}')
        return False

    for b in defaultData.blocks:
        if len(defaultData.blocks[b].entries) != len(newData.blocks[b].entries):
            logger.error( f'Numbers of entries in block {defaultData.blocks[b].name} differ' )
            return False
        keys = defaultData.blocks[b].keys()
        bkeys = newData.blocks[b].keys()
        if keys != bkeys:
            logger.error(f"Keys of blocks {b} differ!")
            return False
        for k in keys:
            ei = defaultData.blocks[b].entries[k]
            ej = newData.blocks[b].entries[k]
            if type(ei) == float and type(ej) == float:
                denom = ei + ej
                if denom == 0.:
                    denom = 1e-6
                de = 2. * abs(ei - ej) / denom
                if de > allowedRelDiff:
                    logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                    return False
            elif ei != ej:
                logger.error(f'Entries in block differ: {ei}!={ej} {type(ei)}')
                return False
    return True


def importModule(filename):
    """ import a module, but giving the filename """
    if sys.version_info[0] == 2:
        import imp
        ## python2, use imp
        with open(filename, 'rb') as fp:  # imports file with dots in name
            output_module = imp.load_module("output", fp, filename,
                    ('.py', 'rb', imp.PY_SOURCE))
        return output_module.smodelsOutput
    ### python3, use importlib
    import importlib
    spec = importlib.util.spec_from_file_location("output", filename)
    output_module = importlib.util.module_from_spec(spec)
    # In case the module contains infs
    output_module.inf = float('inf')
    spec.loader.exec_module(output_module)
    try:
        output = output_module.smodelsOutput
    except AttributeError:
        output = output_module.smodelsOutputDefault

    return output

def runMain(filename, timeout=0, suppressStdout=True, development=False,
            inifile="testParameters_v2.ini", overridedatabase=None):
    """ run SModelS proper
    :param filename: slha file
    :param timeout: timeout for the operation, given in seconds
    :param suppressStdout: if True, then redirect stdout and stderr to /dev/null
    :param development: turn on development mode (e.g. no crash report)
    :param inifile: the config file to be used
    :param overridedatabase: if not None, then use the provided database,
           else use databaseLoader.database
    :returns: printer output
    """
    to = None
    oldlevel = getLogLevel()
    # level = 'debug'
    level = 'info'
    if suppressStdout:
        level = 'fatal'
        to = os.devnull
    database = None
    from smodels.base import runtime
    if overridedatabase is not None:
        database = overridedatabase
    else:
        from databaseLoader import database  # to make sure the db exists
    with redirector.stdout_redirected(to=to):
        out = join(iDir(), "unittests/unitTestOutput")
        setLogLevel(level)
        run(filename, parameterFile=join(iDir(), f"unittests/{inifile}" ),
            outputDir=out, db=database, timeout=timeout,
            development=development)
        setLogLevel(oldlevel)
        sfile = join(iDir(), f"unittests/unitTestOutput/{basename(filename)}.py" )
        return sfile


def compareScanSummary(outA, outB, allowedRelDiff):

    # Check header:
    with open(outA, 'r') as f:
        headerA = f.readlines()[:5]
    with open(outB, 'r') as f:
        headerB = f.readlines()[:5]
    for il,lA in enumerate(headerA):
        if lA != headerB[il]:
            logger.error("Headers differ:\n %s\n and\n %s\n" % (lA, headerB[il]))
            return False

    fA = np.genfromtxt(outA, dtype=None, encoding='utf-8',
                       skip_header=5, names=True)

    fB = np.genfromtxt(outB, dtype=None, encoding='utf-8',
                       skip_header=5, names=True)


    if sorted(fA['filename']) != sorted(fB['filename']):
        logger.error("Filenames differ:\n %s\n and\n %s" % (sorted(fA['filename']), sorted(fB['filename'])))
        return False

    for fname in fA['filename']:
        ptA = fA[fA['filename'] == fname][0]
        ptB = fB[fB['filename'] == fname][0]
        for col in fA.dtype.names:
            if ptA[col] == ptB[col]:
                continue
            elif isinstance(ptA[col], (float, int)):
                diff = 2.*abs(ptA[col]-ptB[col])/abs(ptA[col]+ptB[col])
                if diff > allowedRelDiff:
                    logger.error("values for %s differ by %s in %s" % (col, diff, fname))
                    return False
            else:
                logger.error("values for %s differ in %s" % (col, fname))
                return False
    return True


def compareObjs(obj1, obj2, allowedRelDiff=0.05):

    obj1_keys = sorted([k for k in dir(obj1) if not k[0] == '_'])
    obj2_keys = sorted([k for k in dir(obj2) if not k[0] == '_'])

    if obj1_keys != obj2_keys:
        logger.warning('Object attributes differ:\n %s \n and \n %s' % (obj1_keys, obj2_keys))
        return False

    for key in obj1_keys:
        attr1 = getattr(obj1, key)
        attr2 = getattr(obj2, key)
        if type(attr1) != type(attr2):
            logger.warning('Type of attribute differ:\n %s \n and \n %s' % (type(attr1), type(attr2)))
            return False
        elif attr1 != attr2:
            if isinstance(attr1, (float, int)):
                rel_diff = abs(attr1-attr2)/(attr1+attr2)
                if rel_diff > allowedRelDiff:
                    logger.warning('Attribute %s value differ more than %s:\n %s \n and \n %s'
                                   % (key, allowedRelDiff, attr1, attr2))
                    return False
            else:
                logger.warning('Attribute %s value differ:\n %s \n and \n %s' % (key, attr1, attr2))
                return False
    return True

def getNodesIndices(sms):
    '''
    Convenience function to convert the nodes and indices to strings.
    '''
    nodes_and_indices = list(zip(sms.nodes,sms.nodeIndices))
    nodes_and_indices = [(str(node),inode) for node,inode in nodes_and_indices[:]]
    nodes_and_indices = sorted(nodes_and_indices,key = lambda pt: pt[1])

    return nodes_and_indices

def getEdges(sms):
    '''
    Convenience function to convert the edges to strings.
    '''
    edges = sorted([(str(mom),str(daughter))
             for mom,daughter in sms.edges])

    return edges

class Summary():
    """
    Class to access the output given in the summary.txt
    """

    def __init__(self, filename, allowedRelDiff=0.05):
        self._allowedRelDiff = allowedRelDiff
        self._filename = filename
        self._read(filename)

    def __str__(self):
        fn = self._filename.replace("//", "/")
        final = fn.replace(os.getcwd(), ".")
        if len(final) > 50:
            final = "..."+final[-47:]
        return "Summary(%s)" % final

    def __eq__(self, other):
        return compareObjs(self, other, self._allowedRelDiff)

    def _stripComments(self, lineString):
        return re.sub('#.*', '', lineString)

    def _read(self, filename):
        """
        Read input file and store output
        """
        f = open(filename)
        lines = f.readlines()
        f.close()
        blocks = {}
        blockLines = []
        for line in lines:
            if not line.replace('\n', '').strip():  # Skip empty lines
                continue
            if not line.replace('=', '').replace('\n', '').strip():  # Found separator, reset block
                if blockLines:
                    blocks[blockLines[0]] = blockLines
                blockLines = []
            else:
                blockLines.append(line)

        for block in blocks:
            if 'Input status' in block:
                self._getInputStatus(blocks[block])
            elif '#Analysis' in block and " Theory_Value(fb) " in block:
                self._getResultsOutput(blocks[block])
            elif 'The highest r value' in block:
                self._getResultSummary(blocks[block])
            elif 'Combined Analyses' in block:
                self._getCombinedOutput(blocks[block])
            elif 'Total cross-section for missing' in block:
                self._getMissingSummary(blocks[block])
            elif 'missing topologies with the highest' in block:
                self._getMissedTopos(blocks[block], attrLabel='missingTopos')
            elif 'missing topologies with displaced decays' in block:
                self._getMissedTopos(blocks[block], attrLabel='missingToposDisplaced')
            elif 'missing topologies with prompt decays' in block:
                self._getMissedTopos(blocks[block], attrLabel='missingToposPrompt')
            elif 'topologies outside the grid with the highest' in block:
                self._getMissedTopos(blocks[block], attrLabel='outsideTopos')

    def _getInputStatus(self, inputLines):
        """
        Store input status information
        """

        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            if 'Input status' in line:
                self.inputStatus = eval(line.split(':')[1])
            elif 'Decomposition output status' in line:
                self.decompStatus = eval(line.split(':')[1])

    def _getResultsOutput(self, inputLines):
        """
        Store results information
        """

        results = []
        resultLines = []
        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            if '----' in line and not line.replace('-', '').strip():
                results.append(resultLines)
                resultLines = []
            else:
                resultLines.append(line)
        self.results = [ResultOutput(res, self._allowedRelDiff) for res in results]

    def _getResultSummary(self, inputLines):
        """
        Store results summary information
        """

        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            if 'The highest r value is' in line:
                rval = eval(line.split('=')[1].split(' from ')[0])
                anaID = line.split('=')[1].split(' from ')[1].split()[0].strip()
                self.highestR = rval
                self.highestAna = anaID
            elif 'CMS analysis with highest available r_expected' in line:
                line = line.split(':')[1]
                try:
                    anaID, rexp, robs = line.split(',')
                    rexp = eval(rexp.split('=')[1])
                    robs = eval(robs.split('=')[1])
                    self.highestRCMS = robs
                    self.highestRexpCMS = rexp
                    self.highestAnaCMS = anaID
                except (IndexError, NameError, ValueError):
                    pass
            elif 'ATLAS analysis with highest available r_expected' in line:
                line = line.split(':')[1]
                try:
                    anaID, rexp, robs = line.split(',')
                    rexp = eval(rexp.split('=')[1])
                    robs = eval(robs.split('=')[1])
                    self.highestRATLAS = robs
                    self.highestRexpATLAS = rexp
                    self.highestAnaATLAS = anaID
                except (IndexError, NameError, ValueError):
                    pass

    def _getCombinedOutput(self, inputLines):
        """
        Store combined analyses information
        """

        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            if 'Combined Analyses' in line:
                anaList = sorted(line.split(':')[1].split(','))
                self.combinedIDs = anaList
            elif 'Likelihoods' in line:
                lstr = line.split(':')[1].strip()
                lkeys, lvals = lstr.split('=')
                lkeys = lkeys.split(',')
                lvals = [eval(val) for val in lvals.split(',')]
                for i, label in enumerate(lkeys):
                    setattr(self, 'comb'+label.strip(), lvals[i])
            elif 'combined r-value (expected)' in line:
                rexp = eval(line.split(':')[1])
                self.rexpComb = rexp
            elif 'combined r-value' in line:
                r = eval(line.split(':')[1])
                self.rComb = r

    def _getMissingSummary(self, inputLines):
        """
        Store summary information about missing topologies
        """

        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            if 'Total cross-section for missing topologies with displaced decays' in line:
                self.missedDisplacedTotal = eval(line.split(':')[1])
            elif 'Total cross-section for missing topologies with prompt decays' in line:
                self.missedPromptTotal = eval(line.split(':')[1])
            elif 'Total cross-section for missing topologies' in line:
                self.missedTotal = eval(line.split(':')[1])
            elif 'Total cross-section for topologies outside the grid' in line:
                self.outsideTotal = eval(line.split(':')[1])

    def _getMissedTopos(self, inputLines, attrLabel):
        """
        Store information about missing topologies
        """

        if len(inputLines) < 3:
            setattr(self, attrLabel, None)
            return
        inputLines = inputLines[2:]
        setattr(self, attrLabel, [])
        for line in inputLines:
            line = line.replace('\n', '').strip()
            line = self._stripComments(line)
            if not line:  # Remove comment and empty lines
                continue
            getattr(self, attrLabel).append(MissedTopoOutput(line, self._allowedRelDiff))


class ResultOutput(object):

    def __init__(self, resLines, allowedRelDiff):
        anaID, sqrts, cond, tpValue, expLimit, r, rExp = resLines[0].split()
        self.anaID = anaID.strip()
        self.sqrts = eval(sqrts)
        self.cond = eval(cond)
        self.theoryPred = eval(tpValue)
        self.expLimit = eval(expLimit)
        self.r = eval(r)
        self._allowedRelDiff = allowedRelDiff
        try:
            self.r_expected = eval(rExp)
        except NameError:
            self.r_expected = rExp
        for line in resLines:
            if 'Signal Region' in line:
                self.signalRegion = line.split(':')[1].strip()
            elif 'Txnames' in line:
                self.txNames = line.split(':')[1].strip()
            elif 'Likelihoods' in line:
                lstr = line.split(':')[1].strip()
                lkeys, lvals = lstr.split('=')
                lkeys = lkeys.split(',')
                lvals = [eval(val) for val in lvals.split(',')]
                for i, label in enumerate(lkeys):
                    setattr(self, label.strip(), lvals[i])

    def __eq__(self, other):
        return compareObjs(self, other, self._allowedRelDiff)


class MissedTopoOutput(object):

    def __init__(self, line, allowedRelDiff):
        sqrts, weight = line.split()
        self.sqrts = eval(sqrts)
        self.weight = eval(weight)
        self._allowedRelDiff = allowedRelDiff

    def __eq__(self, other):
        return compareObjs(self, other, self._allowedRelDiff)
