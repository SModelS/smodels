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
import numpy as np
import redirector
from smodels.tools.runSModelS import run
from os.path import join, basename
from smodels.installation import installDirectory as iDir
from smodels.tools.smodelsLogging import logger, setLogLevel, getLogLevel

def sortExptRes ( exptRes ):
    """ the experimental results may be in different orders. 
        sort by AnalysisId+datasetid+TxNames """
    exptRes.sort ( key = lambda x: x["AnalysisID"]+str(x["DataSetID"])+str(x["TxNames"] ) )
    return exptRes

def sortSModelSOutput ( smodelsOutput ):
    smodelsOutput["ExptRes"] = sortExptRes ( smodelsOutput["ExptRes"] )
    return smodelsOutput

def equalObjs(obj1, obj2, allowedDiff, ignore=[], where=None, fname=None,
              fname2=None, checkBothOrders=True):
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
    :param checkBothOrders: If True, check if obj1 == obj2 and obj2 == obj1.
    :return: True/False
    """
    if type(fname) == str:
        fname = fname.replace(os.getcwd(), ".")
    if type(obj1) in [float, int] and type(obj2) in [float, int]:
        obj1, obj2 = float(obj1), float(obj2)

    if type(obj1) != type(obj2):
        logger.warning("Data types differ: (%s,%s) <-> (%s,%s) in ''%s'':%s" % (obj1, type(obj1), obj2, type(obj2), where, fname ))
        return False

    if isinstance(obj1, unum.Unum):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        return diff.asNumber() < allowedDiff
    elif isinstance(obj1, float):
        if obj1 == obj2:
            return True
        diff = 2.*abs(obj1-obj2)/abs(obj1+obj2)
        if diff > allowedDiff:
            logger.error("values %s and %s differ by %s in ''%s'': %s != %s" % (obj1, obj2, diff, where, fname, fname2))
        return diff < allowedDiff
    elif isinstance(obj1, str):
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
            if not equalObjs(obj1[key], obj2[key], allowedDiff, ignore=ignore, where=key, fname=fname, fname2=fname2):
                return False
    elif isinstance(obj1, list):
        if len(obj1) != len(obj2):
            logger.warning('Lists differ in length:\n   %i (this run)\n and\n   %i (default)' %
                           (len(obj1), len(obj2)))
            return False
        for ival, val in enumerate(obj1):
            if not equalObjs(val, obj2[ival], allowedDiff, fname=fname, ignore=ignore,
                             fname2=fname2):
                # logger.warning('Lists differ:\n   %s (this run)\n and\n   %s (default)' %\
                #                (str(val),str(obj2[ival])))
                return False
    else:
        return obj1 == obj2

    # Now check for the opposite order of the objects
    if checkBothOrders:
        if not equalObjs(obj2, obj1, allowedDiff, ignore, where,
                         fname2, fname, checkBothOrders=False):
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
    spec.loader.exec_module(output_module)
    return output_module.smodelsOutput


def runMain(filename, timeout=0, suppressStdout=True, development=False,
            inifile="testParameters.ini", overridedatabase=None):
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
    if overridedatabase is not None:
        database = overridedatabase
    else:
        from databaseLoader import database  # to make sure the db exists
    with redirector.stdout_redirected(to=to):
        out = join(iDir(), "test/unitTestOutput")
        setLogLevel(level)
        run(filename, parameterFile=join(iDir(), "test/%s" % inifile),
            outputDir=out, db=database, timeout=timeout,
            development=development)
        setLogLevel(oldlevel)
        sfile = join(iDir(), "test/unitTestOutput/%s.py" % basename(filename))
        return sfile


def compareScanSummary(outA, outB, allowedDiff):

    fA = np.genfromtxt(outA, dtype=None, encoding='utf-8',
                       skip_header=3, names=True)

    fB = np.genfromtxt(outB, dtype=None, encoding='utf-8',
                       skip_header=3, names=True)

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
                if diff > allowedDiff:
                    logger.error("values for %s differ by %s in %s" % (col, diff, fname))
                    return False
            else:
                logger.error("values for %s differ in %s" % (col, fname))
                return False
    return True


def compareObjs(obj1, obj2, allowedDiff=0.05):

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
                if rel_diff > allowedDiff:
                    logger.warning('Attribute %s value differ more than %s:\n %s \n and \n %s'
                                   % (key, allowedDiff, attr1, attr2))
                    return False
            else:
                logger.warning('Attribute %s value differ:\n %s \n and \n %s' % (key, attr1, attr2))
                return False
    return True


class Summary():
    """
    Class to access the output given in the summary.txt
    """

    def __init__(self, filename, allowedDiff=0.05):
        self._allowedDiff = allowedDiff
        self._filename = filename
        self._read(filename)

    def __str__(self):
        fn = self.filename.replace("//", "/")
        final = fn.replace(os.getcwd(), ".")
        if len(final) > 50:
            final = "..."+final[-47:]
        return "Summary(%s)" % final

    def __eq__(self, other):
        return compareObjs(self, other, self._allowedDiff)

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
        self.results = [ResultOutput(res, self._allowedDiff) for res in results]

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
            getattr(self, attrLabel).append(MissedTopoOutput(line, self._allowedDiff))


class ResultOutput(object):

    def __init__(self, resLines, allowedDiff):
        anaID, sqrts, cond, tpValue, expLimit, r, rExp = resLines[0].split()
        self.anaID = anaID.strip()
        self.sqrts = eval(sqrts)
        self.cond = eval(cond)
        self.theoryPred = eval(tpValue)
        self.expLimit = eval(expLimit)
        self.r = eval(r)
        self._allowedDiff = allowedDiff
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
        return compareObjs(self, other, self._allowedDiff)


class MissedTopoOutput(object):

    def __init__(self, line, allowedDiff):
        sqrts, weight = line.split()
        self.sqrts = eval(sqrts)
        self.weight = eval(weight)
        self._allowedDiff = allowedDiff

    def __eq__(self, other):
        return compareObjs(self, other, self._allowedDiff)
