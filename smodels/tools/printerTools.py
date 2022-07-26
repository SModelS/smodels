
"""
.. module:: printer
   :synopsis: Facility used to print decomposition, theorypredictions, missing topologies et al
      in various forms

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>
.. moduleauthor:: Ursula Laa <ursula.laa@lpsc.in2p3.fr>
.. moduleauthor:: Suchita Kulkanri <suchita.kulkarni@gmail.com>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from __future__ import print_function
import sys
import os
import copy
from smodels.theory.topologyDict import TopologyDict
from smodels.theory.theoryPrediction import TheoryPredictionList
from smodels.tools.theoryPredictionsCombiner import TheoryPredictionsCombiner
from smodels.experiment.databaseObj import ExpResultList
from smodels.tools.ioObjects import OutputStatus
from smodels.tools.coverage import Uncovered
from smodels.tools.physicsUnits import GeV, fb, TeV
from smodels.tools.smodelsLogging import logger
import numpy as np
from collections import OrderedDict
from xml.dom import minidom
from xml.etree import ElementTree
import unum
import time


def printScanSummary(outputDict, outputFile):
    """
    Method for creating a simple summary of the results when running SModelS
    over multiple files.

    :param outputDict: A dictionary with filenames as keys and the master printer flush dictionary as values.
    :param outputFile: Path to the summary file to be written.
    """

    # Check available types of printer:
    printerTypes = ['slha', 'python', 'summary']
    # All outputs should have the same format
    out = list(outputDict.values())[0]
    if all([(ptype not in out) for ptype in printerTypes]):
        header = "#In order to build the summary, one of the following types of printer must be available:\n %s \n" % str(printerTypes)
        with open(outputFile, 'w') as f:
            f.write(header)
        return

    # Header:
    header = "#Global results summary (%i files)\n" % len(outputDict)
    header += "#The most constraining analysis corresponds to the one with largest observed r.\n"
    header += "#The most senstive (ATLAS/CMS) analysis corresponds to the one with largest expected r from those analyses for which this information is available.\n"

    # Get summary information:
    summaryList = []
    fnames = list(outputDict.keys())
    fnames.sort()  # we want a canonical order

    for fname in fnames:
        output = outputDict[fname]
        if output is None:
            continue
        # default values (in case of empty results):
        summaryDict = OrderedDict({'filename': fname,
                                 'MostConstrainingAnalysis': 'N/A',
                                 'r_max': -1,
                                 'r_exp': -1,
                                 'MostSensitive(ATLAS)': 'N/A',
                                 'r(ATLAS)': -1,
                                 'r_exp(ATLAS)': -1,
                                 'MostSensitive(CMS)': 'N/A',
                                 'r(CMS)': -1,
                                 'r_exp(CMS)': -1
                                   })

        if 'python' in output:
            sDict = getSummaryFrom(output['python'], 'python')
        elif 'slha' in output:
            sDict = getSummaryFrom(output['slha'], 'slha')
        elif 'summary' in output:
            sDict = getSummaryFrom(output['summary'], 'summary')
        else:
            sDict = {}

        for key in summaryDict:
            if key in sDict:
                summaryDict[key] = sDict[key]

        summaryList.append(summaryDict)

    # Get column labels and widths:
    labels = list(summaryList[0].keys())
    cwidths = []
    fstr = '%s'  # format for strings
    ffloat = '%1.3g'  # format for floats
    for label in labels:
        maxlength = max([len(ffloat % entry[label]) if isinstance(entry[label], (float, int))
                       else len(fstr % entry[label]) for entry in summaryList])
        maxlength = max(maxlength, len(label))
        cwidths.append(maxlength)

    columns = '#'
    columns += '  '.join([label.ljust(cwidths[i])
                          for i, label in enumerate(labels)])
    # Remove one blank space to make labels match values
    columns = columns.replace(' ', '', 1)
    columns += '\n'

    with open(outputFile, 'w') as f:
        f.write(header)
        f.write(columns)

        for entry in summaryList:
            row = '  '.join([(ffloat % entry[label]).ljust(cwidths[j])
                           if isinstance(entry[label], (float, int))
                           else (fstr % entry[label]).ljust(cwidths[j])
                           for j, label in enumerate(labels)])
            f.write(row+'\n')
    return


def getSummaryFrom(output, ptype):
    """
    Retrieves information about the output according to the printer type (slha,python or summary)

    :param output: output (dictionary for ptype=python or string for ptype=slha/summary)
    :param ptype: Printer type (slha, python or summary)

    :return: Dictionary with the output information
    """

    summaryDict = {}
    if ptype == 'python':
        info = getInfoFromPython(output)
    elif ptype == 'slha':
        info = getInfoFromSLHA(output)
    elif ptype == 'summary':
        info = getInfoFromSummary(output)
    else:
        return summaryDict

    if info is None:
        return summaryDict
    else:
        rvals, rexp, anaIDs = info

    # Sort results by r_obs:
    rvalswo = copy.deepcopy(rvals)
    rvalswo[rvalswo is None] = -1
    asort = rvalswo.argsort()[::-1]
    rvals = rvals[asort]
    anaIDs = anaIDs[asort]
    rexp = rexp[asort]
    summaryDict['r_max'] = rvals[0]
    summaryDict['r_exp'] = rexp[0]
    summaryDict['MostConstrainingAnalysis'] = anaIDs[0]

    # Sort results by r_obs:
    rvalswo = copy.deepcopy(rexp)
    rvalswo[rvalswo is None] = -1
    # Sort results by r_exp:
    asort = rvalswo.argsort()[::-1]
    rvals = rvals[asort]
    anaIDs = anaIDs[asort]
    rexp = rexp[asort]
    iATLAS, iCMS = -1, -1
    for i, anaID in enumerate(anaIDs):
        if rexp[i] < 0:
            continue
        if 'ATLAS' in anaID and iATLAS < 0:
            iATLAS = i
        elif 'CMS' in anaID and iCMS < 0:
            iCMS = i

    if iATLAS >= 0:
        summaryDict['r(ATLAS)'] = rvals[iATLAS]
        summaryDict['r_exp(ATLAS)'] = rexp[iATLAS]
        summaryDict['MostSensitive(ATLAS)'] = anaIDs[iATLAS]

    if iCMS >= 0:
        summaryDict['r(CMS)'] = rvals[iCMS]
        summaryDict['r_exp(CMS)'] = rexp[iCMS]
        summaryDict['MostSensitive(CMS)'] = anaIDs[iCMS]

    return summaryDict


def getInfoFromPython(output):
    """
    Retrieves information from the python output

    :param output: output (dictionary)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
    """

    if 'ExptRes' not in output or not output['ExptRes']:
        return None
    rvals = np.array([res['r'] for res in output['ExptRes']])
    rexp = np.array([res['r_expected'] if res['r_expected']
                   else -1 for res in output['ExptRes']])
    anaIDs = np.array([res['AnalysisID'] for res in output['ExptRes']])

    return rvals, rexp, anaIDs


def getInfoFromSLHA(output):
    """
    Retrieves information from the SLHA output

    :param output: output (string)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
    """

    import pyslha
    results = pyslha.readSLHA(output, ignorenomass=True, ignorenobr=True)
    bname = None
    for b in results.blocks.values():
        if b.name.lower() == 'SModelS_Exclusion'.lower():
            bname = b.name
    if bname is None or len(results.blocks[bname]) <= 1:
        return None

    # Get indices for results:
    resI = list(set([k[0] for k in results.blocks[bname].keys() if k[0] > 0]))
    rvals = np.array([results.blocks[bname][(i, 1)] for i in resI])
    rexp = np.array([results.blocks[bname][(i, 2)]
                   if results.blocks[bname][(i, 2)] != 'N/A' else -1 for i in resI])
    anaIDs = np.array([results.blocks[bname][(i, 4)] for i in resI])

    return rvals, rexp, anaIDs


def getInfoFromSummary(output):
    """
    Retrieves information from the summary output

    :param output: output (string)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
    """

    lines = output.splitlines()
    rvals = []
    rexp = []
    anaIDs = []
    for line in lines:
        if 'The highest r value is' in line:
            rmax = line.split('=')[1].strip()
            ff = np.where([((not x.isdigit()) and (x not in ['.', '+', '-']))
                         for x in rmax])[0][0]  # Get when the value ends
            rmax = eval(rmax[:ff])
            anaMax = line.split('from')[1].split()[0].replace(',', '')
            rexpMax = -1
            if 'r_expected' in line and "r_expected not available" not in line:
                rexpMax = line.split('r_expected')[-1]
                rexpMax = rexpMax.split('=')[1]
                ff = np.where([((not x.isdigit()) and (x not in ['.', '+', '-']))
                             for x in rexpMax])[0][0]  # Get when the value ends
                rexpMax = eval(rexpMax[:ff])
            rvals.append(rmax)
            anaIDs.append(anaMax)
            rexp.append(rexpMax)
        elif 'analysis with highest available r_expected' in line:
            # the space is required to have at least one non-digit character after the value
            rAna = line.split('=')[-1] + ' '
            ff = np.where([((not x.isdigit()) and (x not in ['.', '+', '-']))
                         for x in rAna])[0][0]  # Get when the value ends
            rAna = eval(rAna[:ff])
            rexpAna = -1
            if 'r_expected' in line:
                rexpAna = line.split('r_expected')[-1]
                rexpAna = rexpAna.split('=')[1]
                ff = np.where([((not x.isdigit()) and (x not in ['.', '+', '-']))
                             for x in rexpAna])[0][0]  # Get when the value ends
                rexpAna = eval(rexpAna[:ff])
            if 'CMS' in line:
                anaID = 'CMS-'+line.split('CMS-')[1].split(' ')[0].split(',')[0]
            else:
                anaID = 'ATLAS-' + \
                        line.split('ATLAS-')[1].split(' ')[0].split(',')[0]
            anaID = anaID.split()[0].strip().replace(',', '')
            rvals.append(rAna)
            anaIDs.append(anaID)
            rexp.append(rexpAna)

    if not rvals:
        return None
    rvals = np.array(rvals)
    rexp = np.array(rexp)
    anaIDs = np.array(anaIDs)

    return rvals, rexp, anaIDs
