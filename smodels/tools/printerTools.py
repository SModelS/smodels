
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
import numpy as np
import itertools
from collections import OrderedDict


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
        header = f"#In order to build the summary, one of the following types of printer must be available:\n {str(printerTypes)} \n"
        with open(outputFile, 'w') as f:
            f.write(header)
        return

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
                                 'r_max': np.nan,
                                 'r_exp': np.nan,
                                 'MostSensitive(ATLAS)': 'N/A',
                                 'r(ATLAS)': np.nan,
                                 'r_exp(ATLAS)': np.nan,
                                 'MostSensitive(CMS)': 'N/A',
                                 'r(CMS)': np.nan,
                                 'r_exp(CMS)': np.nan,
                                 'r(combined)' : np.nan,
                                 'r_exp(combined)' : np.nan,
                                 'CombinedAnalyses' : 'N/A'
                                   })

        if 'python' in output:
            sDict = getScanSummaryFrom(output['python'], 'python')
        elif 'slha' in output:
            sDict = getScanSummaryFrom(output['slha'], 'slha')
        elif 'summary' in output:
            sDict = getScanSummaryFrom(output['summary'], 'summary')
        else:
            sDict = {}

        for key in summaryDict:
            if key in sDict:
                summaryDict[key] = sDict[key]

        summaryList.append(summaryDict)

    # If there are no combined results, remove its (dummy) entries
    anaList = set([])
    if all(np.isnan(summary['r(combined)']) for summary in summaryList):
        for summary in summaryList:
            summary.pop('r(combined)')
            summary.pop('r_exp(combined)')

    else:
        # Get maximum list of combined analyses
        # and remove info from rows
        for summary in summaryList:
            anas = summary.pop('CombinedAnalyses')
            if anas == 'N/A':
                continue
            anas = anas.replace(' ','').split(',')
            anaList.update(anas)
        anaList = sorted(anaList)

    # Header:
    header = f"#Global results summary ({len(outputDict)} files)\n"
    header += "#The most constraining analysis corresponds to the one with largest observed r.\n"
    header += "#The most sensitive (ATLAS/CMS) analysis corresponds to the one with largest expected r from those analyses for which this information is available.\n"
    if anaList:
        header += f"#Analyses used for combination = {','.join(anaList)}.\n"
        header += "#r(combined) = nan means no analysis from the above combination set produced results.\n"


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

def formatNestedDict(outputDict,ident=0,maxLength=50):
    """
    Convert a nested dictionary to a string
    with identation.

    :param outputDict: Dictionary to be printed
    :param ident: Current identation
    :param maxLength: Maximum length allowed without identation

    return: String with formatted output
    """

    if len(str(outputDict)) <= maxLength:
        return str(outputDict)

    output = '{\n'
    for ik,(key,val) in enumerate(outputDict.items()):
        if isinstance(val,dict):
            valStr = formatNestedDict(val,ident=ident+4,maxLength=maxLength)
        elif isinstance(val,list):
            valStr = formatNestedList(val,ident=ident+4,maxLength=maxLength)
        elif isinstance(val,str):
            valStr = val.replace("'","")
            valStr = "'"+valStr+"'"
        else:
            valStr = str(val)

        if isinstance(key,str):
            keyStr = "'"+key+"'"
        else:
            keyStr = str(key)
        if ik == 0:
            output += ' '*ident+f"{keyStr} : {valStr}"
        else:
            output += f"{' ' * ident + keyStr} : {valStr}"

        if ik < len(outputDict)-1:
            output += ",\n"
        else:
            output += "\n"

    output += ' '*(ident-4)+'}'
    return output

def formatNestedList(outputList,ident=0,maxLength=50):
    """
    Convert a nested list to a string
    with identation.

    :param outputList: List to be formatted
    :param ident: Current identation
    :param maxLength: Maximum length allowed without identation

    return: String with formatted output
    """

    if len(str(outputList)) <= maxLength:
        return str(outputList)

    output = '[\n'
    for iv,val in enumerate(outputList):
        if isinstance(val,dict):
            valStr = formatNestedDict(val,ident=ident+4,maxLength=maxLength)
        elif isinstance(val,list):
            valStr = formatNestedList(val,ident=ident+4,maxLength=maxLength)
        elif isinstance(val,str):
            valStr = val.replace("'","")
            valStr = "'"+valStr+"'"
        else:
            valStr = str(val)

        if iv == 0 :
            output += ' '*ident+f'{valStr}'
        else:
            output += f"{' ' * ident + valStr}"
        if iv < len(outputList)-1:
            output += ',\n'
        else:
            output += '\n'
    output += ' '*(ident-4)+']'
    return output

def getScanSummaryFrom(output, ptype):
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
        rvals, rexp, anaIDs, r_comb, rexp_comb, anaID_comb = info

    # Convert Nones to NaN:
    rvals = np.array(rvals,dtype=float) # Convert None to nan
    rexp = np.array(rexp,dtype=float) # Convert None to nan
    r_comb,rexp_comb = np.array([r_comb,rexp_comb],dtype=float)


    # Sort results by r_obs
    # (replace nan by -1 for sorting only)
    rvals_c = np.nan_to_num(rvals,nan=-1.0)
    imax = rvals_c.argsort()[::-1][0]
    summaryDict['r_max'] = rvals[imax]
    summaryDict['r_exp'] = rexp[imax]
    summaryDict['MostConstrainingAnalysis'] = anaIDs[imax]

    # Now keep only results with r_exp != nan:
    rexp_good = rexp[~np.isnan(rexp)]
    rvals_good = rvals[~np.isnan(rexp)]
    anaIDs_good = anaIDs[~np.isnan(rexp)]
    asort = rexp_good.argsort()[::-1]
    rvals_good = rvals_good[asort]
    anaIDs_good = anaIDs_good[asort]
    rexp_good = rexp_good[asort]
    iATLAS = [i for i,ana in enumerate(anaIDs_good) if 'ATLAS' in ana]
    iCMS = [i for i,ana in enumerate(anaIDs_good) if 'CMS' in ana]
    if len(iATLAS) > 0:
        iATLAS = iATLAS[0]
        summaryDict['r(ATLAS)'] = rvals_good[iATLAS]
        summaryDict['r_exp(ATLAS)'] = rexp_good[iATLAS]
        summaryDict['MostSensitive(ATLAS)'] = anaIDs_good[iATLAS]
    if len(iCMS) > 0:
        iCMS = iCMS[0]
        summaryDict['r(CMS)'] = rvals_good[iCMS]
        summaryDict['r_exp(CMS)'] = rexp_good[iCMS]
        summaryDict['MostSensitive(CMS)'] = anaIDs_good[iCMS]

    summaryDict['r(combined)'] = r_comb
    summaryDict['r_exp(combined)'] = rexp_comb
    summaryDict['CombinedAnalyses'] = anaID_comb

    return summaryDict

def getInfoFromPython(output):
    """
    Retrieves information from the python output

    :param output: output (dictionary)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
             If there are results for combined analyses, returns the largest r-value and the
             corresponding r-expected from the combination.

    """

    if 'ExptRes' not in output or not output['ExptRes']:
        return None
    rvals = np.array([res['r'] for res in output['ExptRes']])
    rexp = np.array([res['r_expected'] if res['r_expected']
                   else np.nan for res in output['ExptRes']])
    anaIDs = np.array([res['AnalysisID'] for res in output['ExptRes']])

    r_comb = np.nan
    rexp_comb = np.nan
    anaID_comb = 'N/A'

    if 'CombinedRes' in output:
        for res in output['CombinedRes']:
            if np.isnan(r_comb) or r_comb < res['r']:
                r_comb = res['r']
                rexp_comb = res['r_expected']
                anaID_comb = res['AnalysisID']

    return rvals, rexp, anaIDs, r_comb, rexp_comb, anaID_comb

def getInfoFromSLHA(output):
    """
    Retrieves information from the SLHA output

    :param output: output (string)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
             If there are results for combined analyses, returns the largest r-value and the
             corresponding r-expected from the combination.
    """

    import pyslha
    results = pyslha.readSLHA(output, ignorenomass=True, ignorenobr=True)
    bname = None
    bcombName = None
    for b in results.blocks.values():
        if b.name.lower() == 'SModelS_Exclusion'.lower():
            bname = b.name
        if b.name.lower() == 'SModelS_CombinedAnas'.lower():
            bcombName = b.name

    if bname is None or len(results.blocks[bname]) <= 1:
        return None
    else:
        # Group results by block index:
        groups = itertools.groupby(results.blocks[bname].items(),
                                  key = lambda k: k[0][0])
        resDict = {i : dict(block) for i,block in groups if  i != 0}
        # Get r values:
        rvals = np.array([resDict[i][(i,1)] for i in resDict])
        rexp = np.array([resDict[i][(i,2)] if resDict[i][(i,2)] != 'N/A' else np.nan
                        for i in resDict])
        anaIDs = np.array([resDict[i][(i,4)] for i in resDict])

    if bcombName is None or len(results.blocks[bcombName]) < 1:
        r_comb = np.nan
        rexp_comb = np.nan
        anaID_comb = 'N/A'
    else:
        # Group combined results by block index:
        groups = itertools.groupby(results.blocks[bcombName].items(),
                                  key = lambda k: k[0][0])
        resDict = {i : dict(block) for i,block in groups if  i != 0}
        rvals_comb = np.array([resDict[i][(i,1)] for i in resDict  if  i != 0])
        rexp_comb = np.array([resDict[i][(i,2)] if resDict[i][(i,2)] != 'N/A' else np.nan
                        for i in resDict if i != 0])
        anaID_comb = np.array([resDict[i][(i,6)] for i in resDict if i != 0])

        r_comb = max(rvals_comb[~np.isnan(rvals_comb)])
        rexp_comb = rexp_comb[np.argmax(rvals_comb)]
        anaID_comb = anaID_comb[np.argmax(rvals_comb)]

    return rvals, rexp, anaIDs, r_comb, rexp_comb, anaID_comb

def getInfoFromSummary(output):
    """
    Retrieves information from the summary output

    :param output: output (string)

    :return: list of r-values,r-expected and analysis IDs. None if no results are found.
             If there are results for combined analyses, returns the largest r-value and the
             corresponding r-expected from the combination.

    """

    lines = output.splitlines()
    rvals = []
    rexp = []
    anaIDs = []
    r_comb = np.nan
    rexp_comb = np.nan
    anaID_comb = 'N/A'

    for line in lines:
        if 'The highest r value is' in line:
            rmax = line.split('=')[1].strip()
            ff = np.where([((not x.isdigit()) and (x not in ['.', '+', '-']))
                         for x in rmax])[0][0]  # Get when the value ends
            rmax = eval(rmax[:ff])
            anaMax = line.split('from')[1].split()[0].replace(',', '')
            rexpMax = np.nan
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
            rexpAna = np.nan
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
        elif 'combined r-value:' in line:
            r_comb = float(line.replace('\n','').split(':')[1])
        elif 'combined r-value (expected):' in line:
            rexp_comb = float(line.replace('\n','').split(':')[1])
        elif 'Combined Analyses:' in line:
            anaID_comb = line.replace('\n','').split(':')[1].replace(' ','')

    if not rvals:
        return None
    rvals = np.array(rvals)
    rexp = np.array(rexp)
    anaIDs = np.array(anaIDs)

    return rvals, rexp, anaIDs, r_comb, rexp_comb, anaID_comb

