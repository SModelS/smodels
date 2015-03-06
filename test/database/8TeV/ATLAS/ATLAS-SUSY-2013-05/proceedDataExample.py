#!/usr/bin/env python

"""
.. module:: convert
   :synopsis: uesed to create info.txt,sms.py,sms.root and newSms.py.

.. moduleauthor:: Michael Traub <michael.traub@gmx.at>

"""   

import sys
sys.path.append('../../../../smodels-utils/smodels_utils/dataPreparation')
from origPlotObjects import OrigPlot
import ROOT



def getEquationSets(txName):
    """read the equationSets stored in the axses field
    :param:txName: txName as string
    :return: list with equationSets
    """
    f = open(txName + '.txt')
    lines = f.readlines()
    f.close()
    for line in lines:
        if not line[:4] == 'axes':
            continue
        line = line.replace('\n','')
        eqSets = line.split(': ')[1]
        return eqSets.split(';')
        
def getUpperLimits(txName):

    """read the upperlimits stored in the upperlimitsField field
    of info txt
    :param:txName: txName as string
    :return: upperlimits
    """
    f = open(txName + '.txt')
    lines = f.readlines()
    f.close()
    for line in lines:
        #print line
        if 'upperLimits' in line:
            line = line.replace('\n','')
            upperLimits = line.split(': ')[1]
            upperLimits = upperLimits.replace('*GeV','')
            upperLimits = upperLimits.replace('*pb','')
            return eval(upperLimits)
        
        
def getHisto(points, txName, eqSet):
    """creates a ROOT.TH2F object
    :param:point: list with xValues, yValues and limits
    :param:txName: txName as string
    :param:eqSet: equationSet as String
    """
    xMin = min(points,key = lambda p: p[0])[0]
    xMax = max(points,key = lambda p: p[0])[0]
    yMin = min(points,key = lambda p: p[1])[1]
    yMax = max(points,key = lambda p: p[1])[1]
  
    histo = ROOT.TH2F(txName +'_'+ eqSet, txName +'_'+ eqSet ,\
    int(xMax-xMin),xMin,xMax, int(yMax-yMin), yMin, yMax)
    for x,y,limit in points:
        histo.Fill(x,y,limit)
    histo.SetDirectory(0)
    return histo
    
def main():
    # chose TxName:
    txName = 'T6bbWWoff'
    # get all available equationSets/massPlanes/plots
    # from info.txt:
    equationSets = getEquationSets(txName)
    # chose equationSet:
    equationSet = equationSets[0]
    # init OrigPlot object with chosen equationSet:
    origPlot = OrigPlot.fromString(equationSet)

    # get a list with xValues, yValues and upperLimits 
    # for the chosen equationSet from sms.py
    points = []
    upperLimits = getUpperLimits(txName)
    print upperLimits
    for entry in upperLimits:
        massArray, limit = entry[0][0], entry[1]
        xy = origPlot.getXYValues(massArray)
        print '%s --> %s: %s' %(massArray, xy, limit )
        if not xy: continue
        points.append([xy[0], xy[1], limit])
        
    # open root.txt
    smsRoot = ROOT.TFile('sms.root')
    # chreate a Histogram
    histo = getHisto(points, txName, str(origPlot))
    #draw everything
    c = ROOT.TCanvas('myCanvas')
    histo.Draw('textsame')
    smsRoot.Get('%s/exclusion_%s' %(txName, str(origPlot))).Draw('same')
    smsRoot.Get('%s/exclusionM1_%s' %(txName, str(origPlot))).Draw('same')
    smsRoot.Get('%s/exclusionP1_%s' %(txName, str(origPlot))).Draw('same')
    smsRoot.Get('%s/expectedExclusion_%s' %(txName, str(origPlot))).Draw('same')
    smsRoot.Get('%s/expectedExclusionM1_%s' %(txName, str(origPlot))).Draw('same')
    smsRoot.Get('%s/expectedExclusionP1_%s' %(txName, str(origPlot))).Draw('same')
    c.Update()
    r = raw_input()
    
if __name__ == "__main__":
    main()

      
   
    

