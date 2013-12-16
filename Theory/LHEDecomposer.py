#!/usr/bin/env python

"""
.. module:: LHEDecomposer
    :synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import LHEReader, topology, CrossSection

def decompose(lhefile,inputXsecs=None,nevts=None,doCompress=False,doInvisible=False,minmassgap=None):
  """ Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param inputXsecs: Dictionary with cross-section objects for the pdgids of the mothers appearing in the LHE file.\
      If None, use information from file
    :param nevts: (maximum) number of events used in the decomposition. If None, all events from \
      file are processed.
    :returns: a TopologyList object 
    """
  
    reader = LHEReader.LHEReader(lhefile,nevts)
    SMSTopList=topology.TopologyList ( )
#get cross-section from file
    if not inputXsecs:  XSectionList = CrossSection.getXsecFromLHEFile(lhefile)
    else: XSeciontList = inputXsecs 

  
  
    for Event in reader:    
        momPDG = tuple(Event.getMom())  # Get mother PDGs
    
    
    
    for k in W.keys():              # Get event weight
      if W[k].has_key(momPDG): weight[k]=W[k][momPDG]
      else:
        print "[LHEDecomposer] Error getting weight for",k,momPDG
        return False

    # Get event element
    newElement = element.fromEvent(Event,weight)
    #Do compression:
    if DoCompress or DoInvisible: compElements = newElement.compressElement(DoCompress,DoInvisible,minmassgap)
    allElements = [newElement] + compElements
    for el in allElements:
        Top = topology.Topology(el)            
        SMSTopList.addList([Top])                       

  return SMSTopList
