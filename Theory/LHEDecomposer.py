#!/usr/bin/env python

"""
.. module:: LHEDecomposer
    :synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

def decompose(lhefile,W=None,nevts=None,doCompress=False,doInvisible=False,minmassgap=None):
  """ Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param W: Dictionary with event weights, the pdgids of the mothers being the \
      keys. If W = None, get a constant event weight from lhefile itself
    :param nevts: (maximum) number of generated events. If None, all events from \
      file are processed.
    :returns: a TopologyList object 
    """
  import LHEReader, TopologyBuilder, SMSDataObjects, types, CrossSection
  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=SMSDataObjects.TopologyList ( )
  for Event in reader:
    ## Get mother PDGs:
    momPDG = tuple(Event.getMom())  
    # Get event weight list:
    weight = {}
    if not W is None:
      for k in W.keys(): 
        if W[k].has_key(momPDG): weight[k]=W[k][momPDG]
        else: print "[LHEDecomposer] Error getting weight for",k,momPDG
    elif CrossSection.XSectionInfo:
      for xsec in CrossSection.XSectionInfo.xsecs:
        weight[xsec.label] = reader.metainfo["totalxsec"]/reader.metainfo["nevents"]
    else:
        print "[LHEDecomposer.decompose] Cross-Section information not found\
         (either weight dictionary or CrossSection.XSectionInfo) must be defined"    
    
    # Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, doCompress, doInvisible, minmassgap)
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
