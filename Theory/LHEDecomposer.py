#!/usr/bin/env python

"""
.. module:: LHEDecomposer
    :synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

def decompose(lhefile,W=None,nevts=None,doCompress=False,doInvisible=False,minmassgap=None):
  """ Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param W: Dictionary with cross-sections, the pdgids of the mothers being the \
      keys. If W=None, use cross-section information from LHE file
    :param nevts: (maximum) number of generated events. If None, all events from \
      file are processed.
    :returns: a TopologyList object 
    """
  import LHEReader, TopologyBuilder, topology, types
  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=topology.TopologyList ( )
  for Event in reader:
    ## Get mother PDGs:
    momPDG = tuple(Event.getMom())
    # Get event weight list:
    weight = {}
    for k in W.keys(): 
      if W[k].has_key(momPDG):
        weight[k]=W[k][momPDG]
      else:
        print "[LHEDecomposer] Error getting weight for",k,momPDG
        # weight[k]=0.
        # return None

    # Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, doCompress, doInvisible, minmassgap)
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
