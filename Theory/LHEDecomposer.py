#!/usr/bin/python

"""
.. module:: LHEDecomposer
    :synopsis: smodels-decomposes LHE events, creating TopologyLists 

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

def decompose(lhefile,W,nevts=None,doCompress=False,doInvisible=False,minmassgap=None):
  """ Do LHE-based decomposition. 

    :param lhefile: LHE file with e.g. pythia events
    :param W: Dictionary with event weights, the pdgids of the mothers being the \
      keys. Alternatively, a single float can also be supplied, in that case 
      the float is a global weight that is attributed to all events.
    :param nevts: (maximum) number of generated events. If None, all events from \
      file are processed.
    :returns: a TopologyList object 
    """
  import LHEReader, TopologyBuilder, SMSDataObjects, types
  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=SMSDataObjects.TopologyList ( )
  for Event in reader:
    ## Get mother PDGs:
    momPDG = tuple(Event.getMom())
    # Get event weight list:
    weight = {}
    if type(W)==types.DictType:
      for k in W.keys(): 
        if W[k].has_key(momPDG):
          weight[k]=W[k][momPDG]
        else:
          print "[LHEDecomposer] Error getting weight for",k,momPDG
          # weight[k]=0.
          # return None
    if type(W)==types.FloatType:
      ## a global weight
      weight=W

    # Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, doCompress, doInvisible, minmassgap)
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
