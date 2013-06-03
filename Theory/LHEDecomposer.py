#!/usr/bin/python

"""
.. module:: LHEDecomposition
    :synopsis: I have no idea ...

.. moduleauthor:: someone <email@example.com>

"""

def decompose(lhefile,W,nevts=None,doCompress=False,doInvisible=False,minmassgap=-1):
  """ Do LHE-based decomposition.  lhefile = LHE file with pythia events, W =
    dictionary with event weights, nevts = (maximum) number of generated events.  
    If nevts = None, process all events from file.
    Output is a TopologyList object """
  import LHEReader, TopologyBuilder, SMSDataObjects
  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=SMSDataObjects.TopologyList ( )
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
        return None

    # Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, doCompress, doInvisible, minmassgap)
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
