def decompose(lhefile,W,nevts=None,DoCompress=False,DoInvisible=False,minmassgap=-1):
  """ Do LHE-based decomposition.  lhefile = LHE file with pythia events, W =
    dictionary with event weights, nevts = (maximum) number of generated events.  
    If nevts = None, process all events from file.
    Output is a TopologyList object """
  import LHEReader, TopologyBuilder
  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=TopologyList ( )
  for Event in reader:
    ## Get mother PDGs:
    momPDG = tuple(Event.getMom())
    PList = Event.particles
    # Get event weight list:
    weight = {}
    for k in W.keys(): 
      if W[k].has_key(momPDG):
        weight[k]=W[k][momPDG]
      else:
        print "LHEdecomp: Error getting weight"
        return False

    # Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, DoCompress, DoInvisible, minmassgap)
  
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
