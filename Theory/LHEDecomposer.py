def decompose(lhefile,W,nevts,DoCompress=False,DoInvisible=False,minmassgap=-1):
  """ Do LHE-based decomposition.
      lhefile = LHE file with pythia events
      W = dictionary with event weights
      nevts = number of generated events """
  import LHEReader, TopologyBuilder

  reader = LHEReader.LHEReader(lhefile,nevts)
  SMSTopList=TopologyList ( )
  for iev in range(nevts):
##Read event  
    Event = reader.next()
##Get mother PDGs:
    momPDG = tuple(Event.getMom())
    PList = Event.particles
#Get event weight list:
    weight = {}
    for k in W.keys(): 
      if W[k].has_key(momPDG):
        weight[k]=W[k][momPDG]
      else:
        print "LHEdecomp: Error getting weight"
        return False

#Get event topology  
    SMSTopListEv = TopologyBuilder.fromEvent(Event, weight, DoCompress, DoInvisible, minmassgap)
  
    SMSTopList.addList ( SMSTopListEv )
  return SMSTopList
