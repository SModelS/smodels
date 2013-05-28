""" A facility to read in lhe files and generate events. 
    an event is essentially a list of particles """

import SMSEvent

class LHEReader:
  """ a class that produces events from lhe files """

  def __init__ ( self, filename, nmax=None ):
    """ initiate reader with 'filename'.
        When using the iterator, then nmax is the maximum number of events to
        be reader, nmax=None means read till the end of the file. """
    self.filename=filename
    self.nmax=nmax
    self.ctr=0
    self.File = open ( filename )

  def next ( self ):
    """ needed for the iterator """
    if self.nmax!=None and self.ctr==self.nmax:
      raise StopIteration
    e=self.event()
    if e==None:
      raise StopIteration
    return e

  def __iter__ ( self ):
    """ iterator, to allow constructs like 'for a in lhereader: print a' """
    return self

  def event ( self ):
    """ reads lhe file 'file', returns an event.
        file is an open filehandle. returns none if no event is left to be read. """
    line = " "
    self.ctr+=1
    ret=SMSEvent.SMSEvent( self.ctr )
    # PartList = []

#Find next event
    while line.find("<event>") == -1:
      if line=='': 
        # print "[LHEReader.py] error demanding more events than are available."
        return None
      line = self.File.readline()
        
#Read event info:
    line = self.File.readline()

#Get particles info:            
    line = self.File.readline()  
    while line.find("</event>") == -1:
        if line.find("#")>-1:
          line=line[:line.find('#')]
        if len(line)==0:
          line=self.File.readline()
          continue
        particle = SMSEvent.Particle()
        linep = [float(x) for x in line.split()]
        if len(linep)<11:
          print "[LHEReader] dont understand the following line >>%s<< in %s" % (line, self.filename )
          line=self.File.readline()
          continue
        particle.pdg = int(linep[0])
        particle.status = int(linep[1])
        particle.moms = [int(linep[2]),int(linep[3])]
        particle.px = linep[6]
        particle.py = linep[7]
        particle.pz = linep[8]
        particle.e = linep[9]
        particle.mass = linep[10]
        
        ret.add(particle)
        line = self.File.readline()  

    return ret
