""" A facility to read in lhe files and generate events. 
    An event is essentially a list of particles """

import SMSEvent

class LHEReader:
  """ a class that produces Events from lhe files """

  def __init__ ( self, filename ):
    self.filename=filename
    self.File = open ( filename )

  def event ( self ):
    """ reads LHE file 'File', returns an Event.
        File is an open filehandle.  """
    line = " "
    ret=SMSEvent.SMSEvent()
    # PartList = []

#Find next event
    while line.find("<event>") == -1:
      if line=='': 
        print "[SMSmethods.py] error demanding more events than are available."
        return ret
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
        particle = SMSEvent.MParticle()
        linep = [float(x) for x in line.split()]
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
