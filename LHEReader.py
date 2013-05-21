""" A facility to read in lhe files and generate events. 
    An event is essentially a list of particles """

class MParticle:
    """ simple helper class to encapsulate one particle """ 
    def __init__(self):
        self.pdg = 0
        self.status = 0
        self.moms = []
        self.px = 0.
        self.py = 0.
        self.pz = 0.
        self.e = 0.
        self.mass = 0.
    def __str__(self):
        return "particle pdg %d p=(%f,%f,%f,m=%f) st %d" % ( self.pdg, self.px, self.py, self.pz, self.mass, self.status )


class Event:
  """ a super-simple event class. Basically, it's a list of MParticles """
  def __init__ ( self ):
    self.particles=[]

  def add ( self, particle ):
    self.particles.append ( particle )

  def getMom(self):
    """ returns the pdgs of the mothers, None if a problem occured """
    momspdg = []
    imom = 0
    for p in self.particles:
      if len(p.moms)>1 and p.moms[0] == 1 or p.moms[1] == 1:
        momspdg.append( p.pdg)
        imom += 1
    if imom != 2:
        print "getMom: Number of mother particles != 2"
        return None
    if momspdg[0] > momspdg[1]:
        momspdg[0], momspdg[1] = momspdg[1], momspdg[0]
    return momspdg


  def __str__ ( self ):
    ret="\nEvent:\n"
    for p in self.particles:
      ret+=p.__str__()+"\n"
    return ret
      

class LHEReader:
  """ a class that reads in an lhe file and produces lists of MParticles """

  def __init__ ( self, filename ):
    self.filename=filename
    self.File = open ( filename )

  def event ( self ):
    """ reads LHE file 'File', returns an Event.
        File is an open filehandle.  """
    line = " "
    ret=Event()
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
        particle = MParticle()
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
