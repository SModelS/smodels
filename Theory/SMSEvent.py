""" Simple classes that encapsulate the information of an event """

class Particle:
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
        # return "particle pdg %d p=(%f,%f,%f,m=%f) st %d" % ( self.pdg, self.px, self.py, self.pz, self.mass, self.status )
        return "particle pdg %d p=(%.1f,%.1f,%.1f,m=%.1f) status %d moms %s" % ( self.pdg, self.px, self.py, self.pz, self.mass, self.status, self.moms )


class SMSEvent:
  """ a super-simple event class. Basically, it's a list of Particles,
      plus some convenience functions """
  def __init__ ( self, eventnr=None ):
    self.particles=[]
    self.eventnr=eventnr

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
    nr=""
    if self.eventnr!=None: nr=" "+str(self.eventnr)
    ret="\nEvent%s:\n" % nr
    for p in self.particles:
      ret+=p.__str__()+"\n"
    return ret
      
