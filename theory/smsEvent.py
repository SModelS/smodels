"""
.. module:: SmsEvent
   :synopsis: Provides a class that encapsulates an LHE or SLHA event.
        
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
        
"""

import logging

logger = logging.getLogger(__name__) # pylint: disable-msg=C0103
        
class SmsEvent:
    """
    Event class featuring a list of particles and some convenience functions.
    
    """
    def __init__(self, eventnr=None):
        self.particles = []
        self.eventnr = eventnr
        self.metainfo = {}
        

    def metaInfo(self, key):
        """
        Return the meta information of 'key', None if info does not exist.
        
        """
        if not key in self.metainfo:
            return None
        return self.metainfo[key]

    def add(self, particle):
        """
        TODO: write docstring
        
        """
        self.particles.append(particle)
        

    def getMom(self):
        """
        Return the pdgs of the mothers, None if a problem occurs.
        
        """
        momspdg = []
        imom = 0
        for p in self.particles:
            if len(p.moms) > 1 and p.moms[0] == 1 or p.moms[1] == 1:
                momspdg.append(p.pdg)
                imom += 1
        if imom != 2:
            logger.error("Number of mother particles %d != 2", imom)
            return None
        if momspdg[0] > momspdg[1]:
            momspdg[0], momspdg[1] = momspdg[1], momspdg[0]
        return momspdg
    

    def __str__(self):
        nr = ""
        if self.eventnr != None:
            nr = " " + str(self.eventnr)
        metainfo = ""
        for(key, value) in self.metainfo.items():
            metainfo += " %s:%s" % (key, value)
        ret = "\nEvent%s:%s\n" % (nr, metainfo)
        for p in self.particles:
            ret += p.__str__() + "\n"
        return ret
    

class Particle:
    """
    An instance of this class represents a particle.
    
    """ 
    def __init__(self):
        self.pdg = 0
        self.status = 0
        # moms is a list of the indices of the mother particles
        self.moms = [] 
        self.px = 0.
        self.py = 0.
        self.pz = 0.
        self.e = 0.
        self.mass = 0.
        # position in the event list of particles
        self.position = None    
        
    def __str__(self):
        return "particle pdg %d p=(%.1f,%.1f,%.1f,m=%.1f) status %d moms %s" \
                % (self.pdg, self.px, self.py, self.pz, self.mass,
                   self.status, self.moms)


            
