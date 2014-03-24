"""
.. module:: lheReader
   :synopsis: A class that creates SMSEvents from LHE files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

from . import smsEvent
from tools.physicsUnits import addunit
import logging

logger = logging.getLogger(__name__)

class lheReader:
    def __init__(self, filename, nmax=None):
        """
        Constructor.
        
        :param filename: LHE file name
        :param nmax: When using the iterator, then nmax is the maximum number
        of events to be reader, nmax=None means read till the end of the file.
        If filename is not a string, assume it is already a File object and do
        not open it.
        
        """
        self.filename=filename
        self.nmax=nmax
        self.ctr=0
        if type(filename) == type('str'):
            self.File = open(filename, 'r')
        else: self.File = filename
        self.metainfo = {"nevents" : None, "totalxsec" : None, "sqrts" : None}
        
        # Get global information from file (cross-section sqrts, total
        # cross-section, total number of events)
        self.File.seek(0)
        line = self.File.readline()
        nevts = 0
        # Exit if reached end of events or file
        while not "</LesHouchesEvents>" in line and line != "":   
            if "<init>" in line:
                line = self.File.readline()
                self.metainfo["sqrts"] = addunit(eval(line.split()[2]) + \
                                                 eval(line.split()[3]),'GeV')
                totxsec = addunit(0.,'pb')
                line = self.File.readline()
                while not "</init>" in line:
                    totxsec += addunit(eval(line.split()[0]),'pb')
                    line = self.File.readline()                    
                self.metainfo["totalxsec"] = totxsec
            elif "<event>" in line: nevts += 1
            line = self.File.readline()        
        self.metainfo["nevents"] = nevts
        # Return file to initial reader position
        self.File.seek(0)  
        
        
    def next ( self ):
        """
        Needed for the iterator.
        
        """
        if self.nmax!=None and self.ctr==self.nmax: raise StopIteration
        e=self.event()
        if e==None: raise StopIteration
        return e

    def __iter__ ( self ):
        """
        Iterator, to allow constructs like 'for a in lhereader: print a'.
        
        """
        return self

    def event ( self ):
        """
        Get next event.
        
        :returns: an smsEvent; None if no event is left to be read.
        
        """        
        line = " "
        self.ctr+=1
        ret=smsEvent.smsEvent(self.ctr)
        # pass metainfo from file to event
        for (key,value) in self.metainfo.items(): ret.metainfo[key]=value
        # Find next event
        while line.find("<event>") == -1:
            if line=='': return None
            line = self.File.readline()
        # Read event info:
        line = self.File.readline()

        # Get particles info:            
        line = self.File.readline()  
        while line.find("</event>") == -1:
            if line.find("#") > -1:
                line = line[:line.find('#')]
            if len(line) == 0:
                line = self.File.readline()
                continue
            particle = smsEvent.Particle()
            linep = [float(x) for x in line.split()]
            if len(linep)<11:
                logger.error("Line >>%s<< in %s cannot be parsed" % \
                             (line, self.filename))
                line = self.File.readline()
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
    
    def close(self):
        """
        Closes the lhe file, if open.
        
        """        
        self.File.close()
