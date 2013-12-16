#!/usr/bin/env python

"""
.. module:: LHEReader
    :synopsis: A class that creates SMSEvents from LHE files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import SMSEvent
from Tools.PhysicsUnits import addunit

class LHEReader:
    def __init__ ( self, filename, nmax=None ):
        """ constructor.
          :param filename: LHE file name
          :param nmax: when using the iterator, then nmax is the maximum number of \ 
          events to be reader, nmax=None means read till the end of the file. \
          if filename is not a string, assume it is already a File object and do not open it
        """
        self.filename=filename
        self.nmax=nmax
        self.ctr=0
        if type(filename) == type('str'):  self.File = open(filename, 'r')
        else: self.File = filename
        self.metainfo = {"nevents" : None, "totalxsec" : None, "sqrts" : None}
        
#Get header information from file (cross-section sqrts, total cross-section, total number of events)
        self.File.rewind()        
        while None in self.metainfo.values():
            line = self.File.readline()
            if line == "": break             #Exit if reached end of file
            if line.find("Number of Events        :") > -1:
                nevts = int(line.split()[-1])
                self.metainfo["nevents"] = nevts       
            elif line.find("Integrated weight (pb)") > -1:
                iwght = addunit(float(line.split()[-1]),'pb')
                self.metainfo["totalxsec"] = iwght
            elif line.find("<init>") > -1:
                line = self.File.readline()
                self.metainfo["sqrts"] = addunit(line.split()[2] + line.split()[3],'GeV')
        
        self.File.rewind()  #Return file to initial reader position
        
        
    def next ( self ):
        """ needed for the iterator """
        if self.nmax!=None and self.ctr==self.nmax: raise StopIteration
        e=self.event()
        if e==None: raise StopIteration
        return e

    def __iter__ ( self ):
        """ iterator, to allow constructs like 'for a in lhereader: print a' """
        return self

    def event ( self ):
        """ get next event.
          :returns: an SMSEvent; None if no event is left to be read.
        """
        
        line = " "
        self.ctr+=1
        ret=SMSEvent.SMSEvent(self.ctr)
        for (key,value) in self.metainfo.items(): ret.metainfo[key]=value  # pass metainfo from file to event
#Find next event
        while line.find("<event>") == -1:
            if line=='': return None
            line = self.File.readline()
#Read event info:
        line = self.File.readline()

#Get particles info:            
        line = self.File.readline()  
        while line.find("</event>") == -1:
            if line.find("#")>-1: line=line[:line.find('#')]
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
    
    def close(self):
        """ Closes the lhe file, if open """
        
        self.File.close()
