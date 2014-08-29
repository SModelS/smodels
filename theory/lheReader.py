#!/usr/bin/env python

"""
.. module:: LheReader
   :synopsis: Provides a class that creates SMSEvents from LHE files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import setPath
from smodels.theory import smsEvent
from smodels.tools.physicsUnits import addunit
import logging

logger = logging.getLogger(__name__)


class LheReader(object):
    """
    An instance of this class represents a reader for LHE files.
    
    """
    def __init__(self, filename, nmax=None):
        """
        Constructor.
        
        :param filename: LHE file name
        :param nmax: When using the iterator, then nmax is the maximum number
        of events to be reader, nmax=None means read till the end of the file.
        If filename is not a string, assume it is already a file object and do
        not open it.
        
        """
        self.filename = filename
        self.nmax = nmax
        self.ctr = 0
        if type(filename) == type('str'):
            self.file = open(filename, 'r')
        else: self.file = filename
        self.metainfo = {"nevents" : None, "totalxsec" : None, "sqrts" : None}

        # Get global information from file (cross-section sqrts, total
        # cross-section, total number of events)
        self.file.seek(0)
        line = self.file.readline()
        nevts = 0
        # Exit if reached end of events or file
        while not "</LesHouchesEvents>" in line and line != "":
            if "<init>" in line:
                line = self.file.readline()
                sqrts = addunit((eval(line.split()[2]) + \
                                 eval(line.split()[3])) / 1000., 'TeV')
                self.metainfo["sqrts"] = sqrts
                totxsec = addunit(0., 'pb')
                line = self.file.readline()
                while not "</init>" in line:
                    totxsec += addunit(eval(line.split()[0]), 'pb')
                    line = self.file.readline()
                self.metainfo["totalxsec"] = totxsec
            elif "<event>" in line:
                nevts += 1
            line = self.file.readline()
        self.metainfo["nevents"] = nevts
        # Return file to initial reader position
        self.file.seek(0)


    def next(self):
        """
        Get next element in iteration.
        
        Needed for the iterator.
        
        """
        if self.nmax != None and self.ctr == self.nmax:
            raise StopIteration
        e = self.event()
        if e == None:
            raise StopIteration
        return e


    def __iter__(self):
        """
        Make class iterable.
        
        Allows iterations like 'for a in lhereader: print a'.
        
        """
        return self


    def event(self):
        """
        Get next event.
        
        :returns: SmsEvent; None if no event is left to be read.
        
        """
        line = " "
        self.ctr += 1
        ret = smsEvent.SmsEvent(self.ctr)
        # Pass metainfo from file to event
        for (key, value) in self.metainfo.items():
            ret.metainfo[key] = value
        # Find next event
        while line.find("<event>") == -1:
            if line == '':
                return None
            line = self.file.readline()
        # Read event info
        line = self.file.readline()

        # Get particles info
        line = self.file.readline()
        while line.find("</event>") == -1:
            if line.find("#") > -1:
                line = line[:line.find('#')]
            if len(line) == 0:
                line = self.file.readline()
                continue
            particle = smsEvent.Particle()
            linep = [float(x) for x in line.split()]
            if len(linep) < 11:
                logger.error("Line >>%s<< in %s cannot be parsed",
                             line, self.filename)
                line = self.file.readline()
                continue
            particle.pdg = int(linep[0])
            particle.status = int(linep[1])
            particle.moms = [int(linep[2]), int(linep[3])]
            particle.px = linep[6]
            particle.py = linep[7]
            particle.pz = linep[8]
            particle.e = linep[9]
            particle.mass = linep[10]

            ret.add(particle)
            line = self.file.readline()
        return ret


    def close(self):
        """
        Close the lhe file, if open.
        
        """
        self.file.close()

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser( "The LHE file reader class." )
    argparser.add_argument('-f', '--filename',help = 'filename of input lhe file')
    args = argparser.parse_args()
    reader = LheReader ( args.filename )
    print "Reading",reader.filename
    for event in reader:
        print event
    reader.close()
