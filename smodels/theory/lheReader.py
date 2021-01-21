#!/usr/bin/env python3

"""
.. module:: lheReader
   :synopsis: Provides a class that creates SMSEvents from LHE files.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
"""

from smodels.tools.physicsUnits import TeV, pb
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
import pyslha

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

        # Get global information from file (cross section sqrts, total
        # cross section, total number of events)
        self.file.seek(0)
        line = self.file.readline()
        nevts = None
        totxsec = None
        sqrts = None
        # Exit if reached end of events or file
        while not "</LesHouchesEvents>" in line and line != "":
            if "<init>" in line:
                line = self.file.readline()
                if line.split()[0] == line.split()[1] == "2212":
                    sqrts = (eval(line.split()[2]) + eval(line.split()[3])) / 1000. * TeV
                    self.metainfo["sqrts"] = sqrts
                else: break
                line = self.file.readline()
                while not "</init>" in line:
                    if totxsec is None: totxsec = 0*pb
                    totxsec += eval(line.split()[0])* pb
                    line = self.file.readline()
                self.metainfo["totalxsec"] = totxsec
            elif "<event>" in line:
                if nevts is None: nevts = 0
                nevts += 1
            line = self.file.readline()
        self.metainfo["nevents"] = nevts
        # Return file to initial reader position
        self.file.seek(0)

    def close(self):
        """ close file handle """
        self.file.close()

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

    def __next__(self):
        """ for python3 """
        return self.next()

    def event(self):
        """
        Get next event.

        :returns: SmsEvent; None if no event is left to be read.

        """
        line = " "
        self.ctr += 1
        ret = SmsEvent(self.ctr)
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
            particle = LHEParticle()
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




class SmsEvent(object):
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
        Add particle to the event.

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
            raise SModelSError()
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


class LHEParticle(object):
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



def getDictionariesFrom(lheFile,nevts=None):
    """
    Reads all events in the LHE file and create mass and BR dictionaries for each branch in an event.

    :param lheFile: LHE file
    :param nevts: (maximum) number of events used in the decomposition. If
                  None, all events from file are processed.

    :returns: BR and mass dictionaries for the particles appearing in the event
    """

    reader = LheReader(lheFile, nevts)
    # Loop over events and decompose
    massDict = {}
    decaysDict = {}
    for event in reader:
        eventMass,eventDecays = getDictionariesFromEvent(event)
        for pdg,mass in eventMass.items():
            if not pdg in massDict:
                massDict[pdg] = mass
            else:
                massDict[pdg] += mass

        for pdg,decays in eventDecays.items():
            if not pdg in decaysDict:
                decaysDict[pdg] = decays
            else:
                decaysDict[pdg] += decays


    #Use averaged mass over all events:
    for pdg,masses in massDict.items():
        massDict[pdg] = sum(masses)/len(masses)
        if not pdg in decaysDict:
            decaysDict[pdg] = [] #Make sure all particles appear in decaysDict

    #Compute the decay dictionaries:
    for pdg,eventDecays in list(decaysDict.items()):
        daughterIDs = []
        for eventDec in eventDecays:
            pids = sorted(eventDec.ids)
            if not pids in daughterIDs:
                daughterIDs.append(pids)
        n = len(eventDecays) #Number of times the particle appears in all events
        combinedDecays = []
        for daughterID in daughterIDs:
            decay = pyslha.Decay(br=0.,nda=len(daughterID),ids=daughterID,parentid=pdg)
            for eventDecay in eventDecays:
                if sorted(eventDecay.ids) != daughterID:
                    continue
                decay.br += 1./float(n) #Count this decay
            combinedDecays.append(decay)
        decaysDict[pdg] = combinedDecays #Stable particles will appear with an empty list

    #Remove anti-particle decays and masses:
    for pdg in list(massDict.keys())[:]:
        if -abs(pdg) in massDict and abs(pdg) in massDict:
            massDict.pop(-abs(pdg))
    for pdg in list(decaysDict.keys())[:]:
        if -abs(pdg) in decaysDict and abs(pdg) in decaysDict:
            decaysDict.pop(-abs(pdg))

    #Add widths to decaysDict using the pyslha.Particle class:
    for pdg,decays in decaysDict.items():
        decaysDict[pdg] = pyslha.Particle(pid=pdg)
        if abs(pdg) in massDict:
            decaysDict[pdg].mass = massDict[abs(pdg)]
        if decays:
            decaysDict[pdg].totalwidth = float('inf')
        else:
            decaysDict[pdg].totalwidth = 0.
        decaysDict[pdg].decays = decays

    reader.close()

    return massDict,decaysDict


def getDictionariesFromEvent(event):
    """
    Create mass and BR dictionaries for each branch in an event.

    :param event: LHE event
    :returns: BR and mass dictionaries for the branches in the event

    """

    particles = event.particles

    # Identify and label to which branch each particle belongs
    #(the branch label is the position of the primary mother)
    masses = {}
    decays = {}
    for ip, particle in enumerate(particles):
        if particle.status == -1: #Ignore incoming particles
            continue

        if not particle.pdg in masses:
            masses[particle.pdg] = [particle.mass]
            decays[particle.pdg] = []

        #Get event decay:
        decay = pyslha.Decay(0, 0, [], particle.pdg)
        #Collect all daughters:
        for ipn,newparticle in enumerate(particles):
            if ipn == ip:
                continue
            if not ip+1 in newparticle.moms: #newparticle is not a daughter of particle
                continue

            # Check if particle has a single parent
            # (as it should)
            newparticle.moms = [mom for mom in newparticle.moms if mom != 0]
            if len(newparticle.moms) != 1:
                raise SModelSError("More than one parent particle found")

            decay.nda += 1
            decay.br = 1.
            decay.ids.append(newparticle.pdg)

        if decay.ids:
            decay.ids = sorted(decay.ids)
            decays[particle.pdg].append(decay)

    return masses,decays
