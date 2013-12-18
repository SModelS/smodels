#!/usr/bin/env python

"""
.. module:: NLLXSec
    :synopsis: This module returns SUSY strong production cross-section using \
    NLLfast for 7 and 8 TeV. The return output is a dictionary containing the \
    cross-sections at NLO and NLL+NLO along with the k factors. For the \
    decoupling limit to be obeyed the ratio of the corrospoding partiecles \
    should be > 10. 

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>, Andre Lessa <lessa.a.p@gmail.com>

"""

import sys, commands, os, pyslha
from Tools.PhysicsUnits import addunit
import numpy

squarks = [1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004]

def getNLLfast(process = "gg", pdf = 'cteq', squarkmass=0., gluinomass=0., Energy = 8, base="./nllfast/", interpolate=True ):
    """ method to get the NLL fast results. 
        This is the shortest way to code, I found. """
    import logging
    log = logging.getLogger(__name__)
    energy = str(int(Energy))+'TeV'
    NoValues={"sqmass_"+energy:0, "mgmass_"+energy:0,"LOcs_"+energy:False,"NLOcs_"+energy:False,"NLL+NLO_"+energy:False,"K_NLO_"+energy:False,"K_NLL_"+energy:False}

    mass = 0
    o=None   
    
    inprocess = process

    wd=os.getcwd()
    
    if Energy==8: nllpath = "%s/nllfast-2.1" % base
    elif Energy==7: nllpath = "%s/nllfast-1.2" % base
    if not os.path.isdir(nllpath):
      print "[NLLXsec.py] error:",nllpath,"does not exist or is not a directory."
      return NoValues
    if not os.path.isfile(nllpath+"/nllfast_"+energy) or not os.access(nllpath+"/nllfast_"+energy, os.X_OK):
      print "[NLLXsec.py] error: %s/nllfast-2.1/nllfast_8TeV does not exist or is not executable." % base
      return NoValues
    
    try:
      os.chdir(nllpath)
      if process=="st":    
        s = "./nllfast_"+energy+" %s %s %s" % ( process, pdf, squarkmass )
      else:
        s = "./nllfast_"+energy+" %s %s %s %s" % ( process, pdf, squarkmass, gluinomass )
      o=commands.getoutput(s)
    except Exception,e:
      print "[NLLXsec.py] caught exception",e
    os.chdir(wd) ## back to where we started

    xmass = squarkmass
    ymass = gluinomass

    if process == 'st' and o[1]=='T': return NoValues
    if process != 'st' and o[1]=='T':   #Check for possible errors/decoupling limit
      if process=='sb' and gluinomass > 500. and 'gluino' in o.lower(): process='sdcpl'  #Gluino mass is too high (use decoupling limit)
      elif process=='gg' and squarkmass > 500. and 'squark' in o.lower(): process='gdcpl'  #Squark mass is too high (use decoupling limit)
      else: return NoValues
      
      try:
        xmass = max(squarkmass,gluinomass)
        ymass = min(squarkmass,gluinomass)  #Use non-decoupled mass
        os.chdir(nllpath)
        s="./nllfast_"+energy+" %s %s %s" % ( process, pdf, ymass )        
        o=commands.getoutput(s)
      except Exception,e:
        print "[NLLXSec.py 2] caught",e

    os.chdir(wd) ## back to where we started
    if o[1]=='T':
#        print "[NLLXSec.py] how can I end up here? A. o=",o,squarkmass,gluinomass
        return NoValues

    lines = o.split()
#Convert from strings to values whenever possible:    
    for iline,line in enumerate(lines):
      try: lines[iline] = eval(line)
      except: continue
      
    if len(lines) <= 38 or len(lines) >= 41:
      log.error ( "when parsing nllfast output, I have only %d lines. I guess something is wrong. Maybe sth is wrong with the nll binary?" % len(lines) )
      log.error ( "the lines are" + str(lines) )
    elif len(lines) == 39: lines.insert(27,0.)  #Make sure all output have the same format

#Try to interpolate if decoupling limit is not satisfied (M/m < 10):
    if 'dcpl' in process and interpolate and max(squarkmass,gluinomass)/min(squarkmass,gluinomass) < 10.:
        lines.insert
        os.chdir(nllpath)
        xmass = max(squarkmass,gluinomass)  #Decoupled mass
        ymass = min(squarkmass,gluinomass)  #Non-decoupled mass
        xpts = [10.*ymass]                          #Use the decoupled value (valid for xmass > 10*ymass) as one of the points for interpolation
        ypts = [lines[29:40]]
        while len(xpts) < 2 and xmass > 500.:
          xmass -= 100.
          if ymass == squarkmass: s = "./nllfast_"+energy+" %s %s %s %s" % ( inprocess, pdf, ymass, xmass )
          else: s = "./nllfast_"+energy+" %s %s %s %s" % ( inprocess, pdf, xmass, ymass )
          o=commands.getoutput(s)
          if o[1] =='T': continue
          xpts.append(xmass)                            #Get points in the non-decoupling regime for interpolation
          ypts.append([eval(x) for x in o.split()[29:40]])   #Use absolute values of x (does not affect xsec and k-factor interpolation)

        if len(xpts) > 1:
          xmass = max(squarkmass,gluinomass)
          coeffs = numpy.matrix.transpose(numpy.polyfit(xpts, ypts, len(xpts)-1))
          for i in range(11):
            newval = 0.
            for ip,coeff in enumerate(coeffs[i]): newval += coeff*xmass**(len(xpts)-1-ip)
            if i < 9 and lines[i+29] > 0.: lines[i+29] = max(0.,newval)   #Only interpolate for entries with positive values (xsecs)
            else: lines[i+29] = max(1.,newval)          #Force k-factors >= 1

    os.chdir( wd ) # make sure we always chdir back
#Add units
    for i in range(29,32): lines[i] = addunit(lines[i],'pb')
    msquark = squarkmass
    mgluino = gluinomass
    if process == 'st' or process == 'sdcpl':  mgluino = 0.
    elif process == 'gdcpl' : msquark = 0.    
    Values={"sqmass_"+energy: msquark, "mgmass_"+energy: mgluino,"LOcs_"+energy:lines[29],"NLOcs_"+energy:lines[30],"NLL+NLO_"+energy:lines[31],"K_NLO_"+energy:lines[38],"K_NLL_"+energy:lines[39]}

    return Values

def getNLLresult(pdgid1,pdgid2,inputfile,base="../nllfast",pdf="cteq"):
    """ obtain the nll xsecs for pdgid1/pdgid2 production.

       :param pdgid1: pdg id of first mother
       :param pdgid2: pdg id of second mother
       :param inputfile: slha file name
       :type inputfile: str
       :param base: base of nllfast installation
       :param pdf: pdf tag
       :returns: an array of two dictionaries, one for 7, one for 8 tev. if for \
       the given mothers no results can be obtained, the dictionary \
         { "xsecs": False, masses=0 } will be returned. \
         the same is true if the masses of the decoupled particles in the \
         input file are too light (ratio to mothers must be < 10 )
    """ 
    
    readfile=pyslha.readSLHAFile(inputfile)
    gluinomass = abs(readfile[0]['MASS'].entries[1000021])
    squarkmass = (abs(readfile[0]['MASS'].entries[1000001])+abs(readfile[0]['MASS'].entries[2000001])+abs(readfile[0]['MASS'].entries[1000002])+abs(readfile[0]['MASS'].entries[2000002])+abs(readfile[0]['MASS'].entries[1000003])+abs(readfile[0]['MASS'].entries[2000003])+abs(readfile[0]['MASS'].entries[1000004])+abs(readfile[0]['MASS'].entries[2000004]))/8

    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
    output=[Values_7TeV,Values_8TeV]
    
    process = None
    if pdgid1 < 0 and abs(pdgid1) in squarks and pdgid2 in squarks:
      process = 'sb'
      squarkmass = (abs(readfile[0]['MASS'].entries[abs(pdgid1)])+abs(readfile[0]['MASS'].entries[pdgid2]))/2
    elif pdgid1 in squarks and pdgid2 in squarks:
      process = 'ss'
      squarkmass = (abs(readfile[0]['MASS'].entries[pdgid1])+abs(readfile[0]['MASS'].entries[pdgid2]))/2
    elif pdgid1 == pdgid2 == 1000021: process = 'gg'
    elif pdgid1 in squarks and pdgid2 == 1000021: process = 'sg'
    elif abs(pdgid1) == pdgid2 == 1000005 or abs(pdgid1) == pdgid2 == 2000005 or abs(pdgid1) == pdgid2 == 1000006 or abs(pdgid1) == pdgid2 == 2000006:
      process = 'st'
      squarkmass = abs(readfile[0]['MASS'].entries[abs(pdgid1)])

    if process:
      Values_7TeV = getNLLfast(process,pdf,squarkmass,gluinomass,7,base=base)
      Values_8TeV = getNLLfast(process,pdf,squarkmass,gluinomass,8,base=base)

    output = [Values_7TeV,Values_8TeV]

    return output

if __name__ == "__main__":

    print getNLLresult(pdgid1=1000021,pdgid2=1000021,inputfile="../addiotionalSLHA/T1.slha")
#    print getNLLresult(pdgid1=1000021,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=1000006,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=-1000001,pdgid2=1000001,inputfile="./AndreSLHA/andrePT4_var.slha")
