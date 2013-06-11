#!/usr/bin/env python

"""
.. module:: NLLXSec
    :synopsis: This module returns SUSY strong production cross-section using \
    NLLfast for 7 and 8 TeV. The return output is a dictionary containing the \
    cross-sections at NLO and NLL+NLO along with the k factors. For the \
    decoupling limit to be obeyed the ratio of the corrospoding partiecles \
    should be > 10. 

.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>

"""

import sys, commands, os, pyslha2
from Tools.PhysicsUnits import addunit

squarks = [1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004]

def getNLLfast(process = "gg", pdf = 'cteq', squarkmass=0., gluinomass=0., Energy = 8, base="./nllfast/" ):
    """ method to get the NLL fast results. 
        This is the shortest way to code, I found. """
    Values={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}

    mass = 0
    o=None

    wd=os.getcwd()

    try:
      if Energy==8:
          if not os.path.isdir ( "%s/nllfast-2.1" % base ):
            print "[NLLXsec.py] error: %s/nllfast-2.1 does not exist or is not a directory." % base
            return Values_8TeV
          if not os.path.isfile ( "%s/nllfast-2.1/nllfast_8TeV" % base ) or not os.access ( "%s/nllfast-2.1/nllfast_8TeV" % base, os.X_OK ):
            print "[NLLXsec.py] error: %s/nllfast-2.1/nllfast_8TeV does not exist or is not executable." % base
            return Values_8TeV
          if process=="st":
              os.chdir("%s/nllfast-2.1" % base )
              s = "./nllfast_8TeV %s %s %s" % ( process, pdf, squarkmass )
          else:
              os.chdir("%s/nllfast-2.1" % base )
              s = "./nllfast_8TeV %s %s %s %s" % ( process, pdf, squarkmass, gluinomass )
      if Energy==7:
          if not os.path.isdir ( "%s/nllfast-1.2" % base ):
            print "[NLLXsec.py] error: %s/nllfast-1.2 does not exist or is not a directory." % base
            return Values_7TeV
          if not os.path.isfile ( "%s/nllfast-1.2/nllfast_7TeV" % base ) or not os.access ( "%s/nllfast-1.2/nllfast_7TeV" % base, os.X_OK ):
            print "[NLLXsec.py] error: %s/nllfast-1.2/nllfast_7TeV does not exist or is not executable." % base
            return Values_7TeV
          if process=="st":
              os.chdir("%s/nllfast-1.2" % base )
              s = "./nllfast_7TeV %s %s %s" % ( process, pdf, squarkmass )
          else:
              os.chdir("%s/nllfast-1.2" % base )
              s = "./nllfast_7TeV %s %s %s %s" % ( process, pdf, squarkmass, gluinomass )
      o=commands.getoutput(s)
    except Exception,e:
      print "[NLLXsec.py] caught exception",e
    os.chdir(wd) ## back to where we started

    if process=='st' and o[1]=='T' and Energy == 7:
        return Values_7TeV
    elif process=='st' and o[1]=='T' and Energy == 8:
        return Values_8TeV

    if process!='st' and o[1]=='T':
        if Energy == 7:
            if process=='gg' and squarkmass/gluinomass > 10:
              process='gdcpl'
              mass = gluinomass
            elif process=='sb' and gluinomass/squarkmass > 10:
              process='sdcpl'
              mass = squarkmass
            else: return Values_7TeV
        if Energy == 8:
            if process=='gg' and squarkmass/gluinomass > 10:
              process='gdcpl'
              mass = gluinomass
            elif process=='sb'and gluinomass/squarkmass > 10:
              process='sdcpl'
              mass = squarkmass
            else: return Values_8TeV
        try:
          if Energy==8:
              if not os.path.isdir ( "%s/nllfast-2.1" % base ):
                print "[NLLXsec.py] error: %s/nllfast-2.1 does not exist or is not a directory." % base
                return Values_8TeV
              if not os.path.isfile ( "%s/nllfast-2.1/nllfast_8TeV" % base ) or not os.access ( "%s/nllfast-2.1/nllfast_8TeV" % base, os.X_OK ):
                print "[NLLXsec.py] error: %s/nllfast-2.1/nllfast_8TeV does not exist or is not executable." % base
                return Values_8TeV
              s="./nllfast_8TeV  %s %s %s" % ( process, pdf, mass )
              os.chdir("%s/nllfast-2.1" % base )
          if Energy==7:
              if not os.path.isdir ( "%s/nllfast-1.2" % base ):
                print "[NLLXsec.py] error: %s/nllfast-1.2 does not exist or is not a directory." % base
                sys.exit(0)
              if not os.path.isfile ( "%s/nllfast-1.2/nllfast_7TeV" % base ) or not os.access ( "%s/nllfast-1.2/nllfast_7TeV" % base, os.X_OK ):
                print "[NLLXsec.py] error: %s/nllfast-1.2/nllfast_7TeV does not exist or is not executable." % base
                sys.exit(0)
              os.chdir("%s/nllfast-1.2" % base )
              s="./nllfast_7TeV %s %s %s" % ( process, pdf, mass )
          o=commands.getoutput(s)
        except Exception,e:
          print "[NLLXSec.py 2] caught",e
        os.chdir( wd ) # make sure we always chdir back

# uncomment this line to see the nll fast output
#   print "nllfast output is:  ", o

    if o[1]=='T' and Energy  == 7:
        print "[NLLXSec.py] how can I end up here? A. o=",o
        return Values_7TeV
    elif o[1]=='T' and Energy  == 8:
        print "[NLLXSec.py] how can I end up here? B. o=",o
        return Values_8TeV

    lines=o.split()

#Convert from strings to values whenever possible:
    for il in range(len(lines)):
        try:
            lines[il] = eval(lines[il])
        except:
            continue
#Add units
    lines[27] = addunit(lines[27],'GeV')
    if process == 'st' or process == 'sdcpl' or process == 'gdcpl':
        lines[28] = addunit(lines[28],'pb')
        lines[29] = addunit(lines[29],'pb')
        lines[30] = addunit(lines[30],'pb')
    else:
        lines[28] = addunit(lines[28],'GeV')
        lines[29] = addunit(lines[29],'pb')
        lines[30] = addunit(lines[30],'pb')
        lines[31] = addunit(lines[31],'pb')


    if Energy == 7:
        if process == 'st' or process == 'sdcpl':
            Values={"sqmass_7TeV":lines[27], "mgmass_7TeV":0.0,"LOcs_7TeV":lines[28],"NLOcs_7TeV":lines[29],"NLL+NLO_7TeV":lines[30],"K_NLO_7TeV":lines[37],"K_NLL_7TeV":lines[38]}
        elif process == 'gdcpl' :
            print "gluinos decoupling limit 7 TeV"
            Values={"sqmass_7TeV":0.0, "mgmass_7TeV":lines[27],"LOcs_7TeV":lines[28],"NLOcs_7TeV":lines[29],"NLL+NLO_7TeV":lines[30],"K_NLO_7TeV":lines[37],"K_NLL_7TeV":lines[38]}
        else:
            Values={"sqmass_7TeV":lines[27], "mgmass_7TeV":lines[28],"LOcs_7TeV":lines[29],"NLOcs_7TeV":lines[30],"NLL+NLO_7TeV":lines[31],"K_NLO_7TeV":lines[38],"K_NLL_7TeV":lines[39]}
    if Energy == 8:
        if process == 'st' or process == 'sdcpl':
            Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":0.0,"LOcs_8TeV":lines[28],"NLOcs_8TeV":lines[29],"NLL+NLO_8TeV":lines[30],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
        elif process == 'gdcpl' :
            print "gluinos decoupling limit 7 TeV"
            Values={"sqmass_8TeV":0.0, "mgmass_8TeV":lines[27],"LOcs_8TeV":lines[28],"NLOcs_8TeV":lines[29],"NLL+NLO_8TeV":lines[30],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
        else:
            Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":lines[28],"LOcs_8TeV":lines[29],"NLOcs_8TeV":lines[30],"NLL+NLO_8TeV":lines[31],"K_NLO_8TeV":lines[38],"K_NLL_8TeV":lines[39]}


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
    
    readfile=pyslha2.readSLHAFile(inputfile)
    gluinomass = abs(readfile[0]['MASS'].entries[1000021])
    squarkmass = (abs(readfile[0]['MASS'].entries[1000001])+abs(readfile[0]['MASS'].entries[2000001])+abs(readfile[0]['MASS'].entries[1000002])+abs(readfile[0]['MASS'].entries[2000002])+abs(readfile[0]['MASS'].entries[1000003])+abs(readfile[0]['MASS'].entries[2000003])+abs(readfile[0]['MASS'].entries[1000004])+abs(readfile[0]['MASS'].entries[2000004]))/8

    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
    output=[Values_7TeV,Values_8TeV]

    if pdgid1 < 0 and abs(pdgid1) in squarks and pdgid2 in squarks:
      squarkmass = (abs(readfile[0]['MASS'].entries[abs(pdgid1)])+abs(readfile[0]['MASS'].entries[pdgid2]))/2
      Values_7TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,7,base=base)
      Values_8TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,8,base=base)
      output = [Values_7TeV,Values_8TeV]

    if pdgid1 in squarks and pdgid2 in squarks:
      squarkmass = (abs(readfile[0]['MASS'].entries[pdgid1])+abs(readfile[0]['MASS'].entries[pdgid2]))/2
      Values_7TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,7,base=base)
      Values_8TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,8,base=base)
      output = [Values_7TeV,Values_8TeV]

    if pdgid1 == pdgid2 == 1000021:
        Values_7TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,7,base=base)
        Values_8TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,8,base=base)
        output = [Values_7TeV,Values_8TeV]

    if pdgid1 in squarks and pdgid2 == 1000021:
      Values_7TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,7,base=base)
      Values_8TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,8,base=base)
      output = [Values_7TeV,Values_8TeV]

    if abs(pdgid1) == pdgid2 == 1000005:
        squarkmass = abs(readfile[0]['MASS'].entries[1000005])
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7,base=base)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8,base=base)
        output = [Values_7TeV,Values_8TeV]

    if abs(pdgid1) == pdgid2 == 2000005:
        squarkmass = abs(readfile[0]['MASS'].entries[2000005])
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7,base=base)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8,base=base)
        output = [Values_7TeV,Values_8TeV]

    if abs(pdgid1) == pdgid2 == 1000006:
        squarkmass = abs(readfile[0]['MASS'].entries[1000006])
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7,base=base)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8,base=base)
        output = [Values_7TeV,Values_8TeV]

    if abs(pdgid1) == pdgid2 == 2000006:
        squarkmass = abs(readfile[0]['MASS'].entries[2000006])
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7,base=base)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8,base=base)
        output = [Values_7TeV,Values_8TeV]

    return output

if __name__ == "__main__":

    print getNLLresult(pdgid1=1000021,pdgid2=1000021,inputfile="../addiotionalSLHA/T1.slha")
#    print getNLLresult(pdgid1=1000021,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=1000006,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=-1000001,pdgid2=1000001,inputfile="./AndreSLHA/andrePT4_var.slha")
