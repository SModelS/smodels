#!/usr/bin/env python

import sys, commands, os, random, numpy
workdir = os.getcwd()    
pyslhadir = workdir + "/pyslha-1.4.3"
sys.path.append(pyslhadir)
import pyslha

#define global variables
pdf = 'cteq'


# method to get the NLL fast results. This is the shortest way to code, I found. 
def getNLLfast(process = "gg", pdf = 'cteq', squarkmass=0., gluinomass=0., Energy = 8 ):
  from string import Template

  Values={}

  if Energy==8:
    if process=="st": 
      os.chdir("./NLL_fast/nllfast-2.1")
      runnll = Template("./nll_fast_8TeV $prodproc $in_pdf $prodsq")
      s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass)
    else: 
      os.chdir("./NLL_fast/nllfast-2.1")
      runnll = Template("./nll_fast_8TeV $prodproc $in_pdf $prodsq $prodgl")
      s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass, prodgl=gluinomass)
      print "run command is  ", s
  if Energy==7:
    if process=="st": 
      os.chdir("./NLL_fast/nllfast-1.2")
      runnll = Template("./nll_fast_7TeV $prodproc $in_pdf $prodsq")
      s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass)
    else:
      os.chdir("./NLL_fast/nllfast-1.2")
      runnll = Template("./nll_fast_7TeV $prodproc $in_pdf $prodsq $prodgl")
      s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass, prodgl=gluinomass)
  o=commands.getoutput(s)
  os.chdir("./../../")


  if process=='st':
    if o[1]=='T': 
      Values={"sqmass":0, "mgmass":0,"LOcs":False,"NLOcs":False,"NLL+NLO":False,"K_NLO":False,"K_NLL":False}
      return Values

# the definition of decoupling limit is ill defined at the moment
  if process!='st':
    if o[1]=='T' and Energy == 7: 
      if process=='gg' and gluinomass >=200 and gluinomass <=2000 and squarkmass/gluinomass > 1.5: 
        process='gdcpl'
        mass = gluinomass
      else: return Values
      if process=='ss' and squarkmass>=200 and squarkmass<=2000 and gluinomass/squarkmass > 1.5: 
        process='sdcpl'
        mass = squarkmass
      else: return Values
    if o[1]=='T' and Energy == 8: 
      if process=='gg' and gluinomass >=200 and gluinomass <=2500 and squarkmass/gluinomass > 1.5: 
        process='gdcpl'
        mass = gluinomass
      else: return Values
      if process=='ss' and squarkmass>=200 and squarkmass<=2500 and gluinomass/squarkmass > 1.5: 
        process='sdcpl'
        mass = squarkmass
      else: return Values

      if Energy==8:
        os.chdir("./NLL_fast/nllfast-2.1")
        runnll = Template("./nll_fast_8TeV $prodproc $in_pdf $prodmass")
        s=runnll.substitute(prodproc=process, in_pdf = pdf, prodmass=mass)
      if Energy==7:
        os.chdir("./NLL_fast/nllfast-1.2")
        runnll = Template("./nll_fast_7TeV $prodproc $in_pdf $prodmass")
        s=runnll.substitute(prodproc=process, in_pdf = pdf, prodmass=mass)
        print "command for NLL fast is  ", s
      o=commands.getoutput(s)
      os.chdir("./../../")
        

  if o[1]=='T': 
    Values={"sqmass":0, "mgmass":0,"LOcs":False,"NLOcs":False,"NLL+NLO":False,"K_NLO":False,"K_NLL":False}
    return Values


  lines=o.split()
  print lines
  if Energy == 7:
    if process == 'st' or process == 'sdcpl':
      Values={"sqmass_7TeV":lines[27], "mgmass_7TeV":0.0,"LOcs_7TeV":lines[28],"NLOcs_7TeV":lines[29],"NLL+NLO_7TeV":lines[30],"K_NLO_7TeV":lines[37],"K_NLL_7TeV":lines[38]}
    elif process == 'gdcpl' :
      Values={"sqmass_7TeV":0.0, "mgmass_7TeV":lines[27],"LOcs_7TeV":lines[29],"NLOcs_7TeV":lines[30],"NLL+NLO_7TeV":lines[31],"K_NLO_7TeV":lines[37],"K_NLL_7TeV":lines[38]}
    else:
      Values={"sqmass_7TeV":lines[27], "mgmass_7TeV":lines[28],"LOcs_7TeV":lines[29],"NLOcs_7TeV":lines[30],"NLL+NLO_7TeV":lines[31],"K_NLO_7TeV":lines[38],"K_NLL_7TeV":lines[39]}
  if Energy == 8:
    if process == 'st' or process == 'sdcpl':
      Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":0.0,"LOcs_8TeV":lines[28],"NLOcs_8TeV":lines[29],"NLL+NLO_8TeV":lines[30],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
    elif process == 'gdcpl' :
      Values={"sqmass_8TeV":0.0, "mgmass_8TeV":lines[27],"LOcs_8TeV":lines[29],"NLOcs_8TeV":lines[30],"NLL+NLO_8TeV":lines[31],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
    else:
      Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":lines[28],"LOcs_8TeV":lines[29],"NLOcs_8TeV":lines[30],"NLL+NLO_8TeV":lines[31],"K_NLO_8TeV":lines[38],"K_NLL_8TeV":lines[39]}


  return Values


#various methods to get the NLL 7 and 8TeV output for specific processes. getggNLL will give you gluino gluino, getsgNLL will give squark gluino etc. 
def getggNLL(squarkmass,gluinomass):
    
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
    
    if gluinomass >=200 and gluinomass <=1500 and squarkmass >=200 and squarkmass <=4500:
        Values_7TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,7)
    elif gluinomass >1500 and gluinomass <=2000 and squarkmass >=200 and squarkmass <=2000:
        Values_7TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,7)
    
    if gluinomass >=200 and gluinomass <=2500 and squarkmass >=200 and squarkmass <=4500:
        Values_8TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,8)
    
    output = [Values_7TeV,Values_8TeV]
    return output


def getsgNLL(squarkmass,gluinomass):
    
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
        
    if gluinomass >=200 and gluinomass <=1500 and squarkmass >=200 and squarkmass <=3500:
        Values_7TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,7)
    elif gluinomass >1500 and gluinomass <=2000 and squarkmass >=200 and squarkmass <=2000:
        Values_7TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,7)
        
    if gluinomass >=200 and gluinomass <=2500 and squarkmass >=200 and squarkmass <=4500:
        Values_8TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,8)
    
    output = [Values_7TeV,Values_8TeV]
    return output


def getssNLL(squarkmass,gluinomass):    
    
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
        
    if gluinomass >=200 and gluinomass <=2000 and squarkmass >=200 and squarkmass <=2000:
        Values_7TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,7)
        
    if gluinomass >=200 and gluinomass <=2500 and squarkmass >=200 and squarkmass <=2500:
        Values_8TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,8)
    
    output = [Values_7TeV,Values_8TeV]
    return output


def getsbNLL(squarkmass,gluinomass):    
   
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False}
        
    if gluinomass >=200 and gluinomass <=2000 and squarkmass >=200 and squarkmass <=2000:
        Values_7TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,7)
        
    if gluinomass >=200 and gluinomass <=2500 and squarkmass >=200 and squarkmass <=2500:
        Values_8TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,8)

    output = [Values_7TeV,Values_8TeV]
    return output


def getstopNLL(stopmass, gluinomass):    
    
    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False} 
    
    if stopmass >=200 and stopmass <=2500 and stopmass <=gluinomass:
        Values_8TeV = getNLLfast('st',pdf,stopmass,gluinomass,8)
    
    if stopmass >=100 and stopmass <=1000:
        Values_7TeV = getNLLfast('st',pdf,stopmass,gluinomass,7)

    output = [Values_7TeV,Values_8TeV]
    return output

# main method, call this from the xsec routine with pdgid1, pdgid2 in order to get the NLL results. The output will be xecs =False and masses = 0 if the combination of ids don't have any results at NLL. 
def getNLLresult(pdgid1,pdgid2,inputfile):
    
    readfile=pyslha.readSLHAFile(inputfile)
    gluinomass = abs(readfile[1][1000021].mass)
    squarkmass = (abs(readfile[1][1000001].mass)+abs(readfile[1][2000001].mass)+abs(readfile[1][1000002].mass)+abs(readfile[1][2000002].mass)+abs(readfile[1][1000003].mass)+abs(readfile[1][2000003].mass)+abs(readfile[1][1000004].mass)+abs(readfile[1][2000004].mass))/8
#    print "gluinomass  =", gluinomass
#    print "squarkmass  =", squarkmass

    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False} 
    output=[Values_7TeV,Values_8TeV]

    
    if(pdgid1 < 0):
        for id1 in squark:
            if abs(pdgid1) == id1:
                for id2 in squark:
                    if pdgid2 == id2:
                        output=getsbNLL(squarkmass,gluinomass)

    for id1 in squark:
        if pdgid1 == id1:
            for id2 in squark:
                if pdgid2 == id2:
                    output=getssNLL(squarkmass,gluinomass)
                            
    if pdgid1 == pdgid2 == 1000021:
        output=getggNLL(squarkmass,gluinomass)
                            
    for id1 in squark:
        if pdgid1 == id1:
            if pdgid2 == 1000021:
                output=getsgNLL(squarkmass,gluinomass)

    
    if pdgid1 == pdgid2 == 1000005:
        sbotmass = abs(readfile[1][1000005].mass)
        output=getstopNLL(sbotmass,gluinomass)

    if pdgid1 == pdgid2 == 2000005:
        sbotmass = abs(readfile[1][2000005].mass)
        output=getstopNLL(sbotmass,gluinomass)
        
    if pdgid1 == pdgid2 == 1000006:
        stopmass = abs(readfile[1][1000006].mass)
        output=getstopNLL(stopmass,gluinomass)

    if pdgid1 == pdgid2 == 2000006:
        stopmass = abs(readfile[1][2000006].mass)
        output=getstopNLL(stopmass,gluinomass)

#    print "pdgid1  ==", pdgid1
#    print "pdgid2  ==", pdgid2
#    print "output =",output
    return output


squark = [1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004]

if __name__ == "__main__":
    getNLLresult(pdgid1=1000006,pdgid2=1000006)
    getNLLresult(pdgid1=1000001,pdgid2=1000003)
    getNLLresult(pdgid1=-1000001,pdgid2=1000003)
    getNLLresult(pdgid1=1000006,pdgid2=1000005)
    getNLLresult(pdgid1=1000005,pdgid2=1000005)
    getNLLresult(pdgid1=1000002,pdgid2=-1000002)
