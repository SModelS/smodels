#!/usr/bin/env python

import sys, commands, os, random, numpy
workdir = os.getcwd()    
pyslhadir = workdir + "/pyslha-1.4.3"
sys.path.append(pyslhadir)
import pyslha
from string import Template
from SMSHelpers import addunit

#define global variables
pdf = 'cteq'


# method to get the NLL fast results. This is the shortest way to code, I found. 
def getNLLfast(process = "gg", pdf = 'cteq', squarkmass=0., gluinomass=0., Energy = 8 ):
    from string import Template
    
    Values={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    mass = 0
    
    if Energy==8:
        if process=="st": 
            os.chdir("./nllfast/nllfast-2.1")
            runnll = Template("./nllfast_8TeV $prodproc $in_pdf $prodsq")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass)
        else: 
            os.chdir("./nllfast/nllfast-2.1")
            runnll = Template("./nllfast_8TeV $prodproc $in_pdf $prodsq $prodgl")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass, prodgl=gluinomass)
    if Energy==7:
        if process=="st": 
            os.chdir("./nllfast/nllfast-1.2")
            runnll = Template("./nllfast_7TeV $prodproc $in_pdf $prodsq")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass)
        else:
            os.chdir("./nllfast/nllfast-1.2")
            runnll = Template("./nllfast_7TeV $prodproc $in_pdf $prodsq $prodgl")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodsq=squarkmass, prodgl=gluinomass)
    o=commands.getoutput(s)
    os.chdir("./../../")
    
    
    if process=='st' and o[1]=='T':
        Values={"sqmass":0, "mgmass":0,"LOcs":False,"NLOcs":False,"NLL+NLO":False,"K_NLO":False,"K_NLL":False}
        return Values
    
    if process!='st' and o[1]=='T':
        if Energy == 7: 
            if process=='gg': 
                process='gdcpl'
                mass = gluinomass
            elif process=='ss': 
                process='sdcpl'
                mass = squarkmass
            else: return Values
        if Energy == 8: 
            if process=='gg': 
                process='gdcpl'
                mass = gluinomass
            elif process=='ss': 
                process='sdcpl'
                mass = squarkmass
            else: return Values
        if Energy==8:
            os.chdir("./nllfast/nllfast-2.1")
            runnll = Template("./nllfast_8TeV $prodproc $in_pdf $prodmass")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodmass=mass)
#            print "decoupling regime command is  ", s
        if Energy==7:
            os.chdir("./nllfast/nllfast-1.2")
            runnll = Template("./nllfast_7TeV $prodproc $in_pdf $prodmass")
            s=runnll.substitute(prodproc=process, in_pdf = pdf, prodmass=mass)
#            print "decoupling regime command is  ", s
        o=commands.getoutput(s)
        os.chdir("./../../")

# uncomment this line to see the nll fast output
        #    print "nll fast output is", o
    
    if o[1]=='T': 
        Values={"sqmass":0, "mgmass":0,"LOcs":False,"NLOcs":False,"NLL+NLO":False,"K_NLO":False,"K_NLL":False}
        return Values
    
    
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
            Values={"sqmass_7TeV":0.0, "mgmass_7TeV":lines[27],"LOcs_7TeV":lines[28],"NLOcs_7TeV":lines[29],"NLL+NLO_7TeV":lines[30],"K_NLO_7TeV":lines[37],"K_NLL_7TeV":lines[38]}
        else:
            Values={"sqmass_7TeV":lines[27], "mgmass_7TeV":lines[28],"LOcs_7TeV":lines[29],"NLOcs_7TeV":lines[30],"NLL+NLO_7TeV":lines[31],"K_NLO_7TeV":lines[38],"K_NLL_7TeV":lines[39]}
    if Energy == 8:
        if process == 'st' or process == 'sdcpl':
            Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":0.0,"LOcs_8TeV":lines[28],"NLOcs_8TeV":lines[29],"NLL+NLO_8TeV":lines[30],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
        elif process == 'gdcpl' :
            Values={"sqmass_8TeV":0.0, "mgmass_8TeV":lines[27],"LOcs_8TeV":lines[28],"NLOcs_8TeV":lines[29],"NLL+NLO_8TeV":lines[30],"K_NLO_8TeV":lines[37],"K_NLL_8TeV":lines[38]}
        else:
            Values={"sqmass_8TeV":lines[27], "mgmass_8TeV":lines[28],"LOcs_8TeV":lines[29],"NLOcs_8TeV":lines[30],"NLL+NLO_8TeV":lines[31],"K_NLO_8TeV":lines[38],"K_NLL_8TeV":lines[39]}
    
    
    return Values




# main method, call this from the xsec routine with pdgid1, pdgid2 in order to get the NLL results. The output will be xecs =False and masses = 0 if the combination of ids don't have any results at NLL. 
def getNLLresult(pdgid1,pdgid2,inputfile):
    
    readfile=pyslha.readSLHAFile(inputfile)
    gluinomass = abs(readfile[1][1000021].mass)
    squarkmass = (abs(readfile[1][1000001].mass)+abs(readfile[1][2000001].mass)+abs(readfile[1][1000002].mass)+abs(readfile[1][2000002].mass)+abs(readfile[1][1000003].mass)+abs(readfile[1][2000003].mass)+abs(readfile[1][1000004].mass)+abs(readfile[1][2000004].mass))/8

    Values_7TeV={"sqmass_7TeV":0, "mgmass_7TeV":0,"LOcs_7TeV":False,"NLOcs_7TeV":False,"NLL+NLO_7TeV":False,"K_NLO_7TeV":False,"K_NLL_7TeV":False}
    Values_8TeV={"sqmass_8TeV":0, "mgmass_8TeV":0,"LOcs_8TeV":False,"NLOcs_8TeV":False,"NLL+NLO_8TeV":False,"K_NLO_8TeV":False,"K_NLL_8TeV":False} 
    output=[Values_7TeV,Values_8TeV]

    
    if(pdgid1 < 0):
        for id1 in squark:
            if abs(pdgid1) == id1:
                for id2 in squark:
                    if pdgid2 == id2:
                        squarkmass = (abs(readfile[1][abs(pdgid1)].mass)+abs(readfile[1][pdgid2].mass))/2
                        Values_7TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,7)
                        Values_8TeV = getNLLfast('sb',pdf,squarkmass,gluinomass,8)
                        output = [Values_7TeV,Values_8TeV]

    
    for id1 in squark:
        if pdgid1 == id1:
            for id2 in squark:
                if pdgid2 == id2:
                    squarkmass = (abs(readfile[1][pdgid1].mass)+abs(readfile[1][pdgid2].mass))/2
                    Values_7TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,7)
                    Values_8TeV = getNLLfast('ss',pdf,squarkmass,gluinomass,8)
                    output = [Values_7TeV,Values_8TeV]

                            
    if pdgid1 == pdgid2 == 1000021:
        Values_7TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,7)
        Values_8TeV = getNLLfast('gg',pdf,squarkmass,gluinomass,8)
        output = [Values_7TeV,Values_8TeV]

                            
    for id1 in squark:
        if pdgid1 == id1:
            if pdgid2 == 1000021:
                Values_7TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,7)
                Values_8TeV = getNLLfast('sg',pdf,squarkmass,gluinomass,8)
                output = [Values_7TeV,Values_8TeV]


    
    if pdgid1 == pdgid2 == 1000005:
        squarkmass = abs(readfile[1][1000005].mass)
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8)
        output = [Values_7TeV,Values_8TeV]


    if pdgid1 == pdgid2 == 2000005:
        squarkmass = abs(readfile[1][2000005].mass)
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8)
        output = [Values_7TeV,Values_8TeV]
        
    if pdgid1 == pdgid2 == 1000006:
        squarkmass = abs(readfile[1][1000006].mass)
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8)
        output = [Values_7TeV,Values_8TeV]

    if pdgid1 == pdgid2 == 2000006:
        squarkmass = abs(readfile[1][2000006].mass)
        Values_7TeV = getNLLfast('st',pdf,squarkmass,gluinomass,7)
        Values_8TeV = getNLLfast('st',pdf,squarkmass,gluinomass,8)
        output = [Values_7TeV,Values_8TeV]

    return output


squark = [1000001,2000001,1000002,2000002,1000003,2000003,1000004,2000004]

if __name__ == "__main__":

    print getNLLresult(pdgid1=1000001,pdgid2=1000002,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=1000021,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=1000006,pdgid2=1000021,inputfile="./AndreSLHA/andrePT4_var.slha")
#    print getNLLresult(pdgid1=-1000001,pdgid2=1000001,inputfile="./AndreSLHA/andrePT4_var.slha")
