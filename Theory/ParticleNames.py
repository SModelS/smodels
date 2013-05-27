""" dictionaries that deal with particle names """

Rodd={ 
1000021 : "gluino", 1000022: "N1", 1000023 : "N2", 1000025 : "N3", 1000035 : "N4", 1000024 : "C1", 1000037 : "C2", 1000039 : "gravitino",  1000001 : "squark", 1000002 : "squark", 1000003 : "squark", 1000004 : "squark", 2000001 : "squark", 2000002 : "squark", 2000003 : "squark", 2000004 : "squark", 1000005 : "sbottom", 2000005 : "sbottom", 1000006 : "stop", 2000006 : "stop", 1000011 : "slepton", 1000013 : "slepton", 1000015 : "stau", 2000011 : "slepton", 2000013 : "slepton", 2000015 : "stau", 1000012 : "sneutrino", 1000014 : "sneutrino", 1000016 : "sneutrino", 2000012 : "sneutrino", 2000014 : "sneutrino", 2000016 : "sneutrino", -1000021 : "gluino", -1000022: "N1", -1000023 : "N2", -1000025 : "N3", -1000035 : "N4", -1000024 : "C1", -1000037 : "C2", -1000039 : "gravitino",  -1000001 : "squark", -1000002 : "squark", -1000003 : "squark", -1000004 : "squark", -2000001 : "squark", -2000002 : "squark", -2000003 : "squark", -2000004 : "squark", -1000005 : "sbottom", -2000005 : "sbottom", -1000006 : "stop", -2000006 : "stop", -1000011 : "slepton", -1000013 : "slepton", -1000015 : "stau", -2000011 : "slepton", -2000013 : "slepton", -2000015 : "stau", -1000012 : "sneutrino", -1000014 : "sneutrino", -1000016 : "sneutrino", -2000012 : "sneutrino", -2000014 : "sneutrino", -2000016 : "sneutrino"  

}

Reven={
 25 : "higgs", -25: "higgs", 35 : "H0", -35 : "H0", 36 : "A0", -36 : "A0", 37 : "H+", -37 : "H-", 23 : "Z", -23 : "Z", 22 : "photon", -22 : "photon", 24 : "W+", -24 : "W-", 16 : "nu", -16 : "nu", 15 : "ta-", -15 : "ta+", 14 : "nu", -14 : "nu", 13 : "mu-", -13 : "mu+", 12 : "nu", -12 : "nu", 11 : "e-", -11 : "e+", 5 : "b", -5 : "b", 6 : "t+", -6 : "t-", 1 : "jet", 2 : "jet", 3 : "jet", 4 : "jet", 21 : "jet", -1 : "jet", -2 : "jet", -3 : "jet", -4 : "jet", -21 : "jet"
 }

PtcDic={ 
"e" : ["e+","e-"], "mu" : ["mu+", "mu-"], "ta" : ["ta+","ta-"], "l+" : ["e+","mu+"],"l-" : ["e-","mu-"],"l" : ["e-","mu-","e+","mu+"], "W" : ["W+","W-"], "t" : ["t+","t-"], "L+" : ["e+","mu+","ta+"], "L-" : ["e-","mu-","ta-"], "L" : ["e+","mu+","ta+","e-","mu-","ta-"]
}

#Converts pdg number to particle name according to the dictionaries Rodd              
# and Reven                                                                           
def ptype(pdg):                                                                       
  p=int(pdg)                                                                          
  if p in Rodd: return Rodd[p]                                                        
  if p in Reven:                                                                      
    return Reven[p]                                                                   
  else:                                                                               
    return False         
