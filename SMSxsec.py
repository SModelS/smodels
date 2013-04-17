#Runs pythia at 7 and 8 TeV and compute weights for each production
#process. Returns a dictionary with weights at 7 and 8 TeV and the
#event file to be used for SMS decomposition.
#If runpythia = False, read old fort_8,7.68 files and assumes 
# total xsec = 1 (only useful for fast debugging)
def pytinit(nevts,slhafile,rpythia = True, donlo = True):

    import SMSmethods, shutil
    from SMSHelpers import addunit
    
    if rpythia:
#run pythia at 8 TeV:
        D = runpythia(slhafile,nevts,8)
        shutil.copy2("./data/fort.68", "./data/fort_8.68")
        total_cs8 = addunit(D["xsecfb"],'fb')
    
#run pythia at 7 TeV:
        D = runpythia(slhafile,nevts,7)
        shutil.copy2("./data/fort.68", "./data/fort_7.68")
        total_cs7 = addunit(D["xsecfb"],'fb')
    else:
        total_cs8 = addunit(1.,'fb')
        total_cs7 = addunit(1.,'fb')

    sigma_8={}
    sigma_7={}

#Get 8 TeV event decomposition:
    inputfile=open("./data/fort_8.68")
    for iev in range(nevts):
        particles = SMSmethods.getNextEvent(inputfile)
        mompdg=SMSmethods.getMom(particles)
        
        getprodcs(mompdg[0], mompdg[1], sigma_8)


    n8_evts = sigma_8.copy()
    for key in sigma_8.keys(): 
        sigma_8[key]=total_cs8/nevts
          

#Get 7 TeV event decomposition:
    inputfile=open("./data/fort_7.68")
    for iev in range(nevts):
        particles = SMSmethods.getNextEvent(inputfile)
        mompdg=SMSmethods.getMom(particles)
        getprodcs(mompdg[0], mompdg[1], sigma_7)

    n7_evts = sigma_7.copy()
    for key in sigma_7.keys(): 
        if not n8_evts.has_key(key):    #If process was not present at 8 TeV, set to zero weight
            sigma_7[key]=addunit(0.,'fb')
        else: 
            sigma_7[key]=(float(sigma_7[key])*total_cs7)/(nevts*n8_evts[key])

#Make sure both dictionaries have the same number of entries
    for key in sigma_7.keys() + sigma_8.keys():
        if not sigma_7.has_key(key): sigma_7.update({key : addunit(0.,'fb')})
        if not sigma_8.has_key(key): sigma_8.update({key : addunit(0.,'fb')})

#Get NLO cross-sections from NLLfast:
    if donlo:
        import NLLxsec,sys
        sigma_8NLO = {}
        sigma_7NLO = {}
        sumlo = 0.
        lonllfast = 0.

        for key in sigma_8.keys():

            nllres = NLLxsec.getNLLresult(key[0],key[1],slhafile)
            k7 = 1.
            k8 = 1.
     
            if nllres[0]['K_NLL_7TeV'] and nllres[0]['K_NLO_7TeV']:
                k7 = nllres[0]['K_NLL_7TeV']*nllres[0]['K_NLO_7TeV']
            elif nllres[0]['K_NLO_7TeV']:
                k7 = nllres[0]['K_NLO_7TeV']
            if nllres[1]['K_NLL_8TeV'] and nllres[1]['K_NLO_8TeV']:
                k8 = nllres[1]['K_NLL_8TeV']*nllres[1]['K_NLO_8TeV']
            elif nllres[1]['K_NLO_8TeV']:
                k8 = nllres[1]['K_NLO_8TeV']

            LO7 = sigma_7[key]
            LO8 = sigma_8[key]
            NLO7 = LO7*k7
            NLO8 = LO8*k8
            sigma_7NLO.update({key : NLO7})
            sigma_8NLO.update({key : NLO8}) 





#Weight dictionary
    if donlo:
        Wdic = {'7 TeV (LO)':sigma_7, '8 TeV (LO)':sigma_8, '7 TeV (NLO)':sigma_7NLO, '8 TeV (NLO)':sigma_8NLO}
    else:
        Wdic = {'7 TeV (LO)':sigma_7, '8 TeV (LO)':sigma_8}
#LHE event file
    lhefile = "./data/fort_8.68"

    return {"Wdic" : Wdic, "lhefile" : lhefile}  
    
    

def getprodcs(pdgm1, pdgm2, sigma):  
    if pdgm1 > pdgm2: 
        pdgm1, pdgm2 = pdgm2 , pdgm1
    newkey=(pdgm1,pdgm2) 
    if not sigma: 
        init={newkey:1}
        sigma.update(init)
        return sigma
    if not sigma.has_key(newkey):
        sigma[newkey]=0
    sigma[newkey]+=1
    
    return sigma

def runpythia(slhafile,nevts = 10000, Sqrts = 7):
    import sys, os

    workdir = os.getcwd()    
    pyslhadir = workdir + "/pyslha-1.4.3"
    sys.path.append(pyslhadir)
    
    #run pythia
    print "running pythia"
    D = runPythiaLHE(nevts,slhafile, sqrts=Sqrts )
    print "done running"
    return D

def runPythiaLHE ( n, slhafile, sqrts=7 ):
  """ run pythia_lhe with n events, at sqrt(s)=sqrts.
      slhafile is inputfile
      datadir is where this all should run,
      installdir is where pythia_lhe is to be found.  """
  import commands
  print "try to run pythia_lhe at sqrt(s)=%d with %d events" % (sqrts,n)
  datadir="data/"
  etcdir="etc/"
  o=commands.getoutput ( "cp %s %s/fort.61" % ( slhafile, datadir ) )
  if len(o)>0:
    print "[SMSxsec.py] runPythiaLHE error",o
  f=open(etcdir+"/external_lhe.template") 
  lines=f.readlines()
  f.close()
  g=open(datadir+"/external_lhe.dat","write")
  for line in lines:
    out=line.replace("NEVENTS",str(n)).replace("SQRTS",str(1000*sqrts))
    g.write ( out )
  g.close()
  o=commands.getoutput ( "cd %s; ../pythia_lhe < external_lhe.dat" % \
     ( datadir ) )
  lines=o.split( "\n" )
  xsecfb=None
  for line in lines:
#    print line
    if line.find("All included subprocesses")>0:
      try:
        xsecfb=float(line[67:78].replace("D","E"))*10**12
      except Exception,e:
        print "[ResultsTables.py] Exception",e,"xsecfb=",line[67:78]
        print "  `-- line=",line
        print "  `-- masterkey=",masterkey
  return { "xsecfb": xsecfb }

