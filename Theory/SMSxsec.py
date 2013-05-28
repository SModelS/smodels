import NLLxsec

#Runs pythia at 7 and 8 TeV and compute weights for each production
#process. Returns a dictionary with weights at 7 and 8 TeV and the
#event file to be used for SMS decomposition.
#If runpythia = False, read old fort_8,7.68 files and assumes 
# total xsec = 1 (only useful for fast debugging)
def pytinit(nevts,slhafile,rpythia = True, donlo = True):

  import shutil, LHEReader
  from Tools.PhysicsUnits import addunit
  
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


#Get 8 TeV event decomposition:
  n8_evts={}
  reader = LHEReader.LHEReader("./data/fort_8.68",nevts)
  for iev in range(nevts):
    Event = reader.next()
    mompdg = Event.getMom()    
    getprodcs(mompdg[0], mompdg[1], n8_evts)

  weight_8 = {}
  sigma_8 = {}
  for key in n8_evts.keys():
    weight_8.update({key : total_cs8/nevts})  #Weight for one event
    sigma_8.update({key : n8_evts[key]*total_cs8/nevts})   #Production cross-section

#Get 7 TeV event decomposition:
  n7_evts={}
  reader = LHEReader.LHEReader("./data/fort_7.68",nevts)
  for iev in range(nevts):
    Event = reader.next()
    mompdg = Event.getMom()
    getprodcs(mompdg[0], mompdg[1], n7_evts)

  weight_7 = {}
  sigma_7 = {}
  for key in n7_evts.keys():
    if n8_evts.has_key(key):
      weight_7.update({key : n7_evts[key]*total_cs7/(nevts*n8_evts[key])}) #Weight for one event
      sigma_7.update({key : n7_evts[key]*total_cs7/nevts})  #Production cross-section

      
#Make sure both dictionaries have the same number of entries
  for key in sigma_7.keys() + sigma_8.keys():
    if not sigma_7.has_key(key): 
      sigma_7.update({key : addunit(0.,'fb')})
      weight_7.update({key : addunit(0.,'fb')})
    if not sigma_8.has_key(key): 
      sigma_8.update({key : addunit(0.,'fb')})
      weight_8.update({key : addunit(0.,'fb')})

#Get NLO cross-sections from NLLfast:
  if donlo:
    sigma_8NLO = {}
    sigma_7NLO = {}
    weight_8NLO = {}
    weight_7NLO = {}

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

      LO7 = weight_7[key]
      LO8 = weight_8[key]
      NLO7 = LO7*k7
      NLO8 = LO8*k8
      weight_7NLO.update({key : NLO7})
      weight_8NLO.update({key : NLO8}) 

      LO7 = sigma_7[key]
      LO8 = sigma_8[key]
      NLO7 = LO7*k7
      NLO8 = LO8*k8
      sigma_7NLO.update({key : NLO7})
      sigma_8NLO.update({key : NLO8})





#Weight and production cross-section dictionaries
  Wdic = {'7 TeV (LO)':weight_7, '8 TeV (LO)':weight_8}
  Xsecdic = {'7 TeV (LO)':sigma_7, '8 TeV (LO)':sigma_8}
  if donlo:
    Wdic.update({'7 TeV (NLO)':weight_7NLO, '8 TeV (NLO)':weight_8NLO})
    Xsecdic.update({'7 TeV (NLO)':sigma_7NLO, '8 TeV (NLO)':sigma_8NLO})

#LHE event file
  lhefile = "./data/fort_8.68"
  
#Weight center of mass energies dictionary
  CMdic = {'7 TeV (LO)': addunit(7.,'TeV'), '8 TeV (LO)': addunit(8.,'TeV')}
  if donlo:
    CMdic.update({'7 TeV (NLO)': addunit(7.,'TeV'), '8 TeV (NLO)': addunit(8.,'TeV')})
  

  return {"Wdic" : Wdic, "lhefile" : lhefile, "Xsecdic" : Xsecdic, "CMdic" : CMdic}  
  
  

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
#        print "  `-- masterkey=",masterkey
  return { "xsecfb": xsecfb }

