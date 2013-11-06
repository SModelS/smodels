#!/usr/bin/env python

"""
.. module:: XSecComputer
    :synopsis: The unit responsible for the computation of reference ("theory") \
      production cross sections
    
.. moduleauthor:: Suchita Kulkarni <suchita.kulkarni@gmail.com>
    
"""
from Tools.PhysicsUnits import addunit, rmvunit


def loFromLHE( lhefile, totalxsec, nevts=None ):
  """ compute LO weights for each production process from LHE file.

  :param lhefile: name of the lhe file
  :type lhefile: str
  :param totalxsec: total cross section, in fb. If the total cross section can be \
    obtained from the lhe file (from the madgraph banner), then set totalxsec=None
    to pick up the lhe cross section.
  :type totalxsec: cross section 
  :param nevts: maximum number of events to read. None means 'read all'.
  :returns: an array of dictionaries: weights, xsecs, nevts
  """
  from Theory import LHEReader
  import  types


  #Get event decomposition:
  n_evts={}
  reader = LHEReader.LHEReader( lhefile, nevts )
  ntrueevents=0
  for event in reader:
    ntrueevents+=1
    mompdg = event.getMom()    
    if not mompdg: continue
    getprodcs(mompdg[0], mompdg[1], n_evts)

  ##print "reader metainfo=",reader.metainfo
  if type(totalxsec)==types.NoneType:
    if not reader.metainfo.has_key ( "totalxsec" ):
      print "[XSecComputer.py] error: totalxsec=None, but no xsec can be picked up from lhe reader"
      return None
    else:
      totalxsec=reader.metainfo["totalxsec"]

  weight, sigma = {}, {}
  # print "[XSecComputer] lhe",lhefile
  for (key,xsec) in n_evts.items():
    weight[key]= totalxsec / float(ntrueevents) # weight for one event
    sigma[key]=xsec * totalxsec / float(ntrueevents) # production cross-section
    #print "[XSecComputer] moms=%s n=%d ntot=%d w=%s" % ( key,xsec,nevts, weight[key] )

  return [ weight, sigma, n_evts ]
    

def compute(nevts,slhafile,rpythia = True, basedir=None,datadir=None, XsecsInfo=None, tmpfiles=False, printLHE=True):
  """ Runs pythia at 7 and 8 TeV and compute weights for each production process. 

  :param basedir: directory base. If None, then the we look for the SMSDecomposition \
    directory in the current working directory. nllfast we expect to be in \
    basedir+"/nllfast", the template config we expect to reside in basedir+"/etc"
  :type basedir: str
  :param runpythia: run pythia. If False, read old fort_8,7.68 files and assumes \
    total xsec = 1 (only useful for fast debugging)
  :param XsecsInfo: optional information about cross-sections to be computed. \
    If not defined and not found in CrossSection.XSectionInfo, use default values.
  :param tmpfile: If True, it will generate temporary lhe files. If False, use standard \
    file name (fort_sqrts.68) and replace existing files
  :param printLHE: If False, will not write the LHE event file to disk. It will keep it in memory and delete it at the end.
  :returns: a dictionary with weights at 7 and 8 TeV and the event file to be used \
     for SMS decomposition.
  
  """

  import shutil
  from Theory import NLLXSec, CrossSection, LHEReader
  import os, sys
  import tempfile
  import copy
  import cStringIO

  if XsecsInfo is None:
    try:
      XsecsInfo = CrossSection.XSectionInfo  #Check if cross-section information has been defined
    except:
      pass
  if not XsecsInfo:
    XsecsInfo = CrossSection.XSecInfoList()   #If not, define default cross-sections
    CrossSection.XSectionInfo = XsecsInfo
    import logging
    log = logging.getLogger(__name__)
    log.warning ( "Cross-section information not found. Using default values" )

  XsecsInfo.sort(reverse=True)  #Sort by sqrts (the highest sqrts at the first entry will be used for LHE decomposition and weight calculation)
  donlo = max([xsec.order for xsec in XsecsInfo.xsecs])   # =0 if only LO cross-sections are required

  if basedir==None:
    basedir=os.getcwd()
    if basedir[-3:]=="bin": basedir=basedir[:-3]
    if basedir[-4:]=="test": basedir=basedir[:-4]
    if basedir[-10:]=="regression": basedir=basedir[:-10]
  if datadir==None:
    datadir=basedir+"/data"
  if not os.path.exists ( datadir ):
    print "[XSecComputer.py] directory",datadir,"does not exist. Please create."
    return None
    
  nllbase=basedir+"/nllfast"
  if donlo and not os.path.isdir ( nllbase ):
    print "[XSecComputer] error: %s does not exist or is not a directory." % nllbase
    sys.exit(0)
  if donlo and not os.path.isfile ( slhafile ):
    print "[XSecComputer] error: %s does not exist or is not a file." % slhafile
    sys.exit(0)
  installdir=basedir
  etcdir="%s/etc/" % basedir
  

  
  #Get sqrts:
  Allsqrts = XsecsInfo.getSqrts()
  #LHE event files
  lhefiles = []  
  for sqrts in Allsqrts:
    if printLHE:
      if tmpfiles:
        lhefiles.append(tempfile.mkstemp(suffix="_"+str(rmvunit(sqrts,'TeV'))+'TeV', dir=datadir)[1])
      else:
        lhefiles.append(datadir+"/fort_"+str(int(rmvunit(sqrts,'TeV')))+".68")
    else:
        lhefiles.append(None)

   
  total_cs = []
  weights = []
  sigmas = []
  n_evts = []
#Compute LO cross-sections:
  for isqrts,sqrts in enumerate(Allsqrts):
    if rpythia:
      D = runPythia( slhafile,nevts,rmvunit(sqrts,'TeV'),datadir=datadir,etcdir=etcdir,installdir=installdir,printLHE=printLHE)
      if lhefiles[isqrts]:
        shutil.copy2("%s/fort.68" % datadir, lhefiles[isqrts])   #Use file on disk as event file (here lhefiles is a string)
      else:
        lhefiles[isqrts] = cStringIO.StringIO(D["LHEevents"])    #Use memory string as event file (here lhefiles is a file unit)
      total_cs.append(addunit(D["xsecfb"],'fb'))
    else:
      total_cs.append(addunit(1.,'fb'))

    weight, sigma, n_evt = loFromLHE ( lhefiles[isqrts], total_cs[isqrts], nevts )
    weights.append(weight)
    sigmas.append(sigma)
    n_evts.append(n_evt)

#Reweight all LO weights by the highest sqrtS event file and remove weights/cross-sections not present in the highest sqrtS event file
  for isqrts,sqrts in enumerate(Allsqrts):    
    for key in n_evts[isqrts].keys():
      if n_evts[0].has_key(key):
        weights[isqrts][key]=n_evts[isqrts][key]*total_cs[isqrts]/(nevts*n_evts[0][key])
      else:
        weights[isqrts].pop(key)
        sigmas[isqrts].pop(key)

  #Make sure all dictionaries have the same number of entries
  allkeys = []
  for sigma in sigmas: allkeys += sigma.keys()
  allkeys = list(set(allkeys))
  for key in allkeys:
    for isqrts,sqrts in enumerate(Allsqrts):
      if not sigmas[isqrts].has_key(key): 
        sigmas[isqrts][key]=addunit(0.,'fb')
        weights[isqrts][key]=addunit(0.,'fb')


  #Compute NLO/NLL cross-sections (if required) and store all required cross-sections in the results dictionaries:
  Wdic = {}
  Xsecdic = {}
  for xsec in XsecsInfo.xsecs:
    iS = 0
    Sigma = {}
    Weight = {}
    for isqrts,sqrts in enumerate(Allsqrts):
      if xsec.sqrts == sqrts: 
        Weight = copy.deepcopy(weights[isqrts])
        Sigma = copy.deepcopy(sigmas[isqrts])        #Get LO cross-section

    if xsec.order > 0:
      for key in Sigma.keys():
        k = 1. 
        nllres = NLLXSec.getNLLresult(key[0],key[1],slhafile,base=nllbase)             
        klabel = str(int(rmvunit(xsec.sqrts,'TeV')))+'TeV'        
        for nll in nllres:
          if nll.has_key('K_NLO_'+klabel) and nll['K_NLO_'+klabel]:  k = k*nll['K_NLO_'+klabel]              #NLO k-factor (or NLL+NLO k-factor if there is no NLL result)
          if xsec.order == 2 and nll.has_key('K_NLL_'+klabel) and nll['K_NLL_'+klabel]: k = k*nll['K_NLL_'+klabel]            #NLL+NLO k-factor (or NLL k-factor if there is no NLO result)

        Weight[key] = Weight[key]*k
        Sigma[key] = Sigma[key]*k
        
    Wdic[xsec.label] = Weight
    Xsecdic[xsec.label] = Sigma

#Save lhe file names associated with each cross-section:
  LHEfiles = {}
  for xsec in XsecsInfo.xsecs:
    for isqrts,sqrts in enumerate(Allsqrts):
      if xsec.sqrts == sqrts:
        if not printLHE and lhefiles[isqrts]:
          lhefiles[isqrts].close()
          lhefiles[isqrts] = None
        LHEfiles[xsec.label] = lhefiles[isqrts]

  
  return CrossSection.CrossSection ( {"Wdic" : Wdic, "lhefile" : lhefiles[0], "lhefiles" : LHEfiles, "Xsecdic" : Xsecdic, "XsecList" : XsecsInfo} )
  
  

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

def runPythia ( slhafile, n, sqrts=7, datadir="./data/", etcdir="./etc/",
                installdir="./", printLHE=True ):
  """ run pythia_lhe with n events, at sqrt(s)=sqrts.

    :param slhafile: inputfile
    :type slhafile: str
    :param datadir: directory where this all should run
    :param etcdir: is where external_lhe.template is to be picked up
    :param installdir: is where pythia_lhe is to be found.
    :param printLHE: choose if LHE event file is written to disk or not. If False, returns the events as a string with key LHEevents
  """
  import commands, os, sys
  # print "[SMSXSecs.runPythia] try to run pythia_lhe at sqrt(s)=%d with %d events" % (sqrts,n)
  o=commands.getoutput ( "cp %s %s/fort.61" % ( slhafile, datadir ) )
  if len(o)>0:
    print "[SMSXSecs.py] runPythia error",o
  f=open(etcdir+"/external_lhe.template") 
  lines=f.readlines()
  f.close()
  g=open(datadir+"/external_lhe.dat","write")
  for line in lines:
    if not printLHE and "MSTP(163)=" in line: line = "MSTP(163)=6\n"  #Switches output to screen, so no file is written to disk
    out=line.replace("NEVENTS",str(n)).replace("SQRTS",str(1000*sqrts))
    g.write ( out )
  g.close()
  if not os.path.isdir( datadir ):
    print "[XSecComputer.py] error: %s does not exist or is not a directory." % datadir
    sys.exit(0)
  executable="%s/pythia_lhe" % installdir
  if not os.path.isfile ( executable ) or not os.access(executable,os.X_OK):
    print "[XSecComputer.py] error: %s does not exist." % executable
    sys.exit(0)
  if not os.path.isfile ( "%s/external_lhe.dat" % datadir ):
    print "[XSecComputer.py] error: %s/external_lhe.dat does not exist." % datadir
    sys.exit(0)
  o=commands.getoutput ( "cd %s; %s < external_lhe.dat" % \
     ( datadir, executable ) )
  try:
    lines = o[:o.find("<LesHouchesEvents")]
    LHEevents = o[o.find("<LesHouchesEvents"):]
  except:
    lines = o
    LHEevents = None
  lines=lines.split( "\n" )
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
  return { "xsecfb": xsecfb, "LHEevents" : LHEevents }

def clean ( datadir ):
  """ simple routine that can help to clean up after having computed everything """
  import os
  for i in os.listdir ( datadir ): 
    try:
      os.unlink ( datadir + "/" + i )
    except Exception,e:
      print "[XSecComputer] error, cannot unlink %s/%s: %s" % ( datadir, i, e )
  try:
    os.rmdir ( datadir )
  except Exception,e:
    print "[XSecComputer] error, cannot rmdir %s: %s" % ( datadir, e )
