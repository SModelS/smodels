

#Runs pythia at 7 and 8 TeV and compute weights for each production
#process. Returns a dictionary with weights at 7 and 8 TeV and the
#event file to be used for SMS decomposition.
#If runpythia = False, read old fort_8,7.68 files and assumes 
# total xsec = 1 (only useful for fast debugging)
def pytinit(nevts,slhafile,rpythia = True):

    import SMSmethods, shutil
    
    if rpythia:
#run pythia at 8 TeV:
        D = runpythia(slhafile,nevts,8)
        shutil.copy2("./data/fort.68", "./data/fort_8.68")
        total_cs8 = D["xsecfb"]   
    
#run pythia at 7 TeV:
        D = runpythia(slhafile,nevts,7)
        shutil.copy2("./data/fort.68", "./data/fort_7.68")
        total_cs7 = D["xsecfb"]
    else:
        total_cs8 = 1.
        total_cs7 = 1.

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
            sigma_7[key]=0
        else: 
            sigma_7[key]=(float(sigma_7[key])*total_cs7)/(nevts*n8_evts[key])

#Make sure both dictionaries have the same number of entries
    for key in sigma_7.keys() + sigma_8.keys():
        if not sigma_7.has_key(key): sigma_7.update({key : 0.})
        if not sigma_8.has_key(key): sigma_8.update({key : 0.})


#Weight dictionary
    Wdic = {'7 TeV':sigma_7, '8 TeV':sigma_8}
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
    import PMSSM



    PMSSM.installdir = workdir
    PMSSM.etcdir= workdir + "/etc"
    PMSSM.logdir= workdir + "/log"

    PMSSM.verbose=True
    PMSSM.basedir = workdir
    PMSSM.datadir = workdir + "/data"
    
#run pythia
    print "running pythia"
    D = PMSSM.runPythiaLHE(nevts,slhafile, " ",sqrts=Sqrts)
    print "done running"
    return D
