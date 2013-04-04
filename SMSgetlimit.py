#Get upper limit on sigma*BR for a specific array of masses from plot
#(I have included the Analysis in the input to get the value of masscomp and in case we want to access some additional info)
def GetPlotLimit(massarray,plot,Analysis):
    import ROOT
    from SMSmethods import EqualMasses, MassAvg
    import SMSHelpers  
 
#Data base directory (will have to be modified/improved later! just an example!!!)
    datadir = SMSHelpers.Base+"./2012"    

#Skip empty plots:
    if len(plot) != 2: return False
    
        

    limits = []
#Sanity check:    
    if len(massarray) != 2 or len(massarray[0]) !=  len(massarray[1]):  #For now only allow for equal branches
        print "GetPlotLimit: Wrong mass input"
        return False
#Get masscomp value for defining similar masses:
    masscomp = Analysis.masscomp
   
#Get mother and LSP masses:
    mLSP = [massarray[0].pop(),massarray[1].pop()]
    if not EqualMasses(mLSP[0],mLSP[1],masscomp):
        print "GetPlotLimit: Different LSP masses"    #For now only allow for equal LSP masses
        return False
    my = MassAvg(mLSP,"harmonic")    
    mMom = [massarray[0][0],massarray[1][0]]   
    if not EqualMasses(mMom[0],mMom[1],masscomp):
        print "GetPlotLimit: Different mother masses"   #For now only allow for equal mother masses
        return False
    mx = MassAvg(mMom,"harmonic")    
#Get intermediate masses:
    massI = []
    for imass in range(1,len(massarray[0])):
        if not EqualMasses(massarray[0][imass],massarray[1][imass],masscomp):
            print "GetPlotLimit: Different Intermediate masses"    #For now only allow for equal intermediate masses
            return False
        massI.append(MassAvg([massarray[0][imass],massarray[1][imass]],"harmonic"))

        
    CMSlabel = plot[0]        #CMS-type label
    CMSanalyses = plot[1]     #CMS list of analyses
#Now loop over analyses/results and get limits:
    for analyses in CMSanalyses:
        rootfilename = datadir + "/" + analyses + "/sms.root"  #Get Analyses filename
	try:
            rootfile = ROOT.TFile(rootfilename)
        except:
            continue
        histname = "limit_" + CMSlabel   #Get histogram name
        hist = rootfile.Get(histname)
        if not hist: continue            
        
        b = hist.FindBin(mx, my)
        limits.append([analyses,hist.GetBinContent(b)])   #Save experimental limit and analyses name

        
    return limits
