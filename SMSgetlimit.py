def GetPlotLimit(massarray,plot,Analysis):
    """ Get upper limit on sigma*BR for a specific array of masses from plot
     (I have included the Analysis in the input to get the value of masscomp and in case we want to access some additional info) """
    import SMSResults  
 
#Run label:
    run = Analysis.run
    if run == "": run = None   #If run has not been defined, use latest run

#Skip empty plots:
    if len(plot) != 2: return False
    
#Skip empty mass arrays:
    if len(massarray) < 1: return False        

    limits = []
#Sanity check:    
    if len(massarray) != 2 or len(massarray[0]) !=  len(massarray[1]):  #For now only allow for equal branches
        print "GetPlotLimit: Wrong mass input"
        return False
   
#Get mother and LSP masses:
    mLSP = [massarray[0].pop(),massarray[1].pop()]
    if mLSP[0] != mLSP[1]:
        print "GetPlotLimit: Different LSP masses"    #For now only allow for equal LSP masses
        return False
    mMom = [massarray[0][0],massarray[1][0]]   
    if mMom[0] != mMom[1]:
        print "GetPlotLimit: Different mother masses"   #For now only allow for equal mother masses
        return False
    mx = mMom[0]
    my = mLSP[0]
    
#Get intermediate masses:
    massI = []
    for imass in range(1,len(massarray[0])):
        if massarray[0][imass] != massarray[1][imass]:
            print "GetPlotLimit: Different Intermediate masses"    #For now only allow for equal intermediate masses
            return False
        massI.append(massarray[0][imass])

        
    CMSlabel = plot[0]        #CMS-type label
    CMSanalyses = plot[1]     #CMS list of analyses
#Now loop over analyses/results and get limits:
    for analyses in CMSanalyses:
        limits.append([analyses,SMSResults.getUpperLimit(analyses,CMSlabel,mx, my, run)])

        
    return limits
