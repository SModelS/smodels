def GetPlotLimit(massarray,plot,Analysis):
    """ Get upper limit on sigma*BR for a specific array of masses from plot """
    import SMSResults, numpy 
 
#Run label:
    run = Analysis.run
    if run == "": run = None   #If run has not been defined, use latest run

#Skip empty plots:
    if len(plot) != 2: return False
    
#Skip empty mass arrays:
    if len(massarray) < 1: return False        

    limits = []
#Sanity check:    
    if len(massarray) != 2: #2 branches only
## or len(massarray[0]) !=  len(massarray[1]):  #For now only allow for equal branches
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
	for branch in range(len(massarray)):
       	    try:
                mI = massarray[branch][1]
           	x = (mI-my)/(mx-my)
           	massI.append([mI, x])
            except:
		continue

#        if massarray[0][imass] != massarray[1][imass]:
#            print "GetPlotLimit: Different Intermediate masses"    #For now only allow for equal intermediate masses
#            return False
#        massI.append(massarray[0][imass])

        
    CMSlabel = plot[0]        #CMS-type label
    CMSanalyses = plot[1]     #CMS list of analyses
#Now loop over analyses/results and get limits:
    for analyses in CMSanalyses:
#        limits.append([analyses,SMSResults.getUpperLimit(analyses,CMSlabel,mx, my, run)])
        xdict = SMSResults.getx(analyses)
        if not xdict or not xdict.has_key(CMSlabel):
             limit = SMSResults.getUpperLimit(analyses,CMSlabel,mx, my)
             limits.append([analyses, limit])
             if massI:
                print '\033[93m'+"Could not find histograms for different x for",analyses,'\033[0m'
             continue
        if xdict.has_key(CMSlabel):
             if massI[0][1]<0.05 or massI[0][1]>0.95:
              limits.append([analyses, "x value out of range"])
              continue
             y = []
             x = []
             for k in xdict[CMSlabel]:
                if k == '050':
                     y.append(SMSResults.getUpperLimit(analyses,CMSlabel,mx, my))
                else:
                     y.append(SMSResults.getUpperLimit(analyses,CMSlabel+k,mx, my))
                x.append(float(k)/100)
             if len(massI)>1:
                print "More than one intermediate mass! Use the one that was found first!"
             p = numpy.polyfit(x, y, len(y)-1)
             limit = numpy.polyval(p, massI[0][1])
             limits.append([analyses, limit])
        
    return limits
