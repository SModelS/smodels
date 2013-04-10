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
   
#Get mother and LSP masses:
    mLSP = [massarray[0].pop(),massarray[1].pop()]
    if mLSP[0] != mLSP[1]:
        print "GetPlotLimit: Different LSP masses"    #For now only allow for equal LSP masses
        return None
    mMom = [massarray[0][0],massarray[1][0]]   
    if mMom[0] != mMom[1]:
        print "GetPlotLimit: Different mother masses"   #For now only allow for equal mother masses
        return None
    mx = mMom[0]
    my = mLSP[0]
    
#Get intermediate masses and corresponding x values:
    massI = []
    
    for ib in range(len(massarray)):
        for imass in range(1,len(massarray[ib])):            
            mI = massarray[ib][imass]
            x = (mI-my)/(mx-my)
            massI.append([mI, x])
            
    if len(massI) > 1:
        print 'GetPlotLimit: Multi-dimensional fit not available'
        return False        

        
    CMSlabel = plot[0]        #CMS-type label
    CMSanalyses = plot[1]     #CMS list of analyses
#Now loop over analyses/results and get limits:
    for analyses in CMSanalyses:        
        xexp = [SMSResults.getx(analyses, CMSlabel, run)]   #experimenal X values 
        if not xexp[0]: xexp = []
            
#Sanity check
        if len(xexp) != len(massI):
            print 'GetPlotLimit: Number of intermediate masses do not match histogram'
            limits.append([analyses, False])
            continue

        if len(massI) == 0:  #No need for interpolation
            limits.append([analyses,SMSResults.getUpperLimit(analyses,CMSlabel,mx, my, run)])

        else:
            x = []    #x values
            y = []    #sigma limits
            for k in xexp[0]:
                if k == '050':
                    ylim = SMSResults.getUpperLimit(analyses,CMSlabel,mx, my)
                else:
                    ylim = SMSResults.getUpperLimit(analyses,CMSlabel+k,mx, my)
                    
                if ylim:       #Do not include "None" results
                    x.append(float(k)/100.)
                    y.append(ylim)
                
#Interpolation checks:
            if len(x) == 1:
                print  'GetPlotLimit: only one interpolation point'
                limits.append([analyses, None])
                continue
            if massI[0][1] < min(x) or massI[0][1] > max(x):
                print  'GetPlotLimit: intermediate mass out of bounds'
                limits.append([analyses, None])
                continue
      
            p = numpy.polyfit(x, y, len(y)-1)
            limit = numpy.polyval(p, massI[0][1])
            limits.append([analyses, limit])
                

        
    return limits
