def GetPlotLimit(inmass,plot,Analysis,complain = False):
    """ Get upper limit on sigma*BR for a specific array of masses from plot
        inmass: array of masses in SModelS graph?
        plot: ['T8ChiSlep', ['ATLAS_CONF_2013_007']]
        Analysis: SMSDataObjects.EAnalysis
        FIXME in the end this should become no more than a convenience function
        to facilitate looking up results for Analysis objects. 
        All the algorithmic code should go into SMSResults plus helper classes.
    """
    import copy
    import SMSResults
    from Tools.PhysicsUnits import rmvunit

    massarray = [copy.deepcopy(inmass[0]),copy.deepcopy(inmass[1])]

#Run label:
    run = Analysis.run
    if run == "": run = None   #If run has not been defined, use latest run

#Skip empty plots:
    if len(plot) != 2:
      if complain: print "[LimitGetter.py] length of plot != 2"
      if len(plot) != 0:
        print "[LimitGetter.py] warning: len of plot is neither 2 nor 0"
      return False

#Skip empty mass arrays:
    if len(massarray) < 1: 
      if complain: print "[LimitGetter.py] length of massarray < 1"
      return False

    limits = []

#Get mother and LSP masses:
    mLSP = [(massarray[0]).pop(),(massarray[1]).pop()]
    if mLSP[0] != mLSP[1]:
        if complain: print "GetPlotLimit: Different LSP masses"    #For now only allow for equal LSP masses
        return None
    mMom = [massarray[0][0],massarray[1][0]]
    if mMom[0] != mMom[1]:
        if complain: print "GetPlotLimit: Different mother masses"   #For now only allow for equal mother masses
        return None
    mx = mMom[0]
    my = mLSP[0]

    mx = rmvunit(mx,'GeV')
    my = rmvunit(my,'GeV')
#Get intermediate masses and corresponding x values:
    massI = []

    for ib in range(len(massarray)):
        if len(massarray[ib]) > 2:   #(LSP mass has been removed already)
            if complain: print 'GetPlotLimit: Multi-dimensional fit not available, cannot yet deal with multi-cascade decays.'
            return False

        for imass in range(1,len(massarray[ib])):
            mI = massarray[ib][imass]
            mI = rmvunit(mI,'GeV')
	    massI.append(mI)	
#            x = float((mI-my))/(mx-my)
#           massI.append([mI, x])

    if len(massI) == 2:
#        if massI[0][0] != massI[1][0]:
	if massI[0] != massI[1]:
            if complain: print "GetPlotLimit: Different intermediate masses"   #For now only allow for equal intermediate masses
            return None
        else:
            massI = massI.pop()   #For the case of equal masses, keep just one


    CMSlabel = plot[0]        #CMS-type label
    CMSanalyses = plot[1]     #CMS list of analyses
#Now loop over analyses/results and get limits:
    for analyses in CMSanalyses:

#        if not SMSResults.exists(analyses, CMSlabel, run):
#            limits.append([analyses, 'no histogram'])
#            continue

#        xexp = [SMSResults.getx(analyses, CMSlabel, run)]   #experimenal X values
#        if not xexp[0]: xexp = []

#Sanity check
#        if len(massI) > 0 and len(xexp) == 0:
#            if complain: print 'GetPlotLimit: Number of intermediate masses do not match histogram in',analyses,' for',CMSlabel
#            limits.append([analyses, False])
#            continue

#        if len(massI) == 0:  #No need for interpolation
#            limits.append([analyses,SMSResults.getUpperLimit(analyses,CMSlabel,mx, my, run)])

#        else:
#            x = []    #x values
#            y = []    #sigma limits
#            for k in xexp[0]:
#                if k == '050':
#                    ylim = SMSResults.getUpperLimit(analyses,CMSlabel,mx, my)
#                else:
#                    ylim = SMSResults.getUpperLimit(analyses,CMSlabel+k,mx, my)
#
#                if type(ylim) != type(addunit(1.,'fb')): continue  #Skip errors
#
#                x.append(float(k)/100.)
#                y.append(rmvunit(ylim,'fb'))

#Interpolation checks:
#            if len(x) <= 1:
#                if complain: print  'GetPlotLimit: one interpolation point or less in',analyses,' for',CMSlabel
#                limits.append([analyses, None])
#                continue

#           p = numpy.polyfit(x, y, len(y)-1)
#           limits.append([analyses,addunit(float(numpy.polyval(p, massI[1])),'fb') ])
	masslist=[mx]
	if massI: masslist.append(massI)
	masslist.append(my)
	limits.append([analyses,SMSResults.getSmartUpperLimit(analyses,CMSlabel,masslist,run,debug=complain)])

    return limits
