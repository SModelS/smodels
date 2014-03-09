#Break down the total cross-section (AFTER compression) in tested and not tested components
def getTestComposition(SMSTopList,ListOfAnalyses,maxcond):
  from theory import CrossSection, ClusterTools
  from tools.PhysicsUnits import addunit

  zeroweight = {}
  for xsec in CrossSection.XSectionInfo.xsecs:
    zeroweight[xsec.label] = addunit(0., 'fb')      

  Long = zeroweight
  noAna = zeroweight
  asyTop = zeroweight
  noLim = zeroweight
  noCond = zeroweight
  Tested = zeroweight
  

  for Top in SMSTopList:
    for EEl in Top.ElList:
      isLong = False
      isAsy = False  
      if max(Top.vertnumb) > 3: isLong = True
      elif Top.vertparts[0] != Top.vertparts[1]: isAsy = True
      else:     
        inAna = False
        hasLimit = False
        goodCond = False
        mass = [EEl.B[0].masses,EEl.B[1].masses]
        iana = 0
        while not goodCond and iana < len(ListOfAnalyses):
          Analysis = ListOfAnalyses[iana]
          iana += 1
          if EEl.isInAnalysis(Analysis,igmass=True):   #Check if element is in Analysis
            inAna = True
            gmass = False
            if Analysis.goodmasses.has_key(str(mass)):   #Check if mass has limit
              gmass = Analysis.goodmasses[str(mass)]              
            elif Analysis.goodmasses.has_key(str([mass[1],mass[0]])):
              gmass = Analysis.goodmasses[str([mass[1],mass[0]])]
            
            if not gmass: continue #Mass does not have limit
            hasLimit = True
            for cluster in Analysis.ResultList:
              if not gmass in cluster.allmasses: continue  # Element mass does not belong to cluster
              cluster_maxcond = cluster.getMaxCondition()
              if cluster_maxcond == "N/A": continue  #At least one of the conditions could not be computed
              if cluster_maxcond is None or cluster_maxcond < maxcond:
                goodCond = True   #No conditions or conditions below cut
                break

      if hasLimit and not goodCond:
        print EEl
        print mass,gmass

      if isLong:
        Long = ClusterTools.sumweights([Long,EEl.weight])   #Add up cross-sections of long elements
      elif isAsy:
        asyTop = ClusterTools.sumweights([asyTop,EEl.weight]) #Add up cross-sections of asymmetric elements
      elif not inAna:
        noAna = ClusterTools.sumweights([noAna,EEl.weight]) #Add up cross-sections of elements without analysis
      elif not hasLimit:
        noLim = ClusterTools.sumweights([noLim,EEl.weight]) #Add up cross-sections of elements outside ranges
      elif not goodCond:
        noCond = ClusterTools.sumweights([noCond,EEl.weight]) #Add up cross-sections of elements violating conditions
      else:
        Tested = ClusterTools.sumweights([Tested,EEl.weight]) #Add up cross-sections of elements being tested

  return {'Tested' : Tested, 'Long_Topology' : Long, 'No_Analysis' : noAna, 'No_Limit' : noLim, 'Bad_Conditions' : noCond, 'Asymmetric_Top' : asyTop}

