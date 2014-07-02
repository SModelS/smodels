from smodels.theory.printer import Printer


class ExptResults:
    """
    A class to store all relevant information for one result
    """
    def __init__(self, aname, topo, sqrts, cond, tval, exptval, topo_description, mmother, mlsp):
        self.aname = aname
        self.topo = topo
        self.sqrts = sqrts
        self.cond = cond
        self.tval = tval
        self.exptval = exptval
        self.rval = tval/exptval
        self.topo_description = topo_description
        self.mmother = mmother
        self.mlsp = mlsp

class ResultList(Printer):
    """
    Class that collects ExptResults objects and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresult = [], bestresultonly = None, describeTopo = None):
        self.outputarray = outputarray
        self.bestresult = bestresult
        self.bestresultonly = bestresultonly
        self.describeTopo = describeTopo

    def addResult(self,res):
        self.outputarray.append(res)
        return

    def findBest(self):
        best = None
        for res in self.outputarray:
            if not best or best.rval<res.rval:
                best = res
        if best: self.bestresult = [best]
        return

    def formatData(self):
        """
        to access printout format
        """
        return self.formatResultsData()

class OutputStatus(Printer):
    """
    Object that holds all status information and has a predefined printout 
    """
    def __init__(self, status, slhastatus, warnings):
        self.status = status
        self.slhastatus = slhastatus
        self.warnings = warnings
        self.statusStrings = {-1: "#could not run the decomposition",
                              -3: "#no cross sections above sigmacut found",
                              -2: "#bad input slha, did not run decomposition",
                               0: "#no matching experimental results",
                               1: "#decomposition was successful"}

    def formatData(self):
        """
        to access printout format
        """
        return self.formatStatusData()

class InputParameters():
    """
    Object holding all input parameters, __init__ sets default values, then use setFromFile to change parameters according to textfile
    """
    def __init__(self):
        self.doSLHAdec =True
        self.addMissingXsecs = False
        self.addnlo = False
        self.addnll = False
        self.nevts = 50000
        self.doInvisible = True
        self.doCompress = True
        self.sigmacut = 0.03
        self.minmassgap = 5.
        self.maxcond = 0.2
        self.printGtop = True
        self.printThEl = None
        self.evaluateResults = True
        self.printResults = True
        self.expandedSummary = True
        self.analyses = None
        self.topologies = None
        self.describeTopo = True
        self.printAnaEl = False
        self.parameters = ['printGtop', 'minmassgap','printResults', 'doCompress', 'describeTopo', 'topologies', 'doInvisible', 'addMissingXsecs', 'maxcond', 'expandedSummary', 'doSLHAdec', 'evaluateResults', 'printThEl', 'printAnaEl', 'nevts', 'analyses', 'sigmacut', 'addnlo', 'adnll']


    def setFromFile(self, filename = "parameters.in"):
        """
        will update input parameters according to text file filename
        """
        import logging
        logger = logging.getLogger(__name__)
        io_dict = {}
        for l in open(filename):
            l = l.replace(" ", "")
            if l.startswith("#") or l == "\n": continue
            l = l.replace("\n","").split("#")[0]
            if len(l.split("=")) == 1:
                logger.error("Input file %s is damaged." %filename)
                return None
            k, v = l.split("=")[0], l.split("=")[1]
            if v == "True":
                v = True
            elif v == "False" or v == "None" or v == "all":
                v = None
            elif "[" in v or "]" in v:
                v = v.replace("[", "").replace("]", "")
                v = v.split(",")
            else:
                v = float(v)
            io_dict[k] = v
        checkNotSet = []
        checkNotAvailable = []
        for i in io_dict.keys():  
            if not i in self.parameters:
                checkNotAvailable.append(i)
                del io_dict[i]
        for j in self.parameters:
            if not j in io_dict: checkNotSet.append(j)
        if checkNotSet: logger.warning("%s not set in the input file, use default value instead" % str(checkNotSet).replace('[','').replace(']','').replace("'",""))
        if checkNotAvailable: logger.warning("%s not valid input parameter" % str(checkNotAvailable).replace('[','').replace(']','').replace("'",""))

        self.__dict__.update(io_dict)
        return True

class MissingTopo():
    """
    Object to describe one missing topology result
    """
    def __init__(self, topo, weight):
        self.topo = topo
        self.weight = weight

class MissingTopoList(Printer):
    """
    Object to find and collect MissingTopo objects, plus printout functionality
    """
    def __init__(self):
        self.topos = []
        
    def formatData(self):
        return self.formatMissingData()

    def addToTopos(self, el):
        name = self.orderbranches(self.generalName(el.__str__()))
        for topo in self.topos:
            if name == topo.topo: #FIXME need to give correct format of el, plus need general name function!
                topo.weight.__add__(el.weight)
                return
        self.topos.append(MissingTopo(name, el.weight))
        return

    def generalName(self,instr):
        from smodels.theory.particleNames import ptcDic
        exch = ["W", "l", "t", "ta"]
        for pn in exch:
            for on in ptcDic[pn]: instr = instr.replace(on, pn)
        return instr

    def orderbranches(self,instr):
        from smodels.theory.element import Element
        li = Element(instr).getParticles()
        li.sort()
        return str(li).replace("'","").replace(" ","")

    def findMissingTopos(self, smstoplist, listOfAnalyses, sigmacut, minmassgap):
        for el in smstoplist:
            for sel in el.elementList:
                 if sel.compressElement(True, True, minmassgap): continue
                 covered = None
                 for ana in listOfAnalyses:
                     if not ana.getEfficiencyFor(sel) == 0: covered = True
                 if not covered: self.addToTopos(sel)
        return
        


