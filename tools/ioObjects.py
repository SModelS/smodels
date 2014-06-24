from smodels.theory.printer import Printer


class ExptResults:
    """
    A class to store all relevant information for one result
    """
    def __init__(self, aname, topo, sqrts, cond, tval, exptval, rval, topo_description, mmother, mlsp):
        self.aname = aname
        self.topo = topo
        self.sqrts = sqrts
        self.cond = cond
        self.tval = tval
        self.exptval = exptval
        self.rval = rval
        self.topo_description = topo_description
        self.mmother = mmother
        self.mlsp = mlsp

class ResultList(Printer):
    """
    Class that collects ExptResults objects and has a predefined printout
    """
    def __init__(self, outputarray = [], bestresult = [], bestresultonly = None):
        self.outputarray = outputarry
        self.bestresult = bestresult
        self.bestresultonly = bestresultonly

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
        self.describe_topo = True
        self.printAnaEl = False
        self.parameters = ['printGtop', 'minmassgap','printResults', 'doCompress', 'describe_topo', 'topologies', 'doInvisible', 'addMissingXsecs', 'maxcond', 'expandedSummary', 'doSLHAdec', 'evaluateResults', 'printThEl', 'printAnaEl', 'nevts', 'analyses', 'sigmacut']


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
    def __init__(self, topo, weight, sqrts):
        self.topo = topo
        self.weight = weight
        self.sqrts = sqrts

class MissingTopoList(Printer):
    """
    Object to find and collect MissingTopo objects, plus printout functionality
    """
    def __init__(self):
        self.topos = []
        
    def formatData(self):
        return self.formatMissingData()

    def findMissingTopos(self):
        return None #FIXME adapt IO to develop here

