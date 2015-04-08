#!/usr/bin/env python


#from smodels.tools.statistics import getUL,getPValue,computeCLInterval
import sys
from smodels.theory.topology import TopologyList
from smodels.theory.element import Element
from smodels.theory.theoryPrediction import TheoryPredictionList
from smodels.experiment.txnameObject import TxName
from smodels.tools.ioObjects import OutputStatus, ResultList
from smodels.tools.missingTopologies import MissingTopoList

obj1 = [OutputStatus(status=[0.,1.], inputFile='bla', parameters='bla2', databaseVersion='bla3', outputfile='bla3')]

printingOrder = [OutputStatus,TopologyList,Element,TheoryPredictionList,ResultList,MissingTopoList]

for obj in printingOrder:
   for iobj,objB in enumerate(obj1):
      if type(objB) == obj:
         print type(objB)



