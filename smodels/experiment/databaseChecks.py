"""
.. module:: databaseChecks
   :synopsis: Holds the classes and methods to perform checks on the database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging, os, sys, re
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.experiment import infoObjects, dataObjects

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)


class TxNameInfoChecker(object):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    def __init__(self,path,txnameInfoObj):
        
        self.testObj = txnameInfoObj
        self.path = path
        self.mandatoryFields = ['constraint','condition']
        
        for field in self.mandatoryFields:
            if not hasattr(self.testObj ,field):
                logger.warning('in file %s\nmissing field %s for txName %s'\
                %(self.path,field, self.testObj.txName))
                
        for attripute in self.testObj.__dict__:
            if not attripute in TxNameInfoChecker.__dict__:
                continue
            value = getattr(self.testObj,attripute)
            if not getattr(self,attripute)(value):
                logger.warning('in file %s\ninvalid entry %s in field %s for txName %s'\
                %(self.path, getattr(self.testObj,attripute), attripute, self.testObj.txName))
                
    def _constrints(self, value):
        """find the constraint strings in a given string, 
        and check the strings
        :param:value: string
        :return: list with constraints as string
        """
        startString = '[[['
        endString = ']]]'
  
        constraints = []
        for i in range(len(value)):
            if value[i:i + len(startString)] == startString:
                start = i
            if value[i:i + len(endString)] == endString:
                end = i + len(endString)
                constraints.append(value[start:end])

        try:
            checkConstraints = [eval(constraint) for constraint in constraints]
            for constraint in checkConstraints:
                for branch in constraint:
                    for vertex in branch:
                        for particle in vertex:
                            if not isinstance(particle,str): return False
        except:
            return False
  
        return constraints

    
    def constraint(self,value):
        
        constraints = self._constrints(value)
        if not constraints: return False
  
        for constraint in constraints:
            value = value.replace(constraint,'1.0')
        try:
            value = eval(value)
            if not isinstance(value,float): return False
        except:
            return False

        return True
        
    def condition(self, value):
  
        if value == 'None': return True
        constraints = self._constrints(value)
        if not constraints: return False

        value = value.replace('~','==')
        for constraint in constraints:
            value = value.replace(constraint,'1.0')
        try:
            value = eval(value)
            if isinstance(value,tuple):
                for supValue in value: 
                    if not isinstance(supValue,bool): return False
            elif not isinstance(value,bool): return False
        except:
            return False
  
        return True
    
    def condition(self,value):
        
        if value == 'None': return True
        constraints = self._constrints(value)
        if not constraints: return False

        for constraint in constraints:
            value = value.replace(constraint,'1.0')
        functions = value.split(';')
        for function in functions:
            function = function.strip()
            if not function[:4] == 'Csim' and not function[:4] == 'Cgtr':
                return False
            if not function[4:5] == '(' or not function[-1:] == ')':
                return False
            arguments = function[5:-1].split(',')
            try:
                arguments = [eval(argument) for argument in arguments]
                for argument in arguments:
                    if not isinstance(argument,float): return False
            except: 
                return False
   
        return True
    
    def branchcondition(self,value):

        if 'equal branches' == value: return True
        return False
    
    def axes(self,value):
        
        values = value.split('_')
        values = [value.strip() for value in values]
        if not value[:2] == 'Eq': return False
        if not value[2:3] == '(' or not value[-1:] == ')':
            return  False
        if not ',' in value[3:-1]: return False
        return True
    
    def txName(self, value):
        
        if not value[:1] == 'T': return False
        return True

    

def globalInfoChecker(object):
    
    
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    
    def __init__(self,path,globalInfoObj):
        
        self.testObj = txnameInfoObj
        self.path = path
        self.mandatoryFields = ['sqrts', 'lumi', 'id', 'url', 'digitaldata']
        
        for field in self.mandatoryFields:
            if not hasattr(self.testObj ,field):
                logger.warning('in file: %s\nmissing field %s'\
                %(self.path,field))
                
        for attripute in self.testObj.__dict__:
            if not attripute in globalInfoChecker.__dict__:
                continue
            value = getattr(self.testObj,attripute)
            if not getattr(self,attripute)(value):
                logger.warning('in file %s\ninvalid entry %s in field %s'\
                %(self.path, getattr(self.testObj,attripute), attripute))
        
    
    def sqrts(self,value):
        
        pass
    
    def lumi(self,lumi):
        
        pass
    
    def id(self,value):
        
        experiments = ['ATLAS','CMS']
        paperNames = ['CONF','SUSY','PAS']
        try:
            values = value.split('-')
            if not value[0] in experiments:
                return False
            if not value[1] in paperNames:
                return False
            if not isinstance(eval(values[-1]),int) \
            or not isinstance(eval(values[-2]),int):
                return False
        except:
            return False
        return True
                
    def url(self,value):
        
        if not value[:8] == 'https://' and not value[:7] == 'http://':
            return False
        return True
    
    def digitaldata(self,value):
        
        if not value == 'False' and not value == 'True':
            return False
        return True
  
    def publication(self,value):
        
        if not value[:8] == 'https://' and not value[:7] == 'http://':
            return False
        return True
    
    def arxiv(self,value):
        
        if not value[:8] == 'https://' and not value[:7] == 'http://':
            return False
        return True
    
    def superseded_by(self,value):
        
        try:
            values = value.split(',')
            values = [value.strip() for value in values]
            for value in values:
                return self.id(value)
        except:
            return False
        
    def supersedes(self, value):
        
        try:
            values = value.split(',')
            values = [value.strip() for value in values]
            for value in values:
                return self.id(value)
        except:
            return False
    
    def checked(self, value):
        
        checkers = ['AL', 'SuK', 'WW', 'UL', 'VM', 'MT', 'SK']
        try:
            values = values.split(',')
            values = [value.strip() for value in values]
            for value in values:
                topo, checker = value.split(':')
                checker = checker.strip()
                if not checker in checkers: return False
        except:
            return False
        return True
    
    def private(self, value):
        
        if not value == 'False' and not value == 'True':
            return False
        return True
    
    def prettyname(self, value):
    
        return True
    
    def contact(self, value):
        
        return True
    
    def validated(self, value):
        
        if not value == 'False' and not value == 'True':
            return False
        return True


def checkInfoFile(infoFileObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    
    pass
    
    
    
def checkDataFile(infoFileObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    
    globalFields = infoFileObj.globalInfo.getattr()
    
 
def checkExpResult(expResultObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """  
    
    pass
    
def checkDataBase(dataBaseObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """       
    
    pass

class Test(object):
    
    def __init__(self):
        self.constraint = "71.*([[['mu+','mu-']],[['l','nu']]] + [[['e+','e-']],[['l','nu']]])"
        self.condition = "Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['ta','nu']],[['b'],['jet','jet']]]);Cgtr([[['b'],['L','nu']],[['b'],['jet','jet']]],3.*[[['b'],['e','nu']],[['b'],['jet','jet']]])"
        self.axes = 'Eq(y,lsp)_Eq(x,mother)'
        self.branchcondition = 'equal branches'
        self.condition = "[[['mu+']],[['mu-']]] > [[['e+']],[['e-']]]"
        self.txName = 'T6WW'