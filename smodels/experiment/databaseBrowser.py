"""
.. module:: databaseBrowserNewStructure
   :synopsis: Centralized facility to access smodels-database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>

"""

import logging
#import setPath
import sys
from smodels.experiment.databaseObjects import DataBase,ExpResult
import numpy, unum


FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)

class Browser(object):
    
    """Browses the database, exits if given path does not point to a valid 
    smodels-database. Browser can be restricted to specified run or experiment. 
    Verbosity can be set to specified level.
    
    :ivar database: DataBase object holding all the database information
    
    """
    def __init__(self, database):
        if isinstance(database,str):
            self.database = DataBase(database)
        elif isinstance(database,DataBase):
            self.database = database
        else:
            logger.error("The input must be the database location or a DataBase object.")
            sys.exit()
            
    def getValuesFor(self,attribute=None,expResult=None):
        """
        Returns a list for the possible values appearing in the database
        for the required attribute (sqrts,id,constraint,...).
        
        :param attribute: name of a field in the database (string). If not defined
                          it will return a dictionary with all fields and their respective
                          values
        :param expResult: if defined, restricts the list to the corresponding expResult.
                          Must be an ExpResult object.
        :return: list of values
        """
        
        
        if expResult and isinstance(expResult,ExpResult):
            fieldDict = expResult.__dict__.items()[:]   #Use only the entries for the expResult
        else: fieldDict = self.__dict__.items()[:]     #Use all entries/expResults
        valuesDict = {}
        while fieldDict:
            for field,value in fieldDict[:]:
                if not '<smodels.experiment' in str(value):
                    if not field in valuesDict: valuesDict[field] = [value]
                    else: valuesDict[field].append(value)              
                else:
                    if isinstance(value,list):
                        for entry in value: fieldDict += entry.__dict__.items()[:]
                    else: fieldDict += value.__dict__.items()[:]
                fieldDict.remove((field,value))                

        #Try to keep only the set of unique values
        for key,val in valuesDict.items():
            try: valuesDict[key] = list(set(val))
            except: pass
        if not attribute: return valuesDict
        elif not attribute in valuesDict:
            logger.warning("Could not find field %s in database" % attribute)
            return False
        else:
            return valuesDict[attribute]
            
            
    def getAttributes(self,showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.
        
        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)
        """
        
        fields = self.getValuesFor().keys()
        fields = list(set(fields))
        
        if not showPrivate:
            for field in fields[:]:
                if "_" == field[0]: fields.remove(field)
               
        return fields
     
  
    def getExpResultsWith(self,restrDict = {}):
        """Returns a list of all experimental results (pair of InfoFile and DataFile)
        with the fields and field values given as a restriction dictionary.
        """
          
        #First check all the selected fields exist and build the corresponding
        #restriction dictionary
        rDict = {}
        allfields = self.getAttributes()
        for tag,value in restrDict.items():
            if tag in allfields: rDict[tag] = value
            else: logger.warning("Field/attribute %s not found (will be ignored)." % tag)
      
        results = self.database.expResultList[:]
        for expRes in results[:]:
            values = self.getValuesFor(attribute=None, expResult=expRes)           
            for tag in restrDict:
                #Check if the restriction tag appears in the experimental result
                if not tag in values:
                    results.remove(expRes)
                    break
                vals = values[tag]                
                #If it does, check if any of the values in expResult match any of the values given
                #as restrictions
                if not isinstance(restrDict[tag],list): rvals = [restrDict[tag]]
                else: rvals = restrDict[tag]                
                #If there is a type mismatch, also remove
                try: intersec = numpy.intersect1d(vals,rvals)
                except unum.IncompatibleUnitsError:
                    logger.warning("Incompatible units, skipping result.")
                    results.remove(expRes)
                    break                    
                if len(intersec) == 0:
                    results.remove(expRes)
                    break                        

              
        return results