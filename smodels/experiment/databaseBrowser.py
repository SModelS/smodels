"""
.. module:: databaseBrowserNewStructure
   :synopsis: Centralized facility to access smodels-database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>

"""

import logging
#import setPath
import sys
from smodels.experiment.databaseObjects import Database,ExpResult
import numpy, unum


FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)

class Browser(object):
    
    """Browses the database, exits if given path does not point to a valid 
    smodels-database. Browser can be restricted to specified run or experiment. 
    Verbosity can be set to specified level.
    
    :ivar database: Database object holding all the database information
    :ivar browserList: list of experimental results loaded in the browser.
                       Can be used to hold a subset of results in the database.
                       By default all results are loaded. 
    
    """
    def __init__(self, database):
        
        self.browserList = []
        if isinstance(database,str):
            self.database = Database(database)            
        elif isinstance(database,Database):
            self.database = database
        else:
            logger.error("The input must be the database location or a Database object.")
            sys.exit()            
        self.loadAllResults()
        
    def __iter__(self):
        return iter(self.browserList)

    def __getitem__(self, index):
        return self.browserList[index]

    def __setitem__(self, index, expRes):
        if not isinstance(expRes,ExpResult):
            logger.error("Input object must be a ExpResult() object")
            sys.exit()
        else:
            self.browserList[index] = expRes

    def __len__(self):
        return len(self.browserList)

    def __str__(self):
        return str([expRes.globalInfo.getInfo('id') for expRes in self.browserList])
    
        
    def loadAllResults(self):
        """
        Saves all the results from database to the browserList.
        Can be used to restore all results to browserList.
        """
        
        self.browserList = self.database.expResultList[:]
            
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
        
        
        fieldDict = []
        if expResult and isinstance(expResult,ExpResult):
            fieldDict = expResult.__dict__.items()[:]   #Use only the entries for the expResult
        else:
            for expResult in self:
                fieldDict += expResult.__dict__.items()[:]     #Use all entries/expResults
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

    def getULFor(self,expid,txname=None,massarray=None,datasetID=None):
        """
        Get an upper limit for the given experimental id, the txname, and the massarray. 
        Interpolation is done, if necessary.
        :param expid: experimental id (string)
        :param txname: txname (string). ONLY required for upper limit results
        :param massarray: list of masses with units, e.g.
                          [[ 400.*GeV, 100.*GeV],[400.*GeV, 100.*GeV]]
                          ONLY required for upper limit results
        :param datasetID: string defining the dataset id, e.g. ANA5-CUT3. 
                          ONLY required for efficiency map results
        :return: upper limit [fb]
        """
        
        #First select the experimental results matching the id and the result type:
        expres = None
        for expResult in self:
            if expResult.getValuesFor('id') != expid:
                continue
            else:
                if 'upperLimit' in expResult.getValuesFor('dataType'):
                    if not txname or not massarray: continue
                    expres = expResult
                    break
                elif 'efficiencyMap' in expResult.getValuesFor('dataType'):
                    if not datasetID: continue
                    expres = expResult
                    break

        if not expres:
            logger.warning ( "browser could not find %s . For upper limit results \
            txname and massarray must be defined, while for efficiency map results only \
            dataset must be defined" % (expid))
            return None
        
        if 'upperLimit' in expres.getValuesFor('datatype'):
            txnames = expres.getTxNames()
            for tx in txnames:
                if not tx.txName == txname: continue                
                return tx.txnameData.getValueFor(massarray)
            
        elif 'efficiencyMap' in expres.getValuesFor('datatype'):
            for dataset in expres.datasets:
                if dataset.dataInfo.dataid != datasetID:
                    continue
                return dataset.getUpperLimit()

        logger.warning ( "browser could not find upper limit.")
        return None
     
  
    def loadExpResultsWith(self,restrDict = {}):
        """
        Loads the list of the experimental results (pair of InfoFile and DataFile)
        satisfying the restrictions to the browserList.
        The restrictions specified as a dictionary.
        
        :param restrDict: dictionary containing the fields and their allowed values.
                          E.g. {'lumi' : [19.4/fb, 20.3/fb], 'txname' : 'T1',....}
                          The dictionary values can be single entries or a list of values.
                          For the fields not listed, all values are assumed to be allowed.
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
              
        self.browserList = results[:]
        
        if not self.browserList: logger.warning("Zero results loaded.")
   
