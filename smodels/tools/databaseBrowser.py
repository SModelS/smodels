"""
.. module:: databaseBrowser
   :synopsis: Centralized facility to access smodels-database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.experiment.databaseObj import Database,ExpResult
import numpy, unum
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
from smodels.tools.smodelsLogging import logger
from smodels.theory.auxiliaryFunctions import getAttributesFrom,getValuesForObj

#logger.setLevel(level=logging.INFO)

class Browser(object):

    """
    Browses the database, exits if given path does not point to a valid
    smodels-database. Browser can be restricted to specified run or experiment.
    """

    def __init__(self, database, force_txt = False ):
        """
        :param force_txt: If True forces loading the text database.
        :param database: Path to the database or Database object
        """

        self._selectedExpResults = []
        load = None
        if force_txt == True:
            load = "txt"
        if isinstance(database,str):
            if database.endswith(".pcl"):
                load = "pcl"
            self.database = Database(database, load )
        elif isinstance(database,Database):
            self.database = database
        else:
            logger.error("The input must be the database location or a Database object.")
            raise SModelSError()
        self.loadAllResults()

    def __iter__(self):
        return iter(self._selectedExpResults)

    def __getitem__(self, index):
        return self._selectedExpResults[index]

    def __setitem__(self, index, expRes):
        if not isinstance(expRes,ExpResult):
            logger.error("Input object must be a ExpResult() object")
            raise SModelSError()
        else:
            self._selectedExpResults[index] = expRes

    def __len__(self):
        return len(self._selectedExpResults)

    def __str__(self):
        return str([expRes.globalInfo.getInfo('id') for expRes in self._selectedExpResults])


    def loadAllResults(self):
        """
        Saves all the results from database to the _selectedExpResults.
        Can be used to restore all results to _selectedExpResults.
        """

        self._selectedExpResults = self.database.expResultList[:]

    def getValuesFor(self,attribute,expResult=None):
        """
        Returns a list for the possible values appearing in the database
        for the required attribute (sqrts,id,constraint,...).

        :param attribute: name of a field in the database (string).
        :param expResult: if defined, restricts the list to the corresponding expResult.
                          Must be an ExpResult object.
        :return: list of values
        """

        if not expResult:
            return getValuesForObj(self,attribute)
        else:
            return getValuesForObj(expResult,attribute)


    def getAttributes(self,showPrivate=False):
        """
        Checks for all the fields/attributes it contains as well as the
        attributes of its objects if they belong to smodels.experiment.

        :param showPrivate: if True, also returns the protected fields (_field)
        :return: list of field names (strings)
        """

        attributes = getAttributesFrom(self)

        if not showPrivate:
            attributes = list(filter(lambda a: a[0] != '_', attributes))

        return attributes

    def getEfficiencyFor(self,expid,dataset,txname,massarray):
        """
        Get an efficiency for the given experimental id,
        the dataset name, the txname, and the massarray.
        Can only be used for EfficiencyMap-type experimental results.
        Interpolation is done, if necessary.

        :param expid: experimental id (string)
        :param dataset: dataset name (string)
        :param txname: txname (string).
        :param massarray: list of masses with units, e.g.
                          [[ 400.*GeV, 100.*GeV],[400.*GeV, 100.*GeV]]
        :return: efficiency
        """
        #First select the experimental results matching the id and the result type:
        expres = None
        for expResult in self:
            if expResult.globalInfo.id != expid:
                continue
            else:
                if 'efficiencyMap' in [ds.dataInfo.dataType for ds in expResult.datasets]:
                    expres = expResult
                    break

        if not expres:
            logger.warning( "Could not find efficiencyMap result %s."\
                   " getEfficiencyForr can only be\
            used for efficiencyMap results." % (expid))
            return None

        return expres.getEfficiencyFor(txname=txname, mass=massarray, dataset=dataset)

    def getULFor(self,expid,txname,massarray, expected=False ):
        """
        Get an upper limit for the given experimental id, the txname,
        and the massarray.
        Can only be used for UL experimental results.
        Interpolation is done, if necessary.

        :param expid: experimental id (string)
        :param txname: txname (string). ONLY required for upper limit results
        :param massarray: list of masses with units, e.g.
                          [[ 400.*GeV, 100.*GeV],[400.*GeV, 100.*GeV]]
        :param expected: If true, return expected upper limit, otherwise
                         return observed upper limit.
        :return: upper limit [fb]
        """

        #First select the experimental results matching the id and the result type:
        expres = None
        for expResult in self:
            if expResult.globalInfo.id != expid:
                continue
            else:
                if 'upperLimit' in [ds.dataInfo.dataType for ds in expResult.datasets]:
                    expres = expResult
                    break

        if not expres:
            logger.warning( "Could not find UL result %s. getULFor can only be\
            used for upper-limit results." % (expid))
            return None

        txnames = expres.getTxNames()
        for tx in txnames:
            if not tx.txName == txname:
                continue
            return tx.getULFor(massarray,expected)

        logger.warning( "Could not find TxName %s ." % (txname))
        return None

    def getULForSR(self,expid,datasetID):
        """
        Get an upper limit for the given experimental id and dataset (signal region).
        Can only be used for efficiency-map results.
        :param expid: experimental id (string)
        :param datasetID: string defining the dataset id, e.g. ANA5-CUT3.
        :return: upper limit [fb]
        """

        #First select the experimental results matching the id and the result type:
        expres = None
        for expResult in self:
            if expResult.globalInfo.id != expid:
                continue
            else:
                if 'efficiencyMap' in [ds.dataInfo.dataType for ds in expResult.datasets]:
                    expres = expResult
                    break

        if not expres:
            logger.warning ("Could not find efficiency-map result %s . getULForSR can only be\
            used for efficiency-map results." % (expid))
            return None

        for dataset in expres.datasets:
            if dataset.getID() != datasetID:
                continue
            return dataset.getSRUpperLimit()

        logger.warning ( "Could not find dataset %s ." % (datasetID))
        return None


    def selectExpResultsWith(self,**restrDict):
        """
        Loads the list of the experimental results (pair of InfoFile and DataFile)
        satisfying the restrictions to the _selectedExpResults.
        The restrictions specified as a dictionary.

        :param restrDict: selection fields and their allowed values.
                          E.g. lumi = [19.4/fb, 20.3/fb], txName = 'T1',....}
                          The values can be single entries or a list of values.
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
            expAttributes = expRes.getAttributes()
            for tag in restrDict:
                #Check if the restriction tag appears in the experimental result
                if not tag in expAttributes:
                    results.remove(expRes)
                    break
                vals = expRes.getValuesFor(tag)
                #If it does, check if any of the values in expResult match any of the values given
                #as restrictions
                if not isinstance(restrDict[tag],list):
                    rvals = [restrDict[tag]]
                else:
                    rvals = restrDict[tag]
                #If there is a type mismatch, also remove
                try: intersec = numpy.intersect1d(vals,rvals)
                except unum.IncompatibleUnitsError:
                    logger.warning("Incompatible units, skipping result.")
                    results.remove(expRes)
                    break
                if len(intersec) == 0:
                    results.remove(expRes)
                    break

        self._selectedExpResults = results[:]

        if not self._selectedExpResults: logger.warning("Zero results loaded.")


def main(args):
    """
    IPython interface for browsing the Database.
    """

    try:
        import IPython
    except ImportError:
        print("IPython is not installed.")
        print("To use this script, please install ipython.")
        import sys
        sys.exit()


    from smodels.tools import databaseBrowser
    from smodels.tools.physicsUnits import fb, pb, GeV, TeV

    def getHeader ():
        from smodels.installation import installDirectory
        header = ""
        with open( installDirectory()+"smodels/share/BANNER") as f:
            lines=f.readlines()
            for line in lines: header+=line

        header += "\n"
        header += "fb, pb, GeV, TeV defined.\n"
        header +=  "\nBrowser loaded for %s \n" %( args.path_to_database )
        header += "Try 'print(browser)' for the list of available results.\n"
        header += "More examples on how to access the database can be found in the SModelS manual.\n"
        header += "\nType 'exit' to exit this session."

        return header

    browser = databaseBrowser.Browser(args.path_to_database, args.text )
    IPython.embed(header=getHeader())
