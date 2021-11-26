#!/usr/bin/env python3

"""
.. module:: signalRegionsCombiner
   :synopsis: a module that deals with combining signal regions from 
              different analyses into a single fake result with an appropriate
              covariance matrix for the corresponding simplified likelihood.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>
.. moduleauthor:: Jamie Yellen <j.yellen.1@research.gla.ac.uk>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.smodelsLogging import logger
from smodels.experiment.datasetObj import DataSet
from smodels.experiment.expResultObj import ExpResult
from smodels.experiment.infoObj import Info
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
import numpy as np

class SignalRegionsCombiner():
    """
    Facility used to combine signal regions from different analyes.
    Creates a fake result with a (usually diagonal) covariance matrix
    Written originally by Jamie Yellen, ported into SModelS proper by WW.
    """
    def __init__ (self):
        self.fakeResult = None


    def selectDatasetsFrom ( self, objList, filter_labels = None ):
        """
        creates a list of datasets using the filter dictionary. The input list can be
        a list of Dataset, ExptResult or TheoryPrediction objects.

        :param objList: a list of Dataset, ExptResult or TheoryPrediction objects.
        :filter_labels: A list of strings with labels for filtering the objList. The labels should be
                      in the format ['datasetID1', 'datasetID2', ...] (e.g. ['SR1', 'SR2', ...])
                      or ['analsyisID:datasetID1',...] (e.g. ['CMS-PAS-EXO-16-036:c000'])
                      if the datasets have non-unique labels. A value of None means that all 
                      retrievable datasets should be retrieved.

        :return: a dictionary with the filtered dataset labels as keys and the corresponding datasets
                 as values.
        """

        if not isinstance(objList,(list,np.array,tuple)):
            logger.error("Input should be an iterable")
            raise SModelSError()

        filtered_datasets = []
        filtered_labels = []

        #Get all datasets appearing in the in the list of objects:
        datasets = []
        for obj in objList:
            if isinstance(obj,DataSet):
                datasets.append(obj)
            elif isinstance(obj,ExpResult):
                datasets += obj.datasets
            elif isinstance(obj,TheoryPrediction):
                datasets.append(obj.dataset)
            else:
                logger.error('Wrong input format. Must be a list of datasets, ExpResults or TheoryPrediction objects')
                raise SModelSError()

        #Loop over datasets and check which ones should be selected:
        for ds in datasets:
            #Get possible labels (datasetId or analysisID:datasetID)
            dsLabels = [ds.dataInfo.dataId]
            longlabel = ds.globalInfo.id+':'+ds.dataInfo.dataId
            dsLabels.append( longlabel )
            #If any of the labels match, add it to the filtered datasets
            for dsLabel in dsLabels:
                if filter_labels is None or dsLabel in filter_labels:
                    filtered_datasets.append(ds)
                    filtered_labels.append( longlabel )
                    continue

        #Make sure there are only unique labels:
        if len(np.unique(filtered_labels)) != len(filtered_labels):
            line = "Could not filter datasets with unique labels"
            logger.error( line )
            raise SModelSError( line )

        #Create a dictionary with {datasetLabel : dataset} entries:
        datasetDict = dict(list(zip(filtered_labels,filtered_datasets)))

        return datasetDict


    def srnameInList ( self, srname, srnameslist ):
        """ check if signal region name or its abbreviation is in list 'srname',
        :param srname: list of signal region names
        :returns: full signal region name, if in list, else False
        """
        # print ( "is", srname,"in", srnameslist )
        if srname in srnameslist:
            return srname
        isIn = False
        for x in srnameslist:
            if x.endswith ( ":"+srname ):
                if isIn is not False:
                    line = f"{srname} appears more than once in list"
                    raise SModelSError ( line )
                isIn = x
        p1 = srname.find(":")

        return isIn

    def getCorrelation ( self, longname, corrs ):
        """ get the correlation 'longname' from corrs """
        if longname in corrs:
            return corrs[longname]
        p1 = longname.find(":")
        shortname = longname[p1+1:]
        if shortname in corrs:
            return corrs[shortname]
        return {}

    def getLongNames ( self, corrs, longnames ):
        """ in the corrs dictionary, replace short names with long names """
        ret = {}
        for cname, line in corrs.items():
            if ":" in cname:
                if type(line)==dict:
                    line = self.getLongNames ( line, longnames )
                ret[cname]=line
                continue
            hasAppeared = False
            for c in longnames:
                if c.endswith ( ":"+cname ):
                    if hasAppeared:
                        logger.error ( f"{cname} is non-unique abbreviation for SR name in corrs. specify" )
                        raise SModelSError()
                    if type(line)==dict:
                        line = self.getLongNames ( line, longnames )
                    ret[c]=line
                    hasAppeared = True
        return ret

    def fromDatasets ( self, datasetDict, corrs: dict = {} ):
        """
        creates a fake experimental result with covariance matrix from a list of datasets

        :param datasetDict: a dictionary with the filtered dataset labels as keys and the corresponding datasets
                            as values.
        :param corrs: dictionary of dictionary of correlations,
                      using the dataset labels as keys,
                      e.g. { "SR1": { "SR2": 0.2, "SR4": 0.1 } }
                      omittance implies no correlation, correlation matrix
                      is by default assumed to be symmetric. The dataset labels must be the same used in the
                      datasetDict.
        """

        datasetOrder = sorted(list(datasetDict.keys())) #Sort for reproducibility
        corrs = self.getLongNames ( corrs, datasetOrder )

        datasets = [datasetDict[dsLabel] for dsLabel in datasetOrder]
        n_datasets = len(datasets)

        #Check if the dataset labels are unique:
        uniqueIds = np.unique(datasetOrder)
        if len(uniqueIds) < len(datasetOrder):
            logger.error ("Recurring dataset label when constructing a fake result: %s.\n Cannot yet handle." %(datasetOrder) )
            raise SModelSError()

        #Construct a covariance matrix:
        covariance_matrix = np.zeros((n_datasets,n_datasets))
        bgErrors = np.array([d_s.dataInfo.bgError for d_s in datasets])
        #Define the diagonal entries (in case they are not defined in the correlations dict)
        np.fill_diagonal(covariance_matrix,bgErrors**2)
        #Now include the correlations:
        for i,sr1 in enumerate(datasetOrder):
            srname = self.srnameInList ( sr1, datasetOrder )
            if not srname or not srname in corrs:
                continue #If datasetId does not appear in correlations, ignore
            sr1_corrs = corrs[srname]
            for j,sr2 in enumerate(datasetOrder):
                srname2 = self.srnameInList ( sr2, sr1_corrs )
                if not srname2:
                    continue #If datasetId does not appear in correlations for sr1, ignore
                corr = sr1_corrs [ srname2 ]
                covariance_matrix[i,j] = corr*bgErrors[i]*bgErrors[j]
                if covariance_matrix[j,i] == 0.:
                    covariance_matrix[j,i] = covariance_matrix[i,j]

        #Check consistency of the input correlations (should be symmetric and with diagonal entries = 1)
        if not np.allclose(covariance_matrix, covariance_matrix.T, rtol=1e-3, atol=1e-8):
            logger.warning("The provided correlations matrix is not symmetric. Symmetry will be enforced.")
        if not np.allclose(np.diag(covariance_matrix),bgErrors**2):
            logger.warning("The provided correlations matrix nas non-unity diagonal entries. Diagonal entries will be ignored.")

        #Make sure the covariance matrix is symmetric and the diagonals are bgErrors**2:
        covariance_matrix = (covariance_matrix + covariance_matrix.T)/2.
        np.fill_diagonal(covariance_matrix,bgErrors**2)

        #Get experimental result Ids:
        ana_ids = [d_s.globalInfo.id for d_s in datasets]

        #logger.debug ( "cov_matrx", covariance_matrix )
        #logger.debug ( "datasets", datasets )
        #logger.debug ( "ana_ids", ana_ids )

        #Construct a fake result with these N datasets and
        #an N x N covariance matrix
        self.fakeResult = ExpResult ( path = None, discard_zeroes = True,
                                      databaseParticles = None )
        self.fakeResult.datasets = datasets
        self.fakeResult.analysisIDs = ana_ids
        self.fakeResult.globalInfo = Info()
        self.fakeResult.globalInfo.id = f"FakeResult[{'+'.join(set(ana_ids))}]"
        self.fakeResult.globalInfo.datasetOrder = datasetOrder
        self.fakeResult.globalInfo.covariance = covariance_matrix

    @property
    def covariance(self):
        """
        Return covariance matrix
        """

        if self.fakeResult is not None:
            return self.fakeResult.globalInfo.covariance
        logger.error ('No fake result defined' )
        return np.array([])


if __name__ == "__main__":
    C = [ 18774.2, -2866.97, -5807.3, -4460.52, -2777.25, -1572.97, -846.653, -442.531,
       -2866.97, 496.273, 900.195, 667.591, 403.92, 222.614, 116.779, 59.5958,
       -5807.3, 900.195, 1799.56, 1376.77, 854.448, 482.435, 258.92, 134.975,
       -4460.52, 667.591, 1376.77, 1063.03, 664.527, 377.714, 203.967, 106.926,
       -2777.25, 403.92, 854.448, 664.527, 417.837, 238.76, 129.55, 68.2075,
       -1572.97, 222.614, 482.435, 377.714, 238.76, 137.151, 74.7665, 39.5247,
       -846.653, 116.779, 258.92, 203.967, 129.55, 74.7665, 40.9423, 21.7285,
       -442.531, 59.5958, 134.975, 106.926, 68.2075, 39.5247, 21.7285, 11.5732]
    nsignal = [ x/100. for x in [47,29.4,21.1,14.3,9.4,7.1,4.7,4.3] ]
    m=Data( observed=[1964,877,354,182,82,36,15,11],
              backgrounds=[2006.4,836.4,350.,147.1,62.0,26.2,11.1,4.7],
              covariance= C,
              #third_moment = [ 0.1, 0.02, 0.1, 0.1, 0.003, 0.0001, 0.0002, 0.0005 ],
              third_moment = [ 0. ] * 8,
              nsignal = nsignal,
              name="CMS-NOTE-2017-001 model" )
    ulComp = UpperLimitComputer(ntoys=500, cl=.95)
    #uls = ulComp.ulSigma ( Data ( 15,17.5,3.2,0.00454755 ) )
    #print ( "uls=", uls )
    ul_old = 131.828*sum(nsignal) #With respect to the older refernece value one must normalize the xsec
    print ( "old ul=", ul_old )
    ul = ulComp.ulSigma ( m )
    print ( "ul (marginalized)", ul )
    ul = ulComp.ulSigma ( m, marginalize=False )
    print ( "ul (profiled)", ul )
