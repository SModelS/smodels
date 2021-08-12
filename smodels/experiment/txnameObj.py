#!/usr/bin/env python3

"""
.. module:: txnameObj
   :synopsis: Holds the classes and methods used to read and store the
              information in the txname.txt files.
              Also contains the interpolation methods.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os,sys
from smodels.tools import physicsUnits
from smodels.tools.physicsUnits import GeV
from smodels.theory.auxiliaryFunctions import (elementsInStr, removeUnits, unscaleWidth,
                                               rescaleWidth, flattenArray, reshapeList,
                                               removeInclusives, addInclusives)
from smodels.tools.stringTools import concatenateLines
from smodels.theory.element import Element
from smodels.theory.topology import TopologyList
from smodels.tools.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.caching import _memoize
from smodels.tools.reweighting import defaultEffReweight,defaultULReweight
from scipy.linalg import svd, LinAlgError
import scipy.spatial.qhull as qhull
import numpy as np
import unum
import math,itertools
from math import floor, log10


#Build a dictionary with defined units. It can be used to evaluate
#expressions containing units.
unitsDict = dict([[varname,varobj] for varname,varobj \
                  in physicsUnits.__dict__.items()
                  if isinstance(varobj,unum.Unum)])


class TxName(object):
    """
    Holds the information related to one txname in the Txname.txt
    file (constraint, condition,...) as well as the data.
    """

    def __init__(self, path, globalObj, infoObj, databaseParticles):
        self.path = path
        self.globalInfo = globalObj
        self._infoObj = infoObj
        self.txnameData = None
        self.txnameDataExp = None ## expected Data
        self._topologyList = TopologyList()
        self.finalState = ['MET','MET'] #default final state
        self.intermediateState = None #default intermediate state

        logger.debug('Creating object based on txname file: %s' %self.path)
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Txname file %s not found" % path)
            raise SModelSError()
        txtFile = open(path,'r')
        txdata = txtFile.read()
        txtFile.close()
        if not "txName" in txdata: raise TypeError
        if not 'upperLimits' in txdata and not 'efficiencyMap' in txdata:
            raise TypeError
        content = concatenateLines(txdata.split("\n"))

        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        data = None
        expectedData = None
        dataType = None
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()
            if tags.count(tag) != 1:
                logger.info("Duplicated field %s found in file %s" \
                             % (tag, self.path))
            if ';' in value: value = value.split(';')
            if tag == 'upperLimits':
                data = value
                dataType = 'upperLimit'
            elif tag == 'expectedUpperLimits':
                expectedData = value
                dataType = 'upperLimit'
            elif tag == 'efficiencyMap':
                data = value
                dataType = 'efficiencyMap'
            else:
                self.addInfo(tag,value)

        ident = self.globalInfo.id+":"+dataType[0]+":"+ str(self._infoObj.dataId)
        ident += ":" + self.txName

        #Get detector size (if not found in self, look for it in datasetInfo or globalInfo).
        #If not defined anywhere, set it to None and default values will be used for reweighting.
        self.Leff_inner = self.fetchAttribute('Leff_inner',fillvalue=None)
        self.Leff_outer = self.fetchAttribute('Leff_outer',fillvalue=None)

        self.txnameData = TxNameData(data, dataType, ident,
                                        Leff_inner=self.Leff_inner,
                                        Leff_outer=self.Leff_outer)
        if expectedData:
            self.txnameDataExp = TxNameData( expectedData, dataType, ident,
                                            Leff_inner=self.Leff_inner,
                                            Leff_outer=self.Leff_outer)

        #Builds up a list of elements appearing in constraints:
        elements = []
        if not databaseParticles:
            raise SModelSError("Database particles is empty. Can not create TxName object.")
        if hasattr(self,'constraint'):
            elements += [Element(el,self.finalState,self.intermediateState,databaseParticles)
                    for el in elementsInStr(str(self.constraint))]

        if any((elA == elB and not elA is elB) for elA in elements for elB in elements):
            txt = "Duplicate elements in constraint: %s in %s" % \
                  ( self.constraint, self.globalInfo.id )
            logger.error( txt )
            raise SModelSError( txt )

        if hasattr(self,'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds,list): conds = [conds]
            for cond in conds:
                for el in elementsInStr(str(cond)):
                    newEl = Element(el,self.finalState,self.intermediateState,
                                    databaseParticles)
                    if not newEl in elements:
                        elements.append(newEl)

        # Builds up TopologyList with all the elements appearing in constraints
        # and conditions:
        for el in elements:
            self._topologyList.addElement(el)

    def hasOnlyZeroes(self):
        ozs = self.txnameData.onlyZeroValues()
        if self.txnameDataExp:
            e_ozs = self.txnameDataExp.onlyZeroValues()
            if ozs and e_ozs:
                return True
            if (ozs and not e_ozs) or (e_ozs and not ozs):
                logger.warning("%s is weird. One of the (expected, observed) results is zeroes-only, the other one isnt." )
                return False
        return ozs

    def fetchAttribute(self,attr,fillvalue=None):
        """
        Auxiliary method to get the attribute from self. If
        not found, look for it in datasetInfo and if still not found
        look for it in globalInfo.
        If not found in either of the above, return fillvalue.

        :param attr: Name of attribute (string)
        :param fillvalue: Value to be returned if attribute is not found.

        :return: Value of the attribute or fillvalue, if attribute was not found.
        """

        if hasattr(self,attr):
            return getattr(self,attr)
        elif  hasattr(self._infoObj,attr):
            return getattr(self._infoObj,attr)
        elif  hasattr(self.globalInfo,attr):
            return getattr(self.globalInfo,attr)
        else:
            return fillvalue

    def __str__(self):
        return self.txName

    def __repr__(self):
        return self.__str__()

    def __lt__ ( self, other ):
        """ sort by txName """
        return self.txName < other.txName

    def getULFor(self,element,expected=False):
        """
        Returns the upper limit (or expected) for element (only for upperLimit-type).
        Includes the lifetime reweighting (ul/reweight).
        If called for efficiencyMap results raises an error.
        If a mass array is given as input, no lifetime reweighting will be applied.

        :param element: Element object or mass array (with units)
        :param expected: look in self.txnameDataExp, not self.txnameData
        """
        if hasattr ( self, "dbClient" ):
            ## we have a databaseClient, so we send the request
            ## over the network
            # query = "obs:ATLAS-SUSY-2013-05:ul:T2bb:[[300,100],[300,100]]"
            query = "obs:"
            if expected:
                query = "exp:"
            query += self.globalInfo.id + ":ul:"
            query += self.txName + ":"
            query += self.getMassVectorFromElement ( element )
            logger.info ( "sending ul query %s to %s:%d" % \
                          ( query, self.dbClient.servername, self.dbClient.port ) )
            from smodels.tools.physicsUnits import fb
            return self.dbClient.query ( query )

        if not self.txnameData.dataType == 'upperLimit':
            logger.error("getULFor method can only be used in UL-type data.")
            raise SModelSError()

        if not expected:
            ul = self.txnameData.getValueFor(element)
        else:
            if not self.txnameDataExp:
                return None
            else:
                ul = self.txnameDataExp.getValueFor(element)

        return ul

    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.

        :param tag: information label (string)
        :param value: value for the field in string format
        """

        if tag == 'constraint' or tag == 'condition':
            if isinstance(value,list):
                value = [val.replace("'","") for val in value]
            else:
                value = value.replace("'","")
            if value == 'None':
                setattr(self,tag, eval(value))
            else:
                setattr(self,tag,value) #Make sure constraints/conditions are not evaluated
        else:
            try:
                setattr(self,tag,eval(value, unitsDict))
            except SyntaxError:
                setattr(self,tag,value)
            except NameError:
                setattr(self,tag,value)
            except TypeError:
                setattr(self,tag,value)

    def getInfo(self, infoLabel):
        """
        Returns the value of info field.

        :param infoLabel: label of the info field (string). It must be an attribute of
                          the TxNameInfo object
        """

        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False

    def hasElementAs(self,element):
        """
        Verify if the conditions or constraint in Txname contains the element.
        Check both branch orderings. If both orderings match, returns the one
        with the highest mass array.

        :param element: Element object
        :return: A copy of the element on the correct branch ordering appearing
                in the Txname constraint or condition.
        """

        #Stores all orderings of elements which matches the txname
        matches = []
        for el in self._topologyList.getElements():
            #Compare branches:
            for branchesA in itertools.permutations(element.branches):
                branchesA = list(branchesA)
                if branchesA == el.branches:
                    newEl = element.copy()
                    newEl.branches = [br.copy() for br in branchesA]
                    matches.append(newEl)

        #No elements matched:
        if not matches:
            return False
        elif len(matches) == 1:
            return matches[0]
        else: #If more than one element ordering matches, return the one with largest mass (relevant for clustering)
            matches = sorted(matches, key = lambda el: el.mass,reverse=True)
            return matches[0]

    def hasLikelihood(self):
        """ can I construct a likelihood for this map?
        True for all efficiency maps, and for upper limits maps
        with expected Values. """
        if self._infoObj.dataType == "efficiencyMap":
            return True
        if self.txnameDataExp != None:
            return True
        return False


    def getMassVectorFromElement(self, element ):
        """
        given element, extract the mass vector for the server query.
        element can be list of masses or "Element"

        :returns: eg [[300,100],[300,100]]
        """
        if type(element)==list:
            return str(element).replace(" [GeV]","").replace(" ","")
        from smodels.theory.clusterTools import AverageElement
        if type(element) == AverageElement:
            # ret += str(element.mass).replace(" [GeV]","").replace(" ","")
            return str(element.mass).replace(" [GeV]","").replace(" ","")
        ret = "["
        for i,br in enumerate(element.branches):
            ret += str(br.mass).replace(" [GeV]","").replace(" ","")
            if i+1 < len(element.branches):
                ret += ","
        ret += "]"
        # print ( "getMassVectorFromElement returning", ret )
        return ret

    def getQueryStringForElement ( self, element ):
        ## we have a databaseClient, so we send the request
        ## over the network
        # query = "obs:ATLAS-SUSY-2013-05:ul:T2bb:[[300,100],[300,100]]"
        query = "obs:"
        #if expected:
        #    query = "exp:"
        query += self.globalInfo.id + ":"
        dId = self._infoObj.dataId
        if dId == None:
            dId = "ul"
        query += dId + ":"
        query += self.txName + ":"
        query += self.getMassVectorFromElement ( element )
        return query

    def getEfficiencyFor(self,element):
        """
        For upper limit results, checks if the input element falls inside the
        upper limit grid and has a non-zero reweigthing factor.
        If it does, returns efficiency = 1, else returns
        efficiency = 0.  For efficiency map results, returns the
        signal efficiency including the lifetime reweighting.
        If a mass array is given as input, no lifetime reweighting will be applied.

        :param element: Element object or mass array with units.
        :return: efficiency (float)
        """

        if self.txnameData.dataType == 'efficiencyMap':
            if hasattr ( self, "dbClient" ):
                query = self.getQueryStringForElement ( element )
                logger.info ( "sending em query %s to %s:%d" % \
                              ( query, self.dbClient.servername, self.dbClient.port ) )
                #print ( "query will be", query )
                #return 0.001
                eff = self.dbClient.query ( query )
            else:
                eff = self.txnameData.getValueFor(element)

            if not eff or math.isnan(eff):
                eff = 0. #Element is outside the grid or has zero efficiency
        elif self.txnameData.dataType == 'upperLimit':
            if hasattr ( self, "dbClient" ):
                query = self.getQueryStringForElement ( element )
                logger.info ( "sending query %s to %s:%d" % \
                              ( query, self.dbClient.servername, self.dbClient.port ) )
                #print ( "query will be", query )
                #return 0.001
                ul = self.dbClient.query ( query )
            else:
                ul = self.txnameData.getValueFor(element)
            if isinstance(element,Element):
                element._upperLimit = ul #Store the upper limit for convenience
            if ul is None:
                eff = 0. #Element is outside the grid or the decays do not correspond to the txname
            else:
                eff = 1.
        else:
            logger.error("Unknown txnameData type: %s" %self.txnameData.dataType)
            raise SModelSError()

        return eff

class TxNameData(object):
    """
    Holds the data for the Txname object.  It holds Upper limit values or efficiencies.
    """
    _keep_values = False ## keep the original values, only for debugging

    def __init__(self,value,dataType,Id,
                    accept_errors_upto=.05,
                    Leff_inner=None,Leff_outer=None):
        """
        :param value: values in string format
        :param dataType: the dataType (upperLimit or efficiencyMap)
        :param Id: an identifier, must be unique for each TxNameData!
        :param _accept_errors_upto: If None, do not allow extrapolations outside of
                convex hull.  If float value given, allow that much relative
                uncertainty on the upper limit / efficiency
                when extrapolating outside convex hull.
                This method can be used to loosen the equal branches assumption.
        :param Leff_inner: is the effective inner radius of the detector, given in meters (used for reweighting prompt decays). If None, default values will be used.
        :param Leff_outer: is the effective outer radius of the detector, given in meters (used for reweighting decays outside the detector). If None, default values will be used.


        """
        self.dataType = dataType
        self._id = Id
        self._accept_errors_upto=accept_errors_upto
        self.Leff_inner = Leff_inner
        self.Leff_outer = Leff_outer
        self._V = None
        self.loadData(value)
        if self._keep_values:
            self.origdata = value

        if self.dataType == 'efficiencyMap':
            self.reweightF = defaultEffReweight
        elif self.dataType == 'upperLimit':
            self.reweightF = defaultULReweight
        else:
            raise SModelSError("Default reweighting function not defined for data type %s" %self.dataType)

    def __str__ ( self ):
        """ a simple unique string identifier, mostly for _memoize """
        return str ( self._id )

    def round_to_n(self, x, n):
        if x==0.0:
            return x
        return round(x, int(-np.sign(x)* int(floor(log10(abs(x)))) + (n - 1)))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__ ( self, other ):
        if type(self) != type(other):
            return False
        return self._id == other._id

    def evaluateString(self, value):
        """
        Evaluate string.

        :param value: String expression.
        """

        if not isinstance(value,str):
            raise SModelSError("Data should be in string format. Format %s found" %type(value))

        try:
            val = eval(value,unitsDict)
        except (NameError,ValueError,SyntaxError):
            raise SModelSError("data string malformed: %s" %value)

        return val

    def getUnits(self, value):
        """
        Get standard units for the input object.
        Uses the units defined in physicsUnits.standardUnits.
        (e.g. [[100*GeV,100.*GeV],3.*pb] -> returns [[GeV,GeV],fb]
        [[100*GeV,3.],[200.*GeV,2.*pb]] -> returns [[GeV,1.],[GeV,fb]] )

        :param value: Object containing units (e.g. [[100*GeV,100.*GeV],3.*pb])

        :return: Object with same structure containing the standard units used to
                 normalize the data.
        """

        stdUnits = physicsUnits.standardUnits
        if isinstance(value,list):
            return [self.getUnits(x) for x in value]
        if isinstance(value,tuple):
            return tuple([self.getUnits(x) for x in value])
        elif isinstance(value,dict):
            return dict([[self.getUnits(x),self.getUnits(y)]
                                  for x,y in value.items()])
        elif isinstance(value,unum.Unum):
            #Check if value has unit or not:
            if not value._unit:
                return 1.
            #Now try to find standard unit which matches:
            for unit in stdUnits:
                y = (value/unit).normalize()
                if not y._unit:
                    return unit
            raise SModelSError("Could not find standard unit which matches %s. Using the standard units: %s"
                               %(str(value),str(stdUnits)))
        else:
            return 1.

    def getDataShape(self,value):
        """
        Stores the data format (mass shape) and store it for future use.
        If there are inclusive objects (mass or branch = None), store their positions.

        :param value: list of data points
        """

        if isinstance(value,list):
            return [self.getDataShape(m) for m in value]
        elif isinstance(value,(float,int,unum.Unum,tuple)):
            return type(value)
        else:
            return value

    def getWidthPosition(self,value):
        """
        Gets the positions of the widths to be used for interpolation.

        :param value: data point

        :return: A list with the position of the widths. A position is a tuple
                 of the form (branch-index,vertex-index).
        """

        widthPositions = [(ibr,im)  for ibr,br in enumerate(value) for im,m in enumerate(br)
                          if isinstance(m,tuple)]

        return widthPositions

    def dataToCoordinates(self,dataPoint,rotMatrix=None,
                             transVector=None):
        """
        Format a dataPoint to the format used for interpolation.
        All the units are removed, the widths are rescaled and the masses
        and widths are combined in a flat array.
        The input can be an Element object or a massAndWidth nested arrays
        (with tuples to store the relevant widths).

        :param dataPoint: Element object from which the mass and width arrays will be extracted or
                          a nested mass array from the database, which contain tuples to include
                          the width values

        :param rotMatrix: Rotation matrix for PCA (e.g. self._V).
                          If None, no rotation is performed.
        :param transVector: Translation vector for PCA (e.g. self.delta_x).
                            If None no translation is performed

        :return: Point (list of floats)
        """

        #Collect the data
        if isinstance(dataPoint,Element):
            masses = dataPoint.mass
            widths = dataPoint.totalwidth
        elif isinstance(dataPoint,list):
            masses = [[mw[0] if isinstance(mw,tuple) else mw for mw in br] for br in dataPoint]
            widths = [[mw[1] if isinstance(mw,tuple) else None for mw in br] for br in dataPoint]
        else:
            logger.error("dataPoint must be an element or a nested array including masses and widths")
            raise SModelSError()

        #Select the required masses (remove entries corresponding to inclusive entries in data)
        masses = removeInclusives(masses, self.dataShape)
        #Select the required widths (remove widths not used in interpolation)
        widths = [[widths[ibr][im] for im,_ in enumerate(br)
                   if (ibr,im) in self.widthPosition]
                   for ibr,br in enumerate(widths)]
        if None in removeUnits ( flattenArray(widths), GeV ):
            logger.error("Error obtaining widths from %s" %str(dataPoint))
            raise SModelSError()

        #Remove units and flatten arrays:
        masses = flattenArray(masses)
        masses = removeUnits(masses,physicsUnits.standardUnits)
        widths = flattenArray(widths)
        widths = removeUnits(widths,physicsUnits.standardUnits)
        #Rescale widths:
        xwidths = [rescaleWidth(w) for w in widths]

        #Combine masses and rescaled widths in a single point
        point = masses + xwidths

        #Now transform to PCA coordinates (if rotMatrix and transVector are defined:
        if transVector is not None:
            point = np.array([point])
            point = ((point - transVector)).tolist()[0] #Translate
        if rotMatrix is not None:
            point = np.dot(point,rotMatrix)  # Rotate
            point = point.tolist()

        return point

    def coordinatesToData(self,point,rotMatrix=None,transVector=None):
        """
        A function that return the original mass and width array (including the widths
        as tuples) for a given point in PCA space (inverse of dataToCoordinates).

        :param point: Point in PCA space (1D list with size equal
                      to self.full_dimensionality or self.dimensionality)

        :param rotMatrix: Rotation matrix for PCA (e.g. self._V).
                          If None, no rotation is performed.
        :param transVector: Translation vector for PCA (e.g. self.delta_x).
                            If None no translation is performed

        :return: nested mass array including the widths as tuples (e.g. [[(200,1e-10),100],[(200,1e-10),100]])
        """

        if len(point) != self.full_dimensionality and len(point) != self.dimensionality:
            logger.error("Wrong point dimensions (%i), it should be %i (reduced dimensions) or %i (full dimensionts)"
                         %(len(point),self.dimensionality, self.full_dimensionality))
        elif len(point) != self.full_dimensionality:
            pointFull = np.array(point[:])
            pointFull = np.append(pointFull,[0.]*(self.full_dimensionality-len(point)))
        else:
            pointFull = np.array(point[:])

        massAndWidths = pointFull
        if rotMatrix is not None:
            massAndWidths = np.dot(rotMatrix,massAndWidths)
        if transVector is not None:
            massAndWidths = massAndWidths + transVector

        massAndWidths = massAndWidths.tolist()
        if type(massAndWidths[0])==list:
            massAndWidths = massAndWidths[0]
        #Extract masses and transformed widths
        masses = massAndWidths[:len(massAndWidths)-len(self.widthPosition)]
        xwidths = massAndWidths[len(massAndWidths)-len(self.widthPosition):]
        #Rescale widths and add unit:
        widths = [unscaleWidth(xw) for xw in xwidths]
        #Add units (make sure it is consistent with standardUnits)
        massUnit = [unit for unit in physicsUnits.standardUnits
                    if not (1*GeV/unit).normalize()._unit][0]
        masses = [m*massUnit for m in masses[:]]
        #Add inclusive entries to mass
        flatShape = flattenArray(self.dataShape)
        if len([x for x in flatShape if str(x) != '*']) != len(masses):
            logger.error("Error trying to add inclusive entries (%s) to flat mass array (%s)."
                         %(flatShape,masses))
            raise SModelSError()
        masses = addInclusives(masses, flatShape)
        #Reshape masses according to dataShape:
        if len(masses) != len(flatShape):
            logger.error("Number of elements in %s do not match the number of entries in %s"
                         %(masses,self.dataShape))
            raise SModelSError()

        massArray = reshapeList(masses, self.dataShape)
        #Add widths to the mass array
        if len(widths) != len(self.widthPosition):
            logger.error("The number of converted widths (%i) is not the expected (%i)"
                         %(len(widths),len(self.widthPosition)))
            raise SModelSError()

        #Combine masses and widths
        massAndWidthArray = []
        for ibr,br in enumerate(massArray):
            if str(br) != '*':
                newBr = [(m,widths.pop(0)) if (ibr,im) in self.widthPosition else m
                             for im,m in enumerate(br)]
            else:
                newBr = br
            massAndWidthArray.append(newBr)

        return massAndWidthArray

    def loadData(self,value):
        """
        Uses the information in value to generate the data grid used for
        interpolation.
        """

        if self._V:
            return

        if isinstance(value,str):
            val = self.evaluateString(value)
        else:
            val = value

        if len(val)==0:
            logger.error ( f"no values for {self._id} found" )
            sys.exit(-1)
        self.units = self.getUnits(val)[0] #Store standard units
        self.dataShape = self.getDataShape(val[0][0]) #Store the data (mass) format (useful if there are inclusives)
        self.widthPosition = self.getWidthPosition(val[0][0]) #Store the position of the required widths
        values = removeUnits(val,physicsUnits.standardUnits) #Remove units and store the normalization units

        if len(values) < 1 or len(values[0]) < 2:
                raise SModelSError("input value not in correct format. expecting sth " \
                               "like [ [ [[ 300.*GeV,100.*GeV], "\
                               "[ 300.*GeV,100.*GeV] ], 10.*fb ], ... ] "\
                               "for upper limits or [ [ [[ 300.*GeV,100.*GeV],"\
                               " [ 300.*GeV,100.*GeV] ], .1 ], ... ] for "\
                               "efficiency maps. Received %s" % values[:80])

        if not isinstance(self.units[-1],unum.Unum) and not isinstance(self.units[-1],float):
            raise SModelSError("Error obtaining units from value: %s " %values[:80])


        self.y_values = np.array(values,dtype=object)[:,1]
        self.computeV(values)

    def getValueFor(self,element):
        """
        Interpolates the value and returns the UL or efficiency for the
        respective element rescaled according to the reweighting function
        self.reweightF. For UL-type data the default rescaling is ul -> ul/(fraction of prompt decays)
        and for EM-type data it is eff -> eff*(fraction of prompt decays).
        If a mass array is given as input, no lifetime reweighting will be applied.

        :param element: Element object or mass array (with units)
        """

        #For backward compatibility:
        if not hasattr(self,'Leff_inner'):
            self.Leff_inner = None
        if not hasattr(self,'Leff_outer'):
            self.Leff_outer = None


        #Compute reweight factor according to lifetime/widths
        #For the widths not used in interpolation we assume that the
        #analysis require prompt decays
        #(width=inf for intermediate particles and width=0 for the last particle)
        if isinstance(element,Element):
            #Replaced the widths to be used for interpolation
            #with "prompt" widths (inf for intermediate particles and zero for final particles).
            #This way the reweight factor is only applied for the widths not used
            #for interpolation (since inf and zero result in no reweighting).
            widths = []
            for ibr,br in enumerate(element.totalwidth):
                widths.append([])
                for iw,w in enumerate(br):
                    if (ibr,iw) in self.widthPosition:
                        if iw != len(br)-1:
                            widths[ibr].append(float('inf')*GeV)
                        else:
                            widths[ibr].append(0.*GeV)
                    else:
                        widths[ibr].append(w)
            reweightFactor = self.reweightF(widths,
                                            Leff_inner=self.Leff_inner,
                                            Leff_outer=self.Leff_outer)
        elif isinstance(element,list):
            reweightFactor = 1.
        else:
            logger.error("Input of getValueFor must be an Element object or a mass array and not %s" %str(type(element)))
            raise SModelSError()

        #Returns None or zero, if reweightFactor is None or zero:
        if not reweightFactor:
            return reweightFactor

        #Extract the mass and width of the element
        #and convert it to the PCA coordinates (len(point) = self.full_dimensionality):
        point = self.dataToCoordinates(element,rotMatrix=self._V,
                                          transVector=self.delta_x)
        val = self.getValueForPoint(point)
        if not isinstance(val,(float,int,unum.Unum)):
            return val

        #Apply reweightFactor (if data has no width or partial width dependence)
        val *= reweightFactor

        return val

    @_memoize
    def getValueForPoint(self,point):
        """
        Returns the UL or efficiency for the point (in coordinates) using interpolation

        :param point: Point in coordinate space (length = self.full_dimensionality)

        :return: Value of UL or efficiency (float) without units
        """

        #Make sure the point is a numpy array
        point = np.array(point)
        self.projected_value = self.interpolate(point[:self.dimensionality])

        #Check if input point has larger dimensionality:
        dp = self.countNonZeros(point)
        if dp > self.dimensionality: ## we have data in different dimensions
            if self._accept_errors_upto == None:
                return None
            logger.debug( "attempting to interpolate outside of convex hull "\
                    "(d=%d,dp=%d,point=%s)" %
                     ( self.dimensionality, dp, str(point) ) )
            val =  self._interpolateOutsideConvexHull(point)
        else:
            val = self._returnProjectedValue()
        return val

    def interpolate(self, point, fill_value=np.nan):
        """
        Returns the interpolated value for the point (in coordinates)

        :param point: Point in coordinate space (length = self.dimensionality)

        :return: Value for point without units
        """


        tol = 1e-6
        # tol = sys.float_info.epsilon * 1e10
        simplex = self.tri.find_simplex(point, tol=tol)
        if simplex==-1: ## not inside any simplex?
            return fill_value

        #Transformation matrix for the simplex:
        simplexTrans = np.take(self.tri.transform, simplex, axis=0)
        #Space dimension:
        d = simplexTrans.shape[-1]
        #Rotation and translation to baryocentric coordinates:
        delta_x = simplexTrans[d,:]
        rot = simplexTrans[:d,:]
        bary = np.dot(rot,point-delta_x) #Point coordinates in the baryocentric system
        #Weights for the vertices:
        wts = np.append(bary, 1. - bary.sum())
        #Vertex indices:
        vertices = np.take(self.tri.simplices, simplex, axis=0)
        #Compute the value:
        values = np.array(self.y_values)
        ret = np.dot(np.take(values, vertices),wts)
        minXsec = min(np.take(values, vertices))
        if ret < minXsec:
            logger.debug('Interpolation below simplex values. Will take the smallest simplex value.')
            ret = minXsec
        return float(ret)

    def _estimateExtrapolationError(self, point):
        """
        When projecting a point from full_dimensionality to self.dimensionality, we
        estimate the expected extrapolation error with the following
        strategy: we compute the gradient at point P, and let alpha be the
        distance between p and P. We then walk one step of length alpha in
        the direction of the greatest ascent, and the opposite direction.
        Whichever relative change is greater is reported as the expected
        extrapolation error.

        :param point: Point in coordinate space (length = self.full_dimensionality)
        """


        #Make sure the point is a numpy array
        point = np.array(point)
        ## how far are we away from the "plane": distance alpha
        alpha = float(np.sqrt(np.dot(point[self.dimensionality:],
                        point[self.dimensionality:])))
        if alpha == 0.:
            ## no distance to the plane, so no extrapolation error
            return 0.
        ## the value of the grid at the point projected to the "plane"

        ## compute gradient
        gradient=[]
        for i in range(self.dimensionality):
            P2 = np.copy(point)
            P2[i] += alpha
            pv = self.interpolate(P2[:self.dimensionality])
            g = float((pv - self.projected_value)/alpha)
            if math.isnan(g):
                ## if we cannot compute a gradient, we return nan
                return float("nan")
            gradient.append(g)
        ## normalize gradient
        C = float(np.sqrt(np.dot(gradient, gradient)))
        if C == 0.:
            ## zero gradient? we return 0.
            return 0.
        for i,grad in enumerate(gradient):
            gradient[i]=grad/C*alpha
        ## walk one alpha along gradient
        P3 = np.copy(point)
        P4 = np.copy(point)
        for grad in gradient:
            P3[i]+= grad
            P4[i]-= grad
        agp=self.interpolate(P3[:self.dimensionality])
        agm=self.interpolate(P4[:self.dimensionality])
        dep,dem=0.,0.
        if self.projected_value == 0.:
            if agp!=0.:
                dep =1.0
            if agm!=0.:
                dem =1.0
        else:
            dep = abs(agp - self.projected_value)/self.projected_value
            dem = abs(agm - self.projected_value)/self.projected_value
        de=dep
        if dem > de:
            de=dem
        return de

    def _interpolateOutsideConvexHull(self, point):
        """
        Experimental routine, meant to check if we can interpolate outside
        convex hull

        :param point: Point in coordinate space (length = self.full_dimensionality)
        """

        #Make sure the point is a numpy array
        point = np.array(point)
        de = self._estimateExtrapolationError(point)

        if de < self._accept_errors_upto:
            return self._returnProjectedValue()

        if not math.isnan(de):
            logger.debug ( "Expected propagation error of %f too large to " \
                           "propagate." % de )
        return None

    def _returnProjectedValue(self):
        """
        Return interpolation result with the appropriate units.
        """

        ## None is returned without units'
        if self.projected_value is None or math.isnan(self.projected_value):
            logger.debug("Projected value is None. Projected point not in convex hull?")
            return None

        #Set value to zero if it is lower than machine precision (avoids fake negative values)
        if abs(self.projected_value) < 100.*sys.float_info.epsilon:
            self.projected_value = 0.

        return self.projected_value*self.units[-1]

    def countNonZeros(self, mp):
        """ count the nonzeros in a vector """
        nz=0
        lim = 10**-4
        for i in mp:
            if abs(i)>lim:
                nz+=1
        return nz

    def onlyZeroValues(self):
        """ check if the map is zeroes only """
        eps = sys.float_info.epsilon
        negative_values = bool ( sum ( [ x < -eps for x in self.y_values ] ) )
        if negative_values:
            for x in self.y_values:
                if x < -eps:
                    logger.error ( "negative error in result: %f, %s" % \
                                   ( x, self._id) )
                    sys.exit()
        if sum(self.y_values) > 0.:
            return False
        return True

    def computeV(self,values):
        """
        Compute rotation matrix _V, and triangulation self.tri

        :param values: Nested array with the data values without units

        """

        if not self._V is None:
            return

        #Convert nested mass arrays (with width tuples) to coordinates
        #(remove entries in mass corresponding to inclusive values,
        #select the required widths and combine masses and widths
        #in a flat array where the widths are the last entries)
        Morig= [self.dataToCoordinates(pt[0]) for pt in values]

        aM = np.array(Morig)
        MT = aM.T.tolist()
        self.delta_x = np.array([[sum(x)/len(Morig) for x in MT ]])
        M = []

        for Mx in Morig:
            m=(np.array([Mx]) - self.delta_x).tolist()[0]
            M.append(m)

        try:
            ## we dont need thousands of points for SVD
            n = int(math.ceil(len(M)/2000.))
            Vt=svd(M[::n])[2]
        except LinAlgError as e:
            raise SModelSError("exception caught when performing singular value decomposition: %s, %s" %(type(e), e))

        V=Vt.T
        self._V= V ## self.round ( V )
        Mp=[]

        ## the dimensionality of the whole mass space, disrespecting equal branches
        ## assumption
        self.full_dimensionality = len(Morig[0])
        self.dimensionality=0
        for m in M:
            mp=np.dot(m,V)
            Mp.append(mp)
            nz = self.countNonZeros(mp)
            if nz > self.dimensionality:
                self.dimensionality = nz
        MpCut=[]
        for i in Mp:
            MpCut.append(i[:self.dimensionality].tolist() )

        if self.dimensionality > 1:
            self.tri = qhull.Delaunay(MpCut)
        else:
            self.tri = Delaunay1D(MpCut)


class Delaunay1D:
    """
    Uses a 1D data array to interpolate the data.
    The attribute simplices is a list of N-1 pair of ints with the indices of the points
    forming the simplices (e.g. [[0,1],[1,2],[3,4],...]).
    """

    def __init__(self,data):

        self.points = None
        self.simplices = None
        self.transform = None
        if data and self.checkData(data):
            self.points = sorted(data)
            #Create simplices as the point intervals (using the sorted data)
            self.simplices = np.array([[data.index(self.points[i+1]),data.index(pt)]
                                       for i,pt in enumerate(self.points[:-1])])
            transform = []
            #Create trivial transformation to the baryocentric coordinates:
            for simplex in self.simplices:
                xmax,xmin = data[simplex[0]][0],data[simplex[1]][0]
                transform.append([[1./(xmax-xmin)],[xmin]])
            self.transform = np.array(transform)

            #Store convex hull (first and last point):
            self.convex_hull = np.array([data.index(self.points[0]),data.index(self.points[-1])])

        else:
            raise SModelSError()

    def find_simplex(self,x,tol=0.):
        """
        Find 1D data interval (simplex) to which x belongs

        :param x: Point (float) without units
        :param tol: Tolerance. If x is outside the data range with distance < tol, extrapolate.

        :return: simplex index (int)
        """

        xi = self.find_index(self.points,x)
        if xi == -1:
            if abs(x-self.points[0]) < tol:
                return 0
            else:
                return -1
        elif xi == len(self.simplices):
            if abs(x-self.points[-1]) < tol:
                return xi-1
            else:
                return -1
        else:
            return xi

    def checkData(self,data):
        """
        Define the simplices according to data. Compute and store
        the transformation matrix and simplices self.point.
        """
        if not isinstance(data,list):
            logger.error("Input data for 1D Delaunay should be a list.")
            return False
        for pt in data:
            if (not isinstance(pt,list)) or len(pt) != 1 or (not isinstance(pt[0],float)):
                logger.error("Input data for 1D Delaunay is in wrong format. It should be [[x1],[x2],..]")
                return False
        return True


    def find_index(self,xlist, x):
        """
        Efficient way to find x in a list.
        Returns the index (i) of xlist such that xlist[i] < x <= xlist[i+1].
        If x > max(xlist), returns the length of the list.
        If x < min(xlist), returns 0.        vertices = np.take(self.tri.simplices, simplex, axis=0)
        temp = np.take(self.tri.transform, simplex, axis=0)
        d=temp.shape[2]
        delta = uvw - temp[:, d]


        :param xlist: List of x-type objects
        :param x: object to be searched for.

        :return: Index of the list such that xlist[i] < x <= xlist[i+1].
        """

        lo = 0
        hi = len(xlist)
        while lo < hi:
            mid = (lo+hi)//2
            if xlist[mid] < x: lo = mid+1
            else: hi = mid
        return lo-1

if __name__ == "__main__":
    import time
    from smodels.tools.physicsUnits import GeV,fb
    data = [ [ [[ 150.*GeV, 50.*GeV], [ 150.*GeV, 50.*GeV] ],  3.*fb ],
         [ [[ 200.*GeV,100.*GeV], [ 200.*GeV,100.*GeV] ],  5.*fb ],
         [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ],
         [ [[ 300.*GeV,150.*GeV], [ 300.*GeV,150.*GeV] ], 13.*fb ],
         [ [[ 300.*GeV,200.*GeV], [ 300.*GeV,200.*GeV] ], 15.*fb ],
         [ [[ 300.*GeV,250.*GeV], [ 300.*GeV,250.*GeV] ], 20.*fb ],
         [ [[ 400.*GeV,100.*GeV], [ 400.*GeV,100.*GeV] ],  8.*fb ],
         [ [[ 400.*GeV,150.*GeV], [ 400.*GeV,150.*GeV] ], 10.*fb ],
         [ [[ 400.*GeV,200.*GeV], [ 400.*GeV,200.*GeV] ], 12.*fb ],
         [ [[ 400.*GeV,250.*GeV], [ 400.*GeV,250.*GeV] ], 15.*fb ],
         [ [[ 400.*GeV,300.*GeV], [ 400.*GeV,300.*GeV] ], 17.*fb ],
         [ [[ 400.*GeV,350.*GeV], [ 400.*GeV,350.*GeV] ], 19.*fb ], ]
    txnameData=TxNameData ( data, "upperLimit",  sys._getframe().f_code.co_name )
    t0=time.time()
    for masses in [ [[ 302.*GeV,123.*GeV], [ 302.*GeV,123.*GeV]],
                    [[ 254.*GeV,171.*GeV], [ 254.*GeV,170.*GeV]] ]:
        result=txnameData.getValueFor(masses)
        sm = "%.1f %.1f" % (masses[0][0].asNumber(GeV), masses[0][1].asNumber(GeV))
        print ( "%s %.3f fb" % (sm, result.asNumber(fb)))
    print ( "%.2f ms" % ((time.time()-t0)*1000.))
