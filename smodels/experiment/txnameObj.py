# !/usr/bin/env python3

"""
.. module:: txnameObj
   :synopsis: Holds the classes and methods used to read and store the
              information in the txname.txt files.
              Also contains the interpolation methods.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
import sys
from smodels.tools import physicsUnits
from smodels.tools.physicsUnits import GeV
from smodels.theory.auxiliaryFunctions import (elementsInStr, removeUnits,
                                               rescaleWidth, unscaleWidth)
from smodels.tools.stringTools import concatenateLines
from smodels.theory.element import Element
from smodels.theory.topology import TopologyDict
from smodels.tools.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.caching import _memoize
from smodels.tools.reweighting import defaultEffReweight, defaultULReweight
from scipy.linalg import svd, LinAlgError
import numpy as np
import unum
import math
from math import floor, log10


# Build a dictionary with defined units. It can be used to evaluate
# expressions containing units.
unitsDict = dict([[varname, varobj] for varname, varobj
                  in physicsUnits.__dict__.items()
                  if isinstance(varobj, unum.Unum)])


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
        self.txnameDataExp = None  # expected Data
        self.dataMap = None
        self.arrayMap = None
        self._topologyDict = TopologyDict()
        self.finalState = ['MET', 'MET']  # default final state
        self.intermediateState = None  # default intermediate state

        logger.debug('Creating object based on txname file: %s' % self.path)
        # Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Txname file %s not found" % path)
            raise SModelSError()
        txtFile = open(path, 'r')
        txdata = txtFile.read()
        txtFile.close()
        if "txName" not in txdata:
            raise TypeError
        if 'upperLimits' not in txdata and 'efficiencyMap' not in txdata:
            raise TypeError
        content = concatenateLines(txdata.split("\n"))

        # Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        data = None
        expectedData = None
        for i, tag in enumerate(tags):
            if not tag:
                continue
            line = content[i]
            value = line.split(':', 1)[1].strip()
            if tags.count(tag) != 1:
                logger.info("Duplicated field %s found in file %s"
                            % (tag, self.path))
            if ';' in value:
                value = value.split(';')
            if tag == 'upperLimits':
                data = value
                self.dataType = 'upperLimit'
            elif tag == 'expectedUpperLimits':
                expectedData = value
                self.dataType = 'upperLimit'
            elif tag == 'efficiencyMap':
                data = value
                self.dataType = 'efficiencyMap'
            else:
                self.addInfo(tag, value)

        if self.dataType == 'efficiencyMap':
            self.reweightF = defaultEffReweight
        elif self.dataType == 'upperLimit':
            self.reweightF = defaultULReweight

        ident = self.globalInfo.id+":"+self.dataType+":" + str(self._infoObj.dataId)
        ident += ":" + self.txName

        # Builds up a list of elements appearing in constraints:
        elements = []
        if not databaseParticles:
            raise SModelSError("Database particles is empty. Can not create TxName object.")
        # Create unsorted elements (in order to make sure its order matches the data grid)
        if hasattr(self, 'constraint'):
            elements += [Element(el, self.finalState,
                                 self.intermediateState,
                                 model=databaseParticles,
                                 sort=False)
                         for el in elementsInStr(str(self.constraint))]

        if any((elA == elB and elA is not elB) for elA in elements for elB in elements):
            txt = "Duplicate elements in constraint: %s in %s" % \
                  (self.constraint, self.globalInfo.id)
            logger.error(txt)
            raise SModelSError(txt)

        if hasattr(self, 'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds, list):
                conds = [conds]
            for cond in conds:
                for el in elementsInStr(str(cond)):
                    newEl = Element(el, self.finalState, self.intermediateState,
                                    databaseParticles)
                    if newEl not in elements:
                        elements.append(newEl)

        #  Builds up topologyDict with all the elements appearing in constraints
        #  and conditions:
        for el in elements:
            self._topologyDict.addElement(el)

        # Get detector size (if not found in self, look for it in datasetInfo or globalInfo).
        # If not defined anywhere, set it to None and default values will be used for reweighting.
        self.Leff_inner = self.fetchAttribute('Leff_inner', fillvalue=None)
        self.Leff_outer = self.fetchAttribute('Leff_outer', fillvalue=None)

        x_values, y_values = self.preProcessData(data)
        self.txnameData = TxNameData(x=x_values, y=y_values, txdataId=ident)
        if expectedData:
            x_values, y_values = self.preProcessData(expectedData)
            self.txnameDataExp = TxNameData(x=x_values, y=y_values, txdataId=ident)

    def __str__(self):
        return self.txName

    def __repr__(self):
        return self.__str__()

    def __lt__(self, other):
        """ sort by txName """
        return self.txName < other.txName

    def preProcessData(self, rawData):
        """
        Convert input data (from the upperLimits, expectedUpperLimits or efficiencyMap fields)
        to a flat array without units. The output is used to construct the TxNameData object,
        which will further process the data and interpolate it.
        It also builds the dictionary for translating Element properties to the flat data array.

        :parameter rawData: Raw data (either string or list)

        :return: Two flat lists of data, one for the model parameters and the other for the y values
                 (UL or efficiency values)
        """

        if isinstance(rawData, str):
            data = self.evaluateString(rawData)
        else:
            data = rawData

        if len(data) == 0:
            logger.error("no data values for %s found" % self)
            raise SModelSError("no data values for %s found" % self)
        elif len(data[0]) < 2:
            logger.error("No valid data found for %s" % self)
            raise SModelSError("No valid data found for %s" % self)

        xDataPoint = data[0][0]
        yDataPoint = data[0][1]
        # Store y-unit:
        self.y_unit = removeUnits(yDataPoint, returnUnit=True)[1]
        if not isinstance(self.y_unit, (unum.Unum, float)):
            raise SModelSError("Error obtaining units from value: %s " % data[0])

        # Define graph->data mapping
        self.getDataMap(xDataPoint)

        # Transform data:
        x_values, y_values = self.transformData(data)

        return x_values, y_values

    def transformData(self, data):
        """
        Uses the information in self.dataMap (or self.arrayMap) to convert data
        to a list of flat and unitless array. The data is split into two lists, one
        with the x-values (masses/widths) and
        another with the y-values (upper limits/efficiencies).

        :parameter data: 2-D array with the data grid ([[x-value,y-value], ...]).
                         The x-value can be a flat list (e.g. [mass1,mass2,mass3,width1])
                         or a nested list (e.g. [[(mass1,width1),mass2],[mass3]])
        """

        dataArray = np.array(data, dtype=object)
        xvalues = dataArray[:, 0]
        yvalues = dataArray[:, 1]

        # For x we must remove units, rescale widths and flatten array:
        xvalues = [self.transformPoint(x) for x in xvalues[:]]

        # For y we must just remove units:
        yvalues = [removeUnits(y) for y in yvalues[:]]

        return np.array(xvalues), np.array(yvalues)

    def transformPoint(self, x):
        """
        Transforms a x point (mass/width values) to a flat, unitless
        list. The widths are rescaled according to rescaleWidth.
        If x is already flat (e.g. [mass1,mass2,mass3,width3]),
        the transformation will use the mapping in self.dataMap.
        However, if x is a nested array (e.g. [[mass1,mass2],[(mass3,width3)]]),
        the transformation will be done according to the mapping defined in
        self.arrayMap.

        :parameter x: A list (or nested list) with mass/width values.

        :return: A flat and unitless list matching sel.dataMap.
        """

        x = np.array(x, dtype=object)
        # Get length of flat array:
        nDim = max([arrayIndex for arrayIndex in self.dataMap.keys()])+1
        xFlat = [None]*nDim
        for arrayIndex in self.dataMap:
            # If x is already flat, retrieve its value directly from x
            if x.ndim == 1:
                _, attr, unit = self.dataMap[arrayIndex]
                xval = x[arrayIndex]
            # If it is a nested bracket use the arrayMap:
            elif x.ndim == 2:
                multiIndex, attr, unit, _ = self.arrayMap[arrayIndex]
                i, j = multiIndex[:2]
                xval = x[i, j]
                if isinstance(xval, tuple):
                    xval = xval[multiIndex[2]]
            else:
                logger.error("x-value for data point has the wrong dimensions %s" % x)
                raise SModelSError()

            # Remove unit:
            if isinstance(xval, unum.Unum):
                xval = xval.asNumber(unit)
            # Rescale width:
            if attr == 'totalwidth':
                xval = rescaleWidth(xval)
            # Store the transformed value:
            xFlat[arrayIndex] = xval

        # Check if the list has been filled:
        if None in xFlat:
            logger.error("Error transforming point %s" % x)
            raise SModelSError()

        return xFlat

    def inverseTransformPoint(self, xFlat):
        """
        Transforms a 1D unitless array to a list of mass/width values.
        If self.arrayMap is defined, use it to convert to a nested
        bracket array foramt (e.g. [[mass1,(mass2,width2)],[mass3,mass4]]),
        otherwise convert it to a flat array (e.g. [mass1,mass2,mass3,mass4,width2])
        using self.dataMap.

        :parameter x: A 1D unitless array containing masses and rescaled widths

        :return: list (or nested list) with mass/width values (with units).
        """

        mLength = max(self.dataMap.keys())+1
        if self.arrayMap is not None:  # Convert to nested bracket
            xDim = max([v[0][0] for v in self.arrayMap.values()])
            yDim = max([v[0][1] for v in self.arrayMap.values()])
            massPoint = np.empty(shape=(xDim+1, yDim+1), dtype=object)
            for index in self.arrayMap:
                ijk, attr, unit, _ = self.arrayMap[index]
                value = xFlat[index]
                if attr == 'totalwidth':
                    value = unscaleWidth(value)
                value = value*unit
                i, j, k = ijk
                storedValue = massPoint[i, j]
                if isinstance(storedValue, tuple):
                    if k == 1:
                        massPoint[i, j] = (storedValue[0], value)
                    elif k == 0:
                        massPoint[i, j] = (value, storedValue[1])
                elif k == 1:
                    massPoint[i, j] = (storedValue, value)
                else:
                    massPoint[i, j] = value

        else:  # Convert to flat list
            massPoint = [None]*mLength
            for index in self.dataMap:
                node, attr, unit = self.dataMap[index]
                value = xFlat[index]
                if attr == 'totalwidth':
                    value = unscaleWidth(value)
                value = value*unit
                massPoint[index] = value

        return massPoint

    def getDataMap(self, massPoint):
        """
        Using the elements in the topology, construct a dictionary
        mapping the node.number, the node attributes and the corresponding
        index in flatten data array.
        If the dataMap has not been defined, construct from the element topology
        and data point format.

        :param dataPoint: A point with the x-values from the data grid
                          (e.g. [[100*GeV,(50*GeV,1e-3*GeV)],[100*GeV,(50*GeV,1e-3*GeV),10*GeV]])

        :return: Dictionary with the data mapping {dataArrayIndex : (nodeNumber,attr,unit)}
                (e.g. {0  : (1,'mass',GeV), 1 : (1, 'totalwidth',GeV),...})
        """

        massPoint = np.array(massPoint, dtype=object)
        # If dataMap has already been defined, check consistency:
        if self.dataMap is not None:

            if massPoint.ndim != 1:
                msgError = "Inconsistent data format."
                msgError += " The x-values (%s) should be a 1D array" % massPoint
                logger.error(msgError)
                raise SModelSError(msgError)
            dataLength = max([k for k in self.dataMap])+1
            if len(massPoint) != dataLength:
                msgError = "Inconsistent data format."
                msgError += " The number of values in x (%s) do not match dataMap" % massPoint
                logger.error(msgError)
                raise SModelSError(msgError)
            if sorted([k for k in self.dataMap]) != list(range(dataLength)):
                msgError = "Inconsistent data map."
                msgError += " The keys (%s) are missing array indices" % str(self.dataMap.keys())
                logger.error(msgError)
                raise SModelSError(msgError)

            return

        # If dataMap was not previously defined, the massPoint
        # should be a nested array (e.g. [[mass1,(mass2,width2)],[mass3,mass4]])
        if massPoint.ndim != 2:
            msgError = "Inconsistent data format."
            msgError += " The x-values (%s) should be a nested array" % massPoint
            logger.error(msgError)
            raise SModelSError(msgError)

        # Check if all elements in the txname share the same topology:
        if len(self._topologyDict) != 1:
            msgError = "Can not construct a data map for elements with distinct topologies"
            msgError += " (%s,%s)" % (self.globalInfo.id, self)
            raise SModelSError(msgError)

        # Since all elements are equivalent, use the first one
        # to define the map:
        el = self._topologyDict.getElements()[0]
        tree = el.tree

        # Get a nested array of nodes corresponding to the data point:
        nodeArray = []
        for mom, daughters in tree.dfs_successors().items():
            if mom == tree.root:  # Ignore PV
                continue
            nodeArray.append(mom)
            for d in daughters:
                if d.isSM:  # Ignore SM particles
                    continue
                if tree.out_degree(d) != 0:  # Ignore unstable daughters (will appear as mom)
                    continue
                nodeArray.append(d)

        try:
            nodeArray = np.reshape(nodeArray, massPoint.shape)
        except ValueError:
            msgError = "Txname element and data grid have inconsistent formats:"
            msgError += "\n %s\n and %s" % (nodeArray, massPoint)
            logger.error(msgError)
            raise SModelSError(msgError)

        # Iterate over the array and construct a map for the nodes
        # and the flat array:
        arrayMap = {}
        massIndex = 0  # Initial index for the masses
        widthIndex = len(nodeArray.flatten())  # Initial index for the widths
        for index in np.ndindex(massPoint.shape):
            node = nodeArray[index]
            if node.isInclusive:
                continue

            arrayValue = massPoint[index]
            mass, massUnit, width, widthUnit = self.getDataEntry(arrayValue)
            # Add entry for mass
            if mass is not None:
                arrayMap[massIndex] = ((*index, 0), 'mass', massUnit, node.node)
                massIndex += 1
            # Add entry for width
            if width is not None:
                arrayMap[widthIndex] = ((*index, 1), 'totalwidth', widthUnit, node.node)
                widthIndex += 1

        # Store the nested bracket <-> flat array map
        self.arrayMap = arrayMap

        # Also store graph <-> flat array map
        self.dataMap = {}
        for key, val in self.arrayMap.items():
            bracketIndex, attr, unit, nodeIndex = val
            self.dataMap[key] = (nodeIndex, attr, unit)

    def getDataEntry(self, arrayValue):
        """
        Given an array value, extract the masses, widths and their
        units from the array.

        :parameter arrayValue: List with masses and/or masses and widths
                               (e.g. [100*GeV, (50*GeV,1e-6*GeV)])
        """

        mass, massUnit = None, None
        width, widthUnit = None, None
        if isinstance(arrayValue, tuple):
            mass, width = arrayValue
            massUnit = unum.Unum(mass._unit)
            widthUnit = unum.Unum(width._unit)
        elif isinstance(arrayValue, unum.Unum):
            mass = arrayValue
            massUnit = unum.Unum(mass._unit)
        else:
            logger.error("Can not convert array value %s " % (arrayValue))
            raise SModelSError()

        return mass, massUnit, width, widthUnit

    def getDataFromElement(self, element):

        dataMap = self.dataMap
        elementData = [None]*(1+max(dataMap.keys()))
        for indexArray, nodeTuple in dataMap.items():
            nodeNumber, attr, unit = nodeTuple
            node = element.tree.getNode(nodeNumber)
            value = getattr(node, attr)
            if isinstance(unit, unum.Unum):
                value = value.asNumber(unit)
            if attr == 'totalwidth':
                value = rescaleWidth(value)
            elementData[indexArray] = value

        return elementData

    def getReweightingFor(self, element):
        """
        Compute the lifetime reweighting for the element (fraction of prompt decays).
        If element is a list, return 1.0.

        :param element: Element object

        :return: Reweighting factor (float)
        """

        if not isinstance(element, Element):
            msgError = "Input of getReweightingFor must be an Element object"
            msgError += " and not %s" % str(type(element))
            logger.error(msgError)
            raise SModelSError()

        # For backward compatibility:
        if not hasattr(self, 'Leff_inner'):
            self.Leff_inner = None
        if not hasattr(self, 'Leff_outer'):
            self.Leff_outer = None

        # Collect the widths which are taken into accound by data:
        widthsInData = []
        for arrayIndex in self.dataMap:
            node, attr, _ = self.dataMap[arrayIndex]
            if 'width' in attr:
                widthsInData.append(node)

        # Get the widths for all unstable particles  and final state
        # particles not appearing in data:
        unstableWidths = []
        stableWidths = []
        tree = element.tree
        for mom, daughters in tree.dfs_successors().items():
            if mom == tree.root:
                continue  # Ignore primary vertex
            if mom.isInclusive:
                continue  # Ignore inclusive nodes
            if mom.node not in widthsInData:
                unstableWidths.append(mom.totalwidth)
            for d in daughters:
                if tree.out_degree(d) != 0:
                    continue   # Skip intermediate states
                if d.isInclusive:
                    continue
                if d.isSM:
                    continue
                if d.node not in widthsInData:
                    stableWidths.append(d.totalwidth)

        # Compute reweight factor according to lifetime/widths
        # For the widths not used in interpolation we assume that the
        # analysis require prompt decays
        # (width=inf for intermediate particles and width=0 for the last particle)
        reweightFactor = self.reweightF(unstableWidths=unstableWidths,
                                        stableWidths=stableWidths,
                                        Leff_inner=self.Leff_inner,
                                        Leff_outer=self.Leff_outer)
        return reweightFactor

    def evaluateString(self, value):
        """
        Evaluate string.

        :param value: String expression.
        """

        if not isinstance(value, str):
            raise SModelSError("Data should be in string format. Format %s found" % type(value))

        try:
            val = eval(value, unitsDict)
        except (NameError, ValueError, SyntaxError):
            raise SModelSError("data string malformed: %s" % value)

        return val

    def hasOnlyZeroes(self):
        ozs = self.txnameData.onlyZeroValues()
        if self.txnameDataExp:
            e_ozs = self.txnameDataExp.onlyZeroValues()
            if ozs and e_ozs:
                return True
            if (ozs and not e_ozs) or (e_ozs and not ozs):
                logger.warning("%s is weird. One of the (expected, observed) results is zeroes-only, the other one isnt.")
                return False
        return ozs

    def fetchAttribute(self, attr, fillvalue=None):
        """
        Auxiliary method to get the attribute from self. If
        not found, look for it in datasetInfo and if still not found
        look for it in globalInfo.
        If not found in either of the above, return fillvalue.

        :param attr: Name of attribute (string)
        :param fillvalue: Value to be returned if attribute is not found.

        :return: Value of the attribute or fillvalue, if attribute was not found.
        """

        if hasattr(self, attr):
            return getattr(self, attr)
        elif hasattr(self._infoObj, attr):
            return getattr(self._infoObj, attr)
        elif hasattr(self.globalInfo, attr):
            return getattr(self.globalInfo, attr)
        else:
            return fillvalue

    def getULFor(self, element, expected=False):
        """
        Returns the upper limit (or expected) for element (only for upperLimit-type).
        Includes the lifetime reweighting (ul/reweight).
        If called for efficiencyMap results raises an error.
        If a mass array is given as input, no lifetime reweighting will be applied.

        :param element: Element object or mass array (with units)
        :param expected: look in self.txnameDataExp, not self.txnameData
        """

        if not self.dataType == 'upperLimit':
            logger.error("getULFor method can only be used in UL-type data.")
            raise SModelSError()

        point = self.getDataFromElement(element)
        if not expected:
            ul = self.txnameData.getValueFor(point)
        else:
            if not self.txnameDataExp:
                return None
            else:
                ul = self.txnameDataExp.getValueFor(point)

        if ul is None:
            return None

        # Compute reweighting factor:
        reweightF = self.getReweightingFor(element)
        if reweightF is None:
            return None

        ul = ul*reweightF*self.y_unit  # Add unit

        return ul

    def getEfficiencyFor(self, element):
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

        # Get flat data from element:
        point = self.getDataFromElement(element)
        if self.dataType == 'efficiencyMap':
            eff = self.txnameData.getValueFor(point)
            if not eff or math.isnan(eff):
                eff = 0.  # Element is outside the grid or has zero efficiency
            # Compute reweighting factor:
            reweightF = self.getReweightingFor(element)
            eff = eff*reweightF*self.y_unit  # (unit should be 1)

        elif self.txnameData.dataType == 'upperLimit':
            ul = self.txnameData.getValueFor(point)
            element._upperLimit = ul  # Store the upper limit for convenience
            if ul is None:
                eff = 0.  # Element is outside the grid or the decays do not correspond to the txname
            else:
                eff = 1.
        else:
            logger.error("Unknown txnameData type: %s" % self.txnameData.dataType)
            raise SModelSError()

        return eff

    def addInfo(self, tag, value):
        """
        Adds the info field labeled by tag with value value to the object.

        :param tag: information label (string)
        :param value: value for the field in string format
        """

        if tag == 'constraint' or tag == 'condition':
            if isinstance(value, list):
                value = [val for val in value]
            else:
                value = value
            if value == 'None':
                setattr(self, tag, eval(value))
            else:
                setattr(self, tag, value)  # Make sure constraints/conditions are not evaluated
        else:
            try:
                setattr(self, tag, eval(value, unitsDict))
            except SyntaxError:
                setattr(self, tag, value)
            except NameError:
                setattr(self, tag, value)
            except TypeError:
                setattr(self, tag, value)

    def hasElementAs(self, element):
        """
        Verify if the conditions or constraint in Txname contains the element.
        Check both branch orderings. If both orderings match, returns the one
        with the highest mass array.

        :param element: Element object
        :return: A copy of the element on the correct branch ordering appearing
                in the Txname constraint or condition.
        """

        cName = element.canonName
        # Check if the canonical name matches any of the
        # elements in self:
        if cName not in self._topologyDict:
            return False

        # Get list of elements with the same canonical name
        elList = self._topologyDict[cName]
        for el in elList:
            # Compare elements:
            cmp, sortedEl = el.compareTo(element)
            if cmp == 0:
                return sortedEl

        # If this point was reached, there were no macthes
        return False

    def hasLikelihood(self):
        """ can I construct a likelihood for this map?
        True for all efficiency maps, and for upper limits maps
        with expected Values. """
        if self._infoObj.dataType == "efficiencyMap":
            return True
        if self.txnameDataExp is not None:
            return True
        return False


class TxNameData(object):
    """
    Holds the data for the Txname object.  It holds Upper limit values or efficiencies.
    """
    _keep_values = False  # keep the original values, only for debugging

    def __init__(self, x, y, txdataId,
                 accept_errors_upto=.05):
        """
        :param x: 2-D list of flat and unitless x-points (e.g. [ [mass1,mass2,mass3,mass4], ...])
        :param y: 1-D list with y-values (upper limits or efficiencies)
        :param _accept_errors_upto: If None, do not allow extrapolations outside of
                convex hull.  If float value given, allow that much relative
                uncertainty on the upper limit / efficiency
                when extrapolating outside convex hull.
                This method can be used to loosen the equal branches assumption.


        """
        self._id = txdataId
        self._accept_errors_upto = accept_errors_upto
        self._V = None
        self.y_values = y[:]
        # Compute PCA transformation:
        self.computeV(x)
        if self._keep_values:
            self.origdata = x

    def __str__(self):
        """ a simple unique string identifier, mostly for _memoize """
        return str(self._id)

    def round_to_n(self, x, n):
        if x == 0.0:
            return x
        return round(x, int(-np.sign(x) * int(floor(log10(abs(x)))) + (n - 1)))

    def __ne__(self, other):
        return not self.__eq__(other)

    def __eq__(self, other):
        if type(self) != type(other):
            return False
        return self._id == other._id

    def PCAtransf(self, point):
        """
        Transform a flat/unitless point with masses/widths to the PCA
        coordinate space.

        :param point: Flat and unitless mass/rescaled width point (e.g. [mass1,mass2,width1]).
                      Its length should be equal to self.full_dimensionality.

        :return: 1D array in coordinate space

        """

        # Transform to PCA coordinates (if rotMatrix and transVector are defined:
        transVector = self.delta_x  # Translation vector
        rotMatrix = self._V  # Rotation matrix

        point = np.array([point])
        point = ((point - transVector)).tolist()[0]  # Translate
        point = np.dot(point, rotMatrix)  # Rotate

        return point

    def inversePCAtransf(self, point):
        """
        Transform a a flat 1D point from coordinate space to flat/unitless point
        with masses/rescaled widths.

        :param point: 1D array in coordinate space

        :return: Flat and unitless mass/rescaled width point (e.g. [mass1,mass2,width1]).

        """

        if len(point) != self.full_dimensionality and len(point) != self.dimensionality:
            msgError = "Wrong point dimensions (%i)," % (len(point))
            msgError += " it should be % i(reduced dimensions)" % self.dimensionality
            msgError += " or %i(full dimensionts)" % self.full_dimensionality
            logger.error(msgError)
            raise SModelSError(msgError)

        elif len(point) != self.full_dimensionality:
            pointFull = np.array(point[:])
            pointFull = np.append(pointFull, [0.]*(self.full_dimensionality-len(point)))
        else:
            pointFull = np.array(point[:])

        # Transform to PCA coordinates (if rotMatrix and transVector are defined:
        transVector = self.delta_x  # Translation vector
        rotMatrix = self._V  # Rotation matrix

        point = np.array(pointFull)
        point = np.dot(rotMatrix, point)  # Rotate
        point = ((point + transVector)).tolist()[0]  # Translate

        return point

    @_memoize
    def getValueFor(self, point):
        """
        Returns the UL or efficiency for the point.

        :param point: Flat and unitless mass/width point (e.g. [mass1,mass2,width1]).
                      Its length should be equal to self.full_dimensionality.

        :return: Interpolated value for the grid (without units)
        """

        # Transform point accordint to PCA transformation:
        point = self.PCAtransf(point)

        self.projected_value = self.interpolate(point[:self.dimensionality])

        # Check if input point has larger dimensionality:
        dp = self.countNonZeros(point)
        if dp > self.dimensionality:  # we have data in different dimensions
            if self._accept_errors_upto is None:
                return None
            logger.debug("attempting to interpolate outside of convex hull "
                         "(d=%d,dp=%d,point=%s)" %
                         (self.dimensionality, dp, str(point)))
            val = self._interpolateOutsideConvexHull(point)
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
        #  tol = sys.float_info.epsilon * 1e10
        simplex = self.tri.find_simplex(point, tol=tol)
        if simplex == -1:  # not inside any simplex?
            return fill_value

        # Transformation matrix for the simplex:
        simplexTrans = np.take(self.tri.transform, simplex, axis=0)
        # Space dimension:
        d = simplexTrans.shape[-1]
        # Rotation and translation to baryocentric coordinates:
        delta_x = simplexTrans[d, :]
        rot = simplexTrans[:d, :]
        bary = np.dot(rot, point-delta_x)  # Point coordinates in the baryocentric system
        # Weights for the vertices:
        wts = np.append(bary, 1. - bary.sum())
        # Vertex indices:
        vertices = np.take(self.tri.simplices, simplex, axis=0)
        # Compute the value:
        values = np.array(self.y_values)
        ret = np.dot(np.take(values, vertices), wts)
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

        # Make sure the point is a numpy array
        point = np.array(point)
        # #  how far are we away from the "plane": distance alpha
        alpha = float(np.sqrt(np.dot(point[self.dimensionality:],
                                     point[self.dimensionality:])))
        if alpha == 0.:
            # #  no distance to the plane, so no extrapolation error
            return 0.
        # #  the value of the grid at the point projected to the "plane"

        # #  compute gradient
        gradient = []
        for i in range(self.dimensionality):
            P2 = np.copy(point)
            P2[i] += alpha
            pv = self.interpolate(P2[:self.dimensionality])
            g = float((pv - self.projected_value)/alpha)
            if math.isnan(g):
                # #  if we cannot compute a gradient, we return nan
                return float("nan")
            gradient.append(g)
        # #  normalize gradient
        C = float(np.sqrt(np.dot(gradient, gradient)))
        if C == 0.:
            # #  zero gradient? we return 0.
            return 0.
        for i, grad in enumerate(gradient):
            gradient[i] = grad/C*alpha
        # #  walk one alpha along gradient
        P3 = np.copy(point)
        P4 = np.copy(point)
        for grad in gradient:
            P3[i] += grad
            P4[i] -= grad
        agp = self.interpolate(P3[:self.dimensionality])
        agm = self.interpolate(P4[:self.dimensionality])
        dep, dem = 0., 0.
        if self.projected_value == 0.:
            if agp != 0.:
                dep = 1.0
            if agm != 0.:
                dem = 1.0
        else:
            dep = abs(agp - self.projected_value)/self.projected_value
            dem = abs(agm - self.projected_value)/self.projected_value
        de = dep
        if dem > de:
            de = dem
        return de

    def _interpolateOutsideConvexHull(self, point):
        """
        Experimental routine, meant to check if we can interpolate outside
        convex hull

        :param point: Point in coordinate space (length = self.full_dimensionality)
        """

        # Make sure the point is a numpy array
        point = np.array(point)
        de = self._estimateExtrapolationError(point)

        if de < self._accept_errors_upto:
            return self._returnProjectedValue()

        if not math.isnan(de):
            logger.debug("Expected propagation error of %f too large to "
                         "propagate." % de)
        return None

    def _returnProjectedValue(self):
        """
        Return interpolation result with the appropriate units.
        """

        # #  None is returned without units'
        if self.projected_value is None or math.isnan(self.projected_value):
            logger.debug("Projected value is None. Projected point not in convex hull?")
            return None

        # Set value to zero if it is lower than machine precision (avoids fake negative values)
        if abs(self.projected_value) < 100.*sys.float_info.epsilon:
            self.projected_value = 0.

        return self.projected_value

    def countNonZeros(self, mp):
        """ count the nonzeros in a vector """
        nz = 0
        lim = 10**-4
        for i in mp:
            if abs(i) > lim:
                nz += 1
        return nz

    def onlyZeroValues(self):
        """ check if the map is zeroes only """
        eps = sys.float_info.epsilon
        negative_values = bool(sum([x < -eps for x in self.y_values]))
        if negative_values:
            for x in self.y_values:
                if x < -eps:
                    logger.error("negative error in result: %f, %s" %
                                 (x, self._id))
                    sys.exit()
        if sum(self.y_values) > 0.:
            return False
        return True

    def computeV(self, x):
        """
        Compute rotation matrix _V, and triangulation self.tri

        :parameter x: 2-D array with the flatten x-points without units
                      (e.g. [ [mass1,mass2,mass3,mass4], [mass1',mass2',mass3',mass4'], ...])

        """

        if self._V is not None:
            return

        # Convert nested mass arrays (with width tuples) to coordinates
        # (remove entries in mass corresponding to inclusive values,
        # select the required widths and combine masses and widths
        # in a flat array where the widths are the last entries)
        Morig = x[:]

        aM = np.array(Morig)
        MT = aM.T.tolist()
        self.delta_x = np.array([[sum(x)/len(Morig) for x in MT]])
        M = []

        for Mx in Morig:
            m = (np.array([Mx]) - self.delta_x).tolist()[0]
            M.append(m)

        try:
            # #  we dont need thousands of points for SVD
            n = int(math.ceil(len(M)/2000.))
            Vt = svd(M[::n])[2]
        except LinAlgError as e:
            raise SModelSError("exception caught when performing singular value decomposition: %s, %s" % (type(e), e))

        V = Vt.T
        self._V = V  # self.round ( V )
        Mp = []

        # #  the dimensionality of the whole mass space, disrespecting equal branches
        # #  assumption
        self.full_dimensionality = len(Morig[0])
        self.dimensionality = 0
        for m in M:
            mp = np.dot(m, V)
            Mp.append(mp)
            nz = self.countNonZeros(mp)
            if nz > self.dimensionality:
                self.dimensionality = nz
        MpCut = []
        for i in Mp:
            MpCut.append(i[:self.dimensionality].tolist())

        if self.dimensionality > 1:
            try:
                from scipy.spatial import Delaunay
            except (ImportError, ModuleNotFoundError):
                from scipy.spatial.qhull import Delaunay
            self.tri = Delaunay(MpCut)
        else:
            self.tri = Delaunay1D(MpCut)


class Delaunay1D:
    """
    Uses a 1D data array to interpolate the data.
    The attribute simplices is a list of N-1 pair of ints with the indices of the points
    forming the simplices (e.g. [[0,1],[1,2],[3,4],...]).
    """

    def __init__(self, data):

        self.points = None
        self.simplices = None
        self.transform = None
        if data and self.checkData(data):
            self.points = sorted(data)
            # Create simplices as the point intervals (using the sorted data)
            self.simplices = np.array([[data.index(self.points[i+1]), data.index(pt)]
                                       for i, pt in enumerate(self.points[:-1])])
            transform = []
            # Create trivial transformation to the baryocentric coordinates:
            for simplex in self.simplices:
                xmax, xmin = data[simplex[0]][0], data[simplex[1]][0]
                transform.append([[1./(xmax-xmin)], [xmin]])
            self.transform = np.array(transform)

            # Store convex hull (first and last point):
            self.convex_hull = np.array([data.index(self.points[0]), data.index(self.points[-1])])

        else:
            raise SModelSError()

    def find_simplex(self, x, tol=0.):
        """
        Find 1D data interval (simplex) to which x belongs

        :param x: Point (float) without units
        :param tol: Tolerance. If x is outside the data range with distance < tol, extrapolate.

        :return: simplex index (int)
        """

        xi = self.find_index(self.points, x)
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

    def checkData(self, data):
        """
        Define the simplices according to data. Compute and store
        the transformation matrix and simplices self.point.
        """
        if not isinstance(data, list):
            logger.error("Input data for 1D Delaunay should be a list.")
            return False
        for pt in data:
            if (not isinstance(pt, list)) or len(pt) != 1 or (not isinstance(pt[0], float)):
                logger.error("Input data for 1D Delaunay is in wrong format. It should be [[x1],[x2],..]")
                return False
        return True

    def find_index(self, xlist, x):
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
            if xlist[mid] < x:
                lo = mid+1
            else:
                hi = mid
        return lo-1


if __name__ == "__main__":
    import time
    data = [[[[150.*GeV, 50.*GeV], [150.*GeV, 50.*GeV]],  3.*fb],
         [[[200.*GeV, 100.*GeV], [200.*GeV, 100.*GeV]],  5.*fb],
         [[[300.*GeV, 100.*GeV], [300.*GeV, 100.*GeV]], 10.*fb],
         [[[300.*GeV, 150.*GeV], [300.*GeV, 150.*GeV]], 13.*fb],
         [[[300.*GeV, 200.*GeV], [300.*GeV, 200.*GeV]], 15.*fb],
         [[[300.*GeV, 250.*GeV], [300.*GeV, 250.*GeV]], 20.*fb],
         [[[400.*GeV, 100.*GeV], [400.*GeV, 100.*GeV]],  8.*fb],
         [[[400.*GeV, 150.*GeV], [400.*GeV, 150.*GeV]], 10.*fb],
         [[[400.*GeV, 200.*GeV], [400.*GeV, 200.*GeV]], 12.*fb],
         [[[400.*GeV, 250.*GeV], [400.*GeV, 250.*GeV]], 15.*fb],
         [[[400.*GeV, 300.*GeV], [400.*GeV, 300.*GeV]], 17.*fb],
         [[[400.*GeV, 350.*GeV], [400.*GeV, 350.*GeV]], 19.*fb], ]
    txnameData = TxNameData(data, "upperLimit",  sys._getframe().f_code.co_name)
    t0 = time.time()
    for masses in [[[302.*GeV, 123.*GeV], [302.*GeV, 123.*GeV]],
                   [[254.*GeV, 171.*GeV], [254.*GeV, 170.*GeV]]]:
        result = txnameData.getValueFor(masses)
        sm = "%.1f %.1f" % (masses[0][0].asNumber(GeV), masses[0][1].asNumber(GeV))
        print("%s %.3f fb" % (sm, result.asNumber(fb)))
    print("%.2f ms" % ((time.time()-t0)*1000.))
    print("%.2f ms" % ((time.time()-t0)*1000.))
    print("%.2f ms" % ((time.time()-t0)*1000.))
