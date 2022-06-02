# !/usr/bin/env python3

"""
.. module:: txnameObj
   :synopsis: Holds the class for storing the simplified model information
              and data  for experimental results.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
from smodels.tools import physicsUnits
from smodels.theory.auxiliaryFunctions import (elementsInStr, removeUnits,
                                               rescaleWidth, unscaleWidth, cSim, cGtr)
from smodels.tools.stringTools import concatenateLines
from smodels.theory.element import Element
from smodels.tools.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.experiment.txnameDataObj import TxNameData
from smodels.tools.reweighting import defaultEffReweight, defaultULReweight
import numpy as np
import unum
import math


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
        self.elementMap = {}  # Stores the elements and their label representaion
        self._constraintFunc = None
        self._conditionsFunc = None
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

        if not databaseParticles:
            raise SModelSError("Database particles is empty. Can not create TxName object.")

        # Process constraint, simplify it so it can be easily evaluated and
        # stores the elements in self.elementMap
        if hasattr(self, 'constraint'):
            self._constraintFunc = self.processExpr(self.constraint,
                                                    databaseParticles,
                                                    checkUnique=True)

        # Process conditions, simplify it so it can be easily evaluated and
        # stores the elements in self.elementMap
        if hasattr(self, 'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds, list):
                conds = [conds]
            conditionFuncs = [self.processExpr(cond, databaseParticles, checkUnique=False)
                              for cond in conds]
            self._conditionsFunc = conditionFuncs[:]

        # Do consistency checks:
        self.checkConsistency()

        # Define canonical name (should be the same for all elements)
        self.canonName = list(self.elementMap.keys())[0].canonName

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
        """
        Sort by txName
        """

        return self.txName < other.txName

    def checkConsistency(self):
        """
        Checks if all the elements in txname have the same structure (topology/canonical name)
        and verify if its constraints and conditions are valid expressions.

        :return: True if txname is consitency. Raises an error otherwise.
        """

        # Check if all elements have the same canonical name:
        elements = list(self.elementMap.keys())
        cName = elements[0].canonName
        for el in elements[1:]:
            if el.canonName != cName:
                msgError = "Txname %s referes to elements with distinct topologies." % self
                logger.error(msgError)
                raise SModelSError(msgError)

        # Check if the constraint can be evaluated:
        elements = []
        for iel, elLabel in enumerate(self.elementMap.values()):
            elTest = Element()
            elTest.weight = 3.0*(iel+1)*physicsUnits.fb  # dumyy xsec
            elTest.txlabel = elLabel
            elements.append(elTest)

        # Dummy constraint eval:
        try:
            res = self.evalConstraintFor(elements)
        except (TypeError, NameError) as e:
            msgError = "Can not evaluate constraint %s for %s." % (self._constraintFunc, self)
            msgError += "\n error: %s" % str(e)
            logger.error(msgError)
            raise SModelSError(msgError)

        if res is not None and not isinstance(res, (float, int, unum.Unum,)):
            msgError = "Constraint for %s returned an invalid value: %s (%s)" % (self, res, type(res))
            logger.error(msgError)
            raise SModelSError(msgError)

        # Dummy condistions eval:
        try:
            res = self.evalConditionsFor(elements)
        except (TypeError, NameError) as e:
            msgError = "Can not evaluate conditions %s for %s." % (self._conditionsFunc, self)
            msgError += "\n error: %s" % str(e)
            logger.error(msgError)
            raise SModelSError(msgError)

        if res is not None and not isinstance(res, list):
            msgError = "Conditions for %s returned an invalid value: %s (%s)" % (self, res, type(res))
            logger.error(msgError)
            raise SModelSError(msgError)

        return True

    def processExpr(self, stringExpr, databaseParticles, checkUnique=True):
        """
        Process a string expression (constraint or condition) for
        element weights. The element appearing in the string are stored in self.elementMap
        along with their string labels which are used to evaluate the expression.

        :param stringExpr: A mathematical expression for elements (e.g. 2*([[['jet']],[['jet']]]))
        :param databaseParticles: A Model object containing all the particle objects for the database.
        :param checkUnique: If True raises an error if the elements appearing in the expression
                            are not unique (relevant for avoiding double counting in the expression).

        :return: A function that takes a dictionary ({stringLabel : element}) and evaluates
                 the expression.
        """

        # Remove quotes, spaces and curly brackets from expression
        exprFunc = stringExpr[:].replace("'", "").replace('"', '')
        exprFunc = exprFunc.replace(" ", "")
        exprFunc = exprFunc.replace("}", "").replace("{", "")  # New format

        # Get the maximum element ID already used
        if self.elementMap:
            nel = max([el.elID for el in self.elementMap])
        else:
            self.elementMap = {}
            nel = 0

        # Get the elements contained in the expression
        for elStr in elementsInStr(str(stringExpr)):
            elObj = Element(elStr, self.finalState,
                            self.intermediateState,
                            model=databaseParticles,
                            sort=False)

            for el in self.elementMap:
                if elObj == el:
                    elObjLabel = self.elementMap[el]
                    if checkUnique:
                        msgError = "Duplicate elements found in: "
                        msgError += "%s in %s" % (stringExpr, self.globalInfo.id)
                        logger.error(msgError)
                        raise SModelSError(msgError)
                    break
            else:
                # If it is a new element, add it to the element -> label map
                nel += 1
                elObj.elID = nel
                elObjLabel = 'el_%i' % elObj.elID
                self.elementMap[elObj] = elObjLabel

            # Replace the element string by its label
            elStr = elStr.replace("'", "").replace(" ", "")
            exprFunc = exprFunc.replace(elStr, '%s' % elObjLabel)

        return exprFunc

    def evalConstraintFor(self, elements):
        """
        Evaluate the constraint function for a list of elements
        which have been matched to the txname elements. The
        elements must have the attribute txlabel assigned to
        the label appearing in the constraint expression.

        :param elements: List of Element objects

        :return: Value for the evaluated constraint, if the constraint
                 has been defined, None otherwise.
        """

        if not self._constraintFunc:
            return None

        # Build dictionary with elements and functions
        # required for evaluating the constraint expression
        localsDict = {"Cgtr": cGtr, "cGtr": cGtr, "cSim": cSim, "Csim": cSim}
        # Add a dictionary entry with zero weights for all element labels
        localsDict.update({elLabel: 0*physicsUnits.fb for elLabel in self.elementMap.values()})
        # Add elements weights
        for el in elements:
            localsDict[el.txlabel] += el.weight

        return eval(self._constraintFunc, localsDict, {})

    def evalConditionsFor(self, elements):
        """
        Evaluate the conditions for a list of elements
        which have been matched to the txname elements. The
        elements must have the attribute txlabel assigned to
        the label appearing in the constraint expression.

        :param elements: List of Element objects

        :return: List of condition values.
        """

        if not self._conditionsFunc:
            return None

        # Build dictionary with elements and functions
        # required for evaluating the constraint expression
        localsDict = {"Cgtr": cGtr, "cGtr": cGtr, "cSim": cSim, "Csim": cSim}
        # Add a dictionary entry with zero weights for all element labels
        localsDict.update({elLabel: 0*physicsUnits.fb for elLabel in self.elementMap.values()})
        # Add elements weights
        for el in elements:
            localsDict[el.txlabel] += el.weight

        conditions = [eval(cond, localsDict, {}) for cond in self._conditionsFunc]

        return conditions

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
        try:
            self.getDataMap(xDataPoint)
        except (IndexError) as e:
            msgError = "Error constructing dataMap for %s" % self.path
            msgError += ": %s" % str(e)
            logger.error(msgError)
            raise SModelSError(msgError)
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

        # Get length of flat array:
        nDim = max([arrayIndex for arrayIndex in self.dataMap.keys()])+1
        xFlat = [None]*nDim
        for arrayIndex in self.dataMap:
            # If arrayMap has not been defined, retrieve its value directly from x
            if self.arrayMap is None:
                _, attr, unit = self.dataMap[arrayIndex]
                xval = x[arrayIndex]
            # If it is a nested bracket use the arrayMap:
            else:
                multiIndex, attr, unit, _ = self.arrayMap[arrayIndex]
                i, j = multiIndex[:2]
                xval = x[i][j]
                if isinstance(xval, tuple):
                    xval = xval[multiIndex[2]]

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
            # Invert array map:
            mapInv = {v[0]: k for k, v in self.arrayMap.items()}
            # Get sorted array indices
            ijk = sorted(mapInv.keys())
            massPoint = []
            for i, j, k in ijk:
                while i >= len(massPoint):
                    massPoint.append([])
                while j >= len(massPoint[i]):
                    massPoint[i].append([])
                index = mapInv[(i, j, k)]
                _, attr, unit, _ = self.arrayMap[index]
                value = xFlat[index]
                if attr == 'totalwidth':
                    value = unscaleWidth(value)
                value = value*unit
                if k == 0:
                    massPoint[i][j] = value
                elif k == 1:
                    massPoint[i][j] = (massPoint[i][j], value)

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

        # If dataMap has already been defined, do nothing:
        if self.dataMap is not None:
            return

        # Since all elements are equivalent, use the first one
        # to define the map:
        el = list(self.elementMap.keys())[0]
        tree = el.tree

        # Get a nested array of nodes corresponding to the data point:
        nodeArray = []
        for mom, daughters in tree.dfs_successors().items():
            if mom.isInclusive:
                continue
            if mom != tree.root:  # Ignore PV
                nodeArray.append(mom)
            for d in daughters:
                if d.isSM:  # Ignore SM particles
                    continue
                if d.isInclusive:
                    continue
                if tree.out_degree(d) != 0:  # Ignore unstable daughters (will appear as mom)
                    continue
                nodeArray.append(d)

        # Iterate over the array and construct a map for the nodes
        # and the flat array:
        arrayMap = {}
        massIndex = 0  # Initial index for the masses
        widthIndex = len(nodeArray)  # Initial index for the widths
        for i, br in enumerate(massPoint):
            if br == '*':
                continue  # Skip inclusive branch
            for j, m in enumerate(br):
                node = nodeArray.pop(0)
                arrayValue = m
                mass, massUnit, width, widthUnit = self.getDataEntry(arrayValue)
                # Add entry for mass
                if mass is not None:
                    arrayMap[massIndex] = ((i, j, 0), 'mass', massUnit, node.node)
                    massIndex += 1
                # Add entry for width
                if width is not None:
                    arrayMap[widthIndex] = ((i, j, 1), 'totalwidth', widthUnit, node.node)
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

        elif self.dataType == 'upperLimit':
            ul = self.txnameData.getValueFor(point)
            element._upperLimit = ul  # Store the upper limit for convenience
            if ul is None:
                eff = 0.  # Element is outside the grid or the decays do not correspond to the txname
            else:
                eff = 1.
        else:
            logger.error("Unknown txnameData type: %s" % self.dataType)
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

        # Get list of elements with matching canonical name
        # (takes into account inclusive names)
        for el, elLabel in self.elementMap.items():
            # Compare elements:
            cmp, sortedEl = el.compareTo(element)
            if cmp == 0:
                # Attach label used for evaluating expressions
                sortedEl.txlabel = elLabel
                return sortedEl

        # If this point was reached, there were no macthes
        return False

    def hasLikelihood(self):
        """
        Can I construct a likelihood for this map?
        True for all efficiency maps, and for upper limits maps
        with expected Values.
        """
        if self.dataType == "efficiencyMap":
            return True
        if self.txnameDataExp is not None:
            return True
        return False
