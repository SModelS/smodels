# !/usr/bin/env python3

"""
.. module:: txnameObj
   :synopsis: Holds the class for storing the simplified model information
              and data  for experimental results.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import os
from smodels.base import physicsUnits
from smodels.base.genericSMS import GenericSMS
from smodels.base.smodelsLogging import logger
from smodels.statistics.basicStats import observed, apriori, aposteriori, NllEvalType

from smodels.experiment.expAuxiliaryFuncs import (smsInStr, removeUnits,
                                               rescaleWidth, unscaleWidth,
                                               concatenateLines, cSim, cGtr,
                                               bracketToProcessStr)
from smodels.experiment.expSMS import ExpSMS
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.experiment.txnameDataObj import TxNameData
from smodels.experiment.reweighting import defaultEffReweight, defaultULReweight
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

    def __init__(self, path=None, globalObj=None, infoObj=None,
                 databaseParticles=None):
        self.path = path
        self.globalInfo = globalObj
        self._infoObj = infoObj
        self.txnameData = None
        self.txnameDataExp = None  # evaluationType Data
        self.dataMap = None
        self.axesMap = None
        self._arrayMap = None
        self._arrayToNodeDict = None
        self.smsMap = {}  # Stores the SMS and their label representaion
        self._constraintFunc = None
        self._conditionsList = []
        self.finalState = ['MET', 'MET']  # default final state
        self.intermediateState = None  # default intermediate state

        if self.path is None:
            return
        logger.debug(f'Creating object based on txname file: {self.path}')
        # Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error(f"Txname file {path} not found")
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
                logger.info(f"Duplicated field {tag} found in file {self.path}")
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

        # If necessary, convert constraints and conditions
        # from old bracket notation:
        self.convertBracketNotation()

        # Process constraint, simplify it so it can be easily evaluated and
        # stores the SMS in self.smsMap
        if hasattr(self, 'constraint'):            
            exprFunc, smsMap = self.processExpr(self.constraint,
                                               databaseParticles,
                                               checkUnique=True)
            self._constraintFunc = exprFunc
            self.smsMap = smsMap

        # Process conditions, simplify it so it can be easily evaluated and
        # stores the expressions and SMS maps in _conditionsList
        if hasattr(self, 'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds, list):
                conds = [conds]
            conditionsList = []
            for cond in conds:
                exprFunc, smsMap = self.processExpr(cond, databaseParticles)
                conditionsList.append({exprFunc: smsMap})
            self._conditionsList = conditionsList

        # Do consistency checks:
        self.checkConsistency()

        # Define canonical name (should be the same for all SMS)
        self.canonName = list(self.smsMap.keys())[0].canonName

        # Get detector size (if not found in self, look for it in datasetInfo or globalInfo).
        # If not defined anywhere, set it to None and default values will be used for reweighting.
        self.Leff_inner = self.fetchAttribute('Leff_inner', fillvalue=None)
        self.Leff_outer = self.fetchAttribute('Leff_outer', fillvalue=None)

        # Convert data to flat unitless arrays and create TxNameData instances
        # to handle this data
        x_values, y_values = self.preProcessData(data)
        self.txnameData = TxNameData(x=x_values, y=y_values, txdataId=ident)
        if expectedData:
            x_values, y_values = self.preProcessData(expectedData)
            self.txnameDataExp = TxNameData(x=x_values, y=y_values, txdataId=ident)

        # Finally convert axes fields used for validation
        self.convertAxes()

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
        Checks if all the SMS in txname have the same structure (topology/canonical name)
        and verify if its constraints and conditions are valid expressions.

        :return: True if txname is consitency. Raises an error otherwise.
        """

        # Check if all SMS have the same canonical name:
        smsList = list(self.smsMap.keys())
        cName = smsList[0].canonName
        for sms in smsList[1:]:
            if sms.canonName != cName:
                msgError = f"Txname {self} referes to SMS with distinct topologies."
                logger.error(msgError)
                raise SModelSError(msgError)

        # Check if the constraint can be evaluated:
        smsList = []
        for smsObj, smsLabel in self.smsMap.items():
            smsTest = smsObj.copy()
            smsTest.weight = 3.0*(len(smsList)+1)*physicsUnits.fb  # dumyy xsec
            smsTest.txlabel = smsLabel
            smsList.append(smsTest)

        # Dummy constraint eval:
        try:
            res = self.evalConstraintFor(smsList)
        except (TypeError, NameError) as e:
            msgError = f"Can not evaluate constraint {self._constraintFunc} for {self}."
            msgError += f"\n error: {str(e)}"
            logger.error(msgError)
            raise SModelSError(msgError)

        if res is not None and not isinstance(res, (float, int, unum.Unum,)):
            msgError = f"Constraint for {self} returned an invalid value: {res} ({type(res)})"
            logger.error(msgError)
            raise SModelSError(msgError)

        # Dummy condistions eval:
        smsList = []
        for cond in self._conditionsList:
            for smsMap in cond.values():
                for smsObj, smsLabel in smsMap.items():
                    smsTest = smsObj.copy()
                    smsTest.weight = 3.0*(len(smsList)+1)*physicsUnits.fb  # dumyy xsec
                    smsList.append(smsTest)
        try:
            res = self.evalConditionsFor(smsList)
        except (TypeError, NameError) as e:
            msgError = f"Can not evaluate conditions {self._conditionsFunc} for {self}."
            msgError += f"\n error: {str(e)}"
            logger.error(msgError)
            raise SModelSError(msgError)

        if res is not None and not isinstance(res, list):
            msgError = f"Conditions for {self} returned an invalid value: {res} ({type(res)})"
            logger.error(msgError)
            raise SModelSError(msgError)

        return True
    
    def convertBracketNotation(self):
        """
        If the old bracket notation has been found in contraints
        or conditions, convert the strings.
        The _arrayToNodeDict is also defined to keep track of matching
        between the original nested indices and the node indices in
        the new format. The old constraint and conditions are stored
        in self._constraint and self._conditions.
        """

        # Check if constraint needs conversion:
        if hasattr(self,'constraint'):
            self._constraint = str(self.constraint)[:]
            if '[' in self.constraint and ']' in self.constraint:                
                # Get the SMS contained in the expression
                # (remove single quotes to avoid ambiguity)   
                constraint = str(self.constraint)[:].replace("'","").replace(" ","")
                for smsStr in smsInStr(constraint):
                    newSMS,arrayMap = bracketToProcessStr(smsStr,
                                                          finalState=self.finalState,
                                                           intermediateState=self.intermediateState,
                                                           returnNodeDict=True)
                    smsStr = smsStr.replace("'","").replace(" ","")
                    constraint = constraint.replace(smsStr,'{%s}' %newSMS,1)
                    if  self._arrayToNodeDict is None:
                        self._arrayToNodeDict = arrayMap
                    elif self._arrayToNodeDict != arrayMap:
                        logger.error("The SMS in the constraint do not all have the same ordering and topology.")
                self.constraint = constraint
            
        # Check if conditions needs conversion:
        if hasattr(self,'condition') and self.condition:
            if not isinstance(self.condition, list):
                self.condition = [self.condition]
            self._condition = self.condition[:]
            for icond,cond in enumerate(self.condition):
                cond = str(cond)[:].replace("'","").replace(" ","")
                if not ('[' in cond and ']' in cond):
                    continue                
                # Get the SMS contained in the expression                
                for smsStr in smsInStr(cond):
                    newSMS = bracketToProcessStr(smsStr,
                                                    finalState=self.finalState,
                                                    intermediateState=self.intermediateState)
                    smsStr = smsStr.replace("'","").replace(" ","")
                    cond = cond.replace(smsStr, '{%s}' %newSMS,1)
                self.condition[icond] = cond

    def convertAxes(self):
        """
        If the axes field attribute has been defined (v2 format)
        convert it to a list of dictionaries with the format:
        {arrayIndex : axesStr},
        where arrayIndex refer to the index in the flat grid array
        (v3 format) and axesStr defines how this index should be mapped to the
        validation axes (e.g. {0 : 'x', 1 : '(x+y)/2', 2 : 'y', ...}).
        The axes attribute is replaced by ._axes.        
        """

        # Do nothing is axesMap has already been defined
        if hasattr(self,'axesMap') and self.axesMap is not None:
            return
        
        # Check if axes needs conversion:
        if hasattr(self,'axes') and self.axes:
            # Define variables so the expressions can be evaluated:
            from sympy import symbols
            x,y,z,w = symbols('x y z w')
            # In case only one axis has been define, convert it to a list
            if not isinstance(self.axes, list):
                self.axes = [self.axes]
            # Replace axes by _axes
            self._axes = self.axes[:]
            delattr(self,'axes')
            # Convert axes in bracket notation to a dictionary
            self.axesMap = []
            for ax in self._axes:
                try:
                    axList = eval(ax, {'x' : x, 'y' : y, 'z' : z, 'w' : w})
                except NameError:
                    logger.error(f"Error evaluating axis for txname: {self} ({ax})")
                    raise SModelSError()
                axMap = {}
                for flatArrayIndex,indexMap in self._arrayMap.items():
                    # Get the (i,j,k) coordinates for the nested bracket
                    # which corresponds to the new flat array format
                    # (k is only used if there are widths)
                    oldCoords = indexMap[0]
                    oldArrayValue = axList[oldCoords[0]][oldCoords[1]]
                    if isinstance(oldArrayValue,tuple):
                        oldArrayValue = oldArrayValue[oldCoords[2]]
                    axMap[flatArrayIndex] = str(oldArrayValue).replace(' ','')
                self.axesMap.append(axMap)   

    def processExpr(self, stringExpr, databaseParticles,
                    checkUnique=False):
        """
        Process a string expression (constraint or condition) for
        the SMS weights.
        Returns a simplified string expression, which can be readily evaluated using
        a dictionary mapping SMS labels to their weights. It also returns an
        SMS map (dictionary) with the SMS objects as keys and their labels
        (appearing in the simplified expression) as values.

        :param stringExpr: A mathematical expression for SMS weights (e.g. 2*([[['jet']],[['jet']]]))
        :param databaseParticles: A Model object containing all the particle objects for the database.
        :param checkUnique: If True raises an error if the SMS appearing in the expression
                            are not unique (relevant for avoiding double counting in the expression).

        :return: simplfied expression (str), smsMap (dict).
        """

        # Remove quotes, spaces and curly brackets from expression
        exprFunc = stringExpr[:].replace("'", "").replace('"', '')
        exprFunc = exprFunc.replace(" ", "")
        exprFunc = exprFunc.replace("}", "").replace("{", "")  # New format

        # Get the maximum SMS ID already used
        smsMap = {}
        nsms = 0

        # Get the SMS contained in the expression
        for smsStr in smsInStr(str(stringExpr)):
            smsObj = ExpSMS.from_string(smsStr,model=databaseParticles,
                                        finalState=self.finalState,
                                        intermediateState=self.intermediateState)

            if checkUnique and any(smsObj == sms for sms in smsMap):
                msgError = "Duplicate SMS found in: "
                msgError += f"{stringExpr} in {self.globalInfo.id}"
                logger.error(msgError)
                raise SModelSError(msgError)

            # Add a new SMS and label to the map:
            nsms += 1
            smsObj.smsID = nsms
            smsObjLabel = 'sms_%i' % smsObj.smsID
            smsMap[smsObj] = smsObjLabel

            # Replace the SMS string by its label
            smsStr = smsStr.replace("'", "").replace(" ", "")
            exprFunc = exprFunc.replace(smsStr, f'{smsObjLabel}')

        return exprFunc, smsMap

    def evalConstraintFor(self, smsList):
        """
        Evaluate the constraint function for a list of SMS
        which have been matched to the txname SMS. The
        SMS must have the attribute txlabel assigned to
        the label appearing in the constraint expression.

        :param smsList: List of SMS objects with txlabel and weight attributes

        :return: Value for the evaluated constraint, if the constraint
                 has been defined, None otherwise.
        """

        if not self._constraintFunc:
            return None

        # Build dictionary with SMS and functions
        # required for evaluating the constraint expression
        localsDict = {"Cgtr": cGtr, "cGtr": cGtr, "cSim": cSim, "Csim": cSim}
        # Add a dictionary entry with zero weights for all SMS labels
        localsDict.update({smsLabel: 0*physicsUnits.fb for smsLabel in self.smsMap.values()})
        # Add SMS weights
        for sms in smsList:
            localsDict[sms.txlabel] += sms.weight

        return eval(self._constraintFunc, localsDict, {})

    def evalConditionsFor(self, smsList):
        """
        Evaluate the conditions for a list of SMS
        which have been matched to the txname SMS. The
        SMS must have the attribute txlabel assigned to
        the label appearing in the constraint expression.

        :param smsList: List of SMS objects with txlabel and weight attributes

        :return: List of condition values.
        """

        if not self._conditionsList:
            return None

        conditions = []
        for cond in self._conditionsList:
            condExpr = list(cond.keys())[0]
            smsMap = list(cond.values())[0]
            weightsMap = {smsLabel: 0.0*physicsUnits.fb for smsLabel in smsMap.values()}
            # Add weights for matching SMS:
            for sms in smsList:
                for smsC, smsLabel in smsMap.items():
                    if smsC == sms:
                        weightsMap[smsLabel] += sms.weight
            # Build dict with needed functions
            localsDict = {"Cgtr": cGtr, "cGtr": cGtr, "cSim": cSim, "Csim": cSim}
            # Add weightMap:
            localsDict.update(weightsMap)
            # Evaluate condition:
            conditions.append(eval(condExpr, localsDict, {}))

        return conditions

    def preProcessData(self, rawData):
        """
        Convert input data (from the upperLimits, evaluationTypeUpperLimits or efficiencyMap fields)
        to a flat array without units. The output is used to construct the TxNameData object,
        which will further process the data and interpolate it.
        It also builds the dictionary for translating SMS properties to the flat data array.

        :parameter rawData: Raw data (either string or list)

        :return: Two flat lists of data, one for the model parameters and the other for the y values
                 (UL or efficiency values)
        """

        if isinstance(rawData, str):
            data = self.evaluateString(rawData)
        else:
            data = rawData

        if len(data) == 0:
            logger.error(f"no data values for {self} found")
            raise SModelSError(f"no data values for {self} found")
        elif len(data[0]) < 2:
            logger.error(f"No valid data found for {self}")
            raise SModelSError(f"No valid data found for {self}")

        xDataPoint = data[0][0]
        yDataPoint = data[0][1]
        # Store y-unit:
        self.y_unit = removeUnits(yDataPoint, returnUnit=True)[1]
        if not isinstance(self.y_unit, (unum.Unum, float)):
            raise SModelSError(f"Error obtaining units from value: {data[0]} ")

        # Define graph->data mapping
        try:
            self.setDataMap(xDataPoint)
        except (IndexError) as e:
            msgError = f"Error constructing dataMap for {self.path}"
            msgError += f": {str(e)}"
            logger.error(msgError)
            raise SModelSError(msgError)
        # Transform data:
        x_values, y_values = self.transformData(data)

        return x_values, y_values

    def transformData(self, data):
        """
        Uses the information in self.dataMap (or self._arrayMap) to convert data
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
        self._arrayMap.

        :parameter x: A list (or nested list) with mass/width values.

        :return: A flat and unitless list matching sel.dataMap.
        """

        # Get length of flat array:
        nDim = max([arrayIndex for arrayIndex in self.dataMap.keys()])+1
        xFlat = [None]*nDim
        for arrayIndex in self.dataMap:
            # If arrayMap has not been defined, retrieve its value directly from x
            if self._arrayMap is None:
                _, attr, unit = self.dataMap[arrayIndex]
                xval = x[arrayIndex]
            # If it is a nested bracket use the arrayMap:
            else:
                multiIndex, attr, unit, _ = self._arrayMap[arrayIndex]
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
            logger.error(f"Error transforming point {x}")
            raise SModelSError()

        return xFlat

    def inverseTransformPoint(self, xFlat):
        """
        Transforms a 1D unitless array to a list of mass/width values.
        If self._arrayMap is defined, use it to convert to a nested
        bracket array foramt (e.g. [[mass1,(mass2,width2)],[mass3,mass4]]),
        otherwise convert it to a flat array (e.g. [mass1,mass2,mass3,mass4,width2])
        using self.dataMap.

        :parameter x: A 1D unitless array containing masses and rescaled widths

        :return: list (or nested list) with mass/width values (with units).
        """

        mLength = max(self.dataMap.keys())+1
        if self._arrayMap is not None:  # Convert to nested bracket
            # Invert array map:
            mapInv = {v[0]: k for k, v in self._arrayMap.items()}
            # Get sorted array indices
            ijk = sorted(mapInv.keys())
            massPoint = []
            for i, j, k in ijk:
                while i >= len(massPoint):
                    massPoint.append([])
                while j >= len(massPoint[i]):
                    massPoint[i].append([])
                index = mapInv[(i, j, k)]
                _, attr, unit, _ = self._arrayMap[index]
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

    def setDataMap(self, dataPoint):
        """
        If self.dataMap has not been defined,
        sets the dataMap using the first sms in self.smsMap.
        The dataMap is a dictionary mapping the node index to a flat data array.

        :param dataPoint: A point with the x-values from the data grid
                          (e.g. [[100*GeV,(50*GeV,1e-3*GeV)],[100*GeV,(50*GeV,1e-3*GeV),10*GeV]])

        :return: Dictionary with the data mapping {dataArrayIndex : (nodeNumber,attr,unit)}
                (e.g. {0  : (1,'mass',GeV), 1 : (1, 'totalwidth',GeV),...})
        """

        # If dataMap has already been defined, do nothing:
        if self.dataMap is not None:
            return
        # If it was not defined, it should be built from the old
        # nested array format information (generated by convertBracketNotation)
        if not hasattr(self,'_arrayToNodeDict') or not self._arrayToNodeDict:
            errorMsg = f"No information about the data has been found for {self.path}."
            errorMsg += " If the constraint was given using the process string"
            errorMsg += " (e.g. PV > ...) format, a dataMap has to be defined."
            logger.error(errorMsg)
            raise SModelSError(errorMsg)


        # Use the predefined map between branchIndex,vertexIndex -> nodeIndex
        # to define the dataMap

        # First identify which entries have mass and mass and width:
        massData = []
        widthData = []
        for i, br in enumerate(dataPoint):
            if br == '*':
                continue  # Skip inclusive branch
            for j, m in enumerate(br):
                nodeIndex = self._arrayToNodeDict[(i,j)]
                arrayValue = m
                mass, massUnit, width, widthUnit = self.getDataEntry(arrayValue)
                # Add entry for mass
                if mass is not None:
                    massData.append(((i, j, 0), 'mass', massUnit, nodeIndex))
                # Add entry for width
                if width is not None:
                    widthData.append(((i, j, 1), 'totalwidth', widthUnit, nodeIndex))

        # Now combine massData and width data in a single list, with
        # the widths at the end:
        massWidthData = massData + widthData

        # Finally generate an arrayMap with the flat array index as keys
        # and the (ijk,attribute,unit,nodeIndex) tuples as values:
        self._arrayMap = {i : nestedTuple for i,nestedTuple in enumerate(massWidthData)}
        # Also store graph <-> flat array map
        self.dataMap = {}
        for key, val in self._arrayMap.items():
            _, attr, unit, nodeIndex = val
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
            logger.error(f"Can not convert array value {arrayValue} ")
            raise SModelSError()

        return mass, massUnit, width, widthUnit

    def getDataFromSMS(self, sms):

        dataMap = self.dataMap
        smsData = [None]*(1+max(dataMap.keys()))
        for indexArray, nodeTuple in dataMap.items():
            nodeIndex, attr, unit = nodeTuple
            node = sms.indexToNode(nodeIndex)
            value = getattr(node, attr)
            if isinstance(unit, unum.Unum):
                value = value.asNumber(unit)
            if attr == 'totalwidth':
                value = rescaleWidth(value)
            smsData[indexArray] = value

        return smsData

    def getReweightingFor(self, sms):
        """
        Compute the lifetime reweighting for the SMS (fraction of prompt decays).
        If sms is a list, return 1.0.

        :param sms: SMS object

        :return: Reweighting factor (float)
        """

        if not isinstance(sms, GenericSMS):
            msgError = "Input of getReweightingFor must be an SMS object"
            msgError += f" and not {str(type(sms))}"
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
        for nodeIndex in sms.dfsIndexIterator(nodeIndex=sms.rootIndex):
            if nodeIndex in widthsInData:
                continue  # Ignore node if its width does not need reweighting

            # Convert to node object:
            node = sms.indexToNode(nodeIndex)
            if node.isSM:
                continue  # Ignore SM particles

            if sms.out_degree(nodeIndex) != 0:
                unstableWidths.append(node.totalwidth)
            else:
                stableWidths.append(node.totalwidth)

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
            raise SModelSError(f"Data should be in string format. Format {type(value)} found")

        try:
            val = eval(value, unitsDict)
        except (NameError, ValueError, SyntaxError):
            raise SModelSError(f"data string malformed: {value}")

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

    def getULFor(self, sms, evaluationType : NllEvalType = observed, mass=None):
        """
        Returns the upper limit (or expected) for SMS (only for upperLimit-type).
        Includes the lifetime reweighting (ul/reweight).
        If called for efficiencyMap results raises an error.
        If SMS is not defined, but mass is given, compute the UL using only the mass array
        (no width reweighting is applied) and the mass format is assumed
        to follow the evaluationType by the data.


        :param sms: SMS object or mass array (with units)
        :param expected: look in self.txnameDataExp, not self.txnameData
        """

        if not self.dataType == 'upperLimit':
            logger.error("getULFor method can only be used in UL-type data.")
            raise SModelSError()

        # If an SMS has been given, extract the data from it,
        # otherwise use the mass array
        if sms is not None:
            point = self.getDataFromSMS(sms)
            reweightF = self.getReweightingFor(sms)
            # If reweight factor is None (eff = 0),
            # return None
            if reweightF is None:
                return None
        else:
            massFlat = np.array(mass,dtype=object).flatten()
            point = [m.asNumber(physicsUnits.GeV) if isinstance(m,unum.Unum)
                     else m for m in massFlat]
            reweightF = 1.0

        if evaluationType == observed:
            ul = self.txnameData.getValueFor(point)
        else:
            if not self.txnameDataExp:
                return None
            else:
                ul = self.txnameDataExp.getValueFor(point)

        if ul is None:
            return None           

        ul = ul*reweightF*self.y_unit  # Add unit

        return ul

    def getEfficiencyFor(self, sms, mass=None):
        """
        For upper limit results, checks if the input SMS falls inside the
        upper limit grid and has a non-zero reweigthing factor.
        If it does, returns efficiency = 1, else returns
        efficiency = 0.  For efficiency map results, returns the
        signal efficiency including the lifetime reweighting.
        If a mass array is given as input, no lifetime reweighting will be applied.

        :param sms: SMS object.
        :return: efficiency (float)
        """

        
        if self.dataType == 'efficiencyMap':
            # If SMS has been given, extract data from it
            if sms is not None:
                # Get flat data from sms:
                point = self.getDataFromSMS(sms)
                # Compute reweighting factor
                reweightF = self.getReweightingFor(sms)
                if not reweightF:
                    return reweightF
            else:
                # If SMS is None, get data from mass array
                massFlat = np.array(mass,dtype=object).flatten()
                point = [m.asNumber(physicsUnits.GeV) if isinstance(m,unum.Unum)
                        else m for m in massFlat]
                reweightF = 1.0

            eff = self.txnameData.getValueFor(point)
            if not eff or math.isnan(eff):
                eff = 0.  # SMS is outside the grid or has zero efficiency
            eff = eff*reweightF*self.y_unit  # (unit should be 1)

        elif self.dataType == 'upperLimit':
            ul = self.getULFor(sms,mass=mass)
            if sms is not None:
                sms._upperLimit = ul  # Store the upper limit for convenience
            if ul is None:
                eff = 0.  # SMS is outside the grid or the decays do not correspond to the txname
            else:
                eff = 1.
        else:
            logger.error(f"Unknown txnameData type: {self.dataType}")
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

    def getInfo(self, infoLabel):
        """
        Returns the value of info field.

        :param infoLabel: label of the info field (string). It must be an attribute of
                          the TxNameInfo object
        """

        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False

    def hasSMSas(self, theorySMS, useLabel=None):
        """
        Verify if any SMS in conditions or constraint matches sms.
        If possible, check for both branch orderings (for two-branch SMS)
        and return the one with the largest data/reweighting factor.

        :param theorySMS: SMS object
        :param useLabel: String specifying the smsLabel to be used. If None, checks
                         for all ExpSMS in self.smsMap.
        :return: A copy of the sms with its nodes sorted according to
                 the matching topology in the TxName. Nodes matching InclusiveNodes
                 or inclusiveLists are replaced.
        """

        if useLabel:
            checkList = [(expSMS, smsLabel) for expSMS, smsLabel in self.smsMap.items() 
                                            if smsLabel == useLabel]
        else:
            checkList = [(expSMS, smsLabel) for expSMS, smsLabel in self.smsMap.items()]
        
        for expSMS, smsLabel in checkList:
            # Compare sms:
            matches = []
            for sms in [theorySMS,theorySMS.switchBranches()]:
                if sms is None: # If None it means we can not switch branches
                    continue
                matchedSMS = expSMS.matchesTo(sms)
                if matchedSMS is None: # If None it means no matches are possible
                    break
                matchedSMS.txlabel = smsLabel
                matches.append(matchedSMS)
            if not matches:
                continue
            elif len(matches) == 1:
                return matches[0]
            else: # sort matches by point and reweigthing factor
                matches = sorted(matches, 
                                 key = lambda sms: (self.getDataFromSMS(sms),
                                                    self.getReweightingFor(sms)),
                                 reverse=True)
                return matches[0]

        # If this point was reached, there were no macthes
        return None

    def hasLikelihood(self):
        """
        Can I construct a likelihood for this map?
        True for all efficiency maps, and for upper limits maps
        with evaluationType Values.
        """
        if self.dataType == "efficiencyMap":
            return True
        if self.txnameDataExp is not None:
            return True
        return False
        return False
