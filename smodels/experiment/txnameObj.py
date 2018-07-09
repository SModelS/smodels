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
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.theory.particleNames import elementsInStr
from smodels.tools.stringTools import concatenateLines
from smodels.theory.element import Element
from smodels.theory.topology import TopologyList
from smodels.tools.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.caching import _memoize
from scipy.linalg import svd
from scipy.interpolate import interp1d
import scipy.spatial.qhull as qhull
import numpy as np
import unum
import copy
import math

from math import floor, log10

#Build a dictionary with defined units. It can be used to evaluate
#expressions containing units.
unitsDict = dict([[varname,varobj] for varname,varobj in physicsUnits.__dict__.items() 
                  if isinstance(varobj,unum.Unum)])


class TxName(object):
    """
    Holds the information related to one txname in the Txname.txt
    file (constraint, condition,...) as well as the data.
    """

    def __init__(self, path, globalObj, infoObj ):
        self.path = path
        self.globalInfo = globalObj
        self._infoObj = infoObj
        self.txnameData = None
        self.txnameDataExp = None ## expected Data
        self._topologyList = TopologyList()

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
            if tag == 'upperLimits' or tag == 'efficiencyMap':
                data = value
                dataType = tag
            elif tag == 'expectedUpperLimits':
                expectedData = value
                dataType = 'upperLimits'
            else:
                self.addInfo(tag,value)

        ident = self.globalInfo.id+":"+dataType[0]+":"+ str(self._infoObj.dataId)
        ident += ":" + self.txName
        self.txnameData = TxNameData(data, dataType, ident )
        if expectedData:
            self.txnameDataExp = TxNameData( expectedData, dataType, ident )

        #Builds up a list of elements appearing in constraints:
        if hasattr(self,'finalState'):
            finalState = self.finalState
        else:
            finalState = ["MET","MET"]        
        elements = []
        if hasattr(self,'constraint'):
            elements += [Element(el,finalState) for el in elementsInStr(str(self.constraint))]
        if hasattr(self,'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds,list): conds = [conds]
            for cond in conds:
                for el in elementsInStr(str(cond)):
                    newEl = Element(el,finalState)
                    if not newEl in elements: elements.append(newEl)

        # Builds up TopologyList with all the elements appearing in constraints
        # and conditions:
        for el in elements:
            self._topologyList.addElement(el)

    def hasOnlyZeroes ( self ):
        ozs = self.txnameData.onlyZeroValues()
        if self.txnameDataExp:
            e_ozs = self.txnameDataExp.onlyZeroValues()
            if ozs and e_ozs:
                return True
            if (ozs and not e_ozs) or (e_ozs and not ozs):
                logger.warning ( "%s is weird. One of the (expected, observed) results is zeroes-only, the other one isnt." )
                return False
        return ozs


    def __str__(self):
        return self.txName

    def __lt__ ( self, other ):
        """ sort by txName """
        return self.txName < other.txName

    def getValueFor(self,massarray,expected=False ):
        """ 
        Access txnameData and txnameDataExp to get value for 
        massarray.

        :param massarray: mass array values (with units), i.e.
                          [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        :param expected: query self.txnameDataExp
        """
        if not expected:
            return self.txnameData.getValueFor( massarray )
        else:
            if not self.txnameDataExp:
                return None
            else:
                return self.txnameDataExp.getValueFor( massarray )


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

        if tag == 'constraint' or tag == 'condition':
            if isinstance(value,list):
                value = [val.replace("'","") for val in value]
            else:
                value = value.replace("'","")
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
        Check both branch orderings.
        
        :param element: Element object
        :return: A copy of the element on the correct branch ordering appearing
                in the Txname constraint or condition.
        """

        for el in self._topologyList.getElements():
            if element.particlesMatch(el,branchOrder=True):
                return element.copy()
            else:
                elementB = element.switchBranches()
                if elementB.particlesMatch(el,branchOrder=True):
                    return elementB
        return False

    def hasLikelihood ( self ):
        """ can I construct a likelihood for this map? 
        True for all efficiency maps, and for upper limits maps
        with expected Values. """
        if self._infoObj.dataType == "efficiencyMap":
            return True
        if self.txnameDataExp != None:
            return True
        return False

    def getEfficiencyFor(self,mass):
        """
        For upper limit results, checks if the input mass falls inside the
        upper limit grid.  If it does, returns efficiency = 1, else returns
        efficiency = 0.  For efficiency map results, checks if the mass falls
        inside the efficiency map grid.  If it does, returns the corresponding
        efficiency value, else returns efficiency = 0.

        :param element: Element object
        :return: efficiency (float)
        """

        #Check if the element appears in Txname:
        val = self.txnameData.getValueFor(mass)
        if isinstance(val,unum.Unum):
            return 1.  #The element has an UL, return 1
        elif val is None or math.isnan(val):
            return 0.  #The element mass is outside the data grid
        elif isinstance(val,float):
            return val  #The element has an eff
        else:
            logger.error("Unknown txnameData value: %s" % (str(type(val))))
            raise SModelSError()

class TxNameData(object):
    """
    Holds the data for the Txname object.  It holds Upper limit values or efficiencies.
    """
    _keep_values = False ## keep the original values, only for debugging

    def __init__(self,value,datatag,Id,accept_errors_upto=.05):
        """
        :param value: values in string format
        :param dataTag: the dataTag (upperLimits or efficiencyMap)
        :param Id: an identifier, must be unique for each TxNameData!
        :param _accept_errors_upto: If None, do not allow extrapolations outside of
                convex hull.  If float value given, allow that much relative
                uncertainty on the upper limit / efficiency
                when extrapolating outside convex hull.
                This method can be used to loosen the equal branches assumption.
        """
        self.dataTag = datatag
        self._id = Id
        self._accept_errors_upto=accept_errors_upto
        self._V = None
        self.loadData(value)
        if self._keep_values:
            self.origdata = value

    def __str__ ( self ):
        """ a simple unique string identifier, mostly for _memoize """
        return str ( self._id )

    def round_to_n ( self, x, n ):
        if x==0.0:
            return x
        return round(x, int( -np.sign(x)* int(floor(log10(abs(x)))) + (n - 1)))

    def __ne__ ( self, other ):
        return not self.__eq__ ( other )

    def __eq__ ( self, other ):
        if type(self) != type ( other ):
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
        except:
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
        elif isinstance(value,dict):
            return dict([[self.getUnits(x),self.getUnits(y)] 
                                  for x,y in value.items()])
        elif isinstance(value,unum.Unum):
            #Check if value has unit or not:
            if not value._unit:
                return 1.
            #Now try to find stadandard unit which matches:
            for unit in stdUnits:
                y = (value/unit).normalize()
                if not y._unit:
                    return unit
            raise SModelSError("Could not find standard unit which matches %s. Using the standard units: %s" 
                               %(str(value),str(stdUnits)))
        else:
            return 1.    

    def removeUnits(self, value):
        """
        Remove units from unum objects. Uses the units defined
        in physicsUnits.standard units to normalize the data.
        
        :param value: Object containing units (e.g. [[100*GeV,100.*GeV],3.*pb])
        
        :return: Object normalized to standard units (e.g. [[100,100],3000])
        """
        
        stdUnits = physicsUnits.standardUnits
        
        if isinstance(value,list):
            return [self.removeUnits(x) for x in value]
        elif isinstance(value,dict):
            return dict([[self.removeUnits(x),self.removeUnits(y)] for x,y in value.items()])
        elif isinstance(value,unum.Unum):
            #Check if value has unit or not:
            if not value._unit:
                return value.asNumber()
            #Now try to normalize it by one of the standard pre-defined units:
            for unit in stdUnits:
                y = (value/unit).normalize()
                if not y._unit:
                    return value.asNumber(unit)
            raise SModelSError("Could not normalize unit value %s using the standard units: %s" 
                               %(str(value),str(stdUnits)))
        else:
            return value
        
    def getDataShape(self,value):
        """
        Stores the data format (mass shape) and store it for future use.
        If there are wildcards (mass or branch = None), store their positions.
        
        :param value: list of data points
        """

        if isinstance(value,list):
            return [self.getDataShape(m) for m in value]
        elif isinstance(value,(float,int,unum.Unum)):
            return type(value)
        else:
            return value

    def formatInput(self,value,shapeArray):
        """
        Format value according to the shape in shapeArray.
        If shapeArray contains entries = *, the corresponding entries
        in value will be ignored.
        
        :param value: Array to be formatted (e.g. [[200.,100.],[200.,100.]])
        :param shapeArray: Array with format info (e.f. ['*',[float,float]])
        
        :return: formatted array [[200.,100.]]
        
        """

        if shapeArray == '*':
            return None
        elif isinstance(value,list):
            if len(shapeArray) != len(value): 
                raise SModelSError("Input value and data shape mismatch (%s,%s)" 
                                   %(len(shapeArray),len(value)))
            return [self.formatInput(xi,shapeArray[i]) for i,xi in enumerate(value) 
                    if not self.formatInput(xi,shapeArray[i]) is None]
        else:
            return value
        
    def removeWildCards(self,value):
        """
        Remove all entries = '*' from value.
        
        :param value:  usually a list containing floats and '*' (e.g. [[200.,'*'],'*'],0.4],..)
        """
        
        if value == "*":
            return None
        elif isinstance(value,list):
            return [self.removeWildCards(v) for v in value if not self.removeWildCards(v) is None]
        else:
            return value
        
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
            
        self.units = self.getUnits(val)[0] #Store standard units
        self.dataShape = self.getDataShape(val[0][0]) #Store the data (mass) format (useful if there are wildcards)        
        self.value = self.removeUnits(val) #Remove units and store the normalization units
        self.value = self.removeWildCards(self.value)


        if len(self.value) < 1 or len(self.value[0]) < 2:
                raise SModelSError("input value not in correct format. expecting sth " \
                               "like [ [ [[ 300.*GeV,100.*GeV], "\
                               "[ 300.*GeV,100.*GeV] ], 10.*fb ], ... ] "\
                               "for upper limits or [ [ [[ 300.*GeV,100.*GeV],"\
                               " [ 300.*GeV,100.*GeV] ], .1 ], ... ] for "\
                               "efficiency maps. Received %s" % self.value[:80])

                
        if not isinstance(self.units[-1],unum.Unum) and not isinstance(self.units[-1],float):
            raise SModelSError("Error obtaining units from value: %s " %self.value[:80])


        self.y_values = np.array(self.value)[:,1]
        self.computeV()
        self.removeExtraZeroes()            
        self.cleanUp()

    @_memoize
    def getValueFor(self,massarray):
        """
        Interpolates the value and returns the UL or efficiency for the
        respective massarray
        
        :param massarray: mass array values (with units), i.e.
                          [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        """
        
        porig = self.removeUnits(massarray)
        porig = self.formatInput(porig,self.dataShape) #Remove entries which match wildcards
        porig = self.flattenArray(porig) ## flatten        
        self.massarray = massarray ## only for bookkeeping and better error msgs
        
        if len(porig) != self.full_dimensionality:
            logger.error ( "dimensional error. I have been asked to compare a "\
                    "%d-dimensional mass vector with %d-dimensional data!" % \
                    ( len(porig), self.full_dimensionality ) )
            return None

        p = ((np.matrix(porig)[0] - self.delta_x )).tolist()[0]
        P = np.dot(p,self._V)  ## rotate
        #Get value for the truncated point:        
        self.projected_value = self.interpolate(P[:self.dimensionality])
        
        #Check if input point has larger dimensionality:
        dp = self.countNonZeros(P)
        if dp > self.dimensionality: ## we have data in different dimensions
            if self._accept_errors_upto == None:
                return None
            logger.debug( "attempting to interpolate outside of convex hull "\
                    "(d=%d,dp=%d,masses=%s)" %
                     ( self.dimensionality, dp, str(massarray) ) )            
            return self._interpolateOutsideConvexHull(massarray)

        return self._returnProjectedValue()

    def flattenArray(self, objList):
        """
        Flatten any nested list to a 1D list.
        
        :param objList: Any list or nested list of objects (e.g. [[[100.,100.],1.],[[200.,200.],2.],..]
        
        :return: 1D list (e.g. [100.,100.,1.,200.,200.,2.,..])
        """
        
        ret = []
        
        for obj in objList:
            if isinstance(obj,list):
                ret.extend(self.flattenArray(obj))
            else:
                ret.append(obj)
        return ret        

    def interpolate(self, point, fill_value=np.nan):
        
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

    def _estimateExtrapolationError(self, massarray):
        """ when projecting a point p from n to the point P in m dimensions, we
            estimate the expected extrapolation error with the following
            strategy: we compute the gradient at point P, and let alpha be the
            distance between p and P. We then walk one step of length alpha in
            the direction of the greatest ascent, and the opposite direction.
            Whichever relative change is greater is reported as the expected
            extrapolation error.
        """
        
        porig = self.removeUnits(massarray)
        porig = self.formatInput(porig,self.dataShape) #Remove entries which match wildcards
        porig = self.flattenArray(porig)
         
        p = ((np.matrix(porig)[0] - self.delta_x)).tolist()[0]
        P=np.dot(p,self._V)                    ## projected point p in n dimensions
        ## P[self.dimensionality:] is project point p in m dimensions
        # m=self.countNonZeros ( P ) ## dimensionality of input
        ## how far are we away from the "plane": distance alpha
        alpha = float(np.sqrt( np.dot(P[self.dimensionality:],
                        P[self.dimensionality:])))
        if alpha == 0.:
            ## no distance to the plane, so no extrapolation error
            return 0.
        ## the value of the grid at the point projected to the "plane"

        ## compute gradient
        gradient=[]
        for i in range(self.dimensionality):
            P2=copy.deepcopy(P)
            P2[i]+=alpha
            pv = self.interpolate(P2[:self.dimensionality])
            g = float((pv - self.projected_value)/alpha)
            if math.isnan ( g ):
                ## if we cannot compute a gradient, we return nan
                return float("nan")
            gradient.append(g)
        ## normalize gradient
        C= float(np.sqrt( np.dot ( gradient, gradient ) ))
        if C == 0.:
            ## zero gradient? we return 0.
            return 0.
        for i,grad in enumerate(gradient):
            gradient[i]=grad/C*alpha
        ## walk one alpha along gradient
        P3=copy.deepcopy(P)
        P4=copy.deepcopy(P)
        for grad in gradient:
            P3[i]+= grad
            P4[i]-= grad
        agp=self.interpolate( P3[:self.dimensionality] )
        agm=self.interpolate( P4[:self.dimensionality] )
        dep,dem=0.,0.
        if self.projected_value == 0.:
            if agp!=0.:
                dep =1.0
            if agm!=0.:
                dem =1.0
        else:
            dep = abs( agp - self.projected_value)/self.projected_value
            dem = abs( agm - self.projected_value)/self.projected_value
        de=dep
        if dem > de: de=dem
        return de

    def _interpolateOutsideConvexHull(self, massarray):
        """ experimental routine, meant to check if we can interpolate outside
            convex hull """
            
        de = self._estimateExtrapolationError(massarray)
        
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
            logger.debug ( "projected value is None. Projected point not in " \
                    "convex hull? original point=%s" % self.massarray )
            return None
        
        #Set value to zero if it is lower than machine precision (avoids fake negative values)
        if abs(self.projected_value) < 100.*sys.float_info.epsilon:
            self.projected_value = 0.
        
        return self.projected_value*self.units[-1]

    def countNonZeros ( self, mp ):
        """ count the nonzeros in a vector """
        nz=0
        for i in mp:
            if abs(i)>10**-4:
                nz+=1
        return nz

    def onlyZeroValues ( self ):
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

    def computeV(self):
        """
        Compute rotation matrix _V, and triangulation self.tri
        
        """
        
        if not self._V is None:
            return
        
        Morig= [self.flattenArray(pt[0]) for pt in self.value]
        
        aM = np.matrix(Morig)
        MT = aM.T.tolist()
        self.delta_x = np.matrix([ sum (x)/len(Morig) for x in MT ])[0]
        M = []

        for Mx in Morig:
            m=(np.matrix(Mx) - self.delta_x).tolist()[0]
            M.append(m)

        try:
            ## we dont need thousands of points for SVD
            n = int(math.ceil(len(M)/2000.))
            Vt=svd(M[::n])[2]
        except Exception as e:
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
            Mp.append ( mp )
            nz=self.countNonZeros(mp)
            if nz>self.dimensionality:
                self.dimensionality=nz
        MpCut=[]
        for i in Mp:
            MpCut.append(i[:self.dimensionality].tolist() )

        if self.dimensionality > 1:
            self.tri = qhull.Delaunay(MpCut)
        else:            
            self.tri = Delaunay1D(MpCut)           
   
    def cleanUp(self):
        if self._keep_values:
            return
        if hasattr(self, "origdata"):
            del self.origdata
        
    def _getMassArrayFrom(self,pt,unit=physicsUnits.GeV):
        """
        Transforms the point pt in the PCA space to the original mass array
        :param pt: point with the dimentions of the data dimensionality (e.g. [x,y])
        :param unit: Unit for returning the mass array. If None, it will be
                     returned unitless
        :return: Mass array (e.g. [[mass1,mass2],[mass3,mass4]])
        """
        
        if self._V is None:
            logger.error("Data has not been loaded")
            return None
        if len(pt) != self.dimensionality:
            logger.error("Wrong point dimensions (%i), it should be %i" 
                         %(len(pt),self.dimensionality))
            return None
        fullpt = np.append(pt,[0.]*(self.full_dimensionality-len(pt)))
        mass = np.dot(self._V,fullpt) + self.delta_x
        mass = mass.reshape(self.massdim).tolist()
        if isinstance(unit,unum.Unum):
            mass = [[m*unit for m in br] for br in mass]

        return mass

     
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
    from smodels.tools.physicsUnits import GeV
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
    txnameData=TxNameData ( data, "upperLimits",  sys._getframe().f_code.co_name )
    t0=time.time()
    for masses in [ [[ 302.*GeV,123.*GeV], [ 302.*GeV,123.*GeV]],
                    [[ 254.*GeV,171.*GeV], [ 254.*GeV,170.*GeV]],
    ]:
        result=txnameData.getValueFor( masses )
        sm = "%.1f %.1f" % ( masses[0][0].asNumber(GeV), masses[0][1].asNumber(GeV) )
        print ( "%s %.3f fb" % ( sm, result.asNumber(fb) ) )
    print ( "%.2f ms" % ( (time.time()-t0)*1000. ) )
