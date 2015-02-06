"""
.. module:: dataObjects
   :synopsis: Holds the classes and methods used to read and store the information in the
              sms.py files. Also contains the interpolation methods

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import logging,os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.theory.particleNames import elementsInStr
from smodels.theory.element import Element
from scipy.interpolate import griddata
from scipy.linalg import svd
import numpy as np
import unum
import copy

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.DEBUG)


class TxName(object):
    """Holds the information related to one txname in the Txname.txt
    file (constraint, condition,...) as well as the data.
    """
    
    def __init__(self, path):
        self.txnameFile = path
        self.txnameData = None
        self._elements = []
        
        logger.debug('Creating object based on txname file: %s' %self.txnameFile)        
 
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Txname file %s not found" % path)
            sys.exit()      
        txfile = open(self.txnameFile)
        content = txfile.readlines()
        txfile.close()
        
        #Get tags in info file:
        tags = [line.split(':', 1)[0].strip() for line in content]
        for i,tag in enumerate(tags):
            if not tag: continue
            line = content[i]
            value = line.split(':',1)[1].strip()            
            if tags.count(tag) == 1:
                if ';' in value: value = value.split(';')
                if tag == 'upperLimits' or tag == 'efficiencyMap':
                    self.txnameData = TxNameData(tag,value)
                else: self.addInfo(tag,value)
            else:
                logger.info("Ignoring unknown field %s found in file %s" % (tag, self.infopath))
                continue
        
        #Builds up a list of _elements appearing in constraints:        
        if hasattr(self,'constraint'):
            self._elements = [Element(el) for el in elementsInStr(self.constraint)]

    def __str__(self):
        return self.txname
        
    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        :param tag: information label (string)
        :param value: value for the field in string format 
        """

        setattr(self,tag,value)
    
    def getInfo(self, infoLabel):
        """Returns the value of info field.
        :param infoLabel: label of the info field (string). It must be an attribute of
                          the TxNameInfo object
        """
        
        if hasattr(self,infoLabel): return getattr(self,infoLabel)
        else: return False

    def getEfficiencyFor(self,element):
        """
        For upper limit results, checks if the input element appears in the constraints
        and falls inside the upper limit grid.
        If it does, returns efficiency = 1, else returns efficiency = 0.
        For efficiency map results, checks if the input element appears in the constraints.
        and falls inside the efficiency map grid.
        If it does, returns the corresponding efficiency value, else returns efficiency = 0.
        
        :param element: Element object
        :return: efficiency (float)
        """
        
        if self.txnameData.type == 'upperLimits':
            for el in self._elements:                
                if element.particlesMatch(el):
                    ul = self.txnameData.getValueFor(el.getMasses())
                    if type(ul) == type(fb): return 1.
            return 0.
        elif self.txnameData.type == 'efficiencyMap':
            for el in self._elements:
                if element.particlesMatch(el):
                    eff = self.txnameData.getValueFor(element.getMasses())
                    if eff: return eff                    
            return 0.
        else:
            logger.error("Unknown data type: %s" % self.txnameData.type)
            sys.exit()
            
class TxNameData(object):
    """Holds the data for the Txname object.  It holds Upper limit values or efficiencies."""
    
   
    def __init__(self,tag,value, accept_errors_upto=.05):
        """
        :param tag: data tag in string format (upperLimits or efficiencyMap)
        :param value: data in string format
        :param accept_errors_upto: If None, do not allow extrapolations outside of convex hull.
            If float value given, allow that much relative uncertainty on the upper limit / efficiency
            when extrapolating outside convex hull.
            This method can be used to loosen the equal branches assumption.
        """
        self.type = tag
        if type(value)==str:
            self.data = eval(value, {'fb' : fb, 'pb' : pb, 'GeV' : GeV, 'TeV' : TeV})
        else: ## the data can also be given as lists, for debugging
            self.data = value
        self.unit = 1.0 ## store the unit so that we can take arbitrary units for the "z" values.
                        ## default is unitless, which we use for efficiency maps
        if len(self.data)<1 or len(self.data[0])<2:
                logger.error ( "input data not in correct format. expecting sth like " \
         " [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], ... ] for upper " \
         " limits or [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .1 ], ... ] for efficiency maps" )
        if type(self.data[0][1])==unum.Unum:
            ## if its a unum, we store 1.0 * unit
            self.unit=self.data[0][1] / ( self.data[0][1].asNumber() )
        self.unit= self.data[0][1] / ( self.data[0][1].asNumber() )
        self.accept_errors_upto=accept_errors_upto
        self.computeV()

       
    def getValueFor(self,massarray):
        """
        Interpolates the data and returns the UL or efficiency for the respective massarray
        :param massarray: mass array values (with units), i.e. [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        """
        m=self.flattenMassArray ( massarray ) ## flatten
        mrot=np.dot(m,self.V)  ## rotate
        dp=self.countNonZeros ( mrot )
        if dp != self.dimensionality: ## we have data in different dimensions
            if self.accept_errors_upto == None:
                return float('nan')*self.unit
            logger.info ( "attempting to interpolate outside of convex hull (d=%d,dp=%d)" %
                     ( self.dimensionality, dp ) )
            return self._interpolateOutsideConvexHull ( massarray )
        r = griddata( self.Mp, self.xsec, mrot[:self.dimensionality], method="linear")
        return r[0]*self.unit
        
    def flattenMassArray ( self, data ):
        """ flatten mass array and remove units """
        ret=[]
        for i in data:
            for j in i:
                ret.append ( j.asNumber(GeV) )
        return ret
    def _estimateExtrapolationError ( self, massarray ):
        """ when projecting a point p from n to the point P in m dimensions,
            we estimate the expected extrapolation error with the following strategy: 
            we compute the gradient at point P, and let alpha be the distance between
            p and P. We then walk one step of length alpha in the direction of the greatest ascent,
            and the opposite direction. Whichever relative change is greater is 
            reported as the expected extrapolation error.
        """
        p=self.flattenMassArray ( massarray ) ## point p in n dimensions
        P=np.dot(p,self.V)                    ## projected point p in n dimensions
        ## P[self.dimensionality:] is project point p in m dimensions
        # m=self.countNonZeros ( P ) ## dimensionality of input
        ## how far are we away from the "plane": distance alpha
        alpha = np.sqrt ( np.dot ( P[self.dimensionality:], P[self.dimensionality:] ) )
        ## the value of the grid at the point projected to the "plane"
        projected_value=griddata( self.Mp, self.xsec, P[:self.dimensionality], method="linear")[0]
        
        ## compute gradient
        gradient=[]
        for i in range ( self.dimensionality ):
            P2=copy.deepcopy(P)
            P2[i]+=alpha
            gradient.append ( ( 
                griddata( self.Mp, self.xsec, P2[:self.dimensionality], method="linear")[0] - projected_value ) / alpha )
        ## normalize gradient
        # print "gradient=",gradient
        C= np.sqrt ( np.dot ( gradient, gradient ) )
        for i,j in enumerate(gradient):
            gradient[i]=gradient[i]/C*alpha
        ## walk one alpha along gradient
        P3=copy.deepcopy(P)
        P4=copy.deepcopy(P)
        for i,j in enumerate(gradient):
            P3[i]+=gradient[i]
            P4[i]-=gradient[i]
        # print "projected value", projected_value
        ag=griddata( self.Mp, self.xsec, P3[:self.dimensionality], method="linear")[0]
        #print "along gradient", ag
        agm=griddata( self.Mp, self.xsec, P4[:self.dimensionality], method="linear")[0]
        #print "along negative gradient",agm
        dep=abs ( ag - projected_value ) / projected_value
        dem=abs ( agm - projected_value ) / projected_value
        de=dep
        if dem > de: de=dem
        return de

    def _interpolateOutsideConvexHull ( self, massarray ):
        """ experimental routine, meant to check if we can interpolate outside convex hull """
        p=self.flattenMassArray ( massarray )
        P=np.dot(p,self.V)
        projected_value=griddata( self.Mp, self.xsec, P[:self.dimensionality], method="linear")[0]
        de = self._estimateExtrapolationError ( massarray ) 
        if de < self.accept_errors_upto:
            return projected_value * self.unit
        logger.info ( "Expected error of %f too large to propagate outside convext hull" % de )
        return float("nan") * self.unit


    def countNonZeros ( self, mp ):
        """ count the nonzeros in a vector """
        nz=0
        for i in mp:
            if abs(i)>10**-5:
                nz+=1
        return nz

    def computeV ( self ):
        """ compute rotation matrix V, rotate and truncate also
            'data' points and store in self.Mp """
        M=[]
        self.xsec=[]
        for x,y in self.data:
            self.xsec.append ( y / self.unit )
            xp = self.flattenMassArray ( x )
            M.append ( xp )
        U,s,Vt=svd(M)
        V=Vt.T
        self.V=V
        Mp=[]

        self.dimensionality=0
        for m in M:
            mp=np.dot(m,V)
            Mp.append ( mp )
            nz=self.countNonZeros ( mp )
            if nz>self.dimensionality:
                self.dimensionality=nz
        self.MpCut=[]
        MpCut=[]
        for i in Mp:
            MpCut.append ( i[:self.dimensionality].tolist() )
        self.Mp=MpCut ## also keep the rotated points, with truncated zeros
        
        
