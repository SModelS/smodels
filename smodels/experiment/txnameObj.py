#!/usr/bin/env python

"""
.. module:: txnameObj
   :synopsis: Holds the classes and methods used to read and store the
              information in the txname.txt files.
              Also contains the interpolation methods.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import logging,os,sys
from smodels.tools.physicsUnits import GeV, fb, TeV, pb
from smodels.theory.particleNames import elementsInStr
from smodels.tools.stringTools import concatenateLines
from smodels.theory.element import Element
from smodels.theory.topology import TopologyList
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.theory.auxiliaryFunctions import _memoize
from scipy.linalg import svd
import scipy.spatial.qhull as qhull
import numpy as np
import unum
import copy
import math
from math import floor, log10

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)
logger.setLevel(level=logging.ERROR)

round_to_n = lambda x, n: round(x, int( -np.sign(x)* int(floor(log10(abs(x)))) + (n - 1)))

class TxName(object):
    """Holds the information related to one txname in the Txname.txt
    file (constraint, condition,...) as well as the _data.
    """


    def __init__(self, path, infoObj):
        self.path = path
        self.globalInfo = infoObj
        self.txnameData = None
        self._topologyList = TopologyList()

        logger.debug('Creating object based on txname file: %s' %self.path)
        #Open the info file and get the information:
        if not os.path.isfile(path):
            logger.error("Txname file %s not found" % path)
            raise SModelSError()
        txtFile = open(path,'r')
        txdata = txtFile.read()
        if not "txName" in txdata: raise TypeError
        if not 'upperLimits' in txdata and not 'efficiencyMap' in txdata:
            raise TypeError
        txfile = open(self.path)
        content = concatenateLines (  txfile.readlines() )
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
                    self.txnameData = TxNameData(value)
                    self.txnameData.dataTag = tag
                else: self.addInfo(tag,value)
            else:
                logger.info("Ignoring unknown field %s found in file %s" \
                             % (tag, self.path))
                continue

        #Builds up a list of elements appearing in constraints:
        elements = []
        if hasattr(self,'constraint'):
            elements += [Element(el) for el in elementsInStr(self.constraint)]
        if hasattr(self,'condition') and self.condition:
            conds = self.condition
            if not isinstance(conds,list): conds = [conds]
            for cond in conds:
                for el in elementsInStr(cond):
                    newEl = Element(el)
                    if not newEl in elements: elements.append(newEl)

        # Builds up TopologyList with all the elements appearing in constraints
        # and conditions:
        for el in elements:
            el.sortBranches()
            self._topologyList.addElement(el)

    def __str__(self):
        return self.txName

    def __lt__ ( self, other ):
        """ sort by txName """
        return self.txName < other.txName

    def addInfo(self,tag,value):
        """
        Adds the info field labeled by tag with value value to the object.
        :param tag: information label (string)
        :param value: value for the field in string format
        """

        if tag == 'constraint' or tag == 'condition':
            if isinstance(value,list):
                value = [val.replace("'","") for val in value]
            else: value = value.replace("'","")

        try:
            setattr(self,tag,eval(value, {'fb' : fb, 'pb' : pb, 'GeV' : GeV, 'TeV' : TeV}))
        except SyntaxError:
            setattr(self,tag,value)
        except NameError:
            setattr(self,tag,value)
        except TypeError:
            setattr(self,tag,value)

    def getInfo(self, infoLabel):
        """Returns the value of info field.
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
        if type(val) == type(fb):
            return 1.  #The element has an UL, return 1
        elif val is None or math.isnan(val):
            return 0.  #The element mass is outside the _data grid
        elif type(val) == type(1.):
            return val  #The element has an eff
        else:
            logger.error("Unknown txnameData value: %s" % (str(type(val))))
            raise SModelSError()

class TxNameData(object):
    """ Holds the _data for the Txname object.  It holds Upper limit values or
        efficiencies."""

    def __init__(self,value,accept_errors_upto=.05):
        """

        :param value: _data in string format
        :param _accept_errors_upto: If None, do not allow extrapolations outside of
                convex hull.  If float value given, allow that much relative
                uncertainty on the upper limit / efficiency
                when extrapolating outside convex hull.
                This method can be used to loosen the equal branches assumption.
        """
        self._accept_errors_upto=accept_errors_upto
        self._store_value = value
        self._V = None
        self._data = None
        self.loadData()

    def __ne__ ( self, other ):
        return not self.__eq__ ( other )

    def __eq__ ( self, other ):
        return self._data == other._data

    def loadData(self):
        """
        Uses the information in _store_value to generate the _data grid used for
        interpolation.
        """

        if self._data:
            return

        if type(self._store_value)==str:
            self._data = eval(self._store_value,
                             {'fb':fb, 'pb':pb, 'GeV':GeV, 'TeV':TeV})
        else: ## the _data can also be given as lists, for debugging
            self._data = self._store_value
        self.unit = 1.0 ## store the unit so that we can take arbitrary units for
                        ## the "z" values.  default is unitless,
                        ## which we use for efficiency maps
        if len(self._data) < 1 or len(self._data[0]) < 2:
                logger.error ( "input _data not in correct format. expecting sth " \
                               "like [ [ [[ 300.*GeV,100.*GeV], "\
                               "[ 300.*GeV,100.*GeV] ], 10.*fb ], ... ] "\
                               "for upper limits or [ [ [[ 300.*GeV,100.*GeV],"\
                               " [ 300.*GeV,100.*GeV] ], .1 ], ... ] for "\
                               "efficiency maps" )
        if type(self._data[0][1])==unum.Unum:
            ## if its a unum, we store 1.0 * unit
            self.unit=self._data[0][1] / ( self._data[0][1].asNumber() )

        self.computeV()

    @_memoize
    def getValueFor(self,massarray):
        """
        Interpolates the _data and returns the UL or efficiency for the
        respective massarray
        :param massarray: mass array values (with units), i.e.
                          [[100*GeV,10*GeV],[100*GeV,10*GeV]]
        """

        self.loadData()

        porig=self.flattenMassArray ( massarray ) ## flatten
        self.massarray = massarray
        if len(porig)!=self.full_dimensionality:
            logger.error ( "dimensional error. I have been asked to compare a "\
                    "%d-dimensional mass vector with %d-dimensional _data!" % \
                    ( len(porig), self.full_dimensionality ) )
            return None
        p= ( (np.matrix(porig)[0] - self.delta_x ) ).tolist()[0]
        P=np.dot(p,self._V)  ## rotate
        dp=self.countNonZeros ( P )
        self.projected_value = self.interpolate( [ P[:self.dimensionality] ] )

        # self.projected_value = griddata( self.Mp, self.xsec, [ P[:self.dimensionality] ], method="linear")[0]
        # self.projected_value = float(self.projected_value)
        if dp != self.dimensionality: ## we have _data in different dimensions
            if self._accept_errors_upto == None:
                return None
            logger.debug ( "attempting to interpolate outside of convex hull "\
                    "(d=%d,dp=%d,masses=%s)" %
                     ( self.dimensionality, dp, str(massarray) ) )
            return self._interpolateOutsideConvexHull ( massarray )

        return self._returnProjectedValue()

    def flattenMassArray ( self, data ):
        """ flatten mass array and remove units """
        ret=[]
        for i in data:
            for j in i:
                ret.append ( j.asNumber(GeV) )
        return ret

    def interpolate(self, uvw, fill_value=np.nan):
        simplex = self.tri.find_simplex(uvw)
        vertices = np.take(self.tri.simplices, simplex, axis=0)
        temp = np.take(self.tri.transform, simplex, axis=0)
        d=temp.shape[2]
        delta = uvw - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
        self.vtx = vertices
        self.wts = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))
        v=self.xsec
        if type (self.xsec[0]) == float:
            values = np.array ( [ float(x) for x in self.xsec ] )
        else:
            values = np.array ( [ x.asNumber() for x in self.xsec ] )
        # print ("values=",(values[:]))
        # print ("v=",(v[:]))
        #print ("wts=",wts.shape)
        ret = np.einsum('nj,nj->n', np.take(values, self.vtx), self.wts)
        #print ("retb=",ret[0])
        with np.errstate(invalid='ignore'):
            ret[np.any(self.wts < -1e-10, axis=1)] = fill_value
        return float(ret[0])


    def _estimateExtrapolationError ( self, massarray ):
        """ when projecting a point p from n to the point P in m dimensions, we
            estimate the expected extrapolation error with the following
            strategy: we compute the gradient at point P, and let alpha be the
            distance between p and P. We then walk one step of length alpha in
            the direction of the greatest ascent, and the opposite direction.
            Whichever relative change is greater is reported as the expected
            extrapolation error.
        """
        #p=self.flattenMassArray ( massarray ) ## point p in n dimensions
        porig=self.flattenMassArray ( massarray ) ## flatten
        p= ( (np.matrix(porig)[0] - self.delta_x ) ).tolist()[0]
        P=np.dot(p,self._V)                    ## projected point p in n dimensions
        ## P[self.dimensionality:] is project point p in m dimensions
        # m=self.countNonZeros ( P ) ## dimensionality of input
        ## how far are we away from the "plane": distance alpha
        alpha = float ( np.sqrt ( np.dot ( P[self.dimensionality:],
                        P[self.dimensionality:] ) ) )
        ## the value of the grid at the point projected to the "plane"

        ## compute gradient
        gradient=[]
        for i in range ( self.dimensionality ):
            P2=copy.deepcopy(P)
            P2[i]+=alpha
            pv = self.interpolate ( [ P2[:self.dimensionality] ] )
            g=float ( ( pv - self.projected_value ) / alpha )
            #g=float ( ( griddata( self.Mp, self.xsec, [ P2[:self.dimensionality]],
            #            method="linear")[0] - self.projected_value ) / alpha )
            if math.isnan ( g ):
                ## if we cannot compute a gradient, we return nan
                return float("nan")
            gradient.append ( g )
        ## normalize gradient
        # print "gradient=",gradient
        C= float ( np.sqrt ( np.dot ( gradient, gradient ) ) )
        if C == 0.:
            ## zero gradient? we return 0.
            return 0.
        for i,j in enumerate(gradient):
            gradient[i]=gradient[i]/C*alpha
            #print "gradient after=",gradient[i]
        #print "^^^^^^"
        ## walk one alpha along gradient
        P3=copy.deepcopy(P)
        P4=copy.deepcopy(P)
        for i,j in enumerate(gradient):
            P3[i]+=gradient[i]
            P4[i]-=gradient[i]
        # print "projected value", projected_value
        agp=self.interpolate ( [ P3[:self.dimensionality] ] )
        #agp=griddata( self.Mp, self.xsec, [ P3[:self.dimensionality] ],
        #              method="linear")[0]
        #print "along gradient", ag
        agm=self.interpolate ( [ P4[:self.dimensionality] ] )
        #agm=griddata( self.Mp, self.xsec, [ P4[:self.dimensionality] ],
        #              method="linear")[0]
        #print "along negative gradient",agm
        dep,dem=0.,0.
        if self.projected_value == 0.:
            if agp!=0.:
                dep =1.0
            if agm!=0.:
                dem =1.0
        else:
            dep=abs ( agp - self.projected_value) / self.projected_value
            dem=abs ( agm - self.projected_value ) / self.projected_value
        de=dep
        if dem > de: de=dem
        return de

    def _interpolateOutsideConvexHull ( self, massarray ):
        """ experimental routine, meant to check if we can interpolate outside
            convex hull """
        porig=self.flattenMassArray ( massarray ) ## flatten
        p= ( (np.matrix(porig)[0] - self.delta_x ) ).tolist()[0]
        P=np.dot(p,self._V)
        de = self._estimateExtrapolationError ( massarray )
        if de < self._accept_errors_upto:
            return self._returnProjectedValue()
        if not math.isnan(de):
            logger.debug ( "Expected propagation error of %f too large to " \
                           "propagate." % de )
        return None

    def _returnProjectedValue ( self ):
        ## None is returned without units'
        if self.projected_value is None or math.isnan(self.projected_value):
            logger.debug ( "projected value is None. Projected point not in " \
                    "convex hull? original point=%s" % self.massarray )
            return None
        return self.projected_value * self.unit

    def countNonZeros ( self, mp ):
        """ count the nonzeros in a vector """
        nz=0
        for i in mp:
            if abs(i)>10**-5:
                nz+=1
        return nz

    def computeV ( self ):
        """ compute rotation matrix _V, rotate and truncate also
            '_data' points and store in self.Mp """
        if self._V!=None:
             return
        Morig=[]
        self.xsec=[]

        for x,y in self._data:
            self.xsec.append ( y / self.unit )
            xp = self.flattenMassArray ( x )
            Morig.append ( xp )
        aM=np.matrix ( Morig )
        MT=aM.T.tolist()
        self.delta_x = np.matrix ( [ sum (x)/len(Morig) for x in MT ] )[0]
        M = []


        for Mx in Morig:
            m=( np.matrix ( Mx ) - self.delta_x ).tolist()[0]
            M.append ( [ round_to_n ( x, 7 ) for x in m ] )

        U,s,Vt=svd(M)
        V=Vt.T
        self._V=V
        Mp=[]

        ## the dimensionality of the whole mass space, disrespecting equal branches
        ## assumption
        self.full_dimensionality = len(xp)
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
        self.tri = qhull.Delaunay( self.Mp )

if __name__ == "__main__":
    import time
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
    txnameData=TxNameData ( data ) ## "upperlimit", _data )
    t0=time.time()
    for masses in [ [[ 302.*GeV,123.*GeV], [ 302.*GeV,123.*GeV]],
                    [[ 254.*GeV,171.*GeV], [ 254.*GeV,170.*GeV]],
    ]:
        result=txnameData.getValueFor( masses )
        sm = "%.1f %.1f" % ( masses[0][0].asNumber(GeV), masses[0][1].asNumber(GeV) )
        print ( "%s %.3f fb" % ( sm, result.asNumber(fb) ) )
    print ( "%.2f ms" % ( (time.time()-t0)*1000. ) )
