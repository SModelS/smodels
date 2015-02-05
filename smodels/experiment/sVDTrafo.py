"""
.. module:: sVDTrafo
   :synopsis: Holds all methods needed to perform the singular value
       decomposition.  Before interpolating the experimental data, we perform a
       principal component analysis on the data, then check in how many dimensions
       we really want to interpolate. E.g. when interpolating two masses times two
       branches = four dimensions, but the branches are always equal, we really
       want to interpolate in two dimensions.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import numpy as np
from scipy.linalg import svd
from scipy.interpolate import griddata
from smodels.tools.physicsUnits import GeV, fb
import unum
import logging
import copy
import math

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.WARNING )

class SVDTrafo:
    def countNonZeros ( self, mp ):
        """ count the nonzeros in a vector """
        nz=0
        for i in mp:
            if abs(i)>10**-5:
                nz+=1
        return nz

    def flattenMassArray ( self, data ):
        """ flatten mass array and remove units """
        ret=[]
        for i in data:
            for j in i:
                ret.append ( j.asNumber(GeV) )
        return ret
                
    def computeV ( self, data ):
        """ compute the rotation matrix V, compute also rotated points """
        M,Mp=[],[]
        self.xsec=[]
        for x,y in data:
            self.xsec.append ( y / self.unit )
            xp = self.flattenMassArray ( x )
            M.append ( xp )
        # print "M=",M
        U,s,Vt=svd(M)
        V=Vt.T
        # print "V=",V
        self.Vt=Vt ## rotation matrix Vt
        self.V=V   ## Transposed rotation matrix V
        self.M=M ## the points in the original space
        nonzeros=0
        for m in M:
            mp=np.dot(m,V)
            Mp.append ( mp )
            nz=self.countNonZeros ( mp )
            if nz > nonzeros:
                nonzeros=nz
        self.dimensionality = nonzeros
        MpCut=[]
        for i in Mp:
            MpCut.append ( i[:self.dimensionality].tolist() )
        self.Mp=MpCut ## also keep the rotated points, with truncated zeros

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
        return float("nan") * self.unit
            
    def getInterpolatedValue ( self, massarray ):
        """
           :param massarray: e.g. [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ]
           :return: interpolated value, in same units as input
        """
        m=self.flattenMassArray ( massarray )
        mrot=np.dot(m,self.V)
        dp=self.countNonZeros ( mrot )
        if dp != self.dimensionality:
            if self.accept_errors_upto == None:
                return float('nan')*self.unit
            logger.info ( "attempting to interpolate outside of convex hull (d=%d,dp=%d)" % 
                     ( self.dimensionality, dp ) )
            return self._interpolateOutsideConvexHull ( massarray )
            #print "mrot=",mrot
            #print "projected upper limit=",griddata( self.Mp, self.xsec, mrot[:self.dimensionality], method="linear") * self.unit
            #max_displacement = max ( abs(mrot[self.dimensionality:]) )
            #print "maximum displacement=",max_displacement
            #mrotp=mrot
            #mrotp[0]+=max_displacement
            #print "projected+md=",griddata( self.Mp, self.xsec, mrotp[:self.dimensionality], method="linear") * self.unit
        r = griddata( self.Mp, self.xsec, mrot[:self.dimensionality], method="linear") 
        return r[0]*self.unit

    def __init__ ( self, data, accept_errors_upto=None ):
        """ Initialise SVDTrafo, giving data in the form:
            [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], ... ]
            for upper limits or 
            [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .1 ], ... ]
            for efficiency maps

            :param accept_errors_upto: If None, do not allow extrapolations outside of convex hull.
                 If float value given, allow that much relative uncertainty on the upper limit / efficiency 
                 when extrapolating outside convex hull.
                 This method can be used to loosen the equal branches assumption.
        """
        self.data = data
        self.accept_errors_upto=accept_errors_upto
        self.unit=1.0 ## store the unit so that we can take arbitrary units for the "z" values.
                   # default is unitless, which wwe use for efficiency maps
        if len(data)<1 or len(data[0])<2:
                logger.error ( "input data not in correct format. expecting sth like " \
         " [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], ... ] for upper " \
         " limits or [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], .1 ], ... ] for efficiency maps" )
        if type(data[0][1])==unum.Unum:
            ## if its a unum, we store 1.0 * unit
            self.unit=data[0][1] / ( data[0][1].asNumber() )
        self.computeV ( data )
