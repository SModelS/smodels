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
import logging

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.WARNING )

class SVDTrafo:
    def countNonZeros ( self, mp ):
        """ count the nonzeros in a vector """
        nz=0
        for i in mp:
            if abs(i)<10**-5:
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
            self.xsec.append ( y / fb )
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
        #print "Mp,xseec=",zip(self.Mp,self.xsec)

    def getUpperLimit ( self, massarray ):
        """
           :param massarray: e.g. [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ]
           :return: upper limit, in fb 
        """
        m=self.flattenMassArray ( massarray )
        mrot=np.dot(m,self.V)
        if self.countNonZeros ( mrot ) != self.dimensionality:
            logger.warning ( "trying to interpolate outside of convex hull" )
            print "mrot=",mrot
            return float('nan')*fb
        r = griddata( self.Mp, self.xsec, mrot[:self.dimensionality], method="linear") 
        return r*fb

    def __init__ ( self, data ):
        """ Initialise SVDTrafo, giving data in the form:
            [ [ [[ 300.*GeV,100.*GeV], [ 300.*GeV,100.*GeV] ], 10.*fb ], ... ]
        """
        self.data = data
        self.computeV ( data )
