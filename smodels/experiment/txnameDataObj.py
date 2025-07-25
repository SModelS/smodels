# !/usr/bin/env python3

"""
.. module:: txnameDataObj
   :synopsis: Holds the classes to interpolate the data grid.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>
.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import sys

from smodels.base.smodelsLogging import logger
from smodels.experiment.exceptions import SModelSExperimentError as SModelSError
from smodels.tools.caching import roundCache
from scipy.linalg import svd, LinAlgError
import numpy as np
import math
from math import floor, log10


class TxNameData(object):
    """
    Holds the pre-processed data for the Txname object.
    It is responsible for computing the PCA transformation and interpolating.
    Only handles pre-processed data (1D unitless arrays, with widths rescaled).
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

    def __hash__(self):
        """ a simple unique identifier, mostly for caching """
        return id(self)

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
            msgError = f"Wrong point dimensions ({len(point)}),"
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

    @roundCache(argname="point",argpos=1,digits=2)
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
        # Round to 6 significant digits to avoid
        # numerical instabilities
        ret = self.round_to_n(float(ret),6)
        return ret

    def _estimateExtrapolationError(self, point):
        """
        When projecting a point from full_dimensionality to self.dimensionality, we
        estimate the evaluationType extrapolation error with the following
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
                    logger.error(f"negative error in result: {x:f}, {self._id}")
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

        aM = np.array(Morig, dtype=object)
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
            raise SModelSError(f"exception caught when performing singular value decomposition: {type(e)}, {e}")

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
