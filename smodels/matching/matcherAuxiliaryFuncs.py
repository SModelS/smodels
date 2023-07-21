"""
.. module:: matcherAuxiliaryFunctions
   :synopsis: A collection of functions used to evaluate fuzzy the conditions.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""


import unum
import re
import numpy as np
try:
    from collections.abc import Iterable
except (ImportError, ModuleNotFoundError):
    from collections import Iterable

from smodels.base.physicsUnits import standardUnits, GeV
from smodels.matching.exceptions import SModelSMatcherError as SModelSError
from smodels.base.smodelsLogging import logger



def roundValue(value, nround=-1):
    """
    Round a value to nround significant digits. If the input value is not
    a float or Unum object, it is zero or infinity, nothing is done.

    :param nround: number of significant digits

    :return: rounded value
    """

    if nround <= 0:
        return value
    elif not isinstance(value, (float, unum.Unum)):
        return value

    # Remove units, if it is a unum object
    if isinstance(value, unum.Unum):
        if not value.asNumber():
            return value
        elif np.isinf(value.asNumber()):
            return value
        unit = unum.Unum(value._unit)
        v = value.asNumber(unit)
    else:
        unit = 1.0
        v = value

    # Round value:
    if v != 0.0:
        v_rounded = round(v, nround-int(np.floor(np.log10(abs(v))))-1)
    else:
        v_rounded = v
    v_rounded *= unit
    return v_rounded


def average(values, weights=None, nround=-1):
    """
    Compute the weighted average of a list of objects.
    All the objects must be of the same type.
    If all objects are equal returns the first entry of the list.
    Only the average of ints, floats and Unum objects or nested lists
    of these can be computed. If the average can not be computed
    returns None.

    :param values: List of objects of the same type
    :param weights: Weights for computing the weighted average. If None it will assume
                    unit weights.
    :param nround: If greater than zero and the returning attibute is numeric, will round it
                      to this number of significant digits.

    """

    if weights is None:
        weights = [1.]*len(values)
    if not isinstance(values, list):
        raise SModelSError("Values must be a list of objects")
    if not isinstance(weights, list):
        raise SModelSError("Weights must be a list of objects")
    if any(not isinstance(w, (float, int)) for w in weights):
        raise SModelSError("Weights must be a list of integers or floats")
    if len(values) != len(weights):
        return SModelSError("Values and weights must have the same length")

    if any(type(v) != type(values[0]) for v in values):
        raise SModelSError("Can not compute average of distinct types of objects")
    if isinstance(values[0], (float, int, unum.Unum)):
        total = values[0]*weights[0]
        for i, v in enumerate(values[1:]):
            total += v*weights[i+1]
        total = total/sum(weights)
        return roundValue(total, nround)
    elif isinstance(values[0], list):
        ndim = len(values[0])
        if any(len(v) != ndim for v in values):
            raise SModelSError("Can not compute average of lists of distinct lengths")
        res = []
        for i in range(ndim):
            try:
                res.append(average([v[i] for v in values], weights, nround))
            except SModelSError:
                return None
        return res
    else:
        if all(values[0] is v for v in values):
            return roundValue(values[0], nround)
        if all(values[0] == v for v in values):
            return roundValue(values[0], nround)
        return None
