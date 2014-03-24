#!/usr/bin/env python

"""
.. module:: physicsUnits
   :synopsis: This introduces physical units (GeV, fb) to the framework.

.. moduleauthor:: Wolfgang Waltenberger <wolfgang.waltenberger@gmail.com>

"""

import unum
import logging

logger = logging.getLogger(__name__)


# description see
# http://home.scarlet.be/be052320/Unum.html
# http://home.scarlet.be/be052320/docs.html

# install it with:
# sudo pip install ez_setup
# sudo easy_install Unum
# or sth similar

# can be switched off with a single switch
useUnits = True

# make sure we define only once
# if len(Unum.getUnitTable())==0:
unum.Unum.reset()
unum.VALUE_FORMAT = "%0.2E"
unum.UNIT_HIDE_EMPTY = True

fb = unum.Unum.unit('fb')
pb = unum.Unum.unit('pb', 1000 * fb)

eV = unum.Unum.unit('eV')
keV = unum.Unum.unit('keV', 10**3*eV)
MeV = unum.Unum.unit('MeV', 10**6*eV)
GeV = unum.Unum.unit('GeV', 10**9*eV)
TeV = unum.Unum.unit('TeV', 10**12*eV)


def addunit(value, unitstring):
    """
    Adds units to values, but also makes it easy to turn this functionality
    off, in case "units" isn't installed.
    
    """
    if value == None:
        return value
    if not useUnits:
        return value
    if useUnits:
        # for convenience, we add units also to tuples, lists, and dictionaries
        if isinstance(value, list):
            return [addunit(x, unitstring) for x in value]
        if isinstance(value, tuple):
            return tuple([addunit(x, unitstring) for x in value])
        if isinstance(value, dict):
            ret = {}
            for (k, v) in value.items():
                ret[k] = addunit(v, unitstring)
            return ret
        if not isinstance(value, float) and not isinstance(value, int):
            return value
        if unitstring == "GeV":
            return value * GeV
        if unitstring == "TeV":
            return value * TeV
        if unitstring == "fb":
            return value * fb
        if unitstring == "pb":
            return value * pb
        if unitstring == "fb-1":
            return value / fb
        logger.warning("Unknown unit: " + unitstring)
    return value


def rmvunit ( value, unitstring ):
    """
    Removes units from values, but also makes it easy to turn this
    functionality off, in case "units" isn't installed.
    
    """
    if not useUnits: return value
    if useUnits: 
        if type(value) != type(1.*GeV):
            return value
        if unitstring == "GeV":
            return value.asNumber(GeV)
        if unitstring == "TeV":
            return value.asNumber(TeV)
        if unitstring == "fb":
            return value.asNumber(fb)
        if unitstring == "pb":
            return value.asNumber(pb)
        if unitstring == "fb-1":
            return value.asNumber(1/fb)
        logger.warning("Unknown unit: " + unitstring)
        return value


if __name__ == "__main__":
    """
    Called as script, will print some physicsUnits.
    
    """
    three = addunit(3.0, "fb")
    print(three, "=", three.asUnit(pb))
    seven = addunit(7., "TeV")
    print(seven, "=", seven.asUnit(GeV))
