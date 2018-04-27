#!/usr/bin/env python

"""
.. module:: flongCalc
    :synopsis: Tool to compute the fraction of long-lived and prompt decays
                for a given particle.

.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

from smodels.tools.physicsUnits import mm, m, MeV,fm
from smodels.theory.exceptions import SModelSTheoryError as SModelSError
import math


def FlongCalculator(pdg,mass,width,
                    l_inner=10.*mm,gb_inner=10,l_outer=10.*m,gb_outer=0.6,**kargs):
    
    """
    Using the given width, pdg and mass computes the fraction of
    prompt (F_prompt) and long-lived decays (F_long) 
    F_prompt = 1 - exp(-width*l_inner/gb_inner)
    F_long = exp(-width*l_outer/gb_outer)
    where l_inner is the inner radius of the detector, l_outer is the outer radius
    and gb_x is the estimate for the kinematical factor gamma*beta for each case.
    We use gb_outer = 10 and gb_inner= 0.5
    
    :param pdg: Particle pdg code
    :param mass: Particle mass with units (e.g. 100*GeV)
    :param width: Particle Width with units (e.g. 1e-5*GeV)
    :param l_inner: Radius of the inner tracker
    :param gb_inner: Effective gamma*beta factor to be used for prompt decays
    :param l_outer: Radius of the outer detector
    :param gb_outer: Effective gamma*beta factor to be used for long-lived decays
    
    :return: Dictionary = {Flong : Flong, Fprompt : Fprompt}
    """
    
    
    
    hc = 197.327*MeV*fm  #hbar * c

    if mass.asNumber() == 0. or width.asNumber() == 0.:
        Flong = 1.
        Fprompt = 0.
    else:
        Fprompt = 1. - math.exp(-width*l_inner/(gb_inner*hc))
        Flong = math.exp(-width*l_outer/(gb_outer*hc))            
        if Flong < 1e-50:
            Flong = 0.
             
    if (Flong+Fprompt) > 1.:
        raise SModelSError("Sum of decay fractions > 1 for "+str(pdg))

    return {'Flong': Flong, 'Fprompt' : Fprompt}
