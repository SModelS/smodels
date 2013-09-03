#!/usr/bin/env python

"""
.. module:: Experiment.experimentExceptions
   :synopsis: This module contains exception classes for Experiment-specific use.
    
.. moduleauthor:: Wolfgang magerl <wolfgang.magerl@gmail.com>
    
"""

class MetaInfoError(Exception):
    """Exception class that is raised when a meta info field cannot be found
    """
    pass