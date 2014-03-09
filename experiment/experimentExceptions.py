#!/usr/bin/env python

"""
.. module:: experiment.experimentExceptions
   :synopsis: This module contains exception classes for experiment-specific use.
    
.. moduleauthor:: Wolfgang magerl <wolfgang.magerl@gmail.com>
    
"""

class MetaInfoError(Exception):
    """Exception class that is raised when a meta info field cannot be found
    """
    pass