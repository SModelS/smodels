"""
.. module:: databaseChecks
   :synopsis: Holds the classes and methods to perform checks on the database.

.. moduleauthor:: Veronika Magerl <v.magerl@gmx.at>
.. moduleauthor:: Andre Lessa <lessa.a.p@gmail.com>

"""

import logging, os, sys
from smodels.experiment import infoObjects, dataObjects, analysisObjects

FORMAT = '%(levelname)s in %(module)s.%(funcName)s() in %(lineno)s: %(message)s'
logging.basicConfig(format=FORMAT)
logger = logging.getLogger(__name__)

logger.setLevel(level=logging.INFO)


def checkTxNameInfo(txnameInfoObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """


def checkGlobalInfo(globalInfoObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """


def checkInfoFile(infoFileObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    
    globalFields = infoFileObj.globalInfo.getattr()
    
    
    
def checkDataFile(infoFileObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """
    
    globalFields = infoFileObj.globalInfo.getattr()
    
 
def checkExpResult(expResultObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """  
    
def checkDataBase(dataBaseObj):
    """
    Checks if the info file object contains all the necessary fields
    and the field values are correctly set.
    """       