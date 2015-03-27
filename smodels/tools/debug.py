"""
.. module:: tools.debug
   :synopsis: Facility used in runSModelS.py to create and read SModelS debug files.

.. moduleauthor:: Wolfgang Magerl <wolfgang.magerl@gmail.com>

"""

import os
import datetime
import platform
import traceback


class Debug(object):
    """
    Class that handles all debug information.
    
    """    
    def __init__(self):        
        timestamp = datetime.datetime.now().strftime('%Y%m%d%H%M%S%f')
        self.timestampHuman = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S.%f')
        self.debugFileName = 'smodels-' + timestamp + '.debug'


    def createDebugFile(self, inputFileName, parameterFileName):
        """
        Create a new SModelS debug file.
        
        A SModelS debug files contains:
        
        - a timestamp
        - SModelS version
        - platform information (CPU architecture, operating system, ...)
        - Python version
        - stack trace
        - input file name
        - input file content
        - parameter file name
        - parameter file content
        
        :param inputFileName: relative location of the input file
        :param parameterFileName: relative location of the parameter file
        
        """
    
        with open('smodels/version', 'r') as versionFile:
            version = versionFile.readline()
    
        with open(inputFileName, 'r') as inputFile:
            inputFileContent = inputFile.read()
    
        with open(parameterFileName, 'r') as parameterFile:
            parameterFileContent = parameterFile.read()
    
        debugFile = open(self.debugFileName, 'w')
        debugFile.write("================================================================================\n")
        debugFile.write("SModelS Debug File\n")
        debugFile.write("================================================================================\n")
        debugFile.write("Timestamp: " + self.timestampHuman + "\n\n")
        debugFile.write("SModelS Version: " + version + "\n")
        debugFile.write("Platform: " + platform.platform() + "\n")
        debugFile.write("Python Version: " + platform.python_version() + "\n\n")
        debugFile.write("================================================================================\n\n")
        debugFile.write("--------------------------------------------------------------------------------\n")
        debugFile.write("* Output\n")
        debugFile.write("--------------------------------------------------------------------------------\n\n")
        debugFile.write(traceback.format_exc() + "\n\n")
        debugFile.write("--------------------------------------------------------------------------------\n")
        debugFile.write("* Input File\n")
        debugFile.write("  " + os.path.basename(inputFileName) + "\n")
        debugFile.write("--------------------------------------------------------------------------------\n\n")
        debugFile.write(inputFileContent + "\n")
        debugFile.write("--------------------------------------------------------------------------------\n")
        debugFile.write("* Parameter File\n")
        debugFile.write("  " + os.path.basename(parameterFileName) + "\n")
        debugFile.write("--------------------------------------------------------------------------------\n\n")
        debugFile.write(parameterFileContent)
        debugFile.close()
        
    
    def createUnknownErrorMessage(self):
        """
        Create a message for an unknown error.
        
        """
        message = ("\n\n\n"
                   "================================================================================\n\n"
                   "SModelS quit unexpectedly due to an unknown error. The error has been written to\n"
                   "" + self.debugFileName + ".\n\n"
                   "Please send this file to smodels-users@lists.oeaw.ac.at and shortly describe\n"
                   "what you did to help making SModelS better!\n\n"
                   "Alternatively, use the '--development' option when running runSModelS.py to\n"
                   "display all error messages.\n\n"
                   "================================================================================")
        return message
    
    
def readDebugFile(debugFileName):
    """
    Read a debug file to use its input and parameter file sections for a SModelS run.
    
    :param debugFileName: relative location of the debug file
    
    """
    with open(debugFileName, 'r') as debugFile:
        debugFileContent = debugFile.readlines()
        
    lineNumber = 0
    inputStartLine = 0
    inputEndLine = 0
    parameterStartLine = 0
    parameterEndLine = 0
    
    for line in debugFileContent:
        if lineNumber == 1:
            if not line.rstrip() == "SModelS Debug File":
                print("ERROR: Not a SModelS debug file!")
                break
            
        if line.rstrip() == "* Input File":
            inputStartLine = lineNumber + 4
            
        if line.rstrip() == "* Parameter File":
            inputEndLine = lineNumber - 2
            parameterStartLine = lineNumber + 4
        
        lineNumber += 1
        
    parameterEndLine = lineNumber
    
    debugInputFileName = 'debug_input'
    debugParameterFileName = 'debug_parameter'
    
    debugInputFile = open(debugInputFileName, 'w')
    debugParameterFile = open(debugParameterFileName, 'w')
    
    for i in range(inputStartLine, inputEndLine):
        debugInputFile.write(debugFileContent[i])
    
    for i in range(parameterStartLine, parameterEndLine):
        debugParameterFile.write(debugFileContent[i])
        
    debugInputFile.close()
    debugParameterFile.close()
    
    return debugInputFileName, debugParameterFileName
    
    
def createStackTrace():
    """
    Return the stack trace.
    
    """
    return traceback.format_exc()

