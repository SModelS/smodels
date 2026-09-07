#!/usr/bin/env python
# coding: utf-8

# # How To: Write your own printer to use within SModelS

# In[1]:


# Set up the path to SModelS installation folder
import sys; sys.path.append("."); import importlib; importlib.import_module("smodels_paths") if importlib.util.find_spec("smodels_paths") else None


# In[2]:


from smodels.tools.printers.basicPrinter import BasicPrinter


# ### The Printer

# In[3]:


class ExamplePrinter(BasicPrinter):
    """ here is our implementation """
    def __init__ ( self, output : str = "file", 
                   filename : str = "my.file" ):
        """ constructor, do what you will.
        In this example, we cachine in a self.toPrint
        object cache, write when .flush is called
        """
        BasicPrinter.__init__ ( self, output, filename )
        self.toPrint = []

    def setOutPutFile( self, filename : os.PathLike, overwrite : bool = True,       
        silent : bool = False ): 
        """ need to implement. Can implement your own logic here 
        :param filename: slha filename
        :param overwrite: does the user want to overwrite?
        :param silent: usually used to comment on removing old files
        """
        print ( f"ExamplePrinter, outputfile: {filename}" )

    def addObj ( self, obj ):
        """ add an object, either do something immediately with it,
        or write to an object cache 
        """
        self.toPrint.append ( obj )

    def flush ( self ):
        """ this method is called at the end of a model point """
        for obj in self.toPrint:
            print ( f"ExamplePrinter, obj: {type(obj).__name__} object: I do with it as I wish" )

# these lines register the printer with smodels, to handle the "example" extension
from smodels.tools.printers.printerRegistry import PrinterRegistry                  
PrinterRegistry.register ( ExamplePrinter, "example" );


# ### Set up modelTester

# In[4]:


from smodels.matching import modelTester
from smodels.experiment.databaseObj import Database


# ### Set up SModelS

# In[5]:


# Set the path to the database
database = Database("official")
# database.getExpResults ()


# In[6]:


parameterFile = "parameters_oor_printer.ini"
parser = modelTester.getParameters( parameterFile )


# In[7]:


modelTester.loadDatabaseResults(parser, database) 


# In[8]:


fileList, inDir = modelTester.getAllInputFiles( "inputFiles/slha/simplyGluino.slha" )


# In[9]:


### Run SModelS, it will call ExamplePrinter


# In[10]:


modelTester.testPoints ( fileList, inDir, "results/", parser, database, timeout=0, development=False, parameterFile = parameterFile ) 


# In[ ]:




