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
        # self.somearg is defined in parameters_custom_printer.ini 
        print ( f"ExamplePrinter.somearg: {self.somearg}" )
        print ( f"ExamplePrinter: but do note that we lowercase all arguments!!" )
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


parameterFile = "parameters_custom_printer.ini"
parser = modelTester.getParameters( parameterFile )


# ### Configure printers

# In[7]:


## parameters_custom_printer.ini contains a section to configure your new printer:
# [example-printer]
# someArg = "this argument was set in parameters_custom_printer.ini"


# In[8]:


## make sure to list your new printer in the list of printers used:
# [printer]
# outputType = example  ; use the example printer


# In[9]:


modelTester.loadDatabaseResults(parser, database) 


# In[10]:


fileList, inDir = modelTester.getAllInputFiles( "inputFiles/slha/simplyGluino.slha" )


# In[11]:


### Run SModelS, it will call ExamplePrinter


# In[12]:


modelTester.testPoints ( fileList, inDir, "results/", parser, database, timeout=0, development=False, parameterFile = parameterFile ) 


# ### Note how the printer can also be used via runSModelS.py directly, by adding to your ini file:

# In[13]:


## custom code to be executed                                        
# [custom-codes]                                                              
# files = ./examplePrinter.py  # comma separated


# In[ ]:




