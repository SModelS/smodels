#include <SModelS.h>
#include <python2.7/Python.h>
#include <sstream>
#include <iostream>

using namespace std;

// https://docs.python.org/2/extending/embedding.html
// http://stackoverflow.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value

SModelS::SModelS( const string & parameterfile ) 
{
  initialize ( parameterfile, "info" );
}

SModelS::SModelS( const string & parameterfile, const string & verbose )
{
  initialize ( parameterfile, verbose );
}

SModelS::~SModelS ()
{
  Py_Finalize();
}

void SModelS::initialize ( const string & parameterfile, const string & verbose )
{
  cout << "[smodels.cpp] initialising!" << endl;
  /// initialise python, import modules
  Py_SetProgramName( (char *) (const char *) "SModelS" ); // huh??
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("sys.path.insert(0,'../smodels/tools')");
  PyRun_SimpleString("import time");
  PyRun_SimpleString("t0=time.time()");
  PyRun_SimpleString("import modelTester");
  PyRun_SimpleString("from smodelsLogging import setLogLevel");
	ostringstream set_verbosity;
	set_verbosity << "setLogLevel ( \"" << verbose << "\" ) ";
  PyRun_SimpleString( set_verbosity.str().c_str() );
	// PyObject* myModuleString = PyString_FromString((char*)"modelTester");
	// PyObject* myModule = PyImport_Import( myModuleString );
  loadDatabase ( parameterfile );
}

void SModelS::loadDatabase ( const string & parameterfile )
{
	ostringstream buffer;
	buffer << "parameterFile='" << parameterfile << "'";
  PyRun_SimpleString( buffer.str().c_str() );
  PyRun_SimpleString( "parser = modelTester.getParameters( parameterFile )" );
  PyRun_SimpleString( "database, databaseVersion = modelTester.loadDatabase(parser, None )" );
  PyRun_SimpleString( "listOfExpRes = modelTester.loadDatabaseResults(parser, database)" );
  PyRun_SimpleString( "print '[smodels.cpp] %d experimental results found.' % len(listOfExpRes) " );
}

int SModelS::run ( const string & inFile )
{
  cout << "[smodels.cpp] now running over " << inFile << endl;
	ostringstream buffer;
	buffer << "inFile='" << inFile << "'";
  PyRun_SimpleString( buffer.str().c_str() );
  PyRun_SimpleString( "fileList = modelTester.getAllInputFiles( inFile )" );
  PyRun_SimpleString( "modelTester.testPoints( fileList, inFile, 'results', parser, databaseVersion, listOfExpRes, 900, False, parameterFile )" );
  return 0;
}