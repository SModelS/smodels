#include <SModelS.h>
#include <Python.h>
#include <sstream>
#include <iostream>

using namespace std;

// https://docs.python.org/2/extending/embedding.html
// http://stackoverflow.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value

SModelS::SModelS( const string & parameterfile, const string & smodelsdir ) 
{
  initialize ( parameterfile, smodelsdir, "info" );
}

SModelS::SModelS( const string & parameterfile, const string & smodelsdir, 
                  const string & verbose )
{
  initialize ( parameterfile, smodelsdir, verbose );
}

SModelS::~SModelS ()
{
  Py_Finalize();
}

void SModelS::initialize ( const string & parameterfile, const string & smodelsdir, 
                           const string & verbose )
{
  cout << "[smodels.cpp] initialising!" << endl;
  /// initialise python, import modules
  // Py_SetProgramName( (char *) (const char *) "SModelS" ); // huh?! python 2.x
  Py_SetProgramName( (wchar_t *) (const wchar_t *) "SModelS" ); // huh?? python 3.x
  Py_Initialize();
  PyRun_SimpleString("import sys");
  PyRun_SimpleString("import os");
	ostringstream check_path;
  check_path << "if not os.path.exists('" << smodelsdir << "/smodels/version'): print "
             << "('WARNING: installation directory " << smodelsdir << " does not point"
             << " to a SModelS installation!')";
  PyRun_SimpleString( check_path.str().c_str() );
	ostringstream set_path;
	set_path << "sys.path.insert(0,'" << smodelsdir << "') ";
  PyRun_SimpleString( set_path.str().c_str() );
  PyRun_SimpleString("import time");
  PyRun_SimpleString("t0=time.time()");
  PyRun_SimpleString("from smodels.tools import modelTester");
  PyRun_SimpleString("from smodels.tools.smodelsLogging import setLogLevel");
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
  PyRun_SimpleString( "print ( '[smodels.cpp] %d experimental results found.' % len(listOfExpRes) ) " );
}

int SModelS::run ( const string & inFile )
{
  cout << "[smodels.cpp] now running over " << inFile << endl;
	ostringstream buffer;
	buffer << "inFile='" << inFile << "'";
  PyRun_SimpleString( buffer.str().c_str() );
  PyRun_SimpleString( "fileList, inDir = modelTester.getAllInputFiles( inFile )" );
  PyRun_SimpleString( "modelTester.testPoints( fileList, inDir, 'results', parser, databaseVersion, listOfExpRes, 900, False, parameterFile )" );
  //PyRun_SimpleString( "print('pyx=',pyx)" );
  // double result = PyFloat_AsDouble(myResult);
  return 0;
}
