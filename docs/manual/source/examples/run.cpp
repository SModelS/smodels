#include <SModelS.h>
#include <iostream>
#include <string>

/** example only */

using namespace std;

int main( int argc, char * argv[] )
{
  /** initialise SModelS, load database. The first argument 
   *  points to the SModelS config file, 
   *  the second argument (``smodelsdir'') must point to the 
   *  top-level directory of the SModelS installation */
  SModelS smodels ( "./parameters.ini", "../" ); 
  /** run over one single file or directory */
  int ret = smodels.run ( "test.slha" );
  /* return value, zero if succesful */
  return ret;
}
