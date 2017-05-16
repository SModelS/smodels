#include <SModelS.h>
#include <iostream>
#include <string>

/** example only */

using namespace std;

int main( int argc, char * argv[] )
{
  /** initialise SModelS, load database, second argument (installdir) must point to smodels 
   *  installation top-level directory */
  SModelS smodels ( "./parameters.ini", "../" ); 
  /** run over one single file or directory */
  int ret = smodels.run ( "test.slha" );
}
