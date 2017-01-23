#include <SModelS.h>

int main( int argc, char * argv[] )
{
  /** initialise SModelS, load database */
  SModelS smodels ( "./parameters.ini" ); 
  /** run over one single file or directory */
  smodels.run ( "test.slha" );
}
