// This is our program to read slha files and compute cross sections

#include "Pythia8/Pythia.h"
#include <sstream>
using namespace Pythia8;

int run( int nevents, int sqrts /* in TeV */ ) {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readFile("pythia8.cfg");
  ostringstream o1, o2;
  o1 << "Main:numberOfEvents=" << nevents;
  pythia.readString( o1.str() );
  o2 << "Beams:eCM=" << sqrts * 1000;
  pythia.readString( o2.str() );
  pythia.init();

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nevents; ++iEvent ) {
    if (!pythia.next()) continue;
  }
  pythia.stat();
  return 0;
}

int main( int argc, const char * argv[] ) {
  return run ( 100, 7 );
}
