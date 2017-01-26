// This is our program to read slha files and compute cross sections

#include "Pythia8/Pythia.h"
using namespace Pythia8;
int main() {
  // Generator. Process selection. LHC initialization. Histogram.
  Pythia pythia;
  pythia.readFile("pythia8.cfg");
  int nEvent   = pythia.mode("Main:numberOfEvents");
  int nAbort   = pythia.mode("Main:timesAllowErrors");
  double eCM   = pythia.parm("Beams:eCM");
  // LHAupFromPYTHIA8 myLHA(&pythia.process, &pythia.info);
  // myLHA.openLHEF("pythia8.lhe");
  pythia.init();
  // myLHA.initLHEF();

  // Begin event loop. Generate event. Skip if error. List first one.
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) continue;
    //myLHA.setEvent(); myLHA.eventLHEF();
  }
  pythia.stat();
  /* cout << "end pythia.stat" << endl;
  myLHA.updateSigma();
  myLHA.closeLHEF(true); */
  return 0;
}
