// Test program for path searching machinery

#include "LHAPDF/Paths.h"
#include <iostream>
using namespace std;

#ifdef HAVE_MPI
#include <mpi.h>
#endif

int main() {

  #ifdef HAVE_MPI
  MPI::Init();
  #endif

  for (const string& p : LHAPDF::paths())
    cout << p << endl;

  cout << "@" << LHAPDF::findFile("lhapdf.conf") << "@" << endl;

  cout << "List of available PDFs:" << endl;
  for (const string& s : LHAPDF::availablePDFSets())
    cout << " " << s << endl;

  return 0;
}
