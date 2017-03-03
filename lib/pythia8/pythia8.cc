// This is our program to read slha files and compute cross sections

#include "Pythia.h"
// #include "Pythia8/Pythia.h"
#include <sstream>
#include <cstdlib>
using namespace Pythia8;
using namespace std;

void checkEntry ( ParticleDataEntry & entry )
{
    cout << endl;
    int id = entry.id();
    // widths greater than 100 GeV for BSM particles. suspicious!
    if ( (abs(id) > 999999) && ( entry.mWidth() > 100. ) )
    {
      cout << "[ParticleDataEntry] entry " << id << " has very large particle width: "
           << entry.mWidth() << "! (Is there no decay table in the slha file?)" 
           << endl;
      entry.setMWidth(1.0);
      entry.clearChannels();
      entry.addChannel ( 1, 1.0, 100, 1, -1 );
    }
    /*
    cout << "[ParticleDataEntry] id=" << entry.id() << " w=" 
         << entry.mWidth() << " mMin=" << entry.mMin() 
         << " mMax=" << entry.mMax() << " fw=" << entry.doForceWidth()
         << " m0=" << entry.m0() << " mayDecay=" << entry.mayDecay() 
         << " uBW=" << entry.useBreitWigner() << " channels=" << entry.sizeChannels()
         << endl;
    for ( int i=0; i< entry.sizeChannels() ; i++ )
    {
      DecayChannel c = entry.channel ( i );
      cout << "[DecayChannel] onM=" << c.onMode() << " m="
           << c.multiplicity() << " bR=" << c.bRatio() << " changed=" 
           << c.hasChanged() << " width=" << c.onShellWidth();
      for ( int p=0; p< c.multiplicity() ; p++ )
      {
        cout << " pdg=" << c.product(p) << " ";
      };
      cout << endl;
    };
    */
}

void check ( ParticleData & data )
{
  int ctr=0;
  int last = 1;
  while ( last != 0 )
  {
    ctr++;
    int now = data.nextId( last );
    if ( now == 0 )
    {
      // cout << "[pythia8.cc] we hit the end of the particle table" << endl;
      return;
    } else {
      ParticleDataEntry * entry = data.particleDataEntryPtr ( now );
      checkEntry ( *entry );
    }
    if ( ctr > 10000 )
    {
      cout << "[pythia8.cc] we have more than 10,000 particles in our database."
              " This cannot be. I stop." << endl;
      exit(-1);
    }
    last = now;
  }
}

int run( int nevents, float sqrts /* in TeV */, const string & slhafile,
         const string & cfgfile, const string & xmldir )
{
  cout << "[pythia8.exe] we run with " << nevents << ", sqrts=" << sqrts << " TeV"
       << ", slhafile=" << slhafile << endl;
  Pythia pythia ( xmldir, false );
  pythia.readFile( cfgfile );
  ostringstream o1, o2, o3;
  o1 << "Main:numberOfEvents=" << nevents;
  pythia.readString( o1.str() );
  o2 << "Beams:eCM=" << sqrts * 1000;
  pythia.readString( o2.str() );
  o3 << "SLHA:file=" << slhafile;
  pythia.readString ( o3.str() );
  pythia.readString ( "SLHA:allowUserOverride = true" );
  pythia.init();
  check ( pythia.particleData );

  // Begin event loop. Generate event. Skip if error.
  for (int iEvent = 0; iEvent < nevents; ++iEvent ) {
    if (!pythia.next()) continue;
  }
  pythia.stat();
  pythia.particleData.listAll();
  return 0;
}

void help( const char * name )
{
  cout << "syntax: " << name << " [-h] [-n <nevents>] [-s <sqrts>] [-f <slhafile>] " 
       << "[-c <cfgfile>]" << endl;
  cout << "        -n <nevents>:  number of events to simulate [10000]." << endl;
  cout << "        -s <sqrts>:    center-of mass energy, in TeV [13]." << endl;
  cout << "        -f <slhafile>: slhafile [test.slha]." << endl;
  cout << "        -c <cfgfile>: the pythia8 config file [./pythia8.cfg]." << endl;
  cout << "        -x <xmldir>: the pythia8 xmldoc dir [../xmldoc]." << endl;
  exit( 0 );
};

int main( int argc, const char * argv[] ) {
  int nevents = 10000;
  float sqrts = 13;
  string slhafile = "test.slha";
  string cfgfile = "./pythia8.cfg";
  string xmlDir = "../xmldoc";
  for ( int i=1; i!=argc ; ++i )
  {
    string s = argv[i];
    if ( s== "-h" )
    {
      help ( argv[0] );
    }

    if ( s== "-n" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      nevents = atoi ( argv[i+1] );
      i++;
      continue;
    }

    if ( s== "-f" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      slhafile = argv[i+1];
      i++;
      continue;
    }

    if ( s== "-c" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      cfgfile = argv[i+1];
      i++;
      continue;
    }

    if ( s== "-x" )
    {
      if ( argc < i+2 ) help ( argv[0] );
      xmlDir = argv[i+1];
      i++;
      continue;
    }

    if ( (s== "-s") )
    {
      if ( argc < i+2 ) help ( argv[0] );
      sqrts = atof ( argv[i+1] );
      if ( sqrts > 360. )
      {
        cout << "Warning. You supplied an extremely high number for sqrts: "
             << sqrts << " TeV. Could it be that you gave the value in GeV?" << endl;
      }
      i++;
      continue;
    }

    cout << "Error. Argument " << argv[i] << " unknown." << endl;
    help ( argv[0] );
  };

  return run ( nevents, sqrts, slhafile, cfgfile, xmlDir );
}
