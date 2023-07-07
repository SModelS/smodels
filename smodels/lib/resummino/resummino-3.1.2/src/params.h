// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Parameter module.

#ifndef PARAMS_H_
#define PARAMS_H_

#include <complex>
#include <map>

using namespace std;

#define nh 6 // number of Higgs
#define nq 6 // number of quarks
#define nv 5 // number of vector bosons
#define nl 6 // number of neutrinos + leptons
#define nSQ 12 // number of squarks
#define nCH 6 // number of charginos
#define nSL 9 // number of sleptons

using namespace std;

// Coupling structure.
struct Coupling {
  complex<double> L; // Left-handed component
  complex<double> R; // Right-handed component
};

class Parameters {
public:
  // Parameters() {};
  Parameters(){};

  ~Parameters(){};

  // Sets process parameters.
  void set_process_params(const int Icoll, const int O1, const int O2, const double En,
                          const double Mi, const double Pt, const double Ff, const double Fr,
                          const int pdfset);

  // Incoming and outgoing particles.
  int in1;  // incoming particle 1
  int in2;  // incoming particle 2
  int out1; // outgoing particle 1
  int out2; // outgoing particle 2

  // Collider properties.
  int ic;     // collider type
  double sh;  // collider cms energy
  double pts; // transverse momentum (pts < 0 for no restriction)
  double mis; // invariant mass (mis < 0 for default)

  // Renormalisation and Factorisation scales.
  double murs, mufs;

  // Pdf parameters.
  double a1min;
  double afit[8][8];
  int set;

  // Minimum invariant mass.
  double Minv_min;
  // Maximum invariant mass.
  double Minv_max;

  // CKM matrix.
  complex<double> ckm[3][3];

  // SM Masses.
  double mh[nh]; // Higgs masses
  double mq[nq]; // Quark masses
  double mv[nv]; // Gauge boson masses
  double Gv[nv]; // Gauge boson widths
  double ml[nl]; // Lepton masses

  double mhs[nh]; // Higgs masses
  double mqs[nq]; // Quark masses
  double mvs[nv]; // Gauge boson masses
  double Gvs[nv]; // Gauge boson widths
  double mls[nl]; // Lepton masses

  // SUSY Masses.
  double mGL;      // Gluino mass
  double mSQ[nSQ]; // Squark masses
  double Gsq[nSQ]; // Squark widths
  double mCH[nCH]; // Neutralino/chargino masses
  double mSL[nSL]; // Charged Sleptons and sneutrino masses

  // SUSY Masses.
  double mGLs;      // Gluino mass
  double mSQs[nSQ]; // Squark masses
  double mCHs[nCH]; // Neutralino/chargino masses
  double mSLs[nSL]; // Charged Sleptons and sneutrino masses

  // susy breanching ratios
  double SQ_CH_q[nSQ][nCH][nq];
  double SQ_GL_q[nSQ][nq];

  // degenerate squarks?
  bool deg_squarks;

  // Weak Couplings.
  // h: Higgs, q: quark, SQ: squark, CH: gaugino, l: lepton, SL: slepton.
  struct Coupling hqq[nh][nq][nq];
  struct Coupling hSQSQ[nh][nSQ][nSQ];
  struct Coupling hCHCH[nh][nCH][nCH];
  struct Coupling vqq[nv][nq][nq];
  struct Coupling vSQSQ[nv][nSQ][nSQ];
  struct Coupling vCHCH[nv][nCH][nCH];
  struct Coupling vll[nv][nl][nl];
  struct Coupling CHSQq[nCH][nSQ][nq];
  struct Coupling CHqSQ[nCH][nq][nSQ];
  struct Coupling vSLSL[nv][nSL][nSL];

  // Strong Couplings.
  // g: gluon, q: quark, GL: gluino, SQ: squark.
  struct Coupling gqq[nq][nq];
  struct Coupling gSQSQ[nSQ][nSQ];
  struct Coupling GLSQq[nSQ][nq];
  struct Coupling GLqSQ[nq][nSQ];
  struct Coupling SQSQSQ[nSQ][nSQ][nSQ];

  // Coupling Counterterms.
  // CH: gaugino, SQ: squark, q: quark
  struct Coupling dCHSQq[nCH][nSQ][nq];
  struct Coupling dCHqSQ[nCH][nq][nSQ];

  // SUSY Mixings.
  complex<double> RA[2][2];   // Neutral Higgs mixing matrix
  complex<double> RB[2][2];   // Charged Higgs mixing matrix
  complex<double> RSD[6][6];  // Down-type squark mixing matrix
  complex<double> RSU[6][6];  // Up-type squark mixing matrix
  complex<double> RN[4][4];   // Neutralino mixing matrix
  complex<double> RU[2][2];   // Chargino- mixing matrix
  complex<double> RV[2][2];   // Chargino+ mixing matrix
  complex<double> SSLM[2][2]; // Stau LR mixing

  complex<double> XU[6][6];        // p.19 hep/9511250
  complex<double> XD[6][6];        // ^ U -> D
  complex<double> XUU[6][6][6][6]; // only the g_3 part 3X_UX_U-X_UX_U from p.47 hep/9511250
  complex<double> XDD[6][6][6][6]; // ^ with D

  // DR ParameterReader (SM-like part).
  double G_F, a_em, a_s; // SM parameters
  double g1, g2, g3;     // Coupling constants
  double sw, cw, xw;     // Weak angle
  double yq[nq];         // Yukawa couplings

  // Integration parameters.

  double precision;             // desired relative precision
  double abs_precision = 1e-14; // desired absolute precision
  int max_iters;                // maximum number of iterations
  FILE *fout;                   // where to store integration information
  int integration_flag;         // Vegas default
  double ptcut = 0.;
  double width = 1e-3;

  void init_couplings();
  void set_couplings();
  void set_params(const char *file1, const char *file2);
  void set_squark_width(const char *file);
  void read_slha(const char *file);
  void read_ZpWp(const char *file);
  void write_log(const char *file);

  map<string, string> read_input_file(string filename);
};

#endif
