// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Reads and sets the parameters.

#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

#include "slhaea.h"

#include "kinematics.h"
#include "params.h"
#include "utils.h"

#include "Fastor/Fastor.h"
#undef RST
#include "tensors.h"

using namespace std;
using namespace Fastor;

// SUSY couplings.
#define SUSY

// Couplings for W' and Z'.
// Careful here, with SSM defined we basically don't read the input file.
//#define SSM
//#define PYTHIA

// Turn off squark mixing.
//#define NO_SQUARK_MIXING

// Use Gauge coupling conventions.
//#define GAUGE_COUPLING
// Use Prospino coupling conventions.
#define PROSPINO_COUPLING
// Use MadGraph coupling conventions. Also recomputes yu3x3, rest zero like LIGHT_YUKAWAS_ZERO.
//#define MADGRAPH_COUPLING
//#define MADGRAPH_COUPLING_OLD
// Set Light Yukawas to zero
#define LIGHT_YUKAWAS_ZERO

// IF Gauge or Pro coupling fails ( missing params in param_card.dat) resort to MG

// Use diagonal CKM matrix.
#define CKM_DIAG
// Load CKM from slha. CKM_DIAG overwrites this.
//#define CKM_LOAD

// Use degenerate squark masses.
// Prospino conventions.
// without stops and sbottom.
//#define DEGENERATE_SQUARK_MASS

// complex i.
#define II complex<double>(0.0, 1.0)

map<string, string> Parameters::read_input_file(string filename) {
  map<string, string> config;
  errno = 0;
  FILE *input_file;

  if (filename == "-") {
    input_file = stdin;
  } else {
    input_file = fopen(filename.c_str(), "r");
  }
  if (input_file == NULL) {
    fprintf(stderr, "error: could not open file '%s'.\n", filename.c_str());
    char buff[FILENAME_MAX];
    getcwd(buff, FILENAME_MAX);
    printf("Current working dir: %s\n", buff);
    printf("Error %d \n", errno);
    exit(1);
  }

  char *variable = (char *)malloc(128 + 1);
  char *value = (char *)malloc(128 + 1);
  ;
  for (;;) {
    char *input;
    size_t n;
    int chars_read;
    int elements_read;

    input = NULL;
    n = 0;
    chars_read = getline(&input, &n, input_file);
    if (chars_read == -1) {
      break;
    }
    if (chars_read < 2) {
      continue;
    }
    elements_read = sscanf(input, "%128[a-zA-Z0-9_-] = %128s\n", variable, value);
    if (elements_read == 0) {
      // Blank line.
      continue;
    } else if (elements_read == 1 && variable[0] == '#') {
      // Comment.
      continue;
    } else if (elements_read == 1) {
      fprintf(stderr,
              "error: while reading line '%s' of the input file at line %li. (param_card.dat as "
              "input.in?)\n",
              variable, n);
      exit(1);
    } else if (variable[0] == '#') {
      // Comment.
      continue;
    }
    config[variable] = value;
  }
  free(variable);
  free(value);
  fclose(input_file);
  return config;
}

void Parameters::init_couplings() {
  // Higgs couplings
  for (int i0 = 0; i0 < nh; i0++) {
    for (int i1 = 0; i1 < nq; i1++) {
      for (int i2 = 0; i2 < nq; i2++) {
        hqq[i0][i1][i2].L = complex<double>(0.0, 0.0);
        hqq[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  for (int i0 = 0; i0 < nh; i0++) {
    for (int i1 = 0; i1 < nSQ; i1++) {
      for (int i2 = 0; i2 < nSQ; i2++) {
        hSQSQ[i0][i1][i2].L = complex<double>(0.0, 0.0);
        hSQSQ[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  for (int i0 = 0; i0 < nh; i0++) {
    for (int i1 = 0; i1 < nCH; i1++) {
      for (int i2 = 0; i2 < nCH; i2++) {
        hCHCH[i0][i1][i2].L = complex<double>(0.0, 0.0);
        hCHCH[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  // Vector couplings to quarks.
  for (int i0 = 0; i0 < nv; i0++) {
    for (int i1 = 0; i1 < nq; i1++) {
      for (int i2 = 0; i2 < nq; i2++) {
        vqq[i0][i1][i2].L = complex<double>(0.0, 0.0);
        vqq[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  // Vector couplings to leptons.
  for (int i0 = 0; i0 < nv; i0++) {
    for (int i1 = 0; i1 < nl; i1++) {
      for (int i2 = 0; i2 < nl; i2++) {
        vll[i0][i1][i2].L = complex<double>(0.0, 0.0);
        vll[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  // Vector couplings to squarks.
  for (int i0 = 0; i0 < nv; i0++) {
    for (int i1 = 0; i1 < nSQ; i1++) {
      for (int i2 = 0; i2 < nSQ; i2++) {
        vSQSQ[i0][i1][i2].L = complex<double>(0.0, 0.0);
        vSQSQ[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  for (int i0 = 0; i0 < nv; i0++) {
    for (int i1 = 0; i1 < nCH; i1++) {
      for (int i2 = 0; i2 < nCH; i2++) {
        vCHCH[i0][i1][i2].L = complex<double>(0.0, 0.0);
        vCHCH[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
  // Gaugino couplings
  for (int i0 = 0; i0 < nCH; i0++) {
    for (int i1 = 0; i1 < nSQ; i1++) {
      for (int i2 = 0; i2 < nq; i2++) {
        CHSQq[i0][i1][i2].L = complex<double>(0.0, 0.0);
        CHSQq[i0][i1][i2].R = complex<double>(0.0, 0.0);
        CHqSQ[i0][i2][i1].L = complex<double>(0.0, 0.0);
        CHqSQ[i0][i2][i1].R = complex<double>(0.0, 0.0);

        dCHSQq[i0][i1][i2].L = complex<double>(0.0, 0.0);
        dCHSQq[i0][i1][i2].R = complex<double>(0.0, 0.0);
        dCHqSQ[i0][i2][i1].L = complex<double>(0.0, 0.0);
        dCHqSQ[i0][i2][i1].R = complex<double>(0.0, 0.0);
      }
    }
  }
  // Strong couplings
  for (int i0 = 0; i0 < nq; i0++) {
    for (int i1 = 0; i1 < nq; i1++) {
      gqq[i0][i1].L = complex<double>(0.0, 0.0);
      gqq[i0][i1].R = complex<double>(0.0, 0.0);
    }
  }
  for (int i0 = 0; i0 < nSQ; i0++) {
    for (int i1 = 0; i1 < nSQ; i1++) {
      gSQSQ[i0][i1].L = complex<double>(0.0, 0.0);
      gSQSQ[i0][i1].R = complex<double>(0.0, 0.0);
    }
  }
  for (int i0 = 0; i0 < nSQ; i0++) {
    for (int i1 = 0; i1 < nq; i1++) {
      GLSQq[i0][i1].L = complex<double>(0.0, 0.0);
      GLSQq[i0][i1].R = complex<double>(0.0, 0.0);
      GLqSQ[i1][i0].L = complex<double>(0.0, 0.0);
      GLqSQ[i1][i0].R = complex<double>(0.0, 0.0);
    }
  }
  for (int i0 = 0; i0 < nSQ; i0++) {
    for (int i1 = 0; i1 < nSQ; i1++) {
      for (int i2 = 0; i2 < nSQ; i2++) {
        SQSQSQ[i0][i1][i2].L = complex<double>(0.0, 0.0);
        SQSQSQ[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }

  // Vector coupling to sleptons.
  for (int i0 = 0; i0 < nv; i0++) {
    for (int i1 = 0; i1 < nSL; i1++) {
      for (int i2 = 0; i2 < nSL; i2++) {
        vSLSL[i0][i1][i2].L = complex<double>(0.0, 0.0);
        vSLSL[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }
}

void Parameters::set_couplings() {

#ifdef PYTHIA
  xw = 0.231;
  sw = sqrt(xw);
  cw = sqrt(1.0 - xw);
  g2 = sqrt(4.0 * M_PI * 8.1675767e-03) / sw;
  g3 = sqrt(4.0 * M_PI * 8.3187118e-02);

  printf("\nalpha_em = %.7e\nalpha_s = %.7e\nxw = %.7e\n\n", (g2 * sw) * (g2 * sw) / 4.0 / M_PI,
         g3 * g3 / 4.0 / M_PI, xw);
#else
  // Sets the weak angle.
  xw = pow(g1, 2) / (pow(g1, 2) + pow(g2, 2));
  sw = sqrt(xw);
  cw = sqrt(1.0 - xw);
#endif

  // Sets quark charges.
  double eq[6], T3q[6];
  for (int i0 = 0; i0 < 3; i0++) {
    eq[i0] = -1.0 / 3.0;
    eq[i0 + 3] = 2.0 / 3.0;
    T3q[i0] = -0.5;
    T3q[i0 + 3] = 0.5;
  }

  init_couplings();

  // Sets electroweak couplings.
  // Notation: [0] = photon, [1] = Z, [2] = W, [3] = Z', [4] = W'.
  for (int i0 = 0; i0 < 6; i0++) {
    vqq[0][i0][i0].L = -g2 * sw * eq[i0];
    vqq[0][i0][i0].R = vqq[0][i0][i0].L;
    vqq[1][i0][i0].L = -g2 / cw * (T3q[i0] - eq[i0] * xw);
    vqq[1][i0][i0].R = g2 / cw * eq[i0] * xw;
#ifdef SSM
    vqq[3][i0][i0].L = -g2 / cw * (T3q[i0] - eq[i0] * xw);
    vqq[3][i0][i0].R = g2 / cw * eq[i0] * xw;
#endif
  }
  for (int i0 = 0; i0 < 3; i0++) {
    for (int i1 = 0; i1 < 3; i1++) {
      vqq[2][i0 + 3][i1].L = -g2 * M_SQRT1_2 * ckm[i0][i1];
      vqq[2][i0][i1 + 3].L = -g2 * M_SQRT1_2 * conj(ckm[i1][i0]);

#ifdef SSM
      vqq[4][i0 + 3][i1].L = -g2 * M_SQRT1_2 * ckm[i0][i1];
      vqq[4][i0][i1 + 3].L = -g2 * M_SQRT1_2 * conj(ckm[i1][i0]);
#endif
    }
  }

  // Sets lepton couplings.
  for (int i0 = 0; i0 <= 2; i0++) {
    vll[0][i0 + 3][i0 + 3].L = g2 * sw;
    vll[0][i0 + 3][i0 + 3].R = vll[0][i0 + 3][i0 + 3].L;
    vll[1][i0][i0].L = -g2 / cw * 0.5;
    vll[1][i0 + 3][i0 + 3].L = -g2 / cw * (-0.5 + xw);
    vll[1][i0 + 3][i0 + 3].R = -g2 / cw * xw;
    vll[2][i0][i0 + 3].L = -g2 * M_SQRT1_2;
    vll[2][i0 + 3][i0].L = -g2 * M_SQRT1_2;
#ifdef SSM
    vll[3][i0][i0].L = -g2 / cw * 0.5;
    vll[3][i0 + 3][i0 + 3].L = -g2 / cw * (-0.5 + xw);
    vll[3][i0 + 3][i0 + 3].R = -g2 / cw * xw;
    vll[4][i0][i0 + 3].L = -g2 * M_SQRT1_2;
    vll[4][i0 + 3][i0].L = -g2 * M_SQRT1_2;
#endif
  }

  // Sets widths.
  Gv[0] = 0.0;    // y width
  Gv[1] = 2.4952; // Z0 width
  Gv[2] = 2.14;   // W width
#ifdef SSM
  Gv[3] = 0.0;
  Gv[4] = 0.0;
#endif

  //////////////////////////////////////////////////////////////
  // Calculating Z' and W' width
  /////////////////////////////////////////////////////////////
  {

    Gv[0] = 0.0; // y width
    Gv[1] = 0.0; // Z0 width
    Gv[2] = 0.0; // W width

    double Gamma = 0.0;

    // Alphas
    double rcor = 1.0 + pow(g3, 2.0) / (4.0 * M_PI * M_PI);
    double alp = pow(g2 * sw, 2.0) / (4.0 * M_PI);

    // Couplings for SSM from PYTHIA manual (in PYTHIA convention)
    double ZpqA[] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
    double ZpqV[] = {0.387, -0.693, 0.387, -0.693, 0.387, -0.693};
    double ZplA[] = {1.0, -1.0, 1.0, -1.0, 1.0, -1.0};
    double ZplV[] = {1.0, -0.08, 1.0, -0.08, 1.0, -0.08};

    // Gauge boson, quark and lepton masses
    double mzpa;
    double mqa[6] = {172.9, 4.0, 1.29, 100.0e-3, 1.7e-3, 0.0};
    double mla[6] = {0.0, 1776.0e-3, 0.0, 405.0e-3, 0.0, 0.0};

    // David
    // printf("\nZ' debugging parameters:\n1/alpha = %.7e\n", 1.0 / alp);

    // Z' width
    // mzpa = mv[3];
    mzpa = 4000;

    for (int i = 0; i < 6; i++) {
      // Quarks
      Gamma = Gamma + 3.0 * rcor * alp * sqrt(pow(mzpa, 4.0) - 4.0 * pow(mzpa * mqa[i], 2.0)) *
                          (pow(ZpqA[i], 2.0) * (pow(mzpa, 2.0) - 4.0 * pow(mqa[i], 2.0)) +
                           pow(ZpqV[i], 2.0) * (pow(mzpa, 2.0) + 2.0 * pow(mqa[i], 2.0))) /
                          (48.0 * (1.0 - xw) * xw * pow(mzpa, 3.0));

      // Leptons
      Gamma = Gamma + alp * sqrt(pow(mzpa, 4.0) - 4.0 * pow(mzpa * mla[i], 2.0)) *
                          (pow(ZplA[i], 2.0) * (pow(mzpa, 2.0) - 4.0 * pow(mla[i], 2.0)) +
                           pow(ZplV[i], 2.0) * (pow(mzpa, 2.0) + 2.0 * pow(mla[i], 2.0))) /
                          (48.0 * (1.0 - xw) * xw * pow(mzpa, 3.0));
    }

    // David
    // printf("Z' width: %.7e\n", Gamma);

    Gv[3] = Gamma; // Z' width

    // W' width: this is actually for Z', but it seems to work
    Gamma = 0.0;
    mzpa = mv[4];

    for (int i = 0; i < 6; i++) {
      // Quarks
      Gamma = Gamma + 3.0 * rcor * alp * sqrt(pow(mzpa, 4.0) - 4.0 * pow(mzpa * mqa[i], 2.0)) *
                          (pow(ZpqA[i], 2.0) * (pow(mzpa, 2.0) - 4.0 * pow(mqa[i], 2.0)) +
                           pow(ZpqV[i], 2.0) * (pow(mzpa, 2.0) + 2.0 * pow(mqa[i], 2.0))) /
                          (48.0 * (1.0 - xw) * xw * pow(mzpa, 3.0));

      // Leptons
      Gamma = Gamma + alp * sqrt(pow(mzpa, 4.0) - 4.0 * pow(mzpa * mla[i], 2.0)) *
                          (pow(ZplA[i], 2.0) * (pow(mzpa, 2.0) - 4.0 * pow(mla[i], 2.0)) +
                           pow(ZplV[i], 2.0) * (pow(mzpa, 2.0) + 2.0 * pow(mla[i], 2.0))) /
                          (48.0 * (1.0 - xw) * xw * pow(mzpa, 3.0));
    }

    // David
    // printf("W' width: %.7e\n\n", Gamma);

    Gv[4] = Gamma; // Z' width
  }
  /////////////////////////////////////////////////////////////////////////

#ifdef SUSY
  // Sets SUSY couplings.
  for (int i0 = 0; i0 < 12; i0++) {
    vSQSQ[0][i0][i0].L = vqq[0][i0 / 2][i0 / 2].L;
    vSQSQ[0][i0][i0].R = vqq[0][i0 / 2][i0 / 2].R;
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        vSQSQ[1][i0][i1].R += vqq[1][i0 / 2][i1 / 2].L * RSD[i0][i2] * conj(RSD[i1][i2]) +
                              vqq[1][i0 / 2][i1 / 2].R * RSD[i0][i2 + 3] * conj(RSD[i1][i2 + 3]);
        vSQSQ[1][i0 + 6][i1 + 6].R +=
            vqq[1][i0 / 2 + 3][i1 / 2 + 3].L * RSU[i0][i2] * conj(RSU[i1][i2]) +
            vqq[1][i0 / 2 + 3][i1 / 2 + 3].R * RSU[i0][i2 + 3] * conj(RSU[i1][i2 + 3]);
        for (int i3 = 0; i3 < 3; i3++) {
          vSQSQ[2][i0 + 6][i1].R += vqq[2][i2 + 3][i3].L * RSU[i0][i2] * conj(RSD[i1][i3]);
          vSQSQ[2][i1][i0 + 6].R += vqq[2][i2][i3 + 3].L * conj(RSU[i0][i2]) * RSD[i1][i3];
        }
      }
      // BENJ FIX: some couplings are undefined
      vSQSQ[1][i0][i1].L = vSQSQ[1][i0][i1].R;
      vSQSQ[1][i0 + 6][i1 + 6].L = vSQSQ[1][i0 + 6][i1 + 6].R;
      vSQSQ[2][i0][i1].L = vSQSQ[2][i0][i1].R;
      vSQSQ[2][i0 + 6][i1 + 6].L = vSQSQ[2][i0 + 6][i1 + 6].R;
    }
  }

  for (int i0 = 0; i0 < 4; i0++) {
    for (int i1 = 0; i1 < 4; i1++) {
      vCHCH[1][i0][i1].L =
          g2 / cw * .5 * (-RN[i0][2] * conj(RN[i1][2]) + RN[i0][3] * conj(RN[i1][3]));
      // no conjugation below missing as i1 <-> i0 swapped cf. 9511250
      vCHCH[1][i1][i0].R = -vCHCH[1][i0][i1].L;
    }
  }

  for (int i0 = 0; i0 < 4; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      vCHCH[2][i0][i1 + 4].L =
          g2 * (RN[i0][1] * conj(RV[i1][0]) - RN[i0][3] * conj(RV[i1][1]) * M_SQRT1_2);
      vCHCH[2][i0][i1 + 4].R =
          g2 * (conj(RN[i0][1]) * RU[i1][0] + conj(RN[i0][2]) * RU[i1][1] * M_SQRT1_2);

      for (int i2 = 0; i2 < 3; i2++) {
        vCHCH[i2][i1 + 4][i0].L = conj(vCHCH[i2][i0][i1 + 4].L);
        vCHCH[i2][i1 + 4][i0].R = conj(vCHCH[i2][i0][i1 + 4].R);
      }
    }
  }

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      double delta = (i0 == i1 ? 1.0 : 0.0);
      vCHCH[0][i0 + 4][i1 + 4].L = -g2 * sw * delta;
      vCHCH[0][i0 + 4][i1 + 4].R = -g2 * sw * delta;
      vCHCH[1][i0 + 4][i1 + 4].L =
          g2 / cw * (-RV[i0][0] * conj(RV[i1][0]) - .5 * RV[i0][1] * conj(RV[i1][1]) + delta * xw);
      vCHCH[1][i0 + 4][i1 + 4].R =
          g2 / cw * (-conj(RU[i0][0]) * RU[i1][0] - .5 * conj(RU[i0][1]) * RU[i1][1] + delta * xw);
    }
  }

  // Slepton couplings.
  for (int i0 = 3; i0 < 9; i0++) {
    vSLSL[0][i0][i0].L = g2 * sw;
    vSLSL[0][i0][i0].R = vSLSL[0][i0][i0].L;
  }

  for (int i0 = 3; i0 < 9; i0++) {
    if (i0 % 2 != 0) {
      vSLSL[1][i0][i0].L = g2 / (2.0 * cw) * (1 - 2.0 * xw); // L-Type
      vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
    }

    if (i0 % 2 == 0) {
      vSLSL[1][i0][i0].L = -g2 / (2.0 * cw) * 2.0 * xw; // R-Type
      vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
    }
  }

  for (int i0 = 7; i0 < 9; i0++) {
    for (int i1 = 7; i1 < 9; i1++) {
      vSLSL[1][i0][i1].L = -g2 / (2.0 * cw) *
                           (2.0 * xw * SSLM[i0 - 7][1] * SSLM[i1 - 7][1] -
                            SSLM[i0 - 7][0] * SSLM[i1 - 7][0] * (1.0 - 2.0 * xw));
      vSLSL[1][i0][i1].R = vSLSL[1][i0][i1].L;
    }
  }

  for (int i0 = 0; i0 < 3; i0++) {
    vSLSL[1][i0][i0].L = -g2 / (2 * cw);
    vSLSL[1][i0][i0].R = vSLSL[1][i0][i0].L;
  }

  for (int i0 = 0; i0 < 2; i0++) {
    vSLSL[2][i0][i0 + (i0 + 3)].L = g2 * (1.0 / sqrt(2));
    vSLSL[2][i0][i0 + (i0 + 3)].R = g2 * (1.0 / sqrt(2));
    vSLSL[2][i0 + (i0 + 3)][i0].L = g2 * (1.0 / sqrt(2));
    vSLSL[2][i0 + (i0 + 3)][i0].R = g2 * (1.0 / sqrt(2));
  }

  vSLSL[2][7][2].L = g2 * (1.0 / sqrt(2)) * SSLM[0][1];
  vSLSL[2][7][2].R = g2 * (1.0 / sqrt(2)) * SSLM[0][1];
  vSLSL[2][8][2].L = g2 * (1.0 / sqrt(2)) * SSLM[1][1];
  vSLSL[2][8][2].R = g2 * (1.0 / sqrt(2)) * SSLM[1][1];
  vSLSL[2][2][7].L = g2 * (1.0 / sqrt(2)) * SSLM[0][1];
  vSLSL[2][2][7].R = g2 * (1.0 / sqrt(2)) * SSLM[0][1];
  vSLSL[2][2][8].L = g2 * (1.0 / sqrt(2)) * SSLM[1][1];
  vSLSL[2][2][8].R = g2 * (1.0 / sqrt(2)) * SSLM[1][1];

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 6; i2++) {
        for (int i3 = 0; i3 < 3; i3++) {
          hSQSQ[i0][i1][i2].R +=
              mv[1] * (RA[i0][0] * RB[0][1] - RA[i0][1] * RB[0][0]) *
              (vqq[1][i1 / 2][i2 / 2].L * RSD[i1][i3] * conj(RSD[i2][i3]) -
               vqq[1][i1 / 2][i2 / 2].R * RSD[i1][i3 + 3] * conj(RSD[i2][i3 + 3]));
          hSQSQ[i0][i1 + 6][i2 + 6].R +=
              mv[1] * (RA[i0][0] * RB[0][1] - RA[i0][1] * RB[0][0]) *
              (vqq[1][i1 / 2 + 3][i2 / 2 + 3].L * RSU[i1][i3] * conj(RSU[i2][i3]) -
               vqq[1][i1 / 2 + 3][i2 / 2 + 3].R * RSU[i1][i3 + 3] * conj(RSU[i2][i3 + 3]));
        }
      }
    }
  }

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 6; i2++) {
        for (int i3 = 0; i3 < 3; i3++) {
          for (int i4 = 0; i4 < 3; i4++) {
            for (int i5 = 0; i5 < 2; i5++) {
              hSQSQ[i0 + 4][i1 + 6][i2].R += mv[2] * vqq[2][i3 + 3][i4].L * RSU[i1][i3] *
                                             RSD[i2][i4] * RB[0][1 - i5] * RB[1][i5];
            }
          }
        }
        hSQSQ[i0 + 4][i2][i1 + 6].R = conj(hSQSQ[i0 + 4][i1 + 6][i2].R);
      }
    }
  }

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 4; i1++) {
      for (int i2 = 0; i2 < 4; i2++) {
        hCHCH[i0][i1][i2].R =
            .5 * g2 / cw *
            ((RA[i0][0] * RN[i1][2] - RA[i0][1] * RN[i1][3]) * (RN[i2][0] * sw - RN[i2][1] * cw) +
             (RA[i0][0] * RN[i2][2] - RA[i0][1] * RN[i2][3]) * (RN[i1][0] * sw - RN[i1][1] * cw));
        hCHCH[i0][i2][i1].L = conj(hCHCH[i0][i1][i2].R);
        hCHCH[i0 + 2][i1][i2].R =
            .5 * g2 / cw * II *
            ((RB[i0][0] * RN[i1][2] - RB[i0][1] * RN[i1][3]) * (RN[i2][0] * sw - RN[i2][1] * cw) +
             (RB[i0][0] * RN[i2][2] - RB[i0][1] * RN[i2][3]) * (RN[i1][0] * sw - RN[i1][1] * cw));
        hCHCH[i0 + 2][i2][i1].L = conj(hCHCH[i0 + 2][i1][i2].L);
      }
    }
  }

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      for (int i2 = 0; i2 < 2; i2++) {
        hCHCH[i0][i1 + 4][i2 + 4].R =
            -g2 * M_SQRT1_2 *
            (RA[i0][0] * RV[i1][0] * RU[i2][1] + RA[i0][1] * RV[i1][1] * RU[i2][0]);
        hCHCH[i0][i2 + 4][i1 + 4].L = conj(hCHCH[i0][i1 + 4][i2 + 4].R);
        hCHCH[i0 + 2][i1 + 4][i2 + 4].R =
            -g2 * M_SQRT1_2 * II *
            (RB[i0][0] * RV[i1][0] * RU[i2][1] + RB[i0][1] * RV[i1][1] * RU[i2][0]);
        hCHCH[i0 + 2][i2 + 4][i1 + 4].L = conj(hCHCH[i0 + 2][i1 + 4][i2 + 4].R);
      }
    }
  }

  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      for (int i2 = 0; i2 < 4; i2++) {
        hCHCH[i0 + 4][i1 + 4][i2].L =
            g2 * RB[i0][0] *
            conj(M_SQRT1_2 * RU[i1][1] * (RN[i2][0] * sw / cw + RN[i2][1]) - RU[i1][0] * RN[i2][2]);
        hCHCH[i0 + 4][i1 + 4][i2].R =
            -g2 * RB[i0][1] *
            (M_SQRT1_2 * RV[i1][1] * (RN[i2][0] * sw / cw + RN[i2][1]) + RV[i1][0] * RN[i2][3]);
        hCHCH[i0 + 4][i2][i1 + 4].L = conj(hCHCH[i0 + 4][i1 + 4][i2].R);
        hCHCH[i0 + 4][i2][i1 + 4].R = conj(hCHCH[i0 + 4][i1 + 4][i2].L);
      }
    }
  }

  // neutralino - squark - quark
  for (int i0 = 0; i0 < 4; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        CHSQq[i0][i1][i2].L =
            -g2 * M_SQRT2 * RSD[i1][i2] *
                ((eq[i2] - T3q[i2]) * sw / cw * conj(RN[i0][0]) + T3q[i2] * conj(RN[i0][1])) -
            yq[i2] * conj(RN[i0][2]) * RSD[i1][i2 + 3];

        CHSQq[i0][i1][i2].R = g2 * M_SQRT2 * RSD[i1][i2 + 3] * eq[i2] * sw / cw * RN[i0][0] -
                              yq[i2] * RN[i0][2] * RSD[i1][i2];

        CHSQq[i0][i1 + 6][i2 + 3].L = -g2 * M_SQRT2 * RSU[i1][i2] *
                                          ((eq[i2 + 3] - T3q[i2 + 3]) * sw / cw * conj(RN[i0][0]) +
                                           T3q[i2 + 3] * conj(RN[i0][1])) -
                                      yq[i2 + 3] * conj(RN[i0][3]) * RSU[i1][i2 + 3];

        CHSQq[i0][i1 + 6][i2 + 3].R =
            g2 * M_SQRT2 * RSU[i1][i2 + 3] * eq[i2 + 3] * sw / cw * RN[i0][0] -
            yq[i2 + 3] * RN[i0][3] * RSU[i1][i2];
      }
    }
  }

  // chargino quark squark
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        for (int i3 = 0; i3 < 3; i3++) {
          CHSQq[i0 + 4][i1 + 6][i2].L += (-g2 * conj(RV[i0][0]) * RSU[i1][i3] +
                                          yq[i3 + 3] * conj(RV[i0][1]) * RSU[i1][i3 + 3]) *
                                         ckm[i3][i2];
          CHSQq[i0 + 4][i1 + 6][i2].R += yq[i2] * RU[i0][1] * RSU[i1][i3] * ckm[i3][i2];

          CHSQq[i0 + 4][i1][i2 + 3].L +=
              (-g2 * conj(RU[i0][0]) * RSD[i1][i3] + yq[i3] * conj(RU[i0][1]) * RSD[i1][i3 + 3]) *
              conj(ckm[i2][i3]);

          CHSQq[i0 + 4][i1][i2 + 3].R += yq[i2 + 3] * RV[i0][1] * RSD[i1][i3] * conj(ckm[i2][i3]);
        }
      }
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 12; i1++) {
      for (int i2 = 0; i2 < 6; i2++) {
        CHqSQ[i0][i2][i1].L = conj(CHSQq[i0][i1][i2].R);
        CHqSQ[i0][i2][i1].R = conj(CHSQq[i0][i1][i2].L);
      }
    }
  }
#endif

  // Strong couplings.
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      gqq[i0][i1].L = (i0 == i1 ? -g3 : 0.0);
      gqq[i0][i1].R = (i0 == i1 ? -g3 : 0.0);
    }
  }

#ifdef SUSY
  for (int i0 = 0; i0 < 12; i0++) {
    for (int i1 = 0; i1 < 12; i1++) {
      gSQSQ[i0][i1].L = (i0 == i1 ? -g3 : 0.0);
      gSQSQ[i0][i1].R = (i0 == i1 ? -g3 : 0.0);
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 3; i1++) {
      GLSQq[i0][i1].L = -g3 * M_SQRT2 * RSD[i0][i1];
      GLSQq[i0][i1].R = g3 * M_SQRT2 * RSD[i0][i1 + 3];
      GLSQq[i0][i1 + 3].L = 0.0;
      GLSQq[i0][i1 + 3].R = 0.0;
      GLSQq[i0 + 6][i1].L = 0.0;
      GLSQq[i0 + 6][i1].R = 0.0;
      GLSQq[i0 + 6][i1 + 3].L = -g3 * M_SQRT2 * RSU[i0][i1];
      GLSQq[i0 + 6][i1 + 3].R = g3 * M_SQRT2 * RSU[i0][i1 + 3];
    }
  }

  for (int i0 = 0; i0 < 12; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      GLqSQ[i1][i0].L = conj(GLSQq[i0][i1].R);
      GLqSQ[i1][i0].R = conj(GLSQq[i0][i1].L);
    }
  }

  // Effective vertex for SQ-SQloop-SQ.
  complex<double> SD[6][6], SU[6][6];
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 3; i2++) {
        SD[i0][i1] += conj(RSD[i0][i2]) * RSD[i1][i2] - conj(RSD[i0][i2 + 3]) * RSD[i1][i2 + 3];
        SU[i0][i1] += conj(RSU[i0][i2]) * RSU[i1][i2] - conj(RSU[i0][i2 + 3]) * RSU[i1][i2 + 3];
      }
    }
  }

  for (int i0 = 0; i0 < 12; i0++) {
    for (int i1 = 0; i1 < 12; i1++) {
      for (int i2 = 0; i2 < 12; i2++) {
        SQSQSQ[i0][i1][i2].L = complex<double>(0.0, 0.0);
        SQSQSQ[i0][i1][i2].R = complex<double>(0.0, 0.0);
      }
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      for (int i2 = 0; i2 < 6; i2++) {
        SQSQSQ[i0][i1][i2].R = -pow(g3, 2) * SD[i1][i0] * conj(SD[i1][i2]);
        SQSQSQ[i0 + 6][i1 + 6][i2 + 6].R = -pow(g3, 2) * SU[i1][i0] * conj(SU[i1][i2]);
      }
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {          // CH
    for (int i1 = 0; i1 < 12; i1++) {       // SQ
      for (int i2 = 0; i2 < 6; i2++) {      // q
        for (int i3 = 0; i3 < 12; i3++) {   // SQ
          for (int i4 = 0; i4 < 12; i4++) { // SQ
            if (abs(mSQ[i1] - mSQ[i3]) / abs(mSQ[i1] + mSQ[i3]) > 1.0e-7) {
              dCHSQq[i0][i1][i2].L += -pow(mSQ[i4], 2) * SQSQSQ[i1][i4][i3].R /
                                      (pow(mSQ[i1], 2) - pow(mSQ[i3], 2)) * CHSQq[i0][i3][i2].L;
              dCHSQq[i0][i1][i2].R += -pow(mSQ[i4], 2) * SQSQSQ[i1][i4][i3].R /
                                      (pow(mSQ[i1], 2) - pow(mSQ[i3], 2)) * CHSQq[i0][i3][i2].R;
            }
          }
        }
      }
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 12; i1++) {
      for (int i2 = 0; i2 < 6; i2++) {
        dCHqSQ[i0][i2][i1].L = conj(dCHSQq[i0][i1][i2].R);
        dCHqSQ[i0][i2][i1].R = conj(dCHSQq[i0][i1][i2].L);
      }
    }
  }
  /**
   * Calculate squark squark squark squark vertex from complete SUSY feynrules
   * paper.
   * arXiv:hep-ph/9511250
   */
  TensorMap<ComplexType, 6, 6> tXU((ComplexType *)XU);
  TensorMap<ComplexType, 6, 6, 6, 6> tXUU((ComplexType *)XUU);
  Tensor<ComplexType, 6, 6> Su((ComplexType *)RSU);
  tXU.eye();
  tXU(all, all) -=
      2 * einsum<Index<I, L>, Index<I, J>>(conj(Su(fseq<0, 3>(), all)), Su(fseq<0, 3>(), all));
  tXUU(all, all, all, all) =
      3 * permute<Index<I, L, K, J>>(einsum<Index<I, L>, Index<K, J>>(tXU, tXU));
  tXUU(all, all, all, all) -=
      permute<Index<I, J, K, L>>(einsum<Index<I, J>, Index<K, L>>(tXU, tXU));

  TensorMap<ComplexType, 6, 6> tXD((ComplexType *)XD);
  TensorMap<ComplexType, 6, 6, 6, 6> tXDD((ComplexType *)XDD);
  Tensor<ComplexType, 6, 6> Sd((ComplexType *)RSD);
  tXD.eye();
  tXD(all, all) -=
      2 * einsum<Index<I, L>, Index<I, J>>(conj(Sd(fseq<0, 3>(), all)), Sd(fseq<0, 3>(), all));
  tXDD(all, all, all, all) =
      3 * permute<Index<I, L, K, J>>(einsum<Index<I, L>, Index<K, J>>(tXD, tXD));
  tXDD(all, all, all, all) -=
      permute<Index<I, J, K, L>>(einsum<Index<I, J>, Index<K, L>>(tXD, tXD));

#endif
}

void Parameters::read_slha(const char *file) {
  using namespace std;
  using namespace SLHAea;

  ifstream *ifs;

  if (strcmp(file, "-") == 0) {
    ifs = new ifstream("/dev/stdin");
  } else {
    ifs = new ifstream(file);
  }

  if (!(*ifs).good()) {
    fprintf(stderr, "error: Failed to open '%s'\n", file);
    exit(1);
  }

  Coll input(*ifs);
  try {
    mv[0] = 0.0; // photon mass
    try {
      mv[1] = to<double>(input.at("SMINPUTS").at("4").at(1)); // Z mass
    } catch (...) {
      mv[1] = to<double>(input.at("MASS").at("23").at(1)); // Z mass
    }
    mvs[1] = pow2(mv[1]);
    //
    mv[2] = to<double>(input.at("MASS").at("24").at(1)); // W mass
    mvs[2] = pow2(mv[2]);
  } catch (...) {
    fprintf(stderr, "error: No vector boson masses defined.\n");
    exit(1);
  }
  try {
    try {
      mq[5] = to<double>(input.at("SMINPUTS").at("6").at(1)); // top mass
    } catch (...) {
      mq[5] = to<double>(input.at("MASS").at("6").at(1)); // top mass
    }
    mqs[5] = pow2(mq[5]);
  } catch (...) {
    fprintf(stderr, "warning: No top mass defined.\n");
  }
  for (int i = 0; i < 5; i++) {
    mq[i] = 0.0;
    mqs[i] = 0.0;
  }
  try {
    mh[1] = to<double>(input.at("MASS").at("25").at(1)); // h mass
    mh[0] = to<double>(input.at("MASS").at("35").at(1)); // H mass
    mh[2] = to<double>(input.at("MASS").at("36").at(1)); // A mass
    mh[4] = to<double>(input.at("MASS").at("37").at(1)); // H+ mass
    for (int i = 0; i <= 4; i++) {
      mhs[i] = pow2(mh[i]);
    }
  } catch (...) {
    fprintf(stderr, "warning: No Higgs masses defined.\n");
  }
  try {
    mGL = to<double>(input.at("MASS").at("1000021").at(1)); // gluino mass
    mGLs = pow2(mGL);
  } catch (...) {
    fprintf(stderr, "warning: No gluino masses defined.\n");
  }
  // Squark masses
  try {
    // down type left-handed; then down typ right-handed squark masses.
    mSQ[0] = to<double>(input.at("MASS").at("1000001").at(1));
    mSQ[1] = to<double>(input.at("MASS").at("1000003").at(1));
    mSQ[2] = to<double>(input.at("MASS").at("1000005").at(1));
    mSQ[3] = to<double>(input.at("MASS").at("2000001").at(1));
    mSQ[4] = to<double>(input.at("MASS").at("2000003").at(1));
    mSQ[5] = to<double>(input.at("MASS").at("2000005").at(1));

    // up type left-handed; then up type right-handed squark masses.
    mSQ[6] = to<double>(input.at("MASS").at("1000002").at(1));
    mSQ[7] = to<double>(input.at("MASS").at("1000004").at(1));
    mSQ[8] = to<double>(input.at("MASS").at("1000006").at(1));
    mSQ[9] = to<double>(input.at("MASS").at("2000002").at(1));
    mSQ[10] = to<double>(input.at("MASS").at("2000004").at(1));
    mSQ[11] = to<double>(input.at("MASS").at("2000006").at(1));

    for (int i = 0; i <= 11; i++) {
      mSQs[i] = pow2(mSQ[i]);
    }
  } catch (...) {
    fprintf(stderr, "warning: No squark masses defined.\n");
  }
  try {
    // Neutralino masses
    mCH[0] = to<double>(input.at("MASS").at("1000022").at(1));
    mCH[1] = to<double>(input.at("MASS").at("1000023").at(1));
    mCH[2] = to<double>(input.at("MASS").at("1000025").at(1));
    mCH[3] = to<double>(input.at("MASS").at("1000035").at(1));
    // Chargino masses
    mCH[4] = to<double>(input.at("MASS").at("1000024").at(1));
    mCH[5] = to<double>(input.at("MASS").at("1000037").at(1));
    for (int i = 0; i <= 5; i++) {
      mCHs[i] = pow2(mCH[i]);
    }
  } catch (...) {
    fprintf(stderr, "warning: No gaugino masses defined.\n");
  }
  try {
    // Sneutrino masses
    mSL[0] = to<double>(input.at("MASS").at("1000012").at(1));
    mSL[1] = to<double>(input.at("MASS").at("1000014").at(1));
    mSL[2] = to<double>(input.at("MASS").at("1000016").at(1));
    // Charged Sleptons masses
    mSL[3] = to<double>(input.at("MASS").at("1000011").at(1));
    mSL[4] = to<double>(input.at("MASS").at("2000011").at(1));
    mSL[5] = to<double>(input.at("MASS").at("1000013").at(1));
    mSL[6] = to<double>(input.at("MASS").at("2000013").at(1));
    mSL[7] = to<double>(input.at("MASS").at("1000015").at(1));
    mSL[8] = to<double>(input.at("MASS").at("2000015").at(1));
    for (int i = 0; i <= 8; i++) {
      mSLs[i] = pow2(mSL[i]);
    }
  } catch (...) {
    fprintf(stderr, "warning: No slepton masses defined.\n");
  }

  // Mixings

  // Neutralino mixing
  try {
    for (int i0 = 0; i0 < 4; i0++) {
      for (int i1 = 0; i1 < 4; i1++) {
        RN[i0][i1] = to<double>(input.at("NMIX").at(i0 + 1, i1 + 1).at(2));
      }
    }

    for (int i0 = 0; i0 < 4; i0++) {
      if (mCH[i0] < 0.0) {
        mCH[i0] *= -1.0;
        for (int i1 = 0; i1 < 4; i1++) {
          RN[i0][i1] *= complex<double>(0.0, 1.0);
        }
      }
    }
  } catch (...) {
    fprintf(stderr, "warning: No neutralino mixing defined.\n");
  }

  // Chargino mixing matrix U
  try {
    for (int i0 = 0; i0 < 2; i0++) {
      for (int i1 = 0; i1 < 2; i1++) {
        RU[i0][i1] = to<double>(input.at("UMIX").at(i0 + 1, i1 + 1).at(2));
      }
    }

    // Chargino mixing matrix V
    for (int i0 = 0; i0 < 2; i0++) {
      for (int i1 = 0; i1 < 2; i1++) {
        RV[i0][i1] = to<double>(input.at("VMIX").at(i0 + 1, i1 + 1).at(2));
      }
    }
  } catch (...) {
    fprintf(stderr, "warning: No chargino mixing defined.\n");
  }

  // Sbottom mixing
  try {
    complex<double> RSB[2][2];
    for (int i0 = 0; i0 < 2; i0++) {
      for (int i1 = 0; i1 < 2; i1++) {
        RSB[i0][i1] = to<double>(input.at("SBOTMIX").at(i0 + 1, i1 + 1).at(2));
        // Set diagonal (off-diagonal elements proportional to m_b = 0 )
        // only sbottom mixing for gaugino pair production implemented
        if (!is_gaugino_gaugino(out1, out2)) {
          if (i0 == i1) {
            RSB[i0][i1] = 1.0;
          } else {
            RSB[i0][i1] = 0.0;
          }
        }
#ifdef NO_SQUARK_MIXING
        if (i0 == i1) {
          RSB[i0][i1] = 1.0;
        } else {
          RSB[i0][i1] = 0.0;
        }

#endif
      }
    }

    for (int i0 = 0; i0 < 6; i0++) {
      for (int i1 = 0; i1 < 6; i1++) {
        if (i0 % 3 == 2 && i1 % 3 == 2) {
          RSD[i0][i1] = RSB[i0 / 3][i1 / 3];
        } else {
          RSD[i0][i1] = (i0 == i1 ? 1.0 : 0.0);
        }
      }
    }
  } catch (...) {
    fprintf(stderr, "warning: No sbottom mixing defined.\n");
  }

  // Stop mixing
  try {
    complex<double> RST[2][2];
    for (int i0 = 0; i0 < 2; i0++) {
      for (int i1 = 0; i1 < 2; i1++) {
        RST[i0][i1] = to<double>(input.at("STOPMIX").at(i0 + 1, i1 + 1).at(2));

#ifdef NO_SQUARK_MIXING
        if (i0 == i1) {
          RST[i0][i1] = 1.0;
        } else {
          RST[i0][i1] = 0.0;
        }
#endif
      }
    }

    // Heavy squark mixing matrix
    for (int i0 = 0; i0 < 6; i0++) {
      for (int i1 = 0; i1 < 6; i1++) {
        if (i0 % 3 == 2 && i1 % 3 == 2) {
          RSU[i0][i1] = RST[i0 / 3][i1 / 3];
        } else {
          RSU[i0][i1] = (i0 == i1 ? 1.0 : 0.0);
        }
      }
    }
  } catch (...) {
    fprintf(stderr, "warning: No stop mixing defined.\n");
  }

#ifdef NO_SQUARK_MIXING
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      if (i0 == i1) {
        RSU[i0][i1] = 1.0;
      } else {
        RSU[i0][i1] = 0.0;
      }
    }
  }

  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      if (i0 == i1) {
        RSD[i0][i1] = 1.0;
      } else {
        RSD[i0][i1] = 0.0;
      }
    }
  }
#endif

  // Alpha
  try {
    double alpha = 0.0;
    for (Block::const_iterator line = input.at("ALPHA").begin(); line != input.at("ALPHA").end();
         ++line) {
      if (!line->is_data_line()) {
        continue;
      }

      alpha += to<double>(line->at(0));
    }

    RA[0][0] = cos(alpha);
    RA[0][1] = sin(alpha);
    RA[1][0] = -sin(alpha);
    RA[1][1] = cos(alpha);
  } catch (...) {
    fprintf(stderr, "warning: No block ALPHA defined.\n");
  }

  // Higgisino mixing
  try {
    double beta;
    beta = to<double>(input.at("HMIX").at("2").at(1));
    beta = atan(beta);
    RB[0][0] = sin(beta);
    RB[0][1] = cos(beta);
    RB[1][0] = -cos(beta);
    RB[1][1] = sin(beta);
  } catch (...) {
    fprintf(stderr, "warning: No Higgsino mixing defined. (HMIX)\n");
  }

  // TODO: add correct CKM Matrix.
  // Use wolfenstein parameters.
  // (do not forget the CP-violating complex phase)

  ComplexType Vckm[3][3];
  Vckm[0][0] = 0.97425;
  Vckm[0][1] = 0.2253;
  Vckm[0][2] = 4.13E-3;
  Vckm[1][0] = 0.225;
  Vckm[1][1] = 0.986;
  Vckm[1][2] = 41.1E-3;
  Vckm[2][0] = 8.4E-3;
  Vckm[2][1] = 40.0E-3;
  Vckm[2][2] = 1.021;

// From arXiv:0801.0045 [hep-ph] 1106.0935 0406184
#ifdef CKM_LOAD
  try {
    double lambda = to<double>(input.at("VCKMIN").at("1"));
    double cap_A = to<double>(input.at("VCKMIN").at("2"));
    double rho_bar = to<double>(input.at("VCKMIN").at("3"));
    double eta_bar = to<double>(input.at("VCKMIN").at("4"));

    double rho = rho_bar / (1. - pow2(lambda) / 2.);
    double eta = eta_bar / (1. - pow2(lambda) / 2.);

    Vckm[0][0] = 1. - pow2(lambda) / 2.;
    Vckm[0][1] = lambda;
    Vckm[0][2] = cap_A * std::pow(lambda, 3) * (rho - i * eta);
    Vckm[1][0] = -lambda;
    Vckm[1][1] = 1. - pow2(lambda) / 2.;
    Vckm[1][2] = cap_A * pow2(lambda);
    Vckm[2][0] = cap_A * std::pow(lambda, 3) * (1 - rho - i * eta);
    Vckm[2][1] = -cap_A * pow2(lambda);
    Vckm[2][2] = 1.;

  } catch (...) {
    fprintf(stderr, "warning: No CKM mixing defined. (VCKMIN)\n");
    Vckm[0][0] = 0.97425;
    Vckm[0][1] = 0.2253;
    Vckm[0][2] = 4.13E-3;
    Vckm[1][0] = 0.225;
    Vckm[1][1] = 0.986;
    Vckm[1][2] = 41.1E-3;
    Vckm[2][0] = 8.4E-3;
    Vckm[2][1] = 40.0E-3;
    Vckm[2][2] = 1.021;
  }
#endif

  for (int i = 0; i <= 2; i++) {
    for (int j = 0; j <= 2; j++) {
      ckm[i][j] = Vckm[i][j];
    }
  }

// Diagonal CKM matrix.
#ifdef CKM_DIAG
  for (int i0 = 0; i0 < 3; i0++) {
    for (int i1 = 0; i1 < 3; i1++) {
      ckm[i0][i1] = (i0 == i1 ? 1.0 : 0.0);
    }
  }
#endif

  // Gauge coupling constants.
  double g[3];

  // sw = sqrt( 1.0 - mw**2/mz**2 )
  // sw2 = sw**2
  // cw2 = 1.0 - sw2
  // cw  = sqrt(cw2)

  double a_s, G_F;

  try {
#ifdef GAUGE_COUPLING
    g[0] = to<double>(input.at("GAUGE").at("1").at(1));
    g[1] = to<double>(input.at("GAUGE").at("2").at(1));
    g[2] = to<double>(input.at("GAUGE").at("3").at(1));
    g1 = g[0];
    g2 = g[1];
    g3 = g[2];
    // readded and needed
    xw = 1.0 - pow(mv[2] / mv[1], 2);
    sw = sqrt(xw);
    cw = sqrt(1.0 - xw);
#endif
    // g3 = 1.0; // remove that later and get rid of g3 after testing.
#ifdef PROSPINO_COUPLING
    //G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
    a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
    G_F = 1.166379e-5;
    xw = 1.0 - pow(mv[2] / mv[1], 2);
    sw = std::sqrt(xw);
    cw = sqrt(1.0 - xw);
    g2 = 2.0 * mv[2] * sqrt(M_SQRT2 * G_F);
    g3 = sqrt(4 * M_PI * a_s);
    g1 = sw / cw * g2; // not used
#endif
#ifdef MADGRAPH_COUPLING_OLD
    G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
    a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
    double a_EM = 1.0 / to<double>(input.at("SMINPUTS").at("1").at(1));
    mv[1] = to<double>(input.at("SMINPUTS").at("4").at(1));
    sw = sqrt(0.5 *
              (1.0 - sqrt(1.0 - (4.0 * (M_PI) / (sqrt(2) * G_F * (1.0 / a_EM) * pow2(mv[1]))))));
    xw = pow2(sw);
    // xw = 1.0 - pow(mv[2] / mv[1], 2);
    // sw = sqrt(xw);
    // G_F = 1.16639e-5;
    cw = sqrt(1.0 - xw);
    // g2 = 2.0 * mv[2] * sqrt(M_SQRT2 * G_F);
    // g3 = sqrt(4 * M_PI * a_s);
    mv[2] = mv[1] * cw;
    g2 = 2.0 * mv[2] * sqrt((M_PI)*a_EM / (xw * pow2(mv[2])));
    g3 = sqrt(4 * (M_PI)*a_s);
    g1 = sw / cw * g2; // not used
#endif
#ifdef MADGRAPH_COUPLING
    G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
    a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
    double a_EM = 1.0 / to<double>(input.at("SMINPUTS").at("1").at(1));
    // mv[1] = to<double>(input.at("SMINPUTS").at("4").at(1));
    // sw = sqrt(0.5 *
    //           (1.0 - sqrt(1.0 - (4.0 * (M_PI) / (sqrt(2) * G_F * (1.0 / a_EM) * pow2(mv[1]))))));
    cw = mv[2] / mv[1];
    sw = sqrt(1. - pow2(cw));
    xw = pow2(sw);
    // xw = 1.0 - pow(mv[2] / mv[1], 2);
    // sw = sqrt(xw);
    // G_F = 1.16639e-5;
    // cw = sqrt(1.0 - xw);
    double ee = 2. * sqrt(a_EM) * sqrt(M_PI);
    // g2 = 2.0 * mv[2] * sqrt(M_SQRT2 * G_F);
    // g3 = sqrt(4 * M_PI * a_s);
    mv[2] = mv[1] * cw;
    // g2 = 2.0 * mv[2] * sqrt((M_PI)*a_EM / (xw * pow2(mv[2])));
    g2 = ee / sw;
    g3 = sqrt(4 * (M_PI)*a_s);
    g1 = sw / cw * g2; // not used

    // Set Yukawa couplings to 0 for (d,s,b,u,c).
    for (int i0 = 0; i0 < 5; i0++) {
      yq[i0] = 0.0;
    }

    double vev = (2. * cw * mv[1]) / g2;
    double vu = (vev * RB[0][0]).real(); // vev* sin(beta)
    yq[5] = mq[5] * sqrt(2.) / vu;
#endif
  } catch (const std::exception &exc) {
    std::cerr << "warning: couplings might not follow slha parameter input " << std::endl
              << exc.what() << std::endl;
    //G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
    a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
    G_F = 1.166379e-5;
    xw = 1.0 - pow(mv[2] / mv[1], 2);
    sw = std::sqrt(xw);
    cw = sqrt(1.0 - xw);
    g2 = 2.0 * mv[2] * sqrt(M_SQRT2 * G_F);
    g3 = sqrt(4 * M_PI * a_s);
    g1 = sw / cw * g2; // not used
  }

  // Stau mixing
  try {
    complex<double> STA[2][2];
    for (int i0 = 0; i0 < 2; i0++) {
      for (int i1 = 0; i1 < 2; i1++) {
        STA[i0][i1] = to<double>(input.at("STAUMIX").at(i0 + 1, i1 + 1).at(2));
      }
    }
    SSLM[0][0] = STA[0][0];
    SSLM[1][0] = STA[1][0];
    SSLM[0][1] = STA[0][1];
    SSLM[1][1] = STA[1][1];
  } catch (...) {
    fprintf(stderr, "warning: No stau mixing defined.\n");
  }

  // Yukawa Yd and Yu
  try {
    yq[0] = to<double>(input.at("YD").at(1, 1).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Yd defined. Set to zero.\n");
    yq[0] = 0.;
  }
  try {
    yq[1] = to<double>(input.at("YD").at(2, 2).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Ys defined. Set to zero.\n");
    yq[1] = 0.;
  }
  try {
    yq[2] = to<double>(input.at("YD").at(3, 3).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Yb defined. Set to 0*5.86008380E-01.\n");
    yq[2] = 0. * 5.86008380E-01;
  }
  try {
    yq[3] = to<double>(input.at("YU").at(1, 1).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Yu defined. Set to zero.\n");
    yq[3] = 0.;
  }
  try {
    yq[4] = to<double>(input.at("YU").at(2, 2).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Yc defined. Set to zero.\n");
    yq[4] = 0.;
  }
  try {
    yq[5] = to<double>(input.at("YU").at(3, 3).at(2));
  } catch (...) {
    fprintf(stderr, "warning: No Yukawa coupling Yt defined. Set to 8.35181765E-01.\n");
    yq[5] = 8.35181765E-01;
  }
#ifdef LIGHT_YUKAWAS_ZERO
  // Set Yukawa couplings to 0 for (d,s,b,u,c).
  for (int i0 = 0; i0 < 5; i0++) {
    yq[i0] = 0.0;
  }
#endif

// Prospino conventions.
// withpout stops and sbottom.
#ifdef DEGENERATE_SQUARK_MASS
  double aSQm = 0.125 * (mSQ[0] + mSQ[1] + mSQ[3] + mSQ[4] + mSQ[6] + mSQ[7] + mSQ[9] + mSQ[10]);

  // only without stops.
  // double aSQm = 0.1 * (mSQ[0] + mSQ[1] + mSQ[2] + mSQ[3] +
  //                          mSQ[4] + mSQ[5] + mSQ[6] + mSQ[7]
  //                          + mSQ[9] + mSQ[10]);

  for (int i0 = 0; i0 < 6; i0++) {
    //    if (i0 % 3 != 2) {
    mSQ[i0] = aSQm;
    mSQ[i0 + 6] = aSQm;
    mSQs[i0] = pow2(aSQm);
    mSQs[i0 + 6] = pow2(aSQm);
  }
  for (int i0 = 0; i0 < 6; i0++) {
    yq[i0] = 0.0;
  }
#endif
}

void Parameters::read_ZpWp(const char *file) {
  using namespace std;
  using namespace SLHAea;

  ifstream ifs(file);
  if (!ifs.good()) {
    fprintf(stderr, "error: failed to open '%s'.\n", file);
    exit(1);
  }

  Coll input(ifs);

  // Neutrino masses
  ml[0] = 0.0;
  ml[1] = 0.0;
  ml[2] = 0.0;

  mv[0] = 0.0;                                            // photon mass
  mv[1] = to<double>(input.at("SMINPUTS").at("4").at(1)); // Z mass
  mv[2] = to<double>(input.at("MASS").at("24").at(1));    // W mass

  mv[3] = to<double>(input.at("MASS").at("32").at(1)); // Z' mass
  mv[4] = to<double>(input.at("MASS").at("34").at(1)); // W' mass

  a_em = pow(to<double>(input.at("SMINPUTS").at("1").at(1)), -1);
  //   G_F = to<double>(input.at("SMINPUTS").at("2").at(1));
  //   G_F = 1.16639e-5;
  a_s = to<double>(input.at("SMINPUTS").at("3").at(1));
  xw = 1.0 - pow(mv[2] / mv[1], 2);
  sw = sqrt(xw);
  cw = sqrt(1.0 - xw);
  g2 = sqrt(4 * M_PI * a_em) / sw;
  g3 = sqrt(4 * M_PI * a_s);
  g1 = sqrt(4 * M_PI * a_em) / cw;
  fprintf(stderr, "Input parameters and gauge couplings read from Z'/W' input file.\n");

#ifndef SSM
  Gv[3] = to<double>(input.at("32").at("DECAY").at(2)); // Z' width
#endif
  Gv[4] = to<double>(input.at("34").at("DECAY").at(2)); // W' width

#ifndef SSM
  // Coupling of Z' to d-type quarks.
  for (Block::const_iterator line = input.at("ZpddL").begin(); line != input.at("ZpddL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[3][to<int>(line->at(0)) - 1][to<int>(line->at(1)) - 1].L = to<double>(line->at(2));
  }
  for (Block::const_iterator line = input.at("ZpddR").begin(); line != input.at("ZpddR").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[3][to<int>(line->at(0)) - 1][to<int>(line->at(1)) - 1].R = to<double>(line->at(2));
  }

  // Coupling of Z' to u-type quarks.
  for (Block::const_iterator line = input.at("ZpuuL").begin(); line != input.at("ZpuuL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[3][to<int>(line->at(0)) + 2][to<int>(line->at(1)) + 2].L = to<double>(line->at(2));
  }
  for (Block::const_iterator line = input.at("ZpuuR").begin(); line != input.at("ZpuuR").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[3][to<int>(line->at(0)) + 2][to<int>(line->at(1)) + 2].R = to<double>(line->at(2));
  }

  // Coupling of Z' to neutrinos.
  for (Block::const_iterator line = input.at("ZpvvL").begin(); line != input.at("ZpvvL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vll[3][to<int>(line->at(0)) - 1][to<int>(line->at(1)) - 1].L = to<double>(line->at(2));
  }

  // Coupling of Z' to charged leptons.
  for (Block::const_iterator line = input.at("ZpllL").begin(); line != input.at("ZpllL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vll[3][to<int>(line->at(0)) + 2][to<int>(line->at(1)) + 2].L = to<double>(line->at(2));
  }
  for (Block::const_iterator line = input.at("ZpllR").begin(); line != input.at("ZpllR").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vll[3][to<int>(line->at(0)) + 2][to<int>(line->at(1)) + 2].R = to<double>(line->at(2));
  }

  // Coupling of W' to quarks.
  for (Block::const_iterator line = input.at("WpudL").begin(); line != input.at("WpudL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[4][to<int>(line->at(0)) + 2][to<int>(line->at(1)) - 1].L = to<double>(line->at(2));
    vqq[4][to<int>(line->at(0)) - 1][to<int>(line->at(1)) + 2].L = to<double>(line->at(2));
  }
  for (Block::const_iterator line = input.at("IMWpudL").begin(); line != input.at("IMWpudL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vqq[4][to<int>(line->at(0)) + 2][to<int>(line->at(1)) - 1].L += II * to<double>(line->at(2));
    vqq[4][to<int>(line->at(0)) - 1][to<int>(line->at(1)) + 2].L -= II * to<double>(line->at(2));
  }

  // Coupling of W' to leptons.
  for (Block::const_iterator line = input.at("WpvlL").begin(); line != input.at("WpvlL").end();
       ++line) {
    if (!line->is_data_line()) {
      continue;
    }
    vll[4][to<int>(line->at(0)) + 2][to<int>(line->at(1)) - 1].L = to<double>(line->at(2));
    vll[4][to<int>(line->at(0)) - 1][to<int>(line->at(1)) + 2].L = to<double>(line->at(2));
  }
#endif
}

// prints the most important parameters in parameter file
// use with e.g. --parameter-file=param.out
// the output will be in json-format
void Parameters::write_log(const char *file) {
  ofstream fout;
  fout.open(file);
  fout.precision(5);

  fout << "{"
       << "'Couplings':{\n"
       << "'g1':" << g1 << ",\n"
       << "'g2':" << g2 << ",\n"
       << "'g3':" << g3 << ",\n"
       << "'xw':" << xw << ",\n"
       << "'sw':" << sw << ",\n"
       << "'cw':" << cw << ",\n"
       << "'Yd':" << yq[0] << ",\n"
       << "'Ys':" << yq[1] << ",\n"
       << "'Yb':" << yq[2] << ",\n"
       << "'Yu':" << yq[3] << ",\n"
       << "'Yc':" << yq[4] << ",\n"
       << "'Yt':" << yq[5] << ",\n"
       << "},\n";

  fout << "\n"
       << "'Higgs sector':{\n"
       << "'mh':";
  fout << "[";
  for (int i0 = 0; i0 < 6; i0++) {
    if (i0 != 0) {
      fout << ",";
    }
    fout << mh[i0];
  }
  fout << "],\n";
  fout << "'Alpha:'[";
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      if (i0 != 0 & i1 != 0) {
        fout << ",";
      }
      fout << RA[i0][i1];
    }
    fout << "]\n";
  }
  fout << "Beta  =\n";
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      fout << "\t" << RB[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Neutralino sector\n"
       << "----------------------\n"
       << "mN =\n";
  for (int i0 = 0; i0 < 4; i0++) {
    fout << "\t" << mCH[i0];
  }
  fout << "\nN  =\n";
  for (int i0 = 0; i0 < 4; i0++) {
    for (int i1 = 0; i1 < 4; i1++) {
      fout << "\t" << RN[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Chargino sector\n"
       << "----------------------\n"
       << "mC =\n";
  for (int i0 = 0; i0 < 2; i0++) {
    fout << "\t" << mCH[i0 + 4];
  }
  fout << "\nU  =\n";
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      fout << "\t" << RU[i0][i1];
    }
    fout << "\n";
  }
  fout << "V  =\n";
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      fout << "\t" << RV[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Gluino sector\n"
       << "----------------------\n"
       << "mG =\t" << mGL << "\n";

  fout << "----------------------\n"
       << "Down-squark sector\n"
       << "----------------------\n"
       << "mSD =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    fout << "\t" << mSQ[i0];
  }
  fout << "\nRSD =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << RSD[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Up-squark sector\n"
       << "----------------------\n"
       << "mSU =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    fout << "\t" << mSQ[i0 + 6];
  }
  fout << "\nRSU =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << RSU[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Sneutrino sector\n"
       << "----------------------\n"
       << "mSN =\n";
  for (int i0 = 0; i0 < 3; i0++) {
    fout << "\t" << mSL[i0];
  }
  fout << "\n";

  fout << "----------------------\n"
       << "Charged Slepton sector\n"
       << "----------------------\n"
       << "mCSL =\n";
  for (int i0 = 3; i0 < 9; i0++) {
    fout << "\t" << mSL[i0];
  }
  fout << "\nSSLM =\n";
  for (int i0 = 0; i0 < 2; i0++) {
    for (int i1 = 0; i1 < 2; i1++) {
      fout << "\t" << SSLM[i0][i1];
    }
    fout << "\n";
  }

  fout << "----------------------\n"
       << "Z'/W' sector\n"
       << "----------------------\n"
       << "mv =\n";
  for (int i0 = 0; i0 < 5; i0++) {
    fout << "\t" << mv[i0];
  }
  fout << "\nGv =\n";
  for (int i0 = 0; i0 < 5; i0++) {
    fout << "\t" << Gv[i0];
  }
  fout << "\nZ'qqR =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vqq[3][i0][i1].R;
    }
    fout << "\n";
  }
  fout << "\nZ'qqL =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vqq[3][i0][i1].L;
    }
    fout << "\n";
  }
  fout << "\nW'qqR =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vqq[4][i0][i1].R;
    }
    fout << "\n";
  }
  fout << "\nW'qqL =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vqq[4][i0][i1].L;
    }
    fout << "\n";
  }
  fout << "\nZ'llR =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vll[3][i0][i1].R;
    }
    fout << "\n";
  }
  fout << "\nZ'llL =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vll[3][i0][i1].L;
    }
    fout << "\n";
  }
  fout << "\nW'llR =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vll[4][i0][i1].R;
    }
    fout << "\n";
  }
  fout << "\nW'llL =\n";
  for (int i0 = 0; i0 < 6; i0++) {
    for (int i1 = 0; i1 < 6; i1++) {
      fout << "\t" << vll[4][i0][i1].L;
    }
    fout << "\n";
  }
  fout << "\nCKM =\n";
  for (int i0 = 0; i0 < 3; i0++) {
    for (int i1 = 0; i1 < 3; i1++) {
      fout << "\t" << ckm[i0][i1];
    }
    fout << "\n";
  }

  fout.close();
}
