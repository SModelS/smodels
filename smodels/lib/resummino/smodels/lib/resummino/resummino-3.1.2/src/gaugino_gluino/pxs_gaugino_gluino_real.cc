// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>

#include "dipoles.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"
// constants.h must be last include
#include "constants.h"

// Different channels.
#define TT
#define UU
#define TU

// Color factors.
#define cA 3.0
#define cF (4.0 / 3.0)
#define TR 0.5

// On-shell subtraction counterterm.
// (Notes: If you notice a discrepancy between resummino
// and Prospino it is mainly due to the onshell remainder;
// try to compare both codes with gluon PDFs turned off!)
#define ONSUB

// Factor for minimal width approximation: mSQ * WIDTH.
// The smaller the width the less off-shell contributions get subtracted.
// However, too small values lead to an unstable numerical integration.
// (increase number of calls if your on-shell contribution is sizable)

//(changed from 1.0E-3 to 1.0E-2 to get more stable results.
#define WIDTH 1.0E-2

// only regularize resonant propagator if needed -> results are less width
// dependent!
// (seems like Prospino is doing sth like this! see Xmatrix_ng_c.f90 lines 68ff
// in prospino code)
// (get similar results as with smaller WIDTH 1.0E-3, but more stable)
// this is only used in the squared resonant diagrams!
// the interference terms between resonant and non-resonant diagrams lead to
// minor contributions
#define PROSPINO

// Partonic cross section for real gluon emission.
double real_gluon_gaugino_gluino(const double S, const double M2, const double PT2, const double TH,
                                 const double PH, const int YS, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  // Set kinematics for 2 -> 3 processes.
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, M2, PT2, TH, PH, YS);

#ifdef TT
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);

      born += ff->MttG_GLGA();
    }
  }
#endif

#ifdef UU
  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);

      born += ff->MuuG_GLGA();
    }
  }
#endif

#ifdef TU
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);
      born += ff->MtuG_GLGA();
    }
  }
#endif
  delete ff;
  return born;
}

// Partonic cross section for real quark emission without on-shell
// contributions.
double real_quark_gaugino_gluino(const double S, const double M2, const double PT2, const double TH,
                                 const double PH, const int YS, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U<->T
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, M2, PT2, TH, PH, YS);

#ifdef TT
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      born += ff->MttQ_GLGA_wos();
    }
  }
#endif

#ifdef UU
  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      born += ff->MuuQ_GLGA_wos();
    }
  }
#endif

#ifdef TU
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      born += ff->MtuQ_GLGA();
    }
  }
#endif
  delete ff;
  return born;
}

// Partonic cross section for real antiquark emission without on-shell
// contributions.
double real_quarkb_gaugino_gluino(const double S, const double M2, const double PT2,
                                  const double TH, const double PH, const int YS,
                                  Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  // Set different kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, M2, PT2, TH, PH, YS);

#ifdef TT
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      born += ff->MttQB_GLGA_wos();
    }
  }
#endif

#ifdef UU
  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      Cw[2] = params->CHSQq[jj][i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      born += ff->MuuQB_GLGA_wos();
    }
  }
#endif

#ifdef TU
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

      ff->SetWCoupling(Cw);
      born += ff->MtuQB_GLGA();
    }
  }
#endif
  delete ff;
  return born;
}

// Partonic cross section for real quark emission where s23 can lead to an
// on-shell resonance.
double real_quark_gaugino_gluino_onshell_23(const double S, const double S2, const double T1,
                                            const double S1, const double PHI, int flag,
                                            Parameters *params) {
  double born = 0.0;

  // Jacobian factor.
  double factor = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;

  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  double m1 = params->mGL;
  double m2 = params->mCH[jj];

  FI *ff = new FI();
  // Set 2->3 kinematics using helicity frame.
  ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);

#ifdef TT
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

#ifdef PROSPINO
      if ((2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m1) && params->mSQ[i0] >= (ff->m2)) ||
          (2.0 * ff->papb >= pow2(params->mSQ[i1] + ff->m1) && params->mSQ[i1] >= (ff->m2))) {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

      } else {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
      }
#endif

#ifndef PROSPINO
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

#endif

      born += ff->MttQres_GLGA();
#ifdef ONSUB
      if (i0 == i1 && 2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m1) &&
          params->mSQ[i0] >= (ff->m2)) {
        // Set Breit Wigner form.
        ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(params->mSQ[i0]),
                           params->mSQ[i0] * WIDTH);
        // On-shell kinematics.
        ff->SetKinematicRES23(pow2(params->mSQ[i0]), factor);
        // Minimal width needed.
        // Result in principle width independent for very small widths.
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

        // Subtraction counter term.
        born -= ff->BreitWigner * factor * ff->MttQres_GLGA();
        // Reset to ordinary kinematics.
        ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);
      }
#endif
    }
  }

#endif

  delete ff;
  return born;
}

// Partonic cross section for real quark emission where s13 can lead to an
// on-shell resonance.
double real_quark_gaugino_gluino_onshell_13(const double S, const double S2, const double T1,
                                            const double S1, const double PHI, int flag,
                                            Parameters *params) {
  double born = 0.0;
  double factor = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  double m1 = params->mGL;
  double m2 = params->mCH[jj];

  FI *ff = new FI();
  ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);

#ifdef UU

  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetWCoupling(Cw);

#ifdef PROSPINO
      if ((2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m2) && params->mSQs[i0] >= (ff->m1s)) ||
          (2.0 * ff->papb >= pow2(params->mSQ[i1] + ff->m2) && params->mSQs[i1] >= (ff->m1s))) {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

      } else {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
      }
#endif

#ifndef PROSPINO
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);
#endif
      born += ff->MuuQres_GLGA();
#ifdef ONSUB
      if (i0 == i1 && 2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m2) &&
          params->mSQs[i0] >= (ff->m1s)) {

        ff->SetBreitWigner(2.0 * ff->p1p3 + ff->m1s, pow2(params->mSQ[i0]),
                           params->mSQ[i0] * WIDTH);
        ff->SetKinematicRES13(pow2(params->mSQ[i0]), factor);
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

        born -= ff->BreitWigner * factor * ff->MuuQres_GLGA();

        ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);
      }
#endif
    }
  }
#endif

  delete ff;
  return born;
}

// Partonic cross section for real antiquark emission where s13 can lead to an
// on-shell resonance.
double real_quarkb_gaugino_gluino_onshell_13(const double S, const double S2, const double T1,
                                             const double S1, const double PHI, int flag,
                                             Parameters *params) {
  double born = 0.0;
  double factor = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  double m1;
  double m2;

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  m1 = params->mGL;
  m2 = params->mCH[jj];

  FI *ff = new FI();
  ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);

#ifdef TT

  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetWCoupling(Cw);

#ifdef PROSPINO
      if ((2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m2) && params->mSQ[i0] >= (ff->m1)) ||
          (2.0 * ff->papb >= pow2(params->mSQ[i1] + ff->m2) && params->mSQ[i1] >= (ff->m1))) {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

      } else {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
      }
#endif

#ifndef PROSPINO
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

#endif
      born += ff->MttQBres_GLGA();

#ifdef ONSUB
      if (i0 == i1 && 2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m2) &&
          params->mSQ[i0] >= (ff->m1)) {

        ff->SetBreitWigner(2.0 * ff->p1p3 + ff->m1s, pow2(params->mSQ[i0]),
                           params->mSQ[i0] * WIDTH);

        ff->SetKinematicRES13(pow2(params->mSQ[i0]), factor);

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

        born -= ff->BreitWigner * factor * ff->MttQBres_GLGA();

        ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);
      }
#endif
    }
  }
#endif

  delete ff;
  return born;
}

// Partonic cross section for real antiquark emission where s23 can lead to an
// on-shell resonance.
double real_quarkb_gaugino_gluino_onshell_23(const double S, const double S2, const double T1,
                                             const double S1, const double PHI, int flag,
                                             Parameters *params) {
  double born = 0.0;
  double factor = 0.0;
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  double m1;
  double m2;

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  m1 = params->mGL;
  m2 = params->mCH[jj];

  FI *ff = new FI();
  ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);

#ifdef UU

  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetWCoupling(Cw);

#ifdef PROSPINO
      if ((2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m1) && params->mSQs[i0] >= (ff->m2s)) ||
          (2.0 * ff->papb >= pow2(params->mSQ[i1] + ff->m1) && params->mSQs[i1] >= (ff->m2s))) {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

      } else {

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
      }
#endif

#ifndef PROSPINO
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                        params->mSQ[i1] * WIDTH);

#endif
      born += ff->MuuQBres_GLGA();
#ifdef ONSUB

      if (i0 == i1 && 2.0 * ff->papb >= pow2(params->mSQ[i0] + ff->m1) &&
          params->mSQs[i0] >= (ff->m2s)) {

        ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(params->mSQ[i0]),
                           params->mSQ[i0] * WIDTH);
        ff->SetKinematicRES23(pow2(params->mSQ[i0]), factor);
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * WIDTH,
                          params->mSQ[i1] * WIDTH);

        born -= ff->BreitWigner * factor * ff->MuuQBres_GLGA();

        ff->SetKinematic_HelicityFrame(m1, m2, S, S2, T1, S1, PHI, flag);
      }
#endif
    }
  }
#endif

  delete ff;
  return born;
}

// Different dipoles dsigma^A for real emission of gluon, quarks and antiquarks.
double Dip_GLGA(DipoleType Emitter_Spectator, PartonType emitter, PartonType spectator,
                PartonType emitted_parton, const double S, const double M2, const double PT2,
                const double TH, const double PH, const int YS, Parameters *params) {

  double born = 0.0;
  double x = 0.0;
  double zi = 0.0;
  double zj = 0.0;
  double zplus = 0.0;
  double zminus = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt, qu;
  // If chargino is particle (instead of antiparticle) U <-> T.
  if (ii == 30) {
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  } else if (jj == 30) {
    qu = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
    qt = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  }

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  // Set different kinematics for 2 -> 3 processes.
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, M2, PT2, TH, PH, YS);
  double Dipole;

  switch (Emitter_Spectator) {
  case INITIAL_INITIAL:
    if (emitter == PARTON_QUARK && spectator == PARTON_ANTIQUARK &&
        emitted_parton == PARTON_GLUON) { // D_qg_q

      // Set dipole kinematics for emitter A and spectator B.
      ff->SetDipKinematicAB(x);
      // Dipole factor.
      Dipole = D_ai_b(PARTON_QUARK, PARTON_GLUON, ff->pap3, ff->papb, ff->pbp3, x, 1,
                      0.5 * (cA - 2.0 * cF));

    } else if (emitter == PARTON_ANTIQUARK && spectator == PARTON_QUARK &&
               emitted_parton == PARTON_GLUON) { // D_qg_q

      ff->SetDipKinematicBA(x);
      Dipole = D_ai_b(PARTON_ANTIQUARK, PARTON_GLUON, ff->pbp3, ff->papb, ff->pap3, x, 1,
                      0.5 * (cA - 2.0 * cF));

    } else if (emitter == PARTON_GLUON && spectator == PARTON_QUARK &&
               emitted_parton == PARTON_QUARK) { // D_gq_q

      ff->SetDipKinematicBA(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_ANTIQUARK, ff->pbp3, ff->papb, ff->pap3, x, 1,
                      0.5 * (cA - 2.0 * cF) / cF * TR);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_ANTIQUARK &&
               emitted_parton == PARTON_ANTIQUARK) { // D_gq_q

      ff->SetDipKinematicAB(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_QUARK, ff->pap3, ff->papb, ff->pbp3, x, 1,
                      0.5 * (cA - 2.0 * cF) / cF * TR);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_GLUON &&
               emitted_parton == PARTON_ANTIQUARK) { //  D_gq_g

      ff->SetDipKinematicAB(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_ANTIQUARK, ff->pap3, ff->papb, ff->pbp3, x, 1, 0);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_GLUON &&
               emitted_parton == PARTON_QUARK) { //  D_gq_g

      ff->SetDipKinematicBA(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_QUARK, ff->pbp3, ff->papb, ff->pap3, x, 1, 0);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_QUARK &&
               emitted_parton == PARTON_GLUON) { //  D_gg_q

      ff->SetDipKinematicBA(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_GLUON, ff->pbp3, ff->papb, ff->pap3, x, 1, 0);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_ANTIQUARK &&
               emitted_parton == PARTON_GLUON) { //  D_gg_q

      ff->SetDipKinematicAB(x);
      Dipole = D_ai_b(PARTON_GLUON, PARTON_GLUON, ff->pap3, ff->papb, ff->pbp3, x, 1, 0);
    }

  case INITIAL_FINAL:
    if (emitter == PARTON_QUARK && spectator == PARTON_GLUINO &&
        emitted_parton == PARTON_GLUON) { // D_qg_gl

      ff->SetDipKinematicA1(x, zi, zj);
      Dipole = D_ai_j(PARTON_QUARK, PARTON_GLUON, zi, zj, ff->pap3, ff->p1p3, x, 1, -0.5 * cA);

    } else if (emitter == PARTON_ANTIQUARK && spectator == PARTON_GLUINO &&
               emitted_parton == PARTON_GLUON) { // D_qg_gl

      ff->SetDipKinematicB1(x, zi, zj);
      Dipole = D_ai_j(PARTON_ANTIQUARK, PARTON_GLUON, zi, zj, ff->pbp3, ff->p1p3, x, 1, -0.5 * cA);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_GLUINO &&
               emitted_parton == PARTON_ANTIQUARK) { // D_gq_gl

      ff->SetDipKinematicA1(x, zi, zj);
      Dipole = D_ai_j(PARTON_GLUON, PARTON_ANTIQUARK, zi, zj, ff->pap3, ff->p1p3, x, 1,
                      -0.5 * cA / cF * TR);

    } else if (emitter == PARTON_GLUON && spectator == PARTON_GLUINO &&
               emitted_parton == PARTON_QUARK) { // D_gq_gl

      ff->SetDipKinematicB1(x, zi, zj);
      Dipole =
          D_ai_j(PARTON_GLUON, PARTON_QUARK, zi, zj, ff->pbp3, ff->p1p3, x, 1, -0.5 * cA / cF * TR);
    }
  case FINAL_INITIAL:
    if (emitter == PARTON_GLUINO && spectator == PARTON_QUARK && emitted_parton == PARTON_GLUON) {

      ff->SetDipKinematic1A(x, zi, zj, zplus, zminus, ff->m1, 0.0, ff->m1);
      Dipole = D_ij_a(PARTON_GLUON, PARTON_GLUINO, 0.0, params->mGL, params->mGL, ff->p1p3, x, 1,
                      zj, params->mGL, zi, zplus, zminus, -0.5 * cA);

    } else if (emitter == PARTON_GLUINO && spectator == PARTON_ANTIQUARK &&
               emitted_parton == PARTON_GLUON) {
      ff->SetDipKinematic1B(x, zi, zj, zplus, zminus, ff->m1, 0.0, ff->m1);
      Dipole = D_ij_a(PARTON_GLUON, PARTON_GLUINO, 0.0, params->mGL, params->mGL, ff->p1p3, x, 1,
                      zj, params->mGL, zi, zplus, zminus, -0.5 * cA);
    }
  }

// Partonic born cross section with dipole kinematics.
#ifdef TT
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6; i1 = squark_type[qt][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHqSQ[jj][bb][i1];
      Cw[3] = params->GLSQq[i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);
      born += ff->Mtt_GLGA();
    }
  }
#endif

#ifdef UU
  for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6; i0 = squark_type[qu][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[1] = params->GLqSQ[bb][i0];
      Cw[0] = params->CHSQq[jj][i0][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      Cw[2] = params->CHSQq[jj][i1][aa];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);
      born += ff->Muu_GLGA();
    }
  }
#endif

#ifdef TU
  for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6; i0 = squark_type[qt][++index0]) {
    for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6; i1 = squark_type[qu][++index1]) {
      struct Coupling Cw[4] = {0, 0, 0, 0};
      Cw[0] = params->CHqSQ[jj][bb][i0];
      Cw[1] = params->GLSQq[i0][aa];
      Cw[2] = params->CHSQq[jj][i1][aa];
      Cw[3] = params->GLqSQ[bb][i1];
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);
      born += 2 * ff->Mtu_GLGA();
    }
  }
#endif

  delete ff;

  return 1. / (4 * M_PI) * Dipole * born;
}
