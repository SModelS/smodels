// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "dipoles.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"

// The correct color factor is factored out in hxs.cc!
#define cF 1.0
#define TR 1.0

// Partonic cross section for real gluon emission.
double real_gluon_sleptons(const double S, const double M2, const double PT2,
                           const double TH, const double PH, const int YS,
                           Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, M2, PT2, TH,
                   PH, YS);

  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
          params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MGssSL();
    }
  }
  delete ff;
  return born;
}

// Partonic cross section for real quark emission.
double real_quark_sleptons(const double S, const double M2, const double PT2,
                           const double TH, const double PH, const int YS,
                           Parameters *params) {
  double born = 0.0;

  const int aa = params->in1;
  const int bb = params->in2;
  const int ii = params->out1;
  const int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, M2, PT2, TH,
                   PH, YS);

  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
          params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MQssSL();
    }
  }
  delete ff;
  return born;
}

// Partonic cross section for real antiquark emission.
double real_quarkb_sleptons(const double S, const double M2, const double PT2,
                            const double TH, const double PH, const int YS,
                            Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, M2, PT2, TH,
                   PH, YS);

  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
          params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MQBssSL();
    }
  }
  delete ff;
  return born;
}

// Corresponding dipoles dsigma^A.
double Dip_SLEPTONS(DipoleType Emitter_Spectator, PartonType emitter,
                    PartonType spectator, PartonType emitted_parton,
                    const double S, const double M2, const double PT2,
                    const double TH, const double PH, const int YS,
                    Parameters *params) {

  double born = 0.0;
  double x = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
  int qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  // Sets kinematics for 2 -> 3 process.
  FI *ff = new FI();

  JET_KINEMATICS
  double Dipole;

  if (emitter == PARTON_QUARK && spectator == PARTON_ANTIQUARK &&
      emitted_parton == PARTON_GLUON) { // D_qg_q

    ff->SetDipKinematicAB(x);
    Dipole = D_ai_b(PARTON_QUARK, PARTON_GLUON, ff->pap3, ff->papb, ff->pbp3, x,
                    1, 0.5 * (-2.0 * cF));

  } else if (emitter == PARTON_ANTIQUARK && spectator == PARTON_QUARK &&
             emitted_parton == PARTON_GLUON) { // D_qg_q

    ff->SetDipKinematicBA(x);
    Dipole = D_ai_b(PARTON_ANTIQUARK, PARTON_GLUON, ff->pbp3, ff->papb,
                    ff->pap3, x, 1, 0.5 * (-2.0 * cF));

  } else if (emitter == PARTON_GLUON && spectator == PARTON_QUARK &&
             emitted_parton == PARTON_QUARK) { // D_gq_q

    ff->SetDipKinematicBA(x);
    Dipole = D_ai_b(PARTON_GLUON, PARTON_ANTIQUARK, ff->pbp3, ff->papb,
                    ff->pap3, x, 1, 0.5 * (-2.0 * cF) / cF * TR);

  } else if (emitter == PARTON_GLUON && spectator == PARTON_ANTIQUARK &&
             emitted_parton == PARTON_ANTIQUARK) { // D_gq_q

    ff->SetDipKinematicAB(x);
    Dipole = D_ai_b(PARTON_GLUON, PARTON_QUARK, ff->pap3, ff->papb, ff->pbp3, x,
                    1, 0.5 * (-2.0 * cF) / cF * TR);
  }

  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vSLSL[i0][ii - 10][jj - 10], params->vqq[i0][bb][aa],
          params->vSLSL[i1][ii - 10][jj - 10], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MBssSL();
    }
  }
  delete ff;
  return 1. / (4 * M_PI) * Dipole * born;
}
