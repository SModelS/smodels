// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section for real corrections of Drell-Yan process.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "kinematics.h"
#include "params.h"
#include "utils.h"

double real_gluon_leptons(const double S, const double M2, const double PT2,
                          const double TH, const double PH, const int YS,
                          Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Sets different kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);

  // Partonic cross section for leptons + gluon.
  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MGss();
    }
  }
  delete ff;

  return born;
}

double real_quark_leptons(const double S, const double M2, const double PT2,
                          const double TH, const double PH, const int YS,
                          Parameters *params) {
  double born = 0.0;

  const int aa = params->in1;
  const int bb = params->in2;
  const int ii = params->out1;
  const int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);

  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MQss();
    }
  }

  delete ff;

  return born;
}

double real_quarkb_leptons(const double S, const double M2, const double PT2,
                           const double TH, const double PH, const int YS,
                           Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // sets kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);

  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MQBss();
    }
  }

  delete ff;

  return born;
}

double DipGA_leptons(const double S, const double M2, const double PT2,
                     const double TH, const double PH, const int YS,
                     Parameters *params) {
  double born = 0.0;
  double x = 0.0;
  double fact = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Set kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);
  ff->SetDipKinematicA(x, fact);

  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MBss();
    }
  }

  delete ff;

  return fact * (2.0 / (1.0 - x) - 1 - x) * born;
}

double DipGB_leptons(const double S, const double M2, const double PT2,
                     const double TH, const double PH, const int YS,
                     Parameters *params) {
  double born = 0.0;
  double x = 0.0;
  double fact = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Sets kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);
  ff->SetDipKinematicB(x, fact);

  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MBss();
    }
  }

  delete ff;

  return fact * (2.0 / (1.0 - x) - 1 - x) * born;
}

double DipQA_leptons(const double S, const double M2, const double PT2,
                     const double TH, const double PH, const int YS,
                     Parameters *params) {
  double born = 0.0;
  double x = 0.0;
  double fact = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Sets kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);
  ff->SetDipKinematicA(x, fact);

  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MBss();
    }
  }

  delete ff;

  return fact * (1.0 - 2.0 * x * (1 - x)) * born;
}

double DipQB_leptons(const double S, const double M2, const double PT2,
                     const double TH, const double PH, const int YS,
                     Parameters *params) {
  double born = 0.0;
  double x = 0.0;
  double fact = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Set kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH,
                   YS);
  ff->SetDipKinematicB(x, fact);

  // channels[0] = {y, Z0, Z'}, channels[1] = {W, W'}
  static const int channels[][4] = {{0, 1, 3, -1}, {2, 4, -1, -1}};

  for (int i0 = 0; channels[qs][i0] > -1; i0++) {
    for (int i1 = 0; channels[qs][i1] > -1; i1++) {
      struct Coupling Cw[4] = {params->vll[channels[qs][i0]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i0]][bb][aa],
                               params->vll[channels[qs][i1]][ii - 20][jj - 20],
                               params->vqq[channels[qs][i1]][bb][aa]};

      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(
          params->mv[channels[qs][i0]], params->mv[channels[qs][i1]],
          params->Gv[channels[qs][i0]], params->Gv[channels[qs][i1]]);
      ff->SetWCoupling(Cw);
      born += ff->MBss();
    }
  }
  delete ff;

  return fact * (1.0 - 2.0 * x * (1.0 - x)) * born;
}
