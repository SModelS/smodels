// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section for Born cross section of Drell-Yan process.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "kinematics.h"
#include "params.h"
#include "utils.h"

double born_leptons(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = qs = iabs(aa / 3 - bb / 3);

  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, T);

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

  return born;
}

double born_leptons_eps1(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, T);

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
      born += ff->MB1ss();
    }
  }

  delete ff;

  return born;
}

double born_leptons_eps2(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.

  // Sets different Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, T);

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
      born += ff->MB2ss();
    }
  }

  delete ff;

  return born;
}
