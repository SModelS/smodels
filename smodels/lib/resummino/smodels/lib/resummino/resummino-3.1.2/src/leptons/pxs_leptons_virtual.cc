// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Virtual corrections to the partonic cross section of Drell-Yan process.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "kinematics.h"
#include "params.h"
#include "pxs.h"
#include "utils.h"

#define IEPS 0

double Vqcd_leptons(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3);

  // Sets different Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, T);

  // Gluon-quark-quark coupling for strong corrections.
  struct Coupling Cs[2] = {params->gqq[bb][bb], params->gqq[aa][aa]};
  if (is_coupling_null(Cs, 2) == 1) {
    return virt;
  }
  ff->SetSCoupling(Cs);
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

      virt += ff->MVstr1s(0.0, 0.0, 2.0 * ff->papb, 0.0, 0.0, 0.0, IEPS);
    }
  }

  delete ff;

  return virt;
}

double Virt_leptons(const double S, const double T, Parameters *params) {
  const double g3s = std::norm(params->gqq[0][0].R);

  double result = 0.0;
  result += Vqcd_leptons(S, T, params);

  return result / g3s;
}
// Actually not needed!
double Virt2_leptons(const double S, const double T, Parameters *params) {
  const double g3s = std::norm(params->gqq[0][0].R);

  double result = 0.0;
  result += Vqcd_leptons(S, T, params);

  return result / g3s;
}

// Integrated dipole.
double DipI_leptons(const double S, const double T, Parameters *params) {
  const double lnmusq = std::log(params->murs / S);
  if (IEPS == 2) {
    return 2.0 * born_leptons(S, T, params);
  } else if (IEPS == 1) {
    return (3.0 + 2.0 * lnmusq) * born_leptons(S, T, params) +
           2.0 * born_leptons_eps1(S, T, params);
  } else if (IEPS == 0) {
    // Mod: jonathan shifted pieces!
    // return ((3.0 + lnmusq) * lnmusq) * born_leptons(S, T, params)
    //     + (3 + 2.0 * lnmusq) * born_leptons_eps1(S, T, params)
    //     + 2.0 * born_leptons_eps2(S, T, params);
    return (10.0 - pow(M_PI, 2) + (3.0 + lnmusq) * lnmusq) *
               born_leptons(S, T, params) +
           (3 + 2.0 * lnmusq) * born_leptons_eps1(S, T, params) +
           2.0 * born_leptons_eps2(S, T, params);
  } else {
    return 0.0;
  }
}
