// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Partonic cross section for slepton pair production

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "kinematics.h"
#include "params.h"
#include "utils.h"


// Partonic cross section for slepton pair production.
double born_sleptons(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1 - 10;
  int jj = params->out2 - 10;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel

  // Sets different Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->mSL[ii], params->mSL[jj], S, T);

  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vSLSL[i0][ii][jj], params->vqq[i0][bb][aa],
          params->vSLSL[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);

      ff->SetWCoupling(Cw);
      born += ff->MBssSL();
      //std::cout << i0 << " " << i1 << " " << born <<  " " << S << " " << T << " "<< ff->ivs3v1 << " " << ff->ivs3v2<< " " <<params->mv[i0] << " " <<  params->mv[i1] << std::endl;
    }
  }
  delete ff;
  return born;
}
