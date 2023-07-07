// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2015 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>

#include "kinematics.h"
#include "params.h"
#include "utils.h"

#define TT
#define UU
#define TU

// Partonic cross section for gaugino-gluino production at LO.
double born_gagl(const double S, const double T, Parameters *params) {

  double born = 0.0;

  // ingoing and outgoing particle id.
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

  // Set born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

#ifdef TT
  if (qt >= 0) {
    for (int index0 = 0, i0 = squark_type[qt][index0]; index0 < 6;
         index0++, i0 = squark_type[qt][index0]) {
      for (int index1 = 0, i1 = squark_type[qt][index1]; index1 < 6;
           index1++, i1 = squark_type[qt][index1]) {
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
  }
#endif

#ifdef UU
  if (qu <= 1) {
    qu = iabs(qu);
    for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
         index0++, i0 = squark_type[qu][index0]) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           index1++, i1 = squark_type[qu][index1]) {

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
        born += ff->Muu_GLGA();
      }
    }
  }
#endif

#ifdef TU
  if (qt >= 0 && qu <= 1) {
    qu = iabs(qu);
    for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
         index0++, i0 = squark_type[qt][index0]) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           index1++, i1 = squark_type[qu][index1]) {
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
  }
#endif
  delete ff;

  return born;
}

// to access sab,sa1 and sb1 somewhere else easily. (for instance in hxs.cc)
// don't get confused: sja := 2 pj.pa (see Catani Seymour massive dipole paper
// p.41)

// sa1
double sja_gagl(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj >= 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  double sja = 2.0 * ff->pap1;
  delete ff;
  return sja;
}

// sb1
double sjb_gagl(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj >= 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  double sjb = 2.0 * ff->pbp1;
  delete ff;
  return sjb;
}

// sb2
double sab_gagl(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj == 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }
  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  double sab = 2.0 * ff->papb;
  delete ff;
  return sab;
}
