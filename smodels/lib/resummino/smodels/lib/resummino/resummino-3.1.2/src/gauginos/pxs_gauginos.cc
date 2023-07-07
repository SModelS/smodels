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
#include <iostream>
#include <stdio.h>

#include "kinematics.h"
#include "params.h"
#include "utils.h"

#define SS
#define TT
#define UU
#define ST
#define SU
#define TU

// Partonic cross section for gauginos.
double born_gauginos(const double S, const double T, Parameters *params) {
  double born = 0.0;
  // initial (aa,bb) and final (ii,jj) state particles
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
  int qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  // abs() because in the propagator is either a W- or a W+
  int qs = iabs(propagator_charge(aa, bb, ii, jj, CHANNEL_S));

  // Sets Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

// sum over different intermediate states
#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                               params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MBss();
    }
  }
#endif

#ifdef TT
  // qt = -1 not possible (would need 4/3 charge)
  if (qt >= 0) {
    for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
         index0++, i0 = squark_type[qt][index0]) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           index1++, i1 = squark_type[qt][index1]) {

        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += ff->MBtt();
      }
    }
  }
#endif

#ifdef UU
  // qu = 1 not possible (would need 4/3 charge)
  if (qu <= 1) {
    // qu -1 are same squarks as +1; since it is used as an index use abs()
    qu = iabs(qu);
    for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
         index0++, i0 = squark_type[qu][index0]) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           index1++, i1 = squark_type[qu][index1]) {

        struct Coupling Cw[4] = {params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += ff->MBuu();
      }
    }
  }
#endif

#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           index1++, i1 = squark_type[qt][index1]) {

        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MBst();
      }
    }
  }
#endif

#ifdef SU
  if (qu <= 1) {
    qu = iabs(qu);
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           index1++, i1 = squark_type[qu][index1]) {

        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MBsu();
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

        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MBtu();
      }
    }
  }
#endif
  delete ff;
  return born;
}

// pxs-part which is prop. to eps (needed for dipole subtraction)
double born_gauginos_eps1(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qt = propagator_charge(aa, bb, ii, jj, CHANNEL_T);
  int qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  // abs() because in the propagator is either a W- or a W+
  int qs = iabs(propagator_charge(aa, bb, ii, jj, CHANNEL_S));

  // Sets different Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                               params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MB1ss();
    }
  }
#endif
#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           index1++, i1 = squark_type[qt][index1]) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MB1st();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    // qu -1 are same squarks as +1; since it is used as an index use abs()
    qu = iabs(qu);
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           index1++, i1 = squark_type[qu][index1]) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MB1su();
      }
    }
  }
#endif

  delete ff;

  return born;
}

// pxs-part which is prop. to eps (needed for dipole subtraction)
double born_gauginos_eps2(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(propagator_charge(aa, bb, ii, jj, CHANNEL_S));

  // Sets different Born kinematics.
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                               params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MB2ss();
    }
  }
#endif

  delete ff;
  return born;
}
