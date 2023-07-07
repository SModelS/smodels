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

#include "debug.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"

#include "pxs_gausq.h"

#define SS
#define UU
#define SU

double born_gasq_apn(const double S, const double T, Parameters *params) {
  FI *ff = new FI();
  {
    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;
    BORN_KINEMATICS;
  }

  double Q2 = 2 * ff->papb;
  double P1K1 = ff->pap1;
  //std::cout << S << " " << T << " " << P1K1 << std::endl;

  BOX_KINEMATIC;
  double g3s = std::norm(params->gqq[0][0].R);
  delete ff;
  return real(ME_us_born(1, 1, 1, 1, 0, Q2, P1K1, params)) * 96 / g3s;
}
double born_gasq(const double S, const double T, Parameters *params) {

  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = propagator_charge(aa, bb, ii, jj, CHANNEL_S);
  int qu = propagator_charge(aa, bb, ii, jj, CHANNEL_U);
  FI *ff = new FI();

  BORN_KINEMATICS;

  if (ii > 30) {
    ii = ii - 31;
    //_DEBUG_MSG1("ii=%d", ii)
  } else {
    jj = jj - 31;
    //_DEBUG_MSG1("jj=%d", jj)
  }

  struct Coupling Cb[2] = {0, 0};

#ifdef SS
  ff->SetPropagator(0, 0, 0, 0);
  // p2 is gaugino -> Squark production
  if (params->out2 <= 10) {
    Cb[0] = params->CHSQq[jj][ii][aa];
    Cb[1] = params->CHSQq[jj][ii][aa];
    ff->SetBCoupling(Cb);
    born += ff->Mss_SQGA1();

    // p1 is gaugino -> Antisquark production
  } else if (params->out1 <= 10) {
    Cb[0] = params->CHqSQ[ii][bb][jj];
    Cb[1] = params->CHqSQ[ii][bb][jj];
    ff->SetBCoupling(Cb);
    born += ff->Mss_SQGA2();
  }

#endif

#ifdef UU

  // p2 is gaugino
  if (params->out2 <= 10) {
    Cb[0] = params->CHSQq[jj][ii][aa];
    Cb[1] = params->CHSQq[jj][ii][aa];

    ff->SetPropagator(params->mSQ[ii], params->mSQ[ii], 0.0, 0.0);

    ff->SetBCoupling(Cb);
    born += ff->Muu_SQGA1();

    // p1 is gaugino
  } else if (params->out1 <= 10) {
    Cb[0] = params->CHqSQ[ii][bb][jj];
    Cb[1] = params->CHqSQ[ii][bb][jj];

    ff->SetPropagator(params->mSQ[jj], params->mSQ[jj], 0.0, 0.0);

    ff->SetBCoupling(Cb);
    born += ff->Muu_SQGA2();
  }
#endif

#ifdef SU

  // p2 is gaugino
  if (params->out2 <= 10) {
    Cb[0] = params->CHSQq[jj][ii][aa];
    Cb[1] = params->CHSQq[jj][ii][aa];
    ff->SetPropagator(0.0, params->mSQ[ii], 0.0, 0 * params->mSQ[ii] * 1.0e-2);
    ff->SetBCoupling(Cb);
    born += 2 * ff->Msu_SQGA1();

    // p1 is gaugino
  } else if (params->out1 <= 10) {
    Cb[0] = params->CHqSQ[ii][bb][jj];
    Cb[1] = params->CHqSQ[ii][bb][jj];
    ff->SetPropagator(0.0, params->mSQ[jj], 0.0, 0 * params->mSQ[jj] * 1.0e-2);
    ff->SetBCoupling(Cb);
    born += 2 * ff->Msu_SQGA2();
  }

#endif

  delete ff;
  return born;
}

// sa1
double sja_gasq(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj > 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  FI *ff = new FI();
  ff->SetKinematic(params->mSQ[ii - 31], params->mCH[jj], S, T);

  double sja = 2.0 * ff->pap1;
  delete ff;
  return sja;
}

// sb1
double sjb_gasq(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj > 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }

  FI *ff = new FI();
  ff->SetKinematic(params->mSQ[ii - 31], params->mCH[jj], S, T);

  double sjb = 2.0 * ff->pbp1;
  delete ff;
  return sjb;
}

// sb2
double sab_gasq(const double S, const double T, Parameters *params) {

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // swap
  if (jj > 30) {
    int ii_temp = ii;
    ii = jj;
    jj = ii_temp;
  }
  FI *ff = new FI();
  ff->SetKinematic(params->mSQ[ii - 31], params->mCH[jj], S, T);

  double sab = 2.0 * ff->papb;
  delete ff;
  return sab;
}
