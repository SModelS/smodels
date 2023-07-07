// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Computes the virtual part of the partonic NLO cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "kinematics.h"
#include "params.h"
#include "pxs.h"
#include "utils.h"

#include "gsl_all.h"

using namespace std;

// IEPS 0 (finite part), 1 (1/eps coefficient) and 2 (1/eps^2 coefficient).
#define IEPS 0

// Finite shift from DREG to DRED.
#define DRBAR

// Decouple squarks, gluinos and top from the aS running.
// (see arxiv 9610490v1 p. 17 above eq. 19)
#define DECOUPLE_HEAVY

// Color factors.
#define cA 3.0
#define cF (4.0 / 3.0)
#define NC 3.0
#define NA 8.0

// Quark and Gluon "mass" parameters.
// Set to non-zero values to regularize the IR divergences in the loops by a
// small mass.
// Useful to check if UV counterterms kill all UV poles.
#define MGS 0.0
#define MQS 0.0

// Different channels: M_virtT M_bornT = TT, etc.
#define TT
#define UU
#define TU
#define UT

// IR finite contributions.
#define QSELF   // external quark self-enery
#define GLSELFA // external gluino self-energy (squark-quark)
#define SQP1    // internal squark propagator correction 1
#define SQP2    // internal squark propagator correction 2

#define TR2A // squark-gluino-quark loop at leg a
#define TR2B // squark-gluino-quark loop at leg b

#define BO2 // box diagram (gluino-squark-squark-quark)
#define TR4 // gluino-squark-gluon loop at leg 1 (gluino leg)

// 1/eps IR poles.
#define GLSELFB // external gluino self energy (gluino-gluon)
#define TR1     // squark gluon quark loop at leg a and b
#define BU1 // quark self energy (gluon quark; only IR divergent contribution)

// 1/eps^2 & 1/eps IR poles.
#define BO1 // box diagram (quark-gluon-squark-quark)
#define BO3 // box diagram (gluino-gluon-quark-squark)
#define TR3 // gluon-quark-gluino loop at leg 1

// The integrated dipole.
#define DIPOLE

// Includes all the qcd and susy-qcd corrections which are not left- and
// right-handed sensitive.
double qcd(const double S, const double T, Parameters *params) {

  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  const double g3s = norm(params->gqq[0][0].R);

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
        if (is_coupling_null(Cw, 4)) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);

#ifdef TR1
        virt += 2 * g3s *
                ff->Mtt_V1a_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m1s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        virt += 2 * g3s *
                ff->Mtt_V1b_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m2s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mtt_V1a_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m1s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
          virt -= 2 * g3s *
                  ff->Mtt_V1b_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m2s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
        }
#endif

#ifdef SQP2
        virt += 2.0 * g3s *
                ff->Mtt_SQP2_GLGA(ff->m1s - 2.0 * ff->pap1, 0.0,
                                  params->mSQs[i0], params->mSQs[i0],
                                  params->mSQs[i0], IEPS, params->mSQs[i0]);
#endif

#ifdef TR3
        virt += 2 * g3s *
                ff->Mtt_V3a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, MQS,
                                 MGS, params->mGLs, IEPS);

        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mtt_V3a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, 100.0,
                                   100.0, params->mGLs, IEPS);
        }
#endif

#ifdef TR4
        virt += 2 * g3s *
                ff->Mtt_V4a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                 params->mGLs, params->mSQs[i0], MGS, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mtt_V4a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                   params->mGLs, params->mSQs[i0], 100.0, IEPS);
        }
#endif

#ifdef BO1
        virt += 2 * g3s *
                ff->Mtt_B1_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                ff->m1s - 2.0 * ff->pap1, MQS, MGS, MQS,
                                params->mSQs[i0], IEPS);
#endif

#ifdef BO3
        virt +=
            2 * g3s *
            ff->Mtt_B3_GLGA(0.0, ff->m1s, 0.0, ff->m2s,
                            ff->m1s - 2.0 * ff->pbp1, ff->m1s - 2.0 * ff->pap1,
                            MQS, MGS, params->mGLs, params->mSQs[i0], IEPS);
#endif

#ifdef BU1
        if (IEPS == 1) {
          virt += 2.0 * cF * 1.0 / (16.0 * pow2(M_PI)) * g3s * ff->Mtt_GLGA();
        }
#endif

#ifdef GLSELFB
        virt += 2 * g3s *
                ff->Mtt_GLR2_GLGA(params->mGLs, params->mGLs, MGS, params->mGL,
                                  IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          // -aS * cF /(8 pi)
          // is the shift for gaugino-quark-squark vertex and
          // as/(4 pi) * (2.0/3.0 * NC - 1.0/2.0 * cF))
          // is the shift for gluino-quark-squark.
          // (see arxiv:hep-ph/0511344v2 & arXiv:hep-ph/9610490).
          virt += 2.0 * g3s / (32.0 * pow2(M_PI)) * (4.0 / 3.0 * NC - 2 * cF) *
                  ff->Mtt_GLGA();
        }
#endif

#ifdef DECOUPLE_HEAVY
        if (IEPS == 0) {
          double murs = params->murs;
          double mgls = params->mGLs;
          double log_squarks = 0.0;
          for (int i = 0; i < 12; i++) {
            log_squarks += log(params->mSQs[i] / murs);
          }
          double mts = params->mqs[5];

          virt += 2.0 * g3s / (16.0 * pow2(M_PI)) *
                  (-log(mgls / murs) - 1.0 / 12.0 * log_squarks -
                   1.0 / 3.0 * log(mts / murs)) *
                  ff->Mtt_GLGA();
        }
#endif
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
        Cw[1] = params->GLqSQ[bb][i0];
        Cw[0] = params->CHSQq[jj][i0][aa];
        Cw[3] = params->GLqSQ[bb][i1];
        Cw[2] = params->CHSQq[jj][i1][aa];
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
#ifdef SQP2
        virt += 2.0 * g3s *
                ff->Muu_SQP2_GLGA(ff->m2s - 2.0 * ff->pap2, 0.0,
                                  params->mSQs[i0], params->mSQs[i0],
                                  params->mSQs[i0], IEPS, params->mSQs[i0]);
#endif

#ifdef TR1
        virt += 2 * g3s *
                ff->Muu_V1a_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m2s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        virt += 2 * g3s *
                ff->Muu_V1b_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Muu_V1a_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m2s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
          virt -= 2 * g3s *
                  ff->Muu_V1b_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
        }
#endif

#ifdef TR3
        virt += 2 * g3s *
                ff->Muu_V3b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s, MQS,
                                 MGS, params->mGLs, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Muu_V3b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s, 100.0,
                                   100.0, params->mGLs, IEPS);
        }
#endif

#ifdef TR4
        virt += 2 * g3s *
                ff->Muu_V4b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                 params->mGLs, params->mSQs[i0], MGS, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Muu_V4b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                   params->mGLs, params->mSQs[i0], 100.0, IEPS);
        }
#endif

#ifdef BO1
        virt += 2 * g3s *
                ff->Muu_B1_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                ff->m2s - 2.0 * ff->pap2, MQS, MGS, MQS,
                                params->mSQs[i0], IEPS);
#endif

#ifdef BO3
        virt +=
            2 * g3s *
            ff->Muu_B3_GLGA(0.0, ff->m1s, 0.0, ff->m2s,
                            ff->m1s - 2.0 * ff->pap1, ff->m1s - 2.0 * ff->pbp1,
                            MQS, MGS, params->mGLs, params->mSQs[i0], IEPS);
#endif

#ifdef BU1
        if (IEPS == 1) {
          virt += 2.0 * cF * 1.0 / (16.0 * pow2(M_PI)) * g3s * ff->Muu_GLGA();
        }
#endif

#ifdef GLSELFB
        virt += 2 * g3s *
                ff->Muu_GLR2_GLGA(params->mGLs, params->mGLs, MGS, params->mGL,
                                  IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt += 2.0 * g3s / (32.0 * pow2(M_PI)) * (4.0 / 3.0 * NC - 2 * cF) *
                  ff->Muu_GLGA();
        }
#endif

#ifdef DECOUPLE_HEAVY
        if (IEPS == 0) {
          double murs = params->murs;
          double mgls = params->mGLs;
          double log_squarks = 0.0;
          for (int i = 0; i < 12; i++) {
            log_squarks += log(params->mSQs[i] / murs);
          }
          double mts = params->mqs[5];

          virt += 2.0 * g3s / (16.0 * pow2(M_PI)) *
                  (-log(mgls / murs) - 1.0 / 12.0 * log_squarks -
                   1.0 / 3.0 * log(mts / murs)) *
                  ff->Muu_GLGA();
        }
#endif
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
        Cw[3] = params->GLqSQ[bb][i1];
        Cw[2] = params->CHSQq[jj][i1][aa];
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
#ifdef SQP2
        virt += 2.0 * g3s *
                ff->Mtu_SQP2_GLGA(ff->m1s - 2.0 * ff->pap1, 0.0,
                                  params->mSQs[i0], params->mSQs[i0],
                                  params->mSQs[i0], IEPS, params->mSQs[i0]);
#endif

#ifdef TR1
        virt += 2.0 * g3s *
                ff->Mtu_V1a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        virt += 2.0 * g3s *
                ff->Mtu_V1b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        if (IEPS == 1) {
          virt -= 2.0 * g3s *
                  ff->Mtu_V1a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, 100,
                                   100.0, params->mSQs[i0], IEPS);
          virt -= 2.0 * g3s *
                  ff->Mtu_V1b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s, 100,
                                   100.0, params->mSQs[i0], IEPS);
        }
#endif

#ifdef TR3
        virt += 2 * g3s *
                ff->Mtu_V3a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, MQS,
                                 MGS, params->mGLs, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mtu_V3a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, 100.0,
                                   100.0, params->mGLs, IEPS);
        }
#endif

#ifdef TR4
        virt += 2 * g3s *
                ff->Mtu_V4a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                 params->mGLs, params->mSQs[i0], MGS, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mtu_V4a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                   params->mGLs, params->mSQs[i0], 100.0, IEPS);
        }
#endif

#ifdef BO1
        virt += 2 * g3s *
                ff->Mtu_B1_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                ff->m1s - 2.0 * ff->pap1, MQS, MGS, MQS,
                                params->mSQs[i0], IEPS);
#endif

#ifdef BO3
        virt +=
            2 * g3s *
            ff->Mtu_B3_GLGA(0.0, ff->m1s, 0.0, ff->m2s,
                            ff->m1s - 2.0 * ff->pbp1, ff->m1s - 2.0 * ff->pap1,
                            MQS, MGS, params->mGLs, params->mSQs[i0], IEPS);
#endif

#ifdef BU1
        if (IEPS == 1) {
          virt += 2.0 * cF * 1.0 / (16.0 * pow2(M_PI)) * g3s * ff->Mtu_GLGA();
        }
#endif

#ifdef GLSELFB
        virt += 2 * g3s *
                ff->Mtu_GLR2_GLGA(params->mGLs, params->mGLs, MGS, params->mGL,
                                  IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt += 2.0 * g3s / (32.0 * pow2(M_PI)) * (4.0 / 3.0 * NC - 2 * cF) *
                  ff->Mtu_GLGA();
        }
#endif

#ifdef DECOUPLE_HEAVY
        if (IEPS == 0) {
          double murs = params->murs;
          double mgls = params->mGLs;
          double log_squarks = 0.0;
          for (int i = 0; i < 12; i++) {
            log_squarks += log(params->mSQs[i] / murs);
          }
          double mts = params->mqs[5];

          virt += 2.0 * g3s / (16.0 * pow2(M_PI)) *
                  (-log(mgls / murs) - 1.0 / 12.0 * log_squarks -
                   1.0 / 3.0 * log(mts / murs)) *
                  ff->Mtu_GLGA();
        }
#endif
      }
    }
  }
#endif

#ifdef UT
  if (qt >= 0 && qu <= 1) {
    qu = iabs(qu);
    for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
         index0++, i0 = squark_type[qu][index0]) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           index1++, i1 = squark_type[qt][index1]) {
        struct Coupling Cw[4] = {0, 0, 0, 0};
        Cw[1] = params->GLqSQ[bb][i0];
        Cw[0] = params->CHSQq[jj][i0][aa];
        Cw[2] = params->CHqSQ[jj][bb][i1];
        Cw[3] = params->GLSQq[i1][aa];
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
#ifdef SQP2
        virt += 2.0 * g3s *
                ff->Mut_SQP2_GLGA(ff->m2s - 2.0 * ff->pap2, 0.0,
                                  params->mSQs[i0], params->mSQs[i0],
                                  params->mSQs[i0], IEPS, params->mSQs[i0]);
#endif

#ifdef TR1
        virt += 2.0 * g3s *
                ff->Mut_V1a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        virt += 2.0 * g3s *
                ff->Mut_V1b_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, MQS,
                                 MGS, params->mSQs[i0], IEPS);
        if (IEPS == 1) {
          virt -= 2.0 * g3s *
                  ff->Mut_V1a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
          virt -= 2.0 * g3s *
                  ff->Mut_V1b_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, 100.0,
                                   100.0, params->mSQs[i0], IEPS);
        }
#endif

#ifdef TR3
        virt += 2 * g3s *
                ff->Mut_V3b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s, MQS,
                                 MGS, params->mGLs, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mut_V3b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s, 100.0,
                                   100.0, params->mGLs, IEPS);
        }
#endif

#ifdef TR4
        virt += 2 * g3s *
                ff->Mut_V4b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                 params->mGLs, params->mSQs[i0], MGS, IEPS);
        if (IEPS == 1) {
          virt -= 2 * g3s *
                  ff->Mut_V4b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                   params->mGLs, params->mSQs[i0], 100.0, IEPS);
        }
#endif

#ifdef BO1
        virt += 2 * g3s *
                ff->Mut_B1_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                ff->m2s - 2.0 * ff->pap2, MQS, MGS, MQS,
                                params->mSQs[i0], IEPS);
#endif

#ifdef BO3
        virt +=
            2 * g3s *
            ff->Mut_B3_GLGA(0.0, ff->m1s, 0.0, ff->m2s,
                            ff->m1s - 2.0 * ff->pap1, ff->m1s - 2.0 * ff->pbp1,
                            MQS, MGS, params->mGLs, params->mSQs[i0], IEPS);
#endif

#ifdef BU1
        if (IEPS == 1) {
          virt += 2.0 * cF * 1.0 / (16.0 * pow2(M_PI)) * g3s * ff->Mtu_GLGA();
        }
#endif

#ifdef GLSELFB
        virt += 2 * g3s *
                ff->Mut_GLR2_GLGA(params->mGLs, params->mGLs, MGS, params->mGL,
                                  IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt += 2.0 * g3s / (32.0 * pow2(M_PI)) * (4.0 / 3.0 * NC - 2 * cF) *
                  ff->Mtu_GLGA();
        }
#endif

#ifdef DECOUPLE_HEAVY
        if (IEPS == 0) {
          double murs = params->murs;
          double mgls = params->mGLs;
          double log_squarks = 0.0;
          for (int i = 0; i < 12; i++) {
            log_squarks += log(params->mSQs[i] / murs);
          }
          double mts = params->mqs[5];

          virt += 2.0 * g3s / (16.0 * pow2(M_PI)) *
                  (-log(mgls / murs) - 1.0 / 12.0 * log_squarks -
                   1.0 / 3.0 * log(mts / murs))

                  * ff->Mtu_GLGA();
        }
#endif
      }
    }
  }
#endif
  delete ff;
  return virt;
}

// External quark correction at leg a.
// (same results as Jonathan's; only the two strong couplings are interchanged!
// See e.g. Vbu1 and MVtbu1t etc. )
double Residue_quark_a(const double S, const double T, Parameters *params) {
  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  // Sum over squarks inside the loop
  for (int index2 = 0, j0 = squark_type[aa / 3][0]; index2 < 6;
       j0 = squark_type[aa / 3][++index2]) {
    struct Coupling Cs[2] = {params->GLSQq[j0][aa], params->GLqSQ[aa][j0]};
    if (is_coupling_null(Cs, 2)) {
      continue;
    }
    ff->SetSCoupling(Cs);

#ifdef TT
    for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
         i0 = squark_type[qt][++index0]) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           i1 = squark_type[qt][++index1]) {
        struct Coupling Cw[4] = {0, 0, 0, 0};
        Cw[0] = params->CHqSQ[jj][bb][i0];
        Cw[1] = params->GLSQq[i0][aa];
        Cw[2] = params->CHqSQ[jj][bb][i1];
        Cw[3] = params->GLSQq[i1][aa];
        if (is_coupling_null(Cw, 4)) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

        ff->SetWCoupling(Cw);
        virt += 2 * ff->Mtt_QR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
      }
    }
#endif

#ifdef UU
    for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
         i0 = squark_type[qu][++index0]) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           i1 = squark_type[qu][++index1]) {
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
        virt += 2 * ff->Muu_QR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
      }
    }
#endif

#ifdef TU
    for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
         i0 = squark_type[qt][++index0]) {
      for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
           i1 = squark_type[qu][++index1]) {
        struct Coupling Cw[4] = {0, 0, 0, 0};
        Cw[0] = params->CHqSQ[jj][bb][i0];
        Cw[1] = params->GLSQq[i0][aa];
        Cw[3] = params->GLqSQ[bb][i1];
        Cw[2] = params->CHSQq[jj][i1][aa];
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        virt += 2 * ff->Mtu_QR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
      }
    }
#endif

#ifdef UT
    for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
         i0 = squark_type[qu][++index0]) {
      for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
           i1 = squark_type[qt][++index1]) {
        struct Coupling Cw[4] = {0, 0, 0, 0};
        Cw[1] = params->GLqSQ[bb][i0];
        Cw[0] = params->CHSQq[jj][i0][aa];
        Cw[2] = params->CHqSQ[jj][bb][i1];
        Cw[3] = params->GLSQq[i1][aa];
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        virt += 2 * ff->Mut_QR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
      }
    }
#endif
  }
  delete ff;
  return virt;
}

// External quark correction at leg b.
double Residue_quark_b(const double S, const double T, Parameters *params) {
  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  for (int index2 = 0, j0 = squark_type[bb / 3][0]; index2 < 6;
       j0 = squark_type[bb / 3][++index2]) {
    struct Coupling Cs[2] = {params->GLSQq[j0][bb], params->GLqSQ[bb][j0]};
    if (is_coupling_null(Cs, 2)) {
      continue;
    }
    ff->SetSCoupling(Cs);

#ifdef TT
    if (qt >= 0) {
      for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
           i0 = squark_type[qt][++index0]) {
        for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
             i1 = squark_type[qt][++index1]) {
          struct Coupling Cw[4] = {0, 0, 0, 0};
          Cw[0] = params->CHqSQ[jj][bb][i0];
          Cw[1] = params->GLSQq[i0][aa];
          Cw[2] = params->CHqSQ[jj][bb][i1];
          Cw[3] = params->GLSQq[i1][aa];
          if (is_coupling_null(Cw, 4)) {
            continue;
          }
          ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

          ff->SetWCoupling(Cw);
          virt +=
              2 * ff->Mtt_QBR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
        }
      }
    }
#endif

#ifdef UU
    if (qu <= 1) {
      qu = iabs(qu);
      for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
           i0 = squark_type[qu][++index0]) {
        for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
             i1 = squark_type[qu][++index1]) {
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
          virt +=
              2 * ff->Muu_QBR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
        }
      }
    }
#endif

#ifdef TU
    if (qt >= 0 && qu <= 1) {
      qu = iabs(qu);
      for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
           i0 = squark_type[qt][++index0]) {
        for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
             i1 = squark_type[qu][++index1]) {
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
          virt +=
              2 * ff->Mtu_QBR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
        }
      }
    }
#endif

#ifdef UT
    if (qt >= 0 && qu <= 1) {
      qu = iabs(qu);
      for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
           i0 = squark_type[qu][++index0]) {
        for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
             i1 = squark_type[qt][++index1]) {
          struct Coupling Cw[4] = {0, 0, 0, 0};
          Cw[0] = params->CHSQq[jj][i0][aa];
          Cw[1] = params->GLqSQ[bb][i0];
          Cw[2] = params->CHqSQ[jj][bb][i1];
          Cw[3] = params->GLSQq[i1][aa];
          if (is_coupling_null(Cw, 4) == 1) {
            continue;
          }
          ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
          ff->SetWCoupling(Cw);
          virt +=
              2 * ff->Mut_QBR_GLGA(0.0, params->mGLs, params->mSQs[j0], IEPS);
        }
      }
    }
#endif
  }
  delete ff;
  return virt;
}

// External correction for gluino self-energy
double Residue_gluino_a(const double S, const double T, Parameters *params) {
  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  // Sum over all possible quarks and squarks
  for (int j0 = 0; j0 < 6; j0++) {
    for (int index2 = 0, j1 = squark_type[j0 / 3][0]; index2 < 6;
         j1 = squark_type[j0 / 3][++index2]) {
      struct Coupling Cs[2] = {params->GLSQq[j1][j0], params->GLqSQ[j0][j1]};
      if (is_coupling_null(Cs, 2)) {
        continue;
      }
      ff->SetSCoupling(Cs);

#ifdef TT
      for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
           i0 = squark_type[qt][++index0]) {
        for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
             i1 = squark_type[qt][++index1]) {
          struct Coupling Cw[4] = {0, 0, 0, 0};
          Cw[0] = params->CHqSQ[jj][bb][i0];
          Cw[1] = params->GLSQq[i0][aa];
          Cw[2] = params->CHqSQ[jj][bb][i1];
          Cw[3] = params->GLSQq[i1][aa];
          if (is_coupling_null(Cw, 4)) {
            continue;
          }
          ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);

          ff->SetWCoupling(Cw);
          virt += 2 *
                  ff->Mtt_GLR1_GLGA(params->mGLs, params->mqs[j0],
                                    params->mSQs[j1], params->mGL, IEPS);
        }
      }
#endif

#ifdef UU
      for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
           i0 = squark_type[qu][++index0]) {
        for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
             i1 = squark_type[qu][++index1]) {
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
          virt += 2 *
                  ff->Muu_GLR1_GLGA(params->mGLs, params->mqs[j0],
                                    params->mSQs[j1], params->mGL, IEPS);
        }
      }
#endif

#ifdef TU
      for (int index0 = 0, i0 = squark_type[qt][0]; index0 < 6;
           i0 = squark_type[qt][++index0]) {
        for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
             i1 = squark_type[qu][++index1]) {
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
          virt += 2 *
                  ff->Mtu_GLR1_GLGA(params->mGLs, params->mqs[j0],
                                    params->mSQs[j1], params->mGL, IEPS);
        }
      }
#endif

#ifdef UT
      if (qt >= 0 && qu <= 1) {
        qu = iabs(qu);
        for (int index0 = 0, i0 = squark_type[qu][0]; index0 < 6;
             i0 = squark_type[qu][++index0]) {
          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               i1 = squark_type[qt][++index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][i0][aa];
            Cw[1] = params->GLqSQ[bb][i0];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mut_GLR1_GLGA(params->mGLs, params->mqs[j0],
                                      params->mSQs[j1], params->mGL, IEPS);
          }
        }
      }
#endif
    }
  }
  delete ff;
  return virt;
}

// Internal squark propagator correction
double SQprop_a(const double S, const double T, Parameters *params) {

  double virt = 0.0;

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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int index0 = 0, j0 = squark_type[qt][0]; index0 < 6;
         index0++, j0 = squark_type[qt][index0]) {
      for (int index1 = 0, j1 = squark_type[qt][0]; index1 < 6;
           index1++, j1 = squark_type[qt][index1]) {
        for (int index2 = 0, j2 = quark_type[qt][0]; index2 < 3;
             index2++, j2 = quark_type[qt][index2]) {
          struct Coupling Cs[2] = {params->GLqSQ[j2][j0],
                                   params->GLSQq[j1][j2]};
          if (is_coupling_null(Cs, 2)) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef TT
          for (int index3 = 0, i1 = squark_type[qt][0]; index3 < 6;
               index3++, i1 = squark_type[qt][index3]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][bb][j0];
            Cw[1] = params->GLSQq[j1][aa];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4)) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], 0.0, 0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtt_SQP1_GLGA(ff->m1s - 2.0 * ff->pap1, params->mGLs,
                                      params->mqs[j2], params->mSQs[j0],
                                      params->mSQs[j1], IEPS, j1, j0,
                                      params->mSQs[j1]);
          }
#endif

#ifdef TU
          if (qu <= 1) {
            qu = iabs(qu);
            for (int index3 = 0, i1 = squark_type[qu][0]; index3 < 6;
                 index3++, i1 = squark_type[qu][index3]) {
              struct Coupling Cw[4] = {0, 0, 0, 0};
              Cw[0] = params->CHqSQ[jj][bb][j0];
              Cw[1] = params->GLSQq[j1][aa];
              Cw[2] = params->CHSQq[jj][i1][aa];
              Cw[3] = params->GLqSQ[bb][i1];
              if (is_coupling_null(Cw, 4)) {
                continue;
              }
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], 0.0, 0.0);

              ff->SetWCoupling(Cw);
              virt += 2 *
                      ff->Mtu_SQP1_GLGA(ff->m1s - 2.0 * ff->pap1, params->mGLs,
                                        params->mqs[j2], params->mSQs[j0],
                                        params->mSQs[j1], IEPS, j1, j0,
                                        params->mSQs[j1]);
            }
          }
#endif
        }
      }
    }
  }
  if (qu <= 1) {
    qu = iabs(qu);
    for (int index0 = 0, j0 = squark_type[qu][0]; index0 < 6;
         index0++, j0 = squark_type[qu][index0]) {
      for (int index1 = 0, j1 = squark_type[qu][0]; index1 < 6;
           index1++, j1 = squark_type[qu][index1]) {
        for (int index2 = 0, j2 = quark_type[qu][0]; index2 < 3;
             index2++, j2 = quark_type[qu][index2]) {
          struct Coupling Cs[2] = {params->GLSQq[j0][j2],
                                   params->GLqSQ[j2][j1]};
          if (is_coupling_null(Cs, 2)) {
            continue;
          }
          ff->SetSCoupling(Cs);
#ifdef UU
          for (int index3 = 0, i1 = squark_type[qu][0]; index3 < 6;
               index3++, i1 = squark_type[qu][index3]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][j1][aa];
            Cw[1] = params->GLqSQ[bb][j0];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4)) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], 0.0, 0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Muu_SQP1_GLGA(ff->m2s - 2.0 * ff->pap2, params->mGLs,
                                      params->mqs[j2], params->mSQs[j0],
                                      params->mSQs[j1], IEPS, j1, j0,
                                      params->mSQs[j1]);
          }

#endif

#ifdef UT
          if (qt >= 0) {
            for (int index3 = 0, i1 = squark_type[qt][0]; index3 < 6;
                 index3++, i1 = squark_type[qt][index3]) {
              struct Coupling Cw[4] = {0, 0, 0, 0};
              Cw[0] = params->CHSQq[jj][j1][aa];
              Cw[1] = params->GLqSQ[bb][j0];
              Cw[2] = params->CHqSQ[jj][bb][i1];
              Cw[3] = params->GLSQq[i1][aa];
              if (is_coupling_null(Cw, 4)) {
                continue;
              }
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], 0.0, 0.0);

              ff->SetWCoupling(Cw);
              virt += 2 *
                      ff->Mut_SQP1_GLGA(ff->m2s - 2.0 * ff->pap2, params->mGLs,
                                        params->mqs[j2], params->mSQs[j0],
                                        params->mSQs[j1], IEPS, j1, j0,
                                        params->mSQs[j1]);
            }
          }

#endif
        }
      }
    }
  }
  delete ff;
  return virt;
}

// Vertex correction at leg a.
double Triangle2a(const double S, const double T, Parameters *params) {

  double virt = 0.0;

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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  for (int j0 = 0; j0 < 12; j0++) {
    for (int j1 = 0; j1 < 12; j1++) {
      for (int j2 = 0; j2 < 6; j2++) {
        struct Coupling Cs[2] = {params->GLSQq[j1][aa], params->GLSQq[j0][j2]};
        if (is_coupling_null(Cs, 2)) {
          continue;
        }

        ff->SetSCoupling(Cs);

#ifdef TT
        if (qt >= 0) {
          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               index1++, i1 = squark_type[qt][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][bb][j0];
            Cw[1] = params->GLqSQ[j2][j1];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4)) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], params->Gsq[j0],
                              0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtt_V2a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Mtt_V2a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }
        }
#endif

#ifdef TU
        if (qt >= 0 && qu <= 1) {
          for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
               index1++, i1 = squark_type[qu][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][bb][j0];
            Cw[1] = params->GLqSQ[j2][j1];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], params->Gsq[j0],
                              0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtu_V2a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Mtu_V2a_GLGA(0.0, ff->m1s - 2. * ff->pap1, ff->m1s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }
        }
#endif
      }
    }
  }

  for (int j0 = 0; j0 < 12; j0++) {
    for (int j1 = 0; j1 < 12; j1++) {
      for (int j2 = 0; j2 < 6; j2++) {
        struct Coupling Cs[2] = {params->GLSQq[j1][aa], params->GLSQq[j0][j2]};
        if (is_coupling_null(Cs, 2)) {
          continue;
        }

        ff->SetSCoupling(Cs);

#ifdef UU
        if (qu <= 1) {
          qu = iabs(qu);
          for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
               index1++, i1 = squark_type[qu][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][j2][j1];
            Cw[1] = params->GLqSQ[bb][j0];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], params->Gsq[j0],
                              0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Muu_V2a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Muu_V2a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }
        }
#endif

#ifdef UT
        if (qt >= 0 && qu <= 1) {
          qu = iabs(qu);
          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               index1++, i1 = squark_type[qt][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][j2][j1];
            Cw[1] = params->GLqSQ[bb][j0];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], params->Gsq[j0],
                              0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mut_V2a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Mut_V2a_GLGA(0.0, ff->m2s - 2. * ff->pap2, ff->m2s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }
        }
#endif
      }
    }
  }
  delete ff;
  return virt;
}

// Vertex correction at leg b.
double Triangle2b(const double S, const double T, Parameters *params) {
  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  //  for (int index2 = 0 , j0 = squark_type[qt][0]; index2 < 6;  j0 =
  //  squark_type[qt][++index2]){
  if (qt >= 0) {
    for (int j0 = 0; j0 < 12; j0++) {
      for (int j1 = 0; j1 < 12; j1++) {
        for (int j2 = 0; j2 < 6; j2++) {
          struct Coupling Cs[2] = {params->GLqSQ[bb][j1],
                                   params->GLqSQ[j2][j0]};
          if (is_coupling_null(Cs, 2)) {
            continue;
          }

          ff->SetSCoupling(Cs);

#ifdef TT

          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               index1++, i1 = squark_type[qt][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][j1][j2];
            Cw[1] = params->GLSQq[j0][aa];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4)) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], 0.0, 0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtt_V2b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Mtt_V2b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }

#endif

#ifdef TU
          if (qu <= 1) {
            qu = iabs(qu);
            for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
                 index1++, i1 = squark_type[qu][index1]) {
              struct Coupling Cw[4] = {0, 0, 0, 0};
              Cw[0] = params->CHSQq[jj][j1][j2];
              Cw[1] = params->GLSQq[j0][aa];
              Cw[2] = params->CHSQq[jj][i1][aa];
              Cw[3] = params->GLqSQ[bb][i1];
              if (is_coupling_null(Cw, 4) == 1) {
                continue;
              }
              ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], 0.0, 0.0);

              ff->SetWCoupling(Cw);
              virt +=
                  2 *
                  ff->Mtu_V2b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s,
                                   params->mSQs[j0], params->mGLs, MQS, IEPS);
              if (IEPS == 1) {
                virt -= 2 *
                        ff->Mtu_V2b_GLGA(0.0, ff->m2s - 2. * ff->pbp2, ff->m2s,
                                         params->mSQs[j0], params->mGLs, 100.0,
                                         IEPS);
              }
            }
          }
#endif
        }
      }
    }
  }

  if (qu <= 1) {
    iabs(qu);
    for (int j0 = 0; j0 < 12; j0++) {
      for (int j1 = 0; j1 < 12; j1++) {
        for (int j2 = 0; j2 < 6; j2++) {

          struct Coupling Cs[2] = {params->GLqSQ[bb][j1],
                                   params->GLqSQ[j2][j0]};
          if (is_coupling_null(Cs, 2)) {
            continue;
          }

          ff->SetSCoupling(Cs);

#ifdef UU
          for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
               index1++, i1 = squark_type[qu][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][j0][aa];
            Cw[1] = params->GLSQq[j1][j2];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], 0.0, 0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Muu_V2b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                     params->mSQs[j0], params->mGLs, MQS, IEPS);
            if (IEPS == 1) {
              virt -=
                  2 *
                  ff->Muu_V2b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                   params->mSQs[j0], params->mGLs, 100.0, IEPS);
            }
          }

#endif

#ifdef UT
          if (qt >= 0 && qu <= 1) {
            qu = iabs(qu);
            for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
                 index1++, i1 = squark_type[qt][index1]) {
              struct Coupling Cw[4] = {0, 0, 0, 0};
              Cw[0] = params->CHSQq[jj][j0][aa];
              Cw[1] = params->GLSQq[j1][j2];
              Cw[2] = params->CHqSQ[jj][bb][i1];
              Cw[3] = params->GLSQq[i1][aa];
              if (is_coupling_null(Cw, 4) == 1) {
                continue;
              }
              ff->SetPropagator(params->mSQ[j0], params->mSQ[i1], 0.0, 0.0);
              ff->SetWCoupling(Cw);

              virt +=
                  2 *
                  ff->Mut_V2b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                   params->mSQs[j0], params->mGLs, MQS, IEPS);
              if (IEPS == 1) {
                virt -= 2 *
                        ff->Mut_V2b_GLGA(0.0, ff->m1s - 2. * ff->pbp1, ff->m1s,
                                         params->mSQs[j0], params->mGLs, 100.0,
                                         IEPS);
              }
            }
          }
#endif
        }
      }
    }
  }
  delete ff;
  return virt;
}

// Box correction.
double Box2(const double S, const double T, Parameters *params) {
  double virt = 0.0;
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

  FI *ff = new FI();
  ff->SetKinematic(params->mGL, params->mCH[jj], S, T);

  for (int j0 = 0; j0 < 12; j0++) {
    for (int j1 = 0; j1 < 12; j1++) {
      for (int j2 = 0; j2 < 6; j2++) {
        struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j1][aa]};
        if (is_coupling_null(Cs, 2)) {
          continue;
        }

        ff->SetSCoupling(Cs);

#ifdef TT
        if (qt >= 0) {
          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               index1++, i1 = squark_type[qt][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][j0][j2];
            Cw[1] = params->GLqSQ[j2][j1];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4)) {
              continue;
            }
            ff->SetPropagator(params->mSQ[0], params->mSQ[i1], 0.0, 0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtt_B2_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                    ff->m1s - 2.0 * ff->pap1, params->mSQs[j0],
                                    params->mGLs, params->mSQs[j1], 0.0, IEPS);
          }
        }
#endif

#ifdef TU
        if (qt >= 0 && qu <= 1) {
          qu = iabs(qu);
          for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
               index1++, i1 = squark_type[qu][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHSQq[jj][j0][j2];
            Cw[1] = params->GLqSQ[j2][j1];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[0], params->mSQ[i1], 0.0, 0.0);

            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mtu_B2_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                    ff->m1s - 2.0 * ff->pap1, params->mSQs[j0],
                                    params->mGLs, params->mSQs[j1], 0.0, IEPS);
          }
        }
#endif
      }
    }
  }

  for (int j0 = 0; j0 < 12; j0++) {
    for (int j1 = 0; j1 < 12; j1++) {
      for (int j2 = 0; j2 < 6; j2++) {
        struct Coupling Cs[2] = {params->GLSQq[j0][aa], params->GLqSQ[bb][j1]};
        if (is_coupling_null(Cs, 2)) {
          continue;
        }

        ff->SetSCoupling(Cs);
#ifdef UU
        if (qu <= 1) {
          qu = iabs(qu);
          for (int index1 = 0, i1 = squark_type[qu][0]; index1 < 6;
               index1++, i1 = squark_type[qu][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][j2][j0];
            Cw[1] = params->GLSQq[j1][j2];
            Cw[2] = params->CHSQq[jj][i1][aa];
            Cw[3] = params->GLqSQ[bb][i1];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[0], params->mSQ[i1], 0.0, 0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Muu_B2_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                    ff->m1s - 2.0 * ff->pbp1, params->mSQs[j0],
                                    params->mGLs, params->mSQs[j1], 0.0, IEPS);
          }
        }
#endif

#ifdef UT
        if (qt >= 0 && qu <= 1) {
          qu = iabs(qu);
          for (int index1 = 0, i1 = squark_type[qt][0]; index1 < 6;
               index1++, i1 = squark_type[qt][index1]) {
            struct Coupling Cw[4] = {0, 0, 0, 0};
            Cw[0] = params->CHqSQ[jj][j2][j0];
            Cw[1] = params->GLSQq[j1][j2];
            Cw[2] = params->CHqSQ[jj][bb][i1];
            Cw[3] = params->GLSQq[i1][aa];
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[0], params->mSQ[i1], 0.0, 0.0);
            ff->SetWCoupling(Cw);
            virt += 2 *
                    ff->Mut_B2_GLGA(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                    ff->m1s - 2.0 * ff->pbp1, params->mSQs[j0],
                                    params->mGLs, params->mSQs[j1], 0.0, IEPS);
          }
        }
#endif
      }
    }
  }
  delete ff;
  return virt;
}

// Function including all virtual contributions.
double Virt_gaugino_gluino(const double S, const double T, Parameters *params) {
  const double g3s = norm(params->gqq[0][0].R);

  double result = 0.0;

#ifdef QSELF
  result += Residue_quark_a(S, T, params) + Residue_quark_b(S, T, params);
#endif
#ifdef GLSELFA
  result += Residue_gluino_a(S, T, params);
#endif
#ifdef SQP1
  result += SQprop_a(S, T, params);
#endif
#ifdef TR2A
  result += Triangle2a(S, T, params);
#endif
#ifdef TR2B
  result += Triangle2b(S, T, params);
#endif
#ifdef BO2
  result += Box2(S, T, params);
#endif

  result += qcd(S, T, params);

  return result / (g3s * g3s);
}

// Integrated Dipole \int dsigma^A.
double DipI_gagl(const double S, const double T, Parameters *params) {

#ifdef DIPOLE
  const double g3s = norm(params->gqq[0][0].R);
  double murs = params->murs;
  double mgl = sqrt(params->mGLs);
  double sab = S;
  double sa1 = sja_gagl(S, T, params);
  double sb1 = sjb_gagl(S, T, params);
  double born0 = born_gagl(S, T, params);
  // matching factor leads to a total factor of alphaS / (2 pi) in IV in hxs.cc.
  double virtual_matching_factor = 1.0 / g3s * 1.0 / (8.0 * pow2(M_PI));

  // 1/eps^2 divergent term.
  if (IEPS == 2) {
    return virtual_matching_factor * (2 * born0 * cF);
  }
  // 1/eps divergent terms.
  else if (IEPS == 1) {
    return virtual_matching_factor *
           ((born0 * (2 * cA + 6 * cF + cA * log(pow(mgl, 2) / sa1) +
                      cA * log(murs / sa1) - 2 * cA * log(murs / sab) +
                      4 * cF * log(murs / sab) + cA * log(pow(mgl, 2) / sb1) +
                      cA * log(murs / sb1))) /
            2.);
  }
  // finite terms.
  else if (IEPS == 0) {
    return virtual_matching_factor *
           ((born0 *
             (2 * (-6 * cF * (-10 + pow(M_PI, 2)) +
                   cA * (18 - 2 * pow(M_PI, 2) -
                         (9 * mgl) / (mgl + sqrt(pow(mgl, 2) + sa1)) -
                         (9 * mgl) / (mgl + sqrt(pow(mgl, 2) + sb1)))) +
              3 * (-4 * cA * gsl_sf_dilog(sa1 / (pow(mgl, 2) + sa1)) -
                   4 * cA * gsl_sf_dilog(sb1 / (pow(mgl, 2) + sb1)) +
                   2 * cA * log(pow(mgl, 2) / murs) -
                   cA * pow(log(pow(mgl, 2) / sa1), 2) +
                   6 * cA * log(murs / sa1) +
                   2 * cA * log(pow(mgl, 2) / sa1) * log(murs / sa1) +
                   cA * pow(log(murs / sa1), 2) -
                   (2 * cA * pow(mgl, 2) *
                    log(pow(mgl, 2) / (pow(mgl, 2) + sa1))) /
                       sa1 +
                   2 * cA * log(sa1 / (pow(mgl, 2) + sa1)) -
                   2 * cA * log(pow(mgl, 2) / sa1) *
                       log(sa1 / (pow(mgl, 2) + sa1)) -
                   2 * cA * log(pow(mgl, 2) / (pow(mgl, 2) + sa1)) *
                       log(sa1 / (pow(mgl, 2) + sa1)) -
                   6 * cA * log(1 - mgl / sqrt(pow(mgl, 2) + sa1)) -
                   6 * cA * log(murs / sab) + 12 * cF * log(murs / sab) -
                   2 * cA * pow(log(murs / sab), 2) +
                   4 * cF * pow(log(murs / sab), 2) -
                   cA * pow(log(pow(mgl, 2) / sb1), 2) +
                   6 * cA * log(murs / sb1) +
                   2 * cA * log(pow(mgl, 2) / sb1) * log(murs / sb1) +
                   cA * pow(log(murs / sb1), 2) -
                   (2 * cA * pow(mgl, 2) *
                    log(pow(mgl, 2) / (pow(mgl, 2) + sb1))) /
                       sb1 +
                   2 * cA * log(sb1 / (pow(mgl, 2) + sb1)) -
                   2 * cA * log(pow(mgl, 2) / sb1) *
                       log(sb1 / (pow(mgl, 2) + sb1)) -
                   2 * cA * log(pow(mgl, 2) / (pow(mgl, 2) + sb1)) *
                       log(sb1 / (pow(mgl, 2) + sb1)) -
                   6 * cA * log(1 - mgl / sqrt(pow(mgl, 2) + sb1))))) /
            12.);

  } else {
    return 0.0;
  }
#endif
  return 0.0;
}
