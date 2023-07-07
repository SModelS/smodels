// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the virtual part of the partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "kinematics.h"
#include "params.h"
#include "pxs.h"
#include "utils.h"

using namespace std;

#define IEPS 0 // get the finite piece of the virtuals
#define DRBAR  // finite shift from DREG to DRBAR

// activate or deactivate the different channels
// (for simplified models with decoupled squarks only SS channel is relevant)
#define SS
#define TT
#define UU
#define ST
#define SU
#define TU

// the integrated dipole and the virtual correction
#define DIPOLE
#define VIR

// different pieces of the virtual correction
// Vqcd
#define BO1
#define TR1S
#define TR1
#define TR2
#define BU3
#define BU3C

// left and right handed sensitive
#define BU1A
#define BU1B

#define BU4
#define BU4C

// Susy standard qcd-vertex Corr. counterpart diagram
#define TR2S
#define TR2A
#define TR2B

#define BO2

// Higgs
#define TR3

// squark tadople (zero for non-mixing squarks)
#define BU5

// returns all the infrared divergent virtual qcd corrections
/* includes:
   - squark self energy with squark-gluon in loop (bu3)
   - the qcd triangle (tr1) in the s channel: gluon, quark, quark
   - and tr1 in u and t channel: gluon, quark, squark
   - and tr2 in u and t channel: gluino, squark, quark
   - bo1 in t and u channel: gluon, quark, quark squark,
*/
double Vqcd(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  double g3s = norm(params->gqq[0][0].R);

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  // gluon-quark-quark coupling for strong corrections
  // don't get confused: it is also used for gSQSQ in MVttr1t
  struct Coupling Cs[2] = {params->gqq[bb][bb], params->gqq[aa][aa]};
  if (is_coupling_null(Cs, 2) == 1) {
    return virt;
  }

  ff->SetSCoupling(Cs);

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                               params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], params->Gv[i0], params->Gv[i1]);
      ff->SetWCoupling(Cw);
#ifdef TR1S
      virt += ff->MVstr1s(0.0, 0.0, 2.0 * ff->papb, 0.0, 0.0, 0.0, IEPS);
#endif
    }
  }
#endif
#ifdef TT
  if (qt >= 0) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
#ifdef BO1
        virt += ff->MVtbo1t(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb, ff->m1s - 2.0 * ff->pap1,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i0], IEPS);
#endif
#ifdef TR1
        virt += ff->MVttr1t(0.0, ff->m1s - 2. * ff->pap1, ff->m1s, params->mqs[aa], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef TR2
        virt += ff->MVttr2t(0.0, ff->m1s - 2. * ff->pap1, ff->m2s, params->mqs[bb], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef BU3

        virt += ff->MVtbu3t(ff->m1s - 2.0 * ff->pap1, params->mSQs[i0], 0.0, IEPS);

#endif
#ifdef BU3C

        virt += ff->MVtct3t(params->mSQs[i0], params->mSQs[i0], 0.0, IEPS);

#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt -= g3s * ff->MBtt(); // DRbar
        }
#endif
      }
    }
  }
#endif

#ifdef UU
  if (qu <= 1) {
    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
#ifdef BO1
        virt += ff->MVubo1u(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb, ff->m2s - 2.0 * ff->pap2,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i0], IEPS);
#endif
#ifdef TR1
        virt += ff->MVutr1u(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mqs[aa], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef TR2
        virt += ff->MVutr2u(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, params->mqs[bb], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef BU3
        virt += ff->MVubu3u(ff->m2s - 2.0 * ff->pap2, params->mSQs[i0], 0.0, IEPS);
#endif
#ifdef BU3C
        virt += ff->MVuct3u(params->mSQs[i0], params->mSQs[i0], 0.0, IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt -= g3s * ff->MBuu(); // DRbar
        }
#endif
      }
    }
  }
#endif

#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                  params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        struct Coupling Cw2[4] = {params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa],
                                  params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
        if (is_coupling_null(Cw1, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw1);
#ifdef TR1S
        virt += ff->MVstr1t(0.0, 0.0, 2.0 * ff->papb, 0.0, 0.0, 0.0, IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt -= g3s * ff->MBst(); // DRbar
        }
#endif

        ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                          params->Gv[i0]);
        ff->SetWCoupling(Cw2);
#ifdef BO1
        virt += ff->MVtbo1s(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb, ff->m1s - 2.0 * ff->pap1,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i1], IEPS);
#endif
#ifdef TR1
        virt += ff->MVttr1s(0.0, ff->m1s - 2.0 * ff->pap1, ff->m1s, params->mqs[aa], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef TR2
        virt += ff->MVttr2s(0.0, ff->m1s - 2. * ff->pap1, ff->m2s, params->mqs[bb], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef BU3
        virt += ff->MVtbu3s(ff->m1s - 2.0 * ff->pap1, params->mSQs[i1], 0.0, IEPS);
#endif
#ifdef BU3C
        virt += ff->MVtct3s(params->mSQs[i1], params->mSQs[i1], 0.0, IEPS);
#endif
      }
    }
  }
#endif

#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                  params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        struct Coupling Cw2[4] = {params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa],
                                  params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
        if (is_coupling_null(Cw1, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw1);
#ifdef TR1S
        virt += ff->MVstr1u(0.0, 0.0, 2.0 * ff->papb, params->mqs[bb], 0.0, params->mqs[aa], IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt -= g3s * ff->MBsu(); // DRbar
        }
#endif

        ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                          params->Gv[i0]);
        ff->SetWCoupling(Cw2);
#ifdef BO1
        virt += ff->MVubo1s(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb, ff->m2s - 2.0 * ff->pap2,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i1], IEPS);
#endif
#ifdef TR1
        virt += ff->MVutr1s(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mqs[aa], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef TR2
        virt += ff->MVutr2s(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, params->mqs[bb], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef BU3
        virt += ff->MVubu3s(ff->m2s - 2.0 * ff->pap2, params->mSQs[i1], 0.0, IEPS);
#endif
#ifdef BU3C
        virt += ff->MVuct3s(params->mSQs[i1], params->mSQs[i1], 0.0, IEPS);
#endif
      }
    }
  }
#endif

#ifdef TU
  if (qt >= 0 && qu <= 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw1[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                  params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        struct Coupling Cw2[4] = {params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa],
                                  params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa]};
        if (is_coupling_null(Cw1, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw1);
#ifdef BO1
        virt += ff->MVtbo1u(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb, ff->m1s - 2.0 * ff->pap1,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i0], IEPS);
#endif
#ifdef TR1
        virt += ff->MVttr1u(0.0, ff->m1s - 2.0 * ff->pap1, ff->m1s, params->mqs[aa], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef TR2
        virt += ff->MVttr2u(0.0, ff->m1s - 2. * ff->pap1, ff->m2s, params->mqs[bb], 0.0,
                            params->mSQs[i0], IEPS);
#endif
#ifdef BU3
        virt += ff->MVtbu3u(ff->m1s - 2.0 * ff->pap1, params->mSQs[i1], 0.0, IEPS);
#endif
#ifdef BU3C
        virt += ff->MVtct3u(params->mSQs[i1], params->mSQs[i1], 0.0, IEPS);
#endif

#ifdef DRBAR
        if (IEPS == 0) {
          virt -= 2.0 * g3s * ff->MBtu(); // DRbar
        }
#endif

        ff->SetPropagator(params->mSQ[i1], params->mSQ[i0], params->mSQ[i1] * 1.0e-2,
                          params->mSQ[i0] * 1.0e-2);
        ff->SetWCoupling(Cw2);
#ifdef BO1
        virt += ff->MVubo1t(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb, ff->m2s - 2.0 * ff->pap2,
                            params->mqs[bb], 0.0, params->mqs[aa], params->mSQs[i1], IEPS);
#endif
#ifdef TR1
        virt += ff->MVutr1t(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mqs[aa], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef TR2
        virt += ff->MVutr2t(0.0, ff->m2s - 2. * ff->pap2, ff->m1s, params->mqs[bb], 0.0,
                            params->mSQs[i1], IEPS);
#endif
#ifdef BU3
        virt += ff->MVubu3t(ff->m2s - 2.0 * ff->pap2, params->mSQs[i1], 0.0, IEPS);
#endif
#ifdef BU3C
        virt += ff->MVuct3t(params->mSQs[i1], params->mSQs[i1], 0.0, IEPS);
#endif
      }
    }
  }
#endif

  delete ff;

  return virt;
}

// quark self energy with gluino-squark loop
double Vbu1(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
    for (int j1 = (aa / 3) * 3; j1 < (aa / 3 + 1) * 3; j1++) {
      struct Coupling Cs[2] = {params->GLqSQ[j1][j0], params->GLSQq[j0][aa]};
      if (is_coupling_null(Cs, 2) == 1) {
        continue;
      }
      ff->SetSCoupling(Cs);

      // gaugino pair production

#ifdef SS
      for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
        for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
          struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][j1],
                                   params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
          if (is_coupling_null(Cw, 4) == 1) {
            continue;
          }
          ff->SetPropagator(params->mv[i0], params->mv[i1], params->Gv[i0], params->Gv[i1]);
          ff->SetWCoupling(Cw);
#ifdef BU1A
          virt += ff->MVsbu1s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
        }
      }
#endif
#ifdef TT
      if (qt >= 0) {
        for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][j1],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                              params->mSQ[i1] * 1.0e-2);
            ff->SetWCoupling(Cw);
#ifdef BU1A
            virt += ff->MVtbu1t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
#ifdef UU
      if (qu <= 1) {
        for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][j1],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                              params->mSQ[i1] * 1.0e-2);
            ff->SetWCoupling(Cw);
#ifdef BU1A
            virt += ff->MVubu1u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
#ifdef ST
      if (qt >= 0) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][j1],
                                      params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][j1],
                                      params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw1, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                              params->mSQ[i1] * 1.0e-2);
            ff->SetWCoupling(Cw1);
#ifdef BU1A
            virt += ff->MVsbu1t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                              params->Gv[i0]);
            ff->SetWCoupling(Cw2);
#ifdef BU1A
            virt += ff->MVtbu1s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
#ifdef SU
      if (qu <= 1) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][j1],
                                      params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][j1],
                                      params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw1, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0], 0);
            ff->SetWCoupling(Cw1);
#ifdef BU1A
            virt += ff->MVsbu1u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                              params->Gv[i0]);
            ff->SetWCoupling(Cw2);
#ifdef BU1A
            virt += ff->MVubu1s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
#ifdef TU
      if (qt >= 0 && qu <= 1) {
        for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw1[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][j1],
                                      params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][j1],
                                      params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa]};
            if (is_coupling_null(Cw1, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                              params->mSQ[i1] * 1.0e-2);
            ff->SetWCoupling(Cw1);
#ifdef BU1A
            virt += ff->MVtbu1u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            ff->SetPropagator(params->mSQ[i1], params->mSQ[i0], params->mSQ[i1] * 1.0e-2,
                              params->mSQ[i0] * 1.0e-2);
            ff->SetWCoupling(Cw2);
#ifdef BU1A
            virt += ff->MVubu1t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
    }
  }

  delete ff;

  return virt;
}

// quark self energy or antiquark; same as in Vbu1
double Vbu2(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
    for (int j1 = (bb / 3) * 3; j1 < (bb / 3 + 1) * 3; j1++) {
      struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j0][j1]};
      if (is_coupling_null(Cs, 2) == 1) {
        continue;
      }
      ff->SetSCoupling(Cs);

// Gaugino pair production.
#ifdef SS
      for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
        for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
          struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][j1][aa],
                                   params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
          if (is_coupling_null(Cw, 4) == 0) {
            ff->SetPropagator(params->mv[i0], params->mv[i1], params->Gv[i0], params->Gv[i1]);
            ff->SetWCoupling(Cw);
#ifdef BU1B
            virt += ff->MVsbu2s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
          }
        }
      }
#endif
#ifdef TT
      if (qt >= 0) {
        for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][j1][i0], params->CHSQq[ii][i0][aa],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef BU1B
              virt += ff->MVtbu2t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef UU
      if (qu <= 1) {
        for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][j1][i0], params->CHSQq[jj][i0][aa],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef BU1B
              virt += ff->MVubu2u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef ST
      if (qt >= 0) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][j1][aa],
                                      params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[jj][j1][i1], params->CHSQq[ii][i1][aa],
                                      params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw1, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw1);
#ifdef BU1B
              virt += ff->MVsbu2t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
              ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                                params->Gv[i0]);
              ff->SetWCoupling(Cw2);
#ifdef BU1B
              virt += ff->MVtbu2s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef SU
      if (qu <= 1) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw1[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][j1][aa],
                                      params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[ii][j1][i1], params->CHSQq[jj][i1][aa],
                                      params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw1, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw1);
#ifdef BU1B
              virt += ff->MVsbu2u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif

              ff->SetPropagator(params->mSQ[i1], params->mv[i0], params->mSQ[i1] * 1.0e-2,
                                params->Gv[i0]);
              ff->SetWCoupling(Cw2);
#ifdef BU1B
              virt += ff->MVubu2s(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef TU
      if (qt >= 0 && qu <= 1) {
        for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw1[4] = {params->CHqSQ[jj][j1][i0], params->CHSQq[ii][i0][aa],
                                      params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            struct Coupling Cw2[4] = {params->CHqSQ[ii][j1][i1], params->CHSQq[jj][i1][aa],
                                      params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa]};
            if (is_coupling_null(Cw1, 4) == 0) {
              ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw1);
#ifdef BU1B
              virt += ff->MVtbu2u(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif

              ff->SetPropagator(params->mSQ[i1], params->mSQ[i0], params->mSQ[i1] * 1.0e-2,
                                params->mSQ[i0] * 1.0e-2);
              ff->SetWCoupling(Cw2);
#ifdef BU1B
              virt += ff->MVubu2t(0.0, params->mGLs, params->mSQs[j0], IEPS);
#endif
            }
          }
        }
      }
#endif
    }
  }

  delete ff;

  return virt;
}

double Vbu4(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = 6 * qt; j0 < 6 * qt + 6; j0++) {
      for (int j1 = 6 * qt; j1 < 6 * qt + 6; j1++) {
        for (int j2 = 3 * qt; j2 < 3 * qt + 3; j2++) {
          struct Coupling Cs[2] = {params->GLSQq[j0][j2], params->GLqSQ[j2][j1]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef TT
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                              params->mSQ[i1] * 1.0e-2);
            ff->SetWCoupling(Cw);
#ifdef BU4
            virt += ff->MVtbu4t(ff->m1s - 2.0 * ff->pap1, params->mqs[j2], params->mGLs,
                                params->mSQs[j0], params->mSQs[j1], IEPS);

#endif
            if (j1 == j0) {
#ifdef BU4C
              virt += ff->MVtct4t(params->mSQs[j1], params->mqs[j2], params->mGLs, params->mSQs[j0],
                                  params->mSQs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef ST
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef BU4
              virt += ff->MVtbu4s(ff->m1s - 2.0 * ff->pap1, params->mqs[j2], params->mGLs,
                                  params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
              if (j1 == j0) {
#ifdef BU4C
                virt += ff->MVtct4s(params->mSQs[j1], params->mqs[j2], params->mGLs,
                                    params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
              }
            }
          }
#endif
#ifdef TU
          if (qu <= 1) {
            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
              struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                       params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef BU4
                virt += ff->MVtbu4u(ff->m1s - 2.0 * ff->pap1, params->mqs[j2], params->mGLs,
                                    params->mSQs[j0], params->mSQs[j1], IEPS);

#endif
                if (j1 == j0) {
#ifdef BU4C
                  virt += ff->MVtct4u(params->mSQs[j1], params->mqs[j2], params->mGLs,
                                      params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
                }
              }
            }
          }
#endif
        }
      }
    }
  }
  if (qu <= 1) {
    for (int j0 = 6 * qu; j0 < 6 * qu + 6; j0++) {
      for (int j1 = 6 * qu; j1 < 6 * qu + 6; j1++) {
        for (int j2 = 3 * qu; j2 < 3 * qu + 3; j2++) {
          double mls[5] = {params->mq[j2], params->mGL, params->mSQ[j0], params->mSQ[j1]};
          struct Coupling Cs[2] = {params->GLSQq[j0][j2], params->GLqSQ[j2][j1]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }

          ff->SetSCoupling(Cs);
#ifdef UU
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef BU4
              virt += ff->MVubu4u(ff->m2s - 2.0 * ff->pap2, params->mqs[j2], params->mGLs,
                                  params->mSQs[j0], params->mSQs[j1], IEPS);

#endif
              if (j1 == j0) {
#ifdef BU4C
                virt += ff->MVuct4u(params->mSQs[j1], params->mqs[j2], params->mGLs,
                                    params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
              }
            }
          }
#endif
#ifdef SU
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef BU4
              virt += ff->MVubu4s(ff->m2s - 2.0 * ff->pap2, params->mqs[j2], params->mGLs,
                                  params->mSQs[j0], params->mSQs[j1], IEPS);
#endif

              if (j1 == j0) {
#ifdef BU4C
                virt += ff->MVuct4s(params->mSQs[j1], params->mqs[j2], params->mGLs,
                                    params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
              }
            }
          }
#endif
#ifdef TU
          if (qt >= 0) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
              struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                       params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i0], params->mSQ[i0], params->mSQ[i0] * 1.0e-2,
                                  params->mSQ[i0] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef BU4
                virt += ff->MVubu4t(ff->m2s - 2.0 * ff->pap2, params->mqs[j2], params->mGLs,
                                    params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
                if (j1 == j0) {
#ifdef BU4C
                  virt += ff->MVuct4t(params->mSQs[j1], params->mqs[j2], params->mGLs,
                                      params->mSQs[j0], params->mSQs[j1], IEPS);
#endif
                }
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

// infrared safe diagrams for internal squark self energy with a squark inside
// the loop
// only non-zero for mixing squarks.
double Vbu5(const double S, const double T, Parameters *params) {

  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = 6 * qt; j0 < 6 * qt + 6; j0++) {
      for (int j1 = 6 * qt; j1 < 6 * qt + 6; j1++) {
        if (j0 == j1) {
          continue; // Counterterm: dm2 = diag(SelfEnergy); The counterterm
                    // exactly cancels the diagram
        }
        for (int j2 = 6 * qt; j2 < 6 * qt + 6; j2++) {
          struct Coupling DumC = {0.0, 1.0};
          struct Coupling Cs[2] = {params->SQSQSQ[j0][j2][j1], DumC};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef TT
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 1) {
              continue;
            }
            ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], 0, 0);
            ff->SetWCoupling(Cw);
            virt += ff->MVtbu5t(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
          }

#endif

#ifdef ST
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
              virt += ff->MVtbu5s(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
            }
          }
#endif

#ifdef TU
          if (qu <= 1) {
            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
              struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j0], params->CHSQq[ii][j1][aa],
                                       params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
                virt += ff->MVtbu5u(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
              }
            }
          }

#endif

        } // j2
      }   // j1
    }     // jo
  }       // qt
          // end of bubble in t-channel

  if (qu <= 1) {
    for (int j0 = 6 * qu; j0 < 6 * qu + 6; j0++) {
      for (int j1 = 6 * qu; j1 < 6 * qu + 6; j1++) {
        if (j0 == j1) {
          continue; // Counterterm: dm2 = diag(SelfEnergy)
        }
        for (int j2 = 6 * qu; j2 < 6 * qu + 6; j2++) {
          struct Coupling DumC = {0.0, 1.0};
          struct Coupling Cs[2] = {params->SQSQSQ[j0][j2][j1], DumC};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef UU

          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};

            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
              virt += ff->MVubu5u(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
            }
          }

#endif

#ifdef SU

          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};

            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
              virt += ff->MVubu5s(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
            }
          }
#endif

#ifdef TU

          if (qt >= 0) {

            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {

              struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j0], params->CHSQq[jj][j1][aa],
                                       params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa]

              };
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i0], params->mSQ[i0], params->mSQ[i0] * 1.0e-2,
                                  params->mSQ[i0] * 1.0e-2);
                ff->SetWCoupling(Cw);
                virt += ff->MVubu5t(params->mSQs[j2], params->mSQs[j0], params->mSQs[j1], IEPS);
              }
            }
          }
#endif
        } // j2
      }   // j1
    }     // j0
  }       // end of buuble in u channel

  // end of loop in u-channcel

  delete ff;
  return virt;
}

double Vstr2(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
    for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
      struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j1][aa]};
      if (is_coupling_null(Cs, 2) == 1) {
        continue;
      }
      ff->SetSCoupling(Cs);

#ifdef SS
      for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
        for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
          struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vSQSQ[i0][j0][j1],
                                   params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
          if (is_coupling_null(Cw, 4) == 0) {
            ff->SetPropagator(params->mv[i0], params->mv[i1], params->Gv[i0], params->Gv[i1]);
            ff->SetWCoupling(Cw);
#ifdef TR2S
            virt += ff->MVstr2s(0.0, 0.0, 2.0 * ff->papb, params->mGLs, params->mSQs[j0],
                                params->mSQs[j1], IEPS);
#endif
          }
        }
      }
#endif
#ifdef ST
      if (qt >= 0) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vSQSQ[i0][j0][j1],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                                params->mSQ[i1] * 1.0e-2);
#ifdef TR2S
              ff->SetWCoupling(Cw);
              virt += ff->MVstr2t(0.0, 0.0, 2.0 * ff->papb, params->mGLs, params->mSQs[j0],
                                  params->mSQs[j1], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef SU
      if (qu <= 1) {
        for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vSQSQ[i0][j0][j1],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR2S
              virt += ff->MVstr2u(0.0, 0.0, 2.0 * ff->papb, params->mGLs, params->mSQs[j0],
                                  params->mSQs[j1], IEPS);
#endif
            }
          }
        }
      }
#endif
    }
  }

  delete ff;

  return virt;
}

double Vstr3(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
    for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
      struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j1][aa]};
      if (is_coupling_null(Cs, 2) == 1) {
        continue;
      }
      ff->SetSCoupling(Cs);

#ifdef ST
      if (qt >= 0) {
        for (int i0 = 4 * qs; i0 < 2 * qs + 4; i0++) {
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->hCHCH[i0][ii][jj], params->hSQSQ[i0][j0][j1],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mh[i0], params->mSQ[i1], params->mh[i0] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR3
              virt += ff->MVstr3t(0.0, 0.0, 2.0 * ff->papb, params->mGLs, params->mSQs[j0],
                                  params->mSQs[j1], IEPS);
#endif
            }
          }
        }
      }
#endif
#ifdef SU
      if (qu <= 1) {
        for (int i0 = 4 * qs; i0 < 2 * qs + 4; i0++) {
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->hCHCH[i0][ii][jj], params->hSQSQ[i0][j0][j1],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mh[i0], params->mSQ[i1], params->mh[i0] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR3
              virt += ff->MVstr3u(0.0, 0.0, 2.0 * ff->papb, params->mGLs, params->mSQs[j0],
                                  params->mSQs[j1], IEPS);
#endif
            }
          }
        }
      }
#endif
    }
  }

  delete ff;

  return virt;
}

double Vttr3(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Vttr3 not neeeded for slepton and lepton pair production
  if (ii >= 10 && ii < 30) {
    return 0;
  }

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
      for (int j1 = qt * 3; j1 < qt * 3 + 3; j1++) {
        for (int j2 = qt * 6; j2 < qt * 6 + 6; j2++) {
          struct Coupling Cs[2] = {params->GLSQq[j2][j1], params->GLSQq[j0][aa]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef TT
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j2], params->CHqSQ[ii][j1][j0],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR2A
              virt += ff->MVttr3t(0.0, ff->m1s - 2.0 * ff->pap1, ff->m1s, params->mGLs,
                                  params->mSQs[j0], 0.0, IEPS);
#endif
            }
          }
#endif
#ifdef ST
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j2], params->CHqSQ[ii][j1][j0],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mv[i0], params->mSQ[j2] * 1.0e-2,
                                params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef TR2A
              virt += ff->MVttr3s(0.0, ff->m1s - 2.0 * ff->pap1, ff->m1s, params->mGLs,
                                  params->mSQs[j0], params->mqs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef TU
          if (qu <= 1) {
            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
              struct Coupling Cw[4] = {params->CHqSQ[jj][bb][j2], params->CHqSQ[ii][j1][j0],
                                       params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef TR2A
                virt += ff->MVttr3u(0.0, ff->m1s - 2.0 * ff->pap1, ff->m1s, params->mGLs,
                                    params->mSQs[j0], params->mqs[j1], IEPS);
#endif
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

double Vttr4(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
      for (int j1 = qt * 3; j1 < qt * 3 + 3; j1++) {
        for (int j2 = qt * 6; j2 < qt * 6 + 6; j2++) {
          struct Coupling Cs[2] = {params->GLqSQ[j1][j2], params->GLqSQ[bb][j0]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

          if (params->out1 >= 20) {

          } else if (params->out1 >= 10) {

          } else {

#ifdef TT
            for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
              struct Coupling Cw[4] = {params->CHSQq[jj][j0][j1], params->CHSQq[ii][j2][aa],
                                       params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef TR2B
                virt += ff->MVttr4t(0.0, ff->m1s - 2.0 * ff->pap1, ff->m2s, params->mGLs,
                                    params->mSQs[j0], params->mqs[j1], IEPS);

#endif
              }
            }
#endif
#ifdef ST
            for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
              struct Coupling Cw[4] = {params->CHSQq[jj][j0][j1], params->CHSQq[ii][j2][aa],
                                       params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[j2], params->mv[i0], params->mSQ[j2] * 1.0e-2,
                                  params->Gv[i0]);
                ff->SetWCoupling(Cw);
#ifdef TR2B
                virt += ff->MVttr4s(0.0, ff->m1s - 2.0 * ff->pap1, ff->m2s, params->mGLs,
                                    params->mSQs[j0], params->mqs[j1], IEPS);
#endif
              }
            }
#endif
#ifdef TU
            if (qu <= 1) {
              for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
                struct Coupling Cw[4] = {params->CHSQq[jj][j0][j1], params->CHSQq[ii][j2][aa],
                                         params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
                if (is_coupling_null(Cw, 4) == 0) {
                  ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                    params->mSQ[i1] * 1.0e-2);
                  ff->SetWCoupling(Cw);
#ifdef TR2B
                  virt += ff->MVttr4u(0.0, ff->m1s - 2.0 * ff->pap1, ff->m2s, params->mGLs,
                                      params->mSQs[j0], params->mqs[j1], IEPS);
#endif
                }
              }
            }
#endif
          }
        }
      }
    }
  }
  delete ff;

  return virt;
}

double Vutr3(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = (aa / 3) * 6; j0 < (aa / 3 + 1) * 6; j0++) {
      for (int j1 = qu * 3; j1 < qu * 3 + 3; j1++) {
        for (int j2 = qu * 6; j2 < qu * 6 + 6; j2++) {
          struct Coupling Cs[2] = {params->GLSQq[j2][j1], params->GLSQq[j0][aa]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef UU
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j2], params->CHqSQ[jj][j1][j0],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR2A
              virt += ff->MVutr3u(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mGLs,
                                  params->mSQs[j0], params->mqs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef SU
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j2], params->CHqSQ[jj][j1][j0],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mv[i0], params->mSQ[j2] * 1.0e-2,
                                params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef TR2A
              virt += ff->MVutr3s(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mGLs,
                                  params->mSQs[j0], params->mqs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef TU
          if (qt >= 0) {
            for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
              struct Coupling Cw[4] = {params->CHqSQ[ii][bb][j2], params->CHqSQ[jj][j1][j0],
                                       params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef TR2A
                virt += ff->MVutr3t(0.0, ff->m2s - 2.0 * ff->pap2, ff->m2s, params->mGLs,
                                    params->mSQs[j0], params->mqs[j1], IEPS);
#endif
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

double Vutr4(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Vutr4 not neeeded for slepton and lepton pair production
  if (ii >= 10 && ii < 30) {
    return 0;
  }

  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qu <= 1) {
    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
      for (int j1 = qu * 3; j1 < qu * 3 + 3; j1++) {
        for (int j2 = qu * 6; j2 < qu * 6 + 6; j2++) {
          struct Coupling Cs[2] = {params->GLqSQ[j1][j2], params->GLqSQ[bb][j0]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef UU
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHSQq[ii][j0][j1], params->CHSQq[jj][j2][aa],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef TR2B
              virt += ff->MVutr4u(0.0, ff->m2s - 2.0 * ff->pap2, ff->m1s, params->mGLs,
                                  params->mSQs[j0], params->mqs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef SU
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->CHSQq[ii][j0][j1], params->CHSQq[jj][j2][aa],
                                     params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[j2], params->mv[i0], params->mSQ[j2] * 1.0e-2,
                                params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef TR2B
              virt += ff->MVutr4s(0.0, ff->m2s - 2.0 * ff->pap2, ff->m1s, params->mGLs,
                                  params->mSQs[j0], params->mqs[j1], IEPS);
#endif
            }
          }
#endif
#ifdef TU
          if (qt >= 0) {
            for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
              struct Coupling Cw[4] = {params->CHSQq[ii][j0][j1], params->CHSQq[jj][j2][aa],
                                       params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[j2], params->mSQ[i1], params->mSQ[j2] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef TR2B
                virt += ff->MVutr4t(0.0, ff->m2s - 2.0 * ff->pap2, ff->m1s, params->mGLs,
                                    params->mSQs[j0], params->mqs[j1], IEPS);
#endif
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

double Vtbo2(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qt >= 0) {
    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
      for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
        for (int j2 = qt * 3; j2 < 3 * qt + 3; j2++) {
          struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j1][aa]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef TT
          for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[ii][j2][j1], params->CHSQq[jj][j0][j2],
                                     params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef BO2
              virt += ff->MVtbo2t(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                  ff->m1s - 2.0 * ff->pap1, params->mSQs[j0], params->mGLs,
                                  params->mSQs[j1], params->mqs[j2], IEPS);
#endif
            }
          }
#endif
#ifdef ST
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                     params->CHqSQ[ii][j2][j1], params->CHSQq[jj][j0][j2]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef BO2
              virt += ff->MVtbo2s(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                  ff->m1s - 2.0 * ff->pap1, params->mSQs[j0], params->mGLs,
                                  params->mSQs[j1], params->mqs[j2], IEPS);
#endif
            }
          }
#endif
#ifdef TU
          if (qu <= 1) {
            for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
              struct Coupling Cw[4] = {params->CHqSQ[ii][j2][j1], params->CHSQq[jj][j0][j2],
                                       params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                  params->mSQ[i1] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef BO2
                virt += ff->MVtbo2u(0.0, 0.0, ff->m1s, ff->m2s, 2.0 * ff->papb,
                                    ff->m1s - 2.0 * ff->pap1, params->mSQs[j0], params->mGLs,
                                    params->mSQs[j1], params->mqs[j2], IEPS);
#endif
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

double Vubo2(const double S, const double T, Parameters *params) {
  double virt = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

  if (qu >= 0) {
    for (int j0 = (bb / 3) * 6; j0 < (bb / 3 + 1) * 6; j0++) {
      for (int j1 = (aa / 3) * 6; j1 < (aa / 3 + 1) * 6; j1++) {
        for (int j2 = qu * 3; j2 < 3 * qu + 3; j2++) {
          struct Coupling Cs[2] = {params->GLqSQ[bb][j0], params->GLSQq[j1][aa]};
          if (is_coupling_null(Cs, 2) == 1) {
            continue;
          }
          ff->SetSCoupling(Cs);

#ifdef UU
          for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
            struct Coupling Cw[4] = {params->CHqSQ[jj][j2][j1], params->CHSQq[ii][j0][j2],
                                     params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mSQ[i1], params->mSQ[i1], params->mSQ[i1] * 1.0e-2,
                                params->mSQ[i1] * 1.0e-2);
              ff->SetWCoupling(Cw);
#ifdef BO2
              virt += ff->MVubo2u(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb,
                                  ff->m2s - 2.0 * ff->pap2, params->mSQs[j0], params->mGLs,
                                  params->mSQs[j1], params->mqs[j2], IEPS);
#endif
            }
          }
#endif
#ifdef SU
          for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
            struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                     params->CHqSQ[jj][j2][j1], params->CHSQq[ii][j0][j2]};
            if (is_coupling_null(Cw, 4) == 0) {
              ff->SetPropagator(params->mv[i0], params->mv[i0], params->Gv[i0], params->Gv[i0]);
              ff->SetWCoupling(Cw);
#ifdef BO2
              virt += ff->MVubo2s(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb,
                                  ff->m2s - 2.0 * ff->pap2, params->mSQs[j0], params->mGLs,
                                  params->mSQs[j1], params->mqs[j2], IEPS);
#endif
            }
          }
#endif
#ifdef TU
          if (qt >= 0) {
            for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
              struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                       params->CHqSQ[jj][j2][j1], params->CHSQq[ii][j0][j2]};
              if (is_coupling_null(Cw, 4) == 0) {
                ff->SetPropagator(params->mSQ[i0], params->mSQ[i0], params->mSQ[i0] * 1.0e-2,
                                  params->mSQ[i0] * 1.0e-2);
                ff->SetWCoupling(Cw);
#ifdef BO2
                virt += ff->MVubo2t(0.0, 0.0, ff->m2s, ff->m1s, 2.0 * ff->papb,
                                    ff->m2s - 2.0 * ff->pap2, params->mSQs[j0], params->mGLs,
                                    params->mSQs[j1], params->mqs[j2], IEPS);
#endif
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

double Ctt3(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef TT
  if (qt >= 0 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->dCHSQq[ii][i0][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtt();
      }
    }
  }
#endif
#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[jj][bb][i1], params->dCHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        if (IEPS == 1) {
          born += ff->MBst();
        } else if (IEPS == 0) {
          born += ff->MB1st();
        }
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->dCHSQq[ii][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double Ctt4(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef TT
  if (qt >= 0 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {params->dCHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtt();
      }
    }
  }
#endif
#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->dCHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        if (IEPS == 1) {
          born += ff->MBst();
        } else if (IEPS == 0) {
          born += ff->MB1st();
        }
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->dCHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double Ctu3(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef UU
  if (qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[ii][bb][i0], params->dCHSQq[jj][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBuu();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->CHqSQ[ii][bb][i1], params->dCHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        if (IEPS == 1) {
          born += ff->MBsu();
        } else if (IEPS == 0) {
          born += ff->MB1su();
        }
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->dCHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double Ctu4(const double S, const double T, Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel

  // Sets different Born kinematics
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);

#ifdef UU
  if (qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->dCHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
                                 params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBuu();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
                                 params->dCHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], params->Gv[i0],
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        if (IEPS == 1) {
          born += ff->MBsu();
        } else if (IEPS == 0) {
          born += ff->MB1su();
        }
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1 && IEPS == 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
                                 params->dCHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], params->mSQ[i0] * 1.0e-2,
                          params->mSQ[i1] * 1.0e-2);
        ff->SetWCoupling(Cw);
        born += ff->MBtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double Virt_gauginos(const double S, const double T, Parameters *params) {
  const double g3s = norm(params->gqq[0][0].R);
  double result = 0.0;
#ifdef VIR
  result += Vqcd(S, T, params); // Include bu3, tr1, tr2(t, u) and bo1(t, u)
  result += 0.5 * Vbu1(S, T, params) + 0.5 * Vbu2(S, T, params) + Vbu4(S, T, params);
  result += Vstr2(S, T, params) + (Vstr3(S, T, params) + Vttr3(S, T, params) + Vttr4(S, T, params) +
                                   Vutr3(S, T, params) + Vutr4(S, T, params));
#ifdef BU5
  result += Vbu5(S, T, params);
#endif
  result += (Vtbo2(S, T, params) + Vubo2(S, T, params));
  // vertex counterterms (generates only finite pieces for inteference channels
  // with S)
  result += (Ctt3(S, T, params) + Ctt4(S, T, params) + Ctu3(S, T, params) + Ctu4(S, T, params));
#endif
  return result / g3s;
}

// this function is only used for resummation in res.cc
double Virt2_gauginos(const double S, const double T, Parameters *params) {
  const double g3s = norm(params->gqq[0][0].R);

  double result = 0.0;
  result += Vqcd(S, T, params); // Includes bu3, tr1, tr2(t, u) and bo1(t, u)

  return result / g3s;
}

double DipI_gauginos(const double S, const double T, Parameters *params) {
#ifdef DIPOLE
  const double lnmusq = log(params->murs / S);
  const double born = born_gauginos(S, T, params);
  const double eps1 = born_gauginos_eps1(S, T, params);

  if (IEPS == 2) {
    return 2.0 * born_gauginos(S, T, params);
  } else if (IEPS == 1) {
    return ((3.0 + 2.0 * lnmusq) * born_gauginos(S, T, params) +
            2.0 * born_gauginos_eps1(S, T, params));
  } else if (IEPS == 0) {
    // old: 10 - pi^2 missing: Jonathan shifted those pieces to the coll.
    // remainder.
    // return (( 3.0 + lnmusq) * lnmusq) * born_gauginos(S, T, params)
    //                       + (3 + 2.0 * lnmusq) * born_gauginos_eps1(S, T,
    //                       params)
    //     + 2.0 * born_gauginos_eps2(S, T, params);
    return (10.0 - pow(M_PI, 2) + (3.0 + lnmusq) * lnmusq) * born_gauginos(S, T, params) +
           (3 + 2.0 * lnmusq) * born_gauginos_eps1(S, T, params) +
           2.0 * born_gauginos_eps2(S, T, params);

  } else {
    return 0.0;
  }
#endif
  return 0;
};
