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

#include "dipoles.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"
#include "constants.h"

#define ONSUB

#define WIDTH 1.0E-2

#define SS
#define TT
#define UU
#define ST
#define SU
#define TU

#define cF 1.0
#define TR 1.0

double real_gluon_gauginos(const double S, const double M2, const double PT2,
                           const double TH, const double PH, const int YS,
                           Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel.
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel.

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, M2, PT2, TH, PH, YS);

// Partonic cross section for gauginos + gluon.

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
          params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }
      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MGss();
    }
  }
#endif
#ifdef TT
  if (qt >= 0) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += ff->MGtt();
      }
    }
  }
#endif
#ifdef UU
  if (qu <= 1) {
    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += ff->MGuu();
      }
    }
  }
#endif
#ifdef ST
  if (qt >= 0) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MGst();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MGsu();
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }
        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1], 0.0, 0.0);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MGtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double real_quark_gauginos(const double S, const double M2, const double PT2,
                           const double TH, const double PH, const int YS,
                           Parameters *params) {
  double born = 0.0;

  const int aa = params->in1;
  const int bb = params->in2;
  const int ii = params->out1;
  const int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel.
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel.

  // Sets different kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, M2, PT2, TH, PH, YS);

// Partonic cross section for gauginos + quark.

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
          params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MQss();
    }
  }
#endif
#ifdef TT
  if (qt >= 0) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += ff->MQtt();
#ifdef ONSUB
        if (i0 == i1) {
          born -= ff->MQttp(); // On-shell substraction
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
        struct Coupling Cw[4] = {
            params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += ff->MQuu();
#ifdef ONSUB
        if (i0 == i1) {
          born -= ff->MQuup(); // On-shell substraction
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
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0,
                          params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQst();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0,
                          params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQsu();
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);

        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double real_quarkb_gauginos(const double S, const double M2, const double PT2,
                            const double TH, const double PH, const int YS,
                            Parameters *params) {
  double born = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  int qs = iabs(aa / 3 - bb / 3); // Vector charge for s-channel.
  int qt = aa / 3 - ii / 4;       // Squark charge for t-channel.
  int qu = aa / 3 + jj / 4;       // Squark charge for u-channel.

  // Set kinematics for 2 -> 3 process.
  FI *ff = new FI();
  ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, M2, PT2, TH, PH, YS);

// Partonic cross section for gauginos + antiquark.

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
          params->vCHCH[i1][ii][jj], params->vqq[i1][bb][aa]};
      if (is_coupling_null(Cw, 4) == 1) {
        continue;
      }

      ff->SetPropagator(params->mv[i0], params->mv[i1], 0.0, 0.0);
      ff->SetWCoupling(Cw);
      born += ff->MQBss();
    }
  }
#endif
#ifdef TT
  if (qt >= 0) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += ff->MQBtt();

#ifdef ONSUB
        if (i0 == i1) {
          born -= ff->MQBttp(); // On-shell substraction
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
        struct Coupling Cw[4] = {
            params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);

        ff->SetWCoupling(Cw);
        born += ff->MQBuu();
#ifdef ONSUB
        if (i0 == i1) {
          born -= ff->MQBuup(); // On-shell substraction
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
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[jj][bb][i1], params->CHSQq[ii][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0,
                          params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQBst();
      }
    }
  }
#endif
#ifdef SU
  if (qu <= 1) {
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mv[i0], params->mSQ[i1], 0.0,
                          params->mSQ[i1] * WIDTH);
        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQBsu();
      }
    }
  }
#endif
#ifdef TU
  if (qt >= 0 && qu <= 1) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
            params->CHqSQ[ii][bb][i1], params->CHSQq[jj][i1][aa]};
        if (is_coupling_null(Cw, 4) == 1) {
          continue;
        }

        ff->SetPropagator(params->mSQ[i0], params->mSQ[i1],
                          params->mSQ[i0] * WIDTH, params->mSQ[i1] * WIDTH);

        ff->SetWCoupling(Cw);
        born += 2.0 * ff->MQBtu();
      }
    }
  }
#endif

  delete ff;

  return born;
}

double Dip_GAUGINOS(DipoleType Emitter_Spectator, PartonType emitter,
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
  int qs = iabs(aa / 3 - bb / 3);
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

#ifdef SS
  for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
    for (int i1 = 2 * qs; i1 < qs + 2; i1++) {
      struct Coupling Cw[4] = {
          params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
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
  if (qt >= 0) {
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
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
  if (qu <= 1) {
    for (int i0 = 6 * qu; i0 < 6 * qu + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[ii][bb][i0], params->CHSQq[jj][i0][aa],
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
      for (int i1 = 6 * qt; i1 < 6 * qt + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
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
    for (int i0 = 2 * qs; i0 < qs + 2; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->vCHCH[i0][ii][jj], params->vqq[i0][bb][aa],
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
    for (int i0 = 6 * qt; i0 < 6 * qt + 6; i0++) {
      for (int i1 = 6 * qu; i1 < 6 * qu + 6; i1++) {
        struct Coupling Cw[4] = {
            params->CHqSQ[jj][bb][i0], params->CHSQq[ii][i0][aa],
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

  return 1. / (4 * M_PI) * Dipole * born;
}
