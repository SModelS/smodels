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
#include <vector> // TODO change to <array> for speedup?

using namespace std;

#include "pdf.h"

#include "debug.h"
#include "dipoles.h"
#include "kinematics.h"
#include "params.h"
#include "utils.h"

#include "M_tot_uu_UXu_pow0_1.hpp"
#include "M_tot_uu_UXu_pow0_2.hpp"
#include "M_tot_uu_UXu_pow1_1.hpp"
#include "M_tot_uu_UXu_pow1_2.hpp"
#include "M_tot_uu_UXu_pow2_1.hpp"
#include "M_tot_uu_UXu_pow2_2.hpp"
#include "M_tot_uu_UXu_pow2_3.hpp"
#include "constants.h"

#include "hxs.h"

// see pxs_gaugino_gluino_real.cc
#define PROSPINO
#undef ONSUB
// comment to disable the actual subtraction only
#define ONSUB

/** minimum quark index of emitted quark in crossed case. */
#define BMIN 0
/** maximum quark index of emitted quark in crossed case. */
#define BMAX 5

/** minimum quark index of emitted quark in non-crossed case. */
#define KMIN 0
/** maximum quark index of emitted quark in non-crossed case. */
#define KMAX 5

// set PROSPINO to boolean value
#ifdef PROSPINO
#undef PROSPINO
#define PROSPINO true
#else
#define PROSPINO false
#endif
/**
 * Computes kinematical quantities from #FI.
 * Here no squark propagators are included.
 * depending on \b crossed the resonant propagators are regulated by #WIDTH
 */
#define ZERO_SQUARK_DENOMS                                                                         \
  double MQ, MX, MG;                                                                               \
  double K2K3, P1K2, P1K3, P1K1, dM2o, s, s1, QK2, P2K2, QK1, P2K3;                                \
  double K1Q = ff->m1s + ff->p1p2 + ff->p1p3;                                                      \
  QK1 = K1Q;                                                                                       \
  MX = ff->m2;                                                                                     \
  MQ = ff->m1;                                                                                     \
  MG = params->mGL;                                                                                \
  dM2o = MX * MX - MQ * MQ;                                                                        \
  s = 2 * ff->papb;                                                                                \
  s1 = 2 * ff->p1p2 + ff->m1s + ff->m2s;                                                           \
  QK2 = (s + s1) / 2 - K1Q;                                                                        \
  P1K2 = ff->pap2;                                                                                 \
  P1K3 = ff->pap3;                                                                                 \
  K2K3 = ff->p2p3;                                                                                 \
  P2K2 = ff->pbp2;                                                                                 \
  P2K3 = ff->pbp3;                                                                                 \
  P1K1 = ff->pap1;                                                                                 \
  complex<double> Dn_k2pk3MU = s - 2 * K1Q;                                                        \
  auto k2pk3MU = Dn_k2pk3MU;                                                                       \
  if (dM2o < 0 && (!PROSPINO || 2.0 * ff->papb >= pow2(MQ + ff->m1)) && !crossed)                  \
    Dn_k2pk3MU += MQ * MQ * WIDTH * III;                                                           \
  complex<double> rDn_k2pk3MU_pow1 = 1.0 / Dn_k2pk3MU;                                             \
  complex<double> rDn_k2pk3MU_pow2 = 1.0 / norm(Dn_k2pk3MU);                                       \
  double Dn_p1mk2_MQ = MX * MX - 2 * P1K2 - MQ * MQ;                                               \
  double Dn_p1mk1_MG = MQ * MQ - 2 * P1K1 - MG * MG;                                               \
  complex<double> Dn_qmk2MG = MQ * MQ + 2 * ff->p1p3 - MG * MG;                                    \
  auto qmk2MG = Dn_qmk2MG;                                                                         \
  if (MQ * MQ < MG * MG && (!PROSPINO || 2.0 * ff->papb >= pow2(MG + ff->m2)) && !crossed)         \
    Dn_qmk2MG -= MG * MG * WIDTH * III; /*flip? + k2pk3MU with MUi1*/                              \
  complex<double> rDn_qmk2MG_pow1 = 1.0 / Dn_qmk2MG;                                               \
  complex<double> rDn_qmk2MG_pow2 = 1.0 / norm(Dn_qmk2MG);                                         \
  double Dn_p2mk3 = -2 * P2K3;
/**
 * Computes kinematical quantities from #FI.
 * Here a single squark propagator is included.
 * depending on \b crossed the resonant propagators are regulated by #WIDTH
 */
#define SINGLE_SQUARK_DENOMS                                                                       \
  double MQi1 = params->mSQ[m];                                                                    \
  complex<double> Dn_p1mk2MUi1 = MX * MX - MQi1 * MQi1 - 2 * P1K2;                                 \
  complex<double> Dn_qmp1mk2MUi1 = MX * MX - MQi1 * MQi1 - 2 * P2K2;                               \
  auto qmp1mk2MUi1 = Dn_qmp1mk2MUi1;                                                               \
  if (MX * MX < MQi1 * MQi1 && (!PROSPINO || -2.0 * ff->pap3 >= pow2(MQi1 + ff->m1)) && crossed)   \
    Dn_qmp1mk2MUi1 += MQi1 * MQi1 * WIDTH * III;                                                   \
  complex<double> rDn_p1mk2MUi1_pow1 = 1.0 / Dn_p1mk2MUi1;                                         \
  complex<double> rDn_qmp1mk2MUi1_pow1 = 1.0 / Dn_qmp1mk2MUi1;                                     \
  complex<double> Dn_k2pk3MUi1 = s - 2 * K1Q + MQ * MQ - MQi1 * MQi1;                              \
  auto k2pk3MUi1 = Dn_k2pk3MUi1;                                                                   \
  if (MX * MX < MQi1 * MQi1 && (!PROSPINO || 2.0 * ff->papb >= pow2(MQi1 + ff->m1)) && !crossed)   \
    Dn_k2pk3MUi1 += MQi1 * MQi1 * WIDTH * III;                                                     \
  complex<double> rDn_k2pk3MUi1_pow1 = 1.0 / Dn_k2pk3MUi1;                                         \
  complex<double> rDn_k2pk3MUi1_pow2 = 1.0 / norm(Dn_k2pk3MUi1);

/**
 * Computes kinematical quantities from #FI.
 * Here two squark propagators are included.
 * depending on \b crossed the resonant propagators are regulated by #WIDTH
 */
#define DOUBLE_SQUARK_DENOMS                                                                       \
  double MQi1 = params->mSQ[m];                                                                    \
  complex<double> Dn_p1mk2MUi1 = MX * MX - MQi1 * MQi1 - 2 * P1K2;                                 \
  complex<double> Dn_qmp1mk2MUi1 = MX * MX - MQi1 * MQi1 - 2 * P2K2;                               \
  auto qmp1mk2MUi1 = Dn_qmp1mk2MUi1;                                                               \
  if (MX * MX < MQi1 * MQi1 && (!PROSPINO || -2.0 * ff->pap3 >= pow2(MQi1 + ff->m1)) && crossed)   \
    Dn_qmp1mk2MUi1 += MQi1 * MQi1 * WIDTH * III;                                                   \
  complex<double> rDn_p1mk2MUi1_pow1 = 1.0 / Dn_p1mk2MUi1;                                         \
  complex<double> rDn_qmp1mk2MUi1_pow1 = 1.0 / Dn_qmp1mk2MUi1;                                     \
  complex<double> Dn_k2pk3MUi1 = s - 2 * K1Q + MQ * MQ - MQi1 * MQi1;                              \
  auto k2pk3MUi1 = Dn_k2pk3MUi1;                                                                   \
  if (MX * MX < MQi1 * MQi1 && (!PROSPINO || 2.0 * ff->papb >= pow2(MQi1 + ff->m1)) && !crossed)   \
    Dn_k2pk3MUi1 += MQi1 * MQi1 * WIDTH * III;                                                     \
  complex<double> rDn_k2pk3MUi1_pow1 = 1.0 / Dn_k2pk3MUi1;                                         \
  complex<double> rDn_k2pk3MUi1_pow2 = 1.0 / norm(Dn_k2pk3MUi1);                                   \
  double MQj1 = MQi1;                                                                              \
  double MQi2 = params->mSQ[n];                                                                    \
  complex<double> Dn_p1mk2MUj1 = MX * MX - MQj1 * MQj1 - 2 * P1K2;                                 \
  complex<double> Dn_qmp1mk2MUj1 = MX * MX - MQj1 * MQj1 - 2 * P2K2;                               \
  auto qmp1mk2MUj1 = Dn_qmp1mk2MUj1;                                                               \
  if (MX * MX < MQj1 * MQj1 &&                                                                     \
      (!PROSPINO || -2.0 * ff->pap3 >= pow2(MQj1 + ff->m1) ||                                      \
       -2.0 * ff->pap3 >= pow2(MQi2 + ff->m1)) &&                                                  \
      crossed)                                                                                     \
    Dn_qmp1mk2MUj1 += MQj1 * MQj1 * WIDTH * III;                                                   \
  complex<double> rDn_p1mk2MUj1_pow1 = 1.0 / Dn_p1mk2MUj1;                                         \
  complex<double> rDn_qmp1mk2MUj1_pow1 = 1.0 / Dn_qmp1mk2MUj1;                                     \
  complex<double> Dn_p1mk2MUi2 = MX * MX - MQi2 * MQi2 - 2 * P1K2;                                 \
  complex<double> Dn_qmp1mk2MUi2 = MX * MX - MQi2 * MQi2 - 2 * P2K2;                               \
  auto qmp1mk2MUi2 = Dn_qmp1mk2MUi2;                                                               \
  if (MX * MX < MQi2 * MQi2 &&                                                                     \
      (!PROSPINO || -2.0 * ff->pap3 >= pow2(MQi2 + ff->m1) ||                                      \
       -2.0 * ff->pap3 >= pow2(MQj1 + ff->m1)) &&                                                  \
      crossed)                                                                                     \
    Dn_qmp1mk2MUi2 -= MQi2 * MQi2 * WIDTH * III;                                                   \
  complex<double> rDn_p1mk2MUi2_pow1 = 1.0 / Dn_p1mk2MUi2;                                         \
  complex<double> rDn_qmp1mk2MUi2_pow1 = 1.0 / Dn_qmp1mk2MUi2;                                     \
  complex<double> Dn_k2pk3MUi2 = s - 2 * K1Q + MQ * MQ - MQi2 * MQi2;                              \
  auto k2pk3MUi2 = Dn_k2pk3MUi2;                                                                   \
  if (MX * MX < MQi2 * MQi2 &&                                                                     \
      (!PROSPINO || 2.0 * ff->papb >= pow2(MQi2 + ff->m1) ||                                       \
       2.0 * ff->papb >= pow2(MQj1 + ff->m1)) &&                                                   \
      !crossed)                                                                                    \
    Dn_k2pk3MUi2 -= MQi2 * MQi2 * WIDTH * III;                                                     \
  complex<double> rDn_k2pk3MUi2_pow1 = 1.0 / Dn_k2pk3MUi2;                                         \
  complex<double> rDn_k2pk3MUi2_pow2 = 1.0 / norm(Dn_k2pk3MUi2);                                   \
  complex<double> Dn_k2pk3MUj1 = s - 2 * K1Q + MQ * MQ - MQj1 * MQj1;                              \
  auto k2pk3MUj1 = Dn_k2pk3MUj1;                                                                   \
  if (MX * MX < MQj1 * MQj1 &&                                                                     \
      (!PROSPINO || 2.0 * ff->papb >= pow2(MQj1 + ff->m1) ||                                       \
       2.0 * ff->papb >= pow2(MQi2 + ff->m1)) &&                                                   \
      !crossed)                                                                                    \
    Dn_k2pk3MUj1 += MQj1 * MQj1 * WIDTH * III;                                                     \
  complex<double> rDn_k2pk3MUj1_pow1 = 1.0 / Dn_k2pk3MUj1;                                         \
  complex<double> rDn_k2pk3MUj1_pow2 = 1.0 / norm(Dn_k2pk3MUj1);

/**
 * \return new Coupling with flipped sign of the right-handed part of \p c.
 *
 * This is needed as SC defines gluino-couplings differently.
 */
Coupling fix_sign(Coupling c) {
  Coupling r;
  r.L = c.L;
  r.R = -c.R;
  return r;
}
/**
 * Adds \p b to \p a elementwise and saves in \p a.
 */
template <std::size_t N>
void add(std::array<std::complex<double>, N> &a, std::array<std::complex<double>, N> &b) {
  for (int i = 0; i < N; ++i) {
    a[i] += b[i];
  }
}
/**
 * \return True if \p a has any real or imaginary non-zero element.
 */
template <std::size_t N> bool has_non_zero(std::array<std::complex<double>, N> a) {
  for (const auto &elem : a) {
    if (real(elem) != 0.0 || imag(elem) != 0.0)
      return true;
  }
  return false;
}

// TODO explain which coupling belongs to which diagram combination sofar look
// at #M_tot_uu_UXu.h

/** used to calculate all couplings only once at the first_run run
 * (quark/antiquark emission) */
bool first_run = true;
/** used to calculate all crossed couplings only once at the first_run run
 * (quark/antiquark emission)*/
bool first_run_cross = true;
/** all these couplings are used in gaugino_squark/generated/ matrix elements */
std::array<std::array<std::vector<std::tuple<std::array<std::complex<double>, 9>>>, 6>, 6>
    coups_pow0_1;
std::array<std::array<std::vector<std::tuple<std::array<std::complex<double>, 6>>>, 6>, 6>
    coups_pow0_2;
std::array<std::array<std::vector<std::tuple<int, std::array<std::complex<double>, 39>>>, 6>, 6>
    coups_pow1_1;
std::array<std::array<std::vector<std::tuple<int, std::array<std::complex<double>, 50>>>, 6>, 6>
    coups_pow1_2;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 36>>>, 6>,
           6>
    coups_pow2_1;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 77>>>, 6>,
           6>
    coups_pow2_2;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 38>>>, 6>,
           6>
    coups_pow2_3;

/** all these crossed couplings are used in gaugino_squark/generated/ matrix
 * elements */
std::array<std::array<std::vector<std::tuple<std::array<std::complex<double>, 9>>>, 6>, 6>
    coups_pow0_1_cross;
std::array<std::array<std::vector<std::tuple<std::array<std::complex<double>, 6>>>, 6>, 6>
    coups_pow0_2_cross;
std::array<std::array<std::vector<std::tuple<int, std::array<std::complex<double>, 39>>>, 6>, 6>
    coups_pow1_1_cross;
std::array<std::array<std::vector<std::tuple<int, std::array<std::complex<double>, 50>>>, 6>, 6>
    coups_pow1_2_cross;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 36>>>, 6>,
           6>
    coups_pow2_1_cross;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 77>>>, 6>,
           6>
    coups_pow2_2_cross;
std::array<std::array<std::vector<std::tuple<int, int, std::array<std::complex<double>, 38>>>, 6>,
           6>
    coups_pow2_3_cross;

/**
 * \returns gaugino squark quark coupling based on \p ch chargino \p sq squark
 * and \p quark
 *
 * \param params Parameters of current computation/process
 */
Coupling CHSQq(Parameters *params, int ch, int sq, int q) {
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // flip fermion flow - though appears to have no effect
  if (ii > 30) {
    return params->CHSQq[ch][sq][q];
  } else {
    return params->CHqSQ[ch][q][sq];
  }
}
/**
 * \returns gluino squark quark coupling based on \p sq squark and \p quark
 *
 * \param params Parameters of current computation/process
 */
Coupling GLSQq(Parameters *params, int sq, int q) {
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // flip fermion flow - though appears to have no effect
  if (ii > 30) {
    return params->GLSQq[sq][q];
  } else {
    return params->GLqSQ[q][sq];
  }
}

/**
 * Sets the non-crossed couplings based on \p j and \p g for the quark/antiquark
 * emission. The coupling is saved for all 6x6 possible incoming quark/antiquark
 * combinations. The iteration over final state quarks/antiquarks is performed
 * as given by #KMIN and #KMAX.
 *
 * \param params Parameters of current computation/process
 * \param j final state squark
 * \param g final state gaugino
 */
void set_couplings(Parameters *params, int j, int g) {
  int K;
  int aa;
  int bb;
  // reset
  for (int ta = 0; ta < 6; ta++) {
    for (int tb = 0; tb < 6; tb++) {
      coups_pow0_1[ta][tb].clear();
      coups_pow0_2[ta][tb].clear();
      coups_pow1_1[ta][tb].clear();
      coups_pow1_2[ta][tb].clear();
      coups_pow2_1[ta][tb].clear();
      coups_pow2_2[ta][tb].clear();
      coups_pow2_3[ta][tb].clear();
    }
  }

  for (aa = 0; aa < 6; aa++) {
    for (bb = 0; bb < 6; bb++) {

      std::array<std::complex<double>, 9> coup_pow0_1 = {};
      std::array<std::complex<double>, 6> coup_pow0_2 = {};
      for (K = KMIN; K < KMAX; ++K) {
        std::complex<double> LCC5nm4um1U2nuU = CHSQq(params, g, j, aa).L;
        std::complex<double> LCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).L;
        std::complex<double> LCC1nm4u3Um2nuU = CHSQq(params, g, j, K).L;
        std::complex<double> LCC2nm4um6U4nuU = CHSQq(params, g, j, K).L;
        std::complex<double> RCC5nm4um1U2nuU = CHSQq(params, g, j, aa).R;
        std::complex<double> RCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).R;
        std::complex<double> RCC1nm4u3Um2nuU = CHSQq(params, g, j, K).R;
        std::complex<double> RCC2nm4um6U4nuU = CHSQq(params, g, j, K).R;

        double RCC1um3um1g1 = aa == bb;
        double RCC2um3um1g1 = aa == bb;
        double RCC8u4um1g2 = 1.0;
        double RCC1u4um6g2 = 1.0;
        double RCC5um3um6g3 = bb == K;
        double RCC8um3um6g1 = bb == K;

        std::array<std::complex<double>, 9> coup_pow0_1_tmp = {
            (std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
             std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC5um3um6g3, 2.0) * norm(LCC5nm4um1U2nuU) +
              std::pow(RCC5um3um6g3, 2.0) * norm(RCC5nm4um1U2nuU)),
            +(std::pow(RCC1u4um6g2, 2.0) * std::pow(RCC1um3um1g1, 2.0) * norm(LCC1nm4u3Um2nuU) +
              std::pow(RCC1u4um6g2, 2.0) * std::pow(RCC1um3um1g1, 2.0) * norm(RCC1nm4u3Um2nuU)),
            +(std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
              std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
              std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC8u4um1g2, 2.0) * std::pow(RCC8um3um6g1, 2.0) * norm(LCC8nm4u3Um2nuU) +
              std::pow(RCC8u4um1g2, 2.0) * std::pow(RCC8um3um6g1, 2.0) * norm(RCC8nm4u3Um2nuU)),
            +(LCC2nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC2um3um1g1 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC2nm4um6U4nuU * RCC2um3um1g1 * conj(RCC1nm4u3Um2nuU)),
            +(LCC8nm4u3Um2nuU * RCC5um3um6g3 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC5nm4um1U2nuU) +
              RCC5um3um6g3 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC5nm4um1U2nuU)),
            +(LCC2nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC2um3um1g1 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC2nm4um6U4nuU * RCC2um3um1g1 * conj(RCC1nm4u3Um2nuU)),
        };
        add(coup_pow0_1, coup_pow0_1_tmp);
        std::array<std::complex<double>, 6> coup_pow0_2_tmp = {
            (LCC5nm4um1U2nuU * RCC2um3um1g1 * RCC5um3um6g3 * conj(LCC2nm4um6U4nuU) +
             RCC2um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC2nm4um6U4nuU)),
            +(LCC5nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC5um3um6g3 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC1nm4u3Um2nuU)),
            +(LCC5nm4um1U2nuU * RCC2um3um1g1 * RCC5um3um6g3 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC2nm4um6U4nuU)),
            +(LCC8nm4u3Um2nuU * RCC2um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC2nm4um6U4nuU)),
            +(LCC8nm4u3Um2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 *
                  conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                  conj(RCC1nm4u3Um2nuU)),
            +(LCC8nm4u3Um2nuU * RCC2um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC2nm4um6U4nuU)),
        };
        add(coup_pow0_2, coup_pow0_2_tmp);
      }
      if (has_non_zero(coup_pow0_1)) {
        coups_pow0_1[aa][bb].push_back(std::make_tuple(coup_pow0_1));
      }
      if (has_non_zero(coup_pow0_2)) {
        coups_pow0_2[aa][bb].push_back(std::make_tuple(coup_pow0_2));
      }

      for (int m = 0; m < 12; m++) {
        std::array<std::complex<double>, 39> coup_pow1_1 = {};
        std::array<std::complex<double>, 50> coup_pow1_2 = {};
        for (K = KMIN; K < KMAX; K++) {
          std::complex<double> LCC5nm4um1U2nuU = CHSQq(params, g, j, aa).L;
          std::complex<double> LCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).L;
          std::complex<double> LCC1nm4u3Um2nuU = CHSQq(params, g, j, K).L;
          std::complex<double> LCC2nm4um6U4nuU = CHSQq(params, g, j, K).L;
          std::complex<double> RCC5nm4um1U2nuU = CHSQq(params, g, j, aa).R;
          std::complex<double> RCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).R;
          std::complex<double> RCC1nm4u3Um2nuU = CHSQq(params, g, j, K).R;
          std::complex<double> RCC2nm4um6U4nuU = CHSQq(params, g, j, K).R;

          double RCC1um3um1g1 = aa == bb;
          double RCC2um3um1g1 = aa == bb;
          double RCC8u4um1g2 = 1.0;
          double RCC1u4um6g2 = 1.0;
          double RCC5um3um6g3 = bb == K;
          double RCC8um3um6g1 = bb == K;

          complex<double> LCC7um3nm4U1nuU = CHSQq(params, g, m, bb).L;
          complex<double> LCC3um3nm4U3nuU = CHSQq(params, g, m, bb).L;
          complex<double> LCC6nm4um1U2nuU = CHSQq(params, g, m, aa).L;
          complex<double> LCC4nm4um6U4nuU = CHSQq(params, g, m, K).L;
          complex<double> RCC7um3nm4U1nuU = CHSQq(params, g, m, bb).R;
          complex<double> RCC3um3nm4U3nuU = CHSQq(params, g, m, bb).R;
          complex<double> RCC6nm4um1U2nuU = CHSQq(params, g, m, aa).R;
          complex<double> RCC4nm4um6U4nuU = CHSQq(params, g, m, K).R;

          complex<double> LCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
          complex<double> LCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).L;
          complex<double> LCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
          complex<double> LCC7G4um1U2GuU = fix_sign(GLSQq(params, m, aa)).L;
          complex<double> LCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
          complex<double> LCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).L;
          complex<double> LCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
          complex<double> LCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).L;
          complex<double> RCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
          complex<double> RCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).R;
          complex<double> RCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
          complex<double> RCC7G4um1U2GuU = fix_sign(GLSQq(params, m, aa)).R;
          complex<double> RCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
          complex<double> RCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).R;
          complex<double> RCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
          complex<double> RCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).R;

          std::array<std::complex<double>, 39> coup_pow1_1_tmp = {
              (LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                   conj(RCC2nm4um6U4nuU) +
               LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC3um3nm4U3nuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G2um6U4GuU *
                    conj(LCC3um3nm4U3nuU) * conj(RCC1nm4u3Um2nuU) +
                LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC3um3nm4U3nuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G2um6U4GuU *
                    conj(LCC3um3nm4U3nuU) * conj(RCC1nm4u3Um2nuU) +
                LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(LCC8nm4u3Um2nuU) +
                RCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(LCC8nm4u3Um2nuU) +
                RCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4nm4um6U4nuU *
                    conj(LCC4um3G2U3GuU) * conj(RCC1nm4u3Um2nuU) +
                LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4nm4um6U4nuU *
                    conj(LCC4um3G2U3GuU) * conj(RCC1nm4u3Um2nuU) +
                LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),

          };
          add(coup_pow1_1, coup_pow1_1_tmp);

          std::array<std::complex<double>, 50> coup_pow1_2_tmp = {
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(RCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(RCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(RCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(RCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(LCC8nm4u3Um2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(LCC8nm4u3Um2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(LCC8nm4u3Um2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC6nm4um1U2nuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU *
                   conj(LCC6um3G4U1GuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC7um3nm4U1nuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),

          };
          add(coup_pow1_2, coup_pow1_2_tmp);
        }
        if (has_non_zero(coup_pow1_1)) {
          coups_pow1_1[aa][bb].push_back(std::make_tuple(m, coup_pow1_1));
        }
        if (has_non_zero(coup_pow1_2)) {
          coups_pow1_2[aa][bb].push_back(std::make_tuple(m, coup_pow1_2));
        }
      }

      for (int m = 0; m < 12; m++) {
        for (int n = 0; n < 12; n++) {
          std::array<std::complex<double>, 77> coup_pow2_2 = {};
          std::array<std::complex<double>, 38> coup_pow2_3 = {};
          std::array<std::complex<double>, 36> coup_pow2_1 = {};
          for (K = KMIN; K < KMAX; K++) {
            {
              std::complex<double> LCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
              std::complex<double> LCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).L;
              std::complex<double> LCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
              std::complex<double> LCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).L;
              std::complex<double> LCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
              std::complex<double> LCC6um3G4U1GuU = fix_sign(GLSQq(params, n, bb)).L;
              std::complex<double> LCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
              std::complex<double> LCC7G4um1U2GuU = fix_sign(GLSQq(params, n, aa)).L;

              std::complex<double> RCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
              std::complex<double> RCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).R;
              std::complex<double> RCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
              std::complex<double> RCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).R;
              std::complex<double> RCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
              std::complex<double> RCC6um3G4U1GuU = fix_sign(GLSQq(params, n, bb)).R;
              std::complex<double> RCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
              std::complex<double> RCC7G4um1U2GuU = fix_sign(GLSQq(params, n, aa)).R;

              std::complex<double> LCC7um3nm4U1nuU = CHSQq(params, g, n, bb).L;
              std::complex<double> LCC3um3nm4U3nuU = CHSQq(params, g, m, bb).L;
              std::complex<double> LCC6nm4um1U2nuU = CHSQq(params, g, n, aa).L;
              std::complex<double> LCC4nm4um6U4nuU = CHSQq(params, g, m, K).L;
              std::complex<double> RCC7um3nm4U1nuU = CHSQq(params, g, n, bb).R;
              std::complex<double> RCC3um3nm4U3nuU = CHSQq(params, g, m, bb).R;
              std::complex<double> RCC6nm4um1U2nuU = CHSQq(params, g, n, aa).R;
              std::complex<double> RCC4nm4um6U4nuU = CHSQq(params, g, m, K).R;

              std::array<std::complex<double>, 77> coup_pow2_2_tmp = {
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
              };
              add(coup_pow2_2, coup_pow2_2_tmp);

              LCC4um3G2U3GuU = fix_sign(GLSQq(params, n, bb)).L;
              RCC4um3G2U3GuU = fix_sign(GLSQq(params, n, bb)).R;
              LCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).L;
              RCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).R;

              LCC6nm4um1U2nuU = CHSQq(params, g, m, aa).L;
              LCC4nm4um6U4nuU = CHSQq(params, g, n, K).L;
              RCC6nm4um1U2nuU = CHSQq(params, g, m, aa).R;
              RCC4nm4um6U4nuU = CHSQq(params, g, n, K).R;

              std::array<std::complex<double>, 38> coup_pow2_3_tmp = {
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),

              };
              add(coup_pow2_3, coup_pow2_3_tmp);
            }

            {
              std::complex<double> LCC3G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).L);
              std::complex<double> LCC4G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).L);
              std::complex<double> LCC6G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).L);
              std::complex<double> LCC7G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).L);
              std::complex<double> RCC3G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).R);
              std::complex<double> RCC4G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).R);
              std::complex<double> RCC6G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).R);
              std::complex<double> RCC7G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).R);

              std::complex<double> LCC3G2um6U4GuU =
                  fix_sign(GLSQq(params, n, K)).L * conj(fix_sign(GLSQq(params, m, K)).L);
              std::complex<double> LCC4um3G2U3GuU =
                  fix_sign(GLSQq(params, m, bb)).L * conj(fix_sign(GLSQq(params, n, bb)).L);
              std::complex<double> LCC6um3G4U1GuU =
                  fix_sign(GLSQq(params, m, bb)).L * conj(fix_sign(GLSQq(params, n, bb)).L);
              std::complex<double> LCC7G4um1U2GuU =
                  fix_sign(GLSQq(params, n, aa)).L * conj(fix_sign(GLSQq(params, m, aa)).L);

              std::complex<double> RCC3G2um6U4GuU =
                  fix_sign(GLSQq(params, n, K)).R * conj(fix_sign(GLSQq(params, m, K)).R);
              std::complex<double> RCC4um3G2U3GuU =
                  fix_sign(GLSQq(params, m, bb)).R * conj(fix_sign(GLSQq(params, n, bb)).R);
              std::complex<double> RCC6um3G4U1GuU =
                  fix_sign(GLSQq(params, m, bb)).R * conj(fix_sign(GLSQq(params, n, bb)).R);
              std::complex<double> RCC7G4um1U2GuU =
                  fix_sign(GLSQq(params, n, aa)).R * conj(fix_sign(GLSQq(params, m, aa)).R);

              std::complex<double> LCC3um3nm4U3nuU =
                  CHSQq(params, g, m, bb).L * conj(CHSQq(params, g, n, bb).L);
              std::complex<double> LCC4nm4um6U4nuU =
                  CHSQq(params, g, n, K).L * conj(CHSQq(params, g, m, K).L);
              std::complex<double> LCC6nm4um1U2nuU =
                  CHSQq(params, g, n, aa).L * conj(CHSQq(params, g, m, aa).L);
              std::complex<double> LCC7um3nm4U1nuU =
                  CHSQq(params, g, m, bb).L * conj(CHSQq(params, g, n, bb).L);
              std::complex<double> RCC3um3nm4U3nuU =
                  CHSQq(params, g, m, bb).R * conj(CHSQq(params, g, n, bb).R);
              std::complex<double> RCC4nm4um6U4nuU =
                  CHSQq(params, g, n, K).R * conj(CHSQq(params, g, m, K).R);
              std::complex<double> RCC6nm4um1U2nuU =
                  CHSQq(params, g, n, aa).R * conj(CHSQq(params, g, m, aa).R);
              std::complex<double> RCC7um3nm4U1nuU =
                  CHSQq(params, g, m, bb).R * conj(CHSQq(params, g, n, bb).R);

              std::array<std::complex<double>, 36> coup_pow2_1_tmp = {
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) +
                   (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) +
                   (LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) +
                   (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) +
                   (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) +
                   (LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) +
                   (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) +
                   (LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) +
                   (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) +
                   (LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) +
                   (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
              };
              add(coup_pow2_1, coup_pow2_1_tmp);
            }
          }
          if (has_non_zero(coup_pow2_1)) {
            coups_pow2_1[aa][bb].push_back(std::make_tuple(m, n, coup_pow2_1));
          }
          if (has_non_zero(coup_pow2_2)) {
            coups_pow2_2[aa][bb].push_back(std::make_tuple(m, n, coup_pow2_2));
          }
          if (has_non_zero(coup_pow2_3)) {
            coups_pow2_3[aa][bb].push_back(std::make_tuple(m, n, coup_pow2_3));
          }
        }
      }
      // std::cout << "Non-zero elements: " << coups_pow1_1.size() << "/"<< 12<<
      // std::endl; std::cout << "Non-zero elements: " << coups_pow1_2.size() <<
      // "/"<< 12<< std::endl; std::cout << "Non-zero elements: " <<
      // coups_pow2_1.size() << "/"<< 12*12<< std::endl; std::cout << "Non-zero
      // elements: " << coups_pow2_2.size() << "/"<< 12*12<< std::endl;
      // std::cout
      // << "Non-zero elements: " << coups_pow2_3.size() << "/"<< 12*12<<
      // std::endl;
    }
  }
}

/**
 * Sets the non-crossed couplings based on \p j and \p g for the quark/antiquark
 * emission. The coupling is saved for all 6x6 possible incoming quark/antiquark
 * combinations. The iteration over final state quarks/antiquarks is performed
 * as given by #BMIN and #BMAX.
 *
 * \param params Parameters of current computation/process
 * \param j final state squark
 * \param g final state gaugino
 */
void set_couplings_cross(Parameters *params, int j, int g) {
  int K;
  int aa;
  int bb;
  // reset
  for (int ta = 0; ta < 6; ta++) {
    for (int tb = 0; tb < 6; tb++) {
      coups_pow0_1_cross[ta][tb].clear();
      coups_pow0_2_cross[ta][tb].clear();
      coups_pow1_1_cross[ta][tb].clear();
      coups_pow1_2_cross[ta][tb].clear();
      coups_pow2_1_cross[ta][tb].clear();
      coups_pow2_2_cross[ta][tb].clear();
      coups_pow2_3_cross[ta][tb].clear();
    }
  }
  for (aa = 0; aa < 6; aa++) {
    for (K = 0; K < 6; K++) {

      std::array<std::complex<double>, 9> coup_pow0_1 = {};
      std::array<std::complex<double>, 6> coup_pow0_2 = {};
      for (bb = BMIN; bb < BMAX; ++bb) {
        std::complex<double> LCC5nm4um1U2nuU = CHSQq(params, g, j, aa).L;
        std::complex<double> LCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).L;
        std::complex<double> LCC1nm4u3Um2nuU = CHSQq(params, g, j, K).L;
        std::complex<double> LCC2nm4um6U4nuU = CHSQq(params, g, j, K).L;
        std::complex<double> RCC5nm4um1U2nuU = CHSQq(params, g, j, aa).R;
        std::complex<double> RCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).R;
        std::complex<double> RCC1nm4u3Um2nuU = CHSQq(params, g, j, K).R;
        std::complex<double> RCC2nm4um6U4nuU = CHSQq(params, g, j, K).R;

        double RCC1um3um1g1 = aa == bb;
        double RCC2um3um1g1 = aa == bb;
        double RCC8u4um1g2 = 1.0;
        double RCC1u4um6g2 = 1.0;
        double RCC5um3um6g3 = bb == K;
        double RCC8um3um6g1 = bb == K;

        std::array<std::complex<double>, 9> coup_pow0_1_tmp = {
            (std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
             std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC5um3um6g3, 2.0) * norm(LCC5nm4um1U2nuU) +
              std::pow(RCC5um3um6g3, 2.0) * norm(RCC5nm4um1U2nuU)),
            +(std::pow(RCC1u4um6g2, 2.0) * std::pow(RCC1um3um1g1, 2.0) * norm(LCC1nm4u3Um2nuU) +
              std::pow(RCC1u4um6g2, 2.0) * std::pow(RCC1um3um1g1, 2.0) * norm(RCC1nm4u3Um2nuU)),
            +(std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
              std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC2um3um1g1, 2.0) * norm(LCC2nm4um6U4nuU) +
              std::pow(RCC2um3um1g1, 2.0) * norm(RCC2nm4um6U4nuU)),
            +(std::pow(RCC8u4um1g2, 2.0) * std::pow(RCC8um3um6g1, 2.0) * norm(LCC8nm4u3Um2nuU) +
              std::pow(RCC8u4um1g2, 2.0) * std::pow(RCC8um3um6g1, 2.0) * norm(RCC8nm4u3Um2nuU)),
            +(LCC2nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC2um3um1g1 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC2nm4um6U4nuU * RCC2um3um1g1 * conj(RCC1nm4u3Um2nuU)),
            +(LCC8nm4u3Um2nuU * RCC5um3um6g3 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC5nm4um1U2nuU) +
              RCC5um3um6g3 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC5nm4um1U2nuU)),
            +(LCC2nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC2um3um1g1 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC2nm4um6U4nuU * RCC2um3um1g1 * conj(RCC1nm4u3Um2nuU)),
        };
        add(coup_pow0_1, coup_pow0_1_tmp);
        std::array<std::complex<double>, 6> coup_pow0_2_tmp = {
            (LCC5nm4um1U2nuU * RCC2um3um1g1 * RCC5um3um6g3 * conj(LCC2nm4um6U4nuU) +
             RCC2um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC2nm4um6U4nuU)),
            +(LCC5nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC5um3um6g3 * conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC1nm4u3Um2nuU)),
            +(LCC5nm4um1U2nuU * RCC2um3um1g1 * RCC5um3um6g3 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC5nm4um1U2nuU * RCC5um3um6g3 * conj(RCC2nm4um6U4nuU)),
            +(LCC8nm4u3Um2nuU * RCC2um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC2nm4um6U4nuU)),
            +(LCC8nm4u3Um2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 *
                  conj(LCC1nm4u3Um2nuU) +
              RCC1u4um6g2 * RCC1um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                  conj(RCC1nm4u3Um2nuU)),
            +(LCC8nm4u3Um2nuU * RCC2um3um1g1 * RCC8u4um1g2 * RCC8um3um6g1 * conj(LCC2nm4um6U4nuU) +
              RCC2um3um1g1 * RCC8nm4u3Um2nuU * RCC8u4um1g2 * RCC8um3um6g1 * conj(RCC2nm4um6U4nuU)),
        };
        add(coup_pow0_2, coup_pow0_2_tmp);
      }
      if (has_non_zero(coup_pow0_1)) {
        coups_pow0_1_cross[aa][K].push_back(std::make_tuple(coup_pow0_1));
      }
      if (has_non_zero(coup_pow0_2)) {
        coups_pow0_2_cross[aa][K].push_back(std::make_tuple(coup_pow0_2));
      }

      for (int m = 0; m < 12; m++) {
        std::array<std::complex<double>, 39> coup_pow1_1 = {};
        std::array<std::complex<double>, 50> coup_pow1_2 = {};
        for (bb = BMIN; bb < BMAX; bb++) {
          std::complex<double> LCC5nm4um1U2nuU = CHSQq(params, g, j, aa).L;
          std::complex<double> LCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).L;
          std::complex<double> LCC1nm4u3Um2nuU = CHSQq(params, g, j, K).L;
          std::complex<double> LCC2nm4um6U4nuU = CHSQq(params, g, j, K).L;
          std::complex<double> RCC5nm4um1U2nuU = CHSQq(params, g, j, aa).R;
          std::complex<double> RCC8nm4u3Um2nuU = CHSQq(params, g, j, aa).R;
          std::complex<double> RCC1nm4u3Um2nuU = CHSQq(params, g, j, K).R;
          std::complex<double> RCC2nm4um6U4nuU = CHSQq(params, g, j, K).R;

          double RCC1um3um1g1 = aa == bb;
          double RCC2um3um1g1 = aa == bb;
          double RCC8u4um1g2 = 1.0;
          double RCC1u4um6g2 = 1.0;
          double RCC5um3um6g3 = bb == K;
          double RCC8um3um6g1 = bb == K;

          complex<double> LCC7um3nm4U1nuU = CHSQq(params, g, m, bb).L;
          complex<double> LCC3um3nm4U3nuU = CHSQq(params, g, m, bb).L;
          complex<double> LCC6nm4um1U2nuU = CHSQq(params, g, m, aa).L;
          complex<double> LCC4nm4um6U4nuU = CHSQq(params, g, m, K).L;
          complex<double> RCC7um3nm4U1nuU = CHSQq(params, g, m, bb).R;
          complex<double> RCC3um3nm4U3nuU = CHSQq(params, g, m, bb).R;
          complex<double> RCC6nm4um1U2nuU = CHSQq(params, g, m, aa).R;
          complex<double> RCC4nm4um6U4nuU = CHSQq(params, g, m, K).R;

          complex<double> LCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
          complex<double> LCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).L;
          complex<double> LCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
          complex<double> LCC7G4um1U2GuU = fix_sign(GLSQq(params, m, aa)).L;
          complex<double> LCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
          complex<double> LCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).L;
          complex<double> LCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
          complex<double> LCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).L;
          complex<double> RCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
          complex<double> RCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).R;
          complex<double> RCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
          complex<double> RCC7G4um1U2GuU = fix_sign(GLSQq(params, m, aa)).R;
          complex<double> RCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
          complex<double> RCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).R;
          complex<double> RCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
          complex<double> RCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).R;

          std::array<std::complex<double>, 39> coup_pow1_1_tmp = {
              (LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                   conj(RCC2nm4um6U4nuU) +
               LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC3um3nm4U3nuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G2um6U4GuU *
                    conj(LCC3um3nm4U3nuU) * conj(RCC1nm4u3Um2nuU) +
                LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC3um3nm4U3nuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC3um3nm4U3nuU) +
                RCC2um3um1g1 * RCC3G1um1Um2GuU * RCC3G2um6U4GuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G2um6U4GuU *
                    conj(LCC3um3nm4U3nuU) * conj(RCC1nm4u3Um2nuU) +
                LCC3G2um6U4GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC3G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU)),
              +(LCC3G1um1Um2GuU * RCC2um3um1g1 * RCC3G2um6U4GuU * conj(LCC3um3nm4U3nuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC3G2um6U4GuU * RCC2um3um1g1 * RCC3G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC3um3nm4U3nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC7um3nm4U1nuU) +
                RCC5um3um6g3 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(LCC8nm4u3Um2nuU) +
                RCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC7G3um6Um2GuU * RCC5um3um6g3 * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC7G4um1U2GuU * RCC5um3um6g3 * RCC7G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(LCC8nm4u3Um2nuU) +
                RCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC7G3um6Um2GuU * RCC7G4um1U2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC7um3nm4U1nuU) * conj(RCC8nm4u3Um2nuU) +
                LCC7G4um1U2GuU * RCC7G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4nm4um6U4nuU *
                    conj(LCC4um3G2U3GuU) * conj(RCC1nm4u3Um2nuU) +
                LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                    conj(LCC1nm4u3Um2nuU) * conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4nm4um6U4nuU *
                    conj(LCC4um3G2U3GuU) * conj(RCC1nm4u3Um2nuU) +
                LCC4nm4um6U4nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU *
                    conj(LCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
                RCC1u4um6g2 * RCC1um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU *
                    conj(RCC1nm4u3Um2nuU) * conj(RCC4um3G2U3GuU)),
              +(LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                    conj(LCC4um3G2U3GuU) +
                LCC4G1um1Um2GuU * RCC2um3um1g1 * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                    conj(RCC2nm4um6U4nuU) +
                LCC4nm4um6U4nuU * RCC2um3um1g1 * RCC4G1um1Um2GuU * conj(LCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU) +
                RCC2um3um1g1 * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC2nm4um6U4nuU) *
                    conj(RCC4um3G2U3GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                    conj(LCC6um3G4U1GuU) +
                LCC6G3um6Um2GuU * RCC5um3um6g3 * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                    conj(RCC5nm4um1U2nuU) +
                LCC6nm4um1U2nuU * RCC5um3um6g3 * RCC6G3um6Um2GuU * conj(LCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU) +
                RCC5um3um6g3 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC5nm4um1U2nuU) *
                    conj(RCC6um3G4U1GuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),
              +(LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(LCC8nm4u3Um2nuU) +
                LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU) +
                LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(LCC8nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
                RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                    conj(RCC6um3G4U1GuU) * conj(RCC8nm4u3Um2nuU)),

          };
          add(coup_pow1_1, coup_pow1_1_tmp);

          std::array<std::complex<double>, 50> coup_pow1_2_tmp = {
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(RCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(RCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(RCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC5um3um6g3 * conj(LCC3um3nm4U3nuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(LCC5nm4um1U2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(RCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(LCC8nm4u3Um2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC5um3um6g3 * conj(LCC5nm4um1U2nuU) *
                   conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC5um3um6g3 * conj(LCC4um3G2U3GuU) *
                   conj(RCC5nm4um1U2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC3G1um1Um2GuU * LCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(LCC8nm4u3Um2nuU) +
               RCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC3G1um1Um2GuU * RCC3G2um6U4GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC3um3nm4U3nuU) +
               LCC3G2um6U4GuU * RCC3G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC3um3nm4U3nuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(LCC8nm4u3Um2nuU) +
               RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(RCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC8nm4u3Um2nuU) * conj(RCC4um3G2U3GuU) +
               LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * RCC8u4um1g2 * RCC8um3um6g1 *
                   conj(LCC4um3G2U3GuU) * conj(RCC8nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC6um3G4U1GuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC6um3G4U1GuU) +
               RCC2um3um1g1 * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU)),
              (LCC6G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC6nm4um1U2nuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC6G3um6Um2GuU *
                   conj(LCC6um3G4U1GuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 *
                   conj(LCC1nm4u3Um2nuU) * conj(LCC7um3nm4U1nuU) +
               RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU *
                   conj(RCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC2um3um1g1 * conj(LCC2nm4um6U4nuU) *
                   conj(LCC7um3nm4U1nuU) +
               RCC2um3um1g1 * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU)),
              (LCC7G3um6Um2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G4um1U2GuU *
                   conj(LCC1nm4u3Um2nuU) * conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC1u4um6g2 * RCC1um3um1g1 * RCC7G3um6Um2GuU *
                   conj(LCC7um3nm4U1nuU) * conj(RCC1nm4u3Um2nuU)),
              (LCC7G3um6Um2GuU * RCC2um3um1g1 * RCC7G4um1U2GuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC7um3nm4U1nuU) +
               LCC7G4um1U2GuU * RCC2um3um1g1 * RCC7G3um6Um2GuU * conj(LCC7um3nm4U1nuU) *
                   conj(RCC2nm4um6U4nuU)),
              (LCC6G3um6Um2GuU * RCC2um3um1g1 * RCC6nm4um1U2nuU * conj(LCC2nm4um6U4nuU) *
                   conj(RCC6um3G4U1GuU) +
               LCC6nm4um1U2nuU * RCC2um3um1g1 * RCC6G3um6Um2GuU * conj(LCC6um3G4U1GuU) *
                   conj(RCC2nm4um6U4nuU)),

          };
          add(coup_pow1_2, coup_pow1_2_tmp);
        }
        if (has_non_zero(coup_pow1_1)) {
          coups_pow1_1_cross[aa][K].push_back(std::make_tuple(m, coup_pow1_1));
        }
        if (has_non_zero(coup_pow1_2)) {
          coups_pow1_2_cross[aa][K].push_back(std::make_tuple(m, coup_pow1_2));
        }
      }

      for (int m = 0; m < 12; m++) {
        for (int n = 0; n < 12; n++) {
          std::array<std::complex<double>, 77> coup_pow2_2 = {};
          std::array<std::complex<double>, 38> coup_pow2_3 = {};
          std::array<std::complex<double>, 36> coup_pow2_1 = {};
          for (bb = BMIN; bb < BMAX; bb++) {
            {
              std::complex<double> LCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
              std::complex<double> LCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).L;
              std::complex<double> LCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).L;
              std::complex<double> LCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).L;
              std::complex<double> LCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
              std::complex<double> LCC6um3G4U1GuU = fix_sign(GLSQq(params, n, bb)).L;
              std::complex<double> LCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).L;
              std::complex<double> LCC7G4um1U2GuU = fix_sign(GLSQq(params, n, aa)).L;

              std::complex<double> RCC3G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
              std::complex<double> RCC3G2um6U4GuU = fix_sign(GLSQq(params, m, K)).R;
              std::complex<double> RCC4G1um1Um2GuU = fix_sign(GLSQq(params, j, aa)).R;
              std::complex<double> RCC4um3G2U3GuU = fix_sign(GLSQq(params, m, bb)).R;
              std::complex<double> RCC6G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
              std::complex<double> RCC6um3G4U1GuU = fix_sign(GLSQq(params, n, bb)).R;
              std::complex<double> RCC7G3um6Um2GuU = fix_sign(GLSQq(params, j, K)).R;
              std::complex<double> RCC7G4um1U2GuU = fix_sign(GLSQq(params, n, aa)).R;

              std::complex<double> LCC7um3nm4U1nuU = CHSQq(params, g, n, bb).L;
              std::complex<double> LCC3um3nm4U3nuU = CHSQq(params, g, m, bb).L;
              std::complex<double> LCC6nm4um1U2nuU = CHSQq(params, g, n, aa).L;
              std::complex<double> LCC4nm4um6U4nuU = CHSQq(params, g, m, K).L;
              std::complex<double> RCC7um3nm4U1nuU = CHSQq(params, g, n, bb).R;
              std::complex<double> RCC3um3nm4U3nuU = CHSQq(params, g, m, bb).R;
              std::complex<double> RCC6nm4um1U2nuU = CHSQq(params, g, n, aa).R;
              std::complex<double> RCC4nm4um6U4nuU = CHSQq(params, g, m, K).R;

              std::array<std::complex<double>, 77> coup_pow2_2_tmp = {
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC6um3G4U1GuU) +
                   RCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4G1um1Um2GuU) +
                   LCC6nm4um1U2nuU * RCC4um3G2U3GuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * RCC4um3G2U3GuU * RCC6nm4um1U2nuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * RCC4um3G2U3GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(LCC6um3G4U1GuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC6um3G4U1GuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC6G3um6Um2GuU * LCC6nm4um1U2nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC6um3G4U1GuU) +
                   LCC3um3nm4U3nuU * LCC6nm4um1U2nuU * RCC6G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC6um3G4U1GuU) * conj(RCC3G2um6U4GuU) +
                   LCC6G3um6Um2GuU * RCC3um3nm4U3nuU * RCC6nm4um1U2nuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC6um3G4U1GuU) +
                   RCC3um3nm4U3nuU * RCC6G3um6Um2GuU * RCC6nm4um1U2nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC6um3G4U1GuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC4G1um1Um2GuU) * conj(RCC4nm4um6U4nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC4um3G2U3GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC7um3nm4U1nuU) +
                   LCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU) +
                   RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G1um1Um2GuU) +
                   LCC3um3nm4U3nuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC3G2um6U4GuU) +
                   LCC7G3um6Um2GuU * RCC3um3nm4U3nuU * RCC7G4um1U2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC7um3nm4U1nuU) +
                   LCC7G4um1U2GuU * RCC3um3nm4U3nuU * RCC7G3um6Um2GuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC4um3G2U3GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(LCC4nm4um6U4nuU) * conj(LCC7um3nm4U1nuU) +
                   LCC4um3G2U3GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC4nm4um6U4nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC4G1um1Um2GuU) +
                   LCC7G4um1U2GuU * RCC4um3G2U3GuU * RCC7G3um6Um2GuU * conj(LCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU) +
                   RCC4um3G2U3GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC4G1um1Um2GuU) *
                       conj(RCC4nm4um6U4nuU) * conj(RCC7um3nm4U1nuU)),
              };
              add(coup_pow2_2, coup_pow2_2_tmp);

              LCC4um3G2U3GuU = fix_sign(GLSQq(params, n, bb)).L;
              RCC4um3G2U3GuU = fix_sign(GLSQq(params, n, bb)).R;
              LCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).L;
              RCC6um3G4U1GuU = fix_sign(GLSQq(params, m, bb)).R;

              LCC6nm4um1U2nuU = CHSQq(params, g, m, aa).L;
              LCC4nm4um6U4nuU = CHSQq(params, g, n, K).L;
              RCC6nm4um1U2nuU = CHSQq(params, g, m, aa).R;
              RCC4nm4um6U4nuU = CHSQq(params, g, n, K).R;

              std::array<std::complex<double>, 38> coup_pow2_3_tmp = {
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(LCC4um3G2U3GuU) +
                   RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(RCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * LCC4nm4um6U4nuU * RCC4G1um1Um2GuU * conj(LCC3G2um6U4GuU) *
                       conj(LCC4um3G2U3GuU) * conj(RCC3G1um1Um2GuU) +
                   LCC4G1um1Um2GuU * RCC3um3nm4U3nuU * RCC4nm4um6U4nuU * conj(LCC3G1um1Um2GuU) *
                       conj(RCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC3um3nm4U3nuU * RCC4G1um1Um2GuU * RCC4nm4um6U4nuU * conj(LCC4um3G2U3GuU) *
                       conj(RCC3G1um1Um2GuU) * conj(RCC3G2um6U4GuU) +
                   LCC4G1um1Um2GuU * LCC4nm4um6U4nuU * RCC3um3nm4U3nuU * conj(LCC3G1um1Um2GuU) *
                       conj(LCC3G2um6U4GuU) * conj(RCC4um3G2U3GuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * LCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(LCC7um3nm4U1nuU) +
                   RCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(RCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G4um1U2GuU * RCC6um3G4U1GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * RCC7G3um6Um2GuU * RCC7G4um1U2GuU * conj(LCC7um3nm4U1nuU) *
                       conj(RCC6G3um6Um2GuU) * conj(RCC6nm4um1U2nuU) +
                   LCC7G3um6Um2GuU * LCC7G4um1U2GuU * RCC6um3G4U1GuU * conj(LCC6G3um6Um2GuU) *
                       conj(LCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),
                  (LCC6um3G4U1GuU * LCC7G4um1U2GuU * RCC7G3um6Um2GuU * conj(LCC6nm4um1U2nuU) *
                       conj(LCC7um3nm4U1nuU) * conj(RCC6G3um6Um2GuU) +
                   LCC7G3um6Um2GuU * RCC6um3G4U1GuU * RCC7G4um1U2GuU * conj(LCC6G3um6Um2GuU) *
                       conj(RCC6nm4um1U2nuU) * conj(RCC7um3nm4U1nuU)),

              };
              add(coup_pow2_3, coup_pow2_3_tmp);
            }

            {
              std::complex<double> LCC3G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).L);
              std::complex<double> LCC6G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).L);
              std::complex<double> LCC4G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).L);
              std::complex<double> LCC7G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).L);
              std::complex<double> RCC3G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).R);
              std::complex<double> RCC4G1um1Um2GuU = norm(fix_sign(GLSQq(params, j, aa)).R);
              std::complex<double> RCC6G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).R);
              std::complex<double> RCC7G3um6Um2GuU = norm(fix_sign(GLSQq(params, j, K)).R);

              std::complex<double> LCC3G2um6U4GuU =
                  fix_sign(GLSQq(params, n, K)).L * conj(fix_sign(GLSQq(params, m, K)).L);
              std::complex<double> LCC4um3G2U3GuU =
                  fix_sign(GLSQq(params, m, bb)).L * conj(fix_sign(GLSQq(params, n, bb)).L);
              std::complex<double> LCC6um3G4U1GuU =
                  fix_sign(GLSQq(params, m, bb)).L * conj(fix_sign(GLSQq(params, n, bb)).L);
              std::complex<double> LCC7G4um1U2GuU =
                  fix_sign(GLSQq(params, n, aa)).L * conj(fix_sign(GLSQq(params, m, aa)).L);
              std::complex<double> RCC3G2um6U4GuU =
                  fix_sign(GLSQq(params, n, K)).R * conj(fix_sign(GLSQq(params, m, K)).R);
              std::complex<double> RCC4um3G2U3GuU =
                  fix_sign(GLSQq(params, m, bb)).R * conj(fix_sign(GLSQq(params, n, bb)).R);
              std::complex<double> RCC6um3G4U1GuU =
                  fix_sign(GLSQq(params, m, bb)).R * conj(fix_sign(GLSQq(params, n, bb)).R);
              std::complex<double> RCC7G4um1U2GuU =
                  fix_sign(GLSQq(params, n, aa)).R * conj(fix_sign(GLSQq(params, m, aa)).R);

              std::complex<double> LCC3um3nm4U3nuU =
                  CHSQq(params, g, m, bb).L * conj(CHSQq(params, g, n, bb).L);
              std::complex<double> LCC4nm4um6U4nuU =
                  CHSQq(params, g, n, K).L * conj(CHSQq(params, g, m, K).L);
              std::complex<double> LCC6nm4um1U2nuU =
                  CHSQq(params, g, n, aa).L * conj(CHSQq(params, g, m, aa).L);
              std::complex<double> LCC7um3nm4U1nuU =
                  CHSQq(params, g, m, bb).L * conj(CHSQq(params, g, n, bb).L);
              std::complex<double> RCC3um3nm4U3nuU =
                  CHSQq(params, g, m, bb).R * conj(CHSQq(params, g, n, bb).R);
              std::complex<double> RCC4nm4um6U4nuU =
                  CHSQq(params, g, n, K).R * conj(CHSQq(params, g, m, K).R);
              std::complex<double> RCC6nm4um1U2nuU =
                  CHSQq(params, g, n, aa).R * conj(CHSQq(params, g, m, aa).R);
              std::complex<double> RCC7um3nm4U1nuU =
                  CHSQq(params, g, m, bb).R * conj(CHSQq(params, g, n, bb).R);

              std::array<std::complex<double>, 36> coup_pow2_1_tmp = {
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) +
                   (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) +
                   (LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) +
                   (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (LCC4um3G2U3GuU) * (RCC4nm4um6U4nuU) +
                   (LCC4nm4um6U4nuU) * (RCC4G1um1Um2GuU) * (RCC4um3G2U3GuU) +
                   (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU)),
                  ((LCC4G1um1Um2GuU) * (LCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU) * (RCC4um3G2U3GuU) +
                   (LCC4nm4um6U4nuU) * (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) +
                   (LCC4um3G2U3GuU) * (RCC4G1um1Um2GuU) * (RCC4nm4um6U4nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) +
                   (LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) +
                   (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) +
                   (LCC3G1um1Um2GuU) * (LCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) +
                   (RCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC3G1um1Um2GuU) * (LCC3um3nm4U3nuU) * (RCC3G2um6U4GuU) +
                   (LCC3G1um1Um2GuU) * (RCC3G2um6U4GuU) * (RCC3um3nm4U3nuU) +
                   (LCC3G2um6U4GuU) * (LCC3um3nm4U3nuU) * (RCC3G1um1Um2GuU) +
                   (LCC3G2um6U4GuU) * (RCC3G1um1Um2GuU) * (RCC3um3nm4U3nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (LCC6um3G4U1GuU) * (RCC6nm4um1U2nuU) +
                   (LCC6nm4um1U2nuU) * (RCC6G3um6Um2GuU) * (RCC6um3G4U1GuU) +
                   (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU)),
                  ((LCC6G3um6Um2GuU) * (LCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU) * (RCC6um3G4U1GuU) +
                   (LCC6nm4um1U2nuU) * (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) +
                   (LCC6um3G4U1GuU) * (RCC6G3um6Um2GuU) * (RCC6nm4um1U2nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) +
                   (LCC7G3um6Um2GuU) * (LCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) +
                   (RCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU)),
                  ((LCC7G3um6Um2GuU) * (LCC7um3nm4U1nuU) * (RCC7G4um1U2GuU) +
                   (LCC7G3um6Um2GuU) * (RCC7G4um1U2GuU) * (RCC7um3nm4U1nuU) +
                   (LCC7G4um1U2GuU) * (LCC7um3nm4U1nuU) * (RCC7G3um6Um2GuU) +
                   (LCC7G4um1U2GuU) * (RCC7G3um6Um2GuU) * (RCC7um3nm4U1nuU)),
              };
              add(coup_pow2_1, coup_pow2_1_tmp);
            }
          }
          if (has_non_zero(coup_pow2_1)) {
            coups_pow2_1_cross[aa][K].push_back(std::make_tuple(m, n, coup_pow2_1));
          }
          if (has_non_zero(coup_pow2_2)) {
            coups_pow2_2_cross[aa][K].push_back(std::make_tuple(m, n, coup_pow2_2));
          }
          if (has_non_zero(coup_pow2_3)) {
            coups_pow2_3_cross[aa][K].push_back(std::make_tuple(m, n, coup_pow2_3));
          }
        }
      }
      // std::cout << "Non-zero elements: " << coups_pow1_1.size() << "/"<< 12<<
      // std::endl; std::cout << "Non-zero elements: " << coups_pow1_2.size() <<
      // "/"<< 12<< std::endl; std::cout << "Non-zero elements: " <<
      // coups_pow2_1.size() << "/"<< 12*12<< std::endl; std::cout << "Non-zero
      // elements: " << coups_pow2_2.size() << "/"<< 12*12<< std::endl;
      // std::cout
      // << "Non-zero elements: " << coups_pow2_3.size() << "/"<< 12*12<<
      // std::endl;
    }
  }
}

/**
 * \returns onshell subtracted matrix element with double squark propagators of
 * the process \f$ q+q\to Q+\chi+q \f$ or \f$ \bar q+\bar q\to \bar Q+\chi+\bar
 * q \f$ depending on \p params.
 *
 * Here the internal squark going onshell is subtracted.
 * The resonance is regualted by #WIDTH.
 * The matrix element is evaluated in the helicity frame (see
 * #FI::SetKinematic_HelicityFrame) defined through \p S \p M2 \p PT2 \p TH \p
 * PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quark_gaugino_squark_uu_UXu_onshell_23(const double S, const double M2,
                                                   const double PT2, const double TH,
                                                   const double PH, int YS, Parameters *params) {
  complex<double> born = 0.0;
  static const complex<double> III(0.0, 1.0);

  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // int m=0,n=0; // internal squarks (looped)
  int g, j;

  double m1, m2;
  SET_MASS;
  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  // Set 2->3 kinematics using helicity frame.
  // JET_KINEMATICS
  ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);

  if (ii > 30) {
    ii = ii - 31;
    j = ii;
    g = jj;

    ////_DEBUG_MSG1("ii=%d", ii);

  } else {
    jj = jj - 31;
    j = jj;
    g = ii;

    //_DEBUG_MSG1("jj=%d", jj);
  }
  if (first_run_cross) {
    set_couplings_cross(params, j, g);
    first_run_cross = false;
  }

  double CSUM, alpha_s;
  alpha_s = 1.0;
  CSUM = 1.0;
  bool crossed = true;
  double g3s = std::norm(params->gqq[0][0].R);
  double Cf = cF;
  double Nc = NC;
  // std::cout << "a " << aa << " b " << bb << std::endl;
  for (const std::tuple<int, int, std::array<std::complex<double>, 77>> it :
       coups_pow2_2_cross[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 77> coupling;
    std::tie(m, n, coupling) = it;
    // std::cout << " Off diag: " << n << " " << m << std::endl;
    ff->Cross_pb_p3();
    {
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
      // doubled (since we iterate all m and n).
      double f = 1.0;
      if (m != n)
        f = 0.5;
      born +=
          f * 1. / g3s * 1. / g3s *
          M_tot_uu_UXu_pow2_2_new_cross_onshell_qmp1mk2MUij(
              coupling, cF, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, QK2, alpha_s,
              dM2o, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_qmk2MG_pow1, rDn_qmp1mk2MUi2_pow1,
              rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
    }
    ff->Cross_pb_p3();
    double MQi1 = params->mSQ[m];
    double MQj1 = MQi1;
    double MQi2 = params->mSQ[n];
#ifdef ONSUB
    if (m == n && 2.0 * ff->papb >= pow2(MQi2 + ff->m1) && MQi2 > ff->m2) {
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(MQi2), MQi2 * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES23(pow2(MQi2), factor);
      ff->Cross_pb_p3();
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born -=
          ff->BreitWigner * factor * 1. / g3s * 1. / g3s *
          M_tot_uu_UXu_pow2_2_new_cross_onshell_qmp1mk2MUij(
              coupling, cF, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, QK2, alpha_s,
              dM2o, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_qmk2MG_pow1, rDn_qmp1mk2MUi2_pow1,
              rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
      ff->Cross_pb_p3();
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }
  for (const std::tuple<int, int, std::array<std::complex<double>, 36>> it :
       coups_pow2_1_cross[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 36> coupling;
    std::tie(m, n, coupling) = it;
    // std::cout << " diag: " << n << " " << m << std::endl;
    ff->Cross_pb_p3();
    {
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      double f = 1.0;
      // if(m!=n)f=2.0;
      born += f * 1. / g3s * 1. / g3s *
              M_tot_uu_UXu_pow2_1_new_cross_onshell_qmp1mk2MUij(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
    }
    ff->Cross_pb_p3();
    double MQi1 = params->mSQ[m];
    double MQj1 = MQi1;
    double MQi2 = params->mSQ[n];
#ifdef ONSUB
    if (m == n && 2.0 * ff->papb >= pow2(MQi2 + ff->m1) && MQi2 > ff->m2) {
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(MQi2), MQi2 * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES23(pow2(MQi2), factor);
      ff->Cross_pb_p3();
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born -= ff->BreitWigner * factor * 1. / g3s * 1. / g3s *
              M_tot_uu_UXu_pow2_1_new_cross_onshell_qmp1mk2MUij(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
      ff->Cross_pb_p3();
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }

  delete ff;
  return real(born);
}
/**
 * \returns matrix element without the double squark propagators of the process
 * \f$ q+q\to Q+\chi+q \f$ or \f$ \bar q+\bar q\to \bar Q+\chi+\bar q \f$
 * depending on \p params.
 *
 * The matrix element is evaluated in the jet kinematic frame (see
 * #FI::SetKinematic) defined through \p S \p M2 \p PT2 \p TH \p PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quark_gaugino_squark_uu_UXu(const double S, const double M2, const double PT2,
                                        const double TH, const double PH, const int YS,
                                        Parameters *params) {
  complex<double> born = 0.0;
  static const complex<double> III(0.0, 1.0);

  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // int m=0,n=0; // internal squarks (looped)
  int g, j;

  FI *ff = new FI();
  JET_KINEMATICS

  if (ii > 30) {
    ii = ii - 31;
    j = ii;
    g = jj;

    //_DEBUG_MSG1("ii=%d", ii);

  } else {
    jj = jj - 31;
    j = jj;
    g = ii;

    //_DEBUG_MSG1("jj=%d", jj);
  }
  if (first_run_cross) {
    set_couplings_cross(params, j, g);
    first_run_cross = false;
  }

  double CSUM, alpha_s;
  alpha_s = 1.0;
  CSUM = 1.0;
  bool crossed = true;
  double g3s = std::norm(params->gqq[0][0].R);
  double Cf = cF;
  double Nc = NC;
  ff->Cross_pb_p3();
  ZERO_SQUARK_DENOMS;

  for (const std::tuple<std::array<std::complex<double>, 9>> it : coups_pow0_1_cross[aa][bb]) {
    std::array<std::complex<double>, 9> coupling;
    std::tie(coupling) = it;
    born += M_tot_uu_UXu_pow0_1_new_cross(coupling, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2,
                                          P1K3, alpha_s, dM2o, rDn_k2pk3MU_pow1, rDn_k2pk3MU_pow2,
                                          s, s1, k2pk3MU);
  }
  for (const std::tuple<std::array<std::complex<double>, 6>> it : coups_pow0_2_cross[aa][bb]) {
    std::array<std::complex<double>, 6> coupling;
    std::tie(coupling) = it;
    born +=
        M_tot_uu_UXu_pow0_2_new_cross(coupling, cF, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2,
                                      P1K3, alpha_s, dM2o, rDn_k2pk3MU_pow1, s, s1, k2pk3MU);
  }
  for (const std::tuple<int, std::array<std::complex<double>, 39>> it :
       coups_pow1_1_cross[aa][bb]) {
    int m;
    std::array<std::complex<double>, 39> coupling;
    std::tie(m, coupling) = it;
    SINGLE_SQUARK_DENOMS;
    born += 1 / g3s *
            M_tot_uu_UXu_pow1_1_new_cross(coupling, cF, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3,
                                          MG, MQ, MQi1, MX, Nc, P1K2, P1K3, QK2, alpha_s, dM2o,
                                          rDn_k2pk3MU_pow1, rDn_k2pk3MUi1_pow1, rDn_p1mk2MUi1_pow1,
                                          rDn_qmk2MG_pow1, rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MU,
                                          qmk2MG, qmp1mk2MUi1, k2pk3MUi1);
  }
  for (const std::tuple<int, std::array<std::complex<double>, 50>> it :
       coups_pow1_2_cross[aa][bb]) {
    int m;
    std::array<std::complex<double>, 50> coupling;
    std::tie(m, coupling) = it;
    SINGLE_SQUARK_DENOMS;
    born += 1 / g3s *
            M_tot_uu_UXu_pow1_2_new_cross(coupling, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MG,
                                          MQ, MQi1, MX, Nc, P1K2, P1K3, QK2, alpha_s, dM2o,
                                          rDn_k2pk3MU_pow1, rDn_k2pk3MUi1_pow1, rDn_p1mk2MUi1_pow1,
                                          rDn_qmk2MG_pow1, rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MUi1,
                                          qmk2MG, k2pk3MU, qmp1mk2MUi1);
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 77>> it :
       coups_pow2_2_cross[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 77> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    double f = 1.0;
    // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
    // doubled (since we iterate all m and n). Squared diagrams.
    if (m != n)
      f = 0.5;
    born += f * 1 / g3s * 1 / g3s *
            M_tot_uu_UXu_pow2_2_new_cross(coupling, cF, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX,
                                          Nc, P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MUj1_pow1,
                                          rDn_p1mk2MUi2_pow1, rDn_qmk2MG_pow1, rDn_qmp1mk2MUi2_pow1,
                                          rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2, qmp1mk2MUj1,
                                          k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }
  for (const std::tuple<int, int, std::array<std::complex<double>, 38>> it :
       coups_pow2_3_cross[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 38> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    double f = 1.0;
    // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
    // doubled (since we iterate all m and n). Squared diagrams.
    if (m != n)
      f = 0.5;
    born += f * 1 / g3s * 1 / g3s *
            M_tot_uu_UXu_pow2_3_new_cross(coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc,
                                          P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MUi2_pow1,
                                          rDn_p1mk2MUj1_pow1, rDn_qmk2MG_pow1, rDn_qmk2MG_pow2,
                                          rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1, s,
                                          qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 36>> it :
       coups_pow2_1_cross[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 36> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    double f = 1.0;
    // if(m!=n)f=2.0;
    born += f * 1 / g3s * 1 / g3s *
            M_tot_uu_UXu_pow2_1_new_cross(
                coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }
  ff->Cross_pb_p3();

  delete ff;
  return real(born);
}
/**
 * \returns matrix element with the double gluino propagators of the process
 * \f$ q+\bar q\to Q+\chi+\bar q \f$ or \f$ q+\bar q\to \bar Q+\chi+ q \f$
 * depending on \p params.
 *
 * Here the internal gluino going onshell is subtracted.
 * The resonance is regualted by #WIDTH.
 * The matrix element is evaluated in the helicity frame (see
 * #FI::SetKinematic_HelicityFrame) defined through \p S \p M2 \p PT2 \p TH \p
 * PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quarkb_gaugino_squark_uub_UXub_onshell_13(const double S, const double M2,
                                                      const double PT2, const double TH,
                                                      const double PH, int YS, Parameters *params) {
  complex<double> born = 0.0;
  static const complex<double> III(0.0, 1.0);

  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // int m=0,n=0; // internal squarks (looped)
  int g, j;

  double m1, m2;
  SET_MASS;
  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  // Set 2->3 kinematics using helicity frame.
  ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);

  if (ii > 30) {
    ii = ii - 31;
    j = ii;
    g = jj;

    //_DEBUG_MSG1("ii=%d", ii);

  } else {
    jj = jj - 31;
    j = jj;
    g = ii;

    //_DEBUG_MSG1("jj=%d", jj);
  }
  if (first_run) {
    set_couplings(params, j, g);
    first_run = false;
  }

  double CSUM, alpha_s;
  alpha_s = 1.0;
  CSUM = 1.0;
  bool crossed = false;
  double g3s = std::norm(params->gqq[0][0].R);
  double Cf = cF;
  double Nc = NC;

  for (const std::tuple<int, int, std::array<std::complex<double>, 38>> it : coups_pow2_3[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 38> coupling;
    std::tie(m, n, coupling) = it;
    {
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
      // doubled (since we iterate all m and n). Squared diagrams.
      double f = 1.0;
      if (m != n)
        f = 0.5;
      born += f * 1. / g3s * 1. / g3s *
              M_tot_uu_UXu_pow2_3_new_onshell_qmk2MG(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, QK2, alpha_s,
                  dM2o, rDn_k2pk3MUi2_pow1, rDn_p1mk2MUj1_pow1, rDn_qmk2MG_pow1, rDn_qmk2MG_pow2,
                  rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2, qmp1mk2MUj1,
                  k2pk3MUi2, k2pk3MUj1, qmk2MG);
    }
    double MQi1 = params->mSQ[m];
    double MQj1 = MQi1;
    double MQi2 = params->mSQ[n];
// if(std::isnan(real(-born)))std::cout << "pre1 " << std::endl;
#ifdef ONSUB
    if (2.0 * ff->papb >= pow2(params->mGL + ff->m2) && params->mGL > (ff->m1)) {
      // std::cout << params->mGL << std::endl;
      ff->SetBreitWigner(2.0 * ff->p1p3 + ff->m1s, pow2(params->mGL), params->mGL * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES13(pow2(params->mGL), factor);
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      double f = 1.0;
      // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
      // doubled (since we iterate all m and n). Squared diagrams.
      if (m != n)
        f = 0.5;
      born -= f * ff->BreitWigner * factor * 1.0 / g3s * 1.0 / g3s *
              M_tot_uu_UXu_pow2_3_new_onshell_qmk2MG(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, QK2, alpha_s,
                  dM2o, rDn_k2pk3MUi2_pow1, rDn_p1mk2MUj1_pow1, rDn_qmk2MG_pow1, rDn_qmk2MG_pow2,
                  rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2, qmp1mk2MUj1,
                  k2pk3MUi2, k2pk3MUj1, qmk2MG);
      // if(std::isnan(real(-born)))std::cout << "1 " << std::endl;
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 36>> it : coups_pow2_1[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 36> coupling;
    std::tie(m, n, coupling) = it;
    {
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born += 1.0 / g3s * 1.0 / g3s *
              M_tot_uu_UXu_pow2_1_new_onshell_qmk2MG(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
      // if(std::isnan(real(-born)))std::cout <<"pre2 " <<
      // (Dn_qmk2MG*Dn_qmk2MG)/(Dn_qmk2MG*Dn_qmk2MG+MG*MG*MG*MG*WIDTH*WIDTH) <<
      // std::endl;
    }
    double MQi1 = params->mSQ[m];
    double MQj1 = MQi1;
    double MQi2 = params->mSQ[n];
#ifdef ONSUB
    if (2.0 * ff->papb >= pow2(params->mGL + ff->m2) && params->mGL > (ff->m1)) {

      ff->SetBreitWigner(2.0 * ff->p1p3 + ff->m1s, pow2(params->mGL), params->mGL * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES13(pow2(params->mGL), factor);
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born -= ff->BreitWigner * factor * 1.0 / g3s * 1.0 / g3s *
              M_tot_uu_UXu_pow2_1_new_onshell_qmk2MG(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
      // if(std::isnan(real(-born))) {
      // std::cout <<"2 " << M_tot_uu_UXu_pow2_1_new_all(coupling, Dn_p1mk1_MG,
      // K2K3,  MG,  MQ,  MQi2,  MQj1,  MX,  Nc,  P1K2,  P1K3,  alpha_s,
      // rDn_k2pk3MUi2_pow1,  rDn_k2pk3MUj1_pow1,  rDn_p1mk2MUi2_pow1,
      // rDn_p1mk2MUj1_pow1,  rDn_qmk2MG_pow1,  rDn_qmk2MG_pow2,
      // rDn_qmp1mk2MUi2_pow1,  rDn_qmp1mk2MUj1_pow1) << " " << born << " " <<
      // factor << " " << ff->BreitWigner << " " << Dn_p1mk1_MG << " "<< K2K3 <<
      // " "<< MG  << " "<< MQ << " "<< MQi2<< " "<<  MQj1<< " "<<  MX<< " "<<
      // Nc<< " "<<   P1K2<< " "<<  P1K3<< " "<<  alpha_s<< " "<<
      // rDn_k2pk3MUi2_pow1<< " "<<  rDn_k2pk3MUj1_pow1<< " "<<
      // rDn_p1mk2MUi2_pow1<< " "<<  rDn_p1mk2MUj1_pow1<< "  !"<<
      // rDn_qmk2MG_pow1<< "! "<<  rDn_qmk2MG_pow2 << " "<< rDn_qmp1mk2MUi2_pow1
      // << " "<<  rDn_qmp1mk2MUj1_pow1 << " "<<  std::endl;
      //}
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }
  // if(!(real(born) < 10e18))std::cout << "WTF " << real(born)<< std::endl;

  delete ff;
  return real(born);
}
/**
 * \returns matrix element with the double squark propagators of the process
 * \f$ q+\bar q\to Q+\chi+\bar q \f$ or \f$ q+\bar q\to \bar Q+\chi+ q \f$
 * depending on \p params.
 *
 * Here the internal squark going onshell is subtracted.
 * The resonance is regualted by #WIDTH.
 * The matrix element is evaluated in the helicity frame (see
 * #FI::SetKinematic_HelicityFrame) defined through \p S \p M2 \p PT2 \p TH \p
 * PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quarkb_gaugino_squark_uub_UXub_onshell_23(const double S, const double M2,
                                                      const double PT2, const double TH,
                                                      const double PH, int YS, Parameters *params) {
  complex<double> born = 0.0;
  static const complex<double> III(0.0, 1.0);

  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // int m=0,n=0; // internal squarks (looped)
  int g, j;

  double m1, m2;
  SET_MASS;
  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  // Set 2->3 kinematics using helicity frame.
  ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);

  if (ii > 30) {
    ii = ii - 31;
    j = ii;
    g = jj;

    //_DEBUG_MSG1("ii=%d", ii);

  } else {
    jj = jj - 31;
    j = jj;
    g = ii;

    //_DEBUG_MSG1("jj=%d", jj);
  }
  if (first_run) {
    set_couplings(params, j, g);
    first_run = false;
  }

  double CSUM, alpha_s;
  alpha_s = 1.0;
  CSUM = 1.0;
  bool crossed = false;
  double g3s = std::norm(params->gqq[0][0].R);
  double Cf = cF;
  double Nc = NC;
  for (const std::tuple<std::array<std::complex<double>, 9>> it : coups_pow0_1[aa][bb]) {
    std::array<std::complex<double>, 9> coupling;
    std::tie(coupling) = it;
    {
      ZERO_SQUARK_DENOMS;
      born += M_tot_uu_UXu_pow0_1_new_onshell_k2pk3MU(
          coupling, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2, P1K3, alpha_s, dM2o,
          rDn_k2pk3MU_pow1, rDn_k2pk3MU_pow2, s, s1, k2pk3MU);
    }
#ifdef ONSUB
    if (2.0 * ff->papb >= pow2(ff->m1 + ff->m1) && ff->m1 > ff->m2) {
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(ff->m1), ff->m1 * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES23(pow2(ff->m1), factor);
      ZERO_SQUARK_DENOMS;
      born -= ff->BreitWigner * factor *
              M_tot_uu_UXu_pow0_1_new_onshell_k2pk3MU(
                  coupling, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2, P1K3, alpha_s, dM2o,
                  rDn_k2pk3MU_pow1, rDn_k2pk3MU_pow2, s, s1, k2pk3MU);
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }
  for (const std::tuple<int, std::array<std::complex<double>, 39>> it : coups_pow1_1[aa][bb]) {
    int m;
    std::array<std::complex<double>, 39> coupling;
    std::tie(m, coupling) = it;
    {
      ZERO_SQUARK_DENOMS;
      SINGLE_SQUARK_DENOMS;
      born += 1 / g3s *
              M_tot_uu_UXu_pow1_1_new_onshell_k2pk3MUi(
                  coupling, cF, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MG, MQ, MQi1, MX, Nc,
                  P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MU_pow1, rDn_k2pk3MUi1_pow1,
                  rDn_p1mk2MUi1_pow1, rDn_qmk2MG_pow1, rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MU, qmk2MG,
                  qmp1mk2MUi1, k2pk3MUi1);
    }
    double MQi1 = params->mSQ[m];
#ifdef ONSUB
    if (m == j && 2.0 * ff->papb >= pow2(ff->m1 + ff->m1) && ff->m1 > ff->m2) {
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(ff->m1), ff->m1 * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES23(pow2(ff->m1), factor);
      ZERO_SQUARK_DENOMS;
      SINGLE_SQUARK_DENOMS;
      born -= ff->BreitWigner * factor * 1 / g3s *
              M_tot_uu_UXu_pow1_1_new_onshell_k2pk3MUi(
                  coupling, cF, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MG, MQ, MQi1, MX, Nc,
                  P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MU_pow1, rDn_k2pk3MUi1_pow1,
                  rDn_p1mk2MUi1_pow1, rDn_qmk2MG_pow1, rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MU, qmk2MG,
                  qmp1mk2MUi1, k2pk3MUi1);
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 36>> it : coups_pow2_1[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 36> coupling;
    std::tie(m, n, coupling) = it;
    {
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born += 1 / g3s * 1 / g3s *
              M_tot_uu_UXu_pow2_1_new_onshell_k2pk3MUij(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
    }
    double MQi1 = params->mSQ[m];
    double MQj1 = MQi1;
    double MQi2 = params->mSQ[n];
#ifdef ONSUB
    if (m == n && 2.0 * ff->papb >= pow2(MQi2 + ff->m1) && MQi2 > ff->m2) {
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(MQi2), MQi2 * WIDTH);
      double factor = 0.0;
      ff->SetKinematicRES23(pow2(MQi2), factor);
      ZERO_SQUARK_DENOMS;
      DOUBLE_SQUARK_DENOMS;
      born -= ff->BreitWigner * factor * 1 / g3s * 1 / g3s *
              M_tot_uu_UXu_pow2_1_new_onshell_k2pk3MUij(
                  coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2, P1K3, alpha_s,
                  rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1,
                  rDn_qmk2MG_pow1, rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                  qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
      ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }

  delete ff;
  return real(born);
}
/**
 * \returns matrix element without the double squark and double gluino
 * propagators of the process \f$ q+\bar q\to Q+\chi+\bar q \f$ or \f$ q+\bar
 * q\to \bar Q+\chi+ q \f$ depending on \p params.
 *
 * The matrix element is evaluated in the jet kinematic frame (see
 * #FI::SetKinematic) defined through \p S \p M2 \p PT2 \p TH \p PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quarkb_gaugino_squark_uub_UXub(const double S, const double M2, const double PT2,
                                           const double TH, const double PH, const int YS,
                                           Parameters *params) {
  complex<double> born = 0.0;
  static const complex<double> III(0.0, 1.0);

  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino

  // int m=0,n=0; // internal squarks (looped)
  int g, j;

  FI *ff = new FI();
  JET_KINEMATICS

  if (ii > 30) {
    ii = ii - 31;
    j = ii;
    g = jj;

    //_DEBUG_MSG1("ii=%d", ii);

  } else {
    jj = jj - 31;
    j = jj;
    g = ii;

    //_DEBUG_MSG1("jj=%d", jj);
  }
  if (first_run) {
    set_couplings(params, j, g);
    first_run = false;
  }

  double g3s = std::norm(params->gqq[0][0].R);
  double Cf = cF;
  bool crossed = false;
  double alpha_s = 1.0;
  double CSUM = 1.0;
  double Nc = NC;

  ZERO_SQUARK_DENOMS;

  _DEBUG_TABLE("m1s", ff->m1s);
  _DEBUG_TABLE("m2s", ff->m2s);
  _DEBUG_TABLE("MG", MG);
  _DEBUG_TABLE("K1Q", K1Q);
  _DEBUG_TABLE("K1P1", ff->pap1);
  _DEBUG_TABLE("K2P1", ff->pap2);
  _DEBUG_TABLE("K3P1", ff->pap3);
  _DEBUG_TABLE("K3P2", ff->pbp3);
  _DEBUG_TABLE("P1K3", ff->pap3);
  _DEBUG_TABLE("P1K2", ff->pap2);
  _DEBUG_TABLE("K2K3", K2K3);

  for (const std::tuple<std::array<std::complex<double>, 9>> it : coups_pow0_1[aa][bb]) {
    std::array<std::complex<double>, 9> coupling;
    std::tie(coupling) = it;
    born +=
        M_tot_uu_UXu_pow0_1_new(coupling, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2, P1K3,
                                alpha_s, dM2o, rDn_k2pk3MU_pow1, rDn_k2pk3MU_pow2, s, s1, k2pk3MU);
  }
  for (const std::tuple<std::array<std::complex<double>, 6>> it : coups_pow0_2[aa][bb]) {
    std::array<std::complex<double>, 6> coupling;
    std::tie(coupling) = it;
    born += M_tot_uu_UXu_pow0_2_new(coupling, cF, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MQ, MX, Nc, P1K2,
                                    P1K3, alpha_s, dM2o, rDn_k2pk3MU_pow1, s, s1, k2pk3MU);
  }
  for (const std::tuple<int, std::array<std::complex<double>, 39>> it : coups_pow1_1[aa][bb]) {
    int m;
    std::array<std::complex<double>, 39> coupling;
    std::tie(m, coupling) = it;
    SINGLE_SQUARK_DENOMS;
    born += 1 / g3s *
            M_tot_uu_UXu_pow1_1_new(coupling, cF, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MG, MQ,
                                    MQi1, MX, Nc, P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MU_pow1,
                                    rDn_k2pk3MUi1_pow1, rDn_p1mk2MUi1_pow1, rDn_qmk2MG_pow1,
                                    rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MU, qmk2MG, qmp1mk2MUi1,
                                    k2pk3MUi1);
  }
  for (const std::tuple<int, std::array<std::complex<double>, 50>> it : coups_pow1_2[aa][bb]) {
    int m;
    std::array<std::complex<double>, 50> coupling;
    std::tie(m, coupling) = it;
    SINGLE_SQUARK_DENOMS;
    born += 1 / g3s *
            M_tot_uu_UXu_pow1_2_new(coupling, Dn_p1mk1_MG, Dn_p1mk2_MQ, Dn_p2mk3, K2K3, MG, MQ,
                                    MQi1, MX, Nc, P1K2, P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MU_pow1,
                                    rDn_k2pk3MUi1_pow1, rDn_p1mk2MUi1_pow1, rDn_qmk2MG_pow1,
                                    rDn_qmp1mk2MUi1_pow1, s, s1, k2pk3MUi1, qmk2MG, k2pk3MU,
                                    qmp1mk2MUi1);
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 77>> it : coups_pow2_2[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 77> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    double f = 1.0;
    // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
    // doubled (since we iterate all m and n). Squared diagrams.
    if (m != n)
      f = 0.5;
    born +=
        f * 1 / g3s * 1 / g3s *
        M_tot_uu_UXu_pow2_2_new(coupling, cF, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2,
                                P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MUj1_pow1, rDn_p1mk2MUi2_pow1,
                                rDn_qmk2MG_pow1, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1, s,
                                qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }
  for (const std::tuple<int, int, std::array<std::complex<double>, 38>> it : coups_pow2_3[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 38> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    double f = 1.0;
    // here the factor 1/2 is needed to adjust for diagrams that shouldn't be
    // doubled (since we iterate all m and n). Squared diagrams.
    if (m != n)
      f = 0.5;
    born += f * 1 / g3s * 1 / g3s *
            M_tot_uu_UXu_pow2_3_new(coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2,
                                    P1K3, QK2, alpha_s, dM2o, rDn_k2pk3MUi2_pow1,
                                    rDn_p1mk2MUj1_pow1, rDn_qmk2MG_pow1, rDn_qmk2MG_pow2,
                                    rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1, s, qmp1mk2MUi2,
                                    qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }

  for (const std::tuple<int, int, std::array<std::complex<double>, 36>> it : coups_pow2_1[aa][bb]) {
    int m, n;
    std::array<std::complex<double>, 36> coupling;
    std::tie(m, n, coupling) = it;
    DOUBLE_SQUARK_DENOMS;
    born += 1 / g3s * 1 / g3s *
            M_tot_uu_UXu_pow2_1_new(coupling, Dn_p1mk1_MG, K2K3, MG, MQ, MQi2, MQj1, MX, Nc, P1K2,
                                    P1K3, alpha_s, rDn_k2pk3MUi2_pow1, rDn_k2pk3MUj1_pow1,
                                    rDn_p1mk2MUi2_pow1, rDn_p1mk2MUj1_pow1, rDn_qmk2MG_pow1,
                                    rDn_qmk2MG_pow2, rDn_qmp1mk2MUi2_pow1, rDn_qmp1mk2MUj1_pow1,
                                    qmp1mk2MUi2, qmp1mk2MUj1, k2pk3MUi2, k2pk3MUj1, qmk2MG);
  }

  delete ff;
  return real(born);
}
/**
 * \returns matrix element with the double squark propagators of the process
 * \f$ g+g\to Q+\chi+\bar q \f$ or \f$ g+g\to \bar Q+\chi+ q \f$ depending on \p
 * params.
 *
 * Here the internal squark going onshell is subtracted.
 * The resonance is regualted by #WIDTH.
 * The matrix element is evaluated in the helicity frame (see
 * #FI::SetKinematic_HelicityFrame) defined through \p S \p M2 \p PT2 \p TH \p
 * PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_quarkb_gaugino_squark_onshell_23(const double S, const double M2, const double PT2,
                                             const double TH, const double PH, int YS,
                                             Parameters *params) {
  double born = 0.0;
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  double m1, m2;
  SET_MASS;
  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  // Set 2->3 kinematics using helicity frame.
  // JET_KINEMATICS
  ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
  // loop over outgoing quark types
  // for(int ab = 0; ab < 5; ab++)
  {
    // TODO check vs newer version
    struct Coupling Cw[4] = {0, 0, 0, 0};
    if (ii > 30) {
      ii = ii - 31;
      // ff->SetKinematic_HelicityFrame(params->mSQ[ii], params->mCH[jj], S, M2,
      // PT2, TH, PH, YS);
      Cw[0] = params->CHSQq[jj][ii][aa];
      Cw[1] = params->gSQSQ[ii][ii];
      Cw[2] = params->CHSQq[jj][ii][aa];
      Cw[3] = params->gSQSQ[ii][ii];
      //_DEBUG_MSG1("ii=%d", ii)
    } else {
      jj = jj - 31;
      // ff->SetKinematic_HelicityFrame(params->mSQ[jj], params->mCH[ii], S, M2,
      // PT2, TH, PH, YS);
      Cw[0] = params->CHSQq[ii][jj][bb];
      Cw[1] = params->gSQSQ[jj][jj];
      Cw[2] = params->CHSQq[ii][jj][bb];
      Cw[3] = params->gSQSQ[jj][jj];
      //_DEBUG_MSG1("jj=%d", jj)
    }
    ff->SetWCoupling(Cw);
    born += ff->MG_SQGA_gg_onshell_23(params);
    double factor = 0.0;
#ifdef ONSUB
    if (2.0 * ff->papb >= pow2(ff->m1 + ff->m1) && ff->m1 > ff->m2) {
      // ps == sij
      ff->SetBreitWigner(2.0 * ff->p2p3 + ff->m2s, pow2(ff->m1), ff->m1 * WIDTH);

      ff->SetKinematicRES23(pow2(ff->m1), factor);
      // cout << factor << endl;
      // return factor;
      // born += ff->MG_SQGA_gg_onshell_23();
      // born += ff->MG_SQGA_gg_onshell_23()*factor;

      born -= ff->BreitWigner * ff->MG_SQGA_gg_onshell_23(params) * factor;

      // cout << m1 << " " << m2 << endl;

      // Reset to ordinary kinematics.
      // ff->SetKinematic_HelicityFrame(m1, m2, S, M2, PT2, TH, PH, YS);
    }
#endif
  }

  delete ff;
  return born;
}
/**
 * \returns matrix element without the double squark propagators of the process
 * \f$ g+g\to Q+\chi+\bar q \f$ or \f$ g+g\to \bar Q+\chi+ q \f$ depending on \p
 * params.
 *
 * The matrix element is evaluated in the jet kinematic frame (see
 * #FI::SetKinematic) defined through \p S \p M2 \p PT2 \p TH \p PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_gluon_gaugino_squark_gg(const double S, const double M2, const double PT2,
                                    const double TH, const double PH, const int YS,
                                    Parameters *params) {
  double born = 0.0;
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  JET_KINEMATICS
  // outgoing quark loop
  // for( int ab = 0; ab < 5; ab++) {

  // TODO check vs newer version
  struct Coupling Cw[4] = {0, 0, 0, 0};
  if (ii > 30) {
    ii = ii - 31;
    Cw[0] = params->CHSQq[jj][ii][aa];
    Cw[1] = params->gSQSQ[ii][ii];
    Cw[2] = params->CHSQq[jj][ii][aa];
    Cw[3] = params->gSQSQ[ii][ii];
    //_DEBUG_MSG1("ii=%d", ii)
  } else {
    jj = jj - 31;
    Cw[0] = params->CHSQq[ii][jj][bb];
    Cw[1] = params->gSQSQ[jj][jj];
    Cw[2] = params->CHSQq[ii][jj][bb];
    Cw[3] = params->gSQSQ[jj][jj];
    //_DEBUG_MSG1("jj=%d", jj)
  }

  ff->SetWCoupling(Cw);
  born += ff->MG_SQGA_gg(params);
  //}
  delete ff;
  return born;
}
/**
 * \returns matrix element of the process
 * \f$ q+g\to Q+\chi+g\f$ or \f$ \bar q+g\to \bar Q+\chi+ g \f$ depending on \p
 * params.
 *
 * The matrix element is evaluated in the jet kinematic frame (see
 * #FI::SetKinematic) defined through \p S \p M2 \p PT2 \p TH \p PH \p YS.
 *
 * \param S partonic com energy
 * \param M2 intermediate mass/com energy
 * \param PT2 transverse momentum squared
 * \param TH theta angle
 * \param PH phi angle
 * \param YS symmetrization flag for the angles
 * \param params Parameters of current computation/process
 *
 * \note some missing prefactors can be found in hxs.cc
 */
double real_gluon_gaugino_squark_qg(const double S, const double M2, const double PT2,
                                    const double TH, const double PH, const int YS,
                                    Parameters *params) {
  double born = 0.0;
  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // Set different kinematics for 2 -> 3 process
  FI *ff = new FI();
  JET_KINEMATICS

  // TODO check vs newer version
  struct Coupling Cw[4] = {0, 0, 0, 0};
  if (ii > 30) {
    ii = ii - 31;
    Cw[0] = params->CHSQq[jj][ii][aa];
    Cw[1] = params->gSQSQ[ii][ii];
    Cw[2] = params->CHSQq[jj][ii][aa];
    Cw[3] = params->gSQSQ[ii][ii];
    //_DEBUG_MSG1("ii=%d", ii)
  } else {
    jj = jj - 31;
    Cw[0] = params->CHSQq[ii][jj][bb];
    Cw[1] = params->gSQSQ[jj][jj];
    Cw[2] = params->CHSQq[ii][jj][bb];
    Cw[3] = params->gSQSQ[jj][jj];
    //_DEBUG_MSG1("jj=%d", jj)
  }

  //_DEBUG_MSG("aa=%d,bb=%d,ii=%d,jj=%d", aa, bb, ii, jj);
  //_DEBUG_MSG("cw[0]=(%f,%f)", Cw[0].L, Cw[0].R);
  //_DEBUG_MSG("cw[1]=(%f,%f)", Cw[1].L, Cw[1].R);
  //_DEBUG_MSG("cw[2]=(%f,%f)", Cw[2].L, Cw[2].R);
  //_DEBUG_MSG("cw[3]=(%f,%f)", Cw[3].L, Cw[3].R);
  ff->SetWCoupling(Cw);
  born += ff->MG_SQGA_qg();
  //_DEBUG_MSG("born=%f", born);

  delete ff;
  return born;
}

/**
 *  Different dipoles dsigma^A for real emission of gluon, quarks and
 * antiquarks.
 */
double Dip_SQGA(DipoleType Emitter_Spectator, PartonType emitter, PartonType spectator,
                PartonType emitted_parton, const double S, const double M2, const double PT2,
                const double TH, const double PH, const int YS, Parameters *params,
                const int Opt_eo) {

  double born = 0.0;
  vector<double> born_pol;
  double x = 0.0;
  double zi = 0.0;
  double zj = 0.0;
  double zplus = 0.0;
  double zminus = 0.0;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;

  // Set different kinematics for 2 -> 3 processes.
  FI *ff = new FI();

  JET_KINEMATICS

  double Dipole;

  switch (Emitter_Spectator) {
  case INITIAL_INITIAL:
    if ((emitter == PARTON_GLUON) && spectator == PARTON_GLUON &&
        emitted_parton == PARTON_ANTIQUARK) { // D_ga_g
      //         cout << "This is dipole 4" << endl;
      // Set dipole kinematics for emitter A and spectator B.
      if (Opt_eo == 0) {
        ff->SetDipKinematicAB(x);
        // Dipole factor.
        Dipole = D_ai_b(PARTON_GLUON, PARTON_ANTIQUARK, ff->pap3, ff->papb, ff->pbp3, x, 1,
                        -0.5 * cA / cF * TR);
      } else if (Opt_eo == 1) {
        ff->SetDipKinematicBA(x);
        // Dipole factor.
        Dipole = D_ai_b(PARTON_GLUON, PARTON_ANTIQUARK, ff->pbp3, ff->papb, ff->pap3, x, 1,
                        -0.5 * cA / cF * TR); // TODO check pap3 at end or pbp3
      }
      // std::cout << Dipole << std::endl;
    }
    if ((emitter == PARTON_QUARK || emitter == PARTON_ANTIQUARK) && spectator == PARTON_GLUON &&
        emitted_parton == PARTON_GLUON) { // D_qg_g
      //         cout << "This is dipole 4" << endl;
      // Set dipole kinematics for emitter A and spectator B.
      ff->SetDipKinematicAB(x);
      // Dipole factor.
      Dipole = D_ai_b(PARTON_QUARK, PARTON_GLUON, ff->pap3, ff->papb, ff->pbp3, x, 1, -0.5 * cA);

    } else if ((emitter == PARTON_GLUON && emitted_parton == PARTON_GLUON) &&
               (spectator == PARTON_QUARK || spectator == PARTON_ANTIQUARK)) { // D_gg_q
      // ff->papb = -ff->papb;

      double papi = 0.0;
      double s1 = 0.0;
      double ztb = 0.0;
      double PbKi = 0.0;
      double PbKt1 = 0.0;

      // Set dipole kinematics for emitter B and spectator A for polarized Bon
      // cross section contraction
      ff->SetDipKinematicAB_pol(x, s1, papi, ztb, PbKi, PbKt1);

      vector<vector<double>> DIPCONNECT = ff->DIPCONNECT_in_in();
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 7; j++) {
          std::string key = "dip5_dipconnect_" + std::to_string(i) + "_" + std::to_string(j);
          _DEBUG_TABLE(key, DIPCONNECT[i][j]);
        }
      }

      vector<double> Dipole_lorentz =
          D_ai_b_vector(emitter, emitted_parton, papi, x, 1, -0.5 * cA, ztb, PbKi);
      struct Coupling Cw[4] = {0, 0, 0, 0};
      if (ii > 30) {
        ii = ii - 31;
        Cw[0] = params->CHSQq[jj][ii][aa];
        Cw[1] = params->gSQSQ[ii][ii];
        Cw[2] = params->CHSQq[jj][ii][aa];
        Cw[3] = params->gSQSQ[ii][ii];
        //_DEBUG_MSG1("ii=%d", ii)
      } else {
        jj = jj - 31;
        Cw[0] = params->CHSQq[ii][jj][bb];
        Cw[1] = params->gSQSQ[jj][jj];
        Cw[2] = params->CHSQq[ii][jj][bb];
        Cw[3] = params->gSQSQ[jj][jj];
        //_DEBUG_MSG1("jj=%d", jj)
      }

      ff->SetWCoupling(Cw);
      vector<double> pol_born = ff->Mborn_pol(s1, PbKt1);

      delete ff;
      auto tmp = einsum_i_ij_j(Dipole_lorentz, DIPCONNECT, pol_born) / pow2(4.0 * (M_PI));
      // if(tmp > 0.1/ pow2(4.0 * (M_PI))) std::cout << "II: ";
      return tmp;
    } else if ((emitter == PARTON_ANTIQUARK && emitted_parton == PARTON_ANTIQUARK) &&
               (spectator == PARTON_QUARK || spectator == PARTON_ANTIQUARK)) { // D_qq_q

      double papi = 0.0;
      double s1 = 0.0;
      double ztb = 0.0;
      double PbKi = 0.0;
      double PbKt1 = 0.0;

      if (Opt_eo == 1)
        ff->Cross_pa_pb();
      // Set dipole kinematics for emitter B and spectator A for polarized Bon
      // cross section contraction
      ff->SetDipKinematicAB_pol(x, s1, papi, ztb, PbKi, PbKt1);

      vector<vector<double>> DIPCONNECT = ff->DIPCONNECT_in_in();

      vector<double> Dipole_lorentz =
          D_ai_b_vector(emitter, emitted_parton, papi, x, 1, -0.5 * cF, ztb, PbKi);
      struct Coupling Cw[4] = {0, 0, 0, 0};
      if (ii > 30) {
        ii = ii - 31;
        Cw[0] = params->CHSQq[jj][ii][Opt_eo == 1 ? bb : aa];
        Cw[1] = params->gSQSQ[ii][ii];
        Cw[2] = params->CHSQq[jj][ii][Opt_eo == 1 ? bb : aa];
        Cw[3] = params->gSQSQ[ii][ii];
        //_DEBUG_MSG1("ii=%d", ii)
      } else {
        jj = jj - 31;
        Cw[0] = params->CHSQq[ii][jj][Opt_eo == 1 ? bb : aa];
        Cw[1] = params->gSQSQ[jj][jj];
        Cw[2] = params->CHSQq[ii][jj][Opt_eo == 1 ? bb : aa];
        Cw[3] = params->gSQSQ[jj][jj];
        //_DEBUG_MSG1("jj=%d", jj)
      }

      ff->SetWCoupling(Cw);
      vector<double> pol_born = ff->Mborn_pol(s1, PbKt1);

      if (Opt_eo == 1)
        ff->Cross_pa_pb();
      delete ff;
      return einsum_i_ij_j(Dipole_lorentz, DIPCONNECT, pol_born) / pow2(4.0 * (M_PI));
    }
  case INITIAL_FINAL:
    if (emitter == PARTON_GLUON && spectator == PARTON_SQUARK &&
        emitted_parton == PARTON_ANTIQUARK) {
      if (Opt_eo == 0) {
        ff->SetDipKinematicA1(x, zi, zj);
        Dipole = D_ai_j(PARTON_GLUON, PARTON_ANTIQUARK, zi, zj, ff->pap3, ff->p1p3, x, 1,
                        0.5 * (cA - 2.0 * cF) / cF * TR);
      } else if (Opt_eo == 1) {
        ff->SetDipKinematicB1(x, zi, zj);
        Dipole = D_ai_j(PARTON_GLUON, PARTON_ANTIQUARK, zi, zj, ff->pbp3, ff->p1p3, x, 1,
                        0.5 * (cA - 2.0 * cF) / cF * TR);
        _DEBUG_TABLE("dip_papi", ff->pbp3 * x);
        _DEBUG_TABLE("xija", x);
        _DEBUG_TABLE("ztj", zj);
        _DEBUG_TABLE("col", 0.5 * (cA - 2.0 * cF));
        _DEBUG_TABLE("dip_dip", Dipole);
      }
    }
    if ((emitter == PARTON_QUARK || emitter == PARTON_ANTIQUARK) && spectator == PARTON_SQUARK &&
        emitted_parton == PARTON_GLUON) { // D_qg_sq
      ff->SetDipKinematicA1(x, zi, zj);
      Dipole = D_ai_j(PARTON_QUARK, PARTON_GLUON, zi, zj, ff->pap3, ff->p1p3, x, 1,
                      0.5 * (cA - 2.0 * cF));
    } else if ((emitter == PARTON_ANTIQUARK && emitted_parton == PARTON_ANTIQUARK &&
                spectator == PARTON_SQUARK)) { // D_gg_sq D_qq_sq

      double x = 0.0;
      double Q2 = 0.0;
      double P1Ki = 0.0;
      double ztj = 0.0;
      double Pt1Kt1 = 0.0;

      if (Opt_eo == 1)
        ff->Cross_pa_pb();

      // Set dipole kinematics for emitter B and spectator 1 for polarized Bon
      // cross section contraction
      ff->SetDipKinematicB1_pol(x, Q2, P1Ki, ztj, Pt1Kt1);

      vector<vector<double>> DIPCONNECT = ff->DIPCONNECT_in_out();

      vector<double> Dipole_lorentz =
          D_ai_j_vector(emitter, emitted_parton, P1Ki, x, 1, -0.5 * cF, ztj, ff->p1p3);
      struct Coupling Cw[4] = {0, 0, 0, 0};
      if (ii > 30) {
        ii = ii - 31;
        Cw[0] = params->CHSQq[jj][ii][Opt_eo == 1 ? bb : aa];
        Cw[1] = params->gSQSQ[ii][ii];
        Cw[2] = params->CHSQq[jj][ii][Opt_eo == 1 ? bb : aa];
        Cw[3] = params->gSQSQ[ii][ii];
        //_DEBUG_MSG1("ii=%d", ii)
      } else {
        jj = jj - 31;
        Cw[0] = params->CHSQq[ii][jj][Opt_eo == 1 ? bb : aa];
        Cw[1] = params->gSQSQ[jj][jj];
        Cw[2] = params->CHSQq[ii][jj][Opt_eo == 1 ? bb : aa];
        Cw[3] = params->gSQSQ[jj][jj];
        //_DEBUG_MSG1("jj=%d", jj)
      }

      ff->SetWCoupling(Cw);
      vector<double> pol_born = ff->Mborn_pol(Q2, Pt1Kt1);

      if (Opt_eo == 1)
        ff->Cross_pa_pb();
      delete ff;
      return einsum_i_ij_j(Dipole_lorentz, DIPCONNECT, pol_born) / pow2(4.0 * (M_PI));
    } else if ((emitter == PARTON_GLUON && spectator == PARTON_SQUARK &&
                emitted_parton == PARTON_GLUON)) { // D_gg_sq D_qq_sq

      // ff->papb = -ff->papb;

      double x = 0.0;
      double Q2 = 0.0;
      double P1Ki = 0.0;
      double ztj = 0.0;
      double Pt1Kt1 = 0.0;

      // Set dipole kinematics for emitter B and spectator 1 for polarized Bon
      // cross section contraction
      ff->SetDipKinematicB1_pol(x, Q2, P1Ki, ztj, Pt1Kt1);

      vector<vector<double>> DIPCONNECT = ff->DIPCONNECT_in_out();

      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 7; j++) {
          std::string key = "dip6_dipconnect_" + std::to_string(i) + "_" + std::to_string(j);
          _DEBUG_TABLE(key, DIPCONNECT[i][j]);
        }
      }

      vector<double> Dipole_lorentz =
          D_ai_j_vector(emitter, emitted_parton, P1Ki, x, 1, -0.5 * cA, ztj, ff->p1p3);
      struct Coupling Cw[4] = {0, 0, 0, 0};
      if (ii > 30) {
        ii = ii - 31;
        Cw[0] = params->CHSQq[jj][ii][aa];
        Cw[1] = params->gSQSQ[ii][ii];
        Cw[2] = params->CHSQq[jj][ii][aa];
        Cw[3] = params->gSQSQ[ii][ii];
        //_DEBUG_MSG1("ii=%d", ii)
      } else {
        jj = jj - 31;
        Cw[0] = params->CHSQq[ii][jj][bb];
        Cw[1] = params->gSQSQ[jj][jj];
        Cw[2] = params->CHSQq[ii][jj][bb];
        Cw[3] = params->gSQSQ[jj][jj];
        //_DEBUG_MSG1("jj=%d", jj)
      }

      ff->SetWCoupling(Cw);
      vector<double> pol_born = ff->Mborn_pol(Q2, Pt1Kt1);

      delete ff;
      auto tmp = einsum_i_ij_j(Dipole_lorentz, DIPCONNECT, pol_born) / pow2(4.0 * (M_PI));
      // if(tmp > 0.1/ pow2(4.0 * (M_PI))) std::cout << "IO: ";
      return tmp;
    }
  case FINAL_INITIAL:
    if (emitter == PARTON_SQUARK && (spectator == PARTON_QUARK || spectator == PARTON_ANTIQUARK) &&
        emitted_parton == PARTON_GLUON) {

      ff->SetDipKinematic1A(x, zi, zj, zplus, zminus, ff->m1, 0.0, ff->m1);
      Dipole = D_ij_a(PARTON_GLUON, PARTON_SQUARK, 0.0, ff->m1, ff->m1, ff->p1p3, x, 1, zj, ff->m1,
                      zi, zplus, zminus, 0.5 * (cA - 2.0 * cF));

    } else if (emitter == PARTON_SQUARK && spectator == PARTON_GLUON &&
               emitted_parton == PARTON_GLUON) {
      // this case we have a squark as particle 1, meaning that we have a quark
      // in the initial state as particle A
      ff->SetDipKinematic1B(x, zi, zj, zplus, zminus, ff->m1, 0.0, ff->m1);
      Dipole = D_ij_a(PARTON_GLUON, PARTON_SQUARK, 0.0, ff->m1, ff->m1, ff->p1p3, x, 1, zj, ff->m1,
                      zi, zplus, zminus, -0.5 * cA);
    }
  }

  ///*
  struct Coupling Cw[4] = {0, 0, 0, 0};
  if (ii > 30) {
    ii = ii - 31;
    Cw[0] = params->CHSQq[jj][ii][aa];
    Cw[1] = params->gSQSQ[ii][ii];
    Cw[2] = params->CHSQq[jj][ii][aa];
    Cw[3] = params->gSQSQ[ii][ii];
    //_DEBUG_MSG1("ii=%d", ii)
  } else {
    jj = jj - 31;
    Cw[0] = params->CHSQq[ii][jj][bb];
    Cw[1] = params->gSQSQ[jj][jj];
    Cw[2] = params->CHSQq[ii][jj][bb];
    Cw[3] = params->gSQSQ[jj][jj];
    //_DEBUG_MSG1("jj=%d", jj)
  }

  //   ff->SetPropagator(params->mSQ[ii], params->mSQ[ii], params->mSQ[ii]
  //   * 1.0e-2, params->mSQ[ii] * 1.0e-2);

  ff->SetWCoupling(Cw);
  //
  double long MQ, MX, dM2o, Dn_p1mk2_MQ, P1K1, Q2, P1K2;

  MX = ff->m2;
  MQ = ff->m1;
  Q2 = (ff->pap1 + ff->pap2) * 2;
  if (Opt_eo == 0)
    P1K1 = ff->pap1;
  else if (Opt_eo == 1)
    P1K1 = ff->pbp1;

  // MQ = self.mass_o[0]; MX = self.mass_o[1]

  dM2o = MX * MX - MQ * MQ;
  P1K2 = Q2 / 2 - P1K1;
  Dn_p1mk2_MQ = dM2o - 2 * P1K2;
  // return the kinematic part of the born-cross-section
  born =
      (1 + (dM2o - 2 * P1K2) / Q2) - 2 * (dM2o / Q2 + MQ * MQ * dM2o / (Dn_p1mk2_MQ * Dn_p1mk2_MQ) +
                                          dM2o / Dn_p1mk2_MQ - dM2o * dM2o / (Q2 * Dn_p1mk2_MQ));
  born *= real(ff->LLLL + ff->RRRR) / 3 / 4;

  if (Emitter_Spectator == INITIAL_FINAL && (emitter == PARTON_GLUON) &&
      spectator == PARTON_SQUARK && emitted_parton == PARTON_ANTIQUARK) { // D_qg_q
    if (Opt_eo == 1) {
      _DEBUG_TABLE("dip_born", born * 96 / (std::norm(params->gqq[0][0].R)));
    }
  }

  //*/
  delete ff;
  return 1. / (4 * (M_PI)) * Dipole * born;
}
double gausq_qg_gluon(double s, double mi2, double pt2, double th, double ph, Parameters *params) {
  return ((real_gluon_gaugino_squark_qg(s, mi2, pt2, th, ph, 1, params)) +
          (real_gluon_gaugino_squark_qg(s, mi2, pt2, th, ph, 0, params)));
}
double gausq_qg_gluon_dip(double s, double mi2, double pt2, double th, double ph,
                          Parameters *params) {
  return (
      // Initial state emitter (quark) and final state
      // spectator (squark)
      +Dip_SQGA(INITIAL_FINAL, PARTON_QUARK, PARTON_SQUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                params) +
      Dip_SQGA(INITIAL_FINAL, PARTON_QUARK, PARTON_SQUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

      // Final state emitter (squark) and initial state
      // spectator (quark)
      + Dip_SQGA(FINAL_INITIAL, PARTON_SQUARK, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                 params) +
      Dip_SQGA(FINAL_INITIAL, PARTON_SQUARK, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

      // Final state emitter (squark) and initial state
      // spectator (gluon)
      + Dip_SQGA(FINAL_INITIAL, PARTON_SQUARK, PARTON_GLUON, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                 params) +
      Dip_SQGA(FINAL_INITIAL, PARTON_SQUARK, PARTON_GLUON, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

      // Initial state emitter (quark) and initial state
      // spectator (gluon)
      + Dip_SQGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_GLUON, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                 params) +
      Dip_SQGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_GLUON, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

      // Initial state emitter (gluon) and final state
      // spectator (squark)
      + Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                 params) +
      Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

      // Initial state emitter (gluon) and initial state
      // spectator (quark)
      + Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 1,
                 params) +
      Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2, th, ph, 0,
               params)

  );
}

double gausq_qg_gluon_minus_dip(double s, double mi2, double pt2, double th, double ph,
                                Parameters *params) {
  return 0.
#ifdef REAL
         + gausq_qg_gluon(s, mi2, pt2, th, ph, params)
#endif
#ifdef DIPOLE
         - gausq_qg_gluon_dip(s, mi2, pt2, th, ph, params)
#endif
      ;
}

double gausq_gg_gluon(double s, double mi2, double pt2, double th, double ph, Parameters *params) {
  return ((real_gluon_gaugino_squark_gg(s, mi2, pt2, th, ph, 1, params)) +
          (real_gluon_gaugino_squark_gg(s, mi2, pt2, th, ph, 0, params)));
}
double gausq_gg_gluon_dip(double s, double mi2, double pt2, double th, double ph,
                          Parameters *params) {
  return (Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 1, params, 0) +
          Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 0, params, 0) +
          Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 1, params, 1) +
          Dip_SQGA(INITIAL_FINAL, PARTON_GLUON, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 0, params, 1) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_GLUON, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 1, params, 0) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_GLUON, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 0, params, 0) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_GLUON, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 1, params, 1) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_GLUON, PARTON_ANTIQUARK, s, mi2, pt2, th,
                   ph, 0, params, 1));
}

double gausq_gg_gluon_minus_dip(double s, double mi2, double pt2, double th, double ph,
                                Parameters *params) {
  return 0.
#ifdef REAL
         + gausq_gg_gluon(s, mi2, pt2, th, ph, params)
#endif
#ifdef DIPOLE
         - gausq_gg_gluon_dip(s, mi2, pt2, th, ph, params)
#endif
      ;
}
double gausq_qqb_antiquark(double s, double mi2, double pt2, double th, double ph,
                           Parameters *params) {
  double g3s = std::norm(params->gqq[0][0].R);
  return (real_quarkb_gaugino_squark_uub_UXub(s, mi2, pt2, th, ph, 0, params) +
          real_quarkb_gaugino_squark_uub_UXub(s, mi2, pt2, th, ph, 1, params)) *
         pow2(aS(params->murs, params->set)) /
         (pow2(4.0 * (M_PI)*aS(params->murs, params->set)) / g3s);
}
double gausq_qqb_antiquark_dip(double s, double mi2, double pt2, double th, double ph,
                               Parameters *params) {
  return (
      // correct number of dipoles even in the loop
      // since the couplings in the dipole-LO get zero
      // for wrong particle type
      Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
               ph, 0, params) +
      Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
               ph, 1, params) +
      Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
               ph, 0, params) +
      Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2, th,
               ph, 1, params)

  );
}

double gausq_qqb_antiquark_minus_dip(double s, double mi2, double pt2, double th, double ph,
                                     Parameters *params) {
  double g3s = std::norm(params->gqq[0][0].R);
  return 0.
#ifdef REAL
         + gausq_qqb_antiquark(s, mi2, pt2, th, ph, params)
#endif
#ifdef DIPOLE
         - gausq_qqb_antiquark_dip(s, mi2, pt2, th, ph, params)
#endif
      ;
}

double gausq_qq_quark(double s, double mi2, double pt2, double th, double ph, Parameters *params) {
  double g3s = std::norm(params->gqq[0][0].R);
  return +(real_quark_gaugino_squark_uu_UXu(s, mi2, pt2, th, ph, 0, params) +
           real_quark_gaugino_squark_uu_UXu(s, mi2, pt2, th, ph, 1, params)) *
         pow2(aS(params->murs, params->set)) /
         (pow2(4.0 * (M_PI)*aS(params->murs, params->set)) / g3s);
}
double gausq_qq_quark_dip(double s, double mi2, double pt2, double th, double ph,
                          Parameters *params) {
  // QUARK == ANTIQUARK (here)
  return (Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 0, params, 0) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 1, params, 0) +
          Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 0, params, 0) +
          Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 1, params, 0) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 0, params, 1) +
          Dip_SQGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 1, params, 1) +
          Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 0, params, 1) +
          Dip_SQGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_SQUARK, PARTON_ANTIQUARK, s, mi2, pt2,
                   th, ph, 1, params, 1));
}
double gausq_qq_quark_minus_dip(double s, double mi2, double pt2, double th, double ph,
                                Parameters *params) {
  return 0.

#ifdef REAL
         + gausq_qq_quark(s, mi2, pt2, th, ph, params)
#endif
#ifdef DIPOLE
         - gausq_qq_quark_dip(s, mi2, pt2, th, ph, params)
#endif
      ;
}