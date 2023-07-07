// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Born matrix elements.
//
// Notes on notation:
// - Squared Born matrix elements in D = 4 - 2 \epsilon dimensions.
// - The incoming particle (antiparticle) has momentum pb (pa) and the
//   outgoing particle (antiparticle) has momentum p2 (p1).
// TODO: Is this still true?
// - Four general electroweak coupling coefficients referring to the four
//   vertices. Distinguishes between left- and right-handed.
// - General vector coupling: \Gamma_i = \gamma_mu (L_i P_L + R_i P_R )
//   e.g. LLLL = L_1 L_2 L_3^* L_4^* (vqq,vll,vqq,vll).
// - General scalar coupling analoge.
// - MBss:
//   - M = squared Matrix element
//   - B = Born level
//   - s = channel 1 is s-channel
//   - s = channel 2 is s-channel
// - Spin and color sum and average are implemented in `hxs.cc`.

#include "kinematics.h"
#include <complex>
#include <iostream>

using namespace std;

// Process: q + \bar{q} -> f + \bar{f}.

// Squared Born matrix elements of O(\epsilon^0).
double FI::MBss() {
  return real(8.0 * ivs3v1 * ivs3v2 *
              (m1 * m2 * papb * (LLRL + LRRR + RLLL + RRLR) +
               2.0 * (pap1 * pbp2 * (LRLR + RLRL) + pap2 * pbp1 * (LLLL + RRRR))));
}

double FI::MBtt() {
  return real(4.0 * ivt1s1 * ivt1s2 * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR));
}
double FI::MBuu() {
  return real(4.0 * ivu2s1 * ivu2s2 * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR));
}
double FI::MBst() {
  return real(4.0 * ivs3v1 * ivt1s2 *
              (2.0 * pap1 * pbp2 * (LRLR + RLRL) + m1 * m2 * papb * (LLRL + RRLR)));
}
double FI::MBsu() {
  return real(-4.0 * ivs3v1 * ivu2s2 *
              (m1 * m2 * papb * (LRLR + RLRL) + 2.0 * pap2 * pbp1 * (LLRL + RRLR)));
}
double FI::MBtu() {
  return real(-2.0 * ivt1s1 * ivu2s2 *
              ((LRLR + RLRL) * m1 * m2 * papb +
               (LLLL + RRRR) * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)));
}

// Squared Born matrix elements of O(\epsilon^1).
double FI::MB1ss() {
  return real(-8.0 * ivs3v1 * ivs3v2 *
              (LLRL * m1 * m2 * papb + LRRR * m1 * m2 * papb + LLLL * p1p2 * papb +
               LRLR * p1p2 * papb + 3.0 * LLLL * pap2 * pbp1 - 3.0 * LRLR * pap2 * pbp1 -
               3.0 * LLLL * pap1 * pbp2 + 3.0 * LRLR * pap1 * pbp2 + m1 * m2 * papb * RLLL +
               p1p2 * papb * RLRL - 3.0 * pap2 * pbp1 * RLRL + 3.0 * pap1 * pbp2 * RLRL +
               m1 * m2 * papb * RRLR + p1p2 * papb * RRRR + 3.0 * pap2 * pbp1 * RRRR -
               3.0 * pap1 * pbp2 * RRRR));
}

double FI::MB1st() {
  return real(-4.0 * ivs3v1 * ivt1s2 *
              (LLRL * m1 * m2 * papb + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LRLR + RLRL) +
               m1 * m2 * papb * RRLR));
}

double FI::MB1su() {
  return real(4.0 * ivs3v1 * ivu2s2 *
              (LRLR * m1 * m2 * papb + LLRL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
               m1 * m2 * papb * RLRL + p1p2 * papb * RRLR + pap2 * pbp1 * RRLR -
               pap1 * pbp2 * RRLR));
}

// Squared Born matrix elements of O(\epsilon^2).
double FI::MB2ss() {
  return real(16.0 * ivs3v1 * ivs3v2 * (pap2 * pbp1 - pap1 * pbp2) * (LLLL - LRLR - RLRL + RRRR));
}

// Squared Born matrix element for slepton pair production.
// Process: q + \bar{q} -> sl + \bar{sl}.
// Squared Born matrix elements of O(\epsilon^0).
double FI::MBssSL() {
  return real(8.0 * ivs3v1 * ivs3v2 * (LRLR + LLLL) *
              (2.0 * pap1 * pap2 - pap1 * m2s - pap2 * m1s));
}

// Squared Born Matrix elements for gaugino-squark production  of O(\epsilon^0).
// Here I have included the average and color factor
/*
double FI::MBss_SQGA() {
    return real(ivs1s1 * ivs1s2 * papb * pbp2 *  8.0 * (LRLR + RLRL));

}

double FI::MBuu_SQGA() {
    return real(ivu2s1 * ivu2s2 * (- 2.0 * pap2 * RRLR * m2s - 2.0 * pap2 * RRLR
* m1s - 2.0 * pap2 * LLRL * m2s - 2.0 * pap2 * LLRL * m1s
                                   + 4.0 * pap2 * pap2 * RRLR + 4.0 * pap2 *
pap2 * LLRL - 4.0 * pap1 * pap2 * RRLR - 4.0 * pap1 * pap2 * LLRL
                                   + 4.0 * p1p2 * pap2 * RRLR + 4.0 * p1p2 *
pap2 * LLRL));


}

double FI::MBsu_SQGA() {
    return real(ivs1s1 * ivu2s2 * (4.0 * pap2 * pbp2 * RLRR + 4.0 * pap2 * pbp2
* LRLL - 2.0 * pap2 * pbp1 * RLRR - 2.0 * pap2 * pbp1 *
                                   LRLL + 4.0 * pap2 * pap2 * RLRR + 4.0 * pap2
* pap2 * LRLL - 2 * pap1 * pbp2 * RLRR - 2.0 * pap1 * pbp2 *
                                   LRLL - 4.0 * pap1 * pap2 * RLRR - 4.0 * pap1
* pap2 * LRLL - 2.0 * papb * RLRR * m2s - 2.0 * papb *
                                   LRLL * m2s + 2.0 * p1p2 * papb * RLRR + 2.0 *
p1p2 * papb * LRLL));




}

// O(\epsilon^1).
double FI::MB1ss_SQGA() {
    return real(- ivs1s1 * ivs1s2 * 8.0 * papb * pbp2 * (LRLR + RLRL));
}

*/
