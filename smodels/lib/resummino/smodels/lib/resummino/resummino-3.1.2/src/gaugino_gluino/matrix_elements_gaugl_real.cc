// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// This File includes all the squared matrix elements for gaugino gluino pair
// production
// with three final state particles.

// NOTE: Squared amplitude for real quark (antiquark) emission is split in
// |M_res|^2 + |M_nonres|^2 + interference
// terms for the onshell subtraction approach.

#include "kinematics.h"
#include "utils.h"
#include <complex>
#include <iostream>
using namespace std;

#define MPIsI (1.0 / pow2(M_PI))
// color factors
#define cA 3.0
#define NC 3.0
#define cF (4.0 / 3.0)

// Real Gluon Emission
// (can be simplified, but does not increase speed significantly)
// FORM output
double FI::MttG_GLGA() {

  complex<double> MttG;
  // squared matrix element proportional to cA^0
  MttG =
      (-4.0 * pow(cF, 2) * NC *
       (2.0 * ivt2s1 * ivt2s2 * ivt3 * p1p3 * pbp2 *
            (LLLL + LRLR + RLRL + RRRR) +
        ivt1s2 *
            (ivt1s1 * pap1 *
                 (2.0 * ivu3 * p2p3 +
                  ivt2s2 *
                      (m2s - p1p2 + pap2 +
                       2.0 * ivu3 * (p2p3 * (papb - pbp1 + 2.0 * pbp2) +
                                     pbp2 * (-p1p3 + pap3 -
                                             2.0 * (papb - pbp1 + pbp2))))) *
                 (LLLL + LRLR + RLRL + RRRR) +
             ivt2s1 *
                 (ivt2s2 * LRLR * m1s * pbp2 - ivt2s2 * LRLR * p1p2 * pbp2 -
                  2.0 * ivu3 * LRLR * pap1 * pbp2 +
                  ivt1s1 * ivt2s2 * LRLR * m1s * pap1 * pbp2 +
                  ivt1s1 * ivt2s2 * LRLR * m2s * pap1 * pbp2 -
                  2.0 * ivt1s1 * ivt2s2 * LRLR * p1p2 * pap1 * pbp2 -
                  2.0 * ivt1s1 * ivt2s2 * LRLR * pow(pap1, 2) * pbp2 +
                  2.0 * ivt1s1 * ivt2s2 * LRLR * pap1 * pap2 * pbp2 -
                  2.0 * ivt1s1 * ivt2s2 * LRLR * pap1 * papb * pbp2 +
                  ivt2s2 * LRLR * pbp1 * pbp2 +
                  2.0 * ivu3 * LRLR * pbp1 * pbp2 +
                  2.0 * ivt1s1 * ivt2s2 * LRLR * pap1 * pbp1 * pbp2 -
                  2.0 * ivt1s1 * ivt2s2 * LRLR * pap1 * pow(pbp2, 2) +
                  LLLL *
                      (-2.0 * ivt3 * (p1p3 - pap1) *
                           (pap2 + 2.0 * ivu3 * p2p3 * papb) -
                       2.0 * ivu3 * p2p3 * pbp1 + ivt2s2 * m1s * pbp2 -
                       p1p2 * (1.0 +
                               ivt2s2 * (pbp2 + 2.0 * ivt1s1 * pap1 * pbp2)) +
                       pbp2 * (2.0 * ivu3 * (-pap1 + pbp1) +
                               ivt2s2 *
                                   (pbp1 +
                                    ivt1s1 * pap1 * (m1s + m2s -
                                                     2.0 * (pap1 - pap2 + papb -
                                                            pbp1 + pbp2))) -
                               2.0 * ivt3 *
                                   (2.0 * ivt2s2 * pow(pap1, 2) +
                                    ivt2s2 * p1p3 * pap2 -
                                    (ivt2s2 + 2.0 * ivu3) * p1p3 * papb +
                                    pap1 * (1.0 + 4.0 * ivu3 * papb +
                                            ivt2s2 * (-2.0 * p1p3 + p2p3 -
                                                      2.0 * pap2 + 2.0 * papb -
                                                      pbp3))))) +
                  ivt2s2 * m1s * pbp2 * RLRL - ivt2s2 * p1p2 * pbp2 * RLRL -
                  2.0 * ivu3 * pap1 * pbp2 * RLRL +
                  ivt1s1 * ivt2s2 * m1s * pap1 * pbp2 * RLRL +
                  ivt1s1 * ivt2s2 * m2s * pap1 * pbp2 * RLRL -
                  2.0 * ivt1s1 * ivt2s2 * p1p2 * pap1 * pbp2 * RLRL -
                  2.0 * ivt1s1 * ivt2s2 * pow(pap1, 2) * pbp2 * RLRL +
                  2.0 * ivt1s1 * ivt2s2 * pap1 * pap2 * pbp2 * RLRL -
                  2.0 * ivt1s1 * ivt2s2 * pap1 * papb * pbp2 * RLRL +
                  ivt2s2 * pbp1 * pbp2 * RLRL +
                  2.0 * ivu3 * pbp1 * pbp2 * RLRL +
                  2.0 * ivt1s1 * ivt2s2 * pap1 * pbp1 * pbp2 * RLRL -
                  2.0 * ivt1s1 * ivt2s2 * pap1 * pow(pbp2, 2) * RLRL +
                  (ivt2s2 * m1s * pbp2 -
                   2.0 * ivu3 * (p2p3 * pbp1 + (pap1 - pbp1) * pbp2) -
                   p1p2 * (1.0 + ivt2s2 * (pbp2 + 2.0 * ivt1s1 * pap1 * pbp2)) +
                   ivt2s2 * pbp2 * (pbp1 +
                                    ivt1s1 * pap1 * (m1s + m2s -
                                                     2.0 * (pap1 - pap2 + papb -
                                                            pbp1 + pbp2)))) *
                      RRRR +
                  2.0 * ivt3 *
                      (LRLR * (-2.0 * ivt2s2 * pow(pap1, 2) * pbp2 +
                               p1p3 * (2.0 * ivu3 * papb +
                                       ivt2s2 * (-pap2 + papb)) *
                                   pbp2 +
                               pap1 * (pap2 +
                                       2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                                       2.0 * ivt2s2 * pap2 * pbp2 +
                                       pbp2 * (-1.0 +
                                               ivt2s2 * (2.0 * p1p3 - p2p3 -
                                                         2.0 * papb + pbp3)))) +
                       p1p3 * (2.0 * ivu3 * papb + ivt2s2 * (-pap2 + papb)) *
                           pbp2 * RLRL -
                       p1p3 * (pap2 + 2.0 * ivu3 * p2p3 * papb +
                               ivt2s2 * pap2 * pbp2 -
                               (ivt2s2 + 2.0 * ivu3) * papb * pbp2) *
                           RRRR -
                       2.0 * ivt2s2 * pow(pap1, 2) * pbp2 * (RLRL + RRRR) +
                       pap1 * (pap2 + 2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                               2.0 * ivt2s2 * pap2 * pbp2 +
                               pbp2 * (-1.0 +
                                       ivt2s2 * (2.0 * p1p3 - p2p3 -
                                                 2.0 * papb + pbp3))) *
                           (RLRL + RRRR))))));

  // squared matrix element proportional to cA^1
  MttG +=
      cA *
      (2.0 * cF * NC *
       (ivt2s2 *
            (2.0 * ivt2s1 * (-1.0 + 2.0 * ivt3 * pap1) * pbp2 *
                 (LLLL + LRLR + RLRL + RRRR) -
             8.0 * pow(iv13, 2) * ivt2s1 * m1s * (pap1 + pap3) * pbp2 *
                 (LLLL + LRLR + RLRL + RRRR) +
             ivt1s1 * (-(LLLL * (pap2 + 2.0 * ivu3 * p2p3 * papb)) +
                       2.0 * ivu3 * LLLL * (pap1 + papb) * pbp2 +
                       2.0 * ivu3 * (pap1 + papb) * pbp2 * (LRLR + RLRL) -
                       (pap2 + 2.0 * ivu3 * p2p3 * papb -
                        2.0 * ivu3 * (pap1 + papb) * pbp2) *
                           RRRR) +
             2.0 * iv13 *
                 (ivt2s1 *
                      (m1s + 2.0 * (pap1 - 2.0 * ivt3 * pow(pap1, 2) + pap3)) *
                      pbp2 * (LLLL + LRLR + RLRL + RRRR) +
                  ivt1s1 *
                      (LLLL * (pap1 + pap3) *
                           (p1p2 + 2.0 * ivu3 * p2p3 * pbp1) +
                       LLLL * (pap1 - 2.0 * ivu3 * (2.0 * pap1 + pap3) * pbp1) *
                           pbp2 +
                       (p1p2 * pap1 + pap1 * pbp2 +
                        2.0 * ivu3 * pbp1 *
                            (p2p3 * pap1 - (2.0 * pap1 + pap3) * pbp2)) *
                           (LRLR + RLRL) +
                       ((pap1 + pap3) * (p1p2 + 2.0 * ivu3 * p2p3 * pbp1) +
                        (pap1 - 2.0 * ivu3 * (2.0 * pap1 + pap3) * pbp1) *
                            pbp2) *
                           RRRR))) +
        ivt1s2 * ivt2s1 *
            (ivt2s2 * LRLR * m1s * pbp2 - ivt2s2 * LRLR * p1p2 * pbp2 -
             2.0 * ivu3 * LRLR * pap1 * pbp2 -
             4.0 * iv13 * ivt2s2 * LRLR * m1s * pap1 * pbp2 +
             4.0 * iv13 * ivt2s2 * LRLR * p1p2 * pap1 * pbp2 +
             2.0 * iv13 * ivt2s2 * LRLR * p2p3 * pap1 * pbp2 +
             4.0 * iv13 * ivt2s2 * LRLR * pow(pap1, 2) * pbp2 -
             ivt2s2 * LRLR * pap2 * pbp2 -
             2.0 * iv13 * ivt2s2 * LRLR * m1s * pap3 * pbp2 +
             2.0 * iv13 * ivt2s2 * LRLR * p1p2 * pap3 * pbp2 +
             4.0 * iv13 * ivt2s2 * LRLR * pap1 * pap3 * pbp2 +
             ivt2s2 * LRLR * papb * pbp2 + ivt2s2 * LRLR * pbp1 * pbp2 +
             2.0 * ivu3 * LRLR * pbp1 * pbp2 -
             4.0 * iv13 * ivt2s2 * LRLR * pap1 * pbp1 * pbp2 -
             2.0 * iv13 * ivt2s2 * LRLR * pap3 * pbp1 * pbp2 -
             2.0 * iv13 * ivt2s2 * LRLR * pap1 * pbp2 * pbp3 +
             LLLL *
                 (-2.0 * ivt3 * (p1p3 - pap1) *
                      (pap2 + 2.0 * ivu3 * p2p3 * papb) -
                  2.0 * ivu3 * p2p3 * pbp1 + ivt2s2 * m1s * pbp2 +
                  p1p2 * (-1.0 +
                          ivt2s2 * (-1.0 + 2.0 * iv13 * (2.0 * pap1 + pap3)) *
                              pbp2) +
                  pbp2 *
                      (2.0 * ivu3 * (-pap1 + pbp1) -
                       2.0 * ivt3 *
                           (2.0 * ivt2s2 * pow(pap1, 2) + ivt2s2 * p1p3 * pap2 -
                            (ivt2s2 + 2.0 * ivu3) * p1p3 * papb +
                            pap1 * (1.0 + 4.0 * ivu3 * papb +
                                    ivt2s2 * (-2.0 * p1p3 + p2p3 - 2.0 * pap2 +
                                              2.0 * papb - pbp3))) +
                       ivt2s2 * (-pap2 + papb + pbp1 -
                                 2.0 * iv13 *
                                     (2.0 * m1s * pap1 - p2p3 * pap1 -
                                      2.0 * pow(pap1, 2) + m1s * pap3 -
                                      2.0 * pap1 * pap3 + 2.0 * pap1 * pbp1 +
                                      pap3 * pbp1 + pap1 * pbp3)))) +
             ivt2s2 * m1s * pbp2 * RLRL - ivt2s2 * p1p2 * pbp2 * RLRL -
             2.0 * ivu3 * pap1 * pbp2 * RLRL -
             4.0 * iv13 * ivt2s2 * m1s * pap1 * pbp2 * RLRL +
             4.0 * iv13 * ivt2s2 * p1p2 * pap1 * pbp2 * RLRL +
             2.0 * iv13 * ivt2s2 * p2p3 * pap1 * pbp2 * RLRL +
             4.0 * iv13 * ivt2s2 * pow(pap1, 2) * pbp2 * RLRL -
             ivt2s2 * pap2 * pbp2 * RLRL -
             2.0 * iv13 * ivt2s2 * m1s * pap3 * pbp2 * RLRL +
             2.0 * iv13 * ivt2s2 * p1p2 * pap3 * pbp2 * RLRL +
             4.0 * iv13 * ivt2s2 * pap1 * pap3 * pbp2 * RLRL +
             ivt2s2 * papb * pbp2 * RLRL + ivt2s2 * pbp1 * pbp2 * RLRL +
             2.0 * ivu3 * pbp1 * pbp2 * RLRL -
             4.0 * iv13 * ivt2s2 * pap1 * pbp1 * pbp2 * RLRL -
             2.0 * iv13 * ivt2s2 * pap3 * pbp1 * pbp2 * RLRL -
             2.0 * iv13 * ivt2s2 * pap1 * pbp2 * pbp3 * RLRL - p1p2 * RRRR -
             2.0 * ivu3 * p2p3 * pbp1 * RRRR + ivt2s2 * m1s * pbp2 * RRRR -
             ivt2s2 * p1p2 * pbp2 * RRRR - 2.0 * ivu3 * pap1 * pbp2 * RRRR -
             4.0 * iv13 * ivt2s2 * m1s * pap1 * pbp2 * RRRR +
             4.0 * iv13 * ivt2s2 * p1p2 * pap1 * pbp2 * RRRR +
             2.0 * iv13 * ivt2s2 * p2p3 * pap1 * pbp2 * RRRR +
             4.0 * iv13 * ivt2s2 * pow(pap1, 2) * pbp2 * RRRR -
             ivt2s2 * pap2 * pbp2 * RRRR -
             2.0 * iv13 * ivt2s2 * m1s * pap3 * pbp2 * RRRR +
             2.0 * iv13 * ivt2s2 * p1p2 * pap3 * pbp2 * RRRR +
             4.0 * iv13 * ivt2s2 * pap1 * pap3 * pbp2 * RRRR +
             ivt2s2 * papb * pbp2 * RRRR + ivt2s2 * pbp1 * pbp2 * RRRR +
             2.0 * ivu3 * pbp1 * pbp2 * RRRR -
             4.0 * iv13 * ivt2s2 * pap1 * pbp1 * pbp2 * RRRR -
             2.0 * iv13 * ivt2s2 * pap3 * pbp1 * pbp2 * RRRR -
             2.0 * iv13 * ivt2s2 * pap1 * pbp2 * pbp3 * RRRR -
             2.0 * ivt3 *
                 (LRLR *
                      (2.0 * ivt2s2 * pow(pap1, 2) * pbp2 +
                       p1p3 * (ivt2s2 * pap2 - (ivt2s2 + 2.0 * ivu3) * papb) *
                           pbp2 +
                       pap1 * (-2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) -
                               pap2 * (1.0 + 2.0 * ivt2s2 * pbp2) +
                               pbp2 * (1.0 +
                                       ivt2s2 * (-2.0 * p1p3 + p2p3 +
                                                 2.0 * papb - pbp3)))) +
                  2.0 * ivt2s2 * pow(pap1, 2) * pbp2 * (RLRL + RRRR) -
                  pap1 * (pap2 + 2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                          2.0 * ivt2s2 * pap2 * pbp2 +
                          pbp2 * (-1.0 +
                                  ivt2s2 * (2.0 * p1p3 - p2p3 - 2.0 * papb +
                                            pbp3))) *
                      (RLRL + RRRR) +
                  p1p3 * (pap2 * RRRR +
                          ivt2s2 * (pap2 - papb) * pbp2 * (RLRL + RRRR) -
                          2.0 * ivu3 * papb *
                              (-(p2p3 * RRRR) + pbp2 * (RLRL + RRRR)))))));

  return (real(MttG));
}

double FI::MuuG_GLGA() {

  complex<double> MuuG;
  MuuG =
      (4.0 * pow(cF, 2) * NC *
       (-(ivu2s1 * ivu2s2 * pap2 *
          (2.0 * ivu3 * p1p3 +
           ivu1s2 * (m1s - p1p2 + pap1 +
                     2.0 * ivu3 *
                         (pbp1 * (-p2p3 + pap3 - 2.0 * (papb + pbp1 - pbp2)) +
                          p1p3 * (papb + 2.0 * pbp1 - pbp2)))) *
          (LLLL + LRLR + RLRL + RRRR)) +
        ivu1s1 *
            (-2.0 * ivt3 * ivu1s2 * p2p3 * pbp1 * (LLLL + LRLR + RLRL + RRRR) +
             ivu2s2 *
                 (-(ivu1s2 * LRLR * m2s * pbp1) + ivu1s2 * LRLR * p1p2 * pbp1 +
                  2.0 * ivu3 * LRLR * pap2 * pbp1 -
                  ivu1s2 * ivu2s1 * LRLR * m1s * pap2 * pbp1 -
                  ivu1s2 * ivu2s1 * LRLR * m2s * pap2 * pbp1 +
                  2.0 * ivu1s2 * ivu2s1 * LRLR * p1p2 * pap2 * pbp1 -
                  2.0 * ivu1s2 * ivu2s1 * LRLR * pap1 * pap2 * pbp1 +
                  2.0 * ivu1s2 * ivu2s1 * LRLR * pow(pap2, 2) * pbp1 +
                  2.0 * ivu1s2 * ivu2s1 * LRLR * pap2 * papb * pbp1 +
                  2.0 * ivu1s2 * ivu2s1 * LRLR * pap2 * pow(pbp1, 2) -
                  ivu1s2 * LRLR * pbp1 * pbp2 -
                  2.0 * ivu3 * LRLR * pbp1 * pbp2 -
                  2.0 * ivu1s2 * ivu2s1 * LRLR * pap2 * pbp1 * pbp2 +
                  LLLL *
                      (p1p2 +
                       2.0 * ivt3 * (p2p3 - pap2) *
                           (pap1 + 2.0 * ivu3 * p1p3 * papb) -
                       ivu1s2 * m2s * pbp1 +
                       ivu1s2 * p1p2 * (1.0 + 2.0 * ivu2s1 * pap2) * pbp1 +
                       pap2 * (2.0 * ivu3 -
                               ivu1s2 * ivu2s1 *
                                   (m1s + m2s + 2.0 * pap1 - 2.0 * pap2 -
                                    2.0 * papb - 2.0 * pbp1)) *
                           pbp1 -
                       (-2.0 * ivu3 * p1p3 +
                        (ivu1s2 + 2.0 * ivu3 + 2.0 * ivu1s2 * ivu2s1 * pap2) *
                            pbp1) *
                           pbp2 +
                       2.0 * ivt3 * pbp1 *
                           (pap2 - 2.0 * ivu3 * p2p3 * papb +
                            4.0 * ivu3 * pap2 * papb +
                            ivu1s2 * (p2p3 * (pap1 - 2.0 * pap2 - papb) +
                                      pap2 * (p1p3 - 2.0 * pap1 +
                                              2.0 * (pap2 + papb) - pbp3)))) -
                  ivu1s2 * m2s * pbp1 * RLRL + ivu1s2 * p1p2 * pbp1 * RLRL +
                  2.0 * ivu3 * pap2 * pbp1 * RLRL -
                  ivu1s2 * ivu2s1 * m1s * pap2 * pbp1 * RLRL -
                  ivu1s2 * ivu2s1 * m2s * pap2 * pbp1 * RLRL +
                  2.0 * ivu1s2 * ivu2s1 * p1p2 * pap2 * pbp1 * RLRL -
                  2.0 * ivu1s2 * ivu2s1 * pap1 * pap2 * pbp1 * RLRL +
                  2.0 * ivu1s2 * ivu2s1 * pow(pap2, 2) * pbp1 * RLRL +
                  2.0 * ivu1s2 * ivu2s1 * pap2 * papb * pbp1 * RLRL +
                  2.0 * ivu1s2 * ivu2s1 * pap2 * pow(pbp1, 2) * RLRL -
                  ivu1s2 * pbp1 * pbp2 * RLRL -
                  2.0 * ivu3 * pbp1 * pbp2 * RLRL -
                  2.0 * ivu1s2 * ivu2s1 * pap2 * pbp1 * pbp2 * RLRL +
                  (p1p2 - ivu1s2 * m2s * pbp1 +
                   ivu1s2 * p1p2 * (1.0 + 2.0 * ivu2s1 * pap2) * pbp1 +
                   ivu1s2 * pbp1 *
                       (ivu2s1 * pap2 *
                            (-m1s - m2s +
                             2.0 * (-pap1 + pap2 + papb + pbp1 - pbp2)) -
                        pbp2) +
                   2.0 * ivu3 * (pap2 * pbp1 + (p1p3 - pbp1) * pbp2)) *
                      RRRR +
                  2.0 * ivt3 *
                      (LRLR * (-2.0 * ivu3 * papb *
                                   (p1p3 * pap2 + (p2p3 - 2.0 * pap2) * pbp1) -
                               pap1 * (pap2 - ivu1s2 * p2p3 * pbp1 +
                                       2.0 * ivu1s2 * pap2 * pbp1) +
                               pbp1 * (2.0 * ivu1s2 * pow(pap2, 2) -
                                       ivu1s2 * p2p3 * papb +
                                       pap2 * (1.0 +
                                               ivu1s2 * (p1p3 - 2.0 * p2p3 +
                                                         2.0 * papb - pbp3)))) +
                       (-2.0 * ivu3 * papb *
                            (p1p3 * pap2 + (p2p3 - 2.0 * pap2) * pbp1) -
                        pap1 * (pap2 - ivu1s2 * p2p3 * pbp1 +
                                2.0 * ivu1s2 * pap2 * pbp1) +
                        pbp1 * (2.0 * ivu1s2 * pow(pap2, 2) -
                                ivu1s2 * p2p3 * papb +
                                pap2 * (1.0 +
                                        ivu1s2 * (p1p3 - 2.0 * p2p3 +
                                                  2.0 * papb - pbp3)))) *
                           RLRL +
                       (p2p3 * (pap1 + 2.0 * ivu3 * papb * (p1p3 - pbp1) +
                                ivu1s2 * pap1 * pbp1 -
                                ivu1s2 * (2.0 * pap2 + papb) * pbp1) +
                        pap2 * (-2.0 * ivu3 * papb * (p1p3 - 2.0 * pbp1) +
                                pbp1 - pap1 * (1.0 + 2.0 * ivu1s2 * pbp1) +
                                ivu1s2 * pbp1 *
                                    (p1p3 + 2.0 * (pap2 + papb) - pbp3))) *
                           RRRR)))));

  MuuG +=
      cA *
      (2.0 * cF * ivu2s2 * NC *
       (-(ivu2s1 * pap2 *
          (2.0 - 4.0 * ivu3 * pbp1 -
           ivu1s2 *
               (m1s - p1p2 + pap1 +
                papb * (1.0 + 2.0 * ivu3 * (p1p3 - 2.0 * pbp1)) +
                2.0 * ivu3 * (2.0 * p1p3 - p2p3 + pap3 - 2.0 * pbp1) * pbp1 +
                (-1.0 - 2.0 * ivu3 * p1p3 + 4.0 * ivu3 * pbp1) * pbp2) +
           8.0 * pow(iv13, 2) * m1s * (pbp1 + pbp3) +
           2.0 * iv13 * (pbp1 * (-2.0 +
                                 ivu1s2 * (-2.0 * p1p2 - p2p3 + 2.0 * pap1 +
                                           pap3 - 2.0 * pbp1) +
                                 4.0 * ivu3 * pbp1) -
                         (2.0 + ivu1s2 * (p1p2 - pap1 + 2.0 * pbp1)) * pbp3 +
                         m1s * (-1.0 + ivu1s2 * (2.0 * pbp1 + pbp3)))) *
          (LLLL + LRLR + RLRL + RRRR)) +
        ivu1s1 *
            (LLLL * (-2.0 * ivt3 * (p2p3 - pap2) *
                         (pap1 + papb + 2.0 * ivu3 * p1p3 * papb) -
                     pbp2 +
                     2.0 * (iv13 * (2.0 * ivt3 * pap1 * (p2p3 - 2.0 * pap2) +
                                    pap2) *
                                pbp1 -
                            ivu3 * (pap2 - 2.0 * ivt3 * p2p3 * papb +
                                    4.0 * ivt3 * pap2 * papb) *
                                pbp1 +
                            ivu3 * (-p1p3 + pbp1) * pbp2 +
                            2.0 * iv13 * ivt3 * pap1 * (p2p3 - pap2) * pbp3) +
                     p1p2 * (-1.0 + 2.0 * iv13 * (pbp1 + pbp3))) +
             2.0 * pbp1 * (iv13 * (p1p2 + pap2) + ivu3 * (-pap2 + pbp2)) *
                 (LRLR + RLRL) +
             2.0 * ivt3 *
                 (papb * (pap2 + 2.0 * ivu3 * p1p3 * pap2 +
                          2.0 * ivu3 * p2p3 * pbp1 - 4.0 * ivu3 * pap2 * pbp1) +
                  pap1 *
                      (pap2 + 2.0 * iv13 * p2p3 * pbp1 -
                       4.0 * iv13 * pap2 * pbp1 - 2.0 * iv13 * pap2 * pbp3)) *
                 (LRLR + RLRL) +
             (2.0 * (iv13 - ivu3) * pap2 * pbp1 +
              (-1.0 - 2.0 * ivu3 * p1p3 + 2.0 * ivu3 * pbp1) * pbp2 +
              p1p2 * (-1.0 + 2.0 * iv13 * (pbp1 + pbp3))) *
                 RRRR +
             2.0 * ivt3 *
                 (pap2 * (pap1 + papb + 2.0 * ivu3 * p1p3 * papb -
                          4.0 * iv13 * pap1 * pbp1 - 4.0 * ivu3 * papb * pbp1 -
                          2.0 * iv13 * pap1 * pbp3) +
                  p2p3 *
                      (papb * (-1.0 - 2.0 * ivu3 * p1p3 + 2.0 * ivu3 * pbp1) +
                       pap1 * (-1.0 + 2.0 * iv13 * (pbp1 + pbp3)))) *
                 RRRR)));

  return (real(MuuG));
}

double FI::MtuG_GLGA() {
  complex<double> MtuG;

  // -(MtuAA + MtuBB + MtuAB + MtuBA)
  MtuG =
      (-4.0 * cF * NC *
       (-(cA *
          (ivt2s1 * ivu2s2 *
               (ivu3 * (LRRL * m1 * m2 * papb +
                        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + p1p3 * pbp2 -
                                pbp1 * pbp2) +
                        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                        pap2 * pbp1 * RRRR + p1p3 * pbp2 * RRRR -
                        pbp1 * pbp2 * RRRR) +
                ivt3 *
                    (LLLL * (p1p3 * (pap2 + 2.0 * ivu3 * p2p3 * papb) -
                             papb * (p1p2 + 2.0 * ivu3 * p1p2 * papb +
                                     2.0 * ivu3 * (p2p3 - pap2) * pbp1) +
                             pap1 * (-pap2 + pbp2 +
                                     2.0 * ivu3 * papb * (-p2p3 + pbp2))) +
                     m1 * m2 * papb * (1.0 + 2.0 * ivu3 * papb) *
                         (LRRL + RLLR) +
                     (p1p3 - pap1) * pap2 * RRRR +
                     (-(p1p2 * papb * (1.0 + 2.0 * ivu3 * papb)) + pap1 * pbp2 +
                      2.0 * ivu3 * papb *
                          (p1p3 * p2p3 + pap2 * pbp1 - p2p3 * (pap1 + pbp1) +
                           pap1 * pbp2)) *
                         RRRR)) +
           ivt1s1 * ivu1s2 *
               (ivu3 * (LRRL * m1 * m2 * papb +
                        LLLL * (-(p1p2 * papb) + p2p3 * pbp1 + pap1 * pbp2 -
                                pbp1 * pbp2) +
                        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                        p2p3 * pbp1 * RRRR + pap1 * pbp2 * RRRR -
                        pbp1 * pbp2 * RRRR) +
                ivt3 *
                    (LLLL * (-(pap1 * pap2) - p1p2 * papb -
                             2.0 * ivu3 * p1p3 * pap2 * papb -
                             2.0 * ivu3 * p1p2 * pow2(papb) +
                             p2p3 * (pap1 + 2.0 * ivu3 * p1p3 * papb) +
                             pap2 * pbp1 + 2.0 * ivu3 * pap2 * papb * pbp1 +
                             2.0 * ivu3 * (-p1p3 + pap1) * papb * pbp2) +
                     m1 * m2 * papb * (1.0 + 2.0 * ivu3 * papb) *
                         (LRRL + RLLR) +
                     pap1 * (p2p3 - pap2) * RRRR +
                     (-(p1p2 * papb * (1.0 + 2.0 * ivu3 * papb)) + pap2 * pbp1 +
                      2.0 * ivu3 * papb *
                          (pap2 * pbp1 + p1p3 * (p2p3 - pap2 - pbp2) +
                           pap1 * pbp2)) *
                         RRRR)))) +
        2.0 * cF *
            (ivt1s1 *
                 (ivt3 * ivu1s2 *
                      (LLLL * (-(pap1 * pap2) - p1p2 * papb -
                               2.0 * ivu3 * p1p3 * pap2 * papb -
                               2.0 * ivu3 * p1p2 * pow2(papb) +
                               p2p3 * (pap1 + 2.0 * ivu3 * p1p3 * papb) +
                               pap2 * pbp1 + 2.0 * ivu3 * pap2 * papb * pbp1 +
                               2.0 * ivu3 * (-p1p3 + pap1) * papb * pbp2) +
                       m1 * m2 * papb * (1.0 + 2.0 * ivu3 * papb) *
                           (LRRL + RLLR) +
                       pap1 * (p2p3 - pap2) * RRRR +
                       (-(p1p2 * papb * (1.0 + 2.0 * ivu3 * papb)) +
                        pap2 * pbp1 +
                        2.0 * ivu3 * papb *
                            (pap2 * pbp1 + p1p3 * (p2p3 - pap2 - pbp2) +
                             pap1 * pbp2)) *
                           RRRR) -
                  ivu3 *
                      (ivu2s2 *
                           (LLLL * (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) +
                            m1 * m2 * pap3 * (LRRL + RLLR) +
                            (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) * RRRR) -
                       ivu1s2 * (LRRL * m1 * m2 * papb +
                                 LLLL * (-(p1p2 * papb) + p2p3 * pbp1 +
                                         pap1 * pbp2 - pbp1 * pbp2) +
                                 m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                                 p2p3 * pbp1 * RRRR + pap1 * pbp2 * RRRR -
                                 pbp1 * pbp2 * RRRR))) +
             ivt2s1 *
                 (ivu2s2 * ivu3 * (LRRL * m1 * m2 * papb +
                                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 +
                                           p1p3 * pbp2 - pbp1 * pbp2) +
                                   m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                                   pap2 * pbp1 * RRRR + p1p3 * pbp2 * RRRR -
                                   pbp1 * pbp2 * RRRR) +
                  ivt3 *
                      (ivu2s2 *
                           (LLLL *
                                (p1p3 * (pap2 + 2.0 * ivu3 * p2p3 * papb) -
                                 papb * (p1p2 + 2.0 * ivu3 * p1p2 * papb +
                                         2.0 * ivu3 * (p2p3 - pap2) * pbp1) +
                                 pap1 * (-pap2 + pbp2 +
                                         2.0 * ivu3 * papb * (-p2p3 + pbp2))) +
                            m1 * m2 * papb * (1.0 + 2.0 * ivu3 * papb) *
                                (LRRL + RLLR) +
                            (p1p3 - pap1) * pap2 * RRRR +
                            (-(p1p2 * papb * (1.0 + 2.0 * ivu3 * papb)) +
                             pap1 * pbp2 +
                             2.0 * ivu3 * papb *
                                 (p1p3 * p2p3 + pap2 * pbp1 -
                                  p2p3 * (pap1 + pbp1) + pap1 * pbp2)) *
                                RRRR) -
                       ivu1s2 *
                           (LLLL * (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) +
                            m1 * m2 * pbp3 * (LRRL + RLLR) +
                            (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) *
                                RRRR))))));

  // -(MtuA1 + MtuB1 + Mtu1A + Mtu1B)
  MtuG +=
      (-2.0 * cA * cF * NC *
       (ivt1s1 * ivu2s2 * (LLLL * pap2 * (-1.0 + 2.0 * ivu3 * pbp1) +
                           ivu3 * m1 * m2 * papb * (LRRL + RLLR) +
                           pap2 * (-1.0 + 2.0 * ivu3 * pbp1) * RRRR) +
        iv13 *
            (ivt1s1 * ivu2s2 *
                 (LRRL * m1 * m2 * (pap1 +
                                    2.0 * (-pap3 + papb +
                                           ivu3 * (pap3 - 2.0 * papb) * pbp1)) +
                  LLLL * (-2.0 * p2p3 * (pap1 + ivu3 * m1s * papb) +
                          4.0 * ivu3 * p2p3 * pap1 * pbp1 -
                          2.0 * (-1.0 + 2.0 * ivu3 * pbp1) *
                              (p1p2 * pap3 - p1p2 * papb + pap2 * pbp1 +
                               pap1 * pbp2) +
                          m1s * (pap2 + 2.0 * ivu3 * pap3 * pbp2)) +
                  m1 * m2 * (pap1 +
                             2.0 * (-pap3 + papb +
                                    ivu3 * (pap3 - 2.0 * papb) * pbp1)) *
                      RLLR +
                  (-2.0 * p2p3 * (pap1 + ivu3 * m1s * papb) +
                   4.0 * ivu3 * p2p3 * pap1 * pbp1 -
                   2.0 * (-1.0 + 2.0 * ivu3 * pbp1) *
                       (p1p2 * pap3 - p1p2 * papb + pap2 * pbp1 + pap1 * pbp2) +
                   m1s * (pap2 + 2.0 * ivu3 * pap3 * pbp2)) *
                      RRRR) +
             ivt2s1 *
                 (ivu1s2 *
                      (LRRL * m1 * m2 *
                           ((2.0 - 4.0 * ivt3 * pap1) * papb + pbp1 +
                            2.0 * (-1.0 + ivt3 * pap1) * pbp3) +
                       LLLL * (-2.0 * p2p3 * pbp1 + 2.0 * pap2 * pbp1 +
                               m1s * pbp2 + 2.0 * pap1 * pbp2 +
                               2.0 * p1p2 * (-1.0 + 2.0 * ivt3 * pap1) *
                                   (papb - pbp3) -
                               2.0 * ivt3 *
                                   (2.0 * pap1 * (-(p2p3 * pbp1) + pap2 * pbp1 +
                                                  pap1 * pbp2) +
                                    m1s * (p2p3 * papb - pap2 * pbp3))) +
                       m1 * m2 * ((2.0 - 4.0 * ivt3 * pap1) * papb + pbp1 +
                                  2.0 * (-1.0 + ivt3 * pap1) * pbp3) *
                           RLLR +
                       (-2.0 * p2p3 * pbp1 + 2.0 * pap2 * pbp1 + m1s * pbp2 +
                        2.0 * pap1 * pbp2 +
                        2.0 * p1p2 * (-1.0 + 2.0 * ivt3 * pap1) *
                            (papb - pbp3) -
                        2.0 * ivt3 * (2.0 * pap1 * (-(p2p3 * pbp1) +
                                                    pap2 * pbp1 + pap1 * pbp2) +
                                      m1s * (p2p3 * papb - pap2 * pbp3))) *
                           RRRR) +
                  ivu2s2 *
                      (LRRL * m1 * m2 * (pap1 - 4.0 * ivt3 * pap1 * papb +
                                         pbp1 + 2.0 * ivu3 * pap3 * pbp1 -
                                         4.0 * ivu3 * papb * pbp1 +
                                         2.0 * ivt3 * pap1 * pbp3) +
                       LLLL *
                           (m1s * (pap2 + 2.0 * (ivt3 + ivu3) * p2p3 * papb +
                                   pbp2 - 2.0 * ivu3 * pap3 * pbp2 -
                                   2.0 * ivt3 * pap2 * pbp3) +
                            2.0 * (pap2 * pbp1 + pap1 * pbp2 -
                                   2.0 * ivu3 * pbp1 *
                                       (-((p1p2 + p2p3) * papb) + pap2 * pbp1 +
                                        (pap1 + pap3) * pbp2) -
                                   2.0 * ivt3 * pap1 *
                                       (-((p1p2 + p2p3) * papb) + pap1 * pbp2 +
                                        pap2 * (pbp1 + pbp3)))) +
                       m1 * m2 * (pap1 - 4.0 * ivt3 * pap1 * papb + pbp1 +
                                  2.0 * ivu3 * pap3 * pbp1 -
                                  4.0 * ivu3 * papb * pbp1 +
                                  2.0 * ivt3 * pap1 * pbp3) *
                           RLLR +
                       (m1s * (pap2 + 2.0 * (ivt3 + ivu3) * p2p3 * papb + pbp2 -
                               2.0 * ivu3 * pap3 * pbp2 -
                               2.0 * ivt3 * pap2 * pbp3) +
                        2.0 * (pap2 * pbp1 + pap1 * pbp2 -
                               2.0 * ivu3 * pbp1 *
                                   (-((p1p2 + p2p3) * papb) + pap2 * pbp1 +
                                    (pap1 + pap3) * pbp2) -
                               2.0 * ivt3 * pap1 *
                                   (-((p1p2 + p2p3) * papb) + pap1 * pbp2 +
                                    pap2 * (pbp1 + pbp3)))) *
                           RRRR))) +
        ivt2s1 *
            (ivu2s2 *
                 (LLLL *
                      (-2.0 * (ivt3 + ivu3) * (p1p2 + p2p3) * papb +
                       (-1.0 + 2.0 * (ivt3 + ivu3) * pap1 + 2.0 * ivu3 * pap3) *
                           pbp2 +
                       pap2 * (-1.0 + 2.0 * (ivt3 + ivu3) * pbp1 +
                               2.0 * ivt3 * pbp3)) +
                  (ivt3 + ivu3) * m1 * m2 * papb * (LRRL + RLLR) +
                  (-2.0 * (ivt3 + ivu3) * (p1p2 + p2p3) * papb +
                   (-1.0 + 2.0 * (ivt3 + ivu3) * pap1 + 2.0 * ivu3 * pap3) *
                       pbp2 +
                   pap2 * (-1.0 + 2.0 * (ivt3 + ivu3) * pbp1 +
                           2.0 * ivt3 * pbp3)) *
                      RRRR) +
             ivu1s2 * (-(pbp2 * (LLLL + RRRR)) +
                       ivt3 * (m1 * m2 * papb * (LRRL + RLLR) +
                               2.0 * pap1 * pbp2 * (LLLL + RRRR))))));

  // -(Mtu11 + Mtu1P + MtuP1);
  MtuG +=
      (cA * cF * ivt2s1 * ivu2s2 * NC *
       ((ivt1s1 + ivu1s2) * (LRRL * m1 * m2 * papb - LLLL * m2s * papb +
                             2.0 * LLLL * pap2 * pbp2 + m1 * m2 * papb * RLLR -
                             m2s * papb * RRRR + 2.0 * pap2 * pbp2 * RRRR) +
        16.0 * pow2(iv13) * m1s *
            (LRRL * m1 * m2 * papb +
             LLLL * (-((p1p2 + p2p3) * papb) + (pap1 + pap3) * pbp2 +
                     pap2 * (pbp1 + pbp3)) +
             m1 * m2 * papb * RLLR - p1p2 * papb * RRRR - p2p3 * papb * RRRR +
             pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR + pap3 * pbp2 * RRRR +
             pap2 * pbp3 * RRRR) +
        2.0 * iv13 *
            (-2.0 * ivu1s2 * LLLL * m1s * p1p2 * papb +
             2.0 * ivu1s2 * LLLL * pow2(p1p2) * papb +
             4.0 * LLLL * p2p3 * papb - ivu1s2 * LLLL * m1s * p2p3 * papb +
             2.0 * ivu1s2 * LLLL * p1p2 * p2p3 * papb -
             2.0 * ivu1s2 * LLLL * p1p2 * pap1 * papb -
             2.0 * ivu1s2 * LLLL * p1p2 * pap3 * papb +
             2.0 * ivu1s2 * LLLL * m1s * pap2 * pbp1 -
             2.0 * ivu1s2 * LLLL * p1p2 * pap2 * pbp1 +
             2.0 * ivu1s2 * LLLL * pap1 * pap2 * pbp1 -
             ivu1s2 * LLLL * m2s * pap3 * pbp1 +
             2.0 * ivu1s2 * LLLL * pap2 * pap3 * pbp1 +
             2.0 * ivu1s2 * LLLL * p1p2 * papb * pbp1 +
             2.0 * ivu1s2 * LLLL * p2p3 * papb * pbp1 -
             2.0 * ivu1s2 * LLLL * pap2 * pow2(pbp1) +
             2.0 * ivu1s2 * LLLL * m1s * pap1 * pbp2 -
             2.0 * ivu1s2 * LLLL * p1p2 * pap1 * pbp2 -
             2.0 * ivu1s2 * LLLL * p2p3 * pap1 * pbp2 +
             2.0 * ivu1s2 * LLLL * pow2(pap1) * pbp2 -
             4.0 * LLLL * pap3 * pbp2 + ivu1s2 * LLLL * m1s * pap3 * pbp2 +
             2.0 * ivu1s2 * LLLL * pap1 * pap3 * pbp2 -
             2.0 * ivu1s2 * LLLL * pap1 * pbp1 * pbp2 -
             2.0 * ivu1s2 * LLLL * pap3 * pbp1 * pbp2 +
             ivu1s2 * LLLL * m2s * pap1 * pbp3 - 4.0 * LLLL * pap2 * pbp3 +
             ivu1s2 * LLLL * m1s * pap2 * pbp3 -
             2.0 * ivu1s2 * LLLL * p1p2 * pap2 * pbp3 -
             2.0 * ivu1s2 * LLLL * pap2 * pbp1 * pbp3 +
             LRRL * m1 * m2 *
                 (-((ivt1s1 - ivu1s2) *
                    (pap3 * (pbp1 - pbp2) + (-pap1 + pap2) * pbp3)) +
                  papb * (4.0 +
                          ivu1s2 * (2.0 * m1s - 2.0 * p1p2 - p2p3 +
                                    2.0 * (pap1 + pap3 - pbp1)) +
                          ivt1s1 * (2.0 * m1s - 2.0 * p1p2 - p2p3 +
                                    2.0 * (-pap1 + pbp1 + pbp3)))) +
             4.0 * m1 * m2 * papb * RLLR +
             2.0 * ivu1s2 * m1 * m1s * m2 * papb * RLLR -
             2.0 * ivu1s2 * m1 * m2 * p1p2 * papb * RLLR -
             ivu1s2 * m1 * m2 * p2p3 * papb * RLLR +
             2.0 * ivu1s2 * m1 * m2 * pap1 * papb * RLLR +
             2.0 * ivu1s2 * m1 * m2 * pap3 * papb * RLLR +
             ivu1s2 * m1 * m2 * pap3 * pbp1 * RLLR -
             2.0 * ivu1s2 * m1 * m2 * papb * pbp1 * RLLR -
             ivu1s2 * m1 * m2 * pap3 * pbp2 * RLLR -
             ivu1s2 * m1 * m2 * pap1 * pbp3 * RLLR +
             ivu1s2 * m1 * m2 * pap2 * pbp3 * RLLR +
             (4.0 * p2p3 * papb - 4.0 * pap3 * pbp2 - 4.0 * pap2 * pbp3 +
              ivu1s2 * (2.0 * pow2(p1p2) * papb + 2.0 * pap1 * pap2 * pbp1 -
                        m2s * pap3 * pbp1 + 2.0 * pap2 * pap3 * pbp1 +
                        2.0 * p2p3 * papb * pbp1 - 2.0 * pap2 * pow2(pbp1) -
                        2.0 * p2p3 * pap1 * pbp2 + 2.0 * pow2(pap1) * pbp2 +
                        2.0 * pap1 * pap3 * pbp2 - 2.0 * pap1 * pbp1 * pbp2 -
                        2.0 * pap3 * pbp1 * pbp2 + m2s * pap1 * pbp3 -
                        2.0 * pap2 * pbp1 * pbp3 +
                        2.0 * p1p2 * (papb * (p2p3 - pap1 - pap3 + pbp1) -
                                      pap1 * pbp2 - pap2 * (pbp1 + pbp3)) +
                        m1s * (-((2.0 * p1p2 + p2p3) * papb) +
                               (2.0 * pap1 + pap3) * pbp2 +
                               pap2 * (2.0 * pbp1 + pbp3)))) *
                 RRRR +
             ivt1s1 *
                 (LLLL * (2.0 * pow2(p1p2) * papb + 2.0 * p2p3 * pap1 * papb -
                          2.0 * p2p3 * pap2 * pbp1 - 2.0 * pap1 * pap2 * pbp1 +
                          m2s * pap3 * pbp1 + 2.0 * pap2 * pow2(pbp1) -
                          2.0 * pow2(pap1) * pbp2 - 2.0 * pap1 * pap3 * pbp2 +
                          2.0 * pap1 * pbp1 * pbp2 - m2s * pap1 * pbp3 -
                          2.0 * pap1 * pap2 * pbp3 + 2.0 * pap2 * pbp1 * pbp3 +
                          2.0 * pap1 * pbp2 * pbp3 +
                          2.0 * p1p2 * (p2p3 * papb - (pap2 + papb) * pbp1 +
                                        pap1 * (papb - pbp2) - pap3 * pbp2 -
                                        papb * pbp3) +
                          m1s * (-((2.0 * p1p2 + p2p3) * papb) +
                                 (2.0 * pap1 + pap3) * pbp2 +
                                 pap2 * (2.0 * pbp1 + pbp3))) +
                  m1 * m2 *
                      (2.0 * m1s * papb -
                       papb * (2.0 * p1p2 + p2p3 + 2.0 * pap1 - 2.0 * pbp1) +
                       pap3 * (-pbp1 + pbp2) +
                       (pap1 - pap2 + 2.0 * papb) * pbp3) *
                      RLLR +
                  (2.0 * pow2(p1p2) * papb + 2.0 * p2p3 * pap1 * papb -
                   2.0 * p2p3 * pap2 * pbp1 - 2.0 * pap1 * pap2 * pbp1 +
                   m2s * pap3 * pbp1 + 2.0 * pap2 * pow2(pbp1) -
                   2.0 * pow2(pap1) * pbp2 - 2.0 * pap1 * pap3 * pbp2 +
                   2.0 * pap1 * pbp1 * pbp2 - m2s * pap1 * pbp3 -
                   2.0 * pap1 * pap2 * pbp3 + 2.0 * pap2 * pbp1 * pbp3 +
                   2.0 * pap1 * pbp2 * pbp3 +
                   2.0 * p1p2 *
                       (p2p3 * papb - (pap2 + papb) * pbp1 +
                        pap1 * (papb - pbp2) - pap3 * pbp2 - papb * pbp3) +
                   m1s * (-((2.0 * p1p2 + p2p3) * papb) +
                          (2.0 * pap1 + pap3) * pbp2 +
                          pap2 * (2.0 * pbp1 + pbp3))) *
                      RRRR))));

  // MtuPP
  MtuG -= -2.0 * (cA - 2.0 * cF) * cF * ivt1s1 * ivt2s1 * ivu1s2 * ivu2s2 * NC *
          (m1s + m2s - 2.0 * p1p2 + 2.0 * papb) *
          (LRRL * m1 * m2 * papb +
           LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
           m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
           pap1 * pbp2 * RRRR);

  // MtuAP
  MtuG -= 2.0 * pow2(cF) * ivt2s1 * ivu1s2 * ivu2s2 * NC *
          (-(LLLL * m2s * pbp1) + LRRL * m1 * m2 * (pbp1 - pbp2) +
           LLLL * m1s * pbp2 - 2.0 * LLLL * pbp1 * pbp2 +
           m1 * m2 * pbp1 * RLLR - m1 * m2 * pbp2 * RLLR - m2s * pbp1 * RRRR +
           m1s * pbp2 * RRRR - 2.0 * pbp1 * pbp2 * RRRR +
           2.0 * ivt3 *
               (LRRL * m1 * m2 *
                    (papb * (p1p3 - p2p3 + 2.0 * (-pap1 + pap2 + papb)) +
                     (pap1 - pap2 - 2.0 * papb) * pbp3) +
                LLLL * (m1s * p2p3 * papb -
                        2.0 * p1p2 * papb * (p1p3 - pap1 + pap2 + papb) +
                        2.0 * p1p3 * pap2 * pbp1 - 2.0 * p2p3 * pap2 * pbp1 -
                        2.0 * pap1 * pap2 * pbp1 + 2.0 * pow2(pap2) * pbp1 -
                        2.0 * p2p3 * papb * pbp1 + 2.0 * pap2 * papb * pbp1 +
                        2.0 * p1p3 * pap1 * pbp2 - 2.0 * pow2(pap1) * pbp2 -
                        2.0 * p1p3 * pap2 * pbp2 + 2.0 * pap1 * pap2 * pbp2 +
                        2.0 * pap1 * papb * pbp2 - m1s * pap2 * pbp3 +
                        2.0 * p1p2 * (pap2 + papb) * pbp3 -
                        2.0 * pap1 * pbp2 * pbp3 +
                        m2s * (p1p3 * papb - pap1 * pbp3)) +
                m1 * m2 * (papb * (p1p3 - p2p3 + 2.0 * (-pap1 + pap2 + papb)) +
                           (pap1 - pap2 - 2.0 * papb) * pbp3) *
                    RLLR +
                (m1s * p2p3 * papb -
                 2.0 * p1p2 * papb * (p1p3 - pap1 + pap2 + papb) +
                 2.0 * p1p3 * pap2 * pbp1 - 2.0 * p2p3 * pap2 * pbp1 -
                 2.0 * pap1 * pap2 * pbp1 + 2.0 * pow2(pap2) * pbp1 -
                 2.0 * p2p3 * papb * pbp1 + 2.0 * pap2 * papb * pbp1 +
                 2.0 * p1p3 * pap1 * pbp2 - 2.0 * pow2(pap1) * pbp2 -
                 2.0 * p1p3 * pap2 * pbp2 + 2.0 * pap1 * pap2 * pbp2 +
                 2.0 * pap1 * papb * pbp2 - m1s * pap2 * pbp3 +
                 2.0 * p1p2 * (pap2 + papb) * pbp3 - 2.0 * pap1 * pbp2 * pbp3 +
                 m2s * (p1p3 * papb - pap1 * pbp3)) *
                    RRRR));

  // MtuPA
  MtuG -=
      (cA - 2.0 * cF) * cF * ivt1s1 * ivt2s1 * ivu1s2 * NC *
      (-(LLLL * m2s * pbp1) + LRRL * m1 * m2 * (pbp1 - pbp2) +
       LLLL * m1s * pbp2 + 2.0 * LLLL * pbp1 * pbp2 + m1 * m2 * pbp1 * RLLR -
       m1 * m2 * pbp2 * RLLR - m2s * pbp1 * RRRR + m1s * pbp2 * RRRR +
       2.0 * pbp1 * pbp2 * RRRR +
       2.0 * ivt3 *
           (LRRL * m1 * m2 *
                (p1p3 * papb - papb * (p2p3 + 2.0 * (pap1 - pap2 + papb)) +
                 (pap1 - pap2 + 2.0 * papb) * pbp3) +
            LLLL * (m2s * (-(p1p3 * papb) + pap1 * pbp3) +
                    m1s * (-(p2p3 * papb) + pap2 * pbp3) +
                    2.0 * (p1p2 * papb * (p2p3 + pap1 - pap2 + papb) +
                           p2p3 * pap1 * pbp1 - p2p3 * pap2 * pbp1 -
                           pap1 * pap2 * pbp1 + pow2(pap2) * pbp1 -
                           pap2 * papb * pbp1 + p1p3 * pap1 * pbp2 -
                           p2p3 * pap1 * pbp2 - pow2(pap1) * pbp2 +
                           pap1 * pap2 * pbp2 + p1p3 * papb * pbp2 -
                           pap1 * papb * pbp2 - p1p2 * (pap1 + papb) * pbp3 +
                           pap2 * pbp1 * pbp3)) +
            m1 * m2 *
                (p1p3 * papb - papb * (p2p3 + 2.0 * (pap1 - pap2 + papb)) +
                 (pap1 - pap2 + 2.0 * papb) * pbp3) *
                RLLR +
            (m2s * (-(p1p3 * papb) + pap1 * pbp3) +
             m1s * (-(p2p3 * papb) + pap2 * pbp3) +
             2.0 *
                 (p1p2 * papb * (p2p3 + pap1 - pap2 + papb) +
                  p2p3 * pap1 * pbp1 - p2p3 * pap2 * pbp1 - pap1 * pap2 * pbp1 +
                  pow2(pap2) * pbp1 - pap2 * papb * pbp1 + p1p3 * pap1 * pbp2 -
                  p2p3 * pap1 * pbp2 - pow2(pap1) * pbp2 + pap1 * pap2 * pbp2 +
                  p1p3 * papb * pbp2 - pap1 * papb * pbp2 -
                  p1p2 * (pap1 + papb) * pbp3 + pap2 * pbp1 * pbp3)) *
                RRRR));

  // MtuBP
  MtuG -=
      (cA - 2.0 * cF) * cF * ivt1s1 * ivu1s2 * ivu2s2 * NC *
      (LRRL * m1 * m2 *
           (pap1 - pap2 +
            2.0 * ivu3 * (p1p3 * papb - p2p3 * papb + 2.0 * pap3 * papb -
                          2.0 * pow2(papb) + pap3 * pbp1 - 2.0 * papb * pbp1 -
                          pap3 * pbp2 + 2.0 * papb * pbp2)) +
       LLLL *
           (m1s * pap2 -
            m2s * (pap1 + 2.0 * ivu3 * p1p3 * papb - 2.0 * ivu3 * pap3 * pbp1) +
            2.0 * (pap1 * (pap2 +
                           2.0 * ivu3 * (p2p3 * pbp1 -
                                         (p2p3 - pap3 + papb + pbp1) * pbp2 +
                                         pow2(pbp2))) +
                   ivu3 * (m1s * (-(p2p3 * papb) + pap3 * pbp2) +
                           2.0 * (pap2 * (p1p3 * (papb + pbp1) -
                                          pbp1 * (p2p3 + papb + pbp1 - pbp2)) +
                                  p1p2 * (p2p3 * papb +
                                          (-pap3 + papb) * (papb + pbp1) -
                                          papb * pbp2))))) +
       m1 * m2 *
           (pap1 - pap2 +
            2.0 * ivu3 * (p1p3 * papb - p2p3 * papb + 2.0 * pap3 * papb -
                          2.0 * pow2(papb) + pap3 * pbp1 - 2.0 * papb * pbp1 -
                          pap3 * pbp2 + 2.0 * papb * pbp2)) *
           RLLR +
       (m1s * pap2 -
        m2s * (pap1 + 2.0 * ivu3 * p1p3 * papb - 2.0 * ivu3 * pap3 * pbp1) +
        2.0 * (pap1 * (pap2 +
                       2.0 * ivu3 *
                           (p2p3 * pbp1 - (p2p3 - pap3 + papb + pbp1) * pbp2 +
                            pow2(pbp2))) +
               ivu3 * (m1s * (-(p2p3 * papb) + pap3 * pbp2) +
                       2.0 * (pap2 * (p1p3 * (papb + pbp1) -
                                      pbp1 * (p2p3 + papb + pbp1 - pbp2)) +
                              p1p2 * (p2p3 * papb +
                                      (-pap3 + papb) * (papb + pbp1) -
                                      papb * pbp2))))) *
           RRRR);

  // MtuPB
  MtuG -=
      2.0 * pow2(cF) * ivt1s1 * ivt2s1 * ivu2s2 * NC *
      (LRRL * m1 * m2 *
           (pap1 - pap2 +
            2.0 * ivu3 * (p1p3 * papb - p2p3 * papb - 2.0 * pap3 * papb +
                          2.0 * pow2(papb) + pap3 * pbp1 - 2.0 * papb * pbp1 -
                          pap3 * pbp2 + 2.0 * papb * pbp2)) +
       LLLL *
           (m1s * pap2 -
            m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb + 2.0 * ivu3 * pap3 * pbp1) +
            2.0 *
                (ivu3 * (-2.0 * p1p2 * p1p3 + m1s * p2p3) * papb +
                 2.0 * ivu3 *
                     (p1p3 * pap2 * pbp1 +
                      (pap3 - papb + pbp1) * (p1p2 * papb - pap2 * pbp1)) -
                 ivu3 * (2.0 * p1p3 * pap2 + m1s * pap3 - 2.0 * p1p2 * pap3 +
                         2.0 * p1p2 * papb - 2.0 * pap2 * pbp1) *
                     pbp2 +
                 pap1 * (-pap2 +
                         2.0 * ivu3 * (-(p2p3 * (papb + pbp2)) +
                                       pbp2 * (p1p3 + papb - pbp1 + pbp2))))) +
       m1 * m2 *
           (pap1 - pap2 +
            2.0 * ivu3 * (p1p3 * papb - p2p3 * papb - 2.0 * pap3 * papb +
                          2.0 * pow2(papb) + pap3 * pbp1 - 2.0 * papb * pbp1 -
                          pap3 * pbp2 + 2.0 * papb * pbp2)) *
           RLLR +
       (m1s * pap2 -
        m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb + 2.0 * ivu3 * pap3 * pbp1) +
        2.0 *
            (ivu3 * (-2.0 * p1p2 * p1p3 + m1s * p2p3) * papb +
             2.0 * ivu3 * (p1p3 * pap2 * pbp1 +
                           (pap3 - papb + pbp1) * (p1p2 * papb - pap2 * pbp1)) -
             ivu3 * (2.0 * p1p3 * pap2 + m1s * pap3 - 2.0 * p1p2 * pap3 +
                     2.0 * p1p2 * papb - 2.0 * pap2 * pbp1) *
                 pbp2 +
             pap1 * (-pap2 +
                     2.0 * ivu3 * (-(p2p3 * (papb + pbp2)) +
                                   pbp2 * (p1p3 + papb - pbp1 + pbp2))))) *
           RRRR);

  return real(MtuG);
}

// real quark emission without on-shell contributions
double FI::MttQ_GLGA_wos() {

  return real(
      4.0 * cF * ivt1s1 * NC *
      (-(cA * iv23s2 *
         (ivb1 *
              (LLLL * (p1p2 * (pap1 - papb) + pap2 * pbp1 -
                       p2p3 * (pap1 + 4.0 * ivu3 * p1p3 * pap1 -
                               2.0 * ivu3 * p1p3 * papb -
                               2.0 * ivu3 * pap1 * pbp1 +
                               2.0 * ivu3 * pap3 * pbp1) +
                       2.0 * ivu3 * (p1p3 * pap1 - p1p3 * papb + pap3 * pbp1) *
                           pbp2) +
               (p1p2 * pap1 -
                p2p3 * (pap1 + 4.0 * ivu3 * p1p3 * pap1 -
                        2.0 * ivu3 * p1p3 * papb - 2.0 * ivu3 * pap1 * pbp1 +
                        2.0 * ivu3 * pap3 * pbp1) +
                2.0 * ivu3 * p1p3 * pap1 * pbp2) *
                   (LRLR + RLRL) +
               (p1p2 * (pap1 - papb) + pap2 * pbp1 -
                p2p3 * (pap1 + 4.0 * ivu3 * p1p3 * pap1 -
                        2.0 * ivu3 * p1p3 * papb - 2.0 * ivu3 * pap1 * pbp1 +
                        2.0 * ivu3 * pap3 * pbp1) +
                2.0 * ivu3 * (p1p3 * pap1 - p1p3 * papb + pap3 * pbp1) * pbp2) *
                   RRRR) +
          ivs *
              (LLLL * (-(p1p2 * papb) + pap2 * (pap1 + pbp1) +
                       p2p3 * (pap1 - 4.0 * ivu3 * pap1 * pap3 +
                               2.0 * ivu3 * p1p3 * papb +
                               2.0 * ivu3 * pap1 * papb -
                               2.0 * ivu3 * pap3 * pbp1) +
                       2.0 * ivu3 * (-(p1p3 * papb) + pap3 * (pap1 + pbp1)) *
                           pbp2) -
               (-(p2p3 *
                  (pap1 - 4.0 * ivu3 * pap1 * pap3 + 2.0 * ivu3 * p1p3 * papb +
                   2.0 * ivu3 * pap1 * papb - 2.0 * ivu3 * pap3 * pbp1)) -
                pap1 * (pap2 + 2.0 * ivu3 * pap3 * pbp2)) *
                   (LRLR + RLRL) +
               (-(p1p2 * papb) + pap2 * (pap1 + pbp1) +
                p2p3 * (pap1 - 4.0 * ivu3 * pap1 * pap3 +
                        2.0 * ivu3 * p1p3 * papb + 2.0 * ivu3 * pap1 * papb -
                        2.0 * ivu3 * pap3 * pbp1) +
                2.0 * ivu3 * (-(p1p3 * papb) + pap3 * (pap1 + pbp1)) * pbp2) *
                   RRRR))) +
       cF * (-2.0 * ivt1s2 * ivu3 * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR) +
             iv23s2 * (ivt1s2 * pap1 *
                           (m2s - p1p2 + pap2 -
                            2.0 * ivu3 * p2p3 * (-2.0 * p1p3 + 2.0 * p2p3 +
                                                 2.0 * pap3 - papb + pbp1) +
                            2.0 * ivu3 * (-p1p3 + 2.0 * p2p3 + pap3) * pbp2) *
                           (LLLL + LRLR + RLRL + RRRR) +
                       2.0 * ivs *
                           (LLLL * (-(p1p2 * papb) + pap2 * (pap1 + pbp1) +
                                    p2p3 * (pap1 - 4.0 * ivu3 * pap1 * pap3 +
                                            2.0 * ivu3 * p1p3 * papb +
                                            2.0 * ivu3 * pap1 * papb -
                                            2.0 * ivu3 * pap3 * pbp1) +
                                    2.0 * ivu3 * (-(p1p3 * papb) +
                                                  pap3 * (pap1 + pbp1)) *
                                        pbp2) -
                            (-(p2p3 * (pap1 - 4.0 * ivu3 * pap1 * pap3 +
                                       2.0 * ivu3 * p1p3 * papb +
                                       2.0 * ivu3 * pap1 * papb -
                                       2.0 * ivu3 * pap3 * pbp1)) -
                             pap1 * (pap2 + 2.0 * ivu3 * pap3 * pbp2)) *
                                (LRLR + RLRL) +
                            (-(p1p2 * papb) + pap2 * (pap1 + pbp1) +
                             p2p3 * (pap1 - 4.0 * ivu3 * pap1 * pap3 +
                                     2.0 * ivu3 * p1p3 * papb +
                                     2.0 * ivu3 * pap1 * papb -
                                     2.0 * ivu3 * pap3 * pbp1) +
                             2.0 * ivu3 *
                                 (-(p1p3 * papb) + pap3 * (pap1 + pbp1)) *
                                 pbp2) *
                                RRRR)))));
}

// real quark emission without on-shell contributions
double FI::MuuQ_GLGA_wos() {
  return real(
      2.0 * cF * ivu2s2 * NC *
      (2.0 * cF *
           (-2.0 * ivu2s1 * ivu3 * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR) +
            iv13s1 * (ivu2s1 * pap2 *
                          (m1s - p1p2 + pap1 +
                           2.0 * ivu3 *
                               (-2.0 * pow(p1p3, 2) + (-p2p3 + pap3) * pbp1 +
                                p1p3 * (2.0 * p2p3 - 2.0 * pap3 + papb +
                                        2.0 * pbp1 - pbp2))) *
                          (LLLL + LRLR + RLRL + RRRR) +
                      2.0 * ivs *
                          (LLLL * (-(p1p2 * papb) + pap1 * (pap2 + pbp2) +
                                   p1p3 * (pap2 - 4.0 * ivu3 * pap2 * pap3 +
                                           2.0 * ivu3 * p2p3 * papb +
                                           2.0 * ivu3 * pap2 * papb -
                                           2.0 * ivu3 * pap3 * pbp2) +
                                   2.0 * ivu3 * pbp1 * (-(p2p3 * papb) +
                                                        pap3 * (pap2 + pbp2))) +
                           (pap2 * (pap1 + 2.0 * ivu3 * pap3 * pbp1) +
                            p1p3 * (pap2 - 4.0 * ivu3 * pap2 * pap3 +
                                    2.0 * ivu3 * p2p3 * papb +
                                    2.0 * ivu3 * pap2 * papb -
                                    2.0 * ivu3 * pap3 * pbp2)) *
                               (LRLR + RLRL) +
                           (-(p1p2 * papb) + pap1 * (pap2 + pbp2) +
                            p1p3 * (pap2 - 4.0 * ivu3 * pap2 * pap3 +
                                    2.0 * ivu3 * p2p3 * papb +
                                    2.0 * ivu3 * pap2 * papb -
                                    2.0 * ivu3 * pap3 * pbp2) +
                            2.0 * ivu3 * pbp1 *
                                (-(p2p3 * papb) + pap3 * (pap2 + pbp2))) *
                               RRRR))) -
       cA *
           (2.0 * ivb1 * ivu2s1 * pap2 *
                (-2.0 * p1p3 + 2.0 * pbp1 + 4.0 * ivu3 * p1p3 * (-p1p3 + pbp1) -
                 4.0 * ivb1 * pbp1 * pbp3 +
                 m1s * (1.0 + 2.0 * iv13s2 * p1p3 + 4.0 * ivb1 * p1p3 -
                        (iv13s2 + 4.0 * ivb1) * pbp3) +
                 iv13s2 * ((-p2p3 + pap3) * pbp1 +
                           p1p3 * (-2.0 * p1p2 + 2.0 * p1p3 + 2.0 * pap1 -
                                   papb + pbp2) +
                           (p1p2 - 2.0 * p1p3 - pap1) * pbp3)) *
                (LLLL + LRLR + RLRL + RRRR) +
            iv13s1 *
                (ivu2s1 * pap2 *
                     (m1s - p1p2 + pap1 +
                      2.0 * ivu3 *
                          (-2.0 * pow(p1p3, 2) + (-p2p3 + pap3) * pbp1 +
                           p1p3 * (2.0 * p2p3 - 2.0 * pap3 + papb + 2.0 * pbp1 -
                                   pbp2))) *
                     (LLLL + LRLR + RLRL + RRRR) +
                 2.0 * ivs *
                     (LLLL * (-(p1p2 * papb) +
                              p1p3 * (-2.0 * ivb1 * p1p2 * papb +
                                      2.0 * ivu3 * p2p3 * papb +
                                      pap2 * (1.0 - 4.0 * ivu3 * pap3 +
                                              2.0 * ivu3 * papb +
                                              2.0 * ivb1 *
                                                  (2.0 * pap1 - papb + pbp1)) +
                                      2.0 * ivb1 * pap1 * pbp2 -
                                      2.0 * ivu3 * pap3 * pbp2) +
                              2.0 * (ivb1 + ivu3) * pbp1 *
                                  (-(p2p3 * papb) + pap3 * (pap2 + pbp2)) +
                              2.0 * ivb1 * p1p2 * papb * pbp3 -
                              pap1 * (pap2 + pbp2) *
                                  (-1.0 + 2.0 * ivb1 * pbp3)) +
                      LRLR * (p1p3 * (-2.0 * ivb1 * p1p2 * papb +
                                      2.0 * ivu3 * p2p3 * papb +
                                      pap2 * (1.0 - 4.0 * ivu3 * pap3 +
                                              2.0 * ivu3 * papb +
                                              2.0 * ivb1 *
                                                  (2.0 * pap1 - papb + pbp1)) +
                                      2.0 * ivb1 * pap1 * pbp2 -
                                      2.0 * ivu3 * pap3 * pbp2) +
                              pap2 * (pap1 + 2.0 * (ivb1 + ivu3) * pap3 * pbp1 -
                                      2.0 * ivb1 * pap1 * pbp3)) +
                      (p1p3 *
                           (-2.0 * ivb1 * p1p2 * papb +
                            2.0 * ivu3 * p2p3 * papb +
                            pap2 *
                                (1.0 - 4.0 * ivu3 * pap3 + 2.0 * ivu3 * papb +
                                 2.0 * ivb1 * (2.0 * pap1 - papb + pbp1)) +
                            2.0 * ivb1 * pap1 * pbp2 -
                            2.0 * ivu3 * pap3 * pbp2) +
                       pap2 * (pap1 + 2.0 * (ivb1 + ivu3) * pap3 * pbp1 -
                               2.0 * ivb1 * pap1 * pbp3)) *
                          RLRL +
                      (-(p1p2 * papb) +
                       p1p3 *
                           (-2.0 * ivb1 * p1p2 * papb +
                            2.0 * ivu3 * p2p3 * papb +
                            pap2 *
                                (1.0 - 4.0 * ivu3 * pap3 + 2.0 * ivu3 * papb +
                                 2.0 * ivb1 * (2.0 * pap1 - papb + pbp1)) +
                            2.0 * ivb1 * pap1 * pbp2 -
                            2.0 * ivu3 * pap3 * pbp2) +
                       2.0 * (ivb1 + ivu3) * pbp1 *
                           (-(p2p3 * papb) + pap3 * (pap2 + pbp2)) +
                       2.0 * ivb1 * p1p2 * papb * pbp3 -
                       pap1 * (pap2 + pbp2) * (-1.0 + 2.0 * ivb1 * pbp3)) *
                          RRRR)))));
}

// real antiquark emission without on-shell contributions
double FI::MttQB_GLGA_wos() {

  return real(
      -2.0 * cF * NC *
      (-2.0 * cF * ivt2s1 *
           (-2.0 * ivt2s2 * ivt3 * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR) +
            iv13s2 *
                (ivt2s2 * pbp2 *
                     (m1s - p1p2 + pbp1 +
                      2.0 * ivt3 * (-2.0 * pow(p1p3, 2) +
                                    p1p3 * (2.0 * p2p3 + 2.0 * pap1 - pap2 +
                                            papb - 2.0 * pbp3) +
                                    pap1 * (-p2p3 + pbp3))) *
                     (LLLL + LRLR + RLRL + RRRR) +
                 2.0 * ivs *
                     (LLLL *
                          (-(p1p2 * papb) + pap2 * pbp1 + (p1p3 + pbp1) * pbp2 +
                           2.0 * ivt3 * (-(p2p3 * pap1 * papb) +
                                         p1p3 * papb * (p2p3 + pbp2) +
                                         pap1 * (pap2 + pbp2) * pbp3 -
                                         p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) +
                      ((p1p3 + pbp1) * pbp2 +
                       2.0 * ivt3 *
                           (p1p3 * papb * (p2p3 + pbp2) + pap1 * pbp2 * pbp3 -
                            p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) *
                          (LRLR + RLRL) +
                      (-(p1p2 * papb) + pap2 * pbp1 + (p1p3 + pbp1) * pbp2 +
                       2.0 * ivt3 * (-(p2p3 * pap1 * papb) +
                                     p1p3 * papb * (p2p3 + pbp2) +
                                     pap1 * (pap2 + pbp2) * pbp3 -
                                     p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) *
                          RRRR))) +
       cA *
           (iv13s2 * ivt2s1 *
                (ivt2s2 * pbp2 *
                     (m1s - p1p2 + pbp1 +
                      2.0 * ivt3 * (-2.0 * pow(p1p3, 2) +
                                    p1p3 * (2.0 * p2p3 + 2.0 * pap1 - pap2 +
                                            papb - 2.0 * pbp3) +
                                    pap1 * (-p2p3 + pbp3))) *
                     (LLLL + LRLR + RLRL + RRRR) +
                 2.0 * ivs *
                     (LLLL *
                          (-(p1p2 * papb) + pap2 * pbp1 + (p1p3 + pbp1) * pbp2 +
                           2.0 * ivt3 * (-(p2p3 * pap1 * papb) +
                                         p1p3 * papb * (p2p3 + pbp2) +
                                         pap1 * (pap2 + pbp2) * pbp3 -
                                         p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) +
                      ((p1p3 + pbp1) * pbp2 +
                       2.0 * ivt3 *
                           (p1p3 * papb * (p2p3 + pbp2) + pap1 * pbp2 * pbp3 -
                            p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) *
                          (LRLR + RLRL) +
                      (-(p1p2 * papb) + pap2 * pbp1 + (p1p3 + pbp1) * pbp2 +
                       2.0 * ivt3 * (-(p2p3 * pap1 * papb) +
                                     p1p3 * papb * (p2p3 + pbp2) +
                                     pap1 * (pap2 + pbp2) * pbp3 -
                                     p1p3 * (pap2 + 2.0 * pbp2) * pbp3)) *
                          RRRR)) +
            2.0 * iva1 * ivt2s2 *
                (ivt2s1 * (m1s - 2.0 * p1p3 + 4.0 * iva1 * m1s * p1p3 -
                           4.0 * ivt3 * pow(p1p3, 2) + 2.0 * pap1 +
                           4.0 * ivt3 * p1p3 * pap1 -
                           4.0 * iva1 * (m1s + pap1) * pap3) *
                     pbp2 * (LLLL + LRLR + RLRL + RRRR) +
                 iv13s1 *
                     (ivt2s1 * pbp2 *
                          (m1s * (2.0 * p1p3 - pap3) +
                           p1p2 * (-2.0 * p1p3 + pap3) +
                           p1p3 * (2.0 * p1p3 + pap2 - 2.0 * pap3 - papb) +
                           (2.0 * p1p3 - pap3) * pbp1 + pap1 * (-p2p3 + pbp3)) *
                          (LLLL + LRLR + RLRL + RRRR) -
                      2.0 * ivs *
                          (LLLL * (p2p3 * pap1 * papb +
                                   p1p2 * (p1p3 - pap3) * papb +
                                   pap2 * (-p1p3 + pap3) * pbp1 +
                                   (p1p3 * (-pap1 + papb - 2.0 * pbp1) +
                                    pap3 * pbp1) *
                                       pbp2 -
                                   pap1 * (pap2 + pbp2) * pbp3) +
                           (p1p2 * p1p3 * papb -
                            p1p3 * (pap2 * pbp1 +
                                    (pap1 - papb + 2.0 * pbp1) * pbp2) +
                            pbp2 * (pap3 * pbp1 - pap1 * pbp3)) *
                               (LRLR + RLRL) +
                           (p2p3 * pap1 * papb + p1p2 * (p1p3 - pap3) * papb +
                            pap2 * (-p1p3 + pap3) * pbp1 +
                            (p1p3 * (-pap1 + papb - 2.0 * pbp1) + pap3 * pbp1) *
                                pbp2 -
                            pap1 * (pap2 + pbp2) * pbp3) *
                               RRRR))))));
}

// real antiquark emission without on-shell contributions
double FI::MuuQB_GLGA_wos() {

  return real(
      -4.0 * cF * ivu1s1 * NC *
      (cA * iv23s2 *
           (ivs * (LLLL * (-(p1p2 * papb) + p2p3 * pbp1 + (pap1 + pbp1) * pbp2 +
                           2.0 * ivt3 * (p1p3 * (p2p3 - pap2) * papb +
                                         p2p3 * papb * pbp1 +
                                         pap2 * (pap1 + pbp1) * pbp3 -
                                         p2p3 * (pap1 + 2.0 * pbp1) * pbp3)) +
                   (pbp1 * (p2p3 + pbp2) +
                    2.0 * ivt3 *
                        (p2p3 * papb * (p1p3 + pbp1) +
                         (pap2 * pbp1 - p2p3 * (pap1 + 2.0 * pbp1)) * pbp3)) *
                       (LRLR + RLRL) +
                   (-(p1p2 * papb) + p2p3 * pbp1 + (pap1 + pbp1) * pbp2 +
                    2.0 * ivt3 *
                        (p1p3 * (p2p3 - pap2) * papb + p2p3 * papb * pbp1 +
                         pap2 * (pap1 + pbp1) * pbp3 -
                         p2p3 * (pap1 + 2.0 * pbp1) * pbp3)) *
                       RRRR) +
            iva1 *
                (LLLL *
                     (-(p2p3 * pbp1) + p1p2 * (-papb + pbp1) + pap1 * pbp2 +
                      2.0 * ivt3 *
                          (p1p3 * (p2p3 * papb - pap2 * papb -
                                   2.0 * p2p3 * pbp1 + pap2 * pbp1) +
                           pap1 * (p2p3 * pbp1 - p2p3 * pbp3 + pap2 * pbp3))) +
                 ((p1p2 - p2p3) * pbp1 +
                  2.0 * ivt3 *
                      (p1p3 * (p2p3 * papb - 2.0 * p2p3 * pbp1 + pap2 * pbp1) +
                       p2p3 * pap1 * (pbp1 - pbp3))) *
                     (LRLR + RLRL) +
                 (-(p2p3 * pbp1) + p1p2 * (-papb + pbp1) + pap1 * pbp2 +
                  2.0 * ivt3 *
                      (p1p3 * (p2p3 * papb - pap2 * papb - 2.0 * p2p3 * pbp1 +
                               pap2 * pbp1) +
                       pap1 * (p2p3 * pbp1 - p2p3 * pbp3 + pap2 * pbp3))) *
                     RRRR)) +
       cF *
           (2.0 * ivt3 * ivu1s2 * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR) +
            iv23s2 *
                (-(ivu1s2 * pbp1 *
                   (m2s - p1p2 + pbp2 +
                    2.0 * ivt3 *
                        (-2.0 * pow(p2p3, 2) + p1p3 * (2.0 * p2p3 - pap2) +
                         p2p3 * (-pap1 + 2.0 * pap2 + papb - 2.0 * pbp3) +
                         pap2 * pbp3)) *
                   (LLLL + LRLR + RLRL + RRRR)) -
                 2.0 * ivs *
                     (LLLL *
                          (-(p1p2 * papb) + p2p3 * pbp1 + (pap1 + pbp1) * pbp2 +
                           2.0 * ivt3 * (p1p3 * (p2p3 - pap2) * papb +
                                         p2p3 * papb * pbp1 +
                                         pap2 * (pap1 + pbp1) * pbp3 -
                                         p2p3 * (pap1 + 2.0 * pbp1) * pbp3)) +
                      (pbp1 * (p2p3 + pbp2) +
                       2.0 * ivt3 *
                           (p2p3 * papb * (p1p3 + pbp1) +
                            (pap2 * pbp1 - p2p3 * (pap1 + 2.0 * pbp1)) *
                                pbp3)) *
                          (LRLR + RLRL) +
                      (-(p1p2 * papb) + p2p3 * pbp1 + (pap1 + pbp1) * pbp2 +
                       2.0 * ivt3 *
                           (p1p3 * (p2p3 - pap2) * papb + p2p3 * papb * pbp1 +
                            pap2 * (pap1 + pbp1) * pbp3 -
                            p2p3 * (pap1 + 2.0 * pbp1) * pbp3)) *
                          RRRR)))));
}

// resonant diagrams
double FI::MttQres_GLGA() {

  return real(
      -4.0 * cF * iv23s1 * iv23s2 * NC * p2p3 *
      (cA * (4.0 * pow(ivb1, 2) * (m1s * (pap1 - papb) - papb * pbp1) +
             ivs * ivt1s1 *
                 (2.0 * pow(pap1, 2) + (-m1s + p1p2 + p1p3) * papb -
                  (pap2 + pap3) * pbp1 -
                  pap1 * (2.0 * pap2 + 2.0 * pap3 - 2.0 * pbp1 + pbp2 + pbp3)) +
             ivb1 * (2.0 * ivs *
                         (-(m1s * papb) + 2.0 * (pap1 - papb) * (pap1 + pbp1)) +
                     ivt1s1 *
                         (2.0 * m1s * pap1 - 2.0 * p1p2 * pap1 -
                          2.0 * p1p3 * pap1 - 2.0 * pow(pap1, 2) - m1s * papb +
                          p1p2 * papb + p1p3 * papb + 2.0 * pap1 * papb -
                          (pap2 + pap3) * pbp1 + pap1 * pbp2 + pap1 * pbp3))) +
       cF *
           (-4.0 * pow(ivs, 2) * papb * pbp1 +
            ivt1s1 * (ivt1s2 * pap1 *
                          (m1s + m2s -
                           2.0 * (p1p2 + p1p3 - p2p3 + pap1 - pap2 - pap3)) +
                      2.0 * ivs * (-2.0 * pow(pap1, 2) + m1s * papb -
                                   (p1p2 + p1p3) * papb + (pap2 + pap3) * pbp1 +
                                   pap1 * (2.0 * pap2 + 2.0 * pap3 -
                                           2.0 * pbp1 + pbp2 + pbp3))))) *
      (LLLL + LRLR + RLRL + RRRR));
}

double FI::MuuQres_GLGA() {

  return real(
      -4.0 * pow(cF, 2) * iv13s1 * iv13s2 * NC * p1p3 *
      (ivu2s1 * ivu2s2 * pap2 *
           (m1s + m2s - 2.0 * (p1p2 - p1p3 + p2p3 - pap1 + pap2 - pap3)) +
       2.0 * ivs *
           (-2.0 * ivs * papb * pbp2 +
            ivu2s2 * (-2.0 * pow(pap2, 2) + (m2s - p1p2 - p2p3) * papb +
                      pap3 * pbp2 + pap1 * (2.0 * pap2 + pbp2) +
                      pap2 * (2.0 * pap3 + pbp1 - 2.0 * pbp2 + pbp3)))) *
      (LLLL + LRLR + RLRL + RRRR));
}

// real antiquark resonant diagrams
double FI::MttQBres_GLGA() {

  return real(
      4.0 * pow(cF, 2) * iv13s1 * iv13s2 * NC * p1p3 *
      (4.0 * pow(ivs, 2) * pap2 * papb +
       ivt2s1 * ivt2s2 * pbp2 *
           (-m1s - m2s + 2.0 * (p1p2 - p1p3 + p2p3 - pbp1 + pbp2 - pbp3)) -
       2.0 * ivs * ivt2s1 *
           (m2s * papb - (p1p2 + p2p3) * papb +
            pap2 * (pbp1 - 2.0 * pbp2 + pbp3) +
            pbp2 * (pap1 + pap3 + 2.0 * (pbp1 - pbp2 + pbp3)))) *
      (LLLL + LRLR + RLRL + RRRR));
}

double FI::MuuQBres_GLGA() {

  return real(
      4.0 * cF * iv23s1 * iv23s2 * NC * p2p3 *
      (cF *
           (4.0 * pow(ivs, 2) * pap1 * papb +
            ivu1s1 * ivu1s2 * pbp1 *
                (-m1s - m2s + 2.0 * (p1p2 + p1p3 - p2p3 + pbp1 - pbp2 - pbp3)) -
            2.0 * ivs * ivu1s1 *
                (m1s * papb - (p1p2 + p1p3) * papb +
                 pap1 * (-2.0 * pbp1 + pbp2 + pbp3) +
                 pbp1 * (pap2 + pap3 + 2.0 * (-pbp1 + pbp2 + pbp3)))) +
       cA * (4.0 * pow(iva1, 2) * ((m1s + pap1) * papb - m1s * pbp1) +
             iva1 * (2.0 * ivs *
                         (m1s * papb + 2.0 * (papb - pbp1) * (pap1 + pbp1)) +
                     ivu1s1 * ((m1s - p1p2 - p1p3) * papb -
                               (2.0 * m1s - 2.0 * p1p2 - 2.0 * p1p3 + pap2 +
                                pap3 + 2.0 * papb) *
                                   pbp1 +
                               2.0 * pow(pbp1, 2) + pap1 * (pbp2 + pbp3))) +
             ivs * ivu1s1 *
                 (m1s * papb - (p1p2 + p1p3) * papb +
                  pap1 * (-2.0 * pbp1 + pbp2 + pbp3) +
                  pbp1 * (pap2 + pap3 + 2.0 * (-pbp1 + pbp2 + pbp3))))) *
      (LLLL + LRLR + RLRL + RRRR));
}

// real quark emission
double FI::MtuQ_GLGA() {

  complex<double> MtuQ;

  // - (MtuGS + MtuG1 + MtuGP +  MtuGB)
  MtuQ =
      2.0 * cA * cF * iv23s1 * ivb1 * NC *
      (-2.0 * iv13s2 * ivs *
           (LLLL * (2.0 * p1p2 * pap1 * pap3 - 2.0 * p1p2 * pap3 * papb +
                    p2p3 * (m1s * papb - 2.0 * (pap1 - papb) * (pap1 + pbp1)) +
                    m1s * pap3 * pbp2 -
                    2.0 * p1p3 * (pap1 - papb) * (pap2 + pbp2) +
                    2.0 * p1p2 * pap1 * pbp3 - m1s * pap2 * pbp3 -
                    2.0 * p1p2 * papb * pbp3) +
            m1 * m2 * (p1p3 * papb + 2.0 * pap3 * papb - pap3 * pbp1 +
                       2.0 * papb * pbp3 - pap1 * (2.0 * pap3 + pbp3)) *
                (LRRL + RLLR) +
            (2.0 * p1p2 * pap1 * pap3 - 2.0 * p1p2 * pap3 * papb +
             p2p3 * (m1s * papb - 2.0 * (pap1 - papb) * (pap1 + pbp1)) +
             m1s * pap3 * pbp2 - 2.0 * p1p3 * (pap1 - papb) * (pap2 + pbp2) +
             2.0 * p1p2 * pap1 * pbp3 - m1s * pap2 * pbp3 -
             2.0 * p1p2 * papb * pbp3) *
                RRRR) +
       ivu2s2 *
           (8.0 * ivb1 * LLLL * m1s * p2p3 * pap1 -
            4.0 * ivu3 * LLLL * p1p3 * p2p3 * pap1 + LLLL * m1s * pap2 -
            2.0 * LLLL * p1p3 * pap2 + 8.0 * ivb1 * LLLL * m1s * p1p3 * pap2 -
            4.0 * ivu3 * LLLL * pow2(p1p3) * pap2 -
            8.0 * ivb1 * LLLL * m1s * p1p2 * pap3 +
            4.0 * ivu3 * LLLL * p1p2 * p1p3 * pap3 -
            8.0 * ivb1 * LLLL * m1s * p2p3 * papb -
            2.0 * ivu3 * LLLL * m1s * p2p3 * papb +
            4.0 * ivu3 * LLLL * p1p3 * p2p3 * papb +
            4.0 * ivu3 * LLLL * p2p3 * pap1 * pbp1 + 2.0 * LLLL * pap2 * pbp1 +
            4.0 * ivu3 * LLLL * p1p3 * pap2 * pbp1 -
            4.0 * ivu3 * LLLL * p1p2 * pap3 * pbp1 -
            8.0 * ivb1 * LLLL * p2p3 * papb * pbp1 -
            4.0 * ivu3 * LLLL * p2p3 * papb * pbp1 +
            8.0 * ivb1 * LLLL * m1s * pap3 * pbp2 +
            2.0 * ivu3 * LLLL * m1s * pap3 * pbp2 -
            4.0 * ivu3 * LLLL * p1p3 * pap3 * pbp2 +
            8.0 * ivb1 * LLLL * pap3 * pbp1 * pbp2 +
            4.0 * ivu3 * LLLL * pap3 * pbp1 * pbp2 -
            8.0 * ivb1 * LLLL * m1s * pap2 * pbp3 -
            8.0 * ivb1 * LLLL * pap2 * pbp1 * pbp3 +
            LRRL * m1 * m2 *
                (2.0 * (4.0 * ivb1 * m1s - 2.0 * ivu3 * p1p3 +
                        iv13s2 * (m1s - p1p2 + p1p3)) *
                     pap3 +
                 2.0 * ivu3 * p1p3 * papb +
                 2.0 * (-4.0 * ivb1 + ivu3) * pap3 * pbp1 -
                 iv13s2 * (p1p3 * papb - p2p3 * papb +
                           pap3 * (2.0 * papb + pbp1 - pbp2) + pap2 * pbp3) +
                 pap1 * (1.0 + iv13s2 * (2.0 * pap3 + pbp3))) +
            m1 * m2 * pap1 * RLLR + 8.0 * ivb1 * m1 * m1s * m2 * pap3 * RLLR -
            4.0 * ivu3 * m1 * m2 * p1p3 * pap3 * RLLR +
            2.0 * ivu3 * m1 * m2 * p1p3 * papb * RLLR -
            8.0 * ivb1 * m1 * m2 * pap3 * pbp1 * RLLR +
            2.0 * ivu3 * m1 * m2 * pap3 * pbp1 * RLLR +
            (pap2 * (m1s - 2.0 * p1p3 + 2.0 * pbp1) +
             ivu3 * (-4.0 * pow2(p1p3) * pap2 +
                     4.0 * p1p3 * (p1p2 * pap3 + p2p3 * (-pap1 + papb) +
                                   pap2 * pbp1 - pap3 * pbp2) +
                     4.0 * pbp1 * (p2p3 * pap1 - p1p2 * pap3 - p2p3 * papb +
                                   pap3 * pbp2) +
                     m1s * (-2.0 * p2p3 * papb + 2.0 * pap3 * pbp2)) +
             8.0 * ivb1 *
                 (m1s * p2p3 * pap1 + m1s * p1p3 * pap2 - m1s * p1p2 * pap3 -
                  m1s * p2p3 * papb - p2p3 * papb * pbp1 + m1s * pap3 * pbp2 +
                  pap3 * pbp1 * pbp2 - pap2 * (m1s + pbp1) * pbp3)) *
                RRRR +
            iv13s2 *
                (LLLL * (2.0 * p1p3 * p2p3 * pap1 + 2.0 * p2p3 * pow2(pap1) +
                         2.0 * pow2(p1p3) * pap2 + 2.0 * p1p3 * pap1 * pap2 +
                         2.0 * pow2(p1p2) * pap3 + m2s * p1p3 * papb -
                         2.0 * p1p3 * p2p3 * papb - 2.0 * p2p3 * pap1 * papb -
                         2.0 * p1p3 * pap2 * papb - 2.0 * p2p3 * pap2 * pbp1 +
                         m2s * pap3 * pbp1 + 2.0 * p2p3 * pap1 * pbp2 +
                         2.0 * p1p3 * pap3 * pbp2 - m2s * pap1 * pbp3 -
                         2.0 * p1p3 * pap2 * pbp3 +
                         m1s * (2.0 * p2p3 * pap1 + 2.0 * p1p3 * pap2 -
                                2.0 * p1p2 * pap3 - p2p3 * papb + pap3 * pbp2 -
                                pap2 * pbp3) -
                         2.0 * p1p2 *
                             (p2p3 * pap1 + p1p3 * (pap2 + pap3) +
                              pap3 * (pap1 - papb + pbp2) - pap2 * pbp3)) +
                 m1 * m2 *
                     (2.0 * m1s * pap3 - 2.0 * p1p2 * pap3 + 2.0 * p1p3 * pap3 +
                      2.0 * pap1 * pap3 - p1p3 * papb + p2p3 * papb -
                      2.0 * pap3 * papb - pap3 * pbp1 + pap3 * pbp2 +
                      pap1 * pbp3 - pap2 * pbp3) *
                     RLLR +
                 (2.0 * p1p3 * p2p3 * pap1 + 2.0 * p2p3 * pow2(pap1) +
                  2.0 * pow2(p1p3) * pap2 + 2.0 * p1p3 * pap1 * pap2 +
                  2.0 * pow2(p1p2) * pap3 + m2s * p1p3 * papb -
                  2.0 * p1p3 * p2p3 * papb - 2.0 * p2p3 * pap1 * papb -
                  2.0 * p1p3 * pap2 * papb - 2.0 * p2p3 * pap2 * pbp1 +
                  m2s * pap3 * pbp1 + 2.0 * p2p3 * pap1 * pbp2 +
                  2.0 * p1p3 * pap3 * pbp2 - m2s * pap1 * pbp3 -
                  2.0 * p1p3 * pap2 * pbp3 +
                  m1s * (2.0 * p2p3 * pap1 + 2.0 * p1p3 * pap2 -
                         2.0 * p1p2 * pap3 - p2p3 * papb + pap3 * pbp2 -
                         pap2 * pbp3) -
                  2.0 * p1p2 * (p2p3 * pap1 + p1p3 * (pap2 + pap3) +
                                pap3 * (pap1 - papb + pbp2) - pap2 * pbp3)) *
                     RRRR)));

  // - (MtuBB + MtuB1 + MtuBS + MtuBP)
  MtuQ +=
      cF * ivt1s1 * NC *
      (-2.0 * cF *
           (-4.0 * ivu2s2 * ivu3 *
                (LRRL * m1 * m2 * papb +
                 LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                 m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                 pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
            iv13s2 *
                (4.0 * ivs *
                     (LLLL *
                          (-(p1p2 * pap3) +
                           p1p3 * (pap2 - 2.0 * ivu3 * pap2 * pap3 +
                                   2.0 * ivu3 * papb * (p2p3 - pbp2)) +
                           pap1 * (pap2 + 2.0 * ivu3 * p2p3 * (-pap3 + papb) +
                                   pbp2) +
                           2.0 * ivu3 * pap3 * (p1p2 * (pap3 - papb) +
                                                pbp1 * (-p2p3 + pap2 + pbp2))) -
                      m1 * m2 * pap3 * (-1.0 + 2.0 * ivu3 * (pap3 - papb)) *
                          (LRRL + RLLR) +
                      (-(p1p2 * pap3) +
                       p1p3 * (pap2 - 2.0 * ivu3 * pap2 * pap3 +
                               2.0 * ivu3 * papb * (p2p3 - pbp2)) +
                       pap1 *
                           (pap2 + 2.0 * ivu3 * p2p3 * (-pap3 + papb) + pbp2) +
                       2.0 * ivu3 * pap3 * (p1p2 * (pap3 - papb) +
                                            pbp1 * (-p2p3 + pap2 + pbp2))) *
                          RRRR) +
                 ivu2s2 *
                     (LRRL * m1 * m2 *
                          (pap1 - pap2 +
                           2.0 * ivu3 * (p2p3 * (2.0 * pap3 - papb) +
                                         p1p3 * (-2.0 * pap3 + papb) +
                                         pap3 * (-2.0 * pap3 + 2.0 * papb +
                                                 pbp1 - pbp2))) +
                      LLLL *
                          (m1s * pap2 -
                           m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb +
                                  2.0 * ivu3 * pap3 * pbp1) +
                           2.0 * (pap1 * pap2 +
                                  ivu3 *
                                      (2.0 * pow2(p2p3) * pap1 -
                                       2.0 * pow2(p1p3) * pap2 +
                                       pap3 * (2.0 * pap2 * pbp1 - m1s * pbp2 +
                                               2.0 * p1p2 *
                                                   (pap3 - papb + pbp2)) +
                                       p2p3 *
                                           (-2.0 * p1p2 * pap3 + m1s * papb -
                                            2.0 * pap1 * (pap3 - papb + pbp2)) +
                                       2.0 * p1p3 *
                                           (p2p3 * (-pap1 + pap2) +
                                            p1p2 * (pap3 - papb) + pap1 * pbp2 -
                                            pap2 * (pap3 - pbp1 + pbp2))))) +
                      m1 * m2 *
                          (pap1 - pap2 +
                           2.0 * ivu3 * (p2p3 * (2.0 * pap3 - papb) +
                                         p1p3 * (-2.0 * pap3 + papb) +
                                         pap3 * (-2.0 * pap3 + 2.0 * papb +
                                                 pbp1 - pbp2))) *
                          RLLR +
                      (m1s * pap2 -
                       m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb +
                              2.0 * ivu3 * pap3 * pbp1) +
                       2.0 *
                           (pap1 * pap2 +
                            ivu3 * (2.0 * pow2(p2p3) * pap1 -
                                    2.0 * pow2(p1p3) * pap2 +
                                    pap3 * (2.0 * pap2 * pbp1 - m1s * pbp2 +
                                            2.0 * p1p2 * (pap3 - papb + pbp2)) +
                                    p2p3 * (-2.0 * p1p2 * pap3 + m1s * papb -
                                            2.0 * pap1 * (pap3 - papb + pbp2)) +
                                    2.0 * p1p3 *
                                        (p2p3 * (-pap1 + pap2) +
                                         p1p2 * (pap3 - papb) + pap1 * pbp2 -
                                         pap2 * (pap3 - pbp1 + pbp2))))) *
                          RRRR))) +
       cA *
           (2.0 * ivb1 * ivu2s2 *
                (LRRL * m1 * m2 *
                     (pap1 +
                      2.0 * (papb + ivu3 * p1p3 * papb +
                             pap3 * (-1.0 - 2.0 * ivu3 * p1p3 + ivu3 * pbp1))) +
                 LLLL * (-2.0 * p2p3 * (pap1 + 2.0 * ivu3 * p1p3 * pap1 -
                                        ivu3 * m1s * papb) -
                         2.0 * (1.0 + 2.0 * ivu3 * p1p3) *
                             (p1p3 * pap2 - p1p2 * pap3 + p1p2 * papb -
                              pap2 * pbp1 - pap1 * pbp2) +
                         m1s * (pap2 - 2.0 * ivu3 * pap3 * pbp2)) +
                 m1 * m2 *
                     (pap1 +
                      2.0 * (papb + ivu3 * p1p3 * papb +
                             pap3 * (-1.0 - 2.0 * ivu3 * p1p3 + ivu3 * pbp1))) *
                     RLLR +
                 (-2.0 * p2p3 *
                      (pap1 + 2.0 * ivu3 * p1p3 * pap1 - ivu3 * m1s * papb) -
                  2.0 * (1.0 + 2.0 * ivu3 * p1p3) *
                      (p1p3 * pap2 - p1p2 * pap3 + p1p2 * papb - pap2 * pbp1 -
                       pap1 * pbp2) +
                  m1s * (pap2 - 2.0 * ivu3 * pap3 * pbp2)) *
                     RRRR) +
            iv13s2 *
                (4.0 * ivs *
                     (LLLL * (-(p1p2 * pap3) +
                              p1p3 * (pap2 - 2.0 * ivu3 * pap2 * pap3 +
                                      2.0 * ivu3 * papb * (p2p3 - pbp2)) +
                              pap1 *
                                  (pap2 + 2.0 * ivu3 * p2p3 * (-pap3 + papb) +
                                   pbp2) +
                              2.0 * ivu3 * pap3 *
                                  (p1p2 * (pap3 - papb) +
                                   pbp1 * (-p2p3 + pap2 + pbp2))) -
                      m1 * m2 * pap3 * (-1.0 + 2.0 * ivu3 * (pap3 - papb)) *
                          (LRRL + RLLR) +
                      (-(p1p2 * pap3) +
                       p1p3 * (pap2 - 2.0 * ivu3 * pap2 * pap3 +
                               2.0 * ivu3 * papb * (p2p3 - pbp2)) +
                       pap1 *
                           (pap2 + 2.0 * ivu3 * p2p3 * (-pap3 + papb) + pbp2) +
                       2.0 * ivu3 * pap3 * (p1p2 * (pap3 - papb) +
                                            pbp1 * (-p2p3 + pap2 + pbp2))) *
                          RRRR) +
                 ivu2s2 *
                     (LRRL * m1 * m2 *
                          (pap1 - pap2 +
                           2.0 * ivu3 * (p2p3 * (2.0 * pap3 - papb) +
                                         p1p3 * (-2.0 * pap3 + papb) +
                                         pap3 * (-2.0 * pap3 + 2.0 * papb +
                                                 pbp1 - pbp2))) +
                      LLLL *
                          (m1s * pap2 -
                           m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb +
                                  2.0 * ivu3 * pap3 * pbp1) +
                           2.0 * (pap1 * pap2 +
                                  ivu3 *
                                      (2.0 * pow2(p2p3) * pap1 -
                                       2.0 * pow2(p1p3) * pap2 +
                                       pap3 * (2.0 * pap2 * pbp1 - m1s * pbp2 +
                                               2.0 * p1p2 *
                                                   (pap3 - papb + pbp2)) +
                                       p2p3 *
                                           (-2.0 * p1p2 * pap3 + m1s * papb -
                                            2.0 * pap1 * (pap3 - papb + pbp2)) +
                                       2.0 * p1p3 *
                                           (p2p3 * (-pap1 + pap2) +
                                            p1p2 * (pap3 - papb) + pap1 * pbp2 -
                                            pap2 * (pap3 - pbp1 + pbp2))))) +
                      m1 * m2 *
                          (pap1 - pap2 +
                           2.0 * ivu3 * (p2p3 * (2.0 * pap3 - papb) +
                                         p1p3 * (-2.0 * pap3 + papb) +
                                         pap3 * (-2.0 * pap3 + 2.0 * papb +
                                                 pbp1 - pbp2))) *
                          RLLR +
                      (m1s * pap2 -
                       m2s * (pap1 - 2.0 * ivu3 * p1p3 * papb +
                              2.0 * ivu3 * pap3 * pbp1) +
                       2.0 *
                           (pap1 * pap2 +
                            ivu3 * (2.0 * pow2(p2p3) * pap1 -
                                    2.0 * pow2(p1p3) * pap2 +
                                    pap3 * (2.0 * pap2 * pbp1 - m1s * pbp2 +
                                            2.0 * p1p2 * (pap3 - papb + pbp2)) +
                                    p2p3 * (-2.0 * p1p2 * pap3 + m1s * papb -
                                            2.0 * pap1 * (pap3 - papb + pbp2)) +
                                    2.0 * p1p3 *
                                        (p2p3 * (-pap1 + pap2) +
                                         p1p2 * (pap3 - papb) + pap1 * pbp2 -
                                         pap2 * (pap3 - pbp1 + pbp2))))) *
                          RRRR))));

  // - (Mtu21 + Mtu2S + Mtu2P + Mtu2B)
  MtuQ +=
      2.0 * cF * iv23s1 * ivt1s1 * NC *
      (cF *
           (ivu2s2 *
                (LLLL * (-4.0 * ivu3 * p1p3 * p2p3 * pap1 +
                         4.0 * ivu3 * pow2(p2p3) * pap1 + m1s * pap2 -
                         4.0 * ivu3 * pow2(p1p3) * pap2 +
                         4.0 * ivu3 * p1p3 * p2p3 * pap2 - 2.0 * pap1 * pap2 +
                         4.0 * ivu3 * p1p2 * p1p3 * pap3 -
                         4.0 * ivu3 * p1p2 * p2p3 * pap3 +
                         4.0 * ivu3 * p2p3 * pap1 * pap3 +
                         4.0 * ivu3 * p1p3 * pap2 * pap3 -
                         4.0 * ivu3 * p1p2 * pow2(pap3) +
                         2.0 * iv13s2 *
                             (-(p2p3 * pap1) - p1p3 * pap2 + p1p2 * pap3) *
                             (m1s - 2.0 * (p1p2 + pap3)) -
                         2.0 * ivu3 * m1s * p2p3 * papb +
                         4.0 * ivu3 * p1p2 * p2p3 * papb -
                         4.0 * ivu3 * p1p3 * pap2 * papb +
                         4.0 * ivu3 * p1p2 * pap3 * papb +
                         4.0 * ivu3 * p2p3 * pap1 * pbp1 +
                         4.0 * ivu3 * p1p3 * pap2 * pbp1 -
                         4.0 * ivu3 * p2p3 * pap2 * pbp1 -
                         4.0 * ivu3 * p1p2 * pap3 * pbp1 -
                         m2s * (pap1 + 2.0 * iv13s2 * p2p3 * pap1 +
                                2.0 * iv13s2 * p1p3 * pap2 -
                                2.0 * iv13s2 * p1p2 * pap3 +
                                2.0 * ivu3 * p1p3 * papb -
                                2.0 * ivu3 * pap3 * pbp1) +
                         2.0 * ivu3 *
                             (m1s * pap3 - 2.0 * pap1 * (p2p3 + pap3)) * pbp2) +
                 LRRL * m1 * m2 *
                     (pap1 - pap2 -
                      2.0 * (iv13s2 * pap3 * (m1s + m2s - 2.0 * (p1p2 + pap3)) +
                             ivu3 * (p1p3 * (2.0 * pap3 - papb) +
                                     p2p3 * (-2.0 * pap3 + papb) +
                                     pap3 * (-2.0 * pap3 + 2.0 * papb - pbp1 +
                                             pbp2)))) +
                 m1 * m2 *
                     (pap1 - pap2 -
                      2.0 * (iv13s2 * pap3 * (m1s + m2s - 2.0 * (p1p2 + pap3)) +
                             ivu3 * (p1p3 * (2.0 * pap3 - papb) +
                                     p2p3 * (-2.0 * pap3 + papb) +
                                     pap3 * (-2.0 * pap3 + 2.0 * papb - pbp1 +
                                             pbp2)))) *
                     RLLR +
                 (-4.0 * ivu3 * p1p3 * p2p3 * pap1 +
                  4.0 * ivu3 * pow2(p2p3) * pap1 + m1s * pap2 -
                  4.0 * ivu3 * pow2(p1p3) * pap2 +
                  4.0 * ivu3 * p1p3 * p2p3 * pap2 - 2.0 * pap1 * pap2 +
                  4.0 * ivu3 * p1p2 * p1p3 * pap3 -
                  4.0 * ivu3 * p1p2 * p2p3 * pap3 +
                  4.0 * ivu3 * p2p3 * pap1 * pap3 +
                  4.0 * ivu3 * p1p3 * pap2 * pap3 -
                  4.0 * ivu3 * p1p2 * pow2(pap3) +
                  2.0 * iv13s2 * (-(p2p3 * pap1) - p1p3 * pap2 + p1p2 * pap3) *
                      (m1s - 2.0 * (p1p2 + pap3)) -
                  2.0 * ivu3 * m1s * p2p3 * papb +
                  4.0 * ivu3 * p1p2 * p2p3 * papb -
                  4.0 * ivu3 * p1p3 * pap2 * papb +
                  4.0 * ivu3 * p1p2 * pap3 * papb +
                  4.0 * ivu3 * p2p3 * pap1 * pbp1 +
                  4.0 * ivu3 * p1p3 * pap2 * pbp1 -
                  4.0 * ivu3 * p2p3 * pap2 * pbp1 -
                  4.0 * ivu3 * p1p2 * pap3 * pbp1 -
                  m2s *
                      (pap1 + 2.0 * iv13s2 * p2p3 * pap1 +
                       2.0 * iv13s2 * p1p3 * pap2 - 2.0 * iv13s2 * p1p2 * pap3 +
                       2.0 * ivu3 * p1p3 * papb - 2.0 * ivu3 * pap3 * pbp1) +
                  2.0 * ivu3 * (m1s * pap3 - 2.0 * pap1 * (p2p3 + pap3)) *
                      pbp2) *
                     RRRR) +
            2.0 * iv13s2 * ivs *
                (LLLL * (pap3 * (m2s * pbp1 +
                                 2.0 * p1p2 * (pap1 - pap2 - pap3 - pbp2) +
                                 m1s * pbp2) +
                         p2p3 * (-2.0 * pow2(pap1) + m1s * papb -
                                 2.0 * p1p3 * papb + 2.0 * pap3 * pbp1 +
                                 2.0 * pap1 * (pap2 + pap3 - pbp1 + pbp2)) -
                         (m2s * pap1 - 2.0 * p1p2 * pap1 + m1s * pap2 +
                          2.0 * p1p2 * pap3) *
                             pbp3 +
                         p1p3 * (-(m2s * papb) - 2.0 * pap1 * (pap2 + pbp2) +
                                 2.0 * pap2 * (pap2 + pap3 + pbp2 + pbp3))) -
                 m1 * m2 * ((-p1p3 + p2p3) * papb +
                            pap3 * (2.0 * pap1 - 2.0 * pap2 - 2.0 * pap3 +
                                    pbp1 - pbp2) +
                            (pap1 - pap2 - 2.0 * pap3) * pbp3) *
                     (LRRL + RLLR) +
                 (pap3 *
                      (m2s * pbp1 + 2.0 * p1p2 * (pap1 - pap2 - pap3 - pbp2) +
                       m1s * pbp2) +
                  p2p3 * (-2.0 * pow2(pap1) + m1s * papb - 2.0 * p1p3 * papb +
                          2.0 * pap3 * pbp1 +
                          2.0 * pap1 * (pap2 + pap3 - pbp1 + pbp2)) -
                  (m2s * pap1 - 2.0 * p1p2 * pap1 + m1s * pap2 +
                   2.0 * p1p2 * pap3) *
                      pbp3 +
                  p1p3 * (-(m2s * papb) - 2.0 * pap1 * (pap2 + pbp2) +
                          2.0 * pap2 * (pap2 + pap3 + pbp2 + pbp3))) *
                     RRRR)) +
       cA * (ivb1 * ivu2s2 *
                 (LLLL * (-2.0 * p1p3 * p2p3 * pap1 - 2.0 * p2p3 * pow2(pap1) -
                          2.0 * pow2(p1p3) * pap2 - 2.0 * p1p3 * pap1 * pap2 +
                          2.0 * pow2(p1p2) * pap3 - m2s * p1p3 * papb +
                          2.0 * p2p3 * pap1 * papb - 2.0 * p2p3 * pap2 * pbp1 +
                          m2s * pap3 * pbp1 + 2.0 * p1p3 * pap2 * pbp2 -
                          2.0 * pap1 * pap3 * pbp2 +
                          (m2s * pap1 + 2.0 * p2p3 * pap1 +
                           2.0 * (p1p3 + pap1) * pap2) *
                              pbp3 +
                          m1s * (2.0 * p2p3 * pap1 + 2.0 * p1p3 * pap2 -
                                 2.0 * p1p2 * pap3 - p2p3 * papb + pap3 * pbp2 -
                                 pap2 * pbp3) -
                          2.0 * p1p2 *
                              (p1p3 * (pap2 - pap3) + p2p3 * (pap1 - papb) +
                               pap3 * (-pap1 + pbp2 + pbp3))) +
                  m1 * m2 * (2.0 * m1s * pap3 - 2.0 * p1p2 * pap3 -
                             2.0 * p1p3 * pap3 - 2.0 * pap1 * pap3 +
                             p1p3 * papb - p2p3 * papb - pap3 * pbp1 +
                             pap3 * pbp2 + (-pap1 + pap2 + 2.0 * pap3) * pbp3) *
                      (LRRL + RLLR) +
                  (-2.0 * p1p3 * p2p3 * pap1 - 2.0 * p2p3 * pow2(pap1) -
                   2.0 * pow2(p1p3) * pap2 - 2.0 * p1p3 * pap1 * pap2 +
                   2.0 * pow2(p1p2) * pap3 - m2s * p1p3 * papb +
                   2.0 * p2p3 * pap1 * papb - 2.0 * p2p3 * pap2 * pbp1 +
                   m2s * pap3 * pbp1 + 2.0 * p1p3 * pap2 * pbp2 -
                   2.0 * pap1 * pap3 * pbp2 +
                   (m2s * pap1 + 2.0 * p2p3 * pap1 +
                    2.0 * (p1p3 + pap1) * pap2) *
                       pbp3 +
                   m1s * (2.0 * p2p3 * pap1 + 2.0 * p1p3 * pap2 -
                          2.0 * p1p2 * pap3 - p2p3 * papb + pap3 * pbp2 -
                          pap2 * pbp3) -
                   2.0 * p1p2 * (p1p3 * (pap2 - pap3) + p2p3 * (pap1 - papb) +
                                 pap3 * (-pap1 + pbp2 + pbp3))) *
                      RRRR) +
             iv13s2 *
                 (ivu2s2 * (m1s + m2s - 2.0 * (p1p2 + pap3)) *
                      (LLLL * (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) +
                       m1 * m2 * pap3 * (LRRL + RLLR) +
                       (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) * RRRR) +
                  ivs *
                      (LLLL * (-2.0 * p1p2 * pap1 * pap3 +
                               2.0 * p1p2 * pap2 * pap3 +
                               2.0 * p1p2 * pow2(pap3) - m2s * pap3 * pbp1 -
                               m1s * pap3 * pbp2 + 2.0 * p1p2 * pap3 * pbp2 +
                               p2p3 *
                                   (2.0 * pow2(pap1) - m1s * papb +
                                    2.0 * p1p3 * papb - 2.0 * pap3 * pbp1 -
                                    2.0 * pap1 * (pap2 + pap3 - pbp1 + pbp2)) +
                               m2s * pap1 * pbp3 - 2.0 * p1p2 * pap1 * pbp3 +
                               m1s * pap2 * pbp3 + 2.0 * p1p2 * pap3 * pbp3 +
                               p1p3 *
                                   (m2s * papb + 2.0 * pap1 * (pap2 + pbp2) -
                                    2.0 * pap2 * (pap2 + pap3 + pbp2 + pbp3))) +
                       m1 * m2 * ((-p1p3 + p2p3) * papb +
                                  pap3 * (2.0 * pap1 - 2.0 * pap2 - 2.0 * pap3 +
                                          pbp1 - pbp2) +
                                  (pap1 - pap2 - 2.0 * pap3) * pbp3) *
                           (LRRL + RLLR) +
                       (-2.0 * p1p2 * pap1 * pap3 + 2.0 * p1p2 * pap2 * pap3 +
                        2.0 * p1p2 * pow2(pap3) - m2s * pap3 * pbp1 -
                        m1s * pap3 * pbp2 + 2.0 * p1p2 * pap3 * pbp2 +
                        p2p3 * (2.0 * pow2(pap1) - m1s * papb +
                                2.0 * p1p3 * papb - 2.0 * pap3 * pbp1 -
                                2.0 * pap1 * (pap2 + pap3 - pbp1 + pbp2)) +
                        m2s * pap1 * pbp3 - 2.0 * p1p2 * pap1 * pbp3 +
                        m1s * pap2 * pbp3 + 2.0 * p1p2 * pap3 * pbp3 +
                        p1p3 * (m2s * papb + 2.0 * pap1 * (pap2 + pbp2) -
                                2.0 * pap2 * (pap2 + pap3 + pbp2 + pbp3))) *
                           RRRR))));

  // - (MtuSS + MtuSP + MtuS1 + MtuSB )
  MtuQ +=
      -4.0 * cF * iv23s1 * ivs * NC *
      (-(cA * ivu2s2 *
         (LRRL * m1 * m2 * pap3 +
          LLLL * (pap1 * pap2 + 2.0 * ivb1 * p1p3 * pap1 * pap2 - p1p2 * pap3 -
                  2.0 * ivb1 * p1p2 * pap1 * pap3 -
                  2.0 * ivu3 * p1p3 * pap2 * pap3 +
                  2.0 * ivu3 * p1p2 * pow2(pap3) +
                  2.0 * ivu3 * p1p3 * pap2 * papb -
                  2.0 * ivu3 * p1p2 * pap3 * papb + pap2 * pbp1 +
                  2.0 * ivb1 * p1p3 * pap2 * pbp1 -
                  2.0 * ivb1 * p1p2 * pap3 * pbp1 +
                  p2p3 * (2.0 * ivb1 * pow2(pap1) +
                          pap1 * (1.0 - 2.0 * ivu3 * pap3 - 2.0 * ivb1 * papb +
                                  2.0 * ivb1 * pbp1) -
                          papb * (ivb1 * m1s - 2.0 * ivu3 * p1p3 +
                                  2.0 * (ivb1 + ivu3) * pbp1)) +
                  ivb1 * m1s * pap3 * pbp2 - 2.0 * ivu3 * p1p3 * pap3 * pbp2 +
                  2.0 * ivb1 * pap1 * pap3 * pbp2 +
                  2.0 * ivu3 * pap1 * pap3 * pbp2 +
                  2.0 * ivb1 * pap3 * pbp1 * pbp2 +
                  2.0 * ivu3 * pap3 * pbp1 * pbp2 -
                  ivb1 * pap2 * (m1s + 2.0 * (pap1 + pbp1)) * pbp3) +
          ivb1 * m1 * m2 *
              (-(p1p3 * papb) + pap3 * pbp1 + pap1 * (2.0 * pap3 + pbp3)) *
              (LRRL + RLLR) +
          m1 * m2 * pap3 * (RLLR - 2.0 * ivu3 * (pap3 - papb) * (LRRL + RLLR)) +
          (pap1 * pap2 - p1p2 * pap3 - 2.0 * ivu3 * p1p3 * pap2 * pap3 +
           2.0 * ivu3 * p1p2 * pow2(pap3) + 2.0 * ivu3 * p1p3 * pap2 * papb -
           2.0 * ivu3 * p1p2 * pap3 * papb +
           p2p3 * (pap1 - 2.0 * ivu3 * pap1 * pap3 +
                   2.0 * ivu3 * papb * (p1p3 - pbp1)) +
           pap2 * pbp1 + 2.0 * ivu3 * pap3 * (-p1p3 + pap1 + pbp1) * pbp2) *
              RRRR +
          ivb1 * (2.0 * p1p3 * pap2 * (pap1 + pbp1) -
                  2.0 * p1p2 * pap3 * (pap1 + pbp1) +
                  p2p3 * (-(m1s * papb) + 2.0 * (pap1 - papb) * (pap1 + pbp1)) +
                  pap3 * (m1s + 2.0 * (pap1 + pbp1)) * pbp2 -
                  pap2 * (m1s + 2.0 * (pap1 + pbp1)) * pbp3) *
              RRRR)) +
       cF * (4.0 * iv13s2 * ivs * papb *
                 (LLLL * (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) +
                  m1 * m2 * pbp3 * (LRRL + RLLR) +
                  (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) * RRRR) +
             ivu2s2 *
                 (2.0 * LRRL * m1 * m2 * pap3 +
                  LLLL *
                      (2.0 * iv13s2 * p1p3 * pow2(pap2) - 2.0 * p1p2 * pap3 -
                       2.0 * iv13s2 * p1p2 * pap2 * pap3 -
                       2.0 * iv13s2 * p1p3 * pap2 * pap3 -
                       4.0 * ivu3 * p1p3 * pap2 * pap3 +
                       2.0 * iv13s2 * p1p2 * pow2(pap3) +
                       4.0 * ivu3 * p1p2 * pow2(pap3) -
                       iv13s2 * m2s * p1p3 * papb +
                       4.0 * ivu3 * p1p3 * pap2 * papb -
                       4.0 * ivu3 * p1p2 * pap3 * papb + 2.0 * pap2 * pbp1 -
                       2.0 * iv13s2 * p1p3 * pap2 * pbp1 -
                       iv13s2 * m2s * pap3 * pbp1 +
                       2.0 * iv13s2 * p1p2 * pap3 * pbp1 +
                       2.0 * iv13s2 * p1p3 * pap2 * pbp2 -
                       iv13s2 * m1s * pap3 * pbp2 -
                       2.0 * iv13s2 * p1p3 * pap3 * pbp2 -
                       4.0 * ivu3 * p1p3 * pap3 * pbp2 +
                       4.0 * ivu3 * pap3 * pbp1 * pbp2 +
                       iv13s2 * m1s * pap2 * pbp3 -
                       2.0 * iv13s2 * p1p2 * pap2 * pbp3 +
                       2.0 * iv13s2 * p1p2 * pap3 * pbp3 +
                       pap1 * ((2.0 - 2.0 * iv13s2 * p1p3) * pap2 +
                               2.0 * iv13s2 * p1p2 * pap3 +
                               4.0 * ivu3 * pap3 * pbp2 + iv13s2 * m2s * pbp3) +
                       p2p3 * (-2.0 * iv13s2 * pow2(pap1) +
                               4.0 * ivu3 * papb * (p1p3 - pbp1) +
                               iv13s2 * (m1s * papb + 2.0 * p1p3 * papb +
                                         2.0 * pap2 * pbp1) -
                               2.0 * pap1 *
                                   (-1.0 + 2.0 * ivu3 * pap3 +
                                    iv13s2 * (-pap2 + pap3 + pbp1 + pbp3)))) +
                  2.0 * m1 * m2 * pap3 *
                      (RLLR - 2.0 * ivu3 * (pap3 - papb) * (LRRL + RLLR)) +
                  2.0 * (pap1 * pap2 - p1p2 * pap3 -
                         2.0 * ivu3 * p1p3 * pap2 * pap3 +
                         2.0 * ivu3 * p1p2 * pow2(pap3) +
                         2.0 * ivu3 * p1p3 * pap2 * papb -
                         2.0 * ivu3 * p1p2 * pap3 * papb +
                         p2p3 * (pap1 - 2.0 * ivu3 * pap1 * pap3 +
                                 2.0 * ivu3 * papb * (p1p3 - pbp1)) +
                         pap2 * pbp1 +
                         2.0 * ivu3 * pap3 * (-p1p3 + pap1 + pbp1) * pbp2) *
                      RRRR -
                  iv13s2 *
                      (m1 * m2 * ((-p1p3 + p2p3) * papb +
                                  pap3 * (2.0 * pap1 - 2.0 * pap2 + 2.0 * pap3 +
                                          pbp1 - pbp2) +
                                  (pap1 - pap2 + 2.0 * pap3) * pbp3) *
                           (LRRL + RLLR) +
                       (pap3 * (m2s * pbp1 -
                                2.0 * p1p2 * (pap1 - pap2 + pap3 + pbp1) +
                                m1s * pbp2) +
                        p1p3 * (m2s * papb +
                                2.0 * pap2 * (pap1 - pap2 + pap3 + pbp1) +
                                2.0 * (-pap2 + pap3) * pbp2) -
                        (m2s * pap1 + m1s * pap2 - 2.0 * p1p2 * pap2 +
                         2.0 * p1p2 * pap3) *
                            pbp3 +
                        p2p3 * (2.0 * pow2(pap1) - (m1s + 2.0 * p1p3) * papb -
                                2.0 * pap2 * pbp1 +
                                2.0 * pap1 * (-pap2 + pap3 + pbp1 + pbp3))) *
                           RRRR))));

  return real(MtuQ);
}

// real antiquark emission

double FI::MtuQB_GLGA() {

  complex<double> MtuQB;

  // - (MtuAA + MtuA2 + MtuAG + MtuAS);
  MtuQB =
      2.0 * cF * ivt2s1 * NC *
      (cA * iv23s2 *
           (2.0 * ivs *
                (LLLL * (pap1 * pbp2 + pbp1 * (p2p3 + pbp2) - p1p2 * pbp3) +
                 m1 * m2 * pbp3 * (LRRL + RLLR) +
                 (pap1 * pbp2 + pbp1 * (p2p3 + pbp2) - p1p2 * pbp3) * RRRR +
                 2.0 * ivt3 *
                     (LLLL * (papb * (-(p2p3 * pap1) + p1p3 * (p2p3 + pbp2)) -
                              (p1p2 * papb + p2p3 * pbp1 -
                               pap2 * (pap1 + pbp1) + p1p3 * (pap2 + pbp2)) *
                                  pbp3 +
                              p1p2 * pow2(pbp3)) +
                      m1 * m2 * (papb - pbp3) * pbp3 * (LRRL + RLLR) +
                      (papb * (-(p2p3 * pap1) + p1p3 * (p2p3 + pbp2)) -
                       (p1p2 * papb + p2p3 * pbp1 - pap2 * (pap1 + pbp1) +
                        p1p3 * (pap2 + pbp2)) *
                           pbp3 +
                       p1p2 * pow2(pbp3)) *
                          RRRR)) +
            iva1 *
                (m1 * m2 * pbp1 * (LRRL + RLLR) +
                 (m1s - 2.0 * p1p3 + 2.0 * pap1) * pbp2 * (LLLL + RRRR) +
                 2.0 * ivt3 * (-(LLLL * m1s * p2p3 * papb) +
                               2.0 * LLLL * p1p3 * p2p3 * papb -
                               2.0 * LLLL * p2p3 * pap1 * papb -
                               2.0 * LLLL * p1p3 * p2p3 * pbp1 +
                               2.0 * LLLL * p2p3 * pap1 * pbp1 -
                               2.0 * LLLL * pow2(p1p3) * pbp2 +
                               2.0 * LLLL * p1p3 * pap1 * pbp2 +
                               2.0 * LLLL * p1p2 * p1p3 * pbp3 -
                               2.0 * LLLL * p1p2 * pap1 * pbp3 +
                               LLLL * m1s * pap2 * pbp3 -
                               2.0 * LLLL * p1p3 * pap2 * pbp3 +
                               2.0 * LLLL * pap1 * pap2 * pbp3 +
                               LRRL * m1 * m2 *
                                   (p1p3 * (papb - 2.0 * pbp3) + pap1 * pbp3) +
                               m1 * m2 * p1p3 * papb * RLLR -
                               2.0 * m1 * m2 * p1p3 * pbp3 * RLLR +
                               m1 * m2 * pap1 * pbp3 * RLLR +
                               (m1s * (-(p2p3 * papb) + pap2 * pbp3) -
                                2.0 * (p1p3 - pap1) *
                                    (-(p2p3 * papb) + p2p3 * pbp1 +
                                     p1p3 * pbp2 - p1p2 * pbp3 + pap2 * pbp3)) *
                                   RRRR))) +
       cF * (iv23s2 *
                 (ivu1s2 * (LRRL * m1 * m2 * (pbp1 - pbp2) +
                            LLLL * (m1s * pbp2 - pbp1 * (m2s + 2.0 * pbp2)) +
                            m1 * m2 * (pbp1 - pbp2) * RLLR +
                            (m1s * pbp2 - pbp1 * (m2s + 2.0 * pbp2)) * RRRR) -
                  4.0 * ivs *
                      (LLLL *
                           (p2p3 * pbp1 + (pap1 + pbp1) * pbp2 - p1p2 * pbp3) +
                       m1 * m2 * pbp3 * (LRRL + RLLR) +
                       (p2p3 * pbp1 + (pap1 + pbp1) * pbp2 - p1p2 * pbp3) *
                           RRRR)) +
             2.0 * ivt3 *
                 (4.0 * iv23s2 * ivs *
                      (LLLL * (-(p1p3 * papb * (p2p3 + pbp2)) +
                               p1p3 * (pap2 + pbp2) * pbp3 +
                               p2p3 * (pap1 * papb + pbp1 * pbp3) -
                               pbp3 * (pap2 * (pap1 + pbp1) +
                                       p1p2 * (-papb + pbp3))) -
                       m1 * m2 * (papb - pbp3) * pbp3 * (LRRL + RLLR) +
                       (-(p1p3 * papb * (p2p3 + pbp2)) +
                        p1p3 * (pap2 + pbp2) * pbp3 +
                        p2p3 * (pap1 * papb + pbp1 * pbp3) -
                        pbp3 * (pap2 * (pap1 + pbp1) + p1p2 * (-papb + pbp3))) *
                           RRRR) +
                  ivu1s2 *
                      (LRRL * m1 * m2 *
                           (papb * (2.0 + iv23s2 * (p1p3 - p2p3 - 2.0 * pbp3)) +
                            iv23s2 * pbp3 * (-2.0 * p1p3 + 2.0 * p2p3 + pap1 -
                                             pap2 + 2.0 * pbp3)) +
                       LLLL *
                           (2.0 * (pap2 * pbp1 + pap1 * pbp2) +
                            iv23s2 *
                                (-((m2s * p1p3 + m1s * p2p3) * papb) +
                                 2.0 * p2p3 * (-p1p3 + p2p3 + pap1 - pap2) *
                                     pbp1 -
                                 2.0 * ((p1p3 - p2p3) * (p1p3 - pap1) +
                                        p1p3 * papb) *
                                     pbp2 +
                                 (m2s * pap1 + m1s * pap2 + 2.0 * p2p3 * pbp1 -
                                  2.0 * pap2 * pbp1 + 2.0 * p1p3 * pbp2) *
                                     pbp3) +
                            2.0 * p1p2 *
                                (iv23s2 * (p1p3 - p2p3 - pap1 - pbp3) * pbp3 +
                                 papb * (-1.0 + iv23s2 * (p2p3 + pbp3)))) +
                       m1 * m2 *
                           (papb * (2.0 + iv23s2 * (p1p3 - p2p3 - 2.0 * pbp3)) +
                            iv23s2 * pbp3 * (-2.0 * p1p3 + 2.0 * p2p3 + pap1 -
                                             pap2 + 2.0 * pbp3)) *
                           RLLR +
                       (2.0 * (pap2 * pbp1 + pap1 * pbp2) +
                        iv23s2 *
                            (-((m2s * p1p3 + m1s * p2p3) * papb) +
                             2.0 * p2p3 * (-p1p3 + p2p3 + pap1 - pap2) * pbp1 -
                             2.0 *
                                 ((p1p3 - p2p3) * (p1p3 - pap1) + p1p3 * papb) *
                                 pbp2 +
                             (m2s * pap1 + m1s * pap2 + 2.0 * p2p3 * pbp1 -
                              2.0 * pap2 * pbp1 + 2.0 * p1p3 * pbp2) *
                                 pbp3) +
                        2.0 * p1p2 *
                            (iv23s2 * (p1p3 - p2p3 - pap1 - pbp3) * pbp3 +
                             papb * (-1.0 + iv23s2 * (p2p3 + pbp3)))) *
                           RRRR))));

  // - (Mtu1A + Mtu12 +Mtu1G + Mtu1S)
  MtuQB +=
      cF * iv13s1 * ivt2s1 * NC *
      (2.0 * cF *
           (-(ivu1s2 *
              (-(LLLL * m2s * pbp1) + LRRL * m1 * m2 * (pbp1 - pbp2) +
               LLLL * m1s * pbp2 + 2.0 * LLLL * pbp1 * pbp2 +
               m1 * m2 * pbp1 * RLLR - m1 * m2 * pbp2 * RLLR -
               m2s * pbp1 * RRRR + m1s * pbp2 * RRRR +
               2.0 * pbp1 * pbp2 * RRRR +
               2.0 * ivt3 *
                   (LRRL * m1 * m2 *
                        (-(p2p3 * papb) + p1p3 * (papb - 2.0 * pbp3) +
                         (2.0 * p2p3 + pap1 - pap2 + 2.0 * papb) * pbp3 -
                         2.0 * pow2(pbp3)) +
                    LLLL *
                        (m2s * p1p3 * papb - 2.0 * p1p2 * p1p3 * papb +
                         m1s * p2p3 * papb - 2.0 * p1p3 * p2p3 * pbp1 +
                         2.0 * pow2(p2p3) * pbp1 + 2.0 * p1p3 * pap2 * pbp1 -
                         2.0 * p2p3 * pap2 * pbp1 + 2.0 * p2p3 * papb * pbp1 -
                         2.0 * pow2(p1p3) * pbp2 + 2.0 * p1p3 * p2p3 * pbp2 +
                         2.0 * p1p3 * pap1 * pbp2 - 2.0 * p1p3 * pap2 * pbp2 +
                         (-(m2s * pap1) - m1s * pap2 +
                          2.0 * p1p2 * (p1p3 - p2p3 + pap2 - papb) -
                          2.0 * p2p3 * pbp1 - 2.0 * p1p3 * pbp2 +
                          2.0 * pap1 * pbp2) *
                             pbp3 +
                         2.0 * p1p2 * pow2(pbp3)) +
                    m1 * m2 * (-(p2p3 * papb) + p1p3 * (papb - 2.0 * pbp3) +
                               (2.0 * p2p3 + pap1 - pap2 + 2.0 * papb) * pbp3 -
                               2.0 * pow2(pbp3)) *
                        RLLR +
                    (m2s * p1p3 * papb - 2.0 * p1p2 * p1p3 * papb +
                     m1s * p2p3 * papb - 2.0 * p1p3 * p2p3 * pbp1 +
                     2.0 * pow2(p2p3) * pbp1 + 2.0 * p1p3 * pap2 * pbp1 -
                     2.0 * p2p3 * pap2 * pbp1 + 2.0 * p2p3 * papb * pbp1 -
                     2.0 * pow2(p1p3) * pbp2 + 2.0 * p1p3 * p2p3 * pbp2 +
                     2.0 * p1p3 * pap1 * pbp2 - 2.0 * p1p3 * pap2 * pbp2 +
                     (-(m2s * pap1) - m1s * pap2 +
                      2.0 * p1p2 * (p1p3 - p2p3 + pap2 - papb) -
                      2.0 * p2p3 * pbp1 - 2.0 * p1p3 * pbp2 +
                      2.0 * pap1 * pbp2) *
                         pbp3 +
                     2.0 * p1p2 * pow2(pbp3)) *
                        RRRR))) -
            2.0 * iv23s2 *
                (ivu1s2 * (m1s + m2s - 2.0 * (p1p2 + pbp3)) *
                     (LLLL * (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) +
                      m1 * m2 * pbp3 * (LRRL + RLLR) +
                      (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) * RRRR) +
                 ivs *
                     (LLLL *
                          (-(m2s * p1p3 * papb) + m1s * p2p3 * papb +
                           2.0 * p1p3 * p2p3 * papb - 2.0 * p2p3 * pap1 * pbp1 +
                           m2s * pap3 * pbp1 - 2.0 * p2p3 * pap3 * pbp1 -
                           2.0 * p2p3 * pow2(pbp1) - 2.0 * p1p3 * pap1 * pbp2 +
                           2.0 * p2p3 * pap1 * pbp2 + 2.0 * p1p3 * pap2 * pbp2 +
                           m1s * pap3 * pbp2 - 2.0 * p1p2 * pap3 * pbp2 -
                           2.0 * p1p3 * pbp1 * pbp2 + 2.0 * p2p3 * pbp1 * pbp2 +
                           2.0 * p1p3 * pow2(pbp2) -
                           (m2s * pap1 + (m1s + 2.0 * p1p3) * pap2 +
                            2.0 * (p2p3 * pbp1 - p1p2 * (pap1 + pap3 + pbp1) +
                                   (p1p2 + p1p3) * pbp2)) *
                               pbp3 +
                           2.0 * p1p2 * pow2(pbp3)) +
                      LRRL * m1 * m2 *
                          (p1p3 * papb - p2p3 * papb +
                           pap3 * (-pbp1 + pbp2 - 2.0 * pbp3) +
                           pbp3 * (-pap1 + pap2 - 2.0 * (pbp1 - pbp2 + pbp3))) +
                      m1 * m2 *
                          (p1p3 * papb - p2p3 * papb +
                           pap3 * (-pbp1 + pbp2 - 2.0 * pbp3) +
                           pbp3 * (-pap1 + pap2 - 2.0 * (pbp1 - pbp2 + pbp3))) *
                          RLLR +
                      (-(m2s * p1p3 * papb) + m1s * p2p3 * papb +
                       2.0 * p1p3 * p2p3 * papb - 2.0 * p2p3 * pap1 * pbp1 +
                       m2s * pap3 * pbp1 - 2.0 * p2p3 * pap3 * pbp1 -
                       2.0 * p2p3 * pow2(pbp1) - 2.0 * p1p3 * pap1 * pbp2 +
                       2.0 * p2p3 * pap1 * pbp2 + 2.0 * p1p3 * pap2 * pbp2 +
                       m1s * pap3 * pbp2 - 2.0 * p1p2 * pap3 * pbp2 -
                       2.0 * p1p3 * pbp1 * pbp2 + 2.0 * p2p3 * pbp1 * pbp2 +
                       2.0 * p1p3 * pow2(pbp2) -
                       (m2s * pap1 + (m1s + 2.0 * p1p3) * pap2 +
                        2.0 * (p2p3 * pbp1 - p1p2 * (pap1 + pap3 + pbp1) +
                               (p1p2 + p1p3) * pbp2)) *
                           pbp3 +
                       2.0 * p1p2 * pow2(pbp3)) *
                          RRRR))) +
       cA * (ivu1s2 *
                 (-(LLLL * m2s * pbp1) + LRRL * m1 * m2 * (pbp1 - pbp2) +
                  LLLL * m1s * pbp2 + 2.0 * LLLL * pbp1 * pbp2 +
                  m1 * m2 * pbp1 * RLLR - m1 * m2 * pbp2 * RLLR -
                  m2s * pbp1 * RRRR + m1s * pbp2 * RRRR +
                  2.0 * pbp1 * pbp2 * RRRR +
                  2.0 * ivt3 *
                      (LRRL * m1 * m2 *
                           (-(p2p3 * papb) + p1p3 * (papb - 2.0 * pbp3) +
                            (2.0 * p2p3 + pap1 - pap2 + 2.0 * papb) * pbp3 -
                            2.0 * pow2(pbp3)) +
                       LLLL *
                           (m2s * p1p3 * papb - 2.0 * p1p2 * p1p3 * papb +
                            m1s * p2p3 * papb - 2.0 * p1p3 * p2p3 * pbp1 +
                            2.0 * pow2(p2p3) * pbp1 + 2.0 * p1p3 * pap2 * pbp1 -
                            2.0 * p2p3 * pap2 * pbp1 +
                            2.0 * p2p3 * papb * pbp1 - 2.0 * pow2(p1p3) * pbp2 +
                            2.0 * p1p3 * p2p3 * pbp2 +
                            2.0 * p1p3 * pap1 * pbp2 -
                            2.0 * p1p3 * pap2 * pbp2 +
                            (-(m2s * pap1) - m1s * pap2 +
                             2.0 * p1p2 * (p1p3 - p2p3 + pap2 - papb) -
                             2.0 * p2p3 * pbp1 - 2.0 * p1p3 * pbp2 +
                             2.0 * pap1 * pbp2) *
                                pbp3 +
                            2.0 * p1p2 * pow2(pbp3)) +
                       m1 * m2 *
                           (-(p2p3 * papb) + p1p3 * (papb - 2.0 * pbp3) +
                            (2.0 * p2p3 + pap1 - pap2 + 2.0 * papb) * pbp3 -
                            2.0 * pow2(pbp3)) *
                           RLLR +
                       (m2s * p1p3 * papb - 2.0 * p1p2 * p1p3 * papb +
                        m1s * p2p3 * papb - 2.0 * p1p3 * p2p3 * pbp1 +
                        2.0 * pow2(p2p3) * pbp1 + 2.0 * p1p3 * pap2 * pbp1 -
                        2.0 * p2p3 * pap2 * pbp1 + 2.0 * p2p3 * papb * pbp1 -
                        2.0 * pow2(p1p3) * pbp2 + 2.0 * p1p3 * p2p3 * pbp2 +
                        2.0 * p1p3 * pap1 * pbp2 - 2.0 * p1p3 * pap2 * pbp2 +
                        (-(m2s * pap1) - m1s * pap2 +
                         2.0 * p1p2 * (p1p3 - p2p3 + pap2 - papb) -
                         2.0 * p2p3 * pbp1 - 2.0 * p1p3 * pbp2 +
                         2.0 * pap1 * pbp2) *
                            pbp3 +
                        2.0 * p1p2 * pow2(pbp3)) *
                           RRRR)) -
             2.0 * iv23s2 *
                 (-(ivu1s2 * (m1s + m2s - 2.0 * (p1p2 + pbp3)) *
                    (LLLL * (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) +
                     m1 * m2 * pbp3 * (LRRL + RLLR) +
                     (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3) * RRRR)) +
                  iva1 *
                      (LRRL * m1 * m2 *
                           (p1p3 * papb - p2p3 * papb - pap3 * pbp1 +
                            pap3 * pbp2 +
                            (-2.0 * m1s + 2.0 * p1p2 - 2.0 * p1p3 + pap1 -
                             pap2 + 2.0 * papb - 2.0 * pbp1) *
                                pbp3) +
                       LLLL *
                           (-(m2s * (p1p3 * papb - pap3 * pbp1 + pap1 * pbp3)) +
                            m1s * (p2p3 * (papb - 2.0 * pbp1) -
                                   2.0 * p1p3 * pbp2 + pap3 * pbp2 +
                                   2.0 * p1p2 * pbp3 - pap2 * pbp3) +
                            2.0 *
                                (-(p2p3 * pbp1 * (pap2 - papb + pbp1)) -
                                 pow2(p1p3) * pbp2 + p2p3 * pap1 * pbp2 -
                                 pow2(p1p2) * pbp3 +
                                 p1p2 * (p2p3 * pbp1 - pap3 * pbp2 +
                                         (pap2 - papb + pbp1) * pbp3) +
                                 p1p3 * (pap3 * pbp2 +
                                         (papb - pbp1) * (p2p3 + pbp2) -
                                         pap2 * pbp3 + p1p2 * (pbp2 + pbp3)))) +
                       m1 * m2 * (p1p3 * papb - p2p3 * papb - pap3 * pbp1 +
                                  pap3 * pbp2 +
                                  (-2.0 * m1s + 2.0 * p1p2 - 2.0 * p1p3 + pap1 -
                                   pap2 + 2.0 * papb - 2.0 * pbp1) *
                                      pbp3) *
                           RLLR +
                       (-(m2s * (p1p3 * papb - pap3 * pbp1 + pap1 * pbp3)) +
                        m1s * (p2p3 * (papb - 2.0 * pbp1) - 2.0 * p1p3 * pbp2 +
                               pap3 * pbp2 + 2.0 * p1p2 * pbp3 - pap2 * pbp3) +
                        2.0 * (-(p2p3 * pbp1 * (pap2 - papb + pbp1)) -
                               pow2(p1p3) * pbp2 + p2p3 * pap1 * pbp2 -
                               pow2(p1p2) * pbp3 +
                               p1p2 * (p2p3 * pbp1 - pap3 * pbp2 +
                                       (pap2 - papb + pbp1) * pbp3) +
                               p1p3 * (pap3 * pbp2 +
                                       (papb - pbp1) * (p2p3 + pbp2) -
                                       pap2 * pbp3 + p1p2 * (pbp2 + pbp3)))) *
                           RRRR))));

  // - (MtuS2+ MtuSS + MtuSP + MtuSA)
  MtuQB +=
      -2.0 * cF * ivs * NC *
      (-2.0 * cF * iv13s1 *
           (-2.0 * ivu1s2 *
                (LLLL * (p1p3 * pbp2 + pbp1 * (pap2 + pbp2) - p1p2 * pbp3) +
                 m1 * m2 * pbp3 * (LRRL + RLLR) +
                 (p1p3 * pbp2 + pbp1 * (pap2 + pbp2) - p1p2 * pbp3) * RRRR +
                 2.0 * ivt3 *
                     (LLLL * (papb * (p1p3 * (p2p3 - pap2) + p2p3 * pbp1) -
                              (p1p2 * papb + p2p3 * (pap1 + pbp1) +
                               p1p3 * pbp2 - pap1 * (pap2 + pbp2)) *
                                  pbp3 +
                              p1p2 * pow2(pbp3)) +
                      m1 * m2 * (papb - pbp3) * pbp3 * (LRRL + RLLR) +
                      (papb * (p1p3 * (p2p3 - pap2) + p2p3 * pbp1) -
                       (p1p2 * papb + p2p3 * (pap1 + pbp1) + p1p3 * pbp2 -
                        pap1 * (pap2 + pbp2)) *
                           pbp3 +
                       p1p2 * pow2(pbp3)) *
                          RRRR)) +
            iv23s2 *
                (-4.0 * ivs * papb *
                     (LLLL * (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) +
                      m1 * m2 * pap3 * (LRRL + RLLR) +
                      (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3) * RRRR) +
                 ivu1s2 *
                     (LRRL * m1 * m2 * (p1p3 * papb - p2p3 * papb +
                                        pap3 * (-pbp1 + pbp2 + 2.0 * pbp3) +
                                        pbp3 * (-pap1 + pap2 +
                                                2.0 * (-pbp1 + pbp2 + pbp3))) +
                      LLLL *
                          (-(m2s * (p1p3 * papb + pap3 * pbp1 - pap1 * pbp3)) +
                           m1s * (p2p3 * papb - pap3 * pbp2 + pap2 * pbp3) +
                           2.0 *
                               (p2p3 * pbp1 * (-pap1 + pap2 - pbp1 + pbp2) +
                                p2p3 * (pap1 + pbp1) * pbp3 +
                                p1p2 * (pap3 * (pbp1 - pbp3) -
                                        pbp3 * (pap2 - pbp1 + pbp2 + pbp3)) +
                                p1p3 * (-(p2p3 * papb) + pap2 * (-pbp1 + pbp2) +
                                        pbp2 * (pap3 - pbp1 + pbp2 + pbp3)))) +
                      m1 * m2 * (p1p3 * papb - p2p3 * papb +
                                 pap3 * (-pbp1 + pbp2 + 2.0 * pbp3) +
                                 pbp3 * (-pap1 + pap2 +
                                         2.0 * (-pbp1 + pbp2 + pbp3))) *
                          RLLR +
                      (-(m2s * (p1p3 * papb + pap3 * pbp1 - pap1 * pbp3)) +
                       m1s * (p2p3 * papb - pap3 * pbp2 + pap2 * pbp3) +
                       2.0 * (p2p3 * pbp1 * (-pap1 + pap2 - pbp1 + pbp2) +
                              p2p3 * (pap1 + pbp1) * pbp3 +
                              p1p2 * (pap3 * (pbp1 - pbp3) -
                                      pbp3 * (pap2 - pbp1 + pbp2 + pbp3)) +
                              p1p3 * (-(p2p3 * papb) + pap2 * (-pbp1 + pbp2) +
                                      pbp2 * (pap3 - pbp1 + pbp2 + pbp3)))) *
                          RRRR))) +
       cA * (-2.0 * iv13s1 * ivu1s2 *
                 (LLLL * (p1p3 * pbp2 + pbp1 * (pap2 + pbp2) - p1p2 * pbp3) +
                  m1 * m2 * pbp3 * (LRRL + RLLR) +
                  (p1p3 * pbp2 + pbp1 * (pap2 + pbp2) - p1p2 * pbp3) * RRRR +
                  2.0 * ivt3 *
                      (LLLL * (papb * (p1p3 * (p2p3 - pap2) + p2p3 * pbp1) -
                               (p1p2 * papb + p2p3 * (pap1 + pbp1) +
                                p1p3 * pbp2 - pap1 * (pap2 + pbp2)) *
                                   pbp3 +
                               p1p2 * pow2(pbp3)) +
                       m1 * m2 * (papb - pbp3) * pbp3 * (LRRL + RLLR) +
                       (papb * (p1p3 * (p2p3 - pap2) + p2p3 * pbp1) -
                        (p1p2 * papb + p2p3 * (pap1 + pbp1) + p1p3 * pbp2 -
                         pap1 * (pap2 + pbp2)) *
                            pbp3 +
                        p1p2 * pow2(pbp3)) *
                           RRRR)) +
             iv23s2 *
                 (2.0 * iva1 * ivt2s1 *
                      (LRRL * m1 * m2 * (p1p3 * papb - pap3 * pbp1 -
                                         (pap1 + 2.0 * pbp1) * pbp3) +
                       LLLL * (m1s * (p2p3 * papb + pap3 * pbp2 - pap2 * pbp3) +
                               2.0 * (pap1 + pbp1) *
                                   (p2p3 * papb - p2p3 * pbp1 - p1p3 * pbp2 +
                                    pap3 * pbp2 + p1p2 * pbp3 - pap2 * pbp3)) +
                       m1 * m2 * (p1p3 * papb - pap3 * pbp1 -
                                  (pap1 + 2.0 * pbp1) * pbp3) *
                           RLLR +
                       (m1s * (p2p3 * papb + pap3 * pbp2 - pap2 * pbp3) +
                        2.0 * (pap1 + pbp1) *
                            (p2p3 * papb - p2p3 * pbp1 - p1p3 * pbp2 +
                             pap3 * pbp2 + p1p2 * pbp3 - pap2 * pbp3)) *
                           RRRR) +
                  iv13s1 * ivu1s2 *
                      (LRRL * m1 * m2 * (p1p3 * papb - p2p3 * papb +
                                         pap3 * (-pbp1 + pbp2 + 2.0 * pbp3) +
                                         pbp3 * (-pap1 + pap2 +
                                                 2.0 * (-pbp1 + pbp2 + pbp3))) +
                       LLLL *
                           (-(m2s * (p1p3 * papb + pap3 * pbp1 - pap1 * pbp3)) +
                            m1s * (p2p3 * papb - pap3 * pbp2 + pap2 * pbp3) +
                            2.0 * (p2p3 * pbp1 * (-pap1 + pap2 - pbp1 + pbp2) +
                                   p2p3 * (pap1 + pbp1) * pbp3 +
                                   p1p2 * (pap3 * (pbp1 - pbp3) -
                                           pbp3 * (pap2 - pbp1 + pbp2 + pbp3)) +
                                   p1p3 *
                                       (-(p2p3 * papb) + pap2 * (-pbp1 + pbp2) +
                                        pbp2 * (pap3 - pbp1 + pbp2 + pbp3)))) +
                       m1 * m2 * (p1p3 * papb - p2p3 * papb +
                                  pap3 * (-pbp1 + pbp2 + 2.0 * pbp3) +
                                  pbp3 * (-pap1 + pap2 +
                                          2.0 * (-pbp1 + pbp2 + pbp3))) *
                           RLLR +
                       (-(m2s * (p1p3 * papb + pap3 * pbp1 - pap1 * pbp3)) +
                        m1s * (p2p3 * papb - pap3 * pbp2 + pap2 * pbp3) +
                        2.0 * (p2p3 * pbp1 * (-pap1 + pap2 - pbp1 + pbp2) +
                               p2p3 * (pap1 + pbp1) * pbp3 +
                               p1p2 * (pap3 * (pbp1 - pbp3) -
                                       pbp3 * (pap2 - pbp1 + pbp2 + pbp3)) +
                               p1p3 * (-(p2p3 * papb) + pap2 * (-pbp1 + pbp2) +
                                       pbp2 * (pap3 - pbp1 + pbp2 + pbp3)))) *
                           RRRR))));

  // - (MtuP2 + MtuPA + MtuPG + MtuGS)
  MtuQB +=
      2.0 * cA * cF * iva1 * NC *
      (-2.0 * iv13s1 * iv23s2 * ivs *
           (LLLL * m1s * p2p3 * papb + 2.0 * LLLL * p2p3 * pap1 * papb +
            2.0 * LLLL * p1p3 * pap2 * papb - 2.0 * LLLL * p1p2 * pap3 * papb -
            2.0 * LLLL * p2p3 * pap1 * pbp1 - 2.0 * LLLL * p1p3 * pap2 * pbp1 +
            2.0 * LLLL * p1p2 * pap3 * pbp1 + 2.0 * LLLL * p2p3 * papb * pbp1 -
            2.0 * LLLL * p2p3 * pow2(pbp1) - LLLL * m1s * pap3 * pbp2 +
            2.0 * LLLL * p1p3 * papb * pbp2 - 2.0 * LLLL * p1p3 * pbp1 * pbp2 +
            LLLL * m1s * pap2 * pbp3 - 2.0 * LLLL * p1p2 * papb * pbp3 +
            2.0 * LLLL * p1p2 * pbp1 * pbp3 +
            LRRL * m1 * m2 *
                (p1p3 * papb + 2.0 * pap3 * papb - pap3 * pbp1 - pap1 * pbp3 +
                 2.0 * papb * pbp3 - 2.0 * pbp1 * pbp3) +
            m1 * m2 * p1p3 * papb * RLLR + 2.0 * m1 * m2 * pap3 * papb * RLLR -
            m1 * m2 * pap3 * pbp1 * RLLR - m1 * m2 * pap1 * pbp3 * RLLR +
            2.0 * m1 * m2 * papb * pbp3 * RLLR -
            2.0 * m1 * m2 * pbp1 * pbp3 * RLLR +
            (m1s * (p2p3 * papb - pap3 * pbp2 + pap2 * pbp3) +
             2.0 * (papb - pbp1) *
                 (p2p3 * (pap1 + pbp1) + p1p3 * (pap2 + pbp2) -
                  p1p2 * (pap3 + pbp3))) *
                RRRR) +
       ivt2s1 *
           (-8.0 * iv23s2 * iva1 *
                (LLLL * (pap1 * (p2p3 * papb + pap3 * pbp2 - pap2 * pbp3) +
                         m1s * (p2p3 * papb - p2p3 * pbp1 - p1p3 * pbp2 +
                                pap3 * pbp2 + p1p2 * pbp3 - pap2 * pbp3)) -
                 m1 * m2 * (m1s - pap1) * pbp3 * (LRRL + RLLR) +
                 (pap1 * (p2p3 * papb + pap3 * pbp2 - pap2 * pbp3) +
                  m1s * (p2p3 * papb - p2p3 * pbp1 - p1p3 * pbp2 + pap3 * pbp2 +
                         p1p2 * pbp3 - pap2 * pbp3)) *
                     RRRR) +
            ivu1s2 *
                (LRRL * m1 * m2 *
                     ((2.0 + 2.0 * ivt3 * p1p3 + iv23s2 * (p1p3 - p2p3)) *
                          papb +
                      pbp1 - iv23s2 * pap3 * pbp1 + iv23s2 * pap3 * pbp2 +
                      (-2.0 - 4.0 * ivt3 * p1p3 + 2.0 * ivt3 * pap1 +
                       iv23s2 * (2.0 * m1s - 2.0 * p1p2 - 2.0 * p1p3 - pap1 +
                                 pap2 + 2.0 * pap3 - 2.0 * pbp1)) *
                          pbp3) +
                 LLLL *
                     (2.0 * ivt3 * m1s * p2p3 * papb - 2.0 * p2p3 * pbp1 -
                      4.0 * ivt3 * p1p3 * p2p3 * pbp1 + 2.0 * pap2 * pbp1 +
                      4.0 * ivt3 * p1p3 * pap2 * pbp1 + m1s * pbp2 -
                      2.0 * p1p3 * pbp2 - 4.0 * ivt3 * pow2(p1p3) * pbp2 +
                      2.0 * pap1 * pbp2 + 4.0 * ivt3 * p1p3 * pap1 * pbp2 -
                      2.0 * p1p2 * (papb + 2.0 * ivt3 * p1p3 * papb -
                                    iv23s2 * p2p3 * papb +
                                    iv23s2 * p2p3 * pbp1 +
                                    iv23s2 * p1p3 * pbp2) +
                      2.0 * iv23s2 * pow2(p1p2) * pbp3 -
                      2.0 * ivt3 * m1s * pap2 * pbp3 +
                      2.0 * p1p2 *
                          (1.0 + 2.0 * ivt3 * p1p3 -
                           iv23s2 * (m1s - p1p3 + pap2 + pap3 - pbp1)) *
                          pbp3 +
                      iv23s2 *
                          (m2s * (-(p1p3 * papb) + pap3 * pbp1 + pap1 * pbp3) +
                           m1s *
                               (-(p2p3 * papb) + 2.0 * p2p3 * pbp1 +
                                2.0 * p1p3 * pbp2 - pap3 * pbp2 + pap2 * pbp3) -
                           2.0 * (p1p3 * p2p3 * pbp1 -
                                  p2p3 * (pap3 + papb - pbp1) * pbp1 +
                                  pow2(p1p3) * pbp2 + p2p3 * pap1 * pbp2 -
                                  p1p3 * (pap2 + pap3 - pbp1) * pbp2 -
                                  pap3 * pbp1 * pbp2 + pap2 * pbp1 * pbp3))) +
                 2.0 * m1 * m2 * papb * RLLR +
                 iv23s2 * m1 * m2 * p1p3 * papb * RLLR +
                 2.0 * ivt3 * m1 * m2 * p1p3 * papb * RLLR -
                 iv23s2 * m1 * m2 * p2p3 * papb * RLLR + m1 * m2 * pbp1 * RLLR -
                 iv23s2 * m1 * m2 * pap3 * pbp1 * RLLR +
                 iv23s2 * m1 * m2 * pap3 * pbp2 * RLLR -
                 2.0 * m1 * m2 * pbp3 * RLLR +
                 2.0 * iv23s2 * m1 * m1s * m2 * pbp3 * RLLR -
                 2.0 * iv23s2 * m1 * m2 * p1p2 * pbp3 * RLLR -
                 2.0 * iv23s2 * m1 * m2 * p1p3 * pbp3 * RLLR -
                 4.0 * ivt3 * m1 * m2 * p1p3 * pbp3 * RLLR -
                 iv23s2 * m1 * m2 * pap1 * pbp3 * RLLR +
                 2.0 * ivt3 * m1 * m2 * pap1 * pbp3 * RLLR +
                 iv23s2 * m1 * m2 * pap2 * pbp3 * RLLR +
                 2.0 * iv23s2 * m1 * m2 * pap3 * pbp3 * RLLR -
                 2.0 * iv23s2 * m1 * m2 * pbp1 * pbp3 * RLLR +
                 (2.0 * ivt3 * m1s * p2p3 * papb - 2.0 * p2p3 * pbp1 -
                  4.0 * ivt3 * p1p3 * p2p3 * pbp1 + 2.0 * pap2 * pbp1 +
                  4.0 * ivt3 * p1p3 * pap2 * pbp1 + m1s * pbp2 -
                  2.0 * p1p3 * pbp2 - 4.0 * ivt3 * pow2(p1p3) * pbp2 +
                  2.0 * pap1 * pbp2 + 4.0 * ivt3 * p1p3 * pap1 * pbp2 -
                  2.0 * p1p2 *
                      (papb + 2.0 * ivt3 * p1p3 * papb - iv23s2 * p2p3 * papb +
                       iv23s2 * p2p3 * pbp1 + iv23s2 * p1p3 * pbp2) +
                  2.0 * iv23s2 * pow2(p1p2) * pbp3 -
                  2.0 * ivt3 * m1s * pap2 * pbp3 +
                  2.0 * p1p2 * (1.0 + 2.0 * ivt3 * p1p3 -
                                iv23s2 * (m1s - p1p3 + pap2 + pap3 - pbp1)) *
                      pbp3 +
                  iv23s2 *
                      (m2s * (-(p1p3 * papb) + pap3 * pbp1 + pap1 * pbp3) +
                       m1s * (-(p2p3 * papb) + 2.0 * p2p3 * pbp1 +
                              2.0 * p1p3 * pbp2 - pap3 * pbp2 + pap2 * pbp3) -
                       2.0 * (p1p3 * p2p3 * pbp1 -
                              p2p3 * (pap3 + papb - pbp1) * pbp1 +
                              pow2(p1p3) * pbp2 + p2p3 * pap1 * pbp2 -
                              p1p3 * (pap2 + pap3 - pbp1) * pbp2 -
                              pap3 * pbp1 * pbp2 + pap2 * pbp1 * pbp3))) *
                     RRRR)));

  return real(MtuQB);
}
