// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Real matrix elements.
//
// Notes on notation:
// - MGss:
//   - M = Squared Matrix element.
//   - G = Third outgoing particle (gluon G, quark Q or antiquark QB).
//   - ss = Channels 1 and 2 are s-channels.
//
// See mbn.cc for more details on notation.

#include "kinematics.h"
#include "utils.h"
#include <complex>

// for gauginos the color factor is factorized
#define cR 1.0
#define NA 1.0
#define I2R 1.0

#define DJAC std::sqrt(zz / zzp)

// Real gluon emission for: q + \bar{q} -> f + \bar{f} + g

double FI::MGss() {
  std::complex<double> mss;

  mss =
      -16.0 * ivs3v1 * ivs3v2 *

      (ivu3 * (((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 *
                   (pap3 - 2.0 * papb * (1.0 + 2.0 * ivt3 * papb)) +
               2.0 * (LRLR + RLRL) * (p2p3 * (pap1 + 2.0 * ivt3 * pap1 * papb) +

                                      (-pap1 + 2.0 * ivt3 * p1p3 * papb -
                                       4.0 * ivt3 * pap1 * papb + pbp1) *
                                          pbp2) +

               2.0 * (LLLL + RRRR) * (p1p3 * (pap2 + 2.0 * ivt3 * pap2 * papb) +
                                      pbp1 * (-pap2 + 2.0 * ivt3 * p2p3 * papb -
                                              4.0 * ivt3 * pap2 * papb + pbp2)))

       +
       ivt3 *
           (2.0 * (LLLL + RRRR) * (pap1 * pap2 + (p2p3 - pap2) * pbp1) +
            2.0 * (LRLR + RLRL) * (pap1 * pap2 + p1p3 * pbp2 - pap1 * pbp2) -
            ((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 * (2.0 * papb - pbp3)));

  return std::real(mss);
}

double FI::MGtt() {
  std::complex<double> mtt;

  mtt =
      -2.0 *
      (ivt2s2 *
           (4.0 * ivt2s1 * ivt3 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3 *
                pbp2 +
            ivt1s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivt2s1 * (m1s + pap1) + ivu3 * (-pap1 + pbp1)) * pbp2 +
                      ivt3 * (-4.0 * ivt2s1 * std::pow(pap1, 2) * pbp2 +
                              2.0 * ivu3 * p1p3 * papb * pbp2 +
                              pap1 * (pap2 +
                                      2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                                      (-1.0 + 4.0 * ivt2s1 * p1p3) * pbp2))) -
                 (LLLL + RRRR) *
                     (p1p2 +
                      2.0 * (-(ivt2s1 * (m1s + pap1) * pbp2) +
                             ivu3 * (p2p3 * pbp1 + (pap1 - pbp1) * pbp2) +
                             ivt3 *
                                 (p1p3 * (pap2 +
                                          2.0 * ivu3 * papb * (p2p3 - pbp2) -
                                          4.0 * ivt2s1 * pap1 * pbp2) +
                                  pap1 *
                                      (-pap2 -
                                       2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                                       pbp2 + 4.0 * ivt2s1 * pap1 * pbp2)))))) +
       ivt1s2 *
           (-2.0 * ivt1s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * pap1 *
                (-2.0 * ivu3 * p2p3 +
                 ivt2s2 *
                     (p1p2 - pap2 - 2.0 * ivu3 * p2p3 * papb +
                      2.0 * ivu3 * p2p3 * pbp1 + pbp2 +
                      2.0 * ivu3 * p1p3 * pbp2 - 2.0 * ivu3 * pap3 * pbp2 +
                      4.0 * ivu3 * papb * pbp2 - 4.0 * ivu3 * pbp1 * pbp2)) +
            ivt2s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivt2s2 * (m1s + pap1) + ivu3 * (-pap1 + pbp1)) * pbp2 -
                      ivt1s1 * pap1 *
                          (p1p2 - pap2 - 2.0 * ivu3 * p2p3 * papb +
                           2.0 * ivu3 * p2p3 * pbp1 + pbp2 -
                           4.0 * ivt2s2 * m1s * pbp2 -
                           4.0 * ivt2s2 * p1p3 * pbp2 +
                           2.0 * ivu3 * p1p3 * pbp2 +
                           8.0 * ivt2s2 * pap1 * pbp2 +
                           4.0 * ivt2s2 * pap3 * pbp2 -
                           2.0 * ivu3 * pap3 * pbp2 + 4.0 * ivu3 * papb * pbp2 -
                           4.0 * ivu3 * pbp1 * pbp2) +
                      ivt3 * (-4.0 * ivt2s2 * std::pow(pap1, 2) * pbp2 +
                              2.0 * ivu3 * p1p3 * papb * pbp2 +
                              pap1 * (pap2 +
                                      2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) +
                                      (-1.0 + 4.0 * ivt2s2 * p1p3) * pbp2))) -
                 (LLLL + RRRR) *
                     (p1p2 + 2.0 * ivt1s1 * p1p2 * pap1 +
                      2.0 *
                          (ivu3 * p2p3 * pbp1 - ivt2s2 * m1s * pbp2 -
                           ivt2s2 * pap1 * pbp2 + ivu3 * pap1 * pbp2 -
                           ivu3 * pbp1 * pbp2 +
                           ivt3 * (p1p3 * (pap2 +
                                           2.0 * ivu3 * papb * (p2p3 - pbp2) -
                                           4.0 * ivt2s2 * pap1 * pbp2) +
                                   pap1 * (-pap2 -
                                           2.0 * ivu3 * papb *
                                               (p2p3 - 2.0 * pbp2) +
                                           pbp2 + 4.0 * ivt2s2 * pap1 * pbp2)) +
                           ivt1s1 * pap1 *
                               (-pap2 +
                                (1.0 -
                                 4.0 * ivt2s2 *
                                     (m1s + p1p3 - 2.0 * pap1 - pap3)) *
                                    pbp2 -
                                2.0 * ivu3 *
                                    (p2p3 * (papb - pbp1) +
                                     (-p1p3 + pap3 - 2.0 * papb + 2.0 * pbp1) *
                                         pbp2)))))));

  return std::real(mtt);
}

double FI::MGuu() {
  std::complex<double> muu;

  muu =
      -2.0 *
      (-(ivu2s2 *
         (-4.0 * ivu2s1 * ivu3 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3 * pap2 +
          ivu1s1 *
              (2.0 * (LRLR + RLRL) *
                   (ivt3 * (-(pap1 * pap2) + pap2 * pbp1 -
                            2.0 * ivu3 * papb * (p1p3 * pap2 + p2p3 * pbp1 -
                                                 2.0 * pap2 * pbp1)) +
                    ivu3 * pbp1 * (pap2 - pbp2) +
                    ivu2s1 * pap2 *
                        (p1p2 - pap1 - 2.0 * ivu3 * p1p3 * papb + pbp1 +
                         2.0 * ivu3 * p2p3 * pbp1 - 2.0 * ivu3 * pap3 * pbp1 +
                         4.0 * ivu3 * papb * pbp1 + 2.0 * ivu3 * p1p3 * pbp2 -
                         4.0 * ivu3 * pbp1 * pbp2)) +
               (LLLL + RRRR) *
                   (p1p2 + 2.0 * ivu2s1 * p1p2 * pap2 +
                    2.0 *
                        (ivt3 * (p2p3 * (pap1 +
                                         2.0 * ivu3 * papb * (p1p3 - pbp1)) +
                                 pap2 * (-pap1 - 2.0 * ivu3 * p1p3 * papb +
                                         pbp1 + 4.0 * ivu3 * papb * pbp1)) +
                         ivu3 * (pap2 * pbp1 + (p1p3 - pbp1) * pbp2) +
                         ivu2s1 * pap2 *
                             (-pap1 + pbp1 -
                              2.0 * ivu3 * (p1p3 * (papb - pbp2) +
                                            pbp1 * (-p2p3 + pap3 - 2.0 * papb +
                                                    2.0 * pbp2)))))))) +
       ivu1s2 *
           (2.0 * ivu1s1 * ((LLLL + RRRR) + (LRLR + RLRL)) *
                (2.0 * ivt3 * p2p3 +
                 ivu2s2 *
                     (m2s +
                      pap2 * (1.0 + 4.0 * ivt3 * p2p3 - 4.0 * ivt3 * pap2))) *
                pbp1 +
            ivu2s1 *
                (2.0 * (LRLR + RLRL) *
                     (ivt3 * (pap1 * pap2 +
                              pap2 * (-1.0 + 4.0 * ivu1s1 * p2p3 -
                                      4.0 * ivu1s1 * pap2) *
                                  pbp1 +
                              2.0 * ivu3 * papb * (p1p3 * pap2 + p2p3 * pbp1 -
                                                   2.0 * pap2 * pbp1)) -
                      ivu2s2 * pap2 *
                          (p1p2 - pap1 - 2.0 * ivu3 * p1p3 * papb + pbp1 -
                           4.0 * ivu1s1 * m2s * pbp1 -
                           4.0 * ivu1s1 * p2p3 * pbp1 +
                           2.0 * ivu3 * p2p3 * pbp1 +
                           8.0 * ivu1s1 * pap2 * pbp1 +
                           4.0 * ivu1s1 * pap3 * pbp1 -
                           2.0 * ivu3 * pap3 * pbp1 + 4.0 * ivu3 * papb * pbp1 +
                           2.0 * ivu3 * p1p3 * pbp2 -
                           4.0 * ivu3 * pbp1 * pbp2) +
                      pbp1 * (ivu1s1 * (m2s + pap2) + ivu3 * (-pap2 + pbp2))) -
                 (LLLL + RRRR) *
                     (p1p2 + 2.0 * ivu2s2 * p1p2 * pap2 +
                      2.0 *
                          (-(ivu1s1 * m2s * pbp1) - ivu1s1 * pap2 * pbp1 +
                           ivu3 * pap2 * pbp1 +
                           ivt3 * (p2p3 * (pap1 +
                                           2.0 * ivu3 * papb * (p1p3 - pbp1) -
                                           4.0 * ivu1s1 * pap2 * pbp1) +
                                   pap2 * (-pap1 -
                                           2.0 * ivu3 * papb *
                                               (p1p3 - 2.0 * pbp1) +
                                           pbp1 + 4.0 * ivu1s1 * pap2 * pbp1)) +
                           ivu3 * p1p3 * pbp2 - ivu3 * pbp1 * pbp2 +
                           ivu2s2 * pap2 *
                               (-pap1 +
                                (1.0 -
                                 4.0 * ivu1s1 *
                                     (m2s + p2p3 - 2.0 * pap2 - pap3)) *
                                    pbp1 -
                                2.0 * ivu3 *
                                    (p1p3 * (papb - pbp2) +
                                     pbp1 * (-p2p3 + pap3 - 2.0 * papb +
                                             2.0 * pbp2))))))));

  return std::real(muu);
}

double FI::MGst() {
  std::complex<double> mst;

  mst =
      -4.0 * ivs3v1 *
      (-2.0 * ivt2s2 *
           (ivu3 * ((LLRL + RRLR) * m1 * m2 * papb +
                    (LRLR + RLRL) * (pap1 - pbp1) * pbp2) +
            ivt3 * (-((LRLR + RLRL) *
                      (pap1 * (pap2 + 2.0 * ivu3 * papb * (p2p3 - 2.0 * pbp2) -
                               pbp2) +
                       2.0 * p1p3 * (1.0 + ivu3 * papb) * pbp2)) +
                    (LLRL + RRLR) * m1 * m2 *
                        (papb + 2.0 * ivu3 * std::pow(papb, 2) - pbp3))) +
       ivt1s2 *
           (-2.0 *
                (ivt3 * ((LLRL + RRLR) * m1 * m2 * papb +
                         (LRLR + RLRL) * pap1 * (-pap2 + pbp2)) +
                 ivu3 * ((LLRL + RRLR) * m1 * m2 *
                             (-pap3 + papb + 2.0 * ivt3 * std::pow(papb, 2)) -
                         (LRLR + RLRL) *
                             (2.0 * p2p3 * (pap1 + ivt3 * pap1 * papb) +
                              (-pap1 + 2.0 * ivt3 * p1p3 * papb -
                               4.0 * ivt3 * pap1 * papb + pbp1) *
                                  pbp2))) +
            ivt2s2 *
                (2.0 * (LRLR + RLRL) *
                     (-(p1p2 * pap1) + m1s * pbp2 -
                      4.0 * ivt3 * std::pow(pap1, 2) * pbp2 +
                      pap1 * (pap2 +
                              2.0 * (2.0 * ivt3 * p1p3 * pbp2 +
                                     ivu3 * (p2p3 * (papb - pbp1) +
                                             (-p1p3 + pap3 - 2.0 * papb +
                                              2.0 * pbp1) *
                                                 pbp2)))) +
                 (LLRL + RRLR) * m1 * m2 *
                     (2.0 * ivt3 * p1p3 * papb - 2.0 * ivu3 * p1p3 * papb +
                      4.0 * ivu3 * pap3 * papb -
                      4.0 * ivu3 * std::pow(papb, 2) + pbp1 -
                      2.0 * ivu3 * pap3 * pbp1 + 4.0 * ivu3 * papb * pbp1 +
                      pap1 * (-1.0 - 4.0 * ivt3 * papb + 2.0 * ivt3 * pbp3)))));

  return std::real(mst);
}

double FI::MGsu() {
  std::complex<double> msu;

  msu =
      4.0 * ivs3v1 *
      (2.0 * ivu2s2 *
           (ivt3 * (-((LRLR + RLRL) * m1 * m2 * papb) +
                    (LLRL + RRLR) * pap2 * (pap1 - pbp1)) +
            ivu3 *
                ((LRLR + RLRL) * m1 * m2 *
                     (pap3 - papb * (1.0 + 2.0 * ivt3 * papb)) +
                 (LLRL + RRLR) * (2.0 * p1p3 * (pap2 + ivt3 * pap2 * papb) +
                                  pbp1 * (-pap2 + 2.0 * ivt3 * p2p3 * papb -
                                          4.0 * ivt3 * pap2 * papb + pbp2)))) +
       ivu1s2 *
           (2.0 *
                (-(ivu3 * ((LRLR + RLRL) * m1 * m2 * papb +
                           (LLRL + RRLR) * pbp1 * (pap2 - pbp2))) +
                 ivt3 * ((LLRL + RRLR) *
                             (pap1 * pap2 + (2.0 * p2p3 - pap2) * pbp1 +
                              2.0 * ivu3 * papb * (p1p3 * pap2 + p2p3 * pbp1 -
                                                   2.0 * pap2 * pbp1)) +
                         (LRLR + RLRL) * m1 * m2 *
                             (-papb - 2.0 * ivu3 * std::pow(papb, 2) + pbp3))) +
            ivu2s2 * (2.0 * (LLRL + RRLR) *
                          (-(p1p2 * pap2) + pap1 * pap2 +
                           2.0 * ivu3 * p1p3 * pap2 * papb + m2s * pbp1 +
                           4.0 * ivt3 * p2p3 * pap2 * pbp1 -
                           2.0 * ivu3 * p2p3 * pap2 * pbp1 -
                           4.0 * ivt3 * std::pow(pap2, 2) * pbp1 +
                           2.0 * ivu3 * pap2 * pap3 * pbp1 -
                           4.0 * ivu3 * pap2 * papb * pbp1 -
                           2.0 * ivu3 * p1p3 * pap2 * pbp2 +
                           4.0 * ivu3 * pap2 * pbp1 * pbp2) +
                      (LRLR + RLRL) * m1 * m2 *
                          (2.0 * ivt3 * p2p3 * papb - 2.0 * ivu3 * p2p3 * papb +
                           4.0 * ivu3 * pap3 * papb -
                           4.0 * ivu3 * std::pow(papb, 2) + pbp2 -
                           2.0 * ivu3 * pap3 * pbp2 + 4.0 * ivu3 * papb * pbp2 +
                           pap2 * (-1.0 - 4.0 * ivt3 * papb +
                                   2.0 * ivt3 * pbp3)))));

  return std::real(msu);
}

double FI::MGtu() {
  std::complex<double> mtu;

  mtu =
      -2.0 *
      (-(ivt2s1 *
         (ivu2s2 *
              (ivu1s2 * ((LRLR + RLRL) * m1 * m2 * (papb + pbp2) +
                         (LLLL + RRRR) * (-(p1p2 * papb) + m2s * pbp1 +
                                          pap2 * pbp1 + pap1 * pbp2)) -
               2.0 * ivu3 * ((LRLR + RLRL) * m1 * m2 * papb +
                             (LLLL + RRRR) * (-(p1p2 * papb) + pap2 * pbp1 +
                                              p1p3 * pbp2 - pbp1 * pbp2))) +
          2.0 * ivt3 *
              (ivu1s2 *
                   ((LRLR + RLRL) * m1 * m2 * pbp3 +
                    (LLLL + RRRR) * (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3)) +
               ivu2s2 * ((LRLR + RLRL) * m1 * m2 *
                             ((-1.0 + ivu1s2 * (p2p3 - 2.0 * pap2)) * papb -
                              2.0 * ivu3 * std::pow(papb, 2) +
                              ivu1s2 * pap2 * pbp3) +
                         (LLLL + RRRR) *
                             (p1p2 * papb + 2.0 * ivu1s2 * p1p2 * pap2 * papb +
                              2.0 * ivu3 * p1p2 * std::pow(papb, 2) +
                              2.0 * ivu1s2 * p2p3 * pap2 * pbp1 -
                              2.0 * ivu1s2 * std::pow(pap2, 2) * pbp1 +
                              2.0 * ivu3 * p2p3 * papb * pbp1 -
                              2.0 * ivu3 * pap2 * papb * pbp1 +
                              p1p3 * (-pap2 - ivu1s2 * m2s * papb -
                                      2.0 * ivu3 * p2p3 * papb +
                                      2.0 * ivu1s2 * pap2 * pbp2) -
                              2.0 * ivu1s2 * p1p2 * pap2 * pbp3 +
                              pap1 * (pap2 + 2.0 * ivu3 * papb * (p2p3 - pbp2) -
                                      pbp2 - 2.0 * ivu1s2 * pap2 * pbp2 +
                                      ivu1s2 * m2s * pbp3)))))) +
       ivt1s1 *
           (-2.0 * ivu2s2 * ivu3 *
                ((LRLR + RLRL) * m1 * m2 * pap3 +
                 (LLLL + RRRR) * (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3)) +
            ivu1s2 *
                (2.0 * (ivt3 * ((LRLR + RLRL) * m1 * m2 * papb *
                                    (1.0 + 2.0 * ivu3 * papb) +
                                (LLLL + RRRR) *
                                    (p2p3 * pap1 - pap1 * pap2 - p1p2 * papb +
                                     2.0 * ivu3 * p1p3 * p2p3 * papb -
                                     2.0 * ivu3 * p1p3 * pap2 * papb -
                                     2.0 * ivu3 * p1p2 * std::pow(papb, 2) +
                                     pap2 * pbp1 +
                                     2.0 * ivu3 * pap2 * papb * pbp1 -
                                     2.0 * ivu3 * p1p3 * papb * pbp2 +
                                     2.0 * ivu3 * pap1 * papb * pbp2)) +
                        ivu3 * ((LRLR + RLRL) * m1 * m2 * papb +
                                (LLLL + RRRR) * (-(p1p2 * papb) + p2p3 * pbp1 +
                                                 pap1 * pbp2 - pbp1 * pbp2))) +
                 ivu2s2 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (pap2 + 4.0 * ivu3 * std::pow(papb, 2) +
                           2.0 * ivu3 * pap3 * pbp2 +
                           papb * (1.0 + 2.0 * ivu3 * p2p3 - 4.0 * ivu3 * pap3 -
                                   4.0 * ivu3 * pbp2)) +
                      (LLLL + RRRR) *
                          (-(p1p2 * papb) - 4.0 * ivu3 * p1p2 * p2p3 * papb -
                           4.0 * ivu3 * p1p3 * pap2 * papb +
                           4.0 * ivu3 * p1p2 * pap3 * papb -
                           4.0 * ivu3 * p1p2 * std::pow(papb, 2) + pap2 * pbp1 +
                           4.0 * ivu3 * p2p3 * pap2 * pbp1 +
                           4.0 * ivu3 * pap2 * papb * pbp1 +
                           m2s * (pap1 + 2.0 * ivu3 * p1p3 * papb -
                                  2.0 * ivu3 * pap3 * pbp1) +
                           4.0 * ivu3 * p1p2 * papb * pbp2 -
                           4.0 * ivu3 * pap2 * pbp1 * pbp2 +
                           pap1 * (-2.0 * pap2 +
                                   (1.0 +
                                    4.0 * ivu3 * (p2p3 - pap3 + papb - pbp2)) *
                                       pbp2)))) +
            ivt2s1 *
                (ivu2s2 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (pap1 + 4.0 * ivu1s2 * pap1 * papb +
                           4.0 * ivu3 * std::pow(papb, 2) +
                           2.0 * ivu3 * pap3 * pbp1 +
                           papb * (1.0 -
                                   2.0 * ivu1s2 * (2.0 * p1p2 + p1p3 + p2p3 -
                                                   2.0 * pap2 - 2.0 * pap3) +
                                   2.0 * ivu3 * (p1p3 - 2.0 * (pap3 + pbp1)))) +
                      (LLLL + RRRR) *
                          ((-1.0 +
                            2.0 * ivu1s2 * (2.0 * p1p2 + p1p3 + p2p3 -
                                            2.0 * pap2 - 2.0 * pap3) -
                            4.0 * ivu3 * (p1p3 - pap3 + papb - pbp1)) *
                               (p1p2 * papb - pap2 * pbp1) +
                           4.0 * ivu1s2 * std::pow(pap1, 2) * pbp2 +
                           m1s * (pap2 + 2.0 * ivu3 * p2p3 * papb -
                                  2.0 * ivu3 * pap3 * pbp2) +
                           pap1 * (-4.0 * ivu3 * p2p3 * papb + pbp2 +
                                   4.0 * ivu3 * p1p3 * pbp2 +
                                   4.0 * ivu3 * papb * pbp2 -
                                   4.0 * ivu3 * pbp1 * pbp2 -
                                   2.0 * ivu1s2 *
                                       ((p1p3 + p2p3 - 2.0 * pap3) * pbp2 +
                                        2.0 * p1p2 * (papb + pbp2)) +
                                   pap2 * (-2.0 +
                                           4.0 * ivu1s2 * (pbp1 + pbp2))))) -
                 ivu1s2 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (papb + 2.0 * ivt3 * p1p3 * papb -
                           4.0 * ivt3 * pap1 * papb + pbp1 +
                           2.0 * ivt3 * pap1 * pbp3) +
                      (LLLL + RRRR) *
                          (pap2 * pbp1 + m1s * pbp2 + pap1 * pbp2 +
                           p1p2 * ((-1.0 + 4.0 * ivt3 * pap1) * papb -
                                   4.0 * ivt3 * pap1 * pbp3) +
                           ivt3 * (4.0 * pap1 * (p2p3 * pbp1 - pap2 * pbp1 +
                                                 p1p3 * pbp2 - pap1 * pbp2) +
                                   m1s * (-2.0 * p2p3 * papb +
                                          2.0 * pap2 * pbp3)))))));

  return std::real(mtu);
}

// Real quark emission for:  q + g -> f + \bar{f} + q

double FI::MQss() {
  std::complex<double> mss;

  mss =
      16.0 * ivs3v1 * ivs3v2 *
      (ivu3 * (-(((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 *
                 (-2.0 * pap3 + 4.0 * ivs * std::pow(pap3, 2) + papb)) +
               2.0 * (LRLR + RLRL) * (p1p3 * p2p3 +
                                      p2p3 * (pap1 - 4.0 * ivs * pap1 * pap3 -
                                              2.0 * ivs * pap3 * pbp1) +
                                      pap1 * (-1.0 + 2.0 * ivs * pap3) * pbp2) +
               2.0 * (LLLL + RRRR) *
                   (pap2 * (-1.0 + 2.0 * ivs * pap3) * pbp1 +
                    p1p3 * (p2p3 + pap2 - 4.0 * ivs * pap2 * pap3 -
                            2.0 * ivs * pap3 * pbp2))) +
       ivs * (2.0 * (LRLR + RLRL) * (pap1 * pap2 + p2p3 * (pap1 + pbp1)) +
              2.0 * (LLLL + RRRR) * (pap1 * pap2 + p1p3 * (pap2 + pbp2)) +
              ((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 * (2.0 * pap3 + pbp3)));

  return std::real(mss);
}

double FI::MQtt() {
  std::complex<double> mtt;

  mtt =
      2.0 *
      (ivt1s2 *
           (-4.0 * ivt1s1 * ivu3 * ((LLLL + RRRR) + (LRLR + RLRL)) * pap1 *
                pbp2 +
            ivs1s1 *
                (2.0 * (LRLR + RLRL) *
                     (pap1 * (ivs * (p2p3 + pap2) +
                              ivt1s1 * (-p1p2 + p2p3 + pap2)) +
                      ivu3 * (p2p3 * (-2.0 * ivs * pap3 * pbp1 +
                                      pap1 * (1.0 - 4.0 * ivs * pap3 -
                                              2.0 * ivt1s1 *
                                                  (2.0 * pap3 - papb + pbp1))) +
                              2.0 * (ivs + ivt1s1) * pap1 * pap3 * pbp2 +
                              p1p3 * (p2p3 + 4.0 * ivt1s1 * p2p3 * pap1 -
                                      2.0 * ivt1s1 * pap1 * pbp2))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivt1s1 * pap1)) +
                      2.0 * (ivs * p2p3 * pap1 + ivt1s1 * p2p3 * pap1 +
                             ivs * pap1 * pap2 + ivt1s1 * pap1 * pap2 +
                             ivs * pap2 * pbp1 +
                             ivu3 * (p2p3 * (-2.0 * ivs * pap3 * pbp1 +
                                             pap1 * (1.0 - 4.0 * ivs * pap3 -
                                                     2.0 * ivt1s1 *
                                                         (2.0 * pap3 - papb +
                                                          pbp1))) +
                                     2.0 * pap3 *
                                         (ivt1s1 * pap1 + ivs * (pap1 + pbp1)) *
                                         pbp2 +
                                     p1p3 * (p2p3 + 4.0 * ivt1s1 * p2p3 * pap1 -
                                             (1.0 + 2.0 * ivt1s1 * pap1) *
                                                 pbp2)))))) +
       ivs1s2 *
           (-2.0 * ivs1s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * p2p3 *
                (-2.0 * ivs * pbp1 +
                 ivt1s2 * (m1s + pap1 * (1.0 - 4.0 * ivs * (pap1 + pbp1)))) +
            ivt1s1 *
                (2.0 * (LRLR + RLRL) *
                     (pap1 * (ivs * (p2p3 + pap2) +
                              ivt1s2 * (-p1p2 + p2p3 + pap2)) -
                      ivs1s1 * p2p3 *
                          (m1s * (1.0 + 4.0 * ivt1s2 * pap1) -
                           pap1 * (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                   4.0 * ivt1s2 * (2.0 * pap1 - papb + pbp1))) +
                      ivu3 * (p2p3 * (-2.0 * ivs * pap3 * pbp1 +
                                      pap1 * (1.0 - 4.0 * ivs * pap3 -
                                              2.0 * ivt1s2 *
                                                  (2.0 * pap3 - papb + pbp1))) +
                              2.0 * (ivs + ivt1s2) * pap1 * pap3 * pbp2 +
                              p1p3 * (p2p3 + 4.0 * ivt1s2 * p2p3 * pap1 -
                                      2.0 * ivt1s2 * pap1 * pbp2))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivt1s2 * pap1)) +
                      2.0 * (ivs * p2p3 * pap1 + ivt1s2 * p2p3 * pap1 +
                             ivs * pap1 * pap2 + ivt1s2 * pap1 * pap2 +
                             ivs * pap2 * pbp1 -
                             ivs1s1 * p2p3 *
                                 (m1s * (1.0 + 4.0 * ivt1s2 * pap1) -
                                  pap1 * (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                          4.0 * ivt1s2 *
                                              (2.0 * pap1 - papb + pbp1))) +
                             ivu3 * (p2p3 * (-2.0 * ivs * pap3 * pbp1 +
                                             pap1 * (1.0 - 4.0 * ivs * pap3 -
                                                     2.0 * ivt1s2 *
                                                         (2.0 * pap3 - papb +
                                                          pbp1))) +
                                     2.0 * pap3 *
                                         (ivt1s2 * pap1 + ivs * (pap1 + pbp1)) *
                                         pbp2 +
                                     p1p3 * (p2p3 + 4.0 * ivt1s2 * p2p3 * pap1 -
                                             (1.0 + 2.0 * ivt1s2 * pap1) *
                                                 pbp2)))))));

  return std::real(mtt);
}

double FI::MQuu() {
  std::complex<double> muu;

  muu =
      2.0 *
      (ivu2s2 *
           (-4.0 * ivu2s1 * ivu3 * ((LLLL + RRRR) + (LRLR + RLRL)) * pap2 *
                pbp1 +
            ivs2s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivs * (p1p3 + pap1) + ivu2s1 * (-p1p2 + p1p3 + pap1)) *
                          pap2 +
                      ivu3 *
                          (2.0 * pap2 * (ivs * pap3 + ivu2s1 * (-p2p3 + pap3)) *
                               pbp1 +
                           p1p3 * (p2p3 + 4.0 * ivu2s1 * p2p3 * pap2 -
                                   2.0 * ivs * pap3 * pbp2 +
                                   pap2 * (1.0 - 4.0 * ivs * pap3 -
                                           2.0 * ivu2s1 *
                                               (2.0 * pap3 - papb + pbp2))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivu2s1 * pap2)) +
                      2.0 *
                          (ivu2s1 * (p1p3 + pap1) * pap2 +
                           ivs * (p1p3 * pap2 + pap1 * (pap2 + pbp2)) +
                           ivu3 *
                               (pbp1 * (-(p2p3 * (1.0 + 2.0 * ivu2s1 * pap2)) +
                                        2.0 * pap3 * (ivu2s1 * pap2 +
                                                      ivs * (pap2 + pbp2))) +
                                p1p3 *
                                    (p2p3 + 4.0 * ivu2s1 * p2p3 * pap2 -
                                     2.0 * ivs * pap3 * pbp2 +
                                     pap2 * (1.0 - 4.0 * ivs * pap3 -
                                             2.0 * ivu2s1 * (2.0 * pap3 - papb +
                                                             pbp2)))))))) +
       ivs2s2 *
           (-2.0 * ivs2s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3 *
                (-2.0 * ivs * pbp2 +
                 ivu2s2 * (m2s + pap2 * (1.0 - 4.0 * ivs * (pap2 + pbp2)))) +
            ivu2s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivs * (p1p3 + pap1) + ivu2s2 * (-p1p2 + p1p3 + pap1)) *
                          pap2 -
                      ivs2s1 * p1p3 *
                          (m2s * (1.0 + 4.0 * ivu2s2 * pap2) -
                           pap2 * (-1.0 + 4.0 * ivs * (pap2 + pbp2) +
                                   4.0 * ivu2s2 * (2.0 * pap2 - papb + pbp2))) +
                      ivu3 *
                          (2.0 * pap2 * (ivs * pap3 + ivu2s2 * (-p2p3 + pap3)) *
                               pbp1 +
                           p1p3 * (p2p3 + 4.0 * ivu2s2 * p2p3 * pap2 -
                                   2.0 * ivs * pap3 * pbp2 +
                                   pap2 * (1.0 - 4.0 * ivs * pap3 -
                                           2.0 * ivu2s2 *
                                               (2.0 * pap3 - papb + pbp2))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivu2s2 * pap2)) +
                      2.0 *
                          (ivs * p1p3 * pap2 + ivu2s2 * p1p3 * pap2 +
                           ivs * pap1 * pap2 + ivu2s2 * pap1 * pap2 +
                           ivs * pap1 * pbp2 -
                           ivs2s1 * p1p3 *
                               (m2s * (1.0 + 4.0 * ivu2s2 * pap2) -
                                pap2 * (-1.0 + 4.0 * ivs * (pap2 + pbp2) +
                                        4.0 * ivu2s2 *
                                            (2.0 * pap2 - papb + pbp2))) +
                           ivu3 *
                               (pbp1 * (-(p2p3 * (1.0 + 2.0 * ivu2s2 * pap2)) +
                                        2.0 * pap3 * (ivu2s2 * pap2 +
                                                      ivs * (pap2 + pbp2))) +
                                p1p3 *
                                    (p2p3 + 4.0 * ivu2s2 * p2p3 * pap2 -
                                     2.0 * ivs * pap3 * pbp2 +
                                     pap2 * (1.0 - 4.0 * ivs * pap3 -
                                             2.0 * ivu2s2 * (2.0 * pap3 - papb +
                                                             pbp2)))))))));

  return std::real(muu);
}

double FI::MQst() {
  std::complex<double> mst;

  mst = 4.0 * ivs3v1 *
        (2.0 * ivt1s2 *
             (ivs * ((LRLR + RLRL) * pap1 * (p2p3 + pap2) +
                     (LLRL + RRLR) * m1 * m2 * pap3) +
              ivu3 *
                  ((LLRL + RRLR) * m1 * m2 *
                       (pap3 - 2.0 * ivs * std::pow(pap3, 2) - papb) +
                   (LRLR + RLRL) * (p1p3 * p2p3 +
                                    p2p3 * (pap1 - 4.0 * ivs * pap1 * pap3 -
                                            2.0 * ivs * pap3 * pbp1) +
                                    2.0 * pap1 * (-1.0 + ivs * pap3) * pbp2))) +
         ivs1s2 *
             (2.0 * (ivu3 * ((LLRL + RRLR) * m1 * m2 * pap3 *
                                 (1.0 - 2.0 * ivs * pap3) +
                             (LRLR + RLRL) *
                                 (p1p3 * p2p3 +
                                  p2p3 * (pap1 - 4.0 * ivs * pap1 * pap3 -
                                          2.0 * ivs * pap3 * pbp1) +
                                  2.0 * ivs * pap1 * pap3 * pbp2)) +
                     ivs * ((LRLR + RLRL) *
                                (pap1 * pap2 + p2p3 * (pap1 + 2.0 * pbp1)) +
                            (LLRL + RRLR) * m1 * m2 * (pap3 + pbp3))) +
              ivt1s2 *
                  (2.0 * (LRLR + RLRL) *
                       (-(m1s * p2p3) +
                        pap1 * (-p1p2 + 4.0 * ivs * p2p3 * pap1 + pap2 +
                                4.0 * ivs * p2p3 * pbp1 +
                                2.0 * ivu3 *
                                    (2.0 * p1p3 * p2p3 - 2.0 * p2p3 * pap3 +
                                     p2p3 * papb - p2p3 * pbp1 - p1p3 * pbp2 +
                                     pap3 * pbp2))) +
                   (LLRL + RRLR) * m1 * m2 *
                       (p1p3 * (-1.0 + 4.0 * ivu3 * pap3 - 2.0 * ivu3 * papb) +
                        2.0 * pap3 * (ivs * pbp1 -
                                      ivu3 * (2.0 * pap3 - 2.0 * papb + pbp1)) +
                        pap1 * (-1.0 + 2.0 * ivs * (2.0 * pap3 + pbp3))))));

  return std::real(mst);
}

double FI::MQsu() {
  std::complex<double> msu;

  msu =
      -4.0 * ivs3v1 *
      (2.0 * ivu2s2 *
           (ivs * ((LLRL + RRLR) * (p1p3 + pap1) * pap2 +
                   (LRLR + RLRL) * m1 * m2 * pap3) +
            ivu3 * ((LRLR + RLRL) * m1 * m2 *
                        (pap3 - 2.0 * ivs * std::pow(pap3, 2) - papb) +
                    (LLRL + RRLR) *
                        (2.0 * pap2 * (-1.0 + ivs * pap3) * pbp1 +
                         p1p3 * (p2p3 + pap2 - 4.0 * ivs * pap2 * pap3 -
                                 2.0 * ivs * pap3 * pbp2)))) +
       ivs2s2 *
           (2.0 * (ivu3 * ((LRLR + RLRL) * m1 * m2 * pap3 *
                               (1.0 - 2.0 * ivs * pap3) +
                           (LLRL + RRLR) *
                               (2.0 * ivs * pap2 * pap3 * pbp1 +
                                p1p3 * (p2p3 + pap2 - 4.0 * ivs * pap2 * pap3 -
                                        2.0 * ivs * pap3 * pbp2))) +
                   ivs * ((LLRL + RRLR) *
                              (p1p3 * pap2 + pap1 * pap2 + 2.0 * p1p3 * pbp2) +
                          (LRLR + RLRL) * m1 * m2 * (pap3 + pbp3))) +
            ivu2s2 *
                ((LLRL + RRLR) *
                     (-2.0 * m2s * p1p3 +
                      2.0 * pap2 *
                          (-p1p2 + pap1 + 4.0 * ivs * p1p3 * pap2 +
                           2.0 * ivu3 * ((-p2p3 + pap3) * pbp1 +
                                         p1p3 * (2.0 * p2p3 - 2.0 * pap3 +
                                                 papb - pbp2)) +
                           4.0 * ivs * p1p3 * pbp2)) +
                 (LRLR + RLRL) * m1 * m2 *
                     (p2p3 * (-1.0 + 4.0 * ivu3 * pap3 - 2.0 * ivu3 * papb) +
                      2.0 * pap3 * (ivs * pbp2 -
                                    ivu3 * (2.0 * pap3 - 2.0 * papb + pbp2)) +
                      pap2 * (-1.0 + 2.0 * ivs * (2.0 * pap3 + pbp3))))));

  return std::real(msu);
}

double FI::MQtu() {
  std::complex<double> mtu;

  mtu =
      -2.0 *
      (ivt1s1 *
           (-2.0 * ivu2s2 * ivu3 *
                ((LRLR + RLRL) * m1 * m2 * papb +
                 (LLLL + RRRR) * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
            ivs2s2 *
                (-(ivu2s2 * (LLLL + RRRR) * m2s * pap1) +
                 ivu2s2 * (LLLL + RRRR) * p2p3 * pap1 -
                 ivu2s2 * (LRLR + RLRL) * m1 * m2 * pap2 +
                 2.0 * ivs * (LLLL + RRRR) * p1p3 * pap2 +
                 ivu2s2 * (LLLL + RRRR) * p1p3 * pap2 +
                 2.0 * ivs * (LLLL + RRRR) * pap1 * pap2 +
                 2.0 * ivu2s2 * (LLLL + RRRR) * pap1 * pap2 +
                 2.0 * ivs * (LRLR + RLRL) * m1 * m2 * pap3 +
                 ivu2s2 * (LRLR + RLRL) * m1 * m2 * pap3 -
                 2.0 * ivs * (LLLL + RRRR) * p1p2 * pap3 -
                 ivu2s2 * (LLLL + RRRR) * p1p2 * pap3 +
                 2.0 * ivs * (LLLL + RRRR) * pap1 * pbp2 +
                 2.0 * ivu3 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (-2.0 * (ivs + ivu2s2) * std::pow(pap3, 2) -
                           ivu2s2 * p2p3 * papb +
                           pap3 * (1.0 +
                                   ivu2s2 * (2.0 * p2p3 + 2.0 * papb - pbp2))) +
                      (LLLL + RRRR) *
                          (2.0 * ivu2s2 * std::pow(p2p3, 2) * pap1 +
                           p1p3 * (p2p3 + 2.0 * ivu2s2 * p2p3 * pap2 -
                                   2.0 * ivs * pap2 * pap3 -
                                   2.0 * ivu2s2 * pap2 * pap3 +
                                   ivu2s2 * m2s * papb - pbp2 -
                                   2.0 * ivu2s2 * pap2 * pbp2) +
                           p2p3 *
                               (-2.0 * pap3 * (ivu2s2 * p1p2 + ivs * pbp1) +
                                pap1 * (1.0 - 2.0 * ivs * pap3 -
                                        2.0 * ivu2s2 * (pap3 - papb + pbp2))) +
                           pap3 * (pbp1 * (-(ivu2s2 * (m2s - 2.0 * pap2)) +
                                           2.0 * ivs * (pap2 + pbp2)) +
                                   p1p2 * (-1.0 + 2.0 * ivs * pap3 +
                                           2.0 * ivu2s2 *
                                               (pap3 - papb + pbp2))))))) +
       ivs1s1 *
           (ivu2s2 *
                (-(ivt1s1 * (LRLR + RLRL) * m1 * m2 * pap1) +
                 2.0 * ivs * (LLLL + RRRR) * p2p3 * pap1 +
                 ivt1s1 * (LLLL + RRRR) * p2p3 * pap1 -
                 ivt1s1 * (LLLL + RRRR) * m1s * pap2 +
                 ivt1s1 * (LLLL + RRRR) * p1p3 * pap2 +
                 2.0 * ivs * (LLLL + RRRR) * pap1 * pap2 +
                 2.0 * ivt1s1 * (LLLL + RRRR) * pap1 * pap2 +
                 2.0 * ivs * (LRLR + RLRL) * m1 * m2 * pap3 +
                 ivt1s1 * (LRLR + RLRL) * m1 * m2 * pap3 -
                 2.0 * ivs * (LLLL + RRRR) * p1p2 * pap3 -
                 ivt1s1 * (LLLL + RRRR) * p1p2 * pap3 +
                 2.0 * ivs * (LLLL + RRRR) * pap2 * pbp1 +
                 2.0 * ivu3 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (-2.0 * (ivs + ivt1s1) * std::pow(pap3, 2) -
                           ivt1s1 * p1p3 * papb +
                           pap3 * (1.0 +
                                   ivt1s1 * (2.0 * p1p3 + 2.0 * papb - pbp1))) +
                      (LLLL + RRRR) *
                          (2.0 * ivt1s1 * std::pow(p1p3, 2) * pap2 -
                           2.0 * ivs * p2p3 * pap1 * pap3 -
                           2.0 * ivt1s1 * p2p3 * pap1 * pap3 +
                           ivt1s1 * m1s * p2p3 * papb - p2p3 * pbp1 -
                           2.0 * ivt1s1 * p2p3 * pap1 * pbp1 +
                           p1p2 * pap3 * (-1.0 + 2.0 * ivs * pap3 +
                                          2.0 * ivt1s1 * (pap3 - papb + pbp1)) -
                           ivt1s1 * m1s * pap3 * pbp2 +
                           2.0 * ivs * pap1 * pap3 * pbp2 +
                           2.0 * ivt1s1 * pap1 * pap3 * pbp2 +
                           2.0 * ivs * pap3 * pbp1 * pbp2 +
                           p1p3 *
                               (p2p3 + 2.0 * ivt1s1 * p2p3 * pap1 +
                                pap2 * (1.0 - 2.0 * ivs * pap3 -
                                        2.0 * ivt1s1 * (pap3 - papb + pbp1)) -
                                2.0 * pap3 * (ivt1s1 * p1p2 + ivs * pbp2))))) +
            ivs2s2 *
                (2.0 * ivs * ((LRLR + RLRL) * m1 * m2 * pbp3 +
                              (LLLL + RRRR) *
                                  (p2p3 * pbp1 + p1p3 * pbp2 - p1p2 * pbp3)) -
                 ivu2s2 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (p2p3 + pap3 - 4.0 * ivs * pap2 * pap3 -
                           2.0 * ivs * pap3 * pbp2 - 2.0 * ivs * pap2 * pbp3) +
                      (LLLL + RRRR) *
                          (p1p3 * pap2 - 4.0 * ivs * p1p3 * std::pow(pap2, 2) -
                           p1p2 * pap3 + 4.0 * ivs * p1p2 * pap2 * pap3 +
                           p2p3 * (pap1 - 4.0 * ivs * pap1 * pap2 -
                                   4.0 * ivs * pap2 * pbp1) -
                           4.0 * ivs * p1p3 * pap2 * pbp2 +
                           4.0 * ivs * p1p2 * pap2 * pbp3 +
                           m2s * (p1p3 + 2.0 * ivs * pap3 * pbp1 -
                                  2.0 * ivs * pap1 * pbp3))) +
                 ivt1s1 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (-p1p3 +
                           pap3 * (-1.0 + 4.0 * ivs * pap1 + 2.0 * ivs * pbp1 +
                                   2.0 * ivu2s2 *
                                       (-2.0 * p1p2 + 2.0 * pap1 + 2.0 * pap2 -
                                        2.0 * papb + pbp1 + pbp2)) +
                           2.0 * ivs * pap1 * pbp3) +
                      (LLLL + RRRR) *
                          (-(p1p3 * pap2) - 4.0 * ivu2s2 * p1p2 * p1p3 * pap2 +
                           4.0 * ivs * p1p3 * pap1 * pap2 +
                           4.0 * ivu2s2 * p1p3 * pap1 * pap2 +
                           4.0 * ivu2s2 * p1p3 * std::pow(pap2, 2) +
                           p1p2 * pap3 +
                           4.0 * ivu2s2 * std::pow(p1p2, 2) * pap3 -
                           4.0 * ivs * p1p2 * pap1 * pap3 -
                           4.0 * ivu2s2 * p1p2 * pap1 * pap3 -
                           4.0 * ivu2s2 * p1p2 * pap2 * pap3 -
                           4.0 * ivu2s2 * p1p3 * pap2 * papb +
                           4.0 * ivu2s2 * p1p2 * pap3 * papb +
                           2.0 * ivu2s2 * p1p3 * pap2 * pbp1 -
                           2.0 * ivu2s2 * p1p2 * pap3 * pbp1 +
                           4.0 * ivs * p1p3 * pap1 * pbp2 +
                           2.0 * ivu2s2 * p1p3 * pap2 * pbp2 -
                           2.0 * ivu2s2 * p1p2 * pap3 * pbp2 +
                           p2p3 * pap1 *
                               (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                2.0 * ivu2s2 *
                                    (-2.0 * p1p2 + 2.0 * pap1 + 2.0 * pap2 -
                                     2.0 * papb + pbp1 + pbp2)) -
                           4.0 * ivs * p1p2 * pap1 * pbp3 -
                           m1s * (p2p3 + 2.0 * ivs * pap3 * pbp2 -
                                  2.0 * ivs * pap2 * pbp3))))));

  return std::real(mtu);
}

double FI::MQttp() {
  std::complex<double> mttp = std::complex<double>(0.0, 0.0);

  const double mQs = m1s - 2.0 * pap1 - 1.0 / ivt1s1;

  const double p2p3p = .5 * (mQs - m2s);
  const double p1p3p = pap3 + pbp3 - p2p3p;
  const double pbp1p = p1p3p - pap1 + m1s + p1p2;

  const double zz =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  const double zzp =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1p) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  if (mQs > m2s && 2.0 * papb > std::pow(m1 + std::sqrt(mQs), 2)) {
    // pap1, pap2, pbp1p, pbp2p
    if (zz > 0.0 && zzp > 0.0) {
      mttp = DJAC *
             (-4.0 * ivs1s1 * ivs1s2 * ((LLLL + RRRR) + (LRLR + RLRL)) * p2p3p *
              (-2.0 * ivs * pbp1p +
               ivt1s2 * (m1s + pap1 * (1.0 - 4.0 * ivs * (pap1 + pbp1p))) +
               ivt1s1 * (m1s * (1.0 + 4.0 * ivt1s2 * pap1) -
                         pap1 * (-1.0 + 4.0 * ivs * (pap1 + pbp1p) +
                                 4.0 * ivt1s2 * (2.0 * pap1 - papb + pbp1p)))));
    }
  }

  // printf("DJAC = %f \n", DJAC);

  return std::real(mttp);
}

double FI::MQuup() {
  std::complex<double> muup = std::complex<double>(0.0, 0.0);

  const double mQs = m2s - 2.0 * pap2 - 1.0 / ivu2s1;

  const double p1p3p = .5 * (mQs - m1s);
  const double pbp1p = p1p3p - pap1 + m1s + p1p2;
  const double pbp2p = papb - pbp1p - pbp3;

  const double zz =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  const double zzp =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1p) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  if (mQs > m1s && 2.0 * papb > std::pow(m2 + std::sqrt(mQs), 2)) {
    // pap1, pap2, pbp1p, pbp2p
    if (zz > 0.0 && zzp > 0.0) {
      muup = DJAC *
             (-4.0 * ivs2s1 * ivs2s2 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3p *
              (-2.0 * ivs * pbp2p +
               ivu2s2 * (m2s + pap2 * (1.0 - 4.0 * ivs * (pap2 + pbp2p))) +
               ivu2s1 * (m2s * (1.0 + 4.0 * ivu2s2 * pap2) -
                         pap2 * (-1.0 + 4.0 * ivs * (pap2 + pbp2p) +
                                 4.0 * ivu2s2 * (2.0 * pap2 - papb + pbp2p)))));
    }
  }

  return std::real(muup);
}

// Real antiquark emission for: g + \bar{q} -> f + \bar{f} + \bar{q}.

double FI::MQBss() {
  std::complex<double> mss;

  mss =
      16.0 * ivs3v1 * ivs3v2 *
      (ivs * (2.0 * (LLLL + RRRR) * (p2p3 * (pap1 + pbp1) + pbp1 * pbp2) +
              2.0 * (LRLR + RLRL) * (pbp1 * pbp2 + p1p3 * (pap2 + pbp2)) +
              ((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 * (pap3 + 2.0 * pbp3)) +
       ivt3 * (-(((LLRL + RRLR) + (LRRR + RLLL)) * m1 * m2 *
                 (papb + 2.0 * pbp3 * (-1.0 + 2.0 * ivs * pbp3))) +
               2.0 * (LLLL + RRRR) *
                   (p1p3 * p2p3 + pap2 * pbp1 * (-1.0 + 2.0 * ivs * pbp3) +
                    p2p3 * (pbp1 - 2.0 * ivs * pap1 * pbp3 -
                            4.0 * ivs * pbp1 * pbp3)) +
               2.0 * (LRLR + RLRL) *
                   (pap1 * pbp2 * (-1.0 + 2.0 * ivs * pbp3) +
                    p1p3 * (p2p3 + pbp2 - 2.0 * ivs * pap2 * pbp3 -
                            4.0 * ivs * pbp2 * pbp3))));

  return std::real(mss);
}

double FI::MQBtt() {
  std::complex<double> mtt;

  mtt =
      2.0 *
      (ivt2s2 *
           (-4.0 * ivt2s1 * ivt3 * ((LLLL + RRRR) + (LRLR + RLRL)) * pap1 *
                pbp2 +
            ivs2s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivs * (p1p3 + pbp1) + ivt2s1 * (-p1p2 + p1p3 + pbp1)) *
                          pbp2 +
                      ivt3 * (2.0 * pap1 * pbp2 *
                                  (ivs * pbp3 + ivt2s1 * (-p2p3 + pbp3)) +
                              p1p3 * (p2p3 + 4.0 * ivt2s1 * p2p3 * pbp2 -
                                      2.0 * ivs * pap2 * pbp3 +
                                      pbp2 * (1.0 - 4.0 * ivs * pbp3 -
                                              2.0 * ivt2s1 * (pap2 - papb +
                                                              2.0 * pbp3))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivt2s1 * pbp2)) +
                      2.0 * (ivt2s1 * (p1p3 + pbp1) * pbp2 +
                             ivs * (pap2 * pbp1 + (p1p3 + pbp1) * pbp2) +
                             ivt3 * (pap1 * (-(p2p3 *
                                               (1.0 + 2.0 * ivt2s1 * pbp2)) +
                                             2.0 * (ivt2s1 * pbp2 +
                                                    ivs * (pap2 + pbp2)) *
                                                 pbp3) +
                                     p1p3 * (p2p3 + 4.0 * ivt2s1 * p2p3 * pbp2 -
                                             2.0 * ivs * pap2 * pbp3 +
                                             pbp2 * (1.0 - 4.0 * ivs * pbp3 -
                                                     2.0 * ivt2s1 *
                                                         (pap2 - papb +
                                                          2.0 * pbp3)))))))) +
       ivs2s2 *
           (-2.0 * ivs2s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3 *
                (-2.0 * ivs * pap2 +
                 ivt2s2 * (m2s + pbp2 * (1.0 - 4.0 * ivs * (pap2 + pbp2)))) +
            ivt2s1 *
                (2.0 * (LRLR + RLRL) *
                     ((ivs * (p1p3 + pbp1) + ivt2s2 * (-p1p2 + p1p3 + pbp1)) *
                          pbp2 -
                      ivs2s1 * p1p3 *
                          (m2s * (1.0 + 4.0 * ivt2s2 * pbp2) -
                           pbp2 * (-1.0 + 4.0 * ivs * (pap2 + pbp2) +
                                   4.0 * ivt2s2 * (pap2 - papb + 2.0 * pbp2))) +
                      ivt3 * (2.0 * pap1 * pbp2 *
                                  (ivs * pbp3 + ivt2s2 * (-p2p3 + pbp3)) +
                              p1p3 * (p2p3 + 4.0 * ivt2s2 * p2p3 * pbp2 -
                                      2.0 * ivs * pap2 * pbp3 +
                                      pbp2 * (1.0 - 4.0 * ivs * pbp3 -
                                              2.0 * ivt2s2 * (pap2 - papb +
                                                              2.0 * pbp3))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivt2s2 * pbp2)) +
                      2.0 * (ivs * pap2 * pbp1 + ivs * p1p3 * pbp2 +
                             ivt2s2 * p1p3 * pbp2 + ivs * pbp1 * pbp2 +
                             ivt2s2 * pbp1 * pbp2 -
                             ivs2s1 * p1p3 *
                                 (m2s * (1.0 + 4.0 * ivt2s2 * pbp2) -
                                  pbp2 * (-1.0 + 4.0 * ivs * (pap2 + pbp2) +
                                          4.0 * ivt2s2 *
                                              (pap2 - papb + 2.0 * pbp2))) +
                             ivt3 * (pap1 * (-(p2p3 *
                                               (1.0 + 2.0 * ivt2s2 * pbp2)) +
                                             2.0 * (ivt2s2 * pbp2 +
                                                    ivs * (pap2 + pbp2)) *
                                                 pbp3) +
                                     p1p3 * (p2p3 + 4.0 * ivt2s2 * p2p3 * pbp2 -
                                             2.0 * ivs * pap2 * pbp3 +
                                             pbp2 * (1.0 - 4.0 * ivs * pbp3 -
                                                     2.0 * ivt2s2 *
                                                         (pap2 - papb +
                                                          2.0 * pbp3)))))))));

  return std::real(mtt);
}

double FI::MQBuu() {
  std::complex<double> muu;

  muu =
      2.0 *
      (ivu1s2 *
           (-4.0 * ivt3 * ivu1s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * pap2 *
                pbp1 +
            ivs1s1 *
                (2.0 * (LRLR + RLRL) *
                     (pbp1 * (ivs * (p2p3 + pbp2) +
                              ivu1s1 * (-p1p2 + p2p3 + pbp2)) +
                      ivt3 * (p1p3 * (p2p3 + 4.0 * ivu1s1 * p2p3 * pbp1 -
                                      2.0 * ivu1s1 * pap2 * pbp1) +
                              2.0 * (ivs + ivu1s1) * pap2 * pbp1 * pbp3 +
                              p2p3 * (-2.0 * ivs * pap1 * pbp3 +
                                      pbp1 * (1.0 - 4.0 * ivs * pbp3 -
                                              2.0 * ivu1s1 * (pap1 - papb +
                                                              2.0 * pbp3))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivu1s1 * pbp1)) +
                      2.0 *
                          (ivs * p2p3 * pbp1 + ivu1s1 * p2p3 * pbp1 +
                           ivs * pap1 * pbp2 + ivs * pbp1 * pbp2 +
                           ivu1s1 * pbp1 * pbp2 +
                           ivt3 * (p1p3 * (p2p3 + 4.0 * ivu1s1 * p2p3 * pbp1 -
                                           pap2 * (1.0 + 2.0 * ivu1s1 * pbp1)) +
                                   2.0 * pap2 *
                                       (ivu1s1 * pbp1 + ivs * (pap1 + pbp1)) *
                                       pbp3 +
                                   p2p3 * (-2.0 * ivs * pap1 * pbp3 +
                                           pbp1 * (1.0 - 4.0 * ivs * pbp3 -
                                                   2.0 * ivu1s1 *
                                                       (pap1 - papb +
                                                        2.0 * pbp3)))))))) +
       ivs1s2 *
           (-2.0 * ivs1s1 * ((LLLL + RRRR) + (LRLR + RLRL)) * p2p3 *
                (-2.0 * ivs * pap1 +
                 ivu1s2 * (m1s + pbp1 * (1.0 - 4.0 * ivs * (pap1 + pbp1)))) +
            ivu1s1 *
                (2.0 * (LRLR + RLRL) *
                     (-(ivs1s1 * p2p3 *
                        (m1s * (1.0 + 4.0 * ivu1s2 * pbp1) -
                         pbp1 * (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                 4.0 * ivu1s2 * (pap1 - papb + 2.0 * pbp1)))) +
                      pbp1 * (ivs * (p2p3 + pbp2) +
                              ivu1s2 * (-p1p2 + p2p3 + pbp2)) +
                      ivt3 * (p1p3 * (p2p3 + 4.0 * ivu1s2 * p2p3 * pbp1 -
                                      2.0 * ivu1s2 * pap2 * pbp1) +
                              2.0 * (ivs + ivu1s2) * pap2 * pbp1 * pbp3 +
                              p2p3 * (-2.0 * ivs * pap1 * pbp3 +
                                      pbp1 * (1.0 - 4.0 * ivs * pbp3 -
                                              2.0 * ivu1s2 * (pap1 - papb +
                                                              2.0 * pbp3))))) +
                 (LLLL + RRRR) *
                     (-(p1p2 * (1.0 + 2.0 * ivu1s2 * pbp1)) +
                      2.0 *
                          (ivs * p2p3 * pbp1 + ivu1s2 * p2p3 * pbp1 -
                           ivs1s1 * p2p3 *
                               (m1s * (1.0 + 4.0 * ivu1s2 * pbp1) -
                                pbp1 * (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                        4.0 * ivu1s2 *
                                            (pap1 - papb + 2.0 * pbp1))) +
                           ivs * pap1 * pbp2 + ivs * pbp1 * pbp2 +
                           ivu1s2 * pbp1 * pbp2 +
                           ivt3 * (p1p3 * (p2p3 + 4.0 * ivu1s2 * p2p3 * pbp1 -
                                           pap2 * (1.0 + 2.0 * ivu1s2 * pbp1)) +
                                   2.0 * pap2 *
                                       (ivu1s2 * pbp1 + ivs * (pap1 + pbp1)) *
                                       pbp3 +
                                   p2p3 * (-2.0 * ivs * pap1 * pbp3 +
                                           pbp1 * (1.0 - 4.0 * ivs * pbp3 -
                                                   2.0 * ivu1s2 *
                                                       (pap1 - papb +
                                                        2.0 * pbp3)))))))));

  return std::real(muu);
}

double FI::MQBst() {
  std::complex<double> mst;

  mst =
      4.0 * ivs3v1 *
      (2.0 * ivt2s2 *
           (ivs * ((LRLR + RLRL) * (p1p3 + pbp1) * pbp2 +
                   (LLRL + RRLR) * m1 * m2 * pbp3) +
            ivt3 * ((LLRL + RRLR) * m1 * m2 *
                        (-papb + pbp3 - 2.0 * ivs * std::pow(pbp3, 2)) +
                    (LRLR + RLRL) *
                        (2.0 * pap1 * pbp2 * (-1.0 + ivs * pbp3) +
                         p1p3 * (p2p3 + pbp2 - 2.0 * ivs * pap2 * pbp3 -
                                 4.0 * ivs * pbp2 * pbp3)))) +
       ivs2s2 *
           (2.0 * (ivs * ((LRLR + RLRL) *
                              (2.0 * p1p3 * pap2 + p1p3 * pbp2 + pbp1 * pbp2) +
                          (LLRL + RRLR) * m1 * m2 * (pap3 + pbp3)) +
                   ivt3 * ((LLRL + RRLR) * m1 * m2 * pbp3 *
                               (1.0 - 2.0 * ivs * pbp3) +
                           (LRLR + RLRL) *
                               (2.0 * ivs * pap1 * pbp2 * pbp3 +
                                p1p3 * (p2p3 + pbp2 - 2.0 * ivs * pap2 * pbp3 -
                                        4.0 * ivs * pbp2 * pbp3)))) +
            ivt2s2 *
                (-((LLRL + RRLR) * m1 * m2 *
                   (p2p3 * (1.0 + 2.0 * ivt3 * (papb - 2.0 * pbp3)) +
                    pbp2 * (1.0 - 2.0 * ivs * (pap3 + 2.0 * pbp3)) +
                    2.0 * pbp3 * (-(ivs * pap2) +
                                  ivt3 * (pap2 - 2.0 * papb + 2.0 * pbp3)))) +
                 (LRLR + RLRL) *
                     (-2.0 * m2s * p1p3 +
                      2.0 * pbp2 * (-p1p2 + 4.0 * ivs * p1p3 * pap2 + pbp1 +
                                    4.0 * ivs * p1p3 * pbp2 +
                                    2.0 * ivt3 * (p1p3 * (2.0 * p2p3 - pap2 +
                                                          papb - 2.0 * pbp3) +
                                                  pap1 * (-p2p3 + pbp3)))))));

  return std::real(mst);
}

double FI::MQBsu() {
  std::complex<double> msu;

  msu =
      4.0 * ivs3v1 *
      (-2.0 * ivu1s2 *
           (ivs * ((LLRL + RRLR) * pbp1 * (p2p3 + pbp2) +
                   (LRLR + RLRL) * m1 * m2 * pbp3) +
            ivt3 * ((LRLR + RLRL) * m1 * m2 *
                        (-papb + pbp3 - 2.0 * ivs * std::pow(pbp3, 2)) +
                    (LLRL + RRLR) *
                        (p1p3 * p2p3 + 2.0 * pap2 * pbp1 * (-1.0 + ivs * pbp3) +
                         p2p3 * (pbp1 - 2.0 * ivs * pap1 * pbp3 -
                                 4.0 * ivs * pbp1 * pbp3)))) +
       ivs1s2 *
           (ivu1s2 *
                ((LRLR + RLRL) * m1 * m2 *
                     (p1p3 * (1.0 + 2.0 * ivt3 * (papb - 2.0 * pbp3)) +
                      pbp1 * (1.0 - 2.0 * ivs * (pap3 + 2.0 * pbp3)) +
                      2.0 * pbp3 * (-(ivs * pap1) +
                                    ivt3 * (pap1 - 2.0 * papb + 2.0 * pbp3))) +
                 2.0 * (LLRL + RRLR) *
                     (m1s * p2p3 -
                      pbp1 *
                          (-p1p2 + 4.0 * ivs * p2p3 * pap1 +
                           4.0 * ivs * p2p3 * pbp1 + pbp2 +
                           2.0 * ivt3 * (2.0 * p1p3 * p2p3 - p2p3 * pap1 -
                                         p1p3 * pap2 + p2p3 * papb -
                                         2.0 * p2p3 * pbp3 + pap2 * pbp3)))) -
            2.0 * (ivs * ((LLRL + RRLR) *
                              (p2p3 * (2.0 * pap1 + pbp1) + pbp1 * pbp2) +
                          (LRLR + RLRL) * m1 * m2 * (pap3 + pbp3)) +
                   ivt3 * ((LRLR + RLRL) * m1 * m2 * pbp3 *
                               (1.0 - 2.0 * ivs * pbp3) +
                           (LLRL + RRLR) *
                               (p1p3 * p2p3 + 2.0 * ivs * pap2 * pbp1 * pbp3 +
                                p2p3 * (pbp1 - 2.0 * ivs * pap1 * pbp3 -
                                        4.0 * ivs * pbp1 * pbp3))))));

  return std::real(msu);
}

double FI::MQBtu() {
  std::complex<double> mtu;

  mtu =
      2.0 *
      (ivs1s2 *
           (ivs2s1 *
                (-2.0 * ivs * ((LRLR + RLRL) * m1 * m2 * pap3 +
                               (LLLL + RRRR) *
                                   (p2p3 * pap1 + p1p3 * pap2 - p1p2 * pap3)) +
                 ivu1s2 *
                     ((LLLL + RRRR) *
                          (-4.0 * ivs * p1p3 * pap2 * pbp1 +
                           4.0 * ivs * p1p2 * pap3 * pbp1 + p1p3 * pbp2 +
                           4.0 * ivt2s1 * p1p2 * p1p3 * pbp2 -
                           2.0 * ivt2s1 * p1p3 * pap1 * pbp2 -
                           2.0 * ivt2s1 * p1p3 * pap2 * pbp2 +
                           4.0 * ivt2s1 * p1p3 * papb * pbp2 -
                           4.0 * ivs * p1p3 * pbp1 * pbp2 -
                           4.0 * ivt2s1 * p1p3 * pbp1 * pbp2 -
                           4.0 * ivt2s1 * p1p3 * std::pow(pbp2, 2) -
                           p2p3 * pbp1 *
                               (-1.0 + 4.0 * ivs * (pap1 + pbp1) +
                                2.0 * ivt2s1 *
                                    (-2.0 * p1p2 + pap1 + pap2 - 2.0 * papb +
                                     2.0 * pbp1 + 2.0 * pbp2)) -
                           p1p2 * pbp3 -
                           4.0 * ivt2s1 * std::pow(p1p2, 2) * pbp3 +
                           2.0 * ivt2s1 * p1p2 * pap1 * pbp3 +
                           2.0 * ivt2s1 * p1p2 * pap2 * pbp3 -
                           4.0 * ivt2s1 * p1p2 * papb * pbp3 +
                           4.0 * ivs * p1p2 * pbp1 * pbp3 +
                           4.0 * ivt2s1 * p1p2 * pbp1 * pbp3 +
                           4.0 * ivt2s1 * p1p2 * pbp2 * pbp3 +
                           m1s * (p2p3 - 2.0 * ivs * pap3 * pbp2 +
                                  2.0 * ivs * pap2 * pbp3)) +
                      (LRLR + RLRL) * m1 * m2 *
                          (p1p3 +
                           (1.0 +
                            ivt2s1 * (4.0 * p1p2 -
                                      2.0 * (pap1 + pap2 +
                                             2.0 * (-papb + pbp1 + pbp2)))) *
                               pbp3 -
                           2.0 * ivs *
                               (pap3 * pbp1 + (pap1 + 2.0 * pbp1) * pbp3))) +
                 ivt2s1 * ((LLLL + RRRR) *
                               (p1p3 * pbp2 - 4.0 * ivs * p1p3 * pap2 * pbp2 +
                                4.0 * ivs * p1p2 * pap3 * pbp2 -
                                4.0 * ivs * p1p3 * std::pow(pbp2, 2) +
                                p2p3 * (pbp1 - 4.0 * ivs * pap1 * pbp2 -
                                        4.0 * ivs * pbp1 * pbp2) -
                                p1p2 * pbp3 + 4.0 * ivs * p1p2 * pbp2 * pbp3 +
                                m2s * (p1p3 - 2.0 * ivs * pap3 * pbp1 +
                                       2.0 * ivs * pap1 * pbp3)) +
                           (LRLR + RLRL) * m1 * m2 *
                               (p2p3 + pbp3 -
                                2.0 * ivs * (pap3 * pbp2 + pap2 * pbp3 +
                                             2.0 * pbp2 * pbp3)))) -
            ivt2s1 *
                (-(ivu1s2 * (LRLR + RLRL) * m1 * m2 * pbp1) +
                 2.0 * ivs * (LLLL + RRRR) * p2p3 * pbp1 +
                 ivu1s2 * (LLLL + RRRR) * p2p3 * pbp1 -
                 ivu1s2 * (LLLL + RRRR) * m1s * pbp2 +
                 ivu1s2 * (LLLL + RRRR) * p1p3 * pbp2 +
                 2.0 * ivs * (LLLL + RRRR) * pap1 * pbp2 +
                 2.0 * ivs * (LLLL + RRRR) * pbp1 * pbp2 +
                 2.0 * ivu1s2 * (LLLL + RRRR) * pbp1 * pbp2 +
                 2.0 * ivs * (LRLR + RLRL) * m1 * m2 * pbp3 +
                 ivu1s2 * (LRLR + RLRL) * m1 * m2 * pbp3 -
                 2.0 * ivs * (LLLL + RRRR) * p1p2 * pbp3 -
                 ivu1s2 * (LLLL + RRRR) * p1p2 * pbp3 +
                 2.0 * ivt3 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (pbp3 - 2.0 * ivs * std::pow(pbp3, 2) -
                           ivu1s2 * (p1p3 * (papb - 2.0 * pbp3) +
                                     pbp3 * (pap1 - 2.0 * papb + 2.0 * pbp3))) +
                      (LLLL + RRRR) *
                          (2.0 * ivu1s2 * std::pow(p1p3, 2) * pbp2 -
                           p2p3 * (pap1 - ivu1s2 * m1s * papb +
                                   2.0 * ivu1s2 * pap1 * pbp1 +
                                   2.0 * ivs * pbp1 * pbp3 +
                                   2.0 * ivu1s2 * pbp1 * pbp3) +
                           p1p3 *
                               (p2p3 + 2.0 * ivu1s2 * p2p3 * pbp1 -
                                2.0 * (ivu1s2 * p1p2 + ivs * pap2) * pbp3 +
                                pbp2 * (1.0 - 2.0 * ivs * pbp3 -
                                        2.0 * ivu1s2 * (pap1 - papb + pbp3))) +
                           pbp3 * (pap2 * (-(ivu1s2 * (m1s - 2.0 * pbp1)) +
                                           2.0 * ivs * (pap1 + pbp1)) +
                                   p1p2 * (-1.0 + 2.0 * ivs * pbp3 +
                                           2.0 * ivu1s2 *
                                               (pap1 - papb + pbp3))))))) +
       ivu1s2 *
           (2.0 * ivt2s1 * ivt3 *
                ((LRLR + RLRL) * m1 * m2 * papb +
                 (LLLL + RRRR) * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) -
            ivs2s1 *
                (-(ivt2s1 * (LLLL + RRRR) * m2s * pbp1) +
                 ivt2s1 * (LLLL + RRRR) * p2p3 * pbp1 +
                 2.0 * ivs * (LLLL + RRRR) * pap2 * pbp1 -
                 ivt2s1 * (LRLR + RLRL) * m1 * m2 * pbp2 +
                 2.0 * ivs * (LLLL + RRRR) * p1p3 * pbp2 +
                 ivt2s1 * (LLLL + RRRR) * p1p3 * pbp2 +
                 2.0 * ivs * (LLLL + RRRR) * pbp1 * pbp2 +
                 2.0 * ivt2s1 * (LLLL + RRRR) * pbp1 * pbp2 +
                 2.0 * ivs * (LRLR + RLRL) * m1 * m2 * pbp3 +
                 ivt2s1 * (LRLR + RLRL) * m1 * m2 * pbp3 -
                 2.0 * ivs * (LLLL + RRRR) * p1p2 * pbp3 -
                 ivt2s1 * (LLLL + RRRR) * p1p2 * pbp3 +
                 2.0 * ivt3 *
                     ((LRLR + RLRL) * m1 * m2 *
                          (pbp3 - 2.0 * ivs * std::pow(pbp3, 2) -
                           ivt2s1 * (p2p3 * (papb - 2.0 * pbp3) +
                                     pbp3 * (pap2 - 2.0 * papb + 2.0 * pbp3))) +
                      (LLLL + RRRR) *
                          (2.0 * ivt2s1 * std::pow(p2p3, 2) * pbp1 +
                           p1p3 * (p2p3 - pap2 + ivt2s1 * m2s * papb +
                                   2.0 * ivt2s1 * p2p3 * pbp2 -
                                   2.0 * ivt2s1 * pap2 * pbp2 -
                                   2.0 * ivs * pbp2 * pbp3 -
                                   2.0 * ivt2s1 * pbp2 * pbp3) +
                           p2p3 *
                               (-2.0 * (ivt2s1 * p1p2 + ivs * pap1) * pbp3 +
                                pbp1 * (1.0 - 2.0 * ivs * pbp3 -
                                        2.0 * ivt2s1 * (pap2 - papb + pbp3))) +
                           pbp3 * (pap1 * (-(ivt2s1 * (m2s - 2.0 * pbp2)) +
                                           2.0 * ivs * (pap2 + pbp2)) +
                                   p1p2 * (-1.0 + 2.0 * ivs * pbp3 +
                                           2.0 * ivt2s1 *
                                               (pap2 - papb + pbp3))))))));

  return std::real(mtu);
}

double FI::MQBttp() {
  std::complex<double> mttp = std::complex<double>(0.0, 0.0);

  const double mQs = m1s - 2.0 * pap1 - 1.0 / ivt1s1;

  const double p1p3p = .5 * (mQs - m1s);
  const double pbp1p = p1p3p - pap1 + m1s + p1p2;
  const double pbp2p = papb - pbp1p - pbp3;

  const double ivt2s1p = 1.0 / (m2s - 2.0 * pbp2p - mQs);
  const double ivt2s2p = ivt2s1p;

  const double zz =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  const double zzp =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1p) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  if (mQs > m1s && 2.0 * papb > std::pow(m2 + std::sqrt(mQs), 2)) {
    // pap1, pap2, pbp1p, pbp2p
    if (zz > 0.0 && zzp > 0.0) {
      mttp =
          DJAC *
          (-4.0 * ivs2s1 * ivs2s2 * ((LLLL + RRRR) + (LRLR + RLRL)) * p1p3p *
           (-2.0 * ivs * pap2 +
            ivt2s2p * (m2s + pbp2p * (1.0 - 4.0 * ivs * (pap2 + pbp2p))) +
            ivt2s1p * (m2s * (1.0 + 4.0 * ivt2s2p * pbp2p) -
                       pbp2p * (-1.0 + 4.0 * ivs * (pap2 + pbp2p) +
                                4.0 * ivt2s2p * (pap2 - papb + 2.0 * pbp2p)))));
    }
  }

  return std::real(mttp);
}

double FI::MQBuup() {
  std::complex<double> muup = std::complex<double>(0.0, 0.0);

  const double mQs = m2s - 2.0 * pap2 - 1.0 / ivu2s1;
  const double p2p3p = .5 * (mQs - m2s);
  const double p1p3p = pap3 + pbp3 - p2p3p;
  const double pbp1p = p1p3p - pap1 + m1s + p1p2;
  const double ivu1s1p = 1.0 / (m1s - 2.0 * pbp1p - mQs);
  const double ivu1s2p = ivu1s1p;

  const double zz =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  const double zzp =
      1.0 +
      std::pow((pap3 - papb) * (-(m1s * papb) + m2s * papb - 2.0 * pap2 * papb +
                                2.0 * (-pap3 + papb) * pbp1p) +
                   2.0 * ((pap3 - papb) * papb + pap1 * (pap3 + papb)) * pbp3,
               2) /
          (8.0 * papb * pbp3 *
           (m1s * pap3 * (papb - pap3) * pap2 +
            pap1 * pap3 *
                (m2s * (-pap3 + papb) - 2.0 * pap2 * (pap3 - papb + pbp3))));

  if (mQs > m2s && 2.0 * papb > std::pow(m1 + std::sqrt(mQs), 2)) {
    // pap1, pap2, pbp1p, pbp2p
    if (zz > 0.0 && zzp > 0.0) {
      muup =
          DJAC *
          (-4.0 * ivs1s1 * ivs1s2 * ((LLLL + RRRR) + (LRLR + RLRL)) * p2p3p *
           (-2.0 * ivs * pap1 +
            ivu1s2p * (m1s + pbp1p * (1.0 - 4.0 * ivs * (pap1 + pbp1p))) +
            ivu1s1p * (m1s * (1.0 + 4.0 * ivu1s2p * pbp1p) -
                       pbp1p * (-1.0 + 4.0 * ivs * (pap1 + pbp1p) +
                                4.0 * ivu1s2p * (pap1 - papb + 2.0 * pbp1p)))));
    }
  }

  // printf("djac = %f \n", DJAC);

  return std::real(muup);
}

// Real gluon emission for: q + \bar{q} -> sl + \bar{sl} + g

double FI::MGssSL() {
  std::complex<double> mss;

  mss =
      ivt3 * ivs3v1 * ivs3v2 *
          (16.0 * LRLR * papb * p1p2 - 8.0 * LRLR * pbp3 * p1p2 -
           8.0 * LRLR * pbp2 * p2p3 + 8.0 * LRLR * pbp2 * p1p3 +
           8.0 * LRLR * pbp1 * p2p3 - 8.0 * LRLR * pbp1 * p1p3 +
           8.0 * LRLR * pap2 * pbp2 - 8.0 * LRLR * pap2 * pbp1 -
           8.0 * LRLR * pap2 * pap2 - 8.0 * LRLR * pap1 * pbp2 +
           8.0 * LRLR * pap1 * pbp1 + 16.0 * LRLR * pap1 * pap2 -
           8.0 * LRLR * pap1 * pap1 + 16.0 * LLLL * papb * p1p2 -
           8.0 * LLLL * pbp3 * p1p2 - 8.0 * LLLL * pbp2 * p2p3 +
           8.0 * LLLL * pbp2 * p1p3 + 8.0 * LLLL * pbp1 * p2p3 -
           8.0 * LLLL * pbp1 * p1p3 + 8.0 * LLLL * pap2 * pbp2 -
           8.0 * LLLL * pap2 * pbp1 - 8.0 * LLLL * pap2 * pap2 -
           8.0 * LLLL * pap1 * pbp2 + 8.0 * LLLL * pap1 * pbp1 +
           16. * LLLL * pap1 * pap2 - 8.0 * LLLL * pap1 * pap1 -
           8.0 * m2 * m2 * LRLR * papb + 4.0 * m2 * m2 * LRLR * pbp3 -
           8.0 * m2 * m2 * LLLL * papb + 4.0 * m2 * m2 * LLLL * pbp3 -
           8.0 * m1 * m1 * LRLR * papb + 4.0 * m1 * m1 * LRLR * pbp3 -
           8.0 * m1 * m1 * LLLL * papb + 4.0 * m1 * m1 * LLLL * pbp3)

      +
      ivu3 * ivs3v1 * ivs3v2 *
          (16.0 * LRLR * papb * p1p2 - 8.0 * LRLR * pbp2 * pbp2 +
           16.0 * LRLR * pbp1 * pbp2 - 8.0 * LRLR * pbp1 * pbp1 -
           8.0 * LRLR * pap3 * p1p2 - 8.0 * LRLR * pap2 * p2p3 +
           8. * LRLR * pap2 * p1p3 + 8.0 * LRLR * pap2 * pbp2 -
           8.0 * LRLR * pap2 * pbp1 + 8.0 * LRLR * pap1 * p2p3 -
           8.0 * LRLR * pap1 * p1p3 - 8.0 * LRLR * pap1 * pbp2 +
           8.0 * LRLR * pap1 * pbp1 + 16. * LLLL * papb * p1p2 -
           8.0 * LLLL * pbp2 * pbp2 + 16.0 * LLLL * pbp1 * pbp2 -
           8.0 * LLLL * pbp1 * pbp1 - 8.0 * LLLL * pap3 * p1p2 -
           8.0 * LLLL * pap2 * p2p3 + 8.0 * LLLL * pap2 * p1p3 +
           8.0 * LLLL * pap2 * pbp2 - 8.0 * LLLL * pap2 * pbp1 +
           8.0 * LLLL * pap1 * p2p3 - 8.0 * LLLL * pap1 * p1p3 -
           8. * LLLL * pap1 * pbp2 + 8.0 * LLLL * pap1 * pbp1 -
           8.0 * m2 * m2 * LRLR * papb + 4.0 * m2 * m2 * LRLR * pap3 -
           8.0 * m2 * m2 * LLLL * papb + 4.0 * m2 * m2 * LLLL * pap3 -
           8.0 * m1 * m1 * LRLR * papb + 4. * m1 * m1 * LRLR * pap3 -
           8.0 * m1 * m1 * LLLL * papb + 4.0 * m1 * m1 * LLLL * pap3)

      +
      ivu3 * ivt3 * ivs3v1 * ivs3v2 *
          (32.0 * LRLR * papb * papb * p1p2 - 16.0 * LRLR * pbp2 * papb * p2p3 +
           16.0 * LRLR * pbp2 * papb * p1p3 + 16.0 * LRLR * pbp1 * papb * p2p3 -
           16.0 * LRLR * pbp1 * papb * p1p3 - 16.0 * LRLR * pap2 * papb * p2p3 +
           16.0 * LRLR * pap2 * papb * p1p3 + 32. * LRLR * pap2 * pbp2 * papb -
           32.0 * LRLR * pap2 * pbp1 * papb + 16.0 * LRLR * pap1 * papb * p2p3 -
           16.0 * LRLR * pap1 * papb * p1p3 - 32.0 * LRLR * pap1 * pbp2 * papb +
           32.0 * LRLR * pap1 * pbp1 * papb + 32.0 * LLLL * papb * papb * p1p2 -
           16.0 * LLLL * pbp2 * papb * p2p3 + 16.0 * LLLL * pbp2 * papb * p1p3 +
           16.0 * LLLL * pbp1 * papb * p2p3 - 16.0 * LLLL * pbp1 * papb * p1p3 -
           16. * LLLL * pap2 * papb * p2p3 + 16.0 * LLLL * pap2 * papb * p1p3 +
           32.0 * LLLL * pap2 * pbp2 * papb - 32.0 * LLLL * pap2 * pbp1 * papb +
           16.0 * LLLL * pap1 * papb * p2p3 - 16.0 * LLLL * pap1 * papb * p1p3 -
           32.0 * LLLL * pap1 * pbp2 * papb + 32.0 * LLLL * pap1 * pbp1 * papb -
           16. * m2 * m2 * LRLR * papb * papb -
           16.0 * m2 * m2 * LLLL * papb * papb -
           16.0 * m1 * m1 * LRLR * papb * papb -
           16. * m1 * m1 * LLLL * papb * papb);

  return std::real(mss);
}

// Real quark emission for: q + g -> sl + \bar{sl} + q

double FI::MQssSL() {
  std::complex<double> mss;

  mss =
      ivs * ivs3v1 * ivs3v2 *
          (8.0 * pap1 * p1p3 * LRLR + 8.0 * pap1 * p1p3 * LLLL -
           8. * pap1 * p2p3 * LRLR - 8.0 * pap1 * p2p3 * LLLL -
           8.0 * pap2 * p1p3 * LRLR - 8.0 * pap2 * p1p3 * LLLL +
           8.0 * pap2 * p2p3 * LRLR + 8.0 * pap2 * p2p3 * LLLL +
           8.0 * pap3 * p1p2 * LRLR + 8.0 * pap3 * p1p2 * LLLL -
           4.0 * pap3 * m2 * m2 * LRLR - 4.0 * pap3 * m2 * m2 * LLLL -
           4. * pap3 * m1 * m1 * LRLR - 4.0 * pap3 * m1 * m1 * LLLL -
           16.0 * pbp1 * pbp2 * LRLR - 16.0 * pbp1 * pbp2 * LLLL +
           8.0 * pbp1 * p1p3 * LRLR + 8.0 * pbp1 * p1p3 * LLLL -
           8.0 * pbp1 * p2p3 * LRLR - 8.0 * pbp1 * p2p3 * LLLL +
           8.0 * pbp1 * pbp1 * LRLR + 8.0 * pbp1 * pbp1 * LLLL -
           8. * pbp2 * p1p3 * LRLR - 8.0 * pbp2 * p1p3 * LLLL +
           8.0 * pbp2 * p2p3 * LRLR + 8.0 * pbp2 * p2p3 * LLLL +
           8.0 * pbp2 * pbp2 * LRLR + 8.0 * pbp2 * pbp2 * LLLL +
           16.0 * pbp3 * p1p2 * LRLR + 16.0 * pbp3 * p1p2 * LLLL -
           8.0 * pbp3 * m2 * m2 * LRLR - 8.0 * pbp3 * m2 * m2 * LLLL -
           8.0 * pbp3 * m1 * m1 * LRLR - 8.0 * pbp3 * m1 * m1 * LLLL)

      +
      ivt3 * ivs3v1 * ivs3v2 *
          (-8.0 * papb * p1p2 * LRLR - 8.0 * papb * p1p2 * LLLL +
           4. * papb * m2 * m2 * LRLR + 4.0 * papb * m2 * m2 * LLLL +
           4.0 * papb * m1 * m1 * LRLR + 4.0 * papb * m1 * m1 * LLLL -
           8.0 * pap1 * pbp1 * LRLR - 8.0 * pap1 * pbp1 * LLLL +
           8.0 * pap1 * pbp2 * LRLR + 8.0 * pap1 * pbp2 * LLLL +
           8.0 * pap2 * pbp1 * LRLR + 8.0 * pap2 * pbp1 * LLLL -
           8.0 * pap2 * pbp2 * LRLR - 8.0 * pap2 * pbp2 * LLLL +
           8.0 * pbp1 * p1p3 * LRLR + 8. * pbp1 * p1p3 * LLLL -
           8.0 * pbp1 * p2p3 * LRLR - 8.0 * pbp1 * p2p3 * LLLL -
           8.0 * pbp2 * p1p3 * LRLR - 8.0 * pbp2 * p1p3 * LLLL +
           8.0 * pbp2 * p2p3 * LRLR + 8.0 * pbp2 * p2p3 * LLLL +
           16.0 * pbp3 * p1p2 * LRLR + 16.0 * pbp3 * p1p2 * LLLL -
           8.0 * pbp3 * m2 * m2 * LRLR - 8.0 * pbp3 * m2 * m2 * LLLL -
           8.0 * pbp3 * m1 * m1 * LRLR - 8.0 * pbp3 * m1 * m1 * LLLL -
           16. * p1p3 * p2p3 * LRLR - 16.0 * p1p3 * p2p3 * LLLL +
           8.0 * p1p3 * p1p3 * LRLR + 8.0 * p1p3 * p1p3 * LLLL +
           8.0 * p2p3 * p2p3 * LRLR + 8.0 * p2p3 * p2p3 * LLLL)

      +
      ivt3 * ivs * ivs3v1 * ivs3v2 *
          (16.0 * pap1 * pbp1 * pbp3 * LRLR + 16.0 * pap1 * pbp1 * pbp3 * LLLL -
           16.0 * pap1 * pbp2 * pbp3 * LRLR - 16.0 * pap1 * pbp2 * pbp3 * LLLL -
           16.0 * pap1 * pbp3 * p1p3 * LRLR - 16.0 * pap1 * pbp3 * p1p3 * LLLL +
           16.0 * pap1 * pbp3 * p2p3 * LRLR + 16.0 * pap1 * pbp3 * p2p3 * LLLL -
           16.0 * pap2 * pbp1 * pbp3 * LRLR - 16.0 * pap2 * pbp1 * pbp3 * LLLL +
           16.0 * pap2 * pbp2 * pbp3 * LRLR + 16. * pap2 * pbp2 * pbp3 * LLLL +
           16.0 * pap2 * pbp3 * p1p3 * LRLR + 16.0 * pap2 * pbp3 * p1p3 * LLLL -
           16.0 * pap2 * pbp3 * p2p3 * LRLR - 16.0 * pap2 * pbp3 * p2p3 * LLLL -
           32.0 * pbp1 * pbp3 * p1p3 * LRLR - 32.0 * pbp1 * pbp3 * p1p3 * LLLL +
           32.0 * pbp1 * pbp3 * p2p3 * LRLR + 32.0 * pbp1 * pbp3 * p2p3 * LLLL +
           32.0 * pbp2 * pbp3 * p1p3 * LRLR + 32.0 * pbp2 * pbp3 * p1p3 * LLLL -
           32.0 * pbp2 * pbp3 * p2p3 * LRLR - 32. * pbp2 * pbp3 * p2p3 * LLLL -
           32.0 * pbp3 * pbp3 * p1p2 * LRLR - 32.0 * pbp3 * pbp3 * p1p2 * LLLL +
           16.0 * pbp3 * pbp3 * m2 * m2 * LRLR +
           16.0 * pbp3 * pbp3 * m2 * m2 * LLLL +
           16.0 * pbp3 * pbp3 * m1 * m1 * LRLR +
           16.0 * pbp3 * pbp3 * m1 * m1 * LLLL);

  return std::real(mss);
}

// Real quark emission for: \bar{q} + g -> sl + \bar{sl} + \bar{q}
double FI::MQBssSL() {
  std::complex<double> mss;

  mss =
      +ivs * ivs3v1 * ivs3v2 *
          (-16.0 * pap1 * pap2 * LRLR - 16.0 * pap1 * pap2 * LLLL +
           8.0 * pap1 * p1p3 * LRLR + 8.0 * pap1 * p1p3 * LLLL -
           8.0 * pap1 * p2p3 * LRLR - 8. * pap1 * p2p3 * LLLL +
           8.0 * pap1 * pap1 * LRLR + 8.0 * pap1 * pap1 * LLLL -
           8.0 * pap2 * p1p3 * LRLR - 8.0 * pap2 * p1p3 * LLLL +
           8.0 * pap2 * p2p3 * LRLR + 8.0 * pap2 * p2p3 * LLLL +
           8.0 * pap2 * pap2 * LRLR + 8.0 * pap2 * pap2 * LLLL +
           16.0 * pap3 * p1p2 * LRLR + 16.0 * pap3 * p1p2 * LLLL -
           8.0 * pap3 * m2 * m2 * LRLR - 8.0 * pap3 * m2 * m2 * LLLL -
           8.0 * pap3 * m1 * m1 * LRLR - 8.0 * pap3 * m1 * m1 * LLLL +
           8.0 * pbp1 * p1p3 * LRLR + 8.0 * pbp1 * p1p3 * LLLL -
           8.0 * pbp1 * p2p3 * LRLR - 8.0 * pbp1 * p2p3 * LLLL -
           8.0 * pbp2 * p1p3 * LRLR - 8. * pbp2 * p1p3 * LLLL +
           8.0 * pbp2 * p2p3 * LRLR + 8.0 * pbp2 * p2p3 * LLLL +
           8.0 * pbp3 * p1p2 * LRLR + 8.0 * pbp3 * p1p2 * LLLL -
           4.0 * pbp3 * m2 * m2 * LRLR - 4.0 * pbp3 * m2 * m2 * LLLL -
           4.0 * pbp3 * m1 * m1 * LRLR - 4.0 * pbp3 * m1 * m1 * LLLL)

      +
      ivu3 * ivs3v1 * ivs3v2 *
          (-8.0 * papb * p1p2 * LRLR - 8.0 * papb * p1p2 * LLLL +
           4. * papb * m2 * m2 * LRLR + 4.0 * papb * m2 * m2 * LLLL +
           4.0 * papb * m1 * m1 * LRLR + 4.0 * papb * m1 * m1 * LLLL -
           8.0 * pap1 * pbp1 * LRLR - 8.0 * pap1 * pbp1 * LLLL +
           8.0 * pap1 * pbp2 * LRLR + 8.0 * pap1 * pbp2 * LLLL +
           8.0 * pap1 * p1p3 * LRLR + 8.0 * pap1 * p1p3 * LLLL -
           8.0 * pap1 * p2p3 * LRLR - 8.0 * pap1 * p2p3 * LLLL +
           8.0 * pap2 * pbp1 * LRLR + 8. * pap2 * pbp1 * LLLL -
           8.0 * pap2 * pbp2 * LRLR - 8.0 * pap2 * pbp2 * LLLL -
           8.0 * pap2 * p1p3 * LRLR - 8.0 * pap2 * p1p3 * LLLL +
           8.0 * pap2 * p2p3 * LRLR + 8.0 * pap2 * p2p3 * LLLL +
           16.0 * pap3 * p1p2 * LRLR + 16.0 * pap3 * p1p2 * LLLL -
           8.0 * pap3 * m2 * m2 * LRLR - 8.0 * pap3 * m2 * m2 * LLLL -
           8.0 * pap3 * m1 * m1 * LRLR - 8.0 * pap3 * m1 * m1 * LLLL -
           16. * p1p3 * p2p3 * LRLR - 16.0 * p1p3 * p2p3 * LLLL +
           8.0 * p1p3 * p1p3 * LRLR + 8.0 * p1p3 * p1p3 * LLLL +
           8.0 * p2p3 * p2p3 * LRLR + 8.0 * p2p3 * p2p3 * LLLL)

      +
      ivu3 * ivs * ivs3v1 * ivs3v2 *
          (16.0 * pap1 * pap3 * pbp1 * LRLR + 16.0 * pap1 * pap3 * pbp1 * LLLL -
           16.0 * pap1 * pap3 * pbp2 * LRLR - 16.0 * pap1 * pap3 * pbp2 * LLLL -
           32.0 * pap1 * pap3 * p1p3 * LRLR - 32.0 * pap1 * pap3 * p1p3 * LLLL +
           32.0 * pap1 * pap3 * p2p3 * LRLR + 32.0 * pap1 * pap3 * p2p3 * LLLL -
           16.0 * pap2 * pap3 * pbp1 * LRLR - 16.0 * pap2 * pap3 * pbp1 * LLLL +
           16.0 * pap2 * pap3 * pbp2 * LRLR + 16. * pap2 * pap3 * pbp2 * LLLL +
           32.0 * pap2 * pap3 * p1p3 * LRLR + 32.0 * pap2 * pap3 * p1p3 * LLLL -
           32.0 * pap2 * pap3 * p2p3 * LRLR - 32.0 * pap2 * pap3 * p2p3 * LLLL -
           16.0 * pap3 * pbp1 * p1p3 * LRLR - 16.0 * pap3 * pbp1 * p1p3 * LLLL +
           16.0 * pap3 * pbp1 * p2p3 * LRLR + 16.0 * pap3 * pbp1 * p2p3 * LLLL +
           16.0 * pap3 * pbp2 * p1p3 * LRLR + 16.0 * pap3 * pbp2 * p1p3 * LLLL -
           16.0 * pap3 * pbp2 * p2p3 * LRLR - 16. * pap3 * pbp2 * p2p3 * LLLL -
           32.0 * pap3 * pap3 * p1p2 * LRLR - 32.0 * pap3 * pap3 * p1p2 * LLLL +
           16.0 * pap3 * pap3 * m2 * m2 * LRLR +
           16.0 * pap3 * pap3 * m2 * m2 * LLLL +
           16.0 * pap3 * pap3 * m1 * m1 * LRLR +
           16.0 * pap3 * pap3 * m1 * m1 * LLLL);

  return std::real(mss);
}
