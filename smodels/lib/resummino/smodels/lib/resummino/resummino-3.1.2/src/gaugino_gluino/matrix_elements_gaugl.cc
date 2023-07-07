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
// with two final state particles.

// NOTE: counterterms are implemented in MS (not MSbar scheme), because
// looptools
//       sets the universal finite part to zero per default.

#include "kinematics.h"
#include "npf.h"
#include "utils.h"
#include <complex>
#include <iostream>

using namespace std;

// mass counterterm for internal squark propagator correction
#define MC
// color factors
#define NC 3.0
#define cA 3.0
#define cF (4.0 / 3.0)
// 1/pi^2 term from i/16 pi^2
#define MPIsI (1.0 / pow2(M_PI))

// Born
/************************************************/
double FI::Mtt_GLGA() {
  return real(12.0 * cF * ivt1s1 * ivt1s2 * pap1 * pbp2 *
              (LLLL + LRLR + RLRL + RRRR));
}

double FI::Muu_GLGA() {
  return real(12.0 * cF * ivu2s1 * ivu2s2 * pap2 * pbp1 *
              (LLLL + LRLR + RLRL + RRRR));
}

double FI::Mtu_GLGA() {
  return real(-6.0 * cF * ivt1s1 * ivu2s2 *
              (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
               m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
               pap1 * pbp2 * RRRR));
}

// Quark Self Energy (squark-gluino loop)
/************************************************/
double FI::Mtt_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
                        (LR * (LLLL + RLRL) + RL * (LRLR + RRRR))) /
                       16.0);
  }

  return real((3.0 * b1 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
               (LR * (LLLL + RLRL) + RL * (LRLR + RRRR))) /
              8.0) +
         counterterm;
}

double FI::Muu_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
                        (LR * (LLLL + LRLR) + RL * (RLRL + RRRR))) /
                       16.0);
  }
  return real((3.0 * b1 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
               (LR * (LLLL + LRLR) + RL * (RLRL + RRRR))) /
              8.0) +
         counterterm;
}

double FI::Mtu_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm =
        real((3.0 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
              (LLLL * LR * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) -
               m1 * m2 * papb * (LRRL * RL + LR * RLLR) +
               (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRR)) /
             32.0);
  }

  return real((3.0 * b1 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
               (LLLL * LR * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) -
                m1 * m2 * papb * (LRRL * RL + LR * RLLR) +
                (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRR)) /
              16.0) +
         counterterm;
}

double FI::Mut_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm =
        real((-3.0 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
              (LR * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
               RL * (m1 * m2 * papb * RLLR +
                     (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
             32.0);
  }

  return real((-3.0 * b1 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
               (LR * (LRRL * m1 * m2 * papb +
                      LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
                RL * (m1 * m2 * papb * RLLR +
                      (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
              16.0) +
         counterterm;
}

// Antiquark self energy
/************************************************/
double FI::Mtt_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
                        ((LLLL + LRLR) * RL + LR * (RLRL + RRRR))) /
                       16.0);
  }

  return real((3.0 * b1 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
               ((LLLL + LRLR) * RL + LR * (RLRL + RRRR))) /
              8.0) +
         counterterm;
}

double FI::Muu_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
                        (RL * (LLLL + RLRL) + LR * (LRLR + RRRR))) /
                       16.0);
  }

  return real((3.0 * b1 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
               (RL * (LLLL + RLRL) + LR * (LRLR + RRRR))) /
              8.0) +
         counterterm;
}

double FI::Mtu_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm =
        real((-3.0 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
              (LRRL * m1 * m2 * papb * RL +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL +
               LR * (m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR))) /
             32.0);
  }

  return real((-3.0 * b1 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
               (LRRL * m1 * m2 * papb * RL +
                LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL +
                LR * (m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                      pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR))) /
              16.0) +
         counterterm;
}

double FI::Mut_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps) {
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm =
        real((-3.0 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
              (RL * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR) +
               LR * (LRRL * m1 * m2 * papb +
                     (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
             32.0);
  }

  return real((-3.0 * b1 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
               (RL * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                      m1 * m2 * papb * RLLR) +
                LR * (LRRL * m1 * m2 * papb +
                      (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
              16.0) +
         counterterm;
}

// Gluino Self Energy
/************************************************/

// Contribution 1
double FI::Mtt_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);

  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * cF * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                       32.0);
  }
  return real((3.0 * cF * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
               (b1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR) +
                2.0 * (db1 * mGLs * (LR + RL) * (LLLL + LRLR + RLRL + RRRR) -
                       2.0 * db0 * mGL * ml1 *
                           (LL * (LLLL + RLRL) + RR * (LRLR + RRRR))))) /
              16.0) +
         counterterm;
}

double FI::Muu_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);

  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                       32.0);
  }

  return real((3.0 * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
               (b1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR) +
                2.0 * (db1 * mGLs * (LR + RL) * (LLLL + LRLR + RLRL + RRRR) -
                       2.0 * db0 * mGL * ml1 *
                           (LL * (LLLL + RLRL) + RR * (LRLR + RRRR))))) /
              16.0) +
         counterterm;
}

double FI::Mtu_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);

  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((-3.0 * cF * ivt1s1 * ivu2s2 * MPIsI * (LR + RL) *
                        (LRRL * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                       64.0);
  }

  return real((-3.0 * cF * ivt1s1 * ivu2s2 * MPIsI *
               (b1 * (LR + RL) * (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                                  LLLL * pap2 * pbp1 + LLLL * pap1 * pbp2 +
                                  m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                                  pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
                2.0 * (db1 * mGLs * (LR + RL) *
                           (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                            LLLL * pap2 * pbp1 + LLLL * pap1 * pbp2 +
                            m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                            pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
                       2.0 * db0 * mGL * ml1 *
                           (LL * (LLLL * p1p2 * papb - LLLL * pap2 * pbp1 -
                                  LLLL * pap1 * pbp2 - m1 * m2 * papb * RLLR) -
                            RR * (LRRL * m1 * m2 * papb - p1p2 * papb * RRRR +
                                  pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR))))) /
              32.0) +
         counterterm;
}

double FI::Mut_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double mGLs = pow2(mGL);
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);

  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((-3.0 * cF * ivt1s2 * ivu2s1 * MPIsI * (LR + RL) *
                        (LRRL * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                       64.0);
  }

  return real((-3.0 * cF * ivt1s2 * ivu2s1 * MPIsI *
               (b1 * (LR + RL) * (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                                  LLLL * pap2 * pbp1 + LLLL * pap1 * pbp2 +
                                  m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                                  pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
                2.0 * (db1 * mGLs * (LR + RL) *
                           (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                            LLLL * pap2 * pbp1 + LLLL * pap1 * pbp2 +
                            m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                            pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
                       2.0 * db0 * mGL * ml1 *
                           (LL * (LLLL * p1p2 * papb - LLLL * pap2 * pbp1 -
                                  LLLL * pap1 * pbp2 - m1 * m2 * papb * RLLR) -
                            RR * (LRRL * m1 * m2 * papb - p1p2 * papb * RRRR +
                                  pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR))))) /
              32.0) +
         counterterm;
}

// Contribution 2
double FI::Mtt_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * cA * cF * ivt1s1 * ivt1s2 * MPIsI * pap1 * pbp2 *
                        (LLLL + LRLR + RLRL + RRRR)) /
                       8.0);
  }

  complex<double> me0 =
      (3.0 * cA * cF * ivt1s1 * ivt1s2 * (b1 + 2.0 * (2.0 * db0 + db1) * mGLs) *
       MPIsI * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (-3.0 * cA * cF * ivt1s1 * ivt1s2 * (b1 + 2.0 * (db0 + db1) * mGLs) *
       MPIsI * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  return real(me0 + me1) + counterterm;
}

double FI::Muu_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((3.0 * cA * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
                        (LLLL + LRLR + RLRL + RRRR)) /
                       8.0);
  }

  complex<double> me0 =
      (3.0 * cA * cF * ivu2s1 * ivu2s2 * (b1 + 2.0 * (2.0 * db0 + db1) * mGLs) *
       MPIsI * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (-3.0 * cA * cF * ivu2s1 * ivu2s2 * (b1 + 2.0 * (db0 + db1) * mGLs) *
       MPIsI * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  return real(me0 + me1) + counterterm;
}

double FI::Mtu_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((-3.0 * cA * cF * ivt1s1 * ivu2s2 * MPIsI *
                        (LRRL * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                       16.0);
  }

  complex<double> me0 = (-3.0 * cA * cF * ivt1s1 * ivu2s2 *
                         (b1 + 2.0 * (2.0 * db0 + db1) * mGLs) * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        8.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * cA * cF * ivt1s1 * ivu2s2 * (b1 + 2.0 * (db0 + db1) * mGLs) *
       MPIsI * (LRRL * m1 * m2 * papb +
                LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
      8.0;

  return real(me0 + me1) + counterterm;
}

double FI::Mut_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                         int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double mGLs = pow2(mGL);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double counterterm = 0.0;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  if (ieps == 1) {
    counterterm = real((-3.0 * cA * cF * ivt1s2 * ivu2s1 * MPIsI *
                        (LRRL * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                       16.0);
  }

  complex<double> me0 = (-3.0 * cA * cF * ivt1s2 * ivu2s1 *
                         (b1 + 2.0 * (2.0 * db0 + db1) * mGLs) * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        8.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * cA * cF * ivt1s2 * ivu2s1 * (b1 + 2.0 * (db0 + db1) * mGLs) *
       MPIsI * (LRRL * m1 * m2 * papb +
                LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
      8.0;

  return real(me0 + me1) + counterterm;
}

// Squark Propagator Corrections and Counterterms
/************************************************/
// Contribution 1
double FI::Mtt_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, int i, int j, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double fieldcounterterm = 0.0;
  complex<double> me0_mass = (0.0, 0.0);
  complex<double> me1_mass = (0.0, 0.0);

  if (ieps == 1 && i == j) {
    fieldcounterterm = real((3.0 * pow2(cF) * pow2(ivt1s1) * ivt1s2 * MPIsI *
                             (m1s - mp1s - 2.0 * pap1) * pap1 * pbp2 *
                             (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                            4.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = (3.0 * pow2(cF) * ivt1s1 * ivt1s2 * ivt2s1 * MPIsI *
                         (4.0 * b00 + (b1 + b11) * (m1s - 2.0 * pap1)) * pap1 *
                         pbp2 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  if (i == j) {
    SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me0_mass = (-3.0 * pow2(cF) * pow2(ivt1s1) * ivt1s2 *
                (4.0 * b00 + (b1 + b11) * mp1s) * MPIsI * pap1 * pbp2 *
                (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
               2.0;
  }
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -3.0 * b00 * pow2(cF) * ivt1s1 * ivt1s2 * ivt2s1 *
                        MPIsI * pap1 * pbp2 * (LR + RL) *
                        (LLLL + LRLR + RLRL + RRRR);

  if (i == j) {

    SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me1_mass = 3.0 * b00 * pow2(cF) * pow2(ivt1s1) * ivt1s2 * MPIsI * pap1 *
               pbp2 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);
  }
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me1 + me1_mass + me0 + me0_mass) + fieldcounterterm;
}

double FI::Muu_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, int i, int j, double mSQs) {

  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double fieldcounterterm = 0.0;
  complex<double> me0_mass = (0.0, 0.0);
  complex<double> me1_mass = (0.0, 0.0);

  if (ieps == 1 && i == j) {
    fieldcounterterm = real((3.0 * pow2(cF) * pow2(ivu1s1) * ivu2s2 * MPIsI *
                             (m2s - mp1s - 2.0 * pap2) * pap2 * pbp1 *
                             (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                            4.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = (3.0 * pow2(cF) * ivu1s1 * ivu2s1 * ivu2s2 * MPIsI *
                         (4.0 * b00 + (b1 + b11) * (m2s - 2.0 * pap2)) * pap2 *
                         pbp1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  if (i == j) {
    SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me0_mass = (-3.0 * pow2(cF) * pow2(ivu1s1) * ivu2s2 *
                (4.0 * b00 + (b1 + b11) * mp1s) * MPIsI * pap2 * pbp1 *
                (LR + RL) * (LLLL + LRLR + RLRL + RRRR)) /
               2.0;
  }
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -3.0 * b00 * pow2(cF) * ivu1s1 * ivu2s1 * ivu2s2 *
                        MPIsI * pap2 * pbp1 * (LR + RL) *
                        (LLLL + LRLR + RLRL + RRRR);
  if (i == j) {

    SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me1_mass = 3.0 * b00 * pow2(cF) * pow2(ivu1s1) * ivu2s2 * MPIsI * pap2 *
               pbp1 * (LR + RL) * (LLLL + LRLR + RLRL + RRRR);
  }
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me1 + me1_mass + me0 + me0_mass) + fieldcounterterm;
}

double FI::Mtu_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, int i, int j, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double fieldcounterterm = 0.0;
  complex<double> me0_mass = (0.0, 0.0);
  complex<double> me1_mass = (0.0, 0.0);

  if (ieps == 1 && i == j) {
    fieldcounterterm =
        real((-3.0 * pow2(cF) * pow2(ivt1s1) * ivu2s2 * MPIsI *
              (m1s - mp1s - 2.0 * pap1) * (LR + RL) *
              (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
               m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
               pap1 * pbp2 * RRRR)) /
             8.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s1 * ivt2s1 * ivu2s2 * MPIsI *
       (4.0 * b00 + (b1 + b11) * (m1s - 2.0 * pap1)) * (LR + RL) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;

  if (i == j) {
    SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me0_mass = (3.0 * pow2(cF) * pow2(ivt1s1) * ivu2s2 *
                (4.0 * b00 + (b1 + b11) * mp1s) * MPIsI * (LR + RL) *
                (LRRL * m1 * m2 * papb +
                 LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                 m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                 pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
               4.0;
  }
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * b00 * pow2(cF) * ivt1s1 * ivt2s1 * ivu2s2 * MPIsI * (LR + RL) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      2.0;
  if (i == j) {

    SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me1_mass =
        (-3.0 * b00 * pow2(cF) * pow2(ivt1s1) * ivu2s2 * MPIsI * (LR + RL) *
         (LRRL * m1 * m2 * papb +
          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
          pap1 * pbp2 * RRRR)) /
        2.0;
  }
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me1 + me1_mass + me0 + me0_mass) + fieldcounterterm;
}

double FI::Mut_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, int i, int j, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double fieldcounterterm = 0.0;
  complex<double> me0_mass = (0.0, 0.0);
  complex<double> me1_mass = (0.0, 0.0);

  if (ieps == 1 && i == j) {
    fieldcounterterm =
        real((-3.0 * pow2(cF) * ivt1s2 * pow2(ivu1s1) * MPIsI *
              (m2s - mp1s - 2.0 * pap2) * (LR + RL) *
              (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
               m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
               pap1 * pbp2 * RRRR)) /
             8.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s2 * ivu1s1 * ivu2s1 * MPIsI *
       (4.0 * b00 + (b1 + b11) * (m2s - 2.0 * pap2)) * (LR + RL) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;

  if (i == j) {
    SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me0_mass = (3.0 * pow2(cF) * ivt1s2 * pow2(ivu1s1) *
                (4.0 * b00 + (b1 + b11) * mp1s) * MPIsI * (LR + RL) *
                (LRRL * m1 * m2 * papb +
                 LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                 m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                 pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
               4.0;
  }
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * b00 * pow2(cF) * ivt1s2 * ivu1s1 * ivu2s1 * MPIsI * (LR + RL) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      2.0;

  if (i == j) {

    SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
         &db11);

    me1_mass =
        (-3.0 * b00 * pow2(cF) * ivt1s2 * pow2(ivu1s1) * MPIsI * (LR + RL) *
         (LRRL * m1 * m2 * papb +
          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
          pap1 * pbp2 * RRRR)) /
        2.0;
  }
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me1 + me1_mass + me0 + me0_mass) + fieldcounterterm;
}

// Contribution 2
double FI::Mtt_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double fieldcounterterm = 0.0;

  if (ieps == 1) {
    fieldcounterterm = real((-3.0 * pow2(cF) * pow2(ivt1s1) * ivt1s2 * MPIsI *
                             (m1s - mp1s - 2.0 * pap1) * pap1 * pbp2 *
                             (LLLL + LRLR + RLRL + RRRR)) /
                            2.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (3.0 * pow2(cF) * ivt1s1 * ivt1s2 * ivt2s1 * MPIsI *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * (m1s - 2.0 * pap1)) * pap1 *
       pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0_mass =
      (-3.0 * pow2(cF) * pow2(ivt1s1) * ivt1s2 *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * mp1s) * MPIsI * pap1 * pbp2 *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = (-3.0 * b00 * pow2(cF) * ivt1s1 * ivt1s2 * ivt2s1 *
                         MPIsI * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1_mass =
      (3.0 * b00 * pow2(cF) * pow2(ivt1s1) * ivt1s2 * MPIsI * pap1 * pbp2 *
       (LLLL + LRLR + RLRL + RRRR)) /
      2.0;
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me0 + me0_mass + me1 + me1_mass) + fieldcounterterm;
}

double FI::Muu_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);

  double fieldcounterterm = 0.0;

  if (ieps == 1) {
    fieldcounterterm = real((-3.0 * pow2(cF) * pow2(ivu1s1) * ivu2s2 * MPIsI *
                             (m2s - mp1s - 2.0 * pap2) * pap2 * pbp1 *
                             (LLLL + LRLR + RLRL + RRRR)) /
                            2.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (3.0 * pow2(cF) * ivu1s1 * ivu2s1 * ivu2s2 * MPIsI *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * (m2s - 2.0 * pap2)) * pap2 *
       pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0_mass =
      (-3.0 * pow2(cF) * pow2(ivu1s1) * ivu2s2 *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * mp1s) * MPIsI * pap2 * pbp1 *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = (-3.0 * b00 * pow2(cF) * ivu1s1 * ivu2s1 * ivu2s2 *
                         MPIsI * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1_mass =
      (3.0 * b00 * pow2(cF) * pow2(ivu1s1) * ivu2s2 * MPIsI * pap2 * pbp1 *
       (LLLL + LRLR + RLRL + RRRR)) /
      2.0;
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me0 + me0_mass + me1 + me1_mass) + fieldcounterterm;
}

double FI::Mtu_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);

  double fieldcounterterm = 0.0;

  if (ieps == 1) {
    fieldcounterterm =
        real(3.0 * pow2(cF) * pow2(ivt1s1) * ivu2s2 * MPIsI *
             (m1s - mp1s - 2.0 * pap1) *
             (LRRL * m1 * m2 * papb +
              LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
              pap1 * pbp2 * RRRR)) /
        4.0;
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s1 * ivt2s1 * ivu2s2 * MPIsI *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * (m1s - 2.0 * pap1)) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0_mass =
      (3.0 * pow2(cF) * pow2(ivt1s1) * ivu2s2 *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * mp1s) * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * b00 * pow2(cF) * ivt1s1 * ivt2s1 * ivu2s2 * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;

  SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1_mass =
      (-3.0 * b00 * pow2(cF) * pow2(ivt1s1) * ivu2s2 * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me0 + me0_mass + me1 + me1_mass) + fieldcounterterm;
}

double FI::Mut_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                         double mp2s, int ieps, double mSQs) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);
  complex<double> b0, b1, b00, b11, db0, db1, db00, db11;
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);

  double fieldcounterterm = 0.0;

  if (ieps == 1) {
    fieldcounterterm =
        real((3.0 * pow2(cF) * ivt1s2 * pow2(ivu1s1) * MPIsI *
              (m2s - mp1s - 2.0 * pap2) *
              (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
               m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
               pap1 * pbp2 * RRRR)) /
             4.0);
  }

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s2 * ivu1s1 * ivu2s1 * MPIsI *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * (m2s - 2.0 * pap2)) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetB(mSQs, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0_mass =
      (3.0 * pow2(cF) * ivt1s2 * pow2(ivu1s1) *
       (4.0 * b00 + (4.0 * (b0 + b1) + b11) * mp1s) * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      (3.0 * b00 * pow2(cF) * ivt1s2 * ivu1s1 * ivu2s1 * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;

  SetB(mSQs, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1_mass =
      (-3.0 * b00 * pow2(cF) * ivt1s2 * pow2(ivu1s1) * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      4.0;
#ifndef MC
  me0_mass = 0.0;
  me1_mass = 0.0;
#endif
  return real(me0 + me0_mass + me1 + me1_mass) + fieldcounterterm;
}

// Triangle 1
/************************************************/
// Correction at leg a (pa)
double FI::Mtt_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * ivt2s1 * MPIsI * pap1 * pbp2 *
       (-4.0 * c00 + 2.0 * (2.0 * c1 + c12) * (pap2 - papb) -
        (2.0 * c2 + c22) * (m2s - 2.0 * (pap2 - papb + pbp2))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * (cA - 2.0 * cF) * cF * ivt1s2 * ivt2s1 *
                         MPIsI * pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Muu_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
       (4.0 * c00 - 2.0 * (2.0 * c1 + c12) * (pap1 - papb) +
        (2.0 * c2 + c22) * (m1s - 2.0 * (pap1 - papb + pbp1))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (-3.0 * c00 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI *
                         pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  return real(me0 + me1);
}

double FI::Mtu_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivt2s1 * ivu2s2 * MPIsI *
       (4.0 * c00 - 2.0 * (2.0 * c1 + c12) * (pap2 - papb) +
        (2.0 * c2 + c22) * (m2s - 2.0 * (pap2 - papb + pbp2))) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      (-3.0 * c00 * (cA - 2.0 * cF) * cF * ivt2s1 * ivu2s2 * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  return real(me0 + me1);
}

double FI::Mut_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
       (4.0 * c00 - 2.0 * (2.0 * c1 + c12) * (pap1 - papb) +
        (2.0 * c2 + c22) * (m1s - 2.0 * (pap1 - papb + pbp1))) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        4.0;

  return real(me0 + me1);
}

// Correction at leg b (pb)
double FI::Mtt_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * pap1 *
       (4.0 * c00 + 2.0 * (2.0 * c1 + c12) * (papb - pbp1) +
        (2.0 * c2 + c22) * (m1s - 2.0 * (pap1 - papb + pbp1))) *
       pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (-3.0 * c00 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI *
                         pap1 * pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
                        2.0;

  return real(me0 + me1);
}

double FI::Muu_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 * pbp1 *
       (-4.0 * c00 - 2.0 * (2.0 * c1 + c12) * (papb - pbp2) -
        (2.0 * c2 + c22) * (m2s - 2.0 * (pap2 - papb + pbp2))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * (cA - 2.0 * cF) * cF * ivu2s1 * ivu2s2 *
                         MPIsI * pap2 * pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Mtu_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
       (4.0 * c00 + 2.0 * (2.0 * c1 + c12) * (papb - pbp1) +
        (2.0 * c2 + c22) * (m1s - 2.0 * (pap1 - papb + pbp1))) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Mut_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * ivu2s1 * MPIsI *
       (4.0 * c00 + 2.0 * (2.0 * c1 + c12) * (papb - pbp2) +
        (2.0 * c2 + c22) * (m2s - 2.0 * (pap2 - papb + pbp2))) *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      (-3.0 * c00 * (cA - 2.0 * cF) * cF * ivt1s2 * ivu2s1 * MPIsI *
       (LRRL * m1 * m2 * papb +
        LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
        m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
        pap1 * pbp2 * RRRR)) /
      8.0;

  return real(me0 + me1);
}

// Triangle 2
/************************************************/
// Correction at leg a
double FI::Mtt_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      ((cA - 2.0 * cF) * cF * ivt1s2 * ivt2s1 * MPIsI * NC * pbp2 *
       (c22 * LLLL * LR * m2s * pap1 - 2.0 * c12 * LLLL * LR * pap1 * pap2 -
        2.0 * c22 * LLLL * LR * pap1 * pap2 +
        2.0 * c12 * LLLL * LR * pap1 * papb +
        2.0 * c22 * LLLL * LR * pap1 * papb +
        c0 * LL * LRLL * m1 * ml2 * (-pap2 + papb) -
        2.0 * c22 * LLLL * LR * pap1 * pbp2 + c22 * LRLR * m2s * pap1 * RL -
        2.0 * c12 * LRLR * pap1 * pap2 * RL -
        2.0 * c22 * LRLR * pap1 * pap2 * RL +
        2.0 * c12 * LRLR * pap1 * papb * RL +
        2.0 * c22 * LRLR * pap1 * papb * RL -
        2.0 * c22 * LRLR * pap1 * pbp2 * RL + c22 * LR * m2s * pap1 * RLRL -
        2.0 * c12 * LR * pap1 * pap2 * RLRL -
        2.0 * c22 * LR * pap1 * pap2 * RLRL +
        2.0 * c12 * LR * pap1 * papb * RLRL +
        2.0 * c22 * LR * pap1 * papb * RLRL -
        2.0 * c22 * LR * pap1 * pbp2 * RLRL - c0 * LLLR * m1 * ml2 * pap2 * RR +
        c0 * LLLR * m1 * ml2 * papb * RR - c0 * m1 * ml2 * pap2 * RLRR * RR +
        c0 * m1 * ml2 * papb * RLRR * RR - c0 * LL * m1 * ml2 * pap2 * RRRL +
        c0 * LL * m1 * ml2 * papb * RRRL +
        pap1 * (2.0 * c12 * (-pap2 + papb) +
                c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            RL * RRRR +
        c2 * (LLLL * LR * pap1 * (m2s - 2.0 * (pap2 - papb + pbp2)) +
              LRLR * m2s * pap1 * RL - 2.0 * LRLR * pap1 * pap2 * RL +
              2.0 * LRLR * pap1 * papb * RL - 2.0 * LRLR * pap1 * pbp2 * RL +
              LR * m2s * pap1 * RLRL - 2.0 * LR * pap1 * pap2 * RLRL +
              2.0 * LR * pap1 * papb * RLRL - 2.0 * LR * pap1 * pbp2 * RLRL -
              LLLR * m1 * ml2 * pap2 * RR + LLLR * m1 * ml2 * papb * RR -
              m1 * ml2 * pap2 * RLRR * RR + m1 * ml2 * papb * RLRR * RR -
              LL * m1 * ml2 * (pap2 - papb) * (LRLL + RRRL) +
              pap1 * (m2s - 2.0 * (pap2 - papb + pbp2)) * RL * RRRR) +
        4.0 * c00 * pap1 * (LR * (LLLL + RLRL) + RL * (LRLR + RRRR)))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * (cA - 2.0 * cF) * cF * ivt1s2 * ivt2s1 * MPIsI * NC * pap1 *
        pbp2 * (LR * (LLLL + RLRL) + RL * (LRLR + RRRR))) /
      4.0;

  return real(me0 + me1);
}

double FI::Muu_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * NC * pbp1 *
       (-(c22 * LLLL * LR * m1s * pap2) - c22 * LR * LRLR * m1s * pap2 +
        2.0 * c12 * LLLL * LR * pap1 * pap2 +
        2.0 * c22 * LLLL * LR * pap1 * pap2 +
        2.0 * c12 * LR * LRLR * pap1 * pap2 +
        2.0 * c22 * LR * LRLR * pap1 * pap2 -
        2.0 * c12 * LLLL * LR * pap2 * papb -
        2.0 * c22 * LLLL * LR * pap2 * papb -
        2.0 * c12 * LR * LRLR * pap2 * papb -
        2.0 * c22 * LR * LRLR * pap2 * papb +
        2.0 * c22 * LLLL * LR * pap2 * pbp1 +
        2.0 * c22 * LR * LRLR * pap2 * pbp1 + c0 * LL * m2 * ml2 * pap1 * RLLL -
        c0 * LL * m2 * ml2 * papb * RLLL - c22 * m1s * pap2 * RL * RLRL +
        2.0 * c12 * pap1 * pap2 * RL * RLRL +
        2.0 * c22 * pap1 * pap2 * RL * RLRL -
        2.0 * c12 * pap2 * papb * RL * RLRL -
        2.0 * c22 * pap2 * papb * RL * RLRL +
        2.0 * c22 * pap2 * pbp1 * RL * RLRL + c0 * LLRL * m2 * ml2 * pap1 * RR +
        c0 * LRRR * m2 * ml2 * pap1 * RR - c0 * LLRL * m2 * ml2 * papb * RR -
        c0 * LRRR * m2 * ml2 * papb * RR + c0 * LL * m2 * ml2 * pap1 * RRLR -
        c0 * LL * m2 * ml2 * papb * RRLR +
        pap2 * (2.0 * c12 * (pap1 - papb) -
                c22 * (m1s - 2.0 * (pap1 - papb + pbp1))) *
            RL * RRRR -
        c2 * (LLLL * LR * pap2 * (m1s - 2.0 * (pap1 - papb + pbp1)) +
              LR * LRLR * pap2 * (m1s - 2.0 * (pap1 - papb + pbp1)) -
              LL * m2 * ml2 * pap1 * RLLL + LL * m2 * ml2 * papb * RLLL +
              m1s * pap2 * RL * RLRL - 2.0 * pap1 * pap2 * RL * RLRL +
              2.0 * pap2 * papb * RL * RLRL - 2.0 * pap2 * pbp1 * RL * RLRL -
              LLRL * m2 * ml2 * pap1 * RR - LRRR * m2 * ml2 * pap1 * RR +
              LLRL * m2 * ml2 * papb * RR + LRRR * m2 * ml2 * papb * RR -
              LL * m2 * ml2 * pap1 * RRLR + LL * m2 * ml2 * papb * RRLR +
              pap2 * (m1s - 2.0 * (pap1 - papb + pbp1)) * RL * RRRR) -
        4.0 * c00 * pap2 * (LR * (LLLL + LRLR) + RL * (RLRL + RRRR)))) /
      4.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (c00 * pow2(cF) * ivu2s1 * ivu2s2 * MPIsI * NC * pap2 *
                         pbp1 * (LR * (LLLL + LRLR) + RL * (RLRL + RRRR))) /
                        2.0;
  return real(me0 + me1);
}

double FI::Mtu_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -((cA - 2.0 * cF) * cF * ivt2s1 * ivu2s2 * MPIsI * NC *
        (c0 * ml2 * (LL * LRLL * m1 * (m2s * papb - 2.0 * pap2 * pbp2) +
                     (LLRL * m2 * (-(p1p2 * papb) - pap2 * pbp1 +
                                   2.0 * papb * pbp1 + pap1 * pbp2) +
                      m1 * (m2s * papb - 2.0 * pap2 * pbp2) * RLRR) *
                         RR +
                     LL * m2 * (-(p1p2 * papb) - pap2 * pbp1 +
                                2.0 * papb * pbp1 + pap1 * pbp2) *
                         RRLR) +
         (4.0 * c00 + 2.0 * c12 * (-pap2 + papb) +
          c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
             (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              LRRL * m1 * m2 * papb * RL + LR * m1 * m2 * papb * RLLR -
              p1p2 * papb * RL * RRRR + pap2 * pbp1 * RL * RRRR +
              pap1 * pbp2 * RL * RRRR) +
         c2 * (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                   (m2s - 2.0 * (pap2 - papb + pbp2)) +
               LRRL * m1 * m2 * m2s * papb * RL -
               2.0 * LRRL * m1 * m2 * pap2 * papb * RL +
               2.0 * LRRL * m1 * m2 * pow2(papb) * RL -
               2.0 * LRRL * m1 * m2 * papb * pbp2 * RL +
               LR * m1 * m2 * m2s * papb * RLLR -
               2.0 * LR * m1 * m2 * pap2 * papb * RLLR +
               2.0 * LR * m1 * m2 * pow2(papb) * RLLR -
               2.0 * LR * m1 * m2 * papb * pbp2 * RLLR -
               LLRL * m2 * ml2 * p1p2 * papb * RR -
               LLRL * m2 * ml2 * pap2 * pbp1 * RR +
               2.0 * LLRL * m2 * ml2 * papb * pbp1 * RR +
               LLRL * m2 * ml2 * pap1 * pbp2 * RR +
               m1 * m2s * ml2 * papb * RLRR * RR -
               2.0 * m1 * ml2 * pap2 * pbp2 * RLRR * RR +
               LL * ml2 * (LRLL * m1 * (m2s * papb - 2.0 * pap2 * pbp2) +
                           m2 * (-(p1p2 * papb) - pap2 * pbp1 +
                                 2.0 * papb * pbp1 + pap1 * pbp2) *
                               RRLR) +
               (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                   (m2s - 2.0 * (pap2 - papb + pbp2)) * RL * RRRR))) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * (cA - 2.0 * cF) * cF * ivt2s1 * ivu2s2 * MPIsI * NC *
        (LLLL * LR * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) -
         m1 * m2 * papb * (LRRL * RL + LR * RLLR) +
         (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRR)) /
      8.0;

  return real(me0 + me1);
}

double FI::Mut_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (pow2(cF) * ivt1s2 * ivu2s1 * MPIsI * NC *
       (c22 * LR * LRRL * m1 * m1s * m2 * papb -
        c22 * LLLL * LR * m1s * p1p2 * papb -
        2.0 * c12 * LR * LRRL * m1 * m2 * pap1 * papb -
        2.0 * c22 * LR * LRRL * m1 * m2 * pap1 * papb +
        2.0 * c12 * LLLL * LR * p1p2 * pap1 * papb +
        2.0 * c22 * LLLL * LR * p1p2 * pap1 * papb +
        2.0 * c12 * LR * LRRL * m1 * m2 * pow2(papb) +
        2.0 * c22 * LR * LRRL * m1 * m2 * pow2(papb) -
        2.0 * c12 * LLLL * LR * p1p2 * pow2(papb) -
        2.0 * c22 * LLLL * LR * p1p2 * pow2(papb) +
        c22 * LLLL * LR * m1s * pap2 * pbp1 -
        2.0 * c12 * LLLL * LR * pap1 * pap2 * pbp1 -
        2.0 * c22 * LLLL * LR * pap1 * pap2 * pbp1 -
        2.0 * c22 * LR * LRRL * m1 * m2 * papb * pbp1 +
        2.0 * c22 * LLLL * LR * p1p2 * papb * pbp1 +
        2.0 * c12 * LLLL * LR * pap2 * papb * pbp1 +
        2.0 * c22 * LLLL * LR * pap2 * papb * pbp1 -
        2.0 * c22 * LLLL * LR * pap2 * pow2(pbp1) +
        c22 * LLLL * LR * m1s * pap1 * pbp2 -
        2.0 * c12 * LLLL * LR * pow2(pap1) * pbp2 -
        2.0 * c22 * LLLL * LR * pow2(pap1) * pbp2 +
        2.0 * c12 * LLLL * LR * pap1 * papb * pbp2 +
        2.0 * c22 * LLLL * LR * pap1 * papb * pbp2 -
        2.0 * c22 * LLLL * LR * pap1 * pbp1 * pbp2 +
        c0 * LL * m1s * m2 * ml2 * papb * RLLL -
        2.0 * c0 * LL * m2 * ml2 * pap1 * pbp1 * RLLL +
        c22 * m1 * m1s * m2 * papb * RL * RLLR -
        2.0 * c12 * m1 * m2 * pap1 * papb * RL * RLLR -
        2.0 * c22 * m1 * m2 * pap1 * papb * RL * RLLR +
        2.0 * c12 * m1 * m2 * pow2(papb) * RL * RLLR +
        2.0 * c22 * m1 * m2 * pow2(papb) * RL * RLLR -
        2.0 * c22 * m1 * m2 * papb * pbp1 * RL * RLLR +
        c0 * LRRR * m1s * m2 * ml2 * papb * RR -
        c0 * LLLR * m1 * ml2 * p1p2 * papb * RR -
        2.0 * c0 * LRRR * m2 * ml2 * pap1 * pbp1 * RR +
        c0 * LLLR * m1 * ml2 * pap2 * pbp1 * RR -
        c0 * LLLR * m1 * ml2 * pap1 * pbp2 * RR +
        2.0 * c0 * LLLR * m1 * ml2 * papb * pbp2 * RR -
        c0 * LL * m1 * ml2 * p1p2 * papb * RRRL +
        c0 * LL * m1 * ml2 * pap2 * pbp1 * RRRL -
        c0 * LL * m1 * ml2 * pap1 * pbp2 * RRRL +
        2.0 * c0 * LL * m1 * ml2 * papb * pbp2 * RRRL +
        (2.0 * c12 * (-pap1 + papb) +
         c22 * (m1s - 2.0 * (pap1 - papb + pbp1))) *
            (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL * RRRR +
        c2 * (LR * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                  (LRRL * m1 * m2 * papb +
                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
              m1 * m1s * m2 * papb * RL * RLLR -
              2.0 * m1 * m2 * pap1 * papb * RL * RLLR +
              2.0 * m1 * m2 * pow2(papb) * RL * RLLR -
              2.0 * m1 * m2 * papb * pbp1 * RL * RLLR +
              LRRR * m1s * m2 * ml2 * papb * RR -
              LLLR * m1 * ml2 * p1p2 * papb * RR -
              2.0 * LRRR * m2 * ml2 * pap1 * pbp1 * RR +
              LLLR * m1 * ml2 * pap2 * pbp1 * RR -
              LLLR * m1 * ml2 * pap1 * pbp2 * RR +
              2.0 * LLLR * m1 * ml2 * papb * pbp2 * RR +
              LL * ml2 * (m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLL +
                          m1 * (-(p1p2 * papb) + pap2 * pbp1 - pap1 * pbp2 +
                                2.0 * papb * pbp2) *
                              RRRL) +
              (m1s - 2.0 * (pap1 - papb + pbp1)) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
        4.0 * c00 *
            (LR * (LRRL * m1 * m2 * papb +
                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
             RL * (m1 * m2 * papb * RLLR +
                   (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR)))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * pow2(cF) * ivt1s2 * ivu2s1 * MPIsI * NC *
        (LR * (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
         RL * (m1 * m2 * papb * RLLR +
               (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
      4.0;

  return real(me0 + me1);
}

// Correction at leg b
double FI::Mtt_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * NC * pap1 *
       (-(c22 * LLLL * LR * m1s * pbp2) - c22 * LR * LRLR * m1s * pbp2 +
        2.0 * c22 * LLLL * LR * pap1 * pbp2 +
        2.0 * c22 * LR * LRLR * pap1 * pbp2 -
        2.0 * c12 * LLLL * LR * papb * pbp2 -
        2.0 * c22 * LLLL * LR * papb * pbp2 -
        2.0 * c12 * LR * LRLR * papb * pbp2 -
        2.0 * c22 * LR * LRLR * papb * pbp2 +
        2.0 * c12 * LLLL * LR * pbp1 * pbp2 +
        2.0 * c22 * LLLL * LR * pbp1 * pbp2 +
        2.0 * c12 * LR * LRLR * pbp1 * pbp2 +
        2.0 * c22 * LR * LRLR * pbp1 * pbp2 - c0 * LL * m2 * ml2 * papb * RLLL +
        c0 * LL * m2 * ml2 * pbp1 * RLLL - c22 * m1s * pbp2 * RL * RLRL +
        2.0 * c22 * pap1 * pbp2 * RL * RLRL -
        2.0 * c12 * papb * pbp2 * RL * RLRL -
        2.0 * c22 * papb * pbp2 * RL * RLRL +
        2.0 * c12 * pbp1 * pbp2 * RL * RLRL +
        2.0 * c22 * pbp1 * pbp2 * RL * RLRL - c0 * LLRL * m2 * ml2 * papb * RR -
        c0 * LRRR * m2 * ml2 * papb * RR + c0 * LLRL * m2 * ml2 * pbp1 * RR +
        c0 * LRRR * m2 * ml2 * pbp1 * RR - c0 * LL * m2 * ml2 * papb * RRLR +
        c0 * LL * m2 * ml2 * pbp1 * RRLR +
        (2.0 * c12 * (-papb + pbp1) -
         c22 * (m1s - 2.0 * (pap1 - papb + pbp1))) *
            pbp2 * RL * RRRR -
        c2 * (LLLL * LR * (m1s - 2.0 * (pap1 - papb + pbp1)) * pbp2 +
              LR * LRLR * (m1s - 2.0 * (pap1 - papb + pbp1)) * pbp2 +
              (m1s - 2.0 * pap1) * pbp2 * RL * RLRL +
              (papb - pbp1) *
                  (2.0 * pbp2 * RL * RLRL +
                   m2 * ml2 * ((LLRL + LRRR) * RR + LL * (RLLL + RRLR))) +
              (m1s - 2.0 * (pap1 - papb + pbp1)) * pbp2 * RL * RRRR) -
        4.0 * c00 * pbp2 * (LR * (LLLL + LRLR) + RL * (RLRL + RRRR)))) /
      4.0;

  // gives only non-zero if sbottom mixing is included; but we don't!
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (c00 * pow2(cF) * ivt1s1 * ivt1s2 * MPIsI * NC * pap1 *
                         pbp2 * (LR * (LLLL + LRLR) + RL * (RLRL + RRRR))) /
                        2.0;

  return real(me0 + me1);
}

double FI::Muu_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      ((cA - 2.0 * cF) * cF * ivu2s1 * ivu2s2 * MPIsI * NC * pap2 *
       (c0 * m1 * ml2 * (papb - pbp2) *
            ((LLLR + RLRR) * RR + LL * (LRLL + RRRL)) +
        c2 * (LLLL * LR * pbp1 * (m2s - 2.0 * (pap2 - papb + pbp2)) +
              LRLR * m2s * pbp1 * RL - 2.0 * LRLR * pap2 * pbp1 * RL +
              2.0 * LRLR * papb * pbp1 * RL - 2.0 * LRLR * pbp1 * pbp2 * RL +
              LR * m2s * pbp1 * RLRL - 2.0 * LR * pap2 * pbp1 * RLRL +
              2.0 * LR * papb * pbp1 * RLRL - 2.0 * LR * pbp1 * pbp2 * RLRL +
              LLLR * m1 * ml2 * papb * RR - LLLR * m1 * ml2 * pbp2 * RR +
              m1 * ml2 * papb * RLRR * RR - m1 * ml2 * pbp2 * RLRR * RR +
              LL * m1 * ml2 * (papb - pbp2) * (LRLL + RRRL) +
              pbp1 * (m2s - 2.0 * (pap2 - papb + pbp2)) * RL * RRRR) +
        pbp1 * (4.0 * c00 + 2.0 * c12 * (papb - pbp2) +
                c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            (LR * (LLLL + RLRL) + RL * (LRLR + RRRR)))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * (cA - 2.0 * cF) * cF * ivu2s1 * ivu2s2 * MPIsI * NC * pap2 *
        pbp1 * (LR * (LLLL + RLRL) + RL * (LRLR + RRRR))) /
      4.0;

  return real(me0 + me1);
}

double FI::Mtu_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (pow2(cF) * ivt1s1 * ivu2s2 * MPIsI * NC *
       (c22 * LR * LRRL * m1 * m1s * m2 * papb -
        c22 * LLLL * LR * m1s * p1p2 * papb -
        2.0 * c22 * LR * LRRL * m1 * m2 * pap1 * papb +
        2.0 * c22 * LLLL * LR * p1p2 * pap1 * papb +
        2.0 * c12 * LR * LRRL * m1 * m2 * pow2(papb) +
        2.0 * c22 * LR * LRRL * m1 * m2 * pow2(papb) -
        2.0 * c12 * LLLL * LR * p1p2 * pow2(papb) -
        2.0 * c22 * LLLL * LR * p1p2 * pow2(papb) +
        c22 * LLLL * LR * m1s * pap2 * pbp1 -
        2.0 * c22 * LLLL * LR * pap1 * pap2 * pbp1 -
        2.0 * c12 * LR * LRRL * m1 * m2 * papb * pbp1 -
        2.0 * c22 * LR * LRRL * m1 * m2 * papb * pbp1 +
        2.0 * c12 * LLLL * LR * p1p2 * papb * pbp1 +
        2.0 * c22 * LLLL * LR * p1p2 * papb * pbp1 +
        2.0 * c12 * LLLL * LR * pap2 * papb * pbp1 +
        2.0 * c22 * LLLL * LR * pap2 * papb * pbp1 -
        2.0 * c12 * LLLL * LR * pap2 * pow2(pbp1) -
        2.0 * c22 * LLLL * LR * pap2 * pow2(pbp1) +
        c22 * LLLL * LR * m1s * pap1 * pbp2 -
        2.0 * c22 * LLLL * LR * pow2(pap1) * pbp2 +
        2.0 * c12 * LLLL * LR * pap1 * papb * pbp2 +
        2.0 * c22 * LLLL * LR * pap1 * papb * pbp2 -
        2.0 * c12 * LLLL * LR * pap1 * pbp1 * pbp2 -
        2.0 * c22 * LLLL * LR * pap1 * pbp1 * pbp2 +
        c0 * LL * m1s * m2 * ml2 * papb * RLLL -
        2.0 * c0 * LL * m2 * ml2 * pap1 * pbp1 * RLLL +
        c22 * m1 * m1s * m2 * papb * RL * RLLR -
        2.0 * c22 * m1 * m2 * pap1 * papb * RL * RLLR +
        2.0 * c12 * m1 * m2 * pow2(papb) * RL * RLLR +
        2.0 * c22 * m1 * m2 * pow2(papb) * RL * RLLR -
        2.0 * c12 * m1 * m2 * papb * pbp1 * RL * RLLR -
        2.0 * c22 * m1 * m2 * papb * pbp1 * RL * RLLR +
        c0 * LRRR * m1s * m2 * ml2 * papb * RR -
        c0 * LLLR * m1 * ml2 * p1p2 * papb * RR +
        2.0 * c0 * LLLR * m1 * ml2 * pap2 * papb * RR -
        2.0 * c0 * LRRR * m2 * ml2 * pap1 * pbp1 * RR -
        c0 * LLLR * m1 * ml2 * pap2 * pbp1 * RR +
        c0 * LLLR * m1 * ml2 * pap1 * pbp2 * RR -
        c0 * LL * m1 * ml2 * p1p2 * papb * RRRL +
        2.0 * c0 * LL * m1 * ml2 * pap2 * papb * RRRL -
        c0 * LL * m1 * ml2 * pap2 * pbp1 * RRRL +
        c0 * LL * m1 * ml2 * pap1 * pbp2 * RRRL +
        (2.0 * c12 * (papb - pbp1) + c22 * (m1s - 2.0 * (pap1 - papb + pbp1))) *
            (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL * RRRR +
        c2 * (LR * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                  (LRRL * m1 * m2 * papb +
                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
              m1 * m1s * m2 * papb * RL * RLLR -
              2.0 * m1 * m2 * pap1 * papb * RL * RLLR +
              2.0 * m1 * m2 * pow2(papb) * RL * RLLR -
              2.0 * m1 * m2 * papb * pbp1 * RL * RLLR +
              LRRR * m1s * m2 * ml2 * papb * RR -
              LLLR * m1 * ml2 * p1p2 * papb * RR +
              2.0 * LLLR * m1 * ml2 * pap2 * papb * RR -
              2.0 * LRRR * m2 * ml2 * pap1 * pbp1 * RR -
              LLLR * m1 * ml2 * pap2 * pbp1 * RR +
              LLLR * m1 * ml2 * pap1 * pbp2 * RR +
              LL * ml2 * (m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLL +
                          m1 * (-(p1p2 * papb) + 2.0 * pap2 * papb -
                                pap2 * pbp1 + pap1 * pbp2) *
                              RRRL) +
              (m1s - 2.0 * (pap1 - papb + pbp1)) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
        4.0 * c00 *
            (LR * (LRRL * m1 * m2 * papb +
                   LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
             RL * (m1 * m2 * papb * RLLR +
                   (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR)))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * pow2(cF) * ivt1s1 * ivu2s2 * MPIsI * NC *
        (LR * (LRRL * m1 * m2 * papb +
               LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
         RL * (m1 * m2 * papb * RLLR +
               (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
      4.0;

  return real(me0 + me1);
}

double FI::Mut_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -((cA - 2.0 * cF) * cF * ivt1s2 * ivu2s1 * MPIsI * NC *
        (c0 * ml2 * (LL * LRLL * m1 * (m2s * papb - 2.0 * pap2 * pbp2) +
                     (LLRL * m2 * (-(p1p2 * papb) + 2.0 * pap1 * papb +
                                   pap2 * pbp1 - pap1 * pbp2) +
                      m1 * (m2s * papb - 2.0 * pap2 * pbp2) * RLRR) *
                         RR +
                     LL * m2 * (-(p1p2 * papb) + 2.0 * pap1 * papb +
                                pap2 * pbp1 - pap1 * pbp2) *
                         RRLR) +
         (4.0 * c00 + 2.0 * c12 * (papb - pbp2) +
          c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
             (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              LRRL * m1 * m2 * papb * RL + LR * m1 * m2 * papb * RLLR -
              p1p2 * papb * RL * RRRR + pap2 * pbp1 * RL * RRRR +
              pap1 * pbp2 * RL * RRRR) +
         c2 * (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                   (m2s - 2.0 * (pap2 - papb + pbp2)) +
               LRRL * m1 * m2 * m2s * papb * RL -
               2.0 * LRRL * m1 * m2 * pap2 * papb * RL +
               2.0 * LRRL * m1 * m2 * pow2(papb) * RL -
               2.0 * LRRL * m1 * m2 * papb * pbp2 * RL +
               LR * m1 * m2 * m2s * papb * RLLR -
               2.0 * LR * m1 * m2 * pap2 * papb * RLLR +
               2.0 * LR * m1 * m2 * pow2(papb) * RLLR -
               2.0 * LR * m1 * m2 * papb * pbp2 * RLLR -
               LLRL * m2 * ml2 * p1p2 * papb * RR +
               2.0 * LLRL * m2 * ml2 * pap1 * papb * RR +
               LLRL * m2 * ml2 * pap2 * pbp1 * RR -
               LLRL * m2 * ml2 * pap1 * pbp2 * RR +
               m1 * m2s * ml2 * papb * RLRR * RR -
               2.0 * m1 * ml2 * pap2 * pbp2 * RLRR * RR +
               LL * ml2 * (LRLL * m1 * (m2s * papb - 2.0 * pap2 * pbp2) +
                           m2 * (-(p1p2 * papb) + 2.0 * pap1 * papb +
                                 pap2 * pbp1 - pap1 * pbp2) *
                               RRLR) +
               (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                   (m2s - 2.0 * (pap2 - papb + pbp2)) * RL * RRRR))) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -(c00 * (cA - 2.0 * cF) * cF * ivt1s2 * ivu2s1 * MPIsI * NC *
        (LLLL * LR * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) -
         m1 * m2 * papb * (LRRL * RL + LR * RLLR) +
         (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRR)) /
      8.0;

  return real(me0 + me1);
}

// Triangle 3
/************************************************/
// terms proportional to eps^2 contain only c00 which is IR finite.
double FI::Mtt_V3a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * cA * cF * ivt1s1 * ivt1s2 * MPIsI * pbp2 *
       (8.0 * c00 * pap1 +
        c2 * (2.0 * m2s * pap1 + (m1 * ml3 - 6.0 * pap1) * (pap2 - papb) -
              4.0 * pap1 * pbp2) +
        2.0 * pap1 * (-((c0 + c1 + 2.0 * c12) * (pap2 - papb)) +
                      c22 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (-3.0 * cA * cF * ivt1s1 * ivt1s2 * MPIsI * pbp2 *
                         (8.0 * c00 * pap1 +
                          pap1 * (2.0 * c12 * (-pap2 + papb) +
                                  c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) +
                          c2 * (m1 * ml3 * (pap2 - papb) +
                                pap1 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
                         (LLLL + LRLR + RLRL + RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Mtu_V3a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * cA * cF * ivt1s1 * ivu2s2 * MPIsI *
       (8.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
        2.0 * (-((c0 + c1 + 2.0 * c12) * (pap2 - papb)) +
               c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            (LRRL * m1 * m2 * papb +
             LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
             m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
             pap1 * pbp2 * RRRR) +
        c2 *
            (-(LLLL * m1 * m2s * ml3 * papb) - 2.0 * LLLL * m2s * p1p2 * papb +
             6.0 * LLLL * p1p2 * pap2 * papb - 6.0 * LLLL * p1p2 * pow2(papb) +
             2.0 * LLLL * m2s * pap2 * pbp1 - 6.0 * LLLL * pow2(pap2) * pbp1 +
             6.0 * LLLL * pap2 * papb * pbp1 + 2.0 * LLLL * m2s * pap1 * pbp2 +
             2.0 * LLLL * m1 * ml3 * pap2 * pbp2 -
             6.0 * LLLL * pap1 * pap2 * pbp2 + 4.0 * LLLL * p1p2 * papb * pbp2 +
             6.0 * LLLL * pap1 * papb * pbp2 - 4.0 * LLLL * pap2 * pbp1 * pbp2 -
             4.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 * (2.0 * m1 * papb *
                              (m2s - 3.0 * pap2 + 3.0 * papb - 2.0 * pbp2) +
                          ml3 * (p1p2 * papb + pap2 * pbp1 - 2.0 * papb * pbp1 -
                                 pap1 * pbp2)) +
             2.0 * m1 * m2 * m2s * papb * RLLR + m2 * ml3 * p1p2 * papb * RLLR -
             6.0 * m1 * m2 * pap2 * papb * RLLR +
             6.0 * m1 * m2 * pow2(papb) * RLLR + m2 * ml3 * pap2 * pbp1 * RLLR -
             2.0 * m2 * ml3 * papb * pbp1 * RLLR -
             m2 * ml3 * pap1 * pbp2 * RLLR -
             4.0 * m1 * m2 * papb * pbp2 * RLLR +
             (2.0 * (m2s - 3.0 * pap2 + 3.0 * papb - 2.0 * pbp2) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              m1 * ml3 * (-(m2s * papb) + 2.0 * pap2 * pbp2)) *
                 RRRR))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      (3.0 * cA * cF * ivt1s1 * ivu2s2 * MPIsI *
       (8.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
        (2.0 * c12 * (-pap2 + papb) +
         c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            (LRRL * m1 * m2 * papb +
             LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
             m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
             pap1 * pbp2 * RRRR) +
        c2 *
            (-(LLLL * m1 * m2s * ml3 * papb) - LLLL * m2s * p1p2 * papb +
             2.0 * LLLL * p1p2 * pap2 * papb - 2.0 * LLLL * p1p2 * pow2(papb) +
             LLLL * m2s * pap2 * pbp1 - 2.0 * LLLL * pow2(pap2) * pbp1 +
             2.0 * LLLL * pap2 * papb * pbp1 + LLLL * m2s * pap1 * pbp2 +
             2.0 * LLLL * m1 * ml3 * pap2 * pbp2 -
             2.0 * LLLL * pap1 * pap2 * pbp2 + 2.0 * LLLL * p1p2 * papb * pbp2 +
             2.0 * LLLL * pap1 * papb * pbp2 - 2.0 * LLLL * pap2 * pbp1 * pbp2 -
             2.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 * (ml3 * (p1p2 * papb + pap2 * pbp1 - 2.0 * papb * pbp1 -
                                 pap1 * pbp2) +
                          m1 * papb * (m2s - 2.0 * (pap2 - papb + pbp2))) +
             m1 * m2 * m2s * papb * RLLR + m2 * ml3 * p1p2 * papb * RLLR -
             2.0 * m1 * m2 * pap2 * papb * RLLR +
             2.0 * m1 * m2 * pow2(papb) * RLLR + m2 * ml3 * pap2 * pbp1 * RLLR -
             2.0 * m2 * ml3 * papb * pbp1 * RLLR -
             m2 * ml3 * pap1 * pbp2 * RLLR -
             2.0 * m1 * m2 * papb * pbp2 * RLLR +
             (m1 * ml3 * (-(m2s * papb) + 2.0 * pap2 * pbp2) +
              (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                  (m2s - 2.0 * (pap2 - papb + pbp2))) *
                 RRRR))) /
      8.0;

  return real(me0 + me1);
}

double FI::Muu_V3b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * cA * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 *
       (c2 * (-2.0 * (m2s - 2.0 * pap2) * pbp1 +
              (m1 * ml3 - 6.0 * pbp1) * (papb - pbp2)) -
        2.0 * pbp1 * (4.0 * c00 + (c0 + c1 + 2.0 * c12) * (papb - pbp2) +
                      c22 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      (3.0 * cA * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 *
       (-((8.0 * c00 + c22 * (m2s - 2.0 * pap2) + 2.0 * (c12 + c22) * papb) *
          pbp1) +
        2.0 * (c12 + c22) * pbp1 * pbp2 +
        c2 * (m1 * ml3 * (papb - pbp2) -
              pbp1 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      4.0;

  return real(me0 + me1);
}

double FI::Mut_V3b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * cA * cF * ivt1s2 * ivu2s1 * MPIsI *
       (8.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
        2.0 * ((c0 + c1 + 2.0 * c12) * (papb - pbp2) +
               c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            (LRRL * m1 * m2 * papb +
             LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
             m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
             pap1 * pbp2 * RRRR) +
        c2 *
            (-(LLLL * m1 * m2s * ml3 * papb) - 2.0 * LLLL * m2s * p1p2 * papb +
             4.0 * LLLL * p1p2 * pap2 * papb - 6.0 * LLLL * p1p2 * pow2(papb) +
             2.0 * LLLL * m2s * pap2 * pbp1 - 4.0 * LLLL * pow2(pap2) * pbp1 +
             6.0 * LLLL * pap2 * papb * pbp1 + 2.0 * LLLL * m2s * pap1 * pbp2 +
             2.0 * LLLL * m1 * ml3 * pap2 * pbp2 -
             4.0 * LLLL * pap1 * pap2 * pbp2 + 6.0 * LLLL * p1p2 * papb * pbp2 +
             6.0 * LLLL * pap1 * papb * pbp2 - 6.0 * LLLL * pap2 * pbp1 * pbp2 -
             6.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 * (2.0 * m1 * papb *
                              (m2s - 2.0 * pap2 + 3.0 * papb - 3.0 * pbp2) +
                          ml3 * (p1p2 * papb - 2.0 * pap1 * papb - pap2 * pbp1 +
                                 pap1 * pbp2)) +
             2.0 * m1 * m2 * m2s * papb * RLLR + m2 * ml3 * p1p2 * papb * RLLR -
             2.0 * m2 * ml3 * pap1 * papb * RLLR -
             4.0 * m1 * m2 * pap2 * papb * RLLR +
             6.0 * m1 * m2 * pow2(papb) * RLLR - m2 * ml3 * pap2 * pbp1 * RLLR +
             m2 * ml3 * pap1 * pbp2 * RLLR -
             6.0 * m1 * m2 * papb * pbp2 * RLLR +
             (2.0 * (m2s - 2.0 * pap2 + 3.0 * papb - 3.0 * pbp2) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              m1 * ml3 * (-(m2s * papb) + 2.0 * pap2 * pbp2)) *
                 RRRR))) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      (3.0 * cA * cF * ivt1s2 * ivu2s1 * MPIsI *
       (8.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
        (2.0 * c12 * (papb - pbp2) + c22 * (m2s - 2.0 * (pap2 - papb + pbp2))) *
            (LRRL * m1 * m2 * papb +
             LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
             m1 * m2 * papb * RLLR - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
             pap1 * pbp2 * RRRR) +
        c2 *
            (-(LLLL * m1 * m2s * ml3 * papb) - LLLL * m2s * p1p2 * papb +
             2.0 * LLLL * p1p2 * pap2 * papb - 2.0 * LLLL * p1p2 * pow2(papb) +
             LLLL * m2s * pap2 * pbp1 - 2.0 * LLLL * pow2(pap2) * pbp1 +
             2.0 * LLLL * pap2 * papb * pbp1 + LLLL * m2s * pap1 * pbp2 +
             2.0 * LLLL * m1 * ml3 * pap2 * pbp2 -
             2.0 * LLLL * pap1 * pap2 * pbp2 + 2.0 * LLLL * p1p2 * papb * pbp2 +
             2.0 * LLLL * pap1 * papb * pbp2 - 2.0 * LLLL * pap2 * pbp1 * pbp2 -
             2.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 * (ml3 * (p1p2 * papb - 2.0 * pap1 * papb - pap2 * pbp1 +
                                 pap1 * pbp2) +
                          m1 * papb * (m2s - 2.0 * (pap2 - papb + pbp2))) +
             m1 * m2 * m2s * papb * RLLR + m2 * ml3 * p1p2 * papb * RLLR -
             2.0 * m2 * ml3 * pap1 * papb * RLLR -
             2.0 * m1 * m2 * pap2 * papb * RLLR +
             2.0 * m1 * m2 * pow2(papb) * RLLR - m2 * ml3 * pap2 * pbp1 * RLLR +
             m2 * ml3 * pap1 * pbp2 * RLLR -
             2.0 * m1 * m2 * papb * pbp2 * RLLR +
             (m1 * ml3 * (-(m2s * papb) + 2.0 * pap2 * pbp2) +
              (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
                  (m2s - 2.0 * (pap2 - papb + pbp2))) *
                 RRRR))) /
      8.0;

  return real(me0 + me1);
}

// Triangle 4
/************************************************/
double FI::Mtt_V4a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * cA * cF * ivt1s2 * ivt2s1 * MPIsI * pbp2 *
       (4.0 * c00 * pap1 + c0 * m1 * ml1 * (-pap2 + papb) +
        c2 * (-(m2s * pap1) + (m1 * ml1 - 2.0 * pap1) * (pap2 - papb) +
              2.0 * pap1 * pbp2) +
        pap1 * (2.0 * c12 * (-pap2 + papb) +
                c22 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (-3.0 * c00 * cA * cF * ivt1s2 * ivt2s1 * MPIsI * pap1 *
                         pbp2 * (LLLL + LRLR + RLRL + RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Mtu_V4a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * cA * cF * ivt2s1 * ivu2s2 * MPIsI *
       (c22 * LRRL * m1 * m2 * m2s * papb + c0 * LLLL * m1 * m2s * ml1 * papb -
        c22 * LLLL * m2s * p1p2 * papb - c0 * LRRL * m2 * ml1 * p1p2 * papb -
        2.0 * c12 * LRRL * m1 * m2 * pap2 * papb -
        2.0 * c22 * LRRL * m1 * m2 * pap2 * papb +
        2.0 * c12 * LLLL * p1p2 * pap2 * papb +
        2.0 * c22 * LLLL * p1p2 * pap2 * papb +
        2.0 * c12 * LRRL * m1 * m2 * pow2(papb) +
        2.0 * c22 * LRRL * m1 * m2 * pow2(papb) -
        2.0 * c12 * LLLL * p1p2 * pow2(papb) -
        2.0 * c22 * LLLL * p1p2 * pow2(papb) + c22 * LLLL * m2s * pap2 * pbp1 -
        c0 * LRRL * m2 * ml1 * pap2 * pbp1 -
        2.0 * c12 * LLLL * pow2(pap2) * pbp1 -
        2.0 * c22 * LLLL * pow2(pap2) * pbp1 +
        2.0 * c0 * LRRL * m2 * ml1 * papb * pbp1 +
        2.0 * c12 * LLLL * pap2 * papb * pbp1 +
        2.0 * c22 * LLLL * pap2 * papb * pbp1 + c22 * LLLL * m2s * pap1 * pbp2 +
        c0 * LRRL * m2 * ml1 * pap1 * pbp2 -
        2.0 * c0 * LLLL * m1 * ml1 * pap2 * pbp2 -
        2.0 * c12 * LLLL * pap1 * pap2 * pbp2 -
        2.0 * c22 * LLLL * pap1 * pap2 * pbp2 -
        2.0 * c22 * LRRL * m1 * m2 * papb * pbp2 +
        2.0 * c22 * LLLL * p1p2 * papb * pbp2 +
        2.0 * c12 * LLLL * pap1 * papb * pbp2 +
        2.0 * c22 * LLLL * pap1 * papb * pbp2 -
        2.0 * c22 * LLLL * pap2 * pbp1 * pbp2 -
        2.0 * c22 * LLLL * pap1 * pow2(pbp2) +
        c22 * m1 * m2 * m2s * papb * RLLR - c0 * m2 * ml1 * p1p2 * papb * RLLR -
        2.0 * c12 * m1 * m2 * pap2 * papb * RLLR -
        2.0 * c22 * m1 * m2 * pap2 * papb * RLLR +
        2.0 * c12 * m1 * m2 * pow2(papb) * RLLR +
        2.0 * c22 * m1 * m2 * pow2(papb) * RLLR -
        c0 * m2 * ml1 * pap2 * pbp1 * RLLR +
        2.0 * c0 * m2 * ml1 * papb * pbp1 * RLLR +
        c0 * m2 * ml1 * pap1 * pbp2 * RLLR -
        2.0 * c22 * m1 * m2 * papb * pbp2 * RLLR +
        (c0 * m1 * ml1 * (m2s * papb - 2.0 * pap2 * pbp2) +
         (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
             (2.0 * c12 * (-pap2 + papb) +
              c22 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
            RRRR +
        4.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) -
        c2 *
            (LLLL * m1 * m2s * ml1 * papb - LLLL * m2s * p1p2 * papb -
             2.0 * LLLL * p1p2 * pap2 * papb + 2.0 * LLLL * p1p2 * pow2(papb) +
             LLLL * m2s * pap2 * pbp1 + 2.0 * LLLL * pow2(pap2) * pbp1 -
             2.0 * LLLL * pap2 * papb * pbp1 + LLLL * m2s * pap1 * pbp2 -
             2.0 * LLLL * m1 * ml1 * pap2 * pbp2 +
             2.0 * LLLL * pap1 * pap2 * pbp2 + 2.0 * LLLL * p1p2 * papb * pbp2 -
             2.0 * LLLL * pap1 * papb * pbp2 - 2.0 * LLLL * pap2 * pbp1 * pbp2 -
             2.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 *
                 (m1 * papb * (m2s + 2.0 * pap2 - 2.0 * papb - 2.0 * pbp2) +
                  ml1 * (-(p1p2 * papb) - pap2 * pbp1 + 2.0 * papb * pbp1 +
                         pap1 * pbp2)) +
             m1 * m2 * m2s * papb * RLLR - m2 * ml1 * p1p2 * papb * RLLR +
             2.0 * m1 * m2 * pap2 * papb * RLLR -
             2.0 * m1 * m2 * pow2(papb) * RLLR - m2 * ml1 * pap2 * pbp1 * RLLR +
             2.0 * m2 * ml1 * papb * pbp1 * RLLR +
             m2 * ml1 * pap1 * pbp2 * RLLR -
             2.0 * m1 * m2 * papb * pbp2 * RLLR +
             ((m2s + 2.0 * pap2 - 2.0 * papb - 2.0 * pbp2) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              m1 * ml1 * (m2s * papb - 2.0 * pap2 * pbp2)) *
                 RRRR))) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * cA * cF * ivt2s1 * ivu2s2 * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        8.0;

  return real(me0 + me1);
}

double FI::Muu_V4b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (3.0 * cA * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 *
       (c0 * m1 * ml1 * (papb - pbp2) +
        pbp1 * (4.0 * c00 + c22 * m2s - 2.0 * c22 * pap2 + 2.0 * c12 * papb +
                2.0 * c22 * papb - 2.0 * (c12 + c22) * pbp2) +
        c2 * (-(pbp1 * (m2s - 2.0 * (pap2 + papb - pbp2))) +
              m1 * ml1 * (-papb + pbp2))) *
       (LLLL + LRLR + RLRL + RRRR)) /
      8.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (-3.0 * c00 * cA * cF * ivu2s1 * ivu2s2 * MPIsI * pap2 *
                         pbp1 * (LLLL + LRLR + RLRL + RRRR)) /
                        4.0;

  return real(me0 + me1);
}

double FI::Mut_V4b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                        double ml2s, double ml3s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0, c00, c1, c2, c11, c12, c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      (-3.0 * cA * cF * ivt1s2 * ivu2s1 * MPIsI *
       (c22 * LRRL * m1 * m2 * m2s * papb + c0 * LLLL * m1 * m2s * ml1 * papb -
        c22 * LLLL * m2s * p1p2 * papb - c0 * LRRL * m2 * ml1 * p1p2 * papb +
        2.0 * c0 * LRRL * m2 * ml1 * pap1 * papb -
        2.0 * c22 * LRRL * m1 * m2 * pap2 * papb +
        2.0 * c22 * LLLL * p1p2 * pap2 * papb +
        2.0 * c12 * LRRL * m1 * m2 * pow2(papb) +
        2.0 * c22 * LRRL * m1 * m2 * pow2(papb) -
        2.0 * c12 * LLLL * p1p2 * pow2(papb) -
        2.0 * c22 * LLLL * p1p2 * pow2(papb) + c22 * LLLL * m2s * pap2 * pbp1 +
        c0 * LRRL * m2 * ml1 * pap2 * pbp1 -
        2.0 * c22 * LLLL * pow2(pap2) * pbp1 +
        2.0 * c12 * LLLL * pap2 * papb * pbp1 +
        2.0 * c22 * LLLL * pap2 * papb * pbp1 + c22 * LLLL * m2s * pap1 * pbp2 -
        c0 * LRRL * m2 * ml1 * pap1 * pbp2 -
        2.0 * c0 * LLLL * m1 * ml1 * pap2 * pbp2 -
        2.0 * c22 * LLLL * pap1 * pap2 * pbp2 -
        2.0 * c12 * LRRL * m1 * m2 * papb * pbp2 -
        2.0 * c22 * LRRL * m1 * m2 * papb * pbp2 +
        2.0 * c12 * LLLL * p1p2 * papb * pbp2 +
        2.0 * c22 * LLLL * p1p2 * papb * pbp2 +
        2.0 * c12 * LLLL * pap1 * papb * pbp2 +
        2.0 * c22 * LLLL * pap1 * papb * pbp2 -
        2.0 * c12 * LLLL * pap2 * pbp1 * pbp2 -
        2.0 * c22 * LLLL * pap2 * pbp1 * pbp2 -
        2.0 * c12 * LLLL * pap1 * pow2(pbp2) -
        2.0 * c22 * LLLL * pap1 * pow2(pbp2) +
        c22 * m1 * m2 * m2s * papb * RLLR - c0 * m2 * ml1 * p1p2 * papb * RLLR +
        2.0 * c0 * m2 * ml1 * pap1 * papb * RLLR -
        2.0 * c22 * m1 * m2 * pap2 * papb * RLLR +
        2.0 * c12 * m1 * m2 * pow2(papb) * RLLR +
        2.0 * c22 * m1 * m2 * pow2(papb) * RLLR +
        c0 * m2 * ml1 * pap2 * pbp1 * RLLR -
        c0 * m2 * ml1 * pap1 * pbp2 * RLLR -
        2.0 * c12 * m1 * m2 * papb * pbp2 * RLLR -
        2.0 * c22 * m1 * m2 * papb * pbp2 * RLLR +
        (c0 * m1 * ml1 * (m2s * papb - 2.0 * pap2 * pbp2) +
         (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) *
             (2.0 * c12 * (papb - pbp2) +
              c22 * (m2s - 2.0 * (pap2 - papb + pbp2)))) *
            RRRR +
        4.0 * c00 * (LRRL * m1 * m2 * papb +
                     LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                     pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) -
        c2 *
            (LLLL * m1 * m2s * ml1 * papb - LLLL * m2s * p1p2 * papb +
             2.0 * LLLL * p1p2 * pap2 * papb + 2.0 * LLLL * p1p2 * pow2(papb) +
             LLLL * m2s * pap2 * pbp1 - 2.0 * LLLL * pow2(pap2) * pbp1 -
             2.0 * LLLL * pap2 * papb * pbp1 + LLLL * m2s * pap1 * pbp2 -
             2.0 * LLLL * m1 * ml1 * pap2 * pbp2 -
             2.0 * LLLL * pap1 * pap2 * pbp2 - 2.0 * LLLL * p1p2 * papb * pbp2 -
             2.0 * LLLL * pap1 * papb * pbp2 + 2.0 * LLLL * pap2 * pbp1 * pbp2 +
             2.0 * LLLL * pap1 * pow2(pbp2) +
             LRRL * m2 * (m1 * papb * (m2s - 2.0 * (pap2 + papb - pbp2)) +
                          ml1 * (-(p1p2 * papb) + 2.0 * pap1 * papb +
                                 pap2 * pbp1 - pap1 * pbp2)) +
             m1 * m2 * m2s * papb * RLLR - m2 * ml1 * p1p2 * papb * RLLR +
             2.0 * m2 * ml1 * pap1 * papb * RLLR -
             2.0 * m1 * m2 * pap2 * papb * RLLR -
             2.0 * m1 * m2 * pow2(papb) * RLLR + m2 * ml1 * pap2 * pbp1 * RLLR -
             m2 * ml1 * pap1 * pbp2 * RLLR +
             2.0 * m1 * m2 * papb * pbp2 * RLLR +
             ((m2s - 2.0 * (pap2 + papb - pbp2)) *
                  (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
              m1 * ml1 * (m2s * papb - 2.0 * pap2 * pbp2)) *
                 RRRR))) /
      16.0;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = (3.0 * c00 * cA * cF * ivt1s2 * ivu2s1 * MPIsI *
                         (LRRL * m1 * m2 * papb +
                          LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                          m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                          pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        8.0;

  return real(me0 + me1);
}

// Box1
/************************************************/
double FI::Mtt_B1_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (-3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * MPIsI *
       (d33 * (LLLL * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
               (m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                                     2.0 * papb * pbp2) *
                   (LRLR + RLRL) +
               (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR) +
        4.0 * d00 * (LLLL * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                     (p1p2 * papb - pap2 * pbp1) * RRRR +
                     pap1 * pbp2 * (LRLR + RLRL + RRRR)) +
        2.0 *
            (d12 * LLLL * p1p2 * pow2(papb) + d13 * LLLL * p1p2 * pow2(papb) +
             d2 * LLLL * p1p2 * pow2(papb) + d22 * LLLL * p1p2 * pow2(papb) -
             d13 * LLLL * p1p2 * papb * pbp1 - d12 * LLLL * pap2 * papb * pbp1 -
             d13 * LLLL * pap2 * papb * pbp1 - d2 * LLLL * pap2 * papb * pbp1 -
             d22 * LLLL * pap2 * papb * pbp1 + d13 * LLLL * pap2 * pow2(pbp1) -
             d13 * LRLR * m1s * papb * pbp2 +
             2.0 * d1 * LLLL * pap1 * papb * pbp2 +
             d12 * LLLL * pap1 * papb * pbp2 + d13 * LLLL * pap1 * papb * pbp2 +
             d2 * LLLL * pap1 * papb * pbp2 + d22 * LLLL * pap1 * papb * pbp2 +
             2.0 * d1 * LRLR * pap1 * papb * pbp2 +
             2.0 * d12 * LRLR * pap1 * papb * pbp2 +
             2.0 * d13 * LRLR * pap1 * papb * pbp2 +
             2.0 * d2 * LRLR * pap1 * papb * pbp2 +
             2.0 * d22 * LRLR * pap1 * papb * pbp2 -
             d13 * LLLL * pap1 * pbp1 * pbp2 - d13 * m1s * papb * pbp2 * RLRL +
             2.0 * d1 * pap1 * papb * pbp2 * RLRL +
             2.0 * d12 * pap1 * papb * pbp2 * RLRL +
             2.0 * d13 * pap1 * papb * pbp2 * RLRL +
             2.0 * d2 * pap1 * papb * pbp2 * RLRL +
             2.0 * d22 * pap1 * papb * pbp2 * RLRL +
             (((d12 + d13 + d2 + d22) * papb - d13 * pbp1) *
                  (p1p2 * papb - pap2 * pbp1) +
              pap1 * ((2.0 * d1 + d12 + d13 + d2 + d22) * papb - d13 * pbp1) *
                  pbp2) *
                 RRRR +
             d23 * (-(LLLL * (pap1 - 2.0 * papb + pbp1) *
                      (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) -
                    (p1p2 * pap1 * papb - pap1 * pap2 * pbp1 +
                     (pow2(pap1) + m1s * papb - 4.0 * pap1 * papb) * pbp2) *
                        (LRLR + RLRL) -
                    (pap1 - 2.0 * papb + pbp1) *
                        (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR) +
             d3 *
                 (-(LLLL * (pap1 - papb + pbp1) * (p1p2 * papb - pap2 * pbp1)) +
                  LLLL * (m1s * papb - pap1 * (pap1 - papb + pbp1)) * pbp2 -
                  pap1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                          2.0 * papb * pbp2) *
                      (LRLR + RLRL) +
                  (-((pap1 - papb + pbp1) * (p1p2 * papb - pap2 * pbp1)) +
                   (m1s * papb - pap1 * (pap1 - papb + pbp1)) * pbp2) *
                      RRRR)))) /
      8.0;

  return real(me0);
}

double FI::Muu_B1_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (-3.0 * (cA - 2.0 * cF) * cF * ivu2s2 * MPIsI *
       (4.0 * d00 * (LLLL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                     pap2 * pbp1 * (LRLR + RLRL) +
                     (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR) +
        d33 * (LLLL * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
               (m1s - 2.0 * pbp1) * (p1p2 * papb - 2.0 * pap2 * papb +
                                     pap2 * pbp1 - pap1 * pbp2) *
                   (LRLR + RLRL) +
               (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR) +
        2.0 *
            (-(d3 * LLLL * p1p2 * pap1 * papb) + d3 * LLLL * m1s * pap2 * papb +
             d12 * LLLL * p1p2 * pow2(papb) + d2 * LLLL * p1p2 * pow2(papb) +
             d22 * LLLL * p1p2 * pow2(papb) + d3 * LLLL * p1p2 * pow2(papb) -
             d3 * LLLL * pap1 * pap2 * pbp1 - d3 * LLLL * p1p2 * papb * pbp1 -
             d3 * LRLR * p1p2 * papb * pbp1 +
             2.0 * d1 * LLLL * pap2 * papb * pbp1 +
             d12 * LLLL * pap2 * papb * pbp1 + d2 * LLLL * pap2 * papb * pbp1 +
             d22 * LLLL * pap2 * papb * pbp1 + d3 * LLLL * pap2 * papb * pbp1 +
             2.0 * d1 * LRLR * pap2 * papb * pbp1 +
             2.0 * d12 * LRLR * pap2 * papb * pbp1 +
             2.0 * d2 * LRLR * pap2 * papb * pbp1 +
             2.0 * d22 * LRLR * pap2 * papb * pbp1 +
             2.0 * d3 * LRLR * pap2 * papb * pbp1 -
             d3 * LLLL * pap2 * pow2(pbp1) - d3 * LRLR * pap2 * pow2(pbp1) +
             d3 * LLLL * pow2(pap1) * pbp2 - d12 * LLLL * pap1 * papb * pbp2 -
             d2 * LLLL * pap1 * papb * pbp2 - d22 * LLLL * pap1 * papb * pbp2 -
             d3 * LLLL * pap1 * papb * pbp2 + d3 * LLLL * pap1 * pbp1 * pbp2 +
             d3 * LRLR * pap1 * pbp1 * pbp2 - d3 * p1p2 * papb * pbp1 * RLRL +
             2.0 * d1 * pap2 * papb * pbp1 * RLRL +
             2.0 * d12 * pap2 * papb * pbp1 * RLRL +
             2.0 * d2 * pap2 * papb * pbp1 * RLRL +
             2.0 * d22 * pap2 * papb * pbp1 * RLRL +
             2.0 * d3 * pap2 * papb * pbp1 * RLRL -
             d3 * pap2 * pow2(pbp1) * RLRL + d3 * pap1 * pbp1 * pbp2 * RLRL +
             (papb * (2.0 * d1 * pap2 * pbp1 +
                      (d12 + d2 + d22) *
                          (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) +
              d3 * (m1s * pap2 * papb + p1p2 * papb * (-pap1 + papb - pbp1) +
                    (pap1 - papb + pbp1) * (-(pap2 * pbp1) + pap1 * pbp2))) *
                 RRRR +
             d13 * (LLLL * (pap1 - papb) *
                        (-(p1p2 * papb) - pap2 * pbp1 + pap1 * pbp2) -
                    pap2 * papb * (m1s - 2.0 * pbp1) * (LRLR + RLRL) +
                    (-pap1 + papb) * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                        RRRR) +
             d23 *
                 (LLLL * (pap1 - 2.0 * papb + pbp1) *
                      (-(p1p2 * papb) - pap2 * pbp1 + pap1 * pbp2) -
                  (m1s * pap2 * papb +
                   pbp1 * (p1p2 * papb - 4.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2)) *
                      (LRLR + RLRL) +
                  (pap1 - 2.0 * papb + pbp1) *
                      (-(p1p2 * papb) - pap2 * pbp1 + pap1 * pbp2) * RRRR)))) /
      8.0;

  return real(me0);
}

double FI::Mtu_B1_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivu2s2 * MPIsI *
       (2.0 * d00 * (LRRL * m1 * m2 * papb -
                     2.0 * LLLL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                     m1 * m2 * papb * RLLR -
                     2.0 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR) +
        d33 * (LLLL * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (-(p1p2 * papb) - pap2 * pbp1 + pap1 * pbp2) -
               2.0 * m1 * m2 * (pap1 - papb) * (papb - pbp1) * (LRRL + RLLR) +
               (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (-(p1p2 * papb) - pap2 * pbp1 + pap1 * pbp2) * RRRR) +
        2.0 *
            (d1 * LRRL * m1 * m2 * pow2(papb) +
             d12 * LRRL * m1 * m2 * pow2(papb) +
             d13 * LRRL * m1 * m2 * pow2(papb) +
             d2 * LRRL * m1 * m2 * pow2(papb) +
             d22 * LRRL * m1 * m2 * pow2(papb) - d1 * LLLL * p1p2 * pow2(papb) -
             d12 * LLLL * p1p2 * pow2(papb) - d13 * LLLL * p1p2 * pow2(papb) -
             d2 * LLLL * p1p2 * pow2(papb) - d22 * LLLL * p1p2 * pow2(papb) -
             d13 * LRRL * m1 * m2 * papb * pbp1 +
             d13 * LLLL * p1p2 * papb * pbp1 + d1 * LLLL * pap2 * papb * pbp1 -
             d12 * LLLL * pap2 * papb * pbp1 - d13 * LLLL * pap2 * papb * pbp1 -
             d2 * LLLL * pap2 * papb * pbp1 - d22 * LLLL * pap2 * papb * pbp1 +
             d13 * LLLL * pap2 * pow2(pbp1) + d1 * LLLL * pap1 * papb * pbp2 +
             d12 * LLLL * pap1 * papb * pbp2 + d13 * LLLL * pap1 * papb * pbp2 +
             d2 * LLLL * pap1 * papb * pbp2 + d22 * LLLL * pap1 * papb * pbp2 -
             d13 * LLLL * pap1 * pbp1 * pbp2 +
             d1 * m1 * m2 * pow2(papb) * RLLR +
             d12 * m1 * m2 * pow2(papb) * RLLR +
             d13 * m1 * m2 * pow2(papb) * RLLR +
             d2 * m1 * m2 * pow2(papb) * RLLR +
             d22 * m1 * m2 * pow2(papb) * RLLR -
             d13 * m1 * m2 * papb * pbp1 * RLLR -
             (d1 * papb * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) +
              ((d12 + d13 + d2 + d22) * papb - d13 * pbp1) *
                  (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
                 RRRR -
             d23 * (pap1 - 2.0 * papb + pbp1) *
                 (LRRL * m1 * m2 * papb -
                  LLLL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                  m1 * m2 * papb * RLLR - p1p2 * papb * RRRR -
                  pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR) +
             d3 * (LRRL * m1 * m2 * papb * (-pap1 + papb) +
                   LLLL *
                       (-(m1s * pap2 * papb) +
                        p1p2 * papb * (pap1 - papb + pbp1) -
                        (pap1 - papb + pbp1) * (-(pap2 * pbp1) + pap1 * pbp2)) +
                   m1 * m2 * papb * (-pap1 + papb) * RLLR +
                   (-(m1s * pap2 * papb) + p1p2 * papb * (pap1 - papb + pbp1) -
                    (pap1 - papb + pbp1) * (-(pap2 * pbp1) + pap1 * pbp2)) *
                       RRRR)))) /
      8.0;

  return real(me0);
}

double FI::Mut_B1_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * MPIsI *
       (2.0 * d00 * (LRRL * m1 * m2 * papb -
                     2.0 * LLLL * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                     m1 * m2 * papb * RLLR -
                     2.0 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR) -
        d33 * (LLLL * (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
               2.0 * m1 * m2 * (pap1 - papb) * (papb - pbp1) * (LRRL + RLLR) +
               (m1s - 2.0 * (pap1 - papb + pbp1)) *
                   (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR) -
        2.0 *
            (-(d3 * LLLL * p1p2 * pap1 * papb) -
             d1 * LRRL * m1 * m2 * pow2(papb) -
             d12 * LRRL * m1 * m2 * pow2(papb) -
             d2 * LRRL * m1 * m2 * pow2(papb) -
             d22 * LRRL * m1 * m2 * pow2(papb) -
             d3 * LRRL * m1 * m2 * pow2(papb) + d1 * LLLL * p1p2 * pow2(papb) +
             d12 * LLLL * p1p2 * pow2(papb) + d2 * LLLL * p1p2 * pow2(papb) +
             d22 * LLLL * p1p2 * pow2(papb) + d3 * LLLL * p1p2 * pow2(papb) +
             d3 * LLLL * pap1 * pap2 * pbp1 +
             d3 * LRRL * m1 * m2 * papb * pbp1 -
             d3 * LLLL * p1p2 * papb * pbp1 - d1 * LLLL * pap2 * papb * pbp1 -
             d12 * LLLL * pap2 * papb * pbp1 - d2 * LLLL * pap2 * papb * pbp1 -
             d22 * LLLL * pap2 * papb * pbp1 - d3 * LLLL * pap2 * papb * pbp1 +
             d3 * LLLL * pap2 * pow2(pbp1) - d3 * LLLL * pow2(pap1) * pbp2 +
             d3 * LLLL * m1s * papb * pbp2 - d1 * LLLL * pap1 * papb * pbp2 +
             d12 * LLLL * pap1 * papb * pbp2 + d2 * LLLL * pap1 * papb * pbp2 +
             d22 * LLLL * pap1 * papb * pbp2 + d3 * LLLL * pap1 * papb * pbp2 -
             d3 * LLLL * pap1 * pbp1 * pbp2 - d1 * m1 * m2 * pow2(papb) * RLLR -
             d12 * m1 * m2 * pow2(papb) * RLLR -
             d2 * m1 * m2 * pow2(papb) * RLLR -
             d22 * m1 * m2 * pow2(papb) * RLLR -
             d3 * m1 * m2 * pow2(papb) * RLLR +
             d3 * m1 * m2 * papb * pbp1 * RLLR +
             (-((p1p2 * papb - pap2 * pbp1) *
                (-((d1 + d12 + d2 + d22) * papb) + d3 * (pap1 - papb + pbp1))) +
              ((-d1 + d12 + d2 + d22) * pap1 * papb +
               d3 * (m1s * papb - pap1 * (pap1 - papb + pbp1))) *
                  pbp2) *
                 RRRR +
             d13 * (pap1 - papb) * (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                                    LLLL * pap2 * pbp1 - LLLL * pap1 * pbp2 +
                                    m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                                    pap2 * pbp1 * RRRR - pap1 * pbp2 * RRRR) +
             d23 * (pap1 - 2.0 * papb + pbp1) *
                 (LRRL * m1 * m2 * papb - LLLL * p1p2 * papb +
                  LLLL * pap2 * pbp1 - LLLL * pap1 * pbp2 +
                  m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                  pap2 * pbp1 * RRRR - pap1 * pbp2 * RRRR)))) /
      8.0;

  return real(me0);
}

// Box2
/************************************************/
// eps part proportional to d00 which is IR fintie and thus gives no
// contribution.
double FI::Mtt_B2_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (-3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * MPIsI *
       (d0 * LL * LRLL * m1 * ml2 * p1p2 * papb +
        2.0 * d23 * LR * LRLR * p1p2 * pap1 * papb -
        d0 * LL * LRLL * m1 * ml2 * pap2 * pbp1 -
        2.0 * d23 * LR * LRLR * pap1 * pap2 * pbp1 +
        4.0 * d00 * LR * LRLR * pap1 * pbp2 +
        d0 * LL * LRLL * m1 * ml2 * pap1 * pbp2 -
        2.0 * d23 * LR * LRLR * pow2(pap1) * pbp2 +
        2.0 * d13 * LR * LRLR * m1s * papb * pbp2 +
        2.0 * d23 * LR * LRLR * m1s * papb * pbp2 -
        2.0 * d0 * LL * LRLL * m1 * ml2 * papb * pbp2 -
        2.0 * d1 * LL * LRLL * m1 * ml2 * papb * pbp2 -
        2.0 * d2 * LL * LRLL * m1 * ml2 * papb * pbp2 -
        4.0 * d13 * LR * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d23 * LR * LRLR * pap1 * pbp1 * pbp2 +
        2.0 * d00 * LRRL * m1 * m2 * papb * RL -
        d0 * LL * m1s * m2 * ml2 * papb * RLLL +
        2.0 * d0 * LL * m2 * ml2 * pap1 * papb * RLLL +
        2.0 * d2 * LL * m2 * ml2 * pap1 * papb * RLLL +
        2.0 * d00 * LR * m1 * m2 * papb * RLLR +
        2.0 * d23 * p1p2 * pap1 * papb * RL * RLRL -
        2.0 * d23 * pap1 * pap2 * pbp1 * RL * RLRL +
        4.0 * d00 * pap1 * pbp2 * RL * RLRL -
        2.0 * d23 * pow2(pap1) * pbp2 * RL * RLRL +
        2.0 * d13 * m1s * papb * pbp2 * RL * RLRL +
        2.0 * d23 * m1s * papb * pbp2 * RL * RLRL -
        4.0 * d13 * pap1 * pbp1 * pbp2 * RL * RLRL -
        4.0 * d23 * pap1 * pbp1 * pbp2 * RL * RLRL +
        d33 *
            (LR * LRLR * (-((m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1)) +
                          (m1s * (pap1 + 2.0 * papb) -
                           2.0 * pap1 * (pap1 + 2.0 * pbp1)) *
                              pbp2) +
             LRRL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL +
             LR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
             (-((m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1)) +
              (m1s * (pap1 + 2.0 * papb) - 2.0 * pap1 * (pap1 + 2.0 * pbp1)) *
                  pbp2) *
                 RL * RLRL) +
        ml2 * (LRRR * m2 * (-(d0 * m1s) + 2.0 * (d0 + d2) * pap1) * papb +
               m1 * (-2.0 * (d1 + d2) * papb * pbp2 +
                     d0 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                           2.0 * papb * pbp2)) *
                   RLRR) *
            RR +
        d3 * (LR * LRLR * (-((m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1)) +
                           (m1s * (pap1 + 2.0 * papb) -
                            2.0 * pap1 * (pap1 + 2.0 * pbp1)) *
                               pbp2) +
              LRRL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL +
              LL * ml2 * (LRLL * m1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                                       2.0 * papb * pbp2) -
                          m2 * (m1s - 2.0 * pap1) * papb * RLLL) +
              LR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
              (-((m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1)) +
               (m1s * (pap1 + 2.0 * papb) - 2.0 * pap1 * (pap1 + 2.0 * pbp1)) *
                   pbp2) *
                  RL * RLRL +
              ml2 * (-(LRRR * m2 * (m1s - 2.0 * pap1) * papb) +
                     m1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                           2.0 * papb * pbp2) *
                         RLRR) *
                  RR))) /
      16.0;

  return real(me0);
}

double FI::Mtu_B2_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivu2s2 * MPIsI *
       (d33 * LR * LRRL * m1 * m1s * m2 * papb -
        d0 * LL * LRLL * m1 * ml2 * p1p2 * papb +
        2.0 * d0 * LL * LRLL * m1 * ml2 * pap2 * papb +
        2.0 * d2 * LL * LRLL * m1 * ml2 * pap2 * papb -
        2.0 * d33 * LR * LRRL * m1 * m2 * pap1 * pbp1 -
        d0 * LL * LRLL * m1 * ml2 * pap2 * pbp1 +
        d0 * LL * LRLL * m1 * ml2 * pap1 * pbp2 -
        d33 * LRLR * m1s * p1p2 * papb * RL +
        2.0 * d23 * LRLR * m1s * pap2 * papb * RL +
        2.0 * d33 * LRLR * m1s * pap2 * papb * RL +
        d33 * LRLR * m1s * pap2 * pbp1 * RL -
        4.0 * d23 * LRLR * pap1 * pap2 * pbp1 * RL -
        4.0 * d33 * LRLR * pap1 * pap2 * pbp1 * RL +
        2.0 * d13 * LRLR * p1p2 * papb * pbp1 * RL +
        2.0 * d23 * LRLR * p1p2 * papb * pbp1 * RL +
        2.0 * d33 * LRLR * p1p2 * papb * pbp1 * RL -
        2.0 * d13 * LRLR * pap2 * pow2(pbp1) * RL -
        2.0 * d23 * LRLR * pap2 * pow2(pbp1) * RL -
        2.0 * d33 * LRLR * pap2 * pow2(pbp1) * RL +
        d33 * LRLR * m1s * pap1 * pbp2 * RL -
        2.0 * d13 * LRLR * pap1 * pbp1 * pbp2 * RL -
        2.0 * d23 * LRLR * pap1 * pbp1 * pbp2 * RL -
        2.0 * d33 * LRLR * pap1 * pbp1 * pbp2 * RL +
        d0 * LL * m1s * m2 * ml2 * papb * RLLL -
        2.0 * d0 * LL * m2 * ml2 * papb * pbp1 * RLLL -
        2.0 * d1 * LL * m2 * ml2 * papb * pbp1 * RLLL -
        2.0 * d2 * LL * m2 * ml2 * papb * pbp1 * RLLL +
        d33 * m1 * m1s * m2 * papb * RL * RLLR -
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RL * RLLR -
        d33 * LR * m1s * p1p2 * papb * RLRL +
        2.0 * d23 * LR * m1s * pap2 * papb * RLRL +
        2.0 * d33 * LR * m1s * pap2 * papb * RLRL +
        d33 * LR * m1s * pap2 * pbp1 * RLRL -
        4.0 * d23 * LR * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d33 * LR * pap1 * pap2 * pbp1 * RLRL +
        2.0 * d13 * LR * p1p2 * papb * pbp1 * RLRL +
        2.0 * d23 * LR * p1p2 * papb * pbp1 * RLRL +
        2.0 * d33 * LR * p1p2 * papb * pbp1 * RLRL -
        2.0 * d13 * LR * pap2 * pow2(pbp1) * RLRL -
        2.0 * d23 * LR * pap2 * pow2(pbp1) * RLRL -
        2.0 * d33 * LR * pap2 * pow2(pbp1) * RLRL +
        d33 * LR * m1s * pap1 * pbp2 * RLRL -
        2.0 * d13 * LR * pap1 * pbp1 * pbp2 * RLRL -
        2.0 * d23 * LR * pap1 * pbp1 * pbp2 * RLRL -
        2.0 * d33 * LR * pap1 * pbp1 * pbp2 * RLRL +
        2.0 * d00 * (2.0 * LRLR * pap2 * pbp1 * RL +
                     m1 * m2 * papb * (LR * LRRL + RL * RLLR) +
                     2.0 * LR * pap2 * pbp1 * RLRL) +
        ml2 * (LRRR * m2 * papb * (d0 * m1s - 2.0 * (d0 + d1 + d2) * pbp1) +
               m1 * (2.0 * d2 * pap2 * papb +
                     d0 * (-(p1p2 * papb) + 2.0 * pap2 * papb - pap2 * pbp1 +
                           pap1 * pbp2)) *
                   RLRR) *
            RR +
        d3 * (LR * LRRL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) +
              LRLR * (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb +
                                     pap2 * pbp1 + pap1 * pbp2) +
                      m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb + pap2 * pbp1 +
                             pap1 * pbp2)) *
                  RL +
              LL * ml2 * (LRLL * m1 * (-(p1p2 * papb) + 2.0 * pap2 * papb -
                                       pap2 * pbp1 + pap1 * pbp2) +
                          m2 * papb * (m1s - 2.0 * pbp1) * RLLL) +
              m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL * RLLR +
              LR * (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb +
                                   pap2 * pbp1 + pap1 * pbp2) +
                    m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb + pap2 * pbp1 +
                           pap1 * pbp2)) *
                  RLRL +
              ml2 * (LRRR * m2 * papb * (m1s - 2.0 * pbp1) +
                     m1 * (-(p1p2 * papb) + 2.0 * pap2 * papb - pap2 * pbp1 +
                           pap1 * pbp2) *
                         RLRR) *
                  RR))) /
      16.0;

  return real(me0);
}

double FI::Muu_B2_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (-3.0 * (cA - 2.0 * cF) * cF * ivu2s2 * MPIsI *
       (d0 * LL * LRLL * m1 * ml2 * p1p2 * papb +
        2.0 * d13 * LR * LRLR * m1s * pap2 * papb +
        2.0 * d23 * LR * LRLR * m1s * pap2 * papb -
        2.0 * d0 * LL * LRLL * m1 * ml2 * pap2 * papb -
        2.0 * d1 * LL * LRLL * m1 * ml2 * pap2 * papb -
        2.0 * d2 * LL * LRLL * m1 * ml2 * pap2 * papb +
        4.0 * d00 * LR * LRLR * pap2 * pbp1 +
        d0 * LL * LRLL * m1 * ml2 * pap2 * pbp1 -
        4.0 * d13 * LR * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d23 * LR * LRLR * pap1 * pap2 * pbp1 +
        2.0 * d23 * LR * LRLR * p1p2 * papb * pbp1 -
        2.0 * d23 * LR * LRLR * pap2 * pow2(pbp1) -
        d0 * LL * LRLL * m1 * ml2 * pap1 * pbp2 -
        2.0 * d23 * LR * LRLR * pap1 * pbp1 * pbp2 +
        2.0 * d00 * LRRL * m1 * m2 * papb * RL -
        d0 * LL * m1s * m2 * ml2 * papb * RLLL +
        2.0 * d0 * LL * m2 * ml2 * papb * pbp1 * RLLL +
        2.0 * d2 * LL * m2 * ml2 * papb * pbp1 * RLLL +
        2.0 * d00 * LR * m1 * m2 * papb * RLLR +
        2.0 * d13 * m1s * pap2 * papb * RL * RLRL +
        2.0 * d23 * m1s * pap2 * papb * RL * RLRL +
        4.0 * d00 * pap2 * pbp1 * RL * RLRL -
        4.0 * d13 * pap1 * pap2 * pbp1 * RL * RLRL -
        4.0 * d23 * pap1 * pap2 * pbp1 * RL * RLRL +
        2.0 * d23 * p1p2 * papb * pbp1 * RL * RLRL -
        2.0 * d23 * pap2 * pow2(pbp1) * RL * RLRL -
        2.0 * d23 * pap1 * pbp1 * pbp2 * RL * RLRL +
        d33 * (LR * LRLR * (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb +
                                           pap2 * pbp1 + pap1 * pbp2) +
                            m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb +
                                   pap2 * pbp1 + pap1 * pbp2)) +
               LRRL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL +
               LR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
               (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb + pap2 * pbp1 +
                               pap1 * pbp2) +
                m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb + pap2 * pbp1 +
                       pap1 * pbp2)) *
                   RL * RLRL) +
        ml2 * (LRRR * m2 * papb * (-(d0 * m1s) + 2.0 * (d0 + d2) * pbp1) +
               m1 * (-2.0 * (d1 + d2) * pap2 * papb +
                     d0 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2)) *
                   RLRR) *
            RR +
        d3 * (LR * LRLR * (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb +
                                          pap2 * pbp1 + pap1 * pbp2) +
                           m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb +
                                  pap2 * pbp1 + pap1 * pbp2)) +
              LRRL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL +
              LL * ml2 * (LRLL * m1 * (p1p2 * papb - 2.0 * pap2 * papb +
                                       pap2 * pbp1 - pap1 * pbp2) -
                          m2 * papb * (m1s - 2.0 * pbp1) * RLLL) +
              LR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
              (-2.0 * pbp1 * (2.0 * pap1 * pap2 - p1p2 * papb + pap2 * pbp1 +
                              pap1 * pbp2) +
               m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb + pap2 * pbp1 +
                      pap1 * pbp2)) *
                  RL * RLRL +
              ml2 * (-(LRRR * m2 * papb * (m1s - 2.0 * pbp1)) +
                     m1 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2) *
                         RLRR) *
                  RR))) /
      16.0;

  return real(me0);
}

double FI::Mut_B2_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (3.0 * (cA - 2.0 * cF) * cF * ivt1s2 * MPIsI *
       (d3 * LR * LRRL * m1 * m1s * m2 * papb +
        d33 * LR * LRRL * m1 * m1s * m2 * papb -
        d0 * LL * LRLL * m1 * ml2 * p1p2 * papb -
        d3 * LL * LRLL * m1 * ml2 * p1p2 * papb -
        2.0 * d3 * LR * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d33 * LR * LRRL * m1 * m2 * pap1 * pbp1 +
        d0 * LL * LRLL * m1 * ml2 * pap2 * pbp1 +
        d3 * LL * LRLL * m1 * ml2 * pap2 * pbp1 -
        d0 * LL * LRLL * m1 * ml2 * pap1 * pbp2 -
        d3 * LL * LRLL * m1 * ml2 * pap1 * pbp2 +
        2.0 * d0 * LL * LRLL * m1 * ml2 * papb * pbp2 +
        2.0 * d2 * LL * LRLL * m1 * ml2 * papb * pbp2 +
        2.0 * d3 * LL * LRLL * m1 * ml2 * papb * pbp2 -
        d3 * LRLR * m1s * p1p2 * papb * RL -
        d33 * LRLR * m1s * p1p2 * papb * RL +
        2.0 * d13 * LRLR * p1p2 * pap1 * papb * RL +
        2.0 * d23 * LRLR * p1p2 * pap1 * papb * RL +
        2.0 * d3 * LRLR * p1p2 * pap1 * papb * RL +
        2.0 * d33 * LRLR * p1p2 * pap1 * papb * RL +
        d3 * LRLR * m1s * pap2 * pbp1 * RL +
        d33 * LRLR * m1s * pap2 * pbp1 * RL -
        2.0 * d13 * LRLR * pap1 * pap2 * pbp1 * RL -
        2.0 * d23 * LRLR * pap1 * pap2 * pbp1 * RL -
        2.0 * d3 * LRLR * pap1 * pap2 * pbp1 * RL -
        2.0 * d33 * LRLR * pap1 * pap2 * pbp1 * RL +
        d3 * LRLR * m1s * pap1 * pbp2 * RL +
        d33 * LRLR * m1s * pap1 * pbp2 * RL -
        2.0 * d13 * LRLR * pow2(pap1) * pbp2 * RL -
        2.0 * d23 * LRLR * pow2(pap1) * pbp2 * RL -
        2.0 * d3 * LRLR * pow2(pap1) * pbp2 * RL -
        2.0 * d33 * LRLR * pow2(pap1) * pbp2 * RL +
        2.0 * d23 * LRLR * m1s * papb * pbp2 * RL +
        2.0 * d3 * LRLR * m1s * papb * pbp2 * RL +
        2.0 * d33 * LRLR * m1s * papb * pbp2 * RL -
        4.0 * d23 * LRLR * pap1 * pbp1 * pbp2 * RL -
        4.0 * d3 * LRLR * pap1 * pbp1 * pbp2 * RL -
        4.0 * d33 * LRLR * pap1 * pbp1 * pbp2 * RL +
        d0 * LL * m1s * m2 * ml2 * papb * RLLL +
        d3 * LL * m1s * m2 * ml2 * papb * RLLL -
        2.0 * d0 * LL * m2 * ml2 * pap1 * papb * RLLL -
        2.0 * d1 * LL * m2 * ml2 * pap1 * papb * RLLL -
        2.0 * d2 * LL * m2 * ml2 * pap1 * papb * RLLL -
        2.0 * d3 * LL * m2 * ml2 * pap1 * papb * RLLL +
        d3 * m1 * m1s * m2 * papb * RL * RLLR +
        d33 * m1 * m1s * m2 * papb * RL * RLLR -
        2.0 * d3 * m1 * m2 * pap1 * pbp1 * RL * RLLR -
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RL * RLLR -
        d3 * LR * m1s * p1p2 * papb * RLRL -
        d33 * LR * m1s * p1p2 * papb * RLRL +
        2.0 * d13 * LR * p1p2 * pap1 * papb * RLRL +
        2.0 * d23 * LR * p1p2 * pap1 * papb * RLRL +
        2.0 * d3 * LR * p1p2 * pap1 * papb * RLRL +
        2.0 * d33 * LR * p1p2 * pap1 * papb * RLRL +
        d3 * LR * m1s * pap2 * pbp1 * RLRL +
        d33 * LR * m1s * pap2 * pbp1 * RLRL -
        2.0 * d13 * LR * pap1 * pap2 * pbp1 * RLRL -
        2.0 * d23 * LR * pap1 * pap2 * pbp1 * RLRL -
        2.0 * d3 * LR * pap1 * pap2 * pbp1 * RLRL -
        2.0 * d33 * LR * pap1 * pap2 * pbp1 * RLRL +
        d3 * LR * m1s * pap1 * pbp2 * RLRL +
        d33 * LR * m1s * pap1 * pbp2 * RLRL -
        2.0 * d13 * LR * pow2(pap1) * pbp2 * RLRL -
        2.0 * d23 * LR * pow2(pap1) * pbp2 * RLRL -
        2.0 * d3 * LR * pow2(pap1) * pbp2 * RLRL -
        2.0 * d33 * LR * pow2(pap1) * pbp2 * RLRL +
        2.0 * d23 * LR * m1s * papb * pbp2 * RLRL +
        2.0 * d3 * LR * m1s * papb * pbp2 * RLRL +
        2.0 * d33 * LR * m1s * papb * pbp2 * RLRL -
        4.0 * d23 * LR * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d3 * LR * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d33 * LR * pap1 * pbp1 * pbp2 * RLRL +
        2.0 * d00 * (2.0 * LRLR * pap1 * pbp2 * RL +
                     m1 * m2 * papb * (LR * LRRL + RL * RLLR) +
                     2.0 * LR * pap1 * pbp2 * RLRL) +
        ml2 *
            (LRRR * m2 * ((d0 + d3) * m1s - 2.0 * (d0 + d1 + d2 + d3) * pap1) *
                 papb +
             m1 * (-((d0 + d3) * (p1p2 * papb - pap2 * pbp1)) +
                   (-((d0 + d3) * pap1) + 2.0 * (d0 + d2 + d3) * papb) * pbp2) *
                 RLRR) *
            RR)) /
      16.0;

  return real(me0);
}

// Box3
/************************************************/
double FI::Mtt_B3_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (cA * cF * ivt1s2 * MPIsI * NC *
       (-(d22 * LLLL * m1s * p1p2 * papb) -
        2.0 * d23 * LLLL * m1s * p1p2 * papb - d3 * LLLL * m1s * p1p2 * papb -
        d33 * LLLL * m1s * p1p2 * papb - d22 * LRLR * m1s * p1p2 * papb -
        2.0 * d23 * LRLR * m1s * p1p2 * papb - d3 * LRLR * m1s * p1p2 * papb -
        d33 * LRLR * m1s * p1p2 * papb + d3 * LLLL * m1 * ml3 * p1p2 * papb +
        d3 * LRLR * m1 * ml3 * p1p2 * papb +
        2.0 * d23 * LLLL * p1p2 * pap1 * papb +
        2.0 * d3 * LLLL * p1p2 * pap1 * papb +
        2.0 * d33 * LLLL * p1p2 * pap1 * papb +
        2.0 * d23 * LRLR * p1p2 * pap1 * papb +
        2.0 * d3 * LRLR * p1p2 * pap1 * papb +
        2.0 * d33 * LRLR * p1p2 * pap1 * papb -
        2.0 * d13 * LLLL * p1p2 * pow2(papb) -
        2.0 * d23 * LLLL * p1p2 * pow2(papb) -
        2.0 * d3 * LLLL * p1p2 * pow2(papb) -
        2.0 * d33 * LLLL * p1p2 * pow2(papb) + d22 * LLLL * m1s * pap2 * pbp1 +
        2.0 * d23 * LLLL * m1s * pap2 * pbp1 + d3 * LLLL * m1s * pap2 * pbp1 +
        d33 * LLLL * m1s * pap2 * pbp1 + d22 * LRLR * m1s * pap2 * pbp1 +
        2.0 * d23 * LRLR * m1s * pap2 * pbp1 + d3 * LRLR * m1s * pap2 * pbp1 +
        d33 * LRLR * m1s * pap2 * pbp1 - d3 * LLLL * m1 * ml3 * pap2 * pbp1 -
        d3 * LRLR * m1 * ml3 * pap2 * pbp1 -
        2.0 * d23 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d3 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d33 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d23 * LRLR * pap1 * pap2 * pbp1 -
        2.0 * d3 * LRLR * pap1 * pap2 * pbp1 -
        2.0 * d33 * LRLR * pap1 * pap2 * pbp1 +
        2.0 * d12 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d13 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d22 * LLLL * p1p2 * papb * pbp1 +
        4.0 * d23 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d3 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d33 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d13 * LLLL * pap2 * papb * pbp1 +
        2.0 * d23 * LLLL * pap2 * papb * pbp1 +
        2.0 * d3 * LLLL * pap2 * papb * pbp1 +
        2.0 * d33 * LLLL * pap2 * papb * pbp1 -
        2.0 * d12 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d13 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d22 * LLLL * pap2 * pow2(pbp1) -
        4.0 * d23 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d3 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d33 * LLLL * pap2 * pow2(pbp1) + d22 * LLLL * m1s * pap1 * pbp2 +
        2.0 * d23 * LLLL * m1s * pap1 * pbp2 + d3 * LLLL * m1s * pap1 * pbp2 +
        d33 * LLLL * m1s * pap1 * pbp2 + d22 * LRLR * m1s * pap1 * pbp2 +
        2.0 * d23 * LRLR * m1s * pap1 * pbp2 + d3 * LRLR * m1s * pap1 * pbp2 +
        d33 * LRLR * m1s * pap1 * pbp2 + d3 * LLLL * m1 * ml3 * pap1 * pbp2 +
        d3 * LRLR * m1 * ml3 * pap1 * pbp2 -
        2.0 * d23 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d3 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d33 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d23 * LRLR * pow2(pap1) * pbp2 -
        2.0 * d3 * LRLR * pow2(pap1) * pbp2 -
        2.0 * d33 * LRLR * pow2(pap1) * pbp2 +
        2.0 * d1 * LLLL * m1s * papb * pbp2 +
        2.0 * d1 * LRLR * m1s * papb * pbp2 +
        2.0 * d12 * LRLR * m1s * papb * pbp2 +
        2.0 * d13 * LRLR * m1s * papb * pbp2 +
        2.0 * d22 * LRLR * m1s * papb * pbp2 +
        4.0 * d23 * LRLR * m1s * papb * pbp2 +
        2.0 * d3 * LRLR * m1s * papb * pbp2 +
        2.0 * d33 * LRLR * m1s * papb * pbp2 -
        2.0 * d1 * LLLL * m1 * ml3 * papb * pbp2 -
        2.0 * d3 * LLLL * m1 * ml3 * papb * pbp2 -
        2.0 * d1 * LRLR * m1 * ml3 * papb * pbp2 -
        2.0 * d3 * LRLR * m1 * ml3 * papb * pbp2 +
        2.0 * d13 * LLLL * pap1 * papb * pbp2 +
        2.0 * d23 * LLLL * pap1 * papb * pbp2 +
        2.0 * d3 * LLLL * pap1 * papb * pbp2 +
        2.0 * d33 * LLLL * pap1 * papb * pbp2 -
        4.0 * d1 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d12 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d13 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d22 * LLLL * pap1 * pbp1 * pbp2 -
        4.0 * d23 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d3 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d33 * LLLL * pap1 * pbp1 * pbp2 -
        4.0 * d1 * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d12 * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d13 * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d22 * LRLR * pap1 * pbp1 * pbp2 -
        8.0 * d23 * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d3 * LRLR * pap1 * pbp1 * pbp2 -
        4.0 * d33 * LRLR * pap1 * pbp1 * pbp2 - d22 * m1s * p1p2 * papb * RLRL -
        2.0 * d23 * m1s * p1p2 * papb * RLRL - d3 * m1s * p1p2 * papb * RLRL -
        d33 * m1s * p1p2 * papb * RLRL + d3 * m1 * ml3 * p1p2 * papb * RLRL +
        2.0 * d23 * p1p2 * pap1 * papb * RLRL +
        2.0 * d3 * p1p2 * pap1 * papb * RLRL +
        2.0 * d33 * p1p2 * pap1 * papb * RLRL + d22 * m1s * pap2 * pbp1 * RLRL +
        2.0 * d23 * m1s * pap2 * pbp1 * RLRL + d3 * m1s * pap2 * pbp1 * RLRL +
        d33 * m1s * pap2 * pbp1 * RLRL - d3 * m1 * ml3 * pap2 * pbp1 * RLRL -
        2.0 * d23 * pap1 * pap2 * pbp1 * RLRL -
        2.0 * d3 * pap1 * pap2 * pbp1 * RLRL -
        2.0 * d33 * pap1 * pap2 * pbp1 * RLRL + d22 * m1s * pap1 * pbp2 * RLRL +
        2.0 * d23 * m1s * pap1 * pbp2 * RLRL + d3 * m1s * pap1 * pbp2 * RLRL +
        d33 * m1s * pap1 * pbp2 * RLRL + d3 * m1 * ml3 * pap1 * pbp2 * RLRL -
        2.0 * d23 * pow2(pap1) * pbp2 * RLRL -
        2.0 * d3 * pow2(pap1) * pbp2 * RLRL -
        2.0 * d33 * pow2(pap1) * pbp2 * RLRL +
        2.0 * d1 * m1s * papb * pbp2 * RLRL +
        2.0 * d12 * m1s * papb * pbp2 * RLRL +
        2.0 * d13 * m1s * papb * pbp2 * RLRL +
        2.0 * d22 * m1s * papb * pbp2 * RLRL +
        4.0 * d23 * m1s * papb * pbp2 * RLRL +
        2.0 * d3 * m1s * papb * pbp2 * RLRL +
        2.0 * d33 * m1s * papb * pbp2 * RLRL -
        2.0 * d1 * m1 * ml3 * papb * pbp2 * RLRL -
        2.0 * d3 * m1 * ml3 * papb * pbp2 * RLRL -
        4.0 * d1 * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d12 * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d13 * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d22 * pap1 * pbp1 * pbp2 * RLRL -
        8.0 * d23 * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d3 * pap1 * pbp1 * pbp2 * RLRL -
        4.0 * d33 * pap1 * pbp1 * pbp2 * RLRL +
        (-(d3 * m1s * p1p2 * papb) - d33 * m1s * p1p2 * papb +
         d3 * m1 * ml3 * p1p2 * papb + 2.0 * d3 * p1p2 * pap1 * papb +
         2.0 * d33 * p1p2 * pap1 * papb - 2.0 * d13 * p1p2 * pow2(papb) -
         2.0 * d3 * p1p2 * pow2(papb) - 2.0 * d33 * p1p2 * pow2(papb) +
         d3 * m1s * pap2 * pbp1 + d33 * m1s * pap2 * pbp1 -
         d3 * m1 * ml3 * pap2 * pbp1 - 2.0 * d3 * pap1 * pap2 * pbp1 -
         2.0 * d33 * pap1 * pap2 * pbp1 + 2.0 * d12 * p1p2 * papb * pbp1 +
         2.0 * d13 * p1p2 * papb * pbp1 + 2.0 * d3 * p1p2 * papb * pbp1 +
         2.0 * d33 * p1p2 * papb * pbp1 + 2.0 * d13 * pap2 * papb * pbp1 +
         2.0 * d3 * pap2 * papb * pbp1 + 2.0 * d33 * pap2 * papb * pbp1 -
         2.0 * d12 * pap2 * pow2(pbp1) - 2.0 * d13 * pap2 * pow2(pbp1) -
         2.0 * d3 * pap2 * pow2(pbp1) - 2.0 * d33 * pap2 * pow2(pbp1) +
         (2.0 * (d1 * (m1s - m1 * ml3) + d13 * pap1) * papb -
          2.0 * (2.0 * d1 + d12 + d13) * pap1 * pbp1 +
          d33 * pap1 * (m1s - 2.0 * (pap1 - papb + pbp1)) +
          d3 * (m1s * pap1 + m1 * ml3 * (pap1 - 2.0 * papb) -
                2.0 * pap1 * (pap1 - papb + pbp1))) *
             pbp2 +
         d22 * (m1s - 2.0 * pbp1) *
             (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
         2.0 * d23 * (m1s - pap1 + papb - 2.0 * pbp1) *
             (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) *
            RRRR +
        d2 * (LLLL * (m1s - m1 * ml3 - 2.0 * pbp1) *
                  (-(p1p2 * papb) + pap2 * pbp1) +
              LLLL * (m1s * pap1 + m1 * ml3 * (pap1 - 2.0 * papb) -
                      2.0 * pap1 * pbp1) *
                  pbp2 -
              ((m1s - m1 * ml3) * (p1p2 * papb - pap2 * pbp1) -
               (m1 * ml3 * (pap1 - 2.0 * papb) + m1s * (pap1 + 2.0 * papb) -
                4.0 * pap1 * pbp1) *
                   pbp2) *
                  (LRLR + RLRL) +
              ((m1s - m1 * ml3 - 2.0 * pbp1) * (-(p1p2 * papb) + pap2 * pbp1) +
               (m1s * pap1 + m1 * ml3 * (pap1 - 2.0 * papb) -
                2.0 * pap1 * pbp1) *
                   pbp2) *
                  RRRR) +
        4.0 * d00 * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) -
                     p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
                     pap1 * pbp2 * (LRLR + RLRL + RRRR)))) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me1 =
      (cA * cF * ivt1s2 * MPIsI * NC *
       ((((d13 + d23 + d3 + d33) * papb -
          (d12 + d13 + d2 + d22 + 2.0 * d23 + d3 + d33) * pbp1) *
             (p1p2 * papb - pap2 * pbp1) +
         (((d12 + d13 + d2 + d22 + 2.0 * d23 + d3 + d33) * m1s -
           (d13 + d23 + d3 + d33) * pap1) *
              papb -
          (d12 + d13 + d2 + d22 + 2.0 * d23 + d3 + d33) * pap1 * pbp1) *
             pbp2) *
            (LLLL - LRLR - RLRL + RRRR) -
        d00 *
            (LLLL * (-3.0 * p1p2 * papb + 3.0 * pap2 * pbp1 + pap1 * pbp2) +
             (3.0 * p1p2 * papb - 3.0 * pap2 * pbp1 + pap1 * pbp2) *
                 (LRLR + RLRL) +
             (-3.0 * p1p2 * papb + 3.0 * pap2 * pbp1 + pap1 * pbp2) * RRRR))) /
      4.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me2 =
      -(cA * cF * d00 * ivt1s2 * MPIsI * NC * (p1p2 * papb - pap2 * pbp1) *
        (LLLL - LRLR - RLRL + RRRR)) /
      2.0;

  return real(me0 + me1 + me2);
}

double FI::Mtu_B3_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (cA * cF * ivu2s2 * MPIsI * NC *
       (-(d22 * LRRL * m1 * m1s * m2 * papb) -
        2.0 * d23 * LRRL * m1 * m1s * m2 * papb -
        d3 * LRRL * m1 * m1s * m2 * papb - d33 * LRRL * m1 * m1s * m2 * papb +
        d3 * LLLL * m1 * ml3 * p1p2 * papb -
        2.0 * d3 * LLLL * m1 * ml3 * pap2 * papb +
        2.0 * d22 * LRRL * m1 * m2 * pap1 * pbp1 +
        4.0 * d23 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d3 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d33 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d3 * LRRL * m2 * ml3 * pap1 * pbp1 -
        2.0 * d22 * LLLL * m1s * pap2 * pbp1 -
        4.0 * d23 * LLLL * m1s * pap2 * pbp1 -
        2.0 * d3 * LLLL * m1s * pap2 * pbp1 -
        2.0 * d33 * LLLL * m1s * pap2 * pbp1 +
        d3 * LLLL * m1 * ml3 * pap2 * pbp1 +
        4.0 * d23 * LLLL * pap1 * pap2 * pbp1 +
        4.0 * d3 * LLLL * pap1 * pap2 * pbp1 +
        4.0 * d33 * LLLL * pap1 * pap2 * pbp1 +
        2.0 * d1 * LRRL * m2 * ml3 * papb * pbp1 +
        2.0 * d3 * LRRL * m2 * ml3 * papb * pbp1 -
        2.0 * d1 * LLLL * p1p2 * papb * pbp1 -
        4.0 * d13 * LLLL * pap2 * papb * pbp1 -
        4.0 * d23 * LLLL * pap2 * papb * pbp1 -
        4.0 * d3 * LLLL * pap2 * papb * pbp1 -
        4.0 * d33 * LLLL * pap2 * papb * pbp1 +
        2.0 * d1 * LLLL * pap2 * pow2(pbp1) +
        4.0 * d12 * LLLL * pap2 * pow2(pbp1) +
        4.0 * d13 * LLLL * pap2 * pow2(pbp1) +
        4.0 * d22 * LLLL * pap2 * pow2(pbp1) +
        8.0 * d23 * LLLL * pap2 * pow2(pbp1) +
        4.0 * d3 * LLLL * pap2 * pow2(pbp1) +
        4.0 * d33 * LLLL * pap2 * pow2(pbp1) -
        d3 * LLLL * m1 * ml3 * pap1 * pbp2 +
        2.0 * d1 * LLLL * pap1 * pbp1 * pbp2 -
        d22 * m1 * m1s * m2 * papb * RLLR -
        2.0 * d23 * m1 * m1s * m2 * papb * RLLR -
        d3 * m1 * m1s * m2 * papb * RLLR - d33 * m1 * m1s * m2 * papb * RLLR +
        2.0 * d22 * m1 * m2 * pap1 * pbp1 * RLLR +
        4.0 * d23 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d3 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d3 * m2 * ml3 * pap1 * pbp1 * RLLR +
        2.0 * d1 * m2 * ml3 * papb * pbp1 * RLLR +
        2.0 * d3 * m2 * ml3 * papb * pbp1 * RLLR +
        (2.0 * pbp1 *
             (-(d1 * p1p2 * papb) +
              pap2 * (-(d22 * m1s) - 2.0 * d23 * m1s - d33 * m1s +
                      2.0 * d23 * pap1 + 2.0 * d33 * pap1 -
                      2.0 * (d13 + d23 + d33) * papb +
                      (d1 + 2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33)) * pbp1) +
              d1 * pap1 * pbp2) +
         d3 * (2.0 * pap2 * pbp1 * (-m1s + 2.0 * (pap1 - papb + pbp1)) +
               m1 * ml3 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2))) *
            RRRR +
        d2 * (LRRL * m2 * (-(m1 * m1s * papb) +
                           2.0 * (m1 * pap1 + ml3 * (-pap1 + papb)) * pbp1) +
              LLLL * (2.0 * pap2 * pbp1 * (-m1s + 2.0 * pbp1) +
                      m1 * ml3 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) +
              m2 * (-(m1 * m1s * papb) +
                    2.0 * (m1 * pap1 + ml3 * (-pap1 + papb)) * pbp1) *
                  RLLR +
              (2.0 * pap2 * pbp1 * (-m1s + 2.0 * pbp1) +
               m1 * ml3 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
                  RRRR) -
        2.0 * d00 * (m1 * m2 * papb * (LRRL + RLLR) +
                     4.0 * pap2 * pbp1 * (LLLL + RRRR)))) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me1 =
      (cA * cF * ivu2s2 * MPIsI * NC *
       (d22 * LRRL * m1 * m1s * m2 * papb +
        2.0 * d23 * LRRL * m1 * m1s * m2 * papb +
        d3 * LRRL * m1 * m1s * m2 * papb + d33 * LRRL * m1 * m1s * m2 * papb -
        d3 * LRRL * m1s * m2 * ml3 * papb + d22 * LLLL * m1s * p1p2 * papb +
        2.0 * d23 * LLLL * m1s * p1p2 * papb + d3 * LLLL * m1s * p1p2 * papb +
        d33 * LLLL * m1s * p1p2 * papb - d3 * LLLL * m1 * ml3 * p1p2 * papb -
        2.0 * d23 * LLLL * m1s * pap2 * papb -
        2.0 * d3 * LLLL * m1s * pap2 * papb -
        2.0 * d33 * LLLL * m1s * pap2 * papb +
        2.0 * d3 * LLLL * m1 * ml3 * pap2 * papb -
        2.0 * d22 * LRRL * m1 * m2 * pap1 * pbp1 -
        4.0 * d23 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d3 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d33 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d3 * LRRL * m2 * ml3 * pap1 * pbp1 +
        d22 * LLLL * m1s * pap2 * pbp1 + 2.0 * d23 * LLLL * m1s * pap2 * pbp1 +
        d3 * LLLL * m1s * pap2 * pbp1 + d33 * LLLL * m1s * pap2 * pbp1 -
        d3 * LLLL * m1 * ml3 * pap2 * pbp1 -
        2.0 * d12 * LLLL * p1p2 * papb * pbp1 -
        2.0 * d13 * LLLL * p1p2 * papb * pbp1 -
        2.0 * d22 * LLLL * p1p2 * papb * pbp1 -
        4.0 * d23 * LLLL * p1p2 * papb * pbp1 -
        2.0 * d3 * LLLL * p1p2 * papb * pbp1 -
        2.0 * d33 * LLLL * p1p2 * papb * pbp1 +
        4.0 * d13 * LLLL * pap2 * papb * pbp1 +
        4.0 * d23 * LLLL * pap2 * papb * pbp1 +
        4.0 * d3 * LLLL * pap2 * papb * pbp1 +
        4.0 * d33 * LLLL * pap2 * papb * pbp1 -
        2.0 * d12 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d13 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d22 * LLLL * pap2 * pow2(pbp1) -
        4.0 * d23 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d3 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d33 * LLLL * pap2 * pow2(pbp1) - d22 * LLLL * m1s * pap1 * pbp2 -
        2.0 * d23 * LLLL * m1s * pap1 * pbp2 - d3 * LLLL * m1s * pap1 * pbp2 -
        d33 * LLLL * m1s * pap1 * pbp2 + d3 * LLLL * m1 * ml3 * pap1 * pbp2 +
        2.0 * d12 * LLLL * pap1 * pbp1 * pbp2 +
        2.0 * d13 * LLLL * pap1 * pbp1 * pbp2 +
        2.0 * d22 * LLLL * pap1 * pbp1 * pbp2 +
        4.0 * d23 * LLLL * pap1 * pbp1 * pbp2 +
        2.0 * d3 * LLLL * pap1 * pbp1 * pbp2 +
        2.0 * d33 * LLLL * pap1 * pbp1 * pbp2 +
        d22 * m1 * m1s * m2 * papb * RLLR +
        2.0 * d23 * m1 * m1s * m2 * papb * RLLR +
        d3 * m1 * m1s * m2 * papb * RLLR + d33 * m1 * m1s * m2 * papb * RLLR -
        d3 * m1s * m2 * ml3 * papb * RLLR -
        2.0 * d22 * m1 * m2 * pap1 * pbp1 * RLLR -
        4.0 * d23 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d3 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d3 * m2 * ml3 * pap1 * pbp1 * RLLR +
        (d3 * m1s * p1p2 * papb + d33 * m1s * p1p2 * papb -
         d3 * m1 * ml3 * p1p2 * papb - 2.0 * d3 * m1s * pap2 * papb -
         2.0 * d33 * m1s * pap2 * papb + 2.0 * d3 * m1 * ml3 * pap2 * papb +
         d3 * m1s * pap2 * pbp1 + d33 * m1s * pap2 * pbp1 -
         d3 * m1 * ml3 * pap2 * pbp1 - 2.0 * d12 * p1p2 * papb * pbp1 -
         2.0 * d13 * p1p2 * papb * pbp1 - 2.0 * d3 * p1p2 * papb * pbp1 -
         2.0 * d33 * p1p2 * papb * pbp1 + 4.0 * d13 * pap2 * papb * pbp1 +
         4.0 * d3 * pap2 * papb * pbp1 + 4.0 * d33 * pap2 * papb * pbp1 -
         2.0 * d12 * pap2 * pow2(pbp1) - 2.0 * d13 * pap2 * pow2(pbp1) -
         2.0 * d3 * pap2 * pow2(pbp1) - 2.0 * d33 * pap2 * pow2(pbp1) +
         pap1 * (-(d3 * m1s) - d33 * m1s + d3 * m1 * ml3 +
                 2.0 * (d12 + d13 + d3 + d33) * pbp1) *
             pbp2 +
         d22 * (m1s - 2.0 * pbp1) * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
         2.0 * d23 * (m1s - 2.0 * pbp1) *
             (p1p2 * papb - pap2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
            RRRR +
        d2 * (LRRL * m2 * (m1 - ml3) * (m1s * papb - 2.0 * pap1 * pbp1) +
              LLLL * (m1s - m1 * ml3 - 2.0 * pbp1) *
                  (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
              m2 * (m1 - ml3) * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
              (m1s - m1 * ml3 - 2.0 * pbp1) *
                  (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR) +
        4.0 * d00 * (m1 * m2 * papb * (LRRL + RLLR) +
                     2.0 * pap2 * pbp1 * (LLLL + RRRR)))) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me2 = -(cA * cF * d00 * ivu2s2 * MPIsI * NC *
                          (LRRL * m1 * m2 * papb +
                           LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                           m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                           pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        4.0;

  return real(me0 + me1 + me2);
}

double FI::Muu_B3_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (cA * cF * ivu2s2 * MPIsI * NC *
       (-(d22 * LLLL * m1s * p1p2 * papb) -
        2.0 * d23 * LLLL * m1s * p1p2 * papb - d3 * LLLL * m1s * p1p2 * papb -
        d33 * LLLL * m1s * p1p2 * papb - d22 * LRLR * m1s * p1p2 * papb -
        2.0 * d23 * LRLR * m1s * p1p2 * papb - d3 * LRLR * m1s * p1p2 * papb -
        d33 * LRLR * m1s * p1p2 * papb + d3 * LLLL * m1 * ml3 * p1p2 * papb +
        d3 * LRLR * m1 * ml3 * p1p2 * papb +
        2.0 * d12 * LLLL * p1p2 * pap1 * papb +
        2.0 * d13 * LLLL * p1p2 * pap1 * papb +
        2.0 * d22 * LLLL * p1p2 * pap1 * papb +
        4.0 * d23 * LLLL * p1p2 * pap1 * papb +
        2.0 * d3 * LLLL * p1p2 * pap1 * papb +
        2.0 * d33 * LLLL * p1p2 * pap1 * papb +
        2.0 * d1 * LLLL * m1s * pap2 * papb +
        2.0 * d1 * LRLR * m1s * pap2 * papb +
        2.0 * d12 * LRLR * m1s * pap2 * papb +
        2.0 * d13 * LRLR * m1s * pap2 * papb +
        2.0 * d22 * LRLR * m1s * pap2 * papb +
        4.0 * d23 * LRLR * m1s * pap2 * papb +
        2.0 * d3 * LRLR * m1s * pap2 * papb +
        2.0 * d33 * LRLR * m1s * pap2 * papb -
        2.0 * d1 * LLLL * m1 * ml3 * pap2 * papb -
        2.0 * d3 * LLLL * m1 * ml3 * pap2 * papb -
        2.0 * d1 * LRLR * m1 * ml3 * pap2 * papb -
        2.0 * d3 * LRLR * m1 * ml3 * pap2 * papb -
        2.0 * d13 * LLLL * p1p2 * pow2(papb) -
        2.0 * d23 * LLLL * p1p2 * pow2(papb) -
        2.0 * d3 * LLLL * p1p2 * pow2(papb) -
        2.0 * d33 * LLLL * p1p2 * pow2(papb) + d22 * LLLL * m1s * pap2 * pbp1 +
        2.0 * d23 * LLLL * m1s * pap2 * pbp1 + d3 * LLLL * m1s * pap2 * pbp1 +
        d33 * LLLL * m1s * pap2 * pbp1 + d22 * LRLR * m1s * pap2 * pbp1 +
        2.0 * d23 * LRLR * m1s * pap2 * pbp1 + d3 * LRLR * m1s * pap2 * pbp1 +
        d33 * LRLR * m1s * pap2 * pbp1 + d3 * LLLL * m1 * ml3 * pap2 * pbp1 +
        d3 * LRLR * m1 * ml3 * pap2 * pbp1 -
        4.0 * d1 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d12 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d13 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d22 * LLLL * pap1 * pap2 * pbp1 -
        4.0 * d23 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d3 * LLLL * pap1 * pap2 * pbp1 -
        2.0 * d33 * LLLL * pap1 * pap2 * pbp1 -
        4.0 * d1 * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d12 * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d13 * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d22 * LRLR * pap1 * pap2 * pbp1 -
        8.0 * d23 * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d3 * LRLR * pap1 * pap2 * pbp1 -
        4.0 * d33 * LRLR * pap1 * pap2 * pbp1 +
        2.0 * d23 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d3 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d33 * LLLL * p1p2 * papb * pbp1 +
        2.0 * d23 * LRLR * p1p2 * papb * pbp1 +
        2.0 * d3 * LRLR * p1p2 * papb * pbp1 +
        2.0 * d33 * LRLR * p1p2 * papb * pbp1 +
        2.0 * d13 * LLLL * pap2 * papb * pbp1 +
        2.0 * d23 * LLLL * pap2 * papb * pbp1 +
        2.0 * d3 * LLLL * pap2 * papb * pbp1 +
        2.0 * d33 * LLLL * pap2 * papb * pbp1 -
        2.0 * d23 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d3 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d33 * LLLL * pap2 * pow2(pbp1) -
        2.0 * d23 * LRLR * pap2 * pow2(pbp1) -
        2.0 * d3 * LRLR * pap2 * pow2(pbp1) -
        2.0 * d33 * LRLR * pap2 * pow2(pbp1) + d22 * LLLL * m1s * pap1 * pbp2 +
        2.0 * d23 * LLLL * m1s * pap1 * pbp2 + d3 * LLLL * m1s * pap1 * pbp2 +
        d33 * LLLL * m1s * pap1 * pbp2 + d22 * LRLR * m1s * pap1 * pbp2 +
        2.0 * d23 * LRLR * m1s * pap1 * pbp2 + d3 * LRLR * m1s * pap1 * pbp2 +
        d33 * LRLR * m1s * pap1 * pbp2 - d3 * LLLL * m1 * ml3 * pap1 * pbp2 -
        d3 * LRLR * m1 * ml3 * pap1 * pbp2 -
        2.0 * d12 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d13 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d22 * LLLL * pow2(pap1) * pbp2 -
        4.0 * d23 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d3 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d33 * LLLL * pow2(pap1) * pbp2 +
        2.0 * d13 * LLLL * pap1 * papb * pbp2 +
        2.0 * d23 * LLLL * pap1 * papb * pbp2 +
        2.0 * d3 * LLLL * pap1 * papb * pbp2 +
        2.0 * d33 * LLLL * pap1 * papb * pbp2 -
        2.0 * d23 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d3 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d33 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d23 * LRLR * pap1 * pbp1 * pbp2 -
        2.0 * d3 * LRLR * pap1 * pbp1 * pbp2 -
        2.0 * d33 * LRLR * pap1 * pbp1 * pbp2 - d22 * m1s * p1p2 * papb * RLRL -
        2.0 * d23 * m1s * p1p2 * papb * RLRL - d3 * m1s * p1p2 * papb * RLRL -
        d33 * m1s * p1p2 * papb * RLRL + d3 * m1 * ml3 * p1p2 * papb * RLRL +
        2.0 * d1 * m1s * pap2 * papb * RLRL +
        2.0 * d12 * m1s * pap2 * papb * RLRL +
        2.0 * d13 * m1s * pap2 * papb * RLRL +
        2.0 * d22 * m1s * pap2 * papb * RLRL +
        4.0 * d23 * m1s * pap2 * papb * RLRL +
        2.0 * d3 * m1s * pap2 * papb * RLRL +
        2.0 * d33 * m1s * pap2 * papb * RLRL -
        2.0 * d1 * m1 * ml3 * pap2 * papb * RLRL -
        2.0 * d3 * m1 * ml3 * pap2 * papb * RLRL +
        d22 * m1s * pap2 * pbp1 * RLRL + 2.0 * d23 * m1s * pap2 * pbp1 * RLRL +
        d3 * m1s * pap2 * pbp1 * RLRL + d33 * m1s * pap2 * pbp1 * RLRL +
        d3 * m1 * ml3 * pap2 * pbp1 * RLRL -
        4.0 * d1 * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d12 * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d13 * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d22 * pap1 * pap2 * pbp1 * RLRL -
        8.0 * d23 * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d3 * pap1 * pap2 * pbp1 * RLRL -
        4.0 * d33 * pap1 * pap2 * pbp1 * RLRL +
        2.0 * d23 * p1p2 * papb * pbp1 * RLRL +
        2.0 * d3 * p1p2 * papb * pbp1 * RLRL +
        2.0 * d33 * p1p2 * papb * pbp1 * RLRL -
        2.0 * d23 * pap2 * pow2(pbp1) * RLRL -
        2.0 * d3 * pap2 * pow2(pbp1) * RLRL -
        2.0 * d33 * pap2 * pow2(pbp1) * RLRL + d22 * m1s * pap1 * pbp2 * RLRL +
        2.0 * d23 * m1s * pap1 * pbp2 * RLRL + d3 * m1s * pap1 * pbp2 * RLRL +
        d33 * m1s * pap1 * pbp2 * RLRL - d3 * m1 * ml3 * pap1 * pbp2 * RLRL -
        2.0 * d23 * pap1 * pbp1 * pbp2 * RLRL -
        2.0 * d3 * pap1 * pbp1 * pbp2 * RLRL -
        2.0 * d33 * pap1 * pbp1 * pbp2 * RLRL +
        (-(d3 * m1s * p1p2 * papb) - d33 * m1s * p1p2 * papb +
         d3 * m1 * ml3 * p1p2 * papb + 2.0 * d12 * p1p2 * pap1 * papb +
         2.0 * d13 * p1p2 * pap1 * papb + 2.0 * d3 * p1p2 * pap1 * papb +
         2.0 * d33 * p1p2 * pap1 * papb + 2.0 * d1 * m1s * pap2 * papb -
         2.0 * d1 * m1 * ml3 * pap2 * papb - 2.0 * d3 * m1 * ml3 * pap2 * papb -
         2.0 * d13 * p1p2 * pow2(papb) - 2.0 * d3 * p1p2 * pow2(papb) -
         2.0 * d33 * p1p2 * pow2(papb) + d3 * m1s * pap2 * pbp1 +
         d33 * m1s * pap2 * pbp1 + d3 * m1 * ml3 * pap2 * pbp1 -
         4.0 * d1 * pap1 * pap2 * pbp1 - 2.0 * d12 * pap1 * pap2 * pbp1 -
         2.0 * d13 * pap1 * pap2 * pbp1 - 2.0 * d3 * pap1 * pap2 * pbp1 -
         2.0 * d33 * pap1 * pap2 * pbp1 + 2.0 * d3 * p1p2 * papb * pbp1 +
         2.0 * d33 * p1p2 * papb * pbp1 + 2.0 * d13 * pap2 * papb * pbp1 +
         2.0 * d3 * pap2 * papb * pbp1 + 2.0 * d33 * pap2 * papb * pbp1 -
         2.0 * d3 * pap2 * pow2(pbp1) - 2.0 * d33 * pap2 * pow2(pbp1) +
         pap1 * (-2.0 * (d12 + d13) * pap1 + 2.0 * d13 * papb +
                 d33 * (m1s - 2.0 * (pap1 - papb + pbp1)) +
                 d3 * (m1s - m1 * ml3 - 2.0 * (pap1 - papb + pbp1))) *
             pbp2 +
         d22 * (m1s - 2.0 * pap1) *
             (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
         2.0 * d23 * (m1s - 2.0 * pap1 + papb - pbp1) *
             (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) *
            RRRR +
        4.0 * d00 * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                     pap2 * pbp1 * (LRLR + RLRL) +
                     (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR) +
        d2 * (LLLL *
                  (m1 * ml3 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                               pap1 * pbp2) +
                   m1s * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) -
                   2.0 * pap1 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
              (-4.0 * pap1 * pap2 * pbp1 +
               m1 * ml3 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2) +
               m1s * (-(p1p2 * papb) + 2.0 * pap2 * papb + pap2 * pbp1 +
                      pap1 * pbp2)) *
                  (LRLR + RLRL) +
              (m1 * ml3 * (p1p2 * papb - 2.0 * pap2 * papb + pap2 * pbp1 -
                           pap1 * pbp2) +
               m1s * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) -
               2.0 * pap1 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) *
                  RRRR))) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me1 =
      (cA * cF * ivu2s2 * MPIsI * NC *
       (-((d12 * p1p2 * pap1 * papb + d2 * p1p2 * pap1 * papb +
           d22 * p1p2 * pap1 * papb + 2.0 * d23 * p1p2 * pap1 * papb +
           d3 * p1p2 * pap1 * papb + d33 * p1p2 * pap1 * papb -
           d12 * m1s * pap2 * papb - d13 * m1s * pap2 * papb -
           d2 * m1s * pap2 * papb - d22 * m1s * pap2 * papb -
           2.0 * d23 * m1s * pap2 * papb - d3 * m1s * pap2 * papb -
           d33 * m1s * pap2 * papb + d13 * p1p2 * (pap1 - papb) * papb -
           d23 * p1p2 * pow2(papb) - d3 * p1p2 * pow2(papb) -
           d33 * p1p2 * pow2(papb) + d12 * pap1 * pap2 * pbp1 +
           d13 * pap1 * pap2 * pbp1 + d2 * pap1 * pap2 * pbp1 +
           d22 * pap1 * pap2 * pbp1 + 2.0 * d23 * pap1 * pap2 * pbp1 +
           d3 * pap1 * pap2 * pbp1 + d33 * pap1 * pap2 * pbp1 +
           d13 * pap2 * papb * pbp1 + d23 * pap2 * papb * pbp1 +
           d3 * pap2 * papb * pbp1 + d33 * pap2 * papb * pbp1 -
           pap1 * ((d12 + d13 + d2 + d22 + 2.0 * d23 + d3 + d33) * pap1 -
                   (d13 + d23 + d3 + d33) * papb) *
               pbp2) *
          (LLLL - LRLR - RLRL + RRRR)) -
        d00 *
            (LLLL * (-3.0 * p1p2 * papb + pap2 * pbp1 + 3.0 * pap1 * pbp2) +
             (3.0 * p1p2 * papb + pap2 * pbp1 - 3.0 * pap1 * pbp2) *
                 (LRLR + RLRL) +
             (-3.0 * p1p2 * papb + pap2 * pbp1 + 3.0 * pap1 * pbp2) * RRRR))) /
      4.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me2 =
      -(cA * cF * d00 * ivu2s2 * MPIsI * NC * (p1p2 * papb - pap1 * pbp2) *
        (LLLL - LRLR - RLRL + RRRR)) /
      2.0;

  return real(me0 + me1 + me2);
}

double FI::Mut_B3_GLGA(double p1s, double p2s, double p3s, double p4s,
                       double pas, double pbs, double ml1s, double ml2s,
                       double ml3s, double ml4s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0, d1, d2, d3, d4, d00, d11, d22, d33, d12, d13, d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me0 =
      (cA * cF * ivt1s2 * MPIsI * NC *
       (-2.0 * d00 * LRRL * m1 * m2 * papb - d2 * LRRL * m1 * m1s * m2 * papb -
        d22 * LRRL * m1 * m1s * m2 * papb -
        2.0 * d23 * LRRL * m1 * m1s * m2 * papb -
        d3 * LRRL * m1 * m1s * m2 * papb - d33 * LRRL * m1 * m1s * m2 * papb +
        d2 * LLLL * m1 * ml3 * p1p2 * papb +
        d3 * LLLL * m1 * ml3 * p1p2 * papb +
        2.0 * d1 * LRRL * m2 * ml3 * pap1 * papb +
        2.0 * d2 * LRRL * m2 * ml3 * pap1 * papb +
        2.0 * d3 * LRRL * m2 * ml3 * pap1 * papb -
        2.0 * d1 * LLLL * p1p2 * pap1 * papb +
        2.0 * d2 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d22 * LRRL * m1 * m2 * pap1 * pbp1 +
        4.0 * d23 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d3 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d33 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d2 * LRRL * m2 * ml3 * pap1 * pbp1 -
        2.0 * d3 * LRRL * m2 * ml3 * pap1 * pbp1 -
        d2 * LLLL * m1 * ml3 * pap2 * pbp1 -
        d3 * LLLL * m1 * ml3 * pap2 * pbp1 +
        2.0 * d1 * LLLL * pap1 * pap2 * pbp1 - 8.0 * d00 * LLLL * pap1 * pbp2 -
        2.0 * d2 * LLLL * m1s * pap1 * pbp2 -
        2.0 * d22 * LLLL * m1s * pap1 * pbp2 -
        4.0 * d23 * LLLL * m1s * pap1 * pbp2 -
        2.0 * d3 * LLLL * m1s * pap1 * pbp2 -
        2.0 * d33 * LLLL * m1s * pap1 * pbp2 +
        d2 * LLLL * m1 * ml3 * pap1 * pbp2 +
        d3 * LLLL * m1 * ml3 * pap1 * pbp2 +
        2.0 * d1 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d12 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d13 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d2 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d22 * LLLL * pow2(pap1) * pbp2 +
        8.0 * d23 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d3 * LLLL * pow2(pap1) * pbp2 +
        4.0 * d33 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d3 * LLLL * m1 * ml3 * papb * pbp2 -
        4.0 * d13 * LLLL * pap1 * papb * pbp2 -
        4.0 * d23 * LLLL * pap1 * papb * pbp2 -
        4.0 * d3 * LLLL * pap1 * papb * pbp2 -
        4.0 * d33 * LLLL * pap1 * papb * pbp2 +
        4.0 * d23 * LLLL * pap1 * pbp1 * pbp2 +
        4.0 * d3 * LLLL * pap1 * pbp1 * pbp2 +
        4.0 * d33 * LLLL * pap1 * pbp1 * pbp2 -
        2.0 * d00 * m1 * m2 * papb * RLLR - d2 * m1 * m1s * m2 * papb * RLLR -
        d22 * m1 * m1s * m2 * papb * RLLR -
        2.0 * d23 * m1 * m1s * m2 * papb * RLLR -
        d3 * m1 * m1s * m2 * papb * RLLR - d33 * m1 * m1s * m2 * papb * RLLR +
        2.0 * d1 * m2 * ml3 * pap1 * papb * RLLR +
        2.0 * d2 * m2 * ml3 * pap1 * papb * RLLR +
        2.0 * d3 * m2 * ml3 * pap1 * papb * RLLR +
        2.0 * d2 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d22 * m1 * m2 * pap1 * pbp1 * RLLR +
        4.0 * d23 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d3 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d2 * m2 * ml3 * pap1 * pbp1 * RLLR -
        2.0 * d3 * m2 * ml3 * pap1 * pbp1 * RLLR +
        (d2 * m1 * ml3 * (p1p2 * papb - pap2 * pbp1) +
         d2 * pap1 * (-2.0 * m1s + m1 * ml3 + 4.0 * pap1) * pbp2 +
         2.0 * pap1 * (-((4.0 * d00 + (d22 + 2.0 * d23 + d33) * m1s) * pbp2) +
                       2.0 * (d12 * pap1 + d13 * pap1 + d22 * pap1 +
                              2.0 * d23 * pap1 + d33 * pap1 - d13 * papb -
                              d23 * papb - d33 * papb + (d23 + d33) * pbp1) *
                           pbp2 +
                       d1 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
         d3 * (-2.0 * m1s * pap1 * pbp2 +
               4.0 * pap1 * (pap1 - papb + pbp1) * pbp2 +
               m1 * ml3 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2 -
                           2.0 * papb * pbp2))) *
            RRRR)) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me1 =
      (cA * cF * ivt1s2 * MPIsI * NC *
       (d22 * LRRL * m1 * m1s * m2 * papb +
        2.0 * d23 * LRRL * m1 * m1s * m2 * papb +
        d3 * LRRL * m1 * m1s * m2 * papb + d33 * LRRL * m1 * m1s * m2 * papb -
        d3 * LRRL * m1s * m2 * ml3 * papb + d22 * LLLL * m1s * p1p2 * papb +
        2.0 * d23 * LLLL * m1s * p1p2 * papb + d3 * LLLL * m1s * p1p2 * papb +
        d33 * LLLL * m1s * p1p2 * papb - d3 * LLLL * m1 * ml3 * p1p2 * papb -
        2.0 * d12 * LLLL * p1p2 * pap1 * papb -
        2.0 * d13 * LLLL * p1p2 * pap1 * papb -
        2.0 * d22 * LLLL * p1p2 * pap1 * papb -
        4.0 * d23 * LLLL * p1p2 * pap1 * papb -
        2.0 * d3 * LLLL * p1p2 * pap1 * papb -
        2.0 * d33 * LLLL * p1p2 * pap1 * papb -
        2.0 * d22 * LRRL * m1 * m2 * pap1 * pbp1 -
        4.0 * d23 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d3 * LRRL * m1 * m2 * pap1 * pbp1 -
        2.0 * d33 * LRRL * m1 * m2 * pap1 * pbp1 +
        2.0 * d3 * LRRL * m2 * ml3 * pap1 * pbp1 -
        d22 * LLLL * m1s * pap2 * pbp1 - 2.0 * d23 * LLLL * m1s * pap2 * pbp1 -
        d3 * LLLL * m1s * pap2 * pbp1 - d33 * LLLL * m1s * pap2 * pbp1 +
        d3 * LLLL * m1 * ml3 * pap2 * pbp1 +
        2.0 * d12 * LLLL * pap1 * pap2 * pbp1 +
        2.0 * d13 * LLLL * pap1 * pap2 * pbp1 +
        2.0 * d22 * LLLL * pap1 * pap2 * pbp1 +
        4.0 * d23 * LLLL * pap1 * pap2 * pbp1 +
        2.0 * d3 * LLLL * pap1 * pap2 * pbp1 +
        2.0 * d33 * LLLL * pap1 * pap2 * pbp1 + d22 * LLLL * m1s * pap1 * pbp2 +
        2.0 * d23 * LLLL * m1s * pap1 * pbp2 + d3 * LLLL * m1s * pap1 * pbp2 +
        d33 * LLLL * m1s * pap1 * pbp2 - d3 * LLLL * m1 * ml3 * pap1 * pbp2 -
        2.0 * d12 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d13 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d22 * LLLL * pow2(pap1) * pbp2 -
        4.0 * d23 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d3 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d33 * LLLL * pow2(pap1) * pbp2 -
        2.0 * d23 * LLLL * m1s * papb * pbp2 -
        2.0 * d3 * LLLL * m1s * papb * pbp2 -
        2.0 * d33 * LLLL * m1s * papb * pbp2 +
        2.0 * d3 * LLLL * m1 * ml3 * papb * pbp2 +
        4.0 * d13 * LLLL * pap1 * papb * pbp2 +
        4.0 * d23 * LLLL * pap1 * papb * pbp2 +
        4.0 * d3 * LLLL * pap1 * papb * pbp2 +
        4.0 * d33 * LLLL * pap1 * papb * pbp2 +
        d22 * m1 * m1s * m2 * papb * RLLR +
        2.0 * d23 * m1 * m1s * m2 * papb * RLLR +
        d3 * m1 * m1s * m2 * papb * RLLR + d33 * m1 * m1s * m2 * papb * RLLR -
        d3 * m1s * m2 * ml3 * papb * RLLR -
        2.0 * d22 * m1 * m2 * pap1 * pbp1 * RLLR -
        4.0 * d23 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d3 * m1 * m2 * pap1 * pbp1 * RLLR -
        2.0 * d33 * m1 * m2 * pap1 * pbp1 * RLLR +
        2.0 * d3 * m2 * ml3 * pap1 * pbp1 * RLLR +
        ((2.0 * d23 * m1s + d3 * m1s + d33 * m1s - d3 * m1 * ml3 -
          2.0 * (d12 + d13 + 2.0 * d23 + d3 + d33) * pap1) *
             (p1p2 * papb - pap2 * pbp1) +
         (pap1 * (2.0 * d23 * m1s + d3 * m1s + d33 * m1s - d3 * m1 * ml3 -
                  2.0 * (d12 + d13 + 2.0 * d23 + d3 + d33) * pap1) +
          2.0 * (-(d23 * m1s) - d3 * m1s - d33 * m1s + d3 * m1 * ml3 +
                 2.0 * (d13 + d23 + d3 + d33) * pap1) *
              papb) *
             pbp2 +
         d22 * (m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) *
            RRRR +
        d2 * (LRRL * m2 * (m1 - ml3) * (m1s * papb - 2.0 * pap1 * pbp1) +
              LLLL * (m1s - m1 * ml3 - 2.0 * pap1) *
                  (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
              m2 * (m1 - ml3) * (m1s * papb - 2.0 * pap1 * pbp1) * RLLR +
              (m1s - m1 * ml3 - 2.0 * pap1) *
                  (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR) +
        4.0 * d00 * (m1 * m2 * papb * (LRRL + RLLR) +
                     2.0 * pap1 * pbp2 * (LLLL + RRRR)))) /
      8.0;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  complex<double> me2 = -(cA * cF * d00 * ivt1s2 * MPIsI * NC *
                          (LRRL * m1 * m2 * papb +
                           LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                           m1 * m2 * papb * RLLR - p1p2 * papb * RRRR +
                           pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR)) /
                        4.0;

  return real(me0 + me1 + me2);
}
