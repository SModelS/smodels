#include <complex>
#include <iostream>
//#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "utils.h"

using namespace std;

double FI::MVtbo1s(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {

  // Constructs general four point function.
  // Npf *dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  std::complex<double> me0 =
      8.0 * ivs3v2 *
      (-2.0 * LL * m1 * m2 *
           ((d23 + d33) * pap1 * papb -
            ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                std::pow(papb, 2) -
            d33 * pap1 * pbp1 + (d13 + d23 + d33) * papb * pbp1) *
           RLLL +
       LL * ((d33 * m1s - 2.0 * (d23 + d33) * pap1) *
                 (p1p2 * papb - pap2 * pbp1) +
             (d33 * m1s * pap1 - 2.0 * (d23 + d33) * std::pow(pap1, 2) -
              2.0 * (d13 + d23 + d33) * m1s * papb +
              4.0 * ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                  pap1 * papb) *
                 pbp2) *
           RLRL +
       (d33 * (2.0 * LRRR * m1 * m2 * pap1 * pbp1 +
               LRLR * m1s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) +
        2.0 *
            (LRRR * m1 * m2 * papb *
                 (((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                      papb -
                  (d13 + d23 + d33) * pbp1) -
             LRLR *
                 ((d13 + d23 + d33) * m1s -
                  2.0 * ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                      pap1) *
                 papb * pbp2 -
             (d23 + d33) * pap1 * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb -
                                   LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2))) *
           RR +
       2.0 * d00 *
           (LL * m1 * m2 * papb * RLLL + 2.0 * LL * pap1 * pbp2 * RLRL +
            LRRR * m1 * m2 * papb * RR + 2.0 * LRLR * pap1 * pbp2 * RR) -
       2.0 * d3 * pap1 * (LL * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                                pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) +
                          (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb -
                           LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2) *
                              RR));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);

  // Part of scalar integrals of O(1/eps^1) times squared matrix element of
  // O(eps^1).
  std::complex<double> me1 =
      -8.0 * ivs3v2 *
      (d33 * (LL * m1 * m2 * (m1s * papb + 2.0 * pap1 * pbp1) * RLLL +
              2.0 * LL * m1s * p1p2 * papb * RLRL +
              (2.0 * LRLR * m1s * p1p2 * papb +
               LRRR * m1 * m2 * (m1s * papb + 2.0 * pap1 * pbp1)) *
                  RR) +
       2.0 * d00 *
           (LL * (3.0 * m1 * m2 * papb * RLLL +
                  (p1p2 * papb + 3.0 * pap2 * pbp1 + pap1 * pbp2) * RLRL) +
            (3.0 * LRRR * m1 * m2 * papb +
             LRLR * (p1p2 * papb + 3.0 * pap2 * pbp1 + pap1 * pbp2)) *
                RR) +
       2.0 *
           (LL * m1 * m2 * papb * (((d1 + d2 + d3) + (d2 + d3) +
                                    2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                                       papb -
                                   2.0 * (d13 + d23 + d33) * pbp1) *
                RLLL +
            LL * (-((d13 + d23 + d33) * p1p2 * papb * pbp1) +
                  2.0 * (d2 + d3) * pap2 * papb * pbp1 -
                  (d13 + d23 + d33) * pap2 * std::pow(pbp1, 2) -
                  (d13 + d23 + d33) * m1s * papb * pbp2 +
                  (d13 + d23 + d33) * pap1 * pbp1 * pbp2 +
                  (d1 + d2 + d3) * papb *
                      (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                  (d12 + d13 + d22 + 2.0 * d23 + d33) * papb *
                      (p1p2 * papb + pap2 * pbp1 + pap1 * pbp2)) *
                RLRL +
            d3 * LL *
                (-(m1 * m2 * papb * (pap1 + pbp1) * RLLL) +
                 (m1s * pap2 * papb - pap2 * std::pow(pbp1, 2) -
                  p1p2 * papb * (2.0 * pap1 + pbp1) + pap1 * pbp1 * pbp2) *
                     RLRL) -
            d3 * (LRRR * m1 * m2 * papb * (pap1 + pbp1) +
                  LRLR * (-(m1s * pap2 * papb) + pap2 * std::pow(pbp1, 2) +
                          p1p2 * papb * (2.0 * pap1 + pbp1) -
                          pap1 * pbp1 * pbp2)) *
                RR +
            (((d12 + d13 + d22 + 2.0 * d23 + d33) * papb -
              (d13 + d23 + d33) * pbp1) *
                 (2.0 * LRRR * m1 * m2 * papb + LRLR * p1p2 * papb +
                  LRLR * pap2 * pbp1) +
             (d2 + d3) * papb *
                 (LRRR * m1 * m2 * papb + 2.0 * LRLR * pap2 * pbp1) +
             LRLR * (-((d13 + d23 + d33) * m1s * papb) +
                     (d12 + d13 + d22 + 2.0 * d23 + d33) * pap1 * papb +
                     (d13 + d23 + d33) * pap1 * pbp1) *
                 pbp2 +
             (d1 + d2 + d3) * papb *
                 (LRRR * m1 * m2 * papb +
                  LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2))) *
                RR -
            2.0 * (d23 + d33) * pap1 * papb *
                (LL * m1 * m2 * RLLL + LL * p1p2 * RLRL + LRRR * m1 * m2 * RR +
                 LRLR * p1p2 * RR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  // Part of scalar integrals of O(1/eps^2) times squared matrix element of
  // O(eps^2).
  std::complex<double> me2 =
      8.0 * ivs3v2 *
      (d00 * (6.0 * LL * m1 * m2 * papb * RLLL + 4.0 * LL * p1p2 * papb * RLRL +
              4.0 * LL * pap2 * pbp1 * RLRL + 6.0 * LRRR * m1 * m2 * papb * RR +
              4.0 * LRLR * p1p2 * papb * RR + 4.0 * LRLR * pap2 * pbp1 * RR) +
       2.0 * papb *
           (-(LL * m1 * m2 *
              ((d23 + d33) * pap1 -
               ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb +
               (d3 + (d13 + d23 + d33)) * pbp1) *
              RLLL) +
            LL * (-((d23 + d33) * m1s * pap2) +
                  ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * p1p2 *
                      papb -
                  2.0 * (d3 + (d13 + d23 + d33)) * p1p2 * pbp1 +
                  (d3 + (d13 + d23 + d33)) * m1s * pbp2 +
                  ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                      (pap2 * pbp1 - pap1 * pbp2)) *
                RLRL +
            (-((d23 + d33) * (LRRR * m1 * m2 * pap1 + LRLR * m1s * pap2)) +
             (d2 + d3) * LRRR * m1 * m2 * papb +
             (d12 + d13 + d22 + 2.0 * d23 + d33) * LRRR * m1 * m2 * papb +
             (d2 + d3) * LRLR * p1p2 * papb +
             (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * p1p2 * papb -
             d3 * LRRR * m1 * m2 * pbp1 -
             (d13 + d23 + d33) * LRRR * m1 * m2 * pbp1 -
             2.0 * d3 * LRLR * p1p2 * pbp1 -
             2.0 * (d13 + d23 + d33) * LRLR * p1p2 * pbp1 +
             (d2 + d3) * LRLR * pap2 * pbp1 +
             (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * pap2 * pbp1 +
             LRLR * ((d3 + (d13 + d23 + d33)) * m1s -
                     ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1) *
                 pbp2) *
                RR) +
       d33 * m1s * (LL * (m1 * m2 * papb * RLLL +
                          (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLRL) +
                    (LRRR * m1 * m2 * papb +
                     LRLR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
                        RR));

  return std::real(me0 + me1 + me2);
}

double FI::MVtbo1t(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      4.0 * ivt1s2 *
      (d33 * m1s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
           (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
       4.0 * d00 * (LLLL * LR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                    LL * pap1 * pbp2 * RLRL + LRLR * pap1 * pbp2 * RR +
                    (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) -
       2.0 *
           (-((d2 + d3) * LLLL * LR * p1p2 * std::pow(papb, 2)) -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * p1p2 *
                std::pow(papb, 2) +
            (d13 + d23 + d33) * LLLL * LR * p1p2 * papb * pbp1 +
            (d2 + d3) * LLLL * LR * pap2 * papb * pbp1 +
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap2 * papb *
                pbp1 -
            (d13 + d23 + d33) * LLLL * LR * pap2 * std::pow(pbp1, 2) -
            2.0 * (d1 + d2 + d3) * LLLL * LR * pap1 * papb * pbp2 +
            (d2 + d3) * LLLL * LR * pap1 * papb * pbp2 -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap1 * papb *
                pbp2 +
            (d13 + d23 + d33) * LLLL * LR * pap1 * pbp1 * pbp2 +
            (d13 + d23 + d33) * LL * m1s * papb * pbp2 * RLRL -
            2.0 * (d1 + d2 + d3) * LL * pap1 * papb * pbp2 * RLRL -
            2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * pap1 * papb *
                pbp2 * RLRL +
            (d13 + d23 + d33) * LRLR * m1s * papb * pbp2 * RR -
            2.0 * (d1 + d2 + d3) * LRLR * pap1 * papb * pbp2 * RR -
            2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * pap1 * papb *
                pbp2 * RR -
            ((((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb -
              (d13 + d23 + d33) * pbp1) *
                 (p1p2 * papb - pap2 * pbp1) +
             pap1 * ((2.0 * (d1 + d2 + d3) - (d2 + d3) +
                      (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                         papb -
                     (d13 + d23 + d33) * pbp1) *
                 pbp2) *
                RL * RRRR +
            (d23 + d33) * pap1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
            d3 * (LLLL * LR * ((pap1 + pbp1) * (p1p2 * papb - pap2 * pbp1) +
                               (-(m1s * papb) + pap1 * (pap1 + pbp1)) * pbp2) +
                  pap1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                      (LL * RLRL + LRLR * RR) +
                  ((pap1 + pbp1) * (p1p2 * papb - pap2 * pbp1) +
                   (-(m1s * papb) + pap1 * (pap1 + pbp1)) * pbp2) *
                      RL * RRRR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me1 =
      -8.0 * ivt1s2 *
      (((((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb -
         (d3 + (d13 + d23 + d33)) * pbp1) *
            (p1p2 * papb - pap2 * pbp1) +
        (((d3 + (d13 + d23 + d33)) * m1s -
          ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1) *
             papb -
         (d3 + (d13 + d23 + d33)) * pap1 * pbp1) *
            pbp2) *
           (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR) +
       d00 *
           (LLLL * LR * (3.0 * p1p2 * papb - 3.0 * pap2 * pbp1 + pap1 * pbp2) -
            (3.0 * p1p2 * papb - 3.0 * pap2 * pbp1 - pap1 * pbp2) *
                (LL * RLRL + LRLR * RR) +
            (3.0 * p1p2 * papb - 3.0 * pap2 * pbp1 + pap1 * pbp2) * RL * RRRR));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me2 = 16.0 * d00 * ivt1s2 * (p1p2 * papb - pap2 * pbp1) *
                             (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR);

  return std::real(me0 + me1 + me2);
}

double FI::MVtbo1u(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // LL<->RL; RR<->LR;
  // Npf * dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      4.0 * ivu2s2 *
      (2.0 * (d23 + d33) * LR * LRLR * m1 * m2 * pap1 * papb -
       2.0 * (d1 + d2 + d3) * LR * LRLR * m1 * m2 * std::pow(papb, 2) -
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LR * LRLR * m1 * m2 *
           std::pow(papb, 2) -
       2.0 * d33 * LR * LRLR * m1 * m2 * pap1 * pbp1 +
       2.0 * (d13 + d23 + d33) * LR * LRLR * m1 * m2 * papb * pbp1 +
       2.0 * (d23 + d33) * m1 * m2 * pap1 * papb * RL * RLRL -
       2.0 * (d1 + d2 + d3) * m1 * m2 * std::pow(papb, 2) * RL * RLRL -
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * m1 * m2 * std::pow(papb, 2) *
           RL * RLRL -
       2.0 * d33 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
       2.0 * (d13 + d23 + d33) * m1 * m2 * papb * pbp1 * RL * RLRL +
       d33 * LLLL * m1s * p1p2 * papb * RR -
       2.0 * (d23 + d33) * LLLL * p1p2 * pap1 * papb * RR +
       2.0 * (d1 + d2 + d3) * LLLL * p1p2 * std::pow(papb, 2) * RR +
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * p1p2 *
           std::pow(papb, 2) * RR +
       d33 * LLLL * m1s * pap2 * pbp1 * RR -
       2.0 * (d23 + d33) * LLLL * pap1 * pap2 * pbp1 * RR -
       2.0 * (d13 + d23 + d33) * LLLL * p1p2 * papb * pbp1 * RR -
       2.0 * (d1 + d2 + d3) * LLLL * pap2 * papb * pbp1 * RR +
       4.0 * (d2 + d3) * LLLL * pap2 * papb * pbp1 * RR +
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * pap2 * papb * pbp1 *
           RR -
       2.0 * (d13 + d23 + d33) * LLLL * pap2 * std::pow(pbp1, 2) * RR -
       d33 * LLLL * m1s * pap1 * pbp2 * RR +
       2.0 * (d23 + d33) * LLLL * std::pow(pap1, 2) * pbp2 * RR -
       2.0 * (d1 + d2 + d3) * LLLL * pap1 * papb * pbp2 * RR -
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * pap1 * papb * pbp2 *
           RR +
       2.0 * (d13 + d23 + d33) * LLLL * pap1 * pbp1 * pbp2 * RR +
       d33 * LL * m1s * p1p2 * papb * RRRR -
       2.0 * (d23 + d33) * LL * p1p2 * pap1 * papb * RRRR +
       2.0 * (d1 + d2 + d3) * LL * p1p2 * std::pow(papb, 2) * RRRR +
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * p1p2 *
           std::pow(papb, 2) * RRRR +
       d33 * LL * m1s * pap2 * pbp1 * RRRR -
       2.0 * (d23 + d33) * LL * pap1 * pap2 * pbp1 * RRRR -
       2.0 * (d13 + d23 + d33) * LL * p1p2 * papb * pbp1 * RRRR -
       2.0 * (d1 + d2 + d3) * LL * pap2 * papb * pbp1 * RRRR +
       4.0 * (d2 + d3) * LL * pap2 * papb * pbp1 * RRRR +
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * pap2 * papb * pbp1 *
           RRRR -
       2.0 * (d13 + d23 + d33) * LL * pap2 * std::pow(pbp1, 2) * RRRR -
       d33 * LL * m1s * pap1 * pbp2 * RRRR +
       2.0 * (d23 + d33) * LL * std::pow(pap1, 2) * pbp2 * RRRR -
       2.0 * (d1 + d2 + d3) * LL * pap1 * papb * pbp2 * RRRR -
       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * pap1 * papb * pbp2 *
           RRRR +
       2.0 * (d13 + d23 + d33) * LL * pap1 * pbp1 * pbp2 * RRRR -
       2.0 * d00 * (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL -
                    2.0 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                        (LLLL * RR + LL * RRRR)) +
       2.0 * d3 * (LR * LRLR * m1 * m2 * pap1 * papb +
                   m1 * m2 * pap1 * papb * RL * RLRL -
                   (p1p2 * papb * (pap1 + pbp1) +
                    pap2 * (-(m1s * papb) + pbp1 * (pap1 + pbp1)) -
                    pap1 * (pap1 + pbp1) * pbp2) *
                       (LLLL * RR + LL * RRRR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me1 =
      -4.0 * ivu2s2 *
      (-2.0 * ((d23 + d33) * m1s * pap2 * papb +
               pbp1 * (d3 * p1p2 * papb + (d13 + d23 + d33) * p1p2 * papb -
                       2.0 * (d2 + d3) * pap2 * papb -
                       2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * pap2 * papb +
                       d3 * pap2 * pbp1 + (d13 + d23 + d33) * pap2 * pbp1 -
                       (d3 + (d13 + d23 + d33)) * pap1 * pbp2)) *
           (LLLL * RR + LL * RRRR) +
       d33 * (LR * LRLR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) +
              m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RL * RLRL +
              m1s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                  (LLLL * RR + LL * RRRR)) +
       2.0 * d00 * (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL +
                    (p1p2 * papb + 3.0 * pap2 * pbp1 - pap1 * pbp2) *
                        (LLLL * RR + LL * RRRR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me2 =
      8.0 * d00 * ivu2s2 *
      (LR * LRLR * m1 * m2 * papb + m1 * m2 * papb * RL * RLRL -
       (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * (LLLL * RR + LL * RRRR));

  return std::real(me0 + me1 + me2);
}

double FI::MVubo1s(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      -8.0 * ivs3v2 *
      (d33 * LL * (m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL +
                   2.0 * m1 * m2 * pap2 * pbp2 * RLRL) +
       d33 * (2.0 * LRLR * m1 * m2 * pap2 * pbp2 +
              LRRR * m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
           RR +
       2.0 *
           (-(LL *
              (((d13 + d23 + d33) * m2s * papb -
                2.0 * pap2 *
                    (d00 +
                     ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                         papb)) *
                   pbp1 +
               (d23 + d33) * pap2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
              RLLL) +
            LL * m1 * m2 * papb *
                (d00 - (d23 + d33) * pap2 +
                 ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb -
                 (d13 + d23 + d33) * pbp2) *
                RLRL +
            (d00 * (LRLR * m1 * m2 * papb + 2.0 * LRRR * pap2 * pbp1) +
             papb * (-((d13 + d23 + d33) * LRRR * m2s * pbp1) +
                     ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                         (LRLR * m1 * m2 * papb + 2.0 * LRRR * pap2 * pbp1) -
                     (d13 + d23 + d33) * LRLR * m1 * m2 * pbp2) -
             (d23 + d33) * pap2 * (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb +
                                   LRRR * pap2 * pbp1 - LRRR * pap1 * pbp2)) *
                RR -
            d3 * pap2 * (LL * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL -
                               pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                         (LRLR * m1 * m2 * papb +
                          LRRR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) *
                             RR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me1 =
      8.0 * ivs3v2 *
      (2.0 * d00 *
           (LL * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL +
                  3.0 * pap1 * pbp2 * RLLL + 3.0 * m1 * m2 * papb * RLRL) +
            (3.0 * LRLR * m1 * m2 * papb +
             LRRR * (p1p2 * papb + pap2 * pbp1 + 3.0 * pap1 * pbp2)) *
                RR) +
       d33 * (2.0 * LL * m2s * p1p2 * papb * RLLL +
              LL * m1 * m2 * (m2s * papb + 2.0 * pap2 * pbp2) * RLRL +
              (2.0 * LRRR * m2s * p1p2 * papb +
               LRLR * m1 * m2 * (m2s * papb + 2.0 * pap2 * pbp2)) *
                  RR) +
       2.0 *
           (LL * (papb *
                      (-((d13 + d23 + d33) * m2s * pbp1) +
                       ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                           (p1p2 * papb + pap2 * pbp1)) -
                  (((d1 + d2 + d3) - 2.0 * (d2 + d3) -
                    (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                       pap1 * papb +
                   (d13 + d23 + d33) * (p1p2 * papb - pap2 * pbp1)) *
                      pbp2 -
                  (d13 + d23 + d33) * pap1 * std::pow(pbp2, 2)) *
                RLLL +
            LL * m1 * m2 * papb * (((d1 + d2 + d3) + (d2 + d3) +
                                    2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                                       papb -
                                   2.0 * (d13 + d23 + d33) * pbp2) *
                RLRL +
            d3 * LL * (m2s * pap1 * papb * RLLL -
                       p1p2 * papb * (2.0 * pap2 + pbp2) * RLLL +
                       pbp2 * (pap2 * pbp1 - pap1 * pbp2) * RLLL -
                       m1 * m2 * papb * (pap2 + pbp2) * RLRL) +
            (papb *
                 ((d2 + d3) * LRLR * m1 * m2 * papb +
                  2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * m1 * m2 *
                      papb +
                  (d12 + d13 + d22 + 2.0 * d23 + d33) * LRRR * p1p2 * papb -
                  (d13 + d23 + d33) * LRRR * m2s * pbp1 +
                  (d12 + d13 + d22 + 2.0 * d23 + d33) * LRRR * pap2 * pbp1 +
                  (d1 + d2 + d3) * (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb +
                                    LRRR * pap2 * pbp1)) -
             (((d1 + d2 + d3) - 2.0 * (d2 + d3) -
               (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                  LRRR * pap1 * papb +
              (d13 + d23 + d33) * (2.0 * LRLR * m1 * m2 * papb +
                                   LRRR * p1p2 * papb - LRRR * pap2 * pbp1)) *
                 pbp2 -
             (d13 + d23 + d33) * LRRR * pap1 * std::pow(pbp2, 2)) *
                RR +
            d3 *
                (-(LRLR * m1 * m2 * papb * (pap2 + pbp2)) +
                 LRRR * (m2s * pap1 * papb - p1p2 * papb * (2.0 * pap2 + pbp2) +
                         pbp2 * (pap2 * pbp1 - pap1 * pbp2))) *
                RR -
            2.0 * (d23 + d33) * pap2 * papb *
                (LL * p1p2 * RLLL + LL * m1 * m2 * RLRL + LRLR * m1 * m2 * RR +
                 LRRR * p1p2 * RR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me2 =
      -8.0 * ivs3v2 *
      (d00 * (4.0 * LL * p1p2 * papb * RLLL + 4.0 * LL * pap1 * pbp2 * RLLL +
              6.0 * LL * m1 * m2 * papb * RLRL +
              6.0 * LRLR * m1 * m2 * papb * RR + 4.0 * LRRR * p1p2 * papb * RR +
              4.0 * LRRR * pap1 * pbp2 * RR) +
       d33 * m2s * (LL * (p1p2 * papb * RLLL - pap2 * pbp1 * RLLL +
                          pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
                    (LRLR * m1 * m2 * papb +
                     LRRR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) *
                        RR) +
       2.0 * papb *
           (LL * ((d3 + (d13 + d23 + d33)) * m2s * pbp1 +
                  ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                      (p1p2 * papb - pap2 * pbp1) +
                  (-2.0 * (d3 + (d13 + d23 + d33)) * p1p2 +
                   ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1) *
                      pbp2) *
                RLLL +
            LL * m1 * m2 *
                (((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb -
                 (d3 + (d13 + d23 + d33)) * pbp2) *
                RLRL +
            ((d3 + (d13 + d23 + d33)) *
                 (LRRR * m2s * pbp1 - LRLR * m1 * m2 * pbp2 -
                  2.0 * LRRR * p1p2 * pbp2) +
             (d12 + d13 + d22 + 2.0 * d23 + d33) *
                 (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb -
                  LRRR * pap2 * pbp1 + LRRR * pap1 * pbp2) +
             (d2 + d3) * (LRLR * m1 * m2 * papb +
                          LRRR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2))) *
                RR -
            (d23 + d33) *
                (LL * m2s * pap1 * RLLL + LL * m1 * m2 * pap2 * RLRL +
                 LRRR * m2s * pap1 * RR + LRLR * m1 * m2 * pap2 * RR)));

  return std::real(me0 + me1 + me2);
}

double FI::MVubo1t(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      -4.0 * ivt1s2 *
      (2.0 * d00 *
           (-2.0 * LLLL * LR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
            m1 * m2 * papb * (LL * RLRL + LRLR * RR) -
            2.0 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) -
       d33 * (LLLL * LR * m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) -
              2.0 * m1 * m2 * pap2 * pbp2 * (LL * RLRL + LRLR * RR) +
              m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
       2.0 *
           (-((d1 + d2 + d3) * LLLL * LR * p1p2 * std::pow(papb, 2)) -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * p1p2 *
                std::pow(papb, 2) +
            (d1 + d2 + d3) * LLLL * LR * pap2 * papb * pbp1 +
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap2 * papb *
                pbp1 +
            (d13 + d23 + d33) * LLLL * LR * p1p2 * papb * pbp2 +
            (d1 + d2 + d3) * LLLL * LR * pap1 * papb * pbp2 -
            2.0 * (d2 + d3) * LLLL * LR * pap1 * papb * pbp2 -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap1 * papb *
                pbp2 -
            (d13 + d23 + d33) * LLLL * LR * pap2 * pbp1 * pbp2 +
            (d13 + d23 + d33) * LLLL * LR * pap1 * std::pow(pbp2, 2) +
            (d23 + d33) * LLLL * LR * pap2 *
                (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
            d3 * LLLL * LR *
                (-(m2s * pap1 * papb) +
                 (pap2 + pbp2) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) -
            d3 * LL * m1 * m2 * pap2 * papb * RLRL -
            (d23 + d33) * LL * m1 * m2 * pap2 * papb * RLRL +
            (d1 + d2 + d3) * LL * m1 * m2 * std::pow(papb, 2) * RLRL +
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * m1 * m2 *
                std::pow(papb, 2) * RLRL -
            (d13 + d23 + d33) * LL * m1 * m2 * papb * pbp2 * RLRL -
            d3 * LRLR * m1 * m2 * pap2 * papb * RR -
            (d23 + d33) * LRLR * m1 * m2 * pap2 * papb * RR +
            (d1 + d2 + d3) * LRLR * m1 * m2 * std::pow(papb, 2) * RR +
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * m1 * m2 *
                std::pow(papb, 2) * RR -
            (d13 + d23 + d33) * LRLR * m1 * m2 * papb * pbp2 * RR +
            (((d23 + d33) * pap2 -
              ((d1 + d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * papb) *
                 (p1p2 * papb - pap2 * pbp1) +
             ((d23 + d33) * pap1 * pap2 + (d13 + d23 + d33) * p1p2 * papb +
              ((d1 + d2 + d3) - 2.0 * (d2 + d3) -
               (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                  pap1 * papb -
              (d13 + d23 + d33) * pap2 * pbp1) *
                 pbp2 +
             (d13 + d23 + d33) * pap1 * std::pow(pbp2, 2) +
             d3 * (-(m2s * pap1 * papb) +
                   (pap2 + pbp2) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2))) *
                RL * RRRR));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me1 =
      -4.0 * ivt1s2 *
      (-2.0 *
           ((d23 + d33) * m2s * pap1 * papb +
            (((d3 + (d13 + d23 + d33)) * p1p2 -
              2.0 * ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1) *
                 papb -
             (d3 + (d13 + d23 + d33)) * pap2 * pbp1) *
                pbp2 +
            (d3 + (d13 + d23 + d33)) * pap1 * std::pow(pbp2, 2)) *
           (LLLL * LR + RL * RRRR) +
       d33 * (LLLL * LR * m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
              m1 * m2 * (m2s * papb - 2.0 * pap2 * pbp2) *
                  (LL * RLRL + LRLR * RR) +
              m2s * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL * RRRR) +
       2.0 * d00 *
           (LLLL * LR * (p1p2 * papb - pap2 * pbp1 + 3.0 * pap1 * pbp2) +
            m1 * m2 * papb * (LL * RLRL + LRLR * RR) +
            (p1p2 * papb - pap2 * pbp1 + 3.0 * pap1 * pbp2) * RL * RRRR));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me2 =
      8.0 * d00 * ivt1s2 *
      (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
       LL * m1 * m2 * papb * RLRL + LRLR * m1 * m2 * papb * RR -
       p1p2 * papb * RL * RRRR + pap2 * pbp1 * RL * RRRR +
       pap1 * pbp2 * RL * RRRR);

  return std::real(me0 + me1 + me2);
}

double FI::MVubo1u(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      4.0 * ivu2s2 *
      (d33 * m2s * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
           (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
       4.0 * d00 * (LLLL * LR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                    LL * pap2 * pbp1 * RLRL + LRLR * pap2 * pbp1 * RR +
                    (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RL * RRRR) -
       2.0 *
           (-((d2 + d3) * LLLL * LR * p1p2 * std::pow(papb, 2)) -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * p1p2 *
                std::pow(papb, 2) -
            2.0 * (d1 + d2 + d3) * LLLL * LR * pap2 * papb * pbp1 +
            (d2 + d3) * LLLL * LR * pap2 * papb * pbp1 -
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap2 * papb *
                pbp1 +
            (d13 + d23 + d33) * LLLL * LR * p1p2 * papb * pbp2 +
            (d2 + d3) * LLLL * LR * pap1 * papb * pbp2 +
            (d12 + d13 + d22 + 2.0 * d23 + d33) * LLLL * LR * pap1 * papb *
                pbp2 +
            (d13 + d23 + d33) * LLLL * LR * pap2 * pbp1 * pbp2 -
            (d13 + d23 + d33) * LLLL * LR * pap1 * std::pow(pbp2, 2) +
            (d13 + d23 + d33) * LL * m2s * papb * pbp1 * RLRL -
            2.0 * (d1 + d2 + d3) * LL * pap2 * papb * pbp1 * RLRL -
            2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LL * pap2 * papb *
                pbp1 * RLRL +
            (d13 + d23 + d33) * LRLR * m2s * papb * pbp1 * RR -
            2.0 * (d1 + d2 + d3) * LRLR * pap2 * papb * pbp1 * RR -
            2.0 * (d12 + d13 + d22 + 2.0 * d23 + d33) * LRLR * pap2 * papb *
                pbp1 * RR -
            (papb * (((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * p1p2 *
                         papb +
                     (2.0 * (d1 + d2 + d3) - (d2 + d3) +
                      (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                         pap2 * pbp1) -
             (((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1 * papb +
              (d13 + d23 + d33) * (p1p2 * papb + pap2 * pbp1)) *
                 pbp2 +
             (d13 + d23 + d33) * pap1 * std::pow(pbp2, 2)) *
                RL * RRRR +
            (d23 + d33) * pap2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                (LLLL * LR + LL * RLRL + LRLR * RR + RL * RRRR) +
            d3 * (LLLL * LR *
                      ((std::pow(pap2, 2) - m2s * papb) * pbp1 +
                       pap2 * (-pap1 + pbp1) * pbp2 - pap1 * std::pow(pbp2, 2) +
                       p1p2 * papb * (pap2 + pbp2)) +
                  pap2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) *
                      (LL * RLRL + LRLR * RR) +
                  ((std::pow(pap2, 2) - m2s * papb) * pbp1 +
                   pap2 * (-pap1 + pbp1) * pbp2 - pap1 * std::pow(pbp2, 2) +
                   p1p2 * papb * (pap2 + pbp2)) *
                      RL * RRRR)));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 1, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me1 =
      -8.0 * ivu2s2 *
      ((papb * ((d3 + (d13 + d23 + d33)) * m2s * pbp1 +
                ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) *
                    (p1p2 * papb - pap2 * pbp1)) -
        (((d3 + (d13 + d23 + d33)) * p1p2 +
          ((d2 + d3) + (d12 + d13 + d22 + 2.0 * d23 + d33)) * pap1) *
             papb +
         (d3 + (d13 + d23 + d33)) * pap2 * pbp1) *
            pbp2 +
        (d3 + (d13 + d23 + d33)) * pap1 * std::pow(pbp2, 2)) *
           (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR) +
       d00 *
           (LLLL * LR * (3.0 * p1p2 * papb + pap2 * pbp1 - 3.0 * pap1 * pbp2) -
            (3.0 * p1p2 * papb - pap2 * pbp1 - 3.0 * pap1 * pbp2) *
                (LL * RLRL + LRLR * RR) +
            (3.0 * p1p2 * papb + pap2 * pbp1 - 3.0 * pap1 * pbp2) * RL * RRRR));

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps + 2, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me2 = 16.0 * d00 * ivu2s2 * (p1p2 * papb - pap1 * pbp2) *
                             (LLLL * LR - LL * RLRL - LRLR * RR + RL * RRRR);

  return std::real(me0 + me1 + me2);
}

double FI::MVtbo2s(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      4.0 * ivs3v2 *
      (d33 * LR * LRLR * m1 * m1s * m2 * papb -
       d33 * LR * LRRR * m1s * p1p2 * papb +
       2.0 * (d23 + d33) * LR * LRRR * m1s * pap2 * papb -
       2.0 * d33 * LR * LRLR * m1 * m2 * pap1 * pbp1 +
       d33 * LR * LRRR * m1s * pap2 * pbp1 -
       4.0 * (d23 + d33) * LR * LRRR * pap1 * pap2 * pbp1 +
       2.0 * (d13 + d23 + d33) * LR * LRRR * p1p2 * papb * pbp1 -
       2.0 * (d13 + d23 + d33) * LR * LRRR * pap2 * std::pow(pbp1, 2) +
       d33 * LR * LRRR * m1s * pap1 * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRRR * pap1 * pbp1 * pbp2 +
       d33 * LRLL * m1 * m1s * m2 * papb * RL -
       d33 * LRRL * m1s * p1p2 * papb * RL +
       2.0 * (d23 + d33) * LRRL * p1p2 * pap1 * papb * RL -
       2.0 * d33 * LRLL * m1 * m2 * pap1 * pbp1 * RL +
       d33 * LRRL * m1s * pap2 * pbp1 * RL -
       2.0 * (d23 + d33) * LRRL * pap1 * pap2 * pbp1 * RL +
       d33 * LRRL * m1s * pap1 * pbp2 * RL -
       2.0 * (d23 + d33) * LRRL * std::pow(pap1, 2) * pbp2 * RL +
       2.0 * (d13 + d23 + d33) * LRRL * m1s * papb * pbp2 * RL -
       4.0 * (d13 + d23 + d33) * LRRL * pap1 * pbp1 * pbp2 * RL -
       d33 * m1s * p1p2 * papb * RL * RLLL +
       2.0 * (d23 + d33) * m1s * pap2 * papb * RL * RLLL +
       d33 * m1s * pap2 * pbp1 * RL * RLLL -
       4.0 * (d23 + d33) * pap1 * pap2 * pbp1 * RL * RLLL +
       2.0 * (d13 + d23 + d33) * p1p2 * papb * pbp1 * RL * RLLL -
       2.0 * (d13 + d23 + d33) * pap2 * std::pow(pbp1, 2) * RL * RLLL +
       d33 * m1s * pap1 * pbp2 * RL * RLLL -
       2.0 * (d13 + d23 + d33) * pap1 * pbp1 * pbp2 * RL * RLLL -
       d33 * LR * m1s * p1p2 * papb * RLLR +
       2.0 * (d23 + d33) * LR * p1p2 * pap1 * papb * RLLR +
       d33 * LR * m1s * pap2 * pbp1 * RLLR -
       2.0 * (d23 + d33) * LR * pap1 * pap2 * pbp1 * RLLR +
       d33 * LR * m1s * pap1 * pbp2 * RLLR -
       2.0 * (d23 + d33) * LR * std::pow(pap1, 2) * pbp2 * RLLR +
       2.0 * (d13 + d23 + d33) * LR * m1s * papb * pbp2 * RLLR -
       4.0 * (d13 + d23 + d33) * LR * pap1 * pbp1 * pbp2 * RLLR +
       d33 * m1 * m1s * m2 * papb * RL * RLRL -
       2.0 * d33 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
       d33 * LR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) * RLRR +
       2.0 * d00 * (RL * (2.0 * (LRRL * pap1 * pbp2 + pap2 * pbp1 * RLLL) +
                          m1 * m2 * papb * (LRLL + RLRL)) +
                    LR * (2.0 * (LRRR * pap2 * pbp1 + pap1 * pbp2 * RLLR) +
                          m1 * m2 * papb * (LRLR + RLRR))) +
       d3 * (RL * (-(LRRL * m1s * p1p2 * papb) - LLLL * m1 * ml4 * p1p2 * papb +
                   2.0 * LRRL * p1p2 * pap1 * papb + LRRL * m1s * pap2 * pbp1 +
                   LLLL * m1 * ml4 * pap2 * pbp1 -
                   2.0 * LRRL * pap1 * pap2 * pbp1 +
                   LRLL * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) +
                   LLRL * m2 * ml4 * (m1s * papb - 2.0 * pap1 * pbp1) +
                   LRRL * m1s * pap1 * pbp2 + LLLL * m1 * ml4 * pap1 * pbp2 -
                   2.0 * LRRL * std::pow(pap1, 2) * pbp2 +
                   2.0 * LRRL * m1s * papb * pbp2 -
                   4.0 * LRRL * pap1 * pbp1 * pbp2 - m1s * p1p2 * papb * RLLL +
                   2.0 * m1s * pap2 * papb * RLLL + m1s * pap2 * pbp1 * RLLL -
                   4.0 * pap1 * pap2 * pbp1 * RLLL +
                   2.0 * p1p2 * papb * pbp1 * RLLL -
                   2.0 * pap2 * std::pow(pbp1, 2) * RLLL +
                   m1s * pap1 * pbp2 * RLLL - 2.0 * pap1 * pbp1 * pbp2 * RLLL +
                   m1 * m1s * m2 * papb * RLRL -
                   2.0 * m1 * m2 * pap1 * pbp1 * RLRL +
                   m1s * m2 * ml4 * papb * RRLL -
                   2.0 * m2 * ml4 * pap1 * pbp1 * RRLL -
                   m1 * ml4 * p1p2 * papb * RRRL +
                   m1 * ml4 * pap2 * pbp1 * RRRL +
                   m1 * ml4 * pap1 * pbp2 * RRRL) +
             LR * (-(LRRR * m1s * p1p2 * papb) - LLLR * m1 * ml4 * p1p2 * papb +
                   2.0 * LRRR * m1s * pap2 * papb + LRRR * m1s * pap2 * pbp1 +
                   LLLR * m1 * ml4 * pap2 * pbp1 -
                   4.0 * LRRR * pap1 * pap2 * pbp1 +
                   2.0 * LRRR * p1p2 * papb * pbp1 -
                   2.0 * LRRR * pap2 * std::pow(pbp1, 2) +
                   LRLR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) +
                   LLRR * m2 * ml4 * (m1s * papb - 2.0 * pap1 * pbp1) +
                   LRRR * m1s * pap1 * pbp2 + LLLR * m1 * ml4 * pap1 * pbp2 -
                   2.0 * LRRR * pap1 * pbp1 * pbp2 - m1s * p1p2 * papb * RLLR +
                   2.0 * p1p2 * pap1 * papb * RLLR + m1s * pap2 * pbp1 * RLLR -
                   2.0 * pap1 * pap2 * pbp1 * RLLR + m1s * pap1 * pbp2 * RLLR -
                   2.0 * std::pow(pap1, 2) * pbp2 * RLLR +
                   2.0 * m1s * papb * pbp2 * RLLR -
                   4.0 * pap1 * pbp1 * pbp2 * RLLR +
                   m1 * m1s * m2 * papb * RLRR -
                   2.0 * m1 * m2 * pap1 * pbp1 * RLRR +
                   m1s * m2 * ml4 * papb * RRLR -
                   2.0 * m2 * ml4 * pap1 * pbp1 * RRLR -
                   m1 * ml4 * p1p2 * papb * RRRR +
                   m1 * ml4 * pap2 * pbp1 * RRRR +
                   m1 * ml4 * pap1 * pbp2 * RRRR)));

  return std::real(me0);
}

double FI::MVtbo2t(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      2.0 * ivt1s2 *
      (d33 * LR * LRLR * m1 * m1s * m2 * papb -
       d0 * LL * LRLL * m1s * m2 * ml2 * papb +
       d0 * LL * LLLL * ml2 * ml4 * p1p2 * papb +
       2.0 * d0 * LL * LRLL * m2 * ml2 * pap1 * papb +
       2.0 * (d2 + d3) * LL * LRLL * m2 * ml2 * pap1 * papb -
       2.0 * d33 * LR * LRLR * m1 * m2 * pap1 * pbp1 -
       d0 * LL * LLLL * ml2 * ml4 * pap2 * pbp1 +
       d0 * LL * LLLL * ml2 * ml4 * pap1 * pbp2 -
       d33 * LRRL * m1s * p1p2 * papb * RL +
       2.0 * (d23 + d33) * LRRL * p1p2 * pap1 * papb * RL +
       d33 * LRRL * m1s * pap2 * pbp1 * RL -
       2.0 * (d23 + d33) * LRRL * pap1 * pap2 * pbp1 * RL +
       d33 * LRRL * m1s * pap1 * pbp2 * RL -
       2.0 * (d23 + d33) * LRRL * std::pow(pap1, 2) * pbp2 * RL +
       2.0 * (d13 + d23 + d33) * LRRL * m1s * papb * pbp2 * RL -
       4.0 * (d13 + d23 + d33) * LRRL * pap1 * pbp1 * pbp2 * RL +
       d0 * LL * m1 * ml2 * p1p2 * papb * RLLL -
       d0 * LL * m1 * ml2 * pap2 * pbp1 * RLLL +
       d0 * LL * m1 * ml2 * pap1 * pbp2 * RLLL -
       2.0 * d0 * LL * m1 * ml2 * papb * pbp2 * RLLL -
       2.0 * (d1 + d2 + d3) * LL * m1 * ml2 * papb * pbp2 * RLLL -
       d33 * LR * m1s * p1p2 * papb * RLLR +
       2.0 * (d23 + d33) * LR * p1p2 * pap1 * papb * RLLR +
       d33 * LR * m1s * pap2 * pbp1 * RLLR -
       2.0 * (d23 + d33) * LR * pap1 * pap2 * pbp1 * RLLR +
       d33 * LR * m1s * pap1 * pbp2 * RLLR -
       2.0 * (d23 + d33) * LR * std::pow(pap1, 2) * pbp2 * RLLR +
       2.0 * (d13 + d23 + d33) * LR * m1s * papb * pbp2 * RLLR -
       4.0 * (d13 + d23 + d33) * LR * pap1 * pbp1 * pbp2 * RLLR +
       d33 * m1 * m1s * m2 * papb * RL * RLRL -
       2.0 * d33 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
       2.0 * d00 *
           (LR * LRLR * m1 * m2 * papb + 2.0 * LRRL * pap1 * pbp2 * RL +
            2.0 * LR * pap1 * pbp2 * RLLR + m1 * m2 * papb * RL * RLRL) -
       d0 * LLRR * m1 * m2 * ml2 * ml4 * papb * RR +
       d0 * LRRR * m1 * ml2 * p1p2 * papb * RR -
       d0 * LRRR * m1 * ml2 * pap2 * pbp1 * RR +
       d0 * LRRR * m1 * ml2 * pap1 * pbp2 * RR -
       2.0 * d0 * LRRR * m1 * ml2 * papb * pbp2 * RR -
       2.0 * (d1 + d2 + d3) * LRRR * m1 * ml2 * papb * pbp2 * RR -
       d0 * m1s * m2 * ml2 * papb * RLRR * RR +
       2.0 * d0 * m2 * ml2 * pap1 * papb * RLRR * RR +
       2.0 * (d2 + d3) * m2 * ml2 * pap1 * papb * RLRR * RR -
       d0 * LL * m1 * m2 * ml2 * ml4 * papb * RRLL +
       d3 *
           (LLRL * m1s * m2 * ml4 * papb * RL - LRRL * m1s * p1p2 * papb * RL +
            2.0 * LRRL * p1p2 * pap1 * papb * RL -
            2.0 * LLRL * m2 * ml4 * pap1 * pbp1 * RL +
            LRRL * m1s * pap2 * pbp1 * RL -
            2.0 * LRRL * pap1 * pap2 * pbp1 * RL +
            LRRL * m1s * pap1 * pbp2 * RL -
            2.0 * LRRL * std::pow(pap1, 2) * pbp2 * RL +
            2.0 * LRRL * m1s * papb * pbp2 * RL -
            4.0 * LRRL * pap1 * pbp1 * pbp2 * RL +
            LL * ml2 * (-(LRLL * m1s * m2 * papb) +
                        m1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLLL) +
            m1 * m1s * m2 * papb * RL * RLRL -
            2.0 * m1 * m2 * pap1 * pbp1 * RL * RLRL +
            LRRR * m1 * ml2 * p1p2 * papb * RR -
            LRRR * m1 * ml2 * pap2 * pbp1 * RR +
            LRRR * m1 * ml2 * pap1 * pbp2 * RR -
            m1s * m2 * ml2 * papb * RLRR * RR +
            LR * (LRLR * m1 * m2 * (m1s * papb - 2.0 * pap1 * pbp1) +
                  LLLR * m1 * ml4 *
                      (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) -
                  ((m1s - 2.0 * pap1) * (p1p2 * papb - pap2 * pbp1) -
                   (m1s * (pap1 + 2.0 * papb) -
                    2.0 * pap1 * (pap1 + 2.0 * pbp1)) *
                       pbp2) *
                      RLLR +
                  m2 * ml4 * (m1s * papb - 2.0 * pap1 * pbp1) * RRLR) -
            m1 * ml4 * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRL) +
       d0 * ml2 * ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RR * RRRR);

  return std::real(me0);
}

double FI::MVtbo2u(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m1s, m2s, 2.0 * papb, m1s - 2.0 * pap1,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      -2.0 * ivu2s2 *
      (-(d33 * LR * LRLR * m1s * p1p2 * papb) +
       2.0 * (d23 + d33) * LR * LRLR * m1s * pap2 * papb +
       4.0 * d00 * LR * LRLR * pap2 * pbp1 +
       d33 * LR * LRLR * m1s * pap2 * pbp1 -
       4.0 * (d23 + d33) * LR * LRLR * pap1 * pap2 * pbp1 -
       2.0 * (d1 + d2 + d3) * LL * LRLL * m2 * ml2 * papb * pbp1 +
       2.0 * (d13 + d23 + d33) * LR * LRLR * p1p2 * papb * pbp1 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap2 * std::pow(pbp1, 2) +
       d33 * LR * LRLR * m1s * pap1 * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap1 * pbp1 * pbp2 +
       2.0 * d00 * LRRL * m1 * m2 * papb * RL +
       d33 * LRRL * m1 * m1s * m2 * papb * RL -
       2.0 * d33 * LRRL * m1 * m2 * pap1 * pbp1 * RL +
       2.0 * (d2 + d3) * LL * m1 * ml2 * pap2 * papb * RLLL +
       2.0 * d00 * LR * m1 * m2 * papb * RLLR +
       d33 * LR * m1 * m1s * m2 * papb * RLLR -
       2.0 * d33 * LR * m1 * m2 * pap1 * pbp1 * RLLR -
       d33 * m1s * p1p2 * papb * RL * RLRL +
       2.0 * (d23 + d33) * m1s * pap2 * papb * RL * RLRL +
       4.0 * d00 * pap2 * pbp1 * RL * RLRL +
       d33 * m1s * pap2 * pbp1 * RL * RLRL -
       4.0 * (d23 + d33) * pap1 * pap2 * pbp1 * RL * RLRL +
       2.0 * (d13 + d23 + d33) * p1p2 * papb * pbp1 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap2 * std::pow(pbp1, 2) * RL * RLRL +
       d33 * m1s * pap1 * pbp2 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap1 * pbp1 * pbp2 * RL * RLRL +
       2.0 * ml2 * papb *
           ((d2 + d3) * LRRR * m1 * pap2 - (d1 + d2 + d3) * m2 * pbp1 * RLRR) *
           RR +
       d3 *
           (-(LR * LRLR * m1s * p1p2 * papb) +
            2.0 * LR * LRLR * m1s * pap2 * papb +
            LR * LRLR * m1s * pap2 * pbp1 -
            4.0 * LR * LRLR * pap1 * pap2 * pbp1 +
            2.0 * LR * LRLR * p1p2 * papb * pbp1 -
            2.0 * LR * LRLR * pap2 * std::pow(pbp1, 2) +
            LLLR * LR * m2 * ml4 * (m1s * papb - 2.0 * pap1 * pbp1) +
            LR * LRLR * m1s * pap1 * pbp2 -
            2.0 * LR * LRLR * pap1 * pbp1 * pbp2 +
            LRRL * m1 * m1s * m2 * papb * RL -
            LLRL * m1 * ml4 * p1p2 * papb * RL -
            2.0 * LRRL * m1 * m2 * pap1 * pbp1 * RL +
            LLRL * m1 * ml4 * pap2 * pbp1 * RL +
            LLRL * m1 * ml4 * pap1 * pbp2 * RL +
            LL * ml2 * (LRLL * m1s * m2 * papb -
                        m1 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL) +
            LR * m1 * m1s * m2 * papb * RLLR -
            2.0 * LR * m1 * m2 * pap1 * pbp1 * RLLR -
            m1s * p1p2 * papb * RL * RLRL +
            2.0 * m1s * pap2 * papb * RL * RLRL +
            m1s * pap2 * pbp1 * RL * RLRL -
            4.0 * pap1 * pap2 * pbp1 * RL * RLRL +
            2.0 * p1p2 * papb * pbp1 * RL * RLRL -
            2.0 * pap2 * std::pow(pbp1, 2) * RL * RLRL +
            m1s * pap1 * pbp2 * RL * RLRL -
            2.0 * pap1 * pbp1 * pbp2 * RL * RLRL -
            LRRR * m1 * ml2 * p1p2 * papb * RR -
            LRRR * m1 * ml2 * pap2 * pbp1 * RR +
            LRRR * m1 * ml2 * pap1 * pbp2 * RR +
            m1s * m2 * ml2 * papb * RLRR * RR -
            LR * m1 * ml4 * p1p2 * papb * RRLR +
            LR * m1 * ml4 * pap2 * pbp1 * RRLR +
            LR * m1 * ml4 * pap1 * pbp2 * RRLR +
            m2 * ml4 * (m1s * papb - 2.0 * pap1 * pbp1) * RL * RRRL) +
       d0 * ml2 *
           (LL * (LRLL * m2 * papb * (m1s - 2.0 * pbp1) -
                  LLLL * ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) +
                  m1 * (-(p1p2 * papb * RLLL) + 2.0 * pap2 * papb * RLLL -
                        pap2 * pbp1 * RLLL + pap1 * pbp2 * RLLL +
                        m2 * ml4 * papb * RRLL)) +
            RR * (LLRR * m1 * m2 * ml4 * papb +
                  LRRR * m1 * (-(p1p2 * papb) + 2.0 * pap2 * papb -
                               pap2 * pbp1 + pap1 * pbp2) +
                  m2 * papb * (m1s - 2.0 * pbp1) * RLRR -
                  ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRRR)));

  return std::real(me0);
}

double FI::MVubo2s(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      -4.0 * ivs3v2 *
      (d33 * LR * LRRR * m1 * m2 * m2s * papb -
       d33 * LR * LRLR * m2s * p1p2 * papb +
       2.0 * (d23 + d33) * LR * LRLR * m2s * pap1 * papb +
       d33 * LR * LRLR * m2s * pap2 * pbp1 +
       d33 * LR * LRLR * m2s * pap1 * pbp2 -
       2.0 * d33 * LR * LRRR * m1 * m2 * pap2 * pbp2 -
       4.0 * (d23 + d33) * LR * LRLR * pap1 * pap2 * pbp2 +
       2.0 * (d13 + d23 + d33) * LR * LRLR * p1p2 * papb * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap2 * pbp1 * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap1 * std::pow(pbp2, 2) +
       d33 * LRRL * m1 * m2 * m2s * papb * RL -
       d33 * LRLL * m2s * p1p2 * papb * RL +
       2.0 * (d23 + d33) * LRLL * p1p2 * pap2 * papb * RL +
       d33 * LRLL * m2s * pap2 * pbp1 * RL -
       2.0 * (d23 + d33) * LRLL * std::pow(pap2, 2) * pbp1 * RL +
       2.0 * (d13 + d23 + d33) * LRLL * m2s * papb * pbp1 * RL +
       d33 * LRLL * m2s * pap1 * pbp2 * RL -
       2.0 * d33 * LRRL * m1 * m2 * pap2 * pbp2 * RL -
       2.0 * (d23 + d33) * LRLL * pap1 * pap2 * pbp2 * RL -
       4.0 * (d13 + d23 + d33) * LRLL * pap2 * pbp1 * pbp2 * RL +
       d33 * m1 * m2 * m2s * papb * RL * RLLL -
       2.0 * d33 * m1 * m2 * pap2 * pbp2 * RL * RLLL +
       d33 * LR * m1 * m2 * m2s * papb * RLLR -
       2.0 * d33 * LR * m1 * m2 * pap2 * pbp2 * RLLR -
       d33 * m2s * p1p2 * papb * RL * RLRL +
       2.0 * (d23 + d33) * m2s * pap1 * papb * RL * RLRL +
       d33 * m2s * pap2 * pbp1 * RL * RLRL +
       d33 * m2s * pap1 * pbp2 * RL * RLRL -
       4.0 * (d23 + d33) * pap1 * pap2 * pbp2 * RL * RLRL +
       2.0 * (d13 + d23 + d33) * p1p2 * papb * pbp2 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap2 * pbp1 * pbp2 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap1 * std::pow(pbp2, 2) * RL * RLRL +
       LR * (d33 * m2s * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
             2.0 * (-((d23 + d33) * pap2 *
                      (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) +
                    (d13 + d23 + d33) * pbp1 *
                        (m2s * papb - 2.0 * pap2 * pbp2))) *
           RLRR +
       2.0 * d00 * (RL * (m1 * m2 * papb * (LRRL + RLLL) +
                          2.0 * (LRLL * pap2 * pbp1 + pap1 * pbp2 * RLRL)) +
                    LR * (m1 * m2 * papb * (LRRR + RLLR) +
                          2.0 * (LRLR * pap1 * pbp2 + pap2 * pbp1 * RLRR))) +
       d3 * (RL * (-(LRLL * m2s * p1p2 * papb) - LLRL * m2 * ml4 * p1p2 * papb +
                   2.0 * LRLL * p1p2 * pap2 * papb + LRLL * m2s * pap2 * pbp1 +
                   LLRL * m2 * ml4 * pap2 * pbp1 -
                   2.0 * LRLL * std::pow(pap2, 2) * pbp1 +
                   2.0 * LRLL * m2s * papb * pbp1 + LRLL * m2s * pap1 * pbp2 +
                   LLRL * m2 * ml4 * pap1 * pbp2 -
                   2.0 * LRLL * pap1 * pap2 * pbp2 -
                   4.0 * LRLL * pap2 * pbp1 * pbp2 +
                   LRRL * m1 * m2 * (m2s * papb - 2.0 * pap2 * pbp2) +
                   LLLL * m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) +
                   m1 * m2 * m2s * papb * RLLL -
                   2.0 * m1 * m2 * pap2 * pbp2 * RLLL -
                   m2s * p1p2 * papb * RLRL + 2.0 * m2s * pap1 * papb * RLRL +
                   m2s * pap2 * pbp1 * RLRL + m2s * pap1 * pbp2 * RLRL -
                   4.0 * pap1 * pap2 * pbp2 * RLRL +
                   2.0 * p1p2 * papb * pbp2 * RLRL -
                   2.0 * pap2 * pbp1 * pbp2 * RLRL -
                   2.0 * pap1 * std::pow(pbp2, 2) * RLRL -
                   m2 * ml4 * p1p2 * papb * RRLL +
                   m2 * ml4 * pap2 * pbp1 * RRLL +
                   m2 * ml4 * pap1 * pbp2 * RRLL +
                   m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) * RRRL) +
             LR * (-(LRLR * m2s * p1p2 * papb) - LLRR * m2 * ml4 * p1p2 * papb +
                   2.0 * LRLR * m2s * pap1 * papb + LRLR * m2s * pap2 * pbp1 +
                   LLRR * m2 * ml4 * pap2 * pbp1 + LRLR * m2s * pap1 * pbp2 +
                   LLRR * m2 * ml4 * pap1 * pbp2 -
                   4.0 * LRLR * pap1 * pap2 * pbp2 +
                   2.0 * LRLR * p1p2 * papb * pbp2 -
                   2.0 * LRLR * pap2 * pbp1 * pbp2 -
                   2.0 * LRLR * pap1 * std::pow(pbp2, 2) +
                   LRRR * m1 * m2 * (m2s * papb - 2.0 * pap2 * pbp2) +
                   LLLR * m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) +
                   m1 * m2 * m2s * papb * RLLR -
                   2.0 * m1 * m2 * pap2 * pbp2 * RLLR -
                   m2s * p1p2 * papb * RLRR + 2.0 * p1p2 * pap2 * papb * RLRR +
                   m2s * pap2 * pbp1 * RLRR -
                   2.0 * std::pow(pap2, 2) * pbp1 * RLRR +
                   2.0 * m2s * papb * pbp1 * RLRR + m2s * pap1 * pbp2 * RLRR -
                   2.0 * pap1 * pap2 * pbp2 * RLRR -
                   4.0 * pap2 * pbp1 * pbp2 * RLRR -
                   m2 * ml4 * p1p2 * papb * RRLR +
                   m2 * ml4 * pap2 * pbp1 * RRLR +
                   m2 * ml4 * pap1 * pbp2 * RRLR +
                   m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) * RRRR)));

  return std::real(me0);
}

double FI::MVubo2t(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      -2.0 * ivt1s2 *
      (-(d33 * LR * LRLR * m2s * p1p2 * papb) +
       2.0 * (d23 + d33) * LR * LRLR * m2s * pap1 * papb +
       d33 * LR * LRLR * m2s * pap2 * pbp1 +
       4.0 * d00 * LR * LRLR * pap1 * pbp2 +
       d33 * LR * LRLR * m2s * pap1 * pbp2 -
       4.0 * (d23 + d33) * LR * LRLR * pap1 * pap2 * pbp2 -
       2.0 * (d1 + d2 + d3) * LL * LRLL * m1 * ml2 * papb * pbp2 +
       2.0 * (d13 + d23 + d33) * LR * LRLR * p1p2 * papb * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap2 * pbp1 * pbp2 -
       2.0 * (d13 + d23 + d33) * LR * LRLR * pap1 * std::pow(pbp2, 2) +
       2.0 * d00 * LRRL * m1 * m2 * papb * RL +
       d33 * LRRL * m1 * m2 * m2s * papb * RL -
       2.0 * d33 * LRRL * m1 * m2 * pap2 * pbp2 * RL +
       2.0 * (d2 + d3) * LL * m2 * ml2 * pap1 * papb * RLLL +
       2.0 * d00 * LR * m1 * m2 * papb * RLLR +
       d33 * LR * m1 * m2 * m2s * papb * RLLR -
       2.0 * d33 * LR * m1 * m2 * pap2 * pbp2 * RLLR -
       d33 * m2s * p1p2 * papb * RL * RLRL +
       2.0 * (d23 + d33) * m2s * pap1 * papb * RL * RLRL +
       d33 * m2s * pap2 * pbp1 * RL * RLRL +
       4.0 * d00 * pap1 * pbp2 * RL * RLRL +
       d33 * m2s * pap1 * pbp2 * RL * RLRL -
       4.0 * (d23 + d33) * pap1 * pap2 * pbp2 * RL * RLRL +
       2.0 * (d13 + d23 + d33) * p1p2 * papb * pbp2 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap2 * pbp1 * pbp2 * RL * RLRL -
       2.0 * (d13 + d23 + d33) * pap1 * std::pow(pbp2, 2) * RL * RLRL +
       2.0 * ml2 * papb *
           ((d2 + d3) * LRRR * m2 * pap1 - (d1 + d2 + d3) * m1 * pbp2 * RLRR) *
           RR +
       d3 *
           (-(LR * LRLR * m2s * p1p2 * papb) +
            2.0 * LR * LRLR * m2s * pap1 * papb +
            LR * LRLR * m2s * pap2 * pbp1 + LR * LRLR * m2s * pap1 * pbp2 -
            4.0 * LR * LRLR * pap1 * pap2 * pbp2 +
            2.0 * LR * LRLR * p1p2 * papb * pbp2 -
            2.0 * LR * LRLR * pap2 * pbp1 * pbp2 -
            2.0 * LR * LRLR * pap1 * std::pow(pbp2, 2) +
            LLLR * LR * m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) +
            LRRL * m1 * m2 * m2s * papb * RL -
            LLRL * m2 * ml4 * p1p2 * papb * RL +
            LLRL * m2 * ml4 * pap2 * pbp1 * RL +
            LLRL * m2 * ml4 * pap1 * pbp2 * RL -
            2.0 * LRRL * m1 * m2 * pap2 * pbp2 * RL +
            LL * ml2 * (LRLL * m1 * m2s * papb -
                        m2 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLLL) +
            LR * m1 * m2 * m2s * papb * RLLR -
            2.0 * LR * m1 * m2 * pap2 * pbp2 * RLLR -
            m2s * p1p2 * papb * RL * RLRL +
            2.0 * m2s * pap1 * papb * RL * RLRL +
            m2s * pap2 * pbp1 * RL * RLRL + m2s * pap1 * pbp2 * RL * RLRL -
            4.0 * pap1 * pap2 * pbp2 * RL * RLRL +
            2.0 * p1p2 * papb * pbp2 * RL * RLRL -
            2.0 * pap2 * pbp1 * pbp2 * RL * RLRL -
            2.0 * pap1 * std::pow(pbp2, 2) * RL * RLRL -
            LRRR * m2 * ml2 * p1p2 * papb * RR +
            LRRR * m2 * ml2 * pap2 * pbp1 * RR -
            LRRR * m2 * ml2 * pap1 * pbp2 * RR +
            m1 * m2s * ml2 * papb * RLRR * RR -
            LR * m2 * ml4 * p1p2 * papb * RRLR +
            LR * m2 * ml4 * pap2 * pbp1 * RRLR +
            LR * m2 * ml4 * pap1 * pbp2 * RRLR +
            m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) * RL * RRRL) +
       d0 * ml2 *
           (LL * (LRLL * m1 * papb * (m2s - 2.0 * pbp2) -
                  LLLL * ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                  m2 * (-(p1p2 * papb * RLLL) + 2.0 * pap1 * papb * RLLL +
                        pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL +
                        m1 * ml4 * papb * RRLL)) +
            RR * (LLRR * m1 * m2 * ml4 * papb +
                  LRRR * m2 * (-(p1p2 * papb) + 2.0 * pap1 * papb +
                               pap2 * pbp1 - pap1 * pbp2) +
                  m1 * papb * (m2s - 2.0 * pbp2) * RLRR -
                  ml4 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRR)));

  return std::real(me0);
}

double FI::MVubo2u(double p1s, double p2s, double p3s, double p4s, double pas,
                   double pbs, double ml1s, double ml2s, double ml3s,
                   double ml4s, int ieps) {
  // Npf * dd = new Npf(0.0, 0.0, m2s, m1s, 2.0 * papb, m2s - 2.0 * pap2,
  // m1l, m2l, m3l, m4l, mul, ieps);

  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  double ml4 = sqrt(ml4s);
  complex<double> d0;
  complex<double> d1;
  complex<double> d2;
  complex<double> d3;

  complex<double> d00;
  complex<double> d11;
  complex<double> d22;
  complex<double> d33;
  complex<double> d12;
  complex<double> d13;
  complex<double> d23;

  SetD(p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s, ieps, &d0, &d1,
       &d2, &d3, &d00, &d11, &d22, &d33, &d12, &d13, &d23);
  std::complex<double> me0 =
      2.0 * ivu2s2 *
      (d33 * LR * LRLR * m1 * m2 * m2s * papb -
       d0 * LL * LRLL * m1 * m2s * ml2 * papb +
       d0 * LL * LLLL * ml2 * ml4 * p1p2 * papb +
       2.0 * d0 * LL * LRLL * m1 * ml2 * pap2 * papb +
       2.0 * (d2 + d3) * LL * LRLL * m1 * ml2 * pap2 * papb +
       d0 * LL * LLLL * ml2 * ml4 * pap2 * pbp1 -
       d0 * LL * LLLL * ml2 * ml4 * pap1 * pbp2 -
       2.0 * d33 * LR * LRLR * m1 * m2 * pap2 * pbp2 -
       d33 * LRRL * m2s * p1p2 * papb * RL +
       2.0 * (d23 + d33) * LRRL * p1p2 * pap2 * papb * RL +
       d33 * LRRL * m2s * pap2 * pbp1 * RL -
       2.0 * (d23 + d33) * LRRL * std::pow(pap2, 2) * pbp1 * RL +
       2.0 * (d13 + d23 + d33) * LRRL * m2s * papb * pbp1 * RL +
       d33 * LRRL * m2s * pap1 * pbp2 * RL -
       2.0 * (d23 + d33) * LRRL * pap1 * pap2 * pbp2 * RL -
       4.0 * (d13 + d23 + d33) * LRRL * pap2 * pbp1 * pbp2 * RL +
       d0 * LL * m2 * ml2 * p1p2 * papb * RLLL +
       d0 * LL * m2 * ml2 * pap2 * pbp1 * RLLL -
       2.0 * d0 * LL * m2 * ml2 * papb * pbp1 * RLLL -
       2.0 * (d1 + d2 + d3) * LL * m2 * ml2 * papb * pbp1 * RLLL -
       d0 * LL * m2 * ml2 * pap1 * pbp2 * RLLL -
       d33 * LR * m2s * p1p2 * papb * RLLR +
       2.0 * (d23 + d33) * LR * p1p2 * pap2 * papb * RLLR +
       d33 * LR * m2s * pap2 * pbp1 * RLLR -
       2.0 * (d23 + d33) * LR * std::pow(pap2, 2) * pbp1 * RLLR +
       2.0 * (d13 + d23 + d33) * LR * m2s * papb * pbp1 * RLLR +
       d33 * LR * m2s * pap1 * pbp2 * RLLR -
       2.0 * (d23 + d33) * LR * pap1 * pap2 * pbp2 * RLLR -
       4.0 * (d13 + d23 + d33) * LR * pap2 * pbp1 * pbp2 * RLLR +
       d33 * m1 * m2 * m2s * papb * RL * RLRL -
       2.0 * d33 * m1 * m2 * pap2 * pbp2 * RL * RLRL +
       2.0 * d00 *
           (LR * LRLR * m1 * m2 * papb + 2.0 * LRRL * pap2 * pbp1 * RL +
            2.0 * LR * pap2 * pbp1 * RLLR + m1 * m2 * papb * RL * RLRL) -
       d0 * LLRR * m1 * m2 * ml2 * ml4 * papb * RR +
       d0 * LRRR * m2 * ml2 * p1p2 * papb * RR +
       d0 * LRRR * m2 * ml2 * pap2 * pbp1 * RR -
       2.0 * d0 * LRRR * m2 * ml2 * papb * pbp1 * RR -
       2.0 * (d1 + d2 + d3) * LRRR * m2 * ml2 * papb * pbp1 * RR -
       d0 * LRRR * m2 * ml2 * pap1 * pbp2 * RR -
       d0 * m1 * m2s * ml2 * papb * RLRR * RR +
       2.0 * d0 * m1 * ml2 * pap2 * papb * RLRR * RR +
       2.0 * (d2 + d3) * m1 * ml2 * pap2 * papb * RLRR * RR -
       d0 * LL * m1 * m2 * ml2 * ml4 * papb * RRLL +
       d3 *
           (LLRL * m1 * m2s * ml4 * papb * RL - LRRL * m2s * p1p2 * papb * RL +
            2.0 * LRRL * p1p2 * pap2 * papb * RL +
            LRRL * m2s * pap2 * pbp1 * RL -
            2.0 * LRRL * std::pow(pap2, 2) * pbp1 * RL +
            2.0 * LRRL * m2s * papb * pbp1 * RL +
            LRRL * m2s * pap1 * pbp2 * RL -
            2.0 * LLRL * m1 * ml4 * pap2 * pbp2 * RL -
            2.0 * LRRL * pap1 * pap2 * pbp2 * RL -
            4.0 * LRRL * pap2 * pbp1 * pbp2 * RL +
            LL * ml2 * (-(LRLL * m1 * m2s * papb) +
                        m2 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RLLL) +
            m1 * m2 * m2s * papb * RL * RLRL -
            2.0 * m1 * m2 * pap2 * pbp2 * RL * RLRL +
            LRRR * m2 * ml2 * p1p2 * papb * RR +
            LRRR * m2 * ml2 * pap2 * pbp1 * RR -
            LRRR * m2 * ml2 * pap1 * pbp2 * RR -
            m1 * m2s * ml2 * papb * RLRR * RR +
            LR * (LLLR * m2 * ml4 *
                      (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                  LRLR * m1 * m2 * (m2s * papb - 2.0 * pap2 * pbp2) -
                  (m2s * (p1p2 * papb - pap2 * pbp1 - 2.0 * papb * pbp1 -
                          pap1 * pbp2) +
                   2.0 * pap2 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2 +
                                 2.0 * pbp1 * pbp2)) *
                      RLLR +
                  m1 * ml4 * (m2s * papb - 2.0 * pap2 * pbp2) * RRLR) -
            m2 * ml4 * (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) * RL * RRRL) +
       d0 * ml2 * ml4 * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RR * RRRR);

  return std::real(me0);
}
