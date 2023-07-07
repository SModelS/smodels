#include <complex>
#include <iostream>
//#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "utils.h"
using namespace std;

#define MPIs pow2(M_PI)

// sleptons
double FI::MVsbu2sSL(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = b1 * (8.0 * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                              (2.0 * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  return real(me0);
}

double FI::MVsbu1sSL(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = b1 * (8.0 * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                              (2.0 * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  return real(me0);
}

// fermions
// Process: q + \bar{q} -> f + \bar{f}.

double FI::MVsbu1s(double p1s, double ml1s, double ml2s, int ieps) {

  // General idea: M_final = M B  (some squared matrix element times a
  // complicated divergent integral)
  // M = M0 + eps M1 + eps^2 M2 + O(eps^3)
  // B = 1/eps^2 B2 + 1/eps B1 + B0 + O(eps)
  // -> M_final = B0 M0 + B1 M1 + B2 M2 + terms going to zero + divergent terms
  // which we absorb
  // Don't get confused. Number refers to the order of epsilon and not to a
  // certain B-function!

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  // Part of scalar integral of O(1/eps^0) times squared matrix element of
  // O(eps^0).
  complex<double> me0 =
      b1 * ivs3v1 * ivs3v2 *
      (16.0 * pap2 * pbp1 * RL * LLLL + 16.0 * pap2 * pbp1 * LR * RRRR +
       16.0 * pap1 * pbp2 * RL * RLRL + 16.0 * pap1 * pbp2 * LR * LRLR +
       8.0 * m1 * m2 * papb * RL * RLLL + 8.0 * m1 * m2 * papb * RL * LLRL +
       8.0 * m1 * m2 * papb * LR * LRRR + 8 * m1 * m2 * papb * LR * RRLR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11); // get two-point function prop. to O(1/epsilon^1).

  // Part of scalar integral of O(1/eps^1) times squared matrix element of
  // O(eps^1).
  complex<double> me1 =
      b1 * ivs3v1 * ivs3v2 *
      (24.0 * pap2 * pbp1 * RL * RLRL - 24.0 * pap2 * pbp1 * RL * LLLL -
       24.0 * pap2 * pbp1 * LR * RRRR + 24.0 * pap2 * pbp1 * LR * LRLR -
       24.0 * pap1 * pbp2 * RL * RLRL + 24.0 * pap1 * pbp2 * RL * LLLL +
       24.0 * pap1 * pbp2 * LR * RRRR - 24.0 * pap1 * pbp2 * LR * LRLR -
       8.0 * papb * p1p2 * RL * RLRL - 8.0 * papb * p1p2 * RL * LLLL -
       8.0 * papb * p1p2 * LR * RRRR - 8.0 * papb * p1p2 * LR * LRLR -
       8.0 * m1 * m2 * papb * RL * RLLL - 8.0 * m1 * m2 * papb * RL * LLRL -
       8.0 * m1 * m2 * papb * LR * LRRR - 8 * m1 * m2 * papb * LR * RRLR);

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  // Part of scalar integral of O(1/eps^1) times squared matrix element of
  // O(eps^2).
  complex<double> me2 =
      b1 * ivs3v1 * ivs3v2 *
      (-16.0 * pap2 * pbp1 * RL * RLRL + 16.0 * pap2 * pbp1 * RL * LLLL +
       16.0 * pap2 * pbp1 * LR * RRRR - 16.0 * pap2 * pbp1 * LR * LRLR +
       16.0 * pap1 * pbp2 * RL * RLRL - 16.0 * pap1 * pbp2 * RL * LLLL -
       16.0 * pap1 * pbp2 * LR * RRRR + 16.0 * pap1 * pbp2 * LR * LRLR);

  return real(
      me0 + me1 +
      me2); // returns the finite result of the whole squared matrix element.
}

double FI::MVsbu1t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v1 * ivt1s2 *
      (8.0 * pap1 * pbp2 * RL * RLRL + 8.0 * pap1 * pbp2 * LR * LRLR +
       4.0 * m1 * m2 * papb * RL * LLRL + 4.0 * m1 * m2 * papb * LR * RRLR);

  me0 += b0 * ivs3v1 * ivt1s2 *
         (-4. * m2 * pbp1 * ml1 * RR * LLRR - 4. * m2 * pbp1 * ml1 * LL * RRLL -
          8. * m1 * pbp2 * ml1 * RR * RLRR - 8. * m1 * pbp2 * ml1 * LL * LRLL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      b1 * ivs3v1 * ivt1s2 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * LR * LRLR -
       4.0 * pap1 * pbp2 * RL * RLRL - 4.0 * pap1 * pbp2 * LR * LRLR -
       4.0 * papb * p1p2 * RL * RLRL - 4.0 * papb * p1p2 * LR * LRLR -
       4.0 * m1 * m2 * papb * RL * LLRL - 4.0 * m1 * m2 * papb * LR * RRLR);

  me1 += b0 * ivs3v1 * ivt1s2 *
         (4. * m2 * pbp1 * ml1 * RR * LLRR + 4. * m2 * pbp1 * ml1 * LL * RRLL +
          4. * m1 * pbp2 * ml1 * RR * RLRR + 4. * m1 * pbp2 * ml1 * LL * LRLL);

  return real(me0 + me1);
}

double FI::MVsbu1u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v1 * ivu2s2 *
      (-8.0 * pap2 * pbp1 * RL * LLRL - 8.0 * pap2 * pbp1 * LR * RRLR -
       4.0 * m1 * m2 * papb * RL * RLRL - 4.0 * m1 * m2 * papb * LR * LRLR);

  me0 += b0 * ivs3v1 * ivu2s2 *
         (8. * m2 * pbp1 * ml1 * RR * LLRR + 8. * m2 * pbp1 * ml1 * LL * RRLL +
          4. * m1 * pbp2 * ml1 * RR * RLRR + 4. * m1 * pbp2 * ml1 * LL * LRLL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      b1 * ivs3v1 * ivu2s2 *
      (4.0 * pap2 * pbp1 * RL * LLRL + 4.0 * pap2 * pbp1 * LR * RRLR -
       4.0 * pap1 * pbp2 * RL * LLRL - 4.0 * pap1 * pbp2 * LR * RRLR +
       4.0 * papb * p1p2 * RL * LLRL + 4.0 * papb * p1p2 * LR * RRLR +
       4.0 * m1 * m2 * papb * RL * RLRL + 4.0 * m1 * m2 * papb * LR * LRLR);
  me1 += b0 * ivs3v1 * ivu2s2 *
         (-4. * m2 * pbp1 * ml1 * RR * LLRR - 4. * m2 * pbp1 * ml1 * LL * RRLL -
          4. * m1 * pbp2 * ml1 * RR * RLRR - 4. * m1 * pbp2 * ml1 * LL * LRLL);

  return real(me0 + me1);
}

double FI::MVsbu2s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v1 * ivs3v2 *
      (16.0 * pap2 * pbp1 * RL * LLLL + 16.0 * pap2 * pbp1 * LR * RRRR +
       16.0 * pap1 * pbp2 * RL * RLRL + 16.0 * pap1 * pbp2 * LR * LRLR +
       8.0 * m1 * m2 * papb * RL * RLLL + 8.0 * m1 * m2 * papb * RL * LLRL +
       8.0 * m1 * m2 * papb * LR * LRRR + 8 * m1 * m2 * papb * LR * RRLR);
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      b1 * ivs3v1 * ivs3v2 *
      (24.0 * pap2 * pbp1 * RL * RLRL - 24.0 * pap2 * pbp1 * RL * LLLL -
       24.0 * pap2 * pbp1 * LR * RRRR + 24.0 * pap2 * pbp1 * LR * LRLR -
       24.0 * pap1 * pbp2 * RL * RLRL + 24.0 * pap1 * pbp2 * RL * LLLL +
       24.0 * pap1 * pbp2 * LR * RRRR - 24.0 * pap1 * pbp2 * LR * LRLR -
       8.0 * papb * p1p2 * RL * RLRL - 8.0 * papb * p1p2 * RL * LLLL -
       8.0 * papb * p1p2 * LR * RRRR - 8.0 * papb * p1p2 * LR * LRLR -
       8.0 * m1 * m2 * papb * RL * RLLL - 8.0 * m1 * m2 * papb * RL * LLRL -
       8.0 * m1 * m2 * papb * LR * LRRR - 8 * m1 * m2 * papb * LR * RRLR);
  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me2 =
      b1 * ivs3v1 * ivs3v2 *
      (-16.0 * pap2 * pbp1 * RL * RLRL + 16.0 * pap2 * pbp1 * RL * LLLL +
       16.0 * pap2 * pbp1 * LR * RRRR - 16.0 * pap2 * pbp1 * LR * LRLR +
       16.0 * pap1 * pbp2 * RL * RLRL - 16.0 * pap1 * pbp2 * RL * LLLL -
       16.0 * pap1 * pbp2 * LR * RRRR + 16.0 * pap1 * pbp2 * LR * LRLR);
  return real(me0 + me1 + me2);
}

double FI::MVsbu2t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v1 * ivt1s2 *
      (8.0 * pap1 * pbp2 * RL * RLRL + 8.0 * pap1 * pbp2 * LR * LRLR +
       4.0 * m1 * m2 * papb * RL * LLRL + 4.0 * m1 * m2 * papb * LR * RRLR);

  me0 += b0 * ivs3v1 * ivt1s2 *
         (-8. * m2 * pap1 * ml1 * RR * LRRR - 8. * m2 * pap1 * ml1 * LL * RLLL -
          4. * m1 * pap2 * ml1 * RR * RRRR - 4. * m1 * pap2 * ml1 * LL * LLLL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 =
      b1 * ivs3v1 * ivt1s2 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * LR * LRLR -
       4.0 * pap1 * pbp2 * RL * RLRL - 4.0 * pap1 * pbp2 * LR * LRLR -
       4.0 * papb * p1p2 * RL * RLRL - 4.0 * papb * p1p2 * LR * LRLR -
       4.0 * m1 * m2 * papb * RL * LLRL - 4.0 * m1 * m2 * papb * LR * RRLR);

  me1 += b0 * ivs3v1 * ivt1s2 *
         (4. * m2 * pap1 * ml1 * RR * LRRR + 4. * m2 * pap1 * ml1 * LL * RLLL +
          4. * m1 * pap2 * ml1 * RR * RRRR + 4. * m1 * pap2 * ml1 * LL * LLLL);
  return real(me0 + me1);
}

double FI::MVsbu2u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v1 * ivu2s2 *
      (-8.0 * pap2 * pbp1 * RL * LLRL - 8.0 * pap2 * pbp1 * LR * RRLR -
       4.0 * m1 * m2 * papb * RL * RLRL - 4.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivs3v1 * ivu2s2 *
         (4. * m2 * pap1 * ml1 * RR * LRRR + 4. * m2 * pap1 * ml1 * LL * RLLL +
          8. * m1 * pap2 * ml1 * RR * RRRR + 8. * m1 * pap2 * ml1 * LL * LLLL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      b1 * ivs3v1 * ivu2s2 *
      (4.0 * pap2 * pbp1 * RL * LLRL + 4.0 * pap2 * pbp1 * LR * RRLR -
       4.0 * pap1 * pbp2 * RL * LLRL - 4.0 * pap1 * pbp2 * LR * RRLR +
       4.0 * papb * p1p2 * RL * LLRL + 4.0 * papb * p1p2 * LR * RRLR +
       4.0 * m1 * m2 * papb * RL * RLRL + 4.0 * m1 * m2 * papb * LR * LRLR);
  me1 += +b0 * ivs3v1 * ivu2s2 *
         (-4. * m2 * pap1 * ml1 * RR * LRRR - 4. * m2 * pap1 * ml1 * LL * RLLL -
          4. * m1 * pap2 * ml1 * RR * RRRR - 4. * m1 * pap2 * ml1 * LL * LLLL);

  return real(me0 + me1);
}

double FI::MVtbu1s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v2 * ivt1s1 *
      (8.0 * pap1 * pbp2 * RL * RLRL + 8.0 * pap1 * pbp2 * LR * LRLR +
       4.0 * m1 * m2 * papb * RL * RLLL + 4.0 * m1 * m2 * papb * LR * LRRR);

  me0 += +b0 * ivs3v2 * ivt1s1 *
         (-4. * m2 * pbp1 * ml1 * RR * LLRR - 4. * m2 * pbp1 * ml1 * LL * RRLL -
          8. * m1 * pbp2 * ml1 * RR * LLLR - 8. * m1 * pbp2 * ml1 * LL * RRRL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      b1 * ivs3v2 * ivt1s1 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * LR * LRLR -
       4.0 * pap1 * pbp2 * RL * RLRL - 4.0 * pap1 * pbp2 * LR * LRLR -
       4.0 * papb * p1p2 * RL * RLRL - 4.0 * papb * p1p2 * LR * LRLR -
       4.0 * m1 * m2 * papb * RL * RLLL - 4.0 * m1 * m2 * papb * LR * LRRR);

  me1 += +b0 * ivs3v2 * ivt1s1 *
         (4. * m2 * pbp1 * ml1 * RR * LLRR + 4. * m2 * pbp1 * ml1 * LL * RRLL +
          4. * m1 * pbp2 * ml1 * RR * LLLR + 4. * m1 * pbp2 * ml1 * LL * RRRL);
  return real(me0 + me1);
}

double FI::MVtbu1t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivt1s1 * ivt1s2 *
      (4.0 * pap1 * pbp2 * RL * RLRL + 4.0 * pap1 * pbp2 * RL * LLLL +
       4.0 * pap1 * pbp2 * LR * RRRR + 4.0 * pap1 * pbp2 * LR * LRLR);

  me0 += +b0 * ivt1s1 * ivt1s2 *
         (-4. * m1 * pbp2 * ml1 * RR * RLRR - 4. * m1 * pbp2 * ml1 * RR * LLLR -
          4. * m1 * pbp2 * ml1 * LL * RRRL - 4. * m1 * pbp2 * ml1 * LL * LRLL);

  return real(me0);
}

double FI::MVtbu1u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivt1s1 * ivu2s2 *
      (-2 * pap2 * pbp1 * RL * LLLL - 2.0 * pap2 * pbp1 * LR * RRRR -
       2.0 * pap1 * pbp2 * RL * LLLL - 2.0 * pap1 * pbp2 * LR * RRRR +
       2.0 * papb * p1p2 * RL * LLLL + 2.0 * papb * p1p2 * LR * RRRR -
       2.0 * m1 * m2 * papb * RL * RLRL - 2.0 * m1 * m2 * papb * LR * LRLR);

  me0 += b0 * ivt1s1 * ivu2s2 *
         (2. * m2 * pbp1 * ml1 * RR * LLLR + 2. * m2 * pbp1 * ml1 * LL * RRRL +
          2. * m1 * pbp2 * ml1 * RR * RLRR + 2. * m1 * pbp2 * ml1 * LL * LRLL);

  return real(me0);
}

double FI::MVtbu2s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      +b1 * ivs3v2 * ivt1s1 *
      (8.0 * pap1 * pbp2 * RL * RLRL + 8.0 * pap1 * pbp2 * LR * LRLR +
       4.0 * m1 * m2 * papb * RL * RLLL + 4.0 * m1 * m2 * papb * LR * LRRR);

  me0 += b0 * ivs3v2 * ivt1s1 *
         (-8. * m2 * pap1 * ml1 * RR * LLRL - 8. * m2 * pap1 * ml1 * LL * RRLR -
          4. * m1 * pap2 * ml1 * RR * LLLL - 4. * m1 * pap2 * ml1 * LL * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      +b1 * ivs3v2 * ivt1s1 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * LR * LRLR -
       4.0 * pap1 * pbp2 * RL * RLRL - 4.0 * pap1 * pbp2 * LR * LRLR -
       4.0 * papb * p1p2 * RL * RLRL - 4.0 * papb * p1p2 * LR * LRLR -
       4.0 * m1 * m2 * papb * RL * RLLL - 4.0 * m1 * m2 * papb * LR * LRRR);

  me1 += b0 * ivs3v2 * ivt1s1 *
         (4. * m2 * pap1 * ml1 * RR * LLRL + 4. * m2 * pap1 * ml1 * LL * RRLR +
          4. * m1 * pap2 * ml1 * RR * LLLL + 4. * m1 * pap2 * ml1 * LL * RRRR);

  return real(me0 + me1);
}

double FI::MVtbu2t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      +b1 * ivt1s1 * ivt1s2 *
      (4.0 * pap1 * pbp2 * RL * RLRL + 4.0 * pap1 * pbp2 * RL * RRRR +
       4.0 * pap1 * pbp2 * LR * LRLR + 4.0 * pap1 * pbp2 * LR * LLLL);

  me0 += b0 * ivt1s1 * ivt1s2 *
         (-4. * m2 * pap1 * ml1 * RR * LRRR - 4. * m2 * pap1 * ml1 * RR * LLRL -
          4. * m2 * pap1 * ml1 * LL * RRLR - 4. * m2 * pap1 * ml1 * LL * RLLL);

  return real(me0);
}

double FI::MVtbu2u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivt1s1 * ivu2s2 *
      (-2.0 * pap2 * pbp1 * RL * RRRR - 2.0 * pap2 * pbp1 * LR * LLLL -
       2.0 * pap1 * pbp2 * RL * RRRR - 2.0 * pap1 * pbp2 * LR * LLLL +
       2.0 * papb * p1p2 * RL * RRRR + 2.0 * papb * p1p2 * LR * LLLL -
       2.0 * m1 * m2 * papb * RL * RLRL - 2.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivt1s1 * ivu2s2 *
         (2. * m2 * pap1 * ml1 * RR * LRRR + 2. * m2 * pap1 * ml1 * LL * RLLL +
          2. * m1 * pap2 * ml1 * RR * LLRL + 2. * m1 * pap2 * ml1 * LL * RRLR);

  return real(me0);
}

double FI::MVubu1s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v2 * ivu2s1 *
      (-8.0 * pap2 * pbp1 * RL * RLLL - 8.0 * pap2 * pbp1 * LR * LRRR -
       4.0 * m1 * m2 * papb * RL * RLRL - 4.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivs3v2 * ivu2s1 *
         (8. * m2 * pbp1 * ml1 * RR * LLRR + 8. * m2 * pbp1 * ml1 * LL * RRLL +
          4. * m1 * pbp2 * ml1 * RR * LLLR + 4. * m1 * pbp2 * ml1 * LL * RRRL);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      +b1 * ivs3v2 * ivu2s1 *
      (4.0 * pap2 * pbp1 * RL * RLLL + 4.0 * pap2 * pbp1 * LR * LRRR -
       4.0 * pap1 * pbp2 * RL * RLLL - 4.0 * pap1 * pbp2 * LR * LRRR +
       4.0 * papb * p1p2 * RL * RLLL + 4.0 * papb * p1p2 * LR * LRRR +
       4.0 * m1 * m2 * papb * RL * RLRL + 4.0 * m1 * m2 * papb * LR * LRLR);

  me1 += +b0 * ivs3v2 * ivu2s1 *
         (-4. * m2 * pbp1 * ml1 * RR * LLRR - 4. * m2 * pbp1 * ml1 * LL * RRLL -
          4. * m1 * pbp2 * ml1 * RR * LLLR - 4. * m1 * pbp2 * ml1 * LL * RRRL);

  return real(me0 + me1);
}

double FI::MVubu1t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivt1s2 * ivu2s1 *
      (-2.0 * pap2 * pbp1 * RL * LLLL - 2.0 * pap2 * pbp1 * LR * RRRR -
       2.0 * pap1 * pbp2 * RL * LLLL - 2.0 * pap1 * pbp2 * LR * RRRR +
       2.0 * papb * p1p2 * RL * LLLL + 2.0 * papb * p1p2 * LR * RRRR -
       2.0 * m1 * m2 * papb * RL * RLRL - 2.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivt1s2 * ivu2s1 *
         (2. * m2 * pbp1 * ml1 * RR * RLRR + 2. * m2 * pbp1 * ml1 * LL * LRLL +
          2. * m1 * pbp2 * ml1 * RR * LLLR + 2. * m1 * pbp2 * ml1 * LL * RRRL);

  return real(me0);
}

double FI::MVubu1u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      +b1 * ivu2s1 * ivu2s2 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * RL * LLLL +
       4.0 * pap2 * pbp1 * LR * RRRR + 4.0 * pap2 * pbp1 * LR * LRLR);

  me0 += +b0 * ivu2s1 * ivu2s2 *
         (-4. * m2 * pbp1 * ml1 * RR * RLRR - 4. * m2 * pbp1 * ml1 * RR * LLLR -
          4. * m2 * pbp1 * ml1 * LL * RRRL - 4. * m2 * pbp1 * ml1 * LL * LRLL);

  return real(me0);
}

double FI::MVubu2s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivs3v2 * ivu2s1 *
      (-8.0 * pap2 * pbp1 * RL * RLLL - 8.0 * pap2 * pbp1 * LR * LRRR -
       4.0 * m1 * m2 * papb * RL * RLRL - 4.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivs3v2 * ivu2s1 *
         (4. * m2 * pap1 * ml1 * RR * LLRL + 4. * m2 * pap1 * ml1 * LL * RRLR +
          8. * m1 * pap2 * ml1 * RR * LLLL + 8. * m1 * pap2 * ml1 * LL * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      +b1 * ivs3v2 * ivu2s1 *
      (4.0 * pap2 * pbp1 * RL * RLLL + 4.0 * pap2 * pbp1 * LR * LRRR -
       4.0 * pap1 * pbp2 * RL * RLLL - 4.0 * pap1 * pbp2 * LR * LRRR +
       4.0 * papb * p1p2 * RL * RLLL + 4.0 * papb * p1p2 * LR * LRRR +
       4.0 * m1 * m2 * papb * RL * RLRL + 4.0 * m1 * m2 * papb * LR * LRLR);

  me1 += +b0 * ivs3v2 * ivu2s1 *
         (-4. * m2 * pap1 * ml1 * RR * LLRL - 4. * m2 * pap1 * ml1 * LL * RRLR -
          4. * m1 * pap2 * ml1 * RR * LLLL - 4. * m1 * pap2 * ml1 * LL * RRRR);
  return real(me0 + me1);
}

double FI::MVubu2t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      +b1 * ivt1s2 * ivu2s1 *
      (-2.0 * pap2 * pbp1 * RL * RRRR - 2.0 * pap2 * pbp1 * LR * LLLL -
       2.0 * pap1 * pbp2 * RL * RRRR - 2.0 * pap1 * pbp2 * LR * LLLL +
       2.0 * papb * p1p2 * RL * RRRR + 2.0 * papb * p1p2 * LR * LLLL -
       2.0 * m1 * m2 * papb * RL * RLRL - 2.0 * m1 * m2 * papb * LR * LRLR);

  me0 += +b0 * ivt1s2 * ivu2s1 *
         (2. * m2 * pap1 * ml1 * RR * LLRL + 2. * m2 * pap1 * ml1 * LL * RRLR +
          2. * m1 * pap2 * ml1 * RR * LRRR + 2. * m1 * pap2 * ml1 * LL * RLLL);

  return real(me0);
}

double FI::MVubu2u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      b1 * ivu2s1 * ivu2s2 *
      (4.0 * pap2 * pbp1 * RL * RLRL + 4.0 * pap2 * pbp1 * RL * RRRR +
       4.0 * pap2 * pbp1 * LR * LRLR + 4.0 * pap2 * pbp1 * LR * LLLL);

  me0 += +b0 * ivu2s1 * ivu2s2 *
         (-4. * m1 * pap2 * ml1 * RR * LRRR - 4. * m1 * pap2 * ml1 * RR * LLRL -
          4. * m1 * pap2 * ml1 * LL * RRLR - 4. * m1 * pap2 * ml1 * LL * RLLL);

  return real(me0);
}

double FI::MVtbu3s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      4.0 * ivs3v2 * pow(ivt1s1, 2) *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m1s - 2.0 * pap1)) *
      (m1 * m2 * papb * (LRRR + RLLL) + 2.0 * pap1 * pbp2 * (LRLR + RLRL)) * RR;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      -4.0 * ivs3v2 * pow(ivt1s1, 2) *
      ((B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m1s - 2.0 * pap1) *
           (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
            LRLR * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
            pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) +
       b00 * (6.0 * LRRR * m1 * m2 * papb + 4.0 * LRLR * p1p2 * papb -
              4.0 * LRLR * pap2 * pbp1 + 8.0 * LRLR * pap1 * pbp2 +
              6.0 * m1 * m2 * papb * RLLL + 4.0 * p1p2 * papb * RLRL -
              4.0 * pap2 * pbp1 * RLRL + 8.0 * pap1 * pbp2 * RLRL)) *
      RR;

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me2 = 8.0 * b00 * ivs3v2 * pow(ivt1s1, 2) *
                        (LRRR * m1 * m2 * papb +
                         LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                         pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) *
                        RR;

  return real(me0 + me1 + me2);
}

double FI::MVtbu3t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      4.0 * pow(ivt1s1, 2) * ivt1s2 *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m1s - 2.0 * pap1)) *
      pap1 * pbp2 * RR * (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -8.0 * b00 * pow(ivt1s1, 2) * ivt1s2 * pap1 * pbp2 *
                        RR * (LLLL + LRLR + RLRL + RRRR);

  return real(me0 + me1);
}

double FI::MVtbu3u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -2.0 * pow(ivt1s1, 2) * ivu2s2 *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m1s - 2.0 * pap1)) *
      RR * (LRLR * m1 * m2 * papb +
            LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
            m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
            pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = 4.0 * b00 * pow(ivt1s1, 2) * ivu2s2 * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return real(me0 + me1);
}

double FI::MVubu3s(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -4.0 * ivs3v2 * pow(ivu2s1, 2) *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m2s - 2.0 * pap2)) *
      (2.0 * pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) * RR;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      4.0 * ivs3v2 * pow(ivu2s1, 2) *
      ((B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m2s - 2.0 * pap2) *
           (LRLR * m1 * m2 * papb +
            (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
            m1 * m2 * papb * RLRL) +
       b00 * (6.0 * LRLR * m1 * m2 * papb +
              4.0 * (p1p2 * papb + 2.0 * pap2 * pbp1 - pap1 * pbp2) *
                  (LRRR + RLLL) +
              6.0 * m1 * m2 * papb * RLRL)) *
      RR;

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me2 =
      -8.0 * b00 * ivs3v2 * pow(ivu2s1, 2) *
      (LRLR * m1 * m2 * papb +
       (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
       m1 * m2 * papb * RLRL) *
      RR;

  return real(me0 + me1 + me2);
}

double FI::MVubu3t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -2.0 * ivt1s2 * pow(ivu2s1, 2) *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m2s - 2.0 * pap2)) *
      RR * (LRLR * m1 * m2 * papb +
            LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
            m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
            pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = 4.0 * b00 * ivt1s2 * pow(ivu2s1, 2) * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return real(me0 + me1);
}

double FI::MVubu3u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      4.0 * pow(ivu2s1, 2) * ivu2s2 *
      (4.0 * b00 +
       (B0(p1s, ml1s, ml2s) - 2.0 * b1 + b11) * (m2s - 2.0 * pap2)) *
      pap2 * pbp1 * RR * (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -8.0 * b00 * pow(ivu2s1, 2) * ivu2s2 * pap2 * pbp1 *
                        RR * (LLLL + LRLR + RLRL + RRRR);

  return real(me0 + me1);
}

double FI::MVtbu4s(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      8.0 * ivs3v2 * ivt1s1 * ivt2s1 *
      (m1 * m2 * papb * (LRRR + RLLL) + 2.0 * pap1 * pbp2 * (LRLR + RLRL)) *
      ((4.0 * b00 + (b1 + b11) * (m1s - 2.0 * pap1)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR));
  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      -8.0 * ivs3v2 * ivt1s1 * ivt2s1 *
      (2.0 * b00 * (LR + RL) *
           (3.0 * LRRR * m1 * m2 * papb + 2.0 * LRLR * p1p2 * papb -
            2.0 * LRLR * pap2 * pbp1 + 4.0 * LRLR * pap1 * pbp2 +
            3.0 * m1 * m2 * papb * RLLL + 2.0 * p1p2 * papb * RLRL -
            2.0 * pap2 * pbp1 * RLRL + 4.0 * pap1 * pbp2 * RLRL) +
       (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
        LRLR * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
        pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) *
           ((b1 + b11) * (m1s - 2.0 * pap1) * (LR + RL) +
            B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)));

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me2 = 16.0 * b00 * ivs3v2 * ivt1s1 * ivt2s1 * (LR + RL) *
                        (LRRR * m1 * m2 * papb +
                         LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                         pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL);

  return real(me0 + me1 + me2);
}

double FI::MVtbu4t(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      8.0 * ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
      ((4.0 * b00 + (b1 + b11) * (m1s - 2.0 * pap1)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)) *
      (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -16.0 * b00 * ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

  return real(me0 + me1);
}

double FI::MVtbu4u(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -4.0 * ivt1s1 * ivt2s1 * ivu2s2 *
      ((4.0 * b00 + (b1 + b11) * (m1s - 2.0 * pap1)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)) *
      (LRLR * m1 * m2 * papb +
       LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
       m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
       pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = 8.0 * b00 * ivt1s1 * ivt2s1 * ivu2s2 * (LR + RL) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return real(me0 + me1);
}

double FI::MVubu4s(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -8.0 * ivs3v2 * ivu1s1 * ivu2s1 *
      (2.0 * pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) *
      ((4.0 * b00 + (b1 + b11) * (m2s - 2.0 * pap2)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR));

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      8.0 * ivs3v2 * ivu1s1 * ivu2s1 *
      (2.0 * b00 * (LR + RL) *
           (3.0 * LRLR * m1 * m2 * papb +
            2.0 * (p1p2 * papb + 2.0 * pap2 * pbp1 - pap1 * pbp2) *
                (LRRR + RLLL) +
            3.0 * m1 * m2 * papb * RLRL) +
       (LRLR * m1 * m2 * papb +
        (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
        m1 * m2 * papb * RLRL) *
           ((b1 + b11) * (m2s - 2.0 * pap2) * (LR + RL) +
            B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)));

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me2 =
      -16.0 * b00 * ivs3v2 * ivu1s1 * ivu2s1 * (LR + RL) *
      (LRLR * m1 * m2 * papb +
       (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
       m1 * m2 * papb * RLRL);

  return real(me0 + me1 + me2);
}

double FI::MVubu4t(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      -4.0 * ivt1s2 * ivu1s1 * ivu2s1 *
      ((4.0 * b00 + (b1 + b11) * (m2s - 2.0 * pap2)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)) *
      (LRLR * m1 * m2 * papb +
       LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
       m1 * m2 * papb * RLRL - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR +
       pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = 8.0 * b00 * ivt1s2 * ivu1s1 * ivu2s1 * (LR + RL) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return real(me0 + me1);
}

double FI::MVubu4u(double p1s, double ml1s, double ml2s, double mp1s,
                   double mp2s, int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 =
      8.0 * ivu1s1 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
      ((4.0 * b00 + (b1 + b11) * (m2s - 2.0 * pap2)) * (LR + RL) +
       B0(p1s, ml1s, ml2s) * ml1 * ml2 * (LL + RR)) *
      (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = -16.0 * b00 * ivu1s1 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

  return real(me0 + me1);
}

double FI::MVtbu5s(double ml1s, double mp1s, double mp2s, int ieps) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  complex<double> me0 =
      A0(ml1s) * ivs3v2 * ivt1s1 * ivt2s1 *
      (8.0 * RR * pap1 * pbp2 * RLRL + 8.0 * RR * pap1 * pbp2 * LRLR +
       4.0 * RR * m1 * m2 * papb * LRRR + 4.0 * RR * m1 * m2 * papb * RLLL);

  SetA(ml1s, ieps + 1, &a0, &a00);

  complex<double> me1 =
      A0(ml1s) * ivs3v2 * ivt1s1 * ivt2s1 *
      (4.0 * RR * pap2 * pbp1 * RLRL + 4.0 * RR * pap2 * pbp1 * LRLR -
       4.0 * RR * pap1 * pbp2 * RLRL - 4.0 * RR * pap1 * pbp2 * LRLR -
       4.0 * RR * papb * p1p2 * RLRL - 4.0 * RR * papb * p1p2 * LRLR -
       4.0 * RR * m1 * m2 * papb * LRRR - 4.0 * RR * m1 * m2 * papb * RLLL);

  return real(me0 + me1);
}

double FI::MVtbu5t(double ml1s, double mp1s, double mp2s, int ieps) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  //   complex<double> me0 = A0(ml1s) * ivt1s1 * ivt2s1 * ivt1s2 *
  //                          (4.0 * RR * pap1 * pbp2 * RLRL + 4.0 * RR * pap1 *
  //                           pbp2 * RRRR + 4.0 * RR * pap1 * pbp2 * LRLR + 4.0
  //                           * RR * pap1 * pbp2 * LLLL);

  complex<double> me0 = real(4.0 * A0(ml1s) * ivt1s1 * ivt2s1 * ivt1s2 * pap1 *
                             pbp2 * RR * (LLLL + LRLR + RLRL + RRRR));

  return real(me0);
}

double FI::MVtbu5u(double ml1s, double mp1s, double mp2s, int ieps) {
  ivt1s1 = 1.0 / (m1s - 2.0 * pap1 - mp1s);
  ivt2s1 = 1.0 / (m1s - 2.0 * pap1 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  complex<double> me0 =
      A0(ml1s) * ivt1s1 * ivt2s1 * ivu2s2 *
      (-2.0 * RR * pap2 * pbp1 * RRRR - 2.0 * RR * pap2 * pbp1 * LLLL -
       2.0 * RR * pap1 * pbp2 * RRRR - 2.0 * RR * pap1 * pbp2 * LLLL +
       2.0 * RR * papb * p1p2 * RRRR + 2.0 * RR * papb * p1p2 * LLLL -
       2.0 * RR * m1 * m2 * papb * RLRL - 2.0 * RR * m1 * m2 * papb * LRLR);

  return real(me0);
}

double FI::MVubu5s(double ml1s, double mp1s, double mp2s, int ieps) {
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  complex<double> me0 =
      A0(ml1s) * ivs3v2 * ivu1s1 * ivu2s1 *
      (-8.0 * RR * pap2 * pbp1 * LRRR - 8.0 * RR * pap2 * pbp1 * RLLL -
       4.0 * RR * m1 * m2 * papb * RLRL - 4.0 * RR * m1 * m2 * papb * LRLR);

  SetA(ml1s, ieps + 1, &a0, &a00);
  complex<double> me1 =
      A0(ml1s) * ivs3v2 * ivu1s1 * ivu2s1 *
      (4.0 * RR * pap2 * pbp1 * LRRR + 4.0 * RR * pap2 * pbp1 * RLLL -
       4.0 * RR * pap1 * pbp2 * LRRR - 4.0 * RR * pap1 * pbp2 * RLLL +
       4.0 * RR * papb * p1p2 * LRRR + 4.0 * RR * papb * p1p2 * RLLL +
       4.0 * RR * m1 * m2 * papb * RLRL + 4.0 * RR * m1 * m2 * papb * LRLR);

  return real(me0 + me1);
}

double FI::MVubu5t(double ml1s, double mp1s, double mp2s, int ieps) {
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  complex<double> me0 =
      A0(ml1s) * ivt1s2 * ivu1s1 * ivu2s1 *
      (-2.0 * RR * pap2 * pbp1 * RRRR - 2.0 * RR * pap2 * pbp1 * LLLL -
       2.0 * RR * pap1 * pbp2 * RRRR - 2.0 * RR * pap1 * pbp2 * LLLL +
       2.0 * RR * papb * p1p2 * RRRR + 2.0 * RR * papb * p1p2 * LLLL -
       2.0 * RR * m1 * m2 * papb * RLRL - 2.0 * RR * m1 * m2 * papb * LRLR);

  return real(me0);
}

double FI::MVubu5u(double ml1s, double mp1s, double mp2s, int ieps) {
  ivu1s1 = 1.0 / (m2s - 2.0 * pap2 - mp1s);
  ivu2s1 = 1.0 / (m2s - 2.0 * pap2 - mp2s);

  complex<double> a0;
  complex<double> a00;

  SetA(ml1s, ieps, &a0, &a00);

  complex<double> me0 =
      A0(ml1s) * ivu1s1 * ivu2s1 * ivu2s2 *
      (4.0 * RR * pap2 * pbp1 * RLRL + 4.0 * RR * pap2 * pbp1 * RRRR +
       4.0 * RR * pap2 * pbp1 * LRLR + 4.0 * RR * pap2 * pbp1 * LLLL);

  return real(me0);
}
