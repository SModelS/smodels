#include <complex>
#include <iostream>
//#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "utils.h"
using namespace std;

double FI::MVtct3s(double p1s, double ml1s, double ml2s, int ieps) {

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
       (B0(p1s, ml1s, ml2s) - 2.0 * B1(p1s, ml1s, ml2s) +
        B11(p1s, ml1s, ml2s)) *
           ml1s) *
      (m1 * m2 * papb * (LRRR + RLLL) + 2.0 * pap1 * pbp2 * (LRLR + RLRL)) * RR;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 =
      -4.0 * ivs3v2 * pow(ivt1s1, 2) *
      ((B0(p1s, ml1s, ml2s) - 2.0 * B1(p1s, ml1s, ml2s) +
        B11(p1s, ml1s, ml2s)) *
           ml1s *
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

  return -real(me0 + me1 + me2);
}

double FI::MVtct3t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = 4.0 * pow(ivt1s1, 2) * ivt1s2 *
                        (4.0 * b00 + (b0 - 2.0 * b1 + b11) * ml1s) * pap1 *
                        pbp2 * RR * (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = -8.0 * b00 * pow(ivt1s1, 2) * ivt1s2 * pap1 * pbp2 *
                        RR * (LLLL + LRLR + RLRL + RRRR);

  return -real(me0 + me1);
}

double FI::MVtct3u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = -2.0 * pow(ivt1s1, 2) * ivu2s2 *
                        (4.0 * b00 + (b0 - 2.0 * b1 + b11) * ml1s) * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = 4.0 * b00 * pow(ivt1s1, 2) * ivu2s2 * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return -real(me0 + me1);
}

double FI::MVuct3s(double p1s, double ml1s, double ml2s, int ieps) {

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
      (4.0 * b00 + (b0 - 2.0 * b1 + b11) * ml1s) *
      (2.0 * pap2 * pbp1 * (LRRR + RLLL) + m1 * m2 * papb * (LRLR + RLRL)) * RR;

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 =
      4.0 * ivs3v2 * pow(ivu2s1, 2) *
      ((b0 - 2.0 * b1 + b11) * ml1s *
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

  return -real(me0 + me1 + me2);
}

double FI::MVuct3t(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = -2.0 * ivt1s2 * pow(ivu2s1, 2) *
                        (4.0 * b00 + (b0 - 2.0 * b1 + b11) * ml1s) * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);

  complex<double> me1 = 4.0 * b00 * ivt1s2 * pow(ivu2s1, 2) * RR *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return -real(me0 + me1);
}

double FI::MVuct3u(double p1s, double ml1s, double ml2s, int ieps) {

  complex<double> b0;
  complex<double> b1;
  complex<double> b00;
  complex<double> b11;
  complex<double> db0;
  complex<double> db1;
  complex<double> db00;
  complex<double> db11;

  SetB(p1s, ml1s, ml2s, ieps, &b0, &b1, &b00, &b11, &db0, &db1, &db00, &db11);

  complex<double> me0 = 4.0 * pow(ivu2s1, 2) * ivu2s2 *
                        (4.0 * b00 + (b0 - 2.0 * b1 + b11) * ml1s) * pap2 *
                        pbp1 * RR * (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = -8.0 * b00 * pow(ivu2s1, 2) * ivu2s2 * pap2 * pbp1 *
                        RR * (LLLL + LRLR + RLRL + RRRR);

  return -real(me0 + me1);
}

double FI::MVtct4s(double p1s, double ml1s, double ml2s, double mp1s,
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
      (m1 * m2 * papb * (LRRR + RLLL) + 2 * pap1 * pbp2 * (LRLR + RLRL)) *
      ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
       b0 * ml1 * ml2 * (LL + RR));

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
           ((b1 + b11) * mp2s * (LR + RL) + b0 * ml1 * ml2 * (LL + RR)));

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me2 = 16.0 * b00 * ivs3v2 * ivt1s1 * ivt2s1 * (LR + RL) *
                        (LRRR * m1 * m2 * papb +
                         LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                         pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL);

  return -real(me0 + me1 + me2);
}

double FI::MVtct4t(double p1s, double ml1s, double ml2s, double mp1s,
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
  complex<double> me0 = 8.0 * ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
                        ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
                         b0 * ml1 * ml2 * (LL + RR)) *
                        (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = -16.0 * b00 * ivt1s1 * ivt2s1 * ivt1s2 * pap1 * pbp2 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

  return -real(me0 + me1);
}

double FI::MVtct4u(double p1s, double ml1s, double ml2s, double mp1s,
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
  complex<double> me0 = -4.0 * ivt1s1 * ivt2s1 * ivu2s2 *
                        ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
                         b0 * ml1 * ml2 * (LL + RR)) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = 8.0 * b00 * ivt1s1 * ivt2s1 * ivu2s2 * (LR + RL) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return -real(me0 + me1);
}

double FI::MVuct4s(double p1s, double ml1s, double ml2s, double mp1s,
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
      ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
       b0 * ml1 * ml2 * (LL + RR));

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
           ((b1 + b11) * mp2s * (LR + RL) + b0 * ml1 * ml2 * (LL + RR)));

  SetB(p1s, ml1s, ml2s, ieps + 2, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me2 =
      -16.0 * b00 * ivs3v2 * ivu1s1 * ivu2s1 * (LR + RL) *
      (LRLR * m1 * m2 * papb +
       (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * (LRRR + RLLL) +
       m1 * m2 * papb * RLRL);

  return -real(me0 + me1 + me2);
}

double FI::MVuct4t(double p1s, double ml1s, double ml2s, double mp1s,
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

  complex<double> me0 = -4.0 * ivt1s2 * ivu1s1 * ivu2s1 *
                        ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
                         b0 * ml1 * ml2 * (LL + RR)) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = 8.0 * b00 * ivt1s2 * ivu1s1 * ivu2s1 * (LR + RL) *
                        (LRLR * m1 * m2 * papb +
                         LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) +
                         m1 * m2 * papb * RLRL - p1p2 * papb * RRRR +
                         pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR);

  return -real(me0 + me1);
}

double FI::MVuct4u(double p1s, double ml1s, double ml2s, double mp1s,
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
  complex<double> me0 = 8.0 * ivu1s1 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                        ((4.0 * b00 + (b1 + b11) * mp2s) * (LR + RL) +
                         b0 * ml1 * ml2 * (LL + RR)) *
                        (LLLL + LRLR + RLRL + RRRR);

  SetB(p1s, ml1s, ml2s, ieps + 1, &b0, &b1, &b00, &b11, &db0, &db1, &db00,
       &db11);
  complex<double> me1 = -16.0 * b00 * ivu1s1 * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                        (LR + RL) * (LLLL + LRLR + RLRL + RRRR);

  return -real(me0 + me1);
}
