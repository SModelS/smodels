#include <complex>
#include <iostream>
//#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "utils.h"
using namespace std;

// Process: q + \bar{q} -> sl + \bar{sl}.

double FI::MVstr1sSL(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                     int ieps) {

  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = (4. * c00 + 4. * ((c1 + c2) + (c12 + c22)) * papb) *
                        (8. * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                         (2. * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 = (-8. * c00 - 4. * (c2 + (c12 + c22)) * papb) *
                        (8. * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                         (2. * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 = (4. * c00) * (8. * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                                      (2. * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  return real(me0 + me1 + me2);
}

double FI::MVstr2sSL(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                     int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //  Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 2.0 * c00 *
                        (8. * ivs3v1 * ivs3v2 * (LRLR * LR + LLLL * RL) *
                         (2. * pap1 * pap2 - pap1 * m2s - pap2 * m1s));

  return real(me0);
}

double FI::MVstr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;

  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      32. * ivs3v1 * ivs3v2 * (c00 + ((c1 + c2) + (c12 + c22)) * papb) *
      (LL * (m1 * m2 * papb * (LLRL + RLLL) + 2. * (LLLL * pap2 * pbp1 + pap1 * pbp2 * RLRL)) +
       RR * (m1 * m2 * papb * (LRRR + RRLR) + 2. * (LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -32. * ivs3v1 * ivs3v2 *
      (c00 * (LL * (3. * LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 7. * LLLL * pap2 * pbp1 -
                    3. * LLLL * pap1 * pbp2 + 3. * m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                    3. * pap2 * pbp1 * RLRL + 7. * pap1 * pbp2 * RLRL) +
              RR * (3. * LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3. * LRLR * pap2 * pbp1 +
                    7. * LRLR * pap1 * pbp2 + 3. * m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                    7. * pap2 * pbp1 * RRRR - 3. * pap1 * pbp2 * RRRR)) +
       papb *
           ((c1 + c2) *
                (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 3. * LLLL * pap2 * pbp1 -
                       3. * LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                       3. * pap2 * pbp1 * RLRL + 3. * pap1 * pbp2 * RLRL) +
                 RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3. * LRLR * pap2 * pbp1 +
                       3. * LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                       3. * pap2 * pbp1 * RRRR - 3. * pap1 * pbp2 * RRRR)) +
            (c12 + c22) *
                (LL * (2. * LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 5. * LLLL * pap2 * pbp1 -
                       3. * LLLL * pap1 * pbp2 + 2. * m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                       3. * pap2 * pbp1 * RLRL + 5. * pap1 * pbp2 * RLRL) +
                 RR * (2. * LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3. * LRLR * pap2 * pbp1 +
                       5. * LRLR * pap1 * pbp2 + 2. * m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                       5. * pap2 * pbp1 * RRRR - 3. * pap1 * pbp2 * RRRR)) +
            c2 * (LL * (m1 * m2 * papb * (LLRL + RLLL) +
                        2. * (LLLL * pap2 * pbp1 + pap1 * pbp2 * RLRL)) +
                  RR * (m1 * m2 * papb * (LRRR + RRLR) +
                        2. * (LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      32. * ivs3v1 * ivs3v2 *
      (c00 * (LL * (3. * LLRL * m1 * m2 * papb +
                    2. * LLLL * (p1p2 * papb + 5. * pap2 * pbp1 - 4. * pap1 * pbp2) +
                    3. * m1 * m2 * papb * RLLL +
                    2. * (p1p2 * papb - 4. * pap2 * pbp1 + 5. * pap1 * pbp2) * RLRL) +
              RR * (3. * LRRR * m1 * m2 * papb +
                    2. * LRLR * (p1p2 * papb - 4. * pap2 * pbp1 + 5. * pap1 * pbp2) +
                    3. * m1 * m2 * papb * RRLR +
                    2. * (p1p2 * papb + 5. * pap2 * pbp1 - 4. * pap1 * pbp2) * RRRR)) +
       papb * (2. * (c1 + c2) * (pap2 * pbp1 - pap1 * pbp2) *
                   (LL * (LLLL - RLRL) + RR * (-LRLR + RRRR)) +
               (c12 + c22) *
                   (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 5. * LLLL * pap2 * pbp1 -
                          5. * LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                          5. * pap2 * pbp1 * RLRL + 5. * pap1 * pbp2 * RLRL) +
                    RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 5. * LRLR * pap2 * pbp1 +
                          5. * LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                          5. * pap2 * pbp1 * RRRR - 5. * pap1 * pbp2 * RRRR)) +
               c2 * (LL * (LLRL * m1 * m2 * papb + LLLL * p1p2 * papb + 3. * LLLL * pap2 * pbp1 -
                           3. * LLLL * pap1 * pbp2 + m1 * m2 * papb * RLLL + p1p2 * papb * RLRL -
                           3. * pap2 * pbp1 * RLRL + 3. * pap1 * pbp2 * RLRL) +
                     RR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3. * LRLR * pap2 * pbp1 +
                           3. * LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
                           3. * pap2 * pbp1 * RRRR - 3. * pap1 * pbp2 * RRRR))));

  // Returns the finite result of the whole squared matrix element.
  return real(me0 + me1 + me2);
}

double FI::MVstr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {

  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 16. * ivs3v1 * ivt1s2 * (c00 + ((c1 + c2) + (c12 + c22)) * papb) *
                        (LL * LLRL * m1 * m2 * papb + 2. * LL * pap1 * pbp2 * RLRL +
                         2. * LRLR * pap1 * pbp2 * RR + m1 * m2 * papb * RR * RRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -16. * ivs3v1 * ivt1s2 *
      (c00 * (LL * (3. * LLRL * m1 * m2 * papb +
                    (p1p2 * papb - pap2 * pbp1 + 5. * pap1 * pbp2) * RLRL) +
              RR * (LRLR * (p1p2 * papb - pap2 * pbp1 + 5. * pap1 * pbp2) +
                    3. * m1 * m2 * papb * RRLR)) +
       papb * (c2 * (LL * LLRL * m1 * m2 * papb + 2. * LL * pap1 * pbp2 * RLRL +
                     2. * LRLR * pap1 * pbp2 * RR + m1 * m2 * papb * RR * RRLR) +
               (c1 + c2) * (LL * (LLRL * m1 * m2 * papb + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                                  pap1 * pbp2 * RLRL) +
                            RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                                  m1 * m2 * papb * RRLR)) +
               (c12 + c22) * (LL * (2. * LLRL * m1 * m2 * papb + p1p2 * papb * RLRL -
                                    pap2 * pbp1 * RLRL + 3. * pap1 * pbp2 * RLRL) +
                              RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 +
                                    3. * LRLR * pap1 * pbp2 + 2. * m1 * m2 * papb * RRLR))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 = 16. * ivs3v1 * ivt1s2 *
                        ((c2 + (c12 + c22)) * papb *
                             (LL * (LLRL * m1 * m2 * papb + p1p2 * papb * RLRL -
                                    pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL) +
                              RR * (LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2 +
                                    m1 * m2 * papb * RRLR)) +
                         c00 * (LL * (3. * LLRL * m1 * m2 * papb + 2. * p1p2 * papb * RLRL -
                                      2. * pap2 * pbp1 * RLRL + 4. * pap1 * pbp2 * RLRL) +
                                RR * (2. * LRLR * p1p2 * papb - 2. * LRLR * pap2 * pbp1 +
                                      4. * LRLR * pap1 * pbp2 + 3. * m1 * m2 * papb * RRLR)));

  return real(me0 + me1 + me2);
}

// sth wrong here? MARCEL
double FI::MVstr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = -16. * ivs3v1 * ivu2s2 * (c00 + ((c1 + c2) + (c12 + c22)) * papb) *
                        (2. * LL * LLRL * pap2 * pbp1 + LL * m1 * m2 * papb * RLRL +
                         LRLR * m1 * m2 * papb * RR + 2. * pap2 * pbp1 * RR * RRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      16. * ivs3v1 * ivu2s2 *
      (c00 * (LL * (LLRL * (p1p2 * papb + 5. * pap2 * pbp1 - pap1 * pbp2) +
                    3. * m1 * m2 * papb * RLRL) +
              RR * (3. * LRLR * m1 * m2 * papb +
                    (p1p2 * papb + 5. * pap2 * pbp1 - pap1 * pbp2) * RRLR)) +
       papb * (c2 * (2. * LL * LLRL * pap2 * pbp1 + LL * m1 * m2 * papb * RLRL +
                     LRLR * m1 * m2 * papb * RR + 2. * pap2 * pbp1 * RR * RRLR) +
               (c12 + c22) * (LL * LLRL * p1p2 * papb + 3. * LL * LLRL * pap2 * pbp1 -
                              LL * LLRL * pap1 * pbp2 + 2. * LL * m1 * m2 * papb * RLRL +
                              2. * LRLR * m1 * m2 * papb * RR +
                              (p1p2 * papb + 3. * pap2 * pbp1 - pap1 * pbp2) * RR * RRLR) +
               (c1 + c2) * (LL * (LLRL * p1p2 * papb + LLRL * pap2 * pbp1 - LLRL * pap1 * pbp2 +
                                  m1 * m2 * papb * RLRL) +
                            RR * (LRLR * m1 * m2 * papb + p1p2 * papb * RRLR + pap2 * pbp1 * RRLR -
                                  pap1 * pbp2 * RRLR))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 = -16. * ivs3v1 * ivu2s2 *
                        (c00 * (LL * (2. * LLRL * p1p2 * papb + 4. * LLRL * pap2 * pbp1 -
                                      2. * LLRL * pap1 * pbp2 + 3. * m1 * m2 * papb * RLRL) +
                                RR * (3. * LRLR * m1 * m2 * papb + 2. * p1p2 * papb * RRLR +
                                      4. * pap2 * pbp1 * RRLR - 2. * pap1 * pbp2 * RRLR)) +
                         (c2 + (c12 + c22)) * papb *
                             (LL * (LLRL * p1p2 * papb + LLRL * pap2 * pbp1 - LLRL * pap1 * pbp2 +
                                    m1 * m2 * papb * RLRL) +
                              RR * (LRLR * m1 * m2 * papb + p1p2 * papb * RRLR +
                                    pap2 * pbp1 * RRLR - pap1 * pbp2 * RRLR)));

  return real(me0 + me1 + me2);
}

double FI::MVstr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //    Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      16. * c00 * ivs3v1 * ivs3v2 *
      (RL * (m1 * m2 * papb * (LRRL + RRLL) + 2. * (LRLL * pap2 * pbp1 + pap1 * pbp2 * RRRL)) +
       LR * (m1 * m2 * papb * (LRRR + RRLR) + 2. * (LRLR * pap1 * pbp2 + pap2 * pbp1 * RRRR)));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -16. * c00 * ivs3v1 * ivs3v2 *
      (RL * (LRRL * m1 * m2 * papb + LRLL * p1p2 * papb + 3. * LRLL * pap2 * pbp1 -
             3. * LRLL * pap1 * pbp2 + m1 * m2 * papb * RRLL + p1p2 * papb * RRRL -
             3. * pap2 * pbp1 * RRRL + 3. * pap1 * pbp2 * RRRL) +
       LR * (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - 3. * LRLR * pap2 * pbp1 +
             3. * LRLR * pap1 * pbp2 + m1 * m2 * papb * RRLR + p1p2 * papb * RRRR +
             3. * pap2 * pbp1 * RRRR - 3. * pap1 * pbp2 * RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 = 32. * c00 * ivs3v1 * ivs3v2 * (pap2 * pbp1 - pap1 * pbp2) *
                        (RL * (LRLL - RRRL) + LR * (-LRLR + RRRR));

  return real(me0 + me1 + me2);
}

double FI::MVstr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);

  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //  Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivs3v1 * ivt1s2 *
                        (2. * c00 *
                             (2. * LR * LRLR * pap1 * pbp2 + LRRL * m1 * m2 * papb * RL +
                              LR * m1 * m2 * papb * RRLR + 2. * pap1 * pbp2 * RL * RRRL) +
                         ml1 * papb *
                             (c0 * LL * LRLL * m1 * pbp2 + 2. * (c1 + c2) * LL * LRLL * m1 * pbp2 -
                              2. * c2 * LL * LRLL * m1 * pbp2 + c0 * LRRR * m2 * pap1 * RR +
                              2. * c2 * LRRR * m2 * pap1 * RR + c0 * LL * m2 * pap1 * RRLL +
                              2. * c2 * LL * m2 * pap1 * RRLL +
                              (c0 + 2. * (c1 + c2) - 2. * c2) * m1 * pbp2 * RR * RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -8. * c00 * ivs3v1 * ivt1s2 *
      (LR * (LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RRLR) +
       RL * (LRRL * m1 * m2 * papb + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRL));

  return real(me0 + me1);
}

double FI::MVstr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = -4. * ivs3v1 * ivu2s2 *
                        (2. * c00 *
                             (LR * LRLR * m1 * m2 * papb + 2. * LRRL * pap2 * pbp1 * RL +
                              2. * LR * pap2 * pbp1 * RRLR + m1 * m2 * papb * RL * RRRL) +
                         ml1 * papb *
                             (c0 * (LL * LRLL * m1 * pap2 + LRRR * m2 * pbp1 * RR +
                                    LL * m2 * pbp1 * RRLL + m1 * pap2 * RR * RRRR) +
                              2. * ((c1 + c2) * m2 * pbp1 * (LRRR * RR + LL * RRLL) +
                                    c2 * (LL * LRLL * m1 * pap2 - LRRR * m2 * pbp1 * RR -
                                          LL * m2 * pbp1 * RRLL + m1 * pap2 * RR * RRRR))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      8. * c00 * ivs3v1 * ivu2s2 *
      (LR * (LRLR * m1 * m2 * papb + (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) * RRLR) +
       RL * (LRRL * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2) + m1 * m2 * papb * RRRL));

  return real(me0 + me1);
}

double FI::MVstr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      c0 * ml1 * ivs3v1 * ivt1s2 *
      (-2. * pap2 * pbp1 * RR * RRRR - 2. * pap2 * pbp1 * LL * LRLL + 2. * pap1 * pbp2 * RR * RRRR +
       2. * pap1 * pbp2 * LL * LRLL + 2. * papb * p1p2 * RR * RRRR + 2. * papb * p1p2 * LL * LRLL -
       2. * m1 * m2 * papb * RR * LRRR - 2. * m1 * m2 * papb * LL * RRLL);

  return real(me0);
}

double FI::MVstr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., 0., 2.*papb, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      c0 * ml1 * ivs3v1 * ivu2s2 *
      (2. * pap2 * pbp1 * RR * RRRR + 2. * pap2 * pbp1 * LL * LRLL - 2. * pap1 * pbp2 * RR * RRRR -
       2. * pap1 * pbp2 * LL * LRLL + 2. * papb * p1p2 * RR * RRRR + 2. * papb * p1p2 * LL * LRLL -
       2. * m1 * m2 * papb * RR * LRRR - 2. * m1 * m2 * papb * LL * RRLL);

  return real(me0);
}

double FI::MVttr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivs3v2 * ivt1s1 *
                        (4. * c00 + (2. * c2 + c22) * m1s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap1) *
                        (LR * m1 * m2 * papb * RLLL + 2. * LR * pap1 * pbp2 * RLRL +
                         LRRR * m1 * m2 * papb * RR + 2. * LRLR * pap1 * pbp2 * RR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -4. * ivs3v2 * ivt1s1 *
      (((2. * c2 + c22) * m1s + 2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap1) *
           (LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                  pap1 * pbp2 * RLRL) +
            (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2) *
                RR) +
       c00 * (LR * (6. * m1 * m2 * papb * RLLL + 4. * p1p2 * papb * RLRL - 4. * pap2 * pbp1 * RLRL +
                    8. * pap1 * pbp2 * RLRL) +
              2. *
                  (3. * LRRR * m1 * m2 * papb + 2. * LRLR * p1p2 * papb - 2. * LRLR * pap2 * pbp1 +
                   4. * LRLR * pap1 * pbp2) *
                  RR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      8. * c00 * ivs3v2 * ivt1s1 *
      (LR * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) +
       (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) * RR);

  return real(me0 + me1 + me2);
}

double FI::MVttr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivt1s1 * ivt1s2 * pap1 *
                        (4. * c00 + (2. * c2 + c22) * m1s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap1) *
                        pbp2 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -8. * c00 * ivt1s1 * ivt1s2 * pap1 * pbp2 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

  return real(me0 + me1);
}

double FI::MVttr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      2. * ivt1s1 * ivu2s2 *
      (4. * c00 + (2. * c2 + c22) * m1s +
       2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap1) *
      ((RRRR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) * RR +
       LR * (LRLR * m1 * m2 * papb + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * LLLL));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -4. * c00 * ivt1s1 * ivu2s2 *
      ((RRRR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) * RR +
       LR * (LRLR * m1 * m2 * papb + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * LLLL));

  return -real(me0 + me1);
}

double FI::MVttr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivs3v2 * ivt1s1 *
                        (4. * c00 + (2. * c2 + c22) * m2s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp2) *
                        (LR * m1 * m2 * papb * RLLL + 2. * LR * pap1 * pbp2 * RLRL +
                         LRRR * m1 * m2 * papb * RR + 2. * LRLR * pap1 * pbp2 * RR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -4. * ivs3v2 * ivt1s1 *
      (((2. * c2 + c22) * m2s + 2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp2) *
           (LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL +
                  pap1 * pbp2 * RLRL) +
            (LRRR * m1 * m2 * papb + LRLR * p1p2 * papb - LRLR * pap2 * pbp1 + LRLR * pap1 * pbp2) *
                RR) +
       c00 * (LR * (6. * m1 * m2 * papb * RLLL + 4. * p1p2 * papb * RLRL - 4. * pap2 * pbp1 * RLRL +
                    8. * pap1 * pbp2 * RLRL) +
              2. *
                  (3. * LRRR * m1 * m2 * papb + 2. * LRLR * p1p2 * papb - 2. * LRLR * pap2 * pbp1 +
                   4. * LRLR * pap1 * pbp2) *
                  RR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);
  complex<double> me2 =
      8. * c00 * ivs3v2 * ivt1s1 *
      (LR * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) +
       (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) * RR);

  return real(me0 + me1 + me2);
}

double FI::MVttr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivt1s1 * ivt1s2 * pap1 * pbp2 *
                        (4. * c00 + (2. * c2 + c22) * m2s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp2) *
                        ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -8. * c00 * ivt1s1 * ivt1s2 * pap1 * pbp2 * ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

  return real(me0 + me1);
}

double FI::MVttr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -2. * ivt1s1 * ivu2s2 *
      (4. * c00 + (2. * c2 + c22) * m2s +
       2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp2) *
      ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
       LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * c00 * ivt1s1 * ivu2s2 *
      ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
       LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

  return real(me0 + me1);
}

double FI::MVutr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = -4. * ivs3v2 * ivu2s1 *
                        (4. * c00 + (2. * c2 + c22) * m2s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap2) *
                        (2. * LR * pap2 * pbp1 * RLLL + LR * m1 * m2 * papb * RLRL +
                         LRLR * m1 * m2 * papb * RR + 2. * LRRR * pap2 * pbp1 * RR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * ivs3v2 * ivu2s1 *
      (c00 * (LR * (4. * p1p2 * papb * RLLL + 8. * pap2 * pbp1 * RLLL - 4. * pap1 * pbp2 * RLLL +
                    6. * m1 * m2 * papb * RLRL) +
              2. *
                  (3. * LRLR * m1 * m2 * papb + 2. * LRRR * p1p2 * papb + 4. * LRRR * pap2 * pbp1 -
                   2. * LRRR * pap1 * pbp2) *
                  RR) +
       ((2. * c2 + c22) * m2s + 2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap2) *
           (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL +
                  m1 * m2 * papb * RLRL) +
            (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1 - LRRR * pap1 * pbp2) *
                RR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      -8. * c00 * ivs3v2 * ivu2s1 *
      (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
       (LRLR * m1 * m2 * papb + LRRR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR);

  return real(me0 + me1 + me2);
}

double FI::MVutr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -2. * ivt1s2 * ivu2s1 *
      (4. * c00 + (2. * c2 + c22) * m2s +
       2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap2) *
      (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + LR * m1 * m2 * papb * RLRL +
       RR * (LRLR * m1 * m2 * papb - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * c00 * ivt1s2 * ivu2s1 *
      (LLLL * LR * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + LR * m1 * m2 * papb * RLRL +
       RR * (LRLR * m1 * m2 * papb - p1p2 * papb * RRRR + pap2 * pbp1 * RRRR + pap1 * pbp2 * RRRR));

  return real(me0 + me1);
}

double FI::MVutr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivu2s1 * ivu2s2 * pap2 *
                        (4. * c00 + (2. * c2 + c22) * m2s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pap2) *
                        pbp1 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -8. * c00 * ivu2s1 * ivu2s2 * pap2 * pbp1 * (LR * (LLLL + RLRL) + RR * (LRLR + RRRR));

  return real(me0 + me1);
}

double FI::MVutr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = -4. * ivs3v2 * ivu2s1 *
                        (4. * c00 + (2. * c2 + c22) * m1s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp1) *
                        (2. * LR * pap2 * pbp1 * RLLL + LR * m1 * m2 * papb * RLRL +
                         LRLR * m1 * m2 * papb * RR + 2. * LRRR * pap2 * pbp1 * RR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * ivs3v2 * ivu2s1 *
      (c00 * (LR * (4. * p1p2 * papb * RLLL + 8. * pap2 * pbp1 * RLLL - 4. * pap1 * pbp2 * RLLL +
                    6. * m1 * m2 * papb * RLRL) +
              2. *
                  (3. * LRLR * m1 * m2 * papb + 2. * LRRR * p1p2 * papb + 4. * LRRR * pap2 * pbp1 -
                   2. * LRRR * pap1 * pbp2) *
                  RR) +
       ((2. * c2 + c22) * m1s + 2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp1) *
           (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL +
                  m1 * m2 * papb * RLRL) +
            (LRLR * m1 * m2 * papb + LRRR * p1p2 * papb + LRRR * pap2 * pbp1 - LRRR * pap1 * pbp2) *
                RR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      -8. * c00 * ivs3v2 * ivu2s1 *
      (LR * (p1p2 * papb * RLLL + pap2 * pbp1 * RLLL - pap1 * pbp2 * RLLL + m1 * m2 * papb * RLRL) +
       (LRLR * m1 * m2 * papb + LRRR * (p1p2 * papb + pap2 * pbp1 - pap1 * pbp2)) * RR);

  return real(me0 + me1 + me2);
}

double FI::MVutr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);
  complex<double> la[7];

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -2. * ivt1s2 * ivu2s1 *
      (4. * c00 + (2. * c2 + c22) * m1s +
       2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp1) *
      ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
       LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * c00 * ivt1s2 * ivu2s1 *
      ((LRLR * m1 * m2 * papb + LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2)) * RR +
       LR * (m1 * m2 * papb * RLRL + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

  return real(me0 + me1);
}

double FI::MVutr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 = 4. * ivu2s1 * ivu2s2 * pap2 * pbp1 *
                        (4. * c00 + (2. * c2 + c22) * m1s +
                         2. * (2. * (c1 + c2) - 2. * c2 - c22 + (c12 + c22)) * pbp1) *
                        ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -8. * c00 * ivu2s1 * ivu2s2 * pap2 * pbp1 * ((LLLL + LRLR) * RR + LR * (RLRL + RRRR));

  return real(me0 + me1);
}

double FI::MVttr3s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);

  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //   Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -4. * ivs3v2 * ivt1s1 *
      (c22 * LR * LRRR * m1 * m1s * m2 * papb -
       2. * (c12 + c22) * LR * LRRR * m1 * m2 * pap1 * papb +
       2. * c22 * LR * LRLR * m1s * pap1 * pbp2 -
       4. * (c12 + c22) * LR * LRLR * pow(pap1, 2) * pbp2 +
       c0 * LL * m1 * m2 * ml1 * ml3 * papb * RLLL + c22 * m1 * m1s * m2 * papb * RL * RLLL -
       2. * (c12 + c22) * m1 * m2 * pap1 * papb * RL * RLLL +
       2. * c0 * LL * ml1 * ml3 * pap1 * pbp2 * RLRL + 2. * c22 * m1s * pap1 * pbp2 * RL * RLRL -
       4. * (c12 + c22) * pow(pap1, 2) * pbp2 * RL * RLRL +
       4. * c00 *
           (LR * LRRR * m1 * m2 * papb + 2. * LR * LRLR * pap1 * pbp2 + m1 * m2 * papb * RL * RLLL +
            2. * pap1 * pbp2 * RL * RLRL) +
       c0 * LLRR * m1s * m2 * ml1 * papb * RR + c0 * LRRR * m1 * m2 * ml1 * ml3 * papb * RR +
       2. * c0 * LLLR * m1 * ml1 * pap1 * pbp2 * RR +
       2. * c0 * LRLR * ml1 * ml3 * pap1 * pbp2 * RR + c0 * LL * m1s * m2 * ml1 * papb * RRLL +
       2. * c0 * LL * m1 * ml1 * pap1 * pbp2 * RRRL +
       c2 * (LR * m2 * (LLRR * m1s * ml3 + LRRR * m1 * (m1s - 2. * pap1)) * papb +
             2. * LR * (LLLR * m1 * ml3 + LRLR * (m1s - 2. * pap1)) * pap1 * pbp2 +
             2. * m1s * pap1 * pbp2 * RL * RLRL - 4. * pow(pap1, 2) * pbp2 * RL * RLRL +
             LLRR * m1s * m2 * ml1 * papb * RR + m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRLL +
             m1 * (m1s * m2 * papb * RL * RLLL +
                   2. * pap1 *
                       (-(m2 * papb * RL * RLLL) + LLLR * ml1 * pbp2 * RR + LL * ml1 * pbp2 * RRRL +
                        ml3 * pbp2 * RL * RRRL))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * ivs3v2 * ivt1s1 *
      (c22 * LR * LRRR * m1 * m1s * m2 * papb + c22 * LR * LRLR * m1s * p1p2 * papb -
       2. * (c12 + c22) * LR * LRRR * m1 * m2 * pap1 * papb -
       2. * (c12 + c22) * LR * LRLR * p1p2 * pap1 * papb - c22 * LR * LRLR * m1s * pap2 * pbp1 +
       2. * (c12 + c22) * LR * LRLR * pap1 * pap2 * pbp1 + c22 * LR * LRLR * m1s * pap1 * pbp2 -
       2. * (c12 + c22) * LR * LRLR * pow(pap1, 2) * pbp2 +
       c0 * LL * m1 * m2 * ml1 * ml3 * papb * RLLL + c22 * m1 * m1s * m2 * papb * RL * RLLL -
       2. * (c12 + c22) * m1 * m2 * pap1 * papb * RL * RLLL +
       c0 * LL * ml1 * ml3 * p1p2 * papb * RLRL - c0 * LL * ml1 * ml3 * pap2 * pbp1 * RLRL +
       c0 * LL * ml1 * ml3 * pap1 * pbp2 * RLRL + c22 * m1s * p1p2 * papb * RL * RLRL -
       2. * (c12 + c22) * p1p2 * pap1 * papb * RL * RLRL - c22 * m1s * pap2 * pbp1 * RL * RLRL +
       2. * (c12 + c22) * pap1 * pap2 * pbp1 * RL * RLRL + c22 * m1s * pap1 * pbp2 * RL * RLRL -
       2. * (c12 + c22) * pow(pap1, 2) * pbp2 * RL * RLRL +
       c00 * (LR * (6. * LRRR * m1 * m2 * papb +
                    4. * LRLR * (p1p2 * papb - pap2 * pbp1 + 2. * pap1 * pbp2)) +
              6. * m1 * m2 * papb * RL * RLLL +
              4. * (p1p2 * papb - pap2 * pbp1 + 2. * pap1 * pbp2) * RL * RLRL) +
       c0 * LLRR * m1s * m2 * ml1 * papb * RR + c0 * LRRR * m1 * m2 * ml1 * ml3 * papb * RR +
       c0 * LLLR * m1 * ml1 * p1p2 * papb * RR + c0 * LRLR * ml1 * ml3 * p1p2 * papb * RR -
       c0 * LLLR * m1 * ml1 * pap2 * pbp1 * RR - c0 * LRLR * ml1 * ml3 * pap2 * pbp1 * RR +
       c0 * LLLR * m1 * ml1 * pap1 * pbp2 * RR + c0 * LRLR * ml1 * ml3 * pap1 * pbp2 * RR +
       c0 * LL * m1s * m2 * ml1 * papb * RRLL +
       c0 * LL * m1 * ml1 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RRRL +
       c2 * (LR * (LLRR * m1s * m2 * ml3 * papb + LRRR * m1 * m2 * (m1s - 2. * pap1) * papb +
                   (LRLR * m1s + LLLR * m1 * ml3 - 2. * LRLR * pap1) *
                       (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) +
             m1s * p1p2 * papb * RL * RLRL - 2. * p1p2 * pap1 * papb * RL * RLRL -
             m1s * pap2 * pbp1 * RL * RLRL + 2. * pap1 * pap2 * pbp1 * RL * RLRL +
             m1s * pap1 * pbp2 * RL * RLRL - 2. * pow(pap1, 2) * pbp2 * RL * RLRL +
             LLRR * m1s * m2 * ml1 * papb * RR + m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRLL +
             m1 * (m1s * m2 * papb * RL * RLLL - 2. * m2 * pap1 * papb * RL * RLLL +
                   (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) *
                       (LLLR * ml1 * RR + LL * ml1 * RRRL + ml3 * RL * RRRL))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      -8. * c00 * ivs3v2 * ivt1s1 *
      (LR * (LRRR * m1 * m2 * papb + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2)) +
       RL * (m1 * m2 * papb * RLLL + (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL));

  return real(me0 + me1 + me2);
}

double FI::MVttr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -4. * ivt1s1 * ivt1s2 * pap1 * pbp2 *
      (c22 * LR * LRLR * m1s - 2. * (c12 + c22) * LR * LRLR * pap1 + c22 * LLLL * m1s * RL -
       2. * (c12 + c22) * LLLL * pap1 * RL + c22 * m1s * RL * RLRL -
       2. * (c12 + c22) * pap1 * RL * RLRL + c0 * LRLR * ml1 * ml3 * RR +
       c0 * m1 * ml1 * RLRR * RR +
       c0 * ml1 * (LL * (LRLL * m1 + ml3 * (LLLL + RLRL)) + LLLR * m1 * RR) +
       c0 * LL * m1 * ml1 * RRRL + c22 * LR * m1s * RRRR - 2. * (c12 + c22) * LR * pap1 * RRRR +
       c0 * ml1 * ml3 * RR * RRRR + 4. * c00 * (RL * (LLLL + RLRL) + LR * (LRLR + RRRR)) +
       c2 * (LLLL * m1s * RL + LRLL * m1 * ml3 * RL - 2. * LLLL * pap1 * RL + m1s * RL * RLRL -
             2. * pap1 * RL * RLRL + LLLR * m1 * ml1 * RR + m1 * ml1 * RLRR * RR +
             m1 * ml3 * RL * RRRL + LL * m1 * ml1 * (LRLL + RRRL) +
             LR * (m1 * ml3 * (LLLR + RLRR) + (m1s - 2. * pap1) * (LRLR + RRRR))));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      8. * c00 * ivt1s1 * ivt1s2 * pap1 * pbp2 * (RL * (LLLL + RLRL) + LR * (LRLR + RRRR));

  return real(me0 + me1);
}

double FI::MVttr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //   Npf *cc = new Npf(0., m1s - 2.*pap1, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      2. * ivt1s1 * ivu2s2 *
      (4. * c00 * LR * LRLR * m1 * m2 * papb + c22 * LR * LRLR * m1 * m1s * m2 * papb -
       c0 * LL * LRLL * m1 * ml1 * p1p2 * papb - c2 * LL * LRLL * m1 * ml1 * p1p2 * papb -
       c0 * LL * LLLL * ml1 * ml3 * p1p2 * papb +
       c2 * LR * LRLR * m1 * m2 * (m1s - 2. * pap1) * papb -
       2. * (c12 + c22) * LR * LRLR * m1 * m2 * pap1 * papb +
       c0 * LL * LRLL * m1 * ml1 * pap2 * pbp1 + c2 * LL * LRLL * m1 * ml1 * pap2 * pbp1 +
       c0 * LL * LLLL * ml1 * ml3 * pap2 * pbp1 + c0 * LL * LRLL * m1 * ml1 * pap1 * pbp2 +
       c2 * LL * LRLL * m1 * ml1 * pap1 * pbp2 + c0 * LL * LLLL * ml1 * ml3 * pap1 * pbp2 -
       c2 * LLLL * m1s * p1p2 * papb * RL - c22 * LLLL * m1s * p1p2 * papb * RL -
       c2 * LRLL * m1 * ml3 * p1p2 * papb * RL + 2. * c2 * LLLL * p1p2 * pap1 * papb * RL +
       2. * (c12 + c22) * LLLL * p1p2 * pap1 * papb * RL + c2 * LLLL * m1s * pap2 * pbp1 * RL +
       c22 * LLLL * m1s * pap2 * pbp1 * RL + c2 * LRLL * m1 * ml3 * pap2 * pbp1 * RL -
       2. * c2 * LLLL * pap1 * pap2 * pbp1 * RL -
       2. * (c12 + c22) * LLLL * pap1 * pap2 * pbp1 * RL + c2 * LLLL * m1s * pap1 * pbp2 * RL +
       c22 * LLLL * m1s * pap1 * pbp2 * RL + c2 * LRLL * m1 * ml3 * pap1 * pbp2 * RL -
       2. * c2 * LLLL * pow(pap1, 2) * pbp2 * RL -
       2. * (c12 + c22) * LLLL * pow(pap1, 2) * pbp2 * RL +
       c0 * LL * m1 * m2 * ml1 * ml3 * papb * RLRL + c2 * m1 * m1s * m2 * papb * RL * RLRL +
       c22 * m1 * m1s * m2 * papb * RL * RLRL - 2. * c2 * m1 * m2 * pap1 * papb * RL * RLRL -
       2. * (c12 + c22) * m1 * m2 * pap1 * papb * RL * RLRL +
       4. * c00 * RL *
           (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) +
       c2 * LR * ml3 *
           (LLLR * m1s * m2 * papb + m1 * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RLRR) +
       c0 * LLLR * m1s * m2 * ml1 * papb * RR + c2 * LLLR * m1s * m2 * ml1 * papb * RR +
       c0 * LRLR * m1 * m2 * ml1 * ml3 * papb * RR - c0 * m1 * ml1 * p1p2 * papb * RLRR * RR -
       c2 * m1 * ml1 * p1p2 * papb * RLRR * RR + c0 * m1 * ml1 * pap2 * pbp1 * RLRR * RR +
       c2 * m1 * ml1 * pap2 * pbp1 * RLRR * RR + c0 * m1 * ml1 * pap1 * pbp2 * RLRR * RR +
       c2 * m1 * ml1 * pap1 * pbp2 * RLRR * RR + c0 * LL * m1s * m2 * ml1 * papb * RRRL +
       c2 * m1s * m2 * papb * (LL * ml1 + ml3 * RL) * RRRL +
       (p1p2 * papb - pap2 * pbp1 - pap1 * pbp2) *
           (-4. * c00 * LR - (c2 + c22) * LR * m1s + 2. * (c2 + (c12 + c22)) * LR * pap1 -
            c0 * ml1 * ml3 * RR) *
           RRRR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      -4. * c00 * ivt1s1 * ivu2s2 *
      (RL * (LLLL * (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) + m1 * m2 * papb * RLRL) +
       LR * (LRLR * m1 * m2 * papb + (-(p1p2 * papb) + pap2 * pbp1 + pap1 * pbp2) * RRRR));

  return real(me0 + me1);
}

double FI::MVttr4s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      -4. * ivs3v2 * ivt1s1 *
      ((4. * c00 + c22 * m2s - 2. * (c12 + c22) * pbp2) *
           (LRRR * m1 * m2 * papb * RL + LR * m1 * m2 * papb * RLLL +
            2. * pap1 * pbp2 * (LRLR * RL + LR * RLRL)) +
       c2 * (LRRR * m1 * m2 * m2s * papb * RL + 2. * LRLR * m2s * pap1 * pbp2 * RL -
             2. * LRRR * m1 * m2 * papb * pbp2 * RL - 4. * LRLR * pap1 * pow(pbp2, 2) * RL +
             LR * m1 * m2 * m2s * papb * RLLL - 2. * LR * m1 * m2 * papb * pbp2 * RLLL +
             2. * LR * m2s * pap1 * pbp2 * RLRL - 4. * LR * pap1 * pow(pbp2, 2) * RLRL +
             LLLL * m1 * m2s * papb * (LR * ml3 + ml1 * RR) +
             2. * LLRL * m2 * pap1 * pbp2 * (LR * ml3 + ml1 * RR) +
             2. * LL * m2 * ml1 * pap1 * pbp2 * RRLR + 2. * m2 * ml3 * pap1 * pbp2 * RL * RRLR +
             m1 * m2s * papb * (LL * ml1 + ml3 * RL) * RRRR) +
       c0 * ml1 *
           ((LLLL * m1 * m2s * papb + 2. * LLRL * m2 * pap1 * pbp2 + m1 * m2 * ml3 * papb * RLLL +
             2. * ml3 * pap1 * pbp2 * RLRL) *
                RR +
            LL * (LRRR * m1 * m2 * ml3 * papb + 2. * pap1 * pbp2 * (LRLR * ml3 + m2 * RRLR) +
                  m1 * m2s * papb * RRRR)));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      4. * ivs3v2 * ivt1s1 *
      (.6 * c00 * LRRR * m1 * m2 * papb * RL + c22 * LRRR * m1 * m2 * m2s * papb * RL +
       4. * c00 * LRLR * p1p2 * papb * RL + c22 * LRLR * m2s * p1p2 * papb * RL -
       4. * c00 * LRLR * pap2 * pbp1 * RL - c22 * LRLR * m2s * pap2 * pbp1 * RL +
       8. * c00 * LRLR * pap1 * pbp2 * RL + c22 * LRLR * m2s * pap1 * pbp2 * RL -
       2. * (c12 + c22) * LRRR * m1 * m2 * papb * pbp2 * RL -
       2. * (c12 + c22) * LRLR * p1p2 * papb * pbp2 * RL +
       2. * (c12 + c22) * LRLR * pap2 * pbp1 * pbp2 * RL -
       2. * (c12 + c22) * LRLR * pap1 * pow(pbp2, 2) * RL + 6. * c00 * LR * m1 * m2 * papb * RLLL +
       c22 * LR * m1 * m2 * m2s * papb * RLLL -
       2. * (c12 + c22) * LR * m1 * m2 * papb * pbp2 * RLLL +
       LR *
           ((c22 * m2s - 2. * (c12 + c22) * pbp2) * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) +
            4. * c00 * (p1p2 * papb - pap2 * pbp1 + 2. * pap1 * pbp2)) *
           RLRL +
       c2 * (LRRR * m1 * m2 * m2s * papb * RL + LRLR * m2s * p1p2 * papb * RL -
             LRLR * m2s * pap2 * pbp1 * RL + LRLR * m2s * pap1 * pbp2 * RL -
             2. * LRRR * m1 * m2 * papb * pbp2 * RL - 2. * LRLR * p1p2 * papb * pbp2 * RL +
             2. * LRLR * pap2 * pbp1 * pbp2 * RL - 2. * LRLR * pap1 * pow(pbp2, 2) * RL +
             LR * m1 * m2 * m2s * papb * RLLL - 2. * LR * m1 * m2 * papb * pbp2 * RLLL +
             LR * m2s * p1p2 * papb * RLRL - LR * m2s * pap2 * pbp1 * RLRL +
             LR * m2s * pap1 * pbp2 * RLRL - 2. * LR * p1p2 * papb * pbp2 * RLRL +
             2. * LR * pap2 * pbp1 * pbp2 * RLRL - 2. * LR * pap1 * pow(pbp2, 2) * RLRL +
             LLLL * m1 * m2s * papb * (LR * ml3 + ml1 * RR) +
             LLRL * m2 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LR * ml3 + ml1 * RR) +
             LL * m2 * ml1 * p1p2 * papb * RRLR - LL * m2 * ml1 * pap2 * pbp1 * RRLR +
             LL * m2 * ml1 * pap1 * pbp2 * RRLR + m2 * ml3 * p1p2 * papb * RL * RRLR -
             m2 * ml3 * pap2 * pbp1 * RL * RRLR + m2 * ml3 * pap1 * pbp2 * RL * RRLR +
             m1 * m2s * papb * (LL * ml1 + ml3 * RL) * RRRR) +
       c0 * ml1 *
           ((LLLL * m1 * m2s * papb +
             m2 * (LLRL * p1p2 * papb - LLRL * pap2 * pbp1 + LLRL * pap1 * pbp2 +
                   m1 * ml3 * papb * RLLL) +
             ml3 * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RLRL) *
                RR +
            LL * (LRRR * m1 * m2 * ml3 * papb +
                  (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * (LRLR * ml3 + m2 * RRLR) +
                  m1 * m2s * papb * RRRR)));

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      -8. * c00 * ivs3v2 * ivt1s1 *
      (LRRR * m1 * m2 * papb * RL + LRLR * (p1p2 * papb - pap2 * pbp1 + pap1 * pbp2) * RL +
       LR * (m1 * m2 * papb * RLLL + p1p2 * papb * RLRL - pap2 * pbp1 * RLRL + pap1 * pbp2 * RLRL));

  return real(me0 + me1 + me2);
}

double FI::MVttr4t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //    Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      +c0 * ivt1s1 * ivt1s2 *
      (-4. * pap1 * pbp2 * ml1 * ml3 * RR * RLRL - 4. * pap1 * pbp2 * ml1 * ml3 * RR * RRRR -
       4. * pap1 * pbp2 * ml1 * ml3 * LL * LRLR - 4. * pap1 * pbp2 * ml1 * ml3 * LL * LLLL -
       4. * m2 * pap1 * pbp2 * ml1 * RR * LRRR - 4. * m2 * pap1 * pbp2 * ml1 * RR * LLRL -
       4. * m2 * pap1 * pbp2 * ml1 * LL * RRLR - 4. * m2 * pap1 * pbp2 * ml1 * LL * RLLL);
  me0 += +c2 * ivt1s1 * ivt1s2 *
         (-4. * pap1 * pbp2 * RL * LRLR * m2s - 4. * pap1 * pbp2 * RL * LLLL * m2s -
          4. * pap1 * pbp2 * LR * RLRL * m2s - 4. * pap1 * pbp2 * LR * RRRR * m2s +
          8. * pap1 * pow(pbp2, 2) * RL * LRLR + 8. * pap1 * pow(pbp2, 2) * RL * LLLL +
          8. * pap1 * pow(pbp2, 2) * LR * RLRL + 8. * pap1 * pow(pbp2, 2) * LR * RRRR -
          4. * m2 * pap1 * pbp2 * ml3 * RL * RRLR - 4. * m2 * pap1 * pbp2 * ml3 * RL * RLLL -
          4. * m2 * pap1 * pbp2 * ml3 * LR * LRRR - 4. * m2 * pap1 * pbp2 * ml3 * LR * LLRL -
          4. * m2 * pap1 * pbp2 * ml1 * RR * LRRR - 4. * m2 * pap1 * pbp2 * ml1 * RR * LLRL -
          4. * m2 * pap1 * pbp2 * ml1 * LL * RRLR - 4. * m2 * pap1 * pbp2 * ml1 * LL * RLLL);
  me0 += +c22 * ivt1s1 * ivt1s2 *
         (-4. * pap1 * pbp2 * RL * LRLR * m2s - 4. * pap1 * pbp2 * RL * LLLL * m2s -
          4. * pap1 * pbp2 * LR * RLRL * m2s - 4. * pap1 * pbp2 * LR * RRRR * m2s);
  me0 += +(c12 + c22) * ivt1s1 * ivt1s2 *
         (8. * pap1 * pow(pbp2, 2) * RL * LRLR + 8. * pap1 * pow(pbp2, 2) * RL * LLLL +
          8. * pap1 * pow(pbp2, 2) * LR * RLRL + 8. * pap1 * pow(pbp2, 2) * LR * RRRR);
  me0 += +c00 * ivt1s1 * ivt1s2 *
         (-16. * pap1 * pbp2 * RL * LRLR - 16. * pap1 * pbp2 * RL * LLLL -
          16. * pap1 * pbp2 * LR * RLRL - 16. * pap1 * pbp2 * LR * RRRR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 = +c00 * ivt1s1 * ivt1s2 *
                        (8. * pap1 * pbp2 * RL * LRLR + 8. * pap1 * pbp2 * RL * LLLL +
                         8. * pap1 * pbp2 * LR * RLRL + 8. * pap1 * pbp2 * LR * RRRR);

  return real(me0 + me1);
}

double FI::MVttr4u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  //  Npf *cc = new Npf(0., m1s - 2.*pap1, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      +c0 * ivt1s1 * ivu2s2 *
      (2. * pap2 * pbp1 * ml1 * ml3 * RR * RRRR + 2. * pap2 * pbp1 * ml1 * ml3 * LL * LLLL +
       2. * pap1 * pbp2 * ml1 * ml3 * RR * RRRR + 2. * pap1 * pbp2 * ml1 * ml3 * LL * LLLL -
       2. * papb * p1p2 * ml1 * ml3 * RR * RRRR - 2. * papb * p1p2 * ml1 * ml3 * LL * LLLL +
       2. * m2 * pap2 * pbp1 * ml1 * RR * LRRR + 2. * m2 * pap2 * pbp1 * ml1 * LL * RLLL +
       2. * m2 * pap1 * pbp2 * ml1 * RR * LRRR + 2. * m2 * pap1 * pbp2 * ml1 * LL * RLLL -
       2. * m2 * papb * p1p2 * ml1 * RR * LRRR - 2. * m2 * papb * p1p2 * ml1 * LL * RLLL +
       2. * m1 * papb * ml1 * RR * LLRL * m2s + 2. * m1 * papb * ml1 * LL * RRLR * m2s +
       2. * m1 * m2 * papb * ml1 * ml3 * RR * RLRL + 2. * m1 * m2 * papb * ml1 * ml3 * LL * LRLR);
  me0 += +c2 * ivt1s1 * ivu2s2 *
         (2. * pap2 * pbp1 * RL * LLLL * m2s + 2. * pap2 * pbp1 * LR * RRRR * m2s -
          4. * pap2 * pbp1 * pbp2 * RL * LLLL - 4. * pap2 * pbp1 * pbp2 * LR * RRRR +
          2. * pap1 * pbp2 * RL * LLLL * m2s + 2. * pap1 * pbp2 * LR * RRRR * m2s -
          4. * pap1 * pow(pbp2, 2) * RL * LLLL - 4. * pap1 * pow(pbp2, 2) * LR * RRRR -
          2. * papb * p1p2 * RL * LLLL * m2s - 2. * papb * p1p2 * LR * RRRR * m2s +
          4. * papb * pbp2 * p1p2 * RL * LLLL + 4. * papb * pbp2 * p1p2 * LR * RRRR +
          2. * m2 * pap2 * pbp1 * ml3 * RL * RLLL + 2. * m2 * pap2 * pbp1 * ml3 * LR * LRRR +
          2. * m2 * pap2 * pbp1 * ml1 * RR * LRRR + 2. * m2 * pap2 * pbp1 * ml1 * LL * RLLL +
          2. * m2 * pap1 * pbp2 * ml3 * RL * RLLL + 2. * m2 * pap1 * pbp2 * ml3 * LR * LRRR +
          2. * m2 * pap1 * pbp2 * ml1 * RR * LRRR + 2. * m2 * pap1 * pbp2 * ml1 * LL * RLLL -
          2. * m2 * papb * p1p2 * ml3 * RL * RLLL - 2 * m2 * papb * p1p2 * ml3 * LR * LRRR -
          2. * m2 * papb * p1p2 * ml1 * RR * LRRR - 2. * m2 * papb * p1p2 * ml1 * LL * RLLL +
          2. * m1 * papb * ml3 * RL * RRLR * m2s + 2. * m1 * papb * ml3 * LR * LLRL * m2s +
          2. * m1 * papb * ml1 * RR * LLRL * m2s + 2. * m1 * papb * ml1 * LL * RRLR * m2s +
          2. * m1 * m2 * papb * RL * LRLR * m2s + 2. * m1 * m2 * papb * LR * RLRL * m2s -
          4. * m1 * m2 * papb * pbp2 * RL * LRLR - 4. * m1 * m2 * papb * pbp2 * LR * RLRL);
  me0 += +c22 * ivt1s1 * ivu2s2 *
         (2. * pap2 * pbp1 * RL * LLLL * m2s + 2. * pap2 * pbp1 * LR * RRRR * m2s +
          2. * pap1 * pbp2 * RL * LLLL * m2s + 2. * pap1 * pbp2 * LR * RRRR * m2s -
          2. * papb * p1p2 * RL * LLLL * m2s - 2. * papb * p1p2 * LR * RRRR * m2s +
          2. * m1 * m2 * papb * RL * LRLR * m2s + 2. * m1 * m2 * papb * LR * RLRL * m2s);
  me0 += +(c12 + c22) * ivt1s1 * ivu2s2 *
         (-4. * pap2 * pbp1 * pbp2 * RL * LLLL - 4. * pap2 * pbp1 * pbp2 * LR * RRRR -
          4. * pap1 * pow(pbp2, 2) * RL * LLLL - 4. * pap1 * pow(pbp2, 2) * LR * RRRR +
          4. * papb * pbp2 * p1p2 * RL * LLLL + 4. * papb * pbp2 * p1p2 * LR * RRRR -
          4. * m1 * m2 * papb * pbp2 * RL * LRLR - 4. * m1 * m2 * papb * pbp2 * LR * RLRL);
  me0 +=
      +c00 * ivt1s1 * ivu2s2 *
      (8. * pap2 * pbp1 * RL * LLLL + 8. * pap2 * pbp1 * LR * RRRR + 8. * pap1 * pbp2 * RL * LLLL +
       8. * pap1 * pbp2 * LR * RRRR - 8. * papb * p1p2 * RL * LLLL - 8. * papb * p1p2 * LR * RRRR +
       8. * m1 * m2 * papb * RL * LRLR + 8. * m1 * m2 * papb * LR * RLRL);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      +c00 * ivt1s1 * ivu2s2 *
      (-4. * pap2 * pbp1 * RL * LLLL - 4. * pap2 * pbp1 * LR * RRRR - 4. * pap1 * pbp2 * RL * LLLL -
       4. * pap1 * pbp2 * LR * RRRR + 4. * papb * p1p2 * RL * LLLL + 4. * papb * p1p2 * LR * RRRR -
       4. * m1 * m2 * papb * RL * LRLR - 4. * m1 * m2 * papb * LR * RLRL);

  return real(me0 + me1);
}

double FI::MVutr3s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;

  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      +ivs3v2 * ivu2s1 *
      (8. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * LRRR +
       8. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * RLLL +
       8. * c0 * m2 * pap2 * pbp1 * ml1 * RR * LLRR + 8. * c0 * m2 * pap2 * pbp1 * ml1 * LL * RRLL +
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * LRLR +
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL +
       4. * c0 * m1 * m2s * papb * ml1 * RR * LLLR + 4. * c0 * m1 * m2s * papb * ml1 * LL * RRRL -
       16. * c2 * pow(pap2, 2) * pbp1 * RL * RLLL - 16. * c2 * pow(pap2, 2) * pbp1 * LR * LRRR +
       8. * c2 * m2 * pap2 * pbp1 * ml3 * RL * RRLL + 8. * c2 * m2 * pap2 * pbp1 * ml3 * LR * LLRR +
       8. * c2 * m2 * pap2 * pbp1 * ml1 * RR * LLRR + 8. * c2 * m2 * pap2 * pbp1 * ml1 * LL * RRLL +
       8. * c2 * m2s * pap2 * pbp1 * RL * RLLL + 8. * c2 * pow(m2, 2) * pap2 * pbp1 * LR * LRRR -
       8. * c2 * m1 * m2 * papb * pap2 * RL * RLRL - 8. * c2 * m1 * m2 * papb * pap2 * LR * LRLR +
       4. * c2 * m1 * m2s * papb * ml3 * RL * RRRL + 4. * c2 * m1 * m2s * papb * ml3 * LR * LLLR +
       4. * c2 * m1 * pow(m2, 2) * papb * ml1 * RR * LLLR +
       4. * c2 * m1 * m2s * papb * ml1 * LL * RRRL + 4. * c2 * m1 * pow(m2, 3) * papb * RL * RLRL +
       4. * c2 * m1 * pow(m2, 3) * papb * LR * LRLR + 8. * c22 * m2s * pap2 * pbp1 * RL * RLLL +
       8. * c22 * pow(m2, 2) * pap2 * pbp1 * LR * LRRR);
  me0 += +ivs3v2 * ivu2s1 *
         (4. * c22 * m1 * pow(m2, 3) * papb * RL * RLRL +
          4. * c22 * m1 * pow(m2, 3) * papb * LR * LRLR -
          16. * (c12 + c22) * pow(pap2, 2) * pbp1 * RL * RLLL -
          16. * (c12 + c22) * pow(pap2, 2) * pbp1 * LR * LRRR -
          8. * (c12 + c22) * m1 * m2 * papb * pap2 * RL * RLRL -
          8. * (c12 + c22) * m1 * m2 * papb * pap2 * LR * LRLR +
          32. * c00 * pap2 * pbp1 * RL * RLLL + 32. * c00 * pap2 * pbp1 * LR * LRRR +
          16. * c00 * m1 * m2 * papb * RL * RLRL + 16. * c00 * m1 * m2 * papb * LR * LRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      +ivs3v2 * ivu2s1 *
      (-4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * LRRR -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * RLLL +
       4. * c0 * pap1 * pbp2 * ml1 * ml3 * RR * LRRR +
       4. * c0 * pap1 * pbp2 * ml1 * ml3 * LL * RLLL -
       4. * c0 * papb * p1p2 * ml1 * ml3 * RR * LRRR -
       4. * c0 * papb * p1p2 * ml1 * ml3 * LL * RLLL -
       4. * c0 * m2 * pap2 * pbp1 * ml1 * RR * LLRR - 4. * c0 * m2 * pap2 * pbp1 * ml1 * LL * RRLL +
       4. * c0 * m2 * pap1 * pbp2 * ml1 * RR * LLRR + 4. * c0 * m2 * pap1 * pbp2 * ml1 * LL * RRLL -
       4. * c0 * m2 * papb * p1p2 * ml1 * RR * LLRR - 4. * c0 * m2 * papb * p1p2 * ml1 * LL * RRLL -
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * LRLR -
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL -
       4. * c0 * m1 * m2s * papb * ml1 * RR * LLLR - 4. * c0 * m1 * m2s * papb * ml1 * LL * RRRL +
       8. * c2 * pow(pap2, 2) * pbp1 * RL * RLLL + 8. * c2 * pow(pap2, 2) * pbp1 * LR * LRRR -
       8. * c2 * pap1 * pap2 * pbp2 * RL * RLLL - 8. * c2 * pap1 * pap2 * pbp2 * LR * LRRR +
       8. * c2 * papb * pap2 * p1p2 * RL * RLLL + 8. * c2 * papb * pap2 * p1p2 * LR * LRRR -
       4. * c2 * m2 * pap2 * pbp1 * ml3 * RL * RRLL - 4. * c2 * m2 * pap2 * pbp1 * ml3 * LR * LLRR -
       4. * c2 * m2 * pap2 * pbp1 * ml1 * RR * LLRR - 4. * c2 * m2 * pap2 * pbp1 * ml1 * LL * RRLL +
       4. * c2 * m2 * pap1 * pbp2 * ml3 * RL * RRLL);
  me1 +=
      +ivs3v2 * ivu2s1 *
      (4. * c2 * m2 * pap1 * pbp2 * ml3 * LR * LLRR + 4. * c2 * m2 * pap1 * pbp2 * ml1 * RR * LLRR +
       4. * c2 * m2 * pap1 * pbp2 * ml1 * LL * RRLL - 4. * c2 * m2 * papb * p1p2 * ml3 * RL * RRLL -
       4. * c2 * m2 * papb * p1p2 * ml3 * LR * LLRR - 4. * c2 * m2 * papb * p1p2 * ml1 * RR * LLRR -
       4. * c2 * m2 * papb * p1p2 * ml1 * LL * RRLL - 4. * c2 * m2s * pap2 * pbp1 * RL * RLLL -
       4. * c2 * m2s * pap2 * pbp1 * LR * LRRR + 4. * c2 * m2s * pap1 * pbp2 * RL * RLLL +
       4. * c2 * m2s * pap1 * pbp2 * LR * LRRR - 4. * c2 * pow(m2, 2) * papb * p1p2 * RL * RLLL -
       4. * c2 * m2s * papb * p1p2 * LR * LRRR + 8. * c2 * m1 * m2 * papb * pap2 * RL * RLRL +
       8. * c2 * m1 * m2 * papb * pap2 * LR * LRLR - 4. * c2 * m1 * m2s * papb * ml3 * RL * RRRL -
       4. * c2 * m1 * pow(m2, 2) * papb * ml3 * LR * LLLR -
       4. * c2 * m1 * m2s * papb * ml1 * RR * LLLR - 4. * c2 * m1 * m2s * papb * ml1 * LL * RRRL -
       4. * c2 * m1 * pow(m2, 3) * papb * RL * RLRL - 4. * c2 * m1 * pow(m2, 3) * papb * LR * LRLR -
       4. * c22 * pow(m2, 2) * pap2 * pbp1 * RL * RLLL - 4. * c22 * m2s * pap2 * pbp1 * LR * LRRR +
       4. * c22 * m2s * pap1 * pbp2 * RL * RLLL + 4. * c22 * m2s * pap1 * pbp2 * LR * LRRR);
  me1 += +ivs3v2 * ivu2s1 *
         (-4. * c22 * m2s * papb * p1p2 * RL * RLLL - 4. * c22 * m2s * papb * p1p2 * LR * LRRR -
          4. * c22 * m1 * pow(m2, 3) * papb * RL * RLRL -
          4. * c22 * m1 * pow(m2, 3) * papb * LR * LRLR +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * RL * RLLL +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * LR * LRRR -
          8. * (c12 + c22) * pap1 * pap2 * pbp2 * RL * RLLL -
          8. * (c12 + c22) * pap1 * pap2 * pbp2 * LR * LRRR +
          8. * (c12 + c22) * papb * pap2 * p1p2 * RL * RLLL +
          8. * (c12 + c22) * papb * pap2 * p1p2 * LR * LRRR +
          8. * (c12 + c22) * m1 * m2 * papb * pap2 * RL * RLRL +
          8. * (c12 + c22) * m1 * m2 * papb * pap2 * LR * LRLR -
          32. * c00 * pap2 * pbp1 * RL * RLLL - 32. * c00 * pap2 * pbp1 * LR * LRRR +
          16. * c00 * pap1 * pbp2 * RL * RLLL + 16. * c00 * pap1 * pbp2 * LR * LRRR -
          16. * c00 * papb * p1p2 * RL * RLLL - 16. * c00 * papb * p1p2 * LR * LRRR -
          24. * c00 * m1 * m2 * papb * RL * RLRL - 24. * c00 * m1 * m2 * papb * LR * LRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me2 =
      +ivs3v2 * ivu2s1 *
      (8. * c00 * pap2 * pbp1 * RL * RLLL + 8. * c00 * pap2 * pbp1 * LR * LRRR -
       8. * c00 * pap1 * pbp2 * RL * RLLL - 8. * c00 * pap1 * pbp2 * LR * LRRR +
       8. * c00 * papb * p1p2 * RL * RLLL + 8. * c00 * papb * p1p2 * LR * LRRR +
       8. * c00 * m1 * m2 * papb * RL * RLRL + 8. * c00 * m1 * m2 * papb * LR * LRLR);

  return real(me0 + me1 + me2);
}

double FI::MVutr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;

  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      +ivt1s2 * ivu2s1 *
      (2. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR +
       2. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL +
       2. * c0 * pap1 * pbp2 * ml1 * ml3 * RR * RRRR +
       2. * c0 * pap1 * pbp2 * ml1 * ml3 * LL * LLLL -
       2. * c0 * papb * p1p2 * ml1 * ml3 * RR * RRRR -
       2. * c0 * papb * p1p2 * ml1 * ml3 * LL * LLLL +
       2. * c0 * m2 * pap2 * pbp1 * ml1 * RR * RLRR + 2. * c0 * m2 * pap2 * pbp1 * ml1 * LL * LRLL +
       2. * c0 * m2 * pap1 * pbp2 * ml1 * RR * RLRR + 2. * c0 * m2 * pap1 * pbp2 * ml1 * LL * LRLL -
       2. * c0 * m2 * papb * p1p2 * ml1 * RR * RLRR - 2. * c0 * m2 * papb * p1p2 * ml1 * LL * LRLL +
       2. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * LRLR +
       2. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * RLRL +
       2. * c0 * m1 * m2s * papb * ml1 * RR * LLLR + 2. * c0 * m1 * m2s * papb * ml1 * LL * RRRL -
       4. * c2 * pow(pap2, 2) * pbp1 * RL * LLLL - 4. * c2 * pow(pap2, 2) * pbp1 * LR * RRRR -
       4. * c2 * pap1 * pap2 * pbp2 * RL * LLLL - 4. * c2 * pap1 * pap2 * pbp2 * LR * RRRR +
       4. * c2 * papb * pap2 * p1p2 * RL * LLLL + 4. * c2 * papb * pap2 * p1p2 * LR * RRRR +
       2. * c2 * m2 * pap2 * pbp1 * ml3 * RL * LRLL + 2. * c2 * m2 * pap2 * pbp1 * ml3 * LR * RLRR +
       2. * c2 * m2 * pap2 * pbp1 * ml1 * RR * RLRR + 2. * c2 * m2 * pap2 * pbp1 * ml1 * LL * LRLL +
       2. * c2 * m2 * pap1 * pbp2 * ml3 * RL * LRLL);
  me0 +=
      +ivt1s2 * ivu2s1 *
      (2. * c2 * m2 * pap1 * pbp2 * ml3 * LR * RLRR + 2. * c2 * m2 * pap1 * pbp2 * ml1 * RR * RLRR +
       2. * c2 * m2 * pap1 * pbp2 * ml1 * LL * LRLL - 2. * c2 * m2 * papb * p1p2 * ml3 * RL * LRLL -
       2. * c2 * m2 * papb * p1p2 * ml3 * LR * RLRR - 2. * c2 * m2 * papb * p1p2 * ml1 * RR * RLRR -
       2. * c2 * m2 * papb * p1p2 * ml1 * LL * LRLL + 2. * c2 * m2s * pap2 * pbp1 * RL * LLLL +
       2. * c2 * m2s * pap2 * pbp1 * LR * RRRR + 2. * c2 * m2s * pap1 * pbp2 * RL * LLLL +
       2. * c2 * m2s * pap1 * pbp2 * LR * RRRR - 2. * c2 * pow(m2, 2) * papb * p1p2 * RL * LLLL -
       2. * c2 * m2s * papb * p1p2 * LR * RRRR - 4. * c2 * m1 * m2 * papb * pap2 * RL * RLRL -
       4. * c2 * m1 * m2 * papb * pap2 * LR * LRLR + 2. * c2 * m1 * m2s * papb * ml3 * RL * RRRL +
       2. * c2 * m1 * pow(m2, 2) * papb * ml3 * LR * LLLR +
       2. * c2 * m1 * m2s * papb * ml1 * RR * LLLR + 2. * c2 * m1 * m2s * papb * ml1 * LL * RRRL +
       2. * c2 * m1 * pow(m2, 3) * papb * RL * RLRL + 2. * c2 * m1 * pow(m2, 3) * papb * LR * LRLR +
       2. * c22 * pow(m2, 2) * pap2 * pbp1 * RL * LLLL + 2. * c22 * m2s * pap2 * pbp1 * LR * RRRR +
       2. * c22 * m2s * pap1 * pbp2 * RL * LLLL + 2. * c22 * m2s * pap1 * pbp2 * LR * RRRR);
  me0 += +ivt1s2 * ivu2s1 *
         (-2. * c22 * m2s * papb * p1p2 * RL * LLLL - 2. * c22 * m2s * papb * p1p2 * LR * RRRR +
          2. * c22 * m1 * pow(m2, 3) * papb * RL * RLRL +
          2. * c22 * m1 * pow(m2, 3) * papb * LR * LRLR -
          4. * (c12 + c22) * pow(pap2, 2) * pbp1 * RL * LLLL -
          4. * (c12 + c22) * pow(pap2, 2) * pbp1 * LR * RRRR -
          4. * (c12 + c22) * pap1 * pap2 * pbp2 * RL * LLLL -
          4. * (c12 + c22) * pap1 * pap2 * pbp2 * LR * RRRR +
          4. * (c12 + c22) * papb * pap2 * p1p2 * RL * LLLL +
          4. * (c12 + c22) * papb * pap2 * p1p2 * LR * RRRR -
          4. * (c12 + c22) * m1 * m2 * papb * pap2 * RL * RLRL -
          4. * (c12 + c22) * m1 * m2 * papb * pap2 * LR * LRLR +
          8. * c00 * pap2 * pbp1 * RL * LLLL + 8. * c00 * pap2 * pbp1 * LR * RRRR +
          8. * c00 * pap1 * pbp2 * RL * LLLL + 8. * c00 * pap1 * pbp2 * LR * RRRR -
          8. * c00 * papb * p1p2 * RL * LLLL - 8. * c00 * papb * p1p2 * LR * RRRR +
          8. * c00 * m1 * m2 * papb * RL * RLRL + 8. * c00 * m1 * m2 * papb * LR * LRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 =
      +ivt1s2 * ivu2s1 *
      (-4. * c00 * pap2 * pbp1 * RL * LLLL - 4. * c00 * pap2 * pbp1 * LR * RRRR -
       4. * c00 * pap1 * pbp2 * RL * LLLL - 4. * c00 * pap1 * pbp2 * LR * RRRR +
       4. * c00 * papb * p1p2 * RL * LLLL + 4. * c00 * papb * p1p2 * LR * RRRR -
       4. * c00 * m1 * m2 * papb * RL * RLRL - 4. * c00 * m1 * m2 * papb * LR * LRLR);

  return real(me0 + me1);
}

double FI::MVutr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m2s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me0 =
      +ivu2s1 * ivu2s2 *
      (-4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * LRLR -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * RLRL -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL -
       4. * c0 * m2 * pap2 * pbp1 * ml1 * RR * RLRR - 4. * c0 * m2 * pap2 * pbp1 * ml1 * RR * LLLR -
       4. * c0 * m2 * pap2 * pbp1 * ml1 * LL * RRRL - 4. * c0 * m2 * pap2 * pbp1 * ml1 * LL * LRLL +
       8. * c2 * pow(pap2, 2) * pbp1 * RL * RLRL + 8. * c2 * pow(pap2, 2) * pbp1 * RL * LLLL +
       8. * c2 * pow(pap2, 2) * pbp1 * LR * RRRR + 8. * c2 * pow(pap2, 2) * pbp1 * LR * LRLR -
       4. * c2 * m2 * pap2 * pbp1 * ml3 * RL * RRRL - 4. * c2 * m2 * pap2 * pbp1 * ml3 * RL * LRLL -
       4. * c2 * m2 * pap2 * pbp1 * ml3 * LR * RLRR - 4. * c2 * m2 * pap2 * pbp1 * ml3 * LR * LLLR -
       4. * c2 * m2 * pap2 * pbp1 * ml1 * RR * RLRR - 4. * c2 * m2 * pap2 * pbp1 * ml1 * RR * LLLR -
       4. * c2 * m2 * pap2 * pbp1 * ml1 * LL * RRRL - 4. * c2 * m2 * pap2 * pbp1 * ml1 * LL * LRLL -
       4. * c2 * m2s * pap2 * pbp1 * RL * RLRL - 4. * c2 * m2s * pap2 * pbp1 * RL * LLLL -
       4. * c2 * m2s * pap2 * pbp1 * LR * RRRR - 4. * c2 * m2s * pap2 * pbp1 * LR * LRLR -
       4. * c22 * pow(m2, 2) * pap2 * pbp1 * RL * RLRL - 4. * c22 * m2s * pap2 * pbp1 * RL * LLLL);

  me0 += ivu2s1 * ivu2s2 *
         (-4. * c22 * m2s * pap2 * pbp1 * LR * RRRR - 4. * c22 * m2s * pap2 * pbp1 * LR * LRLR +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * RL * RLRL +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * RL * LLLL +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * LR * RRRR +
          8. * (c12 + c22) * pow(pap2, 2) * pbp1 * LR * LRLR - 16. * c00 * pap2 * pbp1 * RL * RLRL -
          16. * c00 * pap2 * pbp1 * RL * LLLL - 16. * c00 * pap2 * pbp1 * LR * RRRR -
          16. * c00 * pap2 * pbp1 * LR * LRLR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  complex<double> me1 = +ivu2s1 * ivu2s2 *
                        (8. * c00 * pap2 * pbp1 * RL * RLRL + 8. * c00 * pap2 * pbp1 * RL * LLLL +
                         8. * c00 * pap2 * pbp1 * LR * RRRR + 8. * c00 * pap2 * pbp1 * LR * LRLR);

  return real(me0);
}

double FI::MVutr4s(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  // Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

  complex<double> me0;
  ;
  complex<double> me1;
  ;
  complex<double> me2;
  ;

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me0 =
      +ivs3v2 * ivu2s1 *
      (8. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RLLL +
       8. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LRRR +
       8. * c0 * m1 * pap2 * pbp1 * ml1 * RR * LLLL + 8. * c0 * m1 * pap2 * pbp1 * ml1 * LL * RRRR +
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * RLRL +
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR +
       4. * c0 * m1s * m2 * papb * ml1 * RR * LLRL + 4. * c0 * m1s * m2 * papb * ml1 * LL * RRLR -
       16. * c2 * pap2 * pow(pbp1, 2) * RL * LRRR - 16. * c2 * pap2 * pow(pbp1, 2) * LR * RLLL +
       8. * c2 * m1 * pap2 * pbp1 * ml3 * RL * RRRR + 8. * c2 * m1 * pap2 * pbp1 * ml3 * LR * LLLL +
       8. * c2 * m1 * pap2 * pbp1 * ml1 * RR * LLLL + 8. * c2 * m1 * pap2 * pbp1 * ml1 * LL * RRRR -
       8. * c2 * m1 * m2 * papb * pbp1 * RL * LRLR - 8. * c2 * m1 * m2 * papb * pbp1 * LR * RLRL +
       8. * c2 * m1s * pap2 * pbp1 * RL * LRRR + 8. * c2 * m1s * pap2 * pbp1 * LR * RLLL +
       4. * c2 * m1s * m2 * papb * ml3 * RL * RRLR + 4. * c2 * m1s * m2 * papb * ml3 * LR * LLRL +
       4. * c2 * pow(m1, 2) * m2 * papb * ml1 * RR * LLRL +
       4. * c2 * m1s * m2 * papb * ml1 * LL * RRLR + 4. * c2 * pow(m1, 3) * m2 * papb * RL * LRLR +
       4. * c2 * pow(m1, 3) * m2 * papb * LR * RLRL + 8. * c22 * m1s * pap2 * pbp1 * RL * LRRR +
       8. * c22 * pow(m1, 2) * pap2 * pbp1 * LR * RLLL);
  me0 += +ivs3v2 * ivu2s1 *
         (4. * c22 * pow(m1, 3) * m2 * papb * RL * LRLR +
          4. * c22 * pow(m1, 3) * m2 * papb * LR * RLRL -
          16. * (c12 + c22) * pap2 * pow(pbp1, 2) * RL * LRRR -
          16. * (c12 + c22) * pap2 * pow(pbp1, 2) * LR * RLLL -
          8. * (c12 + c22) * m1 * m2 * papb * pbp1 * RL * LRLR -
          8. * (c12 + c22) * m1 * m2 * papb * pbp1 * LR * RLRL +
          32. * c00 * pap2 * pbp1 * RL * LRRR + 32. * c00 * pap2 * pbp1 * LR * RLLL +
          16. * c00 * m1 * m2 * papb * RL * LRLR + 16. * c00 * m1 * m2 * papb * LR * RLRL);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me1 =
      +ivs3v2 * ivu2s1 *
      (-4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RLLL -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LRRR +
       4. * c0 * pap1 * pbp2 * ml1 * ml3 * RR * RLLL +
       4. * c0 * pap1 * pbp2 * ml1 * ml3 * LL * LRRR -
       4. * c0 * papb * p1p2 * ml1 * ml3 * RR * RLLL -
       4. * c0 * papb * p1p2 * ml1 * ml3 * LL * LRRR -
       4. * c0 * m1 * pap2 * pbp1 * ml1 * RR * LLLL - 4. * c0 * m1 * pap2 * pbp1 * ml1 * LL * RRRR +
       4. * c0 * m1 * pap1 * pbp2 * ml1 * RR * LLLL + 4. * c0 * m1 * pap1 * pbp2 * ml1 * LL * RRRR -
       4. * c0 * m1 * papb * p1p2 * ml1 * RR * LLLL - 4. * c0 * m1 * papb * p1p2 * ml1 * LL * RRRR -
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * RLRL -
       4. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR -
       4. * c0 * m1s * m2 * papb * ml1 * RR * LLRL - 4. * c0 * m1s * m2 * papb * ml1 * LL * RRLR +
       8. * c2 * pap2 * pow(pbp1, 2) * RL * LRRR + 8. * c2 * pap2 * pow(pbp1, 2) * LR * RLLL -
       8. * c2 * pap1 * pbp1 * pbp2 * RL * LRRR - 8. * c2 * pap1 * pbp1 * pbp2 * LR * RLLL +
       8. * c2 * papb * pbp1 * p1p2 * RL * LRRR + 8. * c2 * papb * pbp1 * p1p2 * LR * RLLL -
       4. * c2 * m1 * pap2 * pbp1 * ml3 * RL * RRRR - 4. * c2 * m1 * pap2 * pbp1 * ml3 * LR * LLLL -
       4. * c2 * m1 * pap2 * pbp1 * ml1 * RR * LLLL - 4. * c2 * m1 * pap2 * pbp1 * ml1 * LL * RRRR +
       4. * c2 * m1 * pap1 * pbp2 * ml3 * RL * RRRR);
  me1 +=
      +ivs3v2 * ivu2s1 *
      (4. * c2 * m1 * pap1 * pbp2 * ml3 * LR * LLLL + 4. * c2 * m1 * pap1 * pbp2 * ml1 * RR * LLLL +
       4. * c2 * m1 * pap1 * pbp2 * ml1 * LL * RRRR - 4. * c2 * m1 * papb * p1p2 * ml3 * RL * RRRR -
       4. * c2 * m1 * papb * p1p2 * ml3 * LR * LLLL - 4. * c2 * m1 * papb * p1p2 * ml1 * RR * LLLL -
       4. * c2 * m1 * papb * p1p2 * ml1 * LL * RRRR + 8. * c2 * m1 * m2 * papb * pbp1 * RL * LRLR +
       8. * c2 * m1 * m2 * papb * pbp1 * LR * RLRL - 4. * c2 * m1s * pap2 * pbp1 * RL * LRRR -
       4. * c2 * m1s * pap2 * pbp1 * LR * RLLL + 4. * c2 * m1s * pap1 * pbp2 * RL * LRRR +
       4. * c2 * m1s * pap1 * pbp2 * LR * RLLL - 4. * c2 * pow(m1, 2) * papb * p1p2 * RL * LRRR -
       4. * c2 * m1s * papb * p1p2 * LR * RLLL - 4. * c2 * m1s * m2 * papb * ml3 * RL * RRLR -
       4. * c2 * m1s * m2 * papb * ml3 * LR * LLRL - 4. * c2 * m1s * m2 * papb * ml1 * RR * LLRL -
       4. * c2 * m1s * m2 * papb * ml1 * LL * RRLR - 4. * c2 * pow(m1, 3) * m2 * papb * RL * LRLR -
       4. * c2 * pow(m1, 3) * m2 * papb * LR * RLRL - 4. * c22 * m1s * pap2 * pbp1 * RL * LRRR -
       4. * c22 * m1s * pap2 * pbp1 * LR * RLLL + 4. * c22 * m1s * pap1 * pbp2 * RL * LRRR +
       4. * c22 * m1s * pap1 * pbp2 * LR * RLLL);
  me1 += +ivs3v2 * ivu2s1 *
         (-4. * c22 * m1s * papb * p1p2 * RL * LRRR - 4. * c22 * m1s * papb * p1p2 * LR * RLLL -
          4. * c22 * pow(m1, 3) * m2 * papb * RL * LRLR -
          4. * c22 * pow(m1, 3) * m2 * papb * LR * RLRL +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * RL * LRRR +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * LR * RLLL -
          8. * (c12 + c22) * pap1 * pbp1 * pbp2 * RL * LRRR -
          8. * (c12 + c22) * pap1 * pbp1 * pbp2 * LR * RLLL +
          8. * (c12 + c22) * papb * pbp1 * p1p2 * RL * LRRR +
          8. * (c12 + c22) * papb * pbp1 * p1p2 * LR * RLLL +
          8. * (c12 + c22) * m1 * m2 * papb * pbp1 * RL * LRLR +
          8. * (c12 + c22) * m1 * m2 * papb * pbp1 * LR * RLRL -
          32. * c00 * pap2 * pbp1 * RL * LRRR - 32. * c00 * pap2 * pbp1 * LR * RLLL +
          16. * c00 * pap1 * pbp2 * RL * LRRR + 16. * c00 * pap1 * pbp2 * LR * RLLL -
          16. * c00 * papb * p1p2 * RL * LRRR - 16. * c00 * papb * p1p2 * LR * RLLL -
          24. * c00 * m1 * m2 * papb * RL * LRLR - 24. * c00 * m1 * m2 * papb * LR * RLRL);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me2 = +ivs3v2 * ivu2s1 *
        (8. * c00 * pap2 * pbp1 * RL * LRRR + 8. * c00 * pap2 * pbp1 * LR * RLLL -
         8. * c00 * pap1 * pbp2 * RL * LRRR - 8. * c00 * pap1 * pbp2 * LR * RLLL +
         8. * c00 * papb * p1p2 * RL * LRRR + 8. * c00 * papb * p1p2 * LR * RLLL +
         8. * c00 * m1 * m2 * papb * RL * LRLR + 8. * c00 * m1 * m2 * papb * LR * RLRL);

  return real(me0 + me1 + me2);
}

double FI::MVutr4t(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  complex<double> me0;
  complex<double> me1;

  //    Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me0 =
      +ivt1s2 * ivu2s1 *
      (2. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR +
       2. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL +
       2. * c0 * pap1 * pbp2 * ml1 * ml3 * RR * RRRR +
       2. * c0 * pap1 * pbp2 * ml1 * ml3 * LL * LLLL -
       2. * c0 * papb * p1p2 * ml1 * ml3 * RR * RRRR -
       2. * c0 * papb * p1p2 * ml1 * ml3 * LL * LLLL +
       2. * c0 * m1 * pap2 * pbp1 * ml1 * RR * LRRR + 2. * c0 * m1 * pap2 * pbp1 * ml1 * LL * RLLL +
       2. * c0 * m1 * pap1 * pbp2 * ml1 * RR * LRRR + 2. * c0 * m1 * pap1 * pbp2 * ml1 * LL * RLLL -
       2. * c0 * m1 * papb * p1p2 * ml1 * RR * LRRR - 2. * c0 * m1 * papb * p1p2 * ml1 * LL * RLLL +
       2. * c0 * m1 * m2 * papb * ml1 * ml3 * RR * RLRL +
       2. * c0 * m1 * m2 * papb * ml1 * ml3 * LL * LRLR +
       2. * c0 * m1s * m2 * papb * ml1 * RR * LLRL + 2. * c0 * m1s * m2 * papb * ml1 * LL * RRLR -
       4. * c2 * pap2 * pow(pbp1, 2) * RL * LLLL - 4. * c2 * pap2 * pow(pbp1, 2) * LR * RRRR -
       4. * c2 * pap1 * pbp1 * pbp2 * RL * LLLL - 4. * c2 * pap1 * pbp1 * pbp2 * LR * RRRR +
       4. * c2 * papb * pbp1 * p1p2 * RL * LLLL + 4. * c2 * papb * pbp1 * p1p2 * LR * RRRR +
       2. * c2 * m1 * pap2 * pbp1 * ml3 * RL * RLLL + 2. * c2 * m1 * pap2 * pbp1 * ml3 * LR * LRRR +
       2. * c2 * m1 * pap2 * pbp1 * ml1 * RR * LRRR + 2. * c2 * m1 * pap2 * pbp1 * ml1 * LL * RLLL +
       2. * c2 * m1 * pap1 * pbp2 * ml3 * RL * RLLL);
  me0 +=
      +ivt1s2 * ivu2s1 *
      (2. * c2 * m1 * pap1 * pbp2 * ml3 * LR * LRRR + 2. * c2 * m1 * pap1 * pbp2 * ml1 * RR * LRRR +
       2. * c2 * m1 * pap1 * pbp2 * ml1 * LL * RLLL - 2. * c2 * m1 * papb * p1p2 * ml3 * RL * RLLL -
       2. * c2 * m1 * papb * p1p2 * ml3 * LR * LRRR - 2. * c2 * m1 * papb * p1p2 * ml1 * RR * LRRR -
       2. * c2 * m1 * papb * p1p2 * ml1 * LL * RLLL - 4. * c2 * m1 * m2 * papb * pbp1 * RL * LRLR -
       4. * c2 * m1 * m2 * papb * pbp1 * LR * RLRL + 2. * c2 * m1s * pap2 * pbp1 * RL * LLLL +
       2. * c2 * m1s * pap2 * pbp1 * LR * RRRR + 2. * c2 * m1s * pap1 * pbp2 * RL * LLLL +
       2. * c2 * m1s * pap1 * pbp2 * LR * RRRR - 2. * c2 * pow(m1, 2) * papb * p1p2 * RL * LLLL -
       2. * c2 * m1s * papb * p1p2 * LR * RRRR + 2. * c2 * m1s * m2 * papb * ml3 * RL * RRLR +
       2. * c2 * m1s * m2 * papb * ml3 * LR * LLRL + 2. * c2 * m1s * m2 * papb * ml1 * RR * LLRL +
       2. * c2 * m1s * m2 * papb * ml1 * LL * RRLR + 2. * c2 * pow(m1, 3) * m2 * papb * RL * LRLR +
       2. * c2 * pow(m1, 3) * m2 * papb * LR * RLRL + 2. * c22 * m1s * pap2 * pbp1 * RL * LLLL +
       2. * c22 * m1s * pap2 * pbp1 * LR * RRRR + 2. * c22 * m1s * pap1 * pbp2 * RL * LLLL +
       2. * c22 * m1s * pap1 * pbp2 * LR * RRRR);
  me0 += +ivt1s2 * ivu2s1 *
         (-2. * c22 * m1s * papb * p1p2 * RL * LLLL - 2. * c22 * m1s * papb * p1p2 * LR * RRRR +
          2. * c22 * pow(m1, 3) * m2 * papb * RL * LRLR +
          2. * c22 * pow(m1, 3) * m2 * papb * LR * RLRL -
          4. * (c12 + c22) * pap2 * pow(pbp1, 2) * RL * LLLL -
          4. * (c12 + c22) * pap2 * pow(pbp1, 2) * LR * RRRR -
          4. * (c12 + c22) * pap1 * pbp1 * pbp2 * RL * LLLL -
          4. * (c12 + c22) * pap1 * pbp1 * pbp2 * LR * RRRR +
          4. * (c12 + c22) * papb * pbp1 * p1p2 * RL * LLLL +
          4. * (c12 + c22) * papb * pbp1 * p1p2 * LR * RRRR -
          4. * (c12 + c22) * m1 * m2 * papb * pbp1 * RL * LRLR -
          4. * (c12 + c22) * m1 * m2 * papb * pbp1 * LR * RLRL +
          8. * c00 * pap2 * pbp1 * RL * LLLL + 8. * c00 * pap2 * pbp1 * LR * RRRR +
          8. * c00 * pap1 * pbp2 * RL * LLLL + 8. * c00 * pap1 * pbp2 * LR * RRRR -
          8. * c00 * papb * p1p2 * RL * LLLL - 8. * c00 * papb * p1p2 * LR * RRRR +
          8. * c00 * m1 * m2 * papb * RL * LRLR + 8. * c00 * m1 * m2 * papb * LR * RLRL);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me1 = +ivt1s2 * ivu2s1 *
        (-4. * c00 * pap2 * pbp1 * RL * LLLL - 4. * c00 * pap2 * pbp1 * LR * RRRR -
         4. * c00 * pap1 * pbp2 * RL * LLLL - 4. * c00 * pap1 * pbp2 * LR * RRRR +
         4. * c00 * papb * p1p2 * RL * LLLL + 4. * c00 * papb * p1p2 * LR * RRRR -
         4. * c00 * m1 * m2 * papb * RL * LRLR - 4. * c00 * m1 * m2 * papb * LR * RLRL);

  return real(me0 + me1);
}

double FI::MVutr4u(double p1s, double p2s, double p3s, double ml1s, double ml2s, double ml3s,
                   int ieps) {
  double ml1 = sqrt(ml1s);
  double ml2 = sqrt(ml2s);
  double ml3 = sqrt(ml3s);
  complex<double> c0;
  complex<double> c00;
  complex<double> c1;
  complex<double> c2;
  complex<double> c11;
  complex<double> c12;
  complex<double> c22;
  complex<double> me0;
  complex<double> me1;

  //   Npf *cc = new Npf(0., m2s - 2.*pap2, m1s, m1l, m2l, m3l, mul, ieps);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me0 =
      +ivu2s1 * ivu2s2 *
      (-4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RLRL -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * RR * RRRR -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LRLR -
       4. * c0 * pap2 * pbp1 * ml1 * ml3 * LL * LLLL -
       4. * c0 * m1 * pap2 * pbp1 * ml1 * RR * LRRR - 4. * c0 * m1 * pap2 * pbp1 * ml1 * RR * LLRL -
       4. * c0 * m1 * pap2 * pbp1 * ml1 * LL * RRLR - 4. * c0 * m1 * pap2 * pbp1 * ml1 * LL * RLLL +
       8. * c2 * pap2 * pow(pbp1, 2) * RL * LRLR + 8. * c2 * pap2 * pow(pbp1, 2) * RL * LLLL +
       8. * c2 * pap2 * pow(pbp1, 2) * LR * RLRL + 8. * c2 * pap2 * pow(pbp1, 2) * LR * RRRR -
       4. * c2 * m1 * pap2 * pbp1 * ml3 * RL * RRLR - 4. * c2 * m1 * pap2 * pbp1 * ml3 * RL * RLLL -
       4. * c2 * m1 * pap2 * pbp1 * ml3 * LR * LRRR - 4. * c2 * m1 * pap2 * pbp1 * ml3 * LR * LLRL -
       4. * c2 * m1 * pap2 * pbp1 * ml1 * RR * LRRR - 4. * c2 * m1 * pap2 * pbp1 * ml1 * RR * LLRL -
       4. * c2 * m1 * pap2 * pbp1 * ml1 * LL * RRLR - 4. * c2 * m1 * pap2 * pbp1 * ml1 * LL * RLLL -
       4. * c2 * m1s * pap2 * pbp1 * RL * LRLR - 4. * c2 * m1s * pap2 * pbp1 * RL * LLLL -
       4. * c2 * m1s * pap2 * pbp1 * LR * RLRL - 4. * c2 * m1s * pap2 * pbp1 * LR * RRRR -
       4. * c22 * pow(m1, 2) * pap2 * pbp1 * RL * LRLR - 4. * c22 * m1s * pap2 * pbp1 * RL * LLLL);
  me0 += +ivu2s1 * ivu2s2 *
         (-4. * c22 * m1s * pap2 * pbp1 * LR * RLRL - 4. * c22 * m1s * pap2 * pbp1 * LR * RRRR +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * RL * LRLR +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * RL * LLLL +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * LR * RLRL +
          8. * (c12 + c22) * pap2 * pow(pbp1, 2) * LR * RRRR - 16. * c00 * pap2 * pbp1 * RL * LRLR -
          16. * c00 * pap2 * pbp1 * RL * LLLL - 16. * c00 * pap2 * pbp1 * LR * RLRL -
          16. * c00 * pap2 * pbp1 * LR * RRRR);

  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11, &c12, &c22);

  me1 = +ivu2s1 * ivu2s2 *
        (8. * c00 * pap2 * pbp1 * RL * LRLR + 8. * c00 * pap2 * pbp1 * RL * LLLL +
         8. * c00 * pap2 * pbp1 * LR * RLRL + 8. * c00 * pap2 * pbp1 * LR * RRRR);

  return real(me0 + me1);
}
