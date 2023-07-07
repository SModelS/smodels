#ifndef NPF_H
#define NPF_H
#include "clooptools.h"
#include <iostream>

// sets PV-integrals
static void SetA(double ml1s, int ieps, complex<double> *a0,
                 complex<double> *a00) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);

    *a0 = A0(ml1s);
    *a00 = A00(ml1s);
  } else {
    *a0 = (0.0, 0.0);
    *a00 = (0.0, 0.0);
  }
}

static void SetB(double p1s, double ml1s, double ml2s, int ieps,
                 complex<double> *b0, complex<double> *b1, complex<double> *b00,
                 complex<double> *b11, complex<double> *db0,
                 complex<double> *db1, complex<double> *db00,
                 complex<double> *db11) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);

    *b0 = B0(p1s, ml1s, ml2s);
    *b1 = B1(p1s, ml1s, ml2s);
    *b00 = B00(p1s, ml1s, ml2s);
    *b11 = B11(p1s, ml1s, ml2s);
    *db0 = DB0(p1s, ml1s, ml2s);
    *db1 = DB1(p1s, ml1s, ml2s);
    *db00 = DB00(p1s, ml1s, ml2s);
    *db11 = DB11(p1s, ml1s, ml2s);
  } else {
    *b0 = (0.0, 0.0);
    *b1 = (0.0, 0.0);
    *b00 = (0.0, 0.0);
    *b11 = (0.0, 0.0);
    *db0 = (0.0, 0.0);
    *db1 = (0.0, 0.0);
    *db00 = (0.0, 0.0);
    *db11 = (0.0, 0.0);
  }
}

static void SetC(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps, complex<double> *c0,
                 complex<double> *c00, complex<double> *c1, complex<double> *c2,
                 complex<double> *c11, complex<double> *c12,
                 complex<double> *c22) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);
    *c0 = C0i(cc0, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c00 = C0i(cc00, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c1 = C0i(cc1, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c2 = C0i(cc2, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c11 = C0i(cc11, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c12 = C0i(cc12, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c22 = C0i(cc22, p1s, p2s, p3s, ml1s, ml2s, ml3s);

  } else {
    *c0 = (0.0, 0.0);
    *c00 = (0.0, 0.0);
    *c1 = (0.0, 0.0);
    *c2 = (0.0, 0.0);
    *c11 = (0.0, 0.0);
    *c12 = (0.0, 0.0);
    *c22 = (0.0, 0.0);
  }
}

static void SetC(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps, complex<double> *c0,
                 complex<double> *c00, complex<double> *c1, complex<double> *c2,
                 complex<double> *c11, complex<double> *c12,
                 complex<double> *c22, complex<double> *c001,
                 complex<double> *c002, complex<double> *c111,
                 complex<double> *c112, complex<double> *c122,
                 complex<double> *c222) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);
    *c0 = C0i(cc0, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c00 = C0i(cc00, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c1 = C0i(cc1, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c2 = C0i(cc2, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c11 = C0i(cc11, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c12 = C0i(cc12, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c22 = C0i(cc22, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c001 = C0i(cc001, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c002 = C0i(cc002, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c111 = C0i(cc111, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c112 = C0i(cc112, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c122 = C0i(cc122, p1s, p2s, p3s, ml1s, ml2s, ml3s);
    *c222 = C0i(cc222, p1s, p2s, p3s, ml1s, ml2s, ml3s);
  } else {
    *c0 = (0.0, 0.0);
    *c00 = (0.0, 0.0);
    *c1 = (0.0, 0.0);
    *c2 = (0.0, 0.0);
    *c11 = (0.0, 0.0);
    *c12 = (0.0, 0.0);
    *c22 = (0.0, 0.0);
    *c001 = (0.0, 0.0);
    *c002 = (0.0, 0.0);
    *c111 = (0.0, 0.0);
    *c112 = (0.0, 0.0);
    *c122 = (0.0, 0.0);
    *c222 = (0.0, 0.0);
  }
}

static void SetD(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps, complex<double> *d0, complex<double> *d1,
                 complex<double> *d2, complex<double> *d3, complex<double> *d00,
                 complex<double> *d11, complex<double> *d22,
                 complex<double> *d33, complex<double> *d12,
                 complex<double> *d13, complex<double> *d23) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);

    *d0 = D0i(dd0, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d1 = D0i(dd1, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d2 = D0i(dd2, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d3 = D0i(dd3, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);

    *d00 = D0i(dd00, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d11 = D0i(dd11, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d22 = D0i(dd22, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d33 = D0i(dd33, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d12 = D0i(dd12, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d13 = D0i(dd13, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d23 = D0i(dd23, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);

  } else {
    *d0 = (0.0, 0.0);
    *d1 = (0.0, 0.0);
    *d2 = (0.0, 0.0);
    *d3 = (0.0, 0.0);

    *d00 = (0.0, 0.0);
    *d11 = (0.0, 0.0);
    *d22 = (0.0, 0.0);
    *d33 = (0.0, 0.0);
    *d12 = (0.0, 0.0);
    *d13 = (0.0, 0.0);
    *d23 = (0.0, 0.0);
  }
}

// so many coeffs to set; maybe better to put those coeffs to class FI?
static void
SetD(double p1s, double p2s, double p3s, double p4s, double pas, double pbs,
     double ml1s, double ml2s, double ml3s, double ml4s, int ieps,
     complex<double> *d0, complex<double> *d1, complex<double> *d2,
     complex<double> *d3, complex<double> *d00, complex<double> *d11,
     complex<double> *d22, complex<double> *d33, complex<double> *d12,
     complex<double> *d13, complex<double> *d23, complex<double> *d001,
     complex<double> *d002, complex<double> *d003, complex<double> *d111,
     complex<double> *d112, complex<double> *d113, complex<double> *d122,
     complex<double> *d123, complex<double> *d133, complex<double> *d222,
     complex<double> *d223, complex<double> *d233, complex<double> *d333) {

  if (ieps >= 0 && ieps <= 2) {
    setlambda(-ieps);

    *d0 = D0i(dd0, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d1 = D0i(dd1, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d2 = D0i(dd2, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d3 = D0i(dd3, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);

    *d00 = D0i(dd00, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d11 = D0i(dd11, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d22 = D0i(dd22, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d33 = D0i(dd33, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d12 = D0i(dd12, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d13 = D0i(dd13, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d23 = D0i(dd23, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);

    *d001 = D0i(dd001, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d002 = D0i(dd002, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d003 = D0i(dd003, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d111 = D0i(dd111, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d112 = D0i(dd112, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d113 = D0i(dd113, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d122 = D0i(dd122, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d123 = D0i(dd123, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d133 = D0i(dd133, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d223 = D0i(dd223, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d222 = D0i(dd222, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d223 = D0i(dd223, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d233 = D0i(dd233, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);
    *d333 = D0i(dd333, p1s, p2s, p3s, p4s, pas, pbs, ml1s, ml2s, ml3s, ml4s);

  } else {
    *d0 = (0.0, 0.0);
    *d1 = (0.0, 0.0);
    *d2 = (0.0, 0.0);
    *d3 = (0.0, 0.0);

    *d00 = (0.0, 0.0);
    *d11 = (0.0, 0.0);
    *d22 = (0.0, 0.0);
    *d33 = (0.0, 0.0);
    *d12 = (0.0, 0.0);
    *d13 = (0.0, 0.0);
    *d23 = (0.0, 0.0);

    *d001 = (0.0, 0.0);
    *d002 = (0.0, 0.0);
    *d003 = (0.0, 0.0);
    *d111 = (0.0, 0.0);
    *d112 = (0.0, 0.0);
    *d113 = (0.0, 0.0);
    *d122 = (0.0, 0, 0);
    *d123 = (0.0, 0, 0);
    *d133 = (0.0, 0, 0);
    *d223 = (0.0, 0, 0);
    *d222 = (0.0, 0, 0);
    *d223 = (0.0, 0, 0);
    *d233 = (0.0, 0, 0);
    *d333 = (0.0, 0, 0);
  }
}

#endif