// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

#include "dipoles.h"
#include "gsl_all.h"
#include "integration_method.h"
#include "kinematics.h"
#include "params.h"
#include "pdf.h"
#include "pxs.h"
#include "utils.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "hxs.h"

using namespace std;
#define CMLLN 0.9  // C coeff. of the inverse Mellin transform
#define VBSSL 2.9  // V coeff. of the inverse Bessel transform
#define PT2CUT 0.0 // pt2 cut for the integration (non-zero for joint resummation)

// two-particle phase space used for LO and virtual corrections
double dPS2_dlnM2(double &xa, double &xb, double &t, double *x, Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);

  // at LO s is the invariant mass (pa + pb)^2 = (p1 + p2)^2
  const double s = params->mis;

  // minimum and maximum invariant mass squared
  if (s < pow(m1 + m2, 2)) {
    cout << "dPS2_dlnM2: M2 too small\n";
    exit(0);
  } else if (s > sh) {
    cout << "dPS2_dlnM2: M2 too large\n";
    exit(1);
  }

  // Jacobian initialization
  double djac = 389379304.0;

  // Integration variable xa
  double xamin = s / sh;
  double xamax = 1.0;
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
  }

  // Getting xb
  xb = xamin / xa;

  // Integration variable t
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;

  // Phase Space factor
  // ( dsigma / dxb = (dsigma / dlnMs) (dlnMs / dxb) = dsigma / dlnMs dln(xa xb
  // S) / dxb)
  // -> d sigma / dlnMs = d sigma/ dxb * xb (that is the reason for the xb
  // factor in the jacobian)
  djac *= 0.125 * M_1_PI / s * xb;

  return djac;
}

// Phase space for real collinear emission.
double dPS2C_dlnM2(double &djacdelta, double &djacplus, double &xa, double &xb, double &xc,
                   double &tc, double *x, Parameters *params) {
  const double sh = params->sh;

  // initialize final state masses
  double m1;
  double m2;

  // macro to set masses of the final state particles
  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);

  // the invariant mass squared
  const double mis = params->mis;

  if (mis < pow(m1 + m2, 2)) {
    cout << "dPS2d_dlnM2: M2 too small\n";
    exit(0);
  } else if (mis > sh) {
    cout << "dPS2d_dlnM2: M2 too large\n";
    exit(1);
  }

  // Jacobian initialization.
  double djac = 389379304.0;

  // Integration variable xa and xb.
  double xamin = mis / sh;
  double xamax = 1.0;
  double xcmin = xamin;
  double xcmax = 1.0;

  // mapping to interval [0:1]
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xcmin /= xa;      // -> xcmin = M^2/ (sh * xa) = xb * M^2 / s = xb * xc := z
    djacdelta = djac; // no x[1] contribution
    xc = (xcmax - xcmin) * x[1] + xcmin;
    djac *= xcmax - xcmin;
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xcmin /= xa;
    djacdelta = djac;
    xc = xcmin * pow(xcmax / xcmin, x[1]);
    djac *= xc * log(xcmax / xcmin);
  }

  // we do not integrate over xb.
  // xb chosen in such a way that M^2 is fixed.
  xb = xamin / xa / xc;

  // sc = xc s is the squared invariant mass of the two final state particles
  const double sc = mis;

  // Integration variable tc
  const double tcmin = -.5 * (sc - m1s - m2s + kln(sc, m1s, m2s));
  const double tcmax = tcmin + kln(sc, m1s, m2s);
  tc = (tcmax - tcmin) * x[2] + tcmin;
  djacdelta *= tcmax - tcmin;
  djac *= tcmax - tcmin;

  // Phase Space factor
  // similar to dPS2C in hxs.cc we have three different phase space factors
  djacdelta *= 0.125 * M_1_PI / sc * xb * xc;
  djac *= 0.125 * M_1_PI / sc * xb; // no xc contribution, see above
  djacplus = djac * xc;

  return djac;
}

// Three-particle phase space.
double dPS3_dlnM2(double &xa, double &xb, double &M2, double &pt2, double &Tp, double &Pp,
                  double *x, Parameters *params) {
  const double sh = params->sh;

  // massive final state masses
  double m1;
  double m2;

  // macro
  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);

  // invariant mass squared
  M2 = params->mis;

  // check minimal and maximal M^2 value
  if (M2 < pow(m1 + m2, 2)) {
    cout << "dPS3_dlnM2: M2 too small\n";
    exit(0);
  } else if (M2 > sh) {
    cout << "dPS3_dlnM2: M2 too large\n";
    exit(0);
  }

  // Jacobian initialization.
  // Conversion factor from GeV^-2 to pb.
  double djac = 389379304.0;

  // Integration variable xa, xb.
  double xamin = M2 / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  // Mapping to interval [0:1].
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac *= xbmax - xbmin;
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac *= xb * log(xbmax / xbmin);
  }

  const double s = xa * xb * sh;

  // Integration variable pt2.
  const double pt2min = PT2CUT;
  const double pt2max = .25 * pow(s - M2, 2) / s;
  if (pt2max < pt2min) {
    return 0.0;
  }
  pt2 = (pt2max - pt2min) * x[2] + pt2min;
  djac *= pt2max - pt2min;

  // Integration variable Tp.
  const double Tpmin = 0.0;
  const double Tpmax = M_PI;
  Tp = (Tpmax - Tpmin) * x[3] + Tpmin;
  djac *= Tpmax - Tpmin;

  // Integration variable Pp.
  const double Ppmin = 0.0;
  const double Ppmax = M_PI;
  djac *= 2.0;
  Pp = (Ppmax - Ppmin) * x[4] + Ppmin;
  djac *= Ppmax - Ppmin;

  // Phase Space factor.
  djac *= kln(M2, m1s, m2s) * sin(Tp) / sqrt(pow(s - M2, 2) - 4.0 * s * pt2) / pow(4.0 * M_PI, 4);

  return djac;
}

// Three-particle phase space used for on-shell subtraction for the
// associated production of gauginos and gluinos.
// Here s2 = s23 can be on-shell. s1 is the usual M^2.
// (Reference: Particle kinematics - Byckling, Kajantie)
// (Helicity frame.)
double dPS3_ONSHELL23_dlnM2(double &xa, double &xb, double &s2, double &t1, double &s1, double &phi,
                            double *x, Parameters *params) {

  // Hadronic com energy.
  const double sh = params->sh;

  // Final state masses.
  double m1;
  double m2;

  // Macro which sets m1 and m2 (see utils.h).
  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  // Invariant mass.
  double M2;
  M2 = params->mis;
  s1 = M2; // s1 is the invariant mass squared.

  // Check minimal and maximal M^2 value.
  if (M2 < pow(m1 + m2, 2)) {
    cout << "dPS3_dlnM2: M2 too small\n";
    exit(0);
  } else if (M2 > sh) {
    cout << "dPS3_dlnM2: M2 too large\n";
    exit(0);
  }

  // Initialize jacobian.
  double djac = 1.0;
  // Conversion factor from GeV^-2 to pb.
  djac = djac * 389379304.0;

  // Integration variable xa, xb and s1
  double xamin = M2 / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  // Mapping to interval [0:1].
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac = djac * (xbmax - xbmin);
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac = djac * xb * log(xbmax / xbmin);
  }

  // Partonic com energy.
  const double s = xa * xb * sh;

  // s2 = s23 = (p2 + p3)^2.
  double s2min = m2s;
  double s2max = pow2(sqrt(s) - m1);

  // Check if phase space point is physical.
  if (s2max < s2min) {
    djac = 0.0;
    return djac;
  }
  // Mapping to [0:1] if Breit Mapping is deactivated.
  // if (params->deg_squarks == true) {
  //   double z4p = atan((s2max - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double z4m = atan((s2min - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double y = (z4p - z4m) * x[2] + z4m;
  //   s2 = params->mSQs[0] + params->mSQs[0] * 1.0E-2 * tan(y);
  //   djac = djac * params->mSQs[0] * 1.0E-2 * (z4p - z4m)
  //     * 1.0/pow2(cos(z4m - x[2] * z4m + x[2] * z4p));
  // } else {
  if (s2max == 0.0 || s2min == 0.0) { // actually can never happen :D
    s2 = (s2max - s2min) * x[2] + s2min;
    djac = djac * (s2max - s2min);
  } else {
    // log mapping
    s2 = s2min * pow(s2max / s2min, x[2]);
    djac *= s2 * log(s2max / s2min);
  }
  //}

  // Integration over t1 = (pa - p2)^2
  double t1min;
  double t1max;
  t1min = -.5 * (s - m1s - s2 + kln(s, s2, m1s));
  t1max = t1min + kln(s, s2, m1s);
  t1 = (t1max - t1min) * x[3] + t1min;
  djac = djac * (t1max - t1min);

  // calculating limits for s1 to check if phase space point is physical.
  double s1min = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m1s) * (s2 + m2s) -
                 1.0 / (2.0 * s2) * (kln(s, s2, m1s) * kln(s2, m2s, 0));
  double s1max = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m1s) * (s2 + m2s) +
                 1.0 / (2.0 * s2) * (kln(s, s2, m1s) * kln(s2, m2s, 0));
  if (s1min > M2 || s1max < M2 || Rkln(s, s2, m1s) < 0.0 || Rkln(s2, m2s, 0) < 0.0) {
    djac = 0.0;
    return djac;
  }

  // Integration over azimuthal angle.
  double phimin = 0.0;
  double phimax = 2.0 * M_PI;
  phi = (phimax - phimin) * x[4] + phimin;
  djac = djac * (phimax - phimin);

  // Normalization factor.
  djac = djac / pow((2.0 * M_PI), 5);
  // kln(s,0,0) for the flux factor
  double sqrt_lambda = s;

  // Final jacobian. (note: kln is sqrt(kaellen))
  djac = djac * M_PI / (2.0 * sqrt_lambda) / (4.0 * kln(s, s2, m1s));

  // check if phase space point is physical.
  if (Rkln(s, s2, m1s) < 0.0 || Rkln(s2, m2s, 0) < 0.0) {
    djac = 0.0;
  }

  return djac;
}

// Three-particle phase space used for on-shell subtraction for the
// associated production of gauginos and gluinos.
// Here s2 = s13 can be on-shell. s1 = s12 is the usual M2.
// same as dPS3_ONSHELL23_dlnM2, but p1 <-> p2;
double dPS3_ONSHELL13_dlnM2(double &xa, double &xb, double &s2, double &t1, double &s1, double &phi,
                            double *x, Parameters *params) {

  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow2(m1);
  const double m2s = pow2(m2);

  double M2;
  M2 = params->mis;
  s1 = M2;

  if (M2 < pow(m1 + m2, 2)) {
    cout << "dPS3_dlnM2: M2 too small\n";
    exit(0);
  } else if (M2 > sh) {
    cout << "dPS3_dlnM2: M2 too large\n";
    exit(0);
  }

  double djac = 1.0;
  // Conversion factor from GeV^-2 to pb.
  djac = djac * 389379304.0;

  // Integration variable xa, xb and s1
  double xamin = M2 / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;

  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac = djac * (xbmax - xbmin);
  } else {
    xa = xamin * pow(xamax / xamin, x[0]);
    djac *= xa * log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * pow(xbmax / xbmin, x[1]);
    djac = djac * xb * log(xbmax / xbmin);
  }

  const double s = xa * xb * sh;

  // integration over s2 and t1
  double s2min = m1s;
  double s2max = pow2(sqrt(s) - m2);
  if (s2max < s2min) {
    djac = 0.0;
    return djac;
  }

  // different mapping of s2
  // if (params->deg_squarks == true) {
  //   double z4p = atan((s2max - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double z4m = atan((s2min - params->mSQs[0])/(params->mSQs[0] * 1.0E-2));
  //   double y = (z4p - z4m) * x[2] + z4m;
  //   s2 = params->mSQs[0] + params->mSQs[0] * 1.0E-2 * tan(y);
  //   djac = djac * params->mSQs[0] * 1.0E-2 * (z4p - z4m)
  //     * 1.0/pow2(cos(z4m - x[2] * z4m + x[2] * z4p));
  // } else {
  if (s2max == 0.0 || s2min == 0.0) { // actually can never happen :D
    s2 = (s2max - s2min) * x[2] + s2min;
    djac = djac * (s2max - s2min);
  } else {
    // log mapping
    s2 = s2min * pow(s2max / s2min, x[2]);
    djac *= s2 * log(s2max / s2min);
  }
  //}

  // Integration over t1.
  double t1min;
  double t1max;
  t1min = -.5 * (s - m2s - s2 + kln(s, s2, m2s));
  t1max = t1min + kln(s, s2, m2s);
  t1 = (t1max - t1min) * x[3] + t1min;
  djac = djac * (t1max - t1min);

  // Check if phase space point is physical.
  double s1min = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m2s) * (s2 + m1s) -
                 1.0 / (2.0 * s2) * (kln(s, s2, m2s) * kln(s2, m1s, 0));
  double s1max = m1s + m2s + 1.0 / (2.0 * s2) * (s - s2 - m2s) * (s2 + m1s) +
                 1.0 / (2.0 * s2) * (kln(s, s2, m2s) * kln(s2, m1s, 0));

  if (s1min > M2 || s1max < M2 || Rkln(s, s2, m2s) < 0.0 || Rkln(s2, m1s, 0) < 0.0) {
    djac = 0.0;

    return djac;
  }

  // Integration over azimuthal angle.
  double phimin = 0.0;
  double phimax = 2.0 * M_PI;

  phi = (phimax - phimin) * x[4] + phimin;
  djac = djac * (phimax - phimin);

  // Normalization factor.
  djac = djac / pow((2.0 * M_PI), 5);

  // For the flux factor.
  double sqrt_lambda = s;

  // Final jacobian.
  djac = djac * M_PI / (2.0 * sqrt_lambda) / (4.0 * kln(s, s2, m2s));

  // Check if physical.
  if (Rkln(s, s2, m2s) < 0.0 || Rkln(s2, m1s, 0) < 0.0) {
    djac = 0.0;
  }

  return djac;
}

// Setting up the different integrands.

// Born
double IB_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // creating an instance of the parameter class

  // setting the needed variables xa,xb,t and s for a specific integration point
  // x (array)
  double xa, xb, t;
  const double djac = dPS2_dlnM2(xa, xb, t, x, params);
  const double s = xa * xb * params->sh;

  // Symmetry factor for two neutralinos in the final state.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  // initialize the cross section.
  double sig = 0.0;

  // Sum over all possible initial states.
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state.
      params->in1 = i0;
      params->in2 = i1;

      // Check if charge is conserved.
      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 born_gagl(s, t, params);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 born_sleptons(s, t, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 born_gauginos(s, t, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 born_leptons(s, t, params);
        } else if (is_squark_gaugino(params->out1, params->out2)) {
          if (params->out1 > 30 && i1 == 0) { // Initial quark+gluon
            params->in2 = -1;                 // gluon, fails if used as (anti-)quark
            sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) * born_gasq_apn(s, t, params);
          } else if (params->out2 > 30 && i0 == 0) { // Initial quark+gluon
            params->in1 = -1;                        // gluon, fails if used as (anti-)quark
            sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) * born_gasq_apn(s, t, params);
          }
        }
      }
    }
  }
  if (is_squark_gaugino(params->out1, params->out2)) {
    return 4.0 * M_PI * aS(params->murs, params->set) * sig * djac * 1.0 / (2.0 * s * 96.0);
  } else if (is_gaugino_gluino(params->out1,
                               params->out2)) { // LO with strong coupling.
    double g3s = std::norm(params->gqq[0][0].R);
    return 4.0 * M_PI * aS(params->murs, params->set) / g3s * 1.0 / (24.0 * s * 3.0) * sig * djac;
  } else { // LO without any strong coupling.
    // Flux, symmetry factor, spin and color average. fc = 3; average 1/36;
    // and 1/2s flux -> 1/24 as overall factor
    return dij * sig * djac / (24.0 * s);
  }
}

// Virtual corrections. (completely analog to IB_dlnM2, but different partonic
// cross section)
double IV_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb, t;
  const double djac = dPS2_dlnM2(xa, xb, t, x, params);
  const double s = xa * xb * params->sh;

  // Symmetry factor.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = 0.5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  // initialize cross section.
  double sig = 0.0;

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Check if charge is conserved.
      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_gaugino_gluino(s, t, params) + DipI_gagl(s, t, params));
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_sleptons(s, t, params) + DipI_sleptons(s, t, params));
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_gauginos(s, t, params) + DipI_gauginos(s, t, params));
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (Virt_leptons(s, t, params) + DipI_leptons(s, t, params));
        } else if (is_squark_gaugino(params->out1, params->out2)) {
          if (params->out1 > 30 && i1 == 0) { // Initial quark+gluon
            params->in2 = -1;                 // gluon, fails if used as (anti-)quark
            sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                   (Virt_gaugino_squark(s, t, params) + DipI_gasq(s, t, params));
          } else if (params->out2 > 30 && i0 == 0) { // Initial quark+gluon
            params->in1 = -1;                        // gluon, fails if used as (anti-)quark
            sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                   (Virt_gaugino_squark(s, t, params) + DipI_gasq(s, t, params));
          }
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average.
  if (is_squark_gaugino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    // factor 96 already in dip and virt
    return sig * djac / s / 2. * (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) / g3s / g3s; // a_s a_s
  } else if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) * dij * sig * djac * 1.0 / (2.0 * s) *
           (1.0 / 36.0);
  } else {
    return aS(params->murs, params->set) / (2.0 * M_PI) * dij * sig * djac * 4.0 / (24.0 * s * 3.0);
  }
}

// Collinear emission.
double IC_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double djacdelta, djacplus, xa, xb, xc, tc;
  double djac = dPS2C_dlnM2(djacdelta, djacplus, xa, xb, xc, tc, x, params);
  // sc = xc s = xa xb xc S = M^2.
  const double sc = params->mis;
  // x value for parton out of proton b
  // leading to the same final state invariant mass without any splitting
  // (used for plus and delta distribution part)
  const double xd = xb * xc;

  // xcmin; needed for the correct use of the plus distribution.
  double z = xb * xc;

  // partonic com energy before the splitting.
  const double s = sc / xc;
  double sa1c, sb1c, sa1d, sb1d;
  double sa1 = 0.0;
  double sb1 = 0.0;

  if (is_gaugino_gluino(params->out1, params->out2)) {
    // 2 * pap1, where p1 belongs to the boosted phase space
    sa1c = sja_gagl(sc, tc, params) / xc;
    // 2 * pbp1, where p1 belongs to the boosted phase space
    sb1c = sjb_gagl(sc, tc, params) / xc;

    // Needed for the plus and delta distribution part.
    // 2 * pap1, where pa and p1 belong to the boosted phase space.
    sa1d = sja_gagl(sc, tc, params);
    // 2 * pbp1, where pb and p1 belong to the boosted phase space.
    sb1d = sjb_gagl(sc, tc, params);

  } else if (is_squark_gaugino(params->out1, params->out2)) {
    // 2 * pap1, where p1 belongs to the boosted phase space
    sa1c = sja_gasq(sc, tc, params) / xc;
    // 2 * pbp1, where p1 belongs to the boosted phase space
    sb1c = sjb_gasq(sc, tc, params) / xc;

    // Needed for the plus and delta distribution part.
    // 2 * pap1, where pa and p1 belong to the boosted phase space.
    sa1d = sja_gasq(sc, tc, params);
    // 2 * pbp1, where pb and p1 belong to the boosted phase space.
    sb1d = sjb_gasq(sc, tc, params);
  }
  // Majorana symmetry factor.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);
  double gd, qd[2][6]; // Needed for plus and delta distribution part.
  pdfX(gd, qd, xd, params->mufs);

  // Initialize collinear born cross section and the total sigma.
  double bornc = 0.0;
  double born = 0.0;
  double sig = 0.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Test charge conservation
      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          bornc = born_gagl(sc, tc, params);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          bornc = born_sleptons(sc, tc, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          bornc = born_gauginos(sc, tc, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          bornc = born_leptons(sc, tc, params);
        } else if (is_squark_gaugino(params->out1, params->out2)) {
          bornc = born_gasq(sc, tc, params);
        }
        // Initial state quark and gluon:
        // for associated squark electroweakino pair production
        // proton-proton -> ic = 0; proton-antiproton -> ic = 1;
        if (is_squark_gaugino(params->out1, params->out2)) {
          int pdf, pdf_ic;
          pdf = ((params->out1 > 30 ? 0 : 1)) % 2;
          pdf_ic = (params->ic + (params->out1 > 30 ? 0 : 1)) % 2;
          if (params->out1 > 30) {

            //*/
            if (i1 == 0) {
              params->in2 = -1;
#ifdef SQGA_QG
              //// INSERTION OPERATOR P

              // Diagonal splitting

              // P_bold_app_b_j for the diagram with particle a = q radiating
              // a gluon, particle b = g (use sja)
              sig +=
                  (qa[0][i0] * gb) *
                  P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out1 - 31],
                                        sa1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qa[0][i0] * gd) *
                  (-P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out1 - 31],
                                          sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out1 - 31],
                                         sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);
              //*/

              sig +=
                  (qb[params->ic][i0] * ga) *
                  P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out1 - 31],
                                        sb1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              ///*
              sig +=
                  (qd[params->ic][i0] * ga) *
                  (-P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out1 - 31],
                                          sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out1 - 31],
                                         sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);
              //*/

              // P_bold_app_b_j for the diagram with particle a = q, particle
              // b = g radiating a gluon (use sjb)
              sig +=
                  (qa[0][i0] * gb) *
                  P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out1 - 31],
                                        sb1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qa[0][i0] * gd) *
                  (-P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out1 - 31],
                                          sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out1 - 31],
                                         sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);
              //*/
              sig +=
                  (qb[params->ic][i0] * ga) *
                  P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out1 - 31],
                                        sa1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qd[params->ic][i0] * ga) *
                  (-P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out1 - 31],
                                          sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out1 - 31],
                                         sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc

                  );
              //*/
#endif
              // Off-diagonal splitting
#ifdef SQGA_GG
              // P_bold_app_b_j for the diagram with particle a = g radiating
              // a quark, particle b = g (use sja)
              sig += (ga * gb) *
                     // P_gq (no delta- or plus-distribution part)
                     P_plus_K_bold(PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK, xc, z,
                                   params->mSQ[params->out1 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                   INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                     xc;

              // P_bold_app_b_j for the diagram with particle a = g, particle
              // b = g radiating a quark (use sjb)
              sig += (ga * gb) *
                     P_plus_K_bold(PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK, xc, z,
                                   params->mSQ[params->out1 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                   INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                     xc;
#endif
            }
            ///*
            // Off-diagonal splitting

            // (K+P)_qg
            // for the diagram with particle a = q radiating
            // a quark, particle b = q (use sja)

#ifdef SQGA_UU
            double f = 1.;
            f = 0.5;

            sig += f * (qa[pdf][i1] * qb[pdf_ic][i0] + qb[pdf_ic][i1] * qa[pdf][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out1 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;

            // (K+P)_qg
            // for the diagram with particle a = q,
            // , particle b = q radiating a quark (use sjb)
            sig += f * (qb[pdf_ic][i1] * qa[pdf][i0] + qa[pdf][i1] * qb[pdf_ic][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out1 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
#endif
#ifdef SQGA_UUB

            // (P+K)_qbg
            // for the diagram with particle a = qb,
            // , particle b = q radiating a antiquark (use sjb)
            sig += (qa[1][i1] * qb[params->ic][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out1 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
            sig += (qb[1 - params->ic][i1] * qa[0][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out1 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
#endif
            // this case we have a squark as particle 1, meaning that we
            // have a quark in the initial state
          } else if (params->out2 > 30) {
            if (i0 == 0) {
              params->in1 = -1;

              //// INSERTION OPERATOR P

              // Diagonal splitting

#ifdef SQGA_QG
              sig +=
                  (qa[1][i0] * gb) *
                  P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out2 - 31],
                                        sa1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qa[1][i0] * gd) *
                  (-P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out2 - 31],
                                          sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out2 - 31],
                                         sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);

              sig +=
                  (qb[1 - params->ic][i0] * ga) *
                  P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out2 - 31],
                                        sb1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qd[1 - params->ic][i0] * ga) *
                  (-P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out2 - 31],
                                          sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out2 - 31],
                                         sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);

              // P_bold_app_b_j for the diagram with particle a = q, particle
              // b = g radiating a gluon (use sjb)
              sig +=
                  (qa[1][i0] * gb) *
                  P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out2 - 31],
                                        sb1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qa[1][i0] * gd) *
                  (-P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out2 - 31],
                                          sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out2 - 31],
                                         sb1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc);

              sig +=
                  (qb[1 - params->ic][i0] * ga) *
                  P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                        FUNCTION_ALL_WDELTA, xc, z, params->mSQ[params->out2 - 31],
                                        sa1c, s, params->mufs, 1.0, INITIAL_AND_FINAL) *
                  djac * bornc;
              sig +=
                  (qd[1 - params->ic][i0] * ga) *
                  (-P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                          FUNCTION_PLUS, xc, z, params->mSQ[params->out2 - 31],
                                          sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacplus * bornc +
                   P_plus_K_bold_aap_b_j(PARTON_GLUON, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK,
                                         FUNCTION_DELTA, xc, z, params->mSQ[params->out2 - 31],
                                         sa1d, sc, params->mufs, 1.0, INITIAL_AND_FINAL) *
                       djacdelta * bornc

                  );
#endif
              // Off-diagonal splitting

#ifdef SQGA_GG
              // P_bold_app_b_j for the diagram with particle a = g radiating
              // a quark, particle b = g (use sja)
              sig += (ga * gb) *
                     // P_gq (no delta- or plus-distribution part)
                     P_plus_K_bold(PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK, xc, z,
                                   params->mSQ[params->out2 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                   INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                     xc;

              // P_bold_app_b_j for the diagram with particle a = g, particle
              // b = g radiating a quark (use sjb)
              sig += (ga * gb) *
                     P_plus_K_bold(PARTON_GLUON, PARTON_QUARK, PARTON_GLUON, PARTON_SQUARK, xc, z,
                                   params->mSQ[params->out2 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                   INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                     xc;
#endif
            }
            ///*
            // this case we have a squark as particle 2, meaning that we have
            // a antiquark in the initial state

            // Off-diagonal splitting
#ifdef SQGA_UU
            double f = 1.;
            f = 0.5;

            // (K+P)_qg
            // for the diagram with particle a = q radiating
            // a quark, particle b = q (use sja)
            sig += f * (qa[pdf][i1] * qb[pdf_ic][i0] + qb[pdf_ic][i1] * qa[pdf][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out2 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;

            // (K+P)_qg
            // for the diagram with particle a = q,
            // , particle b = q radiating a quark (use sjb)
            sig += f * (qb[pdf_ic][i1] * qa[pdf][i0] + qa[pdf][i1] * qb[pdf_ic][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out2 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
#endif
#ifdef SQGA_UUB

            // (P+K)_qbg
            // for the diagram with particle a = qb,
            // , particle b = q radiating a antiquark (use sjb)
            sig += (qb[params->ic][i1] * qa[1][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out2 - 31], sb1, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
            sig += (qa[0][i1] * qb[1 - params->ic][i0]) *
                   P_plus_K_bold(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, PARTON_SQUARK, xc, z,
                                 params->mSQ[params->out2 - 31], sa1, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL, born, bornc, djac, djacplus, djacdelta) *
                   xc;
#endif
            //*/
          }
        } else if (is_gaugino_gluino(params->out1, params->out2)) {

          // P_qq and P_qbqb with sa1c
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          // Parton with momentum fraction xd = xb*xc leads without emission to
          // same invariant mass.
          // s-> sc and sa1c -> sa1d.
          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 (

                     -P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                     FUNCTION_PLUS, xc, z, sa1d, sc, params->mufs, 1.0,
                                     INITIAL_AND_FINAL) *
                         djacplus * bornc +
                     P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                    FUNCTION_DELTA, xc, z, sa1d, sc, params->mufs, 1.0,
                                    INITIAL_AND_FINAL) *
                         djacdelta * bornc);
          // P_qq and P_qbqb with sb1c.
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 (-P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                  FUNCTION_PLUS, xc, z, sb1d, sc, params->mufs, 1.0,
                                  INITIAL_AND_FINAL) *
                      djacplus * bornc +
                  P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_DELTA, xc, z, sb1d, sc, params->mufs, 1.0,
                                 INITIAL_AND_FINAL) *
                      djacdelta * bornc);

          // K_qq and K_qbqb with sa1c.
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, params->mGL, sa1c, s, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 (-K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                  FUNCTION_PLUS, xc, z, params->mGL, sa1d, sc, 1.0,
                                  INITIAL_AND_FINAL) *
                      djacplus * bornc +
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_DELTA, xc, z, params->mGL, sa1d, sc, 1.0,
                                 INITIAL_AND_FINAL) *
                      djacdelta * bornc);

          // K_qq and K_qbqb with sb1c
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, params->mGL, sb1c, s, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 (-K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                  FUNCTION_PLUS, xc, z, params->mGL, sb1d, sc, 1.0,
                                  INITIAL_AND_FINAL) *
                      djacplus * bornc +
                  K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_DELTA, xc, z, params->mGL, sb1d, sc, 1.0,
                                 INITIAL_AND_FINAL) *
                      djacdelta * bornc);

          // Pgq and Pgqb (no plus or delta distribution part).
          // sa1c
          sig += (qb[params->ic][i0] * ga + qb[1 - params->ic][i1] * ga) *
                 (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, sa1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          // sb1c
          sig += (qa[0][i0] * gb + qa[1][i1] * gb) *
                 (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, sb1c, s, params->mufs, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          // Kgq and Kgqb
          // sa1c
          sig += (+qb[params->ic][i0] * ga + qb[1 - params->ic][i1] * ga) *
                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, params->mGL, sa1c, s, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

          // sb1c
          sig += (qa[0][i0] * gb + qa[1][i1] * gb) *
                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUINO,
                                 FUNCTION_ALL_WDELTA, xc, z, params->mGL, sb1c, s, 1.0,
                                 INITIAL_AND_FINAL) *
                  djac * bornc);

        }

        else { // cases without final state partons in LO

          // Pqq and Pqbqb
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1] +
                  qb[0][i0] * qa[1 - params->ic][i1] + qa[params->ic][i0] * qb[1][i1]) *
                 (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                 FUNCTION_ALL_WDELTA, xc, z, 0.0, s, params->mufs, 1.0, INITIAL) *
                  djac * bornc);

          // plus distribution part.
          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1] +
                  qd[0][i0] * qa[1 - params->ic][i1] + qa[params->ic][i0] * qd[1][i1]) *
                 (-P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                  FUNCTION_PLUS, xc, z, 0.0, sc, params->mufs, 1.0, INITIAL) *
                  djacplus * bornc);

          // delta distribution part.
          sig += 2.0 * (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 (P_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                 FUNCTION_DELTA, xc, z, 0.0, sc, params->mufs, 1.0, INITIAL) *
                  djacdelta * bornc);

          // Kqq and Kqbqb
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1] +
                  qb[0][i0] * qa[1 - params->ic][i1] + qa[params->ic][i0] * qb[1][i1]) *
                 (K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                 FUNCTION_ALL_WDELTA, xc, z, 0.0, 0.0, s, 1.0, INITIAL) *
                  djac * bornc);

          // plus distribution part.
          sig += (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1] +
                  qd[0][i0] * qa[1 - params->ic][i1] + qa[params->ic][i0] * qd[1][i1]) *
                 (-K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                  FUNCTION_PLUS, xc, z, 0.0, 0.0, sc, 1.0, INITIAL) *
                  djacplus * bornc);

          // delta distribution part.
          sig += 2.0 * (qa[0][i0] * qd[1 - params->ic][i1] + qd[params->ic][i0] * qa[1][i1]) *
                 K_bold_aap_b_j(PARTON_QUARK, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                FUNCTION_DELTA, xc, z, 0.0, 0.0, sc, 1.0, INITIAL) *
                 djacdelta * bornc;

          // Pgq and Pgqb
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga + qa[1][i1] * gb +
                  qb[1 - params->ic][i1] * ga) *
                 (P_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                 FUNCTION_ALL_WDELTA, xc, z, 0.0, s, params->mufs, 1.0, INITIAL) *
                  djac * bornc);

          // //Kgq and Kgqb
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga + qa[1][i1] * gb +
                  qb[1 - params->ic][i1] * ga) *

                 (K_bold_aap_b_j(PARTON_GLUON, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_NONE,
                                 FUNCTION_ALL_WDELTA, xc, z, 0.0, 0.0, s, 1.0, INITIAL) *
                  djac * bornc);
        }
      }
    }
  }
  // Flux, symmetry factor, spin and color average
  if (is_squark_gaugino(params->out1, params->out2)) {
    return 4.0 * M_PI * aS(params->murs, params->set) * sig * aS(params->murs, params->set) * 1.0 /
           (2.0 * s * 96.0);
  } else if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return 1.0 / (36.0 * 2.0 * sc) * sig * aS(params->murs, params->set) *
           (4.0 * M_PI * aS(params->murs, params->set)) / g3s;
  } else {
    return dij * 3.0 / (36.0 * 2.0 * sc) * sig * aS(params->murs, params->set);
  }
}

// Integrand for real gluon emission.
double IG_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;
  if (djac == 0.0) {
    return 0.0;
  }

  // Numerical cutoff for stability.
  // (Some seeds could directly probe the divergent region and then you get a
  // nan - nan which is nan and not zero.)
  // Checked that for small pt values the dipoles completely reproduce the
  // 2->3 real emission diagrams.
  if (pt2 < 1.0E-6) {
    return 0.0;
  }

  // Symmetry factor.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;
  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Check if charge is conserved.
      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        // gluon emission for the associated gaugino gluino production.
        if (is_gaugino_gluino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 // Real gluon emission.
                 (real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 1, params)

                  // Corresponding Catani-Seymour Dipoles dsigma^A.
                  // Initial state emitter (quark) and initial state spectator
                  // (antiquark)
                  - Dip_GLGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s, mi2,
                             pt2, th, ph, 0, params) -
                  Dip_GLGA(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s, mi2,
                           pt2, th, ph, 1, params)

                  // Initial state emitter (antiquark) and initial state
                  // spectator (quark)
                  - Dip_GLGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON, s, mi2,
                             pt2, th, ph, 0, params) -
                  Dip_GLGA(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON, s, mi2,
                           pt2, th, ph, 1, params)

                  // Initial state emitter (quark) and final state spectator
                  // (gluino)
                  - Dip_GLGA(INITIAL_FINAL, PARTON_QUARK, PARTON_GLUINO, PARTON_GLUON, s, mi2, pt2,
                             th, ph, 0, params) -
                  Dip_GLGA(INITIAL_FINAL, PARTON_QUARK, PARTON_GLUINO, PARTON_GLUON, s, mi2, pt2,
                           th, ph, 1, params)

                  // Initial state emitter (antiquark) and final state spectator
                  // (gluino)
                  - Dip_GLGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_GLUINO, PARTON_GLUON, s, mi2,
                             pt2, th, ph, 0, params) -
                  Dip_GLGA(INITIAL_FINAL, PARTON_ANTIQUARK, PARTON_GLUINO, PARTON_GLUON, s, mi2,
                           pt2, th, ph, 1, params)

                  // Final state emitter (gluino) and initial state spectator
                  // (antiquark)
                  - Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_ANTIQUARK, PARTON_GLUON, s, mi2,
                             pt2, th, ph, 0, params) -
                  Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_ANTIQUARK, PARTON_GLUON, s, mi2,
                           pt2, th, ph, 1, params)

                  // Final state emitter (gluino) and initial state spectator
                  // (quark)
                  - Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2,
                             th, ph, 0, params) -
                  Dip_GLGA(FINAL_INITIAL, PARTON_GLUINO, PARTON_QUARK, PARTON_GLUON, s, mi2, pt2,
                           th, ph, 1, params));

          // Real gluon emission for slepton pair production.
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (real_gluon_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_sleptons(s, mi2, pt2, th, ph, 1, params)

                  // quark antiquark
                  - (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s,
                                  mi2, pt2, th, ph, 0, params) +
                     Dip_SLEPTONS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s,
                                  mi2, pt2, th, ph, 1, params)

                     // antiquark quark
                     + Dip_SLEPTONS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON,
                                    s, mi2, pt2, th, ph, 0, params) +
                     Dip_SLEPTONS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON, s,
                                  mi2, pt2, th, ph, 1, params)));
          // real gluon emission for gaugino pair production.
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (real_gluon_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_gauginos(s, mi2, pt2, th, ph, 1, params) -

                  (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s,
                                mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_QUARK, PARTON_ANTIQUARK, PARTON_GLUON, s,
                                mi2, pt2, th, ph, 1, params)

                   + Dip_GAUGINOS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON, s,
                                  mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_ANTIQUARK, PARTON_QUARK, PARTON_GLUON, s,
                                mi2, pt2, th, ph, 1, params)));
          // real gluon emission for lepton pair production (still old dipole
          // implementation)
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
                 (real_gluon_leptons(s, mi2, pt2, th, ph, 0, params) +
                  real_gluon_leptons(s, mi2, pt2, th, ph, 1, params) -
                  DipGA_leptons(s, mi2, pt2, th, ph, 0, params) -
                  DipGA_leptons(s, mi2, pt2, th, ph, 1, params) -
                  DipGB_leptons(s, mi2, pt2, th, ph, 0, params) -
                  DipGB_leptons(s, mi2, pt2, th, ph, 1, params));
        }
      }
    }
  }
  // real gluon emission for the associated squark gaugino production
#ifdef SQGA_QG
  if (is_squark_gaugino(params->out1, params->out2)) {
    if (params->out1 > 30) {
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;
        // if (is_charge_conserved(params->in1, params->in2,
        // params->out1,params->out2))
        {

          // Real gluon emission for gaugino-squark production
          // PDF not SC like without *(xa*xb)
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (gausq_qg_gluon_minus_dip(s, mi2, pt2, th, ph, params));
        }
      }

    } else if (params->out2 > 30) {
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;
        // if (is_charge_conserved(params->in1, params->in2,
        // params->out1,params->out2))
        {

          // Real gluon emission for gaugino-antisquark production
          // PDF not SC like without *(xa*xb)
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (gausq_qg_gluon_minus_dip(s, mi2, pt2, th, ph, params));
        }
      }
    }
  }
  //*/
#endif

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return 0.5 * 1 / (72 * s) * (4 * M_PI * aS(params->murs, params->set)) *
           (4 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;
  } else if (is_squark_gaugino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    auto tmp = (sig * djac / (4 * sqrt_lambda(s, 0.0, 0.0))) *
               (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s;
    return tmp;
  } else {
    return dij * 0.5 * 4 / (72 * s) * (4 * M_PI * aS(params->murs, params->set)) * sig * djac;
  }
}

// Integrand for real gluon emission for joint resummation.
// No dipoles -> cf. addition to joint res.
double IG_dlnM2j(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;
  if (djac == 0.0) {
    return 0.0;
  }

  // Symmetry factor.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double real_gluon0 = 0.0;
  double real_gluon1 = 0.0;

  double sig = 0.0;
  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;
      // Test charge conservation

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          cout << "No joint resummation for gaugino gluino production." << endl;
          exit(1);
          // real_gluon0 = real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 0,
          // params);
          // real_gluon1 = real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 1,
          // params);
        } else if (is_squark_gaugino(params->out1, params->out2)) {
          cout << "No joint resummation for gaugino squark production." << endl;
          exit(1);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          real_gluon0 = real_gluon_sleptons(s, mi2, pt2, th, ph, 0, params);
          real_gluon1 = real_gluon_sleptons(s, mi2, pt2, th, ph, 1, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          real_gluon0 = real_gluon_gauginos(s, mi2, pt2, th, ph, 0, params);
          real_gluon1 = real_gluon_gauginos(s, mi2, pt2, th, ph, 1, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          real_gluon0 = real_gluon_leptons(s, mi2, pt2, th, ph, 0, params);
          real_gluon1 = real_gluon_leptons(s, mi2, pt2, th, ph, 1, params);
        }

        sig += (qa[0][i0] * qb[1 - params->ic][i1] + qb[params->ic][i0] * qa[1][i1]) *
               (real_gluon0 + real_gluon1);
      }
    }
  }
  // Flux, symmetry factor, spin and color average
  return dij * pow(M_PI, 2) / 4.5 / s * sig * djac * aS(params->murs, params->set) / (2 * M_PI);
}

// Integrand for real quark and antiquark emission.
double IQ_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;

  if (djac == 0.0) {
    return 0.0;
  }

  // Symmetry factor.
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Numerical cutoff for stability.
  // Due to the finite width and interference terms
  // of on- and off-shell contributions
  // non stable integration can occure
  // checked that dipole reproduces the real emission contribution
  // completely close to the divergent region
  if (pt2 < 1.0E-6) {
    return 0.0;
  }

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;
  double color_average_born = 1.0 / 3.0 * 1.0 / 3.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {

          // real quark emission.
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (color_average_2to3 * (real_quark_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                                        real_quark_gaugino_gluino(s, mi2, pt2, th, ph, 1, params)

                                            )

                  // Corresponding dipoles.
                  - color_average_born * (Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                                   PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                                          Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK,
                                                   PARTON_QUARK, s, mi2, pt2, th, ph, 1, params)

                                          + Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                                     PARTON_QUARK, s, mi2, pt2, th, ph, 0, params) +
                                          Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO,
                                                   PARTON_QUARK, s, mi2, pt2, th, ph, 1, params))

                 );

          // real antiquark emission.
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (color_average_2to3 * (real_quarkb_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                                        real_quarkb_gaugino_gluino(s, mi2, pt2, th, ph, 1, params))

                  // corresponding dipoles.
                  - color_average_born *
                        (Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                  s, mi2, pt2, th, ph, 0, params) +
                         Dip_GLGA(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                  s, mi2, pt2, th, ph, 1, params)

                         + Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO, PARTON_ANTIQUARK, s,
                                    mi2, pt2, th, ph, 0, params) +
                         Dip_GLGA(INITIAL_FINAL, PARTON_GLUON, PARTON_GLUINO, PARTON_ANTIQUARK, s,
                                  mi2, pt2, th, ph, 1, params))

                 );

          // Real quark and antiquark emission for Drell-Yan like processes.
        } else if (is_slepton_slepton(params->out1, params->out2)) {

          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_sleptons(s, mi2, pt2, th, ph, 1, params) -

                  (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_QUARK, s, mi2,
                                pt2, th, ph, 0, params) +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_QUARK, s, mi2,
                                pt2, th, ph, 1, params)));

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_sleptons(s, mi2, pt2, th, ph, 1, params) -
                  (Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                s, mi2, pt2, th, ph, 0, params) +
                   Dip_SLEPTONS(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                s, mi2, pt2, th, ph, 1, params)));

        } else if (is_gaugino_gaugino(params->out1, params->out2)) {

          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_gauginos(s, mi2, pt2, th, ph, 1, params) -
                  (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_QUARK, s, mi2,
                                pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_QUARK, PARTON_QUARK, s, mi2,
                                pt2, th, ph, 1, params)));

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_gauginos(s, mi2, pt2, th, ph, 1, params) -
                  (Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                s, mi2, pt2, th, ph, 0, params) +
                   Dip_GAUGINOS(INITIAL_INITIAL, PARTON_GLUON, PARTON_ANTIQUARK, PARTON_ANTIQUARK,
                                s, mi2, pt2, th, ph, 1, params)));

        } else if (is_lepton_lepton(params->out1, params->out2)) {

          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_leptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_leptons(s, mi2, pt2, th, ph, 1, params) -
                  DipQB_leptons(s, mi2, pt2, th, ph, 0, params) -
                  DipQB_leptons(s, mi2, pt2, th, ph, 1, params));

          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_leptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_leptons(s, mi2, pt2, th, ph, 1, params) -
                  DipQA_leptons(s, mi2, pt2, th, ph, 0, params) -
                  DipQA_leptons(s, mi2, pt2, th, ph, 1, params));
        }
      }
    }
  }

  if (is_squark_gaugino(params->out1, params->out2)) {
#ifdef SQGA_GG
    if (params->out1 > 30) {
      // Loop over quark outgoing types
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;
        sig += (ga * gb) *
               ((gausq_gg_gluon_minus_dip(s, mi2, pt2, th, ph, params)) +
                (gausq_gg_gluon_minus_dip(s, mi2, pt2, th, ph + M_PI, params))) /
               2.;
      }
    } else if (params->out2 > 30) {
      // Loop over quark outgoing types
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;
        sig += (ga * gb) *
               ((gausq_gg_gluon_minus_dip(s, mi2, pt2, th, ph, params)) +
                (gausq_gg_gluon_minus_dip(s, mi2, pt2, th, ph + M_PI, params))) /
               2.;
      }
    }
#endif
#ifdef SQGA_UUB
    for (int i0 = QA_MIN; i0 < QA_MAX; i0++) {
      for (int i1 = QB_MIN; i1 < QB_MAX; i1++) {
        params->in1 = i0;
        params->in2 = i1;
        double g3s = std::norm(params->gqq[0][0].R);
        // int pdf = (params->ic + (params->out1 > 30 ? 0 :1) )%2; // flip qq
        // vs qbqb pdf
        int pdf = ((params->out1 > 30 ? 0 : 1)) % 2;

        sig += (qa[(1 + pdf) % 2][i1] * qb[(params->ic + pdf) % 2][i0] +
                qb[(1 - params->ic + pdf) % 2][i1] * qa[(0 + pdf) % 2][i0]) *
               (gausq_qqb_antiquark_minus_dip(s, mi2, pt2, th, ph, params) +
                gausq_qqb_antiquark_minus_dip(s, mi2, pt2, th, ph + M_PI, params)) /
               2.;
      }
    }
#endif
#ifdef SQGA_UU
    for (int i0 = QA_MIN; i0 < QA_MAX; i0++) {
      for (int i1 = QB_MIN; i1 < QB_MAX; i1++) {
        params->in1 = i0;
        params->in2 = i1;
        double g3s = std::norm(params->gqq[0][0].R);
        int pdf, pdf_ic;
        pdf = ((params->out1 > 30 ? 0 : 1)) % 2;
        pdf_ic = (params->ic + (params->out1 > 30 ? 0 : 1)) % 2;
        double f = 1.;
        f = 0.5;

        sig += f * (qa[pdf][i1] * qb[pdf_ic][i0] + qb[pdf_ic][i1] * qa[pdf][i0]) *
               (gausq_qq_quark_minus_dip(s, mi2, pt2, th, ph, params) +
                gausq_qq_quark_minus_dip(s, mi2, pt2, th, ph + M_PI, params)) /
               2.;
      }
    }
#endif
  }

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return 0.5 * 1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
           (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;
  } else if (is_squark_gaugino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (sig * djac / (4 * sqrt_lambda(s, 0.0, 0.0))) *
           (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s;
  } else {
    return dij * 0.5 * 4 / (96 * 2 * s) * (4 * M_PI * aS(params->murs, params->set)) * sig * djac;
  }
}

// resonant diagrams for gaugino-gluino production
double IQ23_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb;

  // Different three-particle phase space for on-shell subtraction
  // for the associated gaugino-gluino production.
  // double s2_13, t1_13, s1_13, phi_13, djac13;
  double s2_23, t1_23, s1_23, phi_23, djac23;
  if (is_gaugino_gluino(params->out1, params->out2) ||
      is_squark_gaugino(params->out1, params->out2)) {
    // djac13 = dPS3_ONSHELL13_dlnM2(xa, xb, s2_13, t1_13, s1_13, phi_13, x,
    // params);
    djac23 = dPS3_ONSHELL23_dlnM2(xa, xb, s2_23, t1_23, s1_23, phi_23, x, params);
  }

  const double s = xa * xb * params->sh;
  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sigOS = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {

        // On-shell subtraction for the associated gaugino gluino production.
        if (is_gaugino_gluino(params->out1, params->out2)) {

          if (abs(djac23) > 1E-15) {
            sigOS +=
                (qa[0][i0] * gb + qb[params->ic][i0] * ga) * color_average_2to3 *
                real_quark_gaugino_gluino_onshell_23(s, s2_23, t1_23, s1_23, phi_23, 1, params);

            sigOS +=
                (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) * color_average_2to3 *
                real_quarkb_gaugino_gluino_onshell_23(s, s2_23, t1_23, s1_23, phi_23, 1, params);
          }
        }
      }
    }
  }
#ifdef SQGA_GG
#ifdef ONSUB
  if (params->out1 > 30) {
    // loop over outgoing quark types
    for (int i0 = 0; i0 < 5; i0++) {
      int i1 = 0;
      params->in1 = i0;
      params->in2 = i1;
      if (is_squark_gaugino(params->out1, params->out2)) {
        if (abs(djac23) > 1E-15) {
          sigOS +=
              (ga * gb) *
              (real_quarkb_gaugino_squark_onshell_23(s, s2_23, t1_23, s1_23, phi_23, 1, params) +
               real_quarkb_gaugino_squark_onshell_23(s, s2_23, t1_23, s1_23, phi_23 + M_PI, 1,
                                                     params)) /
              2.;
        }
      }
    }
  } else if (params->out2 > 30) {
    // loop over outgoing quark types
    for (int i1 = 0; i1 < 5; i1++) {
      int i0 = 0;
      params->in1 = i0;
      params->in2 = i1;
      if (is_squark_gaugino(params->out1, params->out2)) {
        if (abs(djac23) > 1E-15) {
          sigOS +=
              (ga * gb) *
              (real_quarkb_gaugino_squark_onshell_23(s, s2_23, t1_23, s1_23, phi_23, 1, params) +
               real_quarkb_gaugino_squark_onshell_23(s, s2_23, t1_23, s1_23, phi_23 + M_PI, 1,
                                                     params)) /
              2.;
        }
      }
    }
  }
#endif
#endif
#ifdef SQGA_UUB
#ifdef ONSUB
  for (int i0 = QA_MIN; i0 < QA_MAX; i0++) {
    for (int i1 = QB_MIN; i1 < QB_MAX; i1++) {
      params->in1 = i0; // q
      params->in2 = i1; // qb
      int pdf = ((params->out1 > 30 ? 0 : 1)) % 2;

      double g3s = std::norm(params->gqq[0][0].R);
      if (is_squark_gaugino(params->out1, params->out2)) {
        if (abs(djac23) > 1E-15) {
          sigOS += (qa[(1 + pdf) % 2][i1] * qb[(params->ic + pdf) % 2][i0] +
                    qb[(1 - params->ic + pdf) % 2][i1] * qa[(0 + pdf) % 2][i0]) *
                   (real_quarkb_gaugino_squark_uub_UXub_onshell_23(s, s2_23, t1_23, s1_23, phi_23,
                                                                   1, params) +
                    real_quarkb_gaugino_squark_uub_UXub_onshell_23(s, s2_23, t1_23, s1_23,
                                                                   phi_23 + M_PI, 1, params)) /
                   2. * pow2(aS(params->murs, params->set)) /
                   (pow2(4.0 * (M_PI)*aS(params->murs, params->set)) / g3s);
        }
      }
    }
  }
#endif
#endif
#ifdef SQGA_UU
#ifdef ONSUB
  for (int i0 = QA_MIN; i0 < QA_MAX; i0++) {
    for (int i1 = QB_MIN; i1 < QB_MAX; i1++) {
      params->in1 = i0; // q
      params->in2 = i1; // q
      double g3s = std::norm(params->gqq[0][0].R);
      int pdf, pdf_ic;
      pdf = ((params->out1 > 30 ? 0 : 1)) % 2;
      pdf_ic = (params->ic + (params->out1 > 30 ? 0 : 1)) % 2;
      double f = 1.;
      f = 0.5;

      if (is_squark_gaugino(params->out1, params->out2)) {
        if (abs(djac23) > 1E-15) {
          sigOS += f * (qa[pdf][i1] * qb[pdf_ic][i0] + qb[pdf_ic][i1] * qa[pdf][i0]) *
                   (real_quark_gaugino_squark_uu_UXu_onshell_23(s, s2_23, t1_23, s1_23, phi_23, 1,
                                                                params) +
                    real_quark_gaugino_squark_uu_UXu_onshell_23(s, s2_23, t1_23, s1_23,
                                                                phi_23 + M_PI, 1, params)) /
                   2. * pow2(aS(params->murs, params->set)) /
                   (pow2(4.0 * (M_PI)*aS(params->murs, params->set)) / g3s);
        }
      }
    }
  }
#endif
#endif
  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) * (4.0 * M_PI * aS(params->murs, params->set)) *
            (4.0 * M_PI * aS(params->murs, params->set)) * sigOS / g3s * params->mis * djac23);
  } else if (is_squark_gaugino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (sigOS * djac23 / (2 * sqrt_lambda(s, 0.0, 0.0))) *
           (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s * params->mis;
  } else {
    return 0.0;
  }
}

// resonant diagrams for gaugino-gluino production
double IQ13_dlnM2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb;

  // Different three-particle phase space for on-shell subtraction
  // for the associated gaugino-gluino production.
  double s2_13, t1_13, s1_13, phi_13, djac13;
  djac13 = dPS3_ONSHELL13_dlnM2(xa, xb, s2_13, t1_13, s1_13, phi_13, x, params);
  const double s = xa * xb * params->sh;

  // PDFs.
  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sigOS = 0.0;

  // Color average for 2->3 and LO.
  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {

        // On-shell subtraction for the associated gaugino gluino production.
        if (is_gaugino_gluino(params->out1, params->out2)) {

          if (abs(djac13) > 1E-15) {
            sigOS +=
                (qa[0][i0] * gb + qb[params->ic][i0] * ga) * color_average_2to3 *
                real_quark_gaugino_gluino_onshell_13(s, s2_13, t1_13, s1_13, phi_13, 2, params);

            sigOS +=
                (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) * color_average_2to3 *
                real_quarkb_gaugino_gluino_onshell_13(s, s2_13, t1_13, s1_13, phi_13, 2, params);
          }
        }
      }
    }
  }

/*
 *  Symmetrization in phi is needed for (good) convergence of gluino onshell
 * resonance.
 */
#ifdef SQGA_UUB
#ifdef ONSUB
  for (int i0 = QA_MIN; i0 < QA_MAX; i0++) {
    for (int i1 = QB_MIN; i1 < QB_MAX; i1++) {
      params->in1 = i0;
      params->in2 = i1;
      int pdf = ((params->out1 > 30 ? 0 : 1)) % 2;
      double g3s = std::norm(params->gqq[0][0].R);
      if (is_squark_gaugino(params->out1, params->out2)) {
        if (abs(djac13) > 1e-15) {
          sigOS += (qa[(1 + pdf) % 2][i1] * qb[(params->ic + pdf) % 2][i0] +
                    qb[(1 - params->ic + pdf) % 2][i1] * qa[(0 + pdf) % 2][i0]) *
                   (real_quarkb_gaugino_squark_uub_UXub_onshell_13(s, s2_13, t1_13, s1_13, phi_13,
                                                                   2, params) +
                    real_quarkb_gaugino_squark_uub_UXub_onshell_13(s, s2_13, t1_13, s1_13,
                                                                   phi_13 + M_PI, 2, params)) /
                   2. * pow2(aS(params->murs, params->set)) /
                   (pow2(4.0 * (M_PI)*aS(params->murs, params->set)) / g3s);
        }
      }
    }
  }
#endif
#endif

  // Flux, symmetry factor, spin and color average.
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) * (4.0 * M_PI * aS(params->murs, params->set)) *
            (4.0 * M_PI * aS(params->murs, params->set)) * sigOS / g3s * params->mis * djac13);
  } else if (is_squark_gaugino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return (sigOS * djac13 / (2 * sqrt_lambda(s, 0.0, 0.0))) *
           (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s * params->mis;
  } else {
    return 0.0;
  }
}

// Real quark and antiquark emission for joint resummation.
// No dipoles -> cf. addition to the joint.
double IQ_dlnM2j(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj;
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dlnM2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;
  if (djac == 0.0) {
    return 0.0;
  }
  double dij = 1.0;
  if (params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;

  double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;
  double color_average_born = 1.0 / 3.0 * 1.0 / 3.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      // Tests charge conservation.

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_gaugino_gluino(params->out1, params->out2)) {
          cout << "No joint resummation for gaugino gluino production." << endl;
          exit(1);
          // sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
          //   ( color_average_2to3 * (real_quark_gaugino_gluino(s, mi2, pt2,
          //   th, ph, 0, params) +
          //                           real_quark_gaugino_gluino(s, mi2, pt2,
          //                           th, ph, 1, params)));

          // sig +=  (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
          //   (  color_average_2to3 * (real_quarkb_gaugino_gluino(s, mi2, pt2,
          //   th, ph, 0, params) +
          //                            real_quarkb_gaugino_gluino(s, mi2, pt2,
          //                            th, ph, 1, params)));
        } else if (is_squark_gaugino(params->out1, params->out2)) {
          cout << "No joint resummation for squark gaugino production." << endl;
          exit(1);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_sleptons(s, mi2, pt2, th, ph, 1, params));
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_sleptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_sleptons(s, mi2, pt2, th, ph, 1, params));
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_gauginos(s, mi2, pt2, th, ph, 1, params));
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_gauginos(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_gauginos(s, mi2, pt2, th, ph, 1, params));
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                 (real_quark_leptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quark_leptons(s, mi2, pt2, th, ph, 1, params));
          sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                 (real_quarkb_leptons(s, mi2, pt2, th, ph, 0, params) +
                  real_quarkb_leptons(s, mi2, pt2, th, ph, 1, params));
        }
      }
    }
  }

  // Flux, symmetry factor, spin and color average
  if (is_gaugino_gluino(params->out1, params->out2)) {
    double g3s = std::norm(params->gqq[0][0].R);
    return 0.5 * 1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
           (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;
  } else {
    return dij * 0.5 * 4 / (96 * 2 * s) * (4 * M_PI * aS(params->murs, params->set)) * sig * djac;
  }
}

// Threshold resummation.
// Collinear improved version for Drell-Yan like processes
// and "ordinary" resummation for the associated gaugino gluino production.
double IR_dlnM2(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;
  double m1;
  double m2;

  // macro to set final state masses.
  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);

  // ivariant mass squared.
  double s = params->mis;

  if (s < pow(m1 + m2, 2)) {
    cout << "IvMlln: M2 too small\n";
    exit(0);
  } else if (s > sh) {
    cout << "IvMlln: M2 too large\n";
    exit(1);
  }

  // Jacobian initialization
  complex<double> djac(389379304.0, 0.0);

  // Inverse Mellin
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm = (s / sh / .09) + .09 - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;

  // Final jacobian.
  djac /= 8.0 * M_PI * s;
  if (is_squark_gaugino(params->out1, params->out2)) {
    // Ordinary NLL.
    return M_1_PI * imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nll_xs(nm, s, t, params));
  }
  if (is_gaugino_gluino(params->out1, params->out2)) { // NLL
    return M_1_PI * imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nll_xs(nm, s, t, params));
  }
  // Collinear improved version for Drell-Yan like processes
  return M_1_PI *
         imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_xs2(nm, s, t, params)); // NLL impr.
}

// Ordinary resummation for Drell-Yan like processes.
double IR_nll_unimproved_dlnM2(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;
  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);
  double s = params->mis;

  if (s < pow(m1 + m2, 2)) {
    cout << "IvMlln: M2 too small\n";
    exit(0);
  } else if (s > sh) {
    cout << "IvMlln: M2 too large\n";
    exit(1);
  }

  // Jacobian initialization
  complex<double> djac(389379304.0, 0.0);

  // Inverse Mellin
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm = (s / sh / .09) + .09 - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;

  djac /= 8.0 * M_PI * s;

  return M_1_PI * imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nll_xs(nm, s, t, params));
}

// Ordinary resummation for Drell-Yan like processes.
double IR_nnll_unimproved_dlnM2(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm;

  const double sh = params->sh;
  double m1;
  double m2;

  SET_MASS;

  const double m1s = pow(m1, 2);
  const double m2s = pow(m2, 2);
  double s = params->mis;

  if (s < pow(m1 + m2, 2)) {
    cout << "IvMlln: M2 too small\n";
    exit(0);
  } else if (s > sh) {
    cout << "IvMlln: M2 too large\n";
    exit(1);
  }

  // Jacobian initialization
  complex<double> djac(389379304.0, 0.0);

  // Inverse Mellin
  const complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  const complex<double> nm = (s / sh / .09) + .09 - params->a1min - ephi * log(x[0]);
  djac *= ephi / x[0];

  // t integration
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  double t = (tmax - tmin) * x[1] + tmin;
  djac *= tmax - tmin;

  djac /= 8.0 * M_PI * s;

  return M_1_PI * imag(djac * pow(s / sh, -nm + 1.0) * Thadronic_nnll_xs(nm, s, t, params));
}

// Integration.
// If a specific integration is not stable increase the number of calls (3rd
// argument of Integration()).
void hadronic_xs_dlnM2(double &res, double &err, double &chi2, int Flag, int Verb,
                       Parameters *params) {

  // Definition of the integral
  switch (Flag) {

  case 0:
    cout << endl;
    cout << "*********************" << endl;
    cout << "* LO invariant mass *" << endl;
    cout << "*********************" << endl;
    Integration(&IB_dlnM2, 2, 10000, params->precision / 3., params->abs_precision, res, err,
                params, 0.1, 5, 5);
    break;
  case 1:
    cout << endl;
    cout << "*******************************************" << endl;
    cout << "* virtual corrections + integrated dipole *" << endl;
    cout << "*******************************************" << endl;
#ifdef VNLO
    Integration(&IV_dlnM2, 2, 800, params->precision, params->abs_precision, res, err, params, 0.2);
#else
    res = 0;
    err = 0;
#endif
    break;
  case 2:
    cout << endl;
    cout << "*****************************" << endl;
    cout << "* collinear remainder (P+K) *" << endl;
    cout << "*****************************" << endl;
#ifdef PK
    Integration(&IC_dlnM2, 3, 15000, params->precision, params->abs_precision, res, err, params);
#else
    res = 0;
    err = 0;
#endif
    break;

  case 3:
    cout << endl;
    cout << "********************************" << endl;
    cout << "* real gluon emission - dipole *" << endl;
    cout << "********************************" << endl;
#ifdef RNLOg
    Integration(&IG_dlnM2, 5, 20000, params->precision, params->abs_precision, res, err, params);
    break;
#else
    res = 0;
    err = 0;
#endif

  case 4:
    cout << endl;
    cout << "**************************************************" << endl;
    cout << "* real light quark emission - dipole - (onshell) *" << endl;
    cout << "**************************************************" << endl;

#ifdef RNLOq
    // on-shell remainder sizable -> increase precision by increasing number of
    // calls
    if (is_gaugino_gluino(params->out1, params->out2)) {
      double res13 = 0.0, res23 = 0.0, err13 = 0.0, err23 = 0.0;
      Integration(&IQ_dlnM2, 5, 50000, params->precision, params->abs_precision, res, err, params,
                  0.25, 5, 5);
#ifndef OS_DR
      Integration(&IQ23_dlnM2, 5, 200000, params->precision, params->abs_precision, res23, err,
                  params, 0.25, 5, 5);
      int calls = 50000;
      for (int i = 0; i < 12; i++) {
        if (i == 8 || i == 11) {
          continue; // no increase of calls needed
        }
        if (params->mGL < params->mSQ[i]) {
          calls = 100000; // on-shell remainder sizable
          break;
        }
      }
      Integration(&IQ13_dlnM2, 5, calls, params->precision, params->abs_precision, res13, err,
                  params, 0.25, 5, 5);
#endif
      res = res + res13 + res23;
      err = sqrt(pow(err, 2) + pow(err13, 2) + pow(err23, 2));
    } else if (is_squark_gaugino(params->out1, params->out2)) {
      double res13 = 0.0, res23 = 0.0, err13 = 0.0, err23 = 0.0;

      Integration(&IQ_dlnM2, 6, 50000, params->precision, params->abs_precision, res, err, params,
                  0.2, 5, 3);

#ifndef OS_DR
      Integration(&IQ23_dlnM2, 6, 100000, params->precision, params->abs_precision, res23, err23,
                  params, 0.4, 5, 3);

      Integration(&IQ13_dlnM2, 6, 100000, params->precision, params->abs_precision, res13, err13,
                  params, 0.2, 5, 3);
#endif

      res = res + res13 + res23;
      err = sqrt(pow2(err) + pow2(err13) + pow2(err23));
    }

    else {
      Integration(&IQ_dlnM2, 5, 40000, params->precision, params->abs_precision, res, err, params,
                  0.25, 5, 5);
    }
#else
    res = 0;
    err = 0;
#endif
    break;
  case 5: // threshold resummation (improved for DY; ordinary for
          // gaugino-gluino)
    cout << endl;
    cout << "************************************************" << endl;
    cout << "* threshold resummation (improved) - expansion *" << endl;
    cout << "************************************************" << endl;
    Integration(&IR_dlnM2, 2, 20000, params->precision, params->abs_precision, res, err, params,
                0.2, 5, 5);
    break;
  case 6: // ordinary threshold resummation (including nnll logs)
    cout << endl;
    cout << "******************************************************" << endl;
    cout << "* threshold resummation (unimproved nll) - expansion *" << endl;
    cout << "******************************************************" << endl;
    Integration(&IR_nll_unimproved_dlnM2, 2, 20000, params->precision, params->abs_precision, res,
                err, params, 0.2, 5, 3);
    break;
  case 7: // ordinary threshold resummation (including nnll logs)
    cout << endl;
    cout << "*******************************************************" << endl;
    cout << "* threshold resummation (unimproved nnll) - expansion *" << endl;
    cout << "*******************************************************" << endl;
    Integration(&IR_nnll_unimproved_dlnM2, 2, 20000, params->precision, params->abs_precision, res,
                err, params, 0.2, 5, 5);
    break;
  default:
    break;
    cout << "hadronic_xs_dlnM2: Flag=" << Flag << endl;
    exit(0);
  }
  return;
}
