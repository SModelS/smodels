// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2015 David R. Lamprea.
// Copyright 2011-2015 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the hadronic cross section transverse-momentum distribution at LO,
// NLO (collinear, virtual, gluon emission and quark emission) and NLL.

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

#define CMLLN 0.9 // C coeff. of the inverse Mellin transform
#define VBSSL 2.9 // V coeff. of the inverse Bessel transform

// only for gaugino-gluino
#define QUARK_EMISSION
#define ANTIQUARK_EMISSION

double dPS3_dPT2(double &xa, double &xb, double &M2, double &pt2, double &Tp, double &Pp, double *x,
                 Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = std::pow(m1, 2);
  const double m2s = std::pow(m2, 2);
  pt2 = params->pts;

  if (pt2 < 0) {
    std::cout << "IvMlln: PT too small\n";
    exit(0);
  } else if (pt2 > std::pow(sh - std::pow(m1 + m2, 2), 2) / sh * .25 &&
             params->out1 < 20) { // not for leptons
    std::cout << "IvMlln: PT too large\n";
    exit(1);
  }

  // Jacobian initialization
  double djac = 389379304.0;

  // Integration variable xa, xb and M2
  double M2min = std::pow(m1 + m2, 2);
  double M2max = sh;
  // If we have leptons, we have zero mass => special case
  if (is_lepton_lepton(params->out1, params->out2)) {
    M2min = std::pow(params->Minv_min, 2.0);
    M2max = std::pow(params->Minv_max, 2.0);
  }
  double xamin = (M2min + (pt2 + std::sqrt(pt2 * (pt2 + M2min))) * 2.0) / sh;
  double xamax = 1.0;
  double xbmin = xamin;
  double xbmax = 1.0;
  if (xamin == 0.0) {
    xa = (xamax - xamin) * x[0] + xamin;
    djac *= xamax - xamin;
    xbmin /= xa;
    xb = (xbmax - xbmin) * x[1] + xbmin;
    djac *= xbmax - xbmin;
    M2max *= xa * xb;
    M2max -= std::sqrt(M2max * pt2) * 2.0;
    M2 = (M2max - M2min) * x[2] + M2min;
    djac *= M2max - M2min;
  } else {
    xa = xamin * std::pow(xamax / xamin, x[0]);
    djac *= xa * std::log(xamax / xamin);
    xbmin /= xa;
    xb = xbmin * std::pow(xbmax / xbmin, x[1]);
    djac *= xb * std::log(xbmax / xbmin);
    M2max *= xa * xb;
    M2max -= std::sqrt(M2max * pt2) * 2.0;
    M2 = M2min * std::pow(M2max / M2min, x[2]);
    djac *= std::log(M2max / M2min);
  }

  const double s = xa * xb * sh;

  // Integration variable Tp
  const double Tpmin = 0.0;
  const double Tpmax = M_PI;
  Tp = (Tpmax - Tpmin) * x[3] + Tpmin;
  djac *= Tpmax - Tpmin;

  // Integration variable Pp
  const double Ppmin = 0.0;
  const double Ppmax = M_PI;
  djac *= 2.0;
  Pp = (Ppmax - Ppmin) * x[4] + Ppmin;
  djac *= Ppmax - Ppmin;

  // Phase Space factor
  djac *= kln(M2, m1s, m2s) * std::sin(Tp) / std::sqrt(std::pow(s - M2, 2) - 4.0 * s * pt2) /
          std::pow(4.0 * M_PI, 4);

  return djac;
}

std::complex<double> dPSN_dPT2(std::complex<double> &nm, double &s, double &t, double *x,
                               Parameters *params) {
  const double sh = params->sh;

  double m1;
  double m2;

  SET_MASS;

  const double m1s = std::pow(m1, 2);
  const double m2s = std::pow(m2, 2);
  const double pt2 = params->pts;

  if (pt2 < 0) {
    std::cout << "IvMlln: PT too small\n";
    exit(0);
  } else if (pt2 > std::pow(sh - std::pow(m1 + m2, 2), 2) / sh * .25 &&
             params->out1 < 20) { // not for leptons
    std::cout << "IvMlln: PT too large\n";
    exit(1);
  }

  // Jacobian initialization
  std::complex<double> djac(389379304.0, 0.0);

  // Inverse Mellin
  const std::complex<double> ephi(-M_SQRT1_2, M_SQRT1_2);
  nm = CMLLN - params->a1min - ephi * std::log(x[1]);
  djac *= ephi / x[1];

  // M2 integration
  double M2min = std::pow(m1 + m2, 2);
  double M2max = sh - 2.0 * std::sqrt(sh * pt2);
  // If we have leptons, we have zero mass => special case
  if (is_lepton_lepton(params->out1, params->out2)) {
    M2min = std::pow(params->Minv_min, 2.0);
    M2max = std::pow(params->Minv_max, 2.0);
  }
  s = M2min * std::pow(M2max / M2min, x[0]);
  djac *= std::log(M2max / M2min);

  // t integration
  const double tmin = -.5 * (s - m1s - m2s + kln(s, m1s, m2s));
  const double tmax = tmin + kln(s, m1s, m2s);
  t = (tmax - tmin) * x[2] + tmin;
  djac *= tmax - tmin;

  djac /= 8.0 * M_PI * s;

  return djac;
}

double IG_dPT2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dPT2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;

  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double real_gluon0 = 0.0;
  double real_gluon1 = 0.0;
  double sig = 0.0;

  if (!is_squark_gaugino(params->out1, params->out2)) {
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
      for (int i1 = 0; i1 < 5; i1++) {
        // Set initial state
        params->in1 = i0;
        params->in2 = i1;

        if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
          if (is_gaugino_gluino(params->out1, params->out2)) {
            real_gluon0 = real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 0, params);
            real_gluon1 = real_gluon_gaugino_gluino(s, mi2, pt2, th, ph, 1, params);
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
  } else {
#ifdef SQGA_QG
    if (params->out1 > 30) {
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;

        sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
               (gausq_qg_gluon(s, mi2, pt2, th, ph, params));
      }

    } else if (params->out2 > 30) {
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;

        sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
               (gausq_qg_gluon(s, mi2, pt2, th, ph, params));
      }
    }
#endif
  }

  if (is_squark_gaugino(params->out1, params->out2)) {
    cout << "no pt distribution for the associated production yet." << endl;
    double g3s = std::norm(params->gqq[0][0].R);
    auto tmp = (sig * djac / (4 * sqrt_lambda(s, 0.0, 0.0))) *
               (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s;
    return tmp;
  } else if (is_gaugino_gluino(params->out1, params->out2)) {
    cout << "no pt distribution for the associated production yet." << endl;
    double g3s = std::norm(params->gqq[0][0].R);
    // do not understand factor of 0.5
    // APN: Factor 0.5 due to the phase space construction dPS3 with rapidity
    // flags.
    return 0.5 * 1 / (72 * s) * (4 * M_PI * aS(params->murs, params->set)) *
           (4 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;

  } else {
    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 4.5 / s * sig * djac * aS(params->murs, params->set) /
           (2 * M_PI);
  }
}

double IQ_dPT2(double *x, size_t dim, void *jj) {
  Parameters *params = (Parameters *)jj; // conversion of type void* into IOS*
  double xa, xb, mi2, pt2, th, ph;
  const double djac = dPS3_dPT2(xa, xb, mi2, pt2, th, ph, x, params);
  const double s = xa * xb * params->sh;

  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  double ga, qa[2][6];
  pdfX(ga, qa, xa, params->mufs);
  double gb, qb[2][6];
  pdfX(gb, qb, xb, params->mufs);

  double sig = 0.0;
  if (!is_squark_gaugino(params->out1, params->out2)) {
    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
      for (int i1 = 0; i1 < 5; i1++) {
        // Set initial state
        params->in1 = i0;
        params->in2 = i1;

        // Tests charge conservation.

        if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
          if (is_gaugino_gluino(params->out1, params->out2)) {
            sig += (qa[0][i0] * gb + qb[params->ic][i0] * ga) *
                   (real_quark_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                    real_quark_gaugino_gluino(s, mi2, pt2, th, ph, 1, params));
            sig += (qa[1][i1] * gb + qb[1 - params->ic][i1] * ga) *
                   (real_quarkb_gaugino_gluino(s, mi2, pt2, th, ph, 0, params) +
                    real_quarkb_gaugino_gluino(s, mi2, pt2, th, ph, 1, params));
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
  } else {
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
        int pdf = ((params->out1 > 30 ? 0 : 1)) % 2;
        double f = 1.;
        sig += f *
               (qa[(1 + pdf) % 2][i1] * qb[(params->ic + pdf) % 2][i0] +
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
        if (i1 == i0)
          f = 0.5;
        if (i0 != i1)
          f = 0.5;
        sig += f * (qa[pdf][i1] * qb[pdf_ic][i0] + qb[pdf_ic][i1] * qa[pdf][i0]) *
               (gausq_qq_quark_minus_dip(s, mi2, pt2, th, ph, params) +
                gausq_qq_quark_minus_dip(s, mi2, pt2, th, ph + M_PI, params)) /
               2.;
      }
    }
#endif
  }
  if (is_squark_gaugino(params->out1, params->out2)) {
    cout << "no pt distribution for the associated production yet." << endl;
    double g3s = std::norm(params->gqq[0][0].R);
    return (sig * djac / (4 * sqrt_lambda(s, 0.0, 0.0))) *
           (pow2(4.0 * (M_PI)*aS(params->murs, params->set))) / g3s;
  } else if (is_gaugino_gluino(params->out1, params->out2)) {
    cout << "no pt distribution for the associated production yet." << endl;
    return 0.0;
    double g3s = std::norm(params->gqq[0][0].R);
    double color_average_2to3 = 1.0 / 3.0 * 1.0 / 8.0;
    return color_average_2to3 * 0.5 * 1.0 / 2.0 * 1.0 / 2.0 * 1.0 / (2.0 * s) *
           (4.0 * M_PI * aS(params->murs, params->set)) *
           (4.0 * M_PI * aS(params->murs, params->set)) * sig * djac / g3s;
  } else {
    // Flux, symmetry factor, spin and color average
    return dij * std::pow(M_PI, 2) / 12.0 / s * sig * djac * aS(params->murs, params->set) /
           (2 * M_PI);
  }
}

double IR_dPT2(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm; // conversion of type void* into IOS*
  std::complex<double> nm;
  double s, t;
  const double sh = params->sh;
  const double pt = std::sqrt(params->pts);
  const std::complex<double> djac = dPSN_dPT2(nm, s, t, x, params);

  // Inverse Bessel: First branch
  std::complex<double> ephi(M_SQRT1_2, M_SQRT1_2);
  std::complex<double> b = -ephi * std::log(x[3]) / pt;
  std::complex<double> bjac = ephi / x[3] / pt;
  // H function for Bessel transform
  const double v = VBSSL;
  const std::complex<double> II(0.0, 1.0);
  std::complex<double> theta = M_PI * (x[4] * (2.0 * II * v - 1.0) - II * v);
  std::complex<double> hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));
  bjac *= M_PI * (2.0 * II * v - 1.0);

  std::complex<double> result = PtXS(nm, b, s, t, params) * bjac * b * hbssl;

  // Inverse Bessel: Second branch
  ephi = std::complex<double>(M_SQRT1_2, -M_SQRT1_2);
  b = -ephi * std::log(x[3]) / pt;
  bjac = ephi / x[3] / pt;
  // H function for Bessel transform
  theta = M_PI * (x[4] * (2.0 * II * v + 1.0) - II * v);
  bjac *= -M_PI * (2.0 * II * v + 1.0);
  hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));

  result += PtXS(nm, b, s, t, params) * bjac * b * hbssl;

  return .25 * M_1_PI * std::imag(result * std::pow(s / sh, -nm + 1.0) * djac);
}

double IJ_dPT2(double *x, size_t dim, void *prm) {
  Parameters *params = (Parameters *)prm; // conversion of type void* into IOS*
  std::complex<double> nm;
  double s, t;
  const double sh = params->sh;
  const double pt = std::sqrt(params->pts);
  const std::complex<double> djac = dPSN_dPT2(nm, s, t, x, params);

  // Inverse Bessel: First branch
  std::complex<double> ephi(M_SQRT1_2, M_SQRT1_2);
  std::complex<double> b = -ephi * std::log(x[3]) / pt;
  std::complex<double> bjac = ephi / x[3] / pt;
  // H function for Bessel transform
  const double v = VBSSL;
  const std::complex<double> II(0.0, 1.0);
  std::complex<double> theta = M_PI * (x[4] * (2.0 * II * v - 1.0) - II * v);
  std::complex<double> hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));
  bjac *= M_PI * (2.0 * II * v - 1.0);

  std::complex<double> result = JtXS(nm, b, s, t, params) * bjac * b * hbssl;

  // Inverse Bessel: Second branch
  ephi = std::complex<double>(M_SQRT1_2, -M_SQRT1_2);
  b = -ephi * std::log(x[3]) / pt;
  bjac = ephi / x[3] / pt;
  // H function for Bessel transform
  theta = M_PI * (x[4] * (2.0 * II * v + 1.0) - II * v);
  bjac *= -M_PI * (2.0 * II * v + 1.0);
  hbssl = -M_1_PI * std::exp(-II * b * pt * std::sin(theta));

  result += JtXS(nm, b, s, t, params) * bjac * b * hbssl;

  return .25 * M_1_PI * std::imag(result * std::pow(s / sh, -nm + 1.0) * djac);
}

void hadronic_xs_dPT2(double &res, double &err, double &chi2, int Flag, int Verb,
                      Parameters *params) {

  switch (Flag) {
  case 0:
    res = 0.0;
    err = 0.0;
    chi2 = 0.0;
    return;
  case 1:
    res = 0.0;
    err = 0.0;
    chi2 = 0.0;
    return;
  case 2:
    res = 0.0;
    err = 0.0;
    chi2 = 0.0;
    return;
  case 3:
    Integration(&IG_dPT2, 5, 20000, params->precision, 1e-12, res, err, params);
    break;
  case 4:
    Integration(&IQ_dPT2, 5, 20000, params->precision, 1e-12, res, err, params);
    break;
  case 5:
    Integration(&IR_dPT2, 5, 60000, params->precision, 1e-12, res, err, params, 0.2, 5, 5);
    break;
  case 6:
    Integration(&IJ_dPT2, 5, 60000, params->precision, 1e-12, res, err, params, 0.2, 5, 5);
    break;
    std::cout << "hadronic_xs_dPT2: Flag=" << Flag << std::endl;
    exit(0);
  }
  return;
}
