// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2015 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Resummation module.

#include "gsl_all.h"
#include "maths.h"
#include "params.h"
#include "pdf.h"
#include "pxs.h"
#include "utils.h"
#include <cmath>
#include <complex>
#include <iostream>
#include <stdio.h>

using namespace std;

#define NFLVR 5 // Number of active flavours [ 5 RECOMMENDED ]
// Color factors and constants.
#define cA 3.0
#define cF (4.0 / 3.0)
#define NC (3.0)

// Flavor constants given in hep-ph/0201036 p.53 eq. C.12 and p.32 eq. 5.91 and
// p.37 eq. 6.17
#define I2R (0.5) // SU(3) normalization factor
#define GAMMA_Q (3.0 / 2.0 * cF)
#define GAMMA_GL (3.0 / 2.0 * cA)
#define GAMMA_SQ (2.0 * cF)
#define GAMMA_G (11.0 / 6.0 * cA - 2.0 / 3.0 * I2R * nf)
#define K_Q (7.0 / 2.0 - pow2(M_PI) / 6.0) * cF
#define K_GL (7.0 / 2.0 - pow2(M_PI) / 6.0) * cA
#define K_SQ (4.0 - pow2(M_PI) / 6.0) * cF
#define K_G ((67.0 / 18.0 - pow2(M_PI) / 6.0) * cA - 10.0 / 9.0 * I2R * nf)

// Turn on to get only the expansion up to NLO of the resummation result.
// (turn it also "on" in main.cc)
//#define EXPANSION

// Threshold resummation.
// G-exponent function for the resummation
complex<double> ThG_nll(complex<double> LL, double S, double DIJ, double Muf2, double Mur2, int set,
                        Parameters *params) {

  // Coeffs. of the beta function.
  const int nf = NFLVR;
  const double dnf = (double)nf;

  // b-coefficients; related to beta-coeffs by for instance 2 pi b0 = beta 0, etc.
  const double b0 = (11.0 * cA - 2.0 * dnf) / (12.0 * M_PI);
  const double b1 = (17.0 * pow2(cA) - 5.0 * cA * dnf - 3.0 * cF * dnf) / (24.0 * pow2(M_PI));

  // Sudakov coeff. for quarks
  static const double zeta2 = gsl_sf_zeta_int(2); //(pi^2/6)
  //   static const double zeta3 = gsl_sf_zeta_int(3); // (1.20206)

  double A1 = 2.0 * cF;
  double A2 = 2.0 * cF * ((67.0 / 18.0 - zeta2) * cA - 5.0 / 9.0 * dnf);
  if (is_squark_gaugino(params->out1, params->out2)) { // APN CHANGE CHECK TODO REVIEW
    A1 = (cA + cF);
    A2 = (cA + cF) * ((67.0 / 18.0 - zeta2) * cA - 5.0 / 9.0 * dnf);
  }

  // Soft anomalous dimension for the associated gaugino-gluino production
  const double D1_ij = DIJ;

  // Second order soft anomalous dimenstion for the associated gaugino-gluino
  // production
  // This result is only valid for absolute threshold resummation
  // (beta-threshold)!
  // const double D2_ij = cA * (-cA * (115.0/36.0 - 0.5 * zeta2 + 0.5 * zeta3) +
  // 11.0/18.0 * dnf);

  // Logarithms
  const complex<double> lmbd = aS(Mur2, set) * b0 * LL;
  const complex<double> lnlmbd = log(1.0 - 2.0 * lmbd);

  // LL Sudakov (g1 = g1_a + g2_b; symmetric for two initial state quarks).
  const complex<double> g1 =
      A1 / (2.0 * M_PI * b0 * lmbd) * (2.0 * lmbd + (1.0 - 2.0 * lmbd) * lnlmbd);
  // NLL Sudakov
  complex<double> g2 = 0;

  // NLL Sudakov for initial state emission
  // used Christoph Borschensky's notation, but A1 = 2 A1_C, A2 = 4 A2_C.
  // So we have another factor 1/2 in the A2 coefficient.
  // g2 = (g2_a + g2_b)
  g2 = A1 * b1 / (2.0 * M_PI * pow(b0, 3)) * (2.0 * lmbd + lnlmbd + 0.5 * pow2(lnlmbd)) -
       A2 / (pow2(2.0 * M_PI * b0)) * (2.0 * lmbd + lnlmbd) +
       A1 / (2.0 * M_PI * b0) * (lnlmbd * log(S / Mur2) + 2.0 * lmbd * log(Muf2 / Mur2));

  complex<double> h2 = 0.0;
  // For soft wide-angle emission of the final state gluino.
  if (is_gaugino_gluino(params->out1, params->out2) ||
      is_squark_gaugino(params->out1, params->out2)) {
    h2 = D1_ij * lnlmbd / (2.0 * M_PI * b0);
  }

  const double alphaS = aS(params->murs, params->set);

  return exp(LL * g1 + g2 + h2);
}

// Threshold resummation: (unintegrated) hadronic XS for ordinary NLL resummation (matching at NLO)
complex<double> Thadronic_nll_xs(complex<double> NM, double S, double T, Parameters *params) {

  // Majorana symmetry factor.
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  //   // set soft anomalous dimension DIJ for associated gaugino-gluino production
  //   double m1, m2;
  //   SET_MASS;
  //   const double U = -S - T + pow2(m1) + pow2(m2);
  //   const double T13 = log((params->mGLs - T) / (params->mGL * sqrt(S))) - 0.5;
  //   const double T23 = log((params->mGLs - U) / (params->mGL * sqrt(S))) - 0.5;
  //   const double DIJ = cA * (T13 + T23);

  // Coeffs. of the beta function.
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;

  double m1, m2;
  SET_MASS;
  const double U = -S - T + pow2(m1) + pow2(m2);
  double DIJ = 0.0;
  // squark mass
  double msq, msqs;

  // set soft anomalous dimension DIJ
  if (is_gaugino_gluino(params->out1, params->out2)) {
    const double T13 = log((params->mGLs - T) / (params->mGL * sqrt(S))) - 0.5;
    const double T23 = log((params->mGLs - U) / (params->mGL * sqrt(S))) - 0.5;
    DIJ = cA * (T13 + T23);
  } else if (is_squark_gaugino(params->out1, params->out2)) {
    if (params->out1 > 30) {
      msq = params->mSQ[params->out1 - 31];
      msqs = params->mSQs[params->out1 - 31];
    } else if (params->out2 > 30) {
      msq = params->mSQ[params->out2 - 31];
      msqs = params->mSQs[params->out2 - 31];
    }
    // eq. (35) in soft.pdf
    DIJ = cF * (2.0 * log((msqs - T) / (msq * sqrt(S))) - 1) +
          cA * log((msqs - U) / (msqs - T)); // alpha_s / 2 / M_PI
  }

  // Getting N-space PDFs
  complex<double> gN, qN[2][6];
  pdfN(gN, qN, NM, params->afit);

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // log(Nbar) = log(NM * exp(M_Euler))
  const complex<double> NMb = NM * exp(M_EULER);
  const complex<double> lnNb = log(NMb);

  const double alphaS = aS(params->murs, params->set);

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);

  // Collinear remainder part of the hard matching coefficient for Drell-Yan
  // like processes.
  double delta_Pqq1 = -3.0 / 2.0 * log(params->mufs / S);

  // Factorization scale
  double mufs = params->mufs;

  // Gluino mass
  double mgl = params->mGL;
  double mgls = pow2(mgl);

  //   // Initialize different parts of the hard matching coefficient.
  //   double H1_coll = 0.0;
  //   double H1 = 0.0;
  //   double H0 = 0.0;
  //   double brn = 0.0;
  //   double vir = 0.0;

  if (is_squark_gaugino(params->out1, params->out2)) {
    // Initialize different parts of the hard matching coefficient.
    complex<double> H1_coll = 0.0;
    complex<double> H1 = 0.0;
    complex<double> H0 = 0.0;
    complex<double> brn = 0.0;
    complex<double> vir = 0.0;

    if (params->out1 > 30) {
      // this case we have a squark as particle 1, meaning that we have a quark in the initial state
      for (int i0 = 0; i0 < 5; i0++) {
        int i1 = 0;
        params->in1 = i0;
        params->in2 = i1;

        double g3s = std::norm(params->gqq[0][0].R);
        brn = (4.0 * alphaS * M_PI) * born_gasq(S, T, params);
        vir = (qN[0][i0] * gN + qN[params->ic][i0] * gN) * (8.0 * pow2(M_PI)) *
              (4.0 * alphaS * M_PI) / g3s / g3s * 96. *
              (Virt_gaugino_squark(S, T, params) + DipI_gasq(S, T, params));

        // Mellin transform of the insertion operators P and K in the large N limit.
        // Only N-independent terms.
        // alpha_s / 2 / M_PI
        // eq. (32) in hard.pdf
        complex<double> P_T_quark = (qN[0][i0] * gN) * (-3.0 / 4.0) *
                                    (2.0 * cF * log(mufs / (msqs - T)) - cA * log(S / (msqs - T)));
        complex<double> P_U_quark = (gN * qN[params->ic][i0]) * (-3.0 / 4.0) *
                                    (2.0 * cF * log(mufs / (msqs - U)) - cA * log(S / (msqs - U)));
        // eq. (33) in hard.pdf
        complex<double> P_T_gluon =
            (gN * qN[params->ic][i0]) * (-beta0 / 2.0) * log(pow2(mufs) / (S * (msqs - T)));
        complex<double> P_U_gluon =
            (qN[0][i0] * gN) * (-beta0 / 2.0) * log(pow2(mufs) / (S * (msqs - U)));

        // eq. (40) in hard.pdf
        double T_quark =
            (msqs - T) / (2.0 * msqs - T) + (3.0 * msq) / (msq + sqrt(2.0 * msqs - T)) +
            log(msqs / (2.0 * msqs - T)) * (1.0 + 2.0 * log(msqs / (msqs - T))) -
            (3.0 / 2.0) * log((3.0 * msqs - T - 2.0 * msq * sqrt(2.0 * msqs - T)) / (msqs - T)) +
            2.0 * gsl_sf_dilog((2.0 * msqs - T) / msqs) - GAMMA_SQ / cF;
        // eq. (39) in hard.pdf
        complex<double> K_T_quark =
            (qN[0][i0] * gN) * (pow2(M_PI) / 2.0 * cF - GAMMA_Q - K_Q + (cF - cA / 2.0) * T_quark);

        // eq. (40) in hard.pdf
        double U_quark =
            (msqs - U) / (2.0 * msqs - U) + (3.0 * msq) / (msq + sqrt(2.0 * msqs - U)) +
            log(msqs / (2.0 * msqs - U)) * (1.0 + 2.0 * log(msqs / (msqs - U))) -
            (3.0 / 2.0) * log((3.0 * msqs - U - 2.0 * msq * sqrt(2.0 * msqs - U)) / (msqs - U)) +
            2.0 * gsl_sf_dilog((2.0 * msqs - U) / msqs) - GAMMA_SQ / cF;
        // eq. (39) in hard.pdf
        complex<double> K_U_quark = (gN * qN[params->ic][i0]) * (pow2(M_PI) / 2.0 * cF - GAMMA_Q -
                                                                 K_Q + (cF - cA / 2.0) * U_quark);

        // eq. (45) in hard.pdf
        double T_gluon =
            (msqs - T) / (2.0 * msqs - T) +
            log(msqs / (2.0 * msqs - T)) * (1.0 + 2.0 * log(msqs / (msqs - T))) +
            2.0 * gsl_sf_dilog((2.0 * msqs - T) / msqs) - GAMMA_SQ / cF +
            beta0 / cA *
                (log((3.0 * msqs - T - 2.0 * msq * sqrt(2.0 * msqs - T)) / (msqs - T)) +
                 (2.0 * msq) / (msq + sqrt(2.0 * msqs - T)));
        // eq. (44) in hard.pdf
        complex<double> K_T_gluon = (gN * qN[params->ic][i0]) * ((cA / 2.0) * pow2(M_PI) - GAMMA_G -
                                                                 K_G + (cA / 2.0) * T_gluon);
        // eq. (45) in hard.pdf
        double U_gluon =
            (msqs - U) / (2.0 * msqs - U) +
            log(msqs / (2.0 * msqs - U)) * (1.0 + 2.0 * log(msqs / (msqs - U))) +
            2.0 * gsl_sf_dilog((2.0 * msqs - U) / msqs) - GAMMA_SQ / cF +
            beta0 / cA *
                (log((3.0 * msqs - U - 2.0 * msq * sqrt(2.0 * msqs - U)) / (msqs - U)) +
                 (2.0 * msq) / (msq + sqrt(2.0 * msqs - U)));
        // eq. (44) in hard.pdf
        complex<double> K_U_gluon =
            (qN[0][i0] * gN) * ((cA / 2.0) * pow2(M_PI) - GAMMA_G - K_G + (cA / 2.0) * U_gluon);

        H1_coll = (P_T_quark + P_U_gluon) + (K_T_quark + K_U_gluon);

        H0 = (qN[0][i0] * gN + qN[params->ic][i0] * gN) * brn;
        H1 = H1_coll * brn + vir;

        if (brn == 0.0) {
          continue;
        }

        // LO and NLO hard matching coefficients.
        // Resummed XS and expansion up to NLO to avoid double counting.
#ifndef EXPANSION
        RES += (H0 + alphaS / (2.0 * M_PI) * H1);
#endif
        EXP += (alphaS / (2.0 * M_PI) * H1 +
                H0 * (1.0 + alphaS / (M_PI)*lnNb *
                                (-DIJ + (cA + cF) * (log(mufs / params->murs) + lnNb -
                                                     log(S / params->murs)))));
      }

      return (RES * ThG_nll(lnNb, S, DIJ, params->mufs, params->murs, params->set, params) - EXP) *
             dij / (2. * 96.0) / S;

    } else if (params->out2 > 30) {
      // this case we have a squark as particle 1, meaning that we have a quark in the initial state
      for (int i1 = 0; i1 < 5; i1++) {
        int i0 = 0;
        params->in1 = i0;
        params->in2 = i1;

        double g3s = std::norm(params->gqq[0][0].R);
        brn = (4.0 * alphaS * M_PI) * born_gasq(S, T, params);
        vir = (qN[1][i1] * gN + qN[1 - params->ic][i1] * gN) * (8.0 * pow2(M_PI)) *
              (4.0 * alphaS * M_PI) / g3s / g3s * 96. *
              (Virt_gaugino_squark(S, T, params) + DipI_gasq(S, T, params));

        // Mellin transform of the insertion operators P and K in the large N limit.
        // Only N-independent terms.
        // eq. (32) in hard.pdf
        complex<double> P_T_quark = (qN[1][i1] * gN) * (-3.0 / 4.0) *
                                    (2.0 * cF * log(mufs / (msqs - T)) - cA * log(S / (msqs - T)));
        complex<double> P_U_quark = (gN * qN[1 - params->ic][i1]) * (-3.0 / 4.0) *
                                    (2.0 * cF * log(mufs / (msqs - U)) - cA * log(S / (msqs - U)));
        // eq. (33) in hard.pdf
        complex<double> P_T_gluon =
            (gN * qN[1 - params->ic][i1]) * (-beta0 / 2.0) * log(pow2(mufs) / (S * (msqs - T)));
        complex<double> P_U_gluon =
            (qN[1][i1] * gN) * (-beta0 / 2.0) * log(pow2(mufs) / (S * (msqs - U)));
        // eq. (40) in hard.pdf
        double T_quark =
            (msqs - T) / (2.0 * msqs - T) + (3.0 * msq) / (msq + sqrt(2.0 * msqs - T)) +
            log(msqs / (2.0 * msqs - T)) * (1.0 + 2.0 * log(msqs / (msqs - T))) -
            (3.0 / 2.0) * log((3.0 * msqs - T - 2.0 * msq * sqrt(2.0 * msqs - T)) / (msqs - T)) +
            2.0 * gsl_sf_dilog((2.0 * msqs - T) / msqs) - GAMMA_SQ / cF;
        // eq. (39) in hard.pdf
        complex<double> K_T_quark =
            (qN[1][i1] * gN) * (pow2(M_PI) / 2.0 * cF - GAMMA_Q - K_Q + (cF - cA / 2.0) * T_quark);
        // eq. (40) in hard.pdf
        double U_quark =
            (msqs - U) / (2.0 * msqs - U) + (3.0 * msq) / (msq + sqrt(2.0 * msqs - U)) +
            log(msqs / (2.0 * msqs - U)) * (1.0 + 2.0 * log(msqs / (msqs - U))) -
            (3.0 / 2.0) * log((3.0 * msqs - U - 2.0 * msq * sqrt(2.0 * msqs - U)) / (msqs - U)) +
            2.0 * gsl_sf_dilog((2.0 * msqs - U) / msqs) - GAMMA_SQ / cF +
            beta0 / cA *
                (log((3.0 * msqs - U - 2.0 * msq * sqrt(2.0 * msqs - U)) / (msqs - U)) +
                 (2.0 * msq) / (msq + sqrt(2.0 * msqs - U)));
        // eq. (39) in hard.pdf
        complex<double> K_U_quark =
            (gN * qN[1 - params->ic][i1]) *
            (pow2(M_PI) / 2.0 * cF - GAMMA_Q - K_Q + (cF - cA / 2.0) * U_quark);
        // eq. (45) in hard.pdf
        double T_gluon =
            (msqs - T) / (2.0 * msqs - T) +
            log(msqs / (2.0 * msqs - T)) * (1.0 + 2.0 * log(msqs / (msqs - T))) +
            2.0 * gsl_sf_dilog((2.0 * msqs - T) / msqs) - GAMMA_SQ / cF +
            beta0 / cA *
                (log((3.0 * msqs - T - 2.0 * msq * sqrt(2.0 * msqs - T)) / (msqs - T)) +
                 (2.0 * msq) / (msq + sqrt(2.0 * msqs - T)));
        ;
        // eq. (44) in hard.pdf
        complex<double> K_T_gluon =
            (gN * qN[1 - params->ic][i1]) *
            ((cA / 2.0) * pow2(M_PI) - GAMMA_G - K_G + (cA / 2.0) * T_gluon);
        // eq. (45) in hard.pdf
        double U_gluon = (msqs - U) / (2.0 * msqs - U) +
                         log(msqs / (2.0 * msqs - U)) * (1.0 + 2.0 * log(msqs / (msqs - U))) +
                         2.0 * gsl_sf_dilog((2.0 * msqs - U) / msqs) - GAMMA_SQ / cF;
        // eq. (44) in hard.pdf
        complex<double> K_U_gluon =
            (qN[1][i1] * gN) * ((cA / 2.0) * pow2(M_PI) - GAMMA_G - K_G + (cA / 2.0) * U_gluon);

        H1_coll = (P_T_quark + P_U_gluon) + (K_T_quark + K_U_gluon);

        H0 = (qN[1][i1] * gN + qN[1 - params->ic][i1] * gN) * brn;
        H1 = H1_coll * brn + vir;

        if (brn == 0.0) {
          continue;
        }

        // LO and NLO hard matching coefficients.
        // Resummed XS and expansion up to NLO to avoid double counting.
#ifndef EXPANSION
        RES += (H0 + alphaS / (2.0 * M_PI) * H1);
#endif
        EXP += (alphaS / (2.0 * M_PI) * H1 +
                H0 * (1.0 + alphaS / (M_PI)*lnNb *
                                (-DIJ + (cA + cF) * (log(mufs / params->murs) + lnNb -
                                                     log(S / params->murs)))));
      }

      return (RES * ThG_nll(lnNb, S, DIJ, params->mufs, params->murs, params->set, params) - EXP) *
             dij / (2. * 96.0) / S;
    }

  } else {

    // Initialize different parts of the hard matching coefficient.
    double H1_coll = 0.0;
    double H1 = 0.0;
    double H0 = 0.0;
    double brn = 0.0;
    double vir = 0.0;

    // Sum over all possible initial states
    for (int i0 = 0; i0 < 5; i0++) {
      for (int i1 = 0; i1 < 5; i1++) {
        // Set initial state
        params->in1 = i0;
        params->in2 = i1;

        if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {

          // Hard matching coeffs
          if (is_gaugino_gluino(params->out1, params->out2)) {
            double g3s = std::norm(params->gqq[0][0].R);
            brn = (4.0 * alphaS * M_PI) / g3s * born_gagl(S, T, params);
            // (8.0 * pow2(M_PI)) to compensate the aS/(2 Pi) factor (see below)
            // leading to (4 aS pi) = g3s (same prefactor as in IV in hxs.cc)
            vir = (8.0 * pow2(M_PI)) * (4.0 * alphaS * M_PI) *
                  (Virt_gaugino_gluino(S, T, params) + DipI_gagl(S, T, params));

            // Mellin transform of the insertion operators P and K in the large N
            // limit.
            // Only N-independent terms.
            double P_T =
                -3.0 / 4.0 * ((2.0 * cF - cA) * log(mufs / S) + cA * log(mufs / (mgls - T)));
            double P_U =
                -3.0 / 4.0 * ((2.0 * cF - cA) * log(mufs / S) + cA * log(mufs / (mgls - U)));

            double K_T =
                pow2(M_PI) / 2.0 * cF - GAMMA_Q - K_Q +
                cA / 4.0 *
                    (1. + 4.0 * gsl_sf_dilog((2.0 * mgls - T) / (mgls)) +
                     (1.0 + 4.0 * log(mgls / (mgls - T)) + 2.0 * (mgls / (mgls - T))) *
                         log((mgls / (2.0 * mgls - T))) +
                     3.0 * log(1.0 + 2.0 * mgl / (mgls - T) * (mgl - sqrt(2.0 * mgls - T))) +
                     6.0 * mgl / (mgl + sqrt(2.0 * mgls - T)) - 3);

            double K_U =
                pow2(M_PI) / 2.0 * cF - GAMMA_Q - K_Q +
                cA / 4.0 *
                    (1. + 4.0 * gsl_sf_dilog((2.0 * mgls - U) / (mgls)) +
                     (1.0 + 4.0 * log(mgls / (mgls - U)) + 2.0 * (mgls / (mgls - U))) *
                         log((mgls / (2.0 * mgls - U))) +
                     3.0 * log(1.0 + 2.0 * mgl / (mgls - U) * (mgl - sqrt(2.0 * mgls - U))) +
                     6.0 * mgl / (mgl + sqrt(2.0 * mgls - U)) - 3);

            H1_coll = (P_T + P_U) + (K_T + K_U);

          } else if (is_slepton_slepton(params->out1, params->out2)) {
            brn = NC * born_sleptons(S, T, params);
            vir = NC * cF * (Virt_sleptons(S, T, params) + DipI_sleptons(S, T, params));
            // actuall H1_coll, but Jonathan shifted terms from integrated Dipole
            // to coll. remainder.
            H1_coll = 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q + cF * delta_Pqq1);
            // H1_coll = 2.0 *(cF * pow2(M_PI)/6.0 + cF * delta_Pqq1);
          } else if (is_gaugino_gaugino(params->out1, params->out2)) {
            brn = NC * born_gauginos(S, T, params);
            vir = NC * cF * (Virt_gauginos(S, T, params) + DipI_gauginos(S, T, params));
            H1_coll = 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q + cF * delta_Pqq1);
            // H1_coll = 2.0 *(cF * pow2(M_PI)/6.0 + cF * delta_Pqq1);
          } else if (is_lepton_lepton(params->out1, params->out2)) {
            brn = NC * born_leptons(S, T, params);
            vir = NC * cF * (Virt_leptons(S, T, params) + DipI_leptons(S, T, params));
            H1_coll = 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q + cF * delta_Pqq1);
            // H1_coll = 2.0 *(cF * pow2(M_PI)/6.0 + cF * delta_Pqq1);
          }

          // LO and NLO hard matching coefficients.
          H0 = brn;
          H1 = H1_coll * brn + vir;

          if (brn == 0.0) {
            continue;
          }

// Resummed XS and expansion up to NLO to avoid double counting.
#ifndef EXPANSION
          RES += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
                 (H0 + alphaS / (2.0 * M_PI) * H1);
#endif
          if (is_gaugino_gluino(params->out1, params->out2)) {
            EXP +=
                (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
                (alphaS / (2.0 * M_PI) * H1 +
                 H0 * (1.0 +
                       alphaS / (M_PI)*lnNb *
                           (-DIJ +
                            2.0 * cF * (log(mufs / params->murs) + lnNb - log(S / params->murs)))));
          } else {
            EXP += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
                   (alphaS / (2.0 * M_PI) * H1 + H0 +
                    H0 * alphaS / (M_PI)*lnNb * 2.0 * cF *
                        (log(mufs / params->murs) + lnNb - log(S / params->murs)));
          }
        }
      }
    }
  }
  // Hard matching function times exponential; expansion subctracted; color and
  // spin average; flux factor.
  return (RES * ThG_nll(lnNb, S, DIJ, params->mufs, params->murs, params->set, params) - EXP) *
         dij / 72.0 / S;
}

// G-exponent function for the resummation with g3 coefficient valid for Drell-Yan like processes
// with no colour flow (i.e. (s)lepton and electroweakinos pair production)
complex<double> ThG_nnll(complex<double> LL, double S, double Muf2, double Mur2, int set,
                         Parameters *params) {

  // Coeffs. of the beta function.
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double beta1 = 25.5 - 9.5 * dnf / 3.0;

  // b-coefficients; related to beta-coeffs by for instance 2 pi b0 = beta 0, etc.
  const double b0 = (11.0 * cA - 2.0 * dnf) / (12.0 * M_PI);
  const double b1 = (17.0 * pow2(cA) - 5.0 * cA * dnf - 3.0 * cF * dnf) / (24.0 * pow2(M_PI));
  const double b2 =
      1.0 / (64.0 * pow(M_PI, 3)) *
      (2857.0 / 54.0 * pow(cA, 3) - 1415.0 / 54.0 * pow2(cA) * dnf - 205.0 / 18.0 * cA * cF * dnf +
       pow2(cF) * dnf + pow2(dnf) / 9.0 * (79.0 / 6.0 * cA + 11.0 * cF));

  // Sudakov coeff. for quarks
  static const double zeta2 = gsl_sf_zeta_int(2); //(pi^2/6)
  static const double zeta3 = gsl_sf_zeta_int(3); // (1.20206)

  const double A1 = 2.0 * cF;
  const double A2 = 2.0 * cF * ((67.0 / 18.0 - zeta2) * cA - 5.0 / 9.0 * dnf);
  const double A3 =
      2.0 * 0.25 * cF *
      (pow2(cA) *
           (245.0 / 24.0 - 67.0 / 9.0 * zeta2 + 11.0 / 6.0 * zeta3 + 11.0 / 5.0 * pow2(zeta2)) +
       cF * dnf * (-55.0 / 24.0 + 2.0 * zeta3) +
       cA * dnf * (-209.0 / 108.0 + 10.0 / 9.0 * zeta2 - 7.0 / 3.0 * zeta3) - pow2(dnf) / 27.0);

  // Only initial state dependent.
  // Needed for NNLL.
  const double D2_i = 2.0 * cF *
                      (cA * (-101.0 / 27.0 + 11.0 / 3.0 * zeta2 + 7.0 / 2.0 * zeta3) +
                       dnf * (14.0 / 27.0 - 2.0 / 3.0 * zeta2));

  // Logarithms
  const complex<double> lmbd = aS(Mur2, set) * b0 * LL;
  const complex<double> lnlmbd = log(1.0 - 2.0 * lmbd);

  // LL Sudakov (g1 = g1_a + g2_b; symmetric for two initial state quarks).
  const complex<double> g1 =
      A1 / (2.0 * M_PI * b0 * lmbd) * (2.0 * lmbd + (1.0 - 2.0 * lmbd) * lnlmbd);
  // NLL Sudakov
  complex<double> g2 = 0;

  // NLL Sudakov for initial state emission
  // used Christoph Borschensky's notation, but A1 = 2 A1_C, A2 = 4 A2_C.
  // So we have another factor 1/2 in the A2 coefficient.
  // g2 = (g2_a + g2_b)
  g2 = A1 * b1 / (2.0 * M_PI * pow(b0, 3)) * (2.0 * lmbd + lnlmbd + 0.5 * pow2(lnlmbd)) -
       A2 / (pow2(2.0 * M_PI * b0)) * (2.0 * lmbd + lnlmbd) +
       A1 / (2.0 * M_PI * b0) * (lnlmbd * log(S / Mur2) + 2.0 * lmbd * log(Muf2 / Mur2));

  // NNLL Sudakov for initial state emission
  complex<double> g3 = 0;

  //   if (!is_gaugino_gluino(params->out1, params->out2)) {
  //     // Christoph's formulas (see appendix of Christoph Borschensky's PhD thesis).
  //     g3 = A1 * pow2(b1) / (2.0 * M_PI * pow(b0, 4)) * 1.0 / (1.0 - 2.0 * lmbd) * (2.0 *
  //     pow2(lmbd) + 2.0 * lmbd * lnlmbd + 0.5 * pow2(lnlmbd)) +
  //          (A1 * b2) / (2.0 * M_PI * pow(b0, 3)) * (2.0 * lmbd + lnlmbd + 2.0 * pow2(lmbd) / (1.0
  //          - 2.0 * lmbd)) + 2.0 * A1 / M_PI * zeta2 * lmbd / (1.0 - 2.0 * lmbd) - A2 * b1 /
  //          (pow2(2.0 * M_PI) * pow(b0, 3)) * 1.0 / (1.0 - 2.0 * lmbd) * (2.0 * lmbd * (lmbd
  //          + 1.0) + lnlmbd) + A3 / (pow(M_PI, 3) * pow2(b0)) * pow2(lmbd) / (1.0 - 2.0 * lmbd) -
  //          D2_i / (2.0 * pow2(M_PI) * b0) * lmbd / (1.0 - 2.0 * lmbd) +
  //          (A1 * b1 / (2.0 * M_PI * pow2(b0)) * (2.0 * lmbd + lnlmbd) / (1.0 - 2.0 * lmbd)) *
  //          log(S / Mur2) + A1 / (2.0 * M_PI) * (lmbd / (1.0 - 2.0 * lmbd) * pow2(log(S / Mur2)) -
  //          lmbd * pow2(log(Muf2 / Mur2))) - A2 / (2.0 * pow2(M_PI) * b0) * (lmbd / (1.0 - 2.0 *
  //          lmbd) * log(S / Mur2) - lmbd * log(Muf2 / Mur2));
  //   }

  g3 = A1 * pow2(b1) / (2.0 * M_PI * pow(b0, 4)) * 1.0 / (1.0 - 2.0 * lmbd) *
           (2.0 * pow2(lmbd) + 2.0 * lmbd * lnlmbd + 0.5 * pow2(lnlmbd)) +
       (A1 * b2) / (2.0 * M_PI * pow(b0, 3)) *
           (2.0 * lmbd + lnlmbd + 2.0 * pow2(lmbd) / (1.0 - 2.0 * lmbd)) +
       2.0 * A1 / M_PI * zeta2 * lmbd / (1.0 - 2.0 * lmbd) -
       A2 * b1 / (pow2(2.0 * M_PI) * pow(b0, 3)) * 1.0 / (1.0 - 2.0 * lmbd) *
           (2.0 * lmbd * (lmbd + 1.0) + lnlmbd) +
       A3 / (pow(M_PI, 3) * pow2(b0)) * pow2(lmbd) / (1.0 - 2.0 * lmbd) -
       D2_i / (2.0 * pow2(M_PI) * b0) * lmbd / (1.0 - 2.0 * lmbd) +
       (A1 * b1 / (2.0 * M_PI * pow2(b0)) * (2.0 * lmbd + lnlmbd) / (1.0 - 2.0 * lmbd)) *
           log(S / Mur2) +
       A1 / (2.0 * M_PI) *
           (lmbd / (1.0 - 2.0 * lmbd) * pow2(log(S / Mur2)) - lmbd * pow2(log(Muf2 / Mur2))) -
       A2 / (2.0 * pow2(M_PI) * b0) *
           (lmbd / (1.0 - 2.0 * lmbd) * log(S / Mur2) - lmbd * log(Muf2 / Mur2));

  const double alphaS = aS(params->murs, params->set);

  return exp(LL * g1 + g2 + alphaS * g3);
}

// Threshold resummation: (unintegrated) hadronic XS for ordinary NNLL resummation (matching at NLO)
// valid for Drell-Yan like processes with no colour flow (i.e. (s)lepton and electroweakinos pair
// production)
complex<double> Thadronic_nnll_xs(complex<double> NM, double S, double T, Parameters *params) {

  // Majorana symmetry factor.
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Coeffs. of the beta function.
  const int nf = NFLVR;
  const double dnf = (double)nf;
  static const double zeta3 = gsl_sf_zeta_int(3); // (1.20206)

  // Getting N-space PDFs
  complex<double> gN, qN[2][6];
  pdfN(gN, qN, NM, params->afit);

  // log(Nbar) = log(NM * exp(M_Euler))
  const complex<double> NMb = NM * exp(M_EULER);
  const complex<double> lnNb = log(NMb);

  const double alphaS = aS(params->murs, params->set);

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);

  // Collinear remainder part of the hard matching coefficient for Drell-Yan
  // like processes.
  double delta_Pqq1 = -3.0 / 2.0 * log(params->mufs / S);

  // Factorization scale
  double mufs = params->mufs;

  // Initialize different parts of the hard matching coefficient.
  //   double H1_coll = 0.0;
  double H2 = 0.0;
  double H1 = 0.0;
  double H0 = 0.0;
  double brn = 0.0;
  //   double vir = 0.0;
  double C1 = 0.0;
  double C2 = 0.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {

        // Hard matching coeffs
        if (is_slepton_slepton(params->out1, params->out2)) {
          brn = NC * born_sleptons(S, T, params);
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          brn = NC * born_gauginos(S, T, params);
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          brn = NC * born_leptons(S, T, params);
        }

        // LO and NLO hard matching coefficients.
        C1 = cF * (4.0 / 3.0 * pow2(M_PI) - 8.0 - 3.0 * log(params->mufs / S));
        C2 = cF / 720 *
             (5.0 * (-4605.0 * cA + 4599.0 * cF + 762.0 * nf) +
              20.0 * pow2(M_PI) * (188.0 * cA - 297.0 * cF - 32.0 * nf) -
              92.0 * pow(M_PI, 4) * (cA - 6.0 * cF) +
              180.0 * (11.0 * cA + 18.0 * cF - 2.0 * nf) * pow2(log(mufs / S)) -
              160.0 * (11.0 * cA - 2.0 * nf) * (6.0 - pow2(M_PI)) * log(params->murs / S) +
              80.0 * (151.0 * cA - 135.0 * cF + 2.0 * nf) * zeta3 +
              20.0 * log(mufs / S) *
                  (-51.0 * cA + 837.0 * cF + 6.0 * nf -
                   4.0 * pow2(M_PI) * (11.0 * cA + 27.0 * cF - 2.0 * nf) +
                   (-198.0 * cA + 36.0 * nf) * log(params->murs / S) +
                   216.0 * (cA - 2.0 * cF) * zeta3));

        complex<double> EXP_alpha = 4.0 * cF * log(mufs / S) * lnNb + 4.0 * cF * pow2(lnNb);

        H0 = brn;
        H1 = C1 * brn;
        H2 = C2 * brn;

        if (brn == 0.0) {
          continue;
        }

// Resummed XS and expansion up to NLO to avoid double counting.
#ifndef EXPANSION
        RES += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
               (H0 + (alphaS / (2.0 * M_PI)) * H1 + pow2(alphaS / (2.0 * M_PI)) * H2);
#endif

        EXP += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
               (H0 + alphaS / (2.0 * M_PI) * (H1 + H0 * EXP_alpha));
      }
    }
  }
  // Hard matching function times exponential; expansion subctracted; color and
  // spin average; flux factor.
  return (RES * ThG_nnll(lnNb, S, params->mufs, params->murs, params->set, params) - EXP) * dij /
         72.0 / S;
}

// Same as ThG, but with the improvement added (subleading non-diagonal
// splitting functions "readded")
// see Eq. 3.51f of Jonathan's PhD thesis.
complex<double> ThG2(const complex<double> LL, double MIS, double MURS, const int set,
                     Parameters *params) {

  // Coeff. of the beta function
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double beta1 = 25.5 - 9.5 * dnf / 3.0;

  static const double zeta2 = gsl_sf_zeta_int(2); //(pi^2/6)
  static const double zeta3 = gsl_sf_zeta_int(3); // (1.20206)

  // Sudakov coeff. for quarks
  // A_a,b
  const double A1 = 2.0 * cF;
  const double B1 = -3.0 * cF;
  const double A2 = 2.0 * cF * ((67.0 / 18.0 - zeta2) * cA - 5.0 / 9.0 * dnf);

  // logarithms
  // Below Eq. 3.39
  const complex<double> lmbd = aS(MURS, set) / (2 * M_PI) * beta0 * LL;
  const complex<double> lnlmbd = log(1.0 - 2.0 * lmbd);

  // LL + NLL Sudakov
  // Eq. 3.51
  const complex<double> g1 = 2.0 * A1 / beta0 * (1.0 + lnlmbd / (2.0 * lmbd));

  // Exactly Eq. 3.52
  const complex<double> g2 =
      (2.0 * lmbd + lnlmbd) *
          (A1 / beta0 * (log(MIS / MURS) + beta1 / pow(beta0, 2)) - A2 / pow2(beta0) + B1 / beta0) -
      B1 / beta0 * 2.0 * lmbd + .5 * A1 / beta0 * beta1 * pow(lnlmbd / beta0, 2);

  return exp(LL * g1 + g2);
}

// Same as Thadronic_xs, but now for the improved resummation version.
complex<double> Thadronic_xs2(complex<double> N, double S, double T, Parameters *params) {
  // logarithm
  const complex<double> nb = N * exp(M_EULER);
  const complex<double> ll = log(nb);

  // Mellin transform of the AP splitting functions.
  const complex<double> Pqq =
      4.0 / 3.0 * (1.5 + 1.0 / N / (N + 1.0) - 2.0 * (Psi(N + 1.0) + M_EULER));
  const complex<double> Pqg = .5 * (2.0 + N + pow(N, 2)) / N / (N + 1.0) / (N + 2.0);

  // Majorana symmetry factor.
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Getting N-space PDFs
  complex<double> gE, qE[2][6], gN, qN[2][6];
  pdfN(gN, qN, N, params->afit);
  pdfN(gE, qE, N, params->afit);

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // Evolving factorization scale (NLL)
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;

  const double alphaS = aS(params->murs, params->set);

  const complex<double> lmbd = alphaS / (2.0 * M_PI) * beta0 * log(params->mufs / S * pow(nb, 2));

  // Evolve fron N to
  pdfEvolve(gE, qE, N, lmbd);

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);

  double brn = 0.0;
  double vir = 0.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_slepton_slepton(params->out1, params->out2)) {
          brn = born_sleptons(S, T, params);
          vir = Virt_sleptons(S, T, params) + DipI_sleptons(S, T, params);
          vir += 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q) / cF * brn;
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          brn = born_gauginos(S, T, params);
          // "small" vir is Eq. 3.50 (confusing, cause it includes parts of the
          // collinear remainder!).
          // Jon has shifted subtracted the 10 and some pi squared terms arising
          // in DipI and coll.rem.
          // vir = pow(M_PI, 2) / 3.0 * brn + Virt_gauginos(S, T, params) +
          // DipI_gauginos(S, T, params);

          // Using a new version of the hard matching coeff, since I reshifted
          // the (10-pi^2) term
          // from the integrated dipole to the collinear remainder, where it
          // belongs to.
          // my vir + (10-pi^2) leads to the old vir.
          vir = Virt_gauginos(S, T, params) + DipI_gauginos(S, T, params);
          vir += 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q) / cF * brn;
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          brn = born_leptons(S, T, params);
          // vir = pow(M_PI, 2) / 3.0 * brn + Virt_leptons(S, T, params) +
          // DipI_leptons(S, T, params);
          vir = Virt_leptons(S, T, params) + DipI_leptons(S, T, params);
          vir += 2.0 * (cF * pow2(M_PI) / 2.0 - GAMMA_Q - K_Q) / cF * brn;
        }

// Hard matching coeff. and evolved PDFs.
#ifndef EXPANSION
        RES += (qE[0][i0] * qE[1 - params->ic][i1] + qE[params->ic][i0] * qE[1][i1]) *
               (brn + vir * alphaS / (2.0 * M_PI) * 4.0 / 3.0);
#endif
        EXP += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
               (brn + vir * 4.0 / 3.0 * alphaS / (2.0 * M_PI) +
                (2.0 * ll * (-2.0 * ll * 4.0 / 3.0 + 4.0) -
                 2.0 * Pqq * (2.0 * ll - log(S / params->mufs))) *
                    alphaS / (2.0 * M_PI) * brn);

        EXP += (qN[0][i0] + qN[1 - params->ic][i1] + qN[params->ic][i0] + qN[1][i1]) * gN *
               (log(S / params->mufs) - 2.0 * ll) * Pqg * alphaS / (2.0 * M_PI) * brn;
      }
    }
  }
  // Flux, symmetry factor, spin and color average. fc = 3; average 1/36; 1/2s
  // for flux; and expansion subtracted
  return (RES * ThG2(ll, S, params->murs, params->set, params) - EXP) * dij / 24.0 / S;
}

// Transverse momentum resummation.
void CFact(complex<double> &g, complex<double> q[2][6], complex<double> N, double aSmu,
           complex<double> aSbb) {
  complex<double> gdum = g;
  complex<double> qdum[2][6];
  for (size_t i0 = 0; i0 < 2; i0++) {
    for (size_t i1 = 0; i1 < 6; i1++) {
      qdum[i0][i1] = q[i0][i1];
    }
  }

  complex<double> Cqq, Cqg;
  Cqq = 1.0 + aSmu * 4.0 / 3.0 / N / (N + 1.0); // aSmu at NLL / aSbb at NNLL
  Cqg = aSbb / (N + 1.0) / (N + 2.0);           // aSbb at NLL / aSbb at NNLL

  for (size_t i0 = 0; i0 < 2; i0++) {
    for (size_t i1 = 0; i1 < 6; i1++) {
      q[i0][i1] = Cqq * qdum[i0][i1] + Cqg * gdum;
    }
  }
}

complex<double> PtG(const complex<double> LL, double MIS, double MURS, const int set) {
  // Coeff. of the beta function
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double beta1 = 25.5 - 9.5 * dnf / 3.0;

  // Sudakov coeff. for quarks
  const double A1 = 8.0 / 3.0 / beta0;
  const double B1 = -4.0 / beta0;
  const double A2 = .5 * A1 * (67.0 / 3.0 - pow(M_PI, 2) - 10.0 / 9.0 * dnf) / beta0;

  // logarithms
  const complex<double> lmbd = aS(MURS, set) / (2 * M_PI) * beta0 * LL;
  const complex<double> lnlmbd = log(1.0 - lmbd);

  // LL + NLL Sudakov
  const complex<double> g1 = A1 * (1.0 + lnlmbd / lmbd);
  const complex<double> g2 =
      (lmbd / (1.0 - lmbd) + lnlmbd) *
          (A1 * (log(MIS / MURS) + beta1 / pow(beta0, 2) * (1.0 + lnlmbd)) - A2 + B1) -
      B1 * lmbd / (1.0 - lmbd) - .5 * A1 * beta1 * pow(lnlmbd / beta0, 2);

  return exp(LL * g1 + g2);
}

complex<double> PtXS(complex<double> N, complex<double> B, double S, double T, Parameters *params) {
  // logarithm
  const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
  const complex<double> ll = log(1.0 + bb * bb);

  const complex<double> Pqq =
      4.0 / 3.0 * (1.5 + 1.0 / N / (N + 1.0) - 2.0 * (Psi(N + 1.0) + M_EULER));
  const complex<double> Pqg = .5 * (2.0 + N + pow(N, 2)) / N / (N + 1.0) / (N + 2.0);

  // Majorana symmetry factor
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Getting N-space PDFs
  complex<double> gE, qE[2][6], gN, qN[2][6];
  pdfN(gN, qN, N, params->afit);
  pdfN(gE, qE, N, params->afit);

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // Evolving factorization scale [ NLL ]
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double aa = aS(params->murs, params->set) / (2 * M_PI);
  const complex<double> lmbd = aa * beta0 * log(params->mufs / S * (1.0 + bb * bb));
  pdfEvolve(gE, qE, N, lmbd);
  CFact(gE, qE, N, aa, aa / (1.0 - lmbd));

  double brn = 0.0;
  double vir = 0.0;

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);
  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_squark_gaugino(params->out1, params->out2)) {
          cout << "No pt resummation for associated production (gaugino squark)." << endl;
          exit(1);
        } else if (is_gaugino_gluino(params->out1, params->out2)) {
          cout << "No pt resummation for associated production (gaugino gluino)." << endl;
          exit(1);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          brn = born_sleptons(S, T, params);
          vir = Virt_sleptons(S, T, params) + DipI_sleptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          brn = born_gauginos(S, T, params);
          vir = Virt_gauginos(S, T, params) + DipI_gauginos(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          brn = born_leptons(S, T, params);
          vir = Virt_leptons(S, T, params) + DipI_leptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        }
#ifndef EXPANSION
        RES += (qE[0][i0] * qE[1 - params->ic][i1] + qE[params->ic][i0] * qE[1][i1]) *
               (brn + vir * aa * 4.0 / 3.0);
#endif
        EXP += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) * ll *
               (-ll * 4.0 / 3.0 + 4.0 - 2.0 * Pqq) * aa * brn;
        EXP += (qN[0][i0] + qN[1 - params->ic][i1] + qN[params->ic][i0] + qN[1][i1]) * gN * (-ll) *
               Pqg * aa * brn;
      }
    }
  }
  // Flux, symmetry factor, spin and color average
  return (RES * PtG(ll, S, params->murs, params->set) - EXP) * dij / 24.0 / S;
}

// Joint resummation
complex<double> JtG(const complex<double> LL, const complex<double> LN, double MIS, double MURS,
                    const int set) {
  // Coeff. of the beta function
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double beta1 = (25.5 - 9.5 * dnf / 3.0) / pow(beta0, 2);

  // Sudakov coeff. for quarks
  const double A1 = 8.0 / 3.0 / beta0;
  const double B1 = -4.0 / beta0;
  const double A2 = .5 * A1 * (67.0 / 3.0 - pow(M_PI, 2) - 10.0 / 9.0 * dnf) / beta0;

  // logarithms
  const complex<double> lmbdN = aS(MURS, set) / (2 * M_PI) * beta0 * LN;
  const complex<double> lmbd = aS(MURS, set) / (2 * M_PI) * beta0 * LL;
  const complex<double> lnlmbd = log(1.0 - lmbd);

  // LL + NLL Sudakov
  const complex<double> g1 = A1 * (1.0 + lnlmbd / lmbd);
  const complex<double> g2 =
      (lmbd / (1.0 - lmbd) * (1.0 - lmbdN) + lnlmbd) * (A1 * log(MIS / MURS) - A2) +
      A1 * beta1 * (.5 * pow(lnlmbd, 2) + (1.0 - lmbdN) / (1.0 - lmbd) * (lmbd + lnlmbd)) +
      B1 * lnlmbd;

  return exp(LL * g1 + g2);
}

complex<double> JtXS(complex<double> N, complex<double> B, double S, double T, Parameters *params) {
  // logarithm
  const double eta = 1.0;
  const complex<double> nb = N * exp(M_EULER);
  const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
  const complex<double> ll = log(bb + nb / (1.0 + bb / nb * eta)) * 2.0;
  const complex<double> ln = log(nb) * 2.0;

  const complex<double> Pqq =
      4.0 / 3.0 * (1.5 + 1.0 / N / (N + 1.0) - 2.0 * (Psi(N + 1.0) + M_EULER));
  const complex<double> Pqg = .5 * (2.0 + N + pow(N, 2)) / N / (N + 1.0) / (N + 2.0);

  // Majorana symmetry factor
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Getting N-space PDFs
  complex<double> gE, qE[2][6], gN, qN[2][6];
  pdfN(gN, qN, N, params->afit);
  pdfN(gE, qE, N, params->afit);

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // Evolving factorization scale [ NLL ]
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double aa = aS(params->murs, params->set) / (2 * M_PI);
  const complex<double> lmbd =
      aa * beta0 * log(params->mufs / S * pow(bb + nb / (1.0 + bb / nb * eta), 2));
  pdfEvolve(gE, qE, N, lmbd);
  CFact(gE, qE, N, aa, aa / (1.0 - lmbd));

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);

  double brn = 0.0;
  double vir = 0.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_squark_gaugino(params->out1, params->out2)) {
          cout << "No joint resummation for associated production (gaugino squark)." << endl;
          exit(1);
        } else if (is_gaugino_gluino(params->out1, params->out2)) {
          cout << "No joint resummation for associated production (gaugino gluino)." << endl;
          exit(1);
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          brn = born_sleptons(S, T, params);
          vir = Virt_sleptons(S, T, params) + DipI_sleptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          brn = born_gauginos(S, T, params);
          vir = Virt_gauginos(S, T, params) + DipI_gauginos(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          brn = born_leptons(S, T, params);
          vir = Virt_leptons(S, T, params) + DipI_leptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        }

#ifndef EXPANSION
        RES += (qE[0][i0] * qE[1 - params->ic][i1] + qE[params->ic][i0] * qE[1][i1]) *
               (brn + vir * aa * 4.0 / 3.0);
#endif
        EXP += (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) * ll *
               (-ll * 4.0 / 3.0 + 4.0 - 2.0 * Pqq) * aa * brn;
        EXP += (qN[0][i0] + qN[1 - params->ic][i1] + qN[params->ic][i0] + qN[1][i1]) * gN * (-ll) *
               Pqg * aa * brn;
      }
    }
  }
  // Flux, symmetry factor, spin and color average
  return (RES * JtG(ll, ln, S, params->murs, params->set) - EXP) * dij / 24.0 / S;
}

complex<double> JtXS2(complex<double> N, complex<double> B, double S, double T,
                      Parameters *params) {
  // logarithm
  const double eta = 1.0;
  const complex<double> nb = N * exp(M_EULER);
  const complex<double> bb = B * sqrt(S) * exp(M_EULER) * .5;
  const complex<double> ll = log(bb + nb / (1.0 + bb / nb * eta)) * 2.0;
  const complex<double> ln = log(nb) * 2.0;

  const complex<double> Pqq =
      4.0 / 3.0 * (1.5 + 1.0 / N / (N + 1.0) - 2.0 * (Psi(N + 1.0) + M_EULER));
  const complex<double> Pqg = .5 * (2.0 + N + pow(N, 2)) / N / (N + 1.0) / (N + 2.0);

  // Majorana symmetry factor
  double dij = 1.0;
  if (params->out1 < 10 && params->out1 == params->out2 && params->out1 / 4 == 0) {
    dij = .5;
  }

  // Getting N-space PDFs
  complex<double> gE, qE[2][6], gN, qN[2][6];
  pdfN(gN, qN, N, params->afit);
  pdfN(gE, qE, N, params->afit);

  // Reset loop tools cache to prevent memory overflow and to get a speed up
  clearcache();

  // Evolving factorization scale [ NLL ]
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;
  const double aa = aS(params->murs, params->set) / (2 * M_PI);
  const complex<double> lmbd =
      aa * beta0 * log(params->mufs / S * pow(bb + nb / (1.0 + bb / nb * eta), 2));
  pdfEvolve(gE, qE, N, lmbd);
  CFact(gE, qE, N, aa, aa / (1.0 - lmbd));

  complex<double> RES(0.0, 0.0);
  complex<double> EXP(0.0, 0.0);

  double brn = 0.0;
  double vir = 0.0;

  // Sum over all possible initial states
  for (int i0 = 0; i0 < 5; i0++) {
    for (int i1 = 0; i1 < 5; i1++) {
      // Set initial state
      params->in1 = i0;
      params->in2 = i1;

      if (is_charge_conserved(params->in1, params->in2, params->out1, params->out2)) {
        if (is_squark_gaugino(params->out1, params->out2)) {
          cout << "No joint resummation for associated production (gaugino squark)." << endl;
          exit(1);
        } else if (is_gaugino_gluino(params->out1, params->out2)) {
          cout << "No joint resummation for associated production (gaugino gluino)." << endl;
          exit(1);
          // not yet implemented
        } else if (is_slepton_slepton(params->out1, params->out2)) {
          brn = born_sleptons(S, T, params);
          vir = Virt_sleptons(S, T, params) + DipI_sleptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_gaugino_gaugino(params->out1, params->out2)) {
          brn = born_gauginos(S, T, params);
          vir = Virt_gauginos(S, T, params) + DipI_gauginos(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        } else if (is_lepton_lepton(params->out1, params->out2)) {
          brn = born_leptons(S, T, params);
          vir = Virt_leptons(S, T, params) + DipI_leptons(S, T, params);
          vir -= (10.0 - pow2(M_PI)) * brn; // to restore old version with old integrated dipole
        }
#ifndef EXPANSION
        RES += (qE[0][i0] * qE[1 - params->ic][i1] + qE[params->ic][i0] * qE[1][i1]) *
               (brn + vir * aa * 4.0 / 3.0);
#endif
        EXP +=
            (qN[0][i0] * qN[1 - params->ic][i1] + qN[params->ic][i0] * qN[1][i1]) *
            (brn + vir * 4.0 / 3.0 * aa +
             (ll * (-ll * 4.0 / 3.0 + 4.0) - 2.0 * Pqq * (ll - log(S / params->mufs))) * aa * brn);
        EXP += (qN[0][i0] + qN[1 - params->ic][i1] + qN[params->ic][i0] + qN[1][i1]) * gN *
               (log(S / params->mufs) - ll) * Pqg * aa * brn;
      }
    }
  }
  // Flux, symmetry factor, spin and color average
  return (RES * JtG(ll, ln, S, params->murs, params->set) - EXP) * dij / 24.0 / S;
}
