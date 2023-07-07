// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Command-line interface.

// File incldudes catani-seymour dipoles
// based on hep-ph/0201036 and hep-ph/9605323.

#include "dipoles.h"
#include "gsl_all.h"
#include "utils.h"
#include <iostream>

#include "constants.h"
#include "debug.h"
#include "tensors.h"

/*
// color factors
#define cA 3.0
#define cF (4.0 / 3.0)
#define I2R (0.5) // SU(3) normalization factor
#define nf 5.0    // number of light flavors

// Flavor constants given in hep-ph/0201036 p.53 eq. C.12 and p.32 eq. 5.91 and
// p.37 eq. 6.17
#define GAMMA_Q (3.0 / 2.0 * cF)
#define GAMMA_GL (3.0 / 2.0 * cA)
#define GAMMA_SQ (2.0 * cF)
#define GAMMA_G (11.0 / 6.0 * cA - 2.0 / 3.0 * I2R * nf)
#define K_Q (7.0 / 2.0 - pow2(M_PI) / 6.0) * cF
#define K_GL (7.0 / 2.0 - pow2(M_PI) / 6.0) * cA
#define K_SQ (4.0 - pow2(M_PI) / 6.0) * cF
#define K_G ((67.0 / 18.0 - pow2(M_PI) / 6.0) * cA - 10.0 / 9.0 * I2R * nf)
*/

// Regulated Altarelli Parisi functions (hep-ph/0201036 p. 31, eq. 5.89).
double P_aap_reg(PartonType a, PartonType ap, double x) {
  if ((a == PARTON_QUARK && ap == PARTON_QUARK) ||
      (a == PARTON_ANTIQUARK && ap == PARTON_ANTIQUARK)) {
    return -cF * (1.0 + x); // P_qq_reg
  }
  if ((a == PARTON_QUARK && ap == PARTON_GLUON) || (a == PARTON_ANTIQUARK && ap == PARTON_GLUON)) {
    return cF * (1.0 + pow2(1.0 - x)) / x; // P_qg_reg = P_qapg_reg
  }
  if ((a == PARTON_GLUON && ap == PARTON_QUARK) || (a == PARTON_GLUON && ap == PARTON_ANTIQUARK)) {
    return I2R * (pow2(x) + pow2(1.0 - x)); // P_gq_reg = P_qqap_reg
  }
  if (a == PARTON_GLUON && ap == PARTON_GLUON) {
    return 2.0 * cA * ((1.0 - x) / x - 1.0 + x * (1.0 - x)); // P_gg_reg
  }
}

// General expression for the Altarelli Parisi splitting function
// (hep-ph/0201036 p. 32, eq. 5.94).
double P_aap(PartonType a, PartonType ap, FunctionType function_type, double x, double z) {

  // for splitting to same flavor a == ap; (Kronecker Delta)
  int delta_aap;
  if (a == ap) {
    delta_aap = 1;
  } else
    delta_aap = 0;

  // Ta * Ta for a == quark or a == gluon (Expectation values of the color
  // operators)
  double color;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    color = cF;
  } else if (a == PARTON_GLUON) {
    color = cA;
  }

  // gamma_a for quark or gluon
  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  } else if (a == PARTON_GLUINO) { // This case never happens, right?
    gamma_a = GAMMA_GL;
  }

  switch (function_type) {
  case FUNCTION_ALL_WDELTA:
    // all parts which are not proportional to delta(1-x)
    return P_aap_reg(a, ap, x) + delta_aap * 2.0 * color * 1.0 / (1.0 - x);
  case FUNCTION_DELTA:
    // part proportional to delta
    // (because we integrate from z to 1 we have to subtract the integral
    // from 0 to z over 2*color*1/(1-x) leading to -2.0 * color * log(1-z) )
    return delta_aap * (gamma_a + 2.0 * color * log(1 - z));
  case FUNCTION_PLUS:
    // the plus distribution part (if there are prefactors
    // which include x and do not belong to the plus-dist. we set them to 1)
    return delta_aap * (2.0 * color * 1.0 / (1.0 - x));
  }
}

// Epsilon dependent parts of d-dimensional AP splitting functions
// (hep-ph/0201036 p. 32 eq 5.93).
double Php_aap(PartonType a, PartonType ap, double x) {
  if ((a == PARTON_QUARK && ap == PARTON_QUARK) ||
      (a == PARTON_ANTIQUARK && ap == PARTON_ANTIQUARK)) {
    return cF * (1 - x);
  } else if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) && ap == PARTON_GLUON) {
    return cF * x;
  } else if (a == PARTON_GLUON && (ap == PARTON_QUARK || ap == PARTON_ANTIQUARK)) {
    return 2.0 * I2R * x * (1.0 - x);
  } else if (a == PARTON_GLUON && ap == PARTON_GLUON) {
    return 0;
  }
}

// General expression for K_ab_bar (hep-ph/0201036 p. 44, eq. 6.56).
double K_aap_bar(PartonType a, PartonType ap, FunctionType function_type, double x, double z) {

  // for splitting to same flavor a == ap;
  int delta_aap;
  if (a == ap) {
    delta_aap = 1.0;
  } else {
    delta_aap = 0.0;
  }

  // Ta * Ta for a == quark or a == gluon
  double color;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    color = cF;
  } else {
    color = cA;
  }

  // gamma_a for quark or gluon
  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  } else if (a == PARTON_GLUINO) {
    gamma_a = GAMMA_GL;
  } else {
    throw "unset gamma_a";
  }

  // K_a for quark or gluon
  double K_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    K_a = K_Q;
  } else if (a == PARTON_GLUON) {
    K_a = K_G;
  } else if (a == PARTON_GLUINO) {
    K_a = K_GL;
  } else {
    throw "unset K_a";
  }

  switch (function_type) {
  case FUNCTION_ALL_WDELTA:
    return P_aap_reg(a, ap, x) * log((1.0 - x) / x) + Php_aap(a, ap, x) +
           delta_aap * color * 2.0 / (1.0 - x) * log((1 - x) / x);
  case FUNCTION_DELTA:
    //-integral cF * 2.0/(1.0 - x) * log((1-x)/x) from 0 to z
    // (solved with mathematica and a substitution 1-x -> y to "help"
    // mathematica)
    return delta_aap * (-color / 3.0 *
                            (pow(M_PI, 2) - // TODO: APN FIX add missing factor 4?!?
                             3. * pow(log(1 - z), 2) - 6. * gsl_sf_dilog(1 - z)) -
                        (gamma_a + K_a - 5.0 / 6.0 * pow2(M_PI) * color)); // explicit delta part
    // TODO: APN FIXED missing factor 4?!?
    // return delta_aap *
    //       (-color / 3.0 *
    //            (4. * pow(M_PI, 2) +
    //             3. * (log(-1 + z) * (-2 * log(-1. + 1. / z) + log(-1 + z) - 2. * log(z))) -
    //             6. * gsl_sf_dilog(1 - z)) -
    //        (gamma_a + K_a - 5.0 / 6.0 * pow2(M_PI) * color)); // explicit delta part
  case FUNCTION_PLUS:
    return delta_aap * color * 2.0 / (1.0 - x) * log((1 - x) / x);
  }
}

// J_a_ij (hep-ph/0201036 p. 27, eq. 5.58)
double J_a_ij(PartonType i, PartonType j, double muj, FunctionType function_type, double x,
              double z) {
  // J_gQ (J_gQb)
  if (i == PARTON_GLUON && (j == PARTON_QUARK || j == PARTON_ANTIQUARK || j == PARTON_GLUINO)) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return (1.0 - x) / (2.0 * pow2(1.0 - x + pow2(muj))) -
             2.0 / (1.0 - x) * (1.0 + log(1.0 - x + pow2(muj))) +
             2.0 / (1.0 - x) * log(2.0 + pow2(muj) - x);
    case FUNCTION_DELTA: // again substitution used in mathematica 1-x -> y to
                         // obtain a "nice" result
      return -(-((pow(muj, 2) * z) / ((1 + pow(muj, 2)) * (1 + pow(muj, 2) - z))) +
               4. * (1 + 2 * log(muj)) * log(1 - z) - 4. * log(1 + pow(muj, 2)) * log(1 - z) -
               log(1 - z / (1. + pow(muj, 2))) + 4. * gsl_sf_dilog(-pow(muj, -2)) -
               4. * gsl_sf_dilog((-1 + z) / pow(muj, 2))) /
             2.;
    case FUNCTION_PLUS:
      return (1.0 - x) / (2.0 * pow2((1.0 - x + pow2(muj)))) -
             2.0 / (1.0 - x) * (1.0 + log(1.0 - x + pow2(muj))) +
             2.0 / (1.0 - x) * log(2.0 + pow2(muj) - 1.0);
      // the last log is not part of the plus distribution, therefore there is a
      // 1 instead of an x
    }
  }
  throw "unimplemented case for J_a_ij";
}

// (hep-ph/0201036 p. 45, eq. 6.57 - 6.60)
// No general expression.
// later u have to add Kappa_ab_sq in C.15 p.54 (TODO: LASSE)
// FUNCTION_ALL_WDELTA contains all terms except the deltas
// FUNCTION_DELTA contains the terms proportional to the delta and also the
// terms in the plus distribution integrated between 0 and z (while their
// coefficients outside of the plus distribution are evaluated for x=1 [APN: and
// get a minus sign!]) FUNCTION_PLUS contains the terms in the plus
// distributions with their coefficients outside of the plus distribution are
// evaluated for x=1
double Kappa_aap_j(PartonType a, PartonType ap, PartonType j, double mj, double sja,
                   FunctionType function_type, double x, double z) {
  int delta_aap;
  if (a == ap) {
    delta_aap = 1.0;
  } else {
    delta_aap = 0.0;
  }

  if (a == PARTON_GLUON && ap == PARTON_QUARK &&
      (j == PARTON_QUARK || j == PARTON_GLUINO ||
       j == PARTON_SQUARK)) { // Kappa_gq_q or Kappa_gq_gl or Kappa_gq_sq
    return 0;
  } else if (a == PARTON_QUARK && ap == PARTON_QUARK &&
             (j == PARTON_QUARK || j == PARTON_GLUINO)) { // Kappa_qq_q or Kappa_qq_gl
                                                          // (6.58)
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return 2.0 * (log(1 - x) / (1 - x) - log(2 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_ALL_WDELTA, x, z) +
             2.0 * 1.0 / (1.0 - x) * log((2.0 - x) * sja / ((2.0 - x) * sja + pow2(mj)));
    case FUNCTION_DELTA:
      return (log(1 - z) * (2 * log(sja / (pow(mj, 2) + sja)) + log(1 - z))) - GAMMA_Q / cF +
             pow2(mj) / sja * log(pow2(mj) / (sja + pow2(mj))) + 0.5 * pow2(mj) / (sja + pow2(mj)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_DELTA, x, z);
    case FUNCTION_PLUS:
      return 2.0 * (log(1 - x) / (1 - x)) +
             J_a_ij(PARTON_GLUON, j, mj / sqrt(sja), FUNCTION_PLUS, x, z) +
             2.0 * 1 / (1 - x) * log(sja / (sja + pow2(mj)));
    }
  } else if (a == PARTON_QUARK && ap == PARTON_GLUON &&
             (j == PARTON_QUARK || j == PARTON_GLUINO)) { // Kappa_qg_q or Kappa_qg_gl
                                                          // (6.59)
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return 2.0 * cF / cA * pow2(mj) / (x * sja) * log(pow2(mj) / ((1 - x) * sja + pow2(mj)));
    case FUNCTION_DELTA:
      return 0.0;
    case FUNCTION_PLUS:
      return 0.0;
    }
  } else if (a == PARTON_GLUON && ap == PARTON_GLUON && (j == PARTON_QUARK || j == PARTON_GLUINO)) {
    // (6.60)
    return Kappa_aap_j(PARTON_QUARK, PARTON_QUARK, PARTON_QUARK, mj, sja, function_type, x, z) +
           cA / cF *
               Kappa_aap_j(PARTON_QUARK, PARTON_GLUON, PARTON_QUARK, mj, sja, function_type, x, z);
  } else if (j == PARTON_SQUARK) {
    // (C.15)
    auto tmp = Kappa_aap_j(a, ap, PARTON_QUARK, mj, sja, function_type, x, z);
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return tmp - delta_aap * ((1 - x) * pow2(sja) / (2.0 * pow2((1 - x) * sja + pow2(mj))));
    case FUNCTION_DELTA:
      return tmp +
             delta_aap *
                 (((-pow2(mj) * z * sja) / (2.0 * (pow2(mj) + sja) * (pow2(mj) + sja * (1 - z)))) +
                  (1.0 / 2.0) * log((pow2(mj) + sja) / (pow2(mj) + sja * (1 - z)))) +
             delta_aap * (-(pow2(mj) / sja * log(pow2(mj) / (sja + pow2(mj))) +
                            0.5 * pow2(mj) / (sja + pow2(mj))) +
                          (GAMMA_Q - GAMMA_SQ) / cF);
    case FUNCTION_PLUS:
      return tmp - delta_aap * ((1 - x) * pow2(sja)) / (2.0 * pow2((1 - x) * sja + pow2(mj)));
    }
  }
  throw "unimplemented Kappa_aap_j";
}

// For the color operators in the formula see hep-ph/9605323 p.97 eq A.4.
// P_aaprime in (hep-ph/0201036 p. 46, eq. 6.67 and 6.53).
double P_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z, double sja, double sab,
                      double mufs, double alphas, EmissionType emission) {

  // some cuttoff to avoid numerical issues
  if ((1.0 - x) < 1.0E-10) {
    return 0.0;
  }

  double color1 = 0.0; // will be Tj * Tap/Tap^2
  double color2 = 0.0; // will be Tb * Tap/Tap^2

  if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_QUARK) {
    color1 = -1.0 / 2.0;
    color2 = -1.0 / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_GLUINO) {
    color1 = -1.0 / cF * cA / 2.0;
    color2 = +1.0 / cF * (cA - 2 * cF) / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_NONE) {
    color1 = 0.0;
    color2 = -1.0;
  } else if (ap == PARTON_QUARK && b == PARTON_GLUON && j == PARTON_SQUARK) {
    color1 = (cA - 2.0 * cF) / (2.0 * cF);
    color2 = -1.0 / cF * cA / 2.0;
  } else if (ap == PARTON_GLUON && b == PARTON_QUARK && j == PARTON_SQUARK) {
    color1 = -1.0 / 2.0;
    color2 = -1.0 / 2.0;
  } else {
    throw "undefined color";
  }

  if (emission == INITIAL_AND_FINAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 *
                 log(mufs / (x * sja)) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 *
                 log(mufs / (x * sab));
    case FUNCTION_DELTA:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 * log(mufs / sja) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 * log(mufs / sab);
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color1 * log(mufs / sja) +
             alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * color2 * log(mufs / sab);
    }

  } else if (emission == INITIAL) {
    // color2 = -1;
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) *
             log(mufs / (x * sab));
    case FUNCTION_DELTA:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * log(mufs / sab);
    case FUNCTION_PLUS:
      return color2 * alphas / (2 * M_PI) * P_aap(a, ap, function_type, x, z) * log(mufs / sab);
    }
  }
}

// p.44 and 46
// eq 6.55 + 6.68
double K_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z, double mj, double sja,
                      double sab, double alphas, EmissionType emission) {

  // cuttoff to avoid numerical instabilities
  if ((1.0 - x) < 1.0E-10) {
    return 0.0;
  }

  double gamma_a = 0.0;
  if (a == PARTON_QUARK || a == PARTON_ANTIQUARK) {
    gamma_a = GAMMA_Q;
  } else if (a == PARTON_GLUON) {
    gamma_a = GAMMA_G;
  }

  int delta_aap;
  if (a == ap) {
    delta_aap = 1.0;
  } else {
    delta_aap = 0.0;
  }

  // Tj Tap/Tap^2, Tj Tap, Tb Tap/Tap^2, Tb Tap
  double color_factor1, color_factor2, color_factor3, color_factor4;
  if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_QUARK) {
    color_factor1 = -1.0 / 2.0;
    color_factor2 = -1.0 / 2.0 * cF;
    color_factor3 = -1.0 / 2.0;
    color_factor4 = -1.0 / 2.0 * cF;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_GLUINO) {
    color_factor1 = -cA / 2.0 * 1.0 / cF;
    color_factor2 = -cA / 2.0;
    color_factor3 = (cA - 2.0 * cF) / 2.0 * 1.0 / cF;
    color_factor4 = (cA - 2.0 * cF) / 2.0;
  } else if (ap == PARTON_QUARK && b == PARTON_ANTIQUARK && j == PARTON_NONE) {
    color_factor1 = 0.0;
    color_factor2 = 0.0;
    color_factor3 = -1.0;
    color_factor4 = -cF;
  } else if (ap == PARTON_QUARK && b == PARTON_GLUON && j == PARTON_SQUARK) {
    color_factor1 = (cA - 2.0 * cF) / (2.0 * cF);
    color_factor2 = (cA / 2.0) - cF;
    color_factor3 = -cA / (2.0 * cF);
    color_factor4 = -cA / 2.0;
  } else if (ap == PARTON_GLUON && b == PARTON_QUARK && j == PARTON_SQUARK) {
    color_factor1 = -1.0 / 2.0;
    color_factor2 = -cA / 2.0;
    color_factor3 = -1.0 / 2.0;
    color_factor4 = -cA / 2.0;
  } else {
    throw "unknown colored particles";
  }

  if (emission == INITIAL_AND_FINAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      // std::cout << "ok2" << std::endl;
      return alphas / (2 * M_PI) *
                 (K_aap_bar(a, ap, function_type, x, z) -
                  color_factor2 * Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z) -
                  color_factor1 *
                      (P_aap_reg(a, ap, x) * log((1 - x) * sja / ((1 - x) * sja + pow2(mj))))) -
             alphas / (2 * M_PI) *
                 (color_factor3 * P_aap_reg(a, ap, x) * log(1 - x) +
                  color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    case FUNCTION_DELTA:
      return -alphas / (2 * M_PI) * color_factor1 * delta_aap * gamma_a *
                 (log((sja - 2 * mj * sqrt(sja + pow2(mj)) + 2 * pow2(mj)) / sja) +
                  2 * mj / (sqrt(sja + pow2(mj)) + mj)) -
             alphas / (2 * M_PI) * color_factor4 * delta_aap *
                 (+pow(log(1 - z), 2) - pow2(M_PI) / 3) +
             alphas / (2 * M_PI) *
                 (K_aap_bar(a, ap, function_type, x, z) -
                  color_factor2 * Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z));
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) *
             (K_aap_bar(a, ap, function_type, x, z) -
              color_factor2 * Kappa_aap_j(a, ap, j, mj, sja, function_type, x, z) -
              color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    }
  }
  if (emission == INITIAL) {
    switch (function_type) {
    case FUNCTION_ALL_WDELTA:
      return alphas / (2 * M_PI) * (K_aap_bar(a, ap, function_type, x, z)) -
             alphas / (2 * M_PI) *
                 (color_factor3 * P_aap_reg(a, ap, x) * log(1 - x) +
                  color_factor4 * delta_aap * 2.0 * log(1 - x) / (1 - x));
    case FUNCTION_DELTA:
      return alphas / (2 * M_PI) * (-color_factor4) * delta_aap *
                 (pow(log(1 - z), 2) - pow2(M_PI) / 3) +
             alphas / (2 * M_PI) * (K_aap_bar(a, ap, function_type, x, z));
    case FUNCTION_PLUS:
      return alphas / (2 * M_PI) *
             (K_aap_bar(a, ap, function_type, x, z) +
              (-color_factor4) * delta_aap * 2.0 * log(1 - x) / (1 - x));
    }
  }
}

double P_plus_K_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                             FunctionType function_type, double x, double z, double mj, double sja,
                             double sab, double mufs, double alphas, EmissionType emission) {
  return P_bold_aap_b_j(a, ap, b, j, function_type, x, z, sja, sab, mufs, alphas, emission) +
         K_bold_aap_b_j(a, ap, b, j, function_type, x, z, mj, sja, sab, alphas, emission);
}
double P_plus_K_bold(PartonType a, PartonType ap, PartonType b, PartonType j, double x, double z,
                     double mj, double sja, double sjac, double sab, double mufs, double alphas,
                     EmissionType emission, double born, double bornc, double djac, double djacplus,
                     double djacdelta) {
  if (a == ap) {
    return (P_plus_K_bold_aap_b_j(a, ap, b, j, FUNCTION_ALL_WDELTA, x, z, mj, sjac, sab, mufs, 1.0,
                                  emission) *
                djac * bornc * 1. / x
            // subtraction term due to plus distribution
            - P_plus_K_bold_aap_b_j(a, ap, b, j, FUNCTION_PLUS, x, z, mj, sja, sab, mufs, 1.0,
                                    emission) *
                  djacplus * born
            // parts proportional to the delta distribution
            + P_plus_K_bold_aap_b_j(a, ap, b, j, FUNCTION_DELTA, x, z, mj, sja, sab, mufs, 1.0,
                                    emission) *
                  djacdelta * born);
  } else {
    return P_plus_K_bold_aap_b_j(a, ap, b, j, FUNCTION_ALL_WDELTA, x, z, mj, sjac, sab, mufs, 1.0,
                                 emission) *
           djac * bornc * 1. / x;
  }
}

// Splitting functions V for initial state emitter ai and specator b
// (V_ai_b) hep-ph/9605323 p.40 eq. 5.145ff
double V_qg_b(double x, double alpha) {
  return 8.0 * M_PI * alpha * (2.0 / (1.0 - x) - (1.0 + x)); // cF
}

// gluon emits quark or antiquark and the counterpart goes into the hard process
double V_gq_b(double x, double alpha) {
  return 8.0 * M_PI * alpha * (1.0 - 2.0 * x * (1.0 - x)); // TR
}

//// this is for a gluon going into the hard process (not needed yet; check and
//// fix)
// double V_qq_b(double x, double alpha, double papi, double papb, double pbpi)
// {
//  return 8.0 * M_PI * alpha *
//         (-4.0 * x + (1 - x) / x * (2.0 * papb / (papi * pbpi))); // cF
//}

// hep-ph/9605323 p.40 eq. 5.147
// This function returns a vector
vector<double> V_qq_b(double x, double alpha, double ztb, double pipb) { // C_F
  vector<double> kernel;
  double zti = 1.0 - ztb;

  kernel.push_back(8.0 * (M_PI)*alpha * (-x));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * ((1.0 - x) / x * zti / ztb * (1.0 / pipb)));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * ((1.0 - x) / x * ztb / zti * (1.0 / pipb)));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * (-(1.0 - x) / x) * (1.0 / pipb));
  return kernel;
}

// hep-ph/9605323 p.40 eq. 5.148
// This function returns a vector
vector<double> V_gg_b(double x, double alpha, double ztb, double pipb) { // C_A
  vector<double> kernel;
  double zti = 1.0 - ztb;

  kernel.push_back(16.0 * (M_PI)*alpha * (-(x / (1.0 - x) + x * (1.0 - x))));
  kernel.push_back(16.0 * (M_PI)*alpha * ((1.0 - x) / x * zti / ztb * (1.0 / pipb)));
  kernel.push_back(16.0 * (M_PI)*alpha * ((1.0 - x) / x * ztb / zti * (1.0 / pipb)));
  kernel.push_back(16.0 * (M_PI)*alpha * (-(1.0 - x) / x) * (1.0 / pipb));
  return kernel;
}

// Splitting functions V for initial state emitter ai and final state spectator
// j (V_ai_j)
// hep-ph/0201036 p. 30 eq 5.81ff
double V_qg_j(double x, double alpha, double zj) {
  return 8.0 * M_PI * alpha * (2.0 / (2.0 - x - zj) - 1.0 - x); // cF
}

double V_gqb_j(double x, double alpha) {
  return 8 * M_PI * alpha * (1.0 - 2.0 * x * (1.0 - x)); // I2R
}

//// or better 5.84 ? (think so, needed for squark prodction) (TODO: LASSE)
//// needed for associated squark gaugino production (fix)
// double V_qq_j(double x, double alpha, double zi, double zj, double pjpi) {
//  return 8 * M_PI * alpha * (-4.0 * x +
//                             (1.0 - x) / x * 2.0 * zi * zj / pjpi *
//                                 (-2.0 * pjpi / (zi * zj))); // cF
//}

// hep-ph/0201036 p. 31 eq 5.83
// This function returns a vector
vector<double> V_qq_j(double x, double alpha, double ztj, double pipj) { // C_F
  vector<double> kernel;
  double zti = 1 - ztj;
  kernel.push_back(8.0 * (M_PI)*alpha * (-x));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * ((1.0 - x) / x * zti / ztj * (1.0 / pipj)));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * ((1.0 - x) / x * ztj / zti * (1.0 / pipj)));
  kernel.push_back(8.0 * (M_PI)*alpha * 2 * (-(1.0 - x) / x * (1.0 / pipj)));
  return kernel;
}

// hep-ph/0201036 p. 31 eq 5.85
// This function returns a vector
vector<double> V_gg_j(double x, double alpha, double ztj, double pipj) { // C_A
  vector<double> kernel;
  double zti = 1 - ztj;
  kernel.push_back(16.0 * (M_PI)*alpha * (-(1.0 / (2.0 - x - ztj) - 1.0 + x * (1.0 - x))));
  kernel.push_back(16.0 * (M_PI)*alpha * ((1.0 - x) / x * zti / ztj * (1.0 / pipj)));
  kernel.push_back(16.0 * (M_PI)*alpha * ((1.0 - x) / x * ztj / zti * (1.0 / pipj)));
  kernel.push_back(16.0 * (M_PI)*alpha * (-(1.0 - x) / x * (1.0 / pipj)));
  return kernel;
}

// Splitting functions V for final state emitter and initial state spectator
// (V_a_ij)
// hep-ph/0201036 p. 26 eq 5.50ff.
double V_a_gQ(double x, double alpha, double zj, double mQ, double pjpi) {
  // (cF for quarks but cA for gluinos)
  return 8. * M_PI * alpha * (2. / (2. - x - zj) - 1. - zj - pow2(mQ) / pjpi);
}

// eq.5.52 expanded in epsilon
double V_a_QQB(double x, double alpha, double zplus, double zminus, double zi) {
  return 8. * M_PI * alpha * (1. - 2. * (zplus - zi) * (zi - zminus)); // TR
}

// hep-ph/0201036 p. 53 eq C.3
double V_a_gSQ(double x, double alpha, double zj, double mQ, double pjpi) {
  return 8.0 * (M_PI)*alpha * (2.0 / (2.0 - x - zj) - 2.0 - pow2(mQ) / pjpi); // cF
}

// Expressions without the color operators
// p. 39 eq. 5.136
// initial emitter and specator
double D_ai_b(PartonType a, PartonType i, double papi, double papb, double pbpi, double x,
              double alpha, double color) {
  if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) && i == PARTON_GLUON) {
    return -1.0 / (2. * papi) * 1.0 / x * color * V_qg_b(x, alpha);
  } else if ((a == PARTON_GLUON && i == PARTON_QUARK) ||
             (a == PARTON_GLUON && i == PARTON_ANTIQUARK)) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_gq_b(x, alpha);
  }
  // else if (a == PARTON_QUARK && i == PARTON_QUARK) {
  //  return -1.0 / (2.0 * papi) * 1.0 / x * color * V_qq_b(x, alpha, papi,
  //  papb, pbpi);
  //}
}

// The vector is given in [g^\mu\nu, p_?^\mu p_?\nu, p_?\mu p_?^\nu, p_X\mu
// p_?^\nu ]
vector<double> D_ai_b_vector(PartonType a, PartonType i, double papi, double x, double alpha,
                             double color, double ztb, double pipb) {
  if (a == PARTON_GLUON && i == PARTON_GLUON) {
    vector<double> Dipole = V_gg_b(x, alpha, ztb, pipb);

    scale(Dipole, -1.0 / (2.0 * papi * x) * 1.0 * color);
    return Dipole;
  }
  if ((a == PARTON_QUARK && i == PARTON_QUARK) ||
      (a == PARTON_ANTIQUARK && i == PARTON_ANTIQUARK)) {
    vector<double> Dipole = V_qq_b(x, alpha, ztb, pipb);
    // std::cout << "DIV? " << ((1-ztb)/ztb)/papi << std::endl;
    scale(Dipole, -1.0 / (2.0 * papi * x) * 1.0 * color);
    return Dipole;
  }
}

// initial emitter and final spectator
// p.29 eq.5.71
double D_ai_j(PartonType a, PartonType i, double zi, double zj, double papi, double pjpi, double x,
              double alpha, double color) {
  if ((a == PARTON_QUARK || a == PARTON_ANTIQUARK) && i == PARTON_GLUON) {
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_qg_j(x, alpha, zj);
  } else if ((a == PARTON_GLUON && i == PARTON_ANTIQUARK) ||
             (a == PARTON_GLUON && i == PARTON_QUARK)) {
    _DEBUG_TABLE("vdip", V_gqb_j(x, alpha));
    return -1.0 / (2.0 * papi) * 1.0 / x * color * V_gqb_j(x, alpha);
  }
  // else if (a == PARTON_QUARK && i == PARTON_QUARK) {
  //  return -1.0 / (2.0 * papi) * 1.0 / x * color * V_qq_j(x, alpha, zi, zj,
  //  pjpi);
  //}
}

vector<double> D_ai_j_vector(PartonType a, PartonType i, double papi, double x, double alpha,
                             double color, double ztj, double pipj) {
  if (a == PARTON_GLUON && i == PARTON_GLUON) {
    vector<double> Dipole = V_gg_j(x, alpha, ztj, pipj);
    scale(Dipole, -1.0 / (2.0 * papi) * 1.0 / x * color);
    return Dipole;
  }
  if ((a == PARTON_ANTIQUARK && i == PARTON_ANTIQUARK) ||
      (a == PARTON_QUARK && i == PARTON_QUARK)) {
    vector<double> Dipole = V_qq_j(x, alpha, ztj, pipj);
    scale(Dipole, -1.0 / (2.0 * papi) * 1.0 / x * color);
    return Dipole;
  }
}

// Final emitter and initial spectator
// p.25 eq 5.40
double D_ij_a(PartonType i, PartonType j, double mi, double mj, double mij, double pjpi, double x,
              double alpha, double zj, double mQ, double zi, double zplus, double zminus,
              double color) {
  if (i == PARTON_GLUON && (j == PARTON_QUARK || j == PARTON_GLUINO)) {
    return -1.0 / (pow2(mi) + pow2(mj) + 2.0 * pjpi - pow2(mij)) * 1.0 / x * color *
           V_a_gQ(x, alpha, zj, mQ, pjpi);
  } else if (i == PARTON_GLUON && j == PARTON_QUARK) {
    // not needed now and maybe wrong (check!) (needed for squark gaugino
    // production)
    return -1.0 / (pow2(mi) + pow2(mj) + 2 * pjpi - pow2(mij)) * 1.0 / x * color *
           V_a_QQB(x, alpha, zplus, zminus, zi);
  } else if (i == PARTON_GLUON && j == PARTON_SQUARK) {
    return -1.0 / (pow2(mi) + pow2(mj) + 2 * pjpi - pow2(mij)) * 1.0 / x * color *
           V_a_gSQ(x, alpha, zj, mQ, pjpi);
  }
  // TODO throw Errors if function end without return is reached!
}

/**
 * Now follow the integrated dipoles I
 */

/**
 * This correpsonds to T_j^2.
 * [(3.6), p.13, arXiv:hep-ph/0201036]
 */
double Casimir_j(PartonType j) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    return cF;
  case PARTON_GLUON:
    return cA;
  case PARTON_SQUARK:
    return cF;
  case PARTON_GLUINO:
    return cA;
  default:
    throw new runtime_error("Unknown Partontype for Casimir_j");
  }
}
/**
 * [(6.28), p.39, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> Gamma_q() {
  // massless case
  return {0., GAMMA_Q, 0., 0.};
}
/**
 * [(6.28), p.39, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> Gamma_q_massive(double mass, Parameters *params) {
  return {0., cF, cF * (0.5 * log(pow2(mass) / params->murs) - 2), 0.};
}
/**
 * [(6.27), p.39, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> Gamma_g(Parameters *params) {
  double tmp = 0;
  for (int i = NFLVR; i < nq; ++i) {
    tmp += 2. / 3. * TR * log(params->mqs[i] / params->murs);
  }
  return {0., GAMMA_G, -tmp, 0.};
}
/**
 * [(C.13), p.54, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> Gamma_sq_massive(double mass, Parameters *params) {
  return {0., cF, cF * (log(pow2(mass) / params->murs) - 2), 0.};
}
/**
 * [(C.11), p.53, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> Gamma_gl_massive(double mass, Parameters *params) {
  return {0., cA, cA * (0.5 * log(pow2(mass) / params->murs) - 2), 0.};
}
Tensor<ComplexType, 4> Gamma_j(PartonType j, double mass, Parameters *params) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    if (mass > 0.)
      return Gamma_q_massive(mass, params);
    else
      return Gamma_q();
  case PARTON_GLUON:
    return Gamma_g(params);
  case PARTON_SQUARK:
    return Gamma_sq_massive(mass, params);
  case PARTON_GLUINO:
    return Gamma_gl_massive(mass, params);
  default:
    throw new runtime_error("Unknown Partontype for Gamma_j");
  }
}
/**
 *  [(5.91), p.32, arXiv:hep-ph/0201036]
 *  [(C.12), p.53, arXiv:hep-ph/0201036]
 */
double gamma_j(PartonType j) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    return GAMMA_Q;
  case PARTON_GLUON:
    return GAMMA_G;
  case PARTON_SQUARK:
    return GAMMA_SQ;
  case PARTON_GLUINO:
    return GAMMA_GL;
  default:
    throw new runtime_error("Unknown Partontype for gamma_j");
  }
}
/**
 *  [(6.17), p.37, arXiv:hep-ph/0201036]
 */
double K_j(PartonType j) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    return K_Q;
  case PARTON_GLUON:
    return K_G;
  case PARTON_SQUARK:
    return K_SQ;
  case PARTON_GLUINO:
    return K_GL;
  default:
    throw new runtime_error("Unknown Partontype for K_j");
  }
}

#define INIT_KERNELS                                                                               \
  auto ms_j = pow2(m_j);                                                                           \
  auto ms_k = pow2(m_k);                                                                           \
  auto Qs_jk = s_jk + ms_j + ms_k;                                                                 \
  auto Q_jk = sqrt(Qs_jk);                                                                         \
  /* [(5.8), p.19, arXiv:hep-ph/0201036]*/                                                         \
  auto v_jk = sqrt(1 - 4. * ms_j * ms_k / pow2(s_jk));                                             \
  /* [(5.30), p.23, arXiv:hep-ph/0201036]*/                                                        \
  auto rho = sqrt((1 - v_jk) / (1 + v_jk));                                                        \
  auto rho_j = sqrt((1 - v_jk + 2 * ms_j / s_jk) / (1 + v_jk + 2 * ms_j / s_jk));                  \
  auto rho_k = sqrt((1 - v_jk + 2 * ms_k / s_jk) / (1 + v_jk + 2 * ms_k / s_jk));

/**
 * Universal for all PartonTypes
 * [(6.20), p.38, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> VS(double s_jk, double m_j, double m_k) {
  Tensor<ComplexType, 4> ret = {0., 0., 0., 0.};
  if (m_j > 0 && m_k > 0) {
    INIT_KERNELS;
    // [(6.20), p.38, arXiv:hep-ph/0201036]
    ret(1) = log(rho) / v_jk;
    ret(2) = (.5 * log(rho) * log(Qs_jk / s_jk) - .25 * pow2(log(pow2(rho_j))) -
              .25 * pow2(log(pow2(rho_k))) - pow2(M_PI) / 6) /
             v_jk;
  } else if (m_j == 0 && m_k == 0) {
    // [(6.20), p.38, arXiv:hep-ph/0201036]
    ret(0) = 1;
    // return {1., 0., 0., 0.};
  } else {
    // [(6.20), p.38, arXiv:hep-ph/0201036]
    // symmetric in j<->k
    if (m_k > 0 && m_j == 0.) {
      m_j = m_k;
      m_k = 0;
    }
    INIT_KERNELS;
    ret(0) = 0.5;
    ret(1) = 0.5 * log(ms_j / s_jk);
    ret(2) = -.25 * pow2(log(ms_j / s_jk)) - pow2(M_PI) / 12 -
             .5 * log(ms_j / s_jk) * log(s_jk / Qs_jk) - .5 * log(ms_j / Qs_jk) * log(s_jk / Qs_jk);
  }
  return ret;
}

/**
 * Gluino and quark case
 * [(6.21)(6.22), p.38, arXiv:hep-ph/0201036]
 */
double VNS_j(PartonType j, double s_jk, double m_j, double m_k) {
  if (j == PARTON_QUARK || j == PARTON_GLUINO || j == PARTON_ANTIQUARK) {
    double T_qT_q = Casimir_j(j);
    INIT_KERNELS;
    if (m_j > 0 && m_k > 0) {
      // [(6.21), p.38, arXiv:hep-ph/0201036]
      return gamma_j(j) / T_qT_q * log(s_jk / Qs_jk) +
             1. / v_jk *
                 (log(pow2(rho)) * log(1 + pow2(rho)) + 2. * gsl_sf_dilog(pow2(rho)) -
                  gsl_sf_dilog(1 - pow2(rho_j)) - gsl_sf_dilog(pow2(rho_k)) - pow2(M_PI) / 6.) +
             log((Q_jk - m_k) / Q_jk) - 2. * log(pow2((Q_jk - m_k) - m_j) / Qs_jk) -
             2. * ms_j / s_jk * log(m_j / (Q_jk - m_k)) - m_k / (Q_jk - m_k) +
             2. * m_k * (2 * m_k - Q_jk) / s_jk + pow2(M_PI) / 2.;
    } else if (m_j == 0 && m_k == 0) {
      // text above  [(6.26), p.39, arXiv:hep-ph/0201036]
      return 0.;
    } else if (m_k == 0) {
      // [(6.22), p.38, arXiv:hep-ph/0201036]
      return gamma_j(j) / T_qT_q * log(s_jk / Qs_jk) + pow2(M_PI) / 6. -
             gsl_sf_dilog(s_jk / Qs_jk) - 2. * log(s_jk / Qs_jk) - ms_j / s_jk * log(ms_j / Qs_jk);
    } else if (m_j == 0) {
      // [(6.23), p.38, arXiv:hep-ph/0201036]
      return gamma_j(j) / T_qT_q *
                 (log(s_jk / Qs_jk) - 2. * log((Q_jk - m_k) / Q_jk) - 2. * m_k / (Q_jk + m_k)) +
             pow2(M_PI) / 6. - gsl_sf_dilog(s_jk / Qs_jk);
    } else {
      throw std::runtime_error("This code (VNS_q) should never be reached!");
    }
  } else {
    throw std::runtime_error("VNS_j is only for quarks and gluinos!");
  }
}
/**
 * [(6.21)(6.22), p.38, arXiv:hep-ph/0201036]
 */
double VNS_q(double s_jk, double m_j, double m_k) { return VNS_j(PARTON_QUARK, s_jk, m_j, m_k); }
/**
 * [(C.8), p.38, arXiv:hep-ph/0201036]
 */
double VNS_gl(double s_jk, double m_j, double m_k) { return VNS_j(PARTON_GLUINO, s_jk, m_j, m_k); }
/**
 * [(6.21)(6.22), p.38, arXiv:hep-ph/0201036]
 */
double VNS_g(double s_jk, double m_j, double m_k, Parameters *params) {
  PartonType j = PARTON_GLUON;
  double T_gT_g = Casimir_j(j);
  INIT_KERNELS;

  if (m_j > 0 && m_k > 0) {
    throw std::runtime_error("This code (VNS_g) should never be reached! (massive Gluon?!?)");
  } else if (m_j == 0 && m_k == 0) {
    //[(6.26), p.39, arXiv:hep-ph/0201036]
    double ret = 0.;
    for (int f = NFLVR; f < nq; ++f) {
      auto m_f = params->mq[f];
      auto ms_f = params->mqs[f];
      auto rho_1 = sqrt(1 - 4. * ms_f / pow2(Q_jk - m_k));
      auto rho_2 = sqrt(1 - 4. * ms_f / (Qs_jk - ms_k));
      if (s_jk > 4. * m_f * (m_f + m_k)) {
        ret +=
            4. / 3. * TR / cA *
            (log((1 + rho_1) / 2.) - rho_1 / 3. * (3. + pow2(rho_1)) - 1. / 2. * log(ms_f / s_jk));
      }
      ret += 2. / 3. * TR / cA * log(ms_f / params->murs);
    }

    return ret;
  } else if (m_k == 0) {
    throw std::runtime_error("This code (VNS_g) should never be reached! (massive Gluon?!?)");

  } else if (m_j == 0) {
    // [(6.24), p.38, arXiv:hep-ph/0201036]
    double ret = gamma_j(j) / T_gT_g *
                     (log(s_jk / Qs_jk) - 2. * log((Q_jk - m_k) / Q_jk) - 2. * m_k / (Q_jk + m_k)) +
                 pow2(M_PI) / 6. - gsl_sf_dilog(s_jk / Qs_jk);
    for (int f = NFLVR; f < nq; ++f) {
      auto m_f = params->mq[f];
      auto ms_f = params->mqs[f];
      auto rho_1 = sqrt(1 - 4. * ms_f / pow2(Q_jk - m_k));
      auto rho_2 = sqrt(1 - 4. * ms_f / (Qs_jk - ms_k));
      if (s_jk > 4. * m_f * (m_f + m_k)) {
        ret +=
            4. / 3. * TR / cA *
            (log((Q_jk - m_k) / Q_jk) + m_k * rho_1 * rho_1 * rho_1 / (Q_jk + m_k) +
             log((1 + rho_1) / 2.) - rho_1 / 3. * (3. + pow2(rho_1)) - 1. / 2. * log(ms_f / Qs_jk));
      }
      ret += 2. / 3. * TR / cA * log(ms_f / params->murs);
    }
    // we choose kappa = 2/3 to neglect the rest of (6.24)
    // ret +=
    return ret;
  } else {
    throw std::runtime_error("This code (VNS_g) should never be reached!");
  }
}
/**
 * [(C.9)(C.10), p.38, arXiv:hep-ph/0201036]
 */
double VNS_sq(double s_jk, double m_j, double m_k) {
  PartonType j = PARTON_SQUARK;
  double T_sqT_sq = Casimir_j(j);
  INIT_KERNELS;
  if (m_j > 0 && m_k > 0) {
    //[(C.9), p.38, arXiv:hep-ph/0201036]
    return gamma_j(j) / T_sqT_sq * log(s_jk / Qs_jk) -
           2. * log((pow2(Q_jk - m_k) - ms_j) / (Qs_jk)) + 4. * m_k * (m_k - Q_jk) / s_jk +
           pow2(M_PI) / 2. +
           1. / v_jk *
               (log(pow2(rho)) * log(1 + pow2(rho)) + 2. * gsl_sf_dilog(pow2(rho)) -
                gsl_sf_dilog(1 - pow2(rho_j)) - gsl_sf_dilog(1 - pow2(rho_k)) - pow2(M_PI) / 6.);
  } else if (m_j == 0 && m_k == 0) {
    throw std::runtime_error("This code (VNS_sq) should never be reached! (massless squark?!?)");
  } else if (m_k == 0) {
    //[(C.10), p.38, arXiv:hep-ph/0201036]
    return gamma_j(j) / T_sqT_sq * log(s_jk / Qs_jk) + pow2(M_PI) / 6. -
           gsl_sf_dilog(s_jk / Qs_jk) - 2. * log(s_jk / Qs_jk);
  } else if (m_j == 0) {
    throw std::runtime_error("This code (VNS_sq) should never be reached! (massless squark?!?)");
  } else {
    throw std::runtime_error("This code (VNS_g) should never be reached!");
  }
}

/**
 * PartonType dependent
 * [(6.20), p.38, arXiv:hep-ph/0201036]
 */
double VNS(PartonType j, double s_jk, double m_j, double m_k, Parameters *params) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    return VNS_q(s_jk, m_j, m_k);
  case PARTON_GLUON:
    return VNS_g(s_jk, m_j, m_k, params);
  case PARTON_SQUARK:
    return VNS_sq(s_jk, m_j, m_k);
  case PARTON_GLUINO:
    return VNS_gl(s_jk, m_j, m_k);
  default:
    throw new runtime_error("Unknown Partontype for VNS");
  }
}
double dr_gamma_j(PartonType j) {
  switch (j) {
  case PARTON_QUARK:
  case PARTON_ANTIQUARK:
    return cF / 2;
  case PARTON_GLUON:
    return cA / 6;
  case PARTON_SQUARK:
    return cF / 2;
  case PARTON_GLUINO:
    return cA / 6;
  default:
    throw new runtime_error("Unknown Partontype for K_j");
  }
}
/**
 * [(6.20), p.38, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> V_j(PartonType j, double s_jk, double m_j, double m_k, Parameters *params) {
  Tensor<ComplexType, 4> finite = {0., 0., 1., 0.};
  return VS(s_jk, m_j, m_k) + finite * (VNS(j, s_jk, m_j, m_k, params) - 0 * dr_gamma_j(j));
}
/**
 * [(6.16), p.37, arXiv:hep-ph/0201036]
 */
Tensor<ComplexType, 4> I_jk(PartonType j, PartonType k, double s_jk, double m_j, double m_k,
                            double T_jT_k, Parameters *params) {
  Tensor<ComplexType, 4> ret = {0, 0, 0, 0};
  Tensor<ComplexType, 4> finite = {0., 0., 1., 0.};
  double Ts_j = Casimir_j(j);
  // eps_mult gives (mu^2/Q^2)^eps expansion
  ret = (
            ///*
            eps_mult({1., log(params->murs / s_jk), 0.5 * pow2(log(params->murs / s_jk))},
                     Ts_j * (V_j(j, s_jk, m_j, m_k, params) - finite * M_PI * M_PI / 3))
            //*/
            ///*
            + Gamma_j(j, m_j, params)
            //*/
            + finite * ((1. + log(params->murs / s_jk)) * gamma_j(j) + K_j(j))
            //*/
            ) *
        pow2(params->g3) / (8 * M_PI * M_PI) * T_jT_k / Ts_j;
  return -ret;
}