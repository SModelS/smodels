#include "pxs_gausq_3.h"

#include "pxs_gausq_2.h"

#define Power std::pow
//#define Conjugate conj

#define SCHEME(v) scheme(pIEPS, v)
ComplexType scheme(POLE pIEPS, ComplexType val) {
  // UVBAR like in rothering_phd (3.41)
#ifdef MSBAR
  if (pIEPS == UV) {
    return val;
  } else {
    return 0;
  }
#endif
  // OS-scheme
  return val;
}

/**
 * General counter terms can be found in e.g. 1907.04898
 */

ComplexType ME_us_se_qq_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;

  ComplexType r = 0;

  r = +2.0 * (
                 // left handed
                 ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 0, 1, 1, 1, 1).real() *
                     SCHEME(Z_q_massless_wave(pIEPS, q, 0, params)) +
                 // right handed
                 ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 0, 1, 1, 1, 1, 1).real() *
                     SCHEME(Z_q_massless_wave(pIEPS, q, 1, params)));
#ifdef MSBAR
  // Residue see Rothering phd (3.42) for quarks
  if (pIEPS == FIN) {
    // r += -1. * ME_us_born(1, 1, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 1, 1).real() *
    //     (Z_g(pIEPS, params));
    r += -2. * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 0, 1).real() *
         (Z_q_massless_wave(pIEPS, q, 0, params));
    r += -2. * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 1, 0).real() *
         (Z_q_massless_wave(pIEPS, q, 1, params));
    r += -2. * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 1, 0).real() *
         (Z_q_massless_wave(pIEPS, q, 0, params));
    r += -2. * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 0, 1).real() *
         (Z_q_massless_wave(pIEPS, q, 1, params));
  }
#endif
  // r(1) = 0.0;
  return -real(r);
}

ComplexType ME_us_se_QQ_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  ComplexType r = 0;

  ///*
  r += 2.0 * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params).real() *
       (Z_Q_massless_wave(pIEPS, sq, params));

  ///*
  r -= (2.0 * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params).real() * params->mSQs[sq] *
        (Z_Q_massless_mass(pIEPS, sq, params)) / ((MUs + MXs - s - t) - MUs));
  //*/

  // r(1) = 0.0;
  return -real(r);
}
ComplexType ME_us_qqg_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  ComplexType r = 0;
  auto dZa = d_alpha(pIEPS, params);
  auto dZg = (Z_g(pIEPS, params));
  auto dZql = SCHEME(Z_q_massless_wave(pIEPS, q, 0, params));
  auto dZqr = SCHEME(Z_q_massless_wave(pIEPS, q, 1, params));
  auto dZ = (dZa * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real() +
             0.5 * dZg * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real() +
             dZql * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 0, 1, 1).real() +
             dZqr * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 0, 1, 1, 1).real());
  r += dZ * 2.0;
  // shift wrong IR to UV
  if (pIEPS == IR1 || pIEPS == IR_UV) {
    r -= 2. * Cf * 2.0 * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real() *
         (pow2(params->g3 / (4 * M_PI)));
  }
  if (pIEPS == UV || pIEPS == IR_UV) {
    r += 2. * Cf * 2.0 * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real() *
         (pow2(params->g3 / (4 * M_PI)));
  }
  return r;
}
ComplexType ME_us_QQg_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  ComplexType r = 0;
  auto dZa = d_alpha(pIEPS, params);
  auto dZg = (Z_g(pIEPS, params));
  auto dZQ = Z_Q_massless_wave(pIEPS, sq, params);
  auto dZ = (dZa + 0.5 * dZg + dZQ);
  r += dZ * 2.0 * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params).real();
  return r;
}
ComplexType ME_ss_qQX_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  ComplexType r = 0;
  auto dZQ = Z_Q_massless_wave(pIEPS, sq, params);
  auto dZql = SCHEME(Z_q_massless_wave(pIEPS, q, 0, params));
  auto dZqr = SCHEME(Z_q_massless_wave(pIEPS, q, 1, params));
  auto dZ =
      0.5 * (dZql * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 1, 0).real() +
             dZqr * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 0, 1).real() +
             dZQ * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real());
#ifdef DRBAR
  if (pIEPS == FIN) {
    auto a_s = pow2(params->g3) / 4. / M_PI;
    // shift gaugino-quark-squark vertex from CDR to DRBAR
    // (see arxiv:hep-ph/0511344v2 & arXiv:hep-ph/9610490).
    dZ -= a_s / (8 * M_PI) * cF * ME_us_born(1, 0, sc, uc, axial, Q2, P1K1, params).real();
  }
#endif
  r += dZ * 2.;
  return r;
}
ComplexType ME_uu_qQX_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  ComplexType r = 0;
  auto dZQ = Z_Q_massless_wave(pIEPS, sq, params);
  auto dZql = SCHEME(Z_q_massless_wave(pIEPS, q, 0, params));
  auto dZqr = SCHEME(Z_q_massless_wave(pIEPS, q, 1, params));
  // auto dZ = 0.5 * (dZq + dZQ);
  auto dZ =
      0.5 * (dZql * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 1, 0).real() +
             dZqr * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params, 1, 1, 1, 1, 0, 1).real() +
             dZQ * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params).real());
#ifdef DRBAR
  if (pIEPS == FIN) {
    auto a_s = pow2(params->g3) / 4. / M_PI;
    // shift gaugino-quark-squark vertex from CDR to DRBAR
    // (see arxiv:hep-ph/0511344v2 & arXiv:hep-ph/9610490).
    dZ -= a_s / (8 * M_PI) * cF * ME_us_born(0, 1, sc, uc, axial, Q2, P1K1, params).real();
  }
#endif
  r += dZ * 2.;
  return r;
}
