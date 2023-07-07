

#include "pxs_gausq_3.h"

#define Power std::pow
//#define Conjugate conj

ComplexType Z_q_massless_wave(POLE pIEPS, int quark, bool right, Parameters *params) {
  BOX_INDEX;
  auto qq = quark;
  // Tensor<ComplexType, 4> r = {0, 0, 0, 0};
  ComplexType r = 0;
  _EPS0(r, pow2(params->g3 / (4 * M_PI)) * cF *
               (2. * B1(0, params->mGLs,
                        params->mSQs[is_up_quark(qq) * 6 + qq + 3 * right - is_up_quark(qq) * 3]) +
                2. * B1(0, 0, 0)));
  return r;
}

ComplexType Z_q_massless_mass(POLE pIEPS, int quark, bool right, Parameters *params) { return 0; }
ComplexType Z_q_massive_wave(POLE pIEPS, int quark, bool right, Parameters *params) {
  throw new runtime_error("unimplemented");
  return 0;
}
ComplexType Z_q_massive_mass(POLE pIEPS, int quark, bool right, Parameters *params) {
  throw new runtime_error("unimplemented");
  return 0;
}

ComplexType Z_g(POLE pIEPS, Parameters *params) {
  ComplexType r = 0;

  // heavy q
  //*
  for (int itq = nf; itq < 6; itq++) {
    auto ms_q = params->mqs[itq];
    _EPS0(r, -2. / 3. * (2. * ms_q * DB0(0, ms_q, ms_q) + B0(0, ms_q, ms_q)));
    _FIN(r, 2. / 3. * 1. / 3.);
  }
  //*/
  //*
  // light q
  for (int itq = 0; itq < nf; itq++) {
    _EPS0(r, -2. / 3. * (B0(0, 0, 0)));
  }
  // gluons + ghosts
  _EPS0(r, (19. / 4. + 1 / 4.) * B0(0, 0, 0));
  //*/
  //*
  // gluino
  auto ms_G = params->mGLs;
  _EPS0(r, -2. / 3. * NC * (B0(0, ms_G, ms_G) + 2 * ms_G * DB0(0, ms_G, ms_G)));
  _FIN(r, 2. / 3. * 1. / 3.);
  // squarks
  for (int itsq = 0; itsq < 12; itsq++) {
    auto ms_Q = params->mSQs[itsq];
    _EPS0(r, -TR / 3. * (B0(0, ms_Q, ms_Q) - 4 * ms_Q * DB0(0, ms_Q, ms_Q)));
    _FIN(r, -TR / 3. * 2. / 3.);
  }
  //*/

  return pow2(params->g3 / (4 * M_PI)) * r;
}

ComplexType Z_Q_massless_mass(POLE pIEPS, int squark, Parameters *params) {
  auto sq = squark;
  auto ms_Q = params->mSQs[squark];
  auto ms_G = params->mGLs;
  ComplexType r = 0;
  _EPS0(r, -pow2(params->g3 / (4 * M_PI)) * 1. / ms_Q * cF *
               (4 * ms_Q * B0(ms_Q, 0, ms_Q) - A0(ms_Q) - A0(ms_Q) + 2. * A0(ms_G) +
                2 * (ms_G - ms_Q) * B0(ms_Q, ms_G, 0)));
  return r;
}
ComplexType Z_Q_massless_wave(POLE pIEPS, int squark, Parameters *params) {
  auto sq = squark;
  auto ms_Q = params->mSQs[squark];
  auto ms_G = params->mGLs;
  ComplexType r = 0;
  _EPS0(r, pow2(params->g3 / (4 * M_PI)) * cF *
               (2. * B0(ms_Q, 0, ms_Q) + 4 * ms_Q * DB0(ms_Q, 0, ms_Q) +
                2 * (ms_G - ms_Q) * DB0(ms_Q, ms_G, 0) - 2. * B0(ms_Q, ms_G, 0)));
  return r;
}

ComplexType d_alpha(POLE pIEPS, Parameters *params) {
  ComplexType r = 0;
  auto a_s = pow2(params->g3) / 4. / M_PI;
  double log_squarks = 0, log_quarks = 0;
  for (int i = 0; i < 12; i++) {
    log_squarks += log(params->mSQs[i] / params->murs);
  }
  for (int i = 0; i < 6; i++) {
    if (params->mqs[i] > 0) {
      log_quarks += log(params->mqs[i] / params->murs);
    }
  }
// finite
#ifdef DECOUPLE_HEAVY
  if (pIEPS == FIN) {
    r += (-a_s / (6. * M_PI) * log_quarks -
          a_s / (6. * M_PI) * NC * log(params->mGLs / params->murs) -
          a_s / (24. * M_PI) * log_squarks);
  }
#endif
  // UV
  if (pIEPS == UV || pIEPS == IR_UV) {
    r += a_s / (2 * M_PI) * (nf / 3. - (11 - 19. / 4. + 19. / 4.) * NC / 6) +
         a_s / (6 * M_PI) * 1. + a_s / (6 * M_PI) * NC * 1. + a_s / (24 * M_PI) * 12 * 1.;
  }
  // from squark vertex renorm
  // r(3) -= a_s / (4 * M_PI) * (NC * 2. + 0 * NC * NC / 2);
  return 0.5 * r;
}