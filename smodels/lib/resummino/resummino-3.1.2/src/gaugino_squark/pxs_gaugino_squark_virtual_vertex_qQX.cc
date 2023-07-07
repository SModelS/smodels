
#include "pxs_gausq_3.h"

#define C0(a, b, c, d, e, f) C0i(cc0, a, b, c, d, e, f)
#define C1(a, b, c, d, e, f) C0i(cc1, a, b, c, d, e, f)
#define C2(a, b, c, d, e, f) C0i(cc2, a, b, c, d, e, f)
#define C00(a, b, c, d, e, f) C0i(cc00, a, b, c, d, e, f)
#define C11(a, b, c, d, e, f) C0i(cc11, a, b, c, d, e, f)
#define C12(a, b, c, d, e, f) C0i(cc12, a, b, c, d, e, f)
#define C22(a, b, c, d, e, f) C0i(cc22, a, b, c, d, e, f)

#define Power std::pow
//#define Conjugate conj

ComplexType ME_uu_qQX_qgQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  int q_I = q;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  double SS = sc, UU = uc, AXG = axial;

  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 B0(MUs + MXs - s - t, 0, MQs)
#define syFC2 C1(MXs, MUs + MXs - s - t, 0, 0, MQs, 0)
#define syFC3 C2(MXs, MUs + MXs - s - t, 0, 0, MQs, 0)

  ComplexType ret = 0;
  _EPS0(
      ret,
      (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) * (+R * Lp + L * Rp) *
          (+Denom(-768 * pow(Pi, 2) * s * MXs * MQs * Nc + 768 * pow(Pi, 2) * s * MXs * MUs * Nc +
                  768 * pow(Pi, 2) * s * pow(MXs, 2) * Nc + 768 * pow(Pi, 2) * s * t * MQs * Nc -
                  768 * pow(Pi, 2) * s * t * MUs * Nc - 1536 * pow(Pi, 2) * s * t * MXs * Nc +
                  768 * pow(Pi, 2) * s * pow(t, 2) * Nc + 768 * pow(Pi, 2) * pow(s, 2) * MQs * Nc -
                  768 * pow(Pi, 2) * pow(s, 2) * MUs * Nc -
                  1536 * pow(Pi, 2) * pow(s, 2) * MXs * Nc +
                  1536 * pow(Pi, 2) * pow(s, 2) * t * Nc + 768 * pow(Pi, 2) * pow(s, 3) * Nc)) *
          (+syFC1 * (1.) + syFC2 * (2 * MXs) + syFC3 * (2 * MXs - 2 * u)) * (+gs) * (+gs) * (+gs) *
          (+gs) *
          (+2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS - 4 * pow(MX, 4) * MUs * SS +
           2 * pow(MX, 4) * MUs * SS * AXG - 4 * pow(MX, 4) * MUs * UU * AXG -
           4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG - 2 * u * pow(MUs, 2) * SS -
           4 * u * MXs * MUs * SS * AXG + 8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
           8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
           2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
           2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
           4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
           4 * pow(u, 2) * MXs * SS * AXG - 8 * pow(u, 2) * MXs * UU * AXG -
           2 * pow(u, 3) * SS * AXG + 4 * pow(u, 3) * UU * AXG + s * MXs * MUs * SS -
           s * MXs * MUs * SS * AXG - 2 * s * MXs * MUs * UU + 2 * s * MXs * MUs * UU * AXG -
           2 * s * pow(MXs, 2) * UU + 2 * s * pow(MXs, 2) * UU * AXG + 2 * s * pow(MX, 4) * UU -
           2 * s * pow(MX, 4) * UU * AXG + s * u * MUs * SS + s * u * MUs * SS * AXG +
           2 * s * u * MUs * UU - 2 * s * u * MUs * UU * AXG - s * u * MXs * SS +
           s * u * MXs * SS * AXG - 2 * s * u * MXs * UU - 2 * s * u * MXs * UU * AXG -
           s * pow(u, 2) * SS - s * pow(u, 2) * SS * AXG + 2 * s * pow(u, 2) * UU +
           2 * s * pow(u, 2) * UU * AXG));
  ;
  return ret.real();
}

ComplexType ME_ss_qQX_qgQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  int q_I = q;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  double SS = sc, UU = uc, AXG = axial;

  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 B0(MXs, 0, MQs)
#define syFC2 C0(MUs, s, MXs, MQs, 0, 0)
#define syFC3 C1(MUs, s, MXs, MQs, 0, 0)
#define syFC4 C2(MUs, s, MXs, MQs, 0, 0)
  ComplexType ret = 0;
  _EPS0(
      ret,
      (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
          (+Denom(-768 * pow(Pi, 2) * s * MXs * Nc + 768 * pow(Pi, 2) * s * t * Nc +
                  768 * pow(Pi, 2) * pow(s, 2) * Nc)) *
          (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 *
               (-2 * pow(MUs, 2) * SS + 2 * MXs * MUs * UU + 2 * pow(MXs, 2) * UU -
                4 * pow(MX, 4) * UU + 2 * pow(MX, 4) * UU * AXG + 4 * u * MUs * SS -
                2 * u * MUs * UU + 2 * u * MXs * UU - 4 * u * MXs * UU * AXG - 2 * pow(u, 2) * SS +
                2 * pow(u, 2) * UU * AXG + 2 * s * MUs * SS + s * MXs * UU - s * MXs * UU * AXG -
                2 * s * u * SS + s * u * UU + s * u * UU * AXG) +
           syFC2 *
               (-4 * pow(MUs, 3) * SS + 2 * MXs * pow(MUs, 2) * SS - 2 * MXs * pow(MUs, 2) * UU +
                6 * pow(MXs, 2) * MUs * SS + 2 * pow(MXs, 2) * MUs * UU + 4 * pow(MXs, 3) * UU -
                6 * pow(MX, 4) * MUs * SS - 2 * pow(MX, 4) * MUs * UU -
                2 * pow(MX, 4) * MUs * UU * AXG - 14 * pow(MX, 4) * MXs * UU +
                4 * pow(MX, 4) * MXs * UU * AXG + 12 * pow(MX, 6) * UU - 6 * pow(MX, 6) * UU * AXG +
                6 * pow(MU, 4) * MUs * SS + 2 * u * pow(MUs, 2) * SS + 2 * u * pow(MUs, 2) * UU -
                4 * u * MXs * MUs * SS + 4 * u * MXs * MUs * UU * AXG - 6 * u * pow(MXs, 2) * SS +
                10 * u * pow(MXs, 2) * UU - 8 * u * pow(MXs, 2) * UU * AXG +
                6 * u * pow(MX, 4) * SS - 12 * u * pow(MX, 4) * UU +
                12 * u * pow(MX, 4) * UU * AXG - 6 * u * pow(MU, 4) * SS +
                2 * pow(u, 2) * MUs * SS - 2 * pow(u, 2) * MUs * UU * AXG +
                2 * pow(u, 2) * MXs * SS - 2 * pow(u, 2) * MXs * UU * AXG +
                4 * s * pow(MUs, 2) * SS - s * MXs * MUs * UU + s * MXs * MUs * UU * AXG +
                2 * s * pow(MXs, 2) * UU - 2 * s * pow(MXs, 2) * UU * AXG -
                5 * s * pow(MX, 4) * UU + 5 * s * pow(MX, 4) * UU * AXG - 6 * s * pow(MU, 4) * SS +
                2 * s * u * MUs * SS - s * u * MUs * UU - s * u * MUs * UU * AXG +
                s * u * MXs * UU - 3 * s * u * MXs * UU * AXG + pow(s, 2) * MXs * UU -
                pow(s, 2) * MXs * UU * AXG) +
           syFC3 *
               (-2 * MXs * pow(MUs, 2) * SS - 6 * MXs * pow(MUs, 2) * UU -
                2 * pow(MXs, 2) * MUs * SS - 6 * pow(MXs, 2) * MUs * UU +
                2 * pow(MX, 4) * MUs * SS + 14 * pow(MX, 4) * MUs * UU -
                6 * pow(MX, 4) * MUs * UU * AXG + 2 * pow(MX, 4) * MXs * UU - 4 * pow(MX, 6) * UU +
                2 * pow(MX, 6) * UU * AXG + 6 * pow(MU, 4) * MUs * SS - 6 * u * pow(MUs, 2) * SS +
                6 * u * pow(MUs, 2) * UU + 4 * u * MXs * MUs * SS - 8 * u * MXs * MUs * UU +
                12 * u * MXs * MUs * UU * AXG + 2 * u * pow(MXs, 2) * SS -
                2 * u * pow(MXs, 2) * UU - 2 * u * pow(MX, 4) * SS + 4 * u * pow(MX, 4) * UU -
                4 * u * pow(MX, 4) * UU * AXG - 6 * u * pow(MU, 4) * SS + 6 * pow(u, 2) * MUs * SS -
                6 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                2 * pow(u, 2) * MXs * UU * AXG - 3 * s * MXs * MUs * UU +
                3 * s * MXs * MUs * UU * AXG + 3 * s * pow(MX, 4) * UU -
                3 * s * pow(MX, 4) * UU * AXG - 6 * s * pow(MU, 4) * SS + 6 * s * u * MUs * SS -
                3 * s * u * MUs * UU - 3 * s * u * MUs * UU * AXG - s * u * MXs * UU +
                3 * s * u * MXs * UU * AXG - pow(s, 2) * MXs * UU + pow(s, 2) * MXs * UU * AXG) +
           syFC4 *
               (-4 * pow(MUs, 3) * SS + 2 * MXs * pow(MUs, 2) * SS - 2 * MXs * pow(MUs, 2) * UU +
                6 * pow(MXs, 2) * MUs * SS + 2 * pow(MXs, 2) * MUs * UU + 4 * pow(MXs, 3) * UU -
                6 * pow(MX, 4) * MUs * SS - 2 * pow(MX, 4) * MUs * UU -
                2 * pow(MX, 4) * MUs * UU * AXG - 14 * pow(MX, 4) * MXs * UU +
                4 * pow(MX, 4) * MXs * UU * AXG + 12 * pow(MX, 6) * UU - 6 * pow(MX, 6) * UU * AXG +
                6 * pow(MU, 4) * MUs * SS + 2 * u * pow(MUs, 2) * SS + 2 * u * pow(MUs, 2) * UU -
                4 * u * MXs * MUs * SS + 4 * u * MXs * MUs * UU * AXG - 6 * u * pow(MXs, 2) * SS +
                10 * u * pow(MXs, 2) * UU - 8 * u * pow(MXs, 2) * UU * AXG +
                6 * u * pow(MX, 4) * SS - 12 * u * pow(MX, 4) * UU +
                12 * u * pow(MX, 4) * UU * AXG - 6 * u * pow(MU, 4) * SS +
                2 * pow(u, 2) * MUs * SS - 2 * pow(u, 2) * MUs * UU * AXG +
                2 * pow(u, 2) * MXs * SS - 2 * pow(u, 2) * MXs * UU * AXG +
                2 * s * pow(MUs, 2) * SS + s * MXs * MUs * UU + s * MXs * MUs * UU * AXG +
                4 * s * pow(MXs, 2) * UU - 2 * s * pow(MXs, 2) * UU * AXG -
                9 * s * pow(MX, 4) * UU + 7 * s * pow(MX, 4) * UU * AXG - 6 * s * pow(MU, 4) * SS +
                6 * s * u * MUs * SS - 3 * s * u * MUs * UU - s * u * MUs * UU * AXG +
                3 * s * u * MXs * UU - 7 * s * u * MXs * UU * AXG - 2 * s * pow(u, 2) * SS +
                2 * s * pow(u, 2) * UU * AXG + 2 * pow(s, 2) * MUs * SS + 2 * pow(s, 2) * MXs * UU -
                2 * pow(s, 2) * MXs * UU * AXG - 2 * pow(s, 2) * u * SS + pow(s, 2) * u * UU +
                pow(s, 2) * u * UU * AXG)));
  return ret.real();
}

ComplexType ME_uu_qQX_QGq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  auto Denom = [](auto a) { return 1. / a; };

  ComplexType ret = 0;
  int itq = q;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int ftq = 0; ftq < 2; ftq++) {
    int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
    // int iq = (itq + ftq * 3) % 6;
    int iq = (sq - is_up_squark(sq) * 6) % 3 + is_up_squark(sq) * 3;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType kL = (params->CHSQq[ch][isq][iq].L);
    ComplexType kR = (params->CHSQq[ch][isq][iq].R);

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    ComplexType iLGp = conj(params->GLSQq[sq][iq].R);
    ComplexType iRGp = conj(params->GLSQq[sq][iq].L);
    ComplexType jLGp = conj(params->GLSQq[isq][q].R);
    ComplexType jRGp = conj(params->GLSQq[isq][q].L);

#define syFC1 B0(MUs + MXs - s - t, MGs, Mqis)
#define syFC2 C0(MXs, 0, MUs + MXs - s - t, Mqis, MQis, MGs)
#define syFC3 C1(MXs, MUs + MXs - s - t, 0, MQis, Mqis, MGs)

    _EPS0(
        ret,
        -(-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
            (+Denom(
                -768 * pow(Pi, 2) * s * MXs * MQs * Nc + 768 * pow(Pi, 2) * s * MXs * MUs * Nc +
                768 * pow(Pi, 2) * s * pow(MXs, 2) * Nc + 768 * pow(Pi, 2) * s * t * MQs * Nc -
                768 * pow(Pi, 2) * s * t * MUs * Nc - 1536 * pow(Pi, 2) * s * t * MXs * Nc +
                768 * pow(Pi, 2) * s * pow(t, 2) * Nc + 768 * pow(Pi, 2) * pow(s, 2) * MQs * Nc -
                768 * pow(Pi, 2) * pow(s, 2) * MUs * Nc - 1536 * pow(Pi, 2) * pow(s, 2) * MXs * Nc +
                1536 * pow(Pi, 2) * pow(s, 2) * t * Nc + 768 * pow(Pi, 2) * pow(s, 3) * Nc)) *
            (+gs) * (+gs) *
            (+2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS - 4 * pow(MX, 4) * MUs * SS +
             2 * pow(MX, 4) * MUs * SS * AXG - 4 * pow(MX, 4) * MUs * UU * AXG -
             4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
             2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
             8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
             8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
             2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
             2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
             4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
             4 * pow(u, 2) * MXs * SS * AXG - 8 * pow(u, 2) * MXs * UU * AXG -
             2 * pow(u, 3) * SS * AXG + 4 * pow(u, 3) * UU * AXG + s * MXs * MUs * SS -
             s * MXs * MUs * SS * AXG - 2 * s * MXs * MUs * UU + 2 * s * MXs * MUs * UU * AXG -
             2 * s * pow(MXs, 2) * UU + 2 * s * pow(MXs, 2) * UU * AXG + 2 * s * pow(MX, 4) * UU -
             2 * s * pow(MX, 4) * UU * AXG + s * u * MUs * SS + s * u * MUs * SS * AXG +
             2 * s * u * MUs * UU - 2 * s * u * MUs * UU * AXG - s * u * MXs * SS +
             s * u * MXs * SS * AXG - 2 * s * u * MXs * UU - 2 * s * u * MXs * UU * AXG -
             s * pow(u, 2) * SS - s * pow(u, 2) * SS * AXG + 2 * s * pow(u, 2) * UU +
             2 * s * pow(u, 2) * UU * AXG) *
            (+syFC1 * (kR * L * iLGp * jRGp + kL * R * iRGp * jLGp) +
             syFC2 * (kR * MQis * L * iLGp * jRGp + kR * MG * Mqi * L * iRGp * jRGp +
                      kR * MG * MX * R * iLGp * jLGp + kL * MQis * R * iRGp * jLGp +
                      kL * MG * Mqi * R * iLGp * jLGp + kL * MG * MX * L * iRGp * jRGp) +
             syFC3 * (kR * MX * Mqi * R * iRGp * jLGp + kR * pow(MX, 2) * L * iLGp * jRGp +
                      kR * MG * MX * R * iLGp * jLGp + kL * MX * Mqi * L * iLGp * jRGp +
                      kL * pow(MX, 2) * R * iRGp * jLGp + kL * MG * MX * L * iRGp * jRGp)));
    //}
  }
  return ret.real();
}

ComplexType ME_ss_qQX_QGq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  auto Denom = [](auto a) { return 1. / a; };
  ComplexType ret = 0;
  int itq = q;
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int ftq = 0; ftq < 2; ftq++) {
    int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
    // int iq = (itq + ftq * 3) % 6;
    int iq = (sq - is_up_squark(sq) * 6) % 3 + is_up_squark(sq) * 3;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType kL = (params->CHSQq[ch][isq][iq].L);
    ComplexType kR = (params->CHSQq[ch][isq][iq].R);

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    ComplexType iLGp = conj(params->GLSQq[sq][iq].R);
    ComplexType iRGp = conj(params->GLSQq[sq][iq].L);
    ComplexType jLGp = conj(params->GLSQq[isq][q].R);
    ComplexType jRGp = conj(params->GLSQq[isq][q].L);

#define syFC1 B0(MXs, Mqis, MQis)
#define syFC2 C0(MUs, s, MXs, Mqis, MGs, MQis)
#define syFC3 C1(MUs, s, MXs, Mqis, MGs, MQis)
#define syFC4 C2(MUs, s, MXs, Mqis, MGs, MQis)

    _EPS0(ret,
          -(-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
              (+Denom(-768 * pow(Pi, 2) * s * MXs * Nc + 768 * pow(Pi, 2) * s * t * Nc +
                      768 * pow(Pi, 2) * pow(s, 2) * Nc)) *
              (+gs) * (+gs) *
              (+syFC1 *
                   (-2. * kL * pow(MUs, 2) * R * iRGp * jLGp * SS +
                    2. * kL * MXs * MUs * R * iRGp * jLGp * UU +
                    2. * kL * pow(MXs, 2) * R * iRGp * jLGp * UU -
                    4. * kL * pow(MX, 4) * R * iRGp * jLGp * UU +
                    2. * kL * pow(MX, 4) * R * iRGp * jLGp * UU * AXG +
                    4. * kL * u * MUs * R * iRGp * jLGp * SS -
                    2. * kL * u * MUs * R * iRGp * jLGp * UU +
                    2. * kL * u * MXs * R * iRGp * jLGp * UU -
                    4. * kL * u * MXs * R * iRGp * jLGp * UU * AXG -
                    2. * kL * pow(u, 2) * R * iRGp * jLGp * SS +
                    2. * kL * pow(u, 2) * R * iRGp * jLGp * UU * AXG +
                    2. * kL * s * MUs * R * iRGp * jLGp * SS + kL * s * MXs * R * iRGp * jLGp * UU -
                    kL * s * MXs * R * iRGp * jLGp * UU * AXG -
                    2. * kL * s * u * R * iRGp * jLGp * SS + kL * s * u * R * iRGp * jLGp * UU +
                    kL * s * u * R * iRGp * jLGp * UU * AXG -
                    2. * kR * pow(MUs, 2) * L * iLGp * jRGp * SS +
                    2. * kR * MXs * MUs * L * iLGp * jRGp * UU +
                    2. * kR * pow(MXs, 2) * L * iLGp * jRGp * UU -
                    4. * kR * pow(MX, 4) * L * iLGp * jRGp * UU +
                    2. * kR * pow(MX, 4) * L * iLGp * jRGp * UU * AXG +
                    4. * kR * u * MUs * L * iLGp * jRGp * SS -
                    2. * kR * u * MUs * L * iLGp * jRGp * UU +
                    2. * kR * u * MXs * L * iLGp * jRGp * UU -
                    4. * kR * u * MXs * L * iLGp * jRGp * UU * AXG -
                    2. * kR * pow(u, 2) * L * iLGp * jRGp * SS +
                    2. * kR * pow(u, 2) * L * iLGp * jRGp * UU * AXG +
                    2. * kR * s * MUs * L * iLGp * jRGp * SS + kR * s * MXs * L * iLGp * jRGp * UU -
                    kR * s * MXs * L * iLGp * jRGp * UU * AXG -
                    2. * kR * s * u * L * iLGp * jRGp * SS + kR * s * u * L * iLGp * jRGp * UU) +
               syFC1 * (kR * s * u * L * iLGp * jRGp * UU * AXG) +
               syFC2 * (2. * kL * pow(MUs, 3) * R * iRGp * jLGp * SS -
                        2. * kL * MXs * pow(MUs, 2) * R * iRGp * jLGp * UU -
                        2. * kL * pow(MXs, 2) * MUs * R * iRGp * jLGp * UU +
                        2. * kL * MX * pow(MUs, 2) * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * MX * MXs * MUs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * MX * pow(MXs, 2) * Mqi * L * iLGp * jRGp * UU +
                        4. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU -
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU * AXG +
                        4. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU * AXG -
                        2. * kL * MGs * pow(MUs, 2) * R * iRGp * jLGp * SS +
                        2. * kL * MGs * MXs * MUs * R * iRGp * jLGp * UU +
                        2. * kL * MGs * pow(MXs, 2) * R * iRGp * jLGp * UU -
                        4. * kL * MGs * pow(MX, 4) * R * iRGp * jLGp * UU +
                        2. * kL * MGs * pow(MX, 4) * R * iRGp * jLGp * UU * AXG -
                        2. * kL * MG * pow(MUs, 2) * Mqi * R * iLGp * jLGp * SS +
                        2. * kL * MG * MXs * MUs * Mqi * R * iLGp * jLGp * UU +
                        2. * kL * MG * pow(MXs, 2) * Mqi * R * iLGp * jLGp * UU -
                        4. * kL * MG * pow(MX, 4) * Mqi * R * iLGp * jLGp * UU +
                        2. * kL * MG * pow(MX, 4) * Mqi * R * iLGp * jLGp * UU * AXG -
                        4. * kL * u * pow(MUs, 2) * R * iRGp * jLGp * SS +
                        2. * kL * u * pow(MUs, 2) * R * iRGp * jLGp * UU -
                        2. * kL * u * MXs * MUs * R * iRGp * jLGp * UU +
                        4. * kL * u * MXs * MUs * R * iRGp * jLGp * UU * AXG -
                        4. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * UU) +
               syFC2 * (-2. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * UU +
                        4. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * UU * AXG +
                        4. * kL * u * MGs * MUs * R * iRGp * jLGp * SS -
                        2. * kL * u * MGs * MUs * R * iRGp * jLGp * UU +
                        2. * kL * u * MGs * MXs * R * iRGp * jLGp * UU -
                        4. * kL * u * MGs * MXs * R * iRGp * jLGp * UU * AXG +
                        4. * kL * u * MG * MUs * Mqi * R * iLGp * jLGp * SS -
                        2. * kL * u * MG * MUs * Mqi * R * iLGp * jLGp * UU +
                        2. * kL * u * MG * MXs * Mqi * R * iLGp * jLGp * UU -
                        4. * kL * u * MG * MXs * Mqi * R * iLGp * jLGp * UU * AXG +
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * SS -
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * UU * AXG -
                        2. * kL * pow(u, 2) * MGs * R * iRGp * jLGp * SS +
                        2. * kL * pow(u, 2) * MGs * R * iRGp * jLGp * UU * AXG -
                        2. * kL * pow(u, 2) * MG * Mqi * R * iLGp * jLGp * SS +
                        2. * kL * pow(u, 2) * MG * Mqi * R * iLGp * jLGp * UU * AXG -
                        2. * kL * s * pow(MUs, 2) * R * iRGp * jLGp * SS -
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU +
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU * AXG -
                        3. * kL * s * MX * MXs * Mqi * L * iLGp * jRGp * UU +
                        3. * kL * s * MX * MXs * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * s * MGs * MUs * R * iRGp * jLGp * SS +
                        kL * s * MGs * MXs * R * iRGp * jLGp * UU -
                        kL * s * MGs * MXs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * s * MG * MUs * Mqi * R * iLGp * jLGp * SS +
                        kL * s * MG * MXs * Mqi * R * iLGp * jLGp * UU) +
               syFC2 * (-kL * s * MG * MXs * Mqi * R * iLGp * jLGp * UU * AXG +
                        2. * kL * s * u * MUs * R * iRGp * jLGp * SS -
                        kL * s * u * MUs * R * iRGp * jLGp * UU -
                        kL * s * u * MUs * R * iRGp * jLGp * UU * AXG +
                        kL * s * u * MX * Mqi * L * iLGp * jRGp * UU -
                        3. * kL * s * u * MX * Mqi * L * iLGp * jRGp * UU * AXG -
                        2. * kL * s * u * MGs * R * iRGp * jLGp * SS +
                        kL * s * u * MGs * R * iRGp * jLGp * UU +
                        kL * s * u * MGs * R * iRGp * jLGp * UU * AXG -
                        2. * kL * s * u * MG * Mqi * R * iLGp * jLGp * SS +
                        kL * s * u * MG * Mqi * R * iLGp * jLGp * UU +
                        kL * s * u * MG * Mqi * R * iLGp * jLGp * UU * AXG +
                        kL * pow(s, 2) * MX * Mqi * L * iLGp * jRGp * UU -
                        kL * pow(s, 2) * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kR * pow(MUs, 3) * L * iLGp * jRGp * SS -
                        2. * kR * MXs * pow(MUs, 2) * L * iLGp * jRGp * UU -
                        2. * kR * pow(MXs, 2) * MUs * L * iLGp * jRGp * UU +
                        2. * kR * MX * pow(MUs, 2) * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * MX * MXs * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * MX * pow(MXs, 2) * Mqi * R * iRGp * jLGp * UU +
                        4. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU -
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU * AXG +
                        4. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU * AXG -
                        2. * kR * MGs * pow(MUs, 2) * L * iLGp * jRGp * SS +
                        2. * kR * MGs * MXs * MUs * L * iLGp * jRGp * UU +
                        2. * kR * MGs * pow(MXs, 2) * L * iLGp * jRGp * UU -
                        4. * kR * MGs * pow(MX, 4) * L * iLGp * jRGp * UU) +
               syFC2 * (2. * kR * MGs * pow(MX, 4) * L * iLGp * jRGp * UU * AXG -
                        2. * kR * MG * pow(MUs, 2) * Mqi * L * iRGp * jRGp * SS +
                        2. * kR * MG * MXs * MUs * Mqi * L * iRGp * jRGp * UU +
                        2. * kR * MG * pow(MXs, 2) * Mqi * L * iRGp * jRGp * UU -
                        4. * kR * MG * pow(MX, 4) * Mqi * L * iRGp * jRGp * UU +
                        2. * kR * MG * pow(MX, 4) * Mqi * L * iRGp * jRGp * UU * AXG -
                        4. * kR * u * pow(MUs, 2) * L * iLGp * jRGp * SS +
                        2. * kR * u * pow(MUs, 2) * L * iLGp * jRGp * UU -
                        2. * kR * u * MXs * MUs * L * iLGp * jRGp * UU +
                        4. * kR * u * MXs * MUs * L * iLGp * jRGp * UU * AXG -
                        4. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * UU +
                        4. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * UU * AXG +
                        4. * kR * u * MGs * MUs * L * iLGp * jRGp * SS -
                        2. * kR * u * MGs * MUs * L * iLGp * jRGp * UU +
                        2. * kR * u * MGs * MXs * L * iLGp * jRGp * UU -
                        4. * kR * u * MGs * MXs * L * iLGp * jRGp * UU * AXG +
                        4. * kR * u * MG * MUs * Mqi * L * iRGp * jRGp * SS -
                        2. * kR * u * MG * MUs * Mqi * L * iRGp * jRGp * UU +
                        2. * kR * u * MG * MXs * Mqi * L * iRGp * jRGp * UU -
                        4. * kR * u * MG * MXs * Mqi * L * iRGp * jRGp * UU * AXG +
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * SS -
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * UU * AXG -
                        2. * kR * pow(u, 2) * MGs * L * iLGp * jRGp * SS) +
               syFC2 * (2. * kR * pow(u, 2) * MGs * L * iLGp * jRGp * UU * AXG -
                        2. * kR * pow(u, 2) * MG * Mqi * L * iRGp * jRGp * SS +
                        2. * kR * pow(u, 2) * MG * Mqi * L * iRGp * jRGp * UU * AXG -
                        2. * kR * s * pow(MUs, 2) * L * iLGp * jRGp * SS -
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU +
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU * AXG -
                        3. * kR * s * MX * MXs * Mqi * R * iRGp * jLGp * UU +
                        3. * kR * s * MX * MXs * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * s * MGs * MUs * L * iLGp * jRGp * SS +
                        kR * s * MGs * MXs * L * iLGp * jRGp * UU -
                        kR * s * MGs * MXs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * s * MG * MUs * Mqi * L * iRGp * jRGp * SS +
                        kR * s * MG * MXs * Mqi * L * iRGp * jRGp * UU -
                        kR * s * MG * MXs * Mqi * L * iRGp * jRGp * UU * AXG +
                        2. * kR * s * u * MUs * L * iLGp * jRGp * SS -
                        kR * s * u * MUs * L * iLGp * jRGp * UU -
                        kR * s * u * MUs * L * iLGp * jRGp * UU * AXG +
                        kR * s * u * MX * Mqi * R * iRGp * jLGp * UU -
                        3. * kR * s * u * MX * Mqi * R * iRGp * jLGp * UU * AXG -
                        2. * kR * s * u * MGs * L * iLGp * jRGp * SS +
                        kR * s * u * MGs * L * iLGp * jRGp * UU +
                        kR * s * u * MGs * L * iLGp * jRGp * UU * AXG -
                        2. * kR * s * u * MG * Mqi * L * iRGp * jRGp * SS +
                        kR * s * u * MG * Mqi * L * iRGp * jRGp * UU +
                        kR * s * u * MG * Mqi * L * iRGp * jRGp * UU * AXG +
                        kR * pow(s, 2) * MX * Mqi * R * iRGp * jLGp * UU -
                        kR * pow(s, 2) * MX * Mqi * R * iRGp * jLGp * UU * AXG) +
               syFC3 * (-2. * kL * pow(MUs, 3) * R * iRGp * jLGp * SS -
                        2. * kL * MXs * pow(MUs, 2) * R * iRGp * jLGp * UU +
                        2. * kL * pow(MXs, 2) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * pow(MXs, 3) * R * iRGp * jLGp * UU +
                        2. * kL * MX * pow(MUs, 2) * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * MX * MXs * MUs * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * pow(MX, 3) * MUs * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * pow(MX, 3) * MUs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 3) * MXs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU -
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU * AXG -
                        6. * kL * pow(MX, 4) * MXs * R * iRGp * jLGp * UU +
                        2. * kL * pow(MX, 4) * MXs * R * iRGp * jLGp * UU * AXG +
                        4. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU * AXG +
                        4. * kL * pow(MX, 6) * R * iRGp * jLGp * UU -
                        2. * kL * pow(MX, 6) * R * iRGp * jLGp * UU * AXG +
                        4. * kL * pow(MU, 4) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * MG * MX * pow(MUs, 2) * L * iRGp * jRGp * SS +
                        2. * kL * MG * MX * MXs * MUs * L * iRGp * jRGp * SS -
                        2. * kL * MG * pow(MX, 3) * MUs * L * iRGp * jRGp * SS -
                        2. * kL * MG * pow(MX, 3) * MUs * L * iRGp * jRGp * UU -
                        2. * kL * MG * pow(MX, 3) * MXs * L * iRGp * jRGp * UU +
                        4. * kL * MG * pow(MX, 5) * L * iRGp * jRGp * UU -
                        2. * kL * MG * pow(MX, 5) * L * iRGp * jRGp * UU * AXG) +
               syFC3 * (2. * kL * u * pow(MUs, 2) * R * iRGp * jLGp * UU -
                        2. * kL * u * MXs * MUs * R * iRGp * jLGp * UU +
                        4. * kL * u * MXs * MUs * R * iRGp * jLGp * UU * AXG -
                        2. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * SS +
                        4. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * UU -
                        4. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * UU +
                        2. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * SS -
                        4. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU +
                        4. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * u * pow(MX, 4) * R * iRGp * jLGp * SS -
                        4. * kL * u * pow(MX, 4) * R * iRGp * jLGp * UU +
                        4. * kL * u * pow(MX, 4) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * u * pow(MU, 4) * R * iRGp * jLGp * SS -
                        4. * kL * u * MG * MX * MUs * L * iRGp * jRGp * SS +
                        2. * kL * u * MG * MX * MUs * L * iRGp * jRGp * UU -
                        2. * kL * u * MG * MX * MXs * L * iRGp * jRGp * SS +
                        2. * kL * u * MG * MX * MXs * L * iRGp * jRGp * UU +
                        2. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * SS -
                        4. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * UU +
                        4. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * UU * AXG +
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * SS -
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * SS) +
               syFC3 * (-2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * pow(u, 2) * MG * MX * L * iRGp * jRGp * SS -
                        2. * kL * pow(u, 2) * MG * MX * L * iRGp * jRGp * UU * AXG +
                        2. * kL * s * pow(MUs, 2) * R * iRGp * jLGp * SS -
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU +
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU * AXG +
                        kL * s * pow(MXs, 2) * R * iRGp * jLGp * UU -
                        kL * s * pow(MXs, 2) * R * iRGp * jLGp * UU * AXG -
                        3. * kL * s * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU +
                        3. * kL * s * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU * AXG -
                        kL * s * pow(MX, 4) * R * iRGp * jLGp * UU +
                        kL * s * pow(MX, 4) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * s * pow(MU, 4) * R * iRGp * jLGp * SS -
                        3. * kL * s * MG * pow(MX, 3) * L * iRGp * jRGp * UU +
                        3. * kL * s * MG * pow(MX, 3) * L * iRGp * jRGp * UU * AXG +
                        2. * kL * s * u * MUs * R * iRGp * jLGp * SS -
                        kL * s * u * MUs * R * iRGp * jLGp * UU -
                        kL * s * u * MUs * R * iRGp * jLGp * UU * AXG +
                        kL * s * u * MX * Mqi * L * iLGp * jRGp * UU -
                        3. * kL * s * u * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        kL * s * u * MG * MX * L * iRGp * jRGp * UU -
                        3. * kL * s * u * MG * MX * L * iRGp * jRGp * UU * AXG +
                        kL * pow(s, 2) * MX * Mqi * L * iLGp * jRGp * UU -
                        kL * pow(s, 2) * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        kL * pow(s, 2) * MG * MX * L * iRGp * jRGp * UU -
                        kL * pow(s, 2) * MG * MX * L * iRGp * jRGp * UU * AXG -
                        2. * kR * pow(MUs, 3) * L * iLGp * jRGp * SS -
                        2. * kR * MXs * pow(MUs, 2) * L * iLGp * jRGp * UU) +
               syFC3 * (2. * kR * pow(MXs, 2) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * pow(MXs, 3) * L * iLGp * jRGp * UU +
                        2. * kR * MX * pow(MUs, 2) * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * MX * MXs * MUs * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(MX, 3) * MUs * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(MX, 3) * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 3) * MXs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU -
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU * AXG -
                        6. * kR * pow(MX, 4) * MXs * L * iLGp * jRGp * UU +
                        2. * kR * pow(MX, 4) * MXs * L * iLGp * jRGp * UU * AXG +
                        4. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU * AXG +
                        4. * kR * pow(MX, 6) * L * iLGp * jRGp * UU -
                        2. * kR * pow(MX, 6) * L * iLGp * jRGp * UU * AXG +
                        4. * kR * pow(MU, 4) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * MG * MX * pow(MUs, 2) * R * iLGp * jLGp * SS +
                        2. * kR * MG * MX * MXs * MUs * R * iLGp * jLGp * SS -
                        2. * kR * MG * pow(MX, 3) * MUs * R * iLGp * jLGp * SS -
                        2. * kR * MG * pow(MX, 3) * MUs * R * iLGp * jLGp * UU -
                        2. * kR * MG * pow(MX, 3) * MXs * R * iLGp * jLGp * UU +
                        4. * kR * MG * pow(MX, 5) * R * iLGp * jLGp * UU -
                        2. * kR * MG * pow(MX, 5) * R * iLGp * jLGp * UU * AXG +
                        2. * kR * u * pow(MUs, 2) * L * iLGp * jRGp * UU -
                        2. * kR * u * MXs * MUs * L * iLGp * jRGp * UU +
                        4. * kR * u * MXs * MUs * L * iLGp * jRGp * UU * AXG) +
               syFC3 * (-2. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * SS +
                        4. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * UU -
                        4. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * UU +
                        2. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * SS -
                        4. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU +
                        4. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * u * pow(MX, 4) * L * iLGp * jRGp * SS -
                        4. * kR * u * pow(MX, 4) * L * iLGp * jRGp * UU +
                        4. * kR * u * pow(MX, 4) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * u * pow(MU, 4) * L * iLGp * jRGp * SS -
                        4. * kR * u * MG * MX * MUs * R * iLGp * jLGp * SS +
                        2. * kR * u * MG * MX * MUs * R * iLGp * jLGp * UU -
                        2. * kR * u * MG * MX * MXs * R * iLGp * jLGp * SS +
                        2. * kR * u * MG * MX * MXs * R * iLGp * jLGp * UU +
                        2. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * SS -
                        4. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * UU +
                        4. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * UU * AXG +
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * SS -
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * pow(u, 2) * MG * MX * R * iLGp * jLGp * SS -
                        2. * kR * pow(u, 2) * MG * MX * R * iLGp * jLGp * UU * AXG) +
               syFC3 * (2. * kR * s * pow(MUs, 2) * L * iLGp * jRGp * SS -
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU +
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU * AXG +
                        kR * s * pow(MXs, 2) * L * iLGp * jRGp * UU -
                        kR * s * pow(MXs, 2) * L * iLGp * jRGp * UU * AXG -
                        3. * kR * s * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU +
                        3. * kR * s * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU * AXG -
                        kR * s * pow(MX, 4) * L * iLGp * jRGp * UU +
                        kR * s * pow(MX, 4) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * s * pow(MU, 4) * L * iLGp * jRGp * SS -
                        3. * kR * s * MG * pow(MX, 3) * R * iLGp * jLGp * UU +
                        3. * kR * s * MG * pow(MX, 3) * R * iLGp * jLGp * UU * AXG +
                        2. * kR * s * u * MUs * L * iLGp * jRGp * SS -
                        kR * s * u * MUs * L * iLGp * jRGp * UU -
                        kR * s * u * MUs * L * iLGp * jRGp * UU * AXG +
                        kR * s * u * MX * Mqi * R * iRGp * jLGp * UU -
                        3. * kR * s * u * MX * Mqi * R * iRGp * jLGp * UU * AXG +
                        kR * s * u * MG * MX * R * iLGp * jLGp * UU -
                        3. * kR * s * u * MG * MX * R * iLGp * jLGp * UU * AXG +
                        kR * pow(s, 2) * MX * Mqi * R * iRGp * jLGp * UU -
                        kR * pow(s, 2) * MX * Mqi * R * iRGp * jLGp * UU * AXG +
                        kR * pow(s, 2) * MG * MX * R * iLGp * jLGp * UU -
                        kR * pow(s, 2) * MG * MX * R * iLGp * jLGp * UU * AXG) +
               syFC4 * (-2. * kL * pow(MUs, 3) * R * iRGp * jLGp * SS -
                        2. * kL * MXs * pow(MUs, 2) * R * iRGp * jLGp * UU +
                        2. * kL * pow(MXs, 2) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * pow(MXs, 3) * R * iRGp * jLGp * UU +
                        2. * kL * MX * pow(MUs, 2) * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * MX * MXs * MUs * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * pow(MX, 3) * MUs * Mqi * L * iLGp * jRGp * SS -
                        2. * kL * pow(MX, 3) * MUs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 3) * MXs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU -
                        2. * kL * pow(MX, 4) * MUs * R * iRGp * jLGp * UU * AXG -
                        6. * kL * pow(MX, 4) * MXs * R * iRGp * jLGp * UU +
                        2. * kL * pow(MX, 4) * MXs * R * iRGp * jLGp * UU * AXG +
                        4. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * pow(MX, 5) * Mqi * L * iLGp * jRGp * UU * AXG +
                        4. * kL * pow(MX, 6) * R * iRGp * jLGp * UU -
                        2. * kL * pow(MX, 6) * R * iRGp * jLGp * UU * AXG +
                        4. * kL * pow(MU, 4) * MUs * R * iRGp * jLGp * SS +
                        2. * kL * MG * MX * pow(MUs, 2) * L * iRGp * jRGp * SS +
                        2. * kL * MG * MX * MXs * MUs * L * iRGp * jRGp * SS -
                        2. * kL * MG * pow(MX, 3) * MUs * L * iRGp * jRGp * SS -
                        2. * kL * MG * pow(MX, 3) * MUs * L * iRGp * jRGp * UU -
                        2. * kL * MG * pow(MX, 3) * MXs * L * iRGp * jRGp * UU +
                        4. * kL * MG * pow(MX, 5) * L * iRGp * jRGp * UU -
                        2. * kL * MG * pow(MX, 5) * L * iRGp * jRGp * UU * AXG) +
               syFC4 * (2. * kL * u * pow(MUs, 2) * R * iRGp * jLGp * UU -
                        2. * kL * u * MXs * MUs * R * iRGp * jLGp * UU +
                        4. * kL * u * MXs * MUs * R * iRGp * jLGp * UU * AXG -
                        2. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * SS +
                        4. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * UU -
                        4. * kL * u * pow(MXs, 2) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * u * MX * MUs * Mqi * L * iLGp * jRGp * UU -
                        2. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * SS +
                        2. * kL * u * MX * MXs * Mqi * L * iLGp * jRGp * UU +
                        2. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * SS -
                        4. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU +
                        4. * kL * u * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * u * pow(MX, 4) * R * iRGp * jLGp * SS -
                        4. * kL * u * pow(MX, 4) * R * iRGp * jLGp * UU +
                        4. * kL * u * pow(MX, 4) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * u * pow(MU, 4) * R * iRGp * jLGp * SS -
                        4. * kL * u * MG * MX * MUs * L * iRGp * jRGp * SS +
                        2. * kL * u * MG * MX * MUs * L * iRGp * jRGp * UU -
                        2. * kL * u * MG * MX * MXs * L * iRGp * jRGp * SS +
                        2. * kL * u * MG * MX * MXs * L * iRGp * jRGp * UU +
                        2. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * SS -
                        4. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * UU +
                        4. * kL * u * MG * pow(MX, 3) * L * iRGp * jRGp * UU * AXG +
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * SS -
                        2. * kL * pow(u, 2) * MUs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * SS) +
               syFC4 * (-2. * kL * pow(u, 2) * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * pow(u, 2) * MG * MX * L * iRGp * jRGp * SS -
                        2. * kL * pow(u, 2) * MG * MX * L * iRGp * jRGp * UU * AXG -
                        2. * kL * s * MXs * MUs * R * iRGp * jLGp * SS +
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU +
                        kL * s * MXs * MUs * R * iRGp * jLGp * UU * AXG +
                        3. * kL * s * pow(MXs, 2) * R * iRGp * jLGp * UU -
                        kL * s * pow(MXs, 2) * R * iRGp * jLGp * UU * AXG -
                        2. * kL * s * MX * MUs * Mqi * L * iLGp * jRGp * SS -
                        kL * s * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU +
                        kL * s * pow(MX, 3) * Mqi * L * iLGp * jRGp * UU * AXG -
                        3. * kL * s * pow(MX, 4) * R * iRGp * jLGp * UU +
                        kL * s * pow(MX, 4) * R * iRGp * jLGp * UU * AXG -
                        4. * kL * s * pow(MU, 4) * R * iRGp * jLGp * SS -
                        2. * kL * s * MG * MX * MUs * L * iRGp * jRGp * SS -
                        kL * s * MG * pow(MX, 3) * L * iRGp * jRGp * UU +
                        kL * s * MG * pow(MX, 3) * L * iRGp * jRGp * UU * AXG +
                        6. * kL * s * u * MUs * R * iRGp * jLGp * SS -
                        3. * kL * s * u * MUs * R * iRGp * jLGp * UU -
                        kL * s * u * MUs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * s * u * MXs * R * iRGp * jLGp * SS -
                        2. * kL * s * u * MXs * R * iRGp * jLGp * UU * AXG +
                        2. * kL * s * u * MX * Mqi * L * iLGp * jRGp * SS -
                        kL * s * u * MX * Mqi * L * iLGp * jRGp * UU -
                        kL * s * u * MX * Mqi * L * iLGp * jRGp * UU * AXG +
                        2. * kL * s * u * MG * MX * L * iRGp * jRGp * SS -
                        kL * s * u * MG * MX * L * iRGp * jRGp * UU -
                        kL * s * u * MG * MX * L * iRGp * jRGp * UU * AXG -
                        2. * kL * s * pow(u, 2) * R * iRGp * jLGp * SS) +
               syFC4 * (2. * kL * s * pow(u, 2) * R * iRGp * jLGp * UU * AXG +
                        2. * kL * pow(s, 2) * MUs * R * iRGp * jLGp * SS -
                        2. * kL * pow(s, 2) * u * R * iRGp * jLGp * SS +
                        kL * pow(s, 2) * u * R * iRGp * jLGp * UU +
                        kL * pow(s, 2) * u * R * iRGp * jLGp * UU * AXG -
                        2. * kR * pow(MUs, 3) * L * iLGp * jRGp * SS -
                        2. * kR * MXs * pow(MUs, 2) * L * iLGp * jRGp * UU +
                        2. * kR * pow(MXs, 2) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * pow(MXs, 3) * L * iLGp * jRGp * UU +
                        2. * kR * MX * pow(MUs, 2) * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * MX * MXs * MUs * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(MX, 3) * MUs * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(MX, 3) * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 3) * MXs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU -
                        2. * kR * pow(MX, 4) * MUs * L * iLGp * jRGp * UU * AXG -
                        6. * kR * pow(MX, 4) * MXs * L * iLGp * jRGp * UU +
                        2. * kR * pow(MX, 4) * MXs * L * iLGp * jRGp * UU * AXG +
                        4. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * pow(MX, 5) * Mqi * R * iRGp * jLGp * UU * AXG +
                        4. * kR * pow(MX, 6) * L * iLGp * jRGp * UU -
                        2. * kR * pow(MX, 6) * L * iLGp * jRGp * UU * AXG +
                        4. * kR * pow(MU, 4) * MUs * L * iLGp * jRGp * SS +
                        2. * kR * MG * MX * pow(MUs, 2) * R * iLGp * jLGp * SS +
                        2. * kR * MG * MX * MXs * MUs * R * iLGp * jLGp * SS -
                        2. * kR * MG * pow(MX, 3) * MUs * R * iLGp * jLGp * SS) +
               syFC4 * (-2. * kR * MG * pow(MX, 3) * MUs * R * iLGp * jLGp * UU -
                        2. * kR * MG * pow(MX, 3) * MXs * R * iLGp * jLGp * UU +
                        4. * kR * MG * pow(MX, 5) * R * iLGp * jLGp * UU -
                        2. * kR * MG * pow(MX, 5) * R * iLGp * jLGp * UU * AXG +
                        2. * kR * u * pow(MUs, 2) * L * iLGp * jRGp * UU -
                        2. * kR * u * MXs * MUs * L * iLGp * jRGp * UU +
                        4. * kR * u * MXs * MUs * L * iLGp * jRGp * UU * AXG -
                        2. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * SS +
                        4. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * UU -
                        4. * kR * u * pow(MXs, 2) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * u * MX * MUs * Mqi * R * iRGp * jLGp * UU -
                        2. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * SS +
                        2. * kR * u * MX * MXs * Mqi * R * iRGp * jLGp * UU +
                        2. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * SS -
                        4. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU +
                        4. * kR * u * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * u * pow(MX, 4) * L * iLGp * jRGp * SS -
                        4. * kR * u * pow(MX, 4) * L * iLGp * jRGp * UU +
                        4. * kR * u * pow(MX, 4) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * u * pow(MU, 4) * L * iLGp * jRGp * SS -
                        4. * kR * u * MG * MX * MUs * R * iLGp * jLGp * SS +
                        2. * kR * u * MG * MX * MUs * R * iLGp * jLGp * UU -
                        2. * kR * u * MG * MX * MXs * R * iLGp * jLGp * SS +
                        2. * kR * u * MG * MX * MXs * R * iLGp * jLGp * UU +
                        2. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * SS -
                        4. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * UU) +
               syFC4 * (4. * kR * u * MG * pow(MX, 3) * R * iLGp * jLGp * UU * AXG +
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * SS -
                        2. * kR * pow(u, 2) * MUs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * SS -
                        2. * kR * pow(u, 2) * MX * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * pow(u, 2) * MG * MX * R * iLGp * jLGp * SS -
                        2. * kR * pow(u, 2) * MG * MX * R * iLGp * jLGp * UU * AXG -
                        2. * kR * s * MXs * MUs * L * iLGp * jRGp * SS +
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU +
                        kR * s * MXs * MUs * L * iLGp * jRGp * UU * AXG +
                        3. * kR * s * pow(MXs, 2) * L * iLGp * jRGp * UU -
                        kR * s * pow(MXs, 2) * L * iLGp * jRGp * UU * AXG -
                        2. * kR * s * MX * MUs * Mqi * R * iRGp * jLGp * SS -
                        kR * s * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU +
                        kR * s * pow(MX, 3) * Mqi * R * iRGp * jLGp * UU * AXG -
                        3. * kR * s * pow(MX, 4) * L * iLGp * jRGp * UU +
                        kR * s * pow(MX, 4) * L * iLGp * jRGp * UU * AXG -
                        4. * kR * s * pow(MU, 4) * L * iLGp * jRGp * SS -
                        2. * kR * s * MG * MX * MUs * R * iLGp * jLGp * SS -
                        kR * s * MG * pow(MX, 3) * R * iLGp * jLGp * UU +
                        kR * s * MG * pow(MX, 3) * R * iLGp * jLGp * UU * AXG +
                        6. * kR * s * u * MUs * L * iLGp * jRGp * SS -
                        3. * kR * s * u * MUs * L * iLGp * jRGp * UU -
                        kR * s * u * MUs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * s * u * MXs * L * iLGp * jRGp * SS -
                        2. * kR * s * u * MXs * L * iLGp * jRGp * UU * AXG +
                        2. * kR * s * u * MX * Mqi * R * iRGp * jLGp * SS -
                        kR * s * u * MX * Mqi * R * iRGp * jLGp * UU) +
               syFC4 * (-kR * s * u * MX * Mqi * R * iRGp * jLGp * UU * AXG +
                        2. * kR * s * u * MG * MX * R * iLGp * jLGp * SS -
                        kR * s * u * MG * MX * R * iLGp * jLGp * UU -
                        kR * s * u * MG * MX * R * iLGp * jLGp * UU * AXG -
                        2. * kR * s * pow(u, 2) * L * iLGp * jRGp * SS +
                        2. * kR * s * pow(u, 2) * L * iLGp * jRGp * UU * AXG +
                        2. * kR * pow(s, 2) * MUs * L * iLGp * jRGp * SS -
                        2. * kR * pow(s, 2) * u * L * iLGp * jRGp * SS +
                        kR * pow(s, 2) * u * L * iLGp * jRGp * UU +
                        kR * pow(s, 2) * u * L * iLGp * jRGp * UU * AXG)));
    //}
  }
  return ret.real();
}
ComplexType ME_us_qQX_qgQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  return ME_ss_qQX_qgQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
         ME_uu_qQX_qgQ(pIEPS, sc, uc, axial, Q2, P1K1, params);
}

ComplexType ME_us_qQX_QGq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  return ME_ss_qQX_QGq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
         ME_uu_qQX_QGq(pIEPS, sc, uc, axial, Q2, P1K1, params);
}

#undef A0
#undef B0
#undef C0
#undef C1
#undef C2
#undef C00
#undef C11
#undef C12
#undef C22