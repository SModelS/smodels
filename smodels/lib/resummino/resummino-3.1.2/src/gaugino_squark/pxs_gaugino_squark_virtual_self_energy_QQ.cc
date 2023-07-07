#include "pxs_gausq_3.h"

ComplexType ME_us_se_QQ_gQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto MQj = MQ;
  auto MQjs = MQs;
  auto MQi = MQ;
  auto MQis = MQs;

  auto Denom = [](auto a) { return 1. / a; };
  //*

#define syFC1 A0(MQis)
#define syFC2 B0(MUs + MXs - s - t, 0, MQis)

  _EPS0(
      ret,
      -(-1) * (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
          (+Denom(
              -768 * pow(Pi, 2) * s * MXs * MQs * MQjs * Nc +
              768 * pow(Pi, 2) * s * MXs * MUs * MQjs * Nc +
              768 * pow(Pi, 2) * s * MXs * MUs * MQs * Nc -
              768 * pow(Pi, 2) * s * MXs * pow(MUs, 2) * Nc +
              768 * pow(Pi, 2) * s * pow(MXs, 2) * MQjs * Nc +
              768 * pow(Pi, 2) * s * pow(MXs, 2) * MQs * Nc -
              1536 * pow(Pi, 2) * s * pow(MXs, 2) * MUs * Nc -
              768 * pow(Pi, 2) * s * pow(MXs, 3) * Nc + 768 * pow(Pi, 2) * s * t * MQs * MQjs * Nc -
              768 * pow(Pi, 2) * s * t * MUs * MQjs * Nc -
              768 * pow(Pi, 2) * s * t * MUs * MQs * Nc +
              768 * pow(Pi, 2) * s * t * pow(MUs, 2) * Nc -
              1536 * pow(Pi, 2) * s * t * MXs * MQjs * Nc -
              1536 * pow(Pi, 2) * s * t * MXs * MQs * Nc +
              3072 * pow(Pi, 2) * s * t * MXs * MUs * Nc +
              2304 * pow(Pi, 2) * s * t * pow(MXs, 2) * Nc +
              768 * pow(Pi, 2) * s * pow(t, 2) * MQjs * Nc +
              768 * pow(Pi, 2) * s * pow(t, 2) * MQs * Nc -
              1536 * pow(Pi, 2) * s * pow(t, 2) * MUs * Nc -
              2304 * pow(Pi, 2) * s * pow(t, 2) * MXs * Nc + 768 * pow(Pi, 2) * s * pow(t, 3) * Nc +
              768 * pow(Pi, 2) * pow(s, 2) * MQs * MQjs * Nc -
              768 * pow(Pi, 2) * pow(s, 2) * MUs * MQjs * Nc -
              768 * pow(Pi, 2) * pow(s, 2) * MUs * MQs * Nc +
              768 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) * Nc -
              1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQjs * Nc -
              1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQs * Nc +
              3072 * pow(Pi, 2) * pow(s, 2) * MXs * MUs * Nc +
              2304 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) * Nc +
              1536 * pow(Pi, 2) * pow(s, 2) * t * MQjs * Nc +
              1536 * pow(Pi, 2) * pow(s, 2) * t * MQs * Nc -
              3072 * pow(Pi, 2) * pow(s, 2) * t * MUs * Nc -
              4608 * pow(Pi, 2) * pow(s, 2) * t * MXs * Nc +
              2304 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) * Nc +
              768 * pow(Pi, 2) * pow(s, 3) * MQjs * Nc + 768 * pow(Pi, 2) * pow(s, 3) * MQs * Nc -
              1536 * pow(Pi, 2) * pow(s, 3) * MUs * Nc - 2304 * pow(Pi, 2) * pow(s, 3) * MXs * Nc +
              2304 * pow(Pi, 2) * pow(s, 3) * t * Nc + 768 * pow(Pi, 2) * pow(s, 4) * Nc)) *
          (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 * (1.) + syFC2 * (-2 * MQis - 2 * MUs - 2 * MXs + 2 * t + 2 * s)) *
          (+2 * MXs * pow(MUs, 2) * SS - 2 * MXs * pow(MUs, 2) * SS * AXG +
           4 * MXs * pow(MUs, 2) * UU * AXG - 2 * pow(MXs, 2) * MUs * SS +
           8 * pow(MXs, 2) * MUs * UU * AXG - 4 * pow(MXs, 3) * SS + 2 * pow(MXs, 3) * SS * AXG +
           4 * pow(MXs, 3) * UU * AXG - 8 * pow(MX, 4) * MUs * UU * AXG +
           4 * pow(MX, 4) * MXs * SS - 2 * pow(MX, 4) * MXs * SS * AXG -
           8 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG - 2 * t * pow(MUs, 2) * SS +
           2 * t * pow(MUs, 2) * SS * AXG - 4 * t * pow(MUs, 2) * UU * AXG +
           4 * t * MXs * MUs * SS * AXG - 8 * t * MXs * MUs * UU * AXG + 6 * t * pow(MXs, 2) * SS -
           2 * t * pow(MXs, 2) * SS * AXG - 4 * t * pow(MXs, 2) * UU * AXG -
           4 * t * pow(MX, 4) * SS + 2 * t * pow(MX, 4) * SS * AXG + 4 * t * pow(MX, 4) * UU * AXG +
           2 * pow(t, 2) * MUs * SS - 4 * pow(t, 2) * MUs * SS * AXG +
           8 * pow(t, 2) * MUs * UU * AXG - 2 * pow(t, 2) * MXs * SS -
           2 * pow(t, 2) * MXs * SS * AXG + 4 * pow(t, 2) * MXs * UU * AXG +
           2 * pow(t, 3) * SS * AXG - 4 * pow(t, 3) * UU * AXG - 2 * s * pow(MUs, 2) * SS +
           2 * s * pow(MUs, 2) * SS * AXG + 4 * s * pow(MUs, 2) * UU -
           4 * s * pow(MUs, 2) * UU * AXG - s * MXs * MUs * SS + 3 * s * MXs * MUs * SS * AXG +
           2 * s * MXs * MUs * UU - 6 * s * MXs * MUs * UU * AXG + 4 * s * pow(MXs, 2) * SS -
           2 * s * pow(MXs, 2) * SS * AXG - 2 * s * pow(MXs, 2) * UU -
           2 * s * pow(MXs, 2) * UU * AXG - 4 * s * pow(MX, 4) * SS +
           2 * s * pow(MX, 4) * SS * AXG + 2 * s * pow(MX, 4) * UU + 2 * s * pow(MX, 4) * UU * AXG +
           5 * s * t * MUs * SS - 7 * s * t * MUs * SS * AXG - 6 * s * t * MUs * UU +
           14 * s * t * MUs * UU * AXG - s * t * MXs * SS - 3 * s * t * MXs * SS * AXG -
           2 * s * t * MXs * UU + 6 * s * t * MXs * UU * AXG - s * pow(t, 2) * SS +
           5 * s * pow(t, 2) * SS * AXG + 2 * s * pow(t, 2) * UU - 10 * s * pow(t, 2) * UU * AXG +
           3 * pow(s, 2) * MUs * SS - 3 * pow(s, 2) * MUs * SS * AXG - 6 * pow(s, 2) * MUs * UU +
           6 * pow(s, 2) * MUs * UU * AXG + pow(s, 2) * MXs * SS - pow(s, 2) * MXs * SS * AXG -
           2 * pow(s, 2) * MXs * UU + 2 * pow(s, 2) * MXs * UU * AXG - 2 * pow(s, 2) * t * SS +
           4 * pow(s, 2) * t * SS * AXG + 4 * pow(s, 2) * t * UU - 8 * pow(s, 2) * t * UU * AXG -
           pow(s, 3) * SS + pow(s, 3) * SS * AXG + 2 * pow(s, 3) * UU - 2 * pow(s, 3) * UU * AXG));

  return ret.real();
}

ComplexType ME_us_se_QQ_Gq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;

  // int itsq = 0;
  // int itq = 0;

  /*
    for (int itsq = 0; itsq < 6; itsq++) {
      for (int itq = 0; itq < 3; itq++) {
        // != acts as XOR flipping up and down for chargino case
        int isq = (is_chargino(ch) != is_up_quark(q)) * 6 + itsq;
        int iq = (is_up_squark(sq)) * 3 + itq;
        */
  int itq = q;
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int itq = 0; itq < 3; itq++) {
    int isq = (is_chargino(ch) != is_up_quark(itq)) * 6 + itsq * 3 + itq - (is_up_quark(itq)) * 3;
    int iq = (itq + is_chargino(ch) * 3) % 6;


    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType Lp = conj(params->CHSQq[ch][isq][q].R);
    ComplexType Rp = conj(params->CHSQq[ch][isq][q].L);

    ComplexType LGp = conj(params->GLSQq[sq][iq].R);
    ComplexType RGp = conj(params->GLSQq[sq][iq].L);
    ComplexType LG = params->GLSQq[isq][iq].L;
    ComplexType RG = params->GLSQq[isq][iq].R;

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    auto Denom = [](auto a) { return 1. / a; };
    //*

#define syFC1 A0(MGs)
#define syFC2 A0(Mqis)
#define syFC3 B0(MUs + MXs - s - t, MGs, Mqis)

    _EPS0(ret,
          (-1) * (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
              (+R * Lp + L * Rp) * (+gs) * (+gs) *
              (+Denom(-768 * pow(Pi, 2) * s * MXs * MQs * MQis * Nc +
                      768 * pow(Pi, 2) * s * MXs * MUs * MQis * Nc +
                      768 * pow(Pi, 2) * s * MXs * MUs * MQs * Nc -
                      768 * pow(Pi, 2) * s * MXs * pow(MUs, 2) * Nc +
                      768 * pow(Pi, 2) * s * pow(MXs, 2) * MQis * Nc +
                      768 * pow(Pi, 2) * s * pow(MXs, 2) * MQs * Nc -
                      1536 * pow(Pi, 2) * s * pow(MXs, 2) * MUs * Nc -
                      768 * pow(Pi, 2) * s * pow(MXs, 3) * Nc +
                      768 * pow(Pi, 2) * s * t * MQs * MQis * Nc -
                      768 * pow(Pi, 2) * s * t * MUs * MQis * Nc -
                      768 * pow(Pi, 2) * s * t * MUs * MQs * Nc +
                      768 * pow(Pi, 2) * s * t * pow(MUs, 2) * Nc -
                      1536 * pow(Pi, 2) * s * t * MXs * MQis * Nc -
                      1536 * pow(Pi, 2) * s * t * MXs * MQs * Nc +
                      3072 * pow(Pi, 2) * s * t * MXs * MUs * Nc +
                      2304 * pow(Pi, 2) * s * t * pow(MXs, 2) * Nc +
                      768 * pow(Pi, 2) * s * pow(t, 2) * MQis * Nc +
                      768 * pow(Pi, 2) * s * pow(t, 2) * MQs * Nc -
                      1536 * pow(Pi, 2) * s * pow(t, 2) * MUs * Nc -
                      2304 * pow(Pi, 2) * s * pow(t, 2) * MXs * Nc +
                      768 * pow(Pi, 2) * s * pow(t, 3) * Nc +
                      768 * pow(Pi, 2) * pow(s, 2) * MQs * MQis * Nc -
                      768 * pow(Pi, 2) * pow(s, 2) * MUs * MQis * Nc -
                      768 * pow(Pi, 2) * pow(s, 2) * MUs * MQs * Nc +
                      768 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) * Nc -
                      1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQis * Nc -
                      1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQs * Nc +
                      3072 * pow(Pi, 2) * pow(s, 2) * MXs * MUs * Nc +
                      2304 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) * Nc +
                      1536 * pow(Pi, 2) * pow(s, 2) * t * MQis * Nc +
                      1536 * pow(Pi, 2) * pow(s, 2) * t * MQs * Nc -
                      3072 * pow(Pi, 2) * pow(s, 2) * t * MUs * Nc -
                      4608 * pow(Pi, 2) * pow(s, 2) * t * MXs * Nc +
                      2304 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) * Nc +
                      768 * pow(Pi, 2) * pow(s, 3) * MQis * Nc +
                      768 * pow(Pi, 2) * pow(s, 3) * MQs * Nc -
                      1536 * pow(Pi, 2) * pow(s, 3) * MUs * Nc -
                      2304 * pow(Pi, 2) * pow(s, 3) * MXs * Nc +
                      2304 * pow(Pi, 2) * pow(s, 3) * t * Nc + 768 * pow(Pi, 2) * pow(s, 4) * Nc)) *
              (+syFC1 * (RG * LGp + LG * RGp) + syFC2 * (RG * LGp + LG * RGp) +
               syFC3 * (Mqis * RG * LGp + Mqis * LG * RGp + MGs * RG * LGp + MGs * LG * RGp -
                        2 * MG * Mqi * RG * RGp - 2 * MG * Mqi * LG * LGp - u * RG * LGp -
                        u * LG * RGp)) *
              (+2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
               4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
               4 * pow(MX, 4) * MUs * UU * AXG - 4 * pow(MX, 4) * MXs * UU * AXG +
               4 * pow(MX, 6) * UU * AXG - 2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
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
               2 * s * pow(u, 2) * UU * AXG));
    //}
  }

  return ret.real();
}

ComplexType ME_us_se_QQ_g(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  return 0;
}

ComplexType ME_us_se_QQ_Q(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  // Tensor<ComplexType, 4> ret;
  // ret.zeros();
  ComplexType ret = 0;

  int itsq = 0;
  int itq = 0;

  for (int itsq = 0; itsq < 6; itsq++) {
    // for(int itq  = 0; itq  < 3;itq++){

    int isq = (is_chargino(ch) != is_up_quark(q)) * 6 + itsq;
    // int iq = is_up_squark(sq) * 3 + itq;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType Lp = conj(params->CHSQq[ch][isq][q].R);
    ComplexType Rp = conj(params->CHSQq[ch][isq][q].L);

    // ComplexType LGp = conj(params->GLSQq[sq][iq].R);
    // ComplexType RGp = conj(params->GLSQq[sq][iq].L);
    // ComplexType LG  = params->GLSQq[isq][iq].L;
    // ComplexType RG  = params->GLSQq[isq][iq].R;

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    // auto Mqi = params->mq[iq];
    // auto Mqis = pow2(Mqi);

    auto Denom = [](auto a) { return 1. / a; };
    //*

#define syFC1 A0(MQis)

    int k = itsq;
    int tsq = (is_up_squark(sq) ? -6 : 0) + sq;
    auto X = (NC * ((is_up_squark(sq) ? params->XUU : params->XDD)[tsq][tsq][k][k]) +
              (is_up_squark(sq) ? params->XUU : params->XDD)[tsq][k][k][tsq]);

    _EPS0(
        ret,
        (1. / 6.) * (-1 + Nc) * (+1 + Nc) * (+TR) * (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) *
            (+gs) *
            (+Denom(
                -768 * pow(Pi, 2) * s * MXs * MQs * MQis + 768 * pow(Pi, 2) * s * MXs * MUs * MQis +
                768 * pow(Pi, 2) * s * MXs * MUs * MQs - 768 * pow(Pi, 2) * s * MXs * pow(MUs, 2) +
                768 * pow(Pi, 2) * s * pow(MXs, 2) * MQis +
                768 * pow(Pi, 2) * s * pow(MXs, 2) * MQs -
                1536 * pow(Pi, 2) * s * pow(MXs, 2) * MUs - 768 * pow(Pi, 2) * s * pow(MXs, 3) +
                768 * pow(Pi, 2) * s * t * MQs * MQis - 768 * pow(Pi, 2) * s * t * MUs * MQis -
                768 * pow(Pi, 2) * s * t * MUs * MQs + 768 * pow(Pi, 2) * s * t * pow(MUs, 2) -
                1536 * pow(Pi, 2) * s * t * MXs * MQis - 1536 * pow(Pi, 2) * s * t * MXs * MQs +
                3072 * pow(Pi, 2) * s * t * MXs * MUs + 2304 * pow(Pi, 2) * s * t * pow(MXs, 2) +
                768 * pow(Pi, 2) * s * pow(t, 2) * MQis + 768 * pow(Pi, 2) * s * pow(t, 2) * MQs -
                1536 * pow(Pi, 2) * s * pow(t, 2) * MUs - 2304 * pow(Pi, 2) * s * pow(t, 2) * MXs +
                768 * pow(Pi, 2) * s * pow(t, 3) + 768 * pow(Pi, 2) * pow(s, 2) * MQs * MQis -
                768 * pow(Pi, 2) * pow(s, 2) * MUs * MQis -
                768 * pow(Pi, 2) * pow(s, 2) * MUs * MQs +
                768 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) -
                1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQis -
                1536 * pow(Pi, 2) * pow(s, 2) * MXs * MQs +
                3072 * pow(Pi, 2) * pow(s, 2) * MXs * MUs +
                2304 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) +
                1536 * pow(Pi, 2) * pow(s, 2) * t * MQis + 1536 * pow(Pi, 2) * pow(s, 2) * t * MQs -
                3072 * pow(Pi, 2) * pow(s, 2) * t * MUs - 4608 * pow(Pi, 2) * pow(s, 2) * t * MXs +
                2304 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) + 768 * pow(Pi, 2) * pow(s, 3) * MQis +
                768 * pow(Pi, 2) * pow(s, 3) * MQs - 1536 * pow(Pi, 2) * pow(s, 3) * MUs -
                2304 * pow(Pi, 2) * pow(s, 3) * MXs + 2304 * pow(Pi, 2) * pow(s, 3) * t +
                768 * pow(Pi, 2) * pow(s, 4))) *
            (+syFC1 * (1.)) *
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
            (+X));

    //}
  }

  return ret.real();
}
