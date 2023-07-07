#include "pxs_gausq_3.h"

//#define A0 GetA<0>
//#define B0 GetB<0>
//#define B1 GetB<1>
//#define C0 GetC<0>
//#define C1 GetC<1>
//#define C2 GetC<2>
//#define C00 GetC<0, 0>
//#define C11 GetC<1, 1>
//#define C12 GetC<1, 2>
//#define C22 GetC<2, 2>
//#define D0 GetD<0>
//#define D1 GetD<1>
//#define D2 GetD<2>
//#define D3 GetD<3>
//#define D00 GetD<0, 0>
//#define D11 GetD<1, 1>
//#define D12 GetD<1, 2>
//#define D13 GetD<1, 3>
//#define D23 GetD<2, 3>
//#define D33 GetD<3, 3>

ComplexType ME_us_se_qq_gq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
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

  auto Denom = [](auto a) { return 1. / a; };
  //*

#define syFC1 B0(s, 0, 0)

  _EPS0(ret, (-1) * (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
                 (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
                 (+Denom(-768 * pow(Pi, 2) * s * MXs * Nc + 768 * pow(Pi, 2) * s * t * Nc +
                         768 * pow(Pi, 2) * pow(s, 2) * Nc)) *
                 (+syFC1 * (1.)) *
                 (-2 * pow(MUs, 2) * UU + 2 * pow(MUs, 2) * UU * AXG + 2 * MXs * MUs * UU -
                  2 * pow(MXs, 2) * SS + 4 * pow(MXs, 2) * UU - 2 * pow(MXs, 2) * UU * AXG -
                  4 * pow(MX, 4) * UU + 2 * pow(MX, 4) * UU * AXG + 2 * t * MUs * UU -
                  4 * t * MUs * UU * AXG + 4 * t * MXs * SS - 2 * t * MXs * UU -
                  2 * pow(t, 2) * SS + 2 * pow(t, 2) * UU * AXG + 3 * s * MUs * UU -
                  3 * s * MUs * UU * AXG + 2 * s * MXs * SS - 2 * s * t * SS - s * t * UU +
                  3 * s * t * UU * AXG - pow(s, 2) * UU + pow(s, 2) * UU * AXG));
  _EPS1(ret, (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) * (+R * Lp + L * Rp) *
                 (+gs) * (+gs) * (+gs) * (+gs) *
                 (+Denom(-768 * pow(Pi, 2) * s * MXs * Nc + 768 * pow(Pi, 2) * s * t * Nc +
                         768 * pow(Pi, 2) * pow(s, 2) * Nc)) *
                 (+syFC1 * (1.)) *
                 (-2 * pow(MUs, 2) * UU + 2 * pow(MUs, 2) * UU * AXG + 2 * MXs * MUs * UU -
                  2 * pow(MXs, 2) * SS + 4 * pow(MXs, 2) * UU - 2 * pow(MXs, 2) * UU * AXG -
                  4 * pow(MX, 4) * UU + 2 * pow(MX, 4) * UU * AXG + 2 * t * MUs * UU -
                  4 * t * MUs * UU * AXG + 4 * t * MXs * SS - 2 * t * MXs * UU -
                  2 * pow(t, 2) * SS + 2 * pow(t, 2) * UU * AXG + 3 * s * MUs * UU -
                  3 * s * MUs * UU * AXG + 2 * s * MXs * SS - 2 * s * t * SS - s * t * UU +
                  3 * s * t * UU * AXG - pow(s, 2) * UU + pow(s, 2) * UU * AXG));

  return ret.real();
}

ComplexType ME_us_se_qq_GQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
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
  int itq = q;
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int itq = 0; itq < 3; itq++) {
    int isq = is_up_quark(itq) * 6 + itsq * 3 + itq - is_up_quark(itq) * 3;
    int iq = itq;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType Lp = conj(params->CHSQq[ch][sq][iq].R);
    ComplexType Rp = conj(params->CHSQq[ch][sq][iq].L);

    ComplexType LGp = conj(params->GLSQq[isq][q].R);
    ComplexType RGp = conj(params->GLSQq[isq][q].L);
    ComplexType LG = params->GLSQq[isq][iq].L;
    ComplexType RG = params->GLSQq[isq][iq].R;

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Denom = [](auto a) { return 1. / a; };
    //*

#define syFC1 A0(MGs)
#define syFC2 A0(MQis)
#define syFC3 B0(s, MGs, MQis)

    _EPS0(
        ret,
        (-1 + Nc) * (-1 + Nc) * (+1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
            (+Denom(-1536 * pow(Pi, 2) * pow(s, 2) * MXs * Nc +
                    1536 * pow(Pi, 2) * pow(s, 2) * t * Nc + 1536 * pow(Pi, 2) * pow(s, 3) * Nc)) *
            (+syFC1 *
                 (-2 * pow(MUs, 2) * R * RG * Lp * LGp * UU +
                  2 * pow(MUs, 2) * R * RG * Lp * LGp * UU * AXG -
                  2 * pow(MUs, 2) * L * LG * Rp * RGp * UU +
                  2 * pow(MUs, 2) * L * LG * Rp * RGp * UU * AXG +
                  2 * MXs * MUs * R * RG * Lp * LGp * UU + 2 * MXs * MUs * L * LG * Rp * RGp * UU -
                  2 * pow(MXs, 2) * R * RG * Lp * LGp * SS +
                  4 * pow(MXs, 2) * R * RG * Lp * LGp * UU -
                  2 * pow(MXs, 2) * R * RG * Lp * LGp * UU * AXG -
                  2 * pow(MXs, 2) * L * LG * Rp * RGp * SS +
                  4 * pow(MXs, 2) * L * LG * Rp * RGp * UU -
                  2 * pow(MXs, 2) * L * LG * Rp * RGp * UU * AXG -
                  4 * pow(MX, 4) * R * RG * Lp * LGp * UU +
                  2 * pow(MX, 4) * R * RG * Lp * LGp * UU * AXG -
                  4 * pow(MX, 4) * L * LG * Rp * RGp * UU +
                  2 * pow(MX, 4) * L * LG * Rp * RGp * UU * AXG +
                  2 * t * MUs * R * RG * Lp * LGp * UU -
                  4 * t * MUs * R * RG * Lp * LGp * UU * AXG +
                  2 * t * MUs * L * LG * Rp * RGp * UU -
                  4 * t * MUs * L * LG * Rp * RGp * UU * AXG +
                  4 * t * MXs * R * RG * Lp * LGp * SS - 2 * t * MXs * R * RG * Lp * LGp * UU +
                  4 * t * MXs * L * LG * Rp * RGp * SS - 2 * t * MXs * L * LG * Rp * RGp * UU -
                  2 * pow(t, 2) * R * RG * Lp * LGp * SS +
                  2 * pow(t, 2) * R * RG * Lp * LGp * UU * AXG -
                  2 * pow(t, 2) * L * LG * Rp * RGp * SS +
                  2 * pow(t, 2) * L * LG * Rp * RGp * UU * AXG +
                  3 * s * MUs * R * RG * Lp * LGp * UU -
                  3 * s * MUs * R * RG * Lp * LGp * UU * AXG +
                  3 * s * MUs * L * LG * Rp * RGp * UU -
                  3 * s * MUs * L * LG * Rp * RGp * UU * AXG +
                  2 * s * MXs * R * RG * Lp * LGp * SS + 2 * s * MXs * L * LG * Rp * RGp * SS -
                  2 * s * t * R * RG * Lp * LGp * SS) +
             syFC1 *
                 (-s * t * R * RG * Lp * LGp * UU + 3 * s * t * R * RG * Lp * LGp * UU * AXG -
                  2 * s * t * L * LG * Rp * RGp * SS - s * t * L * LG * Rp * RGp * UU +
                  3 * s * t * L * LG * Rp * RGp * UU * AXG - pow(s, 2) * R * RG * Lp * LGp * UU +
                  pow(s, 2) * R * RG * Lp * LGp * UU * AXG - pow(s, 2) * L * LG * Rp * RGp * UU +
                  pow(s, 2) * L * LG * Rp * RGp * UU * AXG) +
             syFC2 *
                 (2 * pow(MUs, 2) * R * RG * Lp * LGp * UU -
                  2 * pow(MUs, 2) * R * RG * Lp * LGp * UU * AXG +
                  2 * pow(MUs, 2) * L * LG * Rp * RGp * UU -
                  2 * pow(MUs, 2) * L * LG * Rp * RGp * UU * AXG -
                  2 * MXs * MUs * R * RG * Lp * LGp * UU - 2 * MXs * MUs * L * LG * Rp * RGp * UU +
                  2 * pow(MXs, 2) * R * RG * Lp * LGp * SS -
                  4 * pow(MXs, 2) * R * RG * Lp * LGp * UU +
                  2 * pow(MXs, 2) * R * RG * Lp * LGp * UU * AXG +
                  2 * pow(MXs, 2) * L * LG * Rp * RGp * SS -
                  4 * pow(MXs, 2) * L * LG * Rp * RGp * UU +
                  2 * pow(MXs, 2) * L * LG * Rp * RGp * UU * AXG +
                  4 * pow(MX, 4) * R * RG * Lp * LGp * UU -
                  2 * pow(MX, 4) * R * RG * Lp * LGp * UU * AXG +
                  4 * pow(MX, 4) * L * LG * Rp * RGp * UU -
                  2 * pow(MX, 4) * L * LG * Rp * RGp * UU * AXG -
                  2 * t * MUs * R * RG * Lp * LGp * UU +
                  4 * t * MUs * R * RG * Lp * LGp * UU * AXG -
                  2 * t * MUs * L * LG * Rp * RGp * UU +
                  4 * t * MUs * L * LG * Rp * RGp * UU * AXG -
                  4 * t * MXs * R * RG * Lp * LGp * SS + 2 * t * MXs * R * RG * Lp * LGp * UU -
                  4 * t * MXs * L * LG * Rp * RGp * SS + 2 * t * MXs * L * LG * Rp * RGp * UU +
                  2 * pow(t, 2) * R * RG * Lp * LGp * SS -
                  2 * pow(t, 2) * R * RG * Lp * LGp * UU * AXG +
                  2 * pow(t, 2) * L * LG * Rp * RGp * SS -
                  2 * pow(t, 2) * L * LG * Rp * RGp * UU * AXG -
                  3 * s * MUs * R * RG * Lp * LGp * UU +
                  3 * s * MUs * R * RG * Lp * LGp * UU * AXG -
                  3 * s * MUs * L * LG * Rp * RGp * UU +
                  3 * s * MUs * L * LG * Rp * RGp * UU * AXG -
                  2 * s * MXs * R * RG * Lp * LGp * SS - 2 * s * MXs * L * LG * Rp * RGp * SS +
                  2 * s * t * R * RG * Lp * LGp * SS) +
             syFC2 *
                 (s * t * R * RG * Lp * LGp * UU - 3 * s * t * R * RG * Lp * LGp * UU * AXG +
                  2 * s * t * L * LG * Rp * RGp * SS + s * t * L * LG * Rp * RGp * UU -
                  3 * s * t * L * LG * Rp * RGp * UU * AXG + pow(s, 2) * R * RG * Lp * LGp * UU -
                  pow(s, 2) * R * RG * Lp * LGp * UU * AXG + pow(s, 2) * L * LG * Rp * RGp * UU -
                  pow(s, 2) * L * LG * Rp * RGp * UU * AXG) +
             syFC3 * (-2 * pow(MUs, 2) * MQis * R * RG * Lp * LGp * UU +
                      2 * pow(MUs, 2) * MQis * R * RG * Lp * LGp * UU * AXG -
                      2 * pow(MUs, 2) * MQis * L * LG * Rp * RGp * UU +
                      2 * pow(MUs, 2) * MQis * L * LG * Rp * RGp * UU * AXG +
                      2 * MXs * MUs * MQis * R * RG * Lp * LGp * UU +
                      2 * MXs * MUs * MQis * L * LG * Rp * RGp * UU -
                      2 * pow(MXs, 2) * MQis * R * RG * Lp * LGp * SS +
                      4 * pow(MXs, 2) * MQis * R * RG * Lp * LGp * UU -
                      2 * pow(MXs, 2) * MQis * R * RG * Lp * LGp * UU * AXG -
                      2 * pow(MXs, 2) * MQis * L * LG * Rp * RGp * SS +
                      4 * pow(MXs, 2) * MQis * L * LG * Rp * RGp * UU -
                      2 * pow(MXs, 2) * MQis * L * LG * Rp * RGp * UU * AXG -
                      4 * pow(MX, 4) * MQis * R * RG * Lp * LGp * UU +
                      2 * pow(MX, 4) * MQis * R * RG * Lp * LGp * UU * AXG -
                      4 * pow(MX, 4) * MQis * L * LG * Rp * RGp * UU +
                      2 * pow(MX, 4) * MQis * L * LG * Rp * RGp * UU * AXG +
                      2 * MGs * pow(MUs, 2) * R * RG * Lp * LGp * UU -
                      2 * MGs * pow(MUs, 2) * R * RG * Lp * LGp * UU * AXG +
                      2 * MGs * pow(MUs, 2) * L * LG * Rp * RGp * UU -
                      2 * MGs * pow(MUs, 2) * L * LG * Rp * RGp * UU * AXG -
                      2 * MGs * MXs * MUs * R * RG * Lp * LGp * UU -
                      2 * MGs * MXs * MUs * L * LG * Rp * RGp * UU +
                      2 * MGs * pow(MXs, 2) * R * RG * Lp * LGp * SS -
                      4 * MGs * pow(MXs, 2) * R * RG * Lp * LGp * UU +
                      2 * MGs * pow(MXs, 2) * R * RG * Lp * LGp * UU * AXG +
                      2 * MGs * pow(MXs, 2) * L * LG * Rp * RGp * SS -
                      4 * MGs * pow(MXs, 2) * L * LG * Rp * RGp * UU +
                      2 * MGs * pow(MXs, 2) * L * LG * Rp * RGp * UU * AXG) +
             syFC3 * (4 * MGs * pow(MX, 4) * R * RG * Lp * LGp * UU -
                      2 * MGs * pow(MX, 4) * R * RG * Lp * LGp * UU * AXG +
                      4 * MGs * pow(MX, 4) * L * LG * Rp * RGp * UU -
                      2 * MGs * pow(MX, 4) * L * LG * Rp * RGp * UU * AXG +
                      2 * t * MUs * MQis * R * RG * Lp * LGp * UU -
                      4 * t * MUs * MQis * R * RG * Lp * LGp * UU * AXG +
                      2 * t * MUs * MQis * L * LG * Rp * RGp * UU -
                      4 * t * MUs * MQis * L * LG * Rp * RGp * UU * AXG +
                      4 * t * MXs * MQis * R * RG * Lp * LGp * SS -
                      2 * t * MXs * MQis * R * RG * Lp * LGp * UU +
                      4 * t * MXs * MQis * L * LG * Rp * RGp * SS -
                      2 * t * MXs * MQis * L * LG * Rp * RGp * UU -
                      2 * t * MGs * MUs * R * RG * Lp * LGp * UU +
                      4 * t * MGs * MUs * R * RG * Lp * LGp * UU * AXG -
                      2 * t * MGs * MUs * L * LG * Rp * RGp * UU +
                      4 * t * MGs * MUs * L * LG * Rp * RGp * UU * AXG -
                      4 * t * MGs * MXs * R * RG * Lp * LGp * SS +
                      2 * t * MGs * MXs * R * RG * Lp * LGp * UU -
                      4 * t * MGs * MXs * L * LG * Rp * RGp * SS +
                      2 * t * MGs * MXs * L * LG * Rp * RGp * UU -
                      2 * pow(t, 2) * MQis * R * RG * Lp * LGp * SS +
                      2 * pow(t, 2) * MQis * R * RG * Lp * LGp * UU * AXG -
                      2 * pow(t, 2) * MQis * L * LG * Rp * RGp * SS +
                      2 * pow(t, 2) * MQis * L * LG * Rp * RGp * UU * AXG +
                      2 * pow(t, 2) * MGs * R * RG * Lp * LGp * SS -
                      2 * pow(t, 2) * MGs * R * RG * Lp * LGp * UU * AXG +
                      2 * pow(t, 2) * MGs * L * LG * Rp * RGp * SS -
                      2 * pow(t, 2) * MGs * L * LG * Rp * RGp * UU * AXG +
                      3 * s * MUs * MQis * R * RG * Lp * LGp * UU -
                      3 * s * MUs * MQis * R * RG * Lp * LGp * UU * AXG +
                      3 * s * MUs * MQis * L * LG * Rp * RGp * UU) +
             syFC3 * (-3 * s * MUs * MQis * L * LG * Rp * RGp * UU * AXG +
                      2 * s * pow(MUs, 2) * R * RG * Lp * LGp * UU -
                      2 * s * pow(MUs, 2) * R * RG * Lp * LGp * UU * AXG +
                      2 * s * pow(MUs, 2) * L * LG * Rp * RGp * UU -
                      2 * s * pow(MUs, 2) * L * LG * Rp * RGp * UU * AXG +
                      2 * s * MXs * MQis * R * RG * Lp * LGp * SS +
                      2 * s * MXs * MQis * L * LG * Rp * RGp * SS -
                      2 * s * MXs * MUs * R * RG * Lp * LGp * UU -
                      2 * s * MXs * MUs * L * LG * Rp * RGp * UU +
                      2 * s * pow(MXs, 2) * R * RG * Lp * LGp * SS -
                      4 * s * pow(MXs, 2) * R * RG * Lp * LGp * UU +
                      2 * s * pow(MXs, 2) * R * RG * Lp * LGp * UU * AXG +
                      2 * s * pow(MXs, 2) * L * LG * Rp * RGp * SS -
                      4 * s * pow(MXs, 2) * L * LG * Rp * RGp * UU +
                      2 * s * pow(MXs, 2) * L * LG * Rp * RGp * UU * AXG +
                      4 * s * pow(MX, 4) * R * RG * Lp * LGp * UU -
                      2 * s * pow(MX, 4) * R * RG * Lp * LGp * UU * AXG +
                      4 * s * pow(MX, 4) * L * LG * Rp * RGp * UU -
                      2 * s * pow(MX, 4) * L * LG * Rp * RGp * UU * AXG -
                      3 * s * MGs * MUs * R * RG * Lp * LGp * UU +
                      3 * s * MGs * MUs * R * RG * Lp * LGp * UU * AXG -
                      3 * s * MGs * MUs * L * LG * Rp * RGp * UU +
                      3 * s * MGs * MUs * L * LG * Rp * RGp * UU * AXG -
                      2 * s * MGs * MXs * R * RG * Lp * LGp * SS -
                      2 * s * MGs * MXs * L * LG * Rp * RGp * SS +
                      4 * s * MG * MX * MUs * R * LG * Rp * LGp * UU -
                      4 * s * MG * MX * MUs * R * LG * Rp * LGp * UU * AXG +
                      4 * s * MG * MX * MUs * L * RG * Lp * RGp * UU -
                      4 * s * MG * MX * MUs * L * RG * Lp * RGp * UU * AXG -
                      4 * s * MG * MX * MXs * R * LG * Rp * LGp * SS -
                      4 * s * MG * MX * MXs * L * RG * Lp * RGp * SS) +
             syFC3 *
                 (-2 * s * t * MQis * R * RG * Lp * LGp * SS -
                  s * t * MQis * R * RG * Lp * LGp * UU +
                  3 * s * t * MQis * R * RG * Lp * LGp * UU * AXG -
                  2 * s * t * MQis * L * LG * Rp * RGp * SS -
                  s * t * MQis * L * LG * Rp * RGp * UU +
                  3 * s * t * MQis * L * LG * Rp * RGp * UU * AXG -
                  2 * s * t * MUs * R * RG * Lp * LGp * UU +
                  4 * s * t * MUs * R * RG * Lp * LGp * UU * AXG -
                  2 * s * t * MUs * L * LG * Rp * RGp * UU +
                  4 * s * t * MUs * L * LG * Rp * RGp * UU * AXG -
                  4 * s * t * MXs * R * RG * Lp * LGp * SS +
                  2 * s * t * MXs * R * RG * Lp * LGp * UU -
                  4 * s * t * MXs * L * LG * Rp * RGp * SS +
                  2 * s * t * MXs * L * LG * Rp * RGp * UU +
                  2 * s * t * MGs * R * RG * Lp * LGp * SS + s * t * MGs * R * RG * Lp * LGp * UU -
                  3 * s * t * MGs * R * RG * Lp * LGp * UU * AXG +
                  2 * s * t * MGs * L * LG * Rp * RGp * SS + s * t * MGs * L * LG * Rp * RGp * UU -
                  3 * s * t * MGs * L * LG * Rp * RGp * UU * AXG +
                  4 * s * t * MG * MX * R * LG * Rp * LGp * SS -
                  4 * s * t * MG * MX * R * LG * Rp * LGp * UU +
                  4 * s * t * MG * MX * R * LG * Rp * LGp * UU * AXG +
                  4 * s * t * MG * MX * L * RG * Lp * RGp * SS -
                  4 * s * t * MG * MX * L * RG * Lp * RGp * UU +
                  4 * s * t * MG * MX * L * RG * Lp * RGp * UU * AXG +
                  2 * s * pow(t, 2) * R * RG * Lp * LGp * SS -
                  2 * s * pow(t, 2) * R * RG * Lp * LGp * UU * AXG +
                  2 * s * pow(t, 2) * L * LG * Rp * RGp * SS -
                  2 * s * pow(t, 2) * L * LG * Rp * RGp * UU * AXG -
                  pow(s, 2) * MQis * R * RG * Lp * LGp * UU +
                  pow(s, 2) * MQis * R * RG * Lp * LGp * UU * AXG -
                  pow(s, 2) * MQis * L * LG * Rp * RGp * UU +
                  pow(s, 2) * MQis * L * LG * Rp * RGp * UU * AXG) +
             syFC3 *
                 (-3 * pow(s, 2) * MUs * R * RG * Lp * LGp * UU +
                  3 * pow(s, 2) * MUs * R * RG * Lp * LGp * UU * AXG -
                  3 * pow(s, 2) * MUs * L * LG * Rp * RGp * UU +
                  3 * pow(s, 2) * MUs * L * LG * Rp * RGp * UU * AXG -
                  2 * pow(s, 2) * MXs * R * RG * Lp * LGp * SS -
                  2 * pow(s, 2) * MXs * L * LG * Rp * RGp * SS +
                  pow(s, 2) * MGs * R * RG * Lp * LGp * UU -
                  pow(s, 2) * MGs * R * RG * Lp * LGp * UU * AXG +
                  pow(s, 2) * MGs * L * LG * Rp * RGp * UU -
                  pow(s, 2) * MGs * L * LG * Rp * RGp * UU * AXG +
                  4 * pow(s, 2) * MG * MX * R * LG * Rp * LGp * SS -
                  2 * pow(s, 2) * MG * MX * R * LG * Rp * LGp * UU +
                  2 * pow(s, 2) * MG * MX * R * LG * Rp * LGp * UU * AXG +
                  4 * pow(s, 2) * MG * MX * L * RG * Lp * RGp * SS -
                  2 * pow(s, 2) * MG * MX * L * RG * Lp * RGp * UU +
                  2 * pow(s, 2) * MG * MX * L * RG * Lp * RGp * UU * AXG +
                  2 * pow(s, 2) * t * R * RG * Lp * LGp * SS +
                  pow(s, 2) * t * R * RG * Lp * LGp * UU -
                  3 * pow(s, 2) * t * R * RG * Lp * LGp * UU * AXG +
                  2 * pow(s, 2) * t * L * LG * Rp * RGp * SS +
                  pow(s, 2) * t * L * LG * Rp * RGp * UU -
                  3 * pow(s, 2) * t * L * LG * Rp * RGp * UU * AXG +
                  pow(s, 3) * R * RG * Lp * LGp * UU - pow(s, 3) * R * RG * Lp * LGp * UU * AXG +
                  pow(s, 3) * L * LG * Rp * RGp * UU - pow(s, 3) * L * LG * Rp * RGp * UU * AXG)) *
            (+gs) * (+gs));
    ;
    //}
  }

  return ret.real();
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