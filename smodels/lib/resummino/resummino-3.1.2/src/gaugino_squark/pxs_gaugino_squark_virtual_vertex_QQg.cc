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

ComplexType ME_us_QQg_Qgg(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
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

  auto MQi = MQ;
  auto MQis = MQs;

#define syFC1 A0(MQs)
#define syFC2 B0(MUs, 0, MQs)
#define syFC3 B0(MUs + MXs - s - t, 0, MQs)
#define syFC4 C0(0, MUs + MXs - s - t, MUs, 0, 0, MQs)
#define syFC5 C1(0, MUs, MUs + MXs - s - t, 0, 0, MQs)
#define syFC6 C2(0, MUs, MUs + MXs - s - t, 0, 0, MQs)
#define syFC7 C00(0, MUs, MUs + MXs - s - t, 0, 0, MQs)
#define syFC8 C11(0, MUs, MUs + MXs - s - t, 0, 0, MQs)
#define syFC9 C12(0, MUs, MUs + MXs - s - t, 0, 0, MQs)
#define syFC10 C22(0, MUs, MUs + MXs - s - t, 0, 0, MQs)

  ComplexType ret = 0;
  _EPS0(
      ret,
      -(-1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
          (+Denom(3072 * pow(Pi, 2) * s * MXs * pow(MUs, 2) * MQis -
                  3072 * pow(Pi, 2) * s * MXs * pow(MUs, 3) +
                  3072 * pow(Pi, 2) * s * pow(MXs, 2) * MUs * MQis -
                  6144 * pow(Pi, 2) * s * pow(MXs, 2) * pow(MUs, 2) -
                  3072 * pow(Pi, 2) * s * pow(MXs, 3) * MUs -
                  3072 * pow(Pi, 2) * s * t * pow(MUs, 2) * MQis +
                  3072 * pow(Pi, 2) * s * t * pow(MUs, 3) -
                  6144 * pow(Pi, 2) * s * t * MXs * MUs * MQis +
                  12288 * pow(Pi, 2) * s * t * MXs * pow(MUs, 2) +
                  9216 * pow(Pi, 2) * s * t * pow(MXs, 2) * MUs +
                  3072 * pow(Pi, 2) * s * pow(t, 2) * MUs * MQis -
                  6144 * pow(Pi, 2) * s * pow(t, 2) * pow(MUs, 2) -
                  9216 * pow(Pi, 2) * s * pow(t, 2) * MXs * MUs +
                  3072 * pow(Pi, 2) * s * pow(t, 3) * MUs -
                  3072 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) * MQis +
                  3072 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 3) -
                  6144 * pow(Pi, 2) * pow(s, 2) * MXs * MUs * MQis +
                  12288 * pow(Pi, 2) * pow(s, 2) * MXs * pow(MUs, 2) +
                  9216 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) * MUs +
                  6144 * pow(Pi, 2) * pow(s, 2) * t * MUs * MQis -
                  12288 * pow(Pi, 2) * pow(s, 2) * t * pow(MUs, 2) -
                  18432 * pow(Pi, 2) * pow(s, 2) * t * MXs * MUs +
                  9216 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) * MUs +
                  3072 * pow(Pi, 2) * pow(s, 3) * MUs * MQis -
                  6144 * pow(Pi, 2) * pow(s, 3) * pow(MUs, 2) -
                  9216 * pow(Pi, 2) * pow(s, 3) * MXs * MUs +
                  9216 * pow(Pi, 2) * pow(s, 3) * t * MUs + 3072 * pow(Pi, 2) * pow(s, 4) * MUs)) *
          (+Nc) * (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 *
               (-2 * MXs * pow(MUs, 3) * SS - 2 * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                2 * pow(MXs, 3) * MUs * SS - 2 * pow(MXs, 3) * MUs * SS * AXG +
                2 * pow(MX, 4) * pow(MUs, 2) * SS + 4 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                2 * pow(MX, 4) * MXs * MUs * SS + 2 * pow(MX, 4) * MXs * MUs * SS * AXG +
                8 * pow(MX, 4) * MXs * MUs * UU * AXG + 4 * pow(MX, 4) * pow(MXs, 2) * UU * AXG -
                8 * pow(MX, 6) * MUs * UU * AXG - 8 * pow(MX, 6) * MXs * UU * AXG +
                4 * pow(MX, 8) * UU * AXG + 2 * u * pow(MUs, 3) * SS +
                2 * u * MXs * pow(MUs, 2) * SS + 4 * u * MXs * pow(MUs, 2) * SS * AXG -
                8 * u * MXs * pow(MUs, 2) * UU * AXG + 2 * u * pow(MXs, 2) * MUs * SS +
                4 * u * pow(MXs, 2) * MUs * SS * AXG - 16 * u * pow(MXs, 2) * MUs * UU * AXG +
                2 * u * pow(MXs, 3) * SS - 8 * u * pow(MXs, 3) * UU * AXG -
                6 * u * pow(MX, 4) * MUs * SS + 8 * u * pow(MX, 4) * MUs * UU * AXG -
                6 * u * pow(MX, 4) * MXs * SS + 2 * u * pow(MX, 4) * MXs * SS * AXG +
                8 * u * pow(MX, 4) * MXs * UU * AXG + 4 * u * pow(MX, 6) * SS -
                2 * u * pow(MX, 6) * SS * AXG - 4 * pow(u, 2) * pow(MUs, 2) * SS -
                2 * pow(u, 2) * pow(MUs, 2) * SS * AXG + 4 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                2 * pow(u, 2) * MXs * MUs * SS - 8 * pow(u, 2) * MXs * MUs * SS * AXG +
                16 * pow(u, 2) * MXs * MUs * UU * AXG + 2 * pow(u, 2) * pow(MXs, 2) * SS -
                4 * pow(u, 2) * pow(MXs, 2) * SS * AXG + 12 * pow(u, 2) * pow(MXs, 2) * UU * AXG) +
           syFC1 *
               (2 * pow(u, 2) * pow(MX, 4) * SS * AXG - 8 * pow(u, 2) * pow(MX, 4) * UU * AXG +
                2 * pow(u, 3) * MUs * SS + 4 * pow(u, 3) * MUs * SS * AXG -
                8 * pow(u, 3) * MUs * UU * AXG - 2 * pow(u, 3) * MXs * SS +
                4 * pow(u, 3) * MXs * SS * AXG - 8 * pow(u, 3) * MXs * UU * AXG -
                2 * pow(u, 4) * SS * AXG + 4 * pow(u, 4) * UU * AXG + s * MXs * pow(MUs, 2) * UU -
                s * MXs * pow(MUs, 2) * UU * AXG + 2 * s * pow(MXs, 2) * MUs * UU -
                2 * s * pow(MXs, 2) * MUs * UU * AXG + s * pow(MXs, 3) * UU -
                s * pow(MXs, 3) * UU * AXG - 2 * s * pow(MX, 4) * MUs * UU +
                2 * s * pow(MX, 4) * MUs * UU * AXG - 2 * s * pow(MX, 4) * MXs * UU +
                2 * s * pow(MX, 4) * MXs * UU * AXG + s * pow(MX, 6) * UU -
                s * pow(MX, 6) * UU * AXG - 2 * s * u * pow(MUs, 2) * SS -
                s * u * pow(MUs, 2) * UU + s * u * pow(MUs, 2) * UU * AXG +
                2 * s * u * MXs * MUs * SS - 2 * s * u * MXs * MUs * SS * AXG +
                4 * s * u * MXs * MUs * UU * AXG + s * u * pow(MXs, 2) * UU +
                3 * s * u * pow(MXs, 2) * UU * AXG - s * u * pow(MX, 4) * UU -
                3 * s * u * pow(MX, 4) * UU * AXG + 2 * s * pow(u, 2) * MUs * SS +
                2 * s * pow(u, 2) * MUs * SS * AXG - 4 * s * pow(u, 2) * MUs * UU * AXG -
                2 * s * pow(u, 2) * MXs * SS + 2 * s * pow(u, 2) * MXs * SS * AXG -
                s * pow(u, 2) * MXs * UU - 3 * s * pow(u, 2) * MXs * UU * AXG -
                2 * s * pow(u, 3) * SS * AXG + s * pow(u, 3) * UU + 3 * s * pow(u, 3) * UU * AXG) +
           syFC2 *
               (-2 * u * MXs * pow(MUs, 2) * MQs * SS + 18 * u * MXs * pow(MUs, 3) * SS -
                2 * u * pow(MXs, 2) * MUs * MQs * SS + 18 * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                4 * u * pow(MX, 4) * MUs * MQs * SS - 2 * u * pow(MX, 4) * MUs * MQs * SS * AXG +
                4 * u * pow(MX, 4) * MUs * MQs * UU * AXG - 36 * u * pow(MX, 4) * pow(MUs, 2) * SS +
                18 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                36 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * u * pow(MX, 4) * MXs * MQs * UU * AXG -
                36 * u * pow(MX, 4) * MXs * MUs * UU * AXG - 4 * u * pow(MX, 6) * MQs * UU * AXG +
                36 * u * pow(MX, 6) * MUs * UU * AXG + 2 * pow(u, 2) * pow(MUs, 2) * MQs * SS -
                18 * pow(u, 2) * pow(MUs, 3) * SS + 4 * pow(u, 2) * MXs * MUs * MQs * SS * AXG -
                8 * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                36 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                72 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                2 * pow(u, 2) * pow(MXs, 2) * MQs * SS -
                8 * pow(u, 2) * pow(MXs, 2) * MQs * UU * AXG -
                18 * pow(u, 2) * pow(MXs, 2) * MUs * SS +
                72 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                4 * pow(u, 2) * pow(MX, 4) * MQs * SS +
                2 * pow(u, 2) * pow(MX, 4) * MQs * SS * AXG +
                4 * pow(u, 2) * pow(MX, 4) * MQs * UU * AXG +
                36 * pow(u, 2) * pow(MX, 4) * MUs * SS -
                18 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                36 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG - 2 * pow(u, 3) * MUs * MQs * SS) +
           syFC2 *
               (-2 * pow(u, 3) * MUs * MQs * SS * AXG + 4 * pow(u, 3) * MUs * MQs * UU * AXG +
                18 * pow(u, 3) * pow(MUs, 2) * SS + 18 * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                36 * pow(u, 3) * pow(MUs, 2) * UU * AXG + 2 * pow(u, 3) * MXs * MQs * SS -
                4 * pow(u, 3) * MXs * MQs * SS * AXG + 8 * pow(u, 3) * MXs * MQs * UU * AXG -
                18 * pow(u, 3) * MXs * MUs * SS + 36 * pow(u, 3) * MXs * MUs * SS * AXG -
                72 * pow(u, 3) * MXs * MUs * UU * AXG + 2 * pow(u, 4) * MQs * SS * AXG -
                4 * pow(u, 4) * MQs * UU * AXG - 18 * pow(u, 4) * MUs * SS * AXG +
                36 * pow(u, 4) * MUs * UU * AXG - 2 * s * u * MXs * MUs * MQs * SS +
                2 * s * u * MXs * MUs * MQs * SS * AXG + 3 * s * u * MXs * MUs * MQs * UU -
                3 * s * u * MXs * MUs * MQs * UU * AXG + 4 * s * u * MXs * pow(MUs, 2) * SS -
                4 * s * u * MXs * pow(MUs, 2) * SS * AXG - 13 * s * u * MXs * pow(MUs, 2) * UU +
                13 * s * u * MXs * pow(MUs, 2) * UU * AXG + 3 * s * u * pow(MXs, 2) * MQs * UU -
                3 * s * u * pow(MXs, 2) * MQs * UU * AXG - 13 * s * u * pow(MXs, 2) * MUs * UU +
                13 * s * u * pow(MXs, 2) * MUs * UU * AXG - 3 * s * u * pow(MX, 4) * MQs * UU +
                3 * s * u * pow(MX, 4) * MQs * UU * AXG + 13 * s * u * pow(MX, 4) * MUs * UU -
                13 * s * u * pow(MX, 4) * MUs * UU * AXG -
                2 * s * pow(u, 2) * MUs * MQs * SS * AXG - 3 * s * pow(u, 2) * MUs * MQs * UU +
                3 * s * pow(u, 2) * MUs * MQs * UU * AXG + 14 * s * pow(u, 2) * pow(MUs, 2) * SS +
                4 * s * pow(u, 2) * pow(MUs, 2) * SS * AXG) +
           syFC2 *
               (13 * s * pow(u, 2) * pow(MUs, 2) * UU -
                13 * s * pow(u, 2) * pow(MUs, 2) * UU * AXG + 2 * s * pow(u, 2) * MXs * MQs * SS -
                2 * s * pow(u, 2) * MXs * MQs * SS * AXG + s * pow(u, 2) * MXs * MQs * UU +
                3 * s * pow(u, 2) * MXs * MQs * UU * AXG - 4 * s * pow(u, 2) * MXs * MUs * SS +
                4 * s * pow(u, 2) * MXs * MUs * SS * AXG - 23 * s * pow(u, 2) * MXs * MUs * UU -
                13 * s * pow(u, 2) * MXs * MUs * UU * AXG + 2 * s * pow(u, 3) * MQs * SS * AXG -
                s * pow(u, 3) * MQs * UU - 3 * s * pow(u, 3) * MQs * UU * AXG -
                14 * s * pow(u, 3) * MUs * SS - 4 * s * pow(u, 3) * MUs * SS * AXG +
                23 * s * pow(u, 3) * MUs * UU + 13 * s * pow(u, 3) * MUs * UU * AXG) +
           syFC3 *
               (2 * MXs * pow(MUs, 3) * MQs * SS + 2 * pow(MXs, 2) * pow(MUs, 2) * MQs * SS -
                4 * pow(MX, 4) * pow(MUs, 2) * MQs * SS +
                2 * pow(MX, 4) * pow(MUs, 2) * MQs * SS * AXG -
                4 * pow(MX, 4) * pow(MUs, 2) * MQs * UU * AXG -
                4 * pow(MX, 4) * MXs * MUs * MQs * UU * AXG +
                4 * pow(MX, 6) * MUs * MQs * UU * AXG - 2 * u * pow(MUs, 3) * MQs * SS -
                4 * u * MXs * pow(MUs, 2) * MQs * SS * AXG +
                8 * u * MXs * pow(MUs, 2) * MQs * UU * AXG - 2 * u * MXs * pow(MUs, 3) * SS -
                2 * u * pow(MXs, 2) * MUs * MQs * SS + 8 * u * pow(MXs, 2) * MUs * MQs * UU * AXG -
                2 * u * pow(MXs, 2) * pow(MUs, 2) * SS + 4 * u * pow(MX, 4) * MUs * MQs * SS -
                2 * u * pow(MX, 4) * MUs * MQs * SS * AXG -
                4 * u * pow(MX, 4) * MUs * MQs * UU * AXG + 4 * u * pow(MX, 4) * pow(MUs, 2) * SS -
                2 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                4 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * u * pow(MX, 4) * MXs * MUs * UU * AXG - 4 * u * pow(MX, 6) * MUs * UU * AXG +
                2 * pow(u, 2) * pow(MUs, 2) * MQs * SS +
                2 * pow(u, 2) * pow(MUs, 2) * MQs * SS * AXG -
                4 * pow(u, 2) * pow(MUs, 2) * MQs * UU * AXG + 2 * pow(u, 2) * pow(MUs, 3) * SS -
                2 * pow(u, 2) * MXs * MUs * MQs * SS + 4 * pow(u, 2) * MXs * MUs * MQs * SS * AXG -
                8 * pow(u, 2) * MXs * MUs * MQs * UU * AXG +
                4 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                8 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                2 * pow(u, 2) * pow(MXs, 2) * MUs * SS) +
           syFC3 * (-8 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                    4 * pow(u, 2) * pow(MX, 4) * MUs * SS +
                    2 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG +
                    4 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                    2 * pow(u, 3) * MUs * MQs * SS * AXG + 4 * pow(u, 3) * MUs * MQs * UU * AXG -
                    2 * pow(u, 3) * pow(MUs, 2) * SS - 2 * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                    4 * pow(u, 3) * pow(MUs, 2) * UU * AXG + 2 * pow(u, 3) * MXs * MUs * SS -
                    4 * pow(u, 3) * MXs * MUs * SS * AXG + 8 * pow(u, 3) * MXs * MUs * UU * AXG +
                    2 * pow(u, 4) * MUs * SS * AXG - 4 * pow(u, 4) * MUs * UU * AXG -
                    s * MXs * pow(MUs, 2) * MQs * UU + s * MXs * pow(MUs, 2) * MQs * UU * AXG -
                    s * pow(MXs, 2) * MUs * MQs * UU + s * pow(MXs, 2) * MUs * MQs * UU * AXG +
                    s * pow(MX, 4) * MUs * MQs * UU - s * pow(MX, 4) * MUs * MQs * UU * AXG +
                    2 * s * u * pow(MUs, 2) * MQs * SS + s * u * pow(MUs, 2) * MQs * UU -
                    s * u * pow(MUs, 2) * MQs * UU * AXG - 3 * s * u * MXs * MUs * MQs * UU -
                    s * u * MXs * MUs * MQs * UU * AXG + 2 * s * u * MXs * pow(MUs, 2) * SS -
                    2 * s * u * MXs * pow(MUs, 2) * SS * AXG - s * u * MXs * pow(MUs, 2) * UU +
                    s * u * MXs * pow(MUs, 2) * UU * AXG - s * u * pow(MXs, 2) * MUs * UU +
                    s * u * pow(MXs, 2) * MUs * UU * AXG + s * u * pow(MX, 4) * MUs * UU -
                    s * u * pow(MX, 4) * MUs * UU * AXG - 2 * s * pow(u, 2) * MUs * MQs * SS +
                    3 * s * pow(u, 2) * MUs * MQs * UU + s * pow(u, 2) * MUs * MQs * UU * AXG) +
           syFC3 *
               (-4 * s * pow(u, 2) * pow(MUs, 2) * SS + 2 * s * pow(u, 2) * pow(MUs, 2) * SS * AXG +
                s * pow(u, 2) * pow(MUs, 2) * UU - s * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                2 * s * pow(u, 2) * MXs * MUs * SS + 2 * s * pow(u, 2) * MXs * MUs * SS * AXG +
                5 * s * pow(u, 2) * MXs * MUs * UU - s * pow(u, 2) * MXs * MUs * UU * AXG +
                4 * s * pow(u, 3) * MUs * SS - 2 * s * pow(u, 3) * MUs * SS * AXG -
                5 * s * pow(u, 3) * MUs * UU + s * pow(u, 3) * MUs * UU * AXG) +
           syFC4 *
               (-8 * u * MXs * pow(MUs, 4) * SS - 16 * u * pow(MXs, 2) * pow(MUs, 3) * SS -
                8 * u * pow(MXs, 3) * pow(MUs, 2) * SS + 24 * u * pow(MX, 4) * pow(MUs, 3) * SS -
                8 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG +
                16 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG +
                24 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS -
                8 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG +
                32 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG +
                16 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG -
                16 * u * pow(MX, 6) * pow(MUs, 2) * SS +
                8 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG -
                32 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG -
                32 * u * pow(MX, 6) * MXs * MUs * UU * AXG + 16 * u * pow(MX, 8) * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MUs, 4) * SS + 8 * pow(u, 2) * MXs * pow(MUs, 3) * SS +
                16 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG -
                32 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS +
                16 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                64 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 3) * MUs * SS -
                32 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                24 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS +
                32 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                24 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS +
                8 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG) +
           syFC4 *
               (32 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG +
                16 * pow(u, 2) * pow(MX, 6) * MUs * SS -
                8 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG - 16 * pow(u, 3) * pow(MUs, 3) * SS -
                8 * pow(u, 3) * pow(MUs, 3) * SS * AXG + 16 * pow(u, 3) * pow(MUs, 3) * UU * AXG +
                8 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                32 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG +
                64 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 3) * pow(MXs, 2) * MUs * SS -
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                48 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG +
                8 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG -
                32 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG + 8 * pow(u, 4) * pow(MUs, 2) * SS +
                16 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 32 * pow(u, 4) * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 4) * MXs * MUs * SS + 16 * pow(u, 4) * MXs * MUs * SS * AXG -
                32 * pow(u, 4) * MXs * MUs * UU * AXG - 8 * pow(u, 5) * MUs * SS * AXG +
                16 * pow(u, 5) * MUs * UU * AXG - 4 * s * u * MXs * pow(MUs, 3) * SS +
                4 * s * u * MXs * pow(MUs, 3) * SS * AXG + 8 * s * u * MXs * pow(MUs, 3) * UU -
                8 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                8 * s * u * pow(MXs, 3) * MUs * UU) +
           syFC4 *
               (-8 * s * u * pow(MXs, 3) * MUs * UU * AXG +
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS -
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                16 * s * u * pow(MX, 4) * MXs * MUs * UU +
                16 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 8 * s * u * pow(MX, 6) * MUs * UU -
                8 * s * u * pow(MX, 6) * MUs * UU * AXG - 4 * s * pow(u, 2) * pow(MUs, 3) * SS -
                4 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG - 8 * s * pow(u, 2) * pow(MUs, 3) * UU +
                8 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                24 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS -
                24 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                32 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS -
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                24 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS +
                4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG +
                8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                24 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                8 * s * pow(u, 3) * pow(MUs, 2) * SS + 24 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 3) * pow(MUs, 2) * UU -
                32 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG) +
           syFC4 *
               (-20 * s * pow(u, 3) * MXs * MUs * SS + 20 * s * pow(u, 3) * MXs * MUs * SS * AXG +
                8 * s * pow(u, 3) * MXs * MUs * UU - 24 * s * pow(u, 3) * MXs * MUs * UU * AXG +
                12 * s * pow(u, 4) * MUs * SS - 20 * s * pow(u, 4) * MUs * SS * AXG -
                8 * s * pow(u, 4) * MUs * UU + 24 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC5 *
               (-16 * u * MXs * pow(MUs, 4) * SS - 32 * u * pow(MXs, 2) * pow(MUs, 3) * SS -
                16 * u * pow(MXs, 3) * pow(MUs, 2) * SS + 48 * u * pow(MX, 4) * pow(MUs, 3) * SS -
                16 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG +
                32 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG +
                48 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS -
                16 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG +
                64 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG +
                32 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG -
                32 * u * pow(MX, 6) * pow(MUs, 2) * SS +
                16 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG -
                64 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG -
                64 * u * pow(MX, 6) * MXs * MUs * UU * AXG + 32 * u * pow(MX, 8) * MUs * UU * AXG +
                16 * pow(u, 2) * pow(MUs, 4) * SS + 16 * pow(u, 2) * MXs * pow(MUs, 3) * SS +
                32 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG -
                64 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                16 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS +
                32 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                128 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 2) * pow(MXs, 3) * MUs * SS -
                64 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                48 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS +
                64 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                48 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS +
                16 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG) +
           syFC5 *
               (64 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG +
                32 * pow(u, 2) * pow(MX, 6) * MUs * SS -
                16 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG - 32 * pow(u, 3) * pow(MUs, 3) * SS -
                16 * pow(u, 3) * pow(MUs, 3) * SS * AXG + 32 * pow(u, 3) * pow(MUs, 3) * UU * AXG +
                16 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                64 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG +
                128 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS -
                32 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                96 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG +
                16 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG -
                64 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG + 16 * pow(u, 4) * pow(MUs, 2) * SS +
                32 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 64 * pow(u, 4) * pow(MUs, 2) * UU * AXG -
                16 * pow(u, 4) * MXs * MUs * SS + 32 * pow(u, 4) * MXs * MUs * SS * AXG -
                64 * pow(u, 4) * MXs * MUs * UU * AXG - 16 * pow(u, 5) * MUs * SS * AXG +
                32 * pow(u, 5) * MUs * UU * AXG + 8 * s * u * MXs * pow(MUs, 3) * UU -
                8 * s * u * MXs * pow(MUs, 3) * UU * AXG +
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                8 * s * u * pow(MXs, 3) * MUs * UU - 8 * s * u * pow(MXs, 3) * MUs * UU * AXG -
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                16 * s * u * pow(MX, 4) * MXs * MUs * UU) +
           syFC5 *
               (16 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 8 * s * u * pow(MX, 6) * MUs * UU -
                8 * s * u * pow(MX, 6) * MUs * UU * AXG - 16 * s * pow(u, 2) * pow(MUs, 3) * SS -
                8 * s * pow(u, 2) * pow(MUs, 3) * UU + 8 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                32 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS -
                32 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                48 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                40 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                40 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG +
                32 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 3) * pow(MUs, 2) * UU -
                48 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG - 32 * s * pow(u, 3) * MXs * MUs * SS +
                32 * s * pow(u, 3) * MXs * MUs * SS * AXG + 8 * s * pow(u, 3) * MXs * MUs * UU -
                40 * s * pow(u, 3) * MXs * MUs * UU * AXG + 16 * s * pow(u, 4) * MUs * SS -
                32 * s * pow(u, 4) * MUs * SS * AXG - 8 * s * pow(u, 4) * MUs * UU +
                40 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC6 *
               (-12 * u * MXs * pow(MUs, 4) * SS - 24 * u * pow(MXs, 2) * pow(MUs, 3) * SS -
                12 * u * pow(MXs, 3) * pow(MUs, 2) * SS + 36 * u * pow(MX, 4) * pow(MUs, 3) * SS -
                12 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG +
                24 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG +
                36 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS -
                12 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG +
                48 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG +
                24 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG -
                24 * u * pow(MX, 6) * pow(MUs, 2) * SS +
                12 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG -
                48 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG -
                48 * u * pow(MX, 6) * MXs * MUs * UU * AXG + 24 * u * pow(MX, 8) * MUs * UU * AXG +
                12 * pow(u, 2) * pow(MUs, 4) * SS + 12 * pow(u, 2) * MXs * pow(MUs, 3) * SS +
                24 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG -
                48 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                12 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS +
                24 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                96 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                12 * pow(u, 2) * pow(MXs, 3) * MUs * SS -
                48 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                36 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS +
                48 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                36 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS +
                12 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG) +
           syFC6 *
               (48 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG +
                24 * pow(u, 2) * pow(MX, 6) * MUs * SS -
                12 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG - 24 * pow(u, 3) * pow(MUs, 3) * SS -
                12 * pow(u, 3) * pow(MUs, 3) * SS * AXG + 24 * pow(u, 3) * pow(MUs, 3) * UU * AXG +
                12 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                48 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG +
                96 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG +
                12 * pow(u, 3) * pow(MXs, 2) * MUs * SS -
                24 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                72 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG +
                12 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG -
                48 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG + 12 * pow(u, 4) * pow(MUs, 2) * SS +
                24 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 48 * pow(u, 4) * pow(MUs, 2) * UU * AXG -
                12 * pow(u, 4) * MXs * MUs * SS + 24 * pow(u, 4) * MXs * MUs * SS * AXG -
                48 * pow(u, 4) * MXs * MUs * UU * AXG - 12 * pow(u, 5) * MUs * SS * AXG +
                24 * pow(u, 5) * MUs * UU * AXG + 6 * s * u * MXs * pow(MUs, 3) * UU -
                6 * s * u * MXs * pow(MUs, 3) * UU * AXG +
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                6 * s * u * pow(MXs, 3) * MUs * UU - 6 * s * u * pow(MXs, 3) * MUs * UU * AXG -
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                12 * s * u * pow(MX, 4) * MXs * MUs * UU) +
           syFC6 *
               (12 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 6 * s * u * pow(MX, 6) * MUs * UU -
                6 * s * u * pow(MX, 6) * MUs * UU * AXG - 12 * s * pow(u, 2) * pow(MUs, 3) * SS -
                6 * s * pow(u, 2) * pow(MUs, 3) * UU + 6 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                24 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS -
                24 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                12 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                36 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                6 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                30 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                6 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                30 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG +
                24 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                12 * s * pow(u, 3) * pow(MUs, 2) * UU -
                36 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG - 24 * s * pow(u, 3) * MXs * MUs * SS +
                24 * s * pow(u, 3) * MXs * MUs * SS * AXG + 6 * s * pow(u, 3) * MXs * MUs * UU -
                30 * s * pow(u, 3) * MXs * MUs * UU * AXG + 12 * s * pow(u, 4) * MUs * SS -
                24 * s * pow(u, 4) * MUs * SS * AXG - 6 * s * pow(u, 4) * MUs * UU +
                30 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC7 *
               (-16 * u * MXs * pow(MUs, 3) * SS - 16 * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                32 * u * pow(MX, 4) * pow(MUs, 2) * SS -
                16 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                32 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                32 * u * pow(MX, 4) * MXs * MUs * UU * AXG - 32 * u * pow(MX, 6) * MUs * UU * AXG +
                16 * pow(u, 2) * pow(MUs, 3) * SS + 32 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                64 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 2) * pow(MXs, 2) * MUs * SS -
                64 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                32 * pow(u, 2) * pow(MX, 4) * MUs * SS +
                16 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG - 16 * pow(u, 3) * pow(MUs, 2) * SS -
                16 * pow(u, 3) * pow(MUs, 2) * SS * AXG + 32 * pow(u, 3) * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 3) * MXs * MUs * SS - 32 * pow(u, 3) * MXs * MUs * SS * AXG +
                64 * pow(u, 3) * MXs * MUs * UU * AXG + 16 * pow(u, 4) * MUs * SS * AXG -
                32 * pow(u, 4) * MUs * UU * AXG + 8 * s * u * MXs * pow(MUs, 2) * UU -
                8 * s * u * MXs * pow(MUs, 2) * UU * AXG + 8 * s * u * pow(MXs, 2) * MUs * UU -
                8 * s * u * pow(MXs, 2) * MUs * UU * AXG - 8 * s * u * pow(MX, 4) * MUs * UU +
                8 * s * u * pow(MX, 4) * MUs * UU * AXG - 16 * s * pow(u, 2) * pow(MUs, 2) * SS -
                8 * s * pow(u, 2) * pow(MUs, 2) * UU + 8 * s * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                24 * s * pow(u, 2) * MXs * MUs * UU) +
           syFC7 * (8 * s * pow(u, 2) * MXs * MUs * UU * AXG + 16 * s * pow(u, 3) * MUs * SS -
                    24 * s * pow(u, 3) * MUs * UU - 8 * s * pow(u, 3) * MUs * UU * AXG) +
           syFC8 *
               (-8 * s * u * MXs * pow(MUs, 3) * SS + 8 * s * u * MXs * pow(MUs, 3) * SS * AXG +
                8 * s * u * MXs * pow(MUs, 3) * UU - 8 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                8 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                8 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                16 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                8 * s * u * pow(MXs, 3) * MUs * UU - 8 * s * u * pow(MXs, 3) * MUs * UU * AXG +
                8 * s * u * pow(MX, 4) * pow(MUs, 2) * SS -
                8 * s * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                16 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                16 * s * u * pow(MX, 4) * MXs * MUs * UU +
                16 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 8 * s * u * pow(MX, 6) * MUs * UU -
                8 * s * u * pow(MX, 6) * MUs * UU * AXG + 8 * s * pow(u, 2) * pow(MUs, 3) * SS -
                8 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG - 8 * s * pow(u, 2) * pow(MUs, 3) * UU +
                8 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG) +
           syFC8 * (-8 * s * pow(u, 2) * pow(MX, 4) * MUs * SS +
                    8 * s * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG +
                    8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                    8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                    16 * s * pow(u, 3) * pow(MUs, 2) * SS +
                    16 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                    16 * s * pow(u, 3) * pow(MUs, 2) * UU -
                    16 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG -
                    8 * s * pow(u, 3) * MXs * MUs * SS + 8 * s * pow(u, 3) * MXs * MUs * SS * AXG +
                    8 * s * pow(u, 3) * MXs * MUs * UU - 8 * s * pow(u, 3) * MXs * MUs * UU * AXG +
                    8 * s * pow(u, 4) * MUs * SS - 8 * s * pow(u, 4) * MUs * SS * AXG -
                    8 * s * pow(u, 4) * MUs * UU + 8 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC9 *
               (8 * u * MXs * pow(MUs, 4) * SS + 16 * u * pow(MXs, 2) * pow(MUs, 3) * SS +
                8 * u * pow(MXs, 3) * pow(MUs, 2) * SS - 24 * u * pow(MX, 4) * pow(MUs, 3) * SS +
                8 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG -
                16 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG -
                24 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS +
                8 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG -
                32 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG -
                16 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG +
                16 * u * pow(MX, 6) * pow(MUs, 2) * SS -
                8 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG +
                32 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG +
                32 * u * pow(MX, 6) * MXs * MUs * UU * AXG - 16 * u * pow(MX, 8) * MUs * UU * AXG -
                8 * pow(u, 2) * pow(MUs, 4) * SS - 8 * pow(u, 2) * MXs * pow(MUs, 3) * SS -
                16 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG +
                32 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG -
                8 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS -
                16 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                64 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 2) * pow(MXs, 3) * MUs * SS +
                32 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG +
                24 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS -
                32 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                24 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS -
                8 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG) +
           syFC9 *
               (-32 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG -
                16 * pow(u, 2) * pow(MX, 6) * MUs * SS +
                8 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG + 16 * pow(u, 3) * pow(MUs, 3) * SS +
                8 * pow(u, 3) * pow(MUs, 3) * SS * AXG - 16 * pow(u, 3) * pow(MUs, 3) * UU * AXG -
                8 * pow(u, 3) * MXs * pow(MUs, 2) * SS +
                32 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG -
                64 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 3) * pow(MXs, 2) * MUs * SS +
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG -
                48 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG -
                8 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG +
                32 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG - 8 * pow(u, 4) * pow(MUs, 2) * SS -
                16 * pow(u, 4) * pow(MUs, 2) * SS * AXG + 32 * pow(u, 4) * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 4) * MXs * MUs * SS - 16 * pow(u, 4) * MXs * MUs * SS * AXG +
                32 * pow(u, 4) * MXs * MUs * UU * AXG + 8 * pow(u, 5) * MUs * SS * AXG -
                16 * pow(u, 5) * MUs * UU * AXG - 4 * s * u * MXs * pow(MUs, 3) * UU +
                4 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                8 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU +
                8 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG -
                4 * s * u * pow(MXs, 3) * MUs * UU + 4 * s * u * pow(MXs, 3) * MUs * UU * AXG +
                8 * s * u * pow(MX, 4) * pow(MUs, 2) * UU -
                8 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                8 * s * u * pow(MX, 4) * MXs * MUs * UU) +
           syFC9 *
               (-8 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG - 4 * s * u * pow(MX, 6) * MUs * UU +
                4 * s * u * pow(MX, 6) * MUs * UU * AXG + 8 * s * pow(u, 2) * pow(MUs, 3) * SS +
                4 * s * pow(u, 2) * pow(MUs, 3) * UU - 4 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG -
                24 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                8 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                28 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                12 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                28 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                12 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                32 * s * pow(u, 3) * pow(MUs, 2) * SS +
                16 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG +
                24 * s * pow(u, 3) * pow(MUs, 2) * UU - 8 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG -
                16 * s * pow(u, 3) * MXs * MUs * SS + 16 * s * pow(u, 3) * MXs * MUs * SS * AXG +
                28 * s * pow(u, 3) * MXs * MUs * UU - 12 * s * pow(u, 3) * MXs * MUs * UU * AXG +
                24 * s * pow(u, 4) * MUs * SS - 16 * s * pow(u, 4) * MUs * SS * AXG -
                28 * s * pow(u, 4) * MUs * UU + 12 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC10 *
               (-16 * pow(u, 2) * MXs * pow(MUs, 3) * SS -
                16 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS +
                32 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS -
                16 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                32 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG -
                32 * pow(u, 2) * pow(MX, 6) * MUs * UU * AXG + 16 * pow(u, 3) * pow(MUs, 3) * SS +
                32 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG -
                64 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS -
                64 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG -
                32 * pow(u, 3) * pow(MX, 4) * MUs * SS +
                16 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG +
                32 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG - 16 * pow(u, 4) * pow(MUs, 2) * SS -
                16 * pow(u, 4) * pow(MUs, 2) * SS * AXG + 32 * pow(u, 4) * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 4) * MXs * MUs * SS - 32 * pow(u, 4) * MXs * MUs * SS * AXG +
                64 * pow(u, 4) * MXs * MUs * UU * AXG + 16 * pow(u, 5) * MUs * SS * AXG -
                32 * pow(u, 5) * MUs * UU * AXG + 8 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU -
                8 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU +
                8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG) +
           syFC10 *
               (-16 * s * pow(u, 3) * pow(MUs, 2) * SS - 8 * s * pow(u, 3) * pow(MUs, 2) * UU +
                8 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG + 24 * s * pow(u, 3) * MXs * MUs * UU +
                8 * s * pow(u, 3) * MXs * MUs * UU * AXG + 16 * s * pow(u, 4) * MUs * SS -
                24 * s * pow(u, 4) * MUs * UU - 8 * s * pow(u, 4) * MUs * UU * AXG)));

  return ret.real();
}

ComplexType ME_us_QQg_gQQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
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

  auto MQi = MQ;
  auto MQis = MQs;

#define syFC1 A0(MQs)
#define syFC2 B0(MUs, 0, MQs)
#define syFC3 B0(MUs + MXs - s - t, 0, MQs)
#define syFC4 C0(0, MUs + MXs - s - t, MUs, MQs, MQs, 0)
#define syFC5 C1(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)
#define syFC6 C2(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)
#define syFC7 C00(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)
#define syFC8 C11(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)
#define syFC9 C12(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)
#define syFC10 C22(0, MUs, MUs + MXs - s - t, MQs, MQs, 0)

  ComplexType ret = 0;
  _EPS0(
      ret,
      -(-1) * (-1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
          (+Denom(3072 * pow(Pi, 2) * s * MXs * pow(MUs, 2) * MQis * Nc -
                  3072 * pow(Pi, 2) * s * MXs * pow(MUs, 3) * Nc +
                  3072 * pow(Pi, 2) * s * pow(MXs, 2) * MUs * MQis * Nc -
                  6144 * pow(Pi, 2) * s * pow(MXs, 2) * pow(MUs, 2) * Nc -
                  3072 * pow(Pi, 2) * s * pow(MXs, 3) * MUs * Nc -
                  3072 * pow(Pi, 2) * s * t * pow(MUs, 2) * MQis * Nc +
                  3072 * pow(Pi, 2) * s * t * pow(MUs, 3) * Nc -
                  6144 * pow(Pi, 2) * s * t * MXs * MUs * MQis * Nc +
                  12288 * pow(Pi, 2) * s * t * MXs * pow(MUs, 2) * Nc +
                  9216 * pow(Pi, 2) * s * t * pow(MXs, 2) * MUs * Nc +
                  3072 * pow(Pi, 2) * s * pow(t, 2) * MUs * MQis * Nc -
                  6144 * pow(Pi, 2) * s * pow(t, 2) * pow(MUs, 2) * Nc -
                  9216 * pow(Pi, 2) * s * pow(t, 2) * MXs * MUs * Nc +
                  3072 * pow(Pi, 2) * s * pow(t, 3) * MUs * Nc -
                  3072 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) * MQis * Nc +
                  3072 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 3) * Nc -
                  6144 * pow(Pi, 2) * pow(s, 2) * MXs * MUs * MQis * Nc +
                  12288 * pow(Pi, 2) * pow(s, 2) * MXs * pow(MUs, 2) * Nc +
                  9216 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) * MUs * Nc +
                  6144 * pow(Pi, 2) * pow(s, 2) * t * MUs * MQis * Nc -
                  12288 * pow(Pi, 2) * pow(s, 2) * t * pow(MUs, 2) * Nc -
                  18432 * pow(Pi, 2) * pow(s, 2) * t * MXs * MUs * Nc +
                  9216 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) * MUs * Nc +
                  3072 * pow(Pi, 2) * pow(s, 3) * MUs * MQis * Nc -
                  6144 * pow(Pi, 2) * pow(s, 3) * pow(MUs, 2) * Nc -
                  9216 * pow(Pi, 2) * pow(s, 3) * MXs * MUs * Nc +
                  9216 * pow(Pi, 2) * pow(s, 3) * t * MUs * Nc +
                  3072 * pow(Pi, 2) * pow(s, 4) * MUs * Nc)) *
          (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 *
               (-2 * MXs * pow(MUs, 3) * SS - 4 * pow(MXs, 2) * pow(MUs, 2) * SS -
                2 * pow(MXs, 3) * MUs * SS + 6 * pow(MX, 4) * pow(MUs, 2) * SS -
                2 * pow(MX, 4) * pow(MUs, 2) * SS * AXG + 4 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                6 * pow(MX, 4) * MXs * MUs * SS - 2 * pow(MX, 4) * MXs * MUs * SS * AXG +
                8 * pow(MX, 4) * MXs * MUs * UU * AXG + 4 * pow(MX, 4) * pow(MXs, 2) * UU * AXG -
                4 * pow(MX, 6) * MUs * SS + 2 * pow(MX, 6) * MUs * SS * AXG -
                8 * pow(MX, 6) * MUs * UU * AXG - 8 * pow(MX, 6) * MXs * UU * AXG +
                4 * pow(MX, 8) * UU * AXG + 2 * u * pow(MUs, 3) * SS -
                2 * u * MXs * pow(MUs, 2) * SS + 4 * u * MXs * pow(MUs, 2) * SS * AXG -
                8 * u * MXs * pow(MUs, 2) * UU * AXG - 2 * u * pow(MXs, 2) * MUs * SS +
                4 * u * pow(MXs, 2) * MUs * SS * AXG - 16 * u * pow(MXs, 2) * MUs * UU * AXG +
                2 * u * pow(MXs, 3) * SS - 8 * u * pow(MXs, 3) * UU * AXG +
                2 * u * pow(MX, 4) * MUs * SS - 4 * u * pow(MX, 4) * MUs * SS * AXG +
                16 * u * pow(MX, 4) * MUs * UU * AXG - 6 * u * pow(MX, 4) * MXs * SS +
                2 * u * pow(MX, 4) * MXs * SS * AXG + 16 * u * pow(MX, 4) * MXs * UU * AXG +
                4 * u * pow(MX, 6) * SS - 2 * u * pow(MX, 6) * SS * AXG -
                8 * u * pow(MX, 6) * UU * AXG - 2 * pow(u, 2) * pow(MUs, 2) * SS * AXG +
                4 * pow(u, 2) * pow(MUs, 2) * UU * AXG + 2 * pow(u, 2) * MXs * MUs * SS +
                6 * pow(u, 2) * pow(MXs, 2) * SS - 4 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                4 * pow(u, 2) * pow(MXs, 2) * UU * AXG) +
           syFC1 *
               (-8 * pow(u, 2) * pow(MX, 4) * SS + 6 * pow(u, 2) * pow(MX, 4) * SS * AXG -
                2 * pow(u, 3) * MUs * SS + 2 * pow(u, 3) * MXs * SS -
                4 * pow(u, 3) * MXs * SS * AXG + 8 * pow(u, 3) * MXs * UU * AXG +
                2 * pow(u, 4) * SS * AXG - 4 * pow(u, 4) * UU * AXG + s * MXs * pow(MUs, 2) * UU -
                s * MXs * pow(MUs, 2) * UU * AXG + 2 * s * pow(MXs, 2) * MUs * UU -
                2 * s * pow(MXs, 2) * MUs * UU * AXG + s * pow(MXs, 3) * UU -
                s * pow(MXs, 3) * UU * AXG - 2 * s * pow(MX, 4) * MUs * UU +
                2 * s * pow(MX, 4) * MUs * UU * AXG - 2 * s * pow(MX, 4) * MXs * UU +
                2 * s * pow(MX, 4) * MXs * UU * AXG + s * pow(MX, 6) * UU -
                s * pow(MX, 6) * UU * AXG - 2 * s * u * pow(MUs, 2) * SS -
                s * u * pow(MUs, 2) * UU + s * u * pow(MUs, 2) * UU * AXG -
                2 * s * u * MXs * MUs * SS + 2 * s * u * MXs * MUs * SS * AXG +
                6 * s * u * MXs * MUs * UU - 2 * s * u * MXs * MUs * UU * AXG +
                7 * s * u * pow(MXs, 2) * UU - 3 * s * u * pow(MXs, 2) * UU * AXG -
                7 * s * u * pow(MX, 4) * UU + 3 * s * u * pow(MX, 4) * UU * AXG +
                2 * s * pow(u, 2) * MUs * SS - 2 * s * pow(u, 2) * MUs * SS * AXG -
                6 * s * pow(u, 2) * MUs * UU + 2 * s * pow(u, 2) * MUs * UU * AXG +
                2 * s * pow(u, 2) * MXs * SS - 2 * s * pow(u, 2) * MXs * SS * AXG +
                s * pow(u, 2) * MXs * UU + 3 * s * pow(u, 2) * MXs * UU * AXG +
                2 * s * pow(u, 3) * SS * AXG - s * pow(u, 3) * UU - 3 * s * pow(u, 3) * UU * AXG) +
           syFC2 *
               (2 * u * MXs * pow(MUs, 2) * MQs * SS + 2 * u * MXs * pow(MUs, 3) * SS +
                2 * u * pow(MXs, 2) * MUs * MQs * SS + 2 * u * pow(MXs, 2) * pow(MUs, 2) * SS -
                4 * u * pow(MX, 4) * MUs * MQs * SS + 2 * u * pow(MX, 4) * MUs * MQs * SS * AXG -
                4 * u * pow(MX, 4) * MUs * MQs * UU * AXG - 4 * u * pow(MX, 4) * pow(MUs, 2) * SS +
                2 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                4 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MX, 4) * MXs * MQs * UU * AXG -
                4 * u * pow(MX, 4) * MXs * MUs * UU * AXG + 4 * u * pow(MX, 6) * MQs * UU * AXG +
                4 * u * pow(MX, 6) * MUs * UU * AXG - 2 * pow(u, 2) * pow(MUs, 2) * MQs * SS -
                2 * pow(u, 2) * pow(MUs, 3) * SS - 4 * pow(u, 2) * MXs * MUs * MQs * SS * AXG +
                8 * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                4 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                8 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                2 * pow(u, 2) * pow(MXs, 2) * MQs * SS +
                8 * pow(u, 2) * pow(MXs, 2) * MQs * UU * AXG -
                2 * pow(u, 2) * pow(MXs, 2) * MUs * SS +
                8 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                4 * pow(u, 2) * pow(MX, 4) * MQs * SS -
                2 * pow(u, 2) * pow(MX, 4) * MQs * SS * AXG -
                4 * pow(u, 2) * pow(MX, 4) * MQs * UU * AXG +
                4 * pow(u, 2) * pow(MX, 4) * MUs * SS -
                2 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG + 2 * pow(u, 3) * MUs * MQs * SS +
                2 * pow(u, 3) * MUs * MQs * SS * AXG) +
           syFC2 * (-4 * pow(u, 3) * MUs * MQs * UU * AXG + 2 * pow(u, 3) * pow(MUs, 2) * SS +
                    2 * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                    4 * pow(u, 3) * pow(MUs, 2) * UU * AXG - 2 * pow(u, 3) * MXs * MQs * SS +
                    4 * pow(u, 3) * MXs * MQs * SS * AXG - 8 * pow(u, 3) * MXs * MQs * UU * AXG -
                    2 * pow(u, 3) * MXs * MUs * SS + 4 * pow(u, 3) * MXs * MUs * SS * AXG -
                    8 * pow(u, 3) * MXs * MUs * UU * AXG - 2 * pow(u, 4) * MQs * SS * AXG +
                    4 * pow(u, 4) * MQs * UU * AXG - 2 * pow(u, 4) * MUs * SS * AXG +
                    4 * pow(u, 4) * MUs * UU * AXG + 2 * s * u * MXs * MUs * MQs * SS -
                    2 * s * u * MXs * MUs * MQs * SS * AXG - 3 * s * u * MXs * MUs * MQs * UU +
                    3 * s * u * MXs * MUs * MQs * UU * AXG - s * u * MXs * pow(MUs, 2) * UU +
                    s * u * MXs * pow(MUs, 2) * UU * AXG - 3 * s * u * pow(MXs, 2) * MQs * UU +
                    3 * s * u * pow(MXs, 2) * MQs * UU * AXG - s * u * pow(MXs, 2) * MUs * UU +
                    s * u * pow(MXs, 2) * MUs * UU * AXG + 3 * s * u * pow(MX, 4) * MQs * UU -
                    3 * s * u * pow(MX, 4) * MQs * UU * AXG + s * u * pow(MX, 4) * MUs * UU -
                    s * u * pow(MX, 4) * MUs * UU * AXG + 2 * s * pow(u, 2) * MUs * MQs * SS * AXG +
                    3 * s * pow(u, 2) * MUs * MQs * UU - 3 * s * pow(u, 2) * MUs * MQs * UU * AXG +
                    2 * s * pow(u, 2) * pow(MUs, 2) * SS + s * pow(u, 2) * pow(MUs, 2) * UU -
                    s * pow(u, 2) * pow(MUs, 2) * UU * AXG - 2 * s * pow(u, 2) * MXs * MQs * SS +
                    2 * s * pow(u, 2) * MXs * MQs * SS * AXG - s * pow(u, 2) * MXs * MQs * UU) +
           syFC2 * (-3 * s * pow(u, 2) * MXs * MQs * UU * AXG - 3 * s * pow(u, 2) * MXs * MUs * UU -
                    s * pow(u, 2) * MXs * MUs * UU * AXG - 2 * s * pow(u, 3) * MQs * SS * AXG +
                    s * pow(u, 3) * MQs * UU + 3 * s * pow(u, 3) * MQs * UU * AXG -
                    2 * s * pow(u, 3) * MUs * SS + 3 * s * pow(u, 3) * MUs * UU +
                    s * pow(u, 3) * MUs * UU * AXG) +
           syFC3 *
               (2 * MXs * pow(MUs, 3) * MQs * SS + 2 * pow(MXs, 2) * pow(MUs, 2) * MQs * SS -
                4 * pow(MX, 4) * pow(MUs, 2) * MQs * SS +
                2 * pow(MX, 4) * pow(MUs, 2) * MQs * SS * AXG -
                4 * pow(MX, 4) * pow(MUs, 2) * MQs * UU * AXG -
                4 * pow(MX, 4) * MXs * MUs * MQs * UU * AXG +
                4 * pow(MX, 6) * MUs * MQs * UU * AXG - 2 * u * pow(MUs, 3) * MQs * SS -
                4 * u * MXs * pow(MUs, 2) * MQs * SS * AXG +
                8 * u * MXs * pow(MUs, 2) * MQs * UU * AXG + 2 * u * MXs * pow(MUs, 3) * SS -
                2 * u * pow(MXs, 2) * MUs * MQs * SS + 8 * u * pow(MXs, 2) * MUs * MQs * UU * AXG +
                2 * u * pow(MXs, 2) * pow(MUs, 2) * SS + 4 * u * pow(MX, 4) * MUs * MQs * SS -
                2 * u * pow(MX, 4) * MUs * MQs * SS * AXG -
                4 * u * pow(MX, 4) * MUs * MQs * UU * AXG - 4 * u * pow(MX, 4) * pow(MUs, 2) * SS +
                2 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                4 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MX, 4) * MXs * MUs * UU * AXG + 4 * u * pow(MX, 6) * MUs * UU * AXG +
                2 * pow(u, 2) * pow(MUs, 2) * MQs * SS +
                2 * pow(u, 2) * pow(MUs, 2) * MQs * SS * AXG -
                4 * pow(u, 2) * pow(MUs, 2) * MQs * UU * AXG - 2 * pow(u, 2) * pow(MUs, 3) * SS -
                2 * pow(u, 2) * MXs * MUs * MQs * SS + 4 * pow(u, 2) * MXs * MUs * MQs * SS * AXG -
                8 * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                4 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                8 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                2 * pow(u, 2) * pow(MXs, 2) * MUs * SS) +
           syFC3 * (8 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                    4 * pow(u, 2) * pow(MX, 4) * MUs * SS -
                    2 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                    4 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                    2 * pow(u, 3) * MUs * MQs * SS * AXG + 4 * pow(u, 3) * MUs * MQs * UU * AXG +
                    2 * pow(u, 3) * pow(MUs, 2) * SS + 2 * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                    4 * pow(u, 3) * pow(MUs, 2) * UU * AXG - 2 * pow(u, 3) * MXs * MUs * SS +
                    4 * pow(u, 3) * MXs * MUs * SS * AXG - 8 * pow(u, 3) * MXs * MUs * UU * AXG -
                    2 * pow(u, 4) * MUs * SS * AXG + 4 * pow(u, 4) * MUs * UU * AXG -
                    s * MXs * pow(MUs, 2) * MQs * UU + s * MXs * pow(MUs, 2) * MQs * UU * AXG -
                    s * pow(MXs, 2) * MUs * MQs * UU + s * pow(MXs, 2) * MUs * MQs * UU * AXG +
                    s * pow(MX, 4) * MUs * MQs * UU - s * pow(MX, 4) * MUs * MQs * UU * AXG +
                    2 * s * u * pow(MUs, 2) * MQs * SS + s * u * pow(MUs, 2) * MQs * UU -
                    s * u * pow(MUs, 2) * MQs * UU * AXG - 3 * s * u * MXs * MUs * MQs * UU -
                    s * u * MXs * MUs * MQs * UU * AXG + 2 * s * u * MXs * pow(MUs, 2) * SS -
                    2 * s * u * MXs * pow(MUs, 2) * SS * AXG - 3 * s * u * MXs * pow(MUs, 2) * UU +
                    3 * s * u * MXs * pow(MUs, 2) * UU * AXG - 3 * s * u * pow(MXs, 2) * MUs * UU +
                    3 * s * u * pow(MXs, 2) * MUs * UU * AXG + 3 * s * u * pow(MX, 4) * MUs * UU -
                    3 * s * u * pow(MX, 4) * MUs * UU * AXG - 2 * s * pow(u, 2) * MUs * MQs * SS +
                    3 * s * pow(u, 2) * MUs * MQs * UU + s * pow(u, 2) * MUs * MQs * UU * AXG) +
           syFC3 *
               (2 * s * pow(u, 2) * pow(MUs, 2) * SS * AXG + 3 * s * pow(u, 2) * pow(MUs, 2) * UU -
                3 * s * pow(u, 2) * pow(MUs, 2) * UU * AXG - 2 * s * pow(u, 2) * MXs * MUs * SS +
                2 * s * pow(u, 2) * MXs * MUs * SS * AXG - s * pow(u, 2) * MXs * MUs * UU -
                3 * s * pow(u, 2) * MXs * MUs * UU * AXG - 2 * s * pow(u, 3) * MUs * SS * AXG +
                s * pow(u, 3) * MUs * UU + 3 * s * pow(u, 3) * MUs * UU * AXG) +
           syFC4 *
               (4 * s * u * MXs * pow(MUs, 2) * MQs * SS -
                4 * s * u * MXs * pow(MUs, 2) * MQs * SS * AXG -
                4 * s * u * MXs * pow(MUs, 2) * MQs * UU +
                4 * s * u * MXs * pow(MUs, 2) * MQs * UU * AXG +
                4 * s * u * MXs * pow(MUs, 3) * SS - 4 * s * u * MXs * pow(MUs, 3) * SS * AXG -
                4 * s * u * MXs * pow(MUs, 3) * UU + 4 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * SS +
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * SS * AXG +
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU -
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU * AXG -
                4 * s * pow(u, 2) * pow(MUs, 3) * SS + 4 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG +
                4 * s * pow(u, 2) * pow(MUs, 3) * UU - 4 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG -
                4 * s * pow(u, 2) * MXs * MUs * MQs * SS +
                4 * s * pow(u, 2) * MXs * MUs * MQs * SS * AXG +
                4 * s * pow(u, 2) * MXs * MUs * MQs * UU -
                4 * s * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS +
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU -
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                4 * s * pow(u, 3) * MUs * MQs * SS - 4 * s * pow(u, 3) * MUs * MQs * SS * AXG -
                4 * s * pow(u, 3) * MUs * MQs * UU + 4 * s * pow(u, 3) * MUs * MQs * UU * AXG +
                4 * s * pow(u, 3) * pow(MUs, 2) * SS - 4 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                4 * s * pow(u, 3) * pow(MUs, 2) * UU + 4 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG) +
           syFC5 *
               (8 * s * u * MXs * pow(MUs, 2) * MQs * SS -
                8 * s * u * MXs * pow(MUs, 2) * MQs * SS * AXG -
                8 * s * u * MXs * pow(MUs, 2) * MQs * UU +
                8 * s * u * MXs * pow(MUs, 2) * MQs * UU * AXG +
                12 * s * u * MXs * pow(MUs, 3) * SS - 12 * s * u * MXs * pow(MUs, 3) * SS * AXG -
                12 * s * u * MXs * pow(MUs, 3) * UU + 12 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                8 * s * pow(u, 2) * pow(MUs, 2) * MQs * SS +
                8 * s * pow(u, 2) * pow(MUs, 2) * MQs * SS * AXG +
                8 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU -
                8 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU * AXG -
                12 * s * pow(u, 2) * pow(MUs, 3) * SS +
                12 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG +
                12 * s * pow(u, 2) * pow(MUs, 3) * UU -
                12 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG -
                8 * s * pow(u, 2) * MXs * MUs * MQs * SS +
                8 * s * pow(u, 2) * MXs * MUs * MQs * SS * AXG +
                8 * s * pow(u, 2) * MXs * MUs * MQs * UU -
                8 * s * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(u, 3) * MUs * MQs * SS - 8 * s * pow(u, 3) * MUs * MQs * SS * AXG -
                8 * s * pow(u, 3) * MUs * MQs * UU + 8 * s * pow(u, 3) * MUs * MQs * UU * AXG +
                16 * s * pow(u, 3) * pow(MUs, 2) * SS -
                16 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                16 * s * pow(u, 3) * pow(MUs, 2) * UU) +
           syFC5 * (16 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG +
                    4 * s * pow(u, 3) * MXs * MUs * SS - 4 * s * pow(u, 3) * MXs * MUs * SS * AXG -
                    4 * s * pow(u, 3) * MXs * MUs * UU + 4 * s * pow(u, 3) * MXs * MUs * UU * AXG -
                    4 * s * pow(u, 4) * MUs * SS + 4 * s * pow(u, 4) * MUs * SS * AXG +
                    4 * s * pow(u, 4) * MUs * UU - 4 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC6 *
               (-8 * u * MXs * pow(MUs, 3) * MQs * SS - 8 * u * MXs * pow(MUs, 4) * SS -
                8 * u * pow(MXs, 2) * pow(MUs, 2) * MQs * SS -
                16 * u * pow(MXs, 2) * pow(MUs, 3) * SS - 8 * u * pow(MXs, 3) * pow(MUs, 2) * SS +
                16 * u * pow(MX, 4) * pow(MUs, 2) * MQs * SS -
                8 * u * pow(MX, 4) * pow(MUs, 2) * MQs * SS * AXG +
                16 * u * pow(MX, 4) * pow(MUs, 2) * MQs * UU * AXG +
                24 * u * pow(MX, 4) * pow(MUs, 3) * SS -
                8 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG +
                16 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG +
                16 * u * pow(MX, 4) * MXs * MUs * MQs * UU * AXG +
                24 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS -
                8 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG +
                32 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG +
                16 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG -
                16 * u * pow(MX, 6) * MUs * MQs * UU * AXG -
                16 * u * pow(MX, 6) * pow(MUs, 2) * SS +
                8 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG -
                32 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG -
                32 * u * pow(MX, 6) * MXs * MUs * UU * AXG + 16 * u * pow(MX, 8) * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MUs, 3) * MQs * SS + 8 * pow(u, 2) * pow(MUs, 4) * SS +
                16 * pow(u, 2) * MXs * pow(MUs, 2) * MQs * SS * AXG -
                32 * pow(u, 2) * MXs * pow(MUs, 2) * MQs * UU * AXG +
                16 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG -
                32 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * MUs * MQs * SS) +
           syFC6 *
               (-32 * pow(u, 2) * pow(MXs, 2) * MUs * MQs * UU * AXG +
                16 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                64 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 3) * MUs * SS -
                32 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                16 * pow(u, 2) * pow(MX, 4) * MUs * MQs * SS +
                8 * pow(u, 2) * pow(MX, 4) * MUs * MQs * SS * AXG +
                16 * pow(u, 2) * pow(MX, 4) * MUs * MQs * UU * AXG -
                8 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                48 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                24 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS +
                8 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG +
                48 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG +
                16 * pow(u, 2) * pow(MX, 6) * MUs * SS -
                8 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG -
                16 * pow(u, 2) * pow(MX, 6) * MUs * UU * AXG -
                8 * pow(u, 3) * pow(MUs, 2) * MQs * SS -
                8 * pow(u, 3) * pow(MUs, 2) * MQs * SS * AXG +
                16 * pow(u, 3) * pow(MUs, 2) * MQs * UU * AXG - 8 * pow(u, 3) * pow(MUs, 3) * SS -
                8 * pow(u, 3) * pow(MUs, 3) * SS * AXG + 16 * pow(u, 3) * pow(MUs, 3) * UU * AXG +
                8 * pow(u, 3) * MXs * MUs * MQs * SS - 16 * pow(u, 3) * MXs * MUs * MQs * SS * AXG +
                32 * pow(u, 3) * MXs * MUs * MQs * UU * AXG +
                8 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                16 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG) +
           syFC6 *
               (32 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG +
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS -
                16 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                16 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG -
                16 * pow(u, 3) * pow(MX, 4) * MUs * SS +
                16 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG -
                16 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG +
                8 * pow(u, 4) * MUs * MQs * SS * AXG - 16 * pow(u, 4) * MUs * MQs * UU * AXG +
                8 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 16 * pow(u, 4) * pow(MUs, 2) * UU * AXG +
                4 * s * u * MXs * pow(MUs, 2) * MQs * UU -
                4 * s * u * MXs * pow(MUs, 2) * MQs * UU * AXG -
                2 * s * u * MXs * pow(MUs, 3) * SS + 2 * s * u * MXs * pow(MUs, 3) * SS * AXG +
                6 * s * u * MXs * pow(MUs, 3) * UU - 6 * s * u * MXs * pow(MUs, 3) * UU * AXG +
                4 * s * u * pow(MXs, 2) * MUs * MQs * UU -
                4 * s * u * pow(MXs, 2) * MUs * MQs * UU * AXG -
                2 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                2 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                6 * s * u * pow(MXs, 3) * MUs * UU - 6 * s * u * pow(MXs, 3) * MUs * UU * AXG -
                4 * s * u * pow(MX, 4) * MUs * MQs * UU +
                4 * s * u * pow(MX, 4) * MUs * MQs * UU * AXG +
                2 * s * u * pow(MX, 4) * pow(MUs, 2) * SS -
                2 * s * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG) +
           syFC6 *
               (-12 * s * u * pow(MX, 4) * MXs * MUs * UU +
                12 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 6 * s * u * pow(MX, 6) * MUs * UU -
                6 * s * u * pow(MX, 6) * MUs * UU * AXG -
                8 * s * pow(u, 2) * pow(MUs, 2) * MQs * SS -
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU +
                4 * s * pow(u, 2) * pow(MUs, 2) * MQs * UU * AXG -
                6 * s * pow(u, 2) * pow(MUs, 3) * SS - 2 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG -
                6 * s * pow(u, 2) * pow(MUs, 3) * UU + 6 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG +
                12 * s * pow(u, 2) * MXs * MUs * MQs * UU +
                4 * s * pow(u, 2) * MXs * MUs * MQs * UU * AXG -
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS +
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                2 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS -
                2 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS * AXG +
                22 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU -
                6 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                2 * s * pow(u, 2) * pow(MX, 4) * MUs * SS +
                2 * s * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                22 * s * pow(u, 2) * pow(MX, 4) * MUs * UU +
                6 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(u, 3) * MUs * MQs * SS - 12 * s * pow(u, 3) * MUs * MQs * UU -
                4 * s * pow(u, 3) * MUs * MQs * UU * AXG + 12 * s * pow(u, 3) * pow(MUs, 2) * SS -
                4 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG - 16 * s * pow(u, 3) * pow(MUs, 2) * UU +
                6 * s * pow(u, 3) * MXs * MUs * SS) +
           syFC6 * (-6 * s * pow(u, 3) * MXs * MUs * SS * AXG - 6 * s * pow(u, 3) * MXs * MUs * UU +
                    6 * s * pow(u, 3) * MXs * MUs * UU * AXG - 6 * s * pow(u, 4) * MUs * SS +
                    6 * s * pow(u, 4) * MUs * SS * AXG + 6 * s * pow(u, 4) * MUs * UU -
                    6 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC7 *
               (16 * u * MXs * pow(MUs, 3) * SS + 16 * u * pow(MXs, 2) * pow(MUs, 2) * SS -
                32 * u * pow(MX, 4) * pow(MUs, 2) * SS +
                16 * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                32 * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                32 * u * pow(MX, 4) * MXs * MUs * UU * AXG + 32 * u * pow(MX, 6) * MUs * UU * AXG -
                16 * pow(u, 2) * pow(MUs, 3) * SS - 32 * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                64 * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                16 * pow(u, 2) * pow(MXs, 2) * MUs * SS +
                64 * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                32 * pow(u, 2) * pow(MX, 4) * MUs * SS -
                16 * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                32 * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG + 16 * pow(u, 3) * pow(MUs, 2) * SS +
                16 * pow(u, 3) * pow(MUs, 2) * SS * AXG - 32 * pow(u, 3) * pow(MUs, 2) * UU * AXG -
                16 * pow(u, 3) * MXs * MUs * SS + 32 * pow(u, 3) * MXs * MUs * SS * AXG -
                64 * pow(u, 3) * MXs * MUs * UU * AXG - 16 * pow(u, 4) * MUs * SS * AXG +
                32 * pow(u, 4) * MUs * UU * AXG + 8 * s * u * MXs * pow(MUs, 2) * SS -
                8 * s * u * MXs * pow(MUs, 2) * SS * AXG - 16 * s * u * MXs * pow(MUs, 2) * UU +
                16 * s * u * MXs * pow(MUs, 2) * UU * AXG - 16 * s * u * pow(MXs, 2) * MUs * UU +
                16 * s * u * pow(MXs, 2) * MUs * UU * AXG + 16 * s * u * pow(MX, 4) * MUs * UU -
                16 * s * u * pow(MX, 4) * MUs * UU * AXG + 8 * s * pow(u, 2) * pow(MUs, 2) * SS +
                8 * s * pow(u, 2) * pow(MUs, 2) * SS * AXG) +
           syFC7 *
               (16 * s * pow(u, 2) * pow(MUs, 2) * UU -
                16 * s * pow(u, 2) * pow(MUs, 2) * UU * AXG - 8 * s * pow(u, 2) * MXs * MUs * SS +
                8 * s * pow(u, 2) * MXs * MUs * SS * AXG - 16 * s * pow(u, 2) * MXs * MUs * UU -
                16 * s * pow(u, 2) * MXs * MUs * UU * AXG - 8 * s * pow(u, 3) * MUs * SS -
                8 * s * pow(u, 3) * MUs * SS * AXG + 16 * s * pow(u, 3) * MUs * UU +
                16 * s * pow(u, 3) * MUs * UU * AXG) +
           syFC8 *
               (4 * s * u * MXs * pow(MUs, 3) * SS - 4 * s * u * MXs * pow(MUs, 3) * SS * AXG -
                8 * s * u * MXs * pow(MUs, 3) * UU + 8 * s * u * MXs * pow(MUs, 3) * UU * AXG +
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS -
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU +
                12 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG -
                4 * s * u * pow(MXs, 3) * MUs * UU + 4 * s * u * pow(MXs, 3) * MUs * UU * AXG -
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS +
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU -
                12 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                8 * s * u * pow(MX, 4) * MXs * MUs * UU -
                8 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG - 4 * s * u * pow(MX, 6) * MUs * UU +
                4 * s * u * pow(MX, 6) * MUs * UU * AXG + 4 * s * u * pow(MU, 4) * MXs * MUs * SS -
                4 * s * u * pow(MU, 4) * MXs * MUs * SS * AXG -
                4 * s * pow(u, 2) * pow(MUs, 3) * SS + 4 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG +
                8 * s * pow(u, 2) * pow(MUs, 3) * UU - 8 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU -
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS +
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS * AXG +
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU) +
           syFC8 * (-8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                    4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS -
                    4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG -
                    8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU +
                    8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                    4 * s * pow(u, 2) * pow(MU, 4) * MUs * SS +
                    4 * s * pow(u, 2) * pow(MU, 4) * MUs * SS * AXG +
                    16 * s * pow(u, 3) * pow(MUs, 2) * SS -
                    16 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG -
                    16 * s * pow(u, 3) * pow(MUs, 2) * UU +
                    16 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG +
                    8 * s * pow(u, 3) * MXs * MUs * SS - 8 * s * pow(u, 3) * MXs * MUs * SS * AXG -
                    8 * s * pow(u, 3) * MXs * MUs * UU + 8 * s * pow(u, 3) * MXs * MUs * UU * AXG -
                    8 * s * pow(u, 4) * MUs * SS + 8 * s * pow(u, 4) * MUs * SS * AXG +
                    8 * s * pow(u, 4) * MUs * UU - 8 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC9 *
               (-4 * u * MXs * pow(MUs, 4) * SS - 8 * u * pow(MXs, 2) * pow(MUs, 3) * SS -
                4 * u * pow(MXs, 3) * pow(MUs, 2) * SS + 12 * u * pow(MX, 4) * pow(MUs, 3) * SS -
                4 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG +
                16 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG +
                12 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS -
                4 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG +
                24 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG +
                8 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG -
                8 * u * pow(MX, 6) * pow(MUs, 2) * SS +
                4 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG -
                24 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG -
                16 * u * pow(MX, 6) * MXs * MUs * UU * AXG + 8 * u * pow(MX, 8) * MUs * UU * AXG -
                4 * u * pow(MU, 4) * MXs * pow(MUs, 2) * SS -
                4 * u * pow(MU, 4) * pow(MXs, 2) * MUs * SS +
                8 * u * pow(MU, 4) * pow(MX, 4) * MUs * SS -
                4 * u * pow(MU, 4) * pow(MX, 4) * MUs * SS * AXG +
                4 * pow(u, 2) * pow(MUs, 4) * SS + 12 * pow(u, 2) * MXs * pow(MUs, 3) * SS +
                8 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG -
                32 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                12 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS +
                8 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG -
                48 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                4 * pow(u, 2) * pow(MXs, 3) * MUs * SS -
                16 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                28 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS) +
           syFC9 *
               (8 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                16 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                12 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS +
                4 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG +
                8 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MX, 6) * MUs * SS -
                4 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG +
                8 * pow(u, 2) * pow(MX, 6) * MUs * UU * AXG +
                4 * pow(u, 2) * pow(MU, 4) * pow(MUs, 2) * SS -
                4 * pow(u, 2) * pow(MU, 4) * MXs * MUs * SS +
                8 * pow(u, 2) * pow(MU, 4) * MXs * MUs * SS * AXG -
                16 * pow(u, 3) * pow(MUs, 3) * SS - 4 * pow(u, 3) * pow(MUs, 3) * SS * AXG +
                16 * pow(u, 3) * pow(MUs, 3) * UU * AXG + 8 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                32 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG +
                64 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                40 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG +
                8 * pow(u, 3) * pow(MX, 4) * MUs * SS -
                24 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(u, 3) * pow(MU, 4) * MUs * SS * AXG + 8 * pow(u, 4) * pow(MUs, 2) * SS +
                16 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 32 * pow(u, 4) * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 4) * MXs * MUs * SS + 16 * pow(u, 4) * MXs * MUs * SS * AXG -
                32 * pow(u, 4) * MXs * MUs * UU * AXG - 8 * pow(u, 5) * MUs * SS * AXG) +
           syFC9 *
               (16 * pow(u, 5) * MUs * UU * AXG - 4 * s * u * MXs * pow(MUs, 3) * SS +
                4 * s * u * MXs * pow(MUs, 3) * SS * AXG + 8 * s * u * MXs * pow(MUs, 3) * UU -
                8 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS +
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                14 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU -
                14 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG +
                6 * s * u * pow(MXs, 3) * MUs * UU - 6 * s * u * pow(MXs, 3) * MUs * UU * AXG +
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS -
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                14 * s * u * pow(MX, 4) * pow(MUs, 2) * UU +
                14 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                12 * s * u * pow(MX, 4) * MXs * MUs * UU +
                12 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG + 6 * s * u * pow(MX, 6) * MUs * UU -
                6 * s * u * pow(MX, 6) * MUs * UU * AXG -
                4 * s * pow(u, 2) * pow(MUs, 3) * SS * AXG - 8 * s * pow(u, 2) * pow(MUs, 3) * UU +
                8 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG -
                8 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS +
                8 * s * pow(u, 2) * MXs * pow(MUs, 2) * SS * AXG +
                16 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS -
                4 * s * pow(u, 2) * pow(MXs, 2) * MUs * SS * AXG +
                16 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU -
                8 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS +
                4 * s * pow(u, 2) * pow(MX, 4) * MUs * SS * AXG) +
           syFC9 *
               (-16 * s * pow(u, 2) * pow(MX, 4) * MUs * UU +
                8 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG -
                4 * s * pow(u, 2) * pow(MU, 4) * MUs * SS + 24 * s * pow(u, 3) * pow(MUs, 2) * SS -
                8 * s * pow(u, 3) * pow(MUs, 2) * SS * AXG - 16 * s * pow(u, 3) * pow(MUs, 2) * UU +
                12 * s * pow(u, 3) * MXs * MUs * SS - 12 * s * pow(u, 3) * MXs * MUs * SS * AXG -
                24 * s * pow(u, 3) * MXs * MUs * UU + 8 * s * pow(u, 3) * MXs * MUs * UU * AXG -
                20 * s * pow(u, 4) * MUs * SS + 12 * s * pow(u, 4) * MUs * SS * AXG +
                24 * s * pow(u, 4) * MUs * UU - 8 * s * pow(u, 4) * MUs * UU * AXG) +
           syFC10 *
               (4 * u * MXs * pow(MUs, 4) * SS + 8 * u * pow(MXs, 2) * pow(MUs, 3) * SS +
                4 * u * pow(MXs, 3) * pow(MUs, 2) * SS - 12 * u * pow(MX, 4) * pow(MUs, 3) * SS +
                4 * u * pow(MX, 4) * pow(MUs, 3) * SS * AXG -
                8 * u * pow(MX, 4) * pow(MUs, 3) * UU * AXG -
                12 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS +
                4 * u * pow(MX, 4) * MXs * pow(MUs, 2) * SS * AXG -
                16 * u * pow(MX, 4) * MXs * pow(MUs, 2) * UU * AXG -
                8 * u * pow(MX, 4) * pow(MXs, 2) * MUs * UU * AXG +
                8 * u * pow(MX, 6) * pow(MUs, 2) * SS -
                4 * u * pow(MX, 6) * pow(MUs, 2) * SS * AXG +
                16 * u * pow(MX, 6) * pow(MUs, 2) * UU * AXG +
                16 * u * pow(MX, 6) * MXs * MUs * UU * AXG - 8 * u * pow(MX, 8) * MUs * UU * AXG -
                4 * pow(u, 2) * pow(MUs, 4) * SS + 12 * pow(u, 2) * MXs * pow(MUs, 3) * SS -
                8 * pow(u, 2) * MXs * pow(MUs, 3) * SS * AXG +
                16 * pow(u, 2) * MXs * pow(MUs, 3) * UU * AXG +
                12 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * SS * AXG +
                32 * pow(u, 2) * pow(MXs, 2) * pow(MUs, 2) * UU * AXG -
                4 * pow(u, 2) * pow(MXs, 3) * MUs * SS +
                16 * pow(u, 2) * pow(MXs, 3) * MUs * UU * AXG -
                20 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS +
                16 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                48 * pow(u, 2) * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                12 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS) +
           syFC10 *
               (-4 * pow(u, 2) * pow(MX, 4) * MXs * MUs * SS * AXG -
                48 * pow(u, 2) * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(u, 2) * pow(MX, 6) * MUs * SS +
                4 * pow(u, 2) * pow(MX, 6) * MUs * SS * AXG +
                32 * pow(u, 2) * pow(MX, 6) * MUs * UU * AXG - 8 * pow(u, 3) * pow(MUs, 3) * SS +
                4 * pow(u, 3) * pow(MUs, 3) * SS * AXG - 8 * pow(u, 3) * pow(MUs, 3) * UU * AXG -
                4 * pow(u, 3) * MXs * pow(MUs, 2) * SS -
                16 * pow(u, 3) * MXs * pow(MUs, 2) * SS * AXG +
                32 * pow(u, 3) * MXs * pow(MUs, 2) * UU * AXG -
                20 * pow(u, 3) * pow(MXs, 2) * MUs * SS +
                8 * pow(u, 3) * pow(MXs, 2) * MUs * SS * AXG +
                40 * pow(u, 3) * pow(MXs, 2) * MUs * UU * AXG +
                32 * pow(u, 3) * pow(MX, 4) * MUs * SS -
                20 * pow(u, 3) * pow(MX, 4) * MUs * SS * AXG -
                16 * pow(u, 3) * pow(MX, 4) * MUs * UU * AXG + 12 * pow(u, 4) * pow(MUs, 2) * SS +
                8 * pow(u, 4) * pow(MUs, 2) * SS * AXG - 16 * pow(u, 4) * pow(MUs, 2) * UU * AXG -
                12 * pow(u, 4) * MXs * MUs * SS + 24 * pow(u, 4) * MXs * MUs * SS * AXG -
                48 * pow(u, 4) * MXs * MUs * UU * AXG - 12 * pow(u, 5) * MUs * SS * AXG +
                24 * pow(u, 5) * MUs * UU * AXG - 2 * s * u * MXs * pow(MUs, 3) * UU +
                2 * s * u * MXs * pow(MUs, 3) * UU * AXG -
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU +
                4 * s * u * pow(MXs, 2) * pow(MUs, 2) * UU * AXG -
                2 * s * u * pow(MXs, 3) * MUs * UU + 2 * s * u * pow(MXs, 3) * MUs * UU * AXG) +
           syFC10 *
               (4 * s * u * pow(MX, 4) * pow(MUs, 2) * UU -
                4 * s * u * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * s * u * pow(MX, 4) * MXs * MUs * UU -
                4 * s * u * pow(MX, 4) * MXs * MUs * UU * AXG - 2 * s * u * pow(MX, 6) * MUs * UU +
                2 * s * u * pow(MX, 6) * MUs * UU * AXG + 4 * s * pow(u, 2) * pow(MUs, 3) * SS +
                2 * s * pow(u, 2) * pow(MUs, 3) * UU - 2 * s * pow(u, 2) * pow(MUs, 3) * UU * AXG -
                12 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU +
                4 * s * pow(u, 2) * MXs * pow(MUs, 2) * UU * AXG -
                14 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU +
                6 * s * pow(u, 2) * pow(MXs, 2) * MUs * UU * AXG +
                14 * s * pow(u, 2) * pow(MX, 4) * MUs * UU -
                6 * s * pow(u, 2) * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(u, 3) * pow(MUs, 2) * SS + 12 * s * pow(u, 3) * pow(MUs, 2) * UU -
                4 * s * pow(u, 3) * pow(MUs, 2) * UU * AXG - 18 * s * pow(u, 3) * MXs * MUs * UU -
                6 * s * pow(u, 3) * MXs * MUs * UU * AXG - 12 * s * pow(u, 4) * MUs * SS +
                18 * s * pow(u, 4) * MUs * UU + 6 * s * pow(u, 4) * MUs * UU * AXG)));
  ;
  return ret.real();
}

ComplexType ME_us_QQg_Qg1(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  int q_I = q;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);
  ;

  double SS = sc, UU = uc, AXG = axial;
  auto Denom = [](auto a) { return 1. / a; };

  auto MQi = MQ;
  auto MQis = MQs;

#define syFC1 A0(MQs)
#define syFC2 B0(MUs, 0, MQs)
  ComplexType ret = 0;
  _EPS0(ret,
        -(-1) * (-2 + pow(Nc, 2)) * (-1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
            (+Denom(-3072 * pow(Pi, 2) * s * pow(MU, 2) * MXs * MQis * Nc +
                    3072 * pow(Pi, 2) * s * pow(MU, 2) * MXs * MUs * Nc +
                    3072 * pow(Pi, 2) * s * pow(MU, 2) * pow(MXs, 2) * Nc +
                    3072 * pow(Pi, 2) * s * t * pow(MU, 2) * MQis * Nc -
                    3072 * pow(Pi, 2) * s * t * pow(MU, 2) * MUs * Nc -
                    6144 * pow(Pi, 2) * s * t * pow(MU, 2) * MXs * Nc +
                    3072 * pow(Pi, 2) * s * pow(t, 2) * pow(MU, 2) * Nc +
                    3072 * pow(Pi, 2) * pow(s, 2) * pow(MU, 2) * MQis * Nc -
                    3072 * pow(Pi, 2) * pow(s, 2) * pow(MU, 2) * MUs * Nc -
                    6144 * pow(Pi, 2) * pow(s, 2) * pow(MU, 2) * MXs * Nc +
                    6144 * pow(Pi, 2) * pow(s, 2) * t * pow(MU, 2) * Nc +
                    3072 * pow(Pi, 2) * pow(s, 3) * pow(MU, 2) * Nc)) *
            (+R * Lp + L * Rp) * (+syFC1 * (1.) + syFC2 * (-MQs - 3 * MUs)) * (+gs) * (+gs) *
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
             2 * pow(u, 3) * SS * AXG + 4 * pow(u, 3) * UU * AXG + 2 * s * MXs * MUs * SS -
             2 * s * MXs * MUs * SS * AXG - 3 * s * MXs * MUs * UU + 3 * s * MXs * MUs * UU * AXG -
             3 * s * pow(MXs, 2) * UU + 3 * s * pow(MXs, 2) * UU * AXG + 3 * s * pow(MX, 4) * UU -
             3 * s * pow(MX, 4) * UU * AXG + 2 * s * u * MUs * SS * AXG + 3 * s * u * MUs * UU -
             3 * s * u * MUs * UU * AXG - 2 * s * u * MXs * SS + 2 * s * u * MXs * SS * AXG -
             s * u * MXs * UU - 3 * s * u * MXs * UU * AXG - 2 * s * pow(u, 2) * SS * AXG +
             s * pow(u, 2) * UU + 3 * s * pow(u, 2) * UU * AXG));

  return ret.real();
}

ComplexType ME_us_QQg_Qg2(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
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

  auto MQi = MQ;
  auto MQis = MQs;

#define syFC1 A0(MQs)
#define syFC2 B0(MUs + MXs - s - t, 0, MQs)
  ComplexType ret = 0;
  _EPS0(ret,
        -(-2 + pow(Nc, 2)) * (-1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
            (+Denom(3072 * pow(Pi, 2) * s * MXs * MUs * MQis * Nc -
                    3072 * pow(Pi, 2) * s * MXs * pow(MUs, 2) * Nc +
                    3072 * pow(Pi, 2) * s * pow(MXs, 2) * MQis * Nc -
                    6144 * pow(Pi, 2) * s * pow(MXs, 2) * MUs * Nc -
                    3072 * pow(Pi, 2) * s * pow(MXs, 3) * Nc -
                    3072 * pow(Pi, 2) * s * t * MUs * MQis * Nc +
                    3072 * pow(Pi, 2) * s * t * pow(MUs, 2) * Nc -
                    6144 * pow(Pi, 2) * s * t * MXs * MQis * Nc +
                    12288 * pow(Pi, 2) * s * t * MXs * MUs * Nc +
                    9216 * pow(Pi, 2) * s * t * pow(MXs, 2) * Nc +
                    3072 * pow(Pi, 2) * s * pow(t, 2) * MQis * Nc -
                    6144 * pow(Pi, 2) * s * pow(t, 2) * MUs * Nc -
                    9216 * pow(Pi, 2) * s * pow(t, 2) * MXs * Nc +
                    3072 * pow(Pi, 2) * s * pow(t, 3) * Nc -
                    3072 * pow(Pi, 2) * pow(s, 2) * MUs * MQis * Nc +
                    3072 * pow(Pi, 2) * pow(s, 2) * pow(MUs, 2) * Nc -
                    6144 * pow(Pi, 2) * pow(s, 2) * MXs * MQis * Nc +
                    12288 * pow(Pi, 2) * pow(s, 2) * MXs * MUs * Nc +
                    9216 * pow(Pi, 2) * pow(s, 2) * pow(MXs, 2) * Nc +
                    6144 * pow(Pi, 2) * pow(s, 2) * t * MQis * Nc -
                    12288 * pow(Pi, 2) * pow(s, 2) * t * MUs * Nc -
                    18432 * pow(Pi, 2) * pow(s, 2) * t * MXs * Nc +
                    9216 * pow(Pi, 2) * pow(s, 2) * pow(t, 2) * Nc +
                    3072 * pow(Pi, 2) * pow(s, 3) * MQis * Nc -
                    6144 * pow(Pi, 2) * pow(s, 3) * MUs * Nc -
                    9216 * pow(Pi, 2) * pow(s, 3) * MXs * Nc +
                    9216 * pow(Pi, 2) * pow(s, 3) * t * Nc + 3072 * pow(Pi, 2) * pow(s, 4) * Nc)) *
            (+R * Lp + L * Rp) * (+syFC1 * (1.) + syFC2 * (-MQs - 3 * u)) * (+gs) * (+gs) * (+gs) *
            (+gs) *
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
             2 * pow(u, 3) * SS * AXG + 4 * pow(u, 3) * UU * AXG - s * MXs * MUs * UU +
             s * MXs * MUs * UU * AXG - s * pow(MXs, 2) * UU + s * pow(MXs, 2) * UU * AXG +
             s * pow(MX, 4) * UU - s * pow(MX, 4) * UU * AXG + 2 * s * u * MUs * SS +
             s * u * MUs * UU - s * u * MUs * UU * AXG - 3 * s * u * MXs * UU -
             s * u * MXs * UU * AXG - 2 * s * pow(u, 2) * SS + 3 * s * pow(u, 2) * UU +
             s * pow(u, 2) * UU * AXG));

  return ret.real();
}
ComplexType ME_us_QQg_Qg(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params) {
  return ME_us_QQg_Qg1(pIEPS, sc, uc, axial, Q2, P1K1, params) +
         ME_us_QQg_Qg2(pIEPS, sc, uc, axial, Q2, P1K1, params);
}

ComplexType ME_us_QQg_qGG(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  auto Denom = [](auto a) { return 1. / a; };
  ComplexType ret = 0;
  int itq = q;
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int itq = 0; itq < 3; itq++) {
    int isq = (is_chargino(ch) != is_up_quark(itq)) * 6 + itsq * 3 + itq -
              (is_up_quark(itq)) * 3;
    int iq = (itq + is_chargino(ch) * 3) % 6;

    ComplexType L = (params->CHSQq[ch][sq][q].L);
    ComplexType R = (params->CHSQq[ch][sq][q].R);
    ComplexType Lp = conj(params->CHSQq[ch][isq][q].R);
    ComplexType Rp = conj(params->CHSQq[ch][isq][q].L);

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    ComplexType LG = (params->GLSQq[isq][iq].L);
    ComplexType RG = (params->GLSQq[isq][iq].R);
    ComplexType LGp = conj(params->GLSQq[sq][iq].R);
    ComplexType RGp = conj(params->GLSQq[sq][iq].L);

#define syFC1 A0(MGs)
#define syFC2 A0(Mqis)
#define syFC3 B0(MUs, MGs, Mqis)
#define syFC4 B0(u, MGs, Mqis)
#define syFC5 C0(0, u, MUs, MGs, MGs, Mqis)
#define syFC6 C1(0, MUs, u, MGs, MGs, Mqis)
#define syFC7 C1(0, u, MUs, MGs, MGs, Mqis)
#define syFC8 C2(0, MUs, u, MGs, MGs, Mqis)
#define syFC9 C2(0, u, MUs, MGs, MGs, Mqis)
#define syFC10 C00(0, MUs, u, MGs, MGs, Mqis)
#define syFC11 C11(0, MUs, u, MGs, MGs, Mqis)
#define syFC12 C12(0, MUs, u, MGs, MGs, Mqis)
#define syFC13 C22(0, MUs, u, MGs, MGs, Mqis)

    _EPS0_(
        ret, auto x0 = pow(MUs, 2); auto x1 = u * x0; auto x2 = LG * RGp; auto x3 = SS * x2;
        auto x4 = 2.0 * x3; auto x5 = x1 * x4; auto x6 = pow(MXs, 2); auto x7 = u * x6;
        auto x8 = x4 * x7; auto x9 = LGp * RG; auto x10 = SS * x9; auto x11 = 2.0 * x10;
        auto x12 = x1 * x11; auto x13 = x11 * x7; auto x14 = pow(MX, 4); auto x15 = 4.0 * x14;
        auto x16 = MUs * x3; auto x17 = x15 * x16; auto x18 = x10 * x15; auto x19 = MUs * x18;
        auto x20 = MUs * x6; auto x21 = x20 * x4; auto x22 = pow(u, 2); auto x23 = MUs * x22;
        auto x24 = x23 * x3; auto x25 = 2.0 * x24; auto x26 = MXs * x0; auto x27 = x26 * x4;
        auto x28 = x11 * x20; auto x29 = x10 * x23; auto x30 = 2.0 * x29; auto x31 = x11 * x26;
        auto x32 = UU * x2; auto x33 = 4.0 * x32; auto x34 = pow(MX, 6); auto x35 = AXG * x34;
        auto x36 = x33 * x35; auto x37 = UU * x9; auto x38 = 4.0 * AXG; auto x39 = x37 * x38;
        auto x40 = x34 * x39; auto x41 = x15 * x3; auto x42 = u * x14; auto x43 = AXG * x4;
        auto x44 = AXG * x11; auto x45 = x14 * x32; auto x46 = x38 * x45; auto x47 = MXs * x46;
        auto x48 = UU * x22; auto x49 = x48 * x9; auto x50 = MUs * x49; auto x51 = x38 * x50;
        auto x52 = x14 * x39; auto x53 = MXs * x52; auto x54 = MUs * x14; auto x55 = x43 * x54;
        auto x56 = x44 * x54; auto x57 = AXG * x30; auto x58 = 8.0 * AXG; auto x59 = x58 * x7;
        auto x60 = x32 * x59; auto x61 = x37 * x59; auto x62 = x3 * x38; auto x63 = MUs * u;
        auto x64 = MXs * x63; auto x65 = x62 * x64; auto x66 = x10 * x38; auto x67 = x64 * x66;
        auto x68 = x32 * x64; auto x69 = x58 * x68; auto x70 = x37 * x63; auto x71 = MXs * x70;
        auto x72 = x58 * x71; auto x73 = MUs * x45; auto x74 = x38 * x73; auto x75 = MUs * x37;
        auto x76 = x14 * x75; auto x77 = x38 * x76; auto x78 = x74 + x77;
        auto x79 = -u * x18 - u * x41 + u * x46 + x12 + x13 + x17 + x19 - x21 - x25 - x27 - x28 -
                   x30 - x31 - x36 + x39 * x42 - x40 + x42 * x43 + x42 * x44 + x47 + x5 + x51 +
                   x53 - x55 - x56 - x57 - x60 - x61 + x65 + x67 - x69 - x72 + x78 + x8;
        auto x80 = pow(u, 3); auto x81 = AXG * x80; auto x82 = x4 * x81; auto x83 = x11 * x81;
        auto x84 = MXs * x22; auto x85 = x4 * x84; auto x86 = x11 * x84; auto x87 = 3.0 * s;
        auto x88 = x32 * x6; auto x89 = x87 * x88; auto x90 = x37 * x87; auto x91 = x6 * x90;
        auto x92 = x45 * x87; auto x93 = x14 * x90; auto x94 = MXs * u; auto x95 = x37 * x94;
        auto x96 = s * x95; auto x97 = s * x11; auto x98 = x2 * x48; auto x99 = MUs * x98;
        auto x100 = x38 * x99; auto x101 = MXs * x58; auto x102 = x101 * x98;
        auto x103 = x101 * x49; auto x104 = MUs * x32; auto x105 = MXs * x104;
        auto x106 = x105 * x87; auto x107 = MXs * x75; auto x108 = x107 * x87;
        auto x109 = AXG * x25; auto x110 = MUs * MXs; auto x111 = x110 * x4; auto x112 = s * x111;
        auto x113 = x110 * x97; auto x114 = x62 * x84; auto x115 = x66 * x84; auto x116 = x32 * x63;
        auto x117 = x116 * x87; auto x118 = x63 * x90; auto x119 = s * x110; auto x120 = x119 * x43;
        auto x121 = x119 * x44; auto x122 = s * x63; auto x123 = AXG * x108; auto x124 = x33 * x80;
        auto x125 = AXG * x124; auto x126 = 4.0 * x80; auto x127 = x126 * x37;
        auto x128 = AXG * x127; auto x129 = -x125 - x128;
        auto x130 = -AXG * x106 + AXG * x117 + AXG * x118 - AXG * x89 - AXG * x91 + AXG * x92 +
                    AXG * x93 - s * x44 * x94 + x100 + x102 + x103 + x106 + x108 - x109 - x112 -
                    x113 - x114 - x115 - x117 - x118 + x120 + x121 - x122 * x43 - x122 * x44 -
                    x123 + x129 + x82 + x83 + x85 + x86 + x89 + x91 - x92 - x93 + x94 * x97 + x96;
        auto x131 = pow(MU, 4); auto x132 = 16.0 * x131; auto x133 = x14 * x3;
        auto x134 = x10 * x132; auto x135 = x131 * x6; auto x136 = 8.0 * x135;
        auto x137 = x136 * x3; auto x138 = x10 * x136; auto x139 = x6 * x63; auto x140 = 8.0 * x139;
        auto x141 = 8.0 * x0; auto x142 = x141 * x3; auto x143 = 8.0 * x131; auto x144 = x10 * x143;
        auto x145 = x10 * x141; auto x146 = 16.0 * x0; auto x147 = AXG * x45;
        auto x148 = AXG * x146; auto x149 = x14 * x37; auto x150 = x143 * x3;
        auto x151 = x14 * x150; auto x152 = AXG * x151; auto x153 = x14 * x144;
        auto x154 = AXG * x153; auto x155 = MXs * x150; auto x156 = MXs * x144;
        auto x157 = 16.0 * x35; auto x158 = 16.0 * x63; auto x159 = x10 * x14;
        auto x160 = x133 * x63; auto x161 = 8.0 * x14; auto x162 = x10 * x161 * x63;
        auto x163 = AXG * x73; auto x164 = 16.0 * MXs; auto x165 = x45 * x63;
        auto x166 = 16.0 * AXG; auto x167 = AXG * x76; auto x168 = x14 * x70;
        auto x169 = AXG * x168; auto x170 = x63 * x88; auto x171 = AXG * x170;
        auto x172 = x32 * x94; auto x173 = x0 * x172; auto x174 = AXG * x173; auto x175 = x6 * x70;
        auto x176 = AXG * x175; auto x177 = x0 * x95; auto x178 = AXG * x177; auto x179 = x22 * x3;
        auto x180 = MUs * x33; auto x181 = x180 * x6; auto x182 = s * x181; auto x183 = x0 * x33;
        auto x184 = MXs * x183; auto x185 = s * x184; auto x186 = x6 * x75; auto x187 = s * x186;
        auto x188 = x26 * x37; auto x189 = s * x188; auto x190 = 8.0 * x81; auto x191 = x16 * x190;
        auto x192 = MUs * x10; auto x193 = x190 * x192; auto x194 = 16.0 * x81;
        auto x195 = AXG * x94; auto x196 = x132 * x195; auto x197 = AXG * x99;
        auto x198 = 32.0 * MXs; auto x199 = AXG * x50; auto x200 = AXG * x24; auto x201 = AXG * x29;
        auto x202 = AXG * x22; auto x203 = x150 * x202; auto x204 = x144 * x202;
        auto x205 = -x203 - x204; auto x206 = 8.0 * x24; auto x207 = MXs * x206;
        auto x208 = 8.0 * x29; auto x209 = MXs * x208; auto x210 = -x150 * x94;
        auto x211 = -x144 * x94; auto x212 = x207 + x209 + x210 + x211; auto x213 = SS * pow(MU, 6);
        auto x214 = 4.0 * x213; auto x215 = x2 * x214; auto x216 = x214 * x9;
        auto x217 = x161 * x213; auto x218 = x14 * x215; auto x219 = x14 * x216;
        auto x220 = MXs * x213; auto x221 = x2 * x220; auto x222 = 4.0 * MUs; auto x223 = x220 * x9;
        auto x224 = x143 * x32; auto x225 = x131 * x37; auto x226 = 8.0 * x35; auto x227 = u * x138;
        auto x228 = x1 * x37; auto x229 = AXG * x228; auto x230 = x135 * x37;
        auto x231 = AXG * x230; auto x232 = 16.0 * u; auto x233 = x131 * x45;
        auto x234 = AXG * x233; auto x235 = 8.0 * MXs; auto x236 = x14 * x225;
        auto x237 = AXG * x236; auto x238 = 8.0 * x237; auto x239 = x155 * x63;
        auto x240 = x156 * x63; auto x241 = AXG * x71; auto x242 = MXs * x131;
        auto x243 = x242 * x32; auto x244 = AXG * x243; auto x245 = MXs * x225;
        auto x246 = AXG * x245; auto x247 = x126 * x3; auto x248 = x10 * x126;
        auto x249 = MUs * x247; auto x250 = MXs * x249; auto x251 = MUs * x248;
        auto x252 = MXs * x251; auto x253 = x141 * x81; auto x254 = x14 * x38;
        auto x255 = x166 * x242; auto x256 = x242 * x3; auto x257 = 16.0 * x202;
        auto x258 = x10 * x242; auto x259 = x144 * x23 - x144 * x81 - x144 * x84 + x150 * x23 -
                                            x150 * x84 + x256 * x257 + x257 * x258;
        auto x260 = 2.0 * UU; auto x261 = x2 * x260; auto x262 = s * x261; auto x263 = x131 * x14;
        auto x264 = x260 * x9; auto x265 = s * x264; auto x266 = pow(u, 4); auto x267 = MUs * x266;
        auto x268 = x135 * x262; auto x269 = x135 * x265; auto x270 = x266 * x58;
        auto x271 = -x150 * x81; auto x272 = 6.0 * s; auto x273 = AXG * x131;
        auto x274 = x14 * x273; auto x275 = MUs * x242; auto x276 = x262 * x275;
        auto x277 = x265 * x275; auto x278 = 6.0 * x0; auto x279 = 6.0 * AXG;
        auto x280 = x279 * x32; auto x281 = s * x139; auto x282 = x279 * x37; auto x283 = x0 * x94;
        auto x284 = 8.0 * x80; auto x285 = x16 * x284; auto x286 = x192 * x284;
        auto x287 = x146 * x81; auto x288 = AXG * x104; auto x289 = 16.0 * x266;
        auto x290 = AXG * x75; auto x291 = 16.0 * x14; auto x292 = 4.0 * x0; auto x293 = 32.0 * x81;
        auto x294 = MXs * x194; auto x295 = x283 * x39; auto x296 = MGs * UU; auto x297 = x2 * x296;
        auto x298 = x34 * x38; auto x299 = x296 * x9; auto x300 = 4.0 * MGs; auto x301 = x16 * x300;
        auto x302 = x14 * x301; auto x303 = x192 * x300; auto x304 = x14 * x303;
        auto x305 = MGs * x4; auto x306 = AXG * x305; auto x307 = MGs * x11; auto x308 = AXG * x307;
        auto x309 = x297 * x54; auto x310 = x309 * x38; auto x311 = MXs * x254;
        auto x312 = x299 * x54; auto x313 = x292 * x94; auto x314 = 4.0 * x10;
        auto x315 = -x10 * x313 - x139 * x314 + 16.0 * x174 + 16.0 * x176 + 16.0 * x178 - x3 * x313;
        auto x316 = AXG * x126; auto x317 = x131 * x22; auto x318 = x317 * x62;
        auto x319 = x104 * x190; auto x320 = x190 * x75; auto x321 = AXG * x249;
        auto x322 = AXG * x251; auto x323 = x3 * x300; auto x324 = AXG * x84;
        auto x325 = x10 * x300; auto x326 = x0 * x38; auto x327 = x326 * x49; auto x328 = -x327;
        auto x329 = x0 * x22; auto x330 = x104 * x6; auto x331 = s * x330; auto x332 = x26 * x32;
        auto x333 = s * x332; auto x334 = x202 * x4; auto x335 = x11 * x202; auto x336 = x4 * x94;
        auto x337 = x11 * x131; auto x338 = x180 * x81; auto x339 = MUs * x128;
        auto x340 = x131 * x63; auto x341 = x326 * x98; auto x342 = x131 * x94;
        auto x343 = MXs * x38; auto x344 = x24 * x343; auto x345 = x29 * x343;
        auto x346 = x344 + x345; auto x347 = x228 * x38; auto x348 = LG * LGp; auto x349 = MG * Mqi;
        auto x350 = x348 * x349; auto x351 = SS * x350; auto x352 = x132 * x14;
        auto x353 = RG * RGp; auto x354 = x349 * x353; auto x355 = SS * x354;
        auto x356 = x133 * x300; auto x357 = 4.0 * x258; auto x358 = 8.0 * x246;
        auto x359 = UU * x148; auto x360 = x350 * x359; auto x361 = x354 * x359;
        auto x362 = AXG * x143; auto x363 = x14 * x362; auto x364 = x110 * x143;
        auto x365 = MUs * UU; auto x366 = x157 * x365; auto x367 = x14 * x365;
        auto x368 = x348 * x367; auto x369 = AXG * x368; auto x370 = x164 * x349;
        auto x371 = x353 * x367; auto x372 = AXG * x371; auto x373 = 4.0 * x131;
        auto x374 = 8.0 * x6; auto x375 = x26 * x58; auto x376 = x143 * x63; auto x377 = x143 * x94;
        auto x378 = 2.0 * x213; auto x379 = x2 * x378; auto x380 = x378 * x9;
        auto x381 = AXG * x379; auto x382 = AXG * x380; auto x383 = 2.0 * x221;
        auto x384 = 2.0 * x223; auto x385 = x131 * x33; auto x386 = 4.0 * x225;
        auto x387 = x135 * x300; auto x388 = AXG * x141; auto x389 = x297 * x388;
        auto x390 = x299 * x388; auto x391 = x233 * x38; auto x392 = x236 * x38;
        auto x393 = x256 * x300; auto x394 = x258 * x300; auto x395 = MUs * x226;
        auto x396 = AXG * x161; auto x397 = x110 * x297; auto x398 = x110 * x299;
        auto x399 = x0 * x80; auto x400 = MUs * x190; auto x401 = MUs * SS; auto x402 = x190 * x401;
        auto x403 = SS * x23; auto x404 = x348 * x403; auto x405 = x143 * x202;
        auto x406 = MUs * x48; auto x407 = x348 * x406; auto x408 = AXG * x407;
        auto x409 = x198 * x349; auto x410 = x353 * x406; auto x411 = AXG * x410;
        auto x412 = AXG * x404; auto x413 = x353 * x403; auto x414 = AXG * x413;
        auto x415 = 4.0 * x6; auto x416 = x179 * x26; auto x417 = x22 * x26;
        auto x418 = x213 * x58 * x94; auto x419 = 16.0 * x6; auto x420 = x166 * x26;
        auto x421 = x135 * x32; auto x422 = AXG * x421; auto x423 = u * x131; auto x424 = x3 * x423;
        auto x425 = x10 * x423;
        auto x426 = u * x137 + u * x152 + u * x154 - x291 * x424 - x291 * x425;
        auto x427 = 32.0 * x6; auto x428 = 32.0 * AXG * x26; auto x429 = AXG * MXs;
        auto x430 = 16.0 * x429; auto x431 = 4.0 * x135; auto x432 = x3 * x431;
        auto x433 = x10 * x431; auto x434 = x141 * x147; auto x435 = x149 * x388;
        auto x436 = x263 * x62; auto x437 = x263 * x66; auto x438 = 4.0 * x256;
        auto x439 = MUs * x438; auto x440 = MUs * x357; auto x441 = x104 * x226;
        auto x442 = x226 * x75; auto x443 = 8.0 * x163; auto x444 = MXs * x443;
        auto x445 = 8.0 * x167; auto x446 = MXs * x445; auto x447 = -AXG * x207;
        auto x448 = -AXG * x209; auto x449 = x38 * x42; auto x450 = Mqis * x98;
        auto x451 = MUs * x38; auto x452 = Mqis * x49; auto x453 = x317 * x66;
        auto x454 = 4.0 * MXs; auto x455 = x24 * x454; auto x456 = x29 * x454;
        auto x457 = x455 + x456; auto x458 = x164 * x197 + x164 * x199 - x179 * x292 - x314 * x329 +
                                             x388 * x49 + x388 * x98 - x453 + x457;
        auto x459 = 2.0 * x256; auto x460 = 2.0 * x258; auto x461 = x14 * x63;
        auto x462 = x165 * x38; auto x463 = x462 + x52 * x63; auto x464 = s * x421;
        auto x465 = s * x230; auto x466 = s * x233; auto x467 = s * x236; auto x468 = AXG * x180;
        auto x469 = x104 * x131; auto x470 = MXs * x469; auto x471 = s * x470;
        auto x472 = MUs * x225; auto x473 = MXs * x472; auto x474 = s * x473;
        auto x475 = x20 * x297; auto x476 = 2.0 * s; auto x477 = x26 * x297; auto x478 = x20 * x299;
        auto x479 = x26 * x299; auto x480 = AXG * x475; auto x481 = AXG * x477;
        auto x482 = AXG * x478; auto x483 = AXG * x479; auto x484 = x194 * x365;
        auto x485 = x3 * x373; auto x486 = x321 + x322; auto x487 = AXG * x150;
        auto x488 = x487 * x94; auto x489 = AXG * x144; auto x490 = x489 * x94;
        auto x491 = x447 + x448 + x488 + x490; auto x492 = MGs * syFC9; auto x493 = 3.0 * AXG;
        auto x494 = x493 * x98; auto x495 = x49 * x493; auto x496 = AXG * x336;
        auto x497 = s * (x172 * x493 + x172 + x334 + x335 + x336 - x49 + x493 * x95 - x494 - x495 -
                         x496 - x98);
        auto x498 = 4.0 * x73; auto x499 = 4.0 * x76; auto x500 = 12.0 * x99;
        auto x501 = 12.0 * x50; auto x502 = x33 * x64; auto x503 = AXG * x502 - x51;
        auto x504 = x1 * x32; auto x505 = -4.0 * x228 + x347 + x38 * x504 + x39 * x64 - 4.0 * x504;
        auto x506 = 6.0 * x131; auto x507 = 12.0 * x80; auto x508 = 14.0 * x80;
        auto x509 = x273 * x98; auto x510 = x273 * x49; auto x511 = 6.0 * x80;
        auto x512 = MXs * x99; auto x513 = s * syFC12; auto x514 = x139 * x33;
        auto x515 = -x338 - x339; auto x516 = x150 * x22; auto x517 = x144 * x22;
        auto x518 = -x516 - x517; auto x519 = x172 * x292 - x173 * x38 + 4.0 * x175;
        auto x520 = x292 * x49;
        auto x521 = MXs * x100 + MXs * x51 - x292 * x98 + x327 + x341 - x520;
        auto x522 = 3.0 * Mqis; auto x523 = x37 * x522; auto x524 = Mqis * x493;
        auto x525 = 3.0 * x110; auto x526 = x299 * x493;
        auto x527 = AXG * x181 - x181 + x186 * x38 + x498 + x499 - x74 - x77; auto x528 = s * syFC3;
        auto x529 = AXG * x424; auto x530 = AXG * x425; auto x531 = 4.0 * x64; auto x532 = UU * x1;
        auto x533 = x348 * x532; auto x534 = 8.0 * x349; auto x535 = x353 * x532;
        auto x536 = x3 * x64; auto x537 = AXG * x536; auto x538 = u * x143; auto x539 = x351 * x538;
        auto x540 = x355 * x538; auto x541 = x348 * x64; auto x542 = UU * x541;
        auto x543 = x353 * x64; auto x544 = UU * x543; auto x545 = x38 * x64;
        auto x546 = AXG * x533; auto x547 = AXG * x535; auto x548 = SS * x541;
        auto x549 = SS * x543; auto x550 = x349 * x58; auto x551 = -x539 - x540;
        auto x552 = s * syFC6; auto x553 = 3.0 * x98; auto x554 = 3.0 * x49; auto x555 = 7.0 * MXs;
        auto x556 = MXs * x50; auto x557 = 12.0 * x349; auto x558 = x349 * x38;
        auto x559 = s * syFC8; auto x560 = x16 * x161; auto x561 = MUs * x297;
        auto x562 = AXG * x155; auto x563 = MUs * x299; auto x564 = AXG * x156;
        auto x565 = 8.0 * x350; auto x566 = SS * x565; auto x567 = 8.0 * x354;
        auto x568 = SS * x567; auto x569 = 16.0 * x54; auto x570 = x396 * x401;
        auto x571 = 16.0 * x349; auto x572 = x365 * x6; auto x573 = x348 * x572;
        auto x574 = AXG * x573; auto x575 = 32.0 * x349; auto x576 = UU * x26;
        auto x577 = x348 * x576; auto x578 = AXG * x577; auto x579 = x353 * x572;
        auto x580 = AXG * x579; auto x581 = x353 * x576; auto x582 = AXG * x581;
        auto x583 = -x393 - x394; auto x584 = syFC8 * u; auto x585 = x0 * x300;
        auto x586 = MXs * x301; auto x587 = MXs * x303; auto x588 = SS * x141;
        auto x589 = x110 * x166; auto x590 = MGs * x101; auto x591 = -x438;
        auto x592 = x562 + x564 + x591; auto x593 = -x151 - x153 + x432 + x433 + x436 + x437 + x439;
        auto x594 = x131 * x38; auto x595 = s * syFC11; auto x596 = u * x215; auto x597 = u * x216;
        auto x598 = 4.0 * x50; auto x599 = AXG * x385; auto x600 = AXG * x386;
        auto x601 = -x596 - x597; auto x602 = x0 * x279; auto x603 = 14.0 * MXs;
        auto x604 = x261 * x273; auto x605 = x264 * x273; auto x606 = 3.0 * x6;
        auto x607 = 3.0 * x14; auto x608 = Mqis * x172; auto x609 = Mqis * x95;
        auto x610 = x297 * x493; auto x611 = Mqis * x94; auto x612 = Mqis * x63;
        auto x613 = x297 * x94; auto x614 = x299 * x94; auto x615 = 2.0 * x424;
        auto x616 = 2.0 * x425; auto x617 = 3.0 * x63; auto x618 = x25 + x30 - x615 - x616;
        auto x619 = AXG * x504; auto x620 = AXG * x68; auto x621 = -x197 - x199;
        auto x622 = AXG * x256; auto x623 = AXG * x258; auto x624 = MXs * x143;
        auto x625 = x143 * x429; auto x626 = -x357; auto x627 = AXG * x587; auto x628 = 4.0 * x16;
        auto x629 = MUs * x314; auto x630 = x161 * x192; auto x631 = Mqis * x14;
        auto x632 = x58 * x6; auto x633 = AXG * x330;
        auto x634 = AXG * x221 + AXG * x223 - AXG * x470 - AXG * x473 - x221 - x223 + x230 - x231 -
                    x233 + x234 - x236 + x237 + x421 - x422 + x470 + x473;
        auto x635 = x256 * x38; auto x636 = x258 * x38; auto x637 = s * x492; auto x638 = x63 * x9;
        auto x639 = 2.0 * MGs; auto x640 = x54 * x639; auto x641 = x2 * x640; auto x642 = x640 * x9;
        auto x643 = x2 * x63; auto x644 = 3.0 * x283; auto x645 = 4.0 * x349;
        auto x646 = x348 * x645; auto x647 = x353 * x645; auto x648 = x493 * x6;
        auto x649 = x283 * x493; auto x650 = MUs * x15; auto x651 = x38 * x54;
        auto x652 = x350 * x38; auto x653 = x354 * x38; auto x654 = MXs * x385;
        auto x655 = MUs * x386; auto x656 = MXs * x386; auto x657 = 2.0 * x0;
        auto x658 = x297 * x657; auto x659 = x299 * x657; auto x660 = x131 * x300;
        auto x661 = 6.0 * x110; auto x662 = UU * x292; auto x663 = 2.0 * AXG; auto x664 = UU * x326;
        auto x665 = x48 * x58; auto x666 = SS * x22;
        auto x667 = AXG * x247 + AXG * x248 - x565 * x666; auto x668 = 7.0 * x80;
        auto x669 = 6.0 * MGs; auto x670 = 3.0 * x81; auto x671 = AXG * x639;
        auto x672 = 12.0 * x48; auto x673 = s * syFC7 * x300; auto x674 = x10 * x64;
        auto x675 = SS * x242; auto x676 = x348 * x675; auto x677 = x353 * x675;
        auto x678 = s * syFC5 * x645; auto x679 = SS * x423; auto x680 = x348 * x679;
        auto x681 = x353 * x679;
        ,
        -Nc * TR * TR * gs * gs * (Nc - 1) * (Nc + 1) * (L * Rp + Lp * R) *
            (-MUs * x552 *
                 (x124 + x127 + x129 + x202 * x566 + x202 * x568 - x247 - x248 - x350 * x665 -
                  x354 * x665 + x48 * x565 + x48 * x567 + x667) -
             MUs * x559 *
                 (-x10 * x511 - x179 * x300 - x3 * x511 + x32 * x668 - x32 * x670 + x350 * x672 +
                  x354 * x672 + x37 * x668 - x37 * x670 + x48 * x652 + x48 * x653 + x49 * x669 +
                  x49 * x671 - x567 * x666 + x667 + x669 * x98 + x671 * x98) +
             UU * x559 *
                 (AXG * x641 + AXG * x642 - x14 * x493 * x638 - x2 * x644 + x2 * x649 + x20 * x646 +
                  x20 * x647 - x20 * x652 - x20 * x653 + x26 * x646 + x26 * x647 - x26 * x652 -
                  x26 * x653 - x350 * x650 + x350 * x651 - x354 * x650 + x354 * x651 - x606 * x638 -
                  x606 * x643 + x607 * x638 + x638 * x648 - x641 - x642 + x643 * x648 - x644 * x9 +
                  x649 * x9) +
             s * syFC10 *
                 (-u * x144 - u * x150 - x100 + x206 + x208 - x498 - x499 - x500 - x501 + x503 +
                  x505 + 12.0 * x68 + 12.0 * x71 + x78) +
             s * syFC13 *
                 (-AXG * x514 + MXs * x500 + MXs * x501 - x104 * x507 - x139 * x39 - x15 * x70 -
                  4.0 * x165 + x285 + x286 + x463 - x507 * x75 + x514 + x515 + x518 + x519 + x521) -
             s * syFC4 *
                 (x163 + x167 - x228 + x229 + x241 - 3.0 * x50 - x504 + x618 + x619 + x620 + x621 +
                  3.0 * x68 + 3.0 * x71 - x73 - 3.0 * x99) +
             syFC1 * x130 + syFC1 * x497 + syFC1 * x79 +
             syFC10 *
                 (-AXG * x182 - AXG * x185 - x104 * x194 + x134 * x195 - x141 * x179 - x145 * x22 +
                  x148 * x49 + x148 * x98 + x150 * x63 - x164 * x200 - x164 * x201 + x182 + x185 -
                  x187 * x38 + 4.0 * x187 - x189 * x38 + 4.0 * x189 + x191 + x193 - x194 * x75 +
                  x196 * x3 + x197 * x198 + x198 * x199 + x205 + x212) +
             syFC10 *
                 (AXG * x162 - MUs * x155 - MUs * x156 + x10 * x140 - x104 * x157 + x132 * x133 -
                  x133 * x158 + x134 * x14 - x137 - x138 + x140 * x3 + x142 * x94 + x144 * x63 +
                  x145 * x94 + x146 * x147 + x148 * x149 - x152 - x154 - x157 * x75 - x158 * x159 +
                  x160 * x58 + x163 * x164 + x164 * x167 + x165 * x166 + 16.0 * x169 - 32.0 * x171 -
                  32.0 * x174 - 32.0 * x176 - 32.0 * x178) +
             syFC12 * (-MXs * x191 - MXs * x193 - x0 * x247 - x0 * x248 + x105 * x194 +
                       x107 * x194 - x143 * x197 - x143 * x199 + x161 * x197 + x190 * x225 +
                       x202 * x215 + x202 * x216 + x24 * x254 + x250 + x252 + x253 * x32 +
                       x253 * x37 - x255 * x49 - x255 * x98 + x259) -
             syFC12 *
                 (8.0 * u * x234 + u * x238 - x161 * x199 + x161 * x24 + x161 * x29 + x197 * x419 +
                  x199 * x419 + x2 * x418 + x215 * x63 - x215 * x94 + x216 * x63 - x216 * x94 -
                  x232 * x422 - x24 * x415 - x254 * x29 - x29 * x415 - x314 * x417 - 4.0 * x416 +
                  x418 * x9 + x420 * x49 + x420 * x98 + x426) +
             syFC12 * (AXG * x218 + AXG * x219 - MUs * x238 + x1 * x45 * x58 + x101 * x165 -
                       x116 * x226 - x143 * x163 + x158 * x244 + x158 * x246 + x161 * x229 +
                       x161 * x241 - x2 * x217 + x215 * x6 + x216 * x6 - x217 * x9 + x221 * x222 +
                       x222 * x223 + x224 * x35 + x225 * x226 - x226 * x70 - x227 + x231 * x232 -
                       x234 * x235 - x235 * x237 - x239 - x240) +
             syFC12 * (AXG * x268 + AXG * x269 + AXG * x276 + AXG * x277 - s * x172 * x278 +
                       s * x173 * x279 + s * x282 * x283 - x104 * x270 + x168 * x272 - x170 * x272 -
                       x175 * x272 + x224 * x81 + x262 * x263 - x262 * x274 + x263 * x265 -
                       x265 * x274 + x267 * x62 + x267 * x66 - x268 - x269 - x270 * x75 + x271 -
                       x276 - x277 - x278 * x96 + x280 * x281 + x281 * x282) +
             syFC13 * (MXs * x285 + MXs * x286 - s * x295 + x105 * x293 + x107 * x293 - x142 * x80 -
                       x145 * x80 + x16 * x270 - x16 * x294 + x161 * x200 + x192 * x270 -
                       x192 * x294 + x197 * x291 - x24 * x291 + x259 + x271 + x287 * x32 +
                       x287 * x37 - x288 * x289 - x289 * x290 + x292 * x96) -
             syFC13 * (-16.0 * x1 * x147 - 8.0 * x10 * x417 + x116 * x157 + x157 * x70 -
                       x161 * x201 - x165 * x430 - x168 * x430 + x197 * x427 - x199 * x291 +
                       x199 * x427 - x206 * x6 - x208 * x6 + x227 - x229 * x291 + x239 + x240 +
                       x29 * x291 - 8.0 * x416 + x426 + x428 * x49 + x428 * x98) -
             syFC2 * x130 - syFC2 * x497 - syFC2 * x79 -
             syFC3 * u *
                 (AXG * x586 + Mqis * x18 + Mqis * x41 - Mqis * x46 - Mqis * x52 + x0 * x305 +
                  x0 * x307 - x101 * x561 - x101 * x563 + x131 * x629 + x16 * x373 - x297 * x632 -
                  x299 * x632 + x305 * x6 + x307 * x6 - x43 * x631 - x44 * x631 + x443 + x445 +
                  x54 * x62 + x54 * x66 - x560 + x592 + x6 * x628 + x626 + x627 - x630 -
                  16.0 * x633) -
             syFC3 * (Mqis * x109 + Mqis * x114 + Mqis * x115 + Mqis * x25 + Mqis * x30 +
                      Mqis * x57 - Mqis * x85 - Mqis * x86 - u * x356 - x101 * x450 - x101 * x452 +
                      x297 * x449 + x299 * x449 + x306 * x42 + x308 * x42 - x325 * x42 + x447 +
                      x448 - x450 * x451 - x451 * x452 + x458) +
             syFC3 *
                 (MGs * x21 + MGs * x27 + MGs * x28 + MGs * x31 + Mqis * x12 + Mqis * x13 +
                  Mqis * x5 - Mqis * x60 - Mqis * x61 + Mqis * x65 + Mqis * x67 - Mqis * x69 -
                  Mqis * x72 + Mqis * x8 + x297 * x298 - x297 * x311 + x298 * x299 - x299 * x311 -
                  x302 - x304 + x306 * x54 + x308 * x54 - x310 - x312 * x38 + x315) -
             syFC3 * (-Mqis * x17 - Mqis * x19 + Mqis * x21 + Mqis * x27 + Mqis * x28 + Mqis * x31 +
                      Mqis * x36 + Mqis * x40 - Mqis * x47 - Mqis * x53 + Mqis * x55 + Mqis * x56 -
                      Mqis * x74 - Mqis * x77 + x151 + x153 - x432 - x433 + x434 + x435 - x436 -
                      x437 - x439 - x440 - x441 - x442 + x444 + x446) +
             syFC3 *
                 (-MGs * x102 - MGs * x103 + MGs * x109 + MGs * x25 + MGs * x30 + MGs * x57 -
                  MGs * x82 - MGs * x83 - MGs * x85 - MGs * x86 + Mqis * x108 - Mqis * x112 -
                  Mqis * x113 + Mqis * x120 + Mqis * x121 - Mqis * x123 - Mqis * x125 -
                  Mqis * x128 + Mqis * x82 + Mqis * x83 - x197 * x300 - x199 * x300 + x297 * x316 +
                  x299 * x316 + x318 + x319 + x320 - x321 - x322 + x323 * x324 + x324 * x325) -
             syFC4 * (-MUs * x459 - MUs * x460 + MXs * x74 + MXs * x77 + x0 * x336 + x0 * x46 +
                      x0 * x52 - x11 * x135 + x11 * x139 - x11 * x274 + x11 * x283 + x131 * x18 +
                      x131 * x41 - x135 * x4 + x139 * x4 - x170 * x58 - x172 * x388 - x175 * x58 -
                      x18 * x63 - x180 * x35 - x274 * x4 - x298 * x75 + x337 * x63 - x388 * x95 -
                      x41 * x63 + x43 * x461 + x44 * x461 + x463) +
             syFC4 * (AXG * x187 + AXG * x189 + AXG * x331 + AXG * x333 - MUs * x82 - MUs * x83 -
                      MXs * x25 - MXs * x30 + s * x76 + x11 * x329 + x131 * x334 + x131 * x335 +
                      x131 * x336 - x187 - x189 - x197 * x235 - x199 * x235 + x328 + x329 * x4 -
                      x331 - x333 + x337 * x94 + x338 + x339 - x340 * x4 - x341 - x342 * x62 -
                      x342 * x66 + x346) +
             syFC8 * x22 *
                 (-x10 * x585 + x110 * x568 - x16 * x590 - x192 * x590 - x224 * x429 - x273 * x323 -
                  x273 * x325 + x297 * x589 + x299 * x589 - x3 * x585 - x350 * x588 - x354 * x588 -
                  x358 + x360 + x361 + x381 + x382 + x389 + x390 + x586 + x587 + x592) +
             syFC8 *
                 (MGs * x151 + x1 * x46 - x136 * x351 - x136 * x355 + x14 * x347 + x14 * x360 +
                  x14 * x361 + x165 * x343 - x273 * x356 - x350 * x366 + x351 * x352 - x351 * x363 -
                  x351 * x364 + x352 * x355 - x354 * x366 - x355 * x363 - x355 * x364 - x357 * x63 +
                  x358 * x63 - x36 * x63 + x369 * x370 + x370 * x372 - x40 * x63 + x52 * x64) +
             syFC8 * (-x100 * x131 + x100 * x14 + x109 * x14 - x131 * x51 + x14 * x51 + x14 * x57 -
                      x15 * x24 - x15 * x29 + x196 * x351 + x196 * x355 - x197 * x374 -
                      x199 * x374 + x22 * x27 + x22 * x31 - x22 * x357 + x24 * x373 + x25 * x6 +
                      x29 * x373 + x30 * x6 + x351 * x376 - x351 * x377 + x355 * x376 -
                      x355 * x377 - x375 * x49 - x375 * x98) +
             syFC8 * (-AXG * x250 - AXG * x252 + MGs * x321 + MGs * x322 + x0 * x128 + x105 * x190 +
                      x107 * x190 + x11 * x110 * x80 - x11 * x399 + x111 * x80 + x183 * x81 +
                      x235 * x349 * x404 - x247 * x273 - x248 * x273 - x297 * x400 - x299 * x400 -
                      x351 * x405 + x354 * x402 - x355 * x405 - x370 * x412 - x370 * x414 +
                      x385 * x81 + x386 * x81 - x399 * x4 + x408 * x409 + x409 * x411) +
             syFC8 * (MGs * x153 + MUs * x383 + MUs * x384 - MUs * x392 - MUs * x393 - MUs * x394 -
                      MXs * x391 - MXs * x392 - x10 * x387 - x131 * x74 + x14 * x381 + x14 * x382 +
                      x14 * x389 + x14 * x390 - x159 * x273 * x300 - x218 - x219 - x297 * x395 -
                      x299 * x395 - x3 * x387 + x35 * x385 + x35 * x386 + x379 * x6 + x380 * x6 +
                      x396 * x397 + x396 * x398) -
             syFC8 * (-AXG * x464 - AXG * x465 + AXG * x466 + AXG * x467 - AXG * x471 - AXG * x474 +
                      x266 * x38 * x75 + x266 * x468 - x267 * x43 - x267 * x44 - x350 * x402 +
                      x350 * x484 + x354 * x484 + x464 + x465 - x466 - x467 + x471 + x474 -
                      x475 * x476 - x476 * x477 - x476 * x478 - x476 * x479 + x476 * x480 +
                      x476 * x481 + x476 * x482 + x476 * x483) +
             u * x552 *
                 (-AXG * x215 - AXG * x216 + AXG * x654 + AXG * x655 + AXG * x656 - x131 * x180 +
                  x155 + x156 + x180 * x273 + x215 + x216 + x292 * x297 + x292 * x299 -
                  x297 * x326 - x299 * x326 + x527 - x562 - x564 - x587 + x627 - x654 - x655 -
                  x656) +
             u * x559 *
                 (-AXG * x469 - AXG * x472 + AXG * x658 + AXG * x659 - x10 * x660 - 3.0 * x243 -
                  x244 - 3.0 * x245 - x246 + x297 * x661 + x299 * x661 - x3 * x660 - x350 * x662 +
                  x350 * x664 - x354 * x662 + x354 * x664 + x357 + x379 + x380 + x397 * x663 +
                  x398 * x663 + x438 + x469 + x472 - x493 * x73 - x635 - x636 - x658 - x659 +
                  3.0 * x73) +
             x22 * x528 * (x180 + x297 + x299 - x306 - x308 + x468 + x526 + x610) +
             x22 * x552 *
                 (-AXG * x183 - AXG * x301 - AXG * x303 + MXs * x180 - MXs * x468 - MXs * x628 -
                  MXs * x629 - x107 * x38 + 4.0 * x107 + x110 * x62 + x110 * x66 - x144 - x150 +
                  x183 - x222 * x297 - x222 * x299 + x297 * x451 + x299 * x451 + x301 + x303 +
                  x385 + x386 + x401 * x567 + x487 + x489 - x599 - x600) -
             x23 * x637 * (x261 + x264 + x280 + x282 - x62 - x66) +
             x492 * (x165 * x58 - x189 * x279 + 6.0 * x189 + x314 * x340 - x314 * x342 - x318 -
                     x319 - x320 + x458 + x461 * x62 + x485 * x63 - x485 * x94 + x486 + x491) -
             x492 * (-4.0 * x139 * x3 + 8.0 * x160 + x162 + 16.0 * x171 + x315 - x396 * x70 - x434 -
                     x435 + x440 + x441 + x442 - x444 - x446 - x461 * x66 + x593) +
             x513 * (-x104 * x508 - x132 * x179 - x134 * x22 + x16 * x507 - x191 + x192 * x507 -
                     x193 + x203 + x204 - x279 * x512 + x288 * x511 + x290 * x511 + x49 * x506 +
                     x506 * x98 - x508 * x75 + 2.0 * x509 + 2.0 * x510) -
             x513 * (x165 * x279 - 6.0 * x165 + x172 * x506 + x212 + 6.0 * x225 * x94 -
                     x261 * x340 - x264 * x340 - x278 * x49 - x278 * x98 + x279 * x556 +
                     x282 * x461 + x49 * x602 + x491 - x50 * x603 + x601 + x602 * x98 - x603 * x99 +
                     x604 * x63 + x604 * x94 + x605 * x63 + x605 * x94) +
             x528 * (AXG * x184 - AXG * x459 - AXG * x460 + x105 * x522 - x105 * x524 +
                     x110 * x305 + x110 * x307 - x110 * x308 + x110 * x526 - x14 * x523 +
                     x149 * x524 - x184 - 4.0 * x186 + x188 * x38 - 4.0 * x188 - x299 * x525 -
                     x37 * x524 * x6 - x45 * x522 + x45 * x524 + x459 + x460 + x522 * x88 +
                     x523 * x6 - x524 * x88 + x527) -
             x528 *
                 (-Mqis * x336 + Mqis * x496 - x11 * x611 + x11 * x64 + x110 * x306 - x110 * x610 +
                  x116 * x522 - x116 * x524 + x14 * x526 + x14 * x610 + x297 * x525 + x297 * x606 -
                  x297 * x607 + x299 * x606 - x299 * x607 + x43 * x612 + x44 * x611 + x44 * x612 -
                  x44 * x64 - x493 * x608 - x493 * x609 + x505 + x523 * x63 - x524 * x70 -
                  x526 * x6 - x6 * x610 - x608 - x609 + 4.0 * x71) -
             x528 * (-AXG * x615 - AXG * x616 + MGs * x336 - MGs * x496 - Mqis * x334 -
                     Mqis * x335 + x109 - x297 * x617 - x299 * x617 - x306 * x63 + x307 * x94 -
                     x308 * x63 - x308 * x94 + x4 * x64 - x43 * x64 + x450 * x493 + x450 +
                     x452 * x493 + x452 + x493 * x613 + x493 * x614 + x502 + x503 + x526 * x63 +
                     x57 - x598 + x610 * x63 + x613 + x614 + x618) +
             4.0 * x552 *
                 (-AXG * x312 + x309 + x312 - x475 - x477 - x478 - x479 + x480 + x481 + x482 +
                  x483 + x634) -
             x552 * (x292 * x95 - x295 + x300 * x622 + x300 * x623 + x310 - x351 * x624 +
                     x351 * x625 - x355 * x624 + x355 * x625 - x368 * x534 + x369 * x534 -
                     x371 * x534 + x372 * x534 + x519 + x534 * x573 - x534 * x574 + x534 * x577 -
                     x534 * x578 + x534 * x579 - x534 * x580 + x534 * x581 - x534 * x582 + x583) +
             x552 * (AXG * x539 + AXG * x540 + x297 * x531 - x297 * x545 + x299 * x531 -
                     x299 * x545 - x300 * x424 - x300 * x425 + x300 * x529 + x300 * x530 -
                     x300 * x536 + x300 * x537 + x328 + x520 + x533 * x534 + x534 * x535 +
                     x534 * x542 + x534 * x544 - x534 * x546 - x534 * x547 - x534 * x548 -
                     x534 * x549 - x542 * x550 - x544 * x550 + x548 * x550 + x549 * x550 + x551) +
             x559 * (-x0 * x494 - x0 * x495 + x0 * x553 + x0 * x554 + x131 * x553 + x131 * x554 +
                     x29 * x300 + x318 + x346 + x453 - x455 - x456 - x493 * x512 - x493 * x556 +
                     x50 * x555 + x509 + x510 + x518 + x542 * x557 + x542 * x558 + x544 * x557 +
                     x544 * x558 + x551 + x555 * x99) -
             x584 * (MGs * x630 + MUs * x379 + MUs * x380 - x105 * x362 - x20 * x325 + x221 * x38 +
                     x223 * x38 - 8.0 * x231 - x26 * x323 - x26 * x325 - x301 * x6 - x383 - x384 +
                     x391 + x392 - 8.0 * x422 + 16.0 * x480 + 16.0 * x481 + 16.0 * x482 +
                     16.0 * x483 + x593) +
             x584 * (AXG * x302 + AXG * x304 - MGs * x560 + MGs * x562 + MGs * x564 + x131 * x301 +
                     x131 * x303 + x20 * x566 + x20 * x568 + x26 * x566 + x26 * x568 + x350 * x570 -
                     x351 * x569 + x354 * x570 - x355 * x569 + x369 * x571 + x372 * x571 +
                     x396 * x561 + x396 * x563 - x574 * x575 - x575 * x578 - x575 * x580 -
                     x575 * x582 + x583) +
             4.0 * x595 *
                 (x165 + x168 - x169 - x170 + x171 - x173 + x174 - x175 + x176 - x177 + x178 +
                  x634) -
             x595 * (MUs * x127 + x180 * x80 + x205 - x249 - x251 - x373 * x49 - x373 * x98 + x486 +
                     x49 * x594 + x515 + x516 + x517 + x594 * x98) -
             x595 * (AXG * x596 + AXG * x597 - MXs * x598 + x172 * x373 - x172 * x594 + x210 +
                     x211 - x344 - x345 + x385 * x63 + x386 * x63 + x386 * x94 - x454 * x99 + x457 +
                     x462 + x488 + x490 + x521 - x599 * x63 - x600 * x63 - x600 * x94 + x601) +
             x637 * (-x186 * x279 + 6.0 * x186 + x228 * x279 - 6.0 * x228 + x261 * x64 +
                     x264 * x64 - x279 * x330 - x279 * x332 + x279 * x504 + x279 * x73 +
                     x279 * x76 + x280 * x64 + x282 * x64 + x314 * x64 + 6.0 * x330 + 6.0 * x332 -
                     x38 * x424 - x38 * x425 - 6.0 * x504 + 4.0 * x536 + x591 + x626 + x635 + x636 -
                     x65 - x67 - 6.0 * x73 - 6.0 * x76) +
             x673 * (x200 + x201 - x24 - x29 + x424 + x425 + x50 - x529 - x530 + x621 + x99) -
             x673 * (AXG * x186 + AXG * x188 + AXG * x332 + AXG * x674 - x163 - x167 - x186 - x188 +
                     x228 - x229 - x241 + x256 + x258 - x330 - x332 + x504 - x536 + x537 - x619 -
                     x620 - x622 - x623 + x633 - x674 + x68 + x71 + x73 + x76) -
             x678 * (AXG * x542 - AXG * x680 - AXG * x681 - x404 + x407 - x408 + x410 - x411 +
                     x412 - x413 + x414 - x542 + x680 + x681) +
             x678 * (-AXG * x544 + AXG * x548 + AXG * x549 - AXG * x676 - AXG * x677 + x368 - x369 +
                     x371 - x372 + x533 + x535 + x544 - x546 - x547 - x548 - x549 - x573 + x574 -
                     x577 + x578 - x579 + x580 - x581 + x582 + x676 + x677)) *
            Denom(-1536.0 * s * (MQis * MUs - MQis * u + x22 - x63) * pow(MU, 2) * pow(Pi, 2)));
  //}
}
return ret.real();
}

ComplexType ME_us_QQg_Gqq(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  auto Denom = [](auto a) { return 1. / a; };
  ComplexType ret = 0;
  int itq = q;
  for (int itsq = 0; itsq < 2; itsq++) {
    // for (int itq = 0; itq < 3; itq++) {

    int isq = (is_chargino(ch) != is_up_quark(itq)) * 6 + itsq * 3 + itq -
              ( is_up_quark(itq)) * 3;
    int iq = (itq + is_chargino(ch) * 3) % 6;

    ComplexType L = (params->CHSQq[ch][isq][iq].L);
    ComplexType R = (params->CHSQq[ch][isq][iq].R);
    ComplexType Lp = conj(params->CHSQq[ch][sq][iq].R);
    ComplexType Rp = conj(params->CHSQq[ch][sq][iq].L);

    auto MQi = params->mSQ[isq];
    auto MQis = pow2(MQi);

    auto Mqi = params->mq[iq];
    auto Mqis = pow2(Mqi);

    ComplexType LG = (params->GLSQq[isq][q].L);
    ComplexType RG = (params->GLSQq[isq][q].R);
    ComplexType LGp = conj(params->GLSQq[sq][q].R);
    ComplexType RGp = conj(params->GLSQq[sq][q].L);

#define syFC1 A0(MGs)
#define syFC2 A0(Mqis)
#define syFC3 B0(MUs, MGs, Mqis)
#define syFC4 B0(u, MGs, Mqis)
#define syFC5 C0(0, u, MUs, Mqis, Mqis, MGs)
#define syFC6 C1(0, MUs, u, Mqis, Mqis, MGs)
#define syFC7 C2(0, MUs, u, Mqis, Mqis, MGs)
#define syFC8 C00(0, MUs, u, Mqis, Mqis, MGs)
#define syFC9 C11(0, MUs, u, Mqis, Mqis, MGs)
#define syFC10 C12(0, MUs, u, Mqis, Mqis, MGs)
#define syFC11 C22(0, MUs, u, Mqis, Mqis, MGs)

    _EPS0_(
        ret, auto x0 = LG * RGp; auto x1 = pow(u, 3); auto x2 = UU * x1; auto x3 = x0 * x2;
        auto x4 = AXG * x3; auto x5 = 4. * x4; auto x6 = LGp * RG; auto x7 = x2 * x6;
        auto x8 = AXG * x7; auto x9 = 4. * x8; auto x10 = 3. * s; auto x11 = pow(MX, 4);
        auto x12 = UU * x0; auto x13 = x11 * x12; auto x14 = x10 * x13; auto x15 = UU * x6;
        auto x16 = x11 * x15; auto x17 = x10 * x16; auto x18 = SS * x1; auto x19 = x0 * x18;
        auto x20 = AXG * x19; auto x21 = 2. * x20; auto x22 = x18 * x6; auto x23 = AXG * x22;
        auto x24 = 2. * x23; auto x25 = 2. * x0; auto x26 = pow(u, 2); auto x27 = MXs * SS;
        auto x28 = x26 * x27; auto x29 = x25 * x28; auto x30 = 2. * x6; auto x31 = x28 * x30;
        auto x32 = pow(MXs, 2); auto x33 = x0 * x32; auto x34 = UU * x10; auto x35 = x33 * x34;
        auto x36 = x32 * x6; auto x37 = x34 * x36; auto x38 = MXs * u; auto x39 = x15 * x38;
        auto x40 = 4. * AXG; auto x41 = x28 * x40; auto x42 = x0 * x41; auto x43 = x41 * x6;
        auto x44 = AXG * x35; auto x45 = AXG * x37; auto x46 = MUs * u; auto x47 = x10 * x46;
        auto x48 = x12 * x47; auto x49 = x15 * x47; auto x50 = MUs * SS; auto x51 = x26 * x50;
        auto x52 = x0 * x51; auto x53 = 2. * x52; auto x54 = AXG * x53; auto x55 = -x54;
        auto x56 = MXs * x50; auto x57 = x0 * x56; auto x58 = 2. * s; auto x59 = x56 * x6;
        auto x60 = SS * u; auto x61 = MXs * x60; auto x62 = x30 * x61; auto x63 = s * x62;
        auto x64 = MUs * x12; auto x65 = MXs * x64; auto x66 = MUs * MXs; auto x67 = x15 * x66;
        auto x68 = x10 * x67; auto x69 = UU * x26; auto x70 = x0 * x69; auto x71 = MUs * x70;
        auto x72 = AXG * x71; auto x73 = 4. * x72; auto x74 = AXG * MXs; auto x75 = 8. * x74;
        auto x76 = x70 * x75; auto x77 = x6 * x69; auto x78 = x75 * x77; auto x79 = AXG * x65;
        auto x80 = AXG * x67; auto x81 = x10 * x80; auto x82 = 2. * SS; auto x83 = x0 * x82;
        auto x84 = x46 * x83; auto x85 = AXG * s; auto x86 = x46 * x82; auto x87 = x6 * x86;
        auto x88 = AXG * x57; auto x89 = AXG * x59;
        auto x90 = AXG * x14 + AXG * x17 + AXG * x48 + AXG * x49 - AXG * x63 + s * x39 + x10 * x65 -
                   x10 * x79 - x14 - x17 + x21 + x24 + x29 + x31 + x35 + x37 - x42 - x43 - x44 -
                   x45 - x48 - x49 - x5 + x55 - x57 * x58 - x58 * x59 + x58 * x88 + x58 * x89 +
                   x63 + x68 + x73 + x76 + x78 - x81 - x84 * x85 - x85 * x87 - x9;
        auto x91 = MGs * UU; auto x92 = x0 * x91; auto x93 = x1 * x40; auto x94 = x6 * x91;
        auto x95 = x51 * x6; auto x96 = 2. * x95; auto x97 = MUs * x3; auto x98 = AXG * x97;
        auto x99 = 8. * x98; auto x100 = MUs * x7; auto x101 = AXG * x100; auto x102 = 8. * x101;
        auto x103 = AXG * x95; auto x104 = 2. * x103; auto x105 = 2. * Mqis; auto x106 = x105 * x57;
        auto x107 = x105 * x59; auto x108 = MUs * x77; auto x109 = 4. * x108;
        auto x110 = AXG * x109; auto x111 = x105 * x88; auto x112 = x105 * x89;
        auto x113 = MUs * x19; auto x114 = AXG * x113; auto x115 = 4. * x114; auto x116 = MUs * x22;
        auto x117 = AXG * x116; auto x118 = 4. * x117; auto x119 = 4. * SS; auto x120 = x119 * x26;
        auto x121 = pow(MU, 4); auto x122 = x0 * x121; auto x123 = AXG * x122;
        auto x124 = -x115 - x118 + x120 * x123; auto x125 = 2. * x61; auto x126 = x122 * x125;
        auto x127 = MUs * x11; auto x128 = x127 * x15; auto x129 = MUs * UU; auto x130 = x129 * x33;
        auto x131 = s * x130; auto x132 = pow(MUs, 2); auto x133 = UU * x132; auto x134 = x0 * x133;
        auto x135 = MXs * x134; auto x136 = s * x135; auto x137 = x129 * x36; auto x138 = s * x137;
        auto x139 = x133 * x6; auto x140 = MXs * x139; auto x141 = s * x140; auto x142 = 2. * x113;
        auto x143 = AXG * x142; auto x144 = 2. * x116; auto x145 = AXG * x144;
        auto x146 = x122 * x82; auto x147 = x146 * x46; auto x148 = x121 * x6;
        auto x149 = x125 * x148; auto x150 = 4. * x98; auto x151 = 4. * x101;
        auto x152 = AXG * x131; auto x153 = AXG * x134; auto x154 = MXs * x153;
        auto x155 = s * x154; auto x156 = AXG * x138; auto x157 = AXG * x139;
        auto x158 = MXs * x157; auto x159 = s * x158; auto x160 = 4. * u; auto x161 = x160 * x27;
        auto x162 = x123 * x161; auto x163 = AXG * x160; auto x164 = x148 * x163 * x27;
        auto x165 = x146 * x26; auto x166 = AXG * x165; auto x167 = SS * x148;
        auto x168 = 2. * x167; auto x169 = x168 * x26; auto x170 = AXG * x169;
        auto x171 = x166 + x170; auto x172 = x0 * x132; auto x173 = x26 * x82;
        auto x174 = x172 * x173; auto x175 = x132 * x6; auto x176 = x173 * x175;
        auto x177 = x132 * x70; auto x178 = AXG * x177; auto x179 = 4. * x178;
        auto x180 = x132 * x77; auto x181 = AXG * x180; auto x182 = 4. * x181;
        auto x183 = MXs * x53; auto x184 = MXs * x96; auto x185 = x66 * x70; auto x186 = AXG * x185;
        auto x187 = 8. * x186; auto x188 = x66 * x77; auto x189 = AXG * x188; auto x190 = 8. * x189;
        auto x191 = MXs * x52; auto x192 = x191 * x40; auto x193 = 4. * MXs;
        auto x194 = x103 * x193;
        auto x195 = x174 + x176 - x179 - x182 - x183 - x184 - x187 - x190 + x192 + x194;
        auto x196 = 4. * s; auto x197 = x12 * x127; auto x198 = Mqis * x197;
        auto x199 = Mqis * x128; auto x200 = Mqis * x130; auto x201 = Mqis * x135;
        auto x202 = Mqis * x137; auto x203 = Mqis * x140; auto x204 = MG * Mqi;
        auto x205 = MUs * x204; auto x206 = LG * LGp; auto x207 = x205 * x206; auto x208 = 6. * AXG;
        auto x209 = x2 * x208; auto x210 = RG * RGp; auto x211 = x205 * x210;
        auto x212 = AXG * x200; auto x213 = Mqis * x154; auto x214 = AXG * x202;
        auto x215 = Mqis * x158; auto x216 = AXG * x197; auto x217 = Mqis * x216;
        auto x218 = AXG * x128; auto x219 = Mqis * x218; auto x220 = 8. * AXG;
        auto x221 = x18 * x220; auto x222 = UU * x148; auto x223 = x222 * x32; auto x224 = s * x223;
        auto x225 = x122 * x40; auto x226 = UU * x122; auto x227 = MUs * x226;
        auto x228 = MXs * x227; auto x229 = s * x228; auto x230 = MUs * x222;
        auto x231 = MXs * x230; auto x232 = s * x231; auto x233 = pow(u, 4); auto x234 = 4. * x64;
        auto x235 = AXG * x234; auto x236 = MUs * x15; auto x237 = x25 * x50;
        auto x238 = AXG * x233; auto x239 = x30 * x50;
        auto x240 = -x233 * x235 - x233 * x236 * x40 + x237 * x238 + x238 * x239;
        auto x241 = AXG * x224 + AXG * x229 + AXG * x232 + x2 * x225 - x224 - x229 - x232 + x240;
        auto x242 = x122 * x32; auto x243 = x242 * x60; auto x244 = x32 * x60;
        auto x245 = x148 * x244; auto x246 = x11 * x121; auto x247 = x0 * x246;
        auto x248 = 2. * x60; auto x249 = x246 * x6; auto x250 = AXG * x247; auto x251 = AXG * x249;
        auto x252 = x27 * x46; auto x253 = u * x134; auto x254 = AXG * x253; auto x255 = 2. * x11;
        auto x256 = x254 * x255; auto x257 = u * x139; auto x258 = AXG * x257;
        auto x259 = x255 * x258; auto x260 = pow(MX, 6); auto x261 = AXG * x260;
        auto x262 = 2. * x261; auto x263 = x12 * x262; auto x264 = x263 * x46;
        auto x265 = x15 * x262; auto x266 = x265 * x46; auto x267 = x13 * x46; auto x268 = 2. * x74;
        auto x269 = x267 * x268; auto x270 = x16 * x46; auto x271 = x268 * x270;
        auto x272 = x33 * x50; auto x273 = x26 * x272; auto x274 = x172 * x27;
        auto x275 = x26 * x274; auto x276 = x36 * x51; auto x277 = x175 * x27;
        auto x278 = x26 * x277; auto x279 = 2. * AXG; auto x280 = x127 * x279;
        auto x281 = MXs * x179; auto x282 = x110 * x32; auto x283 = MXs * x182;
        auto x284 = -x103 * x11 + x11 * x96 - x273 - x275 - x276 - x278 - x280 * x77 + x281 + x282 +
                    x283 + x32 * x73;
        auto x285 = pow(MU, 6); auto x286 = SS * x285; auto x287 = x286 * x46;
        auto x288 = x285 * x60; auto x289 = x0 * x288; auto x290 = MXs * x289;
        auto x291 = x288 * x6; auto x292 = MXs * x291; auto x293 = SS * x160; auto x294 = UU * x242;
        auto x295 = AXG * x294; auto x296 = UU * x247; auto x297 = AXG * x296;
        auto x298 = 2. * x297; auto x299 = UU * x249; auto x300 = AXG * x299; auto x301 = 2. * x300;
        auto x302 = 4. * syFC10; auto x303 = x172 * x60; auto x304 = x0 * x244;
        auto x305 = x175 * x60; auto x306 = x244 * x6; auto x307 = x36 * x50; auto x308 = x11 * x60;
        auto x309 = x25 * x308; auto x310 = x30 * x308; auto x311 = AXG * x308;
        auto x312 = x0 * x311; auto x313 = x311 * x6; auto x314 = UU * x33; auto x315 = UU * x36;
        auto x316 = 2. * x216; auto x317 = 2. * x13; auto x318 = AXG * u; auto x319 = 2. * x218;
        auto x320 = AXG * x108; auto x321 = 2. * x16; auto x322 = x15 * x46; auto x323 = MXs * x40;
        auto x324 = x322 * x323; auto x325 = x252 * x30; auto x326 = AXG * x325;
        auto x327 = x11 * x50; auto x328 = x25 * x327; auto x329 = x30 * x327;
        auto x330 = AXG * x327; auto x331 = x0 * x330; auto x332 = x330 * x6;
        auto x333 = x328 + x329 - x331 - x332; auto x334 = x12 * x46; auto x335 = x25 * x252;
        auto x336 = AXG * x335 - x323 * x334;
        auto x337 = -2. * x103 - 2. * x163 * x314 - 2. * x163 * x315 - 2. * x263 - 2. * x265 -
                    2. * x272 - 2. * x274 - 2. * x277 + 2. * x303 + 2. * x304 + 2. * x305 +
                    2. * x306 - 2. * x307 - 2. * x309 - 2. * x310 + 2. * x312 + 2. * x313 +
                    2. * x316 + 2. * x317 * x318 + 2. * x317 * x74 + 2. * x318 * x321 + 2. * x319 +
                    4. * x320 + 2. * x321 * x74 - 2. * x324 + 2. * x326 + 2. * x333 + 2. * x336 -
                    2. * x52 - 2. * x95;
        auto x338 = MXs * x303; auto x339 = MXs * x305; auto x340 = x153 * x38;
        auto x341 = UU * x46; auto x342 = x341 * x36; auto x343 = AXG * x342;
        auto x344 = x157 * x38; auto x345 = Mqis * x160; auto x346 = AXG * x345;
        auto x347 = x255 * x74; auto x348 = MXs * Mqis; auto x349 = x334 * x348;
        auto x350 = x322 * x348; auto x351 = Mqis * x27; auto x352 = x351 * x46;
        auto x353 = AXG * x352; auto x354 = 2. * syFC3; auto x355 = SS * x46;
        auto x356 = x33 * x355; auto x357 = x355 * x36; auto x358 = SS * x250;
        auto x359 = SS * x251; auto x360 = x11 * x6; auto x361 = x11 * x153; auto x362 = x11 * x157;
        auto x363 = AXG * x355; auto x364 = x33 * x341; auto x365 = AXG * x267;
        auto x366 = AXG * x270; auto x367 = SS * x242; auto x368 = x167 * x32;
        auto x369 = x122 * x56; auto x370 = x148 * x56; auto x371 = -x367 - x368 - x369 - x370;
        auto x372 = -MUs * x265 + MXs * x316 + MXs * x319 + x0 * x11 * x363 - x11 * x84 -
                    x154 * x160 - x158 * x160 + x167 * x46 + x247 * x82 + x249 * x82 - x262 * x64 +
                    x338 + x339 - x342 * x40 + x356 + x357 - x358 - x359 + x360 * x363 -
                    x360 * x86 + 2. * x361 + 2. * x362 - x364 * x40 + 2. * x365 + 2. * x366 + x371;
        auto x373 = x286 * x33; auto x374 = x286 * x36; auto x375 = 2. * x286;
        auto x376 = x375 * x6; auto x377 = x11 * x376; auto x378 = x285 * x50;
        auto x379 = x0 * x378; auto x380 = MXs * x379; auto x381 = x378 * x6;
        auto x382 = MXs * x381; auto x383 = 8. * Mqis; auto x384 = SS * x383;
        auto x385 = AXG * x226; auto x386 = 2. * x260; auto x387 = x385 * x386;
        auto x388 = AXG * x222; auto x389 = x386 * x388; auto x390 = 4. * Mqis;
        auto x391 = Mqis * x11; auto x392 = 8. * x391; auto x393 = MUs * x298;
        auto x394 = MXs * x298; auto x395 = MUs * x301; auto x396 = MXs * x301;
        auto x397 = x261 * x64; auto x398 = x236 * x261; auto x399 = MXs * x216;
        auto x400 = MXs * x218; auto x401 = 2. * syFC7; auto x402 = s * x296; auto x403 = s * x299;
        auto x404 = s * x294; auto x405 = x134 * x38; auto x406 = x139 * x38;
        auto x407 = 2. * syFC10; auto x408 = Mqis * x274; auto x409 = Mqis * x277;
        auto x410 = Mqis * x328; auto x411 = Mqis * x329; auto x412 = 4. * x197;
        auto x413 = 4. * x128; auto x414 = x105 * x216; auto x415 = x105 * x74;
        auto x416 = x105 * x218; auto x417 = -x119 * x247 - x119 * x249 + x250 * x82 + x251 * x82;
        auto x418 = Mqis * x52; auto x419 = Mqis * x95; auto x420 = Mqis * x28;
        auto x421 = Mqis * x70; auto x422 = Mqis * x77; auto x423 = MUs * x421;
        auto x424 = AXG * x423; auto x425 = MUs * x422; auto x426 = AXG * x425;
        auto x427 = x255 * x318; auto x428 = AXG * Mqis; auto x429 = AXG * x418 + Mqis * x103;
        auto x430 = x204 * x206; auto x431 = x121 * x430; auto x432 = x119 * x32;
        auto x433 = x204 * x210; auto x434 = x121 * x433; auto x435 = 8. * x430;
        auto x436 = SS * x246; auto x437 = 8. * x433; auto x438 = 6. * Mqis;
        auto x439 = x129 * x261; auto x440 = AXG * x430; auto x441 = x119 * x246;
        auto x442 = AXG * x433; auto x443 = 4. * x56; auto x444 = AXG * x435;
        auto x445 = x11 * x444; auto x446 = x133 * x433; auto x447 = x220 * x446;
        auto x448 = UU * x127; auto x449 = x430 * x448; auto x450 = x433 * x448;
        auto x451 = x0 * x375; auto x452 = AXG * x286; auto x453 = x0 * x452; auto x454 = x452 * x6;
        auto x455 = -x11 * x451 + x11 * x453 + x11 * x454; auto x456 = Mqis * x120;
        auto x457 = x132 * x220; auto x458 = MXs * x103; auto x459 = x121 * x27;
        auto x460 = x430 * x459; auto x461 = x433 * x459; auto x462 = 4. * x355;
        auto x463 = x208 * x66; auto x464 = x220 * x61; auto x465 = AXG * x120;
        auto x466 = x430 * x51; auto x467 = x433 * x51; auto x468 = x463 * x69;
        auto x469 = -x121 * x24; auto x470 = -2. * x123 * x18 + x469; auto x471 = x132 * x19;
        auto x472 = x132 * x22; auto x473 = MXs * x113; auto x474 = MXs * x116;
        auto x475 = 2. * x132; auto x476 = 2. * x473; auto x477 = 2. * x474;
        auto x478 = -AXG * x476 - AXG * x477 + MXs * x150 + MXs * x151 + 2. * x121 * x8 +
                    x4 * x475 - x471 - x472 + x473 + x474 + x475 * x8;
        auto x479 = 2. * x27; auto x480 = x122 * x479; auto x481 = x148 * x479;
        auto x482 = AXG * x223; auto x483 = MXs * x385; auto x484 = 4. * x483;
        auto x485 = MXs * x388; auto x486 = 4. * x485; auto x487 = x11 * x52;
        auto x488 = x121 * x70; auto x489 = AXG * x488; auto x490 = x121 * x77;
        auto x491 = AXG * x490; auto x492 = 2. * MUs; auto x493 = 2. * x28; auto x494 = 2. * x51;
        auto x495 = x122 * x41 - x122 * x493 + x122 * x494 + x148 * x41 - x148 * x493 + x148 * x494;
        auto x496 = x127 * x40; auto x497 = 4. * syFC11; auto x498 = -x192 - x194;
        auto x499 = x162 + x164; auto x500 = -x126 - x149; auto x501 = 4. * syFC8;
        auto x502 = x12 * x38; auto x503 = 3. * AXG; auto x504 = x503 * x70; auto x505 = x503 * x77;
        auto x506 = AXG * x173; auto x507 = x0 * x506; auto x508 = x506 * x6; auto x509 = x25 * x61;
        auto x510 = AXG * x509; auto x511 = s * (x39 * x503 + x502 * x503 + x502 - x504 - x505 +
                                                 x507 + x508 + x509 - x510 - x70 - x77);
        auto x512 = x122 * x248; auto x513 = x148 * x248; auto x514 = x38 * x92;
        auto x515 = x38 * x94; auto x516 = 3. * x46; auto x517 = MGs * x62; auto x518 = MGs * SS;
        auto x519 = AXG * x518; auto x520 = x25 * x519; auto x521 = x30 * x519;
        auto x522 = x503 * x92; auto x523 = x503 * x94; auto x524 = s * syFC3;
        auto x525 = x205 * x69; auto x526 = x206 * x525; auto x527 = x210 * x525;
        auto x528 = 8. * x26; auto x529 = SS * x122; auto x530 = x167 * x40;
        auto x531 = x100 * x503 - 7. * x100 + 6. * x113 + 6. * x116 + x124 + x26 * x530 +
                    3. * x488 + x489 + 3. * x490 + x491 + x503 * x97 - x528 * x529 - 7. * x97;
        auto x532 = s * syFC7; auto x533 = 3. * Mqis; auto x534 = UU * x533; auto x535 = 4. * x130;
        auto x536 = 4. * x137; auto x537 = AXG * x535; auto x538 = Mqis * x503;
        auto x539 = 3. * x66; auto x540 = AXG * x480; auto x541 = AXG * x481; auto x542 = 2. * MGs;
        auto x543 = AXG * x412; auto x544 = AXG * x413;
        auto x545 = -x480 - x481 - x542 * x89 + x543 + x544; auto x546 = 3. * x32;
        auto x547 = 3. * x11; auto x548 = Mqis * x502; auto x549 = Mqis * x39;
        auto x550 = Mqis * x62; auto x551 = -x542 * x88; auto x552 = 3. * MXs;
        auto x553 = x334 * x74; auto x554 = x322 * x74;
        auto x555 = -3. * x108 - x197 + x216 + x218 - x253 + x254 - x257 + x258 - x320 +
                    x322 * x552 + x334 * x552 - x512 - x513 + x53 + x553 + x554 - 3. * x71 - x72 +
                    x96;
        auto x556 = x167 * x528; auto x557 = x504 * x66; auto x558 = x133 * x430;
        auto x559 = x121 * x60; auto x560 = x430 * x559; auto x561 = x433 * x559;
        auto x562 = MXs * x341; auto x563 = x430 * x562; auto x564 = x433 * x562;
        auto x565 = x323 * x341; auto x566 = MXs * x95;
        auto x567 = x177 * x503 - 3. * x177 + x180 * x503 - 3. * x180 - 7. * x185 - 7. * x188 +
                    4. * x191 + x498 + x505 * x66 + 4. * x566;
        auto x568 = x252 * x430; auto x569 = x252 * x433; auto x570 = x341 * x74;
        auto x571 = syFC5 * x196; auto x572 = Mqis * x60; auto x573 = x122 * x572;
        auto x574 = x148 * x572; auto x575 = Mqis * x253; auto x576 = Mqis * x257;
        auto x577 = x0 * x352; auto x578 = x352 * x6; auto x579 = u * x133; auto x580 = x430 * x579;
        auto x581 = x433 * x579; auto x582 = x318 * x558; auto x583 = x318 * x446;
        auto x584 = x431 * x479; auto x585 = x434 * x479; auto x586 = x129 * x32;
        auto x587 = x430 * x586; auto x588 = MXs * x133; auto x589 = x430 * x588;
        auto x590 = x433 * x586; auto x591 = x433 * x588; auto x592 = x558 * x74;
        auto x593 = x446 * x74; auto x594 = x27 * x285; auto x595 = x0 * x594;
        auto x596 = x594 * x6; auto x597 = AXG * x595; auto x598 = AXG * x596;
        auto x599 = x296 - x297 + x595 + x596 - x597 - x598; auto x600 = syFC6 * x196;
        auto x601 = 2. * x200; auto x602 = 2. * x202; auto x603 = -Mqis * x480 - Mqis * x481;
        auto x604 = AXG * x228; auto x605 = AXG * x231;
        auto x606 = -x223 - x228 - x231 - x294 + x295 + x299 - x300 + x482 + x604 + x605;
        auto x607 = -AXG * x364 - x267 - x270 - x340 + x342 - x343 + x364 + x365 + x366 + x405;
        auto x608 = syFC9 * x196; auto x609 = AXG * x289; auto x610 = AXG * x291;
        auto x611 = x226 * x46; auto x612 = x226 * x38; auto x613 = x222 * x46;
        auto x614 = x222 * x38; auto x615 = x385 * x46; auto x616 = x38 * x385;
        auto x617 = x388 * x46; auto x618 = x38 * x388; auto x619 = x123 * x125;
        auto x620 = AXG * x149; auto x621 = -x177 + x178 - x180 + x181 + x186 + x189;
        auto x622 = x616 + x618; auto x623 = s * x407; auto x624 = 4. * x327; auto x625 = x0 * x624;
        auto x626 = x6 * x624; auto x627 = x391 * x6; auto x628 = 2. * x50; auto x629 = x122 * x628;
        auto x630 = x148 * x628; auto x631 = AXG * SS; auto x632 = AXG * x130;
        auto x633 = AXG * x105; auto x634 = x32 * x40; auto x635 = x40 * x66; auto x636 = 4. * x27;
        auto x637 = x122 * x636; auto x638 = AXG * x637; auto x639 = x148 * x636;
        auto x640 = AXG * x639; auto x641 = x638 + x640; auto x642 = 4. * x32;
        auto x643 = x50 * x642; auto x644 = x132 * x636; auto x645 = 4. * x330;
        auto x646 = AXG * x587; auto x647 = AXG * x590; auto x648 = AXG * x227;
        auto x649 = AXG * x230; auto x650 = x119 * x132; auto x651 = -x165; auto x652 = x122 * x351;
        auto x653 = x148 * x351; auto x654 = x248 * x431; auto x655 = x248 * x434;
        auto x656 = x46 * x479; auto x657 = x268 * x341; auto x658 = x172 * x345;
        auto x659 = x175 * x345; auto x660 = 4. * x127; auto x661 = x132 * x193;
        auto x662 = x348 * x46; auto x663 = 2. * x662; auto x664 = x132 * x323;
        auto x665 = x40 * x662; auto x666 = -x227 - x230 + x648 + x649; auto x667 = 2. * x69;
        auto x668 = x430 * x667; auto x669 = x433 * x667; auto x670 = AXG * x146;
        auto x671 = AXG * x168; auto x672 = Mqis * x237; auto x673 = x105 * x64;
        auto x674 = x105 * x236;
        ,
        TR * TR * gs * gs * (Nc - 1) * (Nc + 1) * (L * Rp + Lp * R) *
            (AXG * x526 * x571 -
             MUs * x600 *
                 (-AXG * x668 - AXG * x669 - x173 * x430 - x173 * x433 - x19 + x20 - x22 + x23 +
                  x3 - x4 + x430 * x506 + x433 * x506 + x668 + x669 + x7 - x8) +
             UU * x532 *
                 (AXG * x242 + AXG * x658 + AXG * x659 + x0 * x663 + x0 * x665 - x207 * x634 +
                  x207 * x642 - x211 * x634 + x211 * x642 - x242 + x247 + x249 - x250 - x251 +
                  x430 * x496 - x430 * x660 + x430 * x661 - x430 * x664 + x433 * x496 -
                  x433 * x660 + x433 * x661 - x433 * x664 + x6 * x663 + x6 * x665 - x658 - x659) -
             s * syFC4 * x555 +
             s * x497 *
                 (-3. * x100 - x101 + x142 + x144 - x169 + 3. * x185 + 3. * x188 + x607 + x621 +
                  x651 - 3. * x97 - x98) +
             s * x501 * (-x128 + x555) - syFC1 * x337 - syFC1 * x511 - syFC1 * x90 -
             8. * syFC11 *
                 (x122 * x252 + x148 * x252 + x243 + x245 - x247 * x248 - x248 * x249 + x250 * x60 +
                  x251 * x60 - x256 - x259 + x264 + x266 - x269 - x271 + x284) +
             syFC2 * x337 + syFC2 * x511 + syFC2 * x90 +
             syFC3 * (-MGs * x104 + MGs * x110 + MGs * x21 + MGs * x24 + MGs * x29 + MGs * x31 -
                      MGs * x42 - MGs * x43 - MGs * x53 - MGs * x54 + MGs * x73 + MGs * x76 +
                      MGs * x78 - MGs * x96 - Mqis * x21 - Mqis * x24 + Mqis * x5 - Mqis * x68 +
                      Mqis * x81 + Mqis * x9 + s * x106 + s * x107 - s * x111 - s * x112 + x102 +
                      x124 - x92 * x93 - x93 * x94 + x99) -
             2. * syFC4 * x372 +
             syFC4 * (s * x128 + x126 - x131 - x136 - x138 - x141 - x143 - x145 - x147 + x149 +
                      x150 + x151 + x152 + x155 + x156 + x159 - x162 - x164 + x171 + x195) +
             syFC7 * x160 *
                 (AXG * x410 + AXG * x411 + Mqis * x543 + Mqis * x544 - Mqis * x625 - Mqis * x626 +
                  Mqis * x629 + Mqis * x630 - 8. * x212 + x361 + x362 + x371 - x397 - x398 + x399 +
                  x40 * x652 + x40 * x653 + x400 + 2. * x482 + x603 + 2. * x604 + 2. * x605) +
             syFC7 * (-x196 * x198 - x196 * x199 + x196 * x200 + x196 * x201 + x196 * x202 +
                      x196 * x203 - x196 * x212 - x196 * x213 - x196 * x214 - x196 * x215 +
                      x196 * x217 + x196 * x219 - x207 * x209 + x207 * x221 - x209 * x211 +
                      x211 * x221 + x241) +
             8. * syFC8 * x372 -
             u * x354 *
                 (AXG * x328 + AXG * x329 + x0 * x391 * x631 + x13 * x633 + x16 * x633 -
                  x172 * x518 - x175 * x518 + 2. * x272 - x33 * x518 - x36 * x518 - x391 * x83 +
                  x545 + x551 - x625 - x626 + x627 * x631 - x627 * x82 + x629 + x630 - 8. * x632 +
                  x634 * x92 + x634 * x94 + x635 * x92 + x635 * x94 + x641) -
             u * x401 *
                 (-x220 * x450 - 4. * x295 + x298 + x301 + x327 * x435 + x327 * x437 + x379 + x381 +
                  x417 - x430 * x643 - x430 * x644 - x430 * x645 - x433 * x643 - x433 * x644 -
                  x433 * x645 + 6. * x592 + 6. * x593 - x595 - x596 + 2. * x597 + 2. * x598 +
                  6. * x646 + 6. * x647) -
             u * x532 *
                 (x128 * x503 - 3. * x128 - x130 * x503 + 3. * x130 + 3. * x135 - x137 * x503 +
                  3. * x137 + 3. * x140 - x153 * x552 - x157 * x552 + x167 * x383 + x197 * x503 -
                  3. * x197 + x222 * x552 + x226 * x552 - x376 + x383 * x529 - x40 * x446 +
                  4. * x446 - x451 + x483 + x485 - x637 - x639 + x641 + x666) +
             u * x600 *
                 (AXG * x137 - Mqis * x146 - Mqis * x168 + Mqis * x670 + Mqis * x671 + x105 * x65 +
                  x105 * x67 - x105 * x79 - x105 * x80 - x106 - x107 + x111 + x112 + x128 - x130 -
                  x135 - x137 - x140 + x154 + x158 + x197 - x216 - x218 + x481 + x632 + x666) -
             x26 * x401 *
                 (Mqis * x119 * x123 + Mqis * x530 - x316 - x319 + x333 + x430 * x650 +
                  x433 * x650 - x447 - x453 - x454 + x480 + x481 + x484 + x486 + x537 - x629 -
                  x630 - x638 - x640 + 2. * x648 + 2. * x649) +
             x26 * x524 * (x234 + x235 + x520 + x521 - x522 - x523 - x92 - x94) -
             x26 * x600 *
                 (AXG * x672 - AXG * x673 - AXG * x674 - x134 - x139 + x146 + x153 + x157 + x168 -
                  x222 - x226 + x239 * x428 + x385 + x388 + x57 + x59 - x65 - x67 - x670 - x671 -
                  x672 + x673 + x674 + x79 + x80 - x88 - x89) +
             x302 * (AXG * x487 - x193 * x489 - x193 * x491 + x26 * x453 + x26 * x454 + x280 * x70 +
                     x469 + x478 - x489 * x492 - x491 * x492 + x495) -
             x302 * (u * x298 + u * x301 + x0 * x287 + x11 * x53 - x160 * x295 + 2. * x243 -
                     x247 * x293 + x248 * x250 + x248 * x251 - x249 * x293 + x279 * x290 +
                     x279 * x292 + x284 + x287 * x6 - x290 - x292) +
             x302 * (x160 * x482 - 2. * x245 + x256 + x259 - x264 - x266 + x269 + x271 + x373 +
                     x374 - x377 + x380 + x382 + x387 + x389 - x393 - x394 - x395 - x396 + x455 -
                     x46 * x480 - x46 * x481 + x46 * x484 + x46 * x486) +
             x354 * (-MGs * x309 - MGs * x310 + MGs * x312 + MGs * x313 - x0 * x420 + x170 + x195 +
                     x29 * x428 + x31 * x428 - x323 * x421 - x323 * x422 + x418 + x419 - x420 * x6 -
                     2. * x424 - 2. * x426 + x427 * x92 + x427 * x94 + x429) +
             x354 * (Mqis * x263 + Mqis * x265 + Mqis * x272 + Mqis * x307 + Mqis * x331 +
                     Mqis * x332 - x13 * x415 - x16 * x415 + x234 * x261 - 4. * x361 - 4. * x362 +
                     2. * x367 + 2. * x368 + 2. * x369 + 2. * x370 + 4. * x398 + x408 + x409 -
                     x410 - x411 - x412 * x74 - x413 * x74 - x414 - x416 + x417) -
             x354 *
                 (MGs * x272 + MGs * x274 + MGs * x277 + MGs * x307 - MGs * x328 - MGs * x329 +
                  MGs * x331 + MGs * x332 + Mqis * x303 + Mqis * x304 + Mqis * x305 + Mqis * x306 +
                  x25 * x353 + x262 * x92 + x262 * x94 - x280 * x92 - x280 * x94 + x30 * x353 -
                  x314 * x346 - x315 * x346 + 2. * x338 + 2. * x339 - 8. * x340 - 8. * x343 -
                  8. * x344 - x347 * x92 - x347 * x94 - x349 * x40 - x350 * x40 + x36 * x86) +
             x401 * (Mqis * x115 + Mqis * x118 - x101 * x383 + x132 * x444 * x69 + x193 * x466 +
                     x193 * x467 - x383 * x98 + x430 * x468 - x431 * x465 + x433 * x468 -
                     x434 * x465 - x466 * x75 - x467 * x75 + x470 + x478) +
             x401 *
                 (x11 * x447 + x133 * x445 + x160 * x408 + x160 * x409 - x340 * x438 - x343 * x438 -
                  x344 * x438 + x356 * x390 + x357 * x390 - x431 * x432 - x431 * x443 -
                  x432 * x434 - x434 * x443 + x435 * x436 - x435 * x439 + x436 * x437 -
                  x437 * x439 - x440 * x441 - x441 * x442 + x449 * x75 + x450 * x75 + x455) +
             x401 * (-x160 * x460 - x160 * x461 - x172 * x456 - x175 * x456 + x193 * x418 +
                     x193 * x419 + x273 + x275 + x276 + x278 - x281 - x282 - x283 + x341 * x445 -
                     x383 * x458 - x418 * x75 + x421 * x457 + x421 * x463 + x422 * x457 +
                     x422 * x463 + x431 * x462 + x431 * x464 + x434 * x462 + x434 * x464) -
             x401 * (-x153 * x392 - x157 * x392 - x247 * x384 - x249 * x384 + x358 * x390 +
                     x359 * x390 + x367 * x390 + x368 * x390 + x369 * x390 + x370 * x390 - x373 -
                     x374 + x377 - x380 - x382 + x383 * x397 + x383 * x398 - x383 * x399 -
                     x383 * x400 - x387 - x389 + x393 + x394 + x395 + x396) +
             x407 * (-AXG * x402 - AXG * x403 + AXG * x404 + x10 * x340 + x10 * x344 - x10 * x405 -
                     x10 * x406 + x17 * x46 - x18 * x225 + x241 - x35 * x46 - x37 * x46 + x402 +
                     x403 - x404 + x44 * x46 + x45 * x46) +
             x497 * (MXs * x102 + MXs * x99 - s * x344 + s * x406 + x11 * x54 + x132 * x5 +
                     x132 * x9 + x240 - x40 * x473 - x40 * x474 + x470 - 2. * x471 - 2. * x472 +
                     x476 + x477 - 4. * x487 + x495 + x496 * x70) +
             x501 * (x131 + x136 + x138 + x141 + x143 + x145 + x147 - x150 - x151 - x152 - x155 -
                     x156 - x159 - x166 - x170 - x174 - x176 + x179 + x182 + x183 + x184 + x187 +
                     x190 + x498 + x499 + x500) -
             x524 * (-AXG * x536 - x13 * x533 + x13 * x538 + 4. * x135 + 4. * x140 - 4. * x154 -
                     4. * x158 - x16 * x533 + x16 * x538 - x314 * x538 - x315 * x538 + x33 * x534 +
                     x36 * x534 - x412 - x413 + x523 * x66 + x533 * x65 + x535 + x536 - x537 -
                     x538 * x65 - x539 * x94 + x540 + x541 + x542 * x57 + x542 * x59 + x545) +
             x524 * (AXG * x512 + AXG * x513 - AXG * x517 + MGs * x509 - MGs * x510 - Mqis * x507 -
                     Mqis * x508 - x104 + x109 + x110 - x193 * x334 - x335 + x336 + x421 * x503 +
                     x421 + x422 * x503 + x422 - x46 * x520 - x46 * x521 + x46 * x522 + x46 * x523 +
                     x503 * x514 + x503 * x515 + x512 + x513 + x514 + x515 - x516 * x92 -
                     x516 * x94 + x517 - x53 + x55 - x96) -
             x524 *
                 (-AXG * x550 + Mqis * x509 - Mqis * x510 - x11 * x522 - x11 * x523 - x134 * x160 -
                  x139 * x160 + x153 * x160 + x157 * x160 + x193 * x322 + x32 * x522 + x32 * x523 -
                  x322 * x533 + x322 * x538 + x324 + x325 - x326 - x334 * x533 + x334 * x538 -
                  x428 * x84 - x428 * x87 + x503 * x548 + x503 * x549 + x522 * x66 - x539 * x92 -
                  x546 * x92 - x546 * x94 + x547 * x92 + x547 * x94 + x548 + x549 + x550 + x551) +
             x532 *
                 (-x40 * x526 - x40 * x527 + 8. * x466 + 8. * x467 - 2. * x526 - 2. * x527 + x531) -
             x532 * (x160 * x558 - x163 * x558 + x40 * x423 + x40 * x425 - 8. * x418 - 8. * x419 +
                     2. * x423 + 2. * x425 - x430 * x565 - x433 * x565 + x556 + x557 + 8. * x560 +
                     8. * x561 - 2. * x563 - 2. * x564 + x567) +
             x571 * (-AXG * x652 - AXG * x653 + x198 + x199 - x200 - x201 - x202 - x203 + x212 +
                     x213 + x214 + x215 - x217 - x219 + x450 - x587 - x589 - x590 - x591 + x592 +
                     x593 + x646 + x647 + x652 + x653) -
             x571 * (AXG * x466 + AXG * x467 - AXG * x527 - AXG * x560 - AXG * x561 - AXG * x568 -
                     AXG * x569 - x418 - x419 + x423 - x424 + x425 - x426 + x429 + x430 * x570 +
                     x433 * x570 - x466 - x467 + x526 + x527 + x560 + x561 - x563 - x564 + x568 +
                     x569) -
             x571 * (AXG * x449 + AXG * x450 + AXG * x460 + AXG * x461 - AXG * x573 - AXG * x574 -
                     AXG * x577 - AXG * x578 + Mqis * x254 + Mqis * x258 + Mqis * x553 +
                     Mqis * x554 - x349 - x350 - x449 - x460 - x461 + x573 + x574 - x575 - x576 +
                     x577 + x578 - x580 - x581 + x582 + x583) -
             x600 * (-AXG * x601 - AXG * x602 + Mqis * x540 + Mqis * x541 - x105 * x154 -
                     x105 * x158 - 2. * x198 - 2. * x199 + 2. * x201 + 2. * x203 + x414 + x416 +
                     x601 + x602 + x603 + x606) -
             x600 * (AXG * x584 + AXG * x585 + x105 * x254 + x105 * x258 + x279 * x449 +
                     x279 * x450 - x279 * x587 - x279 * x590 - 2. * x449 - 2. * x450 - 2. * x575 -
                     2. * x576 - x584 - x585 + 2. * x587 + 2. * x589 + 2. * x590 + 2. * x591 -
                     2. * x592 - 2. * x593 + x599) +
             x600 * (AXG * x654 + AXG * x655 + x126 + x289 + x291 + 2. * x419 - x430 * x656 -
                     x430 * x657 - x433 * x656 - x433 * x657 + x440 * x656 + x442 * x656 +
                     2. * x563 + 2. * x564 + 2. * x580 + 2. * x581 - 2. * x582 - 2. * x583 - x609 -
                     x610 - x612 - x614 - x619 - x620 + x622 - x654 - x655) -
             x608 * (-x344 + x406 + x599 + x606 + x607) +
             x608 * (-x100 + x101 + x113 - x114 + x116 - x117 + x171 + x488 - x489 + x490 - x491 +
                     x651 - x97 + x98) -
             x608 * (-AXG * x191 + x169 - x185 - x188 + x191 - x289 - x291 - x458 + x500 + x566 +
                     x609 + x610 + x611 + x612 + x613 + x614 - x615 - x616 - x617 - x618 + x619 +
                     x620 + x621) +
             x623 * (x531 - x556 - x557) -
             x623 * (-x122 * x161 - x148 * x161 + x267 * x503 - 3. * x267 + x270 * x503 -
                     2. * x289 - 2. * x291 + x499 + x567 - x611 + 3. * x612 - x613 + 3. * x614 +
                     x615 + x617 + x622)) *
            Denom(1536 * Nc * s * (-MQis * MUs + MQis * u - x26 + x46) * pow(MU, 2) * pow(Pi, 2)))

        ;
    //}
  }
  return ret.real();
}
