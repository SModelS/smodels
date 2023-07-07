#include "pxs_gausq_3.h"

#define C0(a, b, c, d, e, f) C0i(cc0, a, b, c, d, e, f)
#define C1(a, b, c, d, e, f) C0i(cc1, a, b, c, d, e, f)
#define C2(a, b, c, d, e, f) C0i(cc2, a, b, c, d, e, f)
#define C00(a, b, c, d, e, f) C0i(cc00, a, b, c, d, e, f)
#define C11(a, b, c, d, e, f) C0i(cc11, a, b, c, d, e, f)
#define C12(a, b, c, d, e, f) C0i(cc12, a, b, c, d, e, f)
#define C22(a, b, c, d, e, f) C0i(cc22, a, b, c, d, e, f)

#define D0(a, b, c, d, e, f, g, h, i, j) D0i(dd0, a, b, c, d, e, f, g, h, i, j)
#define D1(a, b, c, d, e, f, g, h, i, j) D0i(dd1, a, b, c, d, e, f, g, h, i, j)
#define D2(a, b, c, d, e, f, g, h, i, j) D0i(dd2, a, b, c, d, e, f, g, h, i, j)
#define D3(a, b, c, d, e, f, g, h, i, j) D0i(dd3, a, b, c, d, e, f, g, h, i, j)
#define D00(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd00, a, b, c, d, e, f, g, h, i, j)
#define D11(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd11, a, b, c, d, e, f, g, h, i, j)
#define D12(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd12, a, b, c, d, e, f, g, h, i, j)
#define D22(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd22, a, b, c, d, e, f, g, h, i, j)
#define D13(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd13, a, b, c, d, e, f, g, h, i, j)
#define D23(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd23, a, b, c, d, e, f, g, h, i, j)
#define D33(a, b, c, d, e, f, g, h, i, j)                                      \
  D0i(dd33, a, b, c, d, e, f, g, h, i, j)

ComplexType ME_us_box_gqqQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2,
                           double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto Denom = [](auto a) { return 1. / a; };
  //*
#define syFC1 C0(0, MUs + MXs - s - u, MUs, 0, 0, MQs)
#define syFC2 C0(MXs, 0, MUs + MXs - s - u, MQs, 0, 0)
#define syFC3 C0(MXs, s, MUs, MQs, 0, 0)
#define syFC4 D0(MXs, 0, 0, MUs, MUs + MXs - s - u, s, MQs, 0, 0, 0)
#define syFC5 C1(MXs, s, MUs, MQs, 0, 0)
#define syFC6 C1(MXs, MUs + MXs - s - u, 0, 0, MQs, 0)
#define syFC7 D1(0, s, MUs, MUs + MXs - s - u, 0, MXs, 0, 0, 0, MQs)
#define syFC8 D1(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC9 C2(MXs, s, MUs, MQs, 0, 0)
#define syFC10 C2(MXs, MUs + MXs - s - u, 0, 0, MQs, 0)
#define syFC11 D2(0, s, MUs, MUs + MXs - s - u, 0, MXs, 0, 0, 0, MQs)
#define syFC12 D2(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC13 D3(0, s, MUs, MUs + MXs - s - u, 0, MXs, 0, 0, 0, MQs)
#define syFC14 D3(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC15 D00(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC16 D11(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC17 D12(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC18 D13(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC19 D22(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC20 D23(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
#define syFC21 D33(MXs, MUs + MXs - s - u, 0, s, 0, MUs, 0, MQs, 0, 0)
  _EPS0(
      ret,
      (-1 + Nc) * (+1 + Nc) *
          (+Denom(768 * pow(Pi, 2) * s * MUs * Nc -
                  768 * pow(Pi, 2) * s * u * Nc)) *
          (+TR) * (+TR) * (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 *
               (4 * s * pow(MUs, 2) * SS - 4 * s * MXs * MUs * SS +
                4 * s * MXs * MUs * SS * AXG - 4 * s * MXs * MUs * UU * AXG -
                4 * s * pow(MXs, 2) * UU * AXG + 4 * s * pow(MX, 4) * UU * AXG -
                4 * s * u * MUs * SS - 4 * s * u * MUs * SS * AXG +
                4 * s * u * MUs * UU * AXG + 4 * s * u * MXs * SS -
                4 * s * u * MXs * SS * AXG + 4 * s * u * MXs * UU * AXG +
                4 * s * pow(u, 2) * SS * AXG - 4 * s * pow(u, 2) * UU * AXG -
                4 * pow(s, 2) * MUs * SS + 4 * pow(s, 2) * u * SS -
                4 * pow(s, 2) * u * UU) +
           syFC2 *
               (-4 * MXs * pow(MUs, 2) * SS - 4 * pow(MXs, 2) * MUs * SS +
                8 * pow(MX, 4) * MUs * SS - 4 * pow(MX, 4) * MUs * SS * AXG +
                8 * pow(MX, 4) * MUs * UU * AXG +
                8 * pow(MX, 4) * MXs * UU * AXG - 8 * pow(MX, 6) * UU * AXG +
                4 * u * pow(MUs, 2) * SS + 8 * u * MXs * MUs * SS * AXG -
                16 * u * MXs * MUs * UU * AXG + 4 * u * pow(MXs, 2) * SS -
                16 * u * pow(MXs, 2) * UU * AXG - 8 * u * pow(MX, 4) * SS +
                4 * u * pow(MX, 4) * SS * AXG + 8 * u * pow(MX, 4) * UU * AXG -
                4 * pow(u, 2) * MUs * SS - 4 * pow(u, 2) * MUs * SS * AXG +
                8 * pow(u, 2) * MUs * UU * AXG + 4 * pow(u, 2) * MXs * SS -
                8 * pow(u, 2) * MXs * SS * AXG +
                16 * pow(u, 2) * MXs * UU * AXG + 4 * pow(u, 3) * SS * AXG -
                8 * pow(u, 3) * UU * AXG + 6 * s * pow(MUs, 2) * SS -
                2 * s * MXs * MUs * SS + 6 * s * MXs * MUs * SS * AXG +
                2 * s * MXs * MUs * UU - 8 * s * MXs * MUs * UU * AXG +
                2 * s * pow(MXs, 2) * UU - 8 * s * pow(MXs, 2) * UU * AXG -
                2 * s * pow(MX, 4) * UU + 8 * s * pow(MX, 4) * UU * AXG -
                10 * s * u * MUs * SS - 6 * s * u * MUs * SS * AXG -
                2 * s * u * MUs * UU + 8 * s * u * MUs * UU * AXG +
                2 * s * u * MXs * SS - 6 * s * u * MXs * SS * AXG +
                2 * s * u * MXs * UU + 12 * s * u * MXs * UU * AXG +
                4 * s * pow(u, 2) * SS + 6 * s * pow(u, 2) * SS * AXG -
                2 * s * pow(u, 2) * UU - 12 * s * pow(u, 2) * UU * AXG -
                6 * pow(s, 2) * MUs * SS + 6 * pow(s, 2) * u * SS -
                4 * pow(s, 2) * u * UU - 2 * pow(s, 2) * u * UU * AXG) +
           syFC3 *
               (2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
                4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
                2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
                8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
                8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
                2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
                2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
                4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                4 * pow(u, 2) * MXs * SS * AXG -
                8 * pow(u, 2) * MXs * UU * AXG - 2 * pow(u, 3) * SS * AXG +
                4 * pow(u, 3) * UU * AXG - 6 * s * pow(MUs, 2) * SS +
                4 * s * MXs * MUs * SS - 6 * s * MXs * MUs * SS * AXG -
                s * MXs * MUs * UU + 7 * s * MXs * MUs * UU * AXG -
                s * pow(MXs, 2) * UU + 7 * s * pow(MXs, 2) * UU * AXG +
                s * pow(MX, 4) * UU - 7 * s * pow(MX, 4) * UU * AXG +
                8 * s * u * MUs * SS + 6 * s * u * MUs * SS * AXG +
                s * u * MUs * UU - 7 * s * u * MUs * UU * AXG -
                4 * s * u * MXs * SS + 6 * s * u * MXs * SS * AXG -
                s * u * MXs * UU - 9 * s * u * MXs * UU * AXG -
                2 * s * pow(u, 2) * SS - 6 * s * pow(u, 2) * SS * AXG +
                s * pow(u, 2) * UU + 9 * s * pow(u, 2) * UU * AXG +
                6 * pow(s, 2) * MUs * SS - 6 * pow(s, 2) * u * SS +
                5 * pow(s, 2) * u * UU + pow(s, 2) * u * UU * AXG) +
           syFC4 *
               (4 * s * MXs * pow(MUs, 2) * SS +
                4 * s * pow(MXs, 2) * MUs * SS - 8 * s * pow(MX, 4) * MUs * SS +
                4 * s * pow(MX, 4) * MUs * SS * AXG -
                8 * s * pow(MX, 4) * MUs * UU * AXG -
                8 * s * pow(MX, 4) * MXs * UU * AXG +
                8 * s * pow(MX, 6) * UU * AXG - 4 * s * u * pow(MUs, 2) * SS -
                8 * s * u * MXs * MUs * SS * AXG +
                16 * s * u * MXs * MUs * UU * AXG -
                4 * s * u * pow(MXs, 2) * SS +
                16 * s * u * pow(MXs, 2) * UU * AXG +
                8 * s * u * pow(MX, 4) * SS -
                4 * s * u * pow(MX, 4) * SS * AXG -
                8 * s * u * pow(MX, 4) * UU * AXG +
                4 * s * pow(u, 2) * MUs * SS +
                4 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU * AXG -
                4 * s * pow(u, 2) * MXs * SS +
                8 * s * pow(u, 2) * MXs * SS * AXG -
                16 * s * pow(u, 2) * MXs * UU * AXG -
                4 * s * pow(u, 3) * SS * AXG + 8 * s * pow(u, 3) * UU * AXG -
                6 * pow(s, 2) * pow(MUs, 2) * SS -
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                6 * pow(s, 2) * MXs * MUs * UU * AXG +
                6 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                6 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                12 * pow(s, 2) * u * MUs * SS +
                4 * pow(s, 2) * u * MUs * SS * AXG -
                6 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS * AXG -
                4 * pow(s, 2) * u * MXs * UU -
                10 * pow(s, 2) * u * MXs * UU * AXG -
                6 * pow(s, 2) * pow(u, 2) * SS -
                4 * pow(s, 2) * pow(u, 2) * SS * AXG +
                4 * pow(s, 2) * pow(u, 2) * UU +
                10 * pow(s, 2) * pow(u, 2) * UU * AXG +
                6 * pow(s, 3) * MUs * SS - 6 * pow(s, 3) * u * SS) +
           syFC4 * (4 * pow(s, 3) * u * UU + 2 * pow(s, 3) * u * UU * AXG) +
           syFC5 *
               (2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
                4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
                2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
                8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
                8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
                2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
                2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
                4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                4 * pow(u, 2) * MXs * SS * AXG -
                8 * pow(u, 2) * MXs * UU * AXG - 2 * pow(u, 3) * SS * AXG +
                4 * pow(u, 3) * UU * AXG - 2 * s * MXs * MUs * SS -
                s * MXs * MUs * UU + s * MXs * MUs * UU * AXG -
                s * pow(MXs, 2) * UU + s * pow(MXs, 2) * UU * AXG +
                s * pow(MX, 4) * UU - s * pow(MX, 4) * UU * AXG +
                2 * s * u * MUs * SS + s * u * MUs * UU -
                s * u * MUs * UU * AXG + 2 * s * u * MXs * SS -
                s * u * MXs * UU - 3 * s * u * MXs * UU * AXG -
                2 * s * pow(u, 2) * SS + s * pow(u, 2) * UU +
                3 * s * pow(u, 2) * UU * AXG - pow(s, 2) * u * UU +
                pow(s, 2) * u * UU * AXG) +
           syFC6 *
               (2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
                4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
                2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
                8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
                8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
                2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
                2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
                4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                4 * pow(u, 2) * MXs * SS * AXG -
                8 * pow(u, 2) * MXs * UU * AXG - 2 * pow(u, 3) * SS * AXG +
                4 * pow(u, 3) * UU * AXG - s * MXs * MUs * UU +
                s * MXs * MUs * UU * AXG - s * pow(MXs, 2) * UU +
                s * pow(MXs, 2) * UU * AXG - s * pow(MX, 4) * UU +
                s * pow(MX, 4) * UU * AXG + 2 * s * u * MUs * SS +
                s * u * MUs * UU - s * u * MUs * UU * AXG + s * u * MXs * UU -
                5 * s * u * MXs * UU * AXG - 2 * s * pow(u, 2) * SS +
                s * pow(u, 2) * UU + 3 * s * pow(u, 2) * UU * AXG +
                pow(s, 2) * MXs * UU - pow(s, 2) * MXs * UU * AXG -
                pow(s, 2) * u * UU + pow(s, 2) * u * UU * AXG) +
           syFC7 * (-2 * pow(s, 2) * MXs * MUs * SS +
                    2 * pow(s, 2) * MXs * MUs * SS * AXG +
                    2 * pow(s, 2) * MXs * MUs * UU -
                    2 * pow(s, 2) * MXs * MUs * UU * AXG +
                    2 * pow(s, 2) * pow(MXs, 2) * UU -
                    2 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                    2 * pow(s, 2) * pow(MX, 4) * UU +
                    2 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                    2 * pow(s, 2) * u * MUs * SS -
                    2 * pow(s, 2) * u * MUs * SS * AXG -
                    2 * pow(s, 2) * u * MUs * UU +
                    2 * pow(s, 2) * u * MUs * UU * AXG +
                    2 * pow(s, 2) * u * MXs * SS -
                    2 * pow(s, 2) * u * MXs * SS * AXG -
                    2 * pow(s, 2) * u * MXs * UU +
                    2 * pow(s, 2) * u * MXs * UU * AXG -
                    2 * pow(s, 2) * pow(u, 2) * SS +
                    2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                    2 * pow(s, 2) * pow(u, 2) * UU -
                    2 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC8 *
               (-4 * MXs * pow(MUs, 3) * SS -
                8 * pow(MXs, 2) * pow(MUs, 2) * SS -
                4 * pow(MXs, 3) * MUs * SS +
                16 * pow(MX, 4) * pow(MUs, 2) * SS -
                4 * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                16 * pow(MX, 4) * MXs * MUs * SS -
                4 * pow(MX, 4) * MXs * MUs * SS * AXG +
                16 * pow(MX, 4) * MXs * MUs * UU * AXG +
                8 * pow(MX, 4) * pow(MXs, 2) * UU * AXG -
                16 * pow(MX, 6) * MUs * SS + 8 * pow(MX, 6) * MUs * SS * AXG -
                24 * pow(MX, 6) * MUs * UU * AXG -
                24 * pow(MX, 6) * MXs * UU * AXG + 16 * pow(MX, 8) * UU * AXG +
                4 * u * pow(MUs, 3) * SS - 4 * u * MXs * pow(MUs, 2) * SS +
                8 * u * MXs * pow(MUs, 2) * SS * AXG -
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MXs, 2) * MUs * SS +
                8 * u * pow(MXs, 2) * MUs * SS * AXG -
                32 * u * pow(MXs, 2) * MUs * UU * AXG +
                4 * u * pow(MXs, 3) * SS - 16 * u * pow(MXs, 3) * UU * AXG -
                12 * u * pow(MX, 4) * MUs * SS * AXG +
                40 * u * pow(MX, 4) * MUs * UU * AXG -
                16 * u * pow(MX, 4) * MXs * SS +
                4 * u * pow(MX, 4) * MXs * SS * AXG +
                40 * u * pow(MX, 4) * MXs * UU * AXG +
                16 * u * pow(MX, 6) * SS - 8 * u * pow(MX, 6) * SS * AXG -
                16 * u * pow(MX, 6) * UU * AXG -
                4 * pow(u, 2) * pow(MUs, 2) * SS -
                4 * pow(u, 2) * pow(MUs, 2) * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 2) * MXs * MUs * SS -
                4 * pow(u, 2) * MXs * MUs * SS * AXG +
                8 * pow(u, 2) * MXs * MUs * UU * AXG) +
           syFC8 *
               (12 * pow(u, 2) * pow(MXs, 2) * SS -
                8 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                16 * pow(u, 2) * pow(MX, 4) * SS * AXG -
                16 * pow(u, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(u, 3) * MUs * SS * AXG -
                8 * pow(u, 3) * MUs * UU * AXG -
                4 * pow(u, 3) * MXs * SS * AXG +
                8 * pow(u, 3) * MXs * UU * AXG +
                2 * s * MXs * pow(MUs, 2) * SS +
                2 * s * MXs * pow(MUs, 2) * UU -
                2 * s * MXs * pow(MUs, 2) * UU * AXG +
                2 * s * pow(MXs, 2) * MUs * SS +
                4 * s * pow(MXs, 2) * MUs * UU -
                4 * s * pow(MXs, 2) * MUs * UU * AXG +
                2 * s * pow(MXs, 3) * UU - 2 * s * pow(MXs, 3) * UU * AXG -
                4 * s * pow(MX, 4) * MUs * SS -
                2 * s * pow(MX, 4) * MUs * SS * AXG -
                2 * s * pow(MX, 4) * MUs * UU +
                4 * s * pow(MX, 4) * MUs * UU * AXG -
                2 * s * pow(MX, 4) * MXs * UU +
                4 * s * pow(MX, 4) * MXs * UU * AXG -
                2 * s * pow(MX, 6) * UU * AXG - 2 * s * u * pow(MUs, 2) * SS -
                2 * s * u * pow(MUs, 2) * UU +
                2 * s * u * pow(MUs, 2) * UU * AXG +
                4 * s * u * MXs * MUs * SS * AXG +
                8 * s * u * MXs * MUs * UU * AXG -
                2 * s * u * pow(MXs, 2) * SS + 2 * s * u * pow(MXs, 2) * UU +
                6 * s * u * pow(MXs, 2) * UU * AXG +
                4 * s * u * pow(MX, 4) * SS +
                2 * s * u * pow(MX, 4) * SS * AXG -
                8 * s * u * pow(MX, 4) * UU -
                10 * s * u * pow(MX, 4) * UU * AXG +
                2 * s * pow(u, 2) * MUs * SS -
                2 * s * pow(u, 2) * MUs * SS * AXG -
                2 * s * pow(u, 2) * MUs * UU -
                8 * s * pow(u, 2) * MUs * UU * AXG -
                2 * s * pow(u, 2) * MXs * SS) +
           syFC8 *
               (-4 * s * pow(u, 2) * MXs * SS * AXG +
                6 * s * pow(u, 2) * MXs * UU +
                2 * s * pow(u, 2) * MXs * UU * AXG +
                2 * s * pow(u, 3) * SS * AXG + 2 * s * pow(u, 3) * UU * AXG +
                2 * pow(s, 2) * MXs * MUs * SS -
                2 * pow(s, 2) * MXs * MUs * UU +
                2 * pow(s, 2) * MXs * MUs * UU * AXG -
                2 * pow(s, 2) * pow(MXs, 2) * UU +
                2 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                2 * pow(s, 2) * pow(MX, 4) * UU -
                2 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                2 * pow(s, 2) * u * MUs * SS + 4 * pow(s, 2) * u * MUs * UU -
                4 * pow(s, 2) * u * MUs * UU * AXG -
                2 * pow(s, 2) * u * MXs * SS + 6 * pow(s, 2) * u * MXs * UU -
                4 * pow(s, 2) * u * MXs * UU * AXG +
                2 * pow(s, 2) * pow(u, 2) * SS -
                4 * pow(s, 2) * pow(u, 2) * UU +
                6 * pow(s, 2) * pow(u, 2) * UU * AXG - 2 * pow(s, 3) * u * UU +
                2 * pow(s, 3) * u * UU * AXG) +
           syFC9 *
               (2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
                4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
                2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
                8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
                8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
                2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
                2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
                4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                4 * pow(u, 2) * MXs * SS * AXG -
                8 * pow(u, 2) * MXs * UU * AXG - 2 * pow(u, 3) * SS * AXG +
                4 * pow(u, 3) * UU * AXG - 2 * s * pow(MUs, 2) * SS -
                2 * s * MXs * MUs * SS * AXG - s * MXs * MUs * UU +
                3 * s * MXs * MUs * UU * AXG - s * pow(MXs, 2) * UU +
                3 * s * pow(MXs, 2) * UU * AXG + s * pow(MX, 4) * UU -
                3 * s * pow(MX, 4) * UU * AXG + 4 * s * u * MUs * SS +
                2 * s * u * MUs * SS * AXG + s * u * MUs * UU -
                3 * s * u * MUs * UU * AXG + 2 * s * u * MXs * SS * AXG -
                s * u * MXs * UU - 5 * s * u * MXs * UU * AXG -
                2 * s * pow(u, 2) * SS - 2 * s * pow(u, 2) * SS * AXG +
                s * pow(u, 2) * UU + 5 * s * pow(u, 2) * UU * AXG +
                2 * pow(s, 2) * MUs * SS - 2 * pow(s, 2) * u * SS +
                pow(s, 2) * u * UU + pow(s, 2) * u * UU * AXG) +
           syFC10 * (-2 * s * MXs * MUs * SS + 2 * s * MXs * MUs * SS * AXG +
                     2 * s * MXs * MUs * UU - 2 * s * MXs * MUs * UU * AXG +
                     2 * s * pow(MXs, 2) * UU - 2 * s * pow(MXs, 2) * UU * AXG -
                     2 * s * pow(MX, 4) * UU + 2 * s * pow(MX, 4) * UU * AXG +
                     2 * s * u * MUs * SS - 2 * s * u * MUs * SS * AXG -
                     2 * s * u * MUs * UU + 2 * s * u * MUs * UU * AXG +
                     2 * s * u * MXs * SS - 2 * s * u * MXs * SS * AXG -
                     2 * s * u * MXs * UU + 2 * s * u * MXs * UU * AXG -
                     2 * s * pow(u, 2) * SS + 2 * s * pow(u, 2) * SS * AXG +
                     2 * s * pow(u, 2) * UU - 2 * s * pow(u, 2) * UU * AXG) +
           syFC11 * (2 * pow(s, 2) * pow(MX, 4) * UU -
                     2 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                     4 * pow(s, 2) * u * MXs * UU +
                     4 * pow(s, 2) * u * MXs * UU * AXG +
                     2 * pow(s, 2) * pow(u, 2) * UU -
                     2 * pow(s, 2) * pow(u, 2) * UU * AXG -
                     pow(s, 3) * MXs * UU + pow(s, 3) * MXs * UU * AXG +
                     pow(s, 3) * u * UU - pow(s, 3) * u * UU * AXG) +
           syFC12 *
               (8 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                12 * s * pow(MX, 4) * MUs * SS +
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                8 * s * pow(MX, 4) * MUs * UU -
                16 * s * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(MX, 4) * MXs * UU -
                16 * s * pow(MX, 4) * MXs * UU * AXG - 4 * s * pow(MX, 6) * UU +
                12 * s * pow(MX, 6) * UU * AXG - 8 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * UU +
                8 * s * u * MXs * MUs * UU * AXG -
                8 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG +
                4 * s * u * pow(MXs, 2) * UU +
                12 * s * u * pow(MXs, 2) * UU * AXG +
                12 * s * u * pow(MX, 4) * SS -
                8 * s * u * pow(MX, 4) * SS * AXG -
                4 * s * u * pow(MX, 4) * UU -
                4 * s * u * pow(MX, 4) * UU * AXG +
                12 * s * pow(u, 2) * MUs * SS -
                4 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU +
                4 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU -
                12 * s * pow(u, 2) * MXs * UU * AXG - 4 * s * pow(u, 3) * SS +
                4 * s * pow(u, 3) * UU + 4 * s * pow(u, 3) * UU * AXG) +
           syFC12 * (-6 * pow(s, 2) * pow(MUs, 2) * SS +
                     2 * pow(s, 2) * MXs * MUs * SS -
                     6 * pow(s, 2) * MXs * MUs * SS * AXG -
                     2 * pow(s, 2) * MXs * MUs * UU +
                     8 * pow(s, 2) * MXs * MUs * UU * AXG -
                     2 * pow(s, 2) * pow(MXs, 2) * UU +
                     8 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     2 * pow(s, 2) * pow(MX, 4) * UU -
                     8 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                     14 * pow(s, 2) * u * MUs * SS +
                     2 * pow(s, 2) * u * MUs * SS * AXG -
                     2 * pow(s, 2) * u * MUs * UU -
                     4 * pow(s, 2) * u * MUs * UU * AXG -
                     2 * pow(s, 2) * u * MXs * SS +
                     6 * pow(s, 2) * u * MXs * SS * AXG -
                     2 * pow(s, 2) * u * MXs * UU -
                     12 * pow(s, 2) * u * MXs * UU * AXG -
                     8 * pow(s, 2) * pow(u, 2) * SS -
                     2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                     6 * pow(s, 2) * pow(u, 2) * UU +
                     8 * pow(s, 2) * pow(u, 2) * UU * AXG +
                     6 * pow(s, 3) * MUs * SS - 6 * pow(s, 3) * u * SS +
                     4 * pow(s, 3) * u * UU + 2 * pow(s, 3) * u * UU * AXG) +
           syFC13 *
               (-2 * s * MXs * pow(MUs, 2) * SS -
                2 * s * pow(MXs, 2) * MUs * SS + 4 * s * pow(MX, 4) * MUs * SS -
                2 * s * pow(MX, 4) * MUs * SS * AXG +
                4 * s * pow(MX, 4) * MUs * UU * AXG +
                4 * s * pow(MX, 4) * MXs * UU * AXG -
                4 * s * pow(MX, 6) * UU * AXG + 2 * s * u * pow(MUs, 2) * SS +
                4 * s * u * MXs * MUs * SS * AXG -
                8 * s * u * MXs * MUs * UU * AXG +
                2 * s * u * pow(MXs, 2) * SS -
                8 * s * u * pow(MXs, 2) * UU * AXG -
                4 * s * u * pow(MX, 4) * SS +
                2 * s * u * pow(MX, 4) * SS * AXG +
                4 * s * u * pow(MX, 4) * UU * AXG -
                2 * s * pow(u, 2) * MUs * SS -
                2 * s * pow(u, 2) * MUs * SS * AXG +
                4 * s * pow(u, 2) * MUs * UU * AXG +
                2 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG +
                8 * s * pow(u, 2) * MXs * UU * AXG +
                2 * s * pow(u, 3) * SS * AXG - 4 * s * pow(u, 3) * UU * AXG -
                2 * pow(s, 2) * MXs * MUs * SS +
                2 * pow(s, 2) * MXs * MUs * SS * AXG +
                3 * pow(s, 2) * MXs * MUs * UU -
                3 * pow(s, 2) * MXs * MUs * UU * AXG +
                3 * pow(s, 2) * pow(MXs, 2) * UU -
                3 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                pow(s, 2) * pow(MX, 4) * UU +
                pow(s, 2) * pow(MX, 4) * UU * AXG -
                2 * pow(s, 2) * u * MUs * SS * AXG -
                3 * pow(s, 2) * u * MUs * UU +
                3 * pow(s, 2) * u * MUs * UU * AXG +
                2 * pow(s, 2) * u * MXs * SS -
                2 * pow(s, 2) * u * MXs * SS * AXG -
                3 * pow(s, 2) * u * MXs * UU +
                7 * pow(s, 2) * u * MXs * UU * AXG +
                2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                pow(s, 2) * pow(u, 2) * UU) +
           syFC13 * (-5 * pow(s, 2) * pow(u, 2) * UU * AXG -
                     pow(s, 3) * MXs * UU + pow(s, 3) * MXs * UU * AXG +
                     pow(s, 3) * u * UU - pow(s, 3) * u * UU * AXG) +
           syFC14 *
               (12 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                12 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                20 * s * pow(MX, 4) * MUs * SS +
                12 * s * pow(MX, 4) * MUs * SS * AXG +
                12 * s * pow(MX, 4) * MUs * UU -
                24 * s * pow(MX, 4) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MXs * UU -
                24 * s * pow(MX, 4) * MXs * UU * AXG - 8 * s * pow(MX, 6) * UU +
                20 * s * pow(MX, 6) * UU * AXG - 12 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS - 8 * s * u * MXs * MUs * SS * AXG +
                24 * s * u * MXs * MUs * UU * AXG -
                12 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG -
                4 * s * u * pow(MXs, 2) * UU +
                28 * s * u * pow(MXs, 2) * UU * AXG +
                20 * s * u * pow(MX, 4) * SS -
                12 * s * u * pow(MX, 4) * SS * AXG -
                12 * s * u * pow(MX, 4) * UU * AXG +
                16 * s * pow(u, 2) * MUs * SS - 4 * s * pow(u, 2) * MUs * UU -
                8 * s * pow(u, 2) * MUs * UU * AXG -
                4 * s * pow(u, 2) * MXs * SS +
                12 * s * pow(u, 2) * MXs * SS * AXG +
                4 * s * pow(u, 2) * MXs * UU -
                28 * s * pow(u, 2) * MXs * UU * AXG - 4 * s * pow(u, 3) * SS -
                4 * s * pow(u, 3) * SS * AXG) +
           syFC14 * (12 * s * pow(u, 3) * UU * AXG -
                     6 * pow(s, 2) * pow(MUs, 2) * SS -
                     2 * pow(s, 2) * MXs * MUs * SS -
                     6 * pow(s, 2) * MXs * MUs * SS * AXG -
                     4 * pow(s, 2) * MXs * MUs * UU +
                     10 * pow(s, 2) * MXs * MUs * UU * AXG -
                     4 * pow(s, 2) * pow(MXs, 2) * UU +
                     10 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     4 * pow(s, 2) * pow(MX, 4) * UU -
                     10 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                     18 * pow(s, 2) * u * MUs * SS +
                     2 * pow(s, 2) * u * MUs * SS * AXG -
                     6 * pow(s, 2) * u * MUs * UU * AXG +
                     2 * pow(s, 2) * u * MXs * SS +
                     6 * pow(s, 2) * u * MXs * SS * AXG -
                     18 * pow(s, 2) * u * MXs * UU * AXG -
                     12 * pow(s, 2) * pow(u, 2) * SS -
                     2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                     4 * pow(s, 2) * pow(u, 2) * UU +
                     14 * pow(s, 2) * pow(u, 2) * UU * AXG +
                     6 * pow(s, 3) * MUs * SS - 6 * pow(s, 3) * u * SS +
                     2 * pow(s, 3) * u * UU + 4 * pow(s, 3) * u * UU * AXG) +
           syFC15 *
               (8 * MXs * pow(MUs, 2) * SS + 8 * pow(MXs, 2) * MUs * SS -
                16 * pow(MX, 4) * MUs * SS + 8 * pow(MX, 4) * MUs * SS * AXG -
                16 * pow(MX, 4) * MUs * UU * AXG -
                16 * pow(MX, 4) * MXs * UU * AXG + 16 * pow(MX, 6) * UU * AXG -
                8 * u * pow(MUs, 2) * SS - 16 * u * MXs * MUs * SS * AXG +
                32 * u * MXs * MUs * UU * AXG - 8 * u * pow(MXs, 2) * SS +
                32 * u * pow(MXs, 2) * UU * AXG + 16 * u * pow(MX, 4) * SS -
                8 * u * pow(MX, 4) * SS * AXG - 16 * u * pow(MX, 4) * UU * AXG +
                8 * pow(u, 2) * MUs * SS + 8 * pow(u, 2) * MUs * SS * AXG -
                16 * pow(u, 2) * MUs * UU * AXG - 8 * pow(u, 2) * MXs * SS +
                16 * pow(u, 2) * MXs * SS * AXG -
                32 * pow(u, 2) * MXs * UU * AXG - 8 * pow(u, 3) * SS * AXG +
                16 * pow(u, 3) * UU * AXG - 8 * s * pow(MUs, 2) * SS -
                8 * s * MXs * MUs * SS * AXG - 4 * s * MXs * MUs * UU +
                12 * s * MXs * MUs * UU * AXG - 4 * s * pow(MXs, 2) * UU +
                12 * s * pow(MXs, 2) * UU * AXG + 8 * s * pow(MX, 4) * UU -
                16 * s * pow(MX, 4) * UU * AXG + 16 * s * u * MUs * SS +
                8 * s * u * MUs * SS * AXG + 4 * s * u * MUs * UU -
                12 * s * u * MUs * UU * AXG + 8 * s * u * MXs * SS * AXG -
                12 * s * u * MXs * UU - 12 * s * u * MXs * UU * AXG -
                8 * s * pow(u, 2) * SS - 8 * s * pow(u, 2) * SS * AXG +
                8 * s * pow(u, 2) * UU + 16 * s * pow(u, 2) * UU * AXG +
                8 * pow(s, 2) * MUs * SS - 2 * pow(s, 2) * MXs * UU +
                2 * pow(s, 2) * MXs * UU * AXG - 8 * pow(s, 2) * u * SS +
                6 * pow(s, 2) * u * UU + 2 * pow(s, 2) * u * UU * AXG) +
           syFC16 *
               (2 * pow(MX, 4) * pow(MUs, 2) * SS +
                2 * pow(MX, 4) * MXs * MUs * SS - 4 * pow(MX, 6) * MUs * SS +
                2 * pow(MX, 6) * MUs * SS * AXG -
                4 * pow(MX, 6) * MUs * UU * AXG -
                4 * pow(MX, 6) * MXs * UU * AXG + 4 * pow(MX, 8) * UU * AXG -
                2 * u * pow(MX, 4) * MUs * SS -
                2 * u * pow(MX, 4) * MUs * SS * AXG +
                4 * u * pow(MX, 4) * MUs * UU * AXG -
                2 * u * pow(MX, 4) * MXs * SS +
                4 * u * pow(MX, 4) * MXs * UU * AXG + 4 * u * pow(MX, 6) * SS -
                2 * u * pow(MX, 6) * SS * AXG -
                2 * pow(u, 2) * pow(MUs, 2) * SS +
                2 * pow(u, 2) * MXs * MUs * SS -
                2 * pow(u, 2) * MXs * MUs * SS * AXG +
                4 * pow(u, 2) * MXs * MUs * UU * AXG +
                4 * pow(u, 2) * pow(MXs, 2) * UU * AXG +
                2 * pow(u, 2) * pow(MX, 4) * SS * AXG -
                8 * pow(u, 2) * pow(MX, 4) * UU * AXG +
                2 * pow(u, 3) * MUs * SS + 2 * pow(u, 3) * MUs * SS * AXG -
                4 * pow(u, 3) * MUs * UU * AXG - 2 * pow(u, 3) * MXs * SS +
                2 * pow(u, 3) * MXs * SS * AXG -
                4 * pow(u, 3) * MXs * UU * AXG - 2 * pow(u, 4) * SS * AXG +
                4 * pow(u, 4) * UU * AXG + 4 * s * pow(MX, 4) * MUs * SS -
                4 * s * pow(MX, 4) * MUs * SS * AXG -
                s * pow(MX, 4) * MUs * UU +
                9 * s * pow(MX, 4) * MUs * UU * AXG -
                s * pow(MX, 4) * MXs * UU +
                9 * s * pow(MX, 4) * MXs * UU * AXG - s * pow(MX, 6) * UU -
                7 * s * pow(MX, 6) * UU * AXG - 2 * s * u * MXs * MUs * SS +
                4 * s * u * MXs * MUs * SS * AXG -
                8 * s * u * MXs * MUs * UU * AXG -
                8 * s * u * pow(MXs, 2) * UU * AXG) +
           syFC16 *
               (-4 * s * u * pow(MX, 4) * SS +
                4 * s * u * pow(MX, 4) * SS * AXG - s * u * pow(MX, 4) * UU -
                3 * s * u * pow(MX, 4) * UU * AXG +
                2 * s * pow(u, 2) * MUs * SS + s * pow(u, 2) * MUs * UU -
                s * pow(u, 2) * MUs * UU * AXG + 2 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG +
                2 * s * pow(u, 2) * MXs * UU +
                6 * s * pow(u, 2) * MXs * UU * AXG - 2 * s * pow(u, 3) * SS +
                s * pow(u, 3) * UU + 3 * s * pow(u, 3) * UU * AXG +
                2 * pow(s, 2) * MXs * MUs * UU -
                2 * pow(s, 2) * MXs * MUs * UU * AXG +
                2 * pow(s, 2) * pow(MXs, 2) * UU -
                2 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                3 * pow(s, 2) * pow(MX, 4) * UU -
                3 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                2 * pow(s, 2) * u * MXs * UU +
                6 * pow(s, 2) * u * MXs * UU * AXG -
                pow(s, 2) * pow(u, 2) * UU + pow(s, 2) * pow(u, 2) * UU * AXG -
                2 * pow(s, 3) * MXs * UU + 2 * pow(s, 3) * MXs * UU * AXG) +
           syFC17 *
               (-4 * MXs * pow(MUs, 3) * SS -
                8 * pow(MXs, 2) * pow(MUs, 2) * SS -
                4 * pow(MXs, 3) * MUs * SS +
                12 * pow(MX, 4) * pow(MUs, 2) * SS -
                4 * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                12 * pow(MX, 4) * MXs * MUs * SS -
                4 * pow(MX, 4) * MXs * MUs * SS * AXG +
                16 * pow(MX, 4) * MXs * MUs * UU * AXG +
                8 * pow(MX, 4) * pow(MXs, 2) * UU * AXG -
                8 * pow(MX, 6) * MUs * SS + 4 * pow(MX, 6) * MUs * SS * AXG -
                16 * pow(MX, 6) * MUs * UU * AXG -
                16 * pow(MX, 6) * MXs * UU * AXG + 8 * pow(MX, 8) * UU * AXG +
                4 * u * pow(MUs, 3) * SS + 4 * u * MXs * pow(MUs, 2) * SS +
                8 * u * MXs * pow(MUs, 2) * SS * AXG -
                16 * u * MXs * pow(MUs, 2) * UU * AXG +
                4 * u * pow(MXs, 2) * MUs * SS +
                8 * u * pow(MXs, 2) * MUs * SS * AXG -
                32 * u * pow(MXs, 2) * MUs * UU * AXG +
                4 * u * pow(MXs, 3) * SS - 16 * u * pow(MXs, 3) * UU * AXG -
                12 * u * pow(MX, 4) * MUs * SS +
                16 * u * pow(MX, 4) * MUs * UU * AXG -
                12 * u * pow(MX, 4) * MXs * SS +
                4 * u * pow(MX, 4) * MXs * SS * AXG +
                16 * u * pow(MX, 4) * MXs * UU * AXG + 8 * u * pow(MX, 6) * SS -
                4 * u * pow(MX, 6) * SS * AXG -
                8 * pow(u, 2) * pow(MUs, 2) * SS -
                4 * pow(u, 2) * pow(MUs, 2) * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                4 * pow(u, 2) * MXs * MUs * SS -
                16 * pow(u, 2) * MXs * MUs * SS * AXG +
                32 * pow(u, 2) * MXs * MUs * UU * AXG +
                4 * pow(u, 2) * pow(MXs, 2) * SS) +
           syFC17 *
               (-8 * pow(u, 2) * pow(MXs, 2) * SS * AXG +
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG +
                4 * pow(u, 2) * pow(MX, 4) * SS * AXG -
                16 * pow(u, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(u, 3) * MUs * SS + 8 * pow(u, 3) * MUs * SS * AXG -
                16 * pow(u, 3) * MUs * UU * AXG - 4 * pow(u, 3) * MXs * SS +
                8 * pow(u, 3) * MXs * SS * AXG -
                16 * pow(u, 3) * MXs * UU * AXG - 4 * pow(u, 4) * SS * AXG +
                8 * pow(u, 4) * UU * AXG + 2 * s * MXs * pow(MUs, 2) * SS +
                2 * s * MXs * pow(MUs, 2) * UU -
                2 * s * MXs * pow(MUs, 2) * UU * AXG +
                2 * s * pow(MXs, 2) * MUs * SS +
                4 * s * pow(MXs, 2) * MUs * UU -
                4 * s * pow(MXs, 2) * MUs * UU * AXG +
                2 * s * pow(MXs, 3) * UU - 2 * s * pow(MXs, 3) * UU * AXG -
                2 * s * pow(MX, 4) * MUs * SS + 2 * s * pow(MX, 4) * MUs * UU +
                2 * s * pow(MX, 4) * MUs * UU * AXG +
                2 * s * pow(MX, 4) * MXs * UU +
                2 * s * pow(MX, 4) * MXs * UU * AXG - 4 * s * pow(MX, 6) * UU -
                6 * s * u * pow(MUs, 2) * SS - 2 * s * u * pow(MUs, 2) * UU +
                2 * s * u * pow(MUs, 2) * UU * AXG -
                4 * s * u * MXs * MUs * SS - 4 * s * u * MXs * MUs * UU +
                12 * s * u * MXs * MUs * UU * AXG -
                2 * s * u * pow(MXs, 2) * SS - 2 * s * u * pow(MXs, 2) * UU +
                10 * s * u * pow(MXs, 2) * UU * AXG +
                2 * s * u * pow(MX, 4) * SS - 4 * s * u * pow(MX, 4) * UU -
                8 * s * u * pow(MX, 4) * UU * AXG +
                12 * s * pow(u, 2) * MUs * SS - 2 * s * pow(u, 2) * MUs * UU -
                10 * s * pow(u, 2) * MUs * UU * AXG +
                2 * s * pow(u, 2) * MXs * SS) +
           syFC17 *
               (2 * s * pow(u, 2) * MXs * UU -
                10 * s * pow(u, 2) * MXs * UU * AXG - 6 * s * pow(u, 3) * SS +
                4 * s * pow(u, 3) * UU + 8 * s * pow(u, 3) * UU * AXG +
                4 * pow(s, 2) * MXs * MUs * SS -
                4 * pow(s, 2) * MXs * MUs * SS * AXG -
                5 * pow(s, 2) * MXs * MUs * UU +
                5 * pow(s, 2) * MXs * MUs * UU * AXG -
                5 * pow(s, 2) * pow(MXs, 2) * UU +
                5 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                7 * pow(s, 2) * pow(MX, 4) * UU -
                7 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                2 * pow(s, 2) * u * MUs * SS + 3 * pow(s, 2) * u * MUs * UU -
                3 * pow(s, 2) * u * MUs * UU * AXG -
                4 * pow(s, 2) * u * MXs * SS +
                4 * pow(s, 2) * u * MXs * SS * AXG +
                9 * pow(s, 2) * u * MXs * UU -
                5 * pow(s, 2) * u * MXs * UU * AXG -
                2 * pow(s, 2) * pow(u, 2) * SS - pow(s, 2) * pow(u, 2) * UU +
                5 * pow(s, 2) * pow(u, 2) * UU * AXG - pow(s, 3) * MXs * UU +
                pow(s, 3) * MXs * UU * AXG - pow(s, 3) * u * UU +
                pow(s, 3) * u * UU * AXG) +
           syFC18 *
               (-4 * MXs * pow(MUs, 3) * SS -
                8 * pow(MXs, 2) * pow(MUs, 2) * SS -
                4 * pow(MXs, 3) * MUs * SS +
                16 * pow(MX, 4) * pow(MUs, 2) * SS -
                4 * pow(MX, 4) * pow(MUs, 2) * SS * AXG +
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                16 * pow(MX, 4) * MXs * MUs * SS -
                4 * pow(MX, 4) * MXs * MUs * SS * AXG +
                16 * pow(MX, 4) * MXs * MUs * UU * AXG +
                8 * pow(MX, 4) * pow(MXs, 2) * UU * AXG -
                16 * pow(MX, 6) * MUs * SS + 8 * pow(MX, 6) * MUs * SS * AXG -
                24 * pow(MX, 6) * MUs * UU * AXG -
                24 * pow(MX, 6) * MXs * UU * AXG + 16 * pow(MX, 8) * UU * AXG +
                4 * u * pow(MUs, 3) * SS - 4 * u * MXs * pow(MUs, 2) * SS +
                8 * u * MXs * pow(MUs, 2) * SS * AXG -
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MXs, 2) * MUs * SS +
                8 * u * pow(MXs, 2) * MUs * SS * AXG -
                32 * u * pow(MXs, 2) * MUs * UU * AXG +
                4 * u * pow(MXs, 3) * SS - 16 * u * pow(MXs, 3) * UU * AXG -
                12 * u * pow(MX, 4) * MUs * SS * AXG +
                40 * u * pow(MX, 4) * MUs * UU * AXG -
                16 * u * pow(MX, 4) * MXs * SS +
                4 * u * pow(MX, 4) * MXs * SS * AXG +
                40 * u * pow(MX, 4) * MXs * UU * AXG +
                16 * u * pow(MX, 6) * SS - 8 * u * pow(MX, 6) * SS * AXG -
                16 * u * pow(MX, 6) * UU * AXG -
                4 * pow(u, 2) * pow(MUs, 2) * SS -
                4 * pow(u, 2) * pow(MUs, 2) * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                8 * pow(u, 2) * MXs * MUs * SS -
                4 * pow(u, 2) * MXs * MUs * SS * AXG +
                8 * pow(u, 2) * MXs * MUs * UU * AXG) +
           syFC18 *
               (12 * pow(u, 2) * pow(MXs, 2) * SS -
                8 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                16 * pow(u, 2) * pow(MX, 4) * SS * AXG -
                16 * pow(u, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(u, 3) * MUs * SS * AXG -
                8 * pow(u, 3) * MUs * UU * AXG -
                4 * pow(u, 3) * MXs * SS * AXG +
                8 * pow(u, 3) * MXs * UU * AXG -
                2 * s * MXs * pow(MUs, 2) * SS +
                2 * s * MXs * pow(MUs, 2) * UU -
                2 * s * MXs * pow(MUs, 2) * UU * AXG -
                2 * s * pow(MXs, 2) * MUs * SS +
                4 * s * pow(MXs, 2) * MUs * UU -
                4 * s * pow(MXs, 2) * MUs * UU * AXG +
                2 * s * pow(MXs, 3) * UU - 2 * s * pow(MXs, 3) * UU * AXG +
                6 * s * pow(MX, 4) * MUs * SS -
                4 * s * pow(MX, 4) * MUs * SS * AXG +
                12 * s * pow(MX, 4) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MXs * UU * AXG - 8 * s * pow(MX, 6) * UU -
                4 * s * pow(MX, 6) * UU * AXG - 2 * s * u * pow(MUs, 2) * SS -
                2 * s * u * pow(MUs, 2) * UU +
                2 * s * u * pow(MUs, 2) * UU * AXG +
                8 * s * u * MXs * MUs * SS * AXG -
                8 * s * u * MXs * MUs * UU * AXG +
                2 * s * u * pow(MXs, 2) * SS + 2 * s * u * pow(MXs, 2) * UU -
                10 * s * u * pow(MXs, 2) * UU * AXG -
                6 * s * u * pow(MX, 4) * SS +
                4 * s * u * pow(MX, 4) * SS * AXG -
                12 * s * u * pow(MX, 4) * UU * AXG +
                4 * s * pow(u, 2) * MUs * SS -
                4 * s * pow(u, 2) * MUs * SS * AXG -
                4 * s * pow(u, 2) * MUs * UU + 2 * s * pow(u, 2) * MXs * SS -
                8 * s * pow(u, 2) * MXs * SS * AXG +
                4 * s * pow(u, 2) * MXs * UU +
                20 * s * pow(u, 2) * MXs * UU * AXG) +
           syFC18 *
               (-2 * s * pow(u, 3) * SS + 4 * s * pow(u, 3) * SS * AXG -
                4 * s * pow(u, 3) * UU * AXG + 4 * pow(s, 2) * MXs * MUs * SS -
                4 * pow(s, 2) * MXs * MUs * SS * AXG -
                3 * pow(s, 2) * MXs * MUs * UU +
                3 * pow(s, 2) * MXs * MUs * UU * AXG -
                3 * pow(s, 2) * pow(MXs, 2) * UU +
                3 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                16 * pow(s, 2) * pow(MX, 4) * UU -
                16 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                2 * pow(s, 2) * u * MUs * SS + pow(s, 2) * u * MUs * UU -
                pow(s, 2) * u * MUs * UU * AXG - 4 * pow(s, 2) * u * MXs * SS +
                4 * pow(s, 2) * u * MXs * SS * AXG - pow(s, 2) * u * MXs * UU +
                13 * pow(s, 2) * u * MXs * UU * AXG +
                2 * pow(s, 2) * pow(u, 2) * SS -
                2 * pow(s, 2) * pow(u, 2) * UU -
                2 * pow(s, 2) * pow(u, 2) * UU * AXG -
                5 * pow(s, 3) * MXs * UU + 5 * pow(s, 3) * MXs * UU * AXG +
                pow(s, 3) * u * UU - pow(s, 3) * u * UU * AXG) +
           syFC19 *
               (4 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                4 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                4 * s * pow(MX, 4) * MUs * SS +
                4 * s * pow(MX, 4) * MUs * SS * AXG +
                8 * s * pow(MX, 4) * MUs * UU -
                8 * s * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(MX, 4) * MXs * UU -
                8 * s * pow(MX, 4) * MXs * UU * AXG - 4 * s * pow(MX, 6) * UU +
                4 * s * pow(MX, 6) * UU * AXG - 4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * SS * AXG +
                8 * s * u * MXs * MUs * UU - 8 * s * u * MXs * MUs * UU * AXG -
                4 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG +
                4 * s * u * pow(MXs, 2) * UU -
                4 * s * u * pow(MXs, 2) * UU * AXG +
                4 * s * u * pow(MX, 4) * SS -
                4 * s * u * pow(MX, 4) * SS * AXG -
                4 * s * u * pow(MX, 4) * UU +
                4 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * pow(u, 2) * MUs * SS -
                8 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU +
                8 * s * pow(u, 2) * MUs * UU * AXG +
                4 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU +
                4 * s * pow(u, 2) * MXs * UU * AXG) +
           syFC19 * (-4 * s * pow(u, 3) * SS + 4 * s * pow(u, 3) * SS * AXG +
                     4 * s * pow(u, 3) * UU - 4 * s * pow(u, 3) * UU * AXG +
                     2 * pow(s, 2) * MXs * MUs * SS -
                     2 * pow(s, 2) * MXs * MUs * SS * AXG -
                     2 * pow(s, 2) * MXs * MUs * UU +
                     2 * pow(s, 2) * MXs * MUs * UU * AXG -
                     2 * pow(s, 2) * pow(MXs, 2) * UU +
                     2 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     2 * pow(s, 2) * pow(MX, 4) * UU -
                     2 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                     2 * pow(s, 2) * u * MUs * SS -
                     2 * pow(s, 2) * u * MUs * SS * AXG -
                     2 * pow(s, 2) * u * MUs * UU +
                     2 * pow(s, 2) * u * MUs * UU * AXG -
                     2 * pow(s, 2) * u * MXs * SS +
                     2 * pow(s, 2) * u * MXs * SS * AXG +
                     2 * pow(s, 2) * u * MXs * UU -
                     2 * pow(s, 2) * u * MXs * UU * AXG -
                     2 * pow(s, 2) * pow(u, 2) * SS +
                     2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                     2 * pow(s, 2) * pow(u, 2) * UU -
                     2 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC20 *
               (8 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * SS * AXG -
                8 * s * MXs * pow(MUs, 2) * UU +
                8 * s * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                16 * s * pow(MXs, 2) * MUs * UU * AXG -
                8 * s * pow(MXs, 3) * UU + 8 * s * pow(MXs, 3) * UU * AXG -
                12 * s * pow(MX, 4) * MUs * SS +
                12 * s * pow(MX, 4) * MUs * SS * AXG +
                24 * s * pow(MX, 4) * MUs * UU -
                24 * s * pow(MX, 4) * MUs * UU * AXG +
                24 * s * pow(MX, 4) * MXs * UU -
                24 * s * pow(MX, 4) * MXs * UU * AXG -
                16 * s * pow(MX, 6) * UU + 16 * s * pow(MX, 6) * UU * AXG -
                8 * s * u * pow(MUs, 2) * SS +
                8 * s * u * pow(MUs, 2) * SS * AXG +
                8 * s * u * pow(MUs, 2) * UU -
                8 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * SS * AXG -
                8 * s * u * pow(MXs, 2) * SS +
                8 * s * u * pow(MXs, 2) * SS * AXG -
                8 * s * u * pow(MXs, 2) * UU +
                8 * s * u * pow(MXs, 2) * UU * AXG +
                12 * s * u * pow(MX, 4) * SS -
                12 * s * u * pow(MX, 4) * SS * AXG +
                12 * s * pow(u, 2) * MUs * SS -
                12 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU +
                8 * s * pow(u, 2) * MUs * UU * AXG +
                8 * s * pow(u, 2) * MXs * UU -
                8 * s * pow(u, 2) * MXs * UU * AXG - 4 * s * pow(u, 3) * SS +
                4 * s * pow(u, 3) * SS * AXG + 8 * pow(s, 2) * MXs * MUs * SS -
                8 * pow(s, 2) * MXs * MUs * SS * AXG -
                10 * pow(s, 2) * MXs * MUs * UU) +
           syFC20 * (10 * pow(s, 2) * MXs * MUs * UU * AXG -
                     10 * pow(s, 2) * pow(MXs, 2) * UU +
                     10 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     12 * pow(s, 2) * pow(MX, 4) * UU -
                     12 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                     2 * pow(s, 2) * u * MUs * UU -
                     2 * pow(s, 2) * u * MUs * UU * AXG -
                     8 * pow(s, 2) * u * MXs * SS +
                     8 * pow(s, 2) * u * MXs * SS * AXG +
                     10 * pow(s, 2) * u * MXs * UU -
                     10 * pow(s, 2) * u * MXs * UU * AXG -
                     4 * pow(s, 2) * pow(u, 2) * UU +
                     4 * pow(s, 2) * pow(u, 2) * UU * AXG -
                     pow(s, 3) * MXs * UU + pow(s, 3) * MXs * UU * AXG -
                     pow(s, 3) * u * UU + pow(s, 3) * u * UU * AXG) +
           syFC21 *
               (4 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                4 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                8 * s * pow(MX, 4) * MUs * SS +
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                16 * s * pow(MX, 4) * MUs * UU -
                16 * s * pow(MX, 4) * MUs * UU * AXG +
                16 * s * pow(MX, 4) * MXs * UU -
                16 * s * pow(MX, 4) * MXs * UU * AXG -
                16 * s * pow(MX, 6) * UU + 16 * s * pow(MX, 6) * UU * AXG -
                4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * UU + 8 * s * u * MXs * MUs * UU * AXG -
                4 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG -
                12 * s * u * pow(MXs, 2) * UU +
                12 * s * u * pow(MXs, 2) * UU * AXG +
                8 * s * u * pow(MX, 4) * SS -
                8 * s * u * pow(MX, 4) * SS * AXG +
                16 * s * u * pow(MX, 4) * UU -
                16 * s * u * pow(MX, 4) * UU * AXG +
                4 * s * pow(u, 2) * MUs * SS -
                4 * s * pow(u, 2) * MUs * SS * AXG -
                4 * s * pow(u, 2) * MXs * SS +
                4 * s * pow(u, 2) * MXs * SS * AXG +
                6 * pow(s, 2) * MXs * MUs * SS -
                6 * pow(s, 2) * MXs * MUs * SS * AXG -
                8 * pow(s, 2) * MXs * MUs * UU +
                8 * pow(s, 2) * MXs * MUs * UU * AXG -
                8 * pow(s, 2) * pow(MXs, 2) * UU) +
           syFC21 * (8 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     16 * pow(s, 2) * pow(MX, 4) * UU -
                     16 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                     2 * pow(s, 2) * u * MUs * SS +
                     2 * pow(s, 2) * u * MUs * SS * AXG +
                     4 * pow(s, 2) * u * MUs * UU -
                     4 * pow(s, 2) * u * MUs * UU * AXG -
                     6 * pow(s, 2) * u * MXs * SS +
                     6 * pow(s, 2) * u * MXs * SS * AXG -
                     4 * pow(s, 2) * u * MXs * UU +
                     4 * pow(s, 2) * u * MXs * UU * AXG +
                     2 * pow(s, 2) * pow(u, 2) * SS -
                     2 * pow(s, 2) * pow(u, 2) * SS * AXG -
                     3 * pow(s, 3) * MXs * UU + 3 * pow(s, 3) * MXs * UU * AXG +
                     pow(s, 3) * u * UU - pow(s, 3) * u * UU * AXG)));

  return ret.real();
}

ComplexType ME_us_box_qgQQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2,
                           double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C0(MUs, 0, u, 0, MQs, MQs)
#define syFC2 D0(MUs, 0, MXs, 0, u, MUs + MXs - s - u, 0, MQs, MQs, 0)
#define syFC3 C1(MUs, u, 0, MQs, 0, MQs)
#define syFC4 D1(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC5 C2(MUs, u, 0, MQs, 0, MQs)
#define syFC6 D2(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC7 D3(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC8 D00(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC9 D11(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC10 D12(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC11 D13(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC12 D22(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC13 D23(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)
#define syFC14 D33(MUs, u, MXs, MUs + MXs - s - u, 0, 0, MQs, 0, MQs, 0)

  _EPS0(
      ret,
      (-1 + Nc) * (+1 + Nc) *
          (+Denom(768 * pow(Pi, 2) * s * MUs * Nc -
                  768 * pow(Pi, 2) * s * u * Nc)) *
          (+TR) * (+TR) * (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 * (-s * MXs * MUs * SS + s * MXs * MUs * SS * AXG +
                     s * MXs * MUs * UU - s * MXs * MUs * UU * AXG +
                     s * pow(MXs, 2) * UU - s * pow(MXs, 2) * UU * AXG -
                     s * pow(MX, 4) * UU + s * pow(MX, 4) * UU * AXG +
                     s * u * MUs * SS - s * u * MUs * SS * AXG -
                     s * u * MUs * UU + s * u * MUs * UU * AXG +
                     s * u * MXs * SS - s * u * MXs * SS * AXG -
                     s * u * MXs * UU + s * u * MXs * UU * AXG -
                     s * pow(u, 2) * SS + s * pow(u, 2) * SS * AXG +
                     s * pow(u, 2) * UU - s * pow(u, 2) * UU * AXG) +
           syFC2 *
               (-2 * s * MXs * pow(MUs, 2) * UU +
                2 * s * MXs * pow(MUs, 2) * UU * AXG -
                2 * s * pow(MXs, 2) * MUs * UU +
                2 * s * pow(MXs, 2) * MUs * UU * AXG +
                2 * s * pow(MX, 4) * MUs * SS -
                2 * s * pow(MX, 4) * MUs * SS * AXG -
                2 * s * pow(MX, 4) * MXs * UU +
                2 * s * pow(MX, 4) * MXs * UU * AXG + 2 * s * pow(MX, 6) * UU -
                2 * s * pow(MX, 6) * UU * AXG + 2 * s * pow(MU, 4) * MXs * SS -
                2 * s * pow(MU, 4) * MXs * SS * AXG +
                2 * s * u * pow(MUs, 2) * UU -
                2 * s * u * pow(MUs, 2) * UU * AXG -
                6 * s * u * MXs * MUs * SS + 6 * s * u * MXs * MUs * SS * AXG +
                6 * s * u * MXs * MUs * UU - 6 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * UU -
                4 * s * u * pow(MXs, 2) * UU * AXG -
                2 * s * u * pow(MX, 4) * SS +
                2 * s * u * pow(MX, 4) * SS * AXG -
                2 * s * u * pow(MX, 4) * UU +
                2 * s * u * pow(MX, 4) * UU * AXG -
                2 * s * u * pow(MU, 4) * SS +
                2 * s * u * pow(MU, 4) * SS * AXG +
                4 * s * pow(u, 2) * MUs * SS -
                4 * s * pow(u, 2) * MUs * SS * AXG -
                4 * s * pow(u, 2) * MUs * UU +
                4 * s * pow(u, 2) * MUs * UU * AXG +
                4 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU +
                4 * s * pow(u, 2) * MXs * UU * AXG - 2 * s * pow(u, 3) * SS +
                2 * s * pow(u, 3) * SS * AXG + 2 * s * pow(u, 3) * UU -
                2 * s * pow(u, 3) * UU * AXG - 2 * pow(s, 2) * MXs * MUs * SS +
                2 * pow(s, 2) * MXs * MUs * SS * AXG +
                2 * pow(s, 2) * MXs * MUs * UU -
                2 * pow(s, 2) * MXs * MUs * UU * AXG +
                2 * pow(s, 2) * pow(MXs, 2) * UU) +
           syFC2 * (-2 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                    2 * pow(s, 2) * pow(MX, 4) * UU +
                    2 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                    2 * pow(s, 2) * u * MUs * SS -
                    2 * pow(s, 2) * u * MUs * SS * AXG -
                    2 * pow(s, 2) * u * MUs * UU +
                    2 * pow(s, 2) * u * MUs * UU * AXG +
                    2 * pow(s, 2) * u * MXs * SS -
                    2 * pow(s, 2) * u * MXs * SS * AXG -
                    2 * pow(s, 2) * u * MXs * UU +
                    2 * pow(s, 2) * u * MXs * UU * AXG -
                    2 * pow(s, 2) * pow(u, 2) * SS +
                    2 * pow(s, 2) * pow(u, 2) * SS * AXG +
                    2 * pow(s, 2) * pow(u, 2) * UU -
                    2 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC3 *
               (-2 * MXs * pow(MUs, 2) * SS - 2 * pow(MXs, 2) * MUs * SS +
                4 * pow(MX, 4) * MUs * SS - 2 * pow(MX, 4) * MUs * SS * AXG +
                4 * pow(MX, 4) * MUs * UU * AXG +
                4 * pow(MX, 4) * MXs * UU * AXG - 4 * pow(MX, 6) * UU * AXG +
                2 * u * pow(MUs, 2) * SS + 4 * u * MXs * MUs * SS * AXG -
                8 * u * MXs * MUs * UU * AXG + 2 * u * pow(MXs, 2) * SS -
                8 * u * pow(MXs, 2) * UU * AXG - 4 * u * pow(MX, 4) * SS +
                2 * u * pow(MX, 4) * SS * AXG + 4 * u * pow(MX, 4) * UU * AXG -
                2 * pow(u, 2) * MUs * SS - 2 * pow(u, 2) * MUs * SS * AXG +
                4 * pow(u, 2) * MUs * UU * AXG + 2 * pow(u, 2) * MXs * SS -
                4 * pow(u, 2) * MXs * SS * AXG +
                8 * pow(u, 2) * MXs * UU * AXG + 2 * pow(u, 3) * SS * AXG -
                4 * pow(u, 3) * UU * AXG - 2 * s * MXs * MUs * SS +
                2 * s * MXs * MUs * SS * AXG + 3 * s * MXs * MUs * UU -
                3 * s * MXs * MUs * UU * AXG + 3 * s * pow(MXs, 2) * UU -
                3 * s * pow(MXs, 2) * UU * AXG - 3 * s * pow(MX, 4) * UU +
                3 * s * pow(MX, 4) * UU * AXG - 2 * s * u * MUs * SS * AXG -
                3 * s * u * MUs * UU + 3 * s * u * MUs * UU * AXG +
                2 * s * u * MXs * SS - 2 * s * u * MXs * SS * AXG +
                s * u * MXs * UU + 3 * s * u * MXs * UU * AXG +
                2 * s * pow(u, 2) * SS * AXG - s * pow(u, 2) * UU -
                3 * s * pow(u, 2) * UU * AXG) +
           syFC4 *
               (4 * pow(MX, 4) * pow(MUs, 2) * SS -
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * pow(MX, 4) * MXs * MUs * SS -
                8 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 6) * MUs * SS + 4 * pow(MX, 6) * MUs * SS * AXG -
                8 * pow(MX, 6) * MXs * UU * AXG + 8 * pow(MX, 8) * UU * AXG +
                4 * pow(MU, 4) * MXs * MUs * SS +
                4 * pow(MU, 4) * pow(MXs, 2) * SS -
                8 * pow(MU, 4) * pow(MX, 4) * SS +
                4 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                12 * u * MXs * pow(MUs, 2) * SS +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                12 * u * pow(MXs, 2) * MUs * SS +
                16 * u * pow(MXs, 2) * MUs * UU * AXG +
                20 * u * pow(MX, 4) * MUs * SS -
                16 * u * pow(MX, 4) * MUs * SS * AXG +
                16 * u * pow(MX, 4) * MUs * UU * AXG -
                4 * u * pow(MX, 4) * MXs * SS +
                24 * u * pow(MX, 4) * MXs * UU * AXG + 8 * u * pow(MX, 6) * SS -
                4 * u * pow(MX, 6) * SS * AXG - 16 * u * pow(MX, 6) * UU * AXG -
                4 * u * pow(MU, 4) * MUs * SS + 4 * u * pow(MU, 4) * MXs * SS -
                8 * u * pow(MU, 4) * MXs * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                20 * pow(u, 2) * MXs * MUs * SS * AXG -
                40 * pow(u, 2) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * SS -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                12 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                4 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                4 * pow(u, 3) * MUs * SS) +
           syFC4 *
               (-8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                12 * pow(u, 3) * MXs * SS * AXG +
                24 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG - 4 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * UU +
                8 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MUs * SS -
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                2 * s * pow(MX, 4) * MUs * UU +
                6 * s * pow(MX, 4) * MUs * UU * AXG -
                6 * s * pow(MX, 4) * MXs * UU +
                14 * s * pow(MX, 4) * MXs * UU * AXG + 6 * s * pow(MX, 6) * UU -
                14 * s * pow(MX, 6) * UU * AXG + 6 * s * pow(MU, 4) * MXs * SS -
                6 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS + 8 * s * u * pow(MUs, 2) * UU -
                8 * s * u * pow(MUs, 2) * UU * AXG -
                10 * s * u * MXs * MUs * SS +
                22 * s * u * MXs * MUs * SS * AXG +
                12 * s * u * MXs * MUs * UU -
                36 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS + 12 * s * u * pow(MXs, 2) * UU -
                28 * s * u * pow(MXs, 2) * UU * AXG -
                12 * s * u * pow(MX, 4) * SS +
                8 * s * u * pow(MX, 4) * SS * AXG -
                14 * s * u * pow(MX, 4) * UU +
                14 * s * u * pow(MX, 4) * UU * AXG -
                2 * s * u * pow(MU, 4) * SS +
                6 * s * u * pow(MU, 4) * SS * AXG -
                2 * s * pow(u, 2) * MUs * SS -
                14 * s * pow(u, 2) * MUs * SS * AXG -
                6 * s * pow(u, 2) * MUs * UU) +
           syFC4 *
               (22 * s * pow(u, 2) * MUs * UU * AXG +
                8 * s * pow(u, 2) * MXs * SS -
                16 * s * pow(u, 2) * MXs * SS * AXG +
                4 * s * pow(u, 2) * MXs * UU +
                28 * s * pow(u, 2) * MXs * UU * AXG +
                8 * s * pow(u, 3) * SS * AXG - 2 * s * pow(u, 3) * UU -
                14 * s * pow(u, 3) * UU * AXG - 4 * pow(s, 2) * MXs * MUs * SS +
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                6 * pow(s, 2) * MXs * MUs * UU -
                6 * pow(s, 2) * MXs * MUs * UU * AXG +
                6 * pow(s, 2) * pow(MXs, 2) * UU -
                6 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                6 * pow(s, 2) * pow(MX, 4) * UU +
                6 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                4 * pow(s, 2) * u * MUs * SS * AXG -
                6 * pow(s, 2) * u * MUs * UU +
                6 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS -
                4 * pow(s, 2) * u * MXs * SS * AXG +
                2 * pow(s, 2) * u * MXs * UU +
                6 * pow(s, 2) * u * MXs * UU * AXG +
                4 * pow(s, 2) * pow(u, 2) * SS * AXG -
                2 * pow(s, 2) * pow(u, 2) * UU -
                6 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC5 * (-2 * s * MXs * MUs * SS + 2 * s * MXs * MUs * SS * AXG +
                    2 * s * MXs * MUs * UU - 2 * s * MXs * MUs * UU * AXG +
                    2 * s * pow(MXs, 2) * UU - 2 * s * pow(MXs, 2) * UU * AXG -
                    2 * s * pow(MX, 4) * UU + 2 * s * pow(MX, 4) * UU * AXG +
                    2 * s * u * MUs * SS - 2 * s * u * MUs * SS * AXG -
                    2 * s * u * MUs * UU + 2 * s * u * MUs * UU * AXG +
                    2 * s * u * MXs * SS - 2 * s * u * MXs * SS * AXG -
                    2 * s * u * MXs * UU + 2 * s * u * MXs * UU * AXG -
                    2 * s * pow(u, 2) * SS + 2 * s * pow(u, 2) * SS * AXG +
                    2 * s * pow(u, 2) * UU - 2 * s * pow(u, 2) * UU * AXG) +
           syFC6 *
               (2 * s * MXs * pow(MUs, 2) * SS -
                2 * s * MXs * pow(MUs, 2) * SS * AXG -
                6 * s * MXs * pow(MUs, 2) * UU +
                6 * s * MXs * pow(MUs, 2) * UU * AXG +
                2 * s * pow(MXs, 2) * MUs * SS -
                2 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                2 * s * pow(MXs, 3) * UU + 2 * s * pow(MXs, 3) * UU * AXG +
                2 * s * pow(MX, 4) * MUs * SS -
                2 * s * pow(MX, 4) * MUs * SS * AXG +
                4 * s * pow(MX, 4) * MUs * UU -
                4 * s * pow(MX, 4) * MUs * UU * AXG + 2 * s * pow(MX, 6) * UU -
                2 * s * pow(MX, 6) * UU * AXG + 4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG -
                2 * s * u * pow(MUs, 2) * SS +
                2 * s * u * pow(MUs, 2) * SS * AXG +
                6 * s * u * pow(MUs, 2) * UU -
                6 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS +
                16 * s * u * MXs * MUs * SS * AXG +
                16 * s * u * MXs * MUs * UU -
                16 * s * u * MXs * MUs * UU * AXG -
                2 * s * u * pow(MXs, 2) * SS +
                2 * s * u * pow(MXs, 2) * SS * AXG +
                10 * s * u * pow(MXs, 2) * UU -
                10 * s * u * pow(MXs, 2) * UU * AXG -
                2 * s * u * pow(MX, 4) * SS +
                2 * s * u * pow(MX, 4) * SS * AXG -
                6 * s * u * pow(MX, 4) * UU +
                6 * s * u * pow(MX, 4) * UU * AXG -
                4 * s * u * pow(MU, 4) * SS +
                4 * s * u * pow(MU, 4) * SS * AXG +
                12 * s * pow(u, 2) * MUs * SS -
                12 * s * pow(u, 2) * MUs * SS * AXG -
                12 * s * pow(u, 2) * MUs * UU +
                12 * s * pow(u, 2) * MUs * UU * AXG +
                10 * s * pow(u, 2) * MXs * SS -
                10 * s * pow(u, 2) * MXs * SS * AXG) +
           syFC6 *
               (-10 * s * pow(u, 2) * MXs * UU +
                10 * s * pow(u, 2) * MXs * UU * AXG - 6 * s * pow(u, 3) * SS +
                6 * s * pow(u, 3) * SS * AXG + 6 * s * pow(u, 3) * UU -
                6 * s * pow(u, 3) * UU * AXG - 4 * pow(s, 2) * MXs * MUs * SS +
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                4 * pow(s, 2) * MXs * MUs * UU -
                4 * pow(s, 2) * MXs * MUs * UU * AXG +
                4 * pow(s, 2) * pow(MXs, 2) * UU -
                4 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                4 * pow(s, 2) * pow(MX, 4) * UU +
                4 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                6 * pow(s, 2) * u * MUs * SS -
                6 * pow(s, 2) * u * MUs * SS * AXG -
                6 * pow(s, 2) * u * MUs * UU +
                6 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS -
                4 * pow(s, 2) * u * MXs * SS * AXG -
                4 * pow(s, 2) * u * MXs * UU +
                4 * pow(s, 2) * u * MXs * UU * AXG -
                6 * pow(s, 2) * pow(u, 2) * SS +
                6 * pow(s, 2) * pow(u, 2) * SS * AXG +
                6 * pow(s, 2) * pow(u, 2) * UU -
                6 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC7 *
               (4 * pow(MX, 4) * pow(MUs, 2) * SS -
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * pow(MX, 4) * MXs * MUs * SS -
                8 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 6) * MUs * SS + 4 * pow(MX, 6) * MUs * SS * AXG -
                8 * pow(MX, 6) * MXs * UU * AXG + 8 * pow(MX, 8) * UU * AXG +
                4 * pow(MU, 4) * MXs * MUs * SS +
                4 * pow(MU, 4) * pow(MXs, 2) * SS -
                8 * pow(MU, 4) * pow(MX, 4) * SS +
                4 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                12 * u * MXs * pow(MUs, 2) * SS +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                12 * u * pow(MXs, 2) * MUs * SS +
                16 * u * pow(MXs, 2) * MUs * UU * AXG +
                20 * u * pow(MX, 4) * MUs * SS -
                16 * u * pow(MX, 4) * MUs * SS * AXG +
                16 * u * pow(MX, 4) * MUs * UU * AXG -
                4 * u * pow(MX, 4) * MXs * SS +
                24 * u * pow(MX, 4) * MXs * UU * AXG + 8 * u * pow(MX, 6) * SS -
                4 * u * pow(MX, 6) * SS * AXG - 16 * u * pow(MX, 6) * UU * AXG -
                4 * u * pow(MU, 4) * MUs * SS + 4 * u * pow(MU, 4) * MXs * SS -
                8 * u * pow(MU, 4) * MXs * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                20 * pow(u, 2) * MXs * MUs * SS * AXG -
                40 * pow(u, 2) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * SS -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                12 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                4 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                4 * pow(u, 3) * MUs * SS) +
           syFC7 *
               (-8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                12 * pow(u, 3) * MXs * SS * AXG +
                24 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG - 4 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * UU +
                8 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG +
                14 * s * pow(MX, 4) * MUs * SS -
                10 * s * pow(MX, 4) * MUs * SS * AXG -
                4 * s * pow(MX, 4) * MUs * UU +
                12 * s * pow(MX, 4) * MUs * UU * AXG -
                8 * s * pow(MX, 4) * MXs * UU +
                16 * s * pow(MX, 4) * MXs * UU * AXG + 4 * s * pow(MX, 6) * UU -
                12 * s * pow(MX, 6) * UU * AXG + 6 * s * pow(MU, 4) * MXs * SS -
                6 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS + 8 * s * u * pow(MUs, 2) * UU -
                8 * s * u * pow(MUs, 2) * UU * AXG -
                14 * s * u * MXs * MUs * SS +
                26 * s * u * MXs * MUs * SS * AXG +
                24 * s * u * MXs * MUs * UU -
                48 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS + 16 * s * u * pow(MXs, 2) * UU -
                32 * s * u * pow(MXs, 2) * UU * AXG -
                14 * s * u * pow(MX, 4) * SS +
                10 * s * u * pow(MX, 4) * SS * AXG -
                4 * s * u * pow(MX, 4) * UU +
                4 * s * u * pow(MX, 4) * UU * AXG -
                2 * s * u * pow(MU, 4) * SS +
                6 * s * u * pow(MU, 4) * SS * AXG -
                16 * s * pow(u, 2) * MUs * SS * AXG -
                12 * s * pow(u, 2) * MUs * UU +
                28 * s * pow(u, 2) * MUs * UU * AXG) +
           syFC7 *
               (12 * s * pow(u, 2) * MXs * SS -
                20 * s * pow(u, 2) * MXs * SS * AXG -
                12 * s * pow(u, 2) * MXs * UU +
                44 * s * pow(u, 2) * MXs * UU * AXG - 2 * s * pow(u, 3) * SS +
                10 * s * pow(u, 3) * SS * AXG + 4 * s * pow(u, 3) * UU -
                20 * s * pow(u, 3) * UU * AXG - 6 * pow(s, 2) * MXs * MUs * SS +
                6 * pow(s, 2) * MXs * MUs * SS * AXG +
                10 * pow(s, 2) * MXs * MUs * UU -
                10 * pow(s, 2) * MXs * MUs * UU * AXG +
                8 * pow(s, 2) * pow(MXs, 2) * UU -
                8 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                2 * pow(s, 2) * pow(MX, 4) * UU +
                2 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                2 * pow(s, 2) * u * MUs * SS -
                6 * pow(s, 2) * u * MUs * SS * AXG -
                10 * pow(s, 2) * u * MUs * UU +
                10 * pow(s, 2) * u * MUs * UU * AXG +
                6 * pow(s, 2) * u * MXs * SS -
                6 * pow(s, 2) * u * MXs * SS * AXG -
                12 * pow(s, 2) * u * MXs * UU +
                20 * pow(s, 2) * u * MXs * UU * AXG -
                2 * pow(s, 2) * pow(u, 2) * SS +
                6 * pow(s, 2) * pow(u, 2) * SS * AXG +
                6 * pow(s, 2) * pow(u, 2) * UU -
                14 * pow(s, 2) * pow(u, 2) * UU * AXG -
                2 * pow(s, 3) * MXs * UU + 2 * pow(s, 3) * MXs * UU * AXG +
                2 * pow(s, 3) * u * UU - 2 * pow(s, 3) * u * UU * AXG) +
           syFC8 *
               (8 * MXs * pow(MUs, 2) * SS + 8 * pow(MXs, 2) * MUs * SS -
                16 * pow(MX, 4) * MUs * SS + 8 * pow(MX, 4) * MUs * SS * AXG -
                16 * pow(MX, 4) * MUs * UU * AXG -
                16 * pow(MX, 4) * MXs * UU * AXG + 16 * pow(MX, 6) * UU * AXG -
                8 * u * pow(MUs, 2) * SS - 16 * u * MXs * MUs * SS * AXG +
                32 * u * MXs * MUs * UU * AXG - 8 * u * pow(MXs, 2) * SS +
                32 * u * pow(MXs, 2) * UU * AXG + 16 * u * pow(MX, 4) * SS -
                8 * u * pow(MX, 4) * SS * AXG - 16 * u * pow(MX, 4) * UU * AXG +
                8 * pow(u, 2) * MUs * SS + 8 * pow(u, 2) * MUs * SS * AXG -
                16 * pow(u, 2) * MUs * UU * AXG - 8 * pow(u, 2) * MXs * SS +
                16 * pow(u, 2) * MXs * SS * AXG -
                32 * pow(u, 2) * MXs * UU * AXG - 8 * pow(u, 3) * SS * AXG +
                16 * pow(u, 3) * UU * AXG - 8 * s * pow(MUs, 2) * SS -
                8 * s * MXs * MUs * SS * AXG - 4 * s * MXs * MUs * UU +
                12 * s * MXs * MUs * UU * AXG - 4 * s * pow(MXs, 2) * UU +
                12 * s * pow(MXs, 2) * UU * AXG + 4 * s * pow(MX, 4) * UU -
                12 * s * pow(MX, 4) * UU * AXG + 16 * s * u * MUs * SS +
                8 * s * u * MUs * SS * AXG + 4 * s * u * MUs * UU -
                12 * s * u * MUs * UU * AXG + 8 * s * u * MXs * SS * AXG -
                4 * s * u * MXs * UU - 20 * s * u * MXs * UU * AXG -
                8 * s * pow(u, 2) * SS - 8 * s * pow(u, 2) * SS * AXG +
                4 * s * pow(u, 2) * UU + 20 * s * pow(u, 2) * UU * AXG +
                8 * pow(s, 2) * MUs * SS - 8 * pow(s, 2) * u * SS +
                4 * pow(s, 2) * u * UU + 4 * pow(s, 2) * u * UU * AXG) +
           syFC9 *
               (-8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                8 * pow(MX, 4) * MXs * MUs * UU * AXG +
                8 * pow(MX, 6) * MUs * UU * AXG +
                4 * pow(MU, 4) * MXs * MUs * SS +
                4 * pow(MU, 4) * pow(MXs, 2) * SS -
                8 * pow(MU, 4) * pow(MX, 4) * SS +
                4 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                4 * u * MXs * pow(MUs, 2) * SS +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MXs, 2) * MUs * SS +
                16 * u * pow(MXs, 2) * MUs * UU * AXG +
                8 * u * pow(MX, 4) * MUs * SS -
                4 * u * pow(MX, 4) * MUs * SS * AXG -
                8 * u * pow(MX, 4) * MUs * UU * AXG -
                4 * u * pow(MU, 4) * MUs * SS + 4 * u * pow(MU, 4) * MXs * SS -
                8 * u * pow(MU, 4) * MXs * SS * AXG +
                4 * pow(u, 2) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                4 * pow(u, 2) * MXs * MUs * SS +
                8 * pow(u, 2) * MXs * MUs * SS * AXG -
                16 * pow(u, 2) * MXs * MUs * UU * AXG +
                4 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                4 * pow(u, 3) * MUs * SS * AXG +
                8 * pow(u, 3) * MUs * UU * AXG -
                6 * s * MXs * pow(MUs, 2) * UU +
                6 * s * MXs * pow(MUs, 2) * UU * AXG -
                6 * s * pow(MXs, 2) * MUs * UU +
                6 * s * pow(MXs, 2) * MUs * UU * AXG +
                6 * s * pow(MX, 4) * MUs * UU -
                6 * s * pow(MX, 4) * MUs * UU * AXG +
                4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG +
                6 * s * u * pow(MUs, 2) * UU -
                6 * s * u * pow(MUs, 2) * UU * AXG -
                4 * s * u * MXs * MUs * SS + 4 * s * u * MXs * MUs * SS * AXG -
                2 * s * u * MXs * MUs * UU - 6 * s * u * MXs * MUs * UU * AXG) +
           syFC9 * (4 * s * u * pow(MU, 4) * SS * AXG -
                    4 * s * pow(u, 2) * MUs * SS * AXG +
                    2 * s * pow(u, 2) * MUs * UU +
                    6 * s * pow(u, 2) * MUs * UU * AXG) +
           syFC10 *
               (4 * MXs * pow(MUs, 3) * SS +
                8 * pow(MXs, 2) * pow(MUs, 2) * SS +
                4 * pow(MXs, 3) * MUs * SS -
                12 * pow(MX, 4) * pow(MUs, 2) * SS +
                4 * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                12 * pow(MX, 4) * MXs * MUs * SS +
                4 * pow(MX, 4) * MXs * MUs * SS * AXG -
                16 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 4) * pow(MXs, 2) * UU * AXG +
                8 * pow(MX, 6) * MUs * SS - 4 * pow(MX, 6) * MUs * SS * AXG +
                16 * pow(MX, 6) * MUs * UU * AXG +
                16 * pow(MX, 6) * MXs * UU * AXG - 8 * pow(MX, 8) * UU * AXG -
                4 * u * pow(MUs, 3) * SS - 4 * u * MXs * pow(MUs, 2) * SS -
                8 * u * MXs * pow(MUs, 2) * SS * AXG +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MXs, 2) * MUs * SS -
                8 * u * pow(MXs, 2) * MUs * SS * AXG +
                32 * u * pow(MXs, 2) * MUs * UU * AXG -
                4 * u * pow(MXs, 3) * SS + 16 * u * pow(MXs, 3) * UU * AXG +
                12 * u * pow(MX, 4) * MUs * SS -
                16 * u * pow(MX, 4) * MUs * UU * AXG +
                12 * u * pow(MX, 4) * MXs * SS -
                4 * u * pow(MX, 4) * MXs * SS * AXG -
                16 * u * pow(MX, 4) * MXs * UU * AXG - 8 * u * pow(MX, 6) * SS +
                4 * u * pow(MX, 6) * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS +
                4 * pow(u, 2) * pow(MUs, 2) * SS * AXG -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                4 * pow(u, 2) * MXs * MUs * SS +
                16 * pow(u, 2) * MXs * MUs * SS * AXG -
                32 * pow(u, 2) * MXs * MUs * UU * AXG -
                4 * pow(u, 2) * pow(MXs, 2) * SS) +
           syFC10 *
               (8 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                4 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                16 * pow(u, 2) * pow(MX, 4) * UU * AXG -
                4 * pow(u, 3) * MUs * SS - 8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                8 * pow(u, 3) * MXs * SS * AXG +
                16 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                10 * s * MXs * pow(MUs, 2) * UU +
                10 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                16 * s * pow(MXs, 2) * MUs * UU * AXG -
                6 * s * pow(MXs, 3) * UU + 6 * s * pow(MXs, 3) * UU * AXG +
                4 * s * pow(MX, 4) * MUs * SS * AXG +
                16 * s * pow(MX, 4) * MUs * UU -
                16 * s * pow(MX, 4) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MXs * UU -
                12 * s * pow(MX, 4) * MXs * UU * AXG - 6 * s * pow(MX, 6) * UU +
                6 * s * pow(MX, 6) * UU * AXG + 4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                10 * s * u * pow(MUs, 2) * UU -
                10 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 16 * s * u * MXs * MUs * SS * AXG +
                8 * s * u * MXs * MUs * UU - 24 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS * AXG -
                2 * s * u * pow(MXs, 2) * UU -
                14 * s * u * pow(MXs, 2) * UU * AXG) +
           syFC10 * (-4 * s * u * pow(MX, 4) * SS * AXG +
                     2 * s * u * pow(MX, 4) * UU +
                     14 * s * u * pow(MX, 4) * UU * AXG -
                     4 * s * u * pow(MU, 4) * SS +
                     4 * s * u * pow(MU, 4) * SS * AXG -
                     16 * s * pow(u, 2) * MUs * SS * AXG -
                     8 * s * pow(u, 2) * MUs * UU +
                     24 * s * pow(u, 2) * MUs * UU * AXG +
                     4 * s * pow(u, 2) * MXs * SS -
                     8 * s * pow(u, 2) * MXs * SS * AXG +
                     2 * s * pow(u, 2) * MXs * UU +
                     14 * s * pow(u, 2) * MXs * UU * AXG +
                     8 * s * pow(u, 3) * SS * AXG - 2 * s * pow(u, 3) * UU -
                     14 * s * pow(u, 3) * UU * AXG -
                     4 * pow(s, 2) * u * MUs * SS * AXG -
                     6 * pow(s, 2) * u * MUs * UU +
                     6 * pow(s, 2) * u * MUs * UU * AXG +
                     4 * pow(s, 2) * pow(u, 2) * SS * AXG -
                     2 * pow(s, 2) * pow(u, 2) * UU -
                     6 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC11 *
               (4 * pow(MX, 4) * pow(MUs, 2) * SS -
                16 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * pow(MX, 4) * MXs * MUs * SS -
                16 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 6) * MUs * SS + 4 * pow(MX, 6) * MUs * SS * AXG +
                8 * pow(MX, 6) * MUs * UU * AXG -
                8 * pow(MX, 6) * MXs * UU * AXG + 8 * pow(MX, 8) * UU * AXG +
                8 * pow(MU, 4) * MXs * MUs * SS +
                8 * pow(MU, 4) * pow(MXs, 2) * SS -
                16 * pow(MU, 4) * pow(MX, 4) * SS +
                8 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                16 * u * MXs * pow(MUs, 2) * SS +
                32 * u * MXs * pow(MUs, 2) * UU * AXG -
                16 * u * pow(MXs, 2) * MUs * SS +
                32 * u * pow(MXs, 2) * MUs * UU * AXG +
                28 * u * pow(MX, 4) * MUs * SS -
                20 * u * pow(MX, 4) * MUs * SS * AXG +
                8 * u * pow(MX, 4) * MUs * UU * AXG -
                4 * u * pow(MX, 4) * MXs * SS +
                24 * u * pow(MX, 4) * MXs * UU * AXG + 8 * u * pow(MX, 6) * SS -
                4 * u * pow(MX, 6) * SS * AXG - 16 * u * pow(MX, 6) * UU * AXG -
                8 * u * pow(MU, 4) * MUs * SS + 8 * u * pow(MU, 4) * MXs * SS -
                16 * u * pow(MU, 4) * MXs * SS * AXG +
                12 * pow(u, 2) * pow(MUs, 2) * SS -
                16 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                4 * pow(u, 2) * MXs * MUs * SS +
                28 * pow(u, 2) * MXs * MUs * SS * AXG -
                56 * pow(u, 2) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * SS -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                12 * pow(u, 2) * pow(MX, 4) * SS * AXG) +
           syFC11 *
               (8 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                4 * pow(u, 3) * MUs * SS - 12 * pow(u, 3) * MUs * SS * AXG +
                24 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                12 * pow(u, 3) * MXs * SS * AXG +
                24 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG - 4 * s * MXs * pow(MUs, 2) * SS -
                12 * s * MXs * pow(MUs, 2) * UU +
                12 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS -
                12 * s * pow(MXs, 2) * MUs * UU +
                12 * s * pow(MXs, 2) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MUs * SS -
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                2 * s * pow(MX, 4) * MUs * UU +
                6 * s * pow(MX, 4) * MUs * UU * AXG -
                6 * s * pow(MX, 4) * MXs * UU +
                14 * s * pow(MX, 4) * MXs * UU * AXG + 6 * s * pow(MX, 6) * UU -
                14 * s * pow(MX, 6) * UU * AXG + 8 * s * pow(MU, 4) * MXs * SS -
                8 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS + 12 * s * u * pow(MUs, 2) * UU -
                12 * s * u * pow(MUs, 2) * UU * AXG -
                12 * s * u * MXs * MUs * SS +
                24 * s * u * MXs * MUs * SS * AXG +
                16 * s * u * MXs * MUs * UU -
                48 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS + 12 * s * u * pow(MXs, 2) * UU -
                28 * s * u * pow(MXs, 2) * UU * AXG -
                12 * s * u * pow(MX, 4) * SS +
                8 * s * u * pow(MX, 4) * SS * AXG -
                14 * s * u * pow(MX, 4) * UU +
                14 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * u * pow(MU, 4) * SS * AXG -
                4 * s * pow(u, 2) * MUs * SS -
                16 * s * pow(u, 2) * MUs * SS * AXG) +
           syFC11 *
               (-6 * s * pow(u, 2) * MUs * UU +
                30 * s * pow(u, 2) * MUs * UU * AXG +
                8 * s * pow(u, 2) * MXs * SS -
                16 * s * pow(u, 2) * MXs * SS * AXG +
                4 * s * pow(u, 2) * MXs * UU +
                28 * s * pow(u, 2) * MXs * UU * AXG +
                8 * s * pow(u, 3) * SS * AXG - 2 * s * pow(u, 3) * UU -
                14 * s * pow(u, 3) * UU * AXG - 4 * pow(s, 2) * MXs * MUs * SS +
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                8 * pow(s, 2) * MXs * MUs * UU -
                8 * pow(s, 2) * MXs * MUs * UU * AXG +
                6 * pow(s, 2) * pow(MXs, 2) * UU -
                6 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                6 * pow(s, 2) * pow(MX, 4) * UU +
                6 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                4 * pow(s, 2) * u * MUs * SS * AXG -
                8 * pow(s, 2) * u * MUs * UU +
                8 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS -
                4 * pow(s, 2) * u * MXs * SS * AXG +
                2 * pow(s, 2) * u * MXs * UU +
                6 * pow(s, 2) * u * MXs * UU * AXG +
                4 * pow(s, 2) * pow(u, 2) * SS * AXG -
                2 * pow(s, 2) * pow(u, 2) * UU -
                6 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC12 *
               (4 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                4 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                8 * s * pow(MXs, 2) * MUs * UU +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                4 * s * pow(MX, 4) * MUs * SS +
                4 * s * pow(MX, 4) * MUs * SS * AXG +
                8 * s * pow(MX, 4) * MUs * UU -
                8 * s * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(MX, 4) * MXs * UU -
                8 * s * pow(MX, 4) * MXs * UU * AXG - 4 * s * pow(MX, 6) * UU +
                4 * s * pow(MX, 6) * UU * AXG - 4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * SS * AXG +
                8 * s * u * MXs * MUs * UU - 8 * s * u * MXs * MUs * UU * AXG -
                4 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG +
                4 * s * u * pow(MXs, 2) * UU -
                4 * s * u * pow(MXs, 2) * UU * AXG +
                4 * s * u * pow(MX, 4) * SS -
                4 * s * u * pow(MX, 4) * SS * AXG -
                4 * s * u * pow(MX, 4) * UU +
                4 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * pow(u, 2) * MUs * SS -
                8 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU +
                8 * s * pow(u, 2) * MUs * UU * AXG +
                4 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU +
                4 * s * pow(u, 2) * MXs * UU * AXG) +
           syFC12 * (-4 * s * pow(u, 3) * SS + 4 * s * pow(u, 3) * SS * AXG +
                     4 * s * pow(u, 3) * UU - 4 * s * pow(u, 3) * UU * AXG +
                     4 * pow(s, 2) * u * MUs * SS -
                     4 * pow(s, 2) * u * MUs * SS * AXG -
                     4 * pow(s, 2) * u * MUs * UU +
                     4 * pow(s, 2) * u * MUs * UU * AXG -
                     4 * pow(s, 2) * pow(u, 2) * SS +
                     4 * pow(s, 2) * pow(u, 2) * SS * AXG +
                     4 * pow(s, 2) * pow(u, 2) * UU -
                     4 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC13 *
               (4 * MXs * pow(MUs, 3) * SS +
                8 * pow(MXs, 2) * pow(MUs, 2) * SS +
                4 * pow(MXs, 3) * MUs * SS -
                12 * pow(MX, 4) * pow(MUs, 2) * SS +
                4 * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                12 * pow(MX, 4) * MXs * MUs * SS +
                4 * pow(MX, 4) * MXs * MUs * SS * AXG -
                16 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 4) * pow(MXs, 2) * UU * AXG +
                8 * pow(MX, 6) * MUs * SS - 4 * pow(MX, 6) * MUs * SS * AXG +
                16 * pow(MX, 6) * MUs * UU * AXG +
                16 * pow(MX, 6) * MXs * UU * AXG - 8 * pow(MX, 8) * UU * AXG -
                4 * u * pow(MUs, 3) * SS - 4 * u * MXs * pow(MUs, 2) * SS -
                8 * u * MXs * pow(MUs, 2) * SS * AXG +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                4 * u * pow(MXs, 2) * MUs * SS -
                8 * u * pow(MXs, 2) * MUs * SS * AXG +
                32 * u * pow(MXs, 2) * MUs * UU * AXG -
                4 * u * pow(MXs, 3) * SS + 16 * u * pow(MXs, 3) * UU * AXG +
                12 * u * pow(MX, 4) * MUs * SS -
                16 * u * pow(MX, 4) * MUs * UU * AXG +
                12 * u * pow(MX, 4) * MXs * SS -
                4 * u * pow(MX, 4) * MXs * SS * AXG -
                16 * u * pow(MX, 4) * MXs * UU * AXG - 8 * u * pow(MX, 6) * SS +
                4 * u * pow(MX, 6) * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS +
                4 * pow(u, 2) * pow(MUs, 2) * SS * AXG -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                4 * pow(u, 2) * MXs * MUs * SS +
                16 * pow(u, 2) * MXs * MUs * SS * AXG -
                32 * pow(u, 2) * MXs * MUs * UU * AXG -
                4 * pow(u, 2) * pow(MXs, 2) * SS) +
           syFC13 *
               (8 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                4 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                16 * pow(u, 2) * pow(MX, 4) * UU * AXG -
                4 * pow(u, 3) * MUs * SS - 8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                8 * pow(u, 3) * MXs * SS * AXG +
                16 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG -
                4 * s * MXs * pow(MUs, 2) * SS * AXG -
                10 * s * MXs * pow(MUs, 2) * UU +
                10 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                16 * s * pow(MXs, 2) * MUs * UU * AXG -
                6 * s * pow(MXs, 3) * UU + 6 * s * pow(MXs, 3) * UU * AXG +
                4 * s * pow(MX, 4) * MUs * SS + 8 * s * pow(MX, 4) * MUs * UU -
                8 * s * pow(MX, 4) * MUs * UU * AXG +
                4 * s * pow(MX, 4) * MXs * UU -
                4 * s * pow(MX, 4) * MXs * UU * AXG + 2 * s * pow(MX, 6) * UU -
                2 * s * pow(MX, 6) * UU * AXG + 4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG +
                10 * s * u * pow(MUs, 2) * UU -
                10 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS +
                24 * s * u * MXs * MUs * SS * AXG +
                24 * s * u * MXs * MUs * UU -
                40 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS * AXG +
                14 * s * u * pow(MXs, 2) * UU -
                30 * s * u * pow(MXs, 2) * UU * AXG) +
           syFC13 *
               (-4 * s * u * pow(MX, 4) * SS - 6 * s * u * pow(MX, 4) * UU +
                22 * s * u * pow(MX, 4) * UU * AXG -
                4 * s * u * pow(MU, 4) * SS +
                4 * s * u * pow(MU, 4) * SS * AXG +
                4 * s * pow(u, 2) * MUs * SS -
                20 * s * pow(u, 2) * MUs * SS * AXG -
                16 * s * pow(u, 2) * MUs * UU +
                32 * s * pow(u, 2) * MUs * UU * AXG +
                12 * s * pow(u, 2) * MXs * SS -
                16 * s * pow(u, 2) * MXs * SS * AXG -
                14 * s * pow(u, 2) * MXs * UU +
                30 * s * pow(u, 2) * MXs * UU * AXG - 4 * s * pow(u, 3) * SS +
                12 * s * pow(u, 3) * SS * AXG + 6 * s * pow(u, 3) * UU -
                22 * s * pow(u, 3) * UU * AXG - 4 * pow(s, 2) * MXs * MUs * SS +
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                6 * pow(s, 2) * MXs * MUs * UU -
                6 * pow(s, 2) * MXs * MUs * UU * AXG +
                6 * pow(s, 2) * pow(MXs, 2) * UU -
                6 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                6 * pow(s, 2) * pow(MX, 4) * UU +
                6 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(s, 2) * u * MUs * SS -
                8 * pow(s, 2) * u * MUs * SS * AXG -
                12 * pow(s, 2) * u * MUs * UU +
                12 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS -
                4 * pow(s, 2) * u * MXs * SS * AXG -
                10 * pow(s, 2) * u * MXs * UU +
                10 * pow(s, 2) * u * MXs * UU * AXG -
                4 * pow(s, 2) * pow(u, 2) * SS +
                8 * pow(s, 2) * pow(u, 2) * SS * AXG +
                8 * pow(s, 2) * pow(u, 2) * UU -
                16 * pow(s, 2) * pow(u, 2) * UU * AXG + 2 * pow(s, 3) * u * UU -
                2 * pow(s, 3) * u * UU * AXG) +
           syFC14 *
               (4 * pow(MX, 4) * pow(MUs, 2) * SS -
                8 * pow(MX, 4) * pow(MUs, 2) * UU * AXG +
                4 * pow(MX, 4) * MXs * MUs * SS -
                8 * pow(MX, 4) * MXs * MUs * UU * AXG -
                8 * pow(MX, 6) * MUs * SS + 4 * pow(MX, 6) * MUs * SS * AXG -
                8 * pow(MX, 6) * MXs * UU * AXG + 8 * pow(MX, 8) * UU * AXG +
                4 * pow(MU, 4) * MXs * MUs * SS +
                4 * pow(MU, 4) * pow(MXs, 2) * SS -
                8 * pow(MU, 4) * pow(MX, 4) * SS +
                4 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                12 * u * MXs * pow(MUs, 2) * SS +
                16 * u * MXs * pow(MUs, 2) * UU * AXG -
                12 * u * pow(MXs, 2) * MUs * SS +
                16 * u * pow(MXs, 2) * MUs * UU * AXG +
                20 * u * pow(MX, 4) * MUs * SS -
                16 * u * pow(MX, 4) * MUs * SS * AXG +
                16 * u * pow(MX, 4) * MUs * UU * AXG -
                4 * u * pow(MX, 4) * MXs * SS +
                24 * u * pow(MX, 4) * MXs * UU * AXG + 8 * u * pow(MX, 6) * SS -
                4 * u * pow(MX, 6) * SS * AXG - 16 * u * pow(MX, 6) * UU * AXG -
                4 * u * pow(MU, 4) * MUs * SS + 4 * u * pow(MU, 4) * MXs * SS -
                8 * u * pow(MU, 4) * MXs * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS -
                8 * pow(u, 2) * pow(MUs, 2) * UU * AXG +
                20 * pow(u, 2) * MXs * MUs * SS * AXG -
                40 * pow(u, 2) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MXs, 2) * SS -
                24 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                16 * pow(u, 2) * pow(MX, 4) * SS +
                12 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                4 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                4 * pow(u, 3) * MUs * SS) +
           syFC14 *
               (-8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG + 4 * pow(u, 3) * MXs * SS -
                12 * pow(u, 3) * MXs * SS * AXG +
                24 * pow(u, 3) * MXs * UU * AXG + 4 * pow(u, 4) * SS * AXG -
                8 * pow(u, 4) * UU * AXG - 4 * s * MXs * pow(MUs, 2) * SS -
                6 * s * MXs * pow(MUs, 2) * UU +
                6 * s * MXs * pow(MUs, 2) * UU * AXG -
                4 * s * pow(MXs, 2) * MUs * SS -
                6 * s * pow(MXs, 2) * MUs * UU +
                6 * s * pow(MXs, 2) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MUs * SS -
                8 * s * pow(MX, 4) * MUs * SS * AXG -
                4 * s * pow(MX, 4) * MUs * UU +
                12 * s * pow(MX, 4) * MUs * UU * AXG -
                6 * s * pow(MX, 4) * MXs * UU +
                14 * s * pow(MX, 4) * MXs * UU * AXG + 2 * s * pow(MX, 6) * UU -
                10 * s * pow(MX, 6) * UU * AXG + 4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG +
                4 * s * u * pow(MUs, 2) * SS + 6 * s * u * pow(MUs, 2) * UU -
                6 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 20 * s * u * MXs * MUs * SS * AXG +
                18 * s * u * MXs * MUs * UU -
                42 * s * u * MXs * MUs * UU * AXG +
                4 * s * u * pow(MXs, 2) * SS + 12 * s * u * pow(MXs, 2) * UU -
                28 * s * u * pow(MXs, 2) * UU * AXG -
                12 * s * u * pow(MX, 4) * SS +
                8 * s * u * pow(MX, 4) * SS * AXG -
                2 * s * u * pow(MX, 4) * UU +
                2 * s * u * pow(MX, 4) * UU * AXG +
                4 * s * u * pow(MU, 4) * SS * AXG -
                4 * s * pow(u, 2) * MUs * SS -
                12 * s * pow(u, 2) * MUs * SS * AXG -
                8 * s * pow(u, 2) * MUs * UU +
                24 * s * pow(u, 2) * MUs * UU * AXG) +
           syFC14 *
               (8 * s * pow(u, 2) * MXs * SS -
                16 * s * pow(u, 2) * MXs * SS * AXG -
                8 * s * pow(u, 2) * MXs * UU +
                40 * s * pow(u, 2) * MXs * UU * AXG +
                8 * s * pow(u, 3) * SS * AXG + 2 * s * pow(u, 3) * UU -
                18 * s * pow(u, 3) * UU * AXG - 4 * pow(s, 2) * MXs * MUs * SS +
                4 * pow(s, 2) * MXs * MUs * SS * AXG +
                8 * pow(s, 2) * MXs * MUs * UU -
                8 * pow(s, 2) * MXs * MUs * UU * AXG +
                6 * pow(s, 2) * pow(MXs, 2) * UU -
                6 * pow(s, 2) * pow(MXs, 2) * UU * AXG -
                4 * pow(s, 2) * u * MUs * SS * AXG -
                8 * pow(s, 2) * u * MUs * UU +
                8 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * SS -
                4 * pow(s, 2) * u * MXs * SS * AXG -
                10 * pow(s, 2) * u * MXs * UU +
                18 * pow(s, 2) * u * MXs * UU * AXG +
                4 * pow(s, 2) * pow(u, 2) * SS * AXG +
                4 * pow(s, 2) * pow(u, 2) * UU -
                12 * pow(s, 2) * pow(u, 2) * UU * AXG -
                2 * pow(s, 3) * MXs * UU + 2 * pow(s, 3) * MXs * UU * AXG +
                2 * pow(s, 3) * u * UU - 2 * pow(s, 3) * u * UU * AXG)));

  return ret.real();
}

ComplexType ME_us_box_qggQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2,
                           double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C0(0, 0, s, 0, 0, 0)
#define syFC2 C0(MUs, 0, u, MQs, 0, 0)
#define syFC3 C0(MUs, s, MXs, MQs, 0, 0)
#define syFC4 D0(MUs, 0, 0, MXs, u, s, MQs, 0, 0, 0)
#define syFC5 C1(0, s, 0, 0, 0, 0)
#define syFC6 C1(0, u, MXs, 0, 0, MQs)
#define syFC7 C1(MUs, s, MXs, MQs, 0, 0)
#define syFC8 C1(MUs, u, 0, 0, MQs, 0)
#define syFC9 D1(0, s, MXs, u, 0, MUs, 0, 0, 0, MQs)
#define syFC10 D1(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC11 C2(0, u, MXs, 0, 0, MQs)
#define syFC12 C2(MUs, s, MXs, MQs, 0, 0)
#define syFC13 C2(MUs, u, 0, 0, MQs, 0)
#define syFC14 D2(0, s, MXs, u, 0, MUs, 0, 0, 0, MQs)
#define syFC15 D2(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC16 D3(0, s, MXs, u, 0, MUs, 0, 0, 0, MQs)
#define syFC17 D3(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC18 D00(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC19 D11(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC20 D12(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC21 D13(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC22 D22(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC23 D23(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)
#define syFC24 D33(MUs, u, 0, s, 0, MXs, 0, MQs, 0, 0)

  _EPS0(
      ret,
      (-1 + Nc) * (+1 + Nc) * (+TR) * (+TR) *
          (+Denom(1536 * pow(Pi, 2) * s * MUs - 1536 * pow(Pi, 2) * s * u)) *
          (+Nc) * (+R * Lp + L * Rp) * (+gs) * (+gs) * (+gs) * (+gs) *
          (+syFC1 *
               (-4 * s * pow(MXs, 2) * UU + 4 * s * pow(MXs, 2) * UU * AXG +
                8 * s * u * MXs * UU - 8 * s * u * MXs * UU * AXG -
                4 * s * pow(u, 2) * UU + 4 * s * pow(u, 2) * UU * AXG +
                2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG -
                2 * pow(s, 2) * u * UU + 2 * pow(s, 2) * u * UU * AXG) +
           syFC2 *
               (-4 * MXs * pow(MUs, 2) * SS - 4 * pow(MXs, 2) * MUs * SS +
                8 * pow(MX, 4) * MUs * SS - 4 * pow(MX, 4) * MUs * SS * AXG +
                8 * pow(MX, 4) * MUs * UU * AXG +
                8 * pow(MX, 4) * MXs * UU * AXG - 8 * pow(MX, 6) * UU * AXG +
                4 * u * pow(MUs, 2) * SS + 8 * u * MXs * MUs * SS * AXG -
                16 * u * MXs * MUs * UU * AXG + 4 * u * pow(MXs, 2) * SS -
                16 * u * pow(MXs, 2) * UU * AXG - 8 * u * pow(MX, 4) * SS +
                4 * u * pow(MX, 4) * SS * AXG + 8 * u * pow(MX, 4) * UU * AXG -
                4 * pow(u, 2) * MUs * SS - 4 * pow(u, 2) * MUs * SS * AXG +
                8 * pow(u, 2) * MUs * UU * AXG + 4 * pow(u, 2) * MXs * SS -
                8 * pow(u, 2) * MXs * SS * AXG +
                16 * pow(u, 2) * MXs * UU * AXG + 4 * pow(u, 3) * SS * AXG -
                8 * pow(u, 3) * UU * AXG - 2 * s * MXs * MUs * SS +
                2 * s * MXs * MUs * SS * AXG + 4 * s * MXs * MUs * UU -
                4 * s * MXs * MUs * UU * AXG + 4 * s * pow(MXs, 2) * UU -
                4 * s * pow(MXs, 2) * UU * AXG - 4 * s * pow(MX, 4) * UU +
                4 * s * pow(MX, 4) * UU * AXG - 2 * s * u * MUs * SS -
                2 * s * u * MUs * SS * AXG - 4 * s * u * MUs * UU +
                4 * s * u * MUs * UU * AXG + 2 * s * u * MXs * SS -
                2 * s * u * MXs * SS * AXG + 4 * s * u * MXs * UU +
                4 * s * u * MXs * UU * AXG + 2 * s * pow(u, 2) * SS +
                2 * s * pow(u, 2) * SS * AXG - 4 * s * pow(u, 2) * UU -
                4 * s * pow(u, 2) * UU * AXG) +
           syFC3 * (4 * s * MXs * MUs * SS - 4 * s * pow(MX, 4) * UU +
                    4 * s * pow(MX, 4) * UU * AXG - 4 * s * u * MXs * SS +
                    4 * s * u * MXs * UU - 4 * s * u * MXs * UU * AXG +
                    2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG) +
           syFC4 *
               (12 * s * MXs * pow(MUs, 2) * SS -
                4 * s * MXs * pow(MUs, 2) * SS * AXG +
                4 * s * MXs * pow(MUs, 2) * UU +
                4 * s * MXs * pow(MUs, 2) * UU * AXG +
                12 * s * pow(MXs, 2) * MUs * SS -
                4 * s * pow(MXs, 2) * MUs * SS * AXG +
                8 * s * pow(MXs, 2) * MUs * UU * AXG -
                4 * s * pow(MXs, 3) * UU + 4 * s * pow(MXs, 3) * UU * AXG -
                4 * s * pow(MX, 4) * MQs * UU +
                4 * s * pow(MX, 4) * MQs * UU * AXG -
                24 * s * pow(MX, 4) * MUs * SS +
                16 * s * pow(MX, 4) * MUs * SS * AXG -
                24 * s * pow(MX, 4) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MXs * UU -
                28 * s * pow(MX, 4) * MXs * UU * AXG - 8 * s * pow(MX, 6) * UU +
                24 * s * pow(MX, 6) * UU * AXG - 8 * s * pow(MU, 4) * MUs * SS +
                4 * s * u * pow(MUs, 2) * SS +
                4 * s * u * pow(MUs, 2) * SS * AXG -
                4 * s * u * pow(MUs, 2) * UU -
                4 * s * u * pow(MUs, 2) * UU * AXG +
                8 * s * u * MXs * MQs * UU - 8 * s * u * MXs * MQs * UU * AXG -
                16 * s * u * MXs * MUs * SS * AXG - 8 * s * u * MXs * MUs * UU +
                24 * s * u * MXs * MUs * UU * AXG -
                12 * s * u * pow(MXs, 2) * SS +
                4 * s * u * pow(MXs, 2) * SS * AXG -
                12 * s * u * pow(MXs, 2) * UU +
                36 * s * u * pow(MXs, 2) * UU * AXG +
                24 * s * u * pow(MX, 4) * SS -
                16 * s * u * pow(MX, 4) * SS * AXG +
                16 * s * u * pow(MX, 4) * UU -
                24 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * u * pow(MU, 4) * SS - 4 * s * pow(u, 2) * MQs * UU +
                4 * s * pow(u, 2) * MQs * UU * AXG -
                12 * s * pow(u, 2) * MUs * SS +
                4 * s * pow(u, 2) * MUs * SS * AXG) +
           syFC4 *
               (8 * s * pow(u, 2) * MUs * UU -
                8 * s * pow(u, 2) * MUs * UU * AXG -
                12 * s * pow(u, 2) * MXs * SS +
                20 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU -
                20 * s * pow(u, 2) * MXs * UU * AXG + 8 * s * pow(u, 3) * SS -
                8 * s * pow(u, 3) * SS * AXG + 8 * s * pow(u, 3) * UU * AXG +
                2 * pow(s, 2) * MXs * MQs * UU -
                2 * pow(s, 2) * MXs * MQs * UU * AXG +
                12 * pow(s, 2) * MXs * MUs * SS -
                12 * pow(s, 2) * MXs * MUs * SS * AXG -
                14 * pow(s, 2) * MXs * MUs * UU +
                14 * pow(s, 2) * MXs * MUs * UU * AXG -
                16 * pow(s, 2) * pow(MXs, 2) * UU +
                16 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                16 * pow(s, 2) * pow(MX, 4) * UU -
                16 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                8 * pow(s, 2) * pow(MU, 4) * SS - 2 * pow(s, 2) * u * MQs * UU +
                2 * pow(s, 2) * u * MQs * UU * AXG -
                16 * pow(s, 2) * u * MUs * SS +
                8 * pow(s, 2) * u * MUs * SS * AXG +
                18 * pow(s, 2) * u * MUs * UU -
                10 * pow(s, 2) * u * MUs * UU * AXG -
                12 * pow(s, 2) * u * MXs * SS +
                12 * pow(s, 2) * u * MXs * SS * AXG -
                4 * pow(s, 2) * u * MXs * UU -
                12 * pow(s, 2) * u * MXs * UU * AXG +
                8 * pow(s, 2) * pow(u, 2) * SS -
                8 * pow(s, 2) * pow(u, 2) * SS * AXG +
                8 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC5 * (4 * s * pow(MUs, 2) * SS - 4 * s * MXs * MUs * UU -
                    4 * s * pow(MXs, 2) * UU + 4 * s * pow(MX, 4) * UU * AXG -
                    8 * s * u * MUs * SS + 4 * s * u * MUs * UU +
                    12 * s * u * MXs * UU - 8 * s * u * MXs * UU * AXG +
                    4 * s * pow(u, 2) * SS - 8 * s * pow(u, 2) * UU +
                    4 * s * pow(u, 2) * UU * AXG - 4 * pow(s, 2) * MUs * SS +
                    2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG +
                    4 * pow(s, 2) * u * SS - 6 * pow(s, 2) * u * UU +
                    2 * pow(s, 2) * u * UU * AXG) +
           syFC6 * (-4 * s * pow(MX, 4) * UU + 4 * s * pow(MX, 4) * UU * AXG +
                    8 * s * u * MXs * UU - 8 * s * u * MXs * UU * AXG -
                    4 * s * pow(u, 2) * UU + 4 * s * pow(u, 2) * UU * AXG +
                    2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG -
                    2 * pow(s, 2) * u * UU + 2 * pow(s, 2) * u * UU * AXG) +
           syFC7 * (4 * s * pow(MUs, 2) * SS + 4 * s * MXs * MUs * SS -
                    4 * s * MXs * MUs * UU - 4 * s * pow(MXs, 2) * UU +
                    4 * s * pow(MX, 4) * UU - 8 * s * u * MUs * SS +
                    4 * s * u * MUs * UU - 4 * s * u * MXs * SS +
                    4 * s * u * MXs * UU * AXG + 4 * s * pow(u, 2) * SS -
                    4 * s * pow(u, 2) * UU * AXG - 4 * pow(s, 2) * MUs * SS +
                    4 * pow(s, 2) * u * SS - 2 * pow(s, 2) * u * UU -
                    2 * pow(s, 2) * u * UU * AXG) +
           syFC8 *
               (2 * MXs * pow(MUs, 2) * SS + 2 * pow(MXs, 2) * MUs * SS -
                4 * pow(MX, 4) * MUs * SS + 2 * pow(MX, 4) * MUs * SS * AXG -
                4 * pow(MX, 4) * MUs * UU * AXG -
                4 * pow(MX, 4) * MXs * UU * AXG + 4 * pow(MX, 6) * UU * AXG -
                2 * u * pow(MUs, 2) * SS - 4 * u * MXs * MUs * SS * AXG +
                8 * u * MXs * MUs * UU * AXG - 2 * u * pow(MXs, 2) * SS +
                8 * u * pow(MXs, 2) * UU * AXG + 4 * u * pow(MX, 4) * SS -
                2 * u * pow(MX, 4) * SS * AXG - 4 * u * pow(MX, 4) * UU * AXG +
                2 * pow(u, 2) * MUs * SS + 2 * pow(u, 2) * MUs * SS * AXG -
                4 * pow(u, 2) * MUs * UU * AXG - 2 * pow(u, 2) * MXs * SS +
                4 * pow(u, 2) * MXs * SS * AXG -
                8 * pow(u, 2) * MXs * UU * AXG - 2 * pow(u, 3) * SS * AXG +
                4 * pow(u, 3) * UU * AXG + 2 * s * MXs * MUs * SS -
                2 * s * MXs * MUs * SS * AXG - 3 * s * MXs * MUs * UU +
                3 * s * MXs * MUs * UU * AXG - 3 * s * pow(MXs, 2) * UU +
                3 * s * pow(MXs, 2) * UU * AXG + 3 * s * pow(MX, 4) * UU -
                3 * s * pow(MX, 4) * UU * AXG + 2 * s * u * MUs * SS * AXG +
                3 * s * u * MUs * UU - 3 * s * u * MUs * UU * AXG -
                2 * s * u * MXs * SS + 2 * s * u * MXs * SS * AXG -
                s * u * MXs * UU - 3 * s * u * MXs * UU * AXG -
                2 * s * pow(u, 2) * SS * AXG + s * pow(u, 2) * UU +
                3 * s * pow(u, 2) * UU * AXG) +
           syFC9 *
               (4 * s * pow(MUs, 2) * MQs * SS - 4 * s * MXs * MUs * MQs * UU +
                4 * s * MXs * pow(MUs, 2) * UU -
                4 * s * pow(MXs, 2) * MQs * UU +
                4 * s * pow(MXs, 2) * MUs * UU + 4 * s * pow(MX, 4) * MQs * UU -
                4 * s * pow(MX, 4) * MUs * UU - 4 * s * pow(MU, 4) * MUs * SS -
                8 * s * u * MUs * MQs * SS + 4 * s * u * MUs * MQs * UU +
                4 * s * u * pow(MUs, 2) * SS - 4 * s * u * pow(MUs, 2) * UU +
                4 * s * u * MXs * MQs * UU - 4 * s * u * MXs * MUs * UU +
                4 * s * u * pow(MU, 4) * SS + 4 * s * pow(u, 2) * MQs * SS -
                4 * s * pow(u, 2) * MQs * UU - 4 * s * pow(u, 2) * MUs * SS +
                4 * s * pow(u, 2) * MUs * UU - 4 * pow(s, 2) * MUs * MQs * SS +
                2 * pow(s, 2) * MXs * MUs * SS -
                2 * pow(s, 2) * MXs * MUs * SS * AXG -
                2 * pow(s, 2) * MXs * MUs * UU +
                2 * pow(s, 2) * MXs * MUs * UU * AXG -
                2 * pow(s, 2) * pow(MXs, 2) * UU +
                2 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                2 * pow(s, 2) * pow(MX, 4) * UU -
                2 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(s, 2) * pow(MU, 4) * SS + 4 * pow(s, 2) * u * MQs * SS -
                4 * pow(s, 2) * u * MQs * UU - 6 * pow(s, 2) * u * MUs * SS +
                2 * pow(s, 2) * u * MUs * SS * AXG +
                6 * pow(s, 2) * u * MUs * UU -
                2 * pow(s, 2) * u * MUs * UU * AXG -
                2 * pow(s, 2) * u * MXs * SS +
                2 * pow(s, 2) * u * MXs * SS * AXG +
                2 * pow(s, 2) * u * MXs * UU -
                2 * pow(s, 2) * u * MXs * UU * AXG +
                2 * pow(s, 2) * pow(u, 2) * SS -
                2 * pow(s, 2) * pow(u, 2) * SS * AXG -
                2 * pow(s, 2) * pow(u, 2) * UU) +
           syFC9 * (2 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC10 *
               (8 * MXs * pow(MUs, 3) * SS +
                16 * pow(MXs, 2) * pow(MUs, 2) * SS +
                8 * pow(MXs, 3) * MUs * SS -
                32 * pow(MX, 4) * pow(MUs, 2) * SS +
                8 * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                16 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                32 * pow(MX, 4) * MXs * MUs * SS +
                8 * pow(MX, 4) * MXs * MUs * SS * AXG -
                32 * pow(MX, 4) * MXs * MUs * UU * AXG -
                16 * pow(MX, 4) * pow(MXs, 2) * UU * AXG +
                32 * pow(MX, 6) * MUs * SS - 16 * pow(MX, 6) * MUs * SS * AXG +
                48 * pow(MX, 6) * MUs * UU * AXG +
                48 * pow(MX, 6) * MXs * UU * AXG - 32 * pow(MX, 8) * UU * AXG -
                8 * u * pow(MUs, 3) * SS + 8 * u * MXs * pow(MUs, 2) * SS -
                16 * u * MXs * pow(MUs, 2) * SS * AXG +
                32 * u * MXs * pow(MUs, 2) * UU * AXG +
                8 * u * pow(MXs, 2) * MUs * SS -
                16 * u * pow(MXs, 2) * MUs * SS * AXG +
                64 * u * pow(MXs, 2) * MUs * UU * AXG -
                8 * u * pow(MXs, 3) * SS + 32 * u * pow(MXs, 3) * UU * AXG +
                24 * u * pow(MX, 4) * MUs * SS * AXG -
                80 * u * pow(MX, 4) * MUs * UU * AXG +
                32 * u * pow(MX, 4) * MXs * SS -
                8 * u * pow(MX, 4) * MXs * SS * AXG -
                80 * u * pow(MX, 4) * MXs * UU * AXG -
                32 * u * pow(MX, 6) * SS + 16 * u * pow(MX, 6) * SS * AXG +
                32 * u * pow(MX, 6) * UU * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS +
                8 * pow(u, 2) * pow(MUs, 2) * SS * AXG -
                16 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                16 * pow(u, 2) * MXs * MUs * SS +
                8 * pow(u, 2) * MXs * MUs * SS * AXG) +
           syFC10 *
               (-16 * pow(u, 2) * MXs * MUs * UU * AXG -
                24 * pow(u, 2) * pow(MXs, 2) * SS +
                16 * pow(u, 2) * pow(MXs, 2) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * SS -
                32 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * UU * AXG -
                8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG +
                8 * pow(u, 3) * MXs * SS * AXG -
                16 * pow(u, 3) * MXs * UU * AXG +
                18 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * SS * AXG -
                8 * s * MXs * pow(MUs, 2) * UU +
                16 * s * MXs * pow(MUs, 2) * UU * AXG +
                18 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * SS * AXG -
                20 * s * pow(MXs, 2) * MUs * UU +
                28 * s * pow(MXs, 2) * MUs * UU * AXG -
                12 * s * pow(MXs, 3) * UU + 12 * s * pow(MXs, 3) * UU * AXG -
                36 * s * pow(MX, 4) * MUs * SS +
                26 * s * pow(MX, 4) * MUs * SS * AXG +
                32 * s * pow(MX, 4) * MUs * UU -
                60 * s * pow(MX, 4) * MUs * UU * AXG +
                36 * s * pow(MX, 4) * MXs * UU -
                56 * s * pow(MX, 4) * MXs * UU * AXG -
                24 * s * pow(MX, 6) * UU + 44 * s * pow(MX, 6) * UU * AXG -
                8 * s * pow(MU, 4) * MUs * SS - 4 * s * pow(MU, 4) * MXs * SS -
                4 * s * pow(MU, 4) * MXs * SS * AXG +
                6 * s * u * pow(MUs, 2) * SS +
                8 * s * u * pow(MUs, 2) * SS * AXG +
                8 * s * u * pow(MUs, 2) * UU -
                16 * s * u * pow(MUs, 2) * UU * AXG +
                4 * s * u * MXs * MUs * SS - 16 * s * u * MXs * MUs * SS * AXG -
                20 * s * u * MXs * MUs * UU +
                28 * s * u * MXs * MUs * UU * AXG) +
           syFC10 *
               (-18 * s * u * pow(MXs, 2) * SS +
                8 * s * u * pow(MXs, 2) * SS * AXG -
                36 * s * u * pow(MXs, 2) * UU +
                52 * s * u * pow(MXs, 2) * UU * AXG +
                36 * s * u * pow(MX, 4) * SS -
                26 * s * u * pow(MX, 4) * SS * AXG +
                40 * s * u * pow(MX, 4) * UU -
                20 * s * u * pow(MX, 4) * UU * AXG +
                4 * s * u * pow(MU, 4) * SS +
                4 * s * u * pow(MU, 4) * SS * AXG -
                10 * s * pow(u, 2) * MUs * SS -
                2 * s * pow(u, 2) * MUs * SS * AXG +
                8 * s * pow(u, 2) * MUs * UU +
                4 * s * pow(u, 2) * MUs * UU * AXG -
                18 * s * pow(u, 2) * MXs * SS +
                28 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU -
                44 * s * pow(u, 2) * MXs * UU * AXG + 8 * s * pow(u, 3) * SS -
                10 * s * pow(u, 3) * SS * AXG + 12 * s * pow(u, 3) * UU * AXG +
                10 * pow(s, 2) * MXs * MUs * SS -
                10 * pow(s, 2) * MXs * MUs * SS * AXG -
                15 * pow(s, 2) * MXs * MUs * UU +
                15 * pow(s, 2) * MXs * MUs * UU * AXG -
                15 * pow(s, 2) * pow(MXs, 2) * UU +
                15 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                15 * pow(s, 2) * pow(MX, 4) * UU -
                15 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                8 * pow(s, 2) * pow(MU, 4) * SS -
                14 * pow(s, 2) * u * MUs * SS +
                8 * pow(s, 2) * u * MUs * SS * AXG +
                17 * pow(s, 2) * u * MUs * UU -
                9 * pow(s, 2) * u * MUs * UU * AXG -
                10 * pow(s, 2) * u * MXs * SS +
                10 * pow(s, 2) * u * MXs * SS * AXG -
                5 * pow(s, 2) * u * MXs * UU -
                15 * pow(s, 2) * u * MXs * UU * AXG +
                6 * pow(s, 2) * pow(u, 2) * SS -
                8 * pow(s, 2) * pow(u, 2) * SS * AXG) +
           syFC10 * (3 * pow(s, 2) * pow(u, 2) * UU +
                     9 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC11 *
               (4 * s * MXs * MUs * SS - 4 * s * pow(MX, 4) * UU +
                4 * s * pow(MX, 4) * UU * AXG - 4 * s * u * MXs * SS +
                4 * s * u * MXs * UU - 4 * s * u * MXs * UU * AXG +
                2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG) +
           syFC12 *
               (4 * s * MXs * MUs * SS - 4 * s * pow(MX, 4) * UU +
                4 * s * pow(MX, 4) * UU * AXG - 4 * s * u * MXs * SS +
                4 * s * u * MXs * UU - 4 * s * u * MXs * UU * AXG +
                2 * pow(s, 2) * MXs * UU - 2 * pow(s, 2) * MXs * UU * AXG) +
           syFC13 * (2 * s * MXs * MUs * SS - 2 * s * MXs * MUs * SS * AXG -
                     2 * s * MXs * MUs * UU + 2 * s * MXs * MUs * UU * AXG -
                     2 * s * pow(MXs, 2) * UU + 2 * s * pow(MXs, 2) * UU * AXG +
                     2 * s * pow(MX, 4) * UU - 2 * s * pow(MX, 4) * UU * AXG -
                     2 * s * u * MUs * SS + 2 * s * u * MUs * SS * AXG +
                     2 * s * u * MUs * UU - 2 * s * u * MUs * UU * AXG -
                     2 * s * u * MXs * SS + 2 * s * u * MXs * SS * AXG +
                     2 * s * u * MXs * UU - 2 * s * u * MXs * UU * AXG +
                     2 * s * pow(u, 2) * SS - 2 * s * pow(u, 2) * SS * AXG -
                     2 * s * pow(u, 2) * UU + 2 * s * pow(u, 2) * UU * AXG) +
           syFC14 *
               (-4 * s * pow(MX, 4) * MQs * UU +
                4 * s * pow(MX, 4) * MQs * UU * AXG +
                4 * s * pow(MX, 4) * MUs * UU -
                4 * s * pow(MX, 4) * MUs * UU * AXG +
                8 * s * u * MXs * MQs * UU - 8 * s * u * MXs * MQs * UU * AXG -
                8 * s * u * MXs * MUs * UU + 8 * s * u * MXs * MUs * UU * AXG -
                4 * s * pow(u, 2) * MQs * UU +
                4 * s * pow(u, 2) * MQs * UU * AXG +
                4 * s * pow(u, 2) * MUs * UU -
                4 * s * pow(u, 2) * MUs * UU * AXG +
                2 * pow(s, 2) * MXs * MQs * UU -
                2 * pow(s, 2) * MXs * MQs * UU * AXG -
                2 * pow(s, 2) * MXs * MUs * UU +
                2 * pow(s, 2) * MXs * MUs * UU * AXG -
                2 * pow(s, 2) * pow(MX, 4) * UU +
                2 * pow(s, 2) * pow(MX, 4) * UU * AXG -
                2 * pow(s, 2) * u * MQs * UU +
                2 * pow(s, 2) * u * MQs * UU * AXG +
                2 * pow(s, 2) * u * MUs * UU -
                2 * pow(s, 2) * u * MUs * UU * AXG +
                4 * pow(s, 2) * u * MXs * UU -
                4 * pow(s, 2) * u * MXs * UU * AXG -
                2 * pow(s, 2) * pow(u, 2) * UU +
                2 * pow(s, 2) * pow(u, 2) * UU * AXG + pow(s, 3) * MXs * UU -
                pow(s, 3) * MXs * UU * AXG - pow(s, 3) * u * UU +
                pow(s, 3) * u * UU * AXG) +
           syFC15 *
               (12 * s * MXs * pow(MUs, 2) * SS -
                12 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                12 * s * MXs * pow(MUs, 2) * UU * AXG +
                12 * s * pow(MXs, 2) * MUs * SS -
                12 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                24 * s * pow(MXs, 2) * MUs * UU * AXG -
                12 * s * pow(MXs, 3) * UU + 12 * s * pow(MXs, 3) * UU * AXG -
                20 * s * pow(MX, 4) * MUs * SS +
                20 * s * pow(MX, 4) * MUs * SS * AXG +
                24 * s * pow(MX, 4) * MUs * UU -
                32 * s * pow(MX, 4) * MUs * UU * AXG +
                32 * s * pow(MX, 4) * MXs * UU -
                32 * s * pow(MX, 4) * MXs * UU * AXG -
                20 * s * pow(MX, 6) * UU + 20 * s * pow(MX, 6) * UU * AXG -
                8 * s * pow(MU, 4) * MUs * SS + 4 * s * u * pow(MUs, 2) * SS +
                12 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                12 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * SS * AXG -
                8 * s * u * MXs * MUs * UU - 8 * s * u * MXs * MUs * UU * AXG -
                12 * s * u * pow(MXs, 2) * SS +
                12 * s * u * pow(MXs, 2) * SS * AXG -
                12 * s * u * pow(MXs, 2) * UU +
                4 * s * u * pow(MXs, 2) * UU * AXG +
                20 * s * u * pow(MX, 4) * SS -
                20 * s * u * pow(MX, 4) * SS * AXG +
                4 * s * u * pow(MX, 4) * UU +
                4 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * u * pow(MU, 4) * SS - 8 * s * pow(u, 2) * MUs * SS -
                16 * s * pow(u, 2) * MUs * SS * AXG +
                16 * s * pow(u, 2) * MUs * UU * AXG -
                4 * s * pow(u, 2) * MXs * SS +
                4 * s * pow(u, 2) * MXs * SS * AXG) +
           syFC15 *
               (12 * s * pow(u, 2) * MXs * UU -
                4 * s * pow(u, 2) * MXs * UU * AXG + 4 * s * pow(u, 3) * SS +
                4 * s * pow(u, 3) * SS * AXG - 4 * s * pow(u, 3) * UU -
                4 * s * pow(u, 3) * UU * AXG + 10 * pow(s, 2) * MXs * MUs * SS -
                10 * pow(s, 2) * MXs * MUs * SS * AXG -
                10 * pow(s, 2) * MXs * MUs * UU +
                10 * pow(s, 2) * MXs * MUs * UU * AXG -
                10 * pow(s, 2) * pow(MXs, 2) * UU +
                10 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                10 * pow(s, 2) * pow(MX, 4) * UU -
                10 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                8 * pow(s, 2) * pow(MU, 4) * SS -
                14 * pow(s, 2) * u * MUs * SS -
                2 * pow(s, 2) * u * MUs * SS * AXG +
                6 * pow(s, 2) * u * MUs * UU +
                2 * pow(s, 2) * u * MUs * UU * AXG -
                10 * pow(s, 2) * u * MXs * SS +
                10 * pow(s, 2) * u * MXs * SS * AXG +
                10 * pow(s, 2) * u * MXs * UU -
                10 * pow(s, 2) * u * MXs * UU * AXG +
                6 * pow(s, 2) * pow(u, 2) * SS +
                2 * pow(s, 2) * pow(u, 2) * SS * AXG -
                6 * pow(s, 2) * pow(u, 2) * UU -
                2 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC16 *
               (-4 * s * MXs * MUs * MQs * SS - 2 * s * MXs * pow(MUs, 2) * SS -
                2 * s * pow(MXs, 2) * MUs * SS + 4 * s * pow(MX, 4) * MUs * SS -
                2 * s * pow(MX, 4) * MUs * SS * AXG +
                4 * s * pow(MX, 4) * MUs * UU * AXG +
                4 * s * pow(MX, 4) * MXs * UU * AXG -
                4 * s * pow(MX, 6) * UU * AXG + 4 * s * pow(MU, 4) * MXs * SS +
                2 * s * u * pow(MUs, 2) * SS + 4 * s * u * MXs * MQs * SS +
                4 * s * u * MXs * MQs * UU - 4 * s * u * MXs * MQs * UU * AXG -
                4 * s * u * MXs * MUs * SS + 4 * s * u * MXs * MUs * SS * AXG -
                4 * s * u * MXs * MUs * UU - 4 * s * u * MXs * MUs * UU * AXG +
                2 * s * u * pow(MXs, 2) * SS -
                8 * s * u * pow(MXs, 2) * UU * AXG -
                4 * s * u * pow(MX, 4) * SS +
                2 * s * u * pow(MX, 4) * SS * AXG +
                4 * s * u * pow(MX, 4) * UU * AXG -
                4 * s * pow(u, 2) * MQs * UU +
                4 * s * pow(u, 2) * MQs * UU * AXG -
                2 * s * pow(u, 2) * MUs * SS -
                2 * s * pow(u, 2) * MUs * SS * AXG +
                4 * s * pow(u, 2) * MUs * UU + 2 * s * pow(u, 2) * MXs * SS -
                4 * s * pow(u, 2) * MXs * SS * AXG +
                8 * s * pow(u, 2) * MXs * UU * AXG +
                2 * s * pow(u, 3) * SS * AXG - 4 * s * pow(u, 3) * UU * AXG +
                pow(s, 2) * MXs * MUs * UU - pow(s, 2) * MXs * MUs * UU * AXG +
                pow(s, 2) * pow(MXs, 2) * UU -
                pow(s, 2) * pow(MXs, 2) * UU * AXG -
                pow(s, 2) * pow(MX, 4) * UU +
                pow(s, 2) * pow(MX, 4) * UU * AXG -
                2 * pow(s, 2) * u * MQs * UU +
                2 * pow(s, 2) * u * MQs * UU * AXG -
                2 * pow(s, 2) * u * MUs * SS + pow(s, 2) * u * MUs * UU -
                pow(s, 2) * u * MUs * UU * AXG) +
           syFC16 *
               (3 * pow(s, 2) * u * MXs * UU + pow(s, 2) * u * MXs * UU * AXG +
                2 * pow(s, 2) * pow(u, 2) * SS -
                3 * pow(s, 2) * pow(u, 2) * UU -
                pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC17 *
               (20 * s * MXs * pow(MUs, 2) * SS -
                12 * s * MXs * pow(MUs, 2) * SS * AXG -
                4 * s * MXs * pow(MUs, 2) * UU +
                12 * s * MXs * pow(MUs, 2) * UU * AXG +
                20 * s * pow(MXs, 2) * MUs * SS -
                12 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                24 * s * pow(MXs, 2) * MUs * UU * AXG -
                12 * s * pow(MXs, 3) * UU + 12 * s * pow(MXs, 3) * UU * AXG -
                40 * s * pow(MX, 4) * MUs * SS +
                32 * s * pow(MX, 4) * MUs * SS * AXG +
                28 * s * pow(MX, 4) * MUs * UU -
                52 * s * pow(MX, 4) * MUs * UU * AXG +
                44 * s * pow(MX, 4) * MXs * UU -
                60 * s * pow(MX, 4) * MXs * UU * AXG -
                40 * s * pow(MX, 6) * UU + 56 * s * pow(MX, 6) * UU * AXG -
                8 * s * pow(MU, 4) * MUs * SS - 4 * s * u * pow(MUs, 2) * SS +
                12 * s * u * pow(MUs, 2) * SS * AXG +
                4 * s * u * pow(MUs, 2) * UU -
                12 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS * AXG -
                16 * s * u * MXs * MUs * UU +
                32 * s * u * MXs * MUs * UU * AXG -
                20 * s * u * pow(MXs, 2) * SS +
                12 * s * u * pow(MXs, 2) * SS * AXG -
                36 * s * u * pow(MXs, 2) * UU +
                60 * s * u * pow(MXs, 2) * UU * AXG +
                40 * s * u * pow(MX, 4) * SS -
                32 * s * u * pow(MX, 4) * SS * AXG +
                48 * s * u * pow(MX, 4) * UU -
                56 * s * u * pow(MX, 4) * UU * AXG +
                8 * s * u * pow(MU, 4) * SS - 4 * s * pow(u, 2) * MUs * SS -
                4 * s * pow(u, 2) * MUs * SS * AXG +
                4 * s * pow(u, 2) * MUs * UU -
                4 * s * pow(u, 2) * MUs * UU * AXG -
                20 * s * pow(u, 2) * MXs * SS +
                28 * s * pow(u, 2) * MXs * SS * AXG) +
           syFC17 *
               (-4 * s * pow(u, 2) * MXs * UU -
                20 * s * pow(u, 2) * MXs * UU * AXG + 8 * s * pow(u, 3) * SS -
                8 * s * pow(u, 3) * SS * AXG + 8 * s * pow(u, 3) * UU * AXG +
                20 * pow(s, 2) * MXs * MUs * SS -
                20 * pow(s, 2) * MXs * MUs * SS * AXG -
                24 * pow(s, 2) * MXs * MUs * UU +
                24 * pow(s, 2) * MXs * MUs * UU * AXG -
                28 * pow(s, 2) * pow(MXs, 2) * UU +
                28 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                42 * pow(s, 2) * pow(MX, 4) * UU -
                42 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                8 * pow(s, 2) * pow(MU, 4) * SS -
                16 * pow(s, 2) * u * MUs * SS +
                8 * pow(s, 2) * u * MUs * SS * AXG +
                20 * pow(s, 2) * u * MUs * UU -
                12 * pow(s, 2) * u * MUs * UU * AXG -
                20 * pow(s, 2) * u * MXs * SS +
                20 * pow(s, 2) * u * MXs * SS * AXG -
                12 * pow(s, 2) * u * MXs * UU -
                4 * pow(s, 2) * u * MXs * UU * AXG +
                8 * pow(s, 2) * pow(u, 2) * SS -
                8 * pow(s, 2) * pow(u, 2) * SS * AXG +
                2 * pow(s, 2) * pow(u, 2) * UU +
                6 * pow(s, 2) * pow(u, 2) * UU * AXG -
                5 * pow(s, 3) * MXs * UU + 5 * pow(s, 3) * MXs * UU * AXG +
                pow(s, 3) * u * UU - pow(s, 3) * u * UU * AXG) +
           syFC18 *
               (16 * MXs * pow(MUs, 2) * SS + 16 * pow(MXs, 2) * MUs * SS -
                32 * pow(MX, 4) * MUs * SS + 16 * pow(MX, 4) * MUs * SS * AXG -
                32 * pow(MX, 4) * MUs * UU * AXG -
                32 * pow(MX, 4) * MXs * UU * AXG + 32 * pow(MX, 6) * UU * AXG -
                16 * u * pow(MUs, 2) * SS - 32 * u * MXs * MUs * SS * AXG +
                64 * u * MXs * MUs * UU * AXG - 16 * u * pow(MXs, 2) * SS +
                64 * u * pow(MXs, 2) * UU * AXG + 32 * u * pow(MX, 4) * SS -
                16 * u * pow(MX, 4) * SS * AXG -
                32 * u * pow(MX, 4) * UU * AXG + 16 * pow(u, 2) * MUs * SS +
                16 * pow(u, 2) * MUs * SS * AXG -
                32 * pow(u, 2) * MUs * UU * AXG - 16 * pow(u, 2) * MXs * SS +
                32 * pow(u, 2) * MXs * SS * AXG -
                64 * pow(u, 2) * MXs * UU * AXG - 16 * pow(u, 3) * SS * AXG +
                32 * pow(u, 3) * UU * AXG - 16 * s * pow(MUs, 2) * SS -
                16 * s * MXs * MUs * SS * AXG - 8 * s * MXs * MUs * UU +
                24 * s * MXs * MUs * UU * AXG - 8 * s * pow(MXs, 2) * UU +
                24 * s * pow(MXs, 2) * UU * AXG + 8 * s * pow(MX, 4) * UU -
                24 * s * pow(MX, 4) * UU * AXG + 32 * s * u * MUs * SS +
                16 * s * u * MUs * SS * AXG + 8 * s * u * MUs * UU -
                24 * s * u * MUs * UU * AXG + 16 * s * u * MXs * SS * AXG -
                8 * s * u * MXs * UU - 40 * s * u * MXs * UU * AXG -
                16 * s * pow(u, 2) * SS - 16 * s * pow(u, 2) * SS * AXG +
                8 * s * pow(u, 2) * UU + 40 * s * pow(u, 2) * UU * AXG +
                16 * pow(s, 2) * MUs * SS - 16 * pow(s, 2) * u * SS +
                8 * pow(s, 2) * u * UU + 8 * pow(s, 2) * u * UU * AXG) +
           syFC19 *
               (-16 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                16 * pow(MX, 4) * MXs * MUs * UU * AXG +
                16 * pow(MX, 6) * MUs * UU * AXG +
                8 * pow(MU, 4) * MXs * MUs * SS +
                8 * pow(MU, 4) * pow(MXs, 2) * SS -
                16 * pow(MU, 4) * pow(MX, 4) * SS +
                8 * pow(MU, 4) * pow(MX, 4) * SS * AXG -
                8 * u * MXs * pow(MUs, 2) * SS +
                32 * u * MXs * pow(MUs, 2) * UU * AXG -
                8 * u * pow(MXs, 2) * MUs * SS +
                32 * u * pow(MXs, 2) * MUs * UU * AXG +
                16 * u * pow(MX, 4) * MUs * SS -
                8 * u * pow(MX, 4) * MUs * SS * AXG -
                16 * u * pow(MX, 4) * MUs * UU * AXG -
                8 * u * pow(MU, 4) * MUs * SS + 8 * u * pow(MU, 4) * MXs * SS -
                16 * u * pow(MU, 4) * MXs * SS * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS -
                16 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 2) * MXs * MUs * SS +
                16 * pow(u, 2) * MXs * MUs * SS * AXG -
                32 * pow(u, 2) * MXs * MUs * UU * AXG +
                8 * pow(u, 2) * pow(MU, 4) * SS * AXG -
                8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG -
                12 * s * MXs * pow(MUs, 2) * UU +
                12 * s * MXs * pow(MUs, 2) * UU * AXG -
                12 * s * pow(MXs, 2) * MUs * UU +
                12 * s * pow(MXs, 2) * MUs * UU * AXG +
                12 * s * pow(MX, 4) * MUs * UU -
                12 * s * pow(MX, 4) * MUs * UU * AXG +
                8 * s * pow(MU, 4) * MXs * SS -
                8 * s * pow(MU, 4) * MXs * SS * AXG +
                12 * s * u * pow(MUs, 2) * UU -
                12 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS + 8 * s * u * MXs * MUs * SS * AXG) +
           syFC19 * (-4 * s * u * MXs * MUs * UU -
                     12 * s * u * MXs * MUs * UU * AXG +
                     8 * s * u * pow(MU, 4) * SS * AXG -
                     8 * s * pow(u, 2) * MUs * SS * AXG +
                     4 * s * pow(u, 2) * MUs * UU +
                     12 * s * pow(u, 2) * MUs * UU * AXG) +
           syFC20 *
               (8 * MXs * pow(MUs, 3) * SS +
                16 * pow(MXs, 2) * pow(MUs, 2) * SS +
                8 * pow(MXs, 3) * MUs * SS -
                24 * pow(MX, 4) * pow(MUs, 2) * SS +
                8 * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                16 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                24 * pow(MX, 4) * MXs * MUs * SS +
                8 * pow(MX, 4) * MXs * MUs * SS * AXG -
                32 * pow(MX, 4) * MXs * MUs * UU * AXG -
                16 * pow(MX, 4) * pow(MXs, 2) * UU * AXG +
                16 * pow(MX, 6) * MUs * SS - 8 * pow(MX, 6) * MUs * SS * AXG +
                32 * pow(MX, 6) * MUs * UU * AXG +
                32 * pow(MX, 6) * MXs * UU * AXG - 16 * pow(MX, 8) * UU * AXG -
                8 * u * pow(MUs, 3) * SS - 8 * u * MXs * pow(MUs, 2) * SS -
                16 * u * MXs * pow(MUs, 2) * SS * AXG +
                32 * u * MXs * pow(MUs, 2) * UU * AXG -
                8 * u * pow(MXs, 2) * MUs * SS -
                16 * u * pow(MXs, 2) * MUs * SS * AXG +
                64 * u * pow(MXs, 2) * MUs * UU * AXG -
                8 * u * pow(MXs, 3) * SS + 32 * u * pow(MXs, 3) * UU * AXG +
                24 * u * pow(MX, 4) * MUs * SS -
                32 * u * pow(MX, 4) * MUs * UU * AXG +
                24 * u * pow(MX, 4) * MXs * SS -
                8 * u * pow(MX, 4) * MXs * SS * AXG -
                32 * u * pow(MX, 4) * MXs * UU * AXG -
                16 * u * pow(MX, 6) * SS + 8 * u * pow(MX, 6) * SS * AXG +
                16 * pow(u, 2) * pow(MUs, 2) * SS +
                8 * pow(u, 2) * pow(MUs, 2) * SS * AXG -
                16 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                8 * pow(u, 2) * MXs * MUs * SS +
                32 * pow(u, 2) * MXs * MUs * SS * AXG -
                64 * pow(u, 2) * MXs * MUs * UU * AXG -
                8 * pow(u, 2) * pow(MXs, 2) * SS) +
           syFC20 *
               (16 * pow(u, 2) * pow(MXs, 2) * SS * AXG -
                48 * pow(u, 2) * pow(MXs, 2) * UU * AXG -
                8 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * UU * AXG -
                8 * pow(u, 3) * MUs * SS - 16 * pow(u, 3) * MUs * SS * AXG +
                32 * pow(u, 3) * MUs * UU * AXG + 8 * pow(u, 3) * MXs * SS -
                16 * pow(u, 3) * MXs * SS * AXG +
                32 * pow(u, 3) * MXs * UU * AXG + 8 * pow(u, 4) * SS * AXG -
                16 * pow(u, 4) * UU * AXG -
                8 * s * MXs * pow(MUs, 2) * SS * AXG -
                20 * s * MXs * pow(MUs, 2) * UU +
                20 * s * MXs * pow(MUs, 2) * UU * AXG -
                8 * s * pow(MXs, 2) * MUs * SS * AXG -
                32 * s * pow(MXs, 2) * MUs * UU +
                32 * s * pow(MXs, 2) * MUs * UU * AXG -
                12 * s * pow(MXs, 3) * UU + 12 * s * pow(MXs, 3) * UU * AXG +
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                32 * s * pow(MX, 4) * MUs * UU -
                32 * s * pow(MX, 4) * MUs * UU * AXG +
                24 * s * pow(MX, 4) * MXs * UU -
                24 * s * pow(MX, 4) * MXs * UU * AXG -
                12 * s * pow(MX, 6) * UU + 12 * s * pow(MX, 6) * UU * AXG +
                8 * s * pow(MU, 4) * MXs * SS -
                8 * s * pow(MU, 4) * MXs * SS * AXG +
                8 * s * u * pow(MUs, 2) * SS +
                8 * s * u * pow(MUs, 2) * SS * AXG +
                20 * s * u * pow(MUs, 2) * UU -
                20 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS +
                32 * s * u * MXs * MUs * SS * AXG +
                16 * s * u * MXs * MUs * UU -
                48 * s * u * MXs * MUs * UU * AXG +
                8 * s * u * pow(MXs, 2) * SS * AXG -
                4 * s * u * pow(MXs, 2) * UU -
                28 * s * u * pow(MXs, 2) * UU * AXG) +
           syFC20 * (-8 * s * u * pow(MX, 4) * SS * AXG +
                     4 * s * u * pow(MX, 4) * UU +
                     28 * s * u * pow(MX, 4) * UU * AXG -
                     8 * s * u * pow(MU, 4) * SS +
                     8 * s * u * pow(MU, 4) * SS * AXG -
                     32 * s * pow(u, 2) * MUs * SS * AXG -
                     16 * s * pow(u, 2) * MUs * UU +
                     48 * s * pow(u, 2) * MUs * UU * AXG +
                     8 * s * pow(u, 2) * MXs * SS -
                     16 * s * pow(u, 2) * MXs * SS * AXG +
                     4 * s * pow(u, 2) * MXs * UU +
                     28 * s * pow(u, 2) * MXs * UU * AXG +
                     16 * s * pow(u, 3) * SS * AXG - 4 * s * pow(u, 3) * UU -
                     28 * s * pow(u, 3) * UU * AXG -
                     8 * pow(s, 2) * u * MUs * SS * AXG -
                     12 * pow(s, 2) * u * MUs * UU +
                     12 * pow(s, 2) * u * MUs * UU * AXG +
                     8 * pow(s, 2) * pow(u, 2) * SS * AXG -
                     4 * pow(s, 2) * pow(u, 2) * UU -
                     12 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC21 *
               (8 * MXs * pow(MUs, 3) * SS +
                16 * pow(MXs, 2) * pow(MUs, 2) * SS +
                8 * pow(MXs, 3) * MUs * SS -
                32 * pow(MX, 4) * pow(MUs, 2) * SS +
                8 * pow(MX, 4) * pow(MUs, 2) * SS * AXG -
                16 * pow(MX, 4) * pow(MUs, 2) * UU * AXG -
                32 * pow(MX, 4) * MXs * MUs * SS +
                8 * pow(MX, 4) * MXs * MUs * SS * AXG -
                32 * pow(MX, 4) * MXs * MUs * UU * AXG -
                16 * pow(MX, 4) * pow(MXs, 2) * UU * AXG +
                32 * pow(MX, 6) * MUs * SS - 16 * pow(MX, 6) * MUs * SS * AXG +
                48 * pow(MX, 6) * MUs * UU * AXG +
                48 * pow(MX, 6) * MXs * UU * AXG - 32 * pow(MX, 8) * UU * AXG -
                8 * u * pow(MUs, 3) * SS + 8 * u * MXs * pow(MUs, 2) * SS -
                16 * u * MXs * pow(MUs, 2) * SS * AXG +
                32 * u * MXs * pow(MUs, 2) * UU * AXG +
                8 * u * pow(MXs, 2) * MUs * SS -
                16 * u * pow(MXs, 2) * MUs * SS * AXG +
                64 * u * pow(MXs, 2) * MUs * UU * AXG -
                8 * u * pow(MXs, 3) * SS + 32 * u * pow(MXs, 3) * UU * AXG +
                24 * u * pow(MX, 4) * MUs * SS * AXG -
                80 * u * pow(MX, 4) * MUs * UU * AXG +
                32 * u * pow(MX, 4) * MXs * SS -
                8 * u * pow(MX, 4) * MXs * SS * AXG -
                80 * u * pow(MX, 4) * MXs * UU * AXG -
                32 * u * pow(MX, 6) * SS + 16 * u * pow(MX, 6) * SS * AXG +
                32 * u * pow(MX, 6) * UU * AXG +
                8 * pow(u, 2) * pow(MUs, 2) * SS +
                8 * pow(u, 2) * pow(MUs, 2) * SS * AXG -
                16 * pow(u, 2) * pow(MUs, 2) * UU * AXG -
                16 * pow(u, 2) * MXs * MUs * SS +
                8 * pow(u, 2) * MXs * MUs * SS * AXG) +
           syFC21 *
               (-16 * pow(u, 2) * MXs * MUs * UU * AXG -
                24 * pow(u, 2) * pow(MXs, 2) * SS +
                16 * pow(u, 2) * pow(MXs, 2) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * SS -
                32 * pow(u, 2) * pow(MX, 4) * SS * AXG +
                32 * pow(u, 2) * pow(MX, 4) * UU * AXG -
                8 * pow(u, 3) * MUs * SS * AXG +
                16 * pow(u, 3) * MUs * UU * AXG +
                8 * pow(u, 3) * MXs * SS * AXG -
                16 * pow(u, 3) * MXs * UU * AXG +
                8 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * SS * AXG -
                20 * s * MXs * pow(MUs, 2) * UU +
                20 * s * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * SS * AXG -
                32 * s * pow(MXs, 2) * MUs * UU +
                32 * s * pow(MXs, 2) * MUs * UU * AXG -
                12 * s * pow(MXs, 3) * UU + 12 * s * pow(MXs, 3) * UU * AXG -
                24 * s * pow(MX, 4) * MUs * SS +
                24 * s * pow(MX, 4) * MUs * SS * AXG +
                52 * s * pow(MX, 4) * MUs * UU -
                68 * s * pow(MX, 4) * MUs * UU * AXG +
                36 * s * pow(MX, 4) * MXs * UU -
                52 * s * pow(MX, 4) * MXs * UU * AXG -
                24 * s * pow(MX, 6) * UU + 40 * s * pow(MX, 6) * UU * AXG +
                8 * s * pow(MU, 4) * MXs * SS -
                8 * s * pow(MU, 4) * MXs * SS * AXG +
                8 * s * u * pow(MUs, 2) * SS * AXG +
                20 * s * u * pow(MUs, 2) * UU -
                20 * s * u * pow(MUs, 2) * UU * AXG -
                8 * s * u * MXs * MUs * SS - 24 * s * u * MXs * MUs * UU +
                24 * s * u * MXs * MUs * UU * AXG -
                8 * s * u * pow(MXs, 2) * SS +
                8 * s * u * pow(MXs, 2) * SS * AXG -
                28 * s * u * pow(MXs, 2) * UU) +
           syFC21 *
               (28 * s * u * pow(MXs, 2) * UU * AXG +
                24 * s * u * pow(MX, 4) * SS -
                24 * s * u * pow(MX, 4) * SS * AXG +
                32 * s * u * pow(MX, 4) * UU - 8 * s * u * pow(MU, 4) * SS +
                8 * s * u * pow(MU, 4) * SS * AXG +
                8 * s * pow(u, 2) * MUs * SS -
                16 * s * pow(u, 2) * MUs * SS * AXG +
                4 * s * pow(u, 2) * MUs * UU +
                12 * s * pow(u, 2) * MUs * UU * AXG -
                8 * s * pow(u, 2) * MXs * SS +
                16 * s * pow(u, 2) * MXs * SS * AXG -
                4 * s * pow(u, 2) * MXs * UU -
                28 * s * pow(u, 2) * MXs * UU * AXG +
                8 * pow(s, 2) * MXs * MUs * SS -
                8 * pow(s, 2) * MXs * MUs * SS * AXG -
                16 * pow(s, 2) * MXs * MUs * UU +
                16 * pow(s, 2) * MXs * MUs * UU * AXG -
                12 * pow(s, 2) * pow(MXs, 2) * UU +
                12 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                12 * pow(s, 2) * pow(MX, 4) * UU -
                12 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(s, 2) * u * MUs * UU -
                4 * pow(s, 2) * u * MUs * UU * AXG -
                8 * pow(s, 2) * u * MXs * SS +
                8 * pow(s, 2) * u * MXs * SS * AXG -
                4 * pow(s, 2) * u * MXs * UU -
                12 * pow(s, 2) * u * MXs * UU * AXG) +
           syFC22 *
               (8 * s * MXs * pow(MUs, 2) * SS -
                8 * s * MXs * pow(MUs, 2) * SS * AXG -
                8 * s * MXs * pow(MUs, 2) * UU +
                8 * s * MXs * pow(MUs, 2) * UU * AXG +
                8 * s * pow(MXs, 2) * MUs * SS -
                8 * s * pow(MXs, 2) * MUs * SS * AXG -
                16 * s * pow(MXs, 2) * MUs * UU +
                16 * s * pow(MXs, 2) * MUs * UU * AXG -
                8 * s * pow(MXs, 3) * UU + 8 * s * pow(MXs, 3) * UU * AXG -
                8 * s * pow(MX, 4) * MUs * SS +
                8 * s * pow(MX, 4) * MUs * SS * AXG +
                16 * s * pow(MX, 4) * MUs * UU -
                16 * s * pow(MX, 4) * MUs * UU * AXG +
                16 * s * pow(MX, 4) * MXs * UU -
                16 * s * pow(MX, 4) * MXs * UU * AXG - 8 * s * pow(MX, 6) * UU +
                8 * s * pow(MX, 6) * UU * AXG - 8 * s * u * pow(MUs, 2) * SS +
                8 * s * u * pow(MUs, 2) * SS * AXG +
                8 * s * u * pow(MUs, 2) * UU -
                8 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS +
                16 * s * u * MXs * MUs * SS * AXG +
                16 * s * u * MXs * MUs * UU -
                16 * s * u * MXs * MUs * UU * AXG -
                8 * s * u * pow(MXs, 2) * SS +
                8 * s * u * pow(MXs, 2) * SS * AXG +
                8 * s * u * pow(MXs, 2) * UU -
                8 * s * u * pow(MXs, 2) * UU * AXG +
                8 * s * u * pow(MX, 4) * SS -
                8 * s * u * pow(MX, 4) * SS * AXG -
                8 * s * u * pow(MX, 4) * UU +
                8 * s * u * pow(MX, 4) * UU * AXG +
                16 * s * pow(u, 2) * MUs * SS -
                16 * s * pow(u, 2) * MUs * SS * AXG -
                16 * s * pow(u, 2) * MUs * UU +
                16 * s * pow(u, 2) * MUs * UU * AXG +
                8 * s * pow(u, 2) * MXs * SS -
                8 * s * pow(u, 2) * MXs * SS * AXG -
                8 * s * pow(u, 2) * MXs * UU +
                8 * s * pow(u, 2) * MXs * UU * AXG) +
           syFC22 * (-8 * s * pow(u, 3) * SS + 8 * s * pow(u, 3) * SS * AXG +
                     8 * s * pow(u, 3) * UU - 8 * s * pow(u, 3) * UU * AXG +
                     8 * pow(s, 2) * u * MUs * SS -
                     8 * pow(s, 2) * u * MUs * SS * AXG -
                     8 * pow(s, 2) * u * MUs * UU +
                     8 * pow(s, 2) * u * MUs * UU * AXG -
                     8 * pow(s, 2) * pow(u, 2) * SS +
                     8 * pow(s, 2) * pow(u, 2) * SS * AXG +
                     8 * pow(s, 2) * pow(u, 2) * UU -
                     8 * pow(s, 2) * pow(u, 2) * UU * AXG) +
           syFC23 *
               (16 * s * MXs * pow(MUs, 2) * SS -
                16 * s * MXs * pow(MUs, 2) * SS * AXG -
                16 * s * MXs * pow(MUs, 2) * UU +
                16 * s * MXs * pow(MUs, 2) * UU * AXG +
                16 * s * pow(MXs, 2) * MUs * SS -
                16 * s * pow(MXs, 2) * MUs * SS * AXG -
                32 * s * pow(MXs, 2) * MUs * UU +
                32 * s * pow(MXs, 2) * MUs * UU * AXG -
                16 * s * pow(MXs, 3) * UU + 16 * s * pow(MXs, 3) * UU * AXG -
                24 * s * pow(MX, 4) * MUs * SS +
                24 * s * pow(MX, 4) * MUs * SS * AXG +
                48 * s * pow(MX, 4) * MUs * UU -
                48 * s * pow(MX, 4) * MUs * UU * AXG +
                48 * s * pow(MX, 4) * MXs * UU -
                48 * s * pow(MX, 4) * MXs * UU * AXG -
                32 * s * pow(MX, 6) * UU + 32 * s * pow(MX, 6) * UU * AXG -
                16 * s * u * pow(MUs, 2) * SS +
                16 * s * u * pow(MUs, 2) * SS * AXG +
                16 * s * u * pow(MUs, 2) * UU -
                16 * s * u * pow(MUs, 2) * UU * AXG -
                16 * s * u * MXs * MUs * SS +
                16 * s * u * MXs * MUs * SS * AXG -
                16 * s * u * pow(MXs, 2) * SS +
                16 * s * u * pow(MXs, 2) * SS * AXG -
                16 * s * u * pow(MXs, 2) * UU +
                16 * s * u * pow(MXs, 2) * UU * AXG +
                24 * s * u * pow(MX, 4) * SS -
                24 * s * u * pow(MX, 4) * SS * AXG +
                24 * s * pow(u, 2) * MUs * SS -
                24 * s * pow(u, 2) * MUs * SS * AXG -
                16 * s * pow(u, 2) * MUs * UU +
                16 * s * pow(u, 2) * MUs * UU * AXG +
                16 * s * pow(u, 2) * MXs * UU -
                16 * s * pow(u, 2) * MXs * UU * AXG - 8 * s * pow(u, 3) * SS +
                8 * s * pow(u, 3) * SS * AXG + 8 * pow(s, 2) * MXs * MUs * SS -
                8 * pow(s, 2) * MXs * MUs * SS * AXG) +
           syFC23 * (-12 * pow(s, 2) * MXs * MUs * UU +
                     12 * pow(s, 2) * MXs * MUs * UU * AXG -
                     12 * pow(s, 2) * pow(MXs, 2) * UU +
                     12 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                     12 * pow(s, 2) * pow(MX, 4) * UU -
                     12 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                     8 * pow(s, 2) * u * MUs * SS -
                     8 * pow(s, 2) * u * MUs * SS * AXG -
                     4 * pow(s, 2) * u * MUs * UU +
                     4 * pow(s, 2) * u * MUs * UU * AXG -
                     8 * pow(s, 2) * u * MXs * SS +
                     8 * pow(s, 2) * u * MXs * SS * AXG +
                     20 * pow(s, 2) * u * MXs * UU -
                     20 * pow(s, 2) * u * MXs * UU * AXG -
                     8 * pow(s, 2) * pow(u, 2) * SS +
                     8 * pow(s, 2) * pow(u, 2) * SS * AXG -
                     4 * pow(s, 2) * pow(u, 2) * UU +
                     4 * pow(s, 2) * pow(u, 2) * UU * AXG -
                     4 * pow(s, 3) * u * UU + 4 * pow(s, 3) * u * UU * AXG) +
           syFC24 * (8 * s * MXs * pow(MUs, 2) * SS -
                     8 * s * MXs * pow(MUs, 2) * SS * AXG -
                     8 * s * MXs * pow(MUs, 2) * UU +
                     8 * s * MXs * pow(MUs, 2) * UU * AXG +
                     8 * s * pow(MXs, 2) * MUs * SS -
                     8 * s * pow(MXs, 2) * MUs * SS * AXG -
                     16 * s * pow(MXs, 2) * MUs * UU +
                     16 * s * pow(MXs, 2) * MUs * UU * AXG -
                     8 * s * pow(MXs, 3) * UU + 8 * s * pow(MXs, 3) * UU * AXG -
                     16 * s * pow(MX, 4) * MUs * SS +
                     16 * s * pow(MX, 4) * MUs * SS * AXG +
                     32 * s * pow(MX, 4) * MUs * UU -
                     32 * s * pow(MX, 4) * MUs * UU * AXG +
                     32 * s * pow(MX, 4) * MXs * UU -
                     32 * s * pow(MX, 4) * MXs * UU * AXG -
                     32 * s * pow(MX, 6) * UU + 32 * s * pow(MX, 6) * UU * AXG -
                     8 * s * u * pow(MUs, 2) * SS +
                     8 * s * u * pow(MUs, 2) * SS * AXG +
                     8 * s * u * pow(MUs, 2) * UU -
                     8 * s * u * pow(MUs, 2) * UU * AXG -
                     16 * s * u * MXs * MUs * UU +
                     16 * s * u * MXs * MUs * UU * AXG -
                     8 * s * u * pow(MXs, 2) * SS +
                     8 * s * u * pow(MXs, 2) * SS * AXG -
                     24 * s * u * pow(MXs, 2) * UU +
                     24 * s * u * pow(MXs, 2) * UU * AXG +
                     16 * s * u * pow(MX, 4) * SS -
                     16 * s * u * pow(MX, 4) * SS * AXG +
                     32 * s * u * pow(MX, 4) * UU -
                     32 * s * u * pow(MX, 4) * UU * AXG +
                     8 * s * pow(u, 2) * MUs * SS -
                     8 * s * pow(u, 2) * MUs * SS * AXG -
                     8 * s * pow(u, 2) * MXs * SS +
                     8 * s * pow(u, 2) * MXs * SS * AXG +
                     8 * pow(s, 2) * MXs * MUs * SS -
                     8 * pow(s, 2) * MXs * MUs * SS * AXG -
                     12 * pow(s, 2) * MXs * MUs * UU +
                     12 * pow(s, 2) * MXs * MUs * UU * AXG -
                     12 * pow(s, 2) * pow(MXs, 2) * UU) +
           syFC24 *
               (12 * pow(s, 2) * pow(MXs, 2) * UU * AXG +
                24 * pow(s, 2) * pow(MX, 4) * UU -
                24 * pow(s, 2) * pow(MX, 4) * UU * AXG +
                4 * pow(s, 2) * u * MUs * UU -
                4 * pow(s, 2) * u * MUs * UU * AXG -
                8 * pow(s, 2) * u * MXs * SS +
                8 * pow(s, 2) * u * MXs * SS * AXG -
                4 * pow(s, 2) * u * MXs * UU +
                4 * pow(s, 2) * u * MXs * UU * AXG - 4 * pow(s, 3) * MXs * UU +
                4 * pow(s, 3) * MXs * UU * AXG)));

  return ret.real();
}

ComplexType ME_us_box_qgQ(POLE pIEPS, bool sc, bool uc, bool axial, double Q2,
                          double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc, UU = uc, AXG = axial;
  ComplexType ret = 0;
  // int itsq = 0;
  // int itq  = 0;
  // for(int invert_color_flow= 0; invert_color_flow< 2;invert_color_flow++){
  // for(int itsq = 0; itsq < 6;itsq++){
  // for(int itq  = 0; itq  < 3;itq++){

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto Denom = [](auto a) { return 1. / a; };

#define syFC1 C1(0, MUs + MXs - s - t, MXs, 0, 0, MQs)
#define syFC2 C2(0, MUs + MXs - s - t, MXs, 0, 0, MQs)

  _EPS0(ret,
        -(gs * gs * gs * gs * (2 - 3 * Nc * Nc + Nc * Nc * Nc * Nc) *
          (Lp * R + L * Rp) * TR * TR *
          ((-1 + AXG) * syFC1 * (-MUs + s + t) * (-2 * MUs + s + 2 * t) * UU -
           MXs * syFC2 *
               (2 * MXs * SS - 2 * SS * t + 2 * (-1 + AXG) * (MUs - t) * UU +
                s * (-2 * SS + UU - AXG * UU)))) /
            (768 * Nc * Pi * Pi * (-MXs + s + t)));

  return ret.real();
}
