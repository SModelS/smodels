#include "pxs_gausq.h"

#define Power std::pow

ComplexType ME_us_born(bool sc1, bool uc1, bool sc2, bool uc2, bool axial, double Q2, double P1K1,
                       Parameters *params, bool qqleft, bool qqright, bool qqgleft, bool qqgright,
                       bool qQXleft, bool qQXright) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  double SS = sc2, UU = uc2, AXG = axial;
  double Sp = sc1;
  double Up = uc1;
  ComplexType ret = 0;

  double dqqright = qqright;
  double dqqleft = qqleft;
  double dqqgright = qqgright;
  double dqqgleft = qqgleft;
  double dqQXright = qQXright;
  double dqQXleft = qQXleft;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R) * dqQXleft;
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L) * dqQXright;

  auto Denom = [](auto a) { return 1. / a; };
  //*

  ret +=
      Up * ((Power(gs, 2) * (-1 + Power(Nc, 2)) * (Lp * R + L * Rp) * TR *
             (MUs * SS *
                  (2 * (-1 + AXG) * Power(MUs, 2) + 2 * (-2 + AXG) * Power(MX, 4) +
                   4 * Power(MXs, 2) - 2 * AXG * Power(MXs, 2) - Power(s, 2) + AXG * Power(s, 2) -
                   2 * MXs * t - s * t + 3 * AXG * s * t + 2 * AXG * Power(t, 2) +
                   MUs * (2 * MXs + 3 * s - 3 * AXG * s + 2 * t - 4 * AXG * t)) +
              4 * AXG * Power(MX, 6) * UU -
              2 * Power(MX, 4) *
                  ((-2 + AXG) * SS * (MUs + MXs - s - t) + (-1 + 3 * AXG) * s * UU +
                   2 * AXG * (2 * MUs + 2 * MXs - 2 * s - t) * UU) +
              MXs * ((-3 + AXG) * s * SS * (MUs + MXs - s - t) -
                     2 * SS * (MUs + MXs - s - t) * ((3 - 2 * AXG) * (MUs + MXs - s - t) + t) +
                     2 * (-1 + AXG) * Power(s, 2) * UU +
                     4 * AXG * (MUs + MXs - s - t) * (MUs + MXs - s + t) * UU +
                     2 * s * ((-3 + 5 * AXG) * (MUs + MXs - s - t) + (-1 + AXG) * t) * UU) +
              (MUs + MXs - s - t) *
                  (-2 * (-1 + AXG) * Power(s, 2) * UU -
                   (-1 + AXG) * s *
                       (SS * (MUs + MXs - s - t) + 2 * (2 * MUs + 2 * MXs - 2 * s - t) * UU) +
                   2 * (MUs + MXs - s - t) *
                       (SS * (MUs + MXs - s - AXG * (MUs + MXs - s - t)) - 2 * AXG * t * UU)))) /
            (96. * s * (-MXs + s + t) * (MQs - MUs - MXs + s + t))) +
      Sp * (-(Power(gs, 2) * (-1 + Power(Nc, 2)) *
              (Lp * R * dqqgleft * dqqright + L * Rp * dqqgright * dqqleft) * TR *
              (2 * SS * Power(t, 2) -
               2 *
                   ((-1 + AXG) * Power(MUs, 2) + (-2 + AXG) * Power(MX, 4) + AXG * Power(t, 2) +
                    MUs * (t - 2 * AXG * t)) *
                   UU +
               2 * Power(MXs, 2) * (SS + (-2 + AXG) * UU) + Power(s, 2) * (UU - AXG * UU) -
               2 * MXs * (s * SS + 2 * SS * t + (MUs - t) * UU) +
               s * (2 * SS * t + (3 * (-1 + AXG) * MUs + t - 3 * AXG * t) * UU))) /
            (96. * s * (-MXs + s + t)));

  // std::cout << (Lp * R + L * Rp) << " " << MUs << " " << MXs << " " << s << " " << t <<
  // std::endl;
  // std::cout << MUs << " " << MXs << " " << s << " " << t << std::endl;
  // std::cout << ret << std::endl;
  // exit(1);
  return ret;
}

// unused
ComplexType ME_us_born_eps1(double Q2, double P1K1, Parameters *params) {
  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  // auto MQ = MU;
  // auto MQs = MUs;
  ComplexType ret = 0;

  ComplexType L = (params->CHSQq[ch][sq][q].L);
  ComplexType R = (params->CHSQq[ch][sq][q].R);
  ComplexType Lp = conj(params->CHSQq[ch][sq][q].R);
  ComplexType Rp = conj(params->CHSQq[ch][sq][q].L);

  auto Denom = [](auto a) { return 1. / a; };
  //*

  ret += -Power(gs, 2) / 4. * 8. * (Lp * R + L * Rp) * (1 + (MXs - MUs - 2 * P1K2) / s);

  return ret;
}
