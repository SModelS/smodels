// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Computes the virtual part of the partonic cross section.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "clooptools.h"
#include "dipoles.h"
#include "kinematics.h"
#include "npf.h"
#include "options.h"
#include "params.h"
#include "pdf.h"
#include "pxs.h"
#include "pxs_gausq.h"
#include "pxs_gausq_2.h"
#include "tensors.h"
#include "utils.h"

using namespace Fastor;

#include "tensors.h"

#include "constants.h"

ComplexType ME_virt_tot(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                        Parameters *params) {
  ComplexType ret =

///*
#ifdef SQGA_VIRT_BOX
      ME_us_box_gqqQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_QGGq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_qgQQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_GQQq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_qggQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_QGqq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_box_qgQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_uug
      ME_us_qqg_qgg(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_qqg_qqg(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_qqg_QGG(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_qqg_QQG(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_qqg_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_UUg
      ME_us_QQg_Qgg(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_QQg_gQQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_QQg_Qg(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_QQg_qGG(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_QQg_Gqq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_QQg_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_UXu
      ME_uu_qQX_QGq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_ss_qQX_QGq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_uu_qQX_qgQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_ss_qQX_qgQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_ss_qQX_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_uu_qQX_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_gGX
      ME_us_gGX_Qqq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_gGX_qQQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_uu
      ME_us_se_qq_gq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_qq_GQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_qq_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
#ifdef SQGA_VIRT_UU
      ME_us_se_QQ_Gq(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_QQ_gQ(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_QQ_Q(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_QQ_ct(pIEPS, sc, uc, axial, Q2, P1K1, params) +
      ME_us_se_QQ_g(pIEPS, sc, uc, axial, Q2, P1K1, params) +
#endif
      0.;
  return ret;
}

double Virt_gaugino_squark(double S, double T, Parameters *params, POLE pIEPS) {
  int aa = params->in1;  // quark
  int bb = params->in2;  // antiquark
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino
  FI *ff = new FI();
  BORN_KINEMATICS;

  double Q2 = 2 * ff->papb;
  double P1K1 = ff->pap1;

  auto ret = ME_virt_tot(pIEPS, 1, 1, 1, Q2, P1K1, params).real();

  delete ff;
  return ret;
}
double Virt_gaugino_squark(double S, double T, Parameters *params) {
  return Virt_gaugino_squark(S, T, params, static_cast<POLE>(IEPS));
}

double DipI_gasq(const double S, const double T, Parameters *params, POLE pIEPS) {
  FI *ff = new FI();
  {
    int aa = params->in1;
    int bb = params->in2;
    int ii = params->out1;
    int jj = params->out2;
    BORN_KINEMATICS;
  }

  double Q2 = 2 * ff->papb;
  double P1K1 = ff->pap1;

  BOX_KINEMATIC;
  BOX_INDEX;
  BOX_BASE;
  double P2K1 = (-2. * P1K1 + MUs - MXs + Q2) / 2.;
  Tensor<ComplexType, 4> dip = {0., 0., 0., 0.};
  dip += I_jk(PARTON_SQUARK, PARTON_QUARK, 2. * P1K1, MU, 0., TR / NC, params);
  dip += I_jk(PARTON_QUARK, PARTON_SQUARK, 2. * P1K1, 0., MU, TR / NC, params);
  dip += I_jk(PARTON_QUARK, PARTON_GLUON, Q2, 0., 0., -TR * NC, params);
  dip += I_jk(PARTON_GLUON, PARTON_QUARK, Q2, 0., 0., -TR * NC, params);
  dip += I_jk(PARTON_GLUON, PARTON_SQUARK, 2. * P2K1, 0., MU, -TR * NC, params);
  dip += I_jk(PARTON_SQUARK, PARTON_GLUON, 2. * P2K1, MU, 0., -TR * NC, params);
  auto ME_B = (ME_us_born(1, 1, 1, 1, 1, Q2, P1K1, params, 1, 1, 1, 1, 1, 1));
  double ret = dip(2 - pIEPS).real() * (ME_B.real());
  delete ff;
#ifdef SQGA_VINTDIP
  return ret;
#else
  return 0.0;
#endif
}
double DipI_gasq(const double S, const double T, Parameters *params) {
  return DipI_gasq(S, T, params, FIN);
}

using BoxComp = Tensor<ComplexType, 11, 2, 4>;

Tensor<ComplexType, 3> born_ss_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1) {
  double MQ = Mass_o0;
  double MX = Mass_o1;
  double dM2o = MX * MX - MQ * MQ;
  double P1K2 = Q2 / 2 - P1K1;
  double Dn_p1mk2_MQ = dM2o - 2 * P1K2;
  Tensor<ComplexType, 3> dmin2_2 = {1, -DRBAR, 0};
  Tensor<ComplexType, 3> unit = {1, 0, 0};
  return +dmin2_2 * (1 + (dM2o - 2 * P1K2) / Q2);
}

Tensor<ComplexType, 3> born_uu_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1) {
  double MQ = Mass_o0;
  double MX = Mass_o1;
  double dM2o = MX * MX - MQ * MQ;
  double P1K2 = Q2 / 2 - P1K1;
  double Dn_p1mk2_MQ = dM2o - 2 * P1K2;
  Tensor<ComplexType, 3> dmin2_2 = {1, -DRBAR, 0};
  Tensor<ComplexType, 3> unit = {1, 0, 0};
  return +unit *
         (1 - 2 * MQ * MQ * dM2o / Dn_p1mk2_MQ / Dn_p1mk2_MQ - (MX * MX - 3 * MQ * MQ) / Dn_p1mk2_MQ
#ifdef PHYS_GAUGE
          + 1 - 2 * (dM2o + 2 * P1K2) / Q2 - dM2o / Dn_p1mk2_MQ +
          2 * dM2o * dM2o / (Q2 * Dn_p1mk2_MQ)
#endif
         );
}

Tensor<ComplexType, 3> born_su_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1) {
  double MQ = Mass_o0;
  double MX = Mass_o1;
  double dM2o = MX * MX - MQ * MQ;
  double P1K2 = Q2 / 2 - P1K1;
  double Dn_p1mk2_MQ = dM2o - 2 * P1K2;
  Tensor<ComplexType, 3> dmin2_2 = {1, -DRBAR, 0};
  Tensor<ComplexType, 3> unit = {1, 0, 0};
  return
      // not PHYS_GAUGE
      -unit * 1. / 2. *
      ((MX * MX + MQ * MQ) / Dn_p1mk2_MQ - 2 * dM2o * dM2o / (Dn_p1mk2_MQ * Q2) + 2 * dM2o / Q2 +
       1) *
      2.;
}

double born_prefactor(double Q2, double P1K1, double Mass_o0, double Mass_o1, Parameters *params) {
  double S = Q2;
  double T = Mass_o0 * Mass_o0 - 2 * P1K1;

  int aa = params->in1;
  int bb = params->in2;
  int ii = params->out1;
  int jj = params->out2;
  FI fi; // no delete ff needed
  FI *ff = &fi;
  BORN_KINEMATICS;
  struct Coupling Cw[4] = {0, 0, 0, 0};
  if (ii > 30) {
    ii = ii - 31;
    Cw[0] = params->CHSQq[jj][ii][aa];
    Cw[1] = params->gSQSQ[ii][ii];
    Cw[2] = params->CHSQq[jj][ii][aa];
    Cw[3] = params->gSQSQ[ii][ii];
  } else {
    jj = jj - 31;
    Cw[0] = params->CHSQq[ii][jj][bb];
    Cw[1] = params->gSQSQ[jj][jj];
    Cw[2] = params->CHSQq[ii][jj][bb];
    Cw[3] = params->gSQSQ[jj][jj];
  }

  //   ff->SetPropagator(params->mSQ[ii], params->mSQ[ii], params->mSQ[ii]
  //   * 1.0e-2, params->mSQ[ii] * 1.0e-2);

  ff->SetWCoupling(Cw);
  // std::cout << "born" <<
  // real(ff->LLLL+ff->RRRR)/(std::norm(params->gqq[0][0].R)) << std::endl;
  // std::norm(params->gqq[0][0].R)
  return real(ff->LLLL + ff->RRRR) * M_PI / NC / (std::norm(params->gqq[0][0].R)); // alpha_s
}

/**
 * Chiral coupling constant
 */
Tensor<ComplexType, 2> cc_chiral(Parameters *params) {
  int aa = params->in1;  // quark oder gluon
  int bb = params->in2;  // quark order gluon
  int ii = params->out1; // squark if >30 else gaugino
  int jj = params->out2; // antisquark if ii<30 else gaugino
  // std::cout << aa << " " << bb << " " << ii << " " << jj << std::endl;
  if (ii > 30) {
    ii = ii - 31;
    return {params->CHSQq[jj][ii][aa].L, params->CHSQq[jj][ii][aa].R};
  } else {
    jj = jj - 31;
    return {params->CHSQq[ii][jj][bb].L, params->CHSQq[ii][jj][bb].R};
  }
}
