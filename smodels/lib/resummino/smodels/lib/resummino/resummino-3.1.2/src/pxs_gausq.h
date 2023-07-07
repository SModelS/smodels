#ifndef PXS_GAUSQ
#define PXS_GAUSQ

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "params.h"
#include "pdf.h"
#include "pxs.h"
#include "pxs_gausq.h"
#include "tensors.h"
#include "utils.h"

using namespace Fastor;

#include "constants.h"
#include "tensors.h"

#include "pxs_gausq_2.h"
#include "pxs_gausq_3.h"

#define BOX_INDEX                                                                                  \
  int aa = params->in1;  /* quark oder gluon                 */                                    \
  int bb = params->in2;  /* quark order gluon                */                                    \
  int ii = params->out1; /* squark if >30 else gaugino       */                                    \
  int jj = params->out2; /* antisquark if ii<30 else gaugino */                                    \
  int sq, q, ch;                                                                                   \
  if (ii > 30) {                                                                                   \
    ii = ii - 31;                                                                                  \
    sq = ii;                                                                                       \
    q = aa;                                                                                        \
    ch = jj;                                                                                       \
  } else {                                                                                         \
    jj = jj - 31;                                                                                  \
    sq = jj;                                                                                       \
    q = bb;                                                                                        \
    ch = ii;                                                                                       \
  }

#define BOX_MASS                                                                                   \
  double m1, m2;                                                                                   \
  SET_MASS;                                                                                        \
  double MQ = m1;                                                                                  \
  double MX = m2;                                                                                  \
  double MG = params->mGL;                                                                         \
  double dM2o = MX * MX - MQ * MQ;                                                                 \
  double dM2_GX = MG * MG - MX * MX;

#define BOX_KINEMATIC                                                                              \
  BOX_MASS;                                                                                        \
  double P1K2 = Q2 / 2 - P1K1;                                                                     \
  double Dn_p1mk2_MQ = dM2o - 2 * P1K2;                                                            \
  double Dn_p2mk2_MG = -2 * P1K1 - dM2_GX;

#define BOX_BASE                                                                                   \
  auto MQs = MQ * MQ;                                                                              \
  auto MUs = MQ * MQ;                                                                              \
  auto MU = MQ;                                                                                    \
  auto MXs = MX * MX;                                                                              \
  auto MGs = MG * MG;                                                                              \
  auto u = MX * MX - 2 * P1K2;                                                                     \
  auto U = MX * MX - 2 * P1K2;                                                                     \
  auto s = Q2;                                                                                     \
  auto t = MU * MU - 2 * P1K1;                                                                     \
  auto pi_ = M_PI;                                                                                 \
  auto gs = (params->g3);                                                                          \
  auto Ca = 3;                                                                                     \
  auto Cf = 4. / 3.;                                                                               \
  auto Nc = Ca;                                                                                    \
  auto dim = 4;                                                                                    \
  auto Pi = M_PI;

/** used to calculate all couplings only once at the first_run run
 * (quark/antiquark emission) */
extern bool first_run;
/** used to calculate all crossed couplings only once at the first_run run
 * (quark/antiquark emission)*/
extern bool first_run_cross;

using BoxComp = Tensor<ComplexType, 11, 2, 4>;
;
Tensor<ComplexType, 2> cc_chiral(Parameters *params);

// self energies ;
Tensor<ComplexType, 4, 4> SEq_ml_SM(double Q2, int q_I, Parameters *params);
Tensor<ComplexType, 4, 4> SEq_ml_MSSM(double Q2, int q_I, Parameters *params);
Tensor<ComplexType, 4> SEQ_gluon_squark(const double Q2, int Q_I, Parameters *params);
Tensor<ComplexType, 4> SEQ_quark_gluino(const double Q2, int Q_I, Parameters *params);
Tensor<ComplexType, 4> SEQ_squark(const double Q2, int Q_I, Parameters *params);
;
void boxcomp_SEq(Tensor<ComplexType, 2> &Lp_IJk,
                 Tensor<ComplexType, 4, 4> (*SE)(double, int, Parameters *), double Q2, double P1K1,
                 BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);
void boxcomp_SEQ(Tensor<ComplexType, 2> &Lp_IJk,
                 Tensor<ComplexType, 4> (*SE)(double, int, Parameters *), double Q2, double P1K1,
                 BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);

// vertices
Tensor<ComplexType, 2, 2, 4> VFqQx_squark_quark_gluon(double Q2, double K12, Parameters *params);
Tensor<ComplexType, 2, 2, 4> VFqQx_squark_quark_gluino(double Q2, double K12, Parameters *params);
void boxcomp_s_VFqQX(Tensor<ComplexType, 2> &Lp_IJk,
                     Tensor<ComplexType, 2, 2, 4> (*VF)(double, double, Parameters *), double Q2,
                     double P1K1, BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);
void boxcomp_u_VFqQX(Tensor<ComplexType, 2> &Lp_IJk,
                     Tensor<ComplexType, 2, 2, 4> (*VF)(double, double, Parameters *), double Q2,
                     double P1K1, BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);

Tensor<ComplexType, 2, 2, 4> VFQQg_squark_gluon(double U, Parameters *params);
Tensor<ComplexType, 2, 2, 4> VFQQg_squark_gluon_gluon(double U, Parameters *params);
Tensor<ComplexType, 2, 2, 4> VFQQg_squark_squark_gluon(double U, Parameters *params);
Tensor<ComplexType, 2, 2, 4> VFQQg_gluino_quark_quark(double U, Parameters *params);
Tensor<ComplexType, 2, 2, 4> VFQQg_gluino_gluino_quark(double U, Parameters *params);
void boxcomp_u_VFQQg(Tensor<ComplexType, 2> &Lp_IJk,
                     Tensor<ComplexType, 2, 2, 4> (*VF)(double, Parameters *), double Q2,
                     double P1K1, BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);

Tensor<ComplexType, 4, 2, 4> VFqqg_quark_gluon_gluon(double Q2, Parameters *params);
Tensor<ComplexType, 4, 2, 4> VFqqg_quark_quark_gluon(double Q2, Parameters *params);
Tensor<ComplexType, 4, 2, 4> VFqqg_squark_squark_gluino(double Q2, Parameters *params);
Tensor<ComplexType, 4, 2, 4> VFqqg_squark_gluino_gluino(double Q2, Parameters *params);
void boxcomp_s_VFqqg(Tensor<ComplexType, 2> &Lp_IJk,
                     Tensor<ComplexType, 4, 2, 4> (*VF)(double, Parameters *), double Q2,
                     double P1K1, BoxComp &Bx, double Mass_o0, double Mass_o1, Parameters *params);

// ComplexType ME_us_born(bool SS, bool UU, bool SS2, bool UU2, bool axialgauge, double Q2,
//                       double P1K1, Parameters *params, bool qqleft = 1, bool qqright = 1,
//                       bool qqgleft = 1, bool qqgright = 1, bool qQgleft = 1, bool qQgright = 1);

ComplexType ME_us_born_eps1(double Q2, double P1K1, Parameters *params);

ComplexType ME_virt_tot(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                        Parameters *params);
// SE
ComplexType ME_us_se_qq_gq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_se_qq_GQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);

ComplexType ME_us_se_QQ_gQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_se_QQ_Gq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_se_QQ_Q(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_se_QQ_g(POLE pIEPS, bool SS, bool uc, bool axial, double Q2, double P1K1,
                          Parameters *params);

// Vert

ComplexType ME_us_qqg_qgg(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_qqg_qqg(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_qqg_QGG(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_qqg_QQG(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_QQg_gQQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_QQg_Qgg(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_QQg_Qg1(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_QQg_Qg2(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_QQg_Qg(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                         Parameters *params);

ComplexType ME_us_QQg_qGG(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_QQg_Gqq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_uu_qQX_qgQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_ss_qQX_qgQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_qQX_qgQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_uu_qQX_QGq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_ss_qQX_QGq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_qQX_QGq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_gGX_Qqq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_gGX_Qqq_single(POLE pIEPS, bool SS, bool UU, bool axialgauge, int itsq, int itq,
                                 int color_flow, double Q2, double P1K1, Parameters *params);
ComplexType ME_us_gGX_qQQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);
ComplexType ME_us_gGX_qQQ_single(POLE pIEPS, bool SS, bool UU, bool axialgauge, int itsq, int itq,
                                 int color_flow, double Q2, double P1K1, Parameters *params);

// Box

ComplexType ME_us_box_gqqQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_box_qgQQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_box_qggQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);

ComplexType ME_us_box_qgQ(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                          Parameters *params);

ComplexType ME_us_box_GQQq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_box_QGqq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_box_QGGq(POLE pIEPS, bool SS, bool UU, bool axialgauge, double Q2, double P1K1,
                           Parameters *params);

// born
Tensor<ComplexType, 4> MsMBx(Tensor<ComplexType, 2> &CC_chiral, BoxComp &Bx, const double Q2,
                             const double P1K1, double Mass_o0, double Mass_o1, Parameters *params);
Tensor<ComplexType, 4> MuMBx(Tensor<ComplexType, 2> &CC_chiral, BoxComp &Bx, const double Q2,
                             const double P1K1, double Mass_o0, double Mass_o1, Parameters *params);

Tensor<ComplexType, 4> ME_virt(const double, const double, double Mass_o0, double Mass_o1,
                               Parameters *);
Tensor<ComplexType, 3> ME_born(const double, const double, double Mass_o0, double Mass_o1,
                               Parameters *);

Tensor<ComplexType, 3> born_ss_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1);
Tensor<ComplexType, 3> born_uu_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1);
Tensor<ComplexType, 3> born_su_kin_xsec(double Q2, double P1K1, double Mass_o0, double Mass_o1);
double born_prefactor(double Q2, double P1K1, double Mass_o0, double Mass_o1, Parameters *params);
#endif