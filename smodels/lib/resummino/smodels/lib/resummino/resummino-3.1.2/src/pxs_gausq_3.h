#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "clooptools.h"
#include "kinematics.h"
#include "npf.h"
#include "options.h"
#include "pdf.h"
#include "pxs.h"
//#include "pxs_gausq.h"
#include "utils.h"

//#include "pxs_gausq_2.h"

#include "constants.h"

#define NOARG

#define _FIN(res, eq) _FIN_(res, NOARG, eq)
#define _FIN_(res, eq, r)                                                                          \
  if (pIEPS == FIN) {                                                                              \
    eq;                                                                                            \
    res += r;                                                                                      \
  }

#define _EPS0_(res, eq, r)                                                                         \
  if (pIEPS == FIN) {                                                                              \
    setlambda(0);                                                                                  \
  } else if (pIEPS == IR_UV) {                                                                     \
    setuvdiv(1.0);                                                                                 \
    setlambda(-1);                                                                                 \
  } else if (pIEPS == IR1) {                                                                       \
    setuvdiv(0.0);                                                                                 \
    setlambda(-1);                                                                                 \
  } else if (pIEPS == IR2) {                                                                       \
    setuvdiv(0.0);                                                                                 \
    setlambda(-2);                                                                                 \
  } else if (pIEPS == UV) {                                                                        \
    setuvdiv(0.0);                                                                                 \
    setlambda(-1);                                                                                 \
    eq;                                                                                            \
    res -= r;                                                                                      \
    setuvdiv(1.0);                                                                                 \
    setlambda(-1);                                                                                 \
  }                                                                                                \
  eq;                                                                                              \
  res += r;

#define _EPS0(res, eq) _EPS0_(res, NOARG, eq)

#define _EPS1_(res, eq, r)                                                                         \
  if (pIEPS == FIN) {                                                                              \
    setuvdiv(1.0);                                                                                 \
    setlambda(-1);                                                                                 \
    eq;                                                                                            \
    res += r;                                                                                      \
  } else if (pIEPS == IR_UV) {                                                                     \
    setuvdiv(0.0);                                                                                 \
    setlambda(-2);                                                                                 \
    eq;                                                                                            \
    res += r;                                                                                      \
  } else if (pIEPS == IR1) {                                                                       \
    setuvdiv(0.0);                                                                                 \
    setlambda(-2);                                                                                 \
    eq;                                                                                            \
    res += r;                                                                                      \
  } else if (pIEPS == IR2) {                                                                       \
  } else if (pIEPS == UV) {                                                                        \
  }
#define _EPS1(res, eq) _EPS1_(res, NOARG, eq)

#define _EPS2_(res, eq, r)                                                                         \
  if (pIEPS == FIN) {                                                                              \
    setuvdiv(0.0);                                                                                 \
    setlambda(-2);                                                                                 \
    eq;                                                                                            \
    res += r;                                                                                      \
  } else if (pIEPS == IR_UV) {                                                                     \
  } else if (pIEPS == IR1) {                                                                       \
  } else if (pIEPS == IR2) {                                                                       \
  } else if (pIEPS == UV) {                                                                        \
  }
#define _EPS2(res, eq) _EPS2_(res, NOARG, eq)

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
  /*auto DIMD = DIM_D;*/                                                                           \
  auto Pi = M_PI;