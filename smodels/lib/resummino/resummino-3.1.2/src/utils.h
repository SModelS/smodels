// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Utils for other modules.

#ifndef UTILS_H_
#define UTILS_H_

#define RESUMMINO_VERSION "3.1.2"

#include <cassert>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector> // TODO change to <array> for speedup?

#include "clooptools.h"
#include "params.h"

using namespace std;

// Returns x^2.
static inline double pow2(double x) { return x * x; }

// Returns x^2.
static inline complex<double> pow2(complex<double> x) { return x * x; }

// Returns x^4.
static inline double pow4(double x) {
  double y = pow2(x);
  return y * y;
}

enum {
  P_NEUTRALINO_1 = 0,
  P_NEUTRALINO_2,
  P_NEUTRALINO_3,
  P_NEUTRALINO_4,

  P_CHARGINO_1,
  P_CHARGINO_2,

  P_SNEUTRINO_E = 10,
  P_SNEUTRINO_MU,
  P_SNEUTRINO_TAU,

  P_SLEPTON_1,
  P_SLEPTON_2,
  P_SLEPTON_3,
  P_SLEPTON_4,
  P_SLEPTON_5,
  P_SLEPTON_6,

  P_NEUTRINO_E = 20,
  P_NEUTRINO_MU,
  P_NEUTRINO_TAU,

  P_ELECTRON,
  P_MUON,
  P_TAU,

  P_GLUINO = 30,

  P_SQUARK_1 = 31,
  P_SQUARK_2,
  P_SQUARK_3,
  P_SQUARK_4,
  P_SQUARK_5,
  P_SQUARK_6,
  P_SQUARK_7,
  P_SQUARK_8,
  P_SQUARK_9,
  P_SQUARK_10,
  P_SQUARK_11,
  P_SQUARK_12,

  P_GLUON = 50
};

// Returns the charge of the positively charged member of each
// particle/antiparticle pair. And the one of the squarks(not antisquarks) +
// constant (1/3)
static inline int particle_charge(int particle) {
  switch (particle) {
  // Uncharged particles.
  case P_NEUTRALINO_1:
  case P_NEUTRALINO_2:
  case P_NEUTRALINO_3:
  case P_NEUTRALINO_4:

  case P_SNEUTRINO_E:
  case P_SNEUTRINO_MU:
  case P_SNEUTRINO_TAU:

  case P_NEUTRINO_E:
  case P_NEUTRINO_MU:
  case P_NEUTRINO_TAU:

  case P_GLUINO:
  case P_GLUON:

  // down type charge + 1/3
  case P_SQUARK_1:
  case P_SQUARK_2:
  case P_SQUARK_3:
  case P_SQUARK_4:
  case P_SQUARK_5:
  case P_SQUARK_6:
    return 0;
    break;

  // Charged particles.
  case P_CHARGINO_1:
  case P_CHARGINO_2:

  case P_SLEPTON_1:
  case P_SLEPTON_2:
  case P_SLEPTON_3:
  case P_SLEPTON_4:
  case P_SLEPTON_5:
  case P_SLEPTON_6:

  case P_ELECTRON:
  case P_MUON:
  case P_TAU:

  // up type charge +1/3
  case P_SQUARK_7:
  case P_SQUARK_8:
  case P_SQUARK_9:
  case P_SQUARK_10:
  case P_SQUARK_11:
  case P_SQUARK_12:
    return 1;
    break;

  // Charge unknown.
  default:
    fprintf(stderr, "error: Internal error. Particle charge unknown.\n");
    exit(1);
  }
}

// sqrt of kaellen function.
static inline double kln(double x, double y, double z) {
  return sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2) - 2.0 * x * y - 2.0 * y * z - 2.0 * z * x);
}

// kaellen function.
static inline double Rkln(double x, double y, double z) {
  return pow(x, 2) + pow(y, 2) + pow(z, 2) - 2.0 * x * y - 2.0 * y * z - 2.0 * z * x;
}

// universal G-function.
static inline double G(double x, double y, double z, double u, double v, double w) {
  return pow(v, 2) * w + v * pow(w, 2) + u * v * x + u * w * y + pow(x, 2) * y + x * pow(y, 2) +
         pow(u, 2) * z + w * x * z + w * y * z - u * (v + w + x + y) * z + u * pow(z, 2) -
         x * y * (u + v + w + z) - v * w * (u + x + y + z);
}

// Returns the charge of a quark plus an arbitrary constant, which is always the
// same (actually it's 1/3).
// Note: The constant does not matter because this function only appears in
// pairs +/-, so it cancels out, but in this way this is more efficient (and
// this is used at every integration point).
static inline int quark_charge_plus_constant(int quark) { return quark < 3 ? 0 : 1; }

// first row are the down type (s)quarks, second row are the up type ones
static const int squark_type[][7] = {{0, 3, 1, 4, 2, 5, -1}, {6, 9, 7, 10, 8, 11, -1}};
static const int left_right[7] = {0, 1, 0, 1, 0, 1, -1};
static const int down_squarks[] = {0, 3, 1, 4, 2, 5};
static const int down_l_squarks[] = {0, 1, 2};
static const int down_r_squarks[] = {3, 4, 5};
static const int up_squarks[] = {6, 9, 7, 10, 8, 11};
static const int up_l_squarks[] = {6, 7, 8};
static const int up_r_squarks[] = {9, 10, 11};

static const bool is_right_squark(int sq) {
  return sq == 3 || sq == 4 || sq == 5 || 9 == sq || 10 == sq || 11 == sq;
}
static const bool is_left_squark(int sq) { return !is_right_squark(sq); }

static const bool is_up_quark(int q) { return q == 3 || q == 4 || q == 5; }
static const bool is_down_quark(int q) { return q == 0 || q == 1 || q == 2; }
static const bool is_up_squark(int q) {
  return q == 6 || q == 9 || q == 7 || q == 10 || q == 8 || q == 11;
}
static const bool is_down_squark(int q) {
  return q == 0 || q == 1 || q == 2 || q == 3 || q == 4 || q == 5;
}

static const bool is_chargino(int ch) { return ch == 4 || ch == 5; }
// first row are the down type (s)quarks, second row are the up type ones
static const int quark_type[][4] = {{0, 1, 2, -1}, {3, 4, 5, -1}};
// -1 is a placeholder; 0,1 are gamma and Z, 2 is W+-
static const int vboson_type[][2] = {{0, 1}, {2, -1}};

typedef enum { CHANNEL_S = 0, CHANNEL_T, CHANNEL_U } ChannelName;

// this function returns the propagator charge (-1,0,1)
// REMARK: -1 for t-channel not possible
//         +1 for u-channel not possible
static inline int propagator_charge(int in1, int in2, int out1, int out2, ChannelName channel) {
  if (channel == CHANNEL_T) {
    return quark_charge_plus_constant(in1) - particle_charge(out1);
  } else if (channel == CHANNEL_U) {
    return quark_charge_plus_constant(in1) + particle_charge(out2);
  } else if (channel == CHANNEL_S) {
    return quark_charge_plus_constant(in1) - quark_charge_plus_constant(in2);
  }
}

static inline bool is_gaugino_gluino(int out1, int out2) {
  return ((out1 == 30 && out2 < 10) || (out2 == 30 && out1 < 10));
}

static inline bool is_gaugino_gaugino(int out1, int out2) { return (out1 < 10 && out2 < 10); }

static inline bool is_lepton_lepton(int out1, int out2) {
  return (out1 >= 20 && out1 < 30 && out2 >= 20 && out2 < 30);
}

static inline bool is_slepton_slepton(int out1, int out2) {
  return (out1 >= 10 && out1 < 20 && out2 >= 10 && out2 < 20);
}

static inline bool is_squark_gaugino(int out1, int out2) {
  return (out1 >= 31 && out2 < 10 || out2 >= 31 && out1 < 10);
}

static inline bool is_coupling_null(struct Coupling *C, int N) {
  for (int i0 = 0; i0 < N; i0++) {
    if (C[i0].L == 0.0 && C[i0].R == 0.0) {
      return true;
    }
  }
  return false;
}

static inline int iabs(int x) { return x < 0 ? -x : x; }

static inline bool is_pdg_chargino(int pdg_code) {
  switch (pdg_code) {
  case 1000024:  // ~chi_1+
  case 1000037:  // ~chi_2+
  case -1000024: // ~chi_1+
  case -1000037: // ~chi_2+
    return 1;
  default:
    return 0;
  }
}
static inline bool is_pdg_neutralino(int pdg_code) {
  switch (pdg_code) {
  case 1000022:  // ~chi_10
  case 1000023:  // ~chi_20
  case 1000025:  // ~chi_30
  case 1000035:  // ~chi_40
  case -1000022: // ~chi_10
  case -1000023: // ~chi_20
  case -1000025: // ~chi_30
  case -1000035: // ~chi_40
    return 1;
  default:
    return 0;
  }
}
static inline bool is_pdg_gaugino(int pdg_code) {
  return is_pdg_chargino(pdg_code) || is_pdg_neutralino(pdg_code);
}

static inline bool is_pdg_up_type_squark(int pdg_code) {
  switch (pdg_code) {
  case -1000002: // ~u_L*
  case -2000002: // ~u_R*
  case -1000004: // ~c_L*
  case -2000004: // ~c_R*
  case -1000006: // ~t_1*
  case -2000006: // ~t_2*
  case 1000002:  // ~u_L
  case 2000002:  // ~u_R
  case 1000004:  // ~c_L
  case 2000004:  // ~c_R
  case 1000006:  // ~t_1
  case 2000006:  // ~t_2
    return 1;
  default:
    return 0;
  }
}
static inline bool is_pdg_down_type_squark(int pdg_code) {
  switch (pdg_code) {
  case 1000001:  // ~d_L
  case 2000001:  // ~d_R
  case 1000003:  // ~s_L
  case 2000003:  // ~s_R
  case 1000005:  // ~b_1
  case 2000005:  // ~b_2
  case -1000001: // ~d_L*
  case -2000001: // ~d_R*
  case -1000003: // ~s_L*
  case -2000003: // ~s_R*
  case -1000005: // ~b_1*
  case -2000005: // ~b_2*
    return 1;
  default:
    return 0;
  }
}
static inline bool is_pdg_squark(int pdg_code) {
  return is_pdg_up_type_squark(pdg_code) || is_pdg_down_type_squark(pdg_code);
}

static inline int pdg_particle_charge(int pdg_code) {
  switch (pdg_code) {
  case 1000012:  // ~nu_eL
  case 1000014:  // ~nu_muL
  case 1000016:  // ~nu_tau
  case 12:       // nu_eL
  case 14:       // nu_muL
  case 16:       // nu_tauL
  case 1000022:  // ~chi_10
  case 1000023:  // ~chi_20
  case 1000025:  // ~chi_30
  case 1000035:  // ~chi_40
  case 1000021:  // ~g
  case -1000012: // ~nu_eL
  case -1000014: // ~nu_muL
  case -1000016: // ~nu_tau
  case -12:      // nu_eL
  case -14:      // nu_muL
  case -16:      // nu_tauL
  case -1000022: // ~chi_10
  case -1000023: // ~chi_20
  case -1000025: // ~chi_30
  case -1000035: // ~chi_40
  case -1000021: // ~g
  case 1000001:  // ~d_L
  case 2000001:  // ~d_R
  case 1000003:  // ~s_L
  case 2000003:  // ~s_R
  case 1000005:  // ~b_1
  case 2000005:  // ~b_2
  case -1000001: // ~d_L*
  case -2000001: // ~d_R*
  case -1000003: // ~s_L*
  case -2000003: // ~s_R*
  case -1000005: // ~b_1*
  case -2000005: // ~b_2*
    return 0;

  case 1000011:  // ~e_L
  case 2000011:  // ~e_R
  case 1000013:  // ~mu_L
  case 2000013:  // ~mu_R
  case 1000015:  // ~tau_1
  case 2000015:  // ~tau_2
  case 11:       // e
  case 13:       // mu
  case 15:       // tau
  case -1000024: // ~chi_1-
  case -1000037: // ~chi_2-
  case -1000002: // ~u_L*
  case -2000002: // ~u_R*
  case -1000004: // ~c_L*
  case -2000004: // ~c_R*
  case -1000006: // ~t_1*
  case -2000006: // ~t_2*
    return -1;
  case -1000011: // ~e_L
  case -2000011: // ~e_R
  case -1000013: // ~mu_L
  case -2000013: // ~mu_R
  case -1000015: // ~tau_1
  case -2000015: // ~tau_2
  case -11:      // e
  case -13:      // mu
  case -15:      // tau
  case 1000024:  // ~chi_1+
  case 1000037:  // ~chi_2+
  case 1000002:  // ~u_L
  case 2000002:  // ~u_R
  case 1000004:  // ~c_L
  case 2000004:  // ~c_R
  case 1000006:  // ~t_1
  case 2000006:  // ~t_2
    return 1;
  default:
    fprintf(stderr, "error: Wrong PDG code.\n");
    exit(1);
  }
}

// Tests the charge conservation.
static inline bool is_charge_conserved(int in1, int in2, int out1, int out2) {
  // in2 and out2 are the negatively charged counterparts.
  // Notice that the arbitrary constant (which is actually 1/3) in
  // quark_charge_plus_constant cancels out.
  if (is_squark_gaugino(out1, out2)) {
    // Unfortunately resummino does not differ between chi+ and chi- internally
    // so only abs is checked here. The check if the sign is actually correct
    // happens in set_particles in resummino.cc for the squark gaugino process
    if (out1 > 30) {
      return (abs(quark_charge_plus_constant(in1) - particle_charge(out1)) ==
              particle_charge(out2));

    } else {
      return (abs(-quark_charge_plus_constant(in2) + particle_charge(out2)) ==
              particle_charge(out1));
    }
  } else {
    return quark_charge_plus_constant(in1) - quark_charge_plus_constant(in2) ==
           particle_charge(out1) - particle_charge(out2);
  }
}
template <typename T> T scalar_product(const vector<T> &vi, const vector<T> &vj) {
  // check dimensions
  assert(vi.size() == vj.size());
  // sunm over i and j
  T result = T(0);
  for (int i = 0; i < vi.size(); ++i)
    result += vi[i] * vj[i];
  // if(result != 0 && result > 0.1) {
  // for (int i =0; i < vi.size();++i) {
  // 	  std::cout << vi[i]*vj[i]/result << ", ";
  //  }
  //  std::cout << " of " << result << std::endl;
  //}
  return result;
}
// Lorents contraction
template <typename T>
T einsum_i_ij_j(const vector<T> &vi, const vector<vector<T>> &vij, const vector<T> &vj) {
  // check dimensions
  assert(vi.size() == vij.size());
  assert(vj.size() == vij[0].size());
  // sum over i and j
  vector<T> pre_result;
  for (int j = 0; j < vj.size(); ++j) {
    T tmp = T(0);
    for (int i = 0; i < vi.size(); ++i) {
      tmp += vij[i][j] * vi[i];
    }
    pre_result.push_back(tmp);
  }
  auto tmp = scalar_product(pre_result, vj);
  /*
  if(tmp > 0.1) {
        for (int j =0; j < vj.size();++j) {
    std::cout << "\t";
                for (int i =0; i < vi.size();++i) {
      std::cout << std::setw(10) << vi[i] * vij[i][j]*vj[j]/tmp << "\t";
    }
    std::cout << endl;
  }
  std::cout << " of " << tmp << std::endl;
  }
  //*/
  return tmp;
}

// Multiplication of vector components by scalar
template <typename T> void scale(std::vector<T> &to_be_scaled, T factor) {
  for (int i = 0; i < to_be_scaled.size(); i++)
    to_be_scaled[i] *= factor;
}

// Macros
#define SET_MASS                                                                                   \
  if (params->out1 < 10 && params->out2 < 10) {                                                    \
    m1 = params->mCH[params->out1];                                                                \
    m2 = params->mCH[params->out2];                                                                \
  } else if (params->out1 >= 10 && params->out1 < 20) {                                            \
    m1 = params->mSL[params->out1 - 10];                                                           \
    m2 = params->mSL[params->out2 - 10];                                                           \
  } else if (params->out1 < 10 && params->out2 > 30) {                                             \
    m2 = params->mCH[params->out1];                                                                \
    m1 = params->mSQ[params->out2 - 31];                                                           \
  } else if (params->out2 < 10 && params->out1 > 30) {                                             \
    m2 = params->mCH[params->out2];                                                                \
    m1 = params->mSQ[params->out1 - 31];                                                           \
  } else if (params->out1 < 10 && params->out2 == 30) {                                            \
    m2 = params->mCH[params->out1];                                                                \
    m1 = params->mGL;                                                                              \
  } else if (params->out2 < 10 && params->out1 == 30) {                                            \
    m2 = params->mCH[params->out2];                                                                \
    m1 = params->mGL;                                                                              \
  } else if (params->out1 >= 20 && params->out1 < 30) {                                            \
    m1 = params->ml[params->out1 - 20];                                                            \
    m2 = params->ml[params->out2 - 20];                                                            \
  }

#define BORN_KINEMATICS                                                                            \
  if (ii < 10 && jj < 10) {                                                                        \
    ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, T);                                      \
  } else if (ii >= 10 && ii < 20) {                                                                \
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, T);                            \
  } else if (ii >= 20 && ii < 30) {                                                                \
    ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, T);                              \
  } else if (ii == 30 && jj < 10) {                                                                \
    ff->SetKinematic(params->mGL, params->mCH[jj], S, T);                                          \
  } else if (jj == 30 && ii < 10) {                                                                \
    ff->SetKinematic(params->mCH[ii], params->mGL, S, T);                                          \
  } else if (ii >= 31 && jj < 10) {                                                                \
    ff->SetKinematic(params->mSQ[ii - 31], params->mCH[jj], S, T);                                 \
  } else if (jj >= 31 && ii < 10) {                                                                \
    ff->SetKinematic(params->mSQ[jj - 31], params->mCH[ii], S, T);                                 \
  }

#define JET_KINEMATICS                                                                             \
  if (ii < 10 && jj < 10) {                                                                        \
    ff->SetKinematic(params->mCH[ii], params->mCH[jj], S, M2, PT2, TH, PH, YS);                    \
  } else if (ii >= 10 && ii < 20) {                                                                \
    ff->SetKinematic(params->mSL[ii - 10], params->mSL[jj - 10], S, M2, PT2, TH, PH, YS);          \
  } else if (ii >= 20 && ii < 30) {                                                                \
    ff->SetKinematic(params->ml[ii - 20], params->ml[jj - 20], S, M2, PT2, TH, PH, YS);            \
  } else if (ii == 30 && jj < 10) {                                                                \
    ff->SetKinematic(params->mGL, params->mCH[jj], S, M2, PT2, TH, PH, YS);                        \
  } else if (jj == 30 && ii < 10) {                                                                \
    ff->SetKinematic(params->mCH[ii], params->mGL, S, M2, PT2, TH, PH, YS);                        \
  } else if (ii >= 31 && jj < 10) {                                                                \
    ff->SetKinematic(params->mSQ[ii - 31], params->mCH[jj], S, M2, PT2, TH, PH, YS);               \
  } else if (jj >= 31 && ii < 10) {                                                                \
    ff->SetKinematic(params->mSQ[jj - 31], params->mCH[ii], S, M2, PT2, TH, PH, YS);               \
  }

#endif
