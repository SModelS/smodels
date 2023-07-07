// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Header file for dipole functions and collinear remainder.

#include "params.h"
#include "tensors.h"
// TODO change to <array> for speedup?
#include <vector>

// different function types for the collinear remainder
typedef enum {
  FUNCTION_ALL_WDELTA, // all parts except the one proportional to the
                       // delta-disttribution
  FUNCTION_DELTA,      // part proportional to delta-dist.
  FUNCTION_PLUS,       // plus-dist. part
  FUNCTION_XPLUS       // X-plust-dist part
} FunctionType;

// different parton types
typedef enum {
  PARTON_QUARK,
  PARTON_ANTIQUARK,
  PARTON_GLUON,
  PARTON_GLUINO,
  PARTON_SQUARK,
  PARTON_NONE
} PartonType;

// different emission types
typedef enum { INITIAL_AND_FINAL, INITIAL } EmissionType;

// different emitter spectator pairs
typedef enum {
  INITIAL_INITIAL,
  INITIAL_FINAL,
  FINAL_INITIAL,
  FINAL_FINAL
} DipoleType;
double Kappa_aap_j(PartonType a, PartonType ap, PartonType j, double mj,
                   double sja, FunctionType function_type, double x, double z);
double J_a_ij(PartonType i, PartonType j, double muj,
              FunctionType function_type, double x, double z);
// Collinear remainder: Insertion operators K and P (actually only P is
// factorization scale dependent!).
double K_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z, double mj,
                      double sja, double sab, double alphas,
                      EmissionType emission);
double P_bold_aap_b_j(PartonType a, PartonType ap, PartonType b, PartonType j,
                      FunctionType function_type, double x, double z,
                      double sja, double sab, double mufs, double alphas,
                      EmissionType emission);
// combined P + K
double P_plus_K_bold_aap_b_j(PartonType a, PartonType ap, PartonType b,
                             PartonType j, FunctionType function_type, double x,
                             double z, double mj, double sja, double sab,
                             double mufs, double alphas, EmissionType emission);
// combined P + K with plus distributuion evaluation if diagonal else
// FUNCTION_ALL_WDELTA
double P_plus_K_bold(PartonType a, PartonType ap, PartonType b, PartonType j,
                     double x, double z, double mj, double sja, double sjac,
                     double sab, double mufs, double alphas,
                     EmissionType emission, double born, double bornc,
                     double djac, double djacplus, double djacdelta);

// Universal dipoles for 2->3 processes
double D_ai_b(PartonType a, PartonType i, double papi, double papb, double pbpi,
              double x, double alpha, double color);
vector<double> D_ai_b_vector(PartonType a, PartonType i, double papi, double x,
                             double alpha, double color, double ztb,
                             double pipb);
double D_ai_j(PartonType a, PartonType i, double zi, double zj, double papi,
              double pjpi, double x, double alpha, double color);
vector<double> D_ai_j_vector(PartonType a, PartonType i, double papi, double x,
                             double alpha, double color, double ztj,
                             double pipj);
double D_ij_a(PartonType i, PartonType j, double mi, double mj, double mij,
              double pjpi, double x, double alpha, double zj, double mQ,
              double zi, double zplus, double zminus, double color);

// "Dipoles" * born used in pxs_*_real.cc
double Dip_GLGA(DipoleType Emitter_Spectator, PartonType emitter,
                PartonType spectator, PartonType emitted_parton, const double S,
                const double M2, const double PT2, const double TH,
                const double PH, const int YS, Parameters *params);

double Dip_SLEPTONS(DipoleType Emitter_Spectator, PartonType emitter,
                    PartonType spectator, PartonType emitted_parton,
                    const double S, const double M2, const double PT2,
                    const double TH, const double PH, const int YS,
                    Parameters *params);

double Dip_GAUGINOS(DipoleType Emitter_Spectator, PartonType emitter,
                    PartonType spectator, PartonType emitted_parton,
                    const double S, const double M2, const double PT2,
                    const double TH, const double PH, const int YS,
                    Parameters *params);

double Dip_SQGA(DipoleType Emitter_Spectator, PartonType emitter,
                PartonType spectator, PartonType emitted_parton, const double S,
                const double M2, const double PT2, const double TH,
                const double PH, const int YS, Parameters *params,
                const int Opt_eo = 0);

Tensor<ComplexType, 4> I_jk(PartonType j, PartonType k, double s_jk, double m_j,
                            double m_k, double T_jT_k, Parameters *params);