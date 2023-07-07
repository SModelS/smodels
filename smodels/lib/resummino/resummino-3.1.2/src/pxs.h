// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Partonic cross sections at LO (Born), NLO and NLL.
//
// The functions with names ending in eps{1,2} contain the terms proportional to
// \epsilon^{1,2} in the D-dimensional matrix elements, where D=4-2\epsilon.

#ifndef PXS_H_
#define PXS_H_

#include "params.h"
#include <complex>

using namespace std;

enum POLE { FIN, IR_UV, IR2, IR1, UV };

// LO cross section.
double born_gauginos(const double, const double, Parameters *);
double born_gauginos_eps1(const double, const double, Parameters *);
double born_gauginos_eps2(const double, const double, Parameters *);

double born_sleptons(const double, const double, Parameters *);

double born_leptons(const double, const double, Parameters *);
double born_leptons_eps1(const double, const double, Parameters *);
double born_leptons_eps2(const double, const double, Parameters *);

double born_gagl(const double, const double, Parameters *); // Gaugino-gluino.

double born_gasq(const double, const double, Parameters *); // Gaugino-squark.
double born_gasq_apn(const double, const double, Parameters *); // Gaugino-squark.

// Virtual corrections.
double Virt_gauginos(const double, const double, Parameters *);
double Virt2_gauginos(const double, const double, Parameters *);

double Virt_sleptons(const double, const double, Parameters *);
double Virt2_sleptons(const double, const double, Parameters *);

double Virt_leptons(const double, const double, Parameters *);
double Virt2_leptons(const double, const double, Parameters *);

double Virt_gaugino_gluino(const double, const double, Parameters *);

double Virt_gaugino_squark(const double, const double, Parameters *);
double Virt_gaugino_squark(const double, const double, Parameters *,
                           POLE pIEPS);

// Dipole substraction for virtual corrections.
double DipI_gauginos(const double, const double, Parameters *);
double DipI_sleptons(const double, const double, Parameters *);
double DipI_leptons(const double, const double, Parameters *);
double DipI_gagl(const double, const double, Parameters *);
double DipI_gasq(const double, const double, Parameters *);
double DipI_gasq(const double, const double, Parameters *, POLE pIEPS);

// Real gaugino
double real_gluon_gauginos(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quark_gauginos(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quarkb_gauginos(const double, const double, const double,
                            const double, const double, const int,
                            Parameters *);

double real_quark_gaugino_gluino_onshell(const double, const double,
                                         const double, const double,
                                         const double, const int, Parameters *);
double real_quarkb_gaugino_gluino_onshell(const double, const double,
                                          const double, const double,
                                          const double, const int,
                                          Parameters *);

double real_quark_gaugino_gluino_onshell_13(const double, const double,
                                            const double, const double,
                                            const double, int, Parameters *);
double real_quarkb_gaugino_gluino_onshell_13(const double, const double,
                                             const double, const double,
                                             const double, int, Parameters *);

double real_quark_gaugino_gluino_onshell_23(const double, const double,
                                            const double, const double,
                                            const double, int, Parameters *);
double real_quarkb_gaugino_gluino_onshell_23(const double, const double,
                                             const double, const double,
                                             const double, int, Parameters *);

// Real gaugino gluino
double real_gluon_gaugino_gluino(const double, const double, const double,
                                 const double, const double, const int,
                                 Parameters *);
double real_quark_gaugino_gluino(const double, const double, const double,
                                 const double, const double, const int,
                                 Parameters *);
double real_quarkb_gaugino_gluino(const double, const double, const double,
                                  const double, const double, const int,
                                  Parameters *);

double real_quark_gaugino_gluino_os(const double, const double, const double,
                                    const double, const double, const int,
                                    Parameters *);
double real_quarkb_gaugino_gluino_os(const double, const double, const double,
                                     const double, const double, const int,
                                     Parameters *);

double DipGAB_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGBA_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGA1_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipGB1_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipG1A_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);
double DipG1B_GLGA(const double, const double, const double, const double,
                   const double, const int, Parameters *);

/**
 * Real correction to gaugino squark production.
 * Initial particles are a (anti-)quark and a gluon.
 * A new gluon appears in the final state.
 */
double real_gluon_gaugino_squark_qg(const double, const double, const double,
                                    const double, const double, const int,
                                    Parameters *);

/**
 * Real correction to gaugino squark production.
 * Initial particles are a gluon and a gluon.
 * A new (anit-)quark appears in the final state.
 */
double real_gluon_gaugino_squark_gg(const double, const double, const double,
                                    const double, const double, const int,
                                    Parameters *);
/**
 * Onshell resonance of the real correction to gaugino squark production.
 * Initial particles are a gluon and a gluon.
 * A new (anit-)quark appears in the final state.
 */
double real_quarkb_gaugino_squark_onshell_23(const double, const double,
                                             const double, const double,
                                             const double, const int,
                                             Parameters *);
// Real squark qqb

double real_quarkb_gaugino_squark_uub_UXub(const double S, const double M2,
                                           const double PT2, const double TH,
                                           const double PH, const int YS,
                                           Parameters *params);
double real_quarkb_gaugino_squark_uub_UXub_onshell_23(
    const double, const double, const double, const double, const double,
    const int, Parameters *);
double real_quarkb_gaugino_squark_uub_UXub_onshell_13(
    const double, const double, const double, const double, const double,
    const int, Parameters *);
double real_quark_gaugino_squark_uu_UXu(const double S, const double M2,
                                        const double PT2, const double TH,
                                        const double PH, const int YS,
                                        Parameters *params);
double real_quark_gaugino_squark_uu_UXu_onshell_23(const double, const double,
                                                   const double, const double,
                                                   const double, const int,
                                                   Parameters *);

double gausq_qg_gluon(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);
double gausq_gg_gluon(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);
double gausq_qqb_antiquark(double s, double mi2, double pt2,
                                     double th, double ph, Parameters *params);
double gausq_qq_quark(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);
double gausq_qg_gluon_minus_dip(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);
double gausq_gg_gluon_minus_dip(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);
double gausq_qqb_antiquark_minus_dip(double s, double mi2, double pt2,
                                     double th, double ph, Parameters *params);
double gausq_qq_quark_minus_dip(double s, double mi2, double pt2, double th,
                                double ph, Parameters *params);

// Real leptons
double real_gluon_leptons(const double, const double, const double,
                          const double, const double, const int, Parameters *);
double real_quark_leptons(const double, const double, const double,
                          const double, const double, const int, Parameters *);
double real_quarkb_leptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double DipGA_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipGB_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipQA_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);
double DipQB_leptons(const double, const double, const double, const double,
                     const double, const int, Parameters *);

// Real sleptons
double real_gluon_sleptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quark_sleptons(const double, const double, const double,
                           const double, const double, const int, Parameters *);
double real_quarkb_sleptons(const double, const double, const double,
                            const double, const double, const int,
                            Parameters *);
double DipGA_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipGB_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipQA_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);
double DipQB_sleptons(const double, const double, const double, const double,
                      const double, const int, Parameters *);

// Resummation
complex<double> Thadronic_nll_xs(const complex<double>, const double,
                                 const double, Parameters *);
complex<double> Thadronic_nnll_xs(const complex<double>, const double,
                                  const double, Parameters *);
complex<double> Thadronic_xs2(const complex<double>, const double, const double,
                              Parameters *);
complex<double> JtXS2(const complex<double>, const complex<double>,
                      const double, const double, Parameters *);
complex<double> PtXS(const complex<double>, const complex<double>, const double,
                     const double, Parameters *);
complex<double> JtXS(const complex<double>, const complex<double>, const double,
                     const double, Parameters *);

// sja
double sja_gagl(const double S, const double T, Parameters *params);
double sjb_gagl(const double S, const double T, Parameters *params);
double sab_gagl(const double S, const double T, Parameters *params);

double sja_gasq(const double S, const double T, Parameters *params);
double sjb_gasq(const double S, const double T, Parameters *params);
double sab_gasq(const double S, const double T, Parameters *params);
#endif
