// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Kinematics module.

// Notation: The incoming particle (antiparticle) has momentum pb (pa) and the
// outgoing particle (antiparticle) has momentum p2 (p1).

#ifndef FI_H_
#define FI_H_

#include "params.h"
#include <vector> // TODO change to <array> for speedup?

using namespace std;

double sqrt_lambda(const double X,const double Y,const double Z);


/**
 * Contains all Kinematic quantities.
 */
class FI {
public:
  FI(){};
  ~FI(){};

  /** 
   * Sets kinematic variables for Born processes using masses of outgoing
   * particles mi(mj) and the common Mandelstam variables.
   */
  void SetKinematic(const double mi, const double mj, const double s,
                    const double t);

  /**
   * Sets 2->3 particle kinematics using the helicity frame.
   */
  void SetKinematic_HelicityFrame(const double mi, const double mj,
                                  const double s, const double s2,
                                  const double t1, const double s1,
                                  const double phi, int flag);

  /** 
   * Sets 2->3 particle kinematics using transverse momentum.
   */
  void SetKinematic(const double mi, const double mj, const double s,
                    const double mi2, const double pt2, const double th,
                    const double ph, const int ys);

  /**
   *  Swaps pa <-> -p3
   */
  void Cross_pa_p3();
  /**
   *  Swaps pb <-> -p3
   */
  void Cross_pb_p3();
  /** 
   *  Swaps pa <-> pb
   */
  void Cross_pa_pb();


  /**
   *  reshuffles the momenta to on-shell kinematics, where the onshell particle decays in particles 2 and 3.
   */
  void SetKinematicRES23(const double mSQs, double &factor);
  /**
   *  reshuffles the momenta to on-shell kinematics, where the onshell particle decays in particles 1 and 3.
   */
  void SetKinematicRES13(const double mSQs, double &factor);

  /** 
   * Dipole kinematics for emitter A and spectator B. (old)
   */
  void SetDipKinematicA(double &x, double &fact);

  /** 
   * Dipole kinematics for emitter B and spectator A. (old)
   */
  void SetDipKinematicB(double &x, double &fact);

  // Dipole kinematics for different emitter & spectator pairs
  void SetDipKinematicAB(double &x); // Emitter A & Spectator B
  void SetDipKinematicBA(double &x); // Emitter B & Spectator A
  void SetDipKinematicA1(double &x, double &zi,
                         double &zj); // Emitter A & Spectator 1
  void SetDipKinematicB1(double &x, double &zi,
                         double &zj); // Emitter B & Spectator 1
  void SetDipKinematic1A(double &x, double &zi, double &zj, double &zplus,
                         double &zminus, double mi, double mj,
                         double mij); // Emitter 1 & Spectator A
  void SetDipKinematic1B(double &x, double &zi, double &zj, double &zplus,
                         double &zminus, double mi, double mj,
                         double mij); // Emitter 1 & Spectator B

  // Sets propagators using two mediator masses and the width.
  void SetPropagator(const double mass1, const double mass2,
                     const double width1, const double width2);

  /**
   *  Sets the breit wigner form used for on-shell subtraction.
   */
  void SetBreitWigner(double ps, double ms, double DecayWidth);

  // Sets strong coupling coefficients for SUSY-QCD. Distinguishes between L-
  // and R-type.
  void SetSCoupling(struct Coupling C[2]);

  // For squark-gaugino production at NLO
  // Simliar to SetSCoupling
  void SetLsRsCoupling(struct Coupling C[2]);

  // Sets general electroweak couplings. Four coupling coefficients for the
  // squared matrix element.
  // Distinguished between couplings to left- and right-handed fermions.
  // Complex conjugated coefficients in M*.
  // (also used for "new" processes; "W" misleading, since strong coupling is
  // involved)
  void SetWCoupling(struct Coupling C[4]);

  // like other set Couplings but with 6 couplings
  //void Set6Coupling(struct Coupling C[6]);

  // Similar to SetSCoupling, but second coupling coefficient is complex
  // conjugated.
  // (needed for instance for gaugino-squark production in LO)
  void SetBCoupling(struct Coupling C[2]);

  // Born matrix elements
  // for gaugino pair prod (or lepton pair)
  double MBss();
  double MBtt();
  double MBuu();
  double MBst();
  double MBsu();
  double MBtu();
  double MB1ss(); // epsilon coefficient of the Born matrix element in (4-2
                  // epsilon)-dim
  double MB1st();
  double MB1su();
  double MB2ss(); // epsilon^2 coefficient of the Born matrix element in (4-2
                  // epsilon)-dim
  double MBss_width();
  double MB1ss_width();
  double MB2ss_width();

  // Matrix elements for gaugino-pair production
  // Real matrix elements
  double MGss();
  double MGtt();
  double MGuu();
  double MGst();
  double MGsu();
  double MGtu();
  double MQss();
  double MQtt();
  double MQuu();
  double MQst();
  double MQsu();
  double MQtu();
  double MQBss();
  double MQBtt();
  double MQBuu();
  double MQBst();
  double MQBsu();
  double MQBtu();

  // resonant diagrams
  double MQttp();
  double MQuup();
  double MQBttp();
  double MQBuup();

  // Virtual matrix elements
  // Bubbles

  // quark self-energy
  double MVsbu1s(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu1t(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu1u(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu2s(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu2t(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu2u(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu1s(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu1t(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu1u(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu2s(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu2t(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu2u(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu1s(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu1t(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu1u(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu2s(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu2t(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu2u(double p1s, double ml1s, double ml2s, int ieps);

  // squark self-energy

  // squark gluon in loop
  double MVtbu3s(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu3t(double p1s, double ml1s, double ml2s, int ieps);
  double MVtbu3u(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu3s(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu3t(double p1s, double ml1s, double ml2s, int ieps);
  double MVubu3u(double p1s, double ml1s, double ml2s, int ieps);

  // quark gluino in loop
  double MVtbu4s(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVtbu4t(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVtbu4u(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVubu4s(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVubu4t(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVubu4u(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);

  // squark in loop (squark tadpole; only off-diagonal mixing terms do not
  // vanish)
  double MVtbu5s(double ml1s, double mp1s, double mp2s, int ieps);
  double MVtbu5t(double ml1s, double mp1s, double mp2s, int ieps);
  double MVtbu5u(double ml1s, double mp1s, double mp2s, int ieps);
  double MVubu5s(double ml1s, double mp1s, double mp2s, int ieps);
  double MVubu5t(double ml1s, double mp1s, double mp2s, int ieps);
  double MVubu5u(double ml1s, double mp1s, double mp2s, int ieps);

  // Triangles
  // gluon, quark, quark triangle (qcd triangle)
  double MVstr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVstr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVstr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // gluino, squark, squark triangle (susy qcd counterpart to qcd triangle)
  double MVstr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVstr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVstr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // type: gluino, squark, squark in loop
  double MVstr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVstr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // gluon, quark, squark triangle
  double MVttr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr1s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr1t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr1u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // gluino, quark, squark triangle
  double MVttr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr2s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr2t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr2u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // also gluino, quark, squark
  double MVttr3s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr3s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr3t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr3u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // also gluino squark quark
  double MVttr4s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr4t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVttr4u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr4s(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr4t(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);
  double MVutr4u(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                 double ml3s, int ieps);

  // Boxes
  // quark, gluon, squark, quark box
  double MVtbo1s(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVtbo1t(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVtbo1u(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo1s(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo1t(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo1u(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);

  // Boxes
  // squark, gluino, squark, quark box
  double MVtbo2s(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVtbo2t(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVtbo2u(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo2s(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo2t(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);
  double MVubo2u(double p1s, double p2s, double p3s, double p4s, double pas,
                 double pbs, double ml1s, double ml2s, double ml3s, double ml4s,
                 int ieps);

  // Counterterms
  // mass counterterms for internal squark self-energy correction.
  // squark and gluon in loop
  double MVtct3s(double p1s, double ml1s, double ml2s, int ieps);
  double MVtct3t(double p1s, double ml1s, double ml2s, int ieps);
  double MVtct3u(double p1s, double ml1s, double ml2s, int ieps);
  double MVuct3s(double p1s, double ml1s, double ml2s, int ieps);
  double MVuct3t(double p1s, double ml1s, double ml2s, int ieps);
  double MVuct3u(double p1s, double ml1s, double ml2s, int ieps);

  // quark and gluino in loop
  double MVtct4s(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVtct4t(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVtct4u(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVuct4s(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVuct4t(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);
  double MVuct4u(double p1s, double ml1s, double ml2s, double mp1s, double mp2s,
                 int ieps);

  // Matrix elements for sleptons

  // Born
  double MBssSL();

  // Real Corrections
  double MGssSL();  // real gluon emission
  double MQssSL();  // real quark emission
  double MQBssSL(); // real antiquark emission

  // Self Energies
  double MVsbu2sSL(double p1s, double ml1s, double ml2s, int ieps);
  double MVsbu1sSL(double p1s, double ml1s, double ml2s, int ieps);

  // Vertex Corrections
  double MVstr1sSL(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                   double ml3s, int ieps);
  double MVstr2sSL(double p1s, double p2s, double p3s, double ml1s, double ml2s,
                   double ml3s, int ieps);

  // Matrix elements for associated gaugino-gluino production.

  // Born
  double Mtt_GLGA();
  double Muu_GLGA();
  double Mtu_GLGA();

  // Quark self-energy
  double Mtt_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Muu_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Mtu_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Mut_QR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Mtt_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Muu_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Mtu_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps);
  double Mut_QBR_GLGA(double p1s, double ml1s, double ml2s, int ieps);

  // Gluino self-energy
  // Contribution 1 (quark-quark)
  double Mtt_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Muu_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Mtu_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Mut_GLR1_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  // Contribution 2 (gluino-gluon)
  double Mtt_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Muu_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Mtu_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);
  double Mut_GLR2_GLGA(double p1s, double ml1s, double ml2s, double mGL,
                       int ieps);

  // Squark propagator correction
  double Mtt_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, int i, int j, double mSQs);
  double Muu_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, int i, int j, double mSQs);
  double Mtu_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, int i, int j, double mSQs);
  double Mut_SQP1_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, int i, int j, double mSQs);

  double Mtt_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, double mSQs);
  double Muu_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, double mSQs);
  double Mtu_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, double mSQs);
  double Mut_SQP2_GLGA(double p1s, double ml1s, double ml2s, double mp1s,
                       double mp2s, int ieps, double mSQs);

  // Triangle 1
  double Mtt_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Muu_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V1a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtt_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Muu_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V1b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  // Triangle 2
  double Mtt_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Muu_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V2a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtt_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Muu_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V2b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  // Triangle 3
  double Mtt_V3a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtt_V3a_GLGA_BENJ(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps);

  double Muu_V3b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V3a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V3b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  // Triangle 4
  double Mtt_V4a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Muu_V4b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mtu_V4a_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  double Mut_V4b_GLGA(double p1s, double p2s, double p3s, double ml1s,
                      double ml2s, double ml3s, int ieps);
  // Box 1
  double Mtt_B1_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Muu_B1_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mtu_B1_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mut_B1_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  // Box 2
  double Mtt_B2_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Muu_B2_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mtu_B2_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mut_B2_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  // Box 3
  double Mtt_B3_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Muu_B3_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mtu_B3_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);
  double Mut_B3_GLGA(double p1s, double p2s, double p3s, double p4s, double pas,
                     double pbs, double ml1s, double ml2s, double ml3s,
                     double ml4s, int ieps);

  // real gluon emission
  double MttG_GLGA();
  double MuuG_GLGA();
  double MtuG_GLGA();
  // real quark emission
  double MttQ_GLGA();
  double MuuQ_GLGA();
  double MtuQ_GLGA();

  // real quark emission without on-shell diagrams
  double MttQ_GLGA_wos();
  double MuuQ_GLGA_wos();

  // resonant on-shell diagrams for real quark emission
  double MttQres_GLGA();
  double MuuQres_GLGA();

  // real antiquark emission
  double MttQB_GLGA();
  double MuuQB_GLGA();
  double MtuQB_GLGA();

  // real antiquark emission without on-shell contributions
  double MuuQB_GLGA_wos();
  double MttQB_GLGA_wos();

  // resonant on-shell diagrams for real antiquark emission
  double MttQBres_GLGA();
  double MuuQBres_GLGA();

  // squark production
  double Mss_SQGA1();
  double Muu_SQGA1();
  double Msu_SQGA1();

  // antisquark production
  double Mss_SQGA2();
  double Muu_SQGA2();
  double Msu_SQGA2();

  // polarized born matrix element
  void SetDipKinematicAB_pol(double &x, double &s1, double &papi, double &ztb, double &PbKi, double &PbKt1); // Emitter A/B & Spectator B/A (special for the case of contraction with polarized Born matrix element)
  void SetDipKinematicB1_pol(double &x, double &Q2, double &P1Ki, double &ztj, double &Pt1Kt1); // Emitter A/B & Spectator 1/2 (special for the case of contraction with polarized Born matrix element)
  void SetDipKinematicA1_pol(double &x, double &Q2, double &P1Ki, double &ztj, double &Pt1Kt1); // Emitter A/B & Spectator 1/2 (special for the case of contraction with polarized Born matrix element)
  vector<double> Mborn_pol(double &Q2, double &PK);
  vector<vector<double>> DIPCONNECT_in_in();
  vector<vector<double>> DIPCONNECT_in_out();
  complex<double> MG_SQGA_uIg(
    double Dn_p1mk3 ,
    double Dn_p2mk1_MQ ,
    double Dn_qmk2_MQ ,
    complex<double> Dn_p1mk2_MQ ,
    double Dn_p2mk3 ,
    double rDn_p1mk2_MQ_pow0 ,
    complex<double> rDn_p1mk2_MQ_pow1 ,
    complex<double> rDn_p1mk2_MQ_pow2 ,
    double K1Q 
  );
  double MG_SQGA_qg();
  double MG_SQGA_gg(Parameters * params);
  double MG_SQGA_gg_onshell_23(Parameters * params);

  // virtual corrections
  double Mss_qqg1_SQGA(double, double, double, double, double, double, int);
  double Mss_qqg2_SQGA(double, double, double, double, double, double, int);
  double Mss_qqg3_SQGA(double, double, double, double, double, double, int);
  double Mss_qqg4_SQGA(double, double, double, double, double, double, int);
  double Mss_qsqga1_SQGA(double, double, double, double, double, double, int);
  double Mss_qsqga2_SQGA(double, double, double, double, double, double, int);
  double Muu_qsqga1_SQGA(double, double, double, double, double, double, int);
  double Muu_qsqga2_SQGA(double, double, double, double, double, double, int);
  double Muu_gsqsq1_SQGA(double, double, double, double, double, double, int);
  double Muu_gsqsq2_SQGA(double, double, double, double, double, double, int);
  double Muu_gsqsq3_SQGA(double, double, double, double, double, double, int);
  double Muu_gsqsq4_SQGA(double, double, double, double, double, double, int);

  // Widths.
  double G_sq_2_q_gau(double m1s, double m2s, double m3s);
  double G_sq_2_q_gl(double m1s, double m2s, double m3s);
  double G_sq_2_sq_ewb(double m1s, double m2s, double m3s);

public:
  // Widths
  double prop_width[2];

  // Kinematics
  double m1, m2, m1s, m2s;
  double papb, pap1, pap2, pap3, pbp1, pbp2, pbp3, p1p2, p1p3, p2p3;

  // Propagators
  double ivs, BreitWigner;
  complex<double> ivs1s1, ivs1s2, ivs2s1, ivs2s2, ivs3v1, ivs3v2;
  double ivu3, ivt3, ivt1s1, ivu2s1, ivt1s2, ivt2s1, ivt2s2, ivu1s1, ivu1s2,
      ivu2s2, iv12;

  //  double ivu3, ivt3, iv12;

  complex<double> iv13s1, iv13s2, iv23s1, iv23s2, ivb1, ivb1s, iva1, iv13, iv23;

  // Couplings
  // Strong Couplings
  complex<double> LL, LR, RR, RL, LsLs, LsRs, RsRs, RsLs;

  // Weak Couplings
  complex<double> LLLL, LLLR, LLRL, LRLL, RLLL, LLRR, LRLR, LRRL;
  complex<double> RRRR, RRRL, RRLR, RLRR, LRRR, RRLL, RLRL, RLLR;

  // 6 Couplings
  //complex<double>   LLLLLL,   LLLLLR,   LLLLRL,   LLLRLL,   LLRLLL,   LRLLLL,   RLLLLL,   LLLLRR,   LLLRLR,   LLRLLR,   LRLLLR,   RLLLLR,   LLLRRL,   LLRLRL,   LRLLRL,   RLLLRL,   LLRRLL,   LRLRLL,   RLLRLL,   LRRLLL,   RLRLLL,   RRLLLL,     LLLRRR,   LLRLRR,   LRLLRR,   RLLLRR,   LLRRLR,   LRLRLR,   RLLRLR,   LRRLLR,   RLRLLR,   RRLLLR,   LLRRRL,   LRLRRL,   RLLRRL,   LRRLRL,   RLRLRL,   RRLLRL,   LRRRLL,   RLRRLL,   RRLRLL,   RRRLLL,   RRRRRR,   RRRRRL,   RRRRLR,   RRRLRR,   RRLRRR,   RLRRRR,   LRRRRR,   RRRRLL,   RRRLRL,   RRLRRL,   RLRRRL,   LRRRRL,   RRRLLR,   RRLRLR,   RLRRLR,   LRRRLR,   RRLLRR,   RLRLRR,   LRRLRR,   RLLRRR,   LRLRRR,   LLRRRR;
  // Loop masses
  double mul, m1l, m2l, m3l, m4l;
  double ml1, ml2, ml3, ml4;
};

#endif
