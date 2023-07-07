/******************************************************************************
 *                      Code generated with sympy 1.6.2                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                This file is part of 'SunderCalc_Resummino'                 *
 ******************************************************************************/
//#include "ddb_UXub_1323.h"
#include <complex>
#include "ddb_UXub.h"

inline double M_ddb_UXub_1323_0(double K2K3, double MDi1, double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K3, double QK2, double s, double s1) {

   double M_ddb_UXub_1323_0_result;
   M_ddb_UXub_1323_0_result = ME3*MG*MX*(8.0*ME1*ME1ME3*std::pow(P1K3, 2.0)*1.0/(s*s1) + 2.0*ME1ME3*1.0/s1*(2.0*K2K3*ME1 - 1.0*ME1*std::pow(MQ, 2.0) + ME1*std::pow(MX, 2.0) - 4.0*ME1*P1K3 - 2.0*ME1*QK2 + ME1*s) + ME2*ME2ME3 - 1.0*1.0/s*(-1.0*std::pow(MDi1, 2.0)*ME2*ME2ME3 - 4.0*ME1*ME1ME3*P1K3 + ME2*ME2ME3*std::pow(MX, 2.0) + 4.0*ME2*ME2ME3*P1K3));
   return M_ddb_UXub_1323_0_result;

}

inline double M_ddb_UXub_1323_1(double K2K3, double MDi1, double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K3, double dM2o, double s, double s1) {

   double M_ddb_UXub_1323_1_result;
   M_ddb_UXub_1323_1_result = -1.0*ME3*(-2.0*ME1ME3*P1K3*1.0/(s*s1)*(-1.0*std::pow(MDi1, 2.0)*ME1*std::pow(MG, 2.0) + 2.0*std::pow(MDi1, 2.0)*ME1*std::pow(MQ, 2.0) - 1.0*std::pow(MDi1, 2.0)*ME1*std::pow(MX, 2.0) - 1.0*ME1*std::pow(MG, 2.0)*std::pow(MQ, 2.0) + 2.0*ME1*std::pow(MG, 2.0)*P1K3 + ME1*std::pow(MG, 4.0) + 2.0*ME1*std::pow(MX, 2.0)*P1K3 + ME1*std::pow(MX, 2.0)*dM2o) + ME1ME3*1.0/s1*(4.0*K2K3*std::pow(MDi1, 2.0)*ME1 - 4.0*K2K3*ME1*std::pow(MX, 2.0) - 8.0*K2K3*ME1*P1K3 + 4.0*std::pow(K2K3, 2.0)*ME1 - 2.0*std::pow(MDi1, 2.0)*ME1*std::pow(MX, 2.0) - 4.0*std::pow(MDi1, 2.0)*ME1*P1K3 + std::pow(MDi1, 4.0)*ME1 + 2.0*ME1*std::pow(MG, 2.0)*P1K3 + 6.0*ME1*std::pow(MX, 2.0)*P1K3 + ME1*std::pow(MX, 4.0) + 4.0*ME1*std::pow(P1K3, 2.0)) + ME2ME3*(2.0*K2K3*ME2 + std::pow(MDi1, 2.0)*ME2 - 4.0*ME2*P1K3 - 1.0*ME2*dM2o) + 1.0/s*(-2.0*K2K3*std::pow(MDi1, 2.0)*ME1*ME1ME3 + 2.0*K2K3*std::pow(MDi1, 2.0)*ME2*ME2ME3 + 2.0*K2K3*ME1*ME1ME3*std::pow(MX, 2.0) - 2.0*K2K3*ME2*ME2ME3*std::pow(MX, 2.0) + 2.0*std::pow(MDi1, 2.0)*ME1*ME1ME3*std::pow(MX, 2.0) + 2.0*std::pow(MDi1, 2.0)*ME1*ME1ME3*P1K3 - 1.0*std::pow(MDi1, 2.0)*ME2*ME2ME3*std::pow(MG, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*ME2ME3*std::pow(MX, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*ME2ME3*P1K3 - 1.0*std::pow(MDi1, 4.0)*ME1*ME1ME3 + std::pow(MDi1, 4.0)*ME2*ME2ME3 - 2.0*ME1*ME1ME3*std::pow(MG, 2.0)*P1K3 - 4.0*ME1*ME1ME3*std::pow(MX, 2.0)*P1K3 - 1.0*ME1*ME1ME3*std::pow(MX, 4.0) + ME2*ME2ME3*std::pow(MG, 2.0)*std::pow(MX, 2.0) + 2.0*ME2*ME2ME3*std::pow(MG, 2.0)*P1K3 + 4.0*ME2*ME2ME3*std::pow(MX, 2.0)*P1K3 + ME2*ME2ME3*std::pow(MX, 4.0)));
   return M_ddb_UXub_1323_1_result;

}

inline double M_ddb_UXub_1323_2(double K2K3, double MDi1, double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K2, double P1K3, double s, double s1) {

   double M_ddb_UXub_1323_2_result;
   M_ddb_UXub_1323_2_result = ME3*(-2.0*ME1ME3*P1K3*1.0/(s*s1)*(ME1*std::pow(MG, 2.0) - 2.0*ME1*std::pow(MQ, 2.0) + ME1*std::pow(MX, 2.0)) + ME1ME3*1.0/s1*(-4.0*K2K3*ME1 - 1.0*std::pow(MDi1, 2.0)*ME1 + ME1*std::pow(MX, 2.0) + 2.0*ME1*P1K2 + 4.0*ME1*P1K3) - 1.0*ME2*ME2ME3 + 1.0/s*(2.0*K2K3*ME1*ME1ME3 - 2.0*K2K3*ME2*ME2ME3 + std::pow(MDi1, 2.0)*ME1*ME1ME3 - 1.0*std::pow(MDi1, 2.0)*ME2*ME2ME3 - 1.0*ME1*ME1ME3*std::pow(MX, 2.0) - 2.0*ME1*ME1ME3*P1K2 - 2.0*ME1*ME1ME3*P1K3 + ME2*ME2ME3*std::pow(MG, 2.0) + ME2*ME2ME3*std::pow(MX, 2.0) + 2.0*ME2*ME2ME3*P1K2 + 2.0*ME2*ME2ME3*P1K3));
   return M_ddb_UXub_1323_2_result;

}

inline double M_ddb_UXub_1323_3(double MDi1, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K3, double s) {

   double M_ddb_UXub_1323_3_result;
   M_ddb_UXub_1323_3_result = -1.0*ME2ME3*ME3*MG*MX*(-1.0*ME2*std::pow(MQ, 2.0) + ME2*std::pow(MX, 2.0) + 2.0*ME2*P1K3 + 1.0/s*(-1.0*std::pow(MDi1, 2.0)*ME2*std::pow(MG, 2.0) + std::pow(MDi1, 2.0)*ME2*std::pow(MX, 2.0) + 4.0*std::pow(MDi1, 2.0)*ME2*P1K3 + ME2*std::pow(MG, 2.0)*std::pow(MX, 2.0) + 2.0*ME2*std::pow(MG, 2.0)*P1K3 - 6.0*ME2*std::pow(MX, 2.0)*P1K3 - 1.0*ME2*std::pow(MX, 4.0) - 8.0*ME2*std::pow(P1K3, 2.0)));
   return M_ddb_UXub_1323_3_result;

}

inline double M_ddb_UXub_1323_4(double MDi1, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_1323_4_result;
   M_ddb_UXub_1323_4_result = -1.0*ME2ME3*ME3*(-2.0*std::pow(MDi1, 2.0)*ME2*P1K3 - 1.0*std::pow(MDi1, 2.0)*ME2*dM2o + 4.0*ME2*P1K3*dM2o + 4.0*ME2*std::pow(P1K3, 2.0) + ME2*std::pow(dM2o, 2.0) - 1.0*1.0/s*(std::pow(MDi1, 2.0)*ME2*std::pow(MG, 2.0)*std::pow(MQ, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*std::pow(MG, 2.0)*std::pow(MX, 2.0) - 4.0*std::pow(MDi1, 2.0)*ME2*std::pow(MG, 2.0)*P1K3 + 3.0*std::pow(MDi1, 2.0)*ME2*std::pow(MQ, 2.0)*std::pow(MX, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*std::pow(MX, 4.0) - 4.0*std::pow(MDi1, 2.0)*ME2*P1K3*dM2o + std::pow(MDi1, 4.0)*ME2*std::pow(MG, 2.0) - 2.0*std::pow(MDi1, 4.0)*ME2*std::pow(MQ, 2.0) + std::pow(MDi1, 4.0)*ME2*std::pow(MX, 2.0) - 2.0*ME2*std::pow(MG, 2.0)*std::pow(MQ, 2.0)*P1K3 + 4.0*ME2*std::pow(MG, 2.0)*std::pow(MX, 2.0)*P1K3 + ME2*std::pow(MG, 2.0)*std::pow(MX, 2.0)*dM2o + 4.0*ME2*std::pow(MG, 2.0)*std::pow(P1K3, 2.0) - 2.0*ME2*std::pow(MQ, 2.0)*std::pow(MX, 2.0)*P1K3 + 4.0*ME2*std::pow(MX, 2.0)*std::pow(P1K3, 2.0) + 4.0*ME2*std::pow(MX, 4.0)*P1K3 + ME2*std::pow(MX, 4.0)*dM2o));
   return M_ddb_UXub_1323_4_result;

}

inline double M_ddb_UXub_1323_5(double MDi1, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double P1K2, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_1323_5_result;
   M_ddb_UXub_1323_5_result = ME2ME3*ME3*(2.0*ME2*P1K3 + ME2*dM2o + 1.0/s*(std::pow(MDi1, 2.0)*ME2*std::pow(MG, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*std::pow(MQ, 2.0) + std::pow(MDi1, 2.0)*ME2*std::pow(MX, 2.0) - 2.0*ME2*std::pow(MG, 2.0)*P1K2 - 4.0*ME2*std::pow(MG, 2.0)*P1K3 - 1.0*ME2*std::pow(MG, 2.0)*dM2o + 4.0*ME2*std::pow(MQ, 2.0)*P1K2 - 2.0*ME2*std::pow(MX, 2.0)*P1K2 - 1.0*ME2*std::pow(MX, 2.0)*dM2o - 4.0*ME2*P1K3*dM2o));
   return M_ddb_UXub_1323_5_result;

}

inline double M_ddb_UXub_1323_6(double ME2, double ME2ME3, double ME3, double MG, double MX, double P1K3, double s) {

   double M_ddb_UXub_1323_6_result;
   M_ddb_UXub_1323_6_result = ME2ME3*ME3*MG*MX*1.0/s*(ME2*std::pow(MG, 2.0) - 1.0*ME2*std::pow(MX, 2.0) - 4.0*ME2*P1K3);
   return M_ddb_UXub_1323_6_result;

}

inline double M_ddb_UXub_1323_7(double ME2, double ME2ME3, double ME3, double MG, double MX, double s) {

   double M_ddb_UXub_1323_7_result;
   M_ddb_UXub_1323_7_result = ME2*ME2ME3*ME3*MG*MX*1.0/s;
   return M_ddb_UXub_1323_7_result;

}

inline double M_ddb_UXub_1323_8(double K2K3, double MDi1, double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double MG, double MX, double P1K3, double s, double s1) {

   double M_ddb_UXub_1323_8_result;
   M_ddb_UXub_1323_8_result = -1.0*ME3*1.0/s*(-2.0*ME1ME3*P1K3*1.0/s1*(-2.0*K2K3*ME1 - 1.0*std::pow(MDi1, 2.0)*ME1 + ME1*std::pow(MG, 2.0) + 2.0*ME1*P1K3) + ME2ME3*(-1.0*std::pow(MDi1, 2.0)*ME2 + ME2*std::pow(MX, 2.0) + 2.0*ME2*P1K3));
   return M_ddb_UXub_1323_8_result;

}

inline double M_ddb_UXub_1323_9(double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double P1K3, double s, double s1) {

   double M_ddb_UXub_1323_9_result;
   M_ddb_UXub_1323_9_result = ME3*1.0/s*(-2.0*ME1*ME1ME3*P1K3*1.0/s1 + ME2*ME2ME3);
   return M_ddb_UXub_1323_9_result;

}

inline double M_ddb_UXub_1323_10(double MDi1, double ME2, double ME2ME3, double ME3, double MG, double MX, double P1K3, double s) {

   double M_ddb_UXub_1323_10_result;
   M_ddb_UXub_1323_10_result = -1.0*ME2ME3*ME3*MG*MX*1.0/s*(-1.0*std::pow(MDi1, 2.0)*ME2 + ME2*std::pow(MX, 2.0) + 2.0*ME2*P1K3);
   return M_ddb_UXub_1323_10_result;

}

inline double M_ddb_UXub_1323_11(double MDi1, double ME2, double ME2ME3, double ME3, double MQ, double MX, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_1323_11_result;
   M_ddb_UXub_1323_11_result = ME2ME3*ME3*1.0/s*(std::pow(MDi1, 2.0)*ME2*std::pow(MQ, 2.0) - 2.0*std::pow(MDi1, 2.0)*ME2*std::pow(MX, 2.0) - 4.0*std::pow(MDi1, 2.0)*ME2*P1K3 + std::pow(MDi1, 4.0)*ME2 - 2.0*ME2*std::pow(MQ, 2.0)*P1K3 + 4.0*ME2*std::pow(MX, 2.0)*P1K3 + ME2*std::pow(MX, 2.0)*dM2o + 4.0*ME2*std::pow(P1K3, 2.0));
   return M_ddb_UXub_1323_11_result;

}

inline double M_ddb_UXub_1323_12(double ME2, double ME2ME3, double ME3, double MG, double MX, double s) {

   double M_ddb_UXub_1323_12_result;
   M_ddb_UXub_1323_12_result = ME2*ME2ME3*ME3*MG*MX*1.0/s;
   return M_ddb_UXub_1323_12_result;

}

inline double M_ddb_UXub_1323_13(double MDi1, double ME2, double ME2ME3, double ME3, double P1K2, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_1323_13_result;
   M_ddb_UXub_1323_13_result = ME2ME3*ME3*1.0/s*(std::pow(MDi1, 2.0)*ME2 - 2.0*ME2*P1K2 - 4.0*ME2*P1K3 - 1.0*ME2*dM2o);
   return M_ddb_UXub_1323_13_result;

}

inline complex<double> M_ddb_UXub_1323(double CSUM, double K2K3, complex<double> LCCGuU, complex<double> LCCnuU, double MDi1, double ME1, double ME1ME3, double ME2, double ME2ME3, double ME3, double MG, double MQ, double MX, double Nc, double P1K2, double P1K3, double QK2, complex<double> RCCGdD, complex<double> RCCGuU, complex<double> RCCndD, complex<double> RCCnuU, double alpha_s, double dM2o, complex<double> rDn_k2pk3MU_pow1, complex<double> rDn_p1mk2MD_pow1, complex<double> rDn_qmk2MG_pow1, double s, double s1) {

#ifdef TEST_SC_OSSUB13
return 0.0;
#endif
   complex<double> M_ddb_UXub_1323_result;
   M_ddb_UXub_1323_result = 2.0*std::pow(M_PI, 2.0)*CSUM*std::pow(Nc, -2.0)*RCCndD*std::pow(alpha_s, 2.0)*(1.0 - 1.0*std::pow(Nc, 2.0))*(LCCGuU*rDn_k2pk3MU_pow1*rDn_p1mk2MD_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_3(MDi1, ME2, ME2ME3, ME3, MG, MQ, MX, P1K3, s)*conj(LCCnuU) + LCCGuU*rDn_k2pk3MU_pow1*rDn_p1mk2MD_pow1*M_ddb_UXub_1323_10(MDi1, ME2, ME2ME3, ME3, MG, MX, P1K3, s)*conj(LCCnuU) + LCCGuU*rDn_k2pk3MU_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_6(ME2, ME2ME3, ME3, MG, MX, P1K3, s)*conj(LCCnuU) + LCCGuU*rDn_k2pk3MU_pow1*M_ddb_UXub_1323_12(ME2, ME2ME3, ME3, MG, MX, s)*conj(LCCnuU) + LCCGuU*rDn_p1mk2MD_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_0(K2K3, MDi1, ME1, ME1ME3, ME2, ME2ME3, ME3, MG, MQ, MX, P1K3, QK2, s, s1)*conj(LCCnuU) + LCCGuU*rDn_qmk2MG_pow1*M_ddb_UXub_1323_7(ME2, ME2ME3, ME3, MG, MX, s)*conj(LCCnuU) + RCCGuU*rDn_k2pk3MU_pow1*rDn_p1mk2MD_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_4(MDi1, ME2, ME2ME3, ME3, MG, MQ, MX, P1K3, dM2o, s)*conj(RCCnuU) + RCCGuU*rDn_k2pk3MU_pow1*rDn_p1mk2MD_pow1*M_ddb_UXub_1323_11(MDi1, ME2, ME2ME3, ME3, MQ, MX, P1K3, dM2o, s)*conj(RCCnuU) + RCCGuU*rDn_k2pk3MU_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_5(MDi1, ME2, ME2ME3, ME3, MG, MQ, MX, P1K2, P1K3, dM2o, s)*conj(RCCnuU) + RCCGuU*rDn_k2pk3MU_pow1*M_ddb_UXub_1323_13(MDi1, ME2, ME2ME3, ME3, P1K2, P1K3, dM2o, s)*conj(RCCnuU) + RCCGuU*rDn_p1mk2MD_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_1323_1(K2K3, MDi1, ME1, ME1ME3, ME2, ME2ME3, ME3, MG, MQ, MX, P1K3, dM2o, s, s1)*conj(RCCnuU) + RCCGuU*rDn_p1mk2MD_pow1*M_ddb_UXub_1323_8(K2K3, MDi1, ME1, ME1ME3, ME2, ME2ME3, ME3, MG, MX, P1K3, s, s1)*conj(RCCnuU) + RCCGuU*rDn_qmk2MG_pow1*M_ddb_UXub_1323_2(K2K3, MDi1, ME1, ME1ME3, ME2, ME2ME3, ME3, MG, MQ, MX, P1K2, P1K3, s, s1)*conj(RCCnuU) + RCCGuU*M_ddb_UXub_1323_9(ME1, ME1ME3, ME2, ME2ME3, ME3, P1K3, s, s1)*conj(RCCnuU))*conj(RCCGdD);
   return M_ddb_UXub_1323_result;

}
