/******************************************************************************
 *                      Code generated with sympy 1.6.2                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                This file is part of 'SunderCalc_Resummino'                 *
 ******************************************************************************/
//#include "ddb_UXub_34.h"
#include <complex>
#include "ddb_UXub.h"

inline double M_ddb_UXub_34_0(double K2K3, double MDi2, double MDj1, double ME3, double ME4, double MG, double MQ, double MX, double P1K3) {

   double M_ddb_UXub_34_0_result;
   M_ddb_UXub_34_0_result = ME3*ME4*MG*MX*(-2.0*K2K3*std::pow(MDi2, 2.0) + 2.0*K2K3*std::pow(MG, 2.0) - 2.0*std::pow(MDi2, 2.0)*std::pow(MG, 2.0) + 2.0*std::pow(MDi2, 2.0)*std::pow(MQ, 2.0) + 2.0*std::pow(MDi2, 2.0)*P1K3 - 1.0*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + std::pow(MDj1, 2.0)*std::pow(MQ, 2.0) + 2.0*std::pow(MDj1, 2.0)*P1K3 - 2.0*std::pow(MG, 2.0)*std::pow(MQ, 2.0) + std::pow(MG, 2.0)*std::pow(MX, 2.0) - 4.0*std::pow(MG, 2.0)*P1K3 + 2.0*std::pow(MG, 4.0) - 1.0*std::pow(MQ, 2.0)*std::pow(MX, 2.0));
   return M_ddb_UXub_34_0_result;

}

inline double M_ddb_UXub_34_1(double K2K3, double MDi2, double MDj1, double ME3, double ME4, double MG, double MQ, double MX, double P1K3) {

   double M_ddb_UXub_34_1_result;
   M_ddb_UXub_34_1_result = ME3*ME4*MG*MX*(2.0*K2K3*std::pow(MDi2, 2.0) - 2.0*K2K3*std::pow(MG, 2.0) - 2.0*std::pow(MDi2, 2.0)*P1K3 - 1.0*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + std::pow(MDj1, 2.0)*std::pow(MQ, 2.0) - 2.0*std::pow(MDj1, 2.0)*P1K3 + std::pow(MG, 2.0)*std::pow(MX, 2.0) + 4.0*std::pow(MG, 2.0)*P1K3 - 1.0*std::pow(MQ, 2.0)*std::pow(MX, 2.0));
   return M_ddb_UXub_34_1_result;

}

inline double M_ddb_UXub_34_2(double MDi2, double MDj1, double ME3, double ME4, double MG, double MQ, double MX) {

   double M_ddb_UXub_34_2_result;
   M_ddb_UXub_34_2_result = -2.0*ME3*ME4*(-1.0*std::pow(MDi2, 2.0)*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + std::pow(MDi2, 2.0)*std::pow(MDj1, 2.0)*std::pow(MQ, 2.0) - 1.0*std::pow(MG, 2.0)*std::pow(MQ, 2.0)*std::pow(MX, 2.0) + std::pow(MG, 4.0)*std::pow(MX, 2.0));
   return M_ddb_UXub_34_2_result;

}

inline double M_ddb_UXub_34_3(double K2K3, double ME3, double ME4, double MG, double MQ, double MX, double P1K3, double QK2, double s) {

   double M_ddb_UXub_34_3_result;
   M_ddb_UXub_34_3_result = -2.0*ME3*ME4*MG*MX*(K2K3 - 1.0*std::pow(MQ, 2.0) + std::pow(MX, 2.0) - 1.0*P1K3 - 2.0*QK2 + s);
   return M_ddb_UXub_34_3_result;

}

inline double M_ddb_UXub_34_4(double K2K3, double ME3, double ME4, double MG, double MX, double P1K3) {

   double M_ddb_UXub_34_4_result;
   M_ddb_UXub_34_4_result = 2.0*ME3*ME4*MG*MX*(K2K3 - 1.0*P1K3);
   return M_ddb_UXub_34_4_result;

}

inline double M_ddb_UXub_34_5(double K2K3, double MDj1, double ME3, double ME4, double MG, double MQ, double MX, double QK2, double s) {

   double M_ddb_UXub_34_5_result;
   M_ddb_UXub_34_5_result = -1.0*ME3*ME4*(2.0*K2K3*std::pow(MG, 2.0) - 2.0*K2K3*std::pow(MX, 2.0) + 4.0*K2K3*QK2 - 2.0*K2K3*s - 1.0*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + 2.0*std::pow(MDj1, 2.0)*std::pow(MQ, 2.0) - 1.0*std::pow(MDj1, 2.0)*std::pow(MX, 2.0) + 2.0*std::pow(MDj1, 2.0)*QK2 - 1.0*std::pow(MDj1, 2.0)*s + std::pow(MG, 2.0)*std::pow(MX, 2.0) + 2.0*std::pow(MX, 2.0)*QK2 - 1.0*std::pow(MX, 2.0)*s - 1.0*std::pow(MX, 4.0));
   return M_ddb_UXub_34_5_result;

}

inline double M_ddb_UXub_34_6(double ME3, double ME4, double MG, double MQ, double MX, double P1K3, double QK2, double s) {

   double M_ddb_UXub_34_6_result;
   M_ddb_UXub_34_6_result = ME3*ME4*MG*MX*(std::pow(MG, 2.0) + std::pow(MQ, 2.0) - 2.0*std::pow(MX, 2.0) + 2.0*P1K3 + 4.0*QK2 - 2.0*s);
   return M_ddb_UXub_34_6_result;

}

inline double M_ddb_UXub_34_7(double ME3, double ME4, double MG, double MQ, double MX, double P1K3) {

   double M_ddb_UXub_34_7_result;
   M_ddb_UXub_34_7_result = -1.0*ME3*ME4*MG*MX*(std::pow(MG, 2.0) - 1.0*std::pow(MQ, 2.0) + 2.0*P1K3);
   return M_ddb_UXub_34_7_result;

}

inline double M_ddb_UXub_34_8(double K2K3, double MDi2, double MDj1, double ME3, double ME4, double MG, double MX, double P1K2, double QK2, double dM2o, double s) {

   double M_ddb_UXub_34_8_result;
   M_ddb_UXub_34_8_result = ME3*ME4*(-2.0*K2K3*std::pow(MG, 2.0) + 2.0*K2K3*std::pow(MX, 2.0) - 4.0*K2K3*QK2 + 2.0*K2K3*s - 2.0*std::pow(MDi2, 2.0)*QK2 + 2.0*std::pow(MDi2, 2.0)*dM2o + std::pow(MDi2, 2.0)*s - 1.0*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + std::pow(MDj1, 2.0)*std::pow(MX, 2.0) - 2.0*std::pow(MDj1, 2.0)*QK2 + std::pow(MDj1, 2.0)*s - 1.0*std::pow(MG, 2.0)*std::pow(MX, 2.0) + 4.0*std::pow(MG, 2.0)*P1K2 - 1.0*std::pow(MG, 2.0)*s + std::pow(MG, 4.0) - 4.0*std::pow(MX, 2.0)*P1K2 + std::pow(MX, 2.0)*s + 4.0*P1K2*QK2 - 2.0*P1K2*s);
   return M_ddb_UXub_34_8_result;

}

inline double M_ddb_UXub_34_9(double ME3, double ME4, double QK2, double dM2o, double s) {

   double M_ddb_UXub_34_9_result;
   M_ddb_UXub_34_9_result = ME3*ME4*(-2.0*QK2 + 2.0*dM2o + s);
   return M_ddb_UXub_34_9_result;

}

inline double M_ddb_UXub_34_10(double K2K3, double MDi2, double MDj1, double ME3, double ME4, double MG, double MQ, double MX, double P1K3) {

   double M_ddb_UXub_34_10_result;
   M_ddb_UXub_34_10_result = ME3*ME4*MG*MX*(2.0*K2K3 - 2.0*std::pow(MDi2, 2.0) - 1.0*std::pow(MDj1, 2.0) + 4.0*std::pow(MG, 2.0) - 2.0*std::pow(MQ, 2.0) + std::pow(MX, 2.0) - 4.0*P1K3);
   return M_ddb_UXub_34_10_result;

}

inline double M_ddb_UXub_34_11(double K2K3, double MDj1, double ME3, double ME4, double MG, double MX, double P1K3) {

   double M_ddb_UXub_34_11_result;
   M_ddb_UXub_34_11_result = ME3*ME4*MG*MX*(-2.0*K2K3 - 1.0*std::pow(MDj1, 2.0) + std::pow(MX, 2.0) + 4.0*P1K3);
   return M_ddb_UXub_34_11_result;

}

inline double M_ddb_UXub_34_12(double K2K3, double MDi2, double MDj1, double ME3, double ME4, double MG, double MQ, double MX, double P1K3) {

   double M_ddb_UXub_34_12_result;
   M_ddb_UXub_34_12_result = ME3*ME4*(2.0*K2K3*std::pow(MDi2, 2.0) - 2.0*K2K3*std::pow(MG, 2.0) + 2.0*std::pow(MDi2, 2.0)*std::pow(MDj1, 2.0) - 2.0*std::pow(MDi2, 2.0)*P1K3 - 1.0*std::pow(MDj1, 2.0)*std::pow(MG, 2.0) + std::pow(MDj1, 2.0)*std::pow(MQ, 2.0) + 2.0*std::pow(MDj1, 2.0)*P1K3 - 3.0*std::pow(MG, 2.0)*std::pow(MX, 2.0) + std::pow(MQ, 2.0)*std::pow(MX, 2.0));
   return M_ddb_UXub_34_12_result;

}

inline double M_ddb_UXub_34_13(double MDj1, double ME3, double ME4, double MX, double P1K3) {

   double M_ddb_UXub_34_13_result;
   M_ddb_UXub_34_13_result = -1.0*ME3*ME4*(-1.0*std::pow(MDj1, 2.0) + std::pow(MX, 2.0) + 2.0*P1K3);
   return M_ddb_UXub_34_13_result;

}

inline double M_ddb_UXub_34_14(double ME3, double ME4, double MG, double MX) {

   double M_ddb_UXub_34_14_result;
   M_ddb_UXub_34_14_result = ME3*ME4*MG*MX;
   return M_ddb_UXub_34_14_result;

}

inline double M_ddb_UXub_34_15(double ME3, double ME4, double MG, double MX) {

   double M_ddb_UXub_34_15_result;
   M_ddb_UXub_34_15_result = -1.0*ME3*ME4*MG*MX;
   return M_ddb_UXub_34_15_result;

}

inline double M_ddb_UXub_34_16(double K2K3, double MDj1, double ME3, double ME4, double MG, double P1K2, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_34_16_result;
   M_ddb_UXub_34_16_result = ME3*ME4*(-2.0*K2K3 - 1.0*std::pow(MDj1, 2.0) + std::pow(MG, 2.0) + 4.0*P1K2 + 2.0*P1K3 - 1.0*dM2o - 1.0*s);
   return M_ddb_UXub_34_16_result;

}

inline double M_ddb_UXub_34_17(double ME3, double ME4, double MG, double MX) {

   double M_ddb_UXub_34_17_result;
   M_ddb_UXub_34_17_result = 2.0*ME3*ME4*MG*MX;
   return M_ddb_UXub_34_17_result;

}

inline double M_ddb_UXub_34_18(double K2K3, double MDj1, double ME3, double ME4, double MX) {

   double M_ddb_UXub_34_18_result;
   M_ddb_UXub_34_18_result = -1.0*ME3*ME4*(2.0*K2K3 + std::pow(MDj1, 2.0) + std::pow(MX, 2.0));
   return M_ddb_UXub_34_18_result;

}

inline complex<double> M_ddb_UXub_34(double CSUM, double K2K3, complex<double> LCCGdD, complex<double> LCCGuU, complex<double> LCCndD, double MDi2, double MDj1, double ME3, double ME3ME4, double ME4, double MG, double MQ, double MX, double Nc, double P1K2, double P1K3, double QK2, complex<double> RCCGdD, complex<double> RCCGuU, complex<double> RCCndD, double alpha_s, double dM2o, complex<double> rDn_p1mk2MDj1_pow1, complex<double> rDn_qmk2MG_pow1, complex<double> rDn_qmk2MG_pow2, complex<double> rDn_qmp1mk2MDi2_pow1, double s) {

   complex<double> M_ddb_UXub_34_result;
   M_ddb_UXub_34_result = std::pow(M_PI, 2.0)*CSUM*ME3ME4*std::pow(Nc, -2.0)*std::pow(alpha_s, 2.0)*(1.0 - 1.0*std::pow(Nc, 2.0))*(
      #ifndef TEST_SC_OSSUB13
      + LCCGdD*LCCGuU*RCCGdD*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_12(K2K3, MDi2, MDj1, ME3, ME4, MG, MQ, MX, P1K3)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      + LCCGdD*LCCGuU*RCCGdD*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow1*M_ddb_UXub_34_13(MDj1, ME3, ME4, MX, P1K3)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      + LCCGdD*LCCGuU*RCCGdD*rDn_p1mk2MDj1_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_18(K2K3, MDj1, ME3, ME4, MX)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      + LCCGdD*LCCGuU*RCCGdD*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_16(K2K3, MDj1, ME3, ME4, MG, P1K2, P1K3, dM2o, s)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_10(K2K3, MDi2, MDj1, ME3, ME4, MG, MQ, MX, P1K3)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_p1mk2MDj1_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_17(ME3, ME4, MG, MX)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_14(ME3, ME4, MG, MX)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + std::pow(RCCGdD, 2.0)*RCCGuU*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_11(K2K3, MDj1, ME3, ME4, MG, MX, P1K3)*conj(RCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + std::pow(RCCGdD, 2.0)*RCCGuU*rDn_qmk2MG_pow1*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_15(ME3, ME4, MG, MX)*conj(RCCGuU)*std::pow(conj(RCCndD), 2.0) 
      #endif
      //* checked 2
      + LCCGdD*LCCGuU*RCCGdD*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_2(MDi2, MDj1, ME3, ME4, MG, MQ, MX)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_0(K2K3, MDi2, MDj1, ME3, ME4, MG, MQ, MX, P1K3)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + std::pow(RCCGdD, 2.0)*RCCGuU*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_1(K2K3, MDi2, MDj1, ME3, ME4, MG, MQ, MX, P1K3)*conj(RCCGuU)*std::pow(conj(RCCndD), 2.0) 
      //*/

      //* 4.1
      + LCCGdD*LCCGuU*RCCGdD*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*M_ddb_UXub_34_5(K2K3, MDj1, ME3, ME4, MG, MQ, MX, QK2, s)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
   /**/   + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*M_ddb_UXub_34_3(K2K3, ME3, ME4, MG, MQ, MX, P1K3, QK2, s)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + std::pow(RCCGdD, 2.0)*RCCGuU*rDn_p1mk2MDj1_pow1*rDn_qmk2MG_pow2*M_ddb_UXub_34_4(K2K3, ME3, ME4, MG, MX, P1K3)*conj(RCCGuU)*std::pow(conj(RCCndD), 2.0) 
      //*/

      //* 4.2
      + LCCGdD*LCCGuU*RCCGdD*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_8(K2K3, MDi2, MDj1, ME3, ME4, MG, MX, P1K2, QK2, dM2o, s)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
   /**/   + LCCGuU*std::pow(RCCGdD, 2.0)*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_6(ME3, ME4, MG, MQ, MX, P1K3, QK2, s)*conj(LCCGuU)*std::pow(conj(RCCndD), 2.0) 
      + std::pow(RCCGdD, 2.0)*RCCGuU*rDn_qmk2MG_pow2*rDn_qmp1mk2MDi2_pow1*M_ddb_UXub_34_7(ME3, ME4, MG, MQ, MX, P1K3)*conj(RCCGuU)*std::pow(conj(RCCndD), 2.0)
      //*/

      //* checked 1 (no contrib)
      + LCCGdD*LCCGuU*RCCGdD*rDn_qmk2MG_pow2*M_ddb_UXub_34_9(ME3, ME4, QK2, dM2o, s)*conj(LCCGuU)*conj(LCCndD)*conj(RCCndD) 
      //*/
      );
   return M_ddb_UXub_34_result;

}
