/******************************************************************************
 *                      Code generated with sympy 1.6.2                       *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                This file is part of 'SunderCalc_Resummino'                 *
 ******************************************************************************/
//#include "ddb_UXub_111222.h"
#include <complex>
#include "../../debug.h"
#include "ddb_UXub.h"

//double M_ddb_UXub_111222_0(double K2K3, double ME1, double ME1ME1, double ME1ME2, double ME2, double ME2ME2, double P1K2, double P1K3, double dM2o, double s, double s1) {
inline double M_ddb_UXub_111222_0(double K2K3, double ME1, double ME1ME1, double ME1ME2, double ME2, double ME2ME2, double P1K2, double P1K3, double dM2o, double s, double s1) {

   double M_ddb_UXub_111222_0_result;
   M_ddb_UXub_111222_0_result = -4.0*std::pow(ME1, 2.0)*ME1ME1*P1K3*dM2o*1.0/(s*std::pow(s1, 2.0)) + 8.0*ME1ME1*dM2o*std::pow(ME1*P1K3, 2.0)*1.0/(std::pow(s, 2.0)*std::pow(s1, 2.0)) + ME1ME1*dM2o*std::pow(ME1*1.0/s1, 2.0) - 2.0*P1K3*1.0/(std::pow(s, 2.0)*s1)*(-2.0*K2K3*ME1*ME1ME2*ME2 + 2.0*K2K3*std::pow(ME1, 2.0)*ME1ME1 + 4.0*ME1*ME1ME2*ME2*P1K2 + 4.0*ME1*ME1ME2*ME2*P1K3 + ME1*ME1ME2*ME2*dM2o - 4.0*std::pow(ME1, 2.0)*ME1ME1*P1K2 - 4.0*std::pow(ME1, 2.0)*ME1ME1*P1K3 - 1.0*std::pow(ME1, 2.0)*ME1ME1*dM2o) - 1.0*1.0/s*(-2.0*ME1*ME1ME2*ME2 + std::pow(ME1, 2.0)*ME1ME1 + 2.0*std::pow(ME2, 2.0)*ME2ME2) - 1.0*1.0/s1*(ME1*ME1ME2*ME2 - 1.0*std::pow(ME1, 2.0)*ME1ME1) - 1.0*1.0/(s*s1)*(2.0*K2K3*ME1*ME1ME2*ME2 - 6.0*ME1*ME1ME2*ME2*P1K2 - 6.0*ME1*ME1ME2*ME2*P1K3 - 2.0*ME1*ME1ME2*ME2*dM2o + 2.0*std::pow(ME1, 2.0)*ME1ME1*P1K2 + 4.0*std::pow(ME1, 2.0)*ME1ME1*P1K3 + std::pow(ME1, 2.0)*ME1ME1*dM2o) + 2.0*(-3.0*ME1*ME1ME2*ME2*P1K2 - 3.0*ME1*ME1ME2*ME2*P1K3 + std::pow(ME1, 2.0)*ME1ME1*P1K2 + std::pow(ME1, 2.0)*ME1ME1*P1K3 + 2.0*std::pow(ME2, 2.0)*ME2ME2*P1K2 + 2.0*std::pow(ME2, 2.0)*ME2ME2*P1K3)*1.0/(std::pow(s, 2.0));
   //M_ddb_UXub_111222_0_result = _DEBUG_TABLE("M_ddb_UXub_111222_0",real((dM2o/(s1*s1)-4.0*P1K3*dM2o/(s*s1*s1)+(dM2o-2.0*K2K3+2.0*P1K3+4.0*P1K2)/(s*s1)+8.0*P1K3*P1K3*dM2o/(s*s*s1*s1)-1.0/s));
   _DEBUG_TABLE("M_ddb_UXub_111222_0",real(M_ddb_UXub_111222_0_result)*4);
   return M_ddb_UXub_111222_0_result;
}

inline double M_ddb_UXub_111222_1(double ME1, double ME1ME2, double ME2, double ME2ME2, double MQ, double MX, double P1K2, double P1K3, double dM2o, double s, double s1) {

   double M_ddb_UXub_111222_1_result;
   M_ddb_UXub_111222_1_result = 16.0*ME1ME2*ME2*P1K3*dM2o*1.0/(std::pow(s, 2.0)*s1)*(ME1*P1K2 + ME1*P1K3) - 1.0*ME1ME2*ME2*1.0/s1*(ME1*std::pow(MQ, 2.0) + ME1*std::pow(MX, 2.0) - 2.0*ME1*P1K2) - 2.0*ME1ME2*ME2*1.0/(s*s1)*(4.0*ME1*P1K2*P1K3 + 2.0*ME1*P1K2*dM2o + 4.0*ME1*std::pow(P1K2, 2.0) + 4.0*ME1*P1K3*dM2o + ME1*std::pow(dM2o, 2.0)) + 2.0*1.0/s*(ME1*ME1ME2*ME2*std::pow(MQ, 2.0) - 1.0*ME1*ME1ME2*ME2*P1K2 - 1.0*ME1*ME1ME2*ME2*P1K3 - 2.0*std::pow(ME2, 2.0)*ME2ME2*std::pow(MQ, 2.0) + std::pow(ME2, 2.0)*ME2ME2*std::pow(MX, 2.0) + 2.0*std::pow(ME2, 2.0)*ME2ME2*P1K2 + 2.0*std::pow(ME2, 2.0)*ME2ME2*P1K3) - 4.0*(-4.0*ME1*ME1ME2*ME2*P1K2*P1K3 - 1.0*ME1*ME1ME2*ME2*P1K2*dM2o - 2.0*ME1*ME1ME2*ME2*std::pow(P1K2, 2.0) - 1.0*ME1*ME1ME2*ME2*P1K3*dM2o - 2.0*ME1*ME1ME2*ME2*std::pow(P1K3, 2.0) + 4.0*std::pow(ME2, 2.0)*ME2ME2*P1K2*P1K3 + std::pow(ME2, 2.0)*ME2ME2*P1K2*dM2o + 2.0*std::pow(ME2, 2.0)*ME2ME2*std::pow(P1K2, 2.0) + std::pow(ME2, 2.0)*ME2ME2*P1K3*dM2o + 2.0*std::pow(ME2, 2.0)*ME2ME2*std::pow(P1K3, 2.0))*1.0/(std::pow(s, 2.0));
   _DEBUG_TABLE("M_ddb_UXub_111222_1",real(M_ddb_UXub_111222_1_result)*4);
   return M_ddb_UXub_111222_1_result;

}

inline double M_ddb_UXub_111222_2(double ME2, double ME2ME2, double MQ, double P1K2, double P1K3, double dM2o, double s) {

   double M_ddb_UXub_111222_2_result;
   M_ddb_UXub_111222_2_result = 2.0*ME2ME2*dM2o*(1.0/s*(std::pow(ME2, 2.0)*std::pow(MQ, 2.0) - 2.0*std::pow(ME2, 2.0)*P1K2 - 2.0*std::pow(ME2, 2.0)*P1K3) + 4.0*(2.0*std::pow(ME2, 2.0)*P1K2*P1K3 + std::pow(ME2, 2.0)*std::pow(P1K2, 2.0) + std::pow(ME2, 2.0)*std::pow(P1K3, 2.0))*1.0/(std::pow(s, 2.0)));
   _DEBUG_TABLE("M_ddb_UXub_111222_2",real(M_ddb_UXub_111222_2_result)*4);
   return M_ddb_UXub_111222_2_result;

}

inline complex<double> M_ddb_UXub_111222(double CSUM, double K2K3, complex<double> LCCnuU, double ME1, double ME1ME1, double ME1ME2, double ME2, double ME2ME2, double MQ, double MX, double Nc, double P1K2, double P1K3, double alpha_s, double dM2o, complex<double> rDn_k2pk3MU_pow1, complex<double> rDn_k2pk3MU_pow2, double s, double s1) {

#ifdef TEST_SC_OSSUB13
return 0.0;
#endif
   complex<double> M_ddb_UXub_111222_result;
   M_ddb_UXub_111222_result = -4.0*std::pow(M_PI, 2.0)*CSUM*std::pow(Nc, -2.0)*std::pow(alpha_s, 2.0)*(1.0 - 1.0*std::pow(Nc, 2.0))*(rDn_k2pk3MU_pow1*M_ddb_UXub_111222_1(ME1, ME1ME2, ME2, ME2ME2, MQ, MX, P1K2, P1K3, dM2o, s, s1) + rDn_k2pk3MU_pow2*M_ddb_UXub_111222_2(ME2, ME2ME2, MQ, P1K2, P1K3, dM2o, s) + M_ddb_UXub_111222_0(K2K3, ME1, ME1ME1, ME1ME2, ME2, ME2ME2, P1K2, P1K3, dM2o, s, s1))*norm(LCCnuU);
   //_DEBUG_MSG("%e",(M_ddb_UXub_111222_0(K2K3, ME1, ME1ME1, ME1ME2, ME2, ME2ME2, P1K2, P1K3, dM2o, s, s1)));
   //_DEBUG_MSG("%e",(norm(LCCnuU)));
   return M_ddb_UXub_111222_result;

}
