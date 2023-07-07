// This File will include all the Matrix Elements for gaugino squark pair
// production

#include "kinematics.h"
#include "npf.h"
#include "utils.h"
#include <complex>
#include <iostream>
#include "constants.h"

#include "debug.h"
using namespace std;



// Born
double FI::Mss_SQGA1() {
  return real(32.0 * pow2(ivs) * papb * pbp2 * (LL + RR));
}

double FI::Muu_SQGA1() {
  return real(-8.0 * ivu2s1 * ivu2s2 * pap2 *
              (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR));
}

double FI::Msu_SQGA1() {
  return real(-8.0 * ivs * ivu2s2 *
              (-2.0 * pow2(pap2) + (m2s - p1p2) * papb +
               pap2 * (pbp1 - 2.0 * pbp2) + pap1 * (2.0 * pap2 + pbp2)) *
              (LL + RR));
}

double FI::Mss_SQGA2() {
  return real(32.0 * pow2(ivs) * pap1 * papb * (LL + RR));
}

double FI::Muu_SQGA2() {
  return real(-8.0 * ivu1s1 * ivu1s2 * pbp1 *
              (m1s + m2s - 2.0 * (p1p2 + pbp1 - pbp2)) * (LL + RR));
}

double FI::Msu_SQGA2() {
  return real(-8.0 * ivs * ivu1s2 *
              (m1s * papb - p1p2 * papb - 2.0 * pap1 * pbp1 + pap2 * pbp1 -
               2.0 * pow2(pbp1) + pap1 * pbp2 + 2.0 * pbp1 * pbp2) *
              (LL + RR));
}



vector<double> FI::Mborn_pol(double &Q2, double &PK) {
    
    vector<double>born_pol = {
        (m1s - m2s - 2*PK)/(2*Q2),
        (2*(m1s - m2s)*(m1s - m2s - 2*PK) + 2*PK*Q2 - pow2(Q2))/(2*Q2*pow2(m1s - m2s - 2*PK + Q2)),
        (-4*m1s*PK + 4*PK*(m2s + 2*PK) - 6*PK*Q2 + pow2(Q2))/(2*Q2*pow2(m1s - m2s - 2*PK + Q2)),  //flip?
        -((2*(m1s - m2s))/pow2(m1s - m2s - 2*PK + Q2)),  //flip?
        ((2*(-m1s + m2s))/pow2(Q2) - 1/Q2 + (m1s - m2s)/pow2(m1s - m2s - 2*PK + Q2) + 1/(m1s - m2s - 2*PK + Q2))/2,
        (4*(m1s - m2s - 2*PK)*(m1s - m2s - PK) + 3*(m1s - m2s - 2*PK)*Q2 + pow2(Q2))/(2*Q2*pow2(m1s - m2s - 2*PK + Q2)),  //flip?
        -((4*PK*(-m1s + m2s + 2*PK) + (-m1s + m2s - 6*PK)*Q2 + pow2(Q2))/(2*Q2*pow2(m1s - m2s - 2*PK + Q2)))
    };
    
    scale(born_pol,real(LLLL+RRRR)*(M_PI)/NC);
    
    return born_pol;
}

// Simplified version with minimal calculations
vector<vector<double>> FI::DIPCONNECT_in_in () {
    
    vector<vector<double>> DIPCONNECT = {
        {
            4.0,
            0.0,
            m1s + m2s + 2*p1p2,
            m1s,
            m1s + m2s + 2*p1p2,
            2*(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))),
            2*(m1s + p1p2)
        },
        
        {
            0.0,
            0.0,
            pow2(m1s + m2s + 2*p1p2)/4,
            pow2(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))),
            0.0,
            0.0,
            (m1s + m2s + 2*p1p2)*(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3))))
        },
        
        {
            0.0,
            pow2(pap1 + pap2 - papb),
            (pow2(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*pow2(pap3 + pbp3))/(4*pow2(papb)),
            pow2(p1p3 + ((m1s + p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/((m1s + m2s + 2*p1p2)*papb) - ((p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3))*(m1s*(pap3 - papb + pbp3) + m2s*(pap3 - papb + pbp3) + 2*(-pap2*pap3 + pap3*papb + pow2(papb) - pap2*pbp3 + papb*pbp3 - pap1*(pap3 + pbp3) + p1p2*(pap3 - papb + pbp3))))/(papb*(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3))))),
            -(((pap1 + pap2 - papb)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/papb),            
            2*(-pap1 - pap2 + papb)*(p1p3 + ((m1s + p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/((m1s + m2s + 2*p1p2)*papb) - ((p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3))*(m1s*(pap3 - papb + pbp3) + m2s*(pap3 - papb + pbp3) + 2*(-pap2*pap3 + pap3*papb + pow2(papb) - pap2*pbp3 + papb*pbp3 - pap1*(pap3 + pbp3) + p1p2*(pap3 - papb + pbp3))))/(papb*(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3))))),
            (1/papb)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3)*(p1p3 + ((m1s + p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/((m1s + m2s + 2*p1p2)*papb) - ((p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3))*(m1s*(pap3 - papb + pbp3) + m2s*(pap3 - papb + pbp3) + 2*(-pap2*pap3 + pap3*papb + pow2(papb) - pap2*pbp3 + papb*pbp3 - pap1*(pap3 + pbp3) + p1p2*(pap3 - papb + pbp3))))/(papb*(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))))
        },
        
        {
            -2*(pap1 + pap2 - papb),
            0.0,
            ((m1s + m2s + 2*p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/(2*papb),
            2*(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3))))*(p1p3 + ((m1s + p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/((m1s + m2s + 2*p1p2)*papb) - ((p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3))*(m1s*(pap3 - papb + pbp3) + m2s*(pap3 - papb + pbp3) + 2*(-pap2*pap3 + pap3*papb + pow2(papb) - pap2*pbp3 + papb*pbp3 - pap1*(pap3 + pbp3) + p1p2*(pap3 - papb + pbp3))))/(papb*(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3))))),
            -(m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb),
            2*(-pap1 - pap2 + papb)*(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))),
            (1/papb)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3)*(m1s + p1p2 + pap1 - ((m1s + m2s + 2*(p1p2 + pap1 + pap2))*(p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3)))/(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))) + (m1s + m2s + 2*p1p2)*(p1p3 + ((m1s + p1p2)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + papb))*(pap3 + pbp3))/((m1s + m2s + 2*p1p2)*papb) - ((p1p3*pap3 - pap1*pap3 - p1p3*papb + p1p3*pbp3 - pap1*pbp3 + m1s*(pap3 - 2*papb + pbp3) + p1p2*(pap3 - 2*papb + pbp3))*(m1s*(pap3 - papb + pbp3) + m2s*(pap3 - papb + pbp3) + 2*(-pap2*pap3 + pap3*papb + pow2(papb) - pap2*pbp3 + papb*pbp3 - pap1*(pap3 + pbp3) + p1p2*(pap3 - papb + pbp3))))/(papb*(m1s*(pap3 - 3*papb + pbp3) + m2s*(pap3 - 3*papb + pbp3) - 2*(pap1*pap3 + pap2*pap3 - pap3*papb + pow2(papb) + (pap1 + pap2 - papb)*pbp3 - p1p2*(pap3 - 3*papb + pbp3)))))
        }
        
    };

    
    return DIPCONNECT;
}

// Simplified version with minimal calculations
vector<vector<double>> FI::DIPCONNECT_in_out() {
    
    vector<vector<double>> DIPCONNECT = {
        {
            4.0,
            0.0,
            2*papb*(1 - p1p3/(pbp1 + pbp3)),
            m1s,
            2*papb*(1 - p1p3/(pbp1 + pbp3)),
            -2*pap2 + papb*(2 - (2*p1p3)/(pbp1 + pbp3)),
            m1s - m2s + 2*papb - (2*p1p3*papb)/(pbp1 + pbp3)                                 
        },
        
        {
            m1s,
            pow2(pap1),
            pow2(m1s + p1p2 + p1p3 - (p1p3*(m1s + p1p2 + p1p3 - pap1))/(pbp1 + pbp3)),
            pow2((m1s*(-p1p3 + pbp1 + pbp3) + p1p3*(-p1p2 - p1p3 + pap1 + pbp1 + pbp3))/(pbp1 + pbp3)),
            2*pap1*(m1s + p1p2 + p1p3 - (p1p3*(m1s + p1p2 + p1p3 - pap1))/(pbp1 + pbp3)), // check!
            (2*pap1*(m1s*(-p1p3 + pbp1 + pbp3) + p1p3*(-p1p2 - p1p3 + pap1 + pbp1 + pbp3)))/(pbp1 + pbp3),
            (2*(m1s + p1p2 + p1p3 - (p1p3*(m1s + p1p2 + p1p3 - pap1))/(pbp1 + pbp3))*(m1s*(-p1p3 + pbp1 + pbp3) + p1p3*(-p1p2 - p1p3 + pap1 + pbp1 + pbp3)))/(pbp1 + pbp3) //check!
        },
        
        {
            0.0,
            pow2(pap1 + pap2 - papb),
            pow2((m1s + m2s + 2*p1p2 - 2*papb - (p1p3*(m1s + m2s + 2*p1p2 - 2*pap1 - 2*pap2))/(pbp1 + pbp3))/2),
            pow2((p1p3*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + pbp1 + pbp3)))/(2*(pbp1 + pbp3))),
            (1/(pbp1 + pbp3))*(pap1 + pap2 - papb)*(2*p1p3*(pap1 + pap2) - 2*papb*(pbp1 + pbp3) + m1s*(-p1p3 + pbp1 + pbp3) + m2s*(-p1p3 + pbp1 + pbp3) + 2*p1p2*(-p1p3 + pbp1 + pbp3)),
            -((p1p3*(pap1 + pap2 - papb)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + pbp1 + pbp3)))/(pbp1 + pbp3)),
            -(1/(2*pow2(pbp1 + pbp3)))*p1p3*(2*p1p3*(pap1 + pap2) - 2*papb*(pbp1 + pbp3) + m1s*(-p1p3 + pbp1 + pbp3) + m2s*(-p1p3 + pbp1 + pbp3) + 2*p1p2*(-p1p3 + pbp1 + pbp3))*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + pbp1 + pbp3))
        },
        
        {
            2*p1p3,
            -2*pap1*(pap1 + pap2 - papb),
            (m1s + p1p2 + p1p3 - (p1p3*(m1s + p1p2 + p1p3 - pap1))/(pbp1 + pbp3))*(-m1s - m2s - 2*p1p2 + 2*papb + (p1p3*(m1s + m2s + 2*p1p2 - 2*pap1 - 2*pap2))/(pbp1 + pbp3)),
            (p1p3*(m1s*(-p1p3 + pbp1 + pbp3) + p1p3*(-p1p2 - p1p3 + pap1 + pbp1 + pbp3))*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + pbp1 + pbp3)))/pow2(pbp1 + pbp3),
            (1/(pbp1 + pbp3))*(m2s*pap1*(p1p3 - pbp1 - pbp3) + m1s*(3*pap1 + 2*pap2 - 2*papb)*(p1p3 - pbp1 - pbp3) + 2*(pow2(p1p3)*(pap1 + pap2 - papb) + p1p2*(2*pap1 + pap2 - papb)*(p1p3 - pbp1 - pbp3) + pap1*papb*(pbp1 + pbp3) - p1p3*(2*pow2(pap1) + (pap2 - papb)*(pbp1 + pbp3) + pap1*(2*pap2 - papb + pbp1 + pbp3)))),
            (1/(pbp1 + pbp3))*(m2s*p1p3*pap1 + m1s*p1p3*(3*pap1 + 2*pap2 - 2*papb) - 2*m1s*(pap1 + pap2 - papb)*(pbp1 + pbp3) + 2*p1p3*(-2*pow2(pap1) - 2*pap1*pap2 + p1p3*(pap1 + pap2 - papb) + p1p2*(2*pap1 + pap2 - papb) + pap1*papb - pap2*pbp1 + papb*pbp1 - pap2*pbp3 + papb*pbp3)),
            -(1/(pow2(pbp1 + pbp3)))*((m1s*(p1p3 - pbp1 - pbp3) + p1p3*(p1p2 + p1p3 - pap1 - pbp1 - pbp3))*(-2*p1p3*(pap1 + pap2) + m1s*(p1p3 - pbp1 - pbp3) + m2s*(p1p3 - pbp1 - pbp3) + 2*papb*(pbp1 + pbp3) - 2*p1p2*(-p1p3 + pbp1 + pbp3)) + p1p3*(m1s*p1p3 + p1p2*p1p3 + pow2(p1p3) - p1p3*pap1 - m1s*pbp1 - p1p2*pbp1 - p1p3*pbp1 - (m1s + p1p2 + p1p3)*pbp3)*(m1s + m2s + 2*(p1p2 - pap1 - pap2 + pbp1 + pbp3)))
        }
        
    };
    
    return DIPCONNECT;
    
}

double FI::MG_SQGA_gg_onshell_23(Parameters* params) {
    this->Cross_pa_p3();
    static const complex<double> III(0.0, 1.0);
    // Calculate all the derived kinematical quantities
    double Dn_p1mk3 = 2.*(pap2 + pap1 - papb);
    double Dn_p2mk1_MQ = 2.*(pap1 - m1s - p1p2 - p1p3);
    double Dn_qmk2_MQ = 2.*p1p3;
    complex<double> Dn_p1mk2_MQ = m2s - m1s - 2.*pap2;
    // Resonant propagator
    if(m2s-m1s < 0)Dn_p1mk2_MQ+=m1s*WIDTH*III;
    double Dn_p2mk3 = 2.*(p1p2 - pap2 - pap1) + m1s + m2s;
    complex<double> rDn_p1mk2_MQ_pow1 = 1./Dn_p1mk2_MQ;
    complex<double> rDn_p1mk2_MQ_pow2 = 1./norm(Dn_p1mk2_MQ);

    _DEBUG_TABLE("irxG_OS",WIDTH);
   
    // In the case of squark, m1 is the squark mass and m2 is the neutralino mass
    // MX neutralino
    // MQ squark
    // Adjusted the macro in utils.h, I get that always m1 is the (anti)squark mass and m2 is the neutralino mass
    
     // From ug_UXg_1.cpp (shall be gluon contribution)
    double K1Q = m1s + p1p2 + p1p3;
    double value =  -1/(2*cF)*real(MG_SQGA_uIg(
        Dn_p1mk3 ,
        Dn_p2mk1_MQ ,
        Dn_qmk2_MQ ,
        Dn_p1mk2_MQ ,
        Dn_p2mk3 ,
        0.0,
        0.0 ,
        rDn_p1mk2_MQ_pow2 ,
        K1Q 
    ));
    this->Cross_pa_p3();
    return value;

}

double FI::MG_SQGA_gg(Parameters * params) {
    this->Cross_pa_p3();
    static const complex<double> III(0.0, 1.0);
    // Calculate all the derived kinematical quantities
    double Dn_p1mk3 = 2.*(pap2 + pap1 - papb);
    double Dn_p2mk1_MQ = 2.*(pap1 - m1s - p1p2 - p1p3);
    double Dn_qmk2_MQ = 2.*p1p3;
    complex<double> Dn_p1mk2_MQ = m2s - m1s - 2.*pap2;
    // Resonant propagator
    if(m2s-m1s < 0)Dn_p1mk2_MQ+=m1s*WIDTH*III;
    double Dn_p2mk3 = 2.*(p1p2 - pap2 - pap1) + m1s + m2s;
    complex<double> rDn_p1mk2_MQ_pow1 = 1./Dn_p1mk2_MQ;
    complex<double> rDn_p1mk2_MQ_pow2 = 1./norm(Dn_p1mk2_MQ);

    _DEBUG_TABLE("irxG_OS",WIDTH);
   
    // In the case of squark, m1 is the squark mass and m2 is the neutralino mass
    // MX neutralino
    // MQ squark
    // Adjusted the macro in utils.h, I get that always m1 is the (anti)squark mass and m2 is the neutralino mass
    
     // From ug_UXg_1.cpp (shall be gluon contribution)
     double K1Q = m1s + p1p2 + p1p3;

    double value = -1/(2*cF)*real(MG_SQGA_uIg(
        Dn_p1mk3 ,
        Dn_p2mk1_MQ ,
        Dn_qmk2_MQ ,
        Dn_p1mk2_MQ ,
        Dn_p2mk3 ,
        1.0,
        rDn_p1mk2_MQ_pow1 ,
        0.0 ,
        K1Q 
    ));
    this->Cross_pa_p3();
    return value;

}

// Matrix element for real emission diagrams
double FI::MG_SQGA_qg() {
   
    // Calculate all the derived kinematical quantities
    double Dn_p1mk3 = 2.*(pap2 + pap1 - papb);
    double Dn_p2mk1_MQ = 2.*(pap1 - m1s - p1p2 - p1p3);
    double Dn_qmk2_MQ = 2.*p1p3;
    double Dn_p1mk2_MQ = m2s - m1s - 2.*pap2;
    double Dn_p2mk3 = 2.*(p1p2 - pap2 - pap1) + m1s + m2s;
    double rDn_p1mk2_MQ_pow1 = 1./Dn_p1mk2_MQ;
    double rDn_p1mk2_MQ_pow2 = 1./pow2(Dn_p1mk2_MQ);
   
    // In the case of squark, m1 is the squark mass and m2 is the neutralino mass
    // MX neutralino
    // MQ squark
    // Adjusted the macro in utils.h, I get that always m1 is the (anti)squark mass and m2 is the neutralino mass
    
     // From ug_UXg_1.cpp (shall be gluon contribution)
     double K1Q = m1s + p1p2 + p1p3;

   return real(MG_SQGA_uIg(
        Dn_p1mk3 ,
        Dn_p2mk1_MQ ,
        Dn_qmk2_MQ ,
        Dn_p1mk2_MQ ,
        Dn_p2mk3 ,
        1.0,
        rDn_p1mk2_MQ_pow1 ,
        rDn_p1mk2_MQ_pow2 ,
        K1Q 
    ));
}
 
complex<double> FI::MG_SQGA_uIg(
    double Dn_p1mk3 ,
    double Dn_p2mk1_MQ ,
    double Dn_qmk2_MQ ,
    complex<double> Dn_p1mk2_MQ ,
    double Dn_p2mk3 ,
    double rDn_p1mk2_MQ_pow0 ,
    complex<double> rDn_p1mk2_MQ_pow1 ,
    complex<double> rDn_p1mk2_MQ_pow2 ,
    double K1Q 
 ) {
 
    _DEBUG_TABLE("LLLL",real(LLLL));
    _DEBUG_TABLE("RRRR",real(RRRR));
    _DEBUG_TABLE("m1s",m1s);
    _DEBUG_TABLE("m2s",m2s);
    _DEBUG_TABLE("Dn_p1mk3",Dn_p1mk3);
    _DEBUG_TABLE("Dn_p2mk3",Dn_p2mk3);
    _DEBUG_TABLE("Dn_p2mk1_MQ",Dn_p2mk1_MQ);
    _DEBUG_TABLE("Dn_p1mk2_MQ",real(Dn_p1mk2_MQ));
    _DEBUG_TABLE("iDn_p1mk2_MQ",imag(Dn_p1mk2_MQ));
    _DEBUG_TABLE("Dn_qmk2_MQ",Dn_qmk2_MQ);
    _DEBUG_TABLE("K1Q",K1Q);
    _DEBUG_TABLE("K1P1",pap1);      
    _DEBUG_TABLE("K2P1",pap2);      
    _DEBUG_TABLE("K3P1",pap3);      
    _DEBUG_TABLE("K3P2",pbp3);

    if(m2s > 5.6e4 || m1s < 5.6e4) {
        //cout << "WRONG MASS!!!!" << endl;
    }
     
     complex<double> M_tot_1_0 = (pow(m1s,3)*pap1 + 2*pow2(m1s)*((p1p3 + 2*pap1 - 2*pap2)*(pap1 + pap2) + m2s*(3*pap1 + 4*pap2 - 4*papb) + p1p2*(2*pap1 + pap2 - papb) + 3*pap2*papb - pow2(papb)) + m1s*(-3*pow2(m2s)*pap1 - 2*pap1*(pow2(p1p2) - pow2(p1p3) - 2*(p1p2 + p1p3)*pap1 + pow2(pap1)) - 2*(pow2(p1p2) - pow2(p1p3) + 2*(p1p2 + p1p3)*pap1 + pow2(pap1))*pap2 + 4*(-2*(p1p2 + p1p3) + pap1)*pow2(pap2) + 2*(pow2(p1p2) + 6*p1p3*pap2 + pap1*(4*pap1 + pap2) + 2*p1p2*(-p1p3 + pap1 + 3*pap2))*papb - 2*(2*(p1p2 + p1p3) + 3*pap1)*pow2(papb) + 2*m2s*((3*p1p3 - 4*pap1)*(pap1 + pap2) + p1p2*(pap1 + 5*pap2 - 5*papb) + (-8*p1p3 + 5*pap1 + 2*pap2)*papb - 2*pow2(papb))) + 2*(-pow2(m2s)*p1p2*pap1 - 2*(p1p2 + p1p3 - pap1)*pap1*(p1p2*(p1p2 + p1p3) - p1p3*pap1) - 2*p1p2*(pow2(p1p2 + p1p3) + pow2(pap1))*pap2 - 2*(p1p2 + p1p3)*(p1p2 + p1p3 - pap1)*pow2(pap2) + (2*pow(p1p2,3) + p1p2*pap1*(p1p3 + 3*pap1) + p1p2*(5*p1p3 + pap1)*pap2 + pow2(p1p2)*(4*p1p3 - 3*pap1 + pap2) + (2*pow2(p1p3) - p1p3*pap1 - 2*pow2(pap1))*(pap1 + pap2))*papb + (pow2(p1p2) - 2*pow2(p1p3) + p1p3*pap1 + 2*pow2(pap1) - p1p2*(p1p3 + 3*pap1))*pow2(papb) - m2s*(pow2(p1p3)*(pap1 + pap2) - p1p3*(2*pap1 - papb)*(pap1 + pap2 - papb) + pow2(p1p2)*(pap1 - pap2 + papb) - pow2(pap1)*(pap1 - pap2 + papb) + 2*p1p2*(pap1*(pap1 + pap2) + p1p3*papb))))/(2*(m1s + m2s + 2*p1p2)*p1p3*(K1Q - pap1)*(pap1 + pap2 - papb)*papb);
     
     //cout << "M_tot_1_0 = " << M_tot_1_0 << endl;
     
         
     complex<double> M_tot_1_1 = 2/Dn_p1mk3*(-2*pap1 + m2s - m1s) + 1/Dn_p2mk1_MQ*(-10*pap1 + 2*(m2s - m1s) + 2*papb) - 1/Dn_qmk2_MQ*(8*pap1 - 3*(m2s - m1s) + 2*papb) + 32*m1s*(m2s - m1s)/(Dn_p1mk3*Dn_p2mk1_MQ*Dn_qmk2_MQ)*(p1p2 + m1s) + 32*m1s*(m2s - m1s)/(Dn_p2mk1_MQ*Dn_qmk2_MQ*2*papb)*(p1p2 + m1s) + 8*pap1/(Dn_p1mk3*Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(m1s*(m2s - m1s) + 2*m2s*pap1 - 4*pow2(pap1)) + 8*pap1/(Dn_p2mk1_MQ*2*papb*(2*p1p2 + m1s + m2s))*(m1s*(m2s - m1s) + 2*m2s*pap1 - 4*pow2(pap1)) + 4*(m2s - m1s)/(Dn_p1mk3*Dn_p2mk1_MQ)*(p1p2 + 3*m1s + pap1) + 4*(m2s - m1s)/(Dn_qmk2_MQ*2*papb)*(p1p2 + 3*m1s + pap1) + 1/papb*(-2*pap1 + m2s - m1s) - 2/(2*p1p2 + m1s + m2s)*(2*pap1 + m2s - m1s) - 8/(Dn_p1mk3*Dn_qmk2_MQ)*(m1s*pap1 - 2*m1s*(m2s - m1s) - 2*pow2(pap1)) - 4/(Dn_p1mk3*(2*p1p2 + m1s + m2s))*(m1s*(m2s - m1s) + 2*m2s*pap1 - 4*pow2(pap1)) - 8/(Dn_p2mk1_MQ*2*papb)*(m1s*pap1 - 2*m1s*(m2s - m1s) - 2*pow2(pap1)) + 1/(Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(14*m1s*pap1 - 4*m1s*(m2s - m1s) - 6*m2s*pap1 - 2*pap1*2*papb + 20*pow2(pap1) + (m2s - m1s)*2*papb) + 1/(Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(-6*m1s*m2s + 14*m1s*pap1 + 5*pow2(m1s) - 6*m2s*pap1 + pow2(m2s) + 2*pap1*2*papb + 16*pow2(pap1) - (m2s - m1s)*2*papb) - 2/(papb*(2*p1p2 + m1s + m2s))*(m1s*(m2s - m1s) + 2*m2s*pap1 - 4*pow2(pap1));
     
     //cout << "M_tot_1_1 = " << M_tot_1_1 << endl;
         
     complex<double> M_tot_1_2 = 4*(m2s - m1s)*(1/Dn_p2mk1_MQ*(-1*p1p2 + 3*m1s + pap1) + 1/Dn_qmk2_MQ*(-1*p1p2 + 3*m1s + pap1) + 8*pow2(m1s)/(Dn_p2mk1_MQ*Dn_qmk2_MQ));
     
     //cout << "M_tot_1_2 = " << M_tot_1_2 << endl;
     
     // Total contribution from gluon
     complex<double> M_tot_1 = real(-(1.0/(16*pow2(NC)))*(rDn_p1mk2_MQ_pow1*M_tot_1_1 + rDn_p1mk2_MQ_pow2*M_tot_1_2 +rDn_p1mk2_MQ_pow0* M_tot_1_0)*(LLLL + RRRR));

     
     //cout << "M_tot_1 = " << M_tot_1 << endl;
 
     

    //complex<double> M_tot_1 = real(-(1/(64*pow2(NC)))*(4/(m1s + m2s + 2*p1p2) + 4/(m1s + p1p2 + p1p3 - pap1) + (2*(m1s + m2s))/(p1p3*(m1s + p1p2 + p1p3 - pap1)) + (4*pap1)/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)) + (8*(m1s - m2s)*pow2(m1s - p1p2 + pap1))/(p1p3*(m1s + p1p2 + p1p3 - pap1)*pow2(m1s - m2s + 2*pap2)) + (2*(-m1s + m2s + 4*pap1 - (4*m2s*pap1)/(pap1 + pap2 - papb)))/((m1s + m2s + 2*p1p2)*p1p3) - 4/(pap1 + pap2 - papb) - 4/papb - (4*(m1s + p1p2)*(m1s + 2*p1p2))/(p1p3*(m1s + p1p2 + p1p3 - pap1)*papb) + (8*m2s*pap1)/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*papb) + (4*(m1s - 2*pap2))/((m1s + p1p2 + p1p3 - pap1)*papb) + (8*(m1s + pap1 - pap2))/((m1s + m2s + 2*p1p2)*papb) + (4*(m2s - p1p2 + pap2))/(p1p3*papb) - (4*(m1s - m2s)*(m1s - m2s - pap1 + pap2))/((m1s + m2s + 2*p1p2)*p1p3*papb) + (2*(m1s + p1p2)*(5*pow2(m1s) - 3*m1s*m2s + 8*m1s*p1p2 - 2*m2s*p1p2 + 4*pow2(p1p2)))/(p1p3*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)*papb) + (2*(-5*pow2(m1s) + 3*m1s*m2s - 8*m1s*p1p2 + 2*m2s*p1p2 - 4*pow2(p1p2) + 4*(m1s + p1p2)*pap2 - 4*pow2(pap2)))/(p1p3*(pap1 + pap2 - papb)*papb) + (2*(5*pow2(m1s) - 3*m1s*m2s + 8*m1s*p1p2 - 2*m2s*p1p2 + 4*pow2(p1p2) - 4*(m1s + p1p2)*pap2 + 4*pow2(pap2)))/((m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)*papb) + (2*(-5*m1s + m2s - pap1 + 3*pap2 + papb))/(p1p3*(pap1 + pap2 - papb)) + (4*(pow2(m1s) + pow2(m2s) - 2*pap1*(pap1 + pap2 - papb) + m2s*(-2*pap2 + papb) - m1s*(2*m2s - 2*pap2 + papb)))/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + (6*m1s - 6*m2s + 4*(p1p2 - pap1 - 2*pap2 + papb))/((m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + (8*m1s + 12*pap1 - 4*(pap2 + papb))/((m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb)) - (2*(8*m1s*(m1s - m2s)*(m1s + p1p2)*(m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb) - 2*(m1s - m2s)*(m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(3*m1s + p1p2 + pap1)*(pap1 + pap2 - papb) - 2*(m1s + m2s + 2*p1p2)*p1p3*(m1s + p1p2 + p1p3 - pap1)*(m1s - m2s + 2*pap1)*(pap1 + pap2 - papb) + 4*p1p3*(m1s + p1p2 + p1p3 - pap1)*(pow2(m1s) - m1s*m2s - 2*(m2s - 2*pap1)*pap1)*(pap1 + pap2 - papb) + 4*p1p3*pap1*(pow2(m1s) - m1s*m2s - 2*(m2s - 2*pap1)*pap1)*(pap1 + pap2 - papb) + 4*(m1s + m2s + 2*p1p2)*p1p3*(2*m1s*(m1s - m2s) + m1s*pap1 - 2*pow2(pap1))*(pap1 + pap2 - papb) + 8*m1s*(m1s - m2s)*(m1s + p1p2)*(m1s + m2s + 2*p1p2)*papb + 2*(m1s - m2s)*(m1s + m2s + 2*p1p2)*p1p3*(3*m1s + p1p2 + pap1)*papb - 2*(m1s + m2s + 2*p1p2)*p1p3*(m1s + p1p2 + p1p3 - pap1)*(m1s - m2s + 2*pap1)*papb + 4*p1p3*(m1s + p1p2 + p1p3 - pap1)*(pow2(m1s) - m1s*m2s - 2*(m2s - 2*pap1)*pap1)*papb + 4*(m1s + p1p2 + p1p3 - pap1)*pap1*(-pow2(m1s) + m1s*m2s + 2*(m2s - 2*pap1)*pap1)*papb + 4*(m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(2*m1s*(-m1s + m2s) - m1s*pap1 + 2*pow2(pap1))*papb + 4*p1p3*(m1s - m2s - 2*pap1)*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)*papb + 2*(m1s + m2s + 2*p1p2)*p1p3*(m1s - m2s + 5*pap1 - papb)*(pap1 + pap2 - papb)*papb - (m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)*papb*(3*m1s - 3*m2s + 8*pap1 + 2*papb) - 2*p1p3*(pap1 + pap2 - papb)*papb*(2*pow2(m1s) + 2*pap1*(5*pap1 - papb) - m1s*(2*m2s - 7*pap1 + papb) + m2s*(-3*pap1 + papb)) + (m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)*papb*(5*pow2(m1s) + pow2(m2s) - 2*m2s*(3*pap1 + papb) + 4*pap1*(4*pap1 + papb) + 2*m1s*(-3*m2s + 7*pap1 + papb))))/((m1s + m2s + 2*p1p2)*p1p3*(m1s + p1p2 + p1p3 - pap1)*(m1s - m2s + 2*pap2)*(pap1 + pap2 - papb)*papb) + (2*(pow2(m1s) - 4*pow2(p1p2) - m2s*(p1p2 + pap2 - papb) - m1s*(m2s + 3*(p1p2 - pap2 + papb))))/(p1p3*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)))*(LLLL + RRRR));
    
    
     //cout << "M_tot_1 = " << M_tot_1 << endl;
    
     // From ug_UXg_1mNc2.cpp (shall be more contribution from gluon)
     complex<double> M_tot_1mNc2_0 = 2/Dn_p2mk1_MQ + 1/Dn_qmk2_MQ - 8*m1s/(Dn_p1mk3*pow2(Dn_p2mk1_MQ))*(2*pap2 + m2s - m1s) - 8*m1s/(pow2(Dn_qmk2_MQ)*2*papb)*(2*pap2 + m2s - m1s) + 8*(m2s - m1s)/(Dn_p1mk3*pow2(2*p1p2 + m1s + m2s))*(pap1 + pap2) + 8*(m2s - m1s)/(2*papb*pow2(2*p1p2 + m1s + m2s))*(pap1 + pap2) - 8*(m2s - m1s)/pow2(2*p1p2 + m1s + m2s) + 2/papb + 2/(2.*p1p2 + m1s + m2s) - 4/(Dn_p1mk3*Dn_p2mk1_MQ)*(-3*m1s + pap2) + 4/(Dn_p1mk3*(2*p1p2 + m1s + m2s))*(m2s - pap1 - pap2 + 2*(K1Q)) + 8*pow2(1/Dn_p2mk1_MQ*m1) + 2/(Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(3*m1s + m2s + pap1) + 8*pow2(1/Dn_qmk2_MQ*m1) - 4/(Dn_qmk2_MQ*2*papb)*(-3*m1s + pap2) + 1/(Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(7*m1s + m2s + 4*pap1) + 2/(papb*(2*p1p2 + m1s + m2s))*(m1s + pap1 -pap2 - 2*(K1Q)) - 4/(Dn_p1mk3*Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(m1s*pap1 + 5*m1s*pap2 + 3*m2s*pap1 - m2s*pap2 + m2s*(m2s - m1s)) - 4/(Dn_qmk2_MQ*2*papb*(2*p1p2 + m1s + m2s))*(m1s*pap1 + 5*m1s*pap2 + 3*m2s*pap1 - m2s*pap2 + m2s*(m2s - m1s));
     
     //cout << "M_tot_1mNc2_0 = " << M_tot_1mNc2_0 << endl;
     //cout << "DDDD M_tot_1mNc2_0 = " << (M_tot_1mNc2_0+176479.23076965)/176479.23076965 << endl;
     
     complex<double> M_tot_1mNc2_1 = -4/Dn_p1mk3*(m2s - m1s) + 1/Dn_p2mk1_MQ*(16*m1s - 8*m2s - 2*pap1 + 2*papb) - 1/Dn_qmk2_MQ*(-15*m1s + 7*m2s + 2*papb) - 2*(m2s - m1s)/papb - 16*(m2s - m1s)*pow2(1/Dn_p2mk1_MQ*m1) - 16*(m2s - m1s)*pow2(1/Dn_qmk2_MQ*m1) - 2/(2*p1p2 + m1s + m2s)*(2*pap1 + m2s - m1s) + 4/(Dn_p1mk3*Dn_p2mk1_MQ)*(-7*m1s*m2s + 5*pow2(m1s) + 2*pow2(m2s)) + 16/(Dn_p1mk3*pow2(Dn_p2mk1_MQ))*pow2(m1*(m2s - m1s)) + 4/(Dn_p1mk3*(2*p1p2 + m1s + m2s))*(-3*m1s*m2s + 2*pow2(m1s) - 2*m2s*pap1 + pow2(m2s)) + 1/(Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(-2*m1s*pap1 - 6*m2s*pap1 + 4*m2s*(m2s - m1s) - 4*pap1*papb + 4*pow2(pap1) + (m2s - m1s)*2*papb) + 4/(Dn_qmk2_MQ*2*papb)*(-7*m1s*m2s + 5*pow2(m1s) + 2*pow2(m2s)) + 1/(Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(-6*m1s*m2s - 2*m1s*pap1 + pow2(m1s) - 6*m2s*pap1 + 5*pow2(m2s) + 4*pap1*papb - (m2s - m1s)*2*papb) + 16/(pow2(Dn_qmk2_MQ)*2*papb)*pow2(m1*(m2s - m1s)) + 4/(2*papb*(2*p1p2 + m1s + m2s))*(-3*m1s*m2s + 2*pow2(m1s) - 2*m2s*pap1 + pow2(m2s)) - 4/(Dn_p1mk3*Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(-4*m1s*pow2(m2s) + 5*pow2(m1s)*m2s - 2*pow(m1s, 3) - 2*m2s*pap1*(m2s - m1s) + pow(m2s, 3)) - 4/(Dn_qmk2_MQ*2*papb*(2*p1p2 + m1s + m2s))*(-4*m1s*pow2(m2s) + 5*pow2(m1s)*m2s - 2*pow(m1s, 3) - 2*m2s*pap1*(m2s - m1s) + pow(m2s, 3)) + 8;
     
     //cout << "M_tot_1mNc2_1 = " << M_tot_1mNc2_1 << endl;
     
     complex<double> M_tot_1mNc2_2 = -8*(m2s - m1s)*(2/Dn_p2mk1_MQ*m1s + 2/Dn_qmk2_MQ*m1s + 2*pow2(m1s)/pow2(Dn_p2mk1_MQ) + 2*pow2(m1s)/pow2(Dn_qmk2_MQ) + 1);
     
     //cout << "M_tot_1mNc2_2 = " << M_tot_1mNc2_2 << endl;
     complex<double> M_tot_1mNc2 = real(-(1.0/(16*pow2(NC)))*(pow2(NC) - 1)*(rDn_p1mk2_MQ_pow1*M_tot_1mNc2_1 + rDn_p1mk2_MQ_pow2*M_tot_1mNc2_2 + rDn_p1mk2_MQ_pow0* M_tot_1mNc2_0)*(LLLL + RRRR));
     
     //cout << "M_tot_1mNc2 = " << M_tot_1mNc2 << endl;


    // Total contribution from gluon (Nc^2 terms)
//    double M_tot_1mNc2 = real(-(1/(16*pow2(NC)))*(-1 + pow2(NC))*((8*(m1s - m2s))/pow2(m1s + m2s + 2*p1p2) + 2/(m1s + m2s + 2*p1p2) + (2*m1s)/pow2(p1p3) + 1/(2*p1p3) + (2*m1s)/pow2(m1s + p1p2 + p1p3 - pap1) - 1/(m1s + p1p2 + p1p3 - pap1) - (3*m1s + m2s + pap1)/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)) + (7*m1s + m2s + 4*pap1)/(2*m1s*p1p3 + 2*m2s*p1p3 + 4*p1p2*p1p3) - (4*(-m1s + m2s)*(2 + m1s*(2/p1p3 + m1s*(1/pow2(p1p3) + 1/pow2(m1s + p1p2 + p1p3 - pap1)) - 2/(m1s + p1p2 + p1p3 - pap1))))/pow2(m1s - m2s + 2*pap2) + (m1s*(m1s - m2s - 2*pap2))/(pow2(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + (2*(m2s + 2*(m1s + p1p2 + p1p3) - pap1 - pap2))/((m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb)) + (-3*m1s + pap2)/((m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + (4*(-m1s + m2s)*(pap1 + pap2))/(pow2(m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb)) + (m2s*(m2s + 3*pap1 - pap2) + m1s*(-m2s + pap1 + 5*pap2))/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + 2/papb + (m1s*(m1s - m2s - 2*pap2))/(pow2(p1p3)*papb) + (3*m1s - pap2)/(p1p3*papb) - (2*(m1s + 2*p1p2 + 2*p1p3 - pap1 + pap2))/((m1s + m2s + 2*p1p2)*papb) + (4*(-m1s + m2s)*(pap1 + pap2))/(pow2(m1s + m2s + 2*p1p2)*papb) - (m2s*(pow2(m2s) + 3*pap1 - pap2) + m1s*(-pow2(m2s) + pap1 + 5*pap2))/((m1s + m2s + 2*p1p2)*p1p3*papb) + (1/(-m1s + m2s - 2*pap2))*(8 + (4*m1s*(m1s - m2s))/pow2(p1p3) + (2*(m1s - m2s - 2*pap1))/(m1s + m2s + 2*p1p2) + (4*m1s*(m1s - m2s))/pow2(m1s + p1p2 + p1p3 - pap1) + (15*m1s - 7*m2s - 2*papb)/(2*p1p3) + (-8*m1s + 4*m2s + pap1 - papb)/(m1s + p1p2 + p1p3 - pap1) + (2*(m1s - m2s))/(pap1 + pap2 - papb) + (2*(2*pow2(m1s) - 3*m1s*m2s + m2s*(m2s - 2*pap1)))/((m1s + m2s + 2*p1p2)*(pap1 + pap2 - papb)) + (2*m1s*pow2(m1s - m2s))/(pow2(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) - ((5*m1s - 2*m2s)*(m1s - m2s))/((m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) - ((m1s - m2s)*(2*pow2(m1s) - 3*m1s*m2s + m2s*(m2s - 2*pap1)))/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)*(pap1 + pap2 - papb)) + (2*(m1s - m2s))/papb + (2*m1s*pow2(m1s - m2s))/(pow2(p1p3)*papb) + ((5*m1s - 2*m2s)*(m1s - m2s))/(p1p3*papb) + (2*(2*pow2(m1s) - 3*m1s*m2s + m2s*(m2s - 2*pap1)))/((m1s + m2s + 2*p1p2)*papb) + ((m1s - m2s)*(2*pow2(m1s) - 3*m1s*m2s + m2s*(m2s - 2*pap1)))/((m1s + m2s + 2*p1p2)*p1p3*papb) + (-2*pow2(m2s) + m2s*(3*pap1 - papb) + 2*pap1*(-pap1 + papb) + m1s*(2*m2s + pap1 + papb))/((m1s + m2s + 2*p1p2)*(m1s + p1p2 + p1p3 - pap1)) + (pow2(m1s) + 5*pow2(m2s) - 2*m1s*(3*m2s + pap1 - papb) + 4*pap1*papb - 2*m2s*(3*pap1 + papb))/(2*(m1s + m2s + 2*p1p2)*p1p3)))*(LLLL + RRRR));
    
//     cout << "M_tot_1mNc2 = " << M_tot_1mNc2 << endl;
    
    
    // From ug_UXg_Nc2.cpp (shall be more contribution from gluon)
    complex<double> M_tot_Nc2_0 = 15/Dn_p1mk3 - 5/Dn_p2mk1_MQ - 8/Dn_p2mk3 + 6/Dn_qmk2_MQ + 40*(m2s - m1s)*papb/(Dn_p2mk3*pow2(2*p1p2 + m1s + m2s)) - 80*(m2s - m1s)*papb/(pow2(Dn_p2mk3)*(2*p1p2 + m1s + m2s)) + 80*(m2s - m1s)*pow2(papb)/(pow2(Dn_p2mk3)*pow2(2*p1p2 + m1s + m2s)) + 8*(m2s - m1s)/(Dn_p1mk3*pow2(2*p1p2 + m1s + m2s))*(-2*pap1 - 2*pap2 + 2*papb) - 4*(m2s - m1s)/(Dn_qmk2_MQ*2*papb*(2*p1p2 + m1s + m2s))*(pap1 - pap2 + m2s - m1s) + 16*(m2s - m1s)/(pow2(2*p1p2 + m1s + m2s)) + 1/papb + 21/(2*p1p2 + m1s + m2s) + 2/(Dn_p1mk3*Dn_p2mk1_MQ)*(7*p1p2 + 2*m1s + 2*m2s + pap1 - 4*(m1s + p1p2 + p1p3)) - 2/(Dn_p1mk3*Dn_p2mk3)*(10*p1p2 + 6*m1s + 4*m2s - 22*pap1 - 25*pap2 + 6*(m1s + p1p2 + p1p3)) - 2/(Dn_p1mk3*(2*p1p2 + m1s + m2s))*(7*m1s - 3*m2s + 8*pap1 + 3*pap2 - 8*(m1s + p1p2 + p1p3) + 2*papb) + 2/(Dn_p2mk1_MQ*Dn_p2mk3)*(6*p1p2 + 5*m1s + m2s - pap1 - 10*pap2 - 5*(m1s + p1p2 + p1p3)) - 1/(Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(2*m1s + 2*m2s - 6*pap2 - 6*papb) - 2/(Dn_p2mk3*Dn_qmk2_MQ)*(2*p1p2 + m1s + m2s - 2*pap1 - 2*pap2) - 1/(Dn_p2mk3*2*papb)*(2*p1p2 + m1s + m2s - 2*pap1 - 2*pap2) - 2/(Dn_p2mk3*(2*p1p2 + m1s + m2s))*(-20*pap2 - 20*(m1s + p1p2 + p1p3) - 14*(m2s - m1s) + 10*papb) + 1/(Dn_qmk2_MQ*2*papb)*(4*p1p2 + 3*m1s + 7*m2s + 4*pap1) + 1/(Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(-9*m1s + m2s - 4*pap1 + 8*pap2 - 8*papb) + 2/(papb*(2*p1p2 + m1s + m2s))*(pap1 - pap2 - m1s - p1p2 - p1p3 - m2s + m1s) + 1/(Dn_p1mk3*Dn_p2mk1_MQ*Dn_p2mk3)*(-12*p1p2*m2s + 4*p1p2*pap1 + 4*p1p2*pap2 + 16*p1p2*(m1s + p1p2 + p1p3) - 12*pow2(p1p2) - 8*m1s*m2s + 26*m1s*pap1 - 2*m1s*pap2 - 20*m1s*(m1s + p1p2 + p1p3) + 7*pow2(m1s) - 34*m2s*pap1 - 6*m2s*pap2 + 36*m2s*(m1s + p1p2 + p1p3) + pow2(m2s) - 28*pap1*pap2 - 24*pap1*(m1s + p1p2 + p1p3) + 8*pow2(pap1) + 12*pap2*(m1s + p1p2 + p1p3) + 8*pow2(pap2) + 8*pow2(m1s + p1p2 + p1p3)) - 1/(Dn_p1mk3*Dn_p2mk1_MQ*(2*p1p2 + m1s + m2s))*(4*m1s*pap2 - 8*m2s*pap1 - 12*m2s*pap2 + 8*m2s*papb + 20*pap1*(m1s + p1p2 + p1p3) - 20*pap1*papb + 8*pap2*(m1s + p1p2 + p1p3) - 4*pow2(pap2) + 4*(m1s + p1p2 + p1p3)*(m2s - m1s) + 4*(m1s + p1p2 + p1p3)*papb - 8*pow2(m1s + p1p2 + p1p3) + 4*pow2(m2s - m1s) + 4*pow2(papb)) + 4/(Dn_p2mk1_MQ*Dn_p2mk3*(2*p1p2 + m1s + m2s))*(4*pap2*papb + 4*pow2(pap2) + 2*(m2s - m1s)*papb + pow2(m2s - m1s)) + 1/(Dn_p2mk3*Dn_qmk2_MQ*2*papb)*(4*p1p2*m1s - 24*p1p2*m2s + 4*p1p2*pap1 + 4*p1p2*pap2 - 12*pow2(p1p2) - 12*m1s*m2s + 9*pow2(m1s) - 5*pow2(m2s) + 6*pap1*(m2s - m1s) + 8*pow2(pap1) + 6*pap2*(m2s - m1s) + 8*pow2(pap2)) + 4/(Dn_p2mk3*Dn_qmk2_MQ*(2*p1p2 + m1s + m2s))*(-4*pap2*papb + 4*pow2(pap2) - 2*(m2s - m1s)*papb + pow2(m2s - m1s)) - (46*p1p2 + 43*m1s + 3*m2s - 46*pap1 - 46*pap2)/pow2(Dn_p2mk3);
    
//     cout << "M_tot_Nc2_0 = " << M_tot_Nc2_0 << endl;
   
    complex<double> M_tot_Nc2_1 = -16/Dn_p1mk3/Dn_p2mk1_MQ/Dn_p2mk3*m1s*(2*p1p2*(m2s - m1s) - pow2(m1s) + pow2(m2s)) - 1/Dn_p1mk3/Dn_p2mk1_MQ*(16*p1p2*pap1 + 4*p1p2*(m1s + p1p2 + p1p3) + 8*p1p2*(m2s - m1s) - 20*p1p2*papb - 10*m1s*m2s + 22*m1s*pap1 - 20*m1s*(m1s + p1p2 + p1p3) + 2*m1s*papb + 9*pow2(m1s) - 14*m2s*pap1 + 16*m2s*(m1s + p1p2 + p1p3) - 6*m2s*papb + pow2(m2s) - 16*pap1*(m1s + p1p2 + p1p3) + 16*pap1*papb - 12*(m1s + p1p2 + p1p3)*papb + 12*pow2(m1s + p1p2 + p1p3)) + 8/Dn_p1mk3/Dn_p2mk3*(11*p1p2*m1s - 3*p1p2*m2s - 6*p1p2*(m1s + p1p2 + p1p3) + 4*pow2(p1p2) - 3*m1s*m2s - 7*m1s*(m1s + p1p2 + p1p3) + 7*pow2(m1s) + m2s*(m1s + p1p2 + p1p3) + 2*pow2(m1s + p1p2 + p1p3)) + 2/Dn_p1mk3/(2*p1p2 + m1s + m2s)*(20*m1s*m2s - 18*m1s*papb - 12*pow2(m1s) + 8*m2s*pap1 + 10*m2s*papb - 8*pow2(m2s) + 12*pap1*(m1s + p1p2 + p1p3) - 12*pap1*papb - 14*(m1s + p1p2 + p1p3)*(m2s - m1s) + 8*(m1s + p1p2 + p1p3)*papb - 8*pow2(m1s + p1p2 + p1p3)) + 1/Dn_p1mk3*(-18*p1p2 - 21*m1s + 7*m2s - 12*pap1 + 18*(m1s + p1p2 + p1p3) + 2*papb) + 16/Dn_p2mk1_MQ/Dn_p2mk3*m1s*(m2s - m1s)*2*papb/(2*p1p2 + m1s + m2s) - 1/Dn_p2mk1_MQ/(2*p1p2 + m1s + m2s)*(-8*m1s*(m2s - m1s) - 22*m1s*papb + 6*m2s*papb + 20*pap1*(m1s + p1p2 + p1p3) - 20*pap1*papb - 10*(m1s + p1p2 + p1p3)*(m2s - m1s) + 8*(m1s + p1p2 + p1p3)*papb - 8*pow2(m1s + p1p2 + p1p3)) + 1/Dn_p2mk1_MQ*(10*p1p2 + 2*m2s + 4*(m1s + p1p2 + p1p3) - 6*papb) + 40/pow2(Dn_p2mk3)*(m2s - m1s)*2*papb/(2*p1p2 + m1s + m2s)*(2*(m1s + p1p2 + p1p3) + m2s - m1s) - 40/pow2(Dn_p2mk3)*(m2s - m1s)*(-2*p1p2 - 2*m1s + 2*(m1s + p1p2 + p1p3) + 2*papb) - 16/Dn_p2mk3/Dn_qmk2_MQ*m1s*(m2s - m1s)*2*papb/(2*p1p2 + m1s + m2s) + 16/Dn_p2mk3/Dn_qmk2_MQ*m1s*(m2s - m1s) - 16/Dn_p2mk3/Dn_qmk2_MQ*m1s/(2*papb)*(2*p1p2*(m2s - m1s) - pow2(m1s) + pow2(m2s)) - 8/Dn_p2mk3/(2*papb)*(2*p1p2*(m1s + p1p2 + p1p3) + p1p2*(m2s - m1s) + m1s*(m1s + p1p2 + p1p3) + m1s*(m2s - m1s) + m2s*(m1s + p1p2 + p1p3) - 2*pow2(m1s + p1p2 + p1p3)) + 8/Dn_p2mk3/(2*p1p2 + m1s + m2s)*(-5*(m1s + p1p2 + p1p3)*(m2s - m1s) + 10*(m1s + p1p2 + p1p3)*papb - 10*pow2(m1s + p1p2 + p1p3) + 10*(m2s - m1s)*papb - 8*pow2(m2s - m1s)) - 4/Dn_p2mk3*(5*p1p2 + 12*m1s - 7*m2s - 15*(m1s + p1p2 + p1p3) + 10*papb) + 2/Dn_qmk2_MQ/(2*papb)*(p1p2*(m2s - m1s) + 3*m1s*m2s - 4*pow2(m1s) + pow2(m2s) + pap1*(m2s - m1s)) + 1/Dn_qmk2_MQ/(2*p1p2 + m1s + m2s)*(6*m1s*pap1 - 22*m1s*papb + 10*m2s*pap1 + 6*m2s*papb - 12*pap1*papb - 3*pow2(m2s - m1s)) + 1/Dn_qmk2_MQ*(-4*(m2s - m1s) + 6*papb) + 2*(m1s + p1p2 + p1p3)/papb/(2*p1p2 + m1s + m2s)*(-2*pap1 - 4*(m1s + p1p2 + p1p3) + m2s - m1s) + 1/(2*papb)*(2*pap1 + 8*(m1s + p1p2 + p1p3) + 3*(m2s - m1s)) + 1/(2*p1p2 + m1s + m2s)*(11*m1s - 35*m2s + 10*pap1 - 8*(m1s + p1p2 + p1p3) + 16*papb) + 9;
    
//     cout << "M_tot_Nc2_1 = " << M_tot_Nc2_1 << endl;
    
   complex<double>  M_tot_Nc2_2 = -(m2s - m1s)*(2/Dn_p2mk1_MQ*(p1p2 - 3*m1s - 5.0*pap1 + 4.0*(m1s + p1p2 + p1p3)) + 8/Dn_p2mk3*(5*p1p2 + m1s - 5*(m1s + p1p2 + p1p3)) - 2/Dn_qmk2_MQ*(-p1p2 + 3*m1s + pap1) - 80*(2*p1p2*m1s - 2*p1p2*(m1s + p1p2 + p1p3) + pow2(p1p2) - 2*m1s*(m1s + p1p2 + p1p3) + pow2(m1s) + pow2(m1s + p1p2 + p1p3))/(pow2(Dn_p2mk3)) + 3.0);
    
//     cout << "M_tot_Nc2_2 = " << M_tot_Nc2_2 << endl;
    
    // Total contribution from gluon (other Nc^2 terms)
   complex<double>  M_tot_Nc2 = real((1.0/16)*(rDn_p1mk2_MQ_pow1*M_tot_Nc2_1 + rDn_p1mk2_MQ_pow2*M_tot_Nc2_2 + rDn_p1mk2_MQ_pow0* M_tot_Nc2_0)*(LLLL + RRRR));
    
//     cout << "M_tot_Nc2 = " << M_tot_Nc2 << endl;
    
    // From ug_UXg_Nc2m2.cpp (shall be more contribution from gluon)
    complex<double> M_tot_Nc2m2_0 = 1/Dn_p1mk3 + 1/Dn_p2mk1_MQ + 1/Dn_qmk2_MQ - 2*(m1s + p1p2 + p1p3)/(papb*(2*p1p2 + m1s + m2s)) + 3/(2*papb) - 2/(Dn_p1mk3*Dn_p2mk1_MQ)*(-p1p2 - 2*m1s + m2s) + 2/(Dn_p1mk3*(2*p1p2 + m1s + m2s))*(-2*pap1 + 2*(m1s + p1p2 + p1p3) + m2s - m1s) - 2/(Dn_qmk2_MQ*2*papb)*(-p1p2 - 2*m1s + m2s);
    
//     cout << "M_tot_Nc2m2_0 = " << M_tot_Nc2m2_0 << endl;
    
    complex<double> M_tot_Nc2m2_1 = 2/Dn_p1mk3/Dn_p2mk1_MQ*(-p1p2*(m2s - m1s) - 3*m1s*m2s + 2*pow2(m1s) + pow2(m2s) - pap1*(m2s - m1s)) - 2/Dn_p1mk3/(2*p1p2 + m1s + m2s)*(4*pap1*(m1s + p1p2 + p1p3) - 4*pow2(pap1) + 2*(m1s + p1p2 + p1p3)*(m2s - m1s) + pow2(m2s - m1s)) - 1/Dn_p1mk3*(-2*pap1 + m2s - m1s) - 1/Dn_p2mk1_MQ*(-7*m1s + 3*m2s) + 1/Dn_qmk2_MQ/papb*(-p1p2*(m2s - m1s) - 3*m1s*m2s + 2*pow2(m1s) + pow2(m2s) - pap1*(m2s - m1s)) - 1/Dn_qmk2_MQ*(-7*m1s + 3*m2s) + 2*(m1s + p1p2 + p1p3)/papb/(2*p1p2 + m1s + m2s)*(2*pap1 + m2s - m1s) - 1/(2*papb)*(2*pap1 + 3*(m2s - m1s)) + 4/(2*p1p2 + m1s + m2s)*(-2*pap1 + m2s - m1s) + 4;
    
//     cout << "M_tot_Nc2m2_1 = " << M_tot_Nc2m2_1 << endl;
    
    complex<double> M_tot_Nc2m2_2 = -2*(m2s - m1s)*(1/Dn_p2mk1_MQ + 1/Dn_qmk2_MQ)*(-p1p2 + 3*m1s + pap1);
    
//     cout << "M_tot_Nc2m2_2 = " << M_tot_Nc2m2_2 << endl;
    
    // Total contribution from gluon (other Nc^2 terms)
    complex<double> M_tot_Nc2m2 = real((pow2(NC) - 2)/(16*pow2(NC))*(rDn_p1mk2_MQ_pow1*M_tot_Nc2m2_1 + rDn_p1mk2_MQ_pow2*M_tot_Nc2m2_2 + rDn_p1mk2_MQ_pow0* M_tot_Nc2m2_0)*(LLLL + RRRR));
    
//     cout << "M_tot_Nc2m2 = " << M_tot_Nc2m2 << endl;
    
    // From ug_UXg_ghost.cpp (ghost diagram contribution)
    complex<double> M_tot_ghost_0 = 4/Dn_p2mk3 + 8*(m2s - m1s)*papb/(Dn_p2mk3*pow2(2*p1p2 + m1s + m2s)) - 16*(m2s - m1s)*papb/(pow2(Dn_p2mk3)*(2*p1p2 + m1s + m2s)) + 16*(m2s - m1s)*pow2(papb)/(pow2(Dn_p2mk3)*pow2(2*p1p2 + m1s + m2s)) - 2/(Dn_p2mk3*(2*p1p2 + m1s + m2s))*(-4*(m1s + p1p2 + p1p3) + 2*papb) - 2*(6*p1p2 + 5*m1s + m2s - 6*pap1 - 6*pap2)/(pow2(Dn_p2mk3));
    
//     cout << "M_tot_ghost_0 = " << M_tot_ghost_0 << endl;
    
    complex<double> M_tot_ghost_1 = -4/Dn_p2mk3*(p1p2 + m2s - 3*(m1s + p1p2 + p1p3) + 2*papb) + 16*(m2s - m1s)*papb/(pow2(Dn_p2mk3)*(2*p1p2 + m1s + m2s))*(2*(m1s + p1p2 + p1p3) + m2s - m1s) - 8*(m2s - m1s)*(-2*p1p2 - 2*m1s + 2*(m1s + p1p2 + p1p3) + 2*papb)/(pow2(Dn_p2mk3)) + 2/(2*p1p2 + m1s + m2s)*(-2*(m1s + p1p2 + p1p3) + 2*papb) + 8/(Dn_p2mk3*(2*p1p2 + m1s + m2s))*(-(m1s + p1p2 + p1p3)*(m2s - m1s) + 2*(m1s + p1p2 + p1p3)*papb - 2*pow2(m1s + p1p2 + p1p3) + 2*(m2s - m1s)*papb) + 1.0;
    
//     cout << "M_tot_ghost_1 = " << M_tot_ghost_1 << endl;
    
    complex<double> M_tot_ghost_2 = (m2s - m1s)*(8/Dn_p2mk3*(-p1p2 + p1p2 + p1p3) + 16*(2*p1p2*m1s - 2*p1p2*(m1s + p1p2 + p1p3) + pow2(p1p2) - 2*m1s*(m1s + p1p2 + p1p3) + pow2(m1s) + pow2(m1s + p1p2 + p1p3))/(pow2(Dn_p2mk3)) + 1);
    
//     cout << "M_tot_ghost_2 = " << M_tot_ghost_2 << endl;
    
    // Total contribution from ghost
   complex<double> M_tot_ghost = real((-1./32)*(rDn_p1mk2_MQ_pow1*M_tot_ghost_1 + rDn_p1mk2_MQ_pow2*M_tot_ghost_2 + rDn_p1mk2_MQ_pow0* M_tot_ghost_0)*(LLLL + RRRR));
    
//     cout << "M_tot_ghost = " << M_tot_ghost << endl;
    
    // From ug_UXg_antighost.cpp (antighost diagram contribution)
    complex<double> M_tot_antighost_0 = 4/Dn_p2mk3 + 8*(m2s - m1s)*papb/(Dn_p2mk3*pow2(2*p1p2 + m1s + m2s)) - 16*(m2s - m1s)*papb/(pow2(Dn_p2mk3)*(2*p1p2 + m1s + m2s)) + 16*(m2s - m1s)*pow2(papb)/(pow2(Dn_p2mk3)*pow2(2*p1p2 + m1s + m2s)) - 2/(Dn_p2mk3*(2*p1p2 + m1s + m2s))*(-4*(m1s + p1p2 + p1p3) + 2*papb) - 2*(6*p1p2 + 5*m1s + m2s - 6*pap1 - 6*pap2)/(pow2(Dn_p2mk3));
    
//     cout << "M_tot_antighost_0 = " << M_tot_antighost_0 << endl;
    
    complex<double> M_tot_antighost_1 = -4/Dn_p2mk3*(p1p2 + m2s - 3*(m1s + p1p2 + p1p3) + 2*papb) + 16*(m2s - m1s)*papb/(pow2(Dn_p2mk3)*(2*p1p2 + m1s + m2s))*(2*(m1s + p1p2 + p1p3) + m2s - m1s) - 8*(m2s - m1s)*(-2*p1p2 - 2*m1s + 2*(m1s + p1p2 + p1p3) + 2*papb)/(pow2(Dn_p2mk3)) + 2/(2*p1p2 + m1s + m2s)*(-2*(m1s + p1p2 + p1p3) + 2*papb) + 8/(Dn_p2mk3*(2*p1p2 + m1s + m2s))*(-(m1s + p1p2 + p1p3)*(m2s - m1s) + 2*(m1s + p1p2 + p1p3)*papb - 2*pow2(m1s + p1p2 + p1p3) + 2*(m2s - m1s)*papb) + 1;
    
//     cout << "M_tot_antighost_1 = " << M_tot_antighost_1 << endl;
    
    complex<double> M_tot_antighost_2 = (m2s - m1s)*(-8/Dn_p2mk3*(p1p2 + m1s - (m1s + p1p2 + p1p3)) + 16*(2*p1p2*m1s - 2*p1p2*(m1s + p1p2 + p1p3) + pow2(p1p2) - 2*m1s*(m1s + p1p2 + p1p3) + pow2(m1s) + pow2(m1s + p1p2 + p1p3))/(pow2(Dn_p2mk3)) + 1);

//     cout << "M_tot_antighost_2 = " << M_tot_antighost_2 << endl;
    
    // Total contribution from antighost
   complex<double>  M_tot_antighost = real((-1./32)*(rDn_p1mk2_MQ_pow1*M_tot_antighost_1 + rDn_p1mk2_MQ_pow2*M_tot_antighost_2 + rDn_p1mk2_MQ_pow0* M_tot_antighost_0)*(LLLL + RRRR));
    
//     cout << "M_tot_antighost = " << M_tot_antighost << endl;


    // Total contribution from gluon emission diagrams
    complex<double> M_tot = M_tot_1 + M_tot_1mNc2 + M_tot_Nc2 + M_tot_Nc2m2 + M_tot_ghost + M_tot_antighost;
    
//     cout << "M_tot = " << M_tot << endl;
    
    //return real(M_tot_1/((1.0/(16*pow2(NC)))*(LLLL + RRRR)));
    return M_tot;
}


// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// first vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = ml2s = ml3s = 0
double FI::Mss_qqg1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 = -16.0 * papb * (c22 * pap2 * papb - 5.0 * c00 * pbp2 -
                                        2.0 * (-c1 + c12 + c22) * papb * pbp2 +
                                        c2 * papb * (pap2 + pbp2)) *
                        (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      16.0 * papb *
      (c22 * pap2 * papb - 2.0 * (4.0 * c00 + (-c1 + c12 + c22) * papb) * pbp2 +
       c2 * papb * (pap2 + pbp2)) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 48.0 * c00 * papb * pbp2 * (LL + RR);

  return real(pow2(cA) * cF / (32.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// second vertex correction (quark-gluon-gluon loop)
// --> Later (have to sum over squarks?)
// m1ls = msquark, ml2s=m3ls=mgluino
double FI::Mss_qqg2_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -8.0 * papb *
      (LL * LsRs * (6.0 * c00 + c0 * ml2s + 2.0 * c12 * papb) * pbp2 +
       c0 * LsRs * m2 * ml2 * papb * RL +
       (c0 * LR * m2 * ml2 * papb +
        (6.0 * c00 + c0 * ml2s + 2.0 * c12 * papb) * pbp2 * RR) *
           RsLs -
       2.0 * c22 * papb * (pap2 - pbp2) * (LL * LsRs + RR * RsLs) -
       c2 * papb * (2.0 * LL * LsRs * (pap2 + pbp2) + LsRs * m2 * ml2 * RL +
                    (LR * m2 * ml2 + 2.0 * (pap2 + pbp2) * RR) * RsLs));
  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8.0 * papb *
      (8.0 * c00 * pbp2 * (LL * LsRs + RR * RsLs) +
       2.0 * (c12 - c2 + c22) * papb * pbp2 * (LL * LsRs + RR * RsLs) +
       c0 * (LL * LsRs * ml2s * pbp2 + LsRs * m2 * ml2 * papb * RL +
             LR * m2 * ml2 * papb * RsLs + ml2s * pbp2 * RR * RsLs));

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 16.0 * c00 * pbp2 * papb * RR * RsLs -
                        16.0 * c00 * pbp2 * papb * LL * LsRs;

  return real(pow2(cA) * cF / (32.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// third vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = ml2s = ml3s = 0
double FI::Mss_qqg3_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 = -32.0 * papb *
                        ((c2 + c22) * pap2 * papb - 3.0 * c00 * pbp2 +
                         (c0 + c1 - c12 + 2.0 * c2 - c22) * papb * pbp2) *
                        (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      32.0 * papb *
      ((c2 + c22) * pap2 * papb - 7.0 * c00 * pbp2 +
       (c0 + c1 - 2.0 * c12 + 3.0 * c2 - 2.0 * c22) * papb * pbp2) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 =
      32.0 * papb * (5.0 * c00 + (c12 - c2 + c22) * papb) * pbp2 * (LL + RR);

  return real(-cF / (32.0 * pow2(M_PI)) * pow2(ivs) * (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qqg vertex (vertex on the left) correction
// fourth vertex correction (quark-gluon-gluon loop)
// arguments: p1s = pa^2; p2s = pb^2; p3^2 = (pa + pb)^2
// ml1s = mgluino, ml2s=ml3s=msquark
// Later sum over squarks
double FI::Mss_qqg4_SQGA(double p1s, double p2s, double p3s, double ml1s,
                         double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      4.0 * papb *
      (4.0 * c00 * LL * LsRs * pbp2 + c0 * LsRs * m1l * m2 * papb * RL +
       c0 * LR * m1l * m2 * papb * RsLs + 4.0 * c00 * pbp2 * RR * RsLs -
       4.0 * c22 * pap2 * papb * (LL * LsRs + RR * RsLs) +
       2.0 * c2 * papb * (-(LL * LsRs * pap2) + LsRs * m1l * m2 * RL +
                          LR * m1l * m2 * RsLs - pap2 * RR * RsLs));
  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 = -16.0 * c00 * papb * pbp2 * (LL * LsRs + RR * RsLs);
  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(cF / (32.0 * pow2(M_PI)) * pow2(ivs) * (me0 + me1 + me2));
}

// squared matrix element (SS) w/ qsqga vertex (vertex on the right) correction
// first vertex correction (quark-gluon-squark-loop)
// arguments: p1s = k^2; p2s = (-p1)^2; p3^2 = (-p2)^2
// ml1s = mgluino, ml2s=ml3s=msquark
// Later sum over squarks
double FI::Mss_qsqga1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -8.0 * papb * ((2.0 * c1 + c2) * (p1p2 * papb - pap2 * pbp1) +
                     (4.0 * c00 + 2.0 * c2 * m1s + c22 * m1s - 2.0 * c1 * pap1 -
                      2.0 * c12 * pap1 - 3.0 * c2 * pap1 - 2.0 * c22 * pap1 +
                      2.0 * c1 * papb + 2.0 * c11 * papb + 4.0 * c12 * papb +
                      2.0 * c2 * papb + 2.0 * c22 * papb -
                      2.0 * (2.0 * c1 + c12 + 2.0 * c2 + c22) * pbp1) *
                         pbp2) *
      (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8 * papb * ((2.0 * c1 + c2) * (p1p2 * papb - pap2 * pbp1) +
                  (6.0 * c00 + 2.0 * c2 * m1s + c22 * m1s - 2.0 * c1 * pap1 -
                   2.0 * c12 * pap1 - 3.0 * c2 * pap1 - 2.0 * c22 * pap1 +
                   2.0 * c1 * papb + 2.0 * c11 * papb + 4.0 * c12 * papb +
                   2.0 * c2 * papb + 2.0 * c22 * papb -
                   2.0 * (2.0 * c1 + c12 + 2.0 * c2 + c22) * pbp1) *
                      pbp2) *
      (LL + RR);
  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = -16.0 * c00 * papb * pbp2 * (LL + RR);

  return real(-cA * pow2(cF) / (16.0 * pow2(M_PI)) * pow2(ivs) *
              (me0 + me1 + me2));
}

// how to define couplings like L*Lsp or L*L or Rsp*Rsp??

// squared matrix element (SS) w/ qsqga vertex (vertex on the right) correction
// second vertex correction (quark-gluon-squark-loop)
// arguments: p1s = k^2; p2s = (-p1)^2; p3^2 = (-p2)^2
// ml1s = msquark, ml2s=mgluino, ml3s=msquark
// Later sum over squarks and quarks(?)
/*
double FI::Mss_qsqga2_SQGA(double p1s, double p2s, double p3s,
                         double ml1s, double ml2s, double ml3s, int ieps) {

  complex<double> c0,c00,c1,c2,c11,c12,c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me0 = 8.0*papb*(-Lsp*Lsp*m2*ml2*
                                  ((c0 + c1 + c2)*papb - (c0 + c2)*pbp1)* RR +
                                2.0*LLsp*((c0 + c1 + c2)*(p1p2*papb - pap2*pbp1)
+
                                         (-4.0*c00 - c2*m1s - c22*m1s + c0*pap1
+
                                                     c1*pap1 + 2.0*c12*pap1 +
3.0*c2*pap1 +
                                                     2.0*c22*pap1 - 2.0*c0*papb
- 4.0*c1*papb -
                                                     2.0*c11*papb - 4.0*c12*papb
- 4.0*c2*papb -
                                          2.0*c22*papb + 2.0*(c12 + c2 +
c22)*pbp1)* pbp2)*RRsp -
                                  LL*m2*ml2*((c0 + c1 + c2)*papb - (c0 +
c2)*pbp1)*RspRsp;

    // I1
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps+1, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me1 = 8*papb*(Power(Lsp,2)*m2*ml2*((c0 + c1 + c2)*papb - (c0 +
c2)*pbp1)*Power(R,2) - 2*L*Lsp*
                                ((c0 + c1 + c2)*(p1p2*papb - pap2*pbp1) +
(-6*c00 - c2*m1s - c22*m1s + c0*pap1 +
                                                                             c1*pap1
+ 2*c12*pap1 + 3*c2*pap1 +
                                                                             2*c22*pap1
- 2*c0*papb - 4*c1*papb -
                                                                             2*c11*papb
- 4*c12*papb - 4*c2*papb -
                                                                  2*c22*papb +
2*(c12 + c2 + c22)*pbp1)*
                                                                 pbp2)*R*Rsp +
                                                                Power(L,2)*m2*ml2*
                                                                ((c0 + c1 +
c2)*papb - (c0 + c2)*pbp1)*
                                Power(Rsp,2));
   // I2
  SetC(p1s, p2s,  p3s,  ml1s,  ml2s,  ml3s, ieps+2, &c0, &c00, &c1, &c2, &c11,
&c12, &c22);

  complex<double> me2 = - 32.0*c00*pbp2*papb*L*R*Rsp*Lsp;

  return real( cA*pow2(cF)  / (16.0 * pow2(M_PI) ) *  pow2(ivs) * ( me0 + me1 +
me2));

}*/

// UU-channel
// squared matrix element (uu) w/ sqaurk-sqaurk-gaugino vertex (top)
// first correction

double FI::Muu_qsqga1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      2.0 * (-((2.0 * c1 + c2) * m2s * pap1) +
             (4.0 * c00 + c11 * m1s + 2.0 * c12 * m1s + c2 * m1s + c22 * m1s +
              2.0 * c2 * m2s + c22 * m2s + 2.0 * c12 * p1p2 + 4.0 * c2 * p1p2 +
              2.0 * c22 * p1p2 + c1 * (m1s + 4.0 * p1p2) -
              2.0 * (c12 + c2 + c22) * pap1) *
                 pap2 -
             2.0 * (2.0 * c2 + c22) * pow2(pap2)) *
      (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      -4.0 * c00 * pap2 * (m1s + m2s - 2.0 * (p1p2 - pap1 + pap2)) * (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(-cA * pow2(cF) / (16.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}

// squared matrix element (uu) w/ sqaurk-sqaurk-gaugino vertex (top)
// second correction (squark-gluino-quark-loop)
// yet to come -> Couplings

// squared matrix element (uu) w/ gluon-squark-squark vertex (bottom) correction
// first correction
double FI::Muu_gsqsq1_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22);

  complex<double> me0 =
      -2.0 * pap2 *
      (-(c1 * p1p2 * papb) + 2.0 * c11 * p1p2 * papb + 3.0 * c2 * p1p2 * papb -
       2.0 * c22 * p1p2 * papb + c1 * pap1 * papb - 2.0 * c11 * pap1 * papb -
       3.0 * c2 * pap1 * papb + 2.0 * c22 * pap1 * papb +
       2.0 * c12 * pow(papb, 2) + 2.0 * c2 * pow(papb, 2) -
       2.0 * c22 * pow(papb, 2) - c1 * m2s * pbp1 + 2.0 * c11 * m2s * pbp1 +
       3.0 * c2 * m2s * pbp1 - 2.0 * c22 * m2s * pbp1 + 2.0 * c1 * pap2 * pbp1 -
       4.0 * c11 * pap2 * pbp1 - 6.0 * c2 * pap2 * pbp1 +
       4.0 * c22 * pap2 * pbp1 + 2.0 * c12 * papb * pbp1 +
       2.0 * c2 * papb * pbp1 - 2.0 * c22 * papb * pbp1 +
       6.0 * c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb +
                    pbp1 - pbp2) +
       (c1 * p1p2 - 2.0 * c11 * p1p2 - 3.0 * c2 * p1p2 + 2.0 * c22 * p1p2 -
        c1 * pap1 + 2.0 * c11 * pap1 + 3.0 * c2 * pap1 - 2.0 * c22 * pap1 -
        4.0 * c12 * papb - 4.0 * c2 * papb + 4.0 * c22 * papb -
        2.0 * (c12 + c2 - c22) * pbp1) *
           pbp2 +
       2.0 * (c12 + c2 - c22) * pow(pbp2, 2) +
       c0 * (pap1 * papb - m2s * pbp1 + 2.0 * pap2 * pbp1 - pap1 * pbp2 +
             p1p2 * (-papb + pbp2))) *
      (LL + RR);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 1, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me1 =
      8.0 * c00 * pap2 *
      (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb + pbp1 - pbp2) *
      (LL + RR);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps + 2, &c0, &c00, &c1, &c2, &c11,
       &c12, &c22);

  complex<double> me2 = 0;

  return real(cF * pow2(cA) / (32.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}

// squared matrix element (uu) w/ gluon-squark-squark vertex (bottom) correction
// second correction
double FI::Muu_gsqsq2_SQGA(double p1s, double p2s, double p3s, double ml1s,
                           double ml2s, double ml3s, int ieps) {

  complex<double> c0, c00, c1, c2, c11, c12, c22, c001, c002, c003, c111, c112,
      c122, c222;

  // ieps = 0 (default)
  // -> finite term of Laurent Series I0
  // ( 1/eps^2 I2 + 1/eps I1 + I0)
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);
  complex<double> c121 = c112;

  complex<double> me0 =
      -2.0 * pap2 *
      (6.0 * c002 * m2s + c1 * pow(m2s, 2) + 2.0 * c11 * pow(m2s, 2) +
       c111 * pow(m2s, 2) + 2.0 * c112 * pow(m2s, 2) + 4.0 * c12 * pow(m2s, 2) +
       c121 * pow(m2s, 2) + 3.0 * c122 * pow(m2s, 2) + c2 * pow(m2s, 2) +
       2.0 * c22 * pow(m2s, 2) + c222 * pow(m2s, 2) + c1 * m2s * pow(ml2, 2) +
       c2 * m2s * pow(ml2, 2) - 6.0 * c002 * p1p2 - c1 * m2s * p1p2 -
       2.0 * c11 * m2s * p1p2 - c111 * m2s * p1p2 - 2.0 * c112 * m2s * p1p2 -
       4.0 * c12 * m2s * p1p2 - c121 * m2s * p1p2 - 3.0 * c122 * m2s * p1p2 -
       c2 * m2s * p1p2 - 2.0 * c22 * m2s * p1p2 - c222 * m2s * p1p2 -
       c1 * pow(ml2, 2) * p1p2 - c2 * pow(ml2, 2) * p1p2 + 6.0 * c002 * pap1 +
       c1 * m2s * pap1 + 2.0 * c11 * m2s * pap1 + c111 * m2s * pap1 +
       2.0 * c112 * m2s * pap1 + 4.0 * c12 * m2s * pap1 + c121 * m2s * pap1 +
       3.0 * c122 * m2s * pap1 + c2 * m2s * pap1 + 2.0 * c22 * m2s * pap1 +
       c222 * m2s * pap1 + c1 * pow(ml2, 2) * pap1 + c2 * pow(ml2, 2) * pap1 +
       6.0 * c001 * (m2s - p1p2 + pap1 - 2.0 * pap2) - 12.0 * c002 * pap2 -
       4.0 * c1 * m2s * pap2 - 8.0 * c11 * m2s * pap2 -
       4.0 * c111 * m2s * pap2 - 8.0 * c112 * m2s * pap2 -
       16.0 * c12 * m2s * pap2 - 4.0 * c121 * m2s * pap2 -
       12.0 * c122 * m2s * pap2 - 4.0 * c2 * m2s * pap2 -
       8.0 * c22 * m2s * pap2 - 4.0 * c222 * m2s * pap2 -
       2.0 * c1 * pow(ml2, 2) * pap2 - 2.0 * c2 * pow(ml2, 2) * pap2 +
       2.0 * c1 * p1p2 * pap2 + 4.0 * c11 * p1p2 * pap2 +
       2.0 * c111 * p1p2 * pap2 + 4.0 * c112 * p1p2 * pap2 +
       8.0 * c12 * p1p2 * pap2 + 2.0 * c121 * p1p2 * pap2 +
       6.0 * c122 * p1p2 * pap2 + 2.0 * c2 * p1p2 * pap2 +
       4.0 * c22 * p1p2 * pap2 + 2.0 * c222 * p1p2 * pap2 -
       2.0 * c1 * pap1 * pap2 - 4.0 * c11 * pap1 * pap2 -
       2.0 * c111 * pap1 * pap2 - 4.0 * c112 * pap1 * pap2 -
       8.0 * c12 * pap1 * pap2 - 2.0 * c121 * pap1 * pap2 -
       6.0 * c122 * pap1 * pap2 - 2.0 * c2 * pap1 * pap2 -
       4.0 * c22 * pap1 * pap2 - 2.0 * c222 * pap1 * pap2 +
       4.0 * c1 * pow(pap2, 2) + 8.0 * c11 * pow(pap2, 2) +
       4.0 * c111 * pow(pap2, 2) + 8.0 * c112 * pow(pap2, 2) +
       16.0 * c12 * pow(pap2, 2) + 4.0 * c121 * pow(pap2, 2) +
       12.0 * c122 * pow(pap2, 2) + 4.0 * c2 * pow(pap2, 2) +
       8.0 * c22 * pow(pap2, 2) + 4.0 * c222 * pow(pap2, 2) +
       6.0 * c002 * papb + c1 * m2s * papb + c11 * m2s * papb +
       2.0 * c112 * m2s * papb + 6.0 * c12 * m2s * papb + c121 * m2s * papb +
       6.0 * c122 * m2s * papb + 2.0 * c2 * m2s * papb +
       5.0 * c22 * m2s * papb + 3.0 * c222 * m2s * papb +
       c2 * pow(ml2, 2) * papb - c112 * p1p2 * papb - 4.0 * c12 * p1p2 * papb -
       c121 * p1p2 * papb - 4.0 * c122 * p1p2 * papb - 2.0 * c2 * p1p2 * papb -
       4.0 * c22 * p1p2 * papb - 2.0 * c222 * p1p2 * papb + c112 * pap1 * papb +
       4.0 * c12 * pap1 * papb + c121 * pap1 * papb + 4.0 * c122 * pap1 * papb +
       2.0 * c2 * pap1 * papb + 4.0 * c22 * pap1 * papb +
       2.0 * c222 * pap1 * papb - 2.0 * c1 * pap2 * papb -
       2.0 * c11 * pap2 * papb - 4.0 * c112 * pap2 * papb -
       12.0 * c12 * pap2 * papb - 2.0 * c121 * pap2 * papb -
       12.0 * c122 * pap2 * papb - 4.0 * c2 * pap2 * papb -
       10.0 * c22 * pap2 * papb - 6.0 * c222 * pap2 * papb +
       2.0 * c12 * pow(papb, 2) + 2.0 * c122 * pow(papb, 2) +
       2.0 * c22 * pow(papb, 2) + 2.0 * c222 * pow(papb, 2) +
       6.0 * c002 * pbp1 + c1 * m2s * pbp1 + c11 * m2s * pbp1 +
       c112 * m2s * pbp1 + 2.0 * c12 * m2s * pbp1 + 2.0 * c122 * m2s * pbp1 +
       c22 * m2s * pbp1 + c222 * m2s * pbp1 + c2 * pow(ml2, 2) * pbp1 -
       2.0 * c1 * pap2 * pbp1 - 2.0 * c11 * pap2 * pbp1 -
       2.0 * c112 * pap2 * pbp1 - 4.0 * c12 * pap2 * pbp1 -
       4.0 * c122 * pap2 * pbp1 - 2.0 * c22 * pap2 * pbp1 -
       2.0 * c222 * pap2 * pbp1 + 2.0 * c12 * papb * pbp1 +
       2.0 * c122 * papb * pbp1 + 2.0 * c22 * papb * pbp1 +
       2.0 * c222 * papb * pbp1 +
       4.0 * c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb +
                    pbp1 - pbp2) -
       (6.0 * c002 + c1 * m2s + c11 * m2s + 2.0 * c112 * m2s + 6.0 * c12 * m2s +
        c121 * m2s + 6.0 * c122 * m2s + 2.0 * c2 * m2s + 5.0 * c22 * m2s +
        3.0 * c222 * m2s + c2 * pow(ml2, 2) - c112 * p1p2 - 4.0 * c12 * p1p2 -
        c121 * p1p2 - 4.0 * c122 * p1p2 - 2.0 * c2 * p1p2 - 4.0 * c22 * p1p2 -
        2.0 * c222 * p1p2 + c112 * pap1 + 4.0 * c12 * pap1 + c121 * pap1 +
        4.0 * c122 * pap1 + 2.0 * c2 * pap1 + 4.0 * c22 * pap1 +
        2.0 * c222 * pap1 - 2.0 * c1 * pap2 - 2.0 * c11 * pap2 -
        4.0 * c112 * pap2 - 12.0 * c12 * pap2 - 2.0 * c121 * pap2 -
        12.0 * c122 * pap2 - 4.0 * c2 * pap2 - 10.0 * c22 * pap2 -
        6.0 * c222 * pap2 + 4.0 * c12 * papb + 4.0 * c122 * papb +
        4.0 * c22 * papb + 4.0 * c222 * papb +
        2.0 * (c12 + c122 + c22 + c222) * pbp1) *
           pbp2 +
       2.0 * (c12 + c122 + c22 + c222) * pow(pbp2, 2)) *
      (LL + RR) * (LsRs + RsLs);

  // I1
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);
  complex<double> me1 =
      4.0 * pap2 *
      (c001 * (m2s - p1p2 + pap1 - 2.0 * pap2) +
       c00 * (2.0 * m2s - 2.0 * p1p2 + 2.0 * pap1 - 4.0 * pap2 + papb + pbp1 -
              pbp2) +
       c002 * (m2s - p1p2 + pap1 - 2.0 * pap2 + papb + pbp1 - pbp2)) *
      (LL + RR) * (LsRs + RsLs);

  // I2
  SetC(p1s, p2s, p3s, ml1s, ml2s, ml3s, ieps, &c0, &c00, &c1, &c2, &c11, &c12,
       &c22, &c001, &c002, &c111, &c112, &c122, &c222);

  complex<double> me2 = 0;

  return real(cF * pow2(cA) / (32.0 * pow2(M_PI)) * ivu1s1 * ivu1s2 *
              (me0 + me1 + me2));
}
