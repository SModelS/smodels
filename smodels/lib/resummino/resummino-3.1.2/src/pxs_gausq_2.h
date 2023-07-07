#ifndef PXS_GAUSQ_2
#define PXS_GAUSQ_2

ComplexType Z_g(POLE pIEPS, Parameters *params);

ComplexType Z_q_massless_wave(POLE pIEPS, int quark, bool right, Parameters *params);
ComplexType Z_Q_massless_wave(POLE pIEPS, int squark, Parameters *params);
ComplexType Z_Q_massless_mass(POLE pIEPS, int squark, Parameters *params);
ComplexType d_alpha(POLE pIEPS, Parameters *params);

ComplexType ME_us_QQg_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params);
ComplexType ME_ss_qQX_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params);
ComplexType ME_uu_qQX_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params);
ComplexType ME_us_qqg_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                         Parameters *params);
ComplexType ME_us_se_QQ_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params);
ComplexType ME_us_se_qq_ct(POLE pIEPS, bool sc, bool uc, bool axial, double Q2, double P1K1,
                           Parameters *params);

ComplexType ME_us_born(bool SS, bool UU, bool SS2, bool UU2, bool axialgauge, double Q2,
                       double P1K1, Parameters *params, bool qqleft = 1, bool qqright = 1,
                       bool qqgleft = 1, bool qqgright = 1, bool qQgleft = 1, bool qQgright = 1);
#endif