//#ifndef CONSTANTS_H_
//#define CONSTANTS_H_

/** On shell resonance width */
#define WIDTH params->width
// Color factors.
#define cA 3.0
#define cF (4.0 / 3.0)
#define TR 0.5
/** Number of colors */
#define NC 3.0
#define NA 8.0

#define I2R 0.5
#define MPIsI 1.0 / pow2(M_PI)

#define NFLVR 5  // Number of active flavours [ 5 RECOMMENDED ]
#define nf NFLVR // Number of active flavours [ 5 RECOMMENDED ]

// Flavor constants given in hep-ph/0201036 p.53 eq. C.12 and p.32 eq. 5.91 and
// p.37 eq. 6.17
#define I2R (0.5) // SU(3) normalization factor
#define GAMMA_Q (3.0 / 2.0 * cF)
#define GAMMA_GL (3.0 / 2.0 * cA)
#define GAMMA_SQ (2.0 * cF)
#define GAMMA_G (11.0 / 6.0 * cA - 2.0 / 3.0 * I2R * nf)
#define K_Q (7.0 / 2.0 - pow2(M_PI) / 6.0) * cF
#define K_GL (7.0 / 2.0 - pow2(M_PI) / 6.0) * cA
#define K_SQ (4.0 - pow2(M_PI) / 6.0) * cF
#define K_G ((67.0 / 18.0 - pow2(M_PI) / 6.0) * cA - 10.0 / 9.0 * I2R * nf)
