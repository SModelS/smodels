#ifndef OPTIONS_H_
#define OPTIONS_H_

#include "params.h"

// Minimum pt cut
// (useful to check e.g. against madgraph; or to check dipoles)
// (use a pt cut and then check separately the DIPOLE and REAL close to pt->0;
// you should get the same results with a different sign)
#define PTMINCUT
#define PTMIN params->ptcut

// Maximum pt cut
//#define PTMAX 99.01
//#define PTMAXCUT

#define DIPOLE // activate or deactivate the 2->3 dipole
#define REAL   // activate or deactivate the real emission

// turn off if you want to neglect the resonant TT and UU diagrams
// for associated gaugino-gluino production.
// (useful to just check the IR divergent contributions).
#define ONSUB

// remove resonance integrations (only squark gaugino for now)
//#define OS_DR

// enable contributions
#define VNLO
#define PK
#define RNLOg
#define RNLOq

// enable gausq initial states
#define SQGA_QG
#define SQGA_GG
#define SQGA_UU
#define SQGA_UUB

// incoming quark loops for gausq
#define QA_MIN 0
#define QA_MAX 5
#define QB_MIN 0
#define QB_MAX 5

// gausq virtual scheme
// default: on-shell
//#define MSBAR

// IEPS 0 (finite part), 1 (1/eps coefficient) and 2 (1/eps^2 coefficient). IEPS
// -1 UV 1/eps rest is IR here.
#define IEPS 0

// Enable dissable virtual diagrams
#define SQGA_VIRT_uu
#define SQGA_VIRT_UU
#define SQGA_VIRT_uug
#define SQGA_VIRT_UXu
#define SQGA_VIRT_UUg
#define SQGA_VIRT_gGX
#define SQGA_VIRT_BOX
// enable integrated dipoles
#define SQGA_VINTDIP

//// Finite shift from DREG to DRED.
// shift gaugino-quark-squark vertex from CDR to DRBAR
// (see arxiv:hep-ph/0511344v2 & arXiv:hep-ph/9610490).
#define DRBAR

// Decouple squarks, gluinos and top from the aS running.
// (see arxiv 9610490v1 p. 17 above eq. 19)
#define DECOUPLE_HEAVY

// unused outdated legacy
#define PHYS_GAUGE
#define DRBAR 1

#endif