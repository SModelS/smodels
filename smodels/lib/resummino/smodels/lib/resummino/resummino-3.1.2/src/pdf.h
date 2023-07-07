// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Returns the PDFs for gluons and quarks and the corresponding $\alpha_S$.
// Uses LHAPDF to retrieve the PDFs in the usual space, and for conjugate
// spaces, it performs a numerical fit.

#ifndef PDF_H_
#define PDF_H_

#include <cmath>
#include <complex>

using namespace std;

#include "LHAPDF/LHAPDF.h"

// The changes from LHAPDF5 to LHAPDF requires us to define this as global variable.                                                                        
#if LHAPDF_MAJOR_VERSION == 6
extern LHAPDF::PDF *pdf;
#endif
const int pdf_members = 13;
extern bool enable_pdfs[pdf_members];

// Returns the value of $\alpha_S$ at energy $Q=\sqrt{Qsq}$.
#if LHAPDF_MAJOR_VERSION == 6
static inline double aS(const double Qsq, const int set) {
  //return 0.0991;
  return pdf->alphasQ2(Qsq);
}
#else
static inline double aS(const double Qsq, const int set) {
  return LHAPDF::alphasPDF(sqrt(Qsq));
}
#endif

// Returns in g and q the values of the PDF of gluons and quarks at the value
// of parameters x and Q^2 = Q2.
void pdfX(double &g, double q[2][6], double x, double Q2);

// Returns the PDFs in N-space for gluons and quarks in g and q.
void pdfN(complex<double> &g, complex<double> q[2][6], const complex<double> nm,
          double A[8][8]);

// Evolves PDF in N-space.
void pdfEvolve(complex<double> &gg, complex<double> qq[2][6],
               complex<double> nm, complex<double> lmbd);

void pdfFit(double &A1MIN, double A[8][8], double tau, double Q2,
            double weight_valence, double weight_sea, double weight_gluon,
            double xmin);

#endif
