// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Thanks to Michele Papucci for the GSL 2 patch.
//
// Implements the PDF module and the PDF fit (needed for Mellin space PDFs).

#include "pdf.h"
#include "LHAPDF/LHAPDF.h"
#include "gsl_all.h"
#include "maths.h"
#include "utils.h"
#include <cmath>
#include <complex>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#if LHAPDF_MAJOR_VERSION == 6
LHAPDF::PDF *pdf;
#endif

bool enable_pdfs[13];

#define NFLVR 5    // Number of active flavours.
#define NDATA 1000 // Number of data points for the fit

// function used to sample the PDFs.
double sampling(double xmin, double xmax, const size_t n, int i) {
  return pow(xmin, 1.0 - (1.0 + (double)i) / (xmax + (double)n));
  // return xmin + (double)i / (double)n * (xmax - xmin) - xmin;
  // double step_size = (xmax - xmin)/(double)n;
  // return(xmin + (double)i * step_size);
}

void pdfX(double &g, double q[2][6], double x, double Q2) {
  g = 0.0;
  for (int i0 = 0; i0 < 6; i0++) {
    q[0][i0] = 0.0;
    q[1][i0] = 0.0;
  }

  // Sets the PDF values to the output vector.
#if LHAPDF_MAJOR_VERSION == 6
  vector<double> cpdf;
  pdf->xfxQ(x, sqrt(Q2), cpdf);
#else
  vector<double> cpdf = LHAPDF::xfx(x, sqrt(Q2));
#endif
  g = enable_pdfs[6]*cpdf[6] / x;        // gluon
  q[0][0] = enable_pdfs[7]*cpdf[7] / x;  // d
  q[1][0] = enable_pdfs[5]*cpdf[5] / x;  // dbar
  q[0][1] = enable_pdfs[9]*cpdf[9] / x;  // s
  q[1][1] = enable_pdfs[3]*cpdf[3] / x;  // sbar
  q[0][2] = enable_pdfs[11]*cpdf[11] / x; // b
  q[1][2] = enable_pdfs[1]*cpdf[1] / x;  // bbar
  q[0][3] = enable_pdfs[8]*cpdf[8] / x;  // u
  q[1][3] = enable_pdfs[4]*cpdf[4] / x;  // ubar
  q[0][4] = enable_pdfs[10]*cpdf[10] / x; // c
  q[1][4] = enable_pdfs[2]*cpdf[2] / x;  // cbar
  // q[0][5] = enable_pdfs[12]*cpdf[12] / x; // t
  // q[1][5] = enable_pdfs[0]*cpdf[0] / x;  // tbar
}

void pdfN(complex<double> &g, complex<double> q[2][6], const complex<double> nm,
          double A[8][8]) {
  g = 0.0;
  for (size_t i0 = 0; i0 < 6; i0++) {
    q[0][i0] = 0.0;
    q[1][i0] = 0.0;
  }

  complex<double> cpdf[8];

  // Convention: 0 = g, 1 = d valence, 2 = u valence, 3 = d sea, 4 = s sea,
  // 5 = b sea, 6 = u sea, 7 = c sea.
  // f = A0 * x^A1 * (1 - x)^A2 * ( 1 + A3 * x^(1/2) + A4 * x + A5 * x^(3/2)
  // + A6 * x^2 + A7 * x^(5/2) )
  for (size_t i0 = 0; i0 < 8; i0++) {
    complex<double> y = A[i0][2] + 1.0;
    cpdf[i0] = A[i0][0] * Gamma(y) *
               (beta_over_gamma(A[i0][1] + nm, y) +
                A[i0][3] * beta_over_gamma(A[i0][1] + nm + 0.5, y) +
                A[i0][4] * beta_over_gamma(A[i0][1] + nm + 1.0, y) +
                A[i0][5] * beta_over_gamma(A[i0][1] + nm + 1.5, y) +
                A[i0][6] * beta_over_gamma(A[i0][1] + nm + 2.0, y) +
                A[i0][7] * beta_over_gamma(A[i0][1] + nm + 2.5, y));
  }
  // See conventions in pdfX function.
  g = cpdf[0];
  q[0][0] = cpdf[1] + cpdf[3];
  q[1][0] = cpdf[3];
  q[0][1] = cpdf[4];
  q[1][1] = cpdf[4];
  q[0][2] = cpdf[5];
  q[1][2] = cpdf[5];
  q[0][3] = cpdf[2] + cpdf[6];
  q[1][3] = cpdf[6];
  q[0][4] = cpdf[7];
  q[1][4] = cpdf[7];
}

void pdfEvolve(complex<double> &gg, complex<double> qq[2][6],
               complex<double> nm, complex<double> lmbd) {
  const int nf = NFLVR;
  const double dnf = (double)nf;
  const double beta0 = 5.5 - dnf / 3.0;

  // Anomalous dimensions.
  // Conventions: \alpha_S/2\pi
  complex<double> Pqq, Pqg, Pgq, Pgg;
  Pqq = 4.0 / 3.0 *
        (1.5 + 1.0 / nm / (nm + 1.0) - 2.0 * (Psi(nm + 1.0) + M_EULER));
  Pqg = dnf * (2.0 + nm + pow2(nm)) / nm / (nm + 1.0) / (nm + 2.0);
  Pgq = 4.0 / 3.0 * (2.0 + nm + pow2(nm)) / (pow(nm, 3.0) - nm);
  Pgg = 6.0 * (1.0 / (nm + 2.0) / (nm + 1.0) + 1.0 / nm / (nm - 1.0) -
               (Psi(nm + 1.0) + M_EULER)) +
        beta0;

  // Non singlet decomposition.
  complex<double> NS[6], VA[6];
  complex<double> tmp(0.0, 0.0);
  double di0 = 0.0;

  for (int i0 = 0; i0 < 6; i0++) {
    NS[i0] = complex<double>(0.0, 0.0);
    VA[i0] = complex<double>(0.0, 0.0);

    // Flavor number tests.
    if (nf < 6 && i0 == 5) {
      continue;
    }
    if (nf < 5 && i0 == 2) {
      continue;
    }

    // Definitions.
    di0 += 1.0;
    tmp += qq[0][i0] + qq[1][i0];
    NS[i0] = tmp - di0 * (qq[0][i0] + qq[1][i0]);
    VA[i0] = qq[0][i0] - qq[1][i0];

    // Evolution of PDFs.
    VA[i0] *= pow(1.0 - lmbd, Pqq / beta0);
    NS[i0] *= pow(1.0 - lmbd, Pqq / beta0);
  }

  complex<double> rp, rm;
  rp = (Pgg + Pqq + sqrt(pow(Pgg - Pqq, 2) + 4.0 * Pgq * Pqg)) * 0.5;
  rm = (Pgg + Pqq - sqrt(pow(Pgg - Pqq, 2) + 4.0 * Pgq * Pqg)) * 0.5;

  // Singlet decomposition and evolution
  complex<double> SI, GL;
  SI = ((Pqq - rp) * tmp + Pqg * gg) / (rm - rp) * pow(1.0 - lmbd, rm / beta0) +
       ((Pqq - rm) * tmp + Pqg * gg) / (rp - rm) * pow(1.0 - lmbd, rp / beta0);
  GL = (Pgq * tmp + (Pgg - rp) * gg) / (rm - rp) * pow(1.0 - lmbd, rm / beta0) +
       (Pgq * tmp + (Pgg - rm) * gg) / (rp - rm) * pow(1.0 - lmbd, rp / beta0);

  // Recover PDFs
  di0 = dnf;
  tmp = complex<double>(0.0, 0.0);
  for (int i0 = 5; i0 >= 0; i0--) {
    // Flavor number tests
    if (nf < 6 && i0 == 5) {
      continue;
    }
    if (nf < 5 && i0 == 2) {
      continue;
    }

    qq[0][i0] = (1.0 / dnf * SI - 1.0 / di0 * NS[i0] + tmp + VA[i0]) * 0.5;
    qq[1][i0] = (1.0 / dnf * SI - 1.0 / di0 * NS[i0] + tmp - VA[i0]) * 0.5;
    tmp += 1.0 / di0 / (di0 - 1.0) * NS[i0];
    di0 -= 1.0;
  }

  gg = GL;
}

// Struct for GSL fit.
struct data {
  size_t n;
  double *y;
  double *sigma;
  double xr;
};

int expb_f(const gsl_vector *x, void *data, gsl_vector *f) {
  size_t n = ((struct data *)data)->n;
  double *y = ((struct data *)data)->y;
  double *sigma = ((struct data *)data)->sigma;
  double xr = ((struct data *)data)->xr;

  double A[8];
  for (size_t i0 = 0; i0 < 8; i0++) {
    A[i0] = gsl_vector_get(x, i0);
  }

  for (size_t i0 = 0; i0 < n; i0++) {
    // double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
    double t = sampling(xr, 1.0, n, i0);
    double Yi = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2]) *
                (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5) +
                 A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5));
    gsl_vector_set(f, i0, (Yi - y[i0]) / sigma[i0]);
  }

  return GSL_SUCCESS;
}

int expb_df(const gsl_vector *x, void *data, gsl_matrix *J) {
  size_t n = ((struct data *)data)->n;
  double *sigma = ((struct data *)data)->sigma;
  double xr = ((struct data *)data)->xr;

  double A[8];
  for (size_t i0 = 0; i0 < 8; i0++) {
    A[i0] = gsl_vector_get(x, i0);
  }

  for (size_t i0 = 0; i0 < n; i0++) {
    // double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
    double t = sampling(xr, 1.0, n, i0);
    double s = sigma[i0];
    double e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2]) *
               (1.0 + A[3] * sqrt(t) + A[4] * t + A[5] * pow(t, 1.5) +
                A[6] * pow(t, 2.0) + A[7] * pow(t, 2.5)) /
               s;
    gsl_matrix_set(J, i0, 0, e / A[0]);
    gsl_matrix_set(J, i0, 1, e * log(t));
    gsl_matrix_set(J, i0, 2, e * log(1.0 - t));

    e = A[0] * pow(t, A[1]) * pow(1.0 - t, A[2]) / s;
    gsl_matrix_set(J, i0, 3, e * sqrt(t));
    gsl_matrix_set(J, i0, 4, e * t);
    gsl_matrix_set(J, i0, 5, e * pow(t, 1.5));
    gsl_matrix_set(J, i0, 6, e * pow(t, 2.0));
    gsl_matrix_set(J, i0, 7, e * pow(t, 2.5));
  }

  return GSL_SUCCESS;
}

int expb_fdf(const gsl_vector *x, void *data, gsl_vector *f, gsl_matrix *J) {
  expb_f(x, data, f);
  expb_df(x, data, J);

  return GSL_SUCCESS;
}

void Fit(double A[8], double E[8], int flag, double xr, double Q2,
         double weight_valence, double weight_sea, double weight_gluon,
         double xmin) {
  // Defines the function to minimize.
  const size_t n = NDATA;
  const size_t p = 8;

  // we set xr to xmin
  xr = xmin;

  double y[NDATA], sigma[NDATA];
  struct data d = {n, y, sigma, xr};

  gsl_multifit_function_fdf f;
  f.f = &expb_f;
  f.df = &expb_df;
  f.fdf = &expb_fdf;
  f.n = n;
  f.p = p;
  f.params = &d;

  // PDFs to fit.
  double gamma_weight = -1.6;
  for (size_t i0 = 0; i0 < n; i0++) {
    // double t = pow(xr, 1.0 - (1.0 + (double)i0) / (1.0 + (double)n));
    double t = sampling(xr, 1.0, n, i0);
    double q[2][6];
    double g;
    double qq;
    pdfX(g, q, t, Q2);
    switch (flag) {
    case 0:
      qq = g;
      gamma_weight = weight_gluon;
      break;
    case 1:
      qq = q[0][0] - q[1][0];
      gamma_weight = weight_valence;
      break;
    case 2:
      qq = q[0][3] - q[1][3];
      gamma_weight = weight_valence;
      break;
    case 3:
      qq = q[1][0];
      gamma_weight = weight_sea;
      break;
    case 4:
      qq = q[1][1];
      gamma_weight = weight_sea;
      break;
    case 5:
      qq = q[1][2];
      gamma_weight = weight_sea;
      break;
    case 6:
      qq = q[1][3];
      gamma_weight = weight_sea;
      break;
    case 7:
      qq = q[1][4];
      gamma_weight = weight_sea;
      break;
    default:
      fprintf(stderr, "error: while retrieving PDF: unkown flag %d\n", flag);
      exit(1);
    }
    y[i0] = qq;
    sigma[i0] = pow(
        t,
        gamma_weight); // DeFlorian-like pdf-weights t^(-1.6); DeJonathan t^(-1)
  }

  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver *s;
  T = gsl_multifit_fdfsolver_lmsder;
  s = gsl_multifit_fdfsolver_alloc(T, n, p);

  // First guess for the parameters.
  double x_init[8] = {1.0, -1.4, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  if (flag == 1 || flag == 2) {
    x_init[1] = -0.6;
  }

  gsl_vector_view x = gsl_vector_view_array(x_init, p);
  gsl_multifit_fdfsolver_set(s, &f, &x.vector);

  int status = GSL_CONTINUE;
  unsigned int iter = 0;

  // Does the fit.
  while (status == GSL_CONTINUE && iter < 500000) {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(s);
    if (status) {
      break;
    }
    status = gsl_multifit_test_delta(s->dx, s->x, 1.0e-12, 1.0e-5);
  }

  gsl_matrix *covar = gsl_matrix_alloc(p, p);
#if GSL_MAJOR_VERSION > 1
  {
    gsl_matrix *J = gsl_matrix_alloc (n, p);;
    gsl_multifit_fdfsolver_jac(s, J);
    gsl_multifit_covar (J, 0.0, covar);
    gsl_matrix_free(J);
  }
#else
  gsl_multifit_covar(s->J, 0.0, covar);
#endif
  // double c = pow(gsl_blas_dnrm2(s->f), 2) / ((double)(n - p));
  // c = (c > 1.0 ? c : 1.0);
  double chi = gsl_blas_dnrm2(s->f);
  double dof = n - p;
  double c = GSL_MAX_DBL(1, chi / sqrt(dof));

  // Gets parameters.
  for (size_t i0 = 0; i0 < 8; i0++) {
    A[i0] = 0.0;
  }
  for (size_t i0 = 0; i0 < p; i0++) {
    A[i0] = gsl_vector_get(s->x, i0);
  }
  for (size_t i0 = 0; i0 < p; i0++) {
    E[i0] = c * gsl_matrix_get(covar, i0, i0);
  }

  printf("#chisq/dof = %g\n", pow(chi, 2.0) / dof);
  printf("#status = %s\n", gsl_strerror(status));

  gsl_multifit_fdfsolver_free(s);
  gsl_matrix_free(covar);
}

void pdfFit(double &A1MIN, double A[8][8], double tau, double Q2,
            double weight_valence, double weight_sea, double weight_gluon,
            double xmin) {
  const int nf = NFLVR;
  fprintf(stderr, "Performing PDF fit with %d flavors with M^2/S = %g, Q^2 = "
                  "%g\n and weights: valence: x^%g, sea: x^%g, gluon: x^%g and "
                  "xmin = %g \n Fit function: f = A0 * x^A1 * (1 - x)^A2 * ( 1 "
                  "+ A3 * x^(1/2) + A4 * x + A5 * x^(3/2) + A6 * x^2 + A7 * "
                  "x^(5/2) )\n",
          nf, tau, Q2, weight_valence, weight_sea, weight_gluon, xmin);

  for (int i0 = 0; i0 < 8; i0++) {
    double err[8];
    switch (i0) {
    case 0:
      fprintf(stderr, "Fitting gluon PDF...");
      break;
    case 1:
      fprintf(stderr, "Fitting valence down quark PDF...");
      break;
    case 2:
      fprintf(stderr, "Fitting valence up quark PDF...");
      break;
    case 3:
      fprintf(stderr, "Fitting sea down quark PDF...");
      break;
    case 4:
      fprintf(stderr, "Fitting strange quark PDF...");
      break;
    case 5:
      fprintf(stderr, "Fitting bottom quark PDF...");
      break;
    case 6:
      fprintf(stderr, "Fitting sea up quark PDF...");
      break;
    case 7:
      fprintf(stderr, "Fitting charm quark PDF...");
      break;
    }
    fflush(stderr);
    // Fit(A[i0], err, i0, tau, Q2); // tau = M^2/S and Q2 = mu_F^2.
    Fit(A[i0], err, i0, tau, Q2, weight_valence, weight_sea, weight_gluon,
        xmin);
    fprintf(stderr, " done.\n");
    fprintf(stderr, "Fit result:\n");

    for (int i1 = 0; i1 < 8; i1++) {
      fprintf(stderr, "A%i =  %.5f # +-%.5f\n", i1, A[i0][i1], err[i1]);
    }
    if (A[i0][1] < A1MIN) {
      A1MIN = A[i0][1];
    }
  }
}
