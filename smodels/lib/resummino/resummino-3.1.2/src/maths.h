// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

// Mathematical functions used by other modules.

#ifndef MATH_H_
#define MATH_H_

#include <cmath>
#include <complex>
#include <cstdlib>

using namespace std;

// Returns the gamma function evaluated at x.
complex<double> Gamma(const complex<double> x);

// Returns the Psi function evaluated at x.
complex<double> Psi(complex<double> x);

// Returns the Beta function evaluated at (x, y).
complex<double> Beta(const complex<double> x, const complex<double> y);

// Returns beta(x, y) / gamma(y). This exists for performance reasons, since
// beta(x, y) = gamma(x) gamma(y) / gamma(x + y), so we can factor out gamma(y)
// in pdfN.
static inline complex<double> beta_over_gamma(const complex<double> x,
                                              const complex<double> y) {
  return Gamma(x) / Gamma(x + y);
}

#endif
