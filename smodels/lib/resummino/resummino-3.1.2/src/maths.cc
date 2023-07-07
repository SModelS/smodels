// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2014 David R. Lamprea.
// Copyright 2011-2014 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Implements the mathematical module.

#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

using namespace std;

complex<double> Gamma(const complex<double> x) {
  const int intx = (int)real(x);
  const double n = -(double)(intx < 0 ? -intx : intx);
  if (real(x) == n && imag(x) == 0.0) {
    cout << "Gamma(" << x << ") undefined\n";
    exit(0);
  }

  // Works with Re(xx) > 0
  const complex<double> xx = (real(x) > 0.0 ? x : -x);

  // Magic numbers for Gamma function
  const double q0 = 75122.6331530;
  const double q1 = 80916.6278952;
  const double q2 = 36308.2951477;
  const double q3 = 8687.24529705;
  const double q4 = 1168.92649479;
  const double q5 = 83.8676043424;
  const double q6 = 2.50662827511;

  complex<double> gamma =
      exp((xx + .5) * log(xx + 5.5) - xx - 5.5) *
      (q0 + q1 * xx + q2 * pow(xx, 2) + q3 * pow(xx, 3) + q4 * pow(xx, 4) +
       q5 * pow(xx, 5) + q6 * pow(xx, 6)) /
      xx / (xx + 1.0) / (xx + 2.0) / (xx + 3.0) / (xx + 4.0) / (xx + 5.0) /
      (xx + 6.0);

  return (x == xx ? gamma : -M_PI / xx / sin(M_PI * xx) / gamma);
}

complex<double> Psi(complex<double> x) {
  complex<double> result(0.0, 0.0);
  while (real(x) < 10.0) {
    result -= 1.0 / x;
    x += 1.0;
  }
  result +=
      log(x) - .5 / x -
      pow(x, -2) / 2520. * (210.0 + pow(x, -2) * (-21.0 + pow(x, -2) * 10.0));

  return result;
}

complex<double> Beta(const complex<double> x, const complex<double> y) {
  return Gamma(x) * Gamma(y) / Gamma(x + y);
}
