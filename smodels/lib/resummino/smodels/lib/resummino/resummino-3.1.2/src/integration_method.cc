// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

#include "gsl_all.h"
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "params.h"
#include "utils.h"

// This function uses the GSL Monte Vegas integration routine
void Integration(double (*IF)(double *x, size_t dim, void *jj), int ndim, int maxeval,
                 double epsrel, double epsabs, double &res, double &err, void *userdata,
                 double calls_warmupfactor = 0.1, int iterations_warmup = 5,
                 int iterations_main = 3) {

  // gsl_vegas
  // Selects the integrand.
  size_t calls = 0;
  gsl_monte_function I;
  I.f = IF;
  I.dim = ndim;
  I.params = userdata;
  calls = maxeval; // set number of calls
  const size_t dnum = I.dim;

  // Integration limits [0:1]
  double xmin[dnum], xmax[dnum];
  for (size_t i0 = 0; i0 < dnum; i0++) {
    xmin[i0] = 0.0;
    xmax[i0] = 1.0;
  }

  // Initializes the integration routine.
  const gsl_rng_type *T;
  gsl_rng *r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  // seed gsl non default integration
  // T = gsl_rng_mt19937;
  r = gsl_rng_alloc(T);

  gsl_monte_vegas_state *s = gsl_monte_vegas_alloc(dnum);

  // Parameter block
  // make adjustments if you want to change some parameters
  gsl_monte_vegas_params params;
  gsl_monte_vegas_params_get(s, &params);
  params.alpha = 1.5;
  params.mode = GSL_VEGAS_MODE_IMPORTANCE_ONLY;
  gsl_monte_vegas_params_set(s, &params);

  // -1 -> no output; 0 -> summary information; 1 and 2 -> prints much more
  // information
  s->verbose = 0;

  // Integration warm-up.
  // (stage = 0 which begins with a new uniform grid and empty weighted average)
  s->stage = 0;
  s->iterations = iterations_warmup;
  gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls * calls_warmupfactor, r, s, &res, &err);

  // Integrates
  // Calling VEGAS with stage = 1 retains the grid from the previous run but
  // discards the weighted average,
  // so that one can “tune” the grid using a relatively small number of points
  // and then do a large run with stage = 1
  s->stage = 1;
  s->iterations = iterations_main;
  gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls, r, s, &res, &err);

  // Convergence step.
  // perform more iterations (<= maxiters) to improve convergence
  // tries to achieve an appropriate chisquared value
  double prec = abs(err / res);

  for (int counter = 1; counter <= ((Parameters *)userdata)->max_iters > 0; counter++) {
    // if result is zero (e.g. LO pt) relative precision does not make any sense
    if (abs(res) < 1.0e-15 || (prec < epsrel && fabs(gsl_monte_vegas_chisq(s) - 1.0) < 0.5) ||
        prec < epsrel / 2.5 || (abs(err) < epsabs && fabs(gsl_monte_vegas_chisq(s) - 1.0) < 0.5)) {
      break;
    }
    s->iterations = 1;
    s->stage = 3;

    gsl_monte_vegas_integrate(&I, xmin, xmax, dnum, calls, r, s, &res, &err);

    prec = abs(err / res);
  }

  gsl_monte_vegas_free(s);
  gsl_rng_free(r);
}
