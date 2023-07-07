// This file is part of Resummino.
//
// Copyright 2008-2011 Jonathan Debove.
// Copyright 2011-2016 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.

void Integration(double (*IF)(double *x, size_t dim, void *jj), int ndim,
                 int maxeval, double epsrel, double epsabs, double &res,
                 double &err, void *userdata, double calls_warmupfactor = 0.1,
                 int iterations_warmup = 5, int iterations_main = 3);
