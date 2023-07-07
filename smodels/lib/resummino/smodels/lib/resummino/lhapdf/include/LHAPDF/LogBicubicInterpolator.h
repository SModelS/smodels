// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_LogBicubicInterpolator_H
#define LHAPDF_LogBicubicInterpolator_H

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// @brief Implementation of bicubic interpolation
  ///
  /// This class will interpolate in 2D using a bicubic hermite spline.
  class LogBicubicInterpolator : public Interpolator {
  public:

    /// Implementation of (x,Q2) interpolation
    double _interpolateXQ2(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) const;

    /// @brief A single set of cached x-variables
    struct XCache {
      /// Defining params from call (initialised to unphysical values, so first use will set the cache)
      double x = -1;
      //size_t ix;

      /// Cached params
      double logx;
      double dlogx_1;
      double tlogx;
    };

    /// A multi-level x-variable cache for a single subgrid hash
    struct XCaches {

      /// Number of cache levels
      static size_t SIZE;
      /// Cache access strategy
      static int UPDATE_STEP;
      /// Cache access strategy
      static bool UPDATE_ON_HIT;

      /// (Re)define cache size and search strategy
      static void setup(size_t size, int update_step=+1, bool update_on_hit=true);
      /// Initialize a cache on this thread (call with explicit locks to ensure safe initialisation for each thread)
      static void init();

      /// Latest-call index
      size_t ilast = 0;

      /// List of N cached-value sets
      vector<XCache> caches{SIZE};

      /// Get the length of the cache vector
      size_t size() { return caches.size(); }

      /// Access the @a index'th element of the cache
      XCache& operator[] (size_t index) { return caches[index]; }

    };


    /// @brief A single set of cached Q2-variables
    struct Q2Cache {
      /// Defining params from call (initialised to unphysical values, so first use will set the cache)
      double q2 = -1;

      /// Cached params
      double logq2;
      double dlogq_0;
      double dlogq_1;
      double dlogq_2;
      double tlogq;
    };

    /// A multi-level Q2-variable cache for a single subgrid hash
    struct Q2Caches {

      /// Number of cache levels
      static size_t SIZE;
      /// Cache access strategy
      static int UPDATE_STEP;
      /// Cache access strategy
      static bool UPDATE_ON_HIT;

      /// (Re)define cache size and search strategy
      static void setup(size_t size, int update_step=+1, bool update_on_hit=true);
      /// Initialize a cache on this thread (call with explicit locks to ensure safe initialisation for each thread)
      static void init();

      /// Latest-call index
      size_t ilast = 0;

      /// List of N cached-value sets
      vector<Q2Cache> caches{SIZE};

      /// Get the length of the cache vector
      size_t size() { return caches.size(); }

      /// Access the @a index'th element of the cache
      Q2Cache& operator[] (size_t index) { return caches[index]; }

    };


    /// @brief Get and update the current caching structs for interpolation params
    ///
    /// @note Caches are handled separately for x and Q since they can be sampled very differently.
    ///@{
    static XCache& _getCacheX(const KnotArray1F& subgrid, double x, size_t ix);
    static Q2Cache& _getCacheQ2(const KnotArray1F& subgrid, double q2, size_t iq2);
    ///@}

  };


}

#endif
