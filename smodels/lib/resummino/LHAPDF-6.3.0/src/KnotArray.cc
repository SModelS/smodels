// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/KnotArray.h"
#include <functional>

namespace LHAPDF {


  size_t KnotArray1F::_mkhash(const std::vector<double>& xx) const {
    std::hash<double> hasher;
    size_t rtn = 0;
    for (double x : xx) rtn = 31*rtn + hasher(x);
    return rtn + 1;
  }


}
