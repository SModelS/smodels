// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Config.h"
#include "LHAPDF/Version.h"
using namespace std;

namespace LHAPDF {


  Config& Config::get() {
    static Config _cfg; //< Could we use the Info(path) constructor for automatic init-once behaviour?
    // Test for emptiness and only initialise *once*:
    if (_cfg._metadict.empty()) {
      std::string confpath = findFile("lhapdf.conf");
      if (!confpath.empty()) _cfg.load(confpath);
    }
    return _cfg;
  }


  Config::~Config() {
    // Emit citation information at the end of the job, via the Config destructor
    // std::cout << "CONFIG DESTRUCTION" << std::endl;
    if (verbosity() > 0) {
      cout << "Thanks for using LHAPDF " << version() << ". Please make sure to cite the paper:\n";
      cout << "  Eur.Phys.J. C75 (2015) 3, 132  (http://arxiv.org/abs/1412.7420)" << endl;
    }
  }


}
