// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#include "LHAPDF/Info.h"
#include "LHAPDF/PDFIndex.h"
#include "LHAPDF/FileIO.h"

#include "yaml-cpp/yaml.h"
#ifdef YAML_NAMESPACE
#define YAML YAML_NAMESPACE
#endif

namespace LHAPDF {


  void Info::load(const string& filepath) {
    // Complain if the path is empty
    if (filepath.empty()) throw ReadError("Empty PDF file name given to Info::load");

    // But complain if a non-empty path is provided, but it's invalid
    if (!file_exists(filepath)) throw ReadError("PDF data file '" + filepath + "' not found");

    // Read the YAML part of the file into the metadata map
    try {

      // Do the parsing "manually" up to the first doc delimiter
      IFile file(filepath.c_str());
      string docstr, line;
      while (getline(*file, line)) {
        if (line == "---") break;
        docstr += line + "\n";
      }
      YAML::Node doc = YAML::Load(docstr);
      for (YAML::const_iterator it = doc.begin(); it != doc.end(); ++it) {
        const string key = it->first.as<string>();
        YAML::Emitter em;
        em << it->second;
        const string val = em.c_str();
        _metadict[key] = val;
        //
        // const YAML::Node& val = it->second;
        // if (val.IsScalar()) {
        //   // Scalar value
        //   _metadict[key] = val.as<string>();
        // } else {
        //   // Process the sequence entries into a comma-separated string
        //   /// @todo Is there a better way? Can maybe use std::any storage in the metadict?
        //   string seqstr = "";
        //   for (size_t i = 0; i < val.size(); ++i)
        //     seqstr += val[i].as<string>() + ((i < val.size()-1) ? "," : "");
        //   _metadict[key] = seqstr;
        // }
      }

    } catch (const YAML::ParserException& ex) {
      throw ReadError("YAML parse error in " + filepath + " :" + ex.what());
    } catch (const LHAPDF::Exception& ex) {
      throw;
    } catch (const std::exception& ex) {
      throw ReadError("Trouble when reading " + filepath + " :" + ex.what());
    }

  }



  // /// @todo Only support loading via PDF set name and member ID, not explicit paths
  // /// @todo Replace the loading of the set metadata into the member info with set-level Info singletons
  // void Info::loadFull(const string& mempath) { //< @todo Need a better method name!
  //   // Extract the set name from the member data file path
  //   const string memberdata = findFile(mempath);
  //   if (memberdata.empty() || !file_exists(memberdata)) throw ReadError("Could not find PDF data file '" + mempath + "'");
  //   const string memname = basename(memberdata); //< Can use this to alternatively work out the set name...
  //   const string setdir = dirname(memberdata);
  //   const string setname = basename(setdir);
  //   path setinfo = findpdfsetinfopath(setname);
  //   // Load the set info
  //   if (file_exists(setinfo)) load(setinfo);
  //   // Load the member info (possibly overriding the set-level metadata)
  //   load(memberdata);
  // }


}
