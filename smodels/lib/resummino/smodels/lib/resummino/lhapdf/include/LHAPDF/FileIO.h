// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2020 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once
#ifndef LHAPDF_FileIO_H
#define LHAPDF_FileIO_H

// STL includes
#include <fstream>
#include <sstream>
#include <string>

/// Namespace for all LHAPDF functions and classes
namespace LHAPDF {


  /// @brief MPI-safe file I/O interface
  template <class FILETYPE>
  class File {
  public:

    // Constructor
    File(const std::string& name)
      : _name(name), _fileptr(nullptr), _streamptr(nullptr) {
      open();
    }

    /// Destructor
    ~File() { close(); }


    /// Open file
    bool open();

    /// Close file
    bool close();

    /// Forward methods via pointer emulation
    FILETYPE* operator->() const { return _fileptr; }

    /// Forward methods via pointer-dereference emulation
    FILETYPE& operator*() const { return *_fileptr; }

    /// Get the file content
    std::string getContent() const { return _streamptr != nullptr ? _streamptr->str() : ""; }


  protected:

    std::string _name;

    FILETYPE* _fileptr;

    std::stringstream* _streamptr;

  };


  // Convenience aliases (NB. have to use #defines due to how .cc template building works)
  // typedef File<std::ifstream> IFile;
  // typedef File<std::ofstream> OFile;
  #define IFile File<std::ifstream>
  #define OFile File<std::ofstream>


  /// Global function to flush the MPI-safe file cache
  void flushFileCache();


}
#endif
