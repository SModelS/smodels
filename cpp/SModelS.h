#ifndef SModelS_H
#define SModelS_H

#include <string>

class SModelS {
  /** 
   * \class SModelS
   * a simple class that wraps python-based SModelS for use within a C++ environment
   */

  public:
    /*** construction from a parameter.ini file.
     *   \paramname verbosity set the verbosity of SModelS. Values are:
     *   debug, info, warn, error */
    SModelS( const std::string & parameterfile, const std::string & verbosity );
    SModelS( const std::string & parameterfile );
    ~SModelS();
    /** run over a single slha file. */
    int run ( const std::string & slhafile );

  private:
    void loadDatabase ( const std::string & parameterfile );
    void initialize ( const std::string & parameterfile, 
                      const std::string & verbosity );
};

#endif /* SModelS_H */
