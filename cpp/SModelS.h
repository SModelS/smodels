#ifndef SModelS_H
#define SModelS_H

#include <string>

class SModelS {
  /** 
   * \class SModelS
   * a simple class that wraps python-based SModelS for use with a C++ environment
   */

  public:
    SModelS( const std::string & parameterfile );
    ~SModelS();
    int run ( const std::string & file );

  private:
    void loadDatabase( const std::string & parameterfile );
};

#endif /* SModelS_H */
