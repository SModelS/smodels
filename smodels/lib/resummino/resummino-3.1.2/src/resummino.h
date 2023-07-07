#ifndef RESUMMINO
#define RESUMMINO
#include "params.h"
#include <map>
#include <string>

struct Result {
  double res;
  double err;
  Result(double r, double e) : res(r), err(e){};
  Result operator+(const Result &r) const {
    return Result(r.res + res, sqrt(r.err * r.err + err * err));
  }
};

struct Results {
  Result LO;
  Result NLO;
  Result NLL;
  Result aNLLj;
  Results(Result lo, Result nlo, Result nll, Result anll)
      : LO(lo), NLO(nlo), NLL(nll), aNLLj(anll){};
  Results operator+(const Results &r) const {
    return Results(r.LO + LO, r.NLO + NLO, r.NLL + NLL, r.aNLLj + aNLLj);
  }
};
enum ResultType { total, pt, m, ptj };

class Args {
public:
  string input_file = "resummino.in";
  string log_file = "";
  int stop_after_lo = 0;
  int stop_after_nlo = 0;
  int nll_unimproved = 0;
  int nnll_unimproved = 0;
  ResultType result_type;
  Parameters *params;
  map<string, string> arguments;
  map<string, string> config;
  Args() { params = new Parameters(); };
  ~Args() { delete params; };
};

void print_banner();
void print_results(Results &r, Args &args);
Results compute(Args &args);
void done();
void init(Args &args);
void parse_args(int argc, char *argv[], Args &args);
int set_particles(int pdg_p1, int pdg_p2, int *out1, int *out2);
#endif