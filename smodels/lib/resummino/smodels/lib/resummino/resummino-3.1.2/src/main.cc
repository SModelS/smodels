#include "resummino.h"
int main(int argc, char *argv[]) {
  print_banner();
  Args args;
  parse_args(argc, argv, args);
  Results r = compute(args);
  print_results(r, args);
}