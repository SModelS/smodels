// This file is part of Resummino.
//
// Copyright 2008-2010 Jonathan Debove.
// Copyright 2011-2015 David R. Lamprea.
// Copyright 2011-2016 Marcel Rothering.
//
// Licensed under the terms of the EUPL version 1.1 or later.
// See the LICENCE file for more information.
//
// Command-line interface.

#include <cstdlib>
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <cmath>
#include <complex>
#include <iostream>
#include <map>
#include <string>

#include "clooptools.h"
#include "gsl_all.h"
#include "hxs.h"
#include "kinematics.h"
#include "params.h"
#include "pdf.h"
#include "pxs_gausq.h"
#include "resummino.h"
#include "utils.h"

using namespace std;

// will print only the expansion of the resummation up to NLO
//#define EXPANSION

// fixed scale for invariant mass distribution
//#define FIXED_SCALE

#define UNITS (args.result_type == total ? "pb" : "pb/GeV")
// Gets option from command-line if set, from input file otherwise.
#define get_option(name) (args.arguments[name] == "" ? args.config[name] : args.arguments[name])

// Simple implementation of get_current_dir_name for system w/o GNU extensions.
static char *my_get_current_dir_name(void) {
  size_t output_size = 100;
  char *output = (char *)malloc(output_size);
  output = getcwd(output, output_size);
  while (output == NULL && errno == ERANGE) {
    output_size += 100;
    output = getcwd(output, output_size);
  }
  return output;
}

void print_banner() {
  printf("\n"
         "  Welcome to Resummino " RESUMMINO_VERSION "\n"
         "  Copyright 2008-2010 Jonathan Debove.\n"
         "  Copyright 2011-2014 David R. Lamprea.\n"
         "  Copyright 2011-2016 Marcel Rothering.\n"
         "\n"
         "  Licensed under the EUPL 1.1 or later.\n"
         "  See https://resummino.hepforge.org/ for more information.\n"
         "\n"
         "  For general Resummino usage, please cite:\n"
         "  - B. Fuks, M. Klasen, D. R. Lamprea, and M. Rothering, Eur. Phys. J. C 73, 2480 (2013).\n"
         "\n"
         "  For slepton pair production, please cite:\n"
         "  - G. Bozzi, B. Fuks, M. Klasen, PRD 74, 015001 (2006); NPB 777, "
         "157 (2007); NPB 794, 46 (2007).\n"
         "  - J. Fiaschi, and M. Klasen, Phys.Rev.D 98, 5, 055014 (2018).\n"
         "  - J. Fiaschi, M. Klasen, and M. Sunder, JHEP 04, 049 (2020).\n"
         "\n"
         "  For gaugino pair production, please cite:\n"
         "  - J. Debove, B. Fuks, M. Klasen, PLB 688, 208 (2010); NPB 842, 51 "
         "(2011); NPB 849, 64 (2011).\n"
         "  - B. Fuks, M. Klasen, D. R. Lamprea, M. Rothering, JHEP 1210, 081 "
         "(2012).\n"
         "  - J. Fiaschi, and M. Klasen, JHEP 03, 094 (2018); Phys.Rev.D 102, 9, 095021 (2020).\n"
         "\n"
         "  For associated gaugino-gluino production, please cite:\n"
         "  - B. Fuks, M. Klasen, M. Rothering, JHEP 07:053, (2016).\n"
         "\n"
         "  For associated gaugino-squark production, please cite:\n"
         "  - J. Fiaschi, B. Fuks, M. Klasen, A. Neuwirth, JHEP 06, 130 (2022).\n"
         "\n"
         "  For W' and Z' processes, please cite:\n"
         "  - T. Jezo, M. Klasen, D. Lamprea, F. Lyonnet, I. Schienbein, JHEP 12, 092 (2014).\n"
         "\n");
}

static int pdg_to_internal_id(int pdg_code) {
  switch (pdg_code) {
  case 1000012: // ~nu_eL
    return 10;
    break;
  case 1000014: // ~nu_muL
    return 11;
    break;
  case 1000016: // ~nu_tau
    return 12;
    break;
  case 1000011: // ~e_L
    return 13;
    break;
  case 2000011: // ~e_R
    return 14;
    break;
  case 1000013: // ~mu_L
    return 15;
    break;
  case 2000013: // ~mu_R
    return 16;
    break;
  case 1000015: // ~tau_1
    return 17;
    break;
  case 2000015: // ~tau_2
    return 18;
    break;
  case 1000022: // ~chi_10
    return 0;
    break;
  case 1000023: // ~chi_20
    return 1;
    break;
  case 1000025: // ~chi_30
    return 2;
    break;
  case 1000035: // ~chi_40
    return 3;
    break;
  case 1000024: // ~chi_1+
    return 4;
    break;
  case 1000037: // ~chi_2+
    return 5;
    break;
  case 12: // nu_eL
    return 20;
    break;
  case 14: // nu_muL
    return 21;
    break;
  case 16: // nu_tauL
    return 22;
    break;
  case 11: // eL
    return 23;
    break;
  case 13: // muL
    return 24;
    break;
  case 15: // tauL
    return 25;
    break;
  case 1000001: // ~d_L
    return 31;
    break;
  case 2000001: // ~d_R
    return 34;
    break;
  case 1000002: // ~u_L
    return 37;
    break;
  case 2000002: // ~u_R
    return 40;
    break;
  case 1000003: // ~s_L
    return 32;
    break;
  case 2000003: // ~s_R
    return 35;
    break;
  case 1000004: // ~c_L
    return 38;
    break;
  case 2000004: // ~c_R
    return 41;
    break;
  case 1000005: // ~b_1
    return 33;
    break;
  case 2000005: // ~b_2
    return 36;
    break;
  case 1000006: // ~t_1
    return 39;
    break;
  case 2000006: // ~t_2
    return 42;
    break;
  case 1000021: // ~g
    return 30;
    break;
  default:
    fprintf(stderr, "error: Wrong PDG code.\n");
    exit(1);
  }
}

// translates the conventions used in resummino.in
// to code conventions where p1 must be positive/neutral and p2 must be
// negative/neutral
int set_particles(int pdg_p1, int pdg_p2, int *out1, int *out2) {
  int charge_p1 = pdg_particle_charge(pdg_p1);
  int charge_p2 = pdg_particle_charge(pdg_p2);
  int internal_p1 = pdg_to_internal_id(abs(pdg_p1));
  int internal_p2 = pdg_to_internal_id(abs(pdg_p2));

  if ((is_pdg_gaugino(pdg_p1) && is_pdg_squark(pdg_p2)) ||
      (is_pdg_gaugino(pdg_p2) && is_pdg_squark(pdg_p1))) {
    if ((pdg_p1 > 0 && is_pdg_squark(pdg_p1) && is_pdg_neutralino(pdg_p2)) ||
        // anti case is swapped
        (pdg_p2 < 0 && is_pdg_squark(pdg_p2) && is_pdg_neutralino(pdg_p1))) {
      *out1 = internal_p1;
      *out2 = internal_p2;
    } else if ((pdg_p2 > 0 && is_pdg_squark(pdg_p2) && is_pdg_neutralino(pdg_p1)) ||
               // anti case is swapped
               (pdg_p1 < 0 && is_pdg_squark(pdg_p1) && is_pdg_neutralino(pdg_p2))) {
      *out1 = internal_p2;
      *out2 = internal_p1;
    }

    // down squarks chargino+
    // and up squarks chargino-
    else if ((
                 // down squarks
                 pdg_p1 > 0 && is_pdg_down_type_squark(pdg_p1) &&
                 // charginos+
                 pdg_p2 > 0 && is_pdg_chargino(pdg_p2)) ||
             (
                 // anti-down squarks
                 pdg_p2 < 0 && is_pdg_down_type_squark(pdg_p2) &&
                 // charginos-
                 pdg_p1 < 0 && is_pdg_chargino(pdg_p1)) ||
             (
                 // up squarks
                 pdg_p1 > 0 && is_pdg_up_type_squark(pdg_p1) &&
                 // charginos-
                 pdg_p2 < 0 && is_pdg_chargino(pdg_p2)) ||
             (
                 // anti-up squarks
                 pdg_p2 < 0 && is_pdg_up_type_squark(pdg_p2) &&
                 // charginos+
                 pdg_p1 > 0 && is_pdg_chargino(pdg_p1))) {
      *out1 = internal_p1;
      *out2 = internal_p2;
    }
    // anti-down squark chargino-
    // and anti-up squark chargino+
    else if ((
                 // anti-down squarks
                 pdg_p1 < 0 && is_pdg_down_type_squark(pdg_p1) &&
                 // charginos-
                 pdg_p2 < 0 && is_pdg_chargino(pdg_p2)) ||
             (
                 // down squarks
                 pdg_p2 > 0 && is_pdg_down_type_squark(pdg_p2) &&
                 // charginos+
                 pdg_p1 > 0 && is_pdg_chargino(pdg_p1)) ||
             (
                 // anti up squarks
                 pdg_p1 < 0 && is_pdg_up_type_squark(pdg_p1) &&
                 // charginos+
                 pdg_p2 > 0 && is_pdg_chargino(pdg_p2)) ||
             (
                 // up squarks
                 pdg_p2 > 0 && is_pdg_up_type_squark(pdg_p2) &&
                 // charginos-
                 pdg_p1 < 0 && is_pdg_chargino(pdg_p1))) {
      *out1 = internal_p2;
      *out2 = internal_p1;
    } else {
      fprintf(stderr, "error: Process not recognized (4: %d %d).\n", pdg_p1, pdg_p2);
      return (1);
    }
  } else if ((charge_p1 == 1 || charge_p1 == 0) && (charge_p2 == -1 || charge_p2 == 0)) {
    *out1 = internal_p1;
    *out2 = internal_p2;
  } else if ((charge_p2 == 1 || charge_p2 == 0) && (charge_p1 == -1 || charge_p1 == 0)) {
    *out1 = internal_p2;
    *out2 = internal_p1;
  } else {
    fprintf(stderr, "error: Process not recognized (1: %d %d).\n", pdg_p1, pdg_p2);
    return (1);
  }

  if (not(((charge_p1 == 1 || charge_p1 == 0) && (charge_p2 == -1 || charge_p2 == 0)) ||
          (charge_p2 == 1 || charge_p2 == 0) && (charge_p1 == -1 || charge_p1 == 0))) {
    fprintf(stderr, "error: Process not recognized (2: %d %d).\n", pdg_p1, pdg_p2);
    return (1);
  }

  ///*
  if (not(is_squark_gaugino(*out1, *out2) || is_gaugino_gaugino(*out1, *out2) ||
          is_lepton_lepton(*out1, *out2) || is_gaugino_gluino(*out1, *out2) ||
          is_slepton_slepton(*out1, *out2))) {
    fprintf(stderr, "error: Process not implemented (3: %d %d).\n", pdg_p1, pdg_p2);
    return (1);
  }
  //*/
  // std::cout << *out1 << " " << *out2 << std::endl;
  return 0;
}

void parse_args(int argc, char *argv[], Args &args) {
  Parameters *params = args.params;
  // A value different than "" for an argument means that the user has set
  // that argment through the command-line.
  args.arguments["key"] = "";
  args.arguments["M"] = "";
  args.arguments["pt"] = "";
  args.arguments["pdfset_lo"] = "";
  args.arguments["pdfset_nlo"] = "";
  args.arguments["mu_f"] = "";
  args.arguments["mu_r"] = "";
  args.arguments["particle1"] = "";
  args.arguments["particle2"] = "";
  args.arguments["output"] = "";
  args.arguments["paramfile"] = "";
  args.arguments["center_of_mass_energy"] = "";
  optind = 1;
  for (;;) {
    static struct option long_options[] = {{"version", no_argument, 0, 'v'},
                                           {"help", no_argument, 0, 'h'},
                                           {"parameter-log", required_argument, 0, 'p'},
                                           {"lo", no_argument, 0, 'l'},
                                           {"nlo", no_argument, 0, 'n'},
                                           {"nll", no_argument, 0, 's'},
                                           {"nnll", no_argument, 0, 'z'},
                                           {"invariant-mass", required_argument, 0, 'm'},
                                           {"transverse-momentum", required_argument, 0, 't'},
                                           {"pdfset_lo", required_argument, 0, 'a'},
                                           {"pdfset_nlo", required_argument, 0, 'b'},
                                           {"mu_f", required_argument, 0, 'f'},
                                           {"mu_r", required_argument, 0, 'r'},
                                           {"particle1", required_argument, 0, 'c'},
                                           {"particle2", required_argument, 0, 'd'},
                                           {"output", required_argument, 0, 'o'},
                                           {"slha", required_argument, 0, 'i'},
                                           {"center_of_mass_energy", required_argument, 0, 'e'},
                                           {0, 0, 0, 0}};
    int option_index = 0;
    int c;
    c = getopt_long(argc, argv, "vhp:lnm:t:a:b:f:r:o:", long_options, &option_index);
    if (c == -1) {
      break;
    }
    switch (c) {
    case 'v':
      fprintf(stdout, "resummino " RESUMMINO_VERSION "\n");
      exit(0);
      break;
    case 'h':
      fprintf(stdout, "Please see https://resummino.hepforge.org/ for instructions.");
      exit(0);
      break;
    case 'l':
      args.stop_after_lo = 1;
      break;
    case 'n':
      args.stop_after_nlo = 1;
      break;
    case 's':
      args.nll_unimproved = 1;
      break;
    case 'z':
      args.nnll_unimproved = 1;
      break;
    case 'p':
      args.log_file = optarg;
      break;
    case 'm':
      args.arguments["M"] = optarg;
      break;
    case 't':
      args.arguments["pt"] = optarg;
      break;
    case 'a':
      args.arguments["pdfset_lo"] = optarg;
      break;
    case 'b':
      args.arguments["pdfset_nlo"] = optarg;
      break;
    case 'f':
      args.arguments["mu_f"] = optarg;
      break;
    case 'r':
      args.arguments["mu_r"] = optarg;
      break;
    case 'c':
      args.arguments["particle1"] = optarg;
      break;
    case 'd':
      args.arguments["particle2"] = optarg;
      break;
    case 'o':
      args.arguments["output"] = optarg;
      break;
    case 'i':
      args.arguments["slha"] = optarg;
      break;
    case 'e':
      args.arguments["center_of_mass_energy"] = optarg;
      break;
    case '?':
      break;
    default:
      break;
    }
  }

  // There can be either 0 or 1 left arguments (input file name).
  if (optind > argc) {
    fprintf(stderr, "error: Too many arguments.\n");
    exit(1);
  } else if (optind == argc - 1) {
    args.input_file = argv[optind];
  } else {
    fprintf(stderr, "warning: You did not specify an input file: "
                    "Using default 'resummino.in'.\n");
  }

  args.config = params->read_input_file(args.input_file);

  // Changes working directory to input file folder for relative paths for,
  // e.g., LHAPDF file.
  string old_wd(my_get_current_dir_name());
  string real_name(realpath(args.input_file.c_str(), NULL));
  int count = chdir(real_name.substr(0, real_name.rfind("/")).c_str());

  if (args.config["collider_type"] == "proton-proton") {
    params->ic = 0;
  } else if (args.config["collider_type"] == "proton-antiproton" ||
             args.config["collider_type"] == "antiproton-proton") {
    params->ic = 1;
  } else {
    fprintf(stderr, "error: Collider type not recognized.\n");
    exit(1);
  }

  if (args.config["result"] == "total") {
    args.result_type = total;
  } else if (args.config["result"] == "pt") {
    args.result_type = pt;
  } else if (args.config["result"] == "m") {
    args.result_type = m;
  } else if (args.config["result"] == "ptj") {
    args.result_type = ptj;
  } else {
    fprintf(stderr,
            "error: Computation '%s' not recognized. "
            "'total', 'pt' or 'm' expected.\n",
            args.config["result"].c_str());
    exit(1);
  }

  params->in1 = 0;
  params->in2 = 0;
  int pdg_p1 = atoi(get_option("particle1").c_str());
  int pdg_p2 = atoi(get_option("particle2").c_str());
  params->out1 = pdg_to_internal_id(abs(pdg_p1));
  params->out2 = pdg_to_internal_id(abs(pdg_p2));
  int r = set_particles(pdg_p1, pdg_p2, &params->out1, &params->out2);
  if (r)
    exit(1);
  params->sh = pow2(atof(get_option("center_of_mass_energy").c_str()));
  //    printf("out1 = %i \t out2 = %i \n",params->out1,params->out2);

  // no pt or pt resummation for gaugino gluino production
  if ((is_gaugino_gluino(params->out1, params->out2) ||
       is_squark_gaugino(params->out1, params->out2)) &&
      (args.config["result"] == "ptj" || args.config["result"] == "pt")) {
    fprintf(stdout, "No pt-distribution for associated gaugino-gluino production.\n");
    exit(1);
  }

  // These proc. dont have collinear improvement
  if ((is_gaugino_gluino(params->out1, params->out2) ||
       is_squark_gaugino(params->out1, params->out2))) {
    args.nll_unimproved = 1;
  }

  //   if (not(is_lepton_lepton(params->out1, params->out2) ||
  //   is_slepton_slepton(params->out1, params->out2) ||
  //   is_gaugino_gaugino(params->out1, params->out2)) && nnll_unimproved != 0)
  //   {
  if (not(is_lepton_lepton(params->out1, params->out2) ||
          is_slepton_slepton(params->out1, params->out2) ||
          is_gaugino_gaugino(params->out1, params->out2)) &&
      args.nnll_unimproved != 0) {
    args.nnll_unimproved = 0;
    fprintf(stdout, "NNLL unavailable for the selected process. Resummation "
                    "will be calculated for improved NLL.\n");
  }

  params->read_slha(get_option("slha").c_str());

  // Sets SM and SUSY couplings.
  params->set_couplings();

  //   // Reads and sets Z' and W' parameters and couplings.
  //   if (params->out1 >= 20 && params->out1 <= 25 && params->out2 >= 20 &&
  //   params->out2 <= 25) {
  //       params->read_ZpWp(get_option("zpwp").c_str());
  //   }

  // Reads and sets Z' and W' parameters and couplings.
  if (is_lepton_lepton(params->out1, params->out2)) {
    params->read_ZpWp(get_option("zpwp").c_str());
  }

  //   // Reads Minv_min & Minv_max cuts.
  //   if (params->out1 >= 20 && params->out1 <= 25 && params->out2 >= 20 &&
  //   params->out2 <= 25 && params->out1 == params->out2) {
  //     if (get_option("Minv_min") == "auto") {
  //         params->Minv_min = 3.0 * params->mv[3] / 4.0;
  //     } else {
  //         params->Minv_min = atof(get_option("Minv_min").c_str());
  //     }
  //     // Reads Minv_max cut.
  //     if (get_option("Minv_max") == "auto") {
  //         params->Minv_max = 5.0 * params->mv[3] / 4.0;
  //     } else if (get_option("Minv_max") == "none") {
  //         params->Minv_max =
  //         atof(get_option("center_of_mass_energy").c_str());
  //     } else {
  //         params->Minv_max = atof(get_option("Minv_max").c_str());
  //     }
  //
  //   } else if (params->out1 >= 20 && params->out1 <= 25 && params->out2 >= 20
  //   && params->out2 <= 25 && params->out1 != params->out2) {
  //     if (get_option("Minv_min") == "auto") {
  //         params->Minv_min = 3.0 * params->mv[4] / 4.0;
  //     } else {
  //         params->Minv_min = atof(get_option("Minv_min").c_str());
  //     }
  //     // Reads Minv_max cut.
  //     if (get_option("Minv_max") == "auto") {
  //         params->Minv_max = 5.0 * params->mv[4] / 4.0;
  //     } else if (get_option("Minv_max") == "none") {
  //         params->Minv_max =
  //         atof(get_option("center_of_mass_energy").c_str());
  //     } else {
  //         params->Minv_max = atof(get_option("Minv_max").c_str());
  //     }
  //   }

  // Reads Minv_min & Minv_max cuts.
  if (is_lepton_lepton(params->out1, params->out2) && params->out1 == params->out2) {
    if (get_option("Minv_min") == "auto") {
      params->Minv_min = 3.0 * params->mv[3] / 4.0;
    } else {
      params->Minv_min = atof(get_option("Minv_min").c_str());
    }
    // Reads Minv_max cut.
    if (get_option("Minv_max") == "auto") {
      params->Minv_max = 5.0 * params->mv[3] / 4.0;
    } else if (get_option("Minv_max") == "none") {
      params->Minv_max = atof(get_option("center_of_mass_energy").c_str());
    } else {
      params->Minv_max = atof(get_option("Minv_max").c_str());
    }

  } else if (is_lepton_lepton(params->out1, params->out2) && params->out1 != params->out2) {
    if (get_option("Minv_min") == "auto") {
      params->Minv_min = 3.0 * params->mv[4] / 4.0;
    } else {
      params->Minv_min = atof(get_option("Minv_min").c_str());
    }
    // Reads Minv_max cut.
    if (get_option("Minv_max") == "auto") {
      params->Minv_max = 5.0 * params->mv[4] / 4.0;
    } else if (get_option("Minv_max") == "none") {
      params->Minv_max = atof(get_option("center_of_mass_energy").c_str());
    } else {
      params->Minv_max = atof(get_option("Minv_max").c_str());
    }
  }

  if (get_option("pt") == "auto" || get_option("pt") == "") {
    params->pts = -1;
  } else {
    params->pts = pow2(atof(get_option("pt").c_str()));
  }

  // check if squarks (except stop), are degenerate
  // bool deg_squarks used in breit wigner mapping
  double sq_mass_sum;
  for (int i = 0; i < 12; i++) {
    if (i == 8 || i == 11) {
      continue; // no stops
    }
    sq_mass_sum += params->mSQ[i];
  }
  double aSQm = sq_mass_sum / 10.0; // average squark mass
  params->deg_squarks = true;
  for (int i = 0; i < 12; i++) {
    if (i == 8 || i == 11) {
      continue; // no stops
    }
    if (abs(params->mSQ[i] - aSQm) > 0.1) { // if one of the squark masses is different-> break
      params->deg_squarks = false;
      break;
    }
  }

  //   // pow2(m1 + m2) (sum of final states masses)
  //   double Mfs = 0.0;
  //   if (params->out1 < 10 && params->out2 < 10) {
  //     Mfs = pow2(params->mCH[params->out1] + params->mCH[params->out2]);
  //   } else if (params->out1 >= 10 && params->out1 < 20) {
  //     Mfs = pow2(params->mSL[params->out1 - 10] + params->mSL[params->out2 -
  //     10]);
  //   } else if (params->out1 == 30) {
  //     Mfs = pow2(params->mGL + params->mCH[params->out2]);
  //   } else if (params->out2 == 30) {
  //     Mfs = pow2(params->mGL + params->mCH[params->out1]);
  //   } else if (params->out1 > 30 && params->out2 < 10) {
  //     Mfs = pow2(params->mSQ[params->out1 - 31] + params->mCH[params->out2]);
  //   } else if (params->out2 > 30 && params->out1 < 10) {
  //     Mfs = pow2(params->mSQ[params->out2 - 31] + params->mCH[params->out1]);
  //   } else if (params->out1 >= 20 && params->out1 <= 25 && params->out2 >= 20
  //   && params->out2 <= 25 && params->out1 == params->out2) {
  //     Mfs = pow2(2*(params->mv[3]));
  //   } else if (params->out1 >= 20 && params->out1 <= 25 && params->out2 >= 20
  //   && params->out2 <= 25 && params->out1 != params->out2) {
  //     Mfs = pow2(2*(params->mv[4]));
  //   } else {
  //     Mfs= pow2(params->ml[params->out1 - 20] + params->ml[params->out2 -
  //     20]);
  //   }

  // pow2(m1 + m2) (sum of final states masses)
  double Mfs = 0.0;
  if (is_gaugino_gaugino(params->out1, params->out2)) {
    Mfs = pow2(params->mCH[params->out1] + params->mCH[params->out2]);
  } else if (is_slepton_slepton(params->out1, params->out2)) {
    Mfs = pow2(params->mSL[params->out1 - 10] + params->mSL[params->out2 - 10]);
  } else if (is_gaugino_gluino(params->out1, params->out2) && (params->out1 == 30)) {
    Mfs = pow2(params->mGL + params->mCH[params->out2]);
  } else if (is_gaugino_gluino(params->out1, params->out2) && params->out2 == 30) {
    Mfs = pow2(params->mGL + params->mCH[params->out1]);
  } else if (is_squark_gaugino(params->out1, params->out2) && params->out1 > 30 &&
             params->out2 < 10) {
    Mfs = pow2(params->mSQ[params->out1 - 31] + params->mCH[params->out2]);
  } else if (is_squark_gaugino(params->out1, params->out2) && params->out2 > 30 &&
             params->out1 < 10) {
    Mfs = pow2(params->mSQ[params->out2 - 31] + params->mCH[params->out1]);
  } else if (is_lepton_lepton(params->out1, params->out2) && params->out1 == params->out2) {
    Mfs = pow2(2 * (params->mv[3]));
  } else if (is_lepton_lepton(params->out1, params->out2) && params->out1 != params->out2) {
    Mfs = pow2(2 * (params->mv[4]));
  } else {
    Mfs = pow2(params->ml[params->out1 - 20] + params->ml[params->out2 - 20]);
  }

  if (get_option("M") == "auto" || get_option("M") == "") {
    // Sets automatic "invariant mass".
    params->mis = Mfs;
  } else {
    // Sets user-defined invariant mass. (for inv mass distribution)
    params->mis = pow2(atof(get_option("M").c_str()));
  }

  // factorization and renormalization scale values
  if (args.config["result"] == "m") {
    // Central scale is the invariant mass for invariant mass distributions.
    params->mufs = pow2(atof(get_option("mu_f").c_str())) * params->mis;
    params->murs = pow2(atof(get_option("mu_r").c_str())) * params->mis;
  } else {
    // Central scale is average mass of final state particles (for total XS)
    // fprintf(stderr, "warning: scale not set. Using average of final state
    // masses! \n");
    params->mufs = pow2(atof(get_option("mu_f").c_str())) * params->mis * 0.25;
    params->murs = pow2(atof(get_option("mu_r").c_str())) * params->mis * 0.25;
  }

// activate for fix scale in inv mass distribution
#ifdef FIXED_SCALE
  params->mufs = pow2(atof(get_option("mu_f").c_str())) * Mfs * 0.25;
  params->murs = pow2(atof(get_option("mu_r").c_str())) * Mfs * 0.25;
#endif

  // integration parameters
  params->precision = atof(get_option("precision").c_str());
  params->max_iters = atoi(get_option("max_iters").c_str());
  params->fout = stderr;

  // Restores working directory and saves the parameter log.
  count = chdir(old_wd.c_str());
  if (args.log_file != "") {
    params->write_log(args.log_file.c_str());
  }
}

void init(Args &args) {
  Parameters *params = args.params;
// Sets the PDF at LO.
#if LHAPDF_MAJOR_VERSION == 6
  pdf = LHAPDF::mkPDF(get_option("pdf_lo"), atoi(get_option("pdfset_lo").c_str()));
#else
  if (get_option("pdf_format_lo") == "lhpdf" ||
      (get_option("pdf_format_lo") == "" && get_option("pdf_format") == "lhpdf")) {
    LHAPDF::initPDFSet(get_option("pdf_lo"), LHAPDF::LHPDF, atoi(get_option("pdfset_lo").c_str()));
  } else {
    LHAPDF::initPDFSet(get_option("pdf_lo"), LHAPDF::LHGRID, atoi(get_option("pdfset_lo").c_str()));
  }
#endif
  // Limit PDFs to certain q or g
  if (get_option("pdf_only_g") != "" || get_option("pdf_only_u") != "" ||
      get_option("pdf_only_d") != "" || get_option("pdf_only_s") != "" ||
      get_option("pdf_only_c") != "" || get_option("pdf_only_b") != "" ||
      get_option("pdf_only_ubar") != "" || get_option("pdf_only_dbar") != "" ||
      get_option("pdf_only_sbar") != "" || get_option("pdf_only_cbar") != "" ||
      get_option("pdf_only_bbar") != "") {
    fill(enable_pdfs, enable_pdfs + pdf_members, false);
  } else {
    fill(enable_pdfs, enable_pdfs + pdf_members, true);
  }
  if (get_option("pdf_only_u") != "")
    enable_pdfs[8] = true;
  if (get_option("pdf_only_ubar") != "")
    enable_pdfs[4] = true;
  if (get_option("pdf_only_d") != "")
    enable_pdfs[7] = true;
  if (get_option("pdf_only_dbar") != "")
    enable_pdfs[5] = true;
  if (get_option("pdf_only_s") != "")
    enable_pdfs[9] = true;
  if (get_option("pdf_only_sbar") != "")
    enable_pdfs[3] = true;
  if (get_option("pdf_only_c") != "")
    enable_pdfs[10] = true;
  if (get_option("pdf_only_cbar") != "")
    enable_pdfs[2] = true;
  if (get_option("pdf_only_b") != "")
    enable_pdfs[11] = true;
  if (get_option("pdf_only_bbar") != "")
    enable_pdfs[1] = true;
  if (get_option("pdf_only_g") != "")
    enable_pdfs[6] = true;
  // print strong coupling values
  printf(" alpha_s(M_Z^2) = %6.4f \t for LO PDF set \n", aS(pow2(91.188000), params->set));
  printf(" alpha_s(%4.3f^2) = %6.4f \t for LO PDF set \n", std::sqrt(params->murs),
         aS(params->murs, params->set));

  // inits looptools
  ltini();

  setmudim(params->murs);
  // to generate randon numbers. seed with the time
  srand(time(NULL));

  // reset cached couplings
  first_run = true;
  first_run_cross = true;
}

void done() {
  // Looptools exit.
  ltexi();
}

Results compute(Args &args) {
  Parameters *params = args.params;
  init(args);
  double res0 = 0.0;
  double err0 = 0.0;
  double res1 = 0.0;
  double err1 = 0.0;
  double res2 = 0.0;
  double err2 = 0.0;
  double res3 = 0.0;
  double err3 = 0.0;
  double res4 = 0.0;
  double err4 = 0.0;
  double chisq = 0.0;

  // LO calculation
  // set abs precision to rel precision of lo for nlo and nll computation
  // divided by 3 since each integration routine contributes to final error

  // The multliplicatiopn by 2/M transforms dsigma/dlnMs to dsigma/dM
  switch (args.result_type) {
  case total:
    hadronic_xs(res0, err0, 0, params);
    params->abs_precision = res0 * params->precision / 3.;
    break;
  case pt:
    hadronic_xs_dPT2(res0, err0, chisq, 0, -0, params);
    params->abs_precision = res0 * params->precision / 3.;
    res0 *= 2.0 * sqrt(params->pts);
    err0 *= 2.0 * sqrt(params->pts);
    break;
  case ptj:
    hadronic_xs_dPT2(res0, err0, chisq, 0, -0, params);
    params->abs_precision = res0 * params->precision / 3.;
    res0 *= 2.0 * sqrt(params->pts);
    err0 *= 2.0 * sqrt(params->pts);
    break;
  case m:
    hadronic_xs_dlnM2(res0, err0, chisq, 0, -0, params);
    params->abs_precision = res0 * params->precision / 3.;
    res0 *= 2.0 / sqrt(params->mis);
    err0 *= 2.0 / sqrt(params->mis);
    break;
  }
  res1 = res0;
  err1 = pow2(err0);

  printf("LO = (%.7e +- %.7e) %s\n", res1, sqrt(err1), UNITS);

  if (!args.stop_after_lo) {
    // Sets the PDF at NLO.
#if LHAPDF_MAJOR_VERSION == 6
    pdf = LHAPDF::mkPDF(get_option("pdf_nlo"), atoi(get_option("pdfset_nlo").c_str()));
#else
    if (get_option("pdf_format_nlo") == "lhpdf" ||
        (get_option("pdf_format_nlo") == "" && get_option("pdf_format") == "lhpdf")) {
      LHAPDF::initPDFSet(get_option("pdf_nlo"), LHAPDF::LHPDF,
                         atoi(get_option("pdfset_nlo").c_str()));
    } else {
      LHAPDF::initPDFSet(get_option("pdf_nlo"), LHAPDF::LHGRID,
                         atoi(get_option("pdfset_nlo").c_str()));
    }
#endif

    // Print values for alpha
    printf(" alpha_s(M_Z^2) = %6.4f \t for NLO PDF set \n", aS(pow2(91.188000), params->set));
    printf(" alpha_s(%4.3f^2) = %6.4f \t for NLO PDF set \n", std::sqrt(params->murs),
           aS(params->murs, params->set));

    // NLO
    int start = 0;
    for (int i0 = start; i0 < 5; i0++) {
      switch (args.result_type) {
      case total:
        hadronic_xs(res0, err0, i0, params);
        break;
      case pt:
        hadronic_xs_dPT2(res0, err0, chisq, i0, -0, params);
        res0 *= 2.0 * sqrt(params->pts);
        err0 *= 2.0 * sqrt(params->pts);
        break;
      case ptj:
        hadronic_xs_dPT2(res0, err0, chisq, i0, -0, params);
        res0 *= 2.0 * sqrt(params->pts);
        err0 *= 2.0 * sqrt(params->pts);
        break;
      case m:
        hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, params);
        res0 *= 2.0 / sqrt(params->mis);
        err0 *= 2.0 / sqrt(params->mis);
        break;
      }
      if (i0 == start) {
        res2 = res0;
        err2 = pow2(err0);
      } else if (i0 < 5) {
        res2 += res0;
        err2 += pow2(err0);
      }
      if (i0 == 0) {
        printf("LO = (%.7e +- %.7e) %s\n", res0, (err0), UNITS);
      } else if (i0 == 1) {
        printf("vNLO = (%.7e +- %.7e) %s\n", res0, (err0), UNITS);
      } else if (i0 == 2) {
        printf("P+K = (%.7e +- %.7e) %s\n", res0, (err0), UNITS);
      } else if (i0 == 3) {
        printf("rNLOg = (%.7e +- %.7e) %s\n", res0, (err0), UNITS);
      } else if (i0 == 4) {
        printf("rNLOq = (%.7e +- %.7e) %s\n", res0, (err0), UNITS);
      }
    }
    printf("NLO = (%.7e +- %.7e) %s\n", res2, sqrt(err2), UNITS);

    if (!args.stop_after_lo && !args.stop_after_nlo) {

      // NLO+NLL

      // default weights for PDF fitting
      double weight_valence = -1.6;
      double weight_sea = -1.6;
      double weight_gluon = -1.6;
      double xmin = params->mis / params->sh;

      // reads fitting parameters from input file if specified
      if (get_option("xmin") == "" || get_option("xmin") == "auto") {
        xmin = params->mis / params->sh;
      } else {
        xmin = atof(get_option("xmin").c_str());
      }
      if (get_option("weight_valence") == "" || get_option("weight_valence") == "auto") {
        weight_valence = -1.6;
      } else {
        weight_valence = atof(get_option("weight_valence").c_str());
      }
      if (get_option("weight_sea") == "" || get_option("weight_sea") == "auto") {
        weight_sea = -1.6;
      } else {
        weight_sea = atof(get_option("weight_sea").c_str());
      }
      if (get_option("weight_gluon") == "" || get_option("weight_gluon") == "auto") {
        weight_gluon = -1.6;
      } else {
        weight_gluon = atof(get_option("weight_gluon").c_str());
      }

      // PDF fitting procedure
      pdfFit(params->a1min, params->afit, params->mis / params->sh, params->mufs, weight_valence,
             weight_sea, weight_gluon, xmin);
      int i0 = 5;
      switch (args.result_type) {
      case total:
        if ((args.nnll_unimproved == 0) && (args.nll_unimproved == 0)) {
          hadronic_xs(res0, err0, i0, params);
        } else if ((args.nnll_unimproved == 0) && (args.nll_unimproved != 0)) {
          i0 = 6;
          hadronic_xs(res0, err0, i0, params);
        } else if (args.nnll_unimproved != 0) {
          i0 = 7;
          hadronic_xs(res0, err0, i0, params);
        }
        break;
      case pt:
        hadronic_xs_dPT2(res0, err0, chisq, i0, -0, params);
        res0 *= 2.0 * sqrt(params->pts);
        err0 *= 2.0 * sqrt(params->pts);
        break;
      case ptj:
        hadronic_xs_dPT2(res0, err0, chisq, i0, -0, params);
        res0 *= 2.0 * sqrt(params->pts);
        err0 *= 2.0 * sqrt(params->pts);
        break;
      case m:
        if ((args.nnll_unimproved == 0) && (args.nll_unimproved == 0)) {
          hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, params);
        } else if ((args.nnll_unimproved == 0) && (args.nll_unimproved != 0)) {
          i0 = 6;
          hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, params);
        } else if (args.nnll_unimproved != 0) {
          i0 = 7;
          hadronic_xs_dlnM2(res0, err0, chisq, i0, -0, params);
        }
        res0 *= 2.0 / sqrt(params->mis);
        err0 *= 2.0 / sqrt(params->mis);
        break;
      }
#ifndef EXPANSION
      res3 = res2 + res0;
      err3 = err2 + pow2(err0);
#endif

#ifdef EXPANSION
      res3 = res0;
      err3 = pow2(err0);
#endif

      if (args.nnll_unimproved != 0) {
        printf("aNNLO+NNLL = (%.7e +- %.7e) %s\n", res3, sqrt(err3), UNITS);
      } else {
        printf("NLO+NLL = (%.7e +- %.7e) %s\n", res3, sqrt(err3), UNITS);
      }
      // Joint resummation.
      if (args.result_type == ptj) {
        hadronic_xs_dPT2(res0, err0, chisq, 6, -0, params);
        res0 *= 2.0 * sqrt(params->pts);
        err0 *= 2.0 * sqrt(params->pts);
        res4 = res2 + res0;
        err4 = err2 + pow2(err0);
        printf("NLO+NLLj = (%.7e +- %.7e) %s\n", res4, sqrt(err4), UNITS);
      }
    }
  }
  done();
  return Results{{res1, sqrt(err1)}, {res2, sqrt(err2)}, {res3, sqrt(err3)}, {res4, sqrt(err4)}};
}

void print_results(Results &r, Args &args) {
  Parameters *params = args.params;
  double res1 = r.LO.res;
  double err1 = r.LO.err;
  double res2 = r.NLO.res;
  double err2 = r.NLO.err;
  double res3 = r.NLL.res;
  double err3 = r.NLL.err;
  double res4 = r.aNLLj.res;
  double err4 = r.aNLLj.err;
  // Prints again the results.
  if (args.result_type == ptj) {
    printf("\nResults:\n"
           "LO = (%.7e +- %.7e) %s\n"
           "NLO = (%.7e +- %.7e) %s\n"
           "NLO+NLL = (%.7e +- %.7e) %s\n"
           "NLO+NLLj = (%.7e +- %.7e) %s\n",
           res1, (err1), UNITS, res2, (err2), UNITS, res3, (err3), UNITS, res4, (err4), UNITS);
  } else if (args.nnll_unimproved != 0) {
    printf("\nResults:\n"
           "LO = (%.7e +- %.7e) %s\n"
           "NLO = (%.7e +- %.7e) %s\n"
           "aNNLO+NNLL = (%.7e +- %.7e) %s\n",
           res1, (err1), UNITS, res2, (err2), UNITS, res3, (err3), UNITS);
  } else {
    printf("\nResults:\n"
           "LO = (%.7e +- %.7e) %s\n"
           "NLO = (%.7e +- %.7e) %s\n"
           "NLO+NLL = (%.7e +- %.7e) %s\n",
           res1, (err1), UNITS, res2, (err2), UNITS, res3, (err3), UNITS);
  }
  // Prints results to a file in JSON format for a subsequent analysis.
  if (get_option("output") != "") {
    FILE *fout = fopen(get_option("output").c_str(), "w");
    if (!fout) {
      fprintf(stderr, "error: Could not open output file %s.\n", get_option("output").c_str());
      // return 1;
      exit(1);
    }
    fprintf(fout,
            "{\n"
            "\"key\": \"%s\",\n"
            "\"pt\": %s,\n"
            "\"m\": %s,\n"
            "\"pdflo\": \"%s\",\n"
            "\"pdfsetlo\": %s,\n"
            "\"pdfnlo\": \"%s\",\n"
            "\"pdfsetnlo\": %s,\n"
            "\"muf\": %s,\n"
            "\"mur\": %s,\n"
            "\"lo\": %.7e,\n"
            "\"nlo\": %.7e,\n"
            "\"nll\": %.7e,\n"
            "\"nllj\": %.7e,\n"
            "\"units\": \"%s\"\n"
            "}",
            get_option("key").c_str(),
            (get_option("pt") == "auto" ? "-1" : get_option("pt").c_str()),
            (get_option("M") == "auto" ? "-1" : get_option("M").c_str()),
            get_option("pdf_lo").c_str(), get_option("pdfset_lo").c_str(),
            get_option("pdf_nlo").c_str(), get_option("pdfset_nlo").c_str(),
            get_option("mu_f").c_str(), get_option("mu_r").c_str(), res1, res2, res3, res4, UNITS);
  }
}
