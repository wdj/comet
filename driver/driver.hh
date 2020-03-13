//-----------------------------------------------------------------------------
/*!
 * \file   driver.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_driver_hh_
#define _comet_driver_hh_

#include "env.hh"
#include "checksum.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Struct to hold driver options (options not in CEnv)---*/

typedef struct {
  int num_field_local;
  int num_vector_local;
  size_t num_field;
  size_t num_vector;
  size_t num_field_active;
  size_t num_vector_active;
  bool num_field_local_initialized;
  bool num_field_active_initialized;
  bool num_vector_local_initialized;
  bool num_vector_active_initialized;
  int verbosity;
  int stage_min_0based;
  int stage_max_0based;
  int phase_min_0based;
  int phase_max_0based;
  char* input_file_path;
  char* metrics_file_path_stub;
  int problem_type;
  size_t num_incorrect;
  double max_incorrect_diff;
  //double threshold;
  bool checksum;
} DriverOptions;

enum {
  GM_PROBLEM_TYPE_RANDOM = 1,
  GM_PROBLEM_TYPE_ANALYTIC = 2
};

//=============================================================================

//void finish_parsing(int argc, char** argv, DriverOptions* do_, CEnv* env);

void perform_run(int argc, char** argv, const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(comet::Checksum& cksum, int argc, char** argv,
                 const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(comet::Checksum& cksum, const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);


void print_output(bool do_print,
                  Checksum& cksum,
                  CEnv& env,
                  char* metrics_file_path_stub = 0,
                  size_t num_written = 0,
                  double vctime = 0,
                  double mctime = 0,
                  double cktime = 0,
                  double intime = 0,
                  double outtime = 0,
                  double tottime = 0);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_driver_hh_

//-----------------------------------------------------------------------------
