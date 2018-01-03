//-----------------------------------------------------------------------------
/*!
 * \file   driver.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_driver_hh_
#define _gm_driver_hh_

#include "env.hh"
#include "checksums.hh"

//=============================================================================
/*---Struct to hold driver options (options not in GMEnv)---*/

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
  char* output_file_path_stub;
  int problem_type;
  size_t num_incorrect;
  double threshold;
  bool checksum;
} DriverOptions;

enum {
  GM_PROBLEM_TYPE_RANDOM = 1,
  GM_PROBLEM_TYPE_ANALYTIC = 2
};

//=============================================================================

//void finish_parsing(int argc, char** argv, DriverOptions* do_, GMEnv* env);

GMChecksum perform_run(int argc, char** argv, const char* const description);

GMChecksum perform_run(const char* const options);

//=============================================================================

#endif /*---_gm_driver_hh_---*/

//-----------------------------------------------------------------------------
