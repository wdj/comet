//-----------------------------------------------------------------------------
/*!
 * \file   driver.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions, header.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef _COMET_DRIVER_HH_
#define _COMET_DRIVER_HH_

#include "env.hh"
#include "checksum.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Struct to hold driver options (options not in CEnv).

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
  int stage_min;
  int stage_max;
  int phase_min;
  int phase_max;
  char* input_file;
  char* histograms_file;
  char* output_file_stub;
  int problem_type;
  size_t num_incorrect;
  double max_incorrect_diff;
  bool checksum;
} DriverOptions;

enum {
  GM_PROBLEM_TYPE_RANDOM = 1,
  GM_PROBLEM_TYPE_ANALYTIC = 2
};

//=============================================================================

//void finish_parsing(int argc, char** argv, DriverOptions* do_, CEnv* env);

void print_output(bool do_print,
                  Checksum& cksum,
                  CEnv& env,
                  char* output_file_stub = 0,
                  size_t num_written = 0,
                  double vctime = 0,
                  double mctime = 0,
                  double cktime = 0,
                  double intime = 0,
                  double outtime = 0,
                  double tottime = 0);

void perform_run(int argc, char** argv, const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(comet::Checksum& cksum, int argc, char** argv,
                 const char* const description,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

void perform_run(comet::Checksum& cksum, const char* const options,
                 MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_DRIVER_HH_

//-----------------------------------------------------------------------------
