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
#include "vectors.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Struct to hold driver options (options not in CEnv).

class Driver {

  struct Options {
    int num_field_local;
    int num_vector_local;
    size_t num_field_active;
    size_t num_vector_active;
    bool is_inited_num_field_local;
    bool is_inited_num_field_active;
    bool is_inited_num_vector_local;
    bool is_inited_num_vector_active;
    int verbosity;
    int stage_min;
    int stage_max;
    int phase_min;
    int phase_max;
    char* input_file;
    char* histograms_file;
    char* output_file_stub;
    int problem_type;
    bool checksum;

    Options(CEnv& env);
  };

  struct Counters {
    size_t num_incorrect;
    double max_incorrect_diff;
    size_t num_metric_items_computed;
    size_t num_metrics_active;
    size_t num_written;
    double vctime;
    double mctime;
    double cktime;
    double intime;
    double outtime;
    double tottime;
    size_t num_metric_items_local_computed;
    size_t num_metrics_active_local;
    size_t num_local_written;

    Counters();
  };

  Options options_;
  Counters counters_;

public:

  Driver(CEnv& env);

  void finish_parsing(int argc, char** argv);

  void set_vectors(GMVectors& vectors);

  void print_output_sync(Checksum& cksum);

  static void print_output_sync(Checksum& cksum, CEnv& env) {
    Driver driver(env);
    driver.print_output_sync(cksum);
  }

  static void perform_run(int argc, char** argv, MPI_Comm base_comm);

  static void perform_run(const char* const options);

  static void perform_run(const char* const options, MPI_Comm base_comm,
                          CEnv& env);

  static void perform_run(Checksum& cksum, const char* const description,
                          MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

private:

  friend class TestProblem;

  CEnv& env_;

  bool do_print() const {
    return env_.is_proc_active() && env_.proc_num() == 0 && options_.verbosity > 0;
  }

  static void perform_run_(Checksum& cksum, int argc, char** argv,
                           const char* const description,
                           MPI_Comm base_comm = MPI_COMM_WORLD, CEnv* env = 0);

  static void perform_run_(Checksum& cksum, int argc, char** argv,
                           MPI_Comm base_comm, CEnv& env);

  // Helper function to flush output on all ranks.

  void fflush_sync_() {
    env_.synced_time();
    fflush(NULL);
    env_.synced_time();
  }

  // Local class to manage timings.

  class Timer {
  public:
    Timer(CEnv& env) : env_(env), time_begin_(0) {
      start();
    }
    void start() {
      time_begin_ = env_.synced_time();
    }
    double elapsed() const {
      return env_.synced_time() - time_begin_;
    }
    void add_elapsed(double timeval) const {
      timeval += env_.synced_time() - time_begin_;
    }

  private:
    CEnv& env_;
    double time_begin_;
  };
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_DRIVER_HH_

//-----------------------------------------------------------------------------
