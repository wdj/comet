//-----------------------------------------------------------------------------
/*!
 * \file   metrics_io.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O utilities for metrics.
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

#include "cstdio"
#include "cstdint"
#include "string.h"
//#include "errno.h"
#include "typeinfo"

#include "env.hh"
#include "metrics.hh"
#include "decomp_mgr.hh"

//#include "driver.hh"
#include "metrics_io.hh"

//=============================================================================

namespace comet {

//=============================================================================
// SingleMetricIO member definitions.

//-----------------------------------------------------------------------------
/*!
 * \brief Constructor for SingleMetricIO class
 *
 */
SingleMetricIO::SingleMetricIO(FILE* file, GMMetrics& metrics, CEnv& env)
  : env_(env)
  , file_(file)
  , data_type_(env.data_type_metrics())
  , num_way_(env.num_way())
  , num_written_(0) {

  if (stdout == file_)
    return;

  // Check for adequate size for indexes.

  typedef BasicTypes::BigUInt IndexMax_t;
  const auto index_max = static_cast<IndexMax_t>(metrics.num_vector_active) *
    (env_.is_metric_type_bitwise() ? 2 : 1) - 1;
  COMET_INSIST_INTERFACE(&env, index_max ==
     static_cast<IndexMax_t>(static_cast<IntForFile_t>(index_max)) &&
                "Too many vectors for output format.");
  // Check this also just to be sure - this check may be unnecessary.
  COMET_INSIST_INTERFACE(&env, index_max ==
     static_cast<IndexMax_t>(static_cast<NV_t>(index_max)) &&
                "Too many vectors for output format.");
}

//-----------------------------------------------------------------------------
/*!
 * \brief Write a metric value: CZEK 2-way case; helper for other cases.
 *
 */
template<typename Float_t>
void SingleMetricIO::write(NV_t iG, NV_t jG, Float_t value) const {

  bool success = true;
  size_t bytes_written = 0;

  // WARNING: possibly precision-reducing cast:
  const auto value_for_file = static_cast<FloatForFile_t>(value);

  const auto iG_for_file = safe_cast_assert<IntForFile_t>(iG);
  const auto jG_for_file = safe_cast_assert<IntForFile_t>(jG);

  // Do the writes to the file.
  success = success && write_<IntForFile_t>(iG_for_file, bytes_written);
  success = success && write_<IntForFile_t>(jG_for_file, bytes_written);
  success = success && write_<FloatForFile_t>(value_for_file, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
  }

//-----------------------------------------------------------------------------
/*!
 * \brief Write a metric value: CZEK 3-way case; helper for other cases.
 *
 */
template<typename Float_t>
void SingleMetricIO::write(NV_t iG, NV_t jG, NV_t kG, Float_t value) const {

  bool success = true;
  size_t bytes_written = 0;

  // WARNING: possibly precision-reducing cast:
  const auto value_for_file = static_cast<FloatForFile_t>(value);

  const auto iG_for_file = safe_cast_assert<IntForFile_t>(iG);
  const auto jG_for_file = safe_cast_assert<IntForFile_t>(jG);
  const auto kG_for_file = safe_cast_assert<IntForFile_t>(kG);

  // Do the writes to the file.
  success = success && write_<IntForFile_t>(iG_for_file, bytes_written);
  success = success && write_<IntForFile_t>(jG_for_file, bytes_written);
  success = success && write_<IntForFile_t>(kG_for_file, bytes_written);
  success = success && write_<FloatForFile_t>(value_for_file, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Write a metric value: CCC/DUO 2-way case.
 *
 */
template<typename Float_t>
void SingleMetricIO::write(NV_t iG, NV_t jG, int iE, int jE,
                           Float_t value) const {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);

  // Use czek write function as helper.
  write(iE + 2 * iG, jE + 2 * jG, value);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Write a metric value: CCC/DUO 3-way case.
 *
 */
template<typename Float_t>
void SingleMetricIO::write(NV_t iG, NV_t jG, NV_t kG,
                           int iE, int jE, int kE, Float_t value) const {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);

  // Use czek write function as helper.
  write(iE + 2 * iG, jE + 2 * jG, kE + 2 * kG, value);
}

//=============================================================================
// Helper functions for MetricsIO class.

//-----------------------------------------------------------------------------
/*!
 * \brief Write to file or print the metrics.
 *
 */
template<typename Float_t>
static void MetricsIO_write_impl_(GMMetrics& metrics, FILE* file,
                                  NML_t& num_written_, CEnv& env) {
  COMET_INSIST(file);
  COMET_INSIST(env.is_double_prec() == (typeid(Float_t) == typeid(double)));
  COMET_INSIST(!env.is_double_prec() == (typeid(Float_t) == typeid(float)));

  if (! env.is_proc_active())
    return;

  // Due to field-redundancy, only results from some processors are needed.
  if (env.proc_num_field() != 0)
    return;

  const bool is_writing_to_stdout = stdout == file;

  typedef double FloatForStdout_t;

  //--------------------
  if (env.data_type_metrics() == DataTypeId::FLOAT &&
      env.num_way() == NumWay::_2) {
  //--------------------

    SingleMetricIO writer(file, metrics, env);

    COMET_INSIST(!env.thresholds().is_multi());

    // Loop over metric item candidates for write.
    for (NML_t index = 0; index < metrics.num_metrics_local; ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;

      // Directly access the metrics values array.
      const auto value = Metrics_elt_const<Float_t>(metrics, index, env);

      // Do not write if doesn't pass threshold.
      if (!env.thresholds().is_pass(value))
        continue;

      // Output the value.
      if (is_writing_to_stdout)
        fprintf(file, "element (%li,%li): value: %.17e\n",
          iG, jG, static_cast<FloatForStdout_t>(value));
      else
        writer.write(iG, jG, value);

    } // for index

    num_written_ += writer.num_written();

  //--------------------
  } else if (env.data_type_metrics() == DataTypeId::FLOAT) {
          // && env.num_way() == NumWay::_3
  //--------------------

    SingleMetricIO writer(file, metrics, env);

    COMET_INSIST(!env.thresholds().is_multi());

    // Loop over metric item candidates for write.
    for (NML_t index = 0; index < metrics.num_metrics_local; ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env);
      const NV_t kG = Metrics_coords_getG(metrics, index, 2, env);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;

      // Directly access the metrics values array.
      const auto value = Metrics_elt_const<Float_t>(metrics, index, env);

      // Do not write if doesn't pass threshold.
      if (!env.thresholds().is_pass(value))
        continue;

      // Output the value.
      if (is_writing_to_stdout)
        fprintf(file, "element (%li,%li,%li): value: %.17e\n",
          iG, jG, kG, static_cast<FloatForStdout_t>(value));
      else
        writer.write(iG, jG, kG, value);

    } // for index

    num_written_ += writer.num_written();

  //--------------------
  } else if (env.data_type_metrics() == DataTypeId::TALLY2X2) {
  //--------------------

    COMET_INSIST(env.num_way() == NumWay::_2);
    COMET_INSIST(env.is_metric_type_bitwise());

    SingleMetricIO writer(file, metrics, env);

    // Loop over metric item candidates for write.
    for (NML_t index = 0; index < metrics.num_metric_items_local_computed;
         ++index) {
      const MetricItemCoords_t coords = metrics.coords_value(index);
      const NV_t iG = CoordsInfo::getiG(coords, metrics, env);
      const NV_t jG = CoordsInfo::getjG(coords, metrics, env);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;

      // Loop over metric entries in this item.
      for (int entry_num = 0; entry_num < env.num_entries_per_metric_item();
           ++entry_num) {
        const int iE = CoordsInfo::getiE(coords, entry_num, metrics, env);
        const int jE = CoordsInfo::getjE(coords, entry_num, metrics, env);

        // Extract or calculate the metric value.
        const auto value = Metrics_ccc_duo_get_2<Float_t>(metrics, index,
          entry_num, env);

        // For non-shrink case, metric not thresholded yet, so threshold.
        if (!Metrics_is_pass_threshold_noshrink(metrics, index, iE, jE, env))
          continue;

        // Output the value.
        if (is_writing_to_stdout)
          fprintf(file, "element (%li,%li) (%i,%i): value: %.17e\n",
            iG, jG, iE, jE, static_cast<FloatForStdout_t>(value));
        else
          writer.write(iG, jG, iE, jE, value);

      } // for entry_num
    } // for index

    num_written_ += writer.num_written();

  //--------------------
  } else { // if (env.data_type_metrics() == DataTypeId::TALLY4X2)
  //--------------------

    COMET_INSIST(env.data_type_metrics() == DataTypeId::TALLY4X2);
    COMET_INSIST(env.num_way() == NumWay::_3);
    COMET_INSIST(env.is_metric_type_bitwise());

    SingleMetricIO writer(file, metrics, env);

    // Loop over metric item candidates for write.
    for (NML_t index = 0; index < metrics.num_metric_items_local_computed;
         ++index) {
      const MetricItemCoords_t coords = metrics.coords_value(index);
      const NV_t iG = CoordsInfo::getiG(coords, metrics, env);
      const NV_t jG = CoordsInfo::getjG(coords, metrics, env);
      const NV_t kG = CoordsInfo::getkG(coords, metrics, env);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;

      // Loop over metric entries in this item.
      for (int entry_num = 0; entry_num < env.num_entries_per_metric_item();
           ++entry_num) {
        const int iE = CoordsInfo::getiE(coords, entry_num, metrics, env);
        const int jE = CoordsInfo::getjE(coords, entry_num, metrics, env);
        const int kE = CoordsInfo::getkE(coords, entry_num, metrics, env);

        // Extract or calculate the metric value.
        const auto value = Metrics_ccc_duo_get_3<Float_t>(metrics, index,
          entry_num, env);

        // For non-shrink case, metric not thresholded yet, so threshold.
        if (!Metrics_is_pass_threshold_noshrink(metrics, index, iE, jE, kE, env))
          continue;

        // Output the value.
        if (is_writing_to_stdout)
          fprintf(file, "element (%li,%li,%li) (%i,%i,%i): value: %.17e\n",
            iG, jG, kG, iE, jE, kE, static_cast<FloatForStdout_t>(value));
        else
          writer.write(iG, jG, kG, iE, jE, kE, value);

      } // for entry_num
    } // for index

    num_written_ += writer.num_written();

  //----------
  } // if
  //----------
}

//-----------------------------------------------------------------------------
/*!
 * \brief Write to file or print the metrics.
 *
 */
static void MetricsIO_write_(GMMetrics& metrics, FILE* file,
                             NML_t& num_written_, CEnv& env) {
  if (!env.is_double_prec())
    MetricsIO_write_impl_<float>(metrics, file, num_written_, env);
  else
    MetricsIO_write_impl_<double>(metrics, file, num_written_, env);
}

//=============================================================================
// MetricsIO member definitions.

//-----------------------------------------------------------------------------
/*!
 * \brief Constructor for MetricsIO class.
 *
 */
MetricsIO::MetricsIO(const char* path_stub, int verbosity, CEnv& env)
  : env_(env)
  , path_stub_str_(path_stub ? path_stub : "")
  , file_(NULL)
  , verbosity_(verbosity)
  , num_written_(0)
  , num_written_last_write_(0)
  , is_allocated_(false) {

  if (! env_.is_proc_active())
    return;

  if (is_path_stub_() && is_leaving_files_open_()) {
    open_();
    COMET_INSIST(NULL != file_ && "Unable to open file.");
  }

  is_allocated_ = true;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Destructor for MetricsIO class.
 *
 */
MetricsIO::~MetricsIO() {
  deallocate();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Termination function called by destructor.
 *
 */
void MetricsIO::deallocate() {

  if (!is_allocated_)
    return;

  if (is_path_stub_() && is_leaving_files_open_())
    close_();

  is_allocated_ = false;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Write elements of Metrics object to file; if requested also print.
 *
 */
void MetricsIO::write(GMMetrics& metrics) {

  if (! env_.is_proc_active())
    return;

  COMET_INSIST(is_allocated_);

  // Due to field-redundancy, only results from some processors are needed.
  if (env_.proc_num_field() != 0)
    return;

  // Output to stdout if requested

  if (verbosity_ > 1)
    MetricsIO_write_(metrics, stdout, num_written_, env_);

  if (!is_path_stub_())
    return;

  // Output to file

  if (!is_leaving_files_open_())
    open_();

  const size_t num_written_before_write = num_written_;
  MetricsIO_write_(metrics, file_, num_written_, env_);
  num_written_last_write_ = num_written_ - num_written_before_write;

  if (!is_leaving_files_open_())
    close_();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Check the metrics elements most recently written to the file.
 *
 */
template<typename Float_t>
void MetricsIO::check_file_impl_(GMMetrics& metrics) {

  if (! env_.is_proc_active())
    return;

  COMET_INSIST(is_allocated_);

  // Due to field-redundancy, only results from some processors are needed.
  if (env_.proc_num_field() != 0)
    return;

  if (!is_path_stub_())
    return;

  // Close file and reopen file for read.

  if (is_leaving_files_open_())
    close_();
  open_("rb"); 

  // Set fpos to beginning of last metrics object write.

  long offset = is_leaving_files_open_() ? bytes_(num_written_ - num_written_last_write_) : 0;
  const int success1 = fseek(file_, offset, SEEK_SET);
  COMET_INSIST(0 == success1 && "Unable to fseek to correct file position.");

  //----------------------------------------
  // Test 1: loop over metrics stored in file, compare to metrics in memory.
  //----------------------------------------

  size_t num_incorrect = 0;

  for (NML_t index_f = 0; index_f < num_written_last_write_; ++index_f) {

    //--------------------
    if (env_.num_way() == NumWay::_2) {
    //--------------------

      // Read current metric from file.

      MetricForFile<NumWay::_2> metric_f;
      metric_f.read(file_, env_);

      const auto iG_f = metric_f.iG(env_);
      const auto jG_f = metric_f.jG(env_);
      const auto iE_f = metric_f.iE(env_);
      const auto jE_f = metric_f.jE(env_);

      // Access current metric from memory.

      const auto metric_value = env_.is_shrink() ?
        // Extract lone metric entry for index.
        Metrics_ccc_duo_get_2<Float_t>(metrics, index_f, 0, env_) :
        // Extract or compute single table entry from metric.
        GMMetrics_get_2<Float_t>(metrics, iG_f, jG_f, iE_f, jE_f, env_);

      // Infer index represented by coords stored in file.
      // NOTE: throws assertion if these global coords not on this rank.
      // NOTE: may differ from index_f because the latter skips
      // thresholded values and this may not.
      const NML_t index_m = Metrics_index_2(metrics, iG_f, jG_f, env_);

      // Compare indexing.

      bool do_coords_match = true;

      if (!env_.is_shrink()) {
        // NOTE: the following is only valid for non-shrink case.
        MetricItemCoords_t coords = metrics.coords_value(index_m);
        do_coords_match =
          CoordsInfo::getiG(coords, metrics, env_) == iG_f &&
          CoordsInfo::getjG(coords, metrics, env_) == jG_f;
      }

      // NOTE: for is_shrink case, some metric entries may be
      // unavailable to check thresholding (specif., if is_multi),
      // so this call will not entirely check this case.

      const bool pass_threshold =
        Metrics_is_pass_threshold_noshrink<Float_t>(metrics, index_m, iE_f, jE_f, env_);

      const bool is_correct = do_coords_match && pass_threshold && 
        static_cast<FloatForFile_t>(metric_value) == metric_f.value();

      num_incorrect += !is_correct;

      // Output some representative failures.
      if (num_incorrect < 10 && !is_correct)
        fprintf(stderr, "Incorrect metric value in file: element %zu %zu "
          " %i %i actual %.17e expected %.17e pass_threshold %i coords_match %i\n",
          iG_f, jG_f, iE_f, jE_f,
          static_cast<double>(metric_f.value()),
          static_cast<double>(metric_value),
          pass_threshold, do_coords_match);

    //--------------------
    } else { // if (env_.num_way() == NumWay::_3)
    //--------------------

      // Read current metric from file.

      MetricForFile<NumWay::_3> metric_f;
      metric_f.read(file_, env_);

      const auto iG_f = metric_f.iG(env_);
      const auto jG_f = metric_f.jG(env_);
      const auto kG_f = metric_f.kG(env_);
      const auto iE_f = metric_f.iE(env_);
      const auto jE_f = metric_f.jE(env_);
      const auto kE_f = metric_f.kE(env_);

      // Access current metric from memory.

      const auto metric_value = env_.is_shrink() ?
        // Extract lone metric entry for index.
        Metrics_ccc_duo_get_3<Float_t>(metrics, index_f, 0, env_) :
        GMMetrics_get_3<Float_t>(metrics, iG_f, jG_f, kG_f, iE_f, jE_f, kE_f, env_);

      // Infer index represented by coords stored in file.
      // NOTE: throws assertion if these global coords not on this rank.
      // NOTE: may differ from index_f because the latter skips
      // thresholded values and this may not.
      const NML_t index_m = Metrics_index_3(metrics, iG_f, jG_f, kG_f, env_);

      // Compare indexing.

      bool do_coords_match = true;

      if (!env_.is_shrink()) {
        // NOTE: the following is only valid for non-shrink case.
        MetricItemCoords_t coords = metrics.coords_value(index_m);
        do_coords_match =
          CoordsInfo::getiG(coords, metrics, env_) == iG_f &&
          CoordsInfo::getjG(coords, metrics, env_) == jG_f &&
          CoordsInfo::getkG(coords, metrics, env_) == kG_f;
      }

      const bool pass_threshold =
        Metrics_is_pass_threshold_noshrink<Float_t>(metrics, index_m, iE_f, jE_f, kE_f, env_);

      const bool is_correct = do_coords_match && pass_threshold && 
        static_cast<FloatForFile_t>(metric_value) == metric_f.value();

      num_incorrect += !is_correct;

      if (!is_correct && num_incorrect < 10)
        fprintf(stderr, "Incorrect metric value in file: " "element %zu %zu %zu "
          " %i %i %i actual %.17e expected %.17e pass_threshold %i coords_match %i\n",
          iG_f, jG_f, kG_f, iE_f, jE_f, kE_f,
          static_cast<double>(metric_f.value()),
          static_cast<double>(metric_value),
          pass_threshold, do_coords_match);

    //--------------------
    } // if (env_.num_way() == NumWay::_2)
    //--------------------

  } // index_f

  COMET_INSIST(0 == num_incorrect &&
    "Incorrect values detected in output file - partial list is displayed.");

  //----------------------------------------
  // Test 2: count metrics (in memory) passing threshold, compare to num written.
  //----------------------------------------

  NML_t num_passed = 0;

  //--------------------
  if (env_.num_way() == NumWay::_2) {
  //--------------------

    for (NML_t index = 0; index <  metrics.num_metric_items_local_computed;
         ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env_);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env_);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;
      for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
        for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
          num_passed += Metrics_is_pass_threshold_noshrink(metrics, index, iE, jE, env_);
        }
      }
    } // for index

  //--------------------
  } else { // if (env_.num_way() == NumWay::_3)
  //--------------------

    for (NML_t index = 0; index <  metrics.num_metric_items_local_computed;
         ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env_);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env_);
      const NV_t kG = Metrics_coords_getG(metrics, index, 2, env_);
      // Skip if in vector pad.
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;
      for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
        for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
          for (int kE = 0; kE < env_.ijkE_max(); ++kE) {
            num_passed += Metrics_is_pass_threshold_noshrink(metrics, index, iE, jE, kE, env_);
          }
        }
      }
    } // for index

  //--------------------
  } // if (env_.num_way() == NumWay::_2)
  //--------------------

  if (num_passed != num_written_last_write_)
    fprintf(stderr, "Incorrect number of metrics written: "
      "actual %lld expected %lld\n",
      static_cast<long long int>(num_written_last_write_),
      static_cast<long long int>(num_passed));

  COMET_INSIST(num_passed == num_written_last_write_ &&
               "Incorrect number of metrics written.");

  // Close file and reopen file for write.

  close_();
  if (is_leaving_files_open_()) {
    open_("ab"); 

    // Set fpos to end of file for future writes.

    offset = bytes_(num_written_);
    const int success2 = fseek(file_, offset, SEEK_SET);
    COMET_INSIST(0 == success2 && "Unable to fseek to correct file position.");
  }
}

//-----------------------------------------------------------------------------
/*!
 * \brief Check the metrics elements most recently written to the file.
 *
 */
void MetricsIO::check_file(GMMetrics& metrics) {
  if (!env_.is_double_prec())
    check_file_impl_<float>(metrics);
  else
    check_file_impl_<double>(metrics);
}

//-----------------------------------------------------------------------------

static int num_digits_needed(int n) {
  COMET_ASSERT(n >= 1);
  int r = 0;
  for (int tmp = 1; ; tmp*=10, ++r) {
    if (tmp > n)
      break;
  }
  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Open (one-per-rank set of) metrics files.
 *
 */
void MetricsIO::open_(const char* mode) {

  if (! env_.is_proc_active())
    return;

  COMET_INSIST(!file_);

  // Prepare to form filename.

  const int num_digits_proc = num_digits_needed(env_.num_proc());
  const int num_digits_phase = num_digits_needed(env_.num_phase());
  const int num_digits_stage = num_digits_needed(env_.num_stage());

  const char* path_stub = path_stub_str_.c_str();

  const size_t len = strlen(path_stub);
  const size_t len_path_max = len + num_digits_proc + num_digits_phase
    + num_digits_stage + 100;

  char* path = (char*)malloc(len_path_max * sizeof(char));
  COMET_INSIST(path);

  // Create a formt string for forming the filename.
  // It will ensure the embedded numbers in the filenames have enough digits.

  enum {LEN_FORMAT_MAX = 100};
  char format[LEN_FORMAT_MAX];

  if (is_leaving_files_open_()) {

    // Filename only stores proc num.

    sprintf(format, "%%s" "_" "%%0" "%i" "i" ".bin", num_digits_proc);
    sprintf(path, format, path_stub, env_.proc_num());

  } else { // ! is_leaving_files_open_()

    // Filename stores proc num, phase num, stage num.

    sprintf(format, "%%s" "_" "%%0" "%i" "i"
                          "_" "%%0" "%i" "i"
                          "_" "%%0" "%i" "i" ".bin",
            num_digits_proc, num_digits_phase, num_digits_stage);
    sprintf(path, format, path_stub,
            env_.proc_num(), env_.phase_num(), env_.stage_num());

  } // is_leaving_files_open_()

  // Do the open.

  file_ = fopen(path, mode);
  COMET_INSIST(NULL != file_ && "Unable to open file.");

  free(path);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Close previously opened file.
 *
 */
void MetricsIO::close_() {

  if (! env_.is_proc_active())
    return;

  COMET_INSIST(is_allocated_);
  COMET_INSIST(file_);

  fclose(file_);

  file_ = NULL;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Static function to check writability to file system (all ranks).
 *
 */
bool MetricsIO::can_write_file(const char* path_stub, CEnv& env) {

  // NOTE: allowing flexibility to check on all ranks, not just active ranks.

  // Form filename

  const size_t len = strlen(path_stub);
  const size_t len_path_max = len + 100;

  char* path = (char*)malloc(len_path_max * sizeof(char));
  COMET_INSIST(path);

  const int num_digits_proc = num_digits_needed(env.num_proc());

  enum {LEN_FORMAT_MAX = 100};
  char format[LEN_FORMAT_MAX];

  sprintf(format, "%s0%ii.IO_TEST_FILE.bin", "%s_%", num_digits_proc);

  sprintf(path, format, path_stub, env.proc_num());

  // Do the open.

  FILE* const file = fopen(path, "wb");

  bool result = NULL != file;

  // Close and complete.

  if (file) {
    fclose(file);
    result = result && 0 == remove(path);
  }

  free(path);

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
