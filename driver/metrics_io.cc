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

#include "env.hh"
#include "metrics.hh"
#include "decomp_mgr.hh"

#include "driver.hh"
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

  typedef size_t IndexMax_t;
  const auto index_max = static_cast<IndexMax_t>(metrics.num_vector_active) *
    (env_.is_metric_type_bitwise() ? 2 : 1) - 1;
  COMET_INSIST_INTERFACE(&env, index_max ==
     static_cast<IndexMax_t>(static_cast<IntForFile_t>(index_max)) &&
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
static void MetricsIO_write_impl_(
  GMMetrics& metrics, FILE* file, size_t& num_written_, CEnv& env) {
  COMET_INSIST(file);

  if (! env.is_proc_active())
    return;

  // Due to field-redundancy, only results from some processors are needed.
  if (env.proc_num_field() != 0)
    return;

  //--------------------
  if (env.data_type_metrics() == DataTypeId::FLOAT &&
      env.num_way() == NumWay::_2) {
  //--------------------

    SingleMetricIO writer(file, metrics, env);

    COMET_INSIST(!env.thresholds().is_multi());

    for (NML_t index = 0; index < metrics.num_metrics_local; ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;

      const auto value = Metrics_elt_const<Float_t>(metrics, index, env);

      if (!env.thresholds().is_pass(value))
        continue;

      /// Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li): value: %.17e\n",
          iG, jG, static_cast<double>(value));
      else
        writer.write(iG, jG, value);

    } // for index

    num_written_ += writer.num_written();

  //--------------------
  } else if (env.data_type_metrics() == DataTypeId::FLOAT) {
  //--------------------
    COMET_INSIST(env.num_way() == NumWay::_3);

    SingleMetricIO writer(file, metrics, env);

    COMET_INSIST(!env.thresholds().is_multi());

    for (NML_t index = 0; index < metrics.num_metrics_local; ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env);
      const NV_t kG = Metrics_coords_getG(metrics, index, 2, env);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;

      const auto value = Metrics_elt_const<Float_t>(metrics, index, env);

      if (!env.thresholds().is_pass(value))
        continue;

      // Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li,%li): value: %.17e\n",
          iG, jG, kG, value);
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

    // Loop over metric items.
    for (NML_t index = 0; index < metrics.num_metric_items_local_computed;
         ++index) {
      const MetricItemCoords_t coords = metrics.coords_value(index);
      const NV_t iG = CoordsInfo::getiG(coords, metrics, env);
      const NV_t jG = CoordsInfo::getjG(coords, metrics, env);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;
      // Loop over metric entries in this item.
      for (int entry_num = 0; entry_num < env.num_entries_per_metric_item();
           ++entry_num) {
        const int iE = CoordsInfo::getiE(coords, entry_num, metrics, env);
        const int jE = CoordsInfo::getjE(coords, entry_num, metrics, env);

        //const GMFloat value = Metrics_ccc_duo_get_2(metrics,
        //  index, iE, jE, env);
        //TODO: obtain result here in correct precision !!!
        const Float_t value = Metrics_ccc_duo_get_2(metrics,
          index, entry_num, env);

        // For non-shrink case, metric not thresholded yet, so threshold.
        if (!Metrics_is_pass_threshold(metrics, index, iE, jE, env))
          continue;

        // Output the value.
        if (stdout == file)
          fprintf(file, "element (%li,%li) (%i,%i): value: %.17e\n",
            iG, jG, iE, jE, static_cast<double>(value));
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

    // Loop over metric items.
    for (NML_t index = 0; index < metrics.num_metric_items_local_computed;
         ++index) {
      const MetricItemCoords_t coords = metrics.coords_value(index);
      const NV_t iG = CoordsInfo::getiG(coords, metrics, env);
      const NV_t jG = CoordsInfo::getjG(coords, metrics, env);
      const NV_t kG = CoordsInfo::getkG(coords, metrics, env);
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

        //TODO: obtain result here in correct precision !!!
        const Float_t value = Metrics_ccc_duo_get_3(metrics,
          index, entry_num, env);

        // For non-shrink case, metric not thresholded yet, so threshold.
        if (!Metrics_is_pass_threshold(metrics, index, iE, jE, kE, env))
          continue;

        // Output the value.
        if (stdout == file)
          fprintf(file, "element (%li,%li,%li) (%i,%i,%i): value: %.17e\n",
            iG, jG, kG, iE, jE, kE, static_cast<double>(value));
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
static void MetricsIO_write_(
  GMMetrics& metrics, FILE* file, size_t& num_written_, CEnv& env) {

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
void MetricsIO::check_file(GMMetrics& metrics) {

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

  for (size_t index_f = 0; index_f < num_written_last_write_; ++index_f) {

    //--------------------
    if (env_.num_way() == NumWay::_2) {
    //--------------------

      // Read current metric from file.

      MetricForFile<NumWay::_2> metric;
      metric.read(file_, env_);

      const NV_t iG = metric.iG(env_);
      const NV_t jG = metric.jG(env_);

      const int iE = metric.iE(env_);
      const int jE = metric.jE(env_);

      // Access current metric from memory.

      const auto metric_value = env_.is_shrink() ?
        Metrics_ccc_duo_get_2(metrics, index_f, 0, env_) :
        GMMetrics_get_2(metrics, iG, jG, iE, jE, env_);

      // Compare.

      const NML_t index = Metrics_index_2(metrics, iG, jG, env_);

      bool do_coords_match = true;
      if (env_.is_shrink()) {
        MetricItemCoords_t coords = metrics.coords_value(index);
        do_coords_match =
          CoordsInfo::getiG(coords, metrics, env_) == iG &&
          CoordsInfo::getjG(coords, metrics, env_) == jG;
      }

//FIXTHRESHOLD - how to handle not is_shrink case
#if 0
CASE 1: not threshold_tc. then MetricFormat::PACKED_DOUBLE. then need to access full table here (see earlier in file for how to do this). Note this is sort of inefficient, done 4-8 times when only need once.
CASE 2: threshold_tc (MetricFormat::SINGLE), not is_shrink: just get table entry, as above, check for nonzero and whether equal
CASE 3: threshold_tc, is_shrink: file should match metrics exactly, no need to apply threshold

      bool pass_threshold = true;
      if (!env->is_shrink()) {
        if (env->is_threshold_tc()) {
          typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
          const auto ttable = Metrics_elt_const<TTable_t>(metrics, index_f, env);
          pass_threshold = !env->thresholds().is_zero(ttable, iE, jE);
            continue;
        } else { // ! env->is_threshold_tc()
          typedef Tally2x2<MetricFormat::PACKED_DOUBLE> TTable_t;
          const auto ttable = Metrics_elt_const<TTable_t>(metrics, index_f, env);
          pass_threshold = env->thresholds().is_pass(ttable, iE, jE);
        } // if (env->is_threshold_tc())
      } // if (!env->is_shrink())

#endif

      //const bool pass_threshold = env_.pass_threshold(metric_value);
      //const bool pass_threshold = env_.thresholds().is_pass(metric_value);

      const bool pass_threshold = Metrics_is_pass_threshold(metrics, index, iE, jE, env_);

      const bool is_correct =
        static_cast<FloatForFile_t>(metric_value) == metric.value() &&
        pass_threshold && do_coords_match;

      num_incorrect += !is_correct;

      if (num_incorrect < 10 && !is_correct)
        fprintf(stderr, "Incorrect metric value: "
          "element %zu %zu  %i %i actual %.17e expected %.17e pass_threshold %i\n",
          iG, jG, iE, jE, static_cast<double>(metric_value),
          static_cast<double>(metric.value()), pass_threshold);

    //--------------------
    } else { // if (env_.num_way() == NumWay::_3)
    //--------------------

      // Read current metric from file.

      MetricForFile<NumWay::_3> metric;
      metric.read(file_, env_);

      const NV_t iG = metric.iG(env_);
      const NV_t jG = metric.jG(env_);
      const NV_t kG = metric.kG(env_);

      const int iE = metric.iE(env_);
      const int jE = metric.jE(env_);
      const int kE = metric.kE(env_);

      // Access current metric from memory.

      const auto metric_value = env_.is_shrink() ?
        Metrics_ccc_duo_get_3(metrics, index_f, 0, env_) :
        GMMetrics_get_3(metrics, iG, jG, kG, iE, jE, kE, env_);

      // Compare.

      const NML_t index = Metrics_index_3(metrics, iG, jG, kG, env_);

      bool do_coords_match = true;
      if (env_.is_shrink()) {
        MetricItemCoords_t coords = metrics.coords_value(index);
        do_coords_match =
          CoordsInfo::getiG(coords, metrics, env_) == iG &&
          CoordsInfo::getjG(coords, metrics, env_) == jG &&
          CoordsInfo::getkG(coords, metrics, env_) == kG;
      }

//FIXTHRESHOLD - how to handle not is_shrink case - see above for 2-way

#if 0
      bool pass_threshold = true;
      if (!env->is_shrink()) {
        if (env->is_threshold_tc()) {
          typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
          const auto ttable = Metrics_elt_const<TTable_t>(metrics, index_f, env);
          pass_threshold = !env->thresholds().is_zero(ttable, iE, jE, kE);
            continue;
        } else { // ! env->is_threshold_tc()
          typedef Tally2x2<MetricFormat::PACKED_DOUBLE> TTable_t;
          const auto ttable = Metrics_elt_const<TTable_t>(metrics, index_f, env);
          pass_threshold = env->thresholds().is_pass(ttable, iE, jE, kE);
        } // if (env->is_threshold_tc())
      } // if (!env->is_shrink())

#endif

      //const bool pass_threshold = env_.pass_threshold(metric_value);
      //const bool pass_threshold = env_.thresholds().is_pass(metric_value);
      const bool pass_threshold = Metrics_is_pass_threshold(metrics, index, iE, jE, kE, env_);

      const bool is_correct =
        static_cast<FloatForFile_t>(metric_value) == metric.value() &&
        pass_threshold && do_coords_match;

      num_incorrect += !is_correct;

      if (!is_correct && num_incorrect < 10)
        fprintf(stderr, "Incorrect metric value: "
          "element %zu %zu %zu  %i %i %i actual %.17e expected %.17e pass_threshold %i\n",
          iG, jG, kG, iE, jE, kE, static_cast<double>(metric_value),
          static_cast<double>(metric.value()), pass_threshold);

    //--------------------
    } // if (env_.num_way() == NumWay::_2)
    //--------------------

  } // index_f

  COMET_INSIST(0 == num_incorrect &&
    "Incorrect values detected in output file - partial list is displayed.");

  //----------------------------------------
  // Test 2: count metrics (in memory) passing threshold, compare to num written.
  //----------------------------------------

  size_t num_passed = 0;

  //--------------------
  if (env_.num_way() == NumWay::_2) {
  //--------------------

    for (NML_t index = 0; index <  metrics.num_metric_items_local_computed;
         ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env_);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env_);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;
      if (env_.is_shrink()) {
        num_passed += 1;
      } else {
//FIXTHRESHOLD - get metrics table, use for the individual comparisons
        for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
          for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
            //const auto metric_value =
            //  GMMetrics_get_2(metrics, index, iE, jE, env_);

#if 0
            bool pass_threshold = true;
            if (env->is_threshold_tc()) {
              typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
              const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
              pass_threshold = !env->thresholds().is_zero(ttable, iE, jE);
                continue;
            } else { // ! env->is_threshold_tc()
              typedef Tally2x2<MetricFormat::PACKED_DOUBLE> TTable_t;
              const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
              pass_threshold = env->thresholds().is_pass(ttable, iE, jE);
            } // if (env->is_threshold_tc())
#endif

            //num_passed += env_.pass_threshold(metric_value);
            //num_passed += env_.thresholds().is_pass(metric_value);
            num_passed += Metrics_is_pass_threshold(metrics, index, iE, jE, env_);
          }
        }
      } // if is_shrink
    } // for index

  //--------------------
  } else { // if (env_.num_way() == NumWay::_3)
  //--------------------

    for (NML_t index = 0; index <  metrics.num_metric_items_local_computed;
         ++index) {
      const NV_t iG = Metrics_coords_getG(metrics, index, 0, env_);
      const NV_t jG = Metrics_coords_getG(metrics, index, 1, env_);
      const NV_t kG = Metrics_coords_getG(metrics, index, 2, env_);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;
      if (env_.is_shrink()) {
        num_passed += 1;
      } else {
//FIXTHRESHOLD - get metrics table, use for the individual comparisons
        for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
          for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
            for (int kE = 0; kE < env_.ijkE_max(); ++kE) {
              //const auto metric_value =
              //  GMMetrics_get_3(metrics, index, iE, jE, kE, env_);

#if 0
              bool pass_threshold = true;
              if (env->is_threshold_tc()) {
                typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
                const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
                pass_threshold = !env->thresholds().is_zero(ttable, iE, jE, kE);
                  continue;
              } else { // ! env->is_threshold_tc()
                typedef Tally2x2<MetricFormat::PACKED_DOUBLE> TTable_t;
                const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
                pass_threshold = env->thresholds().is_pass(ttable, iE, jE, kE);
              } // if (env->is_threshold_tc())
#endif

              //num_passed += env_.pass_threshold(metric_value);
              //num_passed += env_.thresholds().is_pass(metric_value);
              num_passed += Metrics_is_pass_threshold(metrics, index, iE, jE, kE, env_);
            }
          }
        }
      } // if is_shrink
    } // for index

  //--------------------
  } // if (env_.num_way() == NumWay::_2)
  //--------------------

  if (num_passed != num_written_last_write_)
    fprintf(stderr, "Incorrect number of metrics written: "
      "actual %zu expected %zu\n", num_written_last_write_, num_passed);

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
 * \brief Open (one-per-rank set of) metrics files.
 *
 */
void MetricsIO::open_(const char* mode) {

  if (! env_.is_proc_active())
    return;

  COMET_INSIST(!file_);

  // Prepare to form filename.

  const char* path_stub = path_stub_str_.c_str();

  size_t len = strlen(path_stub);
  char* path = (char*)malloc((len+100) * sizeof(char));
  COMET_INSIST(path);

  int num_digits_proc = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits_proc) {
    if (tmp > env_.num_proc()) {
      break;
    }
  }

  // Create a formt string for forming filename.
  char format[100];

  if (is_leaving_files_open_()) {

    //sprintf(format, "%s0%ii.bin", "%s_%", num_digits);
    sprintf(format, "%%s"
                    "_"
                    "%%0" "%i" "i"
                    ".bin",
            num_digits_proc);
    sprintf(path, format, path_stub, env_.proc_num());

    //sprintf(path, "%s_%*i.bin", path_stub, num_digits_proc, env_.proc_num());

  } else { // ! is_leaving_files_open_()

    int num_digits_phase = 0;
    for (int tmp = 1; ; tmp*=10, ++num_digits_phase) {
      if (tmp > env_.num_phase())
        break;
    }

    int num_digits_stage = 0;
    for (int tmp = 1; ; tmp*=10, ++num_digits_stage) {
      if (tmp > env_.num_stage())
        break;
    }

    sprintf(format, "%%s"
                    "_"
                    "%%0" "%i" "i"
                    "_"
                    "%%0" "%i" "i"
                    "_"
                    "%%0" "%i" "i"
                    ".bin",
            num_digits_proc, num_digits_phase, num_digits_stage);
    sprintf(path, format, path_stub,
            env_.proc_num(), env_.phase_num(), env_.stage_num());

    //sprintf(path, "%s_%*i_%*i_%*i.bin", path_stub,
    //  num_digits_proc, env_.proc_num(),
    //  num_digits_phase, env_.phase_num(),
    //  num_digits_stage, env_.stage_num());

  } // is_leaving_files_open_()

  // Do open

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

  size_t len = strlen(path_stub);
  char* path = (char*)malloc((len+100) * sizeof(char));
  COMET_INSIST(path);

  int num_digits = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits) {
    if (tmp > env.num_proc()) {
      break;
    }
  }

  char format[100];
  sprintf(format, "%s0%ii.IO_TEST_FILE.bin", "%s_%", num_digits);

  sprintf(path, format, path_stub, env.proc_num());

  // Do open

  FILE* const file = fopen(path, "wb");
  //COMET_INSIST(file);

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
