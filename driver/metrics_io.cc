//-----------------------------------------------------------------------------
/*!
 * \file   metrics_io.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O utilities for metrics.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

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
// MetricIO member definitions.

//-----------------------------------------------------------------------------
/// \brief Constructor for MetricIO class

MetricIO::MetricIO(FILE* file, GMMetrics& metrics, CEnv& env)
  : env_(env)
  , file_(file)
  , data_type_(env.data_type_metrics())
  , num_way_(env.num_way())
  , num_written_(0) {

  if (stdout != file_)
    COMET_INSIST_INTERFACE(&env, metrics.num_vector_active ==
                  (uint32_t)metrics.num_vector_active &&
                  "Too many vectors for output format.");
}

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CZEK 2-way case.

void MetricIO::write(size_t iG, size_t jG, GMFloat value) const {

  bool success = true;
  size_t bytes_written = 0;

  success = success && write_<uint32_t>(iG, bytes_written);
  success = success && write_<uint32_t>(jG, bytes_written);
  success = success && write_<GMFp32>(value, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
  }

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CZEK 3-way case.

void MetricIO::write(size_t iG, size_t jG, size_t kG,
                         GMFloat value) const {

  bool success = true;
  size_t bytes_written = 0;

  success = success && write_<uint32_t>(iG, bytes_written);
  success = success && write_<uint32_t>(jG, bytes_written);
  success = success && write_<uint32_t>(kG, bytes_written);
  success = success && write_<GMFp32>(value, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
}

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CCC/DUO 2-way case.

void MetricIO::write(size_t iG, size_t jG, int iE, int jE,
                         GMFloat value) const {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);

  write(iE + 2 * iG, jE + 2 * jG, value);
}

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CCC/DUO 3-way case.

void MetricIO::write(size_t iG, size_t jG, size_t kG,
                         int iE, int jE, int kE, GMFloat value) const {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);

  write(iE + 2 * iG, jE + 2 * jG, kE + 2 * kG, value);
}

//=============================================================================
// Helper functions for MetricsIO class.

//-----------------------------------------------------------------------------
/// \brief Write to file, GMTally2x2 case, implementation.

template<int COUNTED_BITS_PER_ELT>
static void MetricsIO_write_tally2x2_bin_impl_(
  GMMetrics* metrics, FILE* file, size_t& num_written_, CEnv* env) {
  COMET_INSIST(metrics && file && env);
  COMET_INSIST(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);
  COMET_INSIST(env->num_way() == NUM_WAY::_2);
  COMET_INSIST(stdout != file);

  MetricIO writer(file, *metrics, *env);

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, iE, jE) entries (potentially) stored in buffer
  const size_t num_buf = 4 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // iG, jG, iE, jE, and value

  char* const do_out_buf = (char*)malloc(num_buf*sizeof(*do_out_buf));
  int* const iG_buf = (int*)malloc(num_buf*sizeof(*iG_buf));
  int* const jG_buf = (int*)malloc(num_buf*sizeof(*jG_buf));
  int* const ijE_buf = (int*)malloc(num_buf*sizeof(*ijE_buf));
  GMFloat* const value_buf = (GMFloat*)malloc(num_buf*sizeof(*value_buf));

  for (int i=0; i<(int)num_buf; ++i)
    do_out_buf[i] = 0;

  // Process num_buf_ind index values at a time
  for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
       ind_base += num_buf_ind) {
    // Largest index value to visit for this loop trip.
    const size_t ind_max = utils::min(metrics->num_elts_local,
                                      ind_base + num_buf_ind);

    // Fill buffer
#pragma omp parallel for schedule(dynamic,1000)
    for (size_t index = ind_base; index < ind_max; ++index) {
      // Do any of the values exceed the threshold
      if (Metrics_ccc_duo_get_threshold_2<COUNTED_BITS_PER_ELT>(
            *metrics, index, *env)) {
        for (int iE = 0; iE < 2; ++iE) {
          for (int jE = 0; jE < 2; ++jE) {
            const GMFloat value =
              Metrics_ccc_duo_get_2<COUNTED_BITS_PER_ELT>( *metrics, index,
               iE, jE, *env);
            if (env->pass_threshold(value)) {
              const size_t iG =
                GMMetrics_coord_global_from_index(metrics, index, 0, env);
              const size_t jG =
                GMMetrics_coord_global_from_index(metrics, index, 1, env);
              const char do_out = iG < metrics->num_vector_active &&
                                  jG < metrics->num_vector_active;
              const size_t ind_buf = jE + 2*(iE + 2*(index-ind_base));
              COMET_ASSERT(ind_buf < num_buf);
              do_out_buf[ind_buf] = do_out;
              iG_buf[ind_buf] = iG;
              jG_buf[ind_buf] = jG;
              ijE_buf[ind_buf] = iE + 2*jE;
              value_buf[ind_buf] = value;
            } // if
          } // jE
        } // iE
      }
    } // pragma omp / for index

    // Flush buffer

    // Number of buffer entries to visit
    const size_t ind_buf_max = 4 * (ind_max - ind_base);
    // Use 4 byte integer ptr to check 4 chars at a time, for speed
    typedef int multi_t;
    //assert(sizeof(multi_t) == 4 * sizeof(char));
    const multi_t* const do_out_ptr_max = (multi_t*)(do_out_buf + ind_buf_max);
    multi_t* do_out_ptr = (multi_t*)do_out_buf;
    for (; do_out_ptr < do_out_ptr_max;) {
      // Do any of the 4 need to be output
      if (*(do_out_ptr++)) {
        for (int i=0; i<4; ++i) {
          const size_t ind_buf = (do_out_ptr - (multi_t*)do_out_buf - 1)*4 + i;
          COMET_ASSERT(ind_buf < num_buf);
          if (do_out_buf[ind_buf]) {
            const int iE = ijE_buf[ind_buf] % 2;
            const int jE = ijE_buf[ind_buf] / 2;
            const size_t iG = iG_buf[ind_buf];
            const size_t jG = jG_buf[ind_buf];
            writer.write(iG, jG, iE, jE, value_buf[ind_buf]);
            // Reset buffer entry to false
            do_out_buf[ind_buf] = 0;
          }
        }
      }
    } // ind_buf

  } // ind_base

  num_written_ += writer.num_written();

  free(do_out_buf);
  free(iG_buf);
  free(jG_buf);
  free(ijE_buf);
  free(value_buf);
}

//-----------------------------------------------------------------------------
/// \brief Write to file, GMTally2x2 case.

static void MetricsIO_write_tally2x2_bin_(
  GMMetrics* metrics, FILE* file, size_t& num_written_, CEnv* env) {
  COMET_INSIST(metrics && file && env);
  COMET_INSIST(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);
  COMET_INSIST(env->num_way() == NUM_WAY::_2);
  COMET_INSIST(stdout != file);

  if (env->metric_type() == MetricType::CCC)
    MetricsIO_write_tally2x2_bin_impl_<CBPE::CCC>(metrics, file,
      num_written_, env);
  else
    MetricsIO_write_tally2x2_bin_impl_<CBPE::DUO>(metrics, file,
      num_written_, env);
}

//-----------------------------------------------------------------------------
/// \brief Write to file, GMTally4x2 case, implementation.

template<int COUNTED_BITS_PER_ELT>
static void MetricsIO_write_tally4x2_bin_impl_(
  GMMetrics* metrics, FILE* file, size_t& num_written_, CEnv* env) {
  COMET_INSIST(metrics && file && env);
  COMET_INSIST(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  COMET_INSIST(env->num_way() == NUM_WAY::_3);
  COMET_INSIST(stdout != file);
  COMET_INSIST(env->is_metric_type_bitwise());

  MetricIO writer(file, *metrics, *env);

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, iE, jE) entries (potentially) stored in buffer
  const size_t num_buf = 8 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // iG, jG, kG, iE, jE, kE, and value

  char* const do_out_buf = (char*)malloc(num_buf*sizeof(*do_out_buf));
  int* const iG_buf = (int*)malloc(num_buf*sizeof(*iG_buf));
  int* const jG_buf = (int*)malloc(num_buf*sizeof(*jG_buf));
  int* const kG_buf = (int*)malloc(num_buf*sizeof(*kG_buf));
  int* const ijkE_buf = (int*)malloc(num_buf*sizeof(*ijkE_buf));
  GMFloat* const value_buf = (GMFloat*)malloc(num_buf*sizeof(*value_buf));

  for (int i=0; i<(int)num_buf; ++i)
    do_out_buf[i] = 0;

  // Process num_buf_ind index values at a time
  for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
       ind_base += num_buf_ind) {
    const size_t ind_max = utils::min(metrics->num_elts_local,
                                     ind_base + num_buf_ind);
    // Fill buffer
#pragma omp parallel for schedule(dynamic,1000)
    for (size_t index = ind_base; index < ind_max; ++index) {
      // Do any of the values exceed the threshold
      if (Metrics_ccc_duo_get_threshold_3<COUNTED_BITS_PER_ELT>(
             *metrics, index, *env)) {
        for (int iE = 0; iE < 2; ++iE) {
          for (int jE = 0; jE < 2; ++jE) {
            for (int kE = 0; kE < 2; ++kE) {
              const GMFloat value =
                Metrics_ccc_duo_get_3<COUNTED_BITS_PER_ELT>(
                  *metrics, index, iE, jE, kE, *env);
              if (env->pass_threshold(value)) {
                const size_t iG =
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t jG =
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const size_t kG =
                  GMMetrics_coord_global_from_index(metrics, index, 2, env);
                const char do_out = iG < metrics->num_vector_active &&
                                    jG < metrics->num_vector_active &&
                                    kG < metrics->num_vector_active;
                const size_t ind_buf = jE + 2*(iE + 2*(kE +2*(index-ind_base)));
                do_out_buf[ind_buf] = do_out;
                iG_buf[ind_buf] = iG;
                jG_buf[ind_buf] = jG;
                kG_buf[ind_buf] = kG;
                ijkE_buf[ind_buf] = iE + 2*(jE + 2*kE);
                value_buf[ind_buf] = value;
              } // if
            } // kE
          } // jE
        } // iE
      }
    } // for index

    // Flush buffer

    // Number of buffer entries to visit
    const size_t ind_buf_max = 8 * (ind_max - ind_base);
    // Use 4 byte integer ptr to check 4 chars at a time, for speed
    typedef size_t multi_t;
    //assert(sizeof(multi_t) == 8 * sizeof(char));
    const multi_t* const do_out_ptr_max = (multi_t*)(do_out_buf + ind_buf_max);
    multi_t* do_out_ptr = (multi_t*)do_out_buf;
    for (; do_out_ptr < do_out_ptr_max;) {
      // Do any of the 8 need to be output
      if (*(do_out_ptr++)) {
        for (int i=0; i<8; ++i) {
          const size_t ind_buf = (do_out_ptr - (multi_t*)do_out_buf - 1)*8 + i;
          if (do_out_buf[ind_buf]) {
            const int iE = ijkE_buf[ind_buf] % 2;
            const int jE = (ijkE_buf[ind_buf] / 2) % 2;
            const int kE = ijkE_buf[ind_buf] / 4;
            const size_t iG = iG_buf[ind_buf];
            const size_t jG = jG_buf[ind_buf];
            const size_t kG = kG_buf[ind_buf];
            writer.write(iG, jG, kG, iE,jE,kE, value_buf[ind_buf]);
            // Reset buffer entry to false
            do_out_buf[ind_buf] = 0;
          }
        }
      }
    } // ind_buf
  } // ind_base

  num_written_ += writer.num_written();

  free(do_out_buf);
  free(iG_buf);
  free(jG_buf);
  free(kG_buf);
  free(ijkE_buf);
  free(value_buf);
}

//-----------------------------------------------------------------------------
/// \brief Write to file, GMTally4x2 case.

static void MetricsIO_write_tally4x2_bin_(
  GMMetrics* metrics, FILE* file, size_t& num_written_, CEnv* env) {
  COMET_INSIST(metrics && file && env);
  COMET_INSIST(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  COMET_INSIST(env->num_way() == NUM_WAY::_3);
  COMET_INSIST(stdout != file);

  if (env->metric_type() == MetricType::CCC)
    MetricsIO_write_tally4x2_bin_impl_<CBPE::CCC>(metrics, file,
      num_written_, env);
  else
    MetricsIO_write_tally4x2_bin_impl_<CBPE::DUO>(metrics, file,
      num_written_, env);
}

//-----------------------------------------------------------------------------
/// \brief Write to file or print the metrics.

static void MetricsIO_write_(
  GMMetrics* metrics, FILE* file, size_t& num_written_, CEnv* env) {
  COMET_INSIST(metrics && file && env);

  if (! env->is_proc_active())
    return;

  // Due to redundancy, only results from some processors are needed.
  if (env->proc_num_field() != 0)
    return;

  //----------
  if (env->data_type_metrics() == GM_DATA_TYPE_FLOAT &&
      env->num_way() == NUM_WAY::_2) {
  //----------

    MetricIO writer(file, *metrics, *env);

    for (size_t index = 0; index < metrics->num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t jG =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      if (iG >= metrics->num_vector_active ||
          jG >= metrics->num_vector_active)
        continue;
      const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);

      if (!env->pass_threshold(value))
        continue;
      /// Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li): value: %.17e\n",
          iG, jG, value);
      else
        writer.write(iG, jG, value);
    } // for index

    num_written_ += writer.num_written();

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_FLOAT &&
           env->num_way() == NUM_WAY::_3) {
  //----------

    MetricIO writer(file, *metrics, *env);

    for (size_t index = 0; index < metrics->num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t jG =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      const size_t kG =
        GMMetrics_coord_global_from_index(metrics, index, 2, env);
      if (iG >= metrics->num_vector_active ||
          jG >= metrics->num_vector_active ||
          kG >= metrics->num_vector_active)
        continue;
      const auto value = Metrics_elt_const<GMFloat>(*metrics, index, *env);
      if (!env->pass_threshold(value))
        continue;

      // Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li,%li): value: %.17e\n",
          iG, jG, kG, value);
      else
        writer.write(iG, jG, kG, value);
    } // for index

    num_written_ += writer.num_written();

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2 &&
             stdout != file) {
  //----------
    COMET_INSIST(env->num_way() == NUM_WAY::_2);
    COMET_INSIST(env->is_metric_type_bitwise());

    MetricsIO_write_tally2x2_bin_(metrics, file, num_written_, env);

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2) {
  //----------
    COMET_INSIST(env->num_way() == NUM_WAY::_2);
    COMET_INSIST(env->is_metric_type_bitwise());

    MetricIO writer(file, *metrics, *env);

    size_t index = 0;
    for (index = 0; index < metrics->num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t jG =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      if (iG >= metrics->num_vector_active ||
          jG >= metrics->num_vector_active)
        continue;
      int num_out_this_line = 0;
      for (int iE = 0; iE < 2; ++iE) {
        for (int jE = 0; jE < 2; ++jE) {
          const GMFloat value = Metrics_ccc_duo_get_2(*metrics,
            index, iE, jE, *env);
          if (!env->pass_threshold(value))
            continue;

          // Output the value.

          if (num_out_this_line == 0)
            fprintf(file, "element (%li,%li): values:", iG, jG);

          fprintf(file, " %i %i %.17e", iE, jE, value);

          num_out_this_line++;

        } // jE
      } // iE
      if (num_out_this_line > 0)
        fprintf(file, "\n");
    } // for index

    num_written_ += writer.num_written();

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2 &&
             stdout != file) {
  //----------
    COMET_INSIST(env->num_way() == NUM_WAY::_3);
    COMET_INSIST(env->is_metric_type_bitwise());

    MetricsIO_write_tally4x2_bin_(metrics, file, num_written_, env);

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2) {
  //----------
    COMET_INSIST(env->num_way() == NUM_WAY::_3);
    COMET_INSIST(env->is_metric_type_bitwise());

    MetricIO writer(file, *metrics, *env);

    for (size_t index = 0; index < metrics->num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t jG =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      const size_t kG =
        GMMetrics_coord_global_from_index(metrics, index, 2, env);
      if (iG >= metrics->num_vector_active ||
          jG >= metrics->num_vector_active ||
          kG >= metrics->num_vector_active)
        continue;
      int num_out_this_line = 0;
      for (int iE = 0; iE < 2; ++iE) {
        for (int jE = 0; jE < 2; ++jE) {
          for (int kE = 0; kE < 2; ++kE) {
            const GMFloat value = Metrics_ccc_duo_get_3(*metrics,
              index, iE, jE, kE, *env);
            if (!env->pass_threshold(value))
              continue;

            // Output the value.

            if (num_out_this_line == 0)
              fprintf(file, "element (%li,%li,%li): values:",
                iG, jG, kG);

            fprintf(file, " %i %i %i %.17e", iE, jE, kE, value);

            num_out_this_line++;

          } // kE
        } // jE
      } // iE
      if (num_out_this_line > 0)
        fprintf(file, "\n");
    } // for index

    num_written_ += writer.num_written();

  //----------
  } else {
  //----------

      COMET_INSIST(false && "Invalid data type.");

  //----------
  } // if
  //----------
}

//=============================================================================
// MetricsIO member definitions.

//-----------------------------------------------------------------------------
/// \brief Constructor for MetricsIO class.

MetricsIO::MetricsIO(const char* path_stub, int verbosity, CEnv& env)
  : env_(env)
  , path_stub_(path_stub ? path_stub : "")
  , file_(NULL)
  , verbosity_(verbosity)
  , num_written_(0)
  , num_written_last_write_(0) {

  if (! env_.is_proc_active())
    return;

  if (path_stub) {
    file_ = MetricsIO::open(path_stub, env_);
    COMET_INSIST(NULL != file_ && "Unable to open file.");
  }
}

//-----------------------------------------------------------------------------
/// \brief Destructor for MetricsIO class.

MetricsIO::~MetricsIO() {
  if (file_)
    fclose(file_);
}

//-----------------------------------------------------------------------------
/// \brief Write elements of GMMetrics object to file; if requested also print.

void MetricsIO::write(GMMetrics& metrics) {

  if (! env_.is_proc_active())
    return;

  // Output to file

  if (file_) {
    const size_t num_written_hold = num_written_;
    MetricsIO_write_(&metrics, file_, num_written_, &env_);
    num_written_last_write_ = num_written_ - num_written_hold;
  }

  // Output to stdout if requested

  if (verbosity_ > 1)
    MetricsIO_write_(&metrics, stdout, num_written_, &env_);
}

//-----------------------------------------------------------------------------
/// \brief Check the metrics elements most recently written to the file.

void MetricsIO::check_file(GMMetrics& metrics) {

  if (! env_.is_proc_active())
    return;

  // Due to redundancy, only results from some processors are needed.
  if (env_.proc_num_field() != 0)
    return;

  if (! file_)
    return;

  // Close file and reopen file for read.

  fclose(file_);
  file_ = MetricsIO::open(path_stub_.c_str(), env_, "rb");
  COMET_INSIST(NULL != file_ && "Unable to open file.");

  // Set fpos to beginning of last metrics object write.

  long offset = bytes_(num_written_ - num_written_last_write_);
  const int success1 = fseek(file_, offset, SEEK_SET);
  COMET_INSIST(0 == success1);

  // Loop over metrics stored in file, compare to metrics in memory.

  size_t num_incorrect = 0;

  for (size_t i = 0; i < num_written_last_write_; ++i) {

    if (env_.num_way() == NUM_WAY::_2) {

      MetricIO::Metric<NUM_WAY::_2> metric;
      MetricIO::read(metric, file_, env_);

      const size_t iG = metric.iG(env_);
      const size_t jG = metric.jG(env_);

      const int iE = metric.iE(env_);
      const int jE = metric.jE(env_);

      const MetricIO::Float_t metric_value =
        (MetricIO::Float_t)GMMetrics_get_2(metrics, iG, jG, iE, jE, env_);

      const bool is_correct = metric_value == metric.value &&
                              env_.pass_threshold(metric_value);
      num_incorrect += !is_correct;

      if (num_incorrect < 10 && !is_correct) {
        fprintf(stderr, "Incorrect metric value: "
          "element %zu %zu actual %.17e expected %.17e\n",
          iG, jG, (double)metric_value, (double)metric.value);

        fprintf(stderr, "Incorrect metric value: "
          "element %zu %zu actual %.17e expected %.17e\n",
          iG, jG, (double)metric_value, (double)metric.value);
      }

    } else { // if (env_.num_way() == NUM_WAY::_3)

      MetricIO::Metric<NUM_WAY::_3> metric;
      MetricIO::read(metric, file_, env_);

      const size_t iG = metric.iG(env_);
      const size_t jG = metric.jG(env_);
      const size_t kG = metric.kG(env_);

      const int iE = metric.iE(env_);
      const int jE = metric.jE(env_);
      const int kE = metric.kE(env_);

      const MetricIO::Float_t metric_value = (MetricIO::Float_t)
        GMMetrics_get_3(metrics, iG, jG, kG, iE, jE, kE, env_);

      const bool is_correct = metric_value == metric.value &&
                              env_.pass_threshold(metric_value);
      num_incorrect += !is_correct;

      if (!is_correct && num_incorrect < 10) {
        fprintf(stderr, "Incorrect metric value: "
          "element %zu %zu %zu actual %.17e expected %.17e\n",
          iG, jG, kG, (double)metric_value, (double)metric.value);
      }

    } // if (env_.num_way() == NUM_WAY::_2)

  } // i

  COMET_INSIST(0 == num_incorrect &&
    "Incorrect values detected in output file - partial list is displayed.");

  // Count metrics (in memory) passing threshold, compare to num written.

  size_t num_passed = 0;

  if (env_.num_way() == NUM_WAY::_2) {

    for (size_t index = 0; index <  metrics.num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(&metrics, index, 0, &env_);
      const size_t jG =
        GMMetrics_coord_global_from_index(&metrics, index, 1, &env_);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active)
        continue;
      for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
        for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
          const MetricIO::Float_t metric_value =
            (MetricIO::Float_t)GMMetrics_get_2(metrics, index, iE, jE, env_);
          num_passed += env_.pass_threshold(metric_value);
        }
      }
    }

  } else { // if (env_.num_way() == NUM_WAY::_3)

    for (size_t index = 0; index <  metrics.num_elts_local; ++index) {
      const size_t iG =
        GMMetrics_coord_global_from_index(&metrics, index, 0, &env_);
      const size_t jG =
        GMMetrics_coord_global_from_index(&metrics, index, 1, &env_);
      const size_t kG =
        GMMetrics_coord_global_from_index(&metrics, index, 2, &env_);
      if (iG >= metrics.num_vector_active ||
          jG >= metrics.num_vector_active ||
          kG >= metrics.num_vector_active)
        continue;
      for (int iE = 0; iE < env_.ijkE_max(); ++iE) {
        for (int jE = 0; jE < env_.ijkE_max(); ++jE) {
          for (int kE = 0; kE < env_.ijkE_max(); ++kE) {
            const GMFloat metric_value =
              GMMetrics_get_3(metrics, index, iE, jE, kE, env_);
            num_passed += env_.pass_threshold(metric_value);
          }
        }
      }
    }

  } // if (env_.num_way() == NUM_WAY::_2)

  if (num_passed != num_written_last_write_)
    fprintf(stderr, "Incorrect number of metrics written: "
      "actual %zu expected %zu\n", num_written_last_write_, num_passed);

  COMET_INSIST(num_passed == num_written_last_write_ &&
               "Incorrect number of metrics written.");

  // Close file and reopen file for write.

  fclose(file_);
  file_ = MetricsIO::open(path_stub_.c_str(), env_, "ab");
  COMET_INSIST(NULL != file_ && "Unable to open file.");

  // Set fpos to end of file.

  offset = bytes_(num_written_);
  const int success2 = fseek(file_, offset, SEEK_SET);
  COMET_INSIST(0 == success2);
}

//-----------------------------------------------------------------------------
/// \brief Static function to open (one-per-rank set of) metrics files.

FILE* MetricsIO::open(const char* path_stub, CEnv& env, const char* mode) {

  // NOTE: allowing flexibility to open on all ranks, not just active ranks.

  // Form filename

  size_t len = strlen(path_stub);
  char* path = (char*)malloc((len+50) * sizeof(char));
  COMET_INSIST(path);

  int num_digits = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits) {
    if (tmp > env.num_proc()) {
      break;
    }
  }

  char format[100];
  sprintf(format, "%s0%ii.bin", "%s_%", num_digits);

  sprintf(path, format, path_stub, env.proc_num());

  // Do open

  FILE* const file = fopen(path, mode);
  COMET_INSIST(file);
  free(path);

  return file;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
