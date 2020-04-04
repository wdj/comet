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

void MetricIO::write(size_t coord0, size_t coord1, GMFloat value) const {

  bool success = true;
  size_t bytes_written = 0;

  success = success && write_<uint32_t>(coord0, bytes_written);
  success = success && write_<uint32_t>(coord1, bytes_written);
  success = success && write_<GMFp32>(value, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
  }

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CZEK 3-way case.

void MetricIO::write(size_t coord0, size_t coord1, size_t coord2,
                         GMFloat value) const {

  bool success = true;
  size_t bytes_written = 0;

  success = success && write_<uint32_t>(coord0, bytes_written);
  success = success && write_<uint32_t>(coord1, bytes_written);
  success = success && write_<uint32_t>(coord2, bytes_written);
  success = success && write_<GMFp32>(value, bytes_written);

  COMET_INSIST(success && "File write failure.");
  num_written_++;
  COMET_INSIST(num_bytes_written_per_metric() == bytes_written);
}

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CCC/DUO 2-way case.

void MetricIO::write(size_t coord0, size_t coord1, int i0, int i1,
                         GMFloat value) const {
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  write(i0 + 2 * coord0, i1 + 2 * coord1, value);
}

//-----------------------------------------------------------------------------
/// \brief Write a metric value: CCC/DUO 3-way case.

void MetricIO::write(size_t coord0, size_t coord1, size_t coord2,
                         int i0, int i1, int i2, GMFloat value) const {
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  write(i0 + 2 * coord0, i1 + 2 * coord1, i2 + 2 * coord2, value);
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

  //enum {vals_per_index = 4;}

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, i0, i1) entries (potentially) stored in buffer
  const size_t num_buf = 4 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // coord0, coord1, i0, i1, and value

  char* const do_out_buf = (char*)malloc(num_buf*sizeof(*do_out_buf));
  int* const coord0_buf = (int*)malloc(num_buf*sizeof(*coord0_buf));
  int* const coord1_buf = (int*)malloc(num_buf*sizeof(*coord1_buf));
  int* const i01_buf = (int*)malloc(num_buf*sizeof(*i01_buf));
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
      if (GMMetrics_ccc_duo_get_from_index_2_threshold<COUNTED_BITS_PER_ELT>(
            metrics, index, env)) {
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            const GMFloat value =
              GMMetrics_ccc_duo_get_from_index_2<COUNTED_BITS_PER_ELT>(
                metrics, index, i0, i1, env);
            if (env->pass_threshold(value)) {
              const size_t coord0 =
                GMMetrics_coord_global_from_index(metrics, index, 0, env);
              const size_t coord1 =
                GMMetrics_coord_global_from_index(metrics, index, 1, env);
              const char do_out = coord0 < metrics->num_vector_active &&
                                  coord1 < metrics->num_vector_active;
              const size_t ind_buf = i1 + 2*(i0 + 2*(index-ind_base));
              COMET_ASSERT(ind_buf < num_buf);
              do_out_buf[ind_buf] = do_out;
              coord0_buf[ind_buf] = coord0;
              coord1_buf[ind_buf] = coord1;
              i01_buf[ind_buf] = i0 + 2*i1;
              value_buf[ind_buf] = value;
            } /*---if---*/
          } /*---i1---*/
        } /*---i0---*/
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
            const int i0 = i01_buf[ind_buf] % 2;
            const int i1 = i01_buf[ind_buf] / 2;
            const size_t coord0 = coord0_buf[ind_buf];
            const size_t coord1 = coord1_buf[ind_buf];
            writer.write(coord0, coord1, i0, i1, value_buf[ind_buf]);
            // Reset buffer entry to false
            do_out_buf[ind_buf] = 0;
          }
        }
      }
    } /*---ind_buf---*/

  } /*---ind_base---*/

  num_written_ += writer.num_written();

  free(do_out_buf);
  free(coord0_buf);
  free(coord1_buf);
  free(i01_buf);
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

  //enum {vals_per_index = 8;}

  // Number of index values values visited for one pass across buffer
  const size_t num_buf_ind = 1000 * 1000;
  // Number of (index, i0, i1) entries (potentially) stored in buffer
  const size_t num_buf = 8 * num_buf_ind;

  // Each buffer entry contains: whether value is to be written,
  // coord0, coord1, coord2, i0, i1, i2, and value

  char* const do_out_buf = (char*)malloc(num_buf*sizeof(*do_out_buf));
  int* const coord0_buf = (int*)malloc(num_buf*sizeof(*coord0_buf));
  int* const coord1_buf = (int*)malloc(num_buf*sizeof(*coord1_buf));
  int* const coord2_buf = (int*)malloc(num_buf*sizeof(*coord2_buf));
  int* const i012_buf = (int*)malloc(num_buf*sizeof(*i012_buf));
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
      if (GMMetrics_ccc_duo_get_from_index_3_threshold<COUNTED_BITS_PER_ELT>(
             metrics, index, env)) {
        for (int i0 = 0; i0 < 2; ++i0) {
          for (int i1 = 0; i1 < 2; ++i1) {
            for (int i2 = 0; i2 < 2; ++i2) {
              const GMFloat value =
                GMMetrics_ccc_duo_get_from_index_3<COUNTED_BITS_PER_ELT>(
                  metrics, index, i0, i1, i2, env);
              if (env->pass_threshold(value)) {
                const size_t coord0 =
                  GMMetrics_coord_global_from_index(metrics, index, 0, env);
                const size_t coord1 =
                  GMMetrics_coord_global_from_index(metrics, index, 1, env);
                const size_t coord2 =
                  GMMetrics_coord_global_from_index(metrics, index, 2, env);
                const char do_out = coord0 < metrics->num_vector_active &&
                                    coord1 < metrics->num_vector_active &&
                                    coord2 < metrics->num_vector_active;
                const size_t ind_buf = i1 + 2*(i0 + 2*(i2 +2*(index-ind_base)));
                do_out_buf[ind_buf] = do_out;
                coord0_buf[ind_buf] = coord0;
                coord1_buf[ind_buf] = coord1;
                coord2_buf[ind_buf] = coord2;
                i012_buf[ind_buf] = i0 + 2*(i1 + 2*i2);
                value_buf[ind_buf] = value;
              } /*---if---*/
            } /*---i2---*/
          } /*---i1---*/
        } /*---i0---*/
      }
    } /*---for index---*/

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
            const int i0 = i012_buf[ind_buf] % 2;
            const int i1 = (i012_buf[ind_buf] / 2) % 2;
            const int i2 = i012_buf[ind_buf] / 4;
            const size_t coord0 = coord0_buf[ind_buf];
            const size_t coord1 = coord1_buf[ind_buf];
            const size_t coord2 = coord2_buf[ind_buf];
            writer.write(coord0, coord1, coord2, i0,i1,i2, value_buf[ind_buf]);
            // Reset buffer entry to false
            do_out_buf[ind_buf] = 0;
          }
        }
      }
    } /*---ind_buf---*/
  } /*---ind_base---*/

  num_written_ += writer.num_written();

  free(do_out_buf);
  free(coord0_buf);
  free(coord1_buf);
  free(coord2_buf);
  free(i012_buf);
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
      const size_t coord0 =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t coord1 =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      if (coord0 >= metrics->num_vector_active ||
          coord1 >= metrics->num_vector_active)
        continue;
      const auto value = Metrics_get<GMFloat>(*metrics, index, *env);
      if (!env->pass_threshold(value))
        continue;
      /// Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li): value: %.17e\n",
          coord0, coord1, value);
      else
        writer.write(coord0, coord1, value);
    } // for index

    num_written_ += writer.num_written();

  //----------
  } else if (env->data_type_metrics() == GM_DATA_TYPE_FLOAT &&
           env->num_way() == NUM_WAY::_3) {
  //----------

    MetricIO writer(file, *metrics, *env);

    for (size_t index = 0; index < metrics->num_elts_local; ++index) {
      const size_t coord0 =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t coord1 =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      const size_t coord2 =
        GMMetrics_coord_global_from_index(metrics, index, 2, env);
      if (coord0 >= metrics->num_vector_active ||
          coord1 >= metrics->num_vector_active ||
          coord2 >= metrics->num_vector_active)
        continue;
      const auto value = Metrics_get<GMFloat>(*metrics, index, *env);
      if (!env->pass_threshold(value))
        continue;
      // Output the value.
      if (stdout == file)
        fprintf(file, "element (%li,%li,%li): value: %.17e\n",
          coord0, coord1, coord2, value);
      else
        writer.write(coord0, coord1, coord2, value);
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
      const size_t coord0 =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t coord1 =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      if (coord0 >= metrics->num_vector_active ||
          coord1 >= metrics->num_vector_active)
        continue;
      int num_out_this_line = 0;
      for (int i0 = 0; i0 < 2; ++i0) {
        for (int i1 = 0; i1 < 2; ++i1) {
          const GMFloat value = env->metric_type() == MetricType::CCC ?
            GMMetrics_ccc_duo_get_from_index_2<CBPE::CCC>(metrics, index, i0, i1, env) :
            GMMetrics_ccc_duo_get_from_index_2<CBPE::DUO>(metrics, index, i0, i1, env);
          if (!env->pass_threshold(value))
            continue;

          // Output the value.

          if (num_out_this_line == 0)
            fprintf(file, "element (%li,%li): values:", coord0, coord1);

          fprintf(file, " %i %i %.17e", i0, i1, value);

          num_out_this_line++;

        } // i1
      } // i0
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
      const size_t coord0 =
        GMMetrics_coord_global_from_index(metrics, index, 0, env);
      const size_t coord1 =
        GMMetrics_coord_global_from_index(metrics, index, 1, env);
      const size_t coord2 =
        GMMetrics_coord_global_from_index(metrics, index, 2, env);
      if (coord0 >= metrics->num_vector_active ||
          coord1 >= metrics->num_vector_active ||
          coord2 >= metrics->num_vector_active)
        continue;
      int num_out_this_line = 0;
      for (int i0 = 0; i0 < 2; ++i0) {
        for (int i1 = 0; i1 < 2; ++i1) {
          for (int i2 = 0; i2 < 2; ++i2) {
            const GMFloat value = env->metric_type() == MetricType::CCC ?
              GMMetrics_ccc_duo_get_from_index_3<CBPE::CCC>(metrics, index, i0, i1, i2, env) :
              GMMetrics_ccc_duo_get_from_index_3<CBPE::DUO>(metrics, index, i0, i1, i2, env);
            if (!env->pass_threshold(value))
              continue;

            // Output the value.

            if (num_out_this_line == 0)
              fprintf(file, "element (%li,%li,%li): values:",
                coord0, coord1, coord2);

            fprintf(file, " %i %i %i %.17e", i0, i1, i2, value);

            num_out_this_line++;

          } // i2
        } // i1
      } // i0
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

  if (! file_)
    return;

  // CLose file.

  fclose(file_);

  // Reopen file for read.

  MetricsIO::open(path_stub_.c_str(), env_, "rb");

  // Set fpos to beginning of last metrics object write.

  long offset = bytes_(num_written_ - num_written_last_write_);
  int success = fseek(file_, offset, SEEK_SET);
  COMET_INSIST(0 == success);

  // Loop over metrics stored in file.

  GMDecompMgr* const dm = metrics.dm;
  size_t num_incorrect = 0;

  for (size_t i = 0; i < num_written_last_write_; ++i) {

    if (env_.num_way() == NUM_WAY::_2) {

      MetricIO::Metric<NUM_WAY::_2> metric;
      MetricIO::read(metric, file_, env_);

      if (env_.metric_type() == MetricType::CZEK) {

        const size_t i_vag = metric.coord0();
        const size_t j_vag = metric.coord1();
        const size_t i = GMDecompMgr_get_vector_local_from_vector_active(
          dm, i_vag, &env_);
        const size_t j = GMDecompMgr_get_vector_local_from_vector_active(
          dm, j_vag, &env_);
        COMET_ASSERT(i >= 0 && i < dm->num_vector_active);
        COMET_ASSERT(j >= 0 && j < dm->num_vector_active);

        const int i_proc = GMDecompMgr_get_proc_vector_from_vector_active(
          dm, i_vag, &env_);
        const int j_proc = GMDecompMgr_get_proc_vector_from_vector_active(
          dm, j_vag, &env_);
        no_unused_variable_warning(i_proc);
        COMET_ASSERT(env_.proc_num_vector() == i_proc);
        COMET_ASSERT(j_proc >= 0 && j_proc < env_.num_proc_vector());

        const MetricIO::Float_t metric_value = env_.all2all() ?
          GMMetrics_float_get_all2all_2(&metrics, i, j, j_proc, &env_) :
          GMMetrics_float_get_2(&metrics, i, j, &env_);

        const bool is_incorrect = metric_value != metric.value;
        num_incorrect += is_incorrect;

        if (is_incorrect && num_incorrect < 10)
          fprintf(stderr, "Incorrect metric value: "
            "element %zu %zu actual %.17e expected %.17e\n",
            i_vag, j_vag, (double)metric_value, (double)metric.value);

      } else { // if (env_.is_metric_type_bitwise())
      

// TODO !!!


      } //  if (env_.metric_type() == MetricType::CZEK)

    } else { // if (env_.num_way() == NUM_WAY::_3)

      MetricIO::Metric<NUM_WAY::_3> metric;
      MetricIO::read(metric, file_, env_);

      if (env_.metric_type() == MetricType::CZEK) {


// TODO !!!

// NOTE: need to count number of metrics that pass threshold, compare to file
// NOTE: to cleanup, maybe move stuff to metrics_?way_*.hh
// NOTE: check to see if i, j, k may be permuted in file.
// NOTE: need more tests in tester: multi proc, phase, stage, all 4-6 methods.

      } else { // if (env_.is_metric_type_bitwise())
      

// TODO !!!


      } //  if (env_.metric_type() == MetricType::CZEK)

    } // if (env_.num_way() == NUM_WAY::_2)

  } // i

  COMET_INSIST(0 == num_incorrect &&
    "Incorrect values detected in output file - partial list is displayed.");

  // Close file.

  fclose(file_);

  // Reopen file for write.

  MetricsIO::open(path_stub_.c_str(), env_, "ab");

  // Set fpos to end of file.

  offset = bytes_(num_written_);
  success = fseek(file_, offset, SEEK_SET);
  COMET_INSIST(0 == success);
}

//-----------------------------------------------------------------------------
/// \brief Static function to open (one-per-rank set of) metrics files.

FILE* MetricsIO::open(const char* path_stub, CEnv& env, const char* mode) {

  // NOTE: allowing flexibility to open on all ranks, not just active ranks.

  // Form filename

  size_t len = strlen(path_stub);
  char* path = (char*)malloc((len+50) * sizeof(char));

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
  free(path);

  return file;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
