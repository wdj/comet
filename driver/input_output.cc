//-----------------------------------------------------------------------------
/*!
 * \file   input_output.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdio.h"
//#include "stdlib.h"
//#include "stddef.h"
#include "string.h"
//#include "float.h"
//#include "errno.h"

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "driver.hh"
#include "input_output.hh"

//=============================================================================
// Input vectors from files

void set_vectors_from_file_float(GMVectors* vectors, DriverOptions* do_,
                                 GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  FILE* input_file = fopen(do_->input_file_path, "r");
  GMInsist(NULL != input_file && "Unable to open input file.");

  const int fl = 0;
  const size_t field_base = fl +
    vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);

  for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

    const size_t proc_num = GMEnv_proc_num_vector_i(env);
    const size_t vector = vl + vectors->num_vector_local * proc_num;
    //---shuffle.
    //const size_t vector = proc_num + GMEnv_num_proc_vector_i(env) * vl;
    /*---Fill pad vectors with copies of the last vector---*/
    const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                 vector : do_->num_vector_active-1;

    const size_t elt_num = field_base + vectors->num_field * vector_capped;
    const size_t addr_file = elt_num * sizeof(GMFloat);
    int fseek_success = fseek(input_file, addr_file, SEEK_SET);
    fseek_success += 0; /*---Avoid unused var warning---*/
    GMInsist(0 == fseek_success && "File seek failure.");
    GMFloat* const addr_mem = GMVectors_float_ptr(vectors, fl, vl, env);
    /*---NOTE: the following call is ok since has no side effects---*/
    GMInsist((fl+1 >= vectors->num_field_local ||
        GMVectors_float_ptr(vectors, fl+1, vl, env) == addr_mem + 1)
        && "Vector layout is incompatible with operation.");

    size_t num_read = fread(addr_mem, sizeof(GMFloat),
                            vectors->num_field_local, input_file);
    num_read += 0; /*---Avoid unused var warning---*/
    GMInsist((size_t)vectors->num_field_local == (size_t)num_read &&
             "File read failure.");

  } /*---vl---*/

  fclose(input_file);
}

//-----------------------------------------------------------------------------

void set_vectors_from_file_bits2(GMVectors* vectors, DriverOptions* do_,
                                 GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMInsistInterface(env, GMEnv_num_proc_field(env) == 1 &&
                    "CCC file read for this case not yet implemented.");
  const int pvfl = 0;

  typedef char input_t;

  FILE* input_file = fopen(do_->input_file_path, "r");
  GMInsist(NULL != input_file && "Unable to open input file.");

  for (int vl = 0; vl < vectors->num_vector_local; ++vl) {

    const size_t proc_num = GMEnv_proc_num_vector_i(env);
    const size_t vector = vl + vectors->num_vector_local * proc_num;
    /*---Fill pad vectors with copies of the last vector---*/
    const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                 vector : do_->num_vector_active-1;

    const int bits_per_field = 2;
    const int bits_per_byte = 8;
    const size_t bytes_per_vector
      = gm_ceil_i8(vectors->num_field * bits_per_field, bits_per_byte);
    const size_t addr_file = bytes_per_vector * vector_capped;

    int fseek_success = fseek(input_file, addr_file, SEEK_SET);
    fseek_success += 0; /*---Avoid unused var warning---*/
    GMInsist(0 == fseek_success && "File seek failure.");

    input_t* const addr_mem
       = (input_t*)GMVectors_bits2x64_ptr(vectors, pvfl, vl, env);

    size_t num_read = fread(addr_mem, sizeof(input_t),
                            bytes_per_vector, input_file);
    num_read += 0; /*---Avoid unused var warning---*/
    GMInsist(bytes_per_vector == (size_t)num_read && "File read failure.");

  } /*---vl---*/

  fclose(input_file);
}

//-----------------------------------------------------------------------------

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMInsist(vectors && do_ && env && do_->input_file_path);

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      set_vectors_from_file_float(vectors, do_, env);

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

      set_vectors_from_file_bits2(vectors, do_, env);

    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================

void write_vectors_to_file(GMVectors* vectors, char* vectors_file_path, 
                           GMEnv* env) {
  GMInsist(vectors && vectors_file_path && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  FILE* vectors_file = fopen(vectors_file_path, "w");

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      const size_t nval = vectors->dm->num_vector_active_local;
      const size_t nfal = vectors->dm->num_field_active_local;

      size_t num_written_total = 0;
      const size_t num_written_attempted_total = nval * nfal;

      for (size_t vl = 0 ; vl < nval; ++vl) {
        for (size_t fl = 0 ; fl < nfal; ++fl) {

          const GMFloat outv = GMVectors_float_get(vectors, fl, vl, env);
          const size_t num_written = fwrite(&outv, sizeof(outv), 1,
                                            vectors_file);
          num_written_total += num_written;
        }
      }

      GMInsist(num_written_attempted_total == num_written_total &&
               "File write failure.");

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/

      const size_t nval = vectors->dm->num_vector_active_local;
      const size_t npfl = vectors->dm->num_packedfield_local;

      size_t num_written_total = 0;
      size_t num_written_attempted_total = 0;

      for (size_t vl = 0 ; vl < nval; ++vl) {
        for (size_t pfl = 0 ; pfl < npfl; ++pfl) {

          GMBits2x64 outv = GMVectors_bits2x64_get(vectors, pfl, vl, env);

          if (pfl == npfl-1) {
            const size_t num_written = fwrite(&outv, sizeof(outv), 1,
                                              vectors_file);
            num_written_total += num_written;
            num_written_attempted_total += sizeof(outv);
            continue;
          }



//FIX




        }
      }

      GMInsist(num_written_attempted_total == num_written_total &&
               "File write failure.");

    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/

  fclose(vectors_file);
}

//=============================================================================
// Output results metrics to file: implementation

void output_metrics_impl(GMMetrics* metrics, FILE* file,
                         double threshold, size_t& num_written, GMEnv* env) {
  GMInsist(metrics && file && env);

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  MetricWriter writer(file, metrics, env);

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
      /*----------*/
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li): value: %.17e\n" :
              "element (%li,%li): value: %.8e\n", coord0, coord1, value);

          } else {

            writer.write(coord0, coord1, value);

          }
        }
      }

      /*----------*/
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
      /*----------*/
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t coord2 =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (coord0 >= metrics->num_vector_active ||
              coord1 >= metrics->num_vector_active ||
              coord2 >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czek_get_from_index(metrics, index, env);
          if (!(threshold < 0. || value > threshold)) {
            continue;
          }
          /*---Output the value---*/
          if (file == stdout) {

            fprintf(file,
              sizeof(GMFloat) == 8 ?
              "element (%li,%li,%li): value: %.17e\n" :
              "element (%li,%li,%li): value: %.8e\n",
              coord0, coord1, coord2, value);

          } else {

            writer.write(coord0, coord1, coord2, value);

          }
        }
      }

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/

      GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_2);

      if (file != stdout) {

        const size_t num_buf_ind = 1000 * 1000;
        const size_t num_buf = 4 * num_buf_ind;

        char do_out_buf[num_buf];
        int coord0_buf[num_buf];
        int coord1_buf[num_buf];
        int i01_buf[num_buf];
        GMFloat value_buf[num_buf];

        for (int i=0; i<(int)num_buf; ++i) {
          do_out_buf[i] = 0;
        }

        const GMFloat threshold_eff = threshold<0. ? -1e20 : threshold;
  
        for (size_t ind_base = 0; ind_base < metrics->num_elts_local;
             ind_base += num_buf_ind) {
          const size_t ind_max = gm_min_i8(metrics->num_elts_local,
                                           ind_base + num_buf_ind);
#pragma omp parallel for schedule(dynamic,1000)
          for (size_t index = ind_base; index < ind_max; ++index) {
            if (GMMetrics_ccc_get_from_index_2_threshold(
                   metrics, index, threshold_eff, env)) {
              for (int i0 = 0; i0 < 2; ++i0) {
                for (int i1 = 0; i1 < 2; ++i1) {
                  const GMFloat value =
                    GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
                  const bool do_thresh = value>threshold_eff;
                  if (do_thresh) {
                    const size_t coord0 =
                      GMMetrics_coord_global_from_index(metrics, index, 0, env);
                    const size_t coord1 =
                      GMMetrics_coord_global_from_index(metrics, index, 1, env);
                    const char do_out = coord0 < metrics->num_vector_active &&
                                        coord1 < metrics->num_vector_active;
                    const size_t ind_buf = i1 + 2*(i0 + 2*(index-ind_base));
                    do_out_buf[ind_buf] = do_out;
                    coord0_buf[ind_buf] = coord0;
                    coord1_buf[ind_buf] = coord1;
                    i01_buf[ind_buf] = i0 + 2*i1;
                    value_buf[ind_buf] = value;
                  } /*---if---*/
                } /*---i1---*/
              } /*---i0---*/
            }
          } /*---for index---*/

          //assert(sizeof(int) == 4 * sizeof(char));
          const size_t ind_buf_max = 4 * (ind_max - ind_base);
          const int* const do_out_ptr_max = (int*)(do_out_buf + ind_buf_max);
          int* do_out_ptr = (int*)do_out_buf;
          for (; do_out_ptr < do_out_ptr_max;) {
            if (*(do_out_ptr++)) {
              for (int i=0; i<4; ++i) {
                const size_t ind_buf = (do_out_ptr - (int*)do_out_buf - 1)*4 + i;
                if (do_out_buf[ind_buf]) {
                  const int i0 = i01_buf[ind_buf] % 2;
                  const int i1 = i01_buf[ind_buf] / 2;
                  const size_t coord0 = coord0_buf[ind_buf];
                  const size_t coord1 = coord1_buf[ind_buf];
                  writer.write(coord0, coord1, i0, i1, value_buf[ind_buf]);
                  do_out_buf[ind_buf] = 0;
                }
              }
            }
          } /*---ind_buf---*/
        } /*---ind_base---*/

      } else /*---stdout---*/ {

        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          const size_t coord0 =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          if (coord0 >= metrics->num_vector_active) {
            continue;
          }
          const size_t coord1 =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (coord1 >= metrics->num_vector_active) {
            continue;
          }
          int num_out_this_line = 0;
          int i0;
          for (i0 = 0; i0 < 2; ++i0) {
            int i1;
            for (i1 = 0; i1 < 2; ++i1) {
              const GMFloat value
                = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
              if (!(threshold < 0. || value > threshold)) {
                continue;
              }

              /*---Output the value---*/
              if (file == stdout) {

                if (num_out_this_line == 0) {
                  fprintf(file,
                    "element (%li,%li): values:", coord0, coord1);
                }

                fprintf(file,
                  sizeof(GMFloat) == 8 ?
                  " %i %i %.17e" :
                  " %i %i %.8e", i0, i1, value);

              } else {

                writer.write(coord0, coord1, i0, i1, value);

              }

              num_out_this_line++;

            } /*---i1---*/
          } /*---i0---*/
          if (file == stdout) {
            if (num_out_this_line > 0) {
              fprintf(file, "\n");
            }
          }
        } /*---index---*/

      } /*---if---*/

    } break;

    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/

      size_t index = 0;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        const size_t coord0 =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t coord1 =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t coord2 =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (coord0 >= metrics->num_vector_active ||
            coord1 >= metrics->num_vector_active ||
            coord2 >= metrics->num_vector_active) {
          continue;
        }
        int num_out_this_line = 0;
        int i0;
        for (i0 = 0; i0 < 2; ++i0) {
          int i1;
          for (i1 = 0; i1 < 2; ++i1) {
            int i2;
            for (i2 = 0; i2 < 2; ++i2) {
              const GMFloat value
                = GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2,
                                                env);
              if (!(threshold < 0. || value > threshold)) {
                continue;
              }

              /*---Output the value---*/
              if (file == stdout) {

                if (num_out_this_line == 0) {
                  fprintf(file,
                    "element (%li,%li,%li): values:", coord0, coord1, coord2);
                }

                fprintf(file,
                  sizeof(GMFloat) == 8 ?
                  " %i %i %i %.17e" :
                  " %i %i %i %.8e", i0, i1, i2, value);

              } else {

                writer.write(coord0, coord1, coord2, i0, i1, i2, value);

              }

              num_out_this_line++;

            } /*---i2---*/
          } /*---i1---*/
        } /*---i0---*/
        if (file == stdout) {
          if (num_out_this_line > 0) {
            fprintf(file, "\n");
          }
        }
      } /*---index---*/

    } break;
    /*--------------------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/

  num_written += writer.get_num_written();
}

//=============================================================================
// MetricWriter member definitions.

MetricWriter::MetricWriter(FILE* file, GMMetrics* metrics, GMEnv* env) :
  file_(file),
  data_type_(GMEnv_data_type_metrics(env)),
  num_way_(GMEnv_num_way(env)),
  num_written_total_(0) {

  if (stdout != file_) {
    GMInsistInterface(env, metrics->num_vector_active ==
                  (GMUInt32)metrics->num_vector_active &&
                  "Too many vectors for output format.");
  }
}

//-----------------------------------------------------------------------------
// CZEK 2-way

void
MetricWriter::write(size_t coord0, size_t coord1, GMFloat value) {

  bool success = true;

  const GMUInt32 outc0 = coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc1 = coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
  }

//-----------------------------------------------------------------------------
// CZEK 3-way

void
MetricWriter::write(size_t coord0, size_t coord1, size_t coord2,
                     GMFloat value) {

  bool success = true;

  const GMUInt32 outc0 = coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc1 = coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc2 = coord2;
  num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//-----------------------------------------------------------------------------
// CCC 2-way

void
MetricWriter::write(size_t coord0, size_t coord1, int i0, int i1,
                     GMFloat value) {
  GMAssert(i0 >=0 && i0 < 2);
  GMAssert(i1 >=0 && i1 < 2);

  bool success = true;

  const GMUInt32 outc0 = i0 + 2 * coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc1 = i1 + 2 * coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
  num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//-----------------------------------------------------------------------------
// CCC 3-way

void
MetricWriter::write(size_t coord0, size_t coord1, size_t coord2,
              int i0, int i1, int i2, GMFloat value) {
  GMAssert(i0 >=0 && i0 < 2);
  GMAssert(i1 >=0 && i1 < 2);
  GMAssert(i2 >=0 && i2 < 2);

  bool success = true;

  const GMUInt32 outc0 = i0 + 2 * coord0;
  size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc1 = i1 + 2 * coord1;
  num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
  success = success && num_written == 1;

  const GMUInt32 outc2 = i2 + 2 * coord2;
  num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
  success = success && num_written == 1;

  const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
  success = success && num_written == 1;

  num_written_total_ += success ? 1 : 0;
  GMInsist(success && "File write failure.");
}

//=============================================================================
// MetricsFile member definitions.

MetricsFile::MetricsFile(DriverOptions* do_, GMEnv* env)
  : file_(NULL)
  , num_written_(0) {
  GMInsist(do_ && env);

  verbosity_ = do_->verbosity;
  threshold_ = do_->threshold;

  char* stub = do_->metrics_file_path_stub;
  if (! (NULL != stub && GMEnv_is_proc_active(env) &&
      GMEnv_proc_num_field(env) == 0) ) {
    return;
  }

  // Form filename

  size_t len = strlen(stub);
  char* path = (char*)malloc((len+50) * sizeof(char));

  int num_digits = 0;
  for (int tmp = 1; ; tmp*=10, ++num_digits) {
    if (tmp > GMEnv_num_proc(env)) {
      break;
    }
  }

  char format[100];
  sprintf(format, "%s0%ii.bin", "%s_%", num_digits);

  sprintf(path, format, stub, GMEnv_proc_num(env));

  /*---Do open---*/

  file_ = fopen(path, "w");
  free(path);
}

//-----------------------------------------------------------------------------

MetricsFile::~MetricsFile() {
  if (file_) {
    fclose(file_);
  }
}

//-----------------------------------------------------------------------------

void MetricsFile::write(GMMetrics* metrics, GMEnv* env) {

  // Output to file

  if (file_) {
    output_metrics_impl(metrics, file_, threshold_, num_written_, env);
  }

  // Output to stdout if requested

  if (verbosity_ > 1) {
    double threshold = verbosity_ > 2 ? -1. : threshold_;
    output_metrics_impl(metrics, stdout, threshold, num_written_, env);
  }
}

//=============================================================================
