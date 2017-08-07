/*---------------------------------------------------------------------------*/
/*!
 * \file   driver.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>
#include <errno.h>

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "checksums.hh"
#include "compute_metrics.hh"
#include "driver.hh"
#include "test_problems.hh"

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

void finish_parsing(int argc, char** argv, DriverOptions* do_, GMEnv* env) {
  errno = 0;
  int i = 0;
  for (i = 1; i < argc; ++i) {
    /*----------*/
    if (strcmp(argv[i], "--num_field") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_field.");
      const long num_field = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_field >= 0
                    && "Invalid setting for num_field.");
      do_->num_field_active = num_field;
      do_->num_field_active_initialized = true;
      do_->num_field_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_field_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_field_local.");
      const long num_field_local = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && num_field_local >= 0 &&
                    (long)(int)num_field_local == num_field_local &&
                    "Invalid setting for num_field_local.");
      do_->num_field_local = num_field_local;
      do_->num_field_local_initialized = true;
      do_->num_field_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_vector.");
      const long num_vector = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_vector >= 0
                    && "Invalid setting for num_vector.");
      do_->num_vector_active = num_vector;
      do_->num_vector_active_initialized = true;
      do_->num_vector_local_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_vector_local.");
      const long num_vector_local = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && num_vector_local >= 0 &&
                    (long)(int)num_vector_local == num_vector_local &&
                    "Invalid setting for num_vector_local.");
      do_->num_vector_local = num_vector_local;
      do_->num_vector_local_initialized = true;
      do_->num_vector_active_initialized = false;
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for verbosity.");
      const long verbosity = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && verbosity >= 0 &&
                    "Invalid setting for verbosity.");
      do_->verbosity = verbosity;
      /*--------------------*/
    } else if (strcmp(argv[i], "--checksum") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for checksum." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        do_->checksum = true;
      } else if (strcmp(argv[i], "no") == 0) {
        do_->checksum = false;
      } else {
        GMInsist(env, false ? "Invalid setting for checksum." : 0);
      }
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_stage.");
      const long num_stage = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_stage >= 1
                    && (long)(int)num_stage == num_stage
                    && "Invalid setting for num_stage.");
      env->num_stage = num_stage;
      do_->stage_min_1based = 1;
      do_->stage_max_1based = env->num_stage;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for stage_min.");
      const long stage_min_1based = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && stage_min_1based >= 1
                    && (long)(int)stage_min_1based == stage_min_1based
                    && "Invalid setting for stage_min.");
      do_->stage_min_1based = stage_min_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for stage_max.");
      const long stage_max_1based = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && stage_max_1based <= env->num_stage
                    && (long)(int)stage_max_1based == stage_max_1based
                    && "Invalid setting for stage_max.");
      do_->stage_max_1based = stage_max_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--num_phase") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_phase.");
      const long num_phase = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_phase >= 1
                    && (long)(int)num_phase == num_phase
                    && "Invalid setting for num_phase.");
      env->num_phase = num_phase;
      do_->phase_min_1based = 1;
      do_->phase_max_1based = env->num_phase;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_min") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for phase_min.");
      const long phase_min_1based = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && phase_min_1based >= 1
                    && (long)(int)phase_min_1based == phase_min_1based
                    && "Invalid setting for phase_min.");
      do_->phase_min_1based = phase_min_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--phase_max") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for phase_max.");
      const long phase_max_1based = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && phase_max_1based <= env->num_phase
                    && (long)(int)phase_max_1based == phase_max_1based
                    && "Invalid setting for phase_max.");
      do_->phase_max_1based = phase_max_1based;
    /*----------*/
    } else if (strcmp(argv[i], "--input_file") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for input_file." : 0);
      do_->input_file_path = argv[i];
    /*----------*/
    } else if (strcmp(argv[i], "--output_file_stub") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for output_file_stub." : 0);
      do_->output_file_path_stub = argv[i];
      /*--------------------*/
    } else if (strcmp(argv[i], "--problem_type") == 0) {
      /*--------------------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for problem_type." : 0);
      if (strcmp(argv[i], "random") == 0) {
        GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_RANDOM);
      } else if (strcmp(argv[i], "analytic") == 0) {
        GMEnv_set_compute_method(env, GM_PROBLEM_TYPE_ANALYTIC);
      } else {
        GMInsist(env, false ? "Invalid setting for problem_type." : 0);
      }
    /*----------*/
    } else if (strcmp(argv[i], "--threshold") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for threshold." : 0);
      errno = 0;
      const double threshold = strtod(argv[i], NULL);
      GMInsist(env, 0 == errno && "Invalid setting for ccc_param.");
      do_->threshold = threshold;
     /*----------*/
    } else if (strcmp(argv[i], "--metric_type") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--ccc_param") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else if (strcmp(argv[i], "--sparse") == 0) {
      ++i; /*---processed elsewhere by GMEnv---*/
    } else {
    /*----------*/
      if (GMEnv_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      GMInsist(env, false ? "Error: argument not recognized." : 0);
    /*----------*/
    } /*---if/else---*/

  } /*---for i---*/

  GMInsist(env, do_->num_field_local_initialized ||
                do_->num_field_active_initialized
                ? "Error: must set num_field_local or num_field."
                : 0);
  GMInsist(env, do_->num_vector_local_initialized ||
                do_->num_vector_active_initialized
                ? "Error: must set num_vector_local or num_vector."
                : 0);
}

/*---------------------------------------------------------------------------*/

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors && do_ && env);
  GMAssertAlways(do_->input_file_path);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      const int fl = 0;
      const size_t field_base = fl +
        vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
      FILE* input_file = fopen(do_->input_file_path, "r");
      GMAssertAlways(NULL != input_file ? "Unable to open input file." : 0);
      int vl = 0;
      for (vl = 0; vl < vectors->num_vector_local; ++vl) {
        const size_t proc_num = GMEnv_proc_num_vector_i(env);
        const size_t vector = vl + vectors->num_vector_local * proc_num;
        //---shuffle.
        //const size_t vector = proc_num + GMEnv_num_proc_vector_i(env) * vl;
        /*---Fill pad vectors with copies of the last vector---*/
        const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                     vector : do_->num_vector_active-1;
        const size_t addr_file =
          (field_base + vectors->num_field * vector_capped) * sizeof(GMFloat);
        int fseek_success = fseek(input_file, addr_file, SEEK_SET);
        fseek_success += 0; /*---Avoid unused var warning---*/
        GMAssertAlways(0 == fseek_success);
        GMFloat* const addr_mem = GMVectors_float_ptr(vectors, fl, vl, env);
        /*---NOTE: the following call is ok since has no side effects---*/
        GMAssertAlways(fl+1 >= vectors->num_field_local ||
            GMVectors_float_ptr(vectors, fl+1, vl, env) == addr_mem + 1
            ? "Vector layout is incompatible with operation." : 0);
        size_t num_read = fread(addr_mem, sizeof(GMFloat),
                                vectors->num_field_local, input_file);
        num_read += 0; /*---Avoid unused var warning---*/
        GMAssertAlways((size_t)vectors->num_field_local == (size_t)num_read);
      } /*---vl---*/
      fclose(input_file);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
      GMInsist(env, false ? "Not yet implemented." : 0);
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(false ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*---------------------------------------------------------------------------*/

void set_vectors(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors && do_ && env);

  if (do_->input_file_path != NULL) {
    set_vectors_from_file(vectors, do_, env);
  } else if (do_->problem_type == GM_PROBLEM_TYPE_RANDOM) {
    set_vectors_random(vectors, do_, env);
  } else {
    set_vectors_analytic(vectors, do_, env);
  }
}

/*===========================================================================*/
/*---Output the result metrics values to file---*/

class Writer {

  FILE* file_;
  const int data_type_;
  int num_way_;
  GMEnv* env_;

  size_t num_written_total_;

  //typedef unsigned char out_t;
  //const GMFloat out_max = 255;
  //const int buf_size = 4096;
  //int buf_elts = 0;
  //out_t buf[buf_size];

public:

  //---------------------------------------------------------------------------

  Writer(FILE* file, GMMetrics* metrics, GMEnv* env) :
    file_(file),
    data_type_(GMEnv_data_type_metrics(env)),
    num_way_(GMEnv_num_way(env)) {

    if (file_ != stdout) {
      GMInsist(env, metrics->num_vector_active ==
                    (GMUInt32)metrics->num_vector_active &&
                    "Too many vectors for output format.");
    }

  }

  //---------------------------------------------------------------------------

  ~Writer() {

      //if (file != stdout) {
        //size_t num_written = fwrite(&buf, sizeof(out_t), buf_elts, file);
        //num_written_total += num_written;
        //GMAssert(num_written == buf_elts*sizeof(out_t));
        //buf_elts = 0;
        //printf("Wrote %lu elements of %lu from proc %i.\n",
        //       num_written_total, metrics->num_elts_local, GMEnv_proc_num(env));
      //}

  }

  //---------------------------------------------------------------------------
  // CZEK 2-way

  void write(size_t coord0, size_t coord1, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

    //out_t out_v = (out_t)(value * out_max);
    //buf[buf_elts++] = out_v;
    //if (buf_elts == buf_size) {
    //  size_t num_written = fwrite(&buf, sizeof(out_t),
    //                              buf_elts, file);
    //  num_written_total += num_written;
    //  GMAssert(num_written == buf_elts*sizeof(out_t));
    //  buf_elts = 0;
    //}
  }

  //---------------------------------------------------------------------------
  // CZEK 3-way

  void write(size_t coord0, size_t coord1, size_t coord2, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMUInt32 outc2 = coord2;
    num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
    success = success && num_written == 1*sizeof(outc2);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;
  }

  //---------------------------------------------------------------------------
  // CCC 2-way

  void write(size_t coord0, size_t coord1, int i0, int i1, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = i0 + 2 * coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = i1 + 2 * coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

  }

  //---------------------------------------------------------------------------
  // CCC 3-way

  void write(size_t coord0, size_t coord1, size_t coord2,
             int i0, int i1, int i2, GMFloat value) {

    bool success = true;

    const GMUInt32 outc0 = i0 + 2 * coord0;
    size_t num_written = fwrite(&outc0, sizeof(outc0), 1, file_);
    success = success && num_written == 1*sizeof(outc0);

    const GMUInt32 outc1 = i1 + 2 * coord1;
    num_written = fwrite(&outc1, sizeof(outc1), 1, file_);
    success = success && num_written == 1*sizeof(outc1);

    const GMUInt32 outc2 = i2 + 2 * coord2;
    num_written = fwrite(&outc2, sizeof(outc2), 1, file_);
    success = success && num_written == 1*sizeof(outc2);

    const GMFp32 outv = value;
    num_written = fwrite(&outv, sizeof(outv), 1, file_);
    success = success && num_written == 1*sizeof(outv);

    num_written_total_ += success ? 1 : 0;

  }

  //---------------------------------------------------------------------------
};

/*---------------------------------------------------------------------------*/

void output_metrics_impl(GMMetrics* metrics, DriverOptions* do_,
                         FILE* file, double threshold, GMEnv* env) {
  GMAssertAlways(metrics && do_ && file && env);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  Writer writer(file, metrics, env);

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

      GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

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
      GMAssertAlways(false ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

void output_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(metrics && do_ && env);

  char* stub = do_->output_file_path_stub;

  /*---Output to file if requested---*/

  if (NULL != stub && GMEnv_is_proc_active(env) &&
      GMEnv_proc_num_field(env) == 0) {

    /*---Form filename---*/

    size_t len = strlen(stub);
    char* path = (char*)malloc((len+50) * sizeof(char));

    GMAssertAlways(env->num_stage < 1000000);
    GMAssertAlways(env->num_phase < 1000000);
    GMAssertAlways(GMEnv_num_proc(env) < 10000000000);
    sprintf(path, "%s_%06i_%06i_%010i", stub, env->phase_num,
            env->stage_num, GMEnv_proc_num(env));

    /*---Do output---*/

    double threshold = do_->threshold;

    FILE* file = fopen(path, "w");
    output_metrics_impl(metrics, do_, file, threshold, env);
    fclose(file);
    free(path);
  }

  /*---Output to stdout if requested---*/

  if (do_->verbosity > 1) {
    double threshold = do_->verbosity > 2 ? -1. : do_->threshold;
    output_metrics_impl(metrics, do_, stdout, threshold, env);
  }
}

/*===========================================================================*/
/*---Perform a single metrics computation run---*/

GMChecksum perform_run(const char* const options) {
  GMAssertAlways(options);

  /*---Convert options string to args---*/

  size_t len = strlen(options);
  char argstring[len+1];
  char* argv[len+1];
  int argc = 0;
  strcpy(argstring, options);
  gm_create_args(argstring, &argc, argv);

  return perform_run(argc, argv, options);
}

/*---------------------------------------------------------------------------*/

GMChecksum perform_run(int argc, char** argv, const char* const description) {

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();

  GMEnv_create(&env, MPI_COMM_WORLD, argc, argv, description);

  /*---Parse remaining unprocessed arguments---*/

  DriverOptions do_ = {0};
  do_.num_field_local_initialized = false;
  do_.num_field_active_initialized = false;
  do_.num_vector_local_initialized = false;
  do_.num_vector_active_initialized = false;
  do_.verbosity = 1;
  do_.stage_min_1based = 1;
  do_.stage_max_1based = env.num_stage;
  do_.phase_min_1based = 1;
  do_.phase_max_1based = env.num_phase;
  do_.input_file_path = NULL;
  do_.output_file_path_stub = NULL;
  //do_.problem_type = GM_PROBLEM_TYPE_RANDOM;
  do_.problem_type = GM_PROBLEM_TYPE_ANALYTIC;
  do_.threshold = -1.;
  do_.checksum = true;
  do_.num_incorrect = 0;

  finish_parsing(argc, argv, &do_, &env);

  if (do_.num_vector_local_initialized) {
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(&env);
    do_.num_vector_active = do_.num_vector;
  } else {
    /*---Pad up so that every proc has same number of vectors---*/
    do_.num_vector_local = GMVectors_num_local_required(
        do_.num_vector_active, &env);
    do_.num_vector = do_.num_vector_local *
      (size_t)GMEnv_num_proc_vector_i(&env);
  }

  if (do_.num_field_local_initialized) {
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(&env);
    do_.num_field_active = do_.num_field;
  } else {
    /*---Pad up so that every proc has same number of fields---*/
    do_.num_field_local = gm_ceil_i8(
        do_.num_field_active, GMEnv_num_proc_field(&env));
    do_.num_field = do_.num_field_local * (size_t) GMEnv_num_proc_field(&env);
  }

  /*---Initialize vectors---*/

  double vctime = 0;
  double time_beg = GMEnv_get_synced_time(&env);
  GMVectors vectors = GMVectors_null();
  GMVectors_create(&vectors, GMEnv_data_type_vectors(&env),
                   do_.num_field, do_.num_field_active,
                   do_.num_vector_local, &env);
  double time_end = GMEnv_get_synced_time(&env);
  vctime += time_end - time_beg;

  double intime = 0;
  time_beg = GMEnv_get_synced_time(&env);
  set_vectors(&vectors, &do_, &env);
  time_end = GMEnv_get_synced_time(&env);
  intime += time_end - time_beg;

  GMChecksum checksum = GMChecksum_null();
  checksum.computing_checksum = do_.checksum;

  double outtime = 0;
  double mctime = 0;
  double cktime = 0;

  /*---Loops over phases, stages---*/

  size_t num_elts_local_computed = 0;

  for (env.stage_num=do_.stage_min_1based-1;
       env.stage_num<=do_.stage_max_1based-1; ++env.stage_num) {

    for (env.phase_num=do_.phase_min_1based-1;
         env.phase_num<=do_.phase_max_1based-1; ++env.phase_num) {

      /*---Set up metrics container for results---*/

      time_beg = GMEnv_get_synced_time(&env);
      GMMetrics metrics = GMMetrics_null();
      GMMetrics_create(&metrics, GMEnv_data_type_metrics(&env),
                       do_.num_field, do_.num_field_active,
                       do_.num_vector_local,
                       do_.num_vector_active, &env);
      time_end = GMEnv_get_synced_time(&env);
      mctime += time_end - time_beg;

      /*---Calculate metrics---*/

      gm_compute_metrics(&metrics, &vectors, &env);

      num_elts_local_computed += metrics.num_elts_local_computed;

      /*---Output results---*/

      time_beg = GMEnv_get_synced_time(&env);
      output_metrics(&metrics, &do_, &env);
      time_end = GMEnv_get_synced_time(&env);
      outtime += time_end - time_beg;

      /*---Check correctness---*/

      check_metrics(&metrics, &do_, &env);

      /*---Compute checksum---*/

      if (do_.checksum) {
        time_beg = GMEnv_get_synced_time(&env);
        GMMetrics_checksum(&metrics, &checksum, &env);
        time_end = GMEnv_get_synced_time(&env);
        cktime += time_end - time_beg;
      }

      GMMetrics_destroy(&metrics, &env);

    }

  } /*---End loops over phases, stages---*/

  GMVectors_destroy(&vectors, &env);

  /*---Perform some checks---*/

  GMAssertAlways(env.cpu_mem == 0);
  GMAssertAlways(env.gpu_mem == 0);

  if (GMEnv_is_proc_active(&env)) {
  int mpi_code = 0;
    size_t num_elts_computed = 0;
    mpi_code = MPI_Allreduce(&num_elts_local_computed, &num_elts_computed, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             GMEnv_mpi_comm_vector(&env));
    GMAssertAlways(mpi_code == MPI_SUCCESS);

    if (GMEnv_num_way(&env) == GM_NUM_WAY_2 && GMEnv_all2all(&env) &&
        do_.phase_min_1based==1 && do_.phase_max_1based==env.num_phase) {
      GMAssertAlways(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) / 2);
    }

    if (GMEnv_num_way(&env) == GM_NUM_WAY_3 && GMEnv_all2all(&env) &&
        do_.stage_min_1based==1 && do_.stage_max_1based==env.num_stage) {
      GMAssertAlways(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) * (size_t)
                                          (do_.num_vector - 2) / 6);
    }
  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(&env) && GMEnv_proc_num(&env) == 0 &&
      do_.verbosity > 0) {
    //-----
    if (do_.checksum) {
      printf("metrics checksum ");
      int i = 0;
      for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
        printf("%s%li", i == 0 ? "" : "-",
               checksum.data[GM_CHECKSUM_SIZE - 1 - i]);
      }
      if (checksum.is_overflowed) {
        printf("-OVFL");
        printf("-%e", checksum.value_max);
      }
      printf(" ");
    }
    //-----
    printf("time %.6f", env.time);
    //-----
    printf(" ops %e", env.ops);
    if (env.time > 0) {
      printf(" rate %e", env.ops / env.time);
      printf(" rate/proc %e", env.ops / (env.time*GMEnv_num_proc(&env)) );
    }
    //-----
    printf(" cmp %e", env.compares);
    if (env.time > 0) {
      printf(" rate %e", env.compares / env.time);
      printf(" rate/proc %e", env.compares / (env.time*GMEnv_num_proc(&env)) );
    }
    //-----
    printf(" vctime %.6f", vctime);
    printf(" mctime %.6f", mctime);
    if (do_.checksum) {
      printf(" cktime %.6f", cktime);
    }
    if (NULL != do_.input_file_path) {
      printf(" intime %.6f", intime);
    }
    if (NULL != do_.output_file_path_stub) {
      printf(" outtime %.6f", outtime);
    }
    //-----
    printf(" cpumem %e", (double)env.cpu_mem_max);
    printf(" gpumem %e", (double)env.gpu_mem_max);
    //-----
    printf("\n");
  }

  GMAssertAlways(do_.num_incorrect == 0);

  /*---Finalize---*/

  GMEnv_destroy(&env);

  return checksum;
}

/*===========================================================================*/
