/*---------------------------------------------------------------------------*/
/*!
 * \file   driver_utils.cc
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
#include "compute_metrics.hh"
#include "driver_utils.hh"

#ifdef __cplusplus
extern "C" {
#endif

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
      long num_field = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_field >= 0
                    && "Invalid setting for num_field.");
      do_->num_field_active = num_field;
      do_->num_field_active_initialized = GM_BOOL_TRUE;
      do_->num_field_local_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--num_field_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_field_local.");
      long num_field_local = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && num_field_local >= 0 &&
                    (long)(int)num_field_local == num_field_local &&
                    "Invalid setting for num_field_local.");
      do_->num_field_local = num_field_local;
      do_->num_field_local_initialized = GM_BOOL_TRUE;
      do_->num_field_active_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_vector.");
      long num_vector = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_vector >= 0
                    && "Invalid setting for num_vector.");
      do_->num_vector_active = num_vector;
      do_->num_vector_active_initialized = GM_BOOL_TRUE;
      do_->num_vector_local_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_vector_local.");
      long num_vector_local = strtol(argv[i], NULL, 10);
      GMInsist(env, 0 == errno && num_vector_local >= 0 &&
                    (long)(int)num_vector_local == num_vector_local &&
                    "Invalid setting for num_vector_local.");
      do_->num_vector_local = num_vector_local;
      do_->num_vector_local_initialized = GM_BOOL_TRUE;
      do_->num_vector_active_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for verbosity.");
      do_->verbosity = atoi(argv[i]);
      GMInsist(env, do_->verbosity >= 0 && "Invalid setting for verbosity.");
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for num_stage.");
      long num_stage = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && num_stage >= 1
                    && (long)(int)num_stage == num_stage
                    && "Invalid setting for num_stage.");
      env->num_stage = num_stage;
      do_->stage_min = 1;
      do_->stage_max = env->num_stage;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for stage_min.");
      long stage_min = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && stage_min >= 1
                    && (long)(int)stage_min == stage_min
                    && "Invalid setting for stage_min.");
      do_->stage_min = stage_min;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc && "Missing value for stage_max.");
      long stage_max = strtol(argv[i], NULL, 10);
      GMInsist(env, errno == 0 && stage_max <= env->num_stage
                    && (long)(int)stage_max == stage_max
                    && "Invalid setting for stage_max.");
      do_->stage_max = stage_max;
    /*----------*/
    } else if (strcmp(argv[i], "--input_file") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for input_file`." : 0);
      do_->input_file_path = argv[i];
    /*----------*/
    } else if (strcmp(argv[i], "--output_file_stub") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for output_file_stub`." : 0);
      do_->output_file_path_stub = argv[i];
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
    } else {
    /*----------*/
      if (GMEnv_proc_num(env) == 0) {
        fprintf(stderr, "Invalid argument \"%s\".", argv[i]);
      }
      GMInsist(env, GM_BOOL_FALSE ? "Error: argument not recognized." : 0);
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

/*===========================================================================*/
/*---Utility to parse a string to construct arguments---*/

void create_args(char* argstring, int* argc, char** argv) {
  size_t len = strlen(argstring);

  argv[0] = &argstring[0];
  *argc = 1;
  _Bool is_delim_prev = GM_BOOL_TRUE;
  int i = 0;
  for (i = 0; i < (int)len; ++i) {
    const _Bool is_delim = argstring[i] == ' ' || argstring[i] == '\t';
    if (is_delim) {
      argstring[i] = 0;
    }
    if (is_delim_prev && !is_delim) {
      argv[*argc] = &(argstring[i]);
      (*argc)++;
    }
    is_delim_prev = is_delim;
  }
}

/*===========================================================================*/
/*---Set the entries of the vectors---*/

void input_vectors(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
    /*--------------------*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented; design not complete." : 0);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      if (do_->input_file_path != NULL) {
        const int fl = 0;
        const size_t field_base = fl +
          vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
        FILE* input_file = fopen(do_->input_file_path, "r");
        GMAssertAlways(NULL != input_file ? "Unable to open input file." : 0);
        int vl = 0;
        for (vl = 0; vl < vectors->num_vector_local; ++vl) {
          const size_t vector = vl +
              vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
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
      } else { /*--------------------*/
        int vl = 0;
        for (vl = 0; vl < vectors->num_vector_local; ++vl) {
          size_t vector = vl +
              vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
          /*---Fill pad vectors with copies of the last vector---*/
          const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                       vector : do_->num_vector_active-1;
          int fl = 0;
          for (fl = 0; fl < vectors->num_field_local; ++fl) {
            size_t field = fl +
                vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
            if (field >= vectors->num_field_active) {
              continue;
            }
            /*---Compute element unique id---*/
            const size_t uid = field + vectors->num_field_active*vector_capped;
            /*---Generate large random number---*/
            size_t rand1 = uid;
            rand1 = gm_randomize(rand1);
            rand1 = gm_randomize(rand1);
            size_t rand2 = uid;
            rand2 = gm_randomize(rand2);
            rand2 = gm_randomize(rand2);
            rand2 = gm_randomize(rand2);
            const size_t randomize_max = gm_randomize_max();
            size_t rand_value = rand1 + randomize_max * rand2;
            /*---Reduce so that after summing num_field times the integer
                 still exactly representable by floating point type---*/
            const size_t rand_max = randomize_max * randomize_max;
            GMAssertAlways(FLT_RADIX == 2);
            const int mant_dig = sizeof(GMFloat) == 8 ? DBL_MANT_DIG :
                                                        FLT_MANT_DIG;
            const int log2_num_summands_3way_numerator = 2;
            const int shift_amount = gm_log2(log2_num_summands_3way_numerator*
                                             rand_max*vectors->num_field_active)
                                     - mant_dig;
            rand_value >>= shift_amount > 0 ? shift_amount : 0;
            /*---Store---*/
            GMFloat float_value = (GMFloat)rand_value;
            GMAssertAlways((size_t)float_value == rand_value);
            GMAssertAlways(float_value * vectors->num_field_active <=
                           ((size_t)1)<<mant_dig);
            GMVectors_float_set(vectors, fl, vl, float_value, env);
            /*---Print---*/
            if (do_->verbosity > 2) {
              printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                     GMEnv_proc_num_vector_i(env), vl,
                     GMEnv_proc_num_field(env), fl, float_value);
            }
          } /*---field_local---*/
        }   /*---vector_local---*/
      } /*---if---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
      if (do_->input_file_path != NULL) {
        GMInsist(env, NULL == do_->input_file_path
                      ? "File input for this case not yet implemented." : 0);
      } else { /*--------------------*/
        int vl = 0;
        for (vl = 0; vl < vectors->num_vector_local; ++vl) {
          size_t vector = vl +
              vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
          /*---Fill pad vectors with copies of the last vector---*/
          const size_t vector_capped = vector <= do_->num_vector_active-1 ?
                                       vector : do_->num_vector_active-1;
          int fl;
          for (fl = 0; fl < vectors->num_field_local; ++fl) {
            size_t field = fl +
                vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
            if (field >= vectors->num_field_active) {
              continue;
            }
            /*---Compute element unique id---*/
            const size_t uid = field + vectors->num_field_active*vector_capped;
            size_t index = uid;
            /*---Randomize---*/
            index = gm_randomize(index);
            index = gm_randomize(index);
            /*---Calculate random number between 0 and 3---*/
            const float float_rand_value = index / (float)gm_randomize_max();
            /*---Create 2-bit value - make extra sure less than 4---*/
            GMBits2 value = (int)((4. - 1e-5) * float_rand_value);
            /*---Store---*/
            GMVectors_bits2_set(vectors, fl, vl, value, env);
            /*---Print---*/
            if (do_->verbosity > 2) {
              printf("vec_proc %i vec %i "
                     "field_proc %i field %i value %.1i%.1i\n",
                     GMEnv_proc_num_vector_i(env), vl,
                     GMEnv_proc_num_field(env), fl, value / 2, value % 2);
            }
          } /*---fl---*/
        }   /*---vl---*/
      } /*---if---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values to file---*/

void output_metrics_file(GMMetrics* metrics, DriverOptions* do_,
                         FILE* file, double threshold, GMEnv* env) {
  GMAssertAlways(NULL != metrics);
  GMAssertAlways(NULL != do_);
  GMAssertAlways(NULL != file);
  GMAssertAlways(NULL != env);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  /*---Due to redundancy, only results from some processors are needed---*/
  if (GMEnv_proc_num_field(env) != 0) {
    return;
  }

  switch (GMEnv_data_type_metrics(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
    /*--------------------*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/

      if (GMEnv_num_way(env) == GM_NUM_WAY_2) {
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
            = GMMetrics_czekanowski_get_from_index(metrics, index, env);
          if (!(threshold >= 0. && value > threshold)) {
            continue;
          }
          fprintf(file,
            sizeof(GMFloat) == 8 ?
            "element (%li,%li): value: %.17e\n" :
            "element (%li,%li): value: %.8e\n", coord0, coord1, value);
        }
      }

      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
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
            = GMMetrics_czekanowski_get_from_index(metrics, index, env);
          if (!(threshold >= 0. && value > threshold)) {
            continue;
          }
          fprintf(file, sizeof(GMFloat) == 8 ?
                        "element (%li,%li,%li): value: %.17e\n" :
                        "element (%li,%li,%li): value: %.8e\n",
                        coord0, coord1, coord2, value);
        }
      }

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
      size_t index = 0;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        int is_active = GM_BOOL_TRUE;
        int coord_num = 0;
        for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
          const size_t coord = GMMetrics_coord_global_from_index(metrics,
               index, coord_num, env);
          is_active = is_active && coord < metrics->num_vector_active;
        }
        if (is_active) {
          /*---Coordinate---*/
          fprintf(file, "element (");
          int coord_num = 0;
          for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
            if (coord_num > 0) {
              fprintf(file, ",");
            }
            const size_t coord = GMMetrics_coord_global_from_index(metrics,
                index, coord_num, env);
            /*---Present to the user as 1-based---*/
            fprintf(file, "%li", 1 + coord);
          }
          /*---Value---*/
          fprintf(file, "): values:");
          int i;
          for (i = 0; i < 2; ++i) {
            int j;
            for (j = 0; j < 2; ++j) {
              fprintf(file, " %.17e",
                  GMMetrics_ccc_get_from_index_2(metrics, index, i, j, env));
            }
          }
          fprintf(file, "    [from proc %i]\n", GMEnv_proc_num(env));
        }
      } /*---for index---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/
      size_t index = 0;
      for (index = 0; index < metrics->num_elts_local; ++index) {
        int is_active = GM_BOOL_TRUE;
        int coord_num = 0;
        for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
          const size_t coord = GMMetrics_coord_global_from_index(metrics,
              index, coord_num, env);
          is_active = is_active && coord < metrics->num_vector_active;
        }
        if (is_active) {
          /*---Coordinate---*/
          fprintf(file, "element (");
          int coord_num = 0;
          for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
            if (coord_num > 0) {
              fprintf(file, ",");
            }
            const size_t coord = GMMetrics_coord_global_from_index(metrics,
                index, coord_num, env);
            /*---Present to the user as 1-based---*/
            fprintf(file, "%li", 1 + coord);
          }
          /*---Value---*/
          fprintf(file, "): values:");
          int i;
          for (i = 0; i < 2; ++i) {
            int j;
            for (j = 0; j < 2; ++j) {
              int k;
              for (k = 0; k < 2; ++k) {
                fprintf(file, " %.17e", GMMetrics_ccc_get_from_index_3(
                   metrics, index, i, j, k, env));
              }
            }
          }
          fprintf(file, "    [from proc %i]\n", GMEnv_proc_num(env));
        }
      } /*---for index---*/
    } break;
    /*--------------------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

void output_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);

  char* stub = do_->output_file_path_stub;

  if (NULL != stub && GMEnv_is_proc_active(env) &&
      GMEnv_proc_num_field(env) == 0) {

    size_t len = strlen(stub);

    char* path = (char*)malloc((len+50) * sizeof(char));

    GMAssertAlways(env->num_stage < 1000000);
    GMAssertAlways(GMEnv_num_proc(env) < 10000000000);
    sprintf(path, "%s_%06i_%010i", stub, env->stage_num, GMEnv_proc_num(env));

    double threshold = 0.;

    FILE* file = fopen(path, "w");
    output_metrics_file(metrics, do_, file, threshold, env);
    fclose(file);
    free(path);
  }

  if (do_->verbosity > 1) {
    output_metrics_file(metrics, do_, stdout, -1., env);
  }
}

/*===========================================================================*/
/*---Perform a single metrics computation run---*/

GMChecksum perform_run(int argc, char** argv, char const * const description) {

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, argv, description);

  /*---Parse remaining unprocessed arguments---*/

  DriverOptions do_ = {0};
  do_.num_field_local_initialized = GM_BOOL_FALSE;
  do_.num_field_active_initialized = GM_BOOL_FALSE;
  do_.num_vector_local_initialized = GM_BOOL_FALSE;
  do_.num_vector_active_initialized = GM_BOOL_FALSE;
  do_.verbosity = 1;
  do_.stage_min = 1;
  do_.stage_max = env.num_stage;
  do_.input_file_path = NULL;
  do_.output_file_path_stub = NULL;

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
  input_vectors(&vectors, &do_, &env);
  time_end = GMEnv_get_synced_time(&env);
  intime += time_end - time_beg;

  GMChecksum checksum = GMChecksum_null();

  double outtime = 0;
  double mctime = 0;
  double cktime = 0;

  /*---Loop over stages---*/

  for (env.stage_num=do_.stage_min-1;
       env.stage_num<=do_.stage_max-1; ++env.stage_num) {

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

    /*---Output results---*/

    time_beg = GMEnv_get_synced_time(&env);
    output_metrics(&metrics, &do_, &env);
    time_end = GMEnv_get_synced_time(&env);
    outtime += time_end - time_beg;

    /*---Compute checksum---*/

    time_beg = GMEnv_get_synced_time(&env);
    GMMetrics_checksum(&metrics, &checksum, &env);
    time_end = GMEnv_get_synced_time(&env);
    cktime += time_end - time_beg;

    GMMetrics_destroy(&metrics, &env);
  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(&env) && GMEnv_proc_num(&env) == 0 &&
      do_.verbosity > 0) {
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
    printf(" time %.6f", env.time);
    printf(" ops %e", env.ops);
    if (env.time > 0) {
      printf(" rate %e", env.ops / env.time);
    }
    printf(" vctime %.6f", vctime);
    printf(" mctime %.6f", mctime);
    printf(" cktime %.6f", cktime);
    if (NULL != do_.input_file_path) {
      printf(" intime %.6f", intime);
    }
    if (NULL != do_.output_file_path_stub) {
      printf(" outtime %.6f", outtime);
    }
    printf("\n");
  }

  /*---Finalize---*/

  GMVectors_destroy(&vectors, &env);

  GMEnv_destroy(&env);

  return checksum;
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
