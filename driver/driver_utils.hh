/*---------------------------------------------------------------------------*/
/*!
 * \file   driver_utils.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code utilitiy functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_driver_utils_hh_
#define _gm_driver_utils_hh_

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <float.h>

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Struct to hold driver options (options not in GMEnv)---*/

typedef struct {
  int num_field;
  int num_vector_local;
  size_t num_vector;
  size_t num_vector_active;
  _Bool num_field_initialized;
  _Bool num_vector_local_initialized;
  _Bool num_vector_active_initialized;
  int verbosity;
  int stage_min;
  int stage_max;
  char* input_file_path;
} DriverOptions;

/*===========================================================================*/
/*---Parse remaining unprocessed arguments---*/

static void finish_parsing(int argc,
                           char** argv,
                           DriverOptions* driver_options,
                           GMEnv* env) {
  int i = 0;
  for (i = 1; i < argc; ++i) {
    /*----------*/
    if (strcmp(argv[i], "--num_field") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_field." : 0);
      //FIX: safe atoi
      driver_options->num_field = atoi(argv[i]);
      GMInsist(env, driver_options->num_field >= 0
                    ? "Invalid setting for num_field."
                    : 0);
      driver_options->num_field_initialized = GM_BOOL_TRUE;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector_local") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_vector_local." : 0);
      //FIX: safe atoi
      driver_options->num_vector_local = atoi(argv[i]);
      GMInsist(env, driver_options->num_vector_local >= 0
                        ? "Invalid setting for num_vector_local."
                        : 0);
      driver_options->num_vector_local_initialized = GM_BOOL_TRUE;
      driver_options->num_vector_active_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--num_vector") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_vector." : 0);
      //FIX: safe atoi
      driver_options->num_vector_active = atoi(argv[i]);
      GMInsist(env, driver_options->num_vector_active >= 0
                        ? "Invalid setting for num_vector."
                        : 0);
      driver_options->num_vector_active_initialized = GM_BOOL_TRUE;
      driver_options->num_vector_local_initialized = GM_BOOL_FALSE;
    /*----------*/
    } else if (strcmp(argv[i], "--verbosity") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for verbosity." : 0);
      //FIX: safe atoi
      driver_options->verbosity = atoi(argv[i]);
      GMInsist(env, driver_options->verbosity >= 0
                    ? "Invalid setting for verbosity."
                    : 0);
    /*----------*/
    } else if (strcmp(argv[i], "--num_stage") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_stage." : 0);
      //FIX: safe atoi
      env->num_stage = atoi(argv[i]);
      GMInsist(env, env->num_stage >= 1 ? "Invalid setting for num_stage." : 0);
      driver_options->stage_min = 1;
      driver_options->stage_max = env->num_stage;
    /*----------*/
    } else if (strcmp(argv[i], "--stage_min") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for stage_min." : 0);
      //FIX: safe atoi
      driver_options->stage_min = atoi(argv[i]);
      GMInsist(env, driver_options->stage_min >= 1
                    ? "Invalid setting for stage_min."
                    : 0);
    /*----------*/
    } else if (strcmp(argv[i], "--stage_max") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for stage_max`." : 0);
      //FIX: safe atoi
      driver_options->stage_max = atoi(argv[i]);
      GMInsist(env, driver_options->stage_max <= env->num_stage
                    ? "Invalid setting for stage_max."
                    : 0);
    /*----------*/
    } else if (strcmp(argv[i], "--input_file") == 0) {
    /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for input_file`." : 0);
      driver_options->input_file_path = argv[i];
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

  GMInsist(env, driver_options->num_field_initialized
                ? "Error: num_field not set." : 0);
  GMInsist(env, driver_options->num_vector_local_initialized ||
                driver_options->num_vector_active_initialized
                ? "Error: must set num_vector_local or num_vector."
                : 0);
}

/*===========================================================================*/
/*---Utility to parse a string to construct arguments---*/

static void create_args(char* argstring, int* argc, char** argv) {
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

static void input_vectors(GMVectors* vectors,
                          DriverOptions* driver_options, GMEnv* env) {
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
      int vl = 0;
      for (vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        int fl = 0;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          /*---compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
          size_t index = uid;
          /*---randomize---*/
          index = gm_randomize(index);
          /*---Calculate random number between 0 and 1---*/
          GMFloat rand_value = index / (GMFloat)gm_randomize_max();
          /*---Create single bit value---*/
          _Bool value = rand_value < .5 ? GM_BOOL_FALSE : GM_BOOL_TRUE;
          /*---Store---*/
          GMVectors_bit_set(vectors, fl, vl, value, env);
          /*---Print---*/
          if (driver_options->verbosity > 2) {
          }
        } /*---field---*/
      }   /*---vector_local---*/

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      if (driver_options->input_file_path != NULL) {
        const int fl = 0;
        const size_t field_base = fl +
          vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
        FILE* input_file = fopen(driver_options->input_file_path, "r");
        GMAssertAlways(NULL != input_file ? "Unable to open input file." : 0);
        int vl = 0;
        for (vl = 0; vl < vectors->num_vector_local; ++vl) {
          const size_t vector = vl +
              vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
          if (vector < driver_options->num_vector_active) {
            size_t addr_file = field_base + vectors->num_field * vector;
            int fseek_success = fseek(input_file, addr_file, SEEK_SET);
            fseek_success += 0; /*---Avoid unused var warning---*/
            GMAssertAlways(0 == fseek_success);
            GMFloat* const addr_mem = GMVectors_float_ptr(vectors, fl, vl, env);
            /*---NOTE: the following call is ok since has no side effects---*/
            GMAssertAlways(fl+1 >= vectors->num_field_local ||
                        GMVectors_float_ptr(vectors, fl+1, vl, env) == addr_mem + 1
                        ? "Vector layout is inompatible with operation." : 0);
            size_t num_read = fread(addr_mem, sizeof(GMFloat),
                                    vectors->num_field_local, input_file);
            num_read += 0; /*---Avoid unused var warning---*/
            GMAssertAlways(sizeof(GMFloat)*vectors->num_field_local == num_read);
          } else {
            const int vl_prev = vl - 1;
            GMAssertAlways(vl_prev >= 0);
            int fl = 0;
            for (fl = 0; fl < vectors->num_field_local; ++fl) {
              GMFloat pad_value = GMVectors_float_get(vectors, fl, vl_prev, env);;
              GMVectors_float_set(vectors, fl, vl, pad_value, env);
            }
          }
        } /*---vl---*/
        fclose(input_file);
      } else {
        GMAssertAlways(driver_options->num_vector==driver_options->num_vector_active);
        int vl = 0;
        for (vl = 0; vl < vectors->num_vector_local; ++vl) {
          size_t vector = vl +
              vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
          int fl = 0;
          for (fl = 0; fl < vectors->num_field_local; ++fl) {
            size_t field = fl +
                vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
            /*---Compute element unique id---*/
            const size_t uid = field + vectors->num_field * vector;
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
            const int shift_amount = gm_log2(rand_max * vectors->num_field) -
                                     mant_dig;
            rand_value >>= shift_amount > 0 ? shift_amount : 0;
            /*---Store---*/
            GMFloat float_value = (GMFloat)rand_value;
            GMAssertAlways((size_t)float_value == rand_value);
            GMAssertAlways(float_value * vectors->num_field <=
                           ((size_t)1)<<mant_dig);
            GMVectors_float_set(vectors, fl, vl, float_value, env);
            /*---Print---*/
            if (driver_options->verbosity > 2) {
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
      GMInsist(env, NULL == driver_options->input_file_path
                    ? "File input for this case not yet implemented." : 0);
      int vl = 0;
      for (vl = 0; vl < vectors->num_vector_local; ++vl) {
        size_t vector = vl +
            vectors->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env);
        int fl;
        for (fl = 0; fl < vectors->num_field_local; ++fl) {
          size_t field = fl +
              vectors->num_field_local * (size_t)GMEnv_proc_num_field(env);
          /*---Compute element unique id---*/
          const size_t uid = field + vectors->num_field * vector;
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
          if (driver_options->verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, value / 2, value % 2);
          }
        } /*---fl---*/
      }   /*---vl---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Output the result metrics values---*/

static void output_metrics(GMMetrics* metrics, DriverOptions* driver_options,
                           GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);

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
      if (driver_options->verbosity > 1) {
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          int is_active = GM_BOOL_TRUE;
          int coord_num = 0;
          for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
            const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                coord_num, env);
            is_active = is_active && coord < metrics->num_vector_active;
          }
          if (is_active) {
            /*---Coordinate---*/
            printf("element (");
            int coord_num = 0;
            for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
              if (coord_num > 0) {
                printf(",");
              }
              const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                  coord_num, env);
              /*---Present to the user as 1-based---*/
              printf("%li", 1 + coord);
            }
            /*---Value---*/
            printf("): value:");
            printf(" %.17e",
                   GMMetrics_czekanowski_get_from_index(metrics, index, env));
            printf("    [from proc %i]\n", GMEnv_proc_num(env));

          }
        } /*---for index---*/
      } /*---verbosity---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
      if (driver_options->verbosity > 1) {
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          int is_active = GM_BOOL_TRUE;
          int coord_num = 0;
          for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
            const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                coord_num, env);
            is_active = is_active && coord < metrics->num_vector_active;
          }
          if (is_active) {
            /*---Coordinate---*/
            printf("element (");
            int coord_num = 0;
            for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
              if (coord_num > 0) {
                printf(",");
              }
              const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                  coord_num, env);
              /*---Present to the user as 1-based---*/
              printf("%li", 1 + coord);
            }
            /*---Value---*/
            printf("): values:");
            int i;
            for (i = 0; i < 2; ++i) {
              int j;
              for (j = 0; j < 2; ++j) {
                printf(" %.17e",
                       GMMetrics_ccc_get_from_index_2(metrics, index, i, j, env));
              }
            }
            printf("    [from proc %i]\n", GMEnv_proc_num(env));
          }
        } /*---for index---*/
      }
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/
      if (driver_options->verbosity > 1) {
        size_t index = 0;
        for (index = 0; index < metrics->num_elts_local; ++index) {
          int is_active = GM_BOOL_TRUE;
          int coord_num = 0;
          for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
            const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                coord_num, env);
            is_active = is_active && coord < metrics->num_vector_active;
          }
          if (is_active) {
            /*---Coordinate---*/
            printf("element (");
            int coord_num = 0;
            for (coord_num = 0; coord_num < GMEnv_num_way(env); ++coord_num) {
              if (coord_num > 0) {
                printf(",");
              }
              const size_t coord = GMMetrics_coord_global_from_index(metrics, index,
                                                                  coord_num, env);
              /*---Present to the user as 1-based---*/
              printf("%li", 1 + coord);
            }
            /*---Value---*/
            printf("): values:");
            int i;
            for (i = 0; i < 2; ++i) {
              int j;
              for (j = 0; j < 2; ++j) {
                int k;
                for (k = 0; k < 2; ++k) {
                  printf(" %.17e", GMMetrics_ccc_get_from_index_3(metrics, index,
                                                                  i, j, k, env));
                }
              }
            }
            printf("    [from proc %i]\n", GMEnv_proc_num(env));
          }
        } /*---for index---*/
      }
    } break;
    /*--------------------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Perform a single metrics computation run---*/

static GMChecksum perform_run(int argc, char** argv,
                              char const * const description) {

  /*---Initialize environment---*/

  GMEnv env = GMEnv_null();
  GMEnv_create_from_args(&env, argc, argv, description);

  /*---Parse remaining unprocessed arguments---*/

  DriverOptions driver_options = {0};
  driver_options.num_field_initialized = GM_BOOL_FALSE;
  driver_options.num_vector_local_initialized = GM_BOOL_FALSE;
  driver_options.num_vector_active_initialized = GM_BOOL_FALSE;
  driver_options.verbosity = 1;
  driver_options.stage_min = 1;
  driver_options.stage_max = env.num_stage;
  driver_options.input_file_path = NULL;

  finish_parsing(argc, argv, &driver_options, &env);

  if (driver_options.num_vector_local_initialized) {
    driver_options.num_vector = driver_options.num_vector_local *
      (size_t) GMEnv_num_proc_vector_i(&env);
    driver_options.num_vector_active = driver_options.num_vector;
  } else {
    driver_options.num_vector_local = gm_ceil_i8(driver_options.num_vector_active,
                                                 GMEnv_num_proc_vector_i(&env));
    driver_options.num_vector = driver_options.num_vector_local *
      (size_t) GMEnv_num_proc_vector_i(&env);
  }

  /*---Initialize vectors---*/

  GMVectors vectors = GMVectors_null();
  GMVectors_create(&vectors, GMEnv_data_type_vectors(&env),
                   driver_options.num_field,
                   driver_options.num_vector_local, &env);

  input_vectors(&vectors, &driver_options, &env);

  GMChecksum checksum = GMChecksum_null();

  /*---Loop over stages---*/

  for (env.stage_num=driver_options.stage_min-1;
       env.stage_num<=driver_options.stage_max-1; ++env.stage_num) {

    /*---Set up metrics container for results---*/

    GMMetrics metrics = GMMetrics_null();
    GMMetrics_create(&metrics, GMEnv_data_type_metrics(&env),
                     driver_options.num_field, driver_options.num_vector_local,
                     driver_options.num_vector_active, &env);

    /*---Calculate metrics---*/

    gm_compute_metrics(&metrics, &vectors, &env);

    /*---Output results---*/

    output_metrics(&metrics, &driver_options, &env);

    /*---Compute checksum---*/

    GMMetrics_checksum(&metrics, &checksum, &env);

    GMMetrics_destroy(&metrics, &env);
  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(&env) && GMEnv_proc_num(&env) == 0 &&
      driver_options.verbosity > 0) {
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

/*===========================================================================*/

#endif /*---_gm_driver_utils_hh_---*/

/*---------------------------------------------------------------------------*/
