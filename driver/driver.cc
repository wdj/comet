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
#include "compute_metrics.hh"
#include "driver.hh"

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
        GMInsist(env,
                 GM_BOOL_FALSE ? "Invalid setting for problem_type." : 0);
      }
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
/*---Set the entries of the vectors---*/

int mant_dig() {
  GMAssertAlways(FLT_RADIX == 2);
  return sizeof(GMFloat) == 8 ? DBL_MANT_DIG : FLT_MANT_DIG;
}

/*---------------------------------------------------------------------------*/

void set_vectors_random(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(do_ != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(do_->problem_type == GM_PROBLEM_TYPE_RANDOM);

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
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
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
          const int log2_num_summands_3way_numerator = 2;
          const int shift_amount = gm_log2(log2_num_summands_3way_numerator*
                                           rand_max*vectors->num_field_active)
                                   - mant_dig();
          rand_value >>= shift_amount > 0 ? shift_amount : 0;
          /*---Store---*/
          GMFloat float_value = (GMFloat)rand_value;
          GMAssertAlways((size_t)float_value == rand_value);
          GMAssertAlways(float_value * vectors->num_field_active <
                         ((size_t)1)<<mant_dig());
          GMVectors_float_set(vectors, fl, vl, float_value, env);
          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, float_value);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
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
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*---------------------------------------------------------------------------*/

size_t perm_shuffle(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)1)) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  const size_t nhalf = (n+1-key)/2;
  const size_t result = i < nhalf ? 2*i + key : 2*(i-nhalf) + 1 - key;
  GMAssert(result>=0 && result<n);
  return result;
}

/*---------------------------------------------------------------------------*/

enum {NUM_SHUFFLE = 3};

size_t perm(size_t key, size_t i, size_t n) {
  GMAssert((key & (~(size_t)((1<<NUM_SHUFFLE)-1))) == 0);
  GMAssert(i>=0 && i<n);
  GMAssert(n>=0);

  const size_t result0 = perm_shuffle(key&1, i, n);
  const size_t result1 = perm_shuffle((key&2)>>1, result0, n);
  const size_t result2 = perm_shuffle((key&4)>>2, result1, n);
  const size_t result = result2;
  GMAssert(result>=0 && result<n);
  return result;
}

/*---------------------------------------------------------------------------*/

void set_vectors_analytic(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(do_ != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(do_->problem_type == GM_PROBLEM_TYPE_ANALYTIC);

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  const size_t nf = vectors->num_field_active;
  const size_t nv = do_->num_vector_active;

  const size_t ng = 1 << NUM_SHUFFLE;
  const size_t gs = (nf+ng-1) / ng;

  const size_t max_float = ((size_t)1) << mant_dig();

  const size_t value_limit = (max_float - 1) / nf;

  const size_t value_min = 1;
  const size_t value_max = (nv+value_min) < value_limit ?
                           (nv+value_min) : value_limit;

  switch (GMEnv_data_type_vectors(env)) {
    /*--------------------*/
    case GM_DATA_TYPE_BITS1: {
    /*--------------------*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented; design not complete." : 0);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
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
          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nf);
          const size_t g = pf / gs;
          GMAssert(g>=0 && g<ng);

          const size_t pv = perm(g, v, nv);

          const size_t value = value_min + ( pv * value_max ) / (nv+value_min);

          const GMFloat float_value = value;

          //OLD GMFloat float_value = 1 + 0*vector_capped;

          GMAssertAlways(float_value * vectors->num_field_active < max_float);
          GMVectors_float_set(vectors, fl, vl, float_value, env);

          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl, float_value);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
#pragma omp parallel for
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
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
          /*---Create 2-bit value - make extra sure less than 4---*/

          const size_t f = field;
          const size_t v = vector_capped;

          const size_t pf = perm(0, f, nf);
          const size_t g = pf / gs;
          GMAssert(g>=0 && g<ng);

          const size_t pv = perm(g, v, nv);

          const size_t value = value_min + ( pv * value_max ) / (nv+value_min);

          const GMBits2 bval = ((size_t)3) & value;

          //OLD GMBits2 bval = 3 + 0*( vector_capped % 4);

          /*---Store---*/
          GMVectors_bits2_set(vectors, fl, vl, bval, env);
          /*---Print---*/
          if (do_->verbosity > 2) {
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   GMEnv_proc_num_vector_i(env), vl,
                   GMEnv_proc_num_field(env), fl,
                   (int)bval / 2, (int)bval % 2);
          }
        } /*---field_local---*/
      }   /*---vector_local---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*---------------------------------------------------------------------------*/

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(do_ != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(do_->input_file_path != NULL);

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
      GMInsist(env, GM_BOOL_FALSE ? "Not yet implemented." : 0);
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*---------------------------------------------------------------------------*/

void set_vectors(GMVectors* vectors, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(do_ != NULL);
  GMAssertAlways(env != NULL);

  if (do_->input_file_path != NULL) {
    set_vectors_from_file(vectors, do_, env);
  } else if (do_->problem_type == GM_PROBLEM_TYPE_RANDOM) {
    set_vectors_random(vectors, do_, env);
  } else {
    set_vectors_analytic(vectors, do_, env);
  }
}

/*===========================================================================*/
/*---Check correctness of metrics, if possible---*/

void check_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(do_ != NULL);
  GMAssertAlways(env != NULL);

  if (GM_PROBLEM_TYPE_ANALYTIC != do_->problem_type ||
      NULL != do_->input_file_path) {
    return;
  }

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  //OLD const GMFloat num_field_f = metrics->num_field_active;

  const size_t nf = metrics->num_field_active;
  const size_t nv = metrics->num_vector_active;

  const size_t ng = 1 << NUM_SHUFFLE;
  const size_t gs = (nf+ng-1) / ng;

  const size_t max_float = ((size_t)1) << mant_dig();
  const size_t value_limit = (max_float - 1) / nf;

  const size_t value_min = 1;
  const size_t value_max = (nv+value_min) < value_limit ?
                           (nv+value_min) : value_limit;

  size_t num_misses = 0;

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
#pragma omp parallel for reduction(+:num_misses)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          if (vi >= metrics->num_vector_active ||
              vj >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czekanowski_get_from_index(metrics, index, env);

          GMFloat n = 0;
          GMFloat d = 0;

          for (size_t g=0; g<ng; ++g) {

            const size_t pf_min = g * gs;
            const size_t pf_max = gm_min_i8((g+1) * gs, nf);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nv);
            const size_t pvj = perm(g, vj, nv);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nv+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nv+value_min);
            n += gm_min_i8(value_i, value_j) * gs_this;
            d += (value_i + value_j) * gs_this;

          } //---g

          const GMFloat value_expected = (((GMFloat)2) * n) / d;

          //OLD const GMFloat value_expected = 2. * num_field_f /
          //OLD                                (num_field_f + num_field_f);

          num_misses += value_expected != value;
        } //---for index
      } //---if
      if (GMEnv_num_way(env) == GM_NUM_WAY_3) {
#pragma omp parallel for reduction(+:num_misses)
        for (size_t index = 0; index < metrics->num_elts_local; ++index) {
          const size_t vi =
            GMMetrics_coord_global_from_index(metrics, index, 0, env);
          const size_t vj =
            GMMetrics_coord_global_from_index(metrics, index, 1, env);
          const size_t vk =
            GMMetrics_coord_global_from_index(metrics, index, 2, env);
          if (vi >= metrics->num_vector_active ||
              vj >= metrics->num_vector_active ||
              vk >= metrics->num_vector_active) {
            continue;
          }
          const GMFloat value
            = GMMetrics_czekanowski_get_from_index(metrics, index, env);

          GMFloat n = 0;
          GMFloat d = 0;

          for (size_t g=0; g<ng; ++g) {

            const size_t pf_min = g * gs;
            const size_t pf_max = gm_min_i8((g+1) * gs, nf);
            const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

            const size_t pvi = perm(g, vi, nv);
            const size_t pvj = perm(g, vj, nv);
            const size_t pvk = perm(g, vk, nv);

            const size_t value_i = value_min + ( pvi * value_max ) /
                                               (nv+value_min);
            const size_t value_j = value_min + ( pvj * value_max ) /
                                               (nv+value_min);
            const size_t value_k = value_min + ( pvk * value_max ) /
                                               (nv+value_min);

            n += gm_min_i8(value_i, value_j) * gs_this;
            n += gm_min_i8(value_i, value_k) * gs_this;
            n += gm_min_i8(value_j, value_k) * gs_this;

            n -= gm_min_i8(value_i, gm_min_i8(value_j, value_k)) * gs_this;

            d += (value_i + value_j + value_k) * gs_this;

          } //---g

          const GMFloat value_expected = (((GMFloat)1.5) * n) / d;

          //OLD const GMFloat value_expected = ((GMFloat)1.5) *
          //OLD  (num_field_f + num_field_f + num_field_f - num_field_f) /
          //OLD  (num_field_f + num_field_f + num_field_f);

          num_misses += value_expected != value;
        } //---for index
      } //---if
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
#pragma omp parallel for reduction(+:num_misses)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        if (vi >= metrics->num_vector_active ||
            vj >= metrics->num_vector_active) {
          continue;
        }
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            const GMFloat value =
                GMMetrics_ccc_get_from_index_2(metrics, index, i, j, env);

            GMTally1 rij = 0;
            GMTally1 si = 0;
            GMTally1 sj = 0;

            for (size_t g=0; g<ng; ++g) {

              const size_t pf_min = g * gs;
              const size_t pf_max = gm_min_i8((g+1) * gs, nf);
              const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

              const size_t pvi = perm(g, vi, nv);
              const size_t pvj = perm(g, vj, nv);

              const size_t value_i = value_min + ( pvi * value_max ) /
                                                 (nv+value_min);
              const size_t value_j = value_min + ( pvj * value_max ) /
                                                 (nv+value_min);

              const GMBits2 bval_i = ((size_t)3) & value_i;
              const GMBits2 bval_j = ((size_t)3) & value_j;

              const int bval_i_0 = !!(bval_i&1);
              const int bval_i_1 = !!(bval_i&2);
              const int bval_j_0 = !!(bval_j&1);
              const int bval_j_1 = !!(bval_j&2);

              si += ((bval_i_0 == i) + (bval_i_1 == i)) * gs_this;
              sj += ((bval_j_0 == j) + (bval_j_1 == j)) * gs_this;

              rij += (((bval_i_0 == i) && (bval_j_0 == j)) +
                      ((bval_i_0 == i) && (bval_j_1 == j)) +
                      ((bval_i_1 == i) && (bval_j_0 == j)) +
                      ((bval_i_1 == i) && (bval_j_1 == j))) *
                     gs_this;

            } //---g

            //OLD const GMFloat rij = i==1 && j==1 ? (4 * num_field_f) : 0;
            //OLD const GMFloat si = i==1 ? (2 * num_field_f) : 0;
            //OLD const GMFloat sj = j==1 ? (2 * num_field_f) : 0;

            const GMFloat value_expected =
              GMMetrics_ccc_value_2(metrics, rij, si, sj, env);
            num_misses += value_expected != value;
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY4X2: {
    /*--------------------*/
#pragma omp parallel for reduction(+:num_misses)
      for (size_t index = 0; index < metrics->num_elts_local; ++index) {
        const size_t vi =
          GMMetrics_coord_global_from_index(metrics, index, 0, env);
        const size_t vj =
          GMMetrics_coord_global_from_index(metrics, index, 1, env);
        const size_t vk =
          GMMetrics_coord_global_from_index(metrics, index, 2, env);
        if (vi >= metrics->num_vector_active ||
            vj >= metrics->num_vector_active ||
            vk >= metrics->num_vector_active) {
          continue;
        }
        for (int i = 0; i < 2; ++i) {
          for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 2; ++k) {
              const GMFloat value =
               GMMetrics_ccc_get_from_index_3( metrics, index, i, j, k, env);

              GMTally1 rijk = 0;
              GMTally1 si = 0;
              GMTally1 sj = 0;
              GMTally1 sk = 0;

              for (size_t g=0; g<ng; ++g) {

                const size_t pf_min = g * gs;
                const size_t pf_max = gm_min_i8((g+1) * gs, nf);
                const size_t gs_this = pf_max >= pf_min ? pf_max - pf_min : 0;

                const size_t pvi = perm(g, vi, nv);
                const size_t pvj = perm(g, vj, nv);
                const size_t pvk = perm(g, vk, nv);

                const size_t value_i = value_min + ( pvi * value_max ) /
                                                   (nv+value_min);
                const size_t value_j = value_min + ( pvj * value_max ) /
                                                   (nv+value_min);
                const size_t value_k = value_min + ( pvk * value_max ) /
                                                   (nv+value_min);

                const GMBits2 bval_i = ((size_t)3) & value_i;
                const GMBits2 bval_j = ((size_t)3) & value_j;
                const GMBits2 bval_k = ((size_t)3) & value_k;

                const int bval_i_0 = !!(bval_i&1);
                const int bval_i_1 = !!(bval_i&2);
                const int bval_j_0 = !!(bval_j&1);
                const int bval_j_1 = !!(bval_j&2);
                const int bval_k_0 = !!(bval_k&1);
                const int bval_k_1 = !!(bval_k&2);

                si += ((bval_i_0 == i) + (bval_i_1 == i)) * gs_this;
                sj += ((bval_j_0 == j) + (bval_j_1 == j)) * gs_this;
                sk += ((bval_k_0 == k) + (bval_k_1 == k)) * gs_this;

                rijk += (((bval_i_0==i) && (bval_j_0==j) && (bval_k_0==k)) +
                         ((bval_i_1==i) && (bval_j_0==j) && (bval_k_0==k)) +
                         ((bval_i_0==i) && (bval_j_1==j) && (bval_k_0==k)) +
                         ((bval_i_1==i) && (bval_j_1==j) && (bval_k_0==k)) +
                         ((bval_i_0==i) && (bval_j_0==j) && (bval_k_1==k)) +
                         ((bval_i_1==i) && (bval_j_0==j) && (bval_k_1==k)) +
                         ((bval_i_0==i) && (bval_j_1==j) && (bval_k_1==k)) +
                         ((bval_i_1==i) && (bval_j_1==j) && (bval_k_1==k))) *
                        gs_this;

              } //---g

              //OLD const GMFloat rijk = i==1 && j==1 && k==1 ? (8 * num_field_f) : 0;

              //OLD const GMFloat si = i==1 ? (2 * num_field_f) : 0;
              //OLD const GMFloat sj = j==1 ? (2 * num_field_f) : 0;
              //OLD const GMFloat sk = k==1 ? (2 * num_field_f) : 0;

              const GMFloat value_expected =
                GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk, env);
              num_misses += value_expected != value;
            } //---k
          } //---j
        } //---i
      } //---for index
    } break;
    /*--------------------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
  do_->num_misses += num_misses;
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

      typedef unsigned char out_t;
      const GMFloat out_max = 255;

      const int buf_size = 4096;
      int buf_elts = 0;
      out_t buf[buf_size];
      size_t num_written_total = 0;

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
          if (file == stdin) {
            if (threshold < 0. || value > threshold) {
              fprintf(file,
                sizeof(GMFloat) == 8 ?
                "element (%li,%li): value: %.17e\n" :
                "element (%li,%li): value: %.8e\n", coord0, coord1, value);
            }
          } else {
            //if (threshold < 0. || value > threshold) {
            {
              out_t out_v = (out_t)(value * out_max);
              buf[buf_elts++] = out_v;
              if (buf_elts == buf_size) {
                size_t num_written = fwrite(&buf, sizeof(out_t),
                                            buf_elts, file);
                num_written_total += num_written;
                GMAssert(num_written == buf_elts*sizeof(out_t));
                buf_elts = 0;
              }
            }
          }
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

          if (file == stdin) {
            if (threshold < 0. || value > threshold) {
              fprintf(file,
                sizeof(GMFloat) == 8 ?
                "element (%li,%li): value: %.17e\n" :
                "element (%li,%li): value: %.8e\n", coord0, coord1, value);
            }
          } else {
            //if (threshold < 0. || value > threshold) {
            {
              out_t out_v = (out_t)(value * out_max);
              buf[buf_elts++] = out_v;
              if (buf_elts == buf_size) {
                size_t num_written = fwrite(&buf, sizeof(out_t),
                                            buf_elts, file);
                num_written_total += num_written;
                GMAssert(num_written == buf_elts*sizeof(out_t));
                buf_elts = 0;
              }
            }
          }
        }
      }

      if (file != stdin) {
        size_t num_written = fwrite(&buf, sizeof(out_t), buf_elts, file);
        num_written_total += num_written;
        GMAssert(num_written == buf_elts*sizeof(out_t));
        buf_elts = 0;
        printf("Wrote %lu elements of %lu from proc %i.\n",
               num_written_total, metrics->num_elts_local, GMEnv_proc_num(env));
      }

    } break;
    /*--------------------*/
    case GM_DATA_TYPE_TALLY2X2: {
    /*--------------------*/
//TODO: make faster, like above code
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
//TODO: make faster, like above code
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

GMChecksum perform_run(const char* const options) {
  GMAssertAlways(NULL != options);

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

#if 0
{
  env.metric_type_ = GM_METRIC_TYPE_CZEKANOWSKI;
  env.num_way_ = 3;
  env.all2all_ = 1;
  env.compute_method_ = GM_COMPUTE_METHOD_GPU;
  env.num_stage = 36;
  env.stage_num = 35;

  env.num_proc_world_ = 14880;
  env.num_proc_ = 14880;
  env.num_proc_field_ = 1;
  env.num_proc_repl_ = 496;
  env.num_proc_vector_i_ = 30;
  env.num_proc_vector_total_ = 30 * 496;

  env.proc_num_field_ = 1;

  env.is_proc_active_ = 1;

  //for (env.proc_num_=0; env.proc_num_<env.num_proc_; ++env.proc_num_) {
  for (env.proc_num_=6; env.proc_num_<7; ++env.proc_num_) {

    env.proc_num_repl_ = env.proc_num_ % env.num_proc_repl_;
    env.proc_num_vector_i_ = env.proc_num_ / env.num_proc_repl_;
    env.proc_num_vector_ = env.proc_num_;

    GMMetrics metrics = GMMetrics_null();

    //int nvl = 6324;

    GMMetrics_create(&metrics, GMEnv_data_type_metrics(&env),
                     100, 100,
                     6324,
                     6324, &env);
//    GMMetrics_3way_num_elts_local(&metrics, 6324, &env);
    printf("%lu\n", metrics.num_elts_local);
  }

exit(1);
}
#endif

  GMEnv_create(&env, MPI_COMM_WORLD, argc, argv, description);

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
  //do_.problem_type = GM_PROBLEM_TYPE_RANDOM;
  do_.problem_type = GM_PROBLEM_TYPE_ANALYTIC;

  do_.num_misses = 0;

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

  double outtime = 0;
  double mctime = 0;
  double cktime = 0;

  /*---Loop over stages---*/

  size_t num_elts_local_computed = 0;

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

    num_elts_local_computed += metrics.num_elts_local_computed;

    /*---Output results---*/

    time_beg = GMEnv_get_synced_time(&env);
    output_metrics(&metrics, &do_, &env);
    time_end = GMEnv_get_synced_time(&env);
    outtime += time_end - time_beg;

    /*---Check correctness---*/

    check_metrics(&metrics, &do_, &env);

    /*---Compute checksum---*/

    time_beg = GMEnv_get_synced_time(&env);
    GMMetrics_checksum(&metrics, &checksum, &env);
    time_end = GMEnv_get_synced_time(&env);
    cktime += time_end - time_beg;

    GMMetrics_destroy(&metrics, &env);
  }

  GMVectors_destroy(&vectors, &env);

  GMAssertAlways(env.cpu_mem == 0);
  GMAssertAlways(env.gpu_mem == 0);

  if (GMEnv_is_proc_active(&env)) {
  int mpi_code = 0;
    size_t num_elts_computed = 0;
    mpi_code = MPI_Allreduce(&num_elts_local_computed, &num_elts_computed, 1,
                             MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                             GMEnv_mpi_comm_vector(&env));
    GMAssertAlways(mpi_code == MPI_SUCCESS);

    if (GMEnv_num_way(&env) == GM_NUM_WAY_2 && GMEnv_all2all(&env)) {
      GMAssertAlways(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) / 2);
    }

    if (GMEnv_num_way(&env) == GM_NUM_WAY_3 && GMEnv_all2all(&env) &&
        do_.stage_min==1 && do_.stage_max==env.num_stage) {
      GMAssertAlways(num_elts_computed == (do_.num_vector) * (size_t)
                                          (do_.num_vector - 1) * (size_t)
                                          (do_.num_vector - 2) / 6);
    }
  }

  /*---Output run information---*/

  if (GMEnv_is_proc_active(&env) && GMEnv_proc_num(&env) == 0 &&
      do_.verbosity > 0) {
    //-----
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
    printf(" cktime %.6f", cktime);
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

  GMAssertAlways(do_.num_misses == 0);

  /*---Finalize---*/

  GMEnv_destroy(&env);

  return checksum;
}

/*---------------------------------------------------------------------------*/

#if 0
#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
} /*---extern "C"---*/
#endif
#endif

/*---------------------------------------------------------------------------*/
