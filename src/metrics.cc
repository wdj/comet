/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mpi.h"

#include "env.h"
#include "metrics.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Helper: round-robin-pack m values into n bins, give ith bin size---*/

static int rr_pack_(int i, int n, int m) {
  return m/n + (i < m % n ? 1 : 0);
}

/*===========================================================================*/
/*---Null object---*/

GMMetrics GMMetrics_null() {
  GMMetrics result;
  memset((void*)&result, 0, sizeof(GMMetrics));
  return result;
}

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      int num_field,
                      int num_vector_local,
                      GMEnv* env) {
  GMAssert(metrics);
  GMAssert(num_field >= 0);
  GMAssert(num_vector_local >= 0);
  GMAssert(env);

  *metrics = GMMetrics_null();

  if (!Env_is_proc_active(env)) {
    return;
  }

  GMInsist(env,
           num_field % Env_num_proc_field(env) == 0
               ? "num_proc_field must exactly divide the total number of fields"
               : 0);
  GMInsist(env, Env_all2all(env) || (Env_num_proc_vector_j(env) == 1 &&
                                     Env_num_proc_vector_k(env) == 1)
          ? "multidim parallelism only available for all2all case" : 0);
  GMInsist(env, Env_num_way(env) == GM_NUM_WAY_3 ||
                Env_num_proc_vector_k(env) == 1);

  metrics->data_type_id = data_type_id;
  metrics->num_field = num_field;
  metrics->num_field_local = num_field / Env_num_proc_field(env);
  metrics->num_vector_local = num_vector_local;
  metrics->index_offset_0_ = 0;
  metrics->index_offset_01_ = 0;
  metrics->recip_m = ((GMFloat)1) / num_field;

  /*---Compute global values---*/

  const int num_block = Env_num_block_vector(env);

  size_t num_vector_bound = metrics->num_vector_local * (size_t)num_block;
  num_vector_bound *= Env_num_proc_vector_j(env);
  num_vector_bound *= Env_num_proc_vector_k(env);
  GMAssert(num_vector_bound == (size_t)(int)num_vector_bound
               ? "Vector count too large to store in 32-bit int."
               : 0);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, Env_mpi_comm_vector(env));
  GMAssert(mpi_code == MPI_SUCCESS);
  GMAssert((size_t)(metrics->num_vector) == num_vector_bound);
  metrics->num_vector /= Env_num_proc_vector_j(env);
  metrics->num_vector /= Env_num_proc_vector_k(env);

  /*---Assume the following to simplify calculations---*/

  GMInsist(
      env,
      num_vector_local >= Env_num_way(env)
          ? "Currently require number of vecs on a proc to be at least num-way"
          : 0);

  const int i_block = Env_proc_num_vector_i(env);

  const int nchoosek = gm_nchoosek(num_vector_local, Env_num_way(env));

  int i = 0;
  int j = 0;
  int k = 0;

  /*---Compute number of elements etc.---*/

  /*==================================================*/
  if (Env_num_way(env) == GM_NUM_WAY_2 && Env_all2all(env)) {
  /*==================================================*/
    /*---Store the following in this block-row:
        1) strict upper triangular part of main diagonal block
        2) half of the off-diagonal blocks, as a "wrapped rectangle"
      For num_proc_vector_j > 1, map these blocks, starting at the
      main diagonal block, to procs in round-robin fashion.
    ---*/
    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_j = Env_proc_num_vector_j(env);
    const int num_proc_j = Env_num_proc_vector_j(env);
    metrics->num_elts_local = 0;
    /*---Calculate index size part 1: (triangle) i_block==j_block part---*/
    metrics->num_elts_local += proc_j == 0 ? nchoosek : 0;
    metrics->index_offset_0_ = metrics->num_elts_local;
    /*---Calculate index size part 2: (wrapped rect) i_block!=j_block part---*/
    /*---Total computed blocks this block row---*/
    const int num_block_this_slab_2 = num_block % 2 == 0 &&
                                      2 * i_block >= num_block
                                        ? (num_block / 2)
                                        : (num_block / 2) + 1;
    /*---Number stored for this proc_j---*/
    const int num_block_this_proc_2 = rr_pack_(proc_j, num_proc_j,
                                               num_block_this_slab_2);
    /*---Now count offdiag blocks only---*/
    const int num_offdiag_block = proc_j == 0 ? num_block_this_proc_2 - 1
                                              : num_block_this_proc_2;
    /*---Now put it all together---*/
    metrics->num_elts_local +=
        num_offdiag_block * num_vector_local * num_vector_local;
    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*===PART C: SET INDEX===*/
    /*---Set index part 1: (triangle) i_block==j_block part---*/
    size_t index = 0;
    if (proc_j == 0) {
      for (j = 0; j < num_vector_local; ++j) {
        const size_t j_global = j + num_vector_local * i_block;
        for (i = 0; i < j; ++i) {
          const size_t i_global = i + num_vector_local * i_block;
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
    }
    /*---Set index part 2: (wrapped rectangle) i_block!=j_block part---*/
    const size_t beg = (i_block + 1) * num_vector_local;
    const size_t end = (i_block + num_block_this_slab_2) * num_vector_local;
    size_t j_global_unwrapped = 0;
    for (j_global_unwrapped = beg; j_global_unwrapped < end;
         ++j_global_unwrapped) {
      const int j_block_unwrapped = (int)(j_global_unwrapped/num_vector_local);
      if ((j_block_unwrapped-i_block) % num_proc_j != proc_j) {
        continue;
      }
      const size_t j_global = j_global_unwrapped % metrics->num_vector;
      for (i = 0; i < num_vector_local; ++i) {
        const size_t i_global = i + num_vector_local * i_block;
        metrics->coords_global_from_index[index++] =
            i_global + metrics->num_vector * j_global;
      }
    }
    GMAssert(index == metrics->num_elts_local);
  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_3 && Env_all2all(env)) {
  /*==================================================*/
    /*---Make the following assumption to greatly simplify calculations---*/
    GMInsist(env, num_block <= 2 || metrics->num_vector_local % 6 == 0
                      ? "3way all2all case requires num vectors per proc "
                        "divisible by 6."
                      : 0);
    const int nchoosekm1 = gm_nchoosek(num_vector_local, Env_num_way(env) - 1);
    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_j = Env_proc_num_vector_j(env);
    const int num_proc_j = Env_num_proc_vector_j(env);
    metrics->num_elts_local = 0;
    int num_block_this_slab = 0;
    /*---Compute size pt 1: (tetrahedron) i_block==j_block==k_block part---*/

    const int num_block_this_slab_1 = 1;
    num_block_this_slab += num_block_this_slab_1;

    const int num_block_this_proc_1 = rr_pack_(proc_j, num_proc_j,
                                               num_block_this_slab);

    const int num_elts_per_block_1 = nchoosek;

    metrics->num_elts_local += num_block_this_proc_1 * num_elts_per_block_1;

    metrics->index_offset_0_ = metrics->num_elts_local;
    /*---Compute size pt 2: (triang prisms) i_block!=j_block==k_block part---*/

    const int num_block_this_slab_2 = num_block - 1;
    num_block_this_slab += num_block_this_slab_2;

    const int num_block_this_proc_12 = rr_pack_(proc_j, num_proc_j,
                                                num_block_this_slab);

    const int num_block_this_proc_2 = num_block_this_proc_12 -
                                      num_block_this_proc_1;

    const int num_elts_per_block_2 = nchoosekm1 * num_vector_local;

    metrics->num_elts_local += num_block_this_proc_2 * num_elts_per_block_2;

    metrics->index_offset_01_ = metrics->num_elts_local;
    /*---Compute size pt 3: (block sections) i_block!=j_block!=k_block part---*/

    const int num_block_this_slab_3 = (num_block - 1) * (num_block - 2);
    num_block_this_slab += num_block_this_slab_3;

    const int num_block_this_proc_123 = rr_pack_(proc_j, num_proc_j,
                                                 num_block_this_slab);

    const int num_block_this_proc_3 = num_block_this_proc_123 -
                                      num_block_this_proc_12;

    const int num_elts_per_block_3 = num_vector_local *
                               num_vector_local * (num_vector_local / 6);

    metrics->num_elts_local += num_block_this_proc_3 * num_elts_per_block_3;
    GMAssert(num_block_this_slab == (num_block-1) * (num_block-1) + 1);
    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*===PART C: SET INDEX===*/
    int block_this_slab = 0;
    /*---Set index part 1: (tetrahedron) i_block==j_block==k_block part---*/
    size_t index = 0;
    if (block_this_slab % num_proc_j == proc_j) {
      for (k = 0; k < num_vector_local; ++k) {
        const int k_block = i_block;
        const size_t k_global = k + num_vector_local * k_block;
        for (j = 0; j < k; ++j) {
          const int j_block = i_block;
          const size_t j_global = j + num_vector_local * j_block;
          for (i = 0; i < j; ++i) {
            const size_t i_global = i + num_vector_local * i_block;
            GMAssert(index < metrics->num_elts_local);
            metrics->coords_global_from_index[index++] =
                i_global +
                metrics->num_vector *
                    (j_global + metrics->num_vector * (k_global));
          }
        }
      }
    }
    block_this_slab += 1;
    /*---Set index part 2: (triang prisms) i_block!=j_block==k_block part---*/
    int j_block = 0;
    for (j_block = 0; j_block < num_block; ++j_block) {
      if (j_block == i_block) {
        continue;
      }
      if (block_this_slab % num_proc_j == proc_j) {
        const int k_block = j_block;
        for (k = 0; k < num_vector_local; ++k) {
          const size_t k_global = k + num_vector_local * k_block;
          for (j = 0; j < k; ++j) {
            const size_t j_global = j + num_vector_local * j_block;
            for (i = 0; i < num_vector_local; ++i) {
              const size_t i_global = i + num_vector_local * i_block;
              GMAssert(index < metrics->num_elts_local);
              metrics->coords_global_from_index[index++] =
                  i_global +
                  metrics->num_vector *
                      (j_global + metrics->num_vector * (k_global));
            }
          }
        }
      }
      block_this_slab += 1;
    } /*---k_block---*/
    /*---Set index part 3: (block sections) i_block!=j_block!=k_block part---*/
    int k_block = 0;
    for (k_block = 0; k_block < num_block; ++k_block) {
      if (k_block == i_block) {
        continue;
      }
      int j_block = 0;
      for (j_block = 0; j_block < num_block; ++j_block) {
        if (j_block == i_block || j_block == k_block) {
          continue;
        }
        if (block_this_slab % num_proc_j == proc_j) {
          const int section_axis = gm_metrics_3way_section_axis(
            metrics, i_block, j_block, k_block, env);
          const int section_num =
           gm_metrics_3way_section_num(metrics, i_block, j_block, k_block, env);
          const int inc = num_vector_local / 6;
          const int k_min = section_axis == 2 ? (section_num)*inc : 0;
          const int k_max =
              section_axis == 2 ? (section_num + 1) * inc : num_vector_local;
          for (k = k_min; k < k_max; ++k) {
            const size_t k_global = k + num_vector_local * k_block;
            const int j_min = section_axis == 1 ? (section_num)*inc : 0;
            const int j_max =
                section_axis == 1 ? (section_num + 1) * inc : num_vector_local;
            for (j = j_min; j < j_max; ++j) {
              const size_t j_global = j + num_vector_local * j_block;
              const int i_min = section_axis == 0 ? (section_num)*inc : 0;
              const int i_max = section_axis == 0 ? (section_num + 1) * inc
                                                  : num_vector_local;
              for (i = i_min; i < i_max; ++i) {
                const size_t i_global = i + num_vector_local * i_block;
                GMAssert(index < metrics->num_elts_local);
                metrics->coords_global_from_index[index++] =
                    i_global +
                    metrics->num_vector *
                        (j_global + metrics->num_vector * (k_global));
              }
            }
          }
        }
        block_this_slab += 1;
      } /*---j_block---*/
    }   /*---k_block---*/
    GMAssert(index == metrics->num_elts_local);
  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_2 && !Env_all2all(env)) {
  /*==================================================*/
    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*---Need store only strict upper triangular part of matrix---*/
    size_t index = 0;
    for (j = 0; j < num_vector_local; ++j) {
      const size_t j_global = j + num_vector_local * i_block;
      for (i = 0; i < j; ++i) {
        const size_t i_global = i + num_vector_local * i_block;
        metrics->coords_global_from_index[index++] =
            i_global + metrics->num_vector * j_global;
      }
    }
    GMAssert(index == metrics->num_elts_local);
  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_3 && !Env_all2all(env)) {
  /*==================================================*/
    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssert(metrics->coords_global_from_index != NULL);
    /*---Need store only strict interior of tetrahedron---*/
    size_t index = 0;
    for (k = 0; k < num_vector_local; ++k) {
      const size_t k_global = k + num_vector_local * i_block;
      for (j = 0; j < k; ++j) {
        const size_t j_global = j + num_vector_local * i_block;
        for (i = 0; i < j; ++i) {
          const size_t i_global = i + num_vector_local * i_block;
          GMAssert(index < metrics->num_elts_local);
          metrics->coords_global_from_index[index++] =
              i_global +
              metrics->num_vector *
                  (j_global + metrics->num_vector * (k_global));
        }
      }
    }
    GMAssert(index == metrics->num_elts_local);
  /*==================================================*/
  } else {
  /*==================================================*/
    GMInsist(env, 0 == 1 ? "Invalid set of options" : 0);
    /*---LATER: generalize this to N-way---*/
  }

  /*---Allocations---*/

  switch (data_type_id) {
    /*----------*/
    case GM_DATA_TYPE_BITS1: {
      /*---(design not complete)---*/
      const size_t num_floats_needed =
          gm_ceil_i8(metrics->num_elts_local, 8 * sizeof(GMFloat));
      metrics->data = malloc(num_floats_needed * sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      metrics->data_type_num_values = 1;
    } break;
    /*----------*/
    case GM_DATA_TYPE_FLOAT:
      //---TODO: remove thee casts on the mallocs, here and elsewhere.
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMFloat));
      GMAssert(metrics->data != NULL);
      metrics->data_type_num_values = 1;
      break;
    /*----------*/
    case GM_DATA_TYPE_TALLY2X2: {
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMTally2x2));
      GMAssert(metrics->data != NULL);
      metrics->data_M = malloc(metrics->num_elts_local * sizeof(GMFloat2));
      GMAssert(metrics->data_M != NULL);
      metrics->data_type_num_values = 4;
    } break;
    /*----------*/
    case GM_DATA_TYPE_TALLY4X2: {
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMTally4x2));
      GMAssert(metrics->data != NULL);
      metrics->data_M = malloc(metrics->num_elts_local * sizeof(GMFloat3));
      GMAssert(metrics->data_M != NULL);
      metrics->data_type_num_values = 8;
    } break;
    /*----------*/
    default:
      GMAssert(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics);
  GMAssert(env);
  GMAssert(metrics->data != NULL || !Env_is_proc_active(env));

  if (!Env_is_proc_active(env)) {
    return;
  }

  free(metrics->data);
  free(metrics->coords_global_from_index);
  if (metrics->data_M != NULL) {
    free(metrics->data_M);
  }
  *metrics = GMMetrics_null();
}

/*===========================================================================*/
/*---Functions for metrics checksum---*/

/*---------------------------------------------------------------------------*/
/*---Helper function - perform one bubble sort step---*/

static void gm_makegreater(size_t* i, size_t* j, int* ind_i, int* ind_j) {
  if (*i < *j) {
    const size_t tmp = *i;
    *i = *j;
    *j = tmp;
    const int tmp2 = *ind_i;
    *ind_i = *ind_j;
    *ind_j = tmp2;
  }
}

/*---------------------------------------------------------------------------*/
/*---Helper function - left shift that works for any shift amount---*/

static size_t gm_lshift(size_t a, int j) {
  if (j >= 64 || j <= -64) {
    return 0;
  }
  return j > 0 ? a << j : a >> (-j);
}

/*---------------------------------------------------------------------------*/
/*---Metrics checksum---*/

GMChecksum GMMetrics_checksum(GMMetrics* metrics, GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(metrics->data != NULL || !Env_is_proc_active(env));

  /*---Initializations---*/

  GMChecksum result = GMChecksum_null();

  if (!Env_is_proc_active(env)) {
    return result;
  }

  enum { num_way_max = GM_NUM_NUM_WAY + 1 };

  GMAssert(Env_num_way(env) <= num_way_max ? "This num_way not supported." : 0);

  /*---Initializations---*/

  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  UI64 sums_l[16];
  int i = 0;
  for (i = 0; i < 16; ++i) {
    sums_l[i] = 0;
  }
  double sum_d = 0;

  const int w = 30;
  GMAssertAlways(64 - 2 * w >= 4);
  const UI64 one64 = 1;

  const UI64 lomask = (one64 << w) - 1;
  const UI64 lohimask = (one64 << (2 * w)) - 1;

  UI64 coords[num_way_max];
  int ind_coords[num_way_max];
  UI64 index = 0;
  for (index = 0; index < metrics->num_elts_local; ++index) {
    /*---Obtain global coords of metrics elt---*/
    for (i = 0; i < num_way_max; ++i) {
      coords[i] = 0;
      ind_coords[i] = i;
    }
    for (i = 0; i < Env_num_way(env); ++i) {
      coords[i] = GMMetrics_coord_global_from_index(metrics, index, i, env);
    }
    /*---Reflect coords by symmetry to get uniform result -
         sort into descending order---*/
    gm_makegreater(&coords[1], &coords[2], &ind_coords[1], &ind_coords[2]);
    gm_makegreater(&coords[0], &coords[1], &ind_coords[0], &ind_coords[1]);
    gm_makegreater(&coords[1], &coords[2], &ind_coords[1], &ind_coords[2]);

    /*---Loop over data values at this index---*/
    int i_value = 0;
    for (i_value = 0; i_value < metrics->data_type_num_values; ++i_value) {
      /*---Construct global id for metrics data vbalue---*/
      UI64 uid = coords[0];
      for (i = 1; i < Env_num_way(env); ++i) {
        uid = uid * metrics->num_vector + coords[i];
      }
      uid = uid * metrics->data_type_num_values + i_value;
      /*---Randomize---*/
      const UI64 rand1 = gm_randomize(uid + 956158765);
      const UI64 rand2 = gm_randomize(uid + 842467637);
      UI64 rand_value = rand1 + gm_randomize_max() * rand2;
      rand_value &= lohimask;
      /*---Now pick up value of this metrics elt---*/
      GMFloat value = 0;
      switch (metrics->data_type_id) {
        /*--------------------*/
        case GM_DATA_TYPE_BITS1: {
          GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
        } break;
        /*--------------------*/
        case GM_DATA_TYPE_FLOAT: {
          value = GMMetrics_czekanowski_get_from_index(metrics, index, env);
        } break;
        /*--------------------*/
        case GM_DATA_TYPE_TALLY2X2: {
          const int i0_unpermuted = i_value / 2;
          const int i1_unpermuted = i_value % 2;
          const int i0 = ind_coords[0] == 0 ? i0_unpermuted : i1_unpermuted;
          const int i1 = ind_coords[0] == 0 ? i1_unpermuted : i0_unpermuted;
          value = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
        } break;
        /*--------------------*/
        case GM_DATA_TYPE_TALLY4X2: {
          const int i0_unpermuted = i_value / 4;
          const int i1_unpermuted = (i_value / 2) % 2;
          const int i2_unpermuted = i_value % 2;
          // FIX: make sure this permutation direction correct.



#if 0
          const int i0 = ind_coords[0] == 0 ? i0_unpermuted :
                         ind_coords[0] == 1 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i1 = ind_coords[1] == 0 ? i0_unpermuted :
                         ind_coords[1] == 1 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i2 = ind_coords[2] == 0 ? i0_unpermuted :
                         ind_coords[2] == 1 ? i1_unpermuted :
                                              i2_unpermuted;
#endif

          const int i0 = ind_coords[0] == 0 ? i0_unpermuted :
                         ind_coords[1] == 0 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i1 = ind_coords[0] == 1 ? i0_unpermuted :
                         ind_coords[1] == 1 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i2 = ind_coords[0] == 2 ? i0_unpermuted :
                         ind_coords[1] == 2 ? i1_unpermuted :
                                              i2_unpermuted;


          value =
              GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2, env);
        } break;
        /*--------------------*/
        default:
          GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
      } /*---switch---*/
      /*---Convert to integer.  Store only 2*w bits max---*/
      const int log2_value_max = 4;
      GMAssertAlways(value >= 0 && value < (1 << log2_value_max));
      UI64 ivalue = value * (one64 << (2 * w - log2_value_max));
      /*---Multiply the two values---*/
      const UI64 a = rand_value;
      const UI64 alo = a & lomask;
      const UI64 ahi = a >> w;
      const UI64 b = ivalue;
      const UI64 blo = b & lomask;
      const UI64 bhi = b >> w;
      const UI64 cx = alo * bhi + ahi * blo;
      UI64 clo = alo * blo + ((cx & lomask) << w);
      UI64 chi = ahi * bhi + (cx >> w);
      const double value_d =
          ivalue * (double)rand_value / ((double)(one64 << (2 * w)));
      sum_d += value_d;
      /*---(move the carry bits)---*/
      chi += clo >> (2 * w);
      clo &= lohimask;
      /*---Split the product into one-char chunks, accumulate to sums---*/
      for (i = 0; i < 8; ++i) {
        sums_l[0 + i] += (clo << (64 - 8 - 8 * i)) >> (64 - 8);
        sums_l[8 + i] += (chi << (64 - 8 - 8 * i)) >> (64 - 8);
      }
    } /*---for i_value---*/
  }   /*---for index---*/
  /*---Global sum---*/
  UI64 sums_g[16];
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(sums_l, sums_g, 16, MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           Env_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  /*---Combine results---*/

  for (i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    int j = 0;
    for (j = 0; j < 8; ++j) {
      result.data[i] += gm_lshift(sums_g[0 + j], 8 * j - 2 * w * i) & lohimask;
      result.data[i] +=
          gm_lshift(sums_g[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
    }
  }
  /*---(move the carry bits---*/
  result.data[1] += result.data[0] >> (2 * w);
  result.data[0] &= lohimask;
  result.data[2] += result.data[1] >> (2 * w);
  result.data[1] &= lohimask;

  /*---Check against floating point result---*/

  const double tmp = sum_d;
  mpi_code = MPI_Allreduce(&tmp, &sum_d, 1, MPI_DOUBLE, MPI_SUM,
                           Env_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  double result_d = result.data[0] / ((double)(one64 << (2 * w))) +
                    result.data[1] +
                    result.data[2] * ((double)(one64 << (2 * w)));
  result_d = 1 * result_d; /*---Avoid unused variable warning---*/
  GMAssertAlways(fabs(sum_d - result_d) <= sum_d * 1.e-10);

  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
