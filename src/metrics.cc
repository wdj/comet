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
/*---Helper: round-robin-pack m values into n bins, return ith bin size---*/

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
  GMAssertAlways(metrics);
  GMAssertAlways(num_field >= 0);
  GMAssertAlways(num_vector_local >= 0);
  GMAssertAlways(env);

  *metrics = GMMetrics_null();

  if (!Env_is_proc_active(env)) {
    return;
  }

  GMInsist(env,
           num_field % Env_num_proc_field(env) == 0
               ? "num_proc_field must exactly divide the total number of fields"
               : 0);
  GMInsist(env, Env_all2all(env) || Env_num_proc_repl(env) == 1
          ? "multidim parallelism only available for all2all case" : 0);

  int i = 0;
  int j = 0;
  int k = 0;

  metrics->data_type_id = data_type_id;
  metrics->num_field = num_field;
  metrics->num_field_local = num_field / Env_num_proc_field(env);
  metrics->num_vector_local = num_vector_local;
  metrics->nvl6 = num_vector_local / 6;
  metrics->index_offset_0_ = 0;
  metrics->index_offset_01_ = 0;
  metrics->recip_m = ((GMFloat)1) / num_field;
  for (i=0; i<6; ++i) {
    metrics->index_offset_section_pt1_[i] = 0;
    metrics->index_offset_section_pt2_[i] = 0;
    metrics->section_num_valid_pt1_[i] = GM_BOOL_FALSE;
    metrics->section_num_valid_pt2_[i] = GM_BOOL_FALSE;
  }

  /*---Compute global values---*/

  const int num_block = Env_num_block_vector(env);

  size_t num_vector_bound = metrics->num_vector_local * (size_t)num_block;
  num_vector_bound *= Env_num_proc_repl(env);
  GMAssertAlways(num_vector_bound == (size_t)(int)num_vector_bound
               ? "Vector count too large to store in 32-bit int."
               : 0);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, Env_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
  GMAssertAlways((size_t)(metrics->num_vector) == num_vector_bound);
  metrics->num_vector /= Env_num_proc_repl(env);

  /*---Assume the following to simplify calculations---*/

  GMInsist(
      env,
      num_vector_local >= Env_num_way(env)
          ? "Currently require number of vecs on a proc to be at least num-way"
          : 0);

  const int i_block = Env_proc_num_vector_i(env);

  const int nchoosek = gm_nchoosek(num_vector_local, Env_num_way(env));

  /*---Compute number of elements etc.---*/

  /*==================================================*/
  if (Env_num_way(env) == GM_NUM_WAY_2 && Env_all2all(env)) {
  /*==================================================*/
    /*---Store the following in this block-row:
        1) strict upper triangular part of main diagonal block
        2) half of the off-diagonal blocks, as a "wrapped rectangle"
      For num_proc_repl > 1, map these blocks, starting at the
      main diagonal block, to procs in round-robin fashion.
    ---*/
    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_num_r = Env_proc_num_repl(env);
    const int num_proc_r = Env_num_proc_repl(env);
    metrics->num_elts_local = 0;
    /*---Calculate index size part 1: (triangle) i_block==j_block part---*/
    metrics->num_elts_local += proc_num_r == 0 ? nchoosek : 0;
    metrics->index_offset_0_ = metrics->num_elts_local;
    /*---Calculate index size part 2: (wrapped rect) i_block!=j_block part---*/
    /*---Total computed blocks this block row---*/
    const int num_block_this_slab_2 = num_block % 2 == 0 &&
                                      2 * i_block >= num_block
                                        ? (num_block / 2)
                                        : (num_block / 2) + 1;
    /*---Number stored for this proc_num_r---*/
    const int num_block_this_proc_2 = rr_pack_(proc_num_r, num_proc_r,
                                               num_block_this_slab_2);
    /*---Now count offdiag blocks only---*/
    const int num_offdiag_block = proc_num_r == 0 ? num_block_this_proc_2 - 1
                                                  : num_block_this_proc_2;
    /*---Now put it all together---*/
    metrics->num_elts_local +=
        num_offdiag_block * num_vector_local * num_vector_local;
    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssertAlways(metrics->coords_global_from_index != NULL);
    /*===PART C: SET INDEX===*/
    /*---Set index part 1: (triangle) i_block==j_block part---*/
    size_t index = 0;
    if (proc_num_r == 0) {
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
      if ((j_block_unwrapped-i_block) % num_proc_r != proc_num_r) {
        continue;
      }
      const size_t j_global = j_global_unwrapped % metrics->num_vector;
      for (i = 0; i < num_vector_local; ++i) {
        const size_t i_global = i + num_vector_local * i_block;
        metrics->coords_global_from_index[index++] =
            i_global + metrics->num_vector * j_global;
      }
    }
    GMAssertAlways(index == metrics->num_elts_local);
  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_3 && Env_all2all(env)) {
  /*==================================================*/

    /*---Make the following assumption to greatly simplify calculations---*/
    GMInsist(env, num_block <= 2 || metrics->num_vector_local % 6 == 0
                      ? "3way all2all case requires num vectors per proc "
                        "divisible by 6."
                      : 0);
    const int nchoosekm1 = gm_nchoosek(num_vector_local, Env_num_way(env) - 1);
    const int num_section_steps = GMEnv_num_section_steps(env);
    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_num_r = Env_proc_num_repl(env);
    const int num_proc_r = Env_num_proc_repl(env);
    const int nvl6 = metrics->nvl6;
    metrics->num_elts_local = 0;
    int num_block_this_slab = 0;
    int num_block_this_proc = 0;

    int section_block_num = 0;




#if 0
convert this to loop form instead of formulaic use of rr_pack

GMSectionInfo should provide num_section_steps, loop over this
  - if not all2all or num_proc_repl==1 then 1, else 6 (only used if 3-way)
  - for this part and metrics.h, may be ok to always assume 6

if proc_repl == section_step, then assign section to this proc

compute section size formulaically

#endif

    /*---Compute size pt 1: (tetrahedron) i_block==j_block==k_block part---*/

    if (num_section_steps == 1) {
      metrics->index_offset_section_pt1_[0] = metrics->num_elts_local;
      metrics->section_num_valid_pt1_[0] = GM_BOOL_TRUE;
      const int num_block_this_slab_1 = 1;
      num_block_this_slab += num_block_this_slab_1;
      const int num_block_this_proc_1 = rr_pack_(proc_num_r, num_proc_r,
                                                 num_block_this_slab);
      num_block_this_proc = num_block_this_proc_1;
      const int num_elts_per_block_1 = nchoosek;
      const int num_elts_local_pt1 = num_block_this_proc_1 *
                                     num_elts_per_block_1;
      metrics->num_elts_local += num_elts_local_pt1;
      metrics->index_offset_0_ = metrics->num_elts_local;
      ++section_block_num;
    } else {
      int section_step = 0;
      for (section_step=0; section_step<num_section_steps; ++section_step) {
        const int section_num = section_step;
        metrics->index_offset_section_pt1_[section_num]
          = metrics->num_elts_local;
        metrics->section_num_valid_pt1_[section_num] = GM_BOOL_TRUE;
        if (section_block_num % num_proc_r == proc_num_r) {

// set metrics->num_elts_local - slice of trapezoid
// should metrics->index_offset_section_pt1_ have other stuff ...

        }
        ++section_block_num;
      }
    }

    /*---Compute size pt 2: (triang prisms) i_block!=j_block==k_block part---*/

    if (num_section_steps == 1) {
      metrics->index_offset_section_pt2_[0] = metrics->num_elts_local;
      metrics->section_num_valid_pt2_[0] = GM_BOOL_TRUE;
      const int num_block_this_slab_2 = num_block - 1;
      num_block_this_slab += num_block_this_slab_2;
      const int num_block_this_proc_12 = rr_pack_(proc_num_r, num_proc_r,
                                                  num_block_this_slab);
      const int num_block_this_proc_2 = num_block_this_proc_12 -
                                        num_block_this_proc;
      num_block_this_proc = num_block_this_proc_12;
      const int num_elts_per_block_2 = nchoosekm1 * num_vector_local;
      const int num_elts_local_pt2 = num_block_this_proc_2 *
                                     num_elts_per_block_2;
      metrics->num_elts_local += num_elts_local_pt2;;
      metrics->index_offset_01_ = metrics->num_elts_local;
    } else {
      int section_step = 0;
      for (section_step=0; section_step<num_section_steps; ++section_step) {



      }
    }

    /*---Compute size pt 3: (block sections) i_block!=j_block!=k_block part---*/

    if (num_section_steps == 1) {
      const int num_block_this_slab_3 = (num_block - 1) * (num_block - 2);
      num_block_this_slab += num_block_this_slab_3;
      const int num_block_this_proc_123 = rr_pack_(proc_num_r, num_proc_r,
                                                   num_block_this_slab);
      const int num_block_this_proc_3 = num_block_this_proc_123 -
                                        num_block_this_proc;
      num_block_this_proc = num_block_this_proc_123;
      const int num_elts_per_block_3 = num_vector_local *
                                 num_vector_local * nvl6;
      metrics->num_elts_local += num_block_this_proc_3 * num_elts_per_block_3;
    } else {
      int section_step = 0;
      for (section_step=0; section_step<num_section_steps; ++section_step) {



      }
    }

    GMAssertAlways(num_block_this_slab == (num_block-1) * (num_block-1) + 1);

    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssertAlways(metrics->coords_global_from_index != NULL);



    /*===PART C: SET INDEX===*/
    int block_num = 0;
    int block_num_this_r = 0;
    size_t index = 0;



    /*---Set index part 1: (tetrahedron) i_block==j_block==k_block part---*/
    if (block_num % num_proc_r == proc_num_r) {
      for (j = 0; j < num_vector_local; ++j) {
        const int j_block = i_block;
        const size_t j_global = j + num_vector_local * j_block;
        for (k = j+1; k < num_vector_local; ++k) {
          const int k_block = i_block;
          const size_t k_global = k + num_vector_local * k_block;
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
      block_num_this_r += 1;
    } /*---if block_num---*/
    block_num += 1;
    //metrics->block_num_offset_0_ = block_num_this_r;


    /*---Set index part 2: (triang prisms) i_block!=j_block==k_block part---*/
    int j_i_block_delta = 0;
    for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta) {
      const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);
      if (block_num % num_proc_r == proc_num_r) {
        for (j = 0; j < num_vector_local; ++j) {
          const size_t j_global = j + num_vector_local * j_block;
          const int k_block = j_block;
          for (k = j+1; k < num_vector_local; ++k) {
            const size_t k_global = k + num_vector_local * k_block;
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
        block_num_this_r += 1;
      } /*---if block_num---*/
      block_num += 1;
    } /*---k_block---*/
    //metrics->block_num_offset_01_ = block_num_this_r;


    /*---Set index part 3: (block sections) i_block!=j_block!=k_block part---*/
    int k_i_block_delta = 0;
    for (k_i_block_delta = 1; k_i_block_delta < num_block; ++k_i_block_delta) {
      const int k_block = gm_mod_i(i_block + k_i_block_delta, num_block);
      int j_i_block_delta = 0;
      for (j_i_block_delta = 1; j_i_block_delta < num_block; ++j_i_block_delta){
        const int j_block = gm_mod_i(i_block + j_i_block_delta, num_block);
        if (j_block == k_block) {
          continue;
        }
        if (block_num % num_proc_r == proc_num_r) {
          const int section_axis = gm_section_axis_part3(i_block, j_block,
                                                         k_block);
          const int section_num = gm_section_num_part3(i_block, j_block,
                                                       k_block);
          const int j_min = section_axis == 1 ? (section_num)*nvl6 : 0;
          const int j_max =
              section_axis == 1 ? (section_num + 1) * nvl6 : num_vector_local;
          for (j = j_min; j < j_max; ++j) {
            const size_t j_global = j + num_vector_local * j_block;
            const int k_min = section_axis == 2 ? (section_num)*nvl6 : 0;
            const int k_max =
                section_axis == 2 ? (section_num + 1) * nvl6 : num_vector_local;
            for (k = k_min; k < k_max; ++k) {
              const size_t k_global = k + num_vector_local * k_block;
              const int i_min = section_axis == 0 ? (section_num)*nvl6 : 0;
              const int i_max = section_axis == 0 ? (section_num + 1) * nvl6
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
          block_num_this_r += 1;
        } /*---if block_num---*/
        block_num += 1;
      } /*---j_block---*/
    }   /*---k_block---*/

    GMAssertAlways(index == metrics->num_elts_local);
    GMAssertAlways(block_num == (num_block-1) * (num_block-1) + 1);





  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_2 && !Env_all2all(env)) {
  /*==================================================*/
    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssertAlways(metrics->coords_global_from_index != NULL);
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
    GMAssertAlways(index == metrics->num_elts_local);
  /*==================================================*/
  } else if (Env_num_way(env) == GM_NUM_WAY_3 && !Env_all2all(env)) {
  /*==================================================*/
    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)malloc(metrics->num_elts_local * sizeof(size_t));
    GMAssertAlways(metrics->coords_global_from_index != NULL);
    /*---Need store only strict interior of tetrahedron---*/
    size_t index = 0;
    for (j = 0; j < num_vector_local; ++j) {
      const int j_block = i_block;
      const size_t j_global = j + num_vector_local * j_block;
      for (k = j+1; k < num_vector_local; ++k) {
        const int k_block = i_block;
        const size_t k_global = k + num_vector_local * k_block;
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
    GMAssertAlways(index == metrics->num_elts_local);
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
      GMAssertAlways(metrics->data != NULL);
      metrics->data_type_num_values = 1;
    } break;
    /*----------*/
    case GM_DATA_TYPE_FLOAT:
      //---TODO: remove thee casts on the mallocs, here and elsewhere.
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMFloat));
      GMAssertAlways(metrics->data != NULL);
      metrics->data_type_num_values = 1;
      break;
    /*----------*/
    case GM_DATA_TYPE_TALLY2X2: {
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMTally2x2));
      GMAssertAlways(metrics->data != NULL);
      metrics->data_M = malloc(metrics->num_elts_local * sizeof(GMFloat2));
      GMAssertAlways(metrics->data_M != NULL);
      metrics->data_type_num_values = 4;
    } break;
    /*----------*/
    case GM_DATA_TYPE_TALLY4X2: {
      metrics->data = malloc(metrics->num_elts_local * sizeof(GMTally4x2));
      GMAssertAlways(metrics->data != NULL);
      metrics->data_M = malloc(metrics->num_elts_local * sizeof(GMFloat3));
      GMAssertAlways(metrics->data_M != NULL);
      metrics->data_type_num_values = 8;
    } break;
    /*----------*/
    default:
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMAssertAlways(metrics);
  GMAssertAlways(env);
  GMAssertAlways(metrics->data != NULL || !Env_is_proc_active(env));

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
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(metrics->data != NULL || !Env_is_proc_active(env));

  /*---Initializations---*/

  GMChecksum result = GMChecksum_null();

  if (!Env_is_proc_active(env)) {
    return result;
  }

  enum { num_way_max = GM_NUM_NUM_WAY + 1 };

  GMAssertAlways(Env_num_way(env) <= num_way_max ? "This num_way not supported." : 0);

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
      GMAssert(value >= 0);
      GMAssert(value < (1 << log2_value_max));
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
