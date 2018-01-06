//-----------------------------------------------------------------------------
/*!
 * \file   metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "metrics.hh"

//=============================================================================
/*---Helper: round-robin-pack m values into n bins, return ith bin size---*/

static int rr_pack_(int i, int n, int m) {
  return m/n + (i < m % n ? 1 : 0);
}

//=============================================================================
/*---Null object---*/

GMMetrics GMMetrics_null() {
  GMMetrics result;
  memset((void*)&result, 0, sizeof(result));
  return result;
}

//=============================================================================

void GMMetrics_3way_num_elts_local(GMMetrics* metrics, int nvl,
                                   GMEnv* env) {
  GMInsist(metrics && env);
  GMInsist(nvl >= 0);
  GMInsist(GMEnv_num_block_vector(env) <= 2 || nvl % 6 == 0);
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_3);

  metrics->num_elts_local = 0;

  //---Fused counter for section_num and block_num, same across all procs.
  int section_block_num = 0;

  //---Compute size part 1: (tetrahedron) i_block==j_block==k_block part.

  const int num_section_steps_1 = gm_num_section_steps(env, 1);
  for (int section_step=0; section_step<num_section_steps_1; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 1, env);
    const int J_hi = gm_J_hi(section_num, nvl, 1, env);
    const GMInt64 trap_size_lo = gm_trap_size(J_lo, nvl);
    const GMInt64 trap_size_hi = gm_trap_size(J_hi, nvl);
    GMInsist(trap_size_hi >= trap_size_lo);
    //---Absorb size_lo into offset for speed in indexing function.
    metrics->index_offset_section_part1_[section_num]
      = (GMInt64)metrics->num_elts_local - trap_size_lo;
    if (gm_proc_r_active(section_block_num, env)) {
      //---Elements in slice of trapezoid.
      const GMInt64 elts_local = trap_size_hi - trap_size_lo;
      GMInsist(elts_local >= 0);
      metrics->num_elts_local += elts_local;
      metrics->section_num_valid_part1_[section_num] = (elts_local != 0);
    }
    ++section_block_num;
  }
  metrics->index_offset_0_ = metrics->num_elts_local;

  //---Compute size part 2: (triang prisms) i_block!=j_block==k_block part.

  const int num_block = GMEnv_num_block_vector(env);
  const int num_section_steps_2 = gm_num_section_steps(env, 2);
  for (int section_step=0; section_step<num_section_steps_2; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 2, env);
    const int J_hi = gm_J_hi(section_num, nvl, 2, env);
    const GMInt64 triang_size_lo = gm_triang_size(J_lo, nvl);
    const GMInt64 triang_size_hi = gm_triang_size(J_hi, nvl);
    //---Absorb size_lo into offset for speed in indexing function.
    metrics->index_offset_section_part2_[section_num]
      = (GMInt64)metrics->num_elts_local - (GMInt64)nvl*(GMInt64)triang_size_lo;
    metrics->section_size_part2[section_num] = triang_size_hi -
                                               triang_size_lo;
    //---Loop over blocks for part2.
    for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
      if (gm_proc_r_active(section_block_num, env)) {
        //---Elements in slice of triang prism.
        const GMInt64 elts_local = (GMInt64)nvl *
                                   (triang_size_hi - triang_size_lo);
        GMInsist(elts_local >= 0);
        metrics->num_elts_local += elts_local;
        metrics->section_num_valid_part2_[section_num] = (elts_local != 0);
      }
      ++section_block_num;
    }
  }
  metrics->index_offset_01_ = metrics->num_elts_local;

  //---Compute size part 3: (block sections) i_block!=j_block!=k_block part.

  //---Loop over block for part3.
  const int i_block = GMEnv_proc_num_vector_i(env);
  for (int k_i_offset=1; k_i_offset<num_block; ++k_i_offset) {
    const int k_block = gm_mod_i(i_block + k_i_offset, num_block);
    for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
      const int j_block = gm_mod_i(i_block + j_i_offset, num_block);
      if (j_block == k_block) {
        continue;
      }
      //---Get slice bounds.
      const int section_num = gm_section_num_part3(i_block, j_block, k_block);
      const int J_lo = gm_J_lo(section_num, nvl, 3, env);
      const int J_hi = gm_J_hi(section_num, nvl, 3, env);
      if (gm_proc_r_active(section_block_num, env)) {
        //---Elements in slice of block/cube.
        const GMInt64 elts_local = (GMInt64)nvl * (GMInt64)nvl *
                                   (GMInt64)(J_hi - J_lo);
        GMInsist(elts_local >= 0);
        metrics->num_elts_local += elts_local;
      }
      ++section_block_num;
    }
  }
}

//=============================================================================
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      GMDecompMgr* dm,
                      GMEnv* env) {
  GMInsist(metrics && dm && env);

  *metrics = GMMetrics_null();

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  GMInsistInterface(env, (GMEnv_all2all(env) || GMEnv_num_proc_repl(env) == 1)
          && "multidim parallelism only available for all2all case");

  /*---The following less important cases are not yet tested---*/

  GMInsistInterface(env, (GMEnv_all2all(env) ||
                    dm->num_field == dm->num_field_active)
                    && "This case currently not supported.");

  GMInsistInterface(env, (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU ||
                    dm->num_field == dm->num_field_active)
                    && "This case currently not supported.");

  metrics->dm = dm;
  metrics->data_type_id = data_type_id;

  metrics->num_field = dm->num_field;
  metrics->num_field_active = dm->num_field_active;
  metrics->num_field_local = dm->num_field_local;
  metrics->num_vector_local = dm->num_vector_local;
  metrics->num_vector_active = dm->num_vector_active;

  metrics->nvl6 = metrics->num_vector_local / 6;
  metrics->index_offset_0_ = 0;
  metrics->index_offset_01_ = 0;
  metrics->recip_m = ((GMFloat)1) / metrics->num_field_active;
  metrics->block_min = 0;
  for (int i=0; i<6; ++i) {
    metrics->index_offset_section_part1_[i] = 0;
    metrics->index_offset_section_part2_[i] = 0;
    metrics->section_num_valid_part1_[i] = false;
    metrics->section_num_valid_part2_[i] = false;
  }

  metrics->num_vector = dm->num_vector;

  const int num_block = GMEnv_num_block_vector(env);
  int mpi_code = 0;

  GMInsistInterface(
      env,
      metrics->num_vector_local >= GMEnv_num_way(env)
        && "Currently require number of vecs on a proc to be at least num-way");

  const int i_block = GMEnv_proc_num_vector_i(env);

  const size_t nchoosek = gm_nchoosek(metrics->num_vector_local,
                                      GMEnv_num_way(env));
  const int nvl = metrics->num_vector_local;
  const size_t nvlsq = nvl * (size_t)nvl;

  /*---Compute number of elements etc.---*/

  GMInsistInterface(env, env->stage_num >= 0 && env->stage_num < env->num_stage
                && "Invalid stage number specified.");

  GMInsistInterface(env, env->phase_num >= 0 && env->phase_num < env->num_phase
                && "Invalid phase number specified.");

  metrics->num_elts_local_computed = 0;

  /*==================================================*/
  if (GMEnv_num_way(env) == GM_NUM_WAY_2 && GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsistInterface(env, env->num_stage == 1
                      && "Staged computations not allowed for 2-way case.");

    GMInsistInterface(env, env->num_phase <= 1 + num_block / 2
                      && "num_phase must be at most 1 + nproc_vector/2.");

    /*---Store the following in this block-row:
        1) strict upper triangular part of main diagonal block
        2) half of the off-diagonal blocks, as a "wrapped rectangle"
      For num_proc_repl > 1, map these blocks, starting at the
      main diagonal block, to procs in round-robin fashion.
      For num_phase > 1, do all this only for a piece of the block row.
    ---*/

    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_num_r = GMEnv_proc_num_repl(env);
    const int num_proc_r = GMEnv_num_proc_repl(env);
    metrics->num_elts_local = 0;

    /*---PART A.1: (triangle) i_block==j_block part---*/
    const bool have_main_diag = proc_num_r == 0 &&
                                gm_diag_computed_min(env) == 0;
    metrics->num_elts_local += have_main_diag ? nchoosek : 0;
    metrics->index_offset_0_ = have_main_diag ? nchoosek - nvlsq : 0;
    metrics->block_min = (i_block + gm_diag_computed_min(env)) % num_block;


    /*---PART A.2: (wrapped rect) i_block!=j_block part---*/
    const int num_computed_blocks_this_row = gm_computed_blocks_this_row(env);
    const int num_computed_blocks_this_proc = rr_pack_(proc_num_r, num_proc_r,
                                               num_computed_blocks_this_row);
    const int num_computed_offdiag_blocks_this_proc =
      num_computed_blocks_this_proc - (have_main_diag ? 1 : 0);
    metrics->num_elts_local += num_computed_offdiag_blocks_this_proc * nvlsq;

    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
//if (env->proc_num_base_ == 0) printf("Hey1\n");

    /*===PART C: SET INDEX===*/

    /*---PART C.1: (triangle) i_block==j_block part---*/
    size_t index = 0;
    if (have_main_diag) {
      for (int j = 0; j < nvl; ++j) {
        const size_t j_global = j + nvl * i_block;
        for (int i = 0; i < j; ++i) {
          const size_t i_global = i + nvl * i_block;
          GMAssert(GMMetrics_helper2way_maindiag_block_(metrics, i, j, i_block,
                                                        env) == index);
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
    }

    /*---PART C.2: (wrapped rectangle) i_block!=j_block part---*/

    const int beg = gm_diag_computed_min(env);
    const int end = beg + num_computed_blocks_this_row;
    for (int diag=beg; diag<end; ++diag) {
      if (diag == 0 || ! gm_proc_r_active(diag-beg, env)) {
        continue;
      }
      const int j_block_unwrapped = i_block + diag;
      #pragma omp parallel for collapse(2) schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
        const size_t j_global_unwrapped = j + j_block_unwrapped * (size_t)nvl;
        const size_t j_global = j_global_unwrapped % metrics->num_vector;
          const size_t i_global = i + nvl * i_block;
          const size_t index_this = index + i + j * (size_t)nvl;
          GMAssert(index_this>=0 && index_this<metrics->num_elts_local);
          metrics->coords_global_from_index[index_this] =
              i_global + metrics->num_vector * j_global;
        }
      }
      index += nvlsq;
    } /*---for diag---*/

    /*---Final check---*/
    GMInsist(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_3 && GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsistInterface(env, env->num_phase == 1
                      && "Phaseed computations not currently implemented "
                        "for 3-way case.");

    /*---Make the following assumption to greatly simplify calculations---*/
    GMInsistInterface(env, (num_block<=2 || metrics->num_vector_local % 6 == 0)
                      && "3way all2all case requires num vectors per proc "
                        "divisible by 6.");

    /*===PART A: CALCULATE INDEX SIZE===*/

    GMMetrics_3way_num_elts_local(metrics, nvl, env);

    /*---Fused counter for section_num and block_num, same across all procs---*/
    int section_block_num = 0;

    const int num_section_steps_12 = gm_num_section_steps(env, 1);

    /*===PART B: ALLOCATE INDEX===*/

    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);

    /*===PART C: SET INDEX===*/

    section_block_num = 0;
    size_t index = 0;

    /*---Set index part 1: (tetrahedron) i_block==j_block==k_block part---*/

    for (int section_step=0; section_step<num_section_steps_12; ++section_step){
      if (gm_proc_r_active(section_block_num, env)) {
        const int section_num = section_step;
        const int J_lo = gm_J_lo(section_num, nvl, 1, env);
        const int J_hi = gm_J_hi(section_num, nvl, 1, env);
        const int j_min = J_lo;
        const int j_max = J_hi;
        for (int j = j_min; j < j_max; ++j) {
          const int j_block = i_block;
          const size_t j_global = j + nvl * j_block;
          #pragma omp parallel for collapse(2) schedule(dynamic,1000)
          for (int k = j+1; k < nvl; ++k) {
            for (int i = 0; i < j; ++i) {
            const int k_block = i_block;
            const size_t k_global = k + nvl * k_block;
              const size_t i_global = i + nvl * i_block;
              const size_t index_this = index + i + j*(size_t)(k-(j+1));
              GMAssert(index_this>=0 && index_this<metrics->num_elts_local);
              metrics->coords_global_from_index[index_this] =
                  i_global +
                metrics->num_vector *
                      (j_global + metrics->num_vector * (k_global));
            }
          }
          index += j * (size_t)(nvl - (j+1));
        }
      } /*---if block_num---*/
      ++section_block_num;
    }

    /*---Set index part 2: (triang prisms) i_block!=j_block==k_block part---*/

    for (int section_step=0; section_step<num_section_steps_12; ++section_step){
      for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
        const int j_block = gm_mod_i(i_block + j_i_offset, num_block);
        if (gm_proc_r_active(section_block_num, env)) {
          const int section_num = section_step;
          const int J_lo = gm_J_lo(section_num, nvl, 2, env);
          const int J_hi = gm_J_hi(section_num, nvl, 2, env);
          const int j_min = J_lo;
          const int j_max = J_hi;
          for (int j = j_min; j < j_max; ++j) {
            const size_t j_global = j + nvl * j_block;
            #pragma omp parallel for collapse(2) schedule(dynamic,1000)
            for (int k = j+1; k < nvl; ++k) {
              for (int i = 0; i < nvl; ++i) {
              const int k_block = j_block;
              const size_t k_global = k + nvl * k_block;
                const size_t i_global = i + nvl * i_block;
                const size_t index_this = index + i + nvl*(size_t)(k-(j+1));
                GMAssert(index_this>=0 && index_this<metrics->num_elts_local);
                metrics->coords_global_from_index[index_this] =
                    i_global +
                  metrics->num_vector *
                        (j_global + metrics->num_vector * (k_global));
              }
            }
            index += nvl * (size_t)(nvl - (j+1));
          }
        } /*---if block_num---*/
        ++section_block_num;
      } /*---k_block---*/
    }

    /*---Set index part 3: (block sections) i_block!=j_block!=k_block part---*/

    for (int k_i_offset = 1; k_i_offset < num_block; ++k_i_offset) {
      const int k_block = gm_mod_i(i_block + k_i_offset, num_block);
      for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset){
        const int j_block = gm_mod_i(i_block + j_i_offset, num_block);
        if (j_block == k_block) {
          continue;
        }
        if (gm_proc_r_active(section_block_num, env)) {

          const int section_axis = gm_section_axis_part3(i_block, j_block,
                                                         k_block);
          const int section_num = gm_section_num_part3(i_block, j_block,
                                                       k_block);

          const int J_lo = gm_J_lo(section_num, nvl, 3, env);
          const int J_hi = gm_J_hi(section_num, nvl, 3, env);
          metrics->J_lo_part3_[section_num] = J_lo;
          metrics->J_wi_part3_[section_num] = J_hi - J_lo;
          const bool sax0 = section_axis == 0;
          const bool sax1 = section_axis == 1;
          for (int J = J_lo; J < J_hi; ++J) {
            #pragma omp parallel for collapse(2) schedule(dynamic,1000)
            for (int K = 0; K < nvl; ++K) {
              for (int I = 0; I < nvl; ++I) {

                /* clang-format off */
                const int i = sax0 ?   J :
                              sax1 ?   I :
                           /* sax2 ?*/ K;
                const int j = sax0 ?   K :
                              sax1 ?   J :
                           /* sax2 ?*/ I;
                const int k = sax0 ?   I :
                              sax1 ?   K :
                           /* sax2 ?*/ J;
                /* clang-format on */

                const size_t j_global = j + nvl * j_block;
                const size_t k_global = k + nvl * k_block;
                const size_t i_global = i + nvl * i_block;
                GMAssert(i_global>=0 && metrics->num_vector-i_global>0);
                GMAssert(j_global>=0 && metrics->num_vector-j_global>0);
                GMAssert(k_global>=0 && metrics->num_vector-k_global>0);
                const size_t index_this = index + I + K * (size_t)nvl;
                GMAssert(index_this>=0 && index_this<metrics->num_elts_local);
                metrics->coords_global_from_index[index_this] =
                    i_global +
                    metrics->num_vector *
                        (j_global + metrics->num_vector * (k_global));
              }
            }
            index += nvl*(size_t)nvl;
          }
        } /*---if block_num---*/
        ++section_block_num;
      } /*---j_block---*/
    }   /*---k_block---*/

    GMInsist(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_2 && ! GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsistInterface(env, env->num_stage == 1
                      && "Staged computations not allowed for non-all2all case.");

    GMInsistInterface(env, env->num_phase == 1
                      && "Phased computations not allowed for non-all2all case.");

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
    /*---Need store only strict upper triangular part of matrix---*/
    size_t index = 0;
    for (int j = 0; j < nvl; ++j) {
      const size_t j_global = j + nvl * i_block;
      for (int i = 0; i < j; ++i) {
        const size_t i_global = i + nvl * i_block;
        metrics->coords_global_from_index[index++] =
            i_global + metrics->num_vector * j_global;
      }
    }
    GMInsist(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_3 && ! GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsistInterface(env, env->num_stage == 1
                      && "Staged computations not allowed for non-all2all case.");

    GMInsistInterface(env, env->num_phase == 1
                      && "Phased computations not allowed for non-all2all case.");

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
    /*---Need store only strict interior of tetrahedron---*/
    size_t index = 0;
    for (int j = 0; j < nvl; ++j) {
      const int j_block = i_block;
      const size_t j_global = j + nvl * j_block;
      for (int k = j+1; k < nvl; ++k) {
        const int k_block = i_block;
        const size_t k_global = k + nvl * k_block;
        for (int i = 0; i < j; ++i) {
          const size_t i_global = i + nvl * i_block;
          GMAssert(index < metrics->num_elts_local);
          metrics->coords_global_from_index[index++] =
              i_global +
              metrics->num_vector *
                  (j_global + metrics->num_vector * (k_global));
        }
      }
    }
    GMInsist(index == metrics->num_elts_local);

  /*==================================================*/
  } else {
  /*==================================================*/
    GMInsistInterface(env, 0 == 1 && "Invalid set of options");
    /*---LATER: generalize this to N-way---*/
  }

  size_t num_elts = 0;
  mpi_code = MPI_Allreduce(&metrics->num_elts_local, &num_elts, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);

  if (GMEnv_num_way(env) == GM_NUM_WAY_2 && env->num_stage == 1 &&
      env->num_phase == 1 && GMEnv_all2all(env)) {
    GMInsist(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) / 2);
  }

  if (GMEnv_num_way(env) == GM_NUM_WAY_3 && env->num_stage == 1 &&
      GMEnv_all2all(env)) {
    GMInsist(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) * (size_t)
                               (metrics->num_vector - 2) / 6);
  }

  /*---Allocations---*/

  switch (data_type_id) {
    /*----------*/
    case GM_DATA_TYPE_FLOAT:
      metrics->data_size = metrics->num_elts_local * sizeof(GMFloat);
      metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data_type_num_values = 1;
      break;
    /*----------*/
    case GM_DATA_TYPE_TALLY2X2: {
      metrics->data_size = metrics->num_elts_local * sizeof(GMTally2x2);
      metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data_S_size = metrics->num_elts_local * sizeof(GMFloat2);
      metrics->data_S = gm_malloc(metrics->data_S_size, env);
      if (env->sparse) {
        metrics->data_C_size = metrics->num_elts_local * sizeof(GMFloat2);
        metrics->data_C = gm_malloc(metrics->data_C_size, env);
      }
      metrics->data_type_num_values = 4;
    } break;
    /*----------*/
    case GM_DATA_TYPE_TALLY4X2: {
      metrics->data_size = metrics->num_elts_local * sizeof(GMTally4x2);
      metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data_S_size = metrics->num_elts_local * sizeof(GMFloat3);
      metrics->data_S = gm_malloc(metrics->data_S_size, env);
      if (env->sparse) {
        metrics->data_C_size = metrics->num_elts_local * sizeof(GMFloat3);
        metrics->data_C = gm_malloc(metrics->data_C_size, env);
      }
      metrics->data_type_num_values = 8;
    } break;
    /*----------*/
    default:
      GMInsist(false && "Invalid data type.");
  } /*---switch---*/
}

//=============================================================================
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMInsist(metrics && env);
  GMInsist(metrics->data || ! GMEnv_is_proc_active(env));

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  gm_free(metrics->coords_global_from_index,
          metrics->num_elts_local * sizeof(size_t), env);
  gm_free(metrics->data, metrics->data_size, env);
  if (metrics->data_S) {
    gm_free(metrics->data_S, metrics->data_S_size, env);
  }
  if (metrics->data_C) {
    gm_free(metrics->data_C, metrics->data_C_size, env);
  }
  *metrics = GMMetrics_null();
}

//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: generic---*/

int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                      size_t index,
                                      int coord_num,
                                      GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(coord_num >= 0 && coord_num < GMEnv_num_way(env));

  size_t result64 = 0;

  GMAssert(GMEnv_num_way(env) <= GM_NUM_NUM_WAY + 1
               && "this num_way currently not supported");

  switch (GMEnv_num_way(env) + 4 * coord_num) {
    case 2 + 4 * 0: /* 2-way, coord 0 */
      result64 = GMMetrics_coord0_global_from_index_2(metrics, index, env);
      break;
    case 2 + 4 * 1: /* 2-way, coord 1 */
      result64 = GMMetrics_coord1_global_from_index_2(metrics, index, env);
      break;
    case 3 + 4 * 0: /* 3-way, coord 0 */
      result64 = GMMetrics_coord0_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 1: /* 3-way, coord 1 */
      result64 = GMMetrics_coord1_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 2: /* 3-way, coord 2 */
      result64 = GMMetrics_coord2_global_from_index_3(metrics, index, env);
      break;
    default:
      GMInsistInterface(env, false && "Unimplemented.");
  } /*---case---*/

  const int result = (int)result64;
  GMAssert((size_t)result == result64);

  return result;
}

//-----------------------------------------------------------------------------

void gm_metrics_pad_adjust(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                           GMEnv* env) {
  GMInsist(metrics && metrics_buf && env);

  if (! (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
         GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU)) {
    return;
  }

  const int pad_adjustment = 4 * metrics->dm->num_pad_field_local;

  const GMFloat float_pad_adjustment = GMTally1_encode(pad_adjustment, 0);

  #pragma omp parallel for collapse(2) schedule(dynamic,1000)
  for (size_t j = 0; j < metrics_buf->dim1; ++j) {
    for (size_t i = 0; i < metrics_buf->dim0; ++i) {

#ifdef GM_ASSERTIONS_ON
      const GMTally2x2 old = GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
#endif

      GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j).data[0]
         -= float_pad_adjustment;

#ifdef GM_ASSERTIONS_ON
      const GMTally2x2 new_ = GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
      GMAssert(GMTally2x2_get(old, 0, 0) ==
               GMTally2x2_get(new_, 0, 0) + pad_adjustment);
#endif

    } // for j
  }   // for i
}

//-----------------------------------------------------------------------------
