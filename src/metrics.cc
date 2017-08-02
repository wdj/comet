/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.cc
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

#include "env.hh"
#include "metrics.hh"

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
  memset((void*)&result, 0, sizeof(result));
  return result;
}

/*===========================================================================*/

void GMMetrics_3way_num_elts_local(GMMetrics* metrics, int nvl,
                                   GMEnv* env) {
  GMAssertAlways(metrics);
  GMAssertAlways(env);
  GMAssertAlways(nvl >= 0);
  GMAssertAlways(GMEnv_num_block_vector(env) <= 2 || nvl % 6 == 0);
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_3);

  metrics->num_elts_local = 0;

  //---Fused counter for section_num and block_num, same across all procs.
  int section_block_num = 0;

  //---Compute size part 1: (tetrahedron) i_block==j_block==k_block part.

  const int num_section_steps_1 = GMEnv_num_section_steps(env, 1);
  for (int section_step=0; section_step<num_section_steps_1; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 1, env);
    const int J_hi = gm_J_hi(section_num, nvl, 1, env);
    const GMInt64 trap_size_lo = gm_trap_size(J_lo, nvl);
    const GMInt64 trap_size_hi = gm_trap_size(J_hi, nvl);
    GMAssertAlways(trap_size_hi >= trap_size_lo);
    //---Absorb size_lo into offset for speed in indexing function.
    metrics->index_offset_section_part1_[section_num]
      = (GMInt64)metrics->num_elts_local - trap_size_lo;
    if (gm_proc_r_active(section_block_num, env)) {
      //---Elements in slice of trapezoid.
      const GMInt64 elts_local = trap_size_hi - trap_size_lo;
      GMAssertAlways(elts_local >= 0);
      metrics->num_elts_local += elts_local;
      metrics->section_num_valid_part1_[section_num] = (elts_local != 0);
    }
    ++section_block_num;
  }
  metrics->index_offset_0_ = metrics->num_elts_local;

  //---Compute size part 2: (triang prisms) i_block!=j_block==k_block part.

  const int num_block = GMEnv_num_block_vector(env);
  const int num_section_steps_2 = GMEnv_num_section_steps(env, 2);
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
        GMAssertAlways(elts_local >= 0);
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
        GMAssertAlways(elts_local >= 0);
        metrics->num_elts_local += elts_local;
      }
      ++section_block_num;
    }
  }
}

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      int num_field,
                      size_t num_field_active,
                      int num_vector_local,
                      size_t num_vector_active,
                      GMEnv* env) {
  GMAssertAlways(metrics);
  GMAssertAlways(num_field >= 0);
  GMAssertAlways(num_field_active >= 0);
  GMAssertAlways(num_field_active <= (size_t)num_field);
  GMAssertAlways(num_vector_local >= 0);
  GMAssertAlways(env);

  *metrics = GMMetrics_null();

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  GMInsist(env, GMEnv_all2all(env) || GMEnv_num_proc_repl(env) == 1
          ? "multidim parallelism only available for all2all case" : 0);

  GMInsist(env,
           num_field % GMEnv_num_proc_field(env) == 0
               ? "num_proc_field must exactly divide the total number of fields"
               : 0);

  /*---These cases less important, not yet tested---*/

  GMInsist(env, GMEnv_all2all(env) || (size_t)num_field == num_field_active
                ? "This case currently not supported." : 0);

  GMInsist(env, GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU ||
                (size_t)num_field == num_field_active
                ? "This case currently not supported." : 0);

  metrics->data_type_id = data_type_id;
  metrics->num_field = num_field;
  metrics->num_field_active = num_field_active;
  metrics->num_field_local = num_field / GMEnv_num_proc_field(env);
  metrics->num_vector_local = num_vector_local;
  metrics->num_vector_active = num_vector_active;
  metrics->nvl6 = num_vector_local / 6;
  metrics->index_offset_0_ = 0;
  metrics->index_offset_01_ = 0;
  metrics->recip_m = ((GMFloat)1) / num_field_active;
  metrics->block_min = 0;
  for (int i=0; i<6; ++i) {
    metrics->index_offset_section_part1_[i] = 0;
    metrics->index_offset_section_part2_[i] = 0;
    metrics->section_num_valid_part1_[i] = GM_BOOL_FALSE;
    metrics->section_num_valid_part2_[i] = GM_BOOL_FALSE;
  }

  /*---Compute global values---*/

  const int num_block = GMEnv_num_block_vector(env);

  const size_t num_vector_bound = metrics->num_vector_local *
                          (size_t)num_block * (size_t)GMEnv_num_proc_repl(env);
  GMAssertAlways(num_vector_bound == (size_t)(int)num_vector_bound
    ? "Vector count too large to store in 32-bit int; please modify code." : 0);

  int mpi_code = 0;
  mpi_code = MPI_Allreduce(&(metrics->num_vector_local), &(metrics->num_vector),
                           1, MPI_INT, MPI_SUM, GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
  GMAssertAlways((size_t)(metrics->num_vector) == num_vector_bound);
  metrics->num_vector /= GMEnv_num_proc_repl(env);
  GMAssertAlways(metrics->num_vector_active <= (size_t)metrics->num_vector);

  /*---Assume the following to simplify calculations---*/

  GMInsist(
      env,
      num_vector_local >= GMEnv_num_way(env)
          ? "Currently require number of vecs on a proc to be at least num-way"
          : 0);

  const int i_block = GMEnv_proc_num_vector_i(env);

  const size_t nchoosek = gm_nchoosek(num_vector_local, GMEnv_num_way(env));
  const int nvl = num_vector_local;
  const size_t nvlsq = nvl * (size_t)nvl;

  /*---Compute number of elements etc.---*/

  GMInsist(env, env->stage_num >= 0 && env->stage_num < env->num_stage
                ? "Invalid stage number specified."
                : 0);

  GMInsist(env, env->phase_num >= 0 && env->phase_num < env->num_phase
                ? "Invalid phase number specified."
                : 0);

  metrics->num_elts_local_computed = 0;

  /*==================================================*/
  if (GMEnv_num_way(env) == GM_NUM_WAY_2 && GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsist(env, env->num_stage == 1
                      ? "Staged computations not allowed for 2-way case."
                      : 0);

    GMInsist(env, env->num_phase <= 1 + num_block / 2
                      ? "num_phase must be at most 1 + nproc_vector/2."
                      : 0);

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
    const _Bool have_main_diag = proc_num_r == 0 &&
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
    GMAssertAlways(metrics->coords_global_from_index != NULL);
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
      if (diag == 0 || !gm_proc_r_active(diag-beg, env)) {
        continue;
      }
      const int j_block_unwrapped = i_block + diag;
#pragma omp parallel for collapse(2)
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
    GMAssertAlways(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_3 && GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsist(env, env->num_phase == 1
                      ? "Phaseed computations not currently implemented "
                        "for 3-way case."
                      : 0);

    /*---Make the following assumption to greatly simplify calculations---*/
    GMInsist(env, num_block <= 2 || metrics->num_vector_local % 6 == 0
                      ? "3way all2all case requires num vectors per proc "
                        "divisible by 6."
                      : 0);

    /*===PART A: CALCULATE INDEX SIZE===*/

    GMMetrics_3way_num_elts_local(metrics, nvl, env);

    /*---Fused counter for section_num and block_num, same across all procs---*/
    int section_block_num = 0;

    const int num_section_steps_12 = GMEnv_num_section_steps(env, 1);

    /*===PART B: ALLOCATE INDEX===*/

    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
    GMAssertAlways(metrics->coords_global_from_index != NULL);

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
#pragma omp parallel for collapse(2)
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
#pragma omp parallel for collapse(2)
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
          const _Bool sax0 = section_axis == 0;
          const _Bool sax1 = section_axis == 1;
          for (int J = J_lo; J < J_hi; ++J) {
#pragma omp parallel for collapse(2)
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

    GMAssertAlways(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_2 && !GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsist(env, env->num_stage == 1
                      ? "Staged computations not allowed for non-all2all case."
                      : 0);

    GMInsist(env, env->num_phase == 1
                      ? "Phased computations not allowed for non-all2all case."
                      : 0);

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
    GMAssertAlways(metrics->coords_global_from_index != NULL);
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
    GMAssertAlways(index == metrics->num_elts_local);

  /*==================================================*/
  } else if (GMEnv_num_way(env) == GM_NUM_WAY_3 && !GMEnv_all2all(env)) {
  /*==================================================*/

    GMInsist(env, env->num_stage == 1
                      ? "Staged computations not allowed for non-all2all case."
                      : 0);

    GMInsist(env, env->num_phase == 1
                      ? "Phased computations not allowed for non-all2all case."
                      : 0);

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        (size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
    GMAssertAlways(metrics->coords_global_from_index != NULL);
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
    GMAssertAlways(index == metrics->num_elts_local);

  /*==================================================*/
  } else {
  /*==================================================*/
    GMInsist(env, 0 == 1 ? "Invalid set of options" : 0);
    /*---LATER: generalize this to N-way---*/
  }

  size_t num_elts = 0;
  mpi_code = MPI_Allreduce(&metrics->num_elts_local, &num_elts, 1,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_vector(env));

  if (GMEnv_num_way(env) == GM_NUM_WAY_2 && env->num_stage == 1 &&
      env->num_phase == 1 && GMEnv_all2all(env)) {
    GMAssertAlways(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) / 2);
  }

  if (GMEnv_num_way(env) == GM_NUM_WAY_3 && env->num_stage == 1 &&
      GMEnv_all2all(env)) {
    GMAssertAlways(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) * (size_t)
                               (metrics->num_vector - 2) / 6);
  }

  /*---Allocations---*/

  switch (data_type_id) {
    /*----------*/
    case GM_DATA_TYPE_BITS1: {
      /*---(design not complete)---*/
      const size_t num_floats_needed =
          gm_ceil_i8(metrics->num_elts_local, 8 * sizeof(GMFloat));
      metrics->data_size = num_floats_needed * sizeof(GMFloat);
      metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data_type_num_values = 1;
    } break;
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
      GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
  } /*---switch---*/
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  GMAssertAlways(metrics);
  GMAssertAlways(env);
  GMAssertAlways(metrics->data != NULL || !GMEnv_is_proc_active(env));

  if (!GMEnv_is_proc_active(env)) {
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

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: generic---*/

int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                      size_t index,
                                      int coord_num,
                                      GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index+1 >= 1);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(coord_num >= 0);
  GMAssert(coord_num < GMEnv_num_way(env));

  size_t result64 = 0;

  GMAssert(GMEnv_num_way(env) <= GM_NUM_NUM_WAY + 1
               ? "this num_way currently not supported"
               : 0);

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
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

  const int result = (int)result64;
  GMAssert((size_t)result == result64);

  return result;
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

void GMMetrics_checksum(GMMetrics* metrics, GMChecksum* cs, GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(metrics->data != NULL || !GMEnv_is_proc_active(env));

  /*---Initializations---*/

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  enum { num_way_max = GM_NUM_NUM_WAY + 1 };
  GMAssertAlways(GMEnv_num_way(env) <= num_way_max ?
                 "This num_way not supported." : 0);

  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  /*---Check for NaNs if appropriate---*/

  switch (metrics->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat_check((GMFloat*)(metrics->data), metrics->num_elts_local);
    } break;
  }

  /*---Calculate the global largest value---*/

  double value_max_this = cs->value_max;
#pragma omp parallel for reduction(max:value_max_this)
  for (UI64 index = 0; index < metrics->num_elts_local; ++index) {
    _Bool is_active = GM_BOOL_TRUE;
    for (int i = 0; i < GMEnv_num_way(env); ++i) {
      const UI64 coord = GMMetrics_coord_global_from_index(metrics, index,
                                                           i, env);
      is_active = is_active && coord < metrics->num_vector_active;
    }
    /*---Loop over data values at this index---*/
    for (int i_value = 0; i_value < metrics->data_type_num_values; ++i_value) {
      /*---Pick up value of this metrics elt---*/
      double value = 0;
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
          const int i0 = i_value / 2;
          const int i1 = i_value % 2;
          value = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
        } break;
        /*--------------------*/
        case GM_DATA_TYPE_TALLY4X2: {
          const int i0 = i_value / 4;
          const int i1 = (i_value / 2) % 2;
          const int i2 = i_value % 2;
          value =
              GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2, env);
        } break;
        /*--------------------*/
        default:
          GMAssertAlways(GM_BOOL_FALSE ? "Invalid data type." : 0);
      } /*---switch---*/
      if (is_active) {
        value_max_this = value > value_max_this ? value : value_max_this;
      }
    }
  } /*---for index---*/

  int mpi_code = 0;
  mpi_code *= 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Allreduce(&value_max_this, &cs->value_max, 1,
                           MPI_DOUBLE, MPI_MAX, GMEnv_mpi_comm_vector(env));

  /*---The largest we expect any value to be if using "special" inputs---*/
  const int log2_value_max_allowed = 4;
  const double value_max_allowed = 1 << log2_value_max_allowed;

  cs->is_overflowed = cs->is_overflowed && cs->value_max > value_max_allowed;

  const double scaling = value_max_allowed;

  //const double scaling = cs->is_overflowed ? cs->value_max : value_max_allowed;

  /*---Calculate checksum---*/

  //GMMultiprecInt sums_l = {0};

  //UI64 sums_l[16];
  //for (int i = 0; i < 16; ++i) {
  //  sums_l[i] = 0;
  //}
  //double sum_d = 0;

  const int w = 30;
  GMAssertAlways(64 - 2 * w >= 4);
  const UI64 one64 = 1;

  const UI64 lomask = (one64 << w) - 1;
  const UI64 lohimask = (one64 << (2 * w)) - 1;

  GMMultiprecInt sum_this = {0};
  double sum_d_this = 0;

#pragma omp parallel
{
  GMMultiprecInt sum_this_private = {0};
  double sum_d_this_private = 0;
#pragma omp for collapse(2)
  for (UI64 index = 0; index < metrics->num_elts_local; ++index) {
    /*---Loop over data values at this index---*/
    for (int i_value = 0; i_value < metrics->data_type_num_values; ++i_value) {

      /*---Obtain global coords of metrics elt---*/
      UI64 coords[num_way_max];
      int ind_coords[num_way_max];
      for (int i = 0; i < num_way_max; ++i) {
        coords[i] = 0;
        ind_coords[i] = i;
      }
      _Bool is_active = GM_BOOL_TRUE;
      for (int i = 0; i < GMEnv_num_way(env); ++i) {
        const UI64 coord = GMMetrics_coord_global_from_index(metrics, index,
                                                             i, env);
        is_active = is_active && coord < metrics->num_vector_active;
        coords[i] = coord;
      }
      /*---Reflect coords by symmetry to get uniform result -
           sort into descending order---*/
      gm_makegreater(&coords[1], &coords[2], &ind_coords[1], &ind_coords[2]);
      gm_makegreater(&coords[0], &coords[1], &ind_coords[0], &ind_coords[1]);
      gm_makegreater(&coords[1], &coords[2], &ind_coords[1], &ind_coords[2]);

      /*---Pick up value of this metrics elt---*/
      double value = 0;
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
      UI64 ivalue = (value / scaling) * (one64 << (2 * w));
      /*---Construct global id for metrics data vbalue---*/
      UI64 uid = coords[0];
      for (int i = 1; i < GMEnv_num_way(env); ++i) {
        uid = uid * metrics->num_vector_active + coords[i];
      }
      uid = uid * metrics->data_type_num_values + i_value;
      /*---Randomize---*/
      const UI64 rand1 = gm_randomize(uid + 956158765);
      const UI64 rand2 = gm_randomize(uid + 842467637);
      UI64 rand_value = rand1 + gm_randomize_max() * rand2;
      rand_value &= lohimask;
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
      /*---(move the carry bits)---*/
      chi += clo >> (2 * w);
      clo &= lohimask;
      const double value_d =
          ivalue * (double)rand_value / ((double)(one64 << (2 * w)));
      if (is_active) {
        sum_d_this_private += value_d; /*---Reduction---*/
        /*---Split the product into one-char chunks, accumulate to sums---*/
        for (int i = 0; i < 8; ++i) {
          const UI64 value0 = (clo << (64 - 8 - 8 * i)) >> (64 - 8);
          const UI64 value1 = (chi << (64 - 8 - 8 * i)) >> (64 - 8);
          sum_this_private.data[0 + i] += value0; /*---Reduction---*/
          sum_this_private.data[8 + i] += value1; /*---Reduction---*/
        }
      }
    } /*---for i_value---*/
  }   /*---for index---*/

#pragma omp critical
  {
      sum_d_this += sum_d_this_private; /*---Reduction---*/
      for (int i = 0; i < 8; ++i) {
        sum_this.data[0 + i] += sum_this_private.data[0 + i];/*---Reduction---*/
        sum_this.data[8 + i] += sum_this_private.data[8 + i];/*---Reduction---*/
      }
  }
} /*---omp parallel---*/

  /*---Global sum---*/
  GMMultiprecInt sum = {0};
  mpi_code = MPI_Allreduce(sum_this.data, sum.data, GM_MULTIPREC_INT_SIZE,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
  for (int i = 0; i < GM_MULTIPREC_INT_SIZE; ++i) {
    cs->sum.data[i] += sum.data[i];
  }

  /*---Combine results---*/

  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    cs->data[i] = 0;
    for (int j = 0; j < 8; ++j) {
      cs->data[i] +=
          gm_lshift(cs->sum.data[0 + j], 8 * j - 2 * w * i) & lohimask;
      cs->data[i] +=
          gm_lshift(cs->sum.data[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
    }
  }
  /*---(move the carry bits---*/
  cs->data[1] += cs->data[0] >> (2 * w);
  cs->data[0] &= lohimask;
  cs->data[2] += cs->data[1] >> (2 * w);
  cs->data[1] &= lohimask;

  /*---Check against floating point result---*/

  const double tmp = sum_d_this;
  mpi_code = MPI_Allreduce(&tmp, &sum_d_this, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_vector(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
  cs->sum_d += sum_d_this;

  double result_d = cs->data[0] / ((double)(one64 << (2 * w))) +
                    cs->data[1] +
                    cs->data[2] * ((double)(one64 << (2 * w)));
  result_d *= 1; /*---Avoid unused variable warning---*/
  GMAssertAlways(fabs(cs->sum_d - result_d) <= cs->sum_d * 1.e-10);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
