//-----------------------------------------------------------------------------
/*!
 * \file   metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the calculated metrics output from the methods.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdlib"
#include "cstdint"
#include "string.h"
#include "math.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Helper class for memory---*/

GMMetricsMem::GMMetricsMem(GMEnv* env) :
  env_(env),
  data_(0),
  data_size_(0),
  data_S_(0),
  data_S_size_(0),
  data_C_(0),
  data_C_size_(0),
  coords_global_from_index_(0),
  coords_global_from_index_size_(0) {
}

//-----------------------------------------------------------------------------

GMMetricsMem::~GMMetricsMem() {

  if (! env_->is_proc_active()) {
    return;
  }
  if (data_) {
    gm_free(data_, data_size_, env_);
  }
  if (data_S_) {
    gm_free(data_S_, data_S_size_, env_);
  }
  if (data_C_) {
    gm_free(data_C_, data_C_size_, env_);
  }
  if (coords_global_from_index_) {
    gm_free(coords_global_from_index_, coords_global_from_index_size_, env_);
  }
}

//-----------------------------------------------------------------------------

void* GMMetricsMem::malloc_data(size_t data_size) {

  if (! env_->is_proc_active()) {
    return NULL;
  }

  if (!data_ || data_size > data_size_) {
    if (data_) {
      gm_free(data_, data_size_, env_);
    }
    data_ = (size_t*)gm_malloc(data_size, env_);
    data_size_ = data_size;
  }

  return data_;
}

//-----------------------------------------------------------------------------

void* GMMetricsMem::malloc_data_S(size_t data_S_size) {

  if (! env_->is_proc_active()) {
    return NULL;
  }

  if (!data_S_ || data_S_size > data_S_size_) {
    if (data_S_) {
      gm_free(data_S_, data_S_size_, env_);
    }
    data_S_ = (size_t*)gm_malloc(data_S_size, env_);
    data_S_size_ = data_S_size;
  }

  return data_S_;
}

//-----------------------------------------------------------------------------

void* GMMetricsMem::malloc_data_C(size_t data_C_size) {

  if (! env_->is_proc_active()) {
    return NULL;
  }

  if (!data_C_ || data_C_size > data_C_size_) {
    if (data_C_) {
      gm_free(data_C_, data_C_size_, env_);
    }
    data_C_ = (size_t*)gm_malloc(data_C_size, env_);
    data_C_size_ = data_C_size;
  }

  return data_C_;
}

//-----------------------------------------------------------------------------

void* GMMetricsMem::malloc_coords_global_from_index(
  size_t coords_global_from_index_size) {

  if (! env_->is_proc_active()) {
    return NULL;
  }

  if (!coords_global_from_index_ ||
      coords_global_from_index_size > coords_global_from_index_size_) {
    if (coords_global_from_index_) {
      gm_free(coords_global_from_index_, coords_global_from_index_size_, env_);
    }
    coords_global_from_index_
       = (size_t*)gm_malloc(coords_global_from_index_size, env_);
    coords_global_from_index_size_ = coords_global_from_index_size;
  }

  return coords_global_from_index_;
}

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
  COMET_INSIST(metrics && env);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(env->num_block_vector() <= 2 || nvl % 6 == 0);
  COMET_INSIST(env->num_way() == NUM_WAY::_3);

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
    const int64_t trap_size_lo = gm_trap_size(J_lo, nvl);
    const int64_t trap_size_hi = gm_trap_size(J_hi, nvl);
    COMET_INSIST(trap_size_hi >= trap_size_lo && "Error in sizes calculation.");
    //---Absorb size_lo into offset for speed in indexing function.
    metrics->index_offset_section_part1_[section_num]
      = (int64_t)metrics->num_elts_local - trap_size_lo;
    if (gm_is_section_block_in_phase(env, section_block_num)) {
      if (gm_proc_r_active(section_block_num, env)) {
        //---Elements in slice of trapezoid.
        const int64_t elts_local = trap_size_hi - trap_size_lo;
        COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
        metrics->num_elts_local += elts_local;
        metrics->section_num_valid_part1_[section_num] = (elts_local != 0);
      } // if
    } // if
    ++section_block_num;
  } // section_step
  metrics->index_offset_0_ = metrics->num_elts_local;

  //---Compute size part 2: (triang prisms) i_block!=j_block==k_block part.

  const int num_block = env->num_block_vector();
  const int num_section_steps_2 = gm_num_section_steps(env, 2);
  for (int section_step=0; section_step<num_section_steps_2; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 2, env);
    const int J_hi = gm_J_hi(section_num, nvl, 2, env);
    const int64_t triang_size_lo = gm_triang_size(J_lo, nvl);
    const int64_t triang_size_hi = gm_triang_size(J_hi, nvl);
    //---Absorb size_lo into offset for speed in indexing function.
    metrics->index_offset_section_part2_[section_num]
      = (int64_t)metrics->num_elts_local - (int64_t)nvl*(int64_t)triang_size_lo;
    metrics->section_size_part2[section_num] = triang_size_hi -
                                               triang_size_lo;
    int block_num_part2 = 0;
    bool is_phase_block_start_set = false;
    //---Loop over blocks for part2.
    for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
      if (gm_is_section_block_in_phase(env, section_block_num)) {
        if (!is_phase_block_start_set) {
          metrics->phase_block_start_2_[section_num] = block_num_part2;
          is_phase_block_start_set = true;
        }
        if (gm_proc_r_active(section_block_num, env)) {
          //---Elements in slice of triang prism.
          const int64_t elts_local = (int64_t)nvl *
                                     (triang_size_hi - triang_size_lo);
          COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
          metrics->num_elts_local += elts_local;
          metrics->section_num_valid_part2_[section_num] = (elts_local != 0);
        } // if
      } // if
      ++section_block_num;
      ++block_num_part2;
    }
  } // section_step
  metrics->index_offset_01_ = metrics->num_elts_local;

  //---Compute size part 3: (block sections) i_block!=j_block!=k_block part.

  //---Loop over block for part3.
  const int num_section_steps_3 = gm_num_section_steps(env, 3); // = 1
  for (int section_step=0; section_step<num_section_steps_3; ++section_step) {
    const int i_block = env->proc_num_vector();
    int block_num_part3 = 0;
    bool is_phase_block_start_set = false;
    for (int k_i_offset=1; k_i_offset<num_block; ++k_i_offset) {
      const int k_block = utils::mod_i(i_block + k_i_offset, num_block);
      for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
        const int j_block = utils::mod_i(i_block + j_i_offset, num_block);
        if (j_block == k_block) {
          continue;
        }
        //---Get slice bounds.
        const int section_num = gm_section_num_part3(i_block, j_block, k_block);
        const int J_lo = gm_J_lo(section_num, nvl, 3, env);
        const int J_hi = gm_J_hi(section_num, nvl, 3, env);
        if (gm_is_section_block_in_phase(env, section_block_num)) {
          if (!is_phase_block_start_set) {
            metrics->phase_block_start_3_ = block_num_part3;
            is_phase_block_start_set = true;
          }
          if (gm_proc_r_active(section_block_num, env)) {
            //---Elements in slice of block/cube.
            const int64_t elts_local = (int64_t)nvl * (int64_t)nvl *
                                       (int64_t)(J_hi - J_lo);
            COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
            metrics->num_elts_local += elts_local;
          } // if
        } // if
        ++section_block_num;
        ++block_num_part3;
      }
    }
  } // section_step
}

//=============================================================================
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      GMDecompMgr* dm,
                      GMMetricsMem* metrics_mem,
                      GMEnv* env) {
  COMET_INSIST(metrics && dm && env);

  *metrics = GMMetrics_null();

  if (! env->is_proc_active()) {
    return;
  }

  COMET_INSIST_INTERFACE(env, (env->metric_type() != MetricType::DUO ||
                          env->sparse()) && "DUO method requires sparse input.");

  COMET_INSIST_INTERFACE(env, (env->all2all() || env->num_proc_repl() == 1)
          && "Multidim parallelism only available for all2all case");

  /*---The following less important cases are not yet tested---*/

  COMET_INSIST_INTERFACE(env, (env->all2all() ||
                    dm->num_field == dm->num_field_active)
                    && "This case currently not supported.");

  COMET_INSIST_INTERFACE(env, (env->compute_method() == ComputeMethod::GPU ||
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
    metrics->phase_block_start_2_[i] = 0;
  }
  metrics->phase_block_start_3_ = 0;

  metrics->num_vector = dm->num_vector;

  const int num_block = env->num_block_vector();

  COMET_INSIST_INTERFACE(
      env,
      metrics->num_vector_local >= env->num_way()
        && "Currently require number of vecs on a proc to be at least num-way");

  const int i_block = env->proc_num_vector();

  const size_t nchoosek = utils::nchoosek(metrics->num_vector_local,
                                          env->num_way());
  const int nvl = metrics->num_vector_local;
  const size_t nvlsq = nvl * (size_t)nvl;

  /*---Compute number of elements etc.---*/

  COMET_INSIST_INTERFACE(env, env->stage_num() >= 0
    && env->stage_num() < env->num_stage()
    && "Invalid stage number specified.");

  COMET_INSIST_INTERFACE(env, env->phase_num() >= 0 && env->phase_num() < env->num_phase()
                    && "Invalid phase number specified.");

  GMMetrics_ccc_check_size_nofp_2(metrics, env);
  GMMetrics_ccc_check_size_nofp_3(metrics, env);

  metrics->num_elts_local_computed = 0;

  /*==================================================*/
  if (env->num_way() == NUM_WAY::_2 && env->all2all()) {
  /*==================================================*/

    COMET_INSIST_INTERFACE(env, env->num_stage() == 1
                      && "Staged computations not allowed for 2-way case.");

    COMET_INSIST_INTERFACE(env, env->num_phase() <= 1 + num_block / 2
                      && "num_phase must be at most 1 + num_proc_vector/2.");

    /*---Store the following in this block-row:
        1) strict upper triangular part of main diagonal block
        2) half of the off-diagonal blocks, as a "wrapped rectangle"
      For num_proc_repl > 1, map these blocks, starting at the
      main diagonal block, to procs in round-robin fashion.
      For num_phase > 1, do all this only for a piece of the block row.
    ---*/

    /*===PART A: CALCULATE INDEX SIZE===*/
    const int proc_num_r = env->proc_num_repl();
    const int num_proc_r = env->num_proc_repl();
    metrics->num_elts_local = 0;

    /*---PART A.1: (triangle) i_block==j_block part---*/
    const bool have_main_diag = proc_num_r == 0 &&
                                gm_bdiag_computed_min(env) == 0;
    metrics->num_elts_local += have_main_diag ? nchoosek : 0;
    metrics->index_offset_0_ = have_main_diag ? nchoosek - nvlsq : 0;
    metrics->block_min = (i_block + gm_bdiag_computed_min(env)) % num_block;


    /*---PART A.2: (wrapped rect) i_block!=j_block part---*/
    const int num_computed_blocks_this_row = gm_blocks_computed_this_row(env);
    const int num_computed_blocks_this_proc = rr_pack_(proc_num_r, num_proc_r,
                                               num_computed_blocks_this_row);
    const int num_computed_offdiag_blocks_this_proc =
      num_computed_blocks_this_proc - (have_main_diag ? 1 : 0);
    metrics->num_elts_local += num_computed_offdiag_blocks_this_proc * nvlsq;

    /*===PART B: ALLOCATE INDEX===*/
    metrics->coords_global_from_index =
        //(size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
        (size_t*)metrics_mem->malloc_coords_global_from_index(metrics->num_elts_local * sizeof(size_t));

    /*===PART C: SET INDEX===*/

    /*---PART C.1: (triangle) i_block==j_block part---*/
    size_t index = 0;
    if (have_main_diag) {
      COMET_ASSERT(env->proc_num_repl() == 0);
      COMET_ASSERT(gm_proc_r_active(0, env));
      for (int j = 0; j < nvl; ++j) {
        const size_t j_global = j + nvl * i_block;
        for (int i = 0; i < j; ++i) {
          const size_t i_global = i + nvl * i_block;
          COMET_ASSERT(GMMetrics_helper2way_maindiag_block_(metrics, i, j, i_block,
                                                        env) == index);
          metrics->coords_global_from_index[index++] =
              i_global + metrics->num_vector * j_global;
        }
      }
    }

    /*---PART C.2: (wrapped rectangle) i_block!=j_block part---*/

    const int beg = gm_bdiag_computed_min(env);
    const int end = beg + num_computed_blocks_this_row;
    for (int diag=beg; diag<end; ++diag) {
      const int diag_offset = diag - beg;
      if (diag == 0 || ! gm_proc_r_active(diag_offset, env)) {
        continue;
      }
      const int j_block_unwrapped = i_block + diag;
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
        const size_t j_global_unwrapped = j + j_block_unwrapped * (size_t)nvl;
        const size_t j_global = j_global_unwrapped % metrics->num_vector;
          const size_t i_global = i + nvl * i_block;
          const size_t index_this = index + i + j * (size_t)nvl;
          COMET_ASSERT(index_this>=0 && index_this<metrics->num_elts_local);
          metrics->coords_global_from_index[index_this] =
              i_global + metrics->num_vector * j_global;
        }
      }
      index += nvlsq;
    } /*---for diag---*/

    /*---Final check---*/
    COMET_INSIST(index == metrics->num_elts_local && "Error in sizes calculation.");

  /*==================================================*/
  } else if (env->num_way() == NUM_WAY::_3 && env->all2all()) {
  /*==================================================*/

    COMET_INSIST_INTERFACE(env, env->num_phase() <= gm_num_section_blocks(env)
     && "num_phase must be at most (num_proc_vector+1)*(num_proc_vector+2)/2.");

    //---Make the following assumption to greatly simplify calculations.
    COMET_INSIST_INTERFACE(env, (num_block<=2 || metrics->num_vector_local % 6 == 0)
                      && "3-way all2all case requires num vectors per proc "
                        "divisible by 6.");

    /*===PART A: CALCULATE INDEX SIZE===*/

    GMMetrics_3way_num_elts_local(metrics, nvl, env);

    /*---Fused counter for section_num and block_num, same across all procs---*/
    int section_block_num = 0;

    /*===PART B: ALLOCATE INDEX===*/

    metrics->coords_global_from_index =
        //(size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
        (size_t*)metrics_mem->malloc_coords_global_from_index(metrics->num_elts_local * sizeof(size_t));

    /*===PART C: SET INDEX===*/

    section_block_num = 0;
    size_t index = 0;

    /*---Set index part 1: (tetrahedron) i_block==j_block==k_block part---*/

    const int num_section_steps_1 = gm_num_section_steps(env, 1);
    for (int section_step=0; section_step<num_section_steps_1; ++section_step){
      if (gm_is_section_block_in_phase(env, section_block_num)) {
        if (gm_proc_r_active(section_block_num, env)) {
          const int section_num = section_step;
          const int J_lo = gm_J_lo(section_num, nvl, 1, env);
          const int J_hi = gm_J_hi(section_num, nvl, 1, env);
          const int j_min = J_lo;
          const int j_max = J_hi;
          for (int j = j_min; j < j_max; ++j) {
            const int j_block = i_block;
            const size_t j_global = j + nvl * j_block;
            // don't use collapse because of overflow for large sizes
            //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
            #pragma omp parallel for schedule(dynamic,1000)
            for (int k = j+1; k < nvl; ++k) {
              for (int i = 0; i < j; ++i) {
              const int k_block = i_block;
              const size_t k_global = k + nvl * k_block;
                const size_t i_global = i + nvl * i_block;
                const size_t index_this = index + i + j*(size_t)(k-(j+1));
                COMET_ASSERT(index_this>=0 && index_this<metrics->num_elts_local);
                metrics->coords_global_from_index[index_this] =
                    i_global +
                  metrics->num_vector *
                        (j_global + metrics->num_vector * (k_global));
              }
            }
            index += j * (size_t)(nvl - (j+1));
          }
        } // if
      } // if
      ++section_block_num;
    } // section_step

    /*---Set index part 2: (triang prisms) i_block!=j_block==k_block part---*/

    const int num_section_steps_2 = gm_num_section_steps(env, 2);
    for (int section_step=0; section_step<num_section_steps_2; ++section_step){
      for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
        const int j_block = utils::mod_i(i_block + j_i_offset, num_block);
        if (gm_is_section_block_in_phase(env, section_block_num)) {
          if (gm_proc_r_active(section_block_num, env)) {
            const int section_num = section_step;
            const int J_lo = gm_J_lo(section_num, nvl, 2, env);
            const int J_hi = gm_J_hi(section_num, nvl, 2, env);
            const int j_min = J_lo;
            const int j_max = J_hi;
            for (int j = j_min; j < j_max; ++j) {
              const size_t j_global = j + nvl * j_block;
              // don't use collapse because of overflow for large sizes
              //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
              #pragma omp parallel for schedule(dynamic,1000)
              for (int k = j+1; k < nvl; ++k) {
                for (int i = 0; i < nvl; ++i) {
                  const int k_block = j_block;
                  const size_t k_global = k + nvl * k_block;
                  const size_t i_global = i + nvl * i_block;
                  const size_t index_this = index + i + nvl*(size_t)(k-(j+1));
                  COMET_ASSERT(index_this>=0 && index_this<metrics->num_elts_local);
                  metrics->coords_global_from_index[index_this] =
                      i_global +
                    metrics->num_vector *
                          (j_global + metrics->num_vector * (k_global));
                } // for i
              } // for k
              index += nvl * (size_t)(nvl - (j+1));
            } // for j
          } // if
        } // if
        ++section_block_num;
      } // for j_i_offset
    } // section_step

    /*---Set index part 3: (block sections) i_block!=j_block!=k_block part---*/

    const int num_section_steps_3 = gm_num_section_steps(env, 3); // = 1
    for (int section_step=0; section_step<num_section_steps_3; ++section_step) {
      for (int k_i_offset = 1; k_i_offset < num_block; ++k_i_offset) {
        const int k_block = utils::mod_i(i_block + k_i_offset, num_block);
        for (int j_i_offset = 1; j_i_offset < num_block; ++j_i_offset){
          const int j_block = utils::mod_i(i_block + j_i_offset, num_block);
          if (j_block == k_block) {
            continue;
          }
          if (gm_is_section_block_in_phase(env, section_block_num)) {
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
                // don't use collapse because of overflow for large sizes
                //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
                #pragma omp parallel for schedule(dynamic,1000)
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
                    COMET_ASSERT(i_global>=0 && metrics->num_vector-i_global>0);
                    COMET_ASSERT(j_global>=0 && metrics->num_vector-j_global>0);
                    COMET_ASSERT(k_global>=0 && metrics->num_vector-k_global>0);
                    const size_t index_this = index + I + K * (size_t)nvl;
                    COMET_ASSERT(index_this>=0 &&
                             index_this<metrics->num_elts_local);
                    metrics->coords_global_from_index[index_this] =
                        i_global +
                        metrics->num_vector *
                            (j_global + metrics->num_vector * (k_global));
                  } // for I
                } // for K
                index += nvl*(size_t)nvl;
              } // for J
            } // if
          } // if
          ++section_block_num;
        } // j_i_offset
      } // k_i_offset
    } // section_step

    COMET_INSIST(index == metrics->num_elts_local && "Error in sizes calculation.");

  /*==================================================*/
  } else if (env->num_way() == NUM_WAY::_2 && ! env->all2all()) {
  /*==================================================*/

    COMET_INSIST_INTERFACE(env, env->num_stage() == 1 &&
                      "Staged computations not allowed for non-all2all case.");

    COMET_INSIST_INTERFACE(env, env->num_phase() == 1 &&
                      "Phased computations not allowed for non-all2all case.");

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        //(size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
        (size_t*)metrics_mem->malloc_coords_global_from_index(
                                    metrics->num_elts_local * sizeof(size_t));
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
    COMET_INSIST(index == metrics->num_elts_local && "Error in sizes calculation.");

  /*==================================================*/
  } else if (env->num_way() == NUM_WAY::_3 && ! env->all2all()) {
  /*==================================================*/

    COMET_INSIST_INTERFACE(env, env->num_stage() == 1
                      && "Staged computations not allowed for non-all2all case.");

    COMET_INSIST_INTERFACE(env, env->num_phase() == 1
                      && "Phased computations not allowed for non-all2all case.");

    metrics->num_elts_local = nchoosek;
    metrics->coords_global_from_index =
        //(size_t*)gm_malloc(metrics->num_elts_local * sizeof(size_t), env);
        (size_t*)metrics_mem->malloc_coords_global_from_index(metrics->num_elts_local * sizeof(size_t));
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
          COMET_ASSERT(index < metrics->num_elts_local);
          metrics->coords_global_from_index[index++] =
              i_global +
              metrics->num_vector *
                  (j_global + metrics->num_vector * (k_global));
        }
      }
    }
    COMET_INSIST(index == metrics->num_elts_local && "Error in sizes calculation.");

  /*==================================================*/
  } else {
  /*==================================================*/
    COMET_INSIST_INTERFACE(env, 0 == 1 && "Invalid set of options");
    /*---LATER: generalize this to N-way---*/
  }

  size_t num_elts = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metrics->num_elts_local, &num_elts, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));

  if (env->num_way() == NUM_WAY::_2 &&
      env->num_phase() == 1 && env->all2all()) {
    COMET_INSIST(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) / 2);
  }

  if (env->num_way() == NUM_WAY::_3 && env->num_stage() == 1 &&
      env->num_phase() == 1 && env->all2all()) {
    COMET_INSIST(num_elts == (metrics->num_vector) * (size_t)
                               (metrics->num_vector - 1) * (size_t)
                               (metrics->num_vector - 2) / 6);
  }

  /*---Allocations---*/

  switch (data_type_id) {
    /*----------*/
    case GM_DATA_TYPE_FLOAT:
      metrics->data_size = metrics->num_elts_local * sizeof(GMFloat);
      //metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data = metrics_mem->malloc_data(metrics->data_size);
      metrics->data_type_num_values = 1;
      break;
    /*----------*/
    case GM_DATA_TYPE_TALLY2X2: {
      metrics->data_size = metrics->num_elts_local * sizeof(GMTally2x2);
      //metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data = metrics_mem->malloc_data(metrics->data_size);
      metrics->data_S_size = metrics->num_elts_local * sizeof(GMFloat2);
      //metrics->data_S = gm_malloc(metrics->data_S_size, env);
      metrics->data_S = metrics_mem->malloc_data_S(metrics->data_S_size);
      if (env->sparse()) {
        metrics->data_C_size = metrics->num_elts_local * sizeof(GMFloat2);
        //metrics->data_C = gm_malloc(metrics->data_C_size, env);
        metrics->data_C = metrics_mem->malloc_data_C(metrics->data_C_size);
      }
      metrics->data_type_num_values = 4;
    } break;
    /*----------*/
    case GM_DATA_TYPE_TALLY4X2: {
      metrics->data_size = metrics->num_elts_local * sizeof(GMTally4x2);
      //metrics->data = gm_malloc(metrics->data_size, env);
      metrics->data = metrics_mem->malloc_data(metrics->data_size);
      metrics->data_S_size = metrics->num_elts_local * sizeof(GMFloat3);
      //metrics->data_S = gm_malloc(metrics->data_S_size, env);
      metrics->data_S = metrics_mem->malloc_data_S(metrics->data_S_size);
      if (env->sparse()) {
        metrics->data_C_size = metrics->num_elts_local * sizeof(GMFloat3);
        //metrics->data_C = gm_malloc(metrics->data_C_size, env);
        metrics->data_C = metrics_mem->malloc_data_C(metrics->data_C_size);
      }
      metrics->data_type_num_values = 8;
    } break;
    /*----------*/
    default:
      COMET_INSIST(false && "Invalid data_type_id.");
  } /*---switch---*/
}

//=============================================================================
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env) {
  COMET_INSIST(metrics && env);
  COMET_INSIST(metrics->data || ! env->is_proc_active());

  if (! env->is_proc_active()) {
    return;
  }

//  gm_free(metrics->coords_global_from_index,
//          metrics->num_elts_local * sizeof(size_t), env);
//  gm_free(metrics->data, metrics->data_size, env);
//  if (metrics->data_S) {
//    gm_free(metrics->data_S, metrics->data_S_size, env);
//  }
//  if (metrics->data_C) {
//    gm_free(metrics->data_C, metrics->data_C_size, env);
//  }
  *metrics = GMMetrics_null();
}

//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: generic---*/

int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                      size_t index,
                                      int coord_num,
                                      GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(coord_num >= 0 && coord_num < env->num_way());

  size_t result64 = 0;

  switch (env->num_way() + 4 * coord_num) {
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
      COMET_INSIST_INTERFACE(env, false && "Invalid num_way or coord_num.");
  } /*---case---*/

  const int result = (int)result64;
  COMET_ASSERT((size_t)result == result64);

  return result;
}

//-----------------------------------------------------------------------------

void gm_metrics_pad_adjust(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                           GMEnv* env) {
  COMET_INSIST(metrics && metrics_buf && env);
//printf("%i %i\n", (int)metrics_buf->dim0, (int)metrics_buf->dim1);

  if (!(env->is_metric_type_bitwise() && env->is_using_linalg()))
    return;

  // TODO: should more of this be owned by decomp_mgr

  const bool count_2 = env->metric_type() == MetricType::CCC;

  const int pad_adjustment = (count_2 ? 4 : 1) *
    metrics->dm->num_pad_field_local;

  const GMFloat float_pad_adjustment = GMTally1_encode(pad_adjustment, 0);

  // don't use collapse because of overflow for large sizes
  //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
  #pragma omp parallel for schedule(dynamic,1000)
  for (size_t j = 0; j < metrics_buf->dim1; ++j) {
    for (size_t i = 0; i < metrics_buf->dim0; ++i) {

#if 0
        GMTally1 mB00, mB01;
        GMTally1_decode(&mB00, &mB01, metrics_buf->elt<GMTally2x2>(i, j).data[0]);
        GMTally1 mB10, mB11;
        GMTally1_decode(&mB10, &mB11, metrics_buf->elt<GMTally2x2>(i, j).data[1]);
printf("%i  %i %i %i %i\n", env->compute_method(), (int)mB00, (int)mB01, (int)mB10, (int)mB11);
#endif

#ifdef COMET_ASSERTIONS_ON
      const GMTally2x2 old = metrics_buf->elt_const<GMTally2x2>(i, j);
#endif

//printf("%i %zu\n", env->compute_method(), (size_t)metrics_buf->elt<GMTally2x2>(i, j).data[0]);
      metrics_buf->elt<GMTally2x2>(i, j).data[0]
         -= float_pad_adjustment;
//printf("2 %zu\n", (size_t)metrics_buf->elt<GMTally2x2>(i, j).data[0]);

#ifdef COMET_ASSERTIONS_ON
      const GMTally2x2 new_ = metrics_buf->elt_const<GMTally2x2>(i, j);
      COMET_ASSERT(GMTally2x2_get(old, 0, 0) ==
               GMTally2x2_get(new_, 0, 0) + pad_adjustment);
#endif

    } // for j
  }   // for i
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
