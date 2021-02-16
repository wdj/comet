//-----------------------------------------------------------------------------
/*!
 * \file   metrics.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the calculated metrics output from the methods.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

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
// Helper class for memory.

MetricsMem::MetricsMem(CEnv* env) :
  env_(env),
  data_(NULL),
  data_S_(NULL),
  data_C_(NULL),
  data_coords_values_(NULL),
  data_size_(0),
  data_S_size_(0),
  data_C_size_(0),
  data_coords_values_size_(0) {
}

//-----------------------------------------------------------------------------

MetricsMem::~MetricsMem() {

  if (! env_->is_proc_active())
    return;

  if (data_)
    gm_free(data_, data_size_, env_);

  if (data_S_)
    gm_free(data_S_, data_S_size_, env_);

  if (data_C_)
    gm_free(data_C_, data_C_size_, env_);

  if (data_coords_values_)
    gm_free(data_coords_values_, data_coords_values_size_, env_);
}

//-----------------------------------------------------------------------------

void* MetricsMem::malloc_data(size_t data_size) {

  if (! env_->is_proc_active())
    return NULL;

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

void* MetricsMem::malloc_data_S(size_t data_S_size) {

  if (! env_->is_proc_active())
    return NULL;

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

void* MetricsMem::malloc_data_C(size_t data_C_size) {

  if (! env_->is_proc_active())
    return NULL;

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

MetricItemCoords_t* MetricsMem::malloc_data_coords_values(
  size_t data_coords_values_size) {

  if (! env_->is_proc_active())
    return NULL;

  if (!data_coords_values_ ||
      data_coords_values_size > data_coords_values_size_) {
    if (data_coords_values_) {
      gm_free(data_coords_values_, data_coords_values_size_, env_);
    }
    data_coords_values_
       = (size_t*)gm_malloc(data_coords_values_size, env_);
    data_coords_values_size_ = data_coords_values_size;
  }

  return data_coords_values_;
}

//=============================================================================
/// \brief Round-robin-pack m values into n bins, return ith bin size.

static int elt_count_rr_decomp_(int i, int n, int m) {
  return m/n + (i < m % n ? 1 : 0);
}

//=============================================================================
/// \brief Block-pack m values into n bins, return ith bin size.

static int elt_count_block_decomp_(int i, int n, int m) {
  return m / n + (i < m % n);
}

//=============================================================================
// Null object.

GMMetrics GMMetrics_null() {
  GMMetrics result;
  memset((void*)&result, 0, sizeof(result));
  return result;
}

//=============================================================================
/// \brief Calc num metrics to store on this proc (w/o shrink), 2-way case.

void GMMetrics_2way_set_num_metrics_(GMMetrics& metrics, int nvl,
  const CEnv& env) {
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(env.num_way() == NumWay::_2);

  const size_t nchoosek = utils::nchoosek(nvl, env.num_way());

  if (!env.all2all()) {
    // Triangle of values.
    metrics.num_metrics_local = nchoosek;
    return;
  }

  metrics.num_metrics_local = 0;

  const int i_block = env.proc_num_vector();
  const size_t nvlsq = nvl * (size_t)nvl;

  /*---Store the following in this block-row:
      1) strict upper triangular part of main diagonal block
      2) half of the off-diagonal blocks, as a "wrapped rectangle"
    For num_proc_repl > 1, map these blocks, starting at the
    main diagonal block, to procs in round-robin fashion.
    For num_phase > 1, do all this only for a piece of the block row.
  ---*/

  /*===PART A: CALCULATE INDEX SIZE===*/

  const int proc_num_repl = env.proc_num_repl();
  const int num_proc_repl = env.num_proc_repl();

  // PART A.1: (triangle) i_block==j_block part.

  const bool have_main_bdiag = proc_num_repl == 0 &&
                              metrics_bdiag_thisphase_min(env) == 0;

  metrics.num_metrics_local += have_main_bdiag ? nchoosek : 0;
  metrics.index_offset_part2_ = have_main_bdiag ? nchoosek - nvlsq : 0;

  // Block number in this block row where part2 starts (stored
  // off-block-diag blocks).
  metrics.block_min_part2_ = (i_block + metrics_bdiag_thisphase_min(env)) %
    env.num_block_vector();

  // PART A.2: (wrapped rect) i_block!=j_block part.

  // Num blocks computed this block row, this phase, all proc_repl.
  const int num_computed_blocks_thisrow = metrics_num_bdiag_thisphase_thisbrow(env);

//FIXRING - DONE

  // Num blocks computed this block row, this phase, this proc_repl.
  const int num_computed_blocks_thisproc = env.is_comm_ring() ?
//    num_computed_blocks_thisrow / num_proc_repl +
//      (proc_num_repl < num_computed_blocks_thisrow % num_proc_repl) :
   elt_count_block_decomp_(proc_num_repl, num_proc_repl,
     num_computed_blocks_thisrow) :
   elt_count_rr_decomp_(proc_num_repl, num_proc_repl,
     num_computed_blocks_thisrow);

  const int num_computed_offbdiag_blocks_thisproc =
    num_computed_blocks_thisproc - (have_main_bdiag ? 1 : 0);
  metrics.num_metrics_local += num_computed_offbdiag_blocks_thisproc * nvlsq;
}

//=============================================================================

void GMMetrics_3way_set_num_metrics_(GMMetrics& metrics, int nvl, CEnv& env) {
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(env.num_block_vector() <= 2 || nvl % 6 == 0);
  COMET_INSIST(env.num_way() == NumWay::_3);

  if (!env.all2all()) {
    metrics.num_metrics_local = utils::nchoosek(nvl, env.num_way());
    return;
  }

  metrics.num_metrics_local = 0;

  //---Fused counter for section_num and block_num, same across all procs.
  int section_block_num = 0;

  //---Compute size part 1: (tetrahedron) i_block==j_block==k_block part.

  const int num_section_steps_1 = gm_num_section_steps(&env, 1);
  for (int section_step=0; section_step<num_section_steps_1; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 1, &env);
    const int J_hi = gm_J_hi(section_num, nvl, 1, &env);
    const int64_t trap_size_lo = trapezoid_size(J_lo, nvl);
    const int64_t trap_size_hi = trapezoid_size(J_hi, nvl);
    COMET_INSIST(trap_size_hi >= trap_size_lo && "Error in sizes calculation.");
    //---Absorb size_lo into offset for speed in indexing function.
    metrics.index_offset_section_part1_[section_num]
      = (int64_t)metrics.num_metrics_local - trap_size_lo;
    if (gm_is_section_block_in_phase(&env, section_block_num)) {
      if (gm_proc_r_active(section_block_num, &env)) {
        //---Elements in slice of trapezoid.
        const int64_t elts_local = trap_size_hi - trap_size_lo;
        COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
        metrics.num_metrics_local += elts_local;
        metrics.is_section_num_valid_part1_[section_num] = (elts_local != 0);
      } // if
    } // if
    ++section_block_num;
  } // section_step
  metrics.index_offset_part2_ = metrics.num_metrics_local;

  //---Compute size part 2: (triang prisms) i_block!=j_block==k_block part.

  const int num_block = env.num_block_vector();
  const int num_section_steps_2 = gm_num_section_steps(&env, 2);
  for (int section_step=0; section_step<num_section_steps_2; ++section_step) {
    //---Get slice bounds.
    const int section_num = section_step;
    const int J_lo = gm_J_lo(section_num, nvl, 2, &env);
    const int J_hi = gm_J_hi(section_num, nvl, 2, &env);
    const int64_t triang_size_lo = triangle_size(J_lo, nvl);
    const int64_t triang_size_hi = triangle_size(J_hi, nvl);
    //---Absorb size_lo into offset for speed in indexing function.
    metrics.index_offset_section_part2_[section_num]
      = (int64_t)metrics.num_metrics_local - (int64_t)nvl*(int64_t)triang_size_lo;
    metrics.section_size_part2_[section_num] = triang_size_hi - triang_size_lo;
    int block_num_part2 = 0;
    bool is_phase_block_start_set = false;
    //---Loop over blocks for part2.
    for (int j_i_offset=1; j_i_offset<num_block; ++j_i_offset) {
      if (gm_is_section_block_in_phase(&env, section_block_num)) {
        if (!is_phase_block_start_set) {
          metrics.phase_block_start_part2_[section_num] = block_num_part2;
          is_phase_block_start_set = true;
        }
        if (gm_proc_r_active(section_block_num, &env)) {
          //---Elements in slice of triang prism.
          const int64_t elts_local = (int64_t)nvl *
                                     (triang_size_hi - triang_size_lo);
          COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
          metrics.num_metrics_local += elts_local;
          metrics.is_section_num_valid_part2_[section_num] = (elts_local != 0);
        } // if
      } // if
      ++section_block_num;
      ++block_num_part2;
    }
  } // section_step
  metrics.index_offset_part3_ = metrics.num_metrics_local;

  //---Compute size part 3: (block sections) i_block!=j_block!=k_block part.

  //---Loop over block for part3.
  const int num_section_steps_3 = gm_num_section_steps(&env, 3); // = 1
  for (int section_step=0; section_step<num_section_steps_3; ++section_step) {
    const int i_block = env.proc_num_vector();
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
        const int J_lo = gm_J_lo(section_num, nvl, 3, &env);
        const int J_hi = gm_J_hi(section_num, nvl, 3, &env);
        if (gm_is_section_block_in_phase(&env, section_block_num)) {
          if (!is_phase_block_start_set) {
            metrics.phase_block_start_part3_ = block_num_part3;
            is_phase_block_start_set = true;
          }
          if (gm_proc_r_active(section_block_num, &env)) {
            //---Elements in slice of block/cube.
            const int64_t elts_local = (int64_t)nvl * (int64_t)nvl *
                                       (int64_t)(J_hi - J_lo);
            COMET_INSIST(elts_local >= 0 && "Error in sizes calculation.");
            metrics.num_metrics_local += elts_local;
          } // if
        } // if
        ++section_block_num;
        ++block_num_part3;
      }
    }
  } // section_step
}

//=============================================================================

void GMMetrics_set_num_metrics(GMMetrics& metrics, int nvl, CEnv& env) {
  COMET_INSIST(nvl >= 0);

  if (env.num_way() == NumWay::_2)
    GMMetrics_2way_set_num_metrics_(metrics, nvl, env);
  else
    GMMetrics_3way_set_num_metrics_(metrics, nvl, env);
}

//=============================================================================
// Metrics pseudo-constructor.

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      GMDecompMgr* dm,
                      MetricsMem* metrics_mem,
                      CEnv* env) {
  COMET_INSIST(metrics && dm && env);

  if (! env->is_proc_active())
    return;

  //--------------------
  // Perform checks.
  //--------------------

  COMET_INSIST_INTERFACE(env, (env->metric_type() != MetricType::DUO ||
                         env->sparse()) && "DUO method requires sparse input.");

  COMET_INSIST_INTERFACE(env, (env->all2all() || env->num_proc_repl() == 1) &&
          "Multidim parallelism only available for all2all case");

  // (The following less important cases are not yet tested.)

  COMET_INSIST_INTERFACE(env, (env->all2all() ||
                    dm->num_field == dm->num_field_active)
                    && "This case currently not supported.");

  COMET_INSIST_INTERFACE(env, (env->is_compute_method_gpu() ||
                    dm->num_field == dm->num_field_active)
                    && "This case currently not supported.");

  COMET_INSIST_INTERFACE(env, dm->num_vector_local >= (size_t)env->num_way()
        && "Currently require number of vecs on a proc to be at least num-way");

  COMET_INSIST_INTERFACE(env, env->stage_num() >= 0 &&
                              env->stage_num() < env->num_stage() &&
                              "Invalid stage number specified.");

  COMET_INSIST_INTERFACE(env, env->phase_num() >= 0 &&
                              env->phase_num() < env->num_phase() &&
                              "Invalid phase number specified.");

  GMMetrics_ccc_check_size_nofp_2(metrics, env);
  GMMetrics_ccc_check_size_nofp_3(metrics, env);

  COMET_INSIST_INTERFACE(env, (env->num_stage() == 1 ||
                               env->num_way() != NumWay::_2) &&
                    "Staged computations not allowed for 2-way case.");

  COMET_INSIST_INTERFACE(env, (env->num_stage() == 1 || env->all2all()) &&
                    "Staged computations not allowed for non-all2all case.");

  COMET_INSIST_INTERFACE(env, (env->num_phase() == 1 || env->all2all()) &&
                    "Phased computations not allowed for non-all2all case.");

  if (env->num_way() == NumWay::_2)
    COMET_INSIST_INTERFACE(env,
      env->num_phase() <= 1 + env->num_block_vector() / 2 &&
      "num_phase must be at most 1 + num_proc_vector/2.");

  if (env->num_way() == NumWay::_3) {
    COMET_INSIST_INTERFACE(env, env->num_phase() <= gm_num_section_blocks(env)
     && "num_phase must be at most (num_proc_vector+1)*(num_proc_vector+2)/2.");

    //---Make the following assumption to greatly simplify calculations.
    COMET_INSIST_INTERFACE(env,
      (env->num_block_vector()<3 || dm->num_vector_local % 6 == 0) &&
      "3-way all2all case requires num vectors per proc divisible by 6.");
  }

  //--------------------
  // Initializations.
  //--------------------

  *metrics = GMMetrics_null();

  metrics->dm = dm;
  metrics->data_type_id = data_type_id;

  metrics->num_field = dm->num_field;
  metrics->num_field_active = dm->num_field_active;
  metrics->num_field_local = dm->num_field_local;
  metrics->num_vector_local = dm->num_vector_local;
  metrics->num_vector_active = dm->num_vector_active;

  metrics->index_offset_part2_ = 0;
  metrics->index_offset_part3_ = 0;
  metrics->recip_m = 1. / metrics->num_field_active;
  metrics->block_min_part2_ = 0;
  for (int i=0; i<6; ++i) {
    metrics->index_offset_section_part1_[i] = 0;
    metrics->index_offset_section_part2_[i] = 0;
    metrics->is_section_num_valid_part1_[i] = false;
    metrics->is_section_num_valid_part2_[i] = false;
    metrics->phase_block_start_part2_[i] = 0;
  }
  metrics->phase_block_start_part3_ = 0;
  // TODO: make this better.
  metrics->num_steps_2way =
    env->num_way() == NumWay::_2 ? metrics_num_steps_2way(*env) : 0;

  metrics->num_vector = dm->num_vector;

  metrics->num_metrics_local = 0;
  metrics->num_metric_items_local = 0;
  metrics->num_metric_items_local_allocated = 0;
  metrics->num_metric_items_local_computed = 0;

  //--------------------
  // Calcuate number of metrics to be computed.
  //--------------------

  GMMetrics_set_num_metrics(*metrics, metrics->num_vector_local, *env);

  metrics->num_metric_items_local = metrics->num_metrics_local *
                                    env->num_metric_items_per_metric();

  //--------------------
  // Allocations: data.
  //--------------------

  metrics->data_elt_size = env->metric_item_size();

  metrics->num_metric_items_local_allocated =
    env->apply_shrink(metrics->num_metric_items_local);

  const size_t data_size = metrics->num_metric_items_local_allocated *
    env->metric_item_size();
//printf("metrics->num_metric_items_local %zu\n", metrics->num_metric_items_local);
//printf("metrics->num_metric_items_local_allocated %zu\n", metrics->num_metric_items_local_allocated);
//printf("data_size %zu\n", data_size);

  //data_size = metrics->num_metrics_local * metrics->data_elt_size;

  metrics->data = metrics_mem->malloc_data(data_size);

  //--------------------
  // Allocations: coords.
  //--------------------

  metrics->data_coords_values_ = metrics_mem->malloc_data_coords_values(
    metrics->num_metric_items_local_allocated * sizeof(MetricItemCoords_t));

  //--------------------
  // Allocations: data_S, data_C
  //--------------------

  switch (data_type_id) {
    //----------
    case GM_DATA_TYPE_FLOAT:
      break;
    //----------
    case GM_DATA_TYPE_TALLY2X2: {

      if (!env->is_threshold_tc()) {

        metrics->data_S_elt_size = sizeof(GMFloat2);
        const size_t data_S_size = metrics->num_metrics_local *
                               metrics->data_S_elt_size;
        metrics->data_S = metrics_mem->malloc_data_S(data_S_size);

         metrics->data_C_elt_size = sizeof(GMFloat2);
        if (env->sparse()) {
          const size_t data_C_size = metrics->num_metrics_local *
                                    metrics->data_C_elt_size;
          metrics->data_C = metrics_mem->malloc_data_C(data_C_size);
        }

      }

    } break;
    //----------
    case GM_DATA_TYPE_TALLY4X2: {

      if (!env->is_threshold_tc()) {

        metrics->data_S_elt_size = sizeof(GMFloat3);
        const size_t data_S_size = metrics->num_metrics_local *
                                   metrics->data_S_elt_size;
        metrics->data_S = metrics_mem->malloc_data_S(data_S_size);

        metrics->data_C_elt_size = sizeof(GMFloat3);
        if (env->sparse()) {
          const size_t data_C_size = metrics->num_metrics_local *
                                     metrics->data_C_elt_size;
          metrics->data_C = metrics_mem->malloc_data_C(data_C_size);
        }

      }

    } break;
    //----------
    default:
      COMET_INSIST(false && "Invalid data_type_id.");
  } // switch

  //--------------------
  // Set coords.
  //--------------------

  if (env->is_shrink())
    return;

  // TODO: put the following in its own function.

  const int num_block = env->num_block_vector();
  const int i_block = env->proc_num_vector();
  const int nvl = metrics->num_vector_local;
  const size_t nvlsq = nvl * (size_t)nvl;

  /*==================================================*/
  if (env->num_way() == NumWay::_2 && env->all2all()) {
  /*==================================================*/

    const int proc_num_repl = env->proc_num_repl();
    const bool have_main_bdiag = proc_num_repl == 0 &&
//                                gm_bdiag_computed_min(env) == 0;
                                metrics_bdiag_thisphase_min(*env) == 0;

    const int num_computed_blocks_thisrow = metrics_num_bdiag_thisphase_thisbrow(*env);

    // PART C.1: (triangle) i_block==j_block part.

    // Running tally of index into metrics array.
    size_t index = 0;

    if (have_main_bdiag) {
      const int bdiag = 0;
      no_unused_variable_warning(bdiag);
      COMET_ASSERT(env->proc_num_repl() == 0);
//FIXRING - DONE
      COMET_ASSERT(metrics_is_proc_repl_active(*metrics, bdiag, *env));
      // WARNING: no omp pragma here because index++ is sequential.
      for (int j = 0; j < nvl; ++j) {
        const size_t jG = j + nvl * i_block;
        for (int i = 0; i < j; ++i) {
          const size_t iG = i + nvl * i_block;
          COMET_ASSERT(Metrics_index_2_part1(*metrics, i, j, i_block, *env) ==
                       index);
          metrics->data_coords_values_[index++] =
            CoordsInfo::set(iG, jG, *metrics, *env);
        }
      }
    } // if

    // PART C.2: (wrapped rectangle) i_block!=j_block part.

    // NOTE: this loops over all proc_repl cases, not just this
    // one (= "active"), could be inefficient.
    const int bdiag_beg = metrics_bdiag_thisphase_min(*env);
    const int bdiag_end = bdiag_beg + num_computed_blocks_thisrow;
    for (int bdiag = bdiag_beg; bdiag < bdiag_end; ++bdiag) {
      const int bdiag_offset = bdiag - bdiag_beg;
//FIXRING - DONE
      // NOTE bdiag == 0 case already handled (above).
      if (bdiag == 0 ||
          ! metrics_is_proc_repl_active(*metrics, bdiag_offset, *env))
        continue;
      // "unwrapped" is without circulant pattern wrapping.
      const int j_block_unwrapped = i_block + bdiag;
      // don't use collapse because overflow for large sizes; could use size_t j
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
        const size_t jG_unwrapped = j + j_block_unwrapped * (size_t)nvl;
        const size_t jG = jG_unwrapped % metrics->num_vector;
          const size_t iG = i + nvl * i_block;
          const size_t index_this = index + i + j * (size_t)nvl;
          COMET_ASSERT(index_this>=0 && index_this<metrics->num_metrics_local);
          metrics->data_coords_values_[index_this] =
            CoordsInfo::set(iG, jG, *metrics, *env);
        }
      }
      index += nvlsq;
    } // for bdiag

    // Final check.
    COMET_INSIST(index == metrics->num_metrics_local &&
                 "Error in sizes calculation.");

  /*==================================================*/
  } else if (env->num_way() == NumWay::_3 && env->all2all()) {
  /*==================================================*/

    /*===PART C: SET INDEX===*/

    // Fused counter for section_num and block_num, same across all procs.
    int section_block_num = 0;

    section_block_num = 0;
    size_t index = 0;

    // Set index part 1: (tetrahedron) i_block==j_block==k_block part.

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
            const size_t jG = j + nvl * j_block;
            // don't use collapse because of overflow for large sizes
            //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
            #pragma omp parallel for schedule(dynamic,1000)
            for (int k = j+1; k < nvl; ++k) {
              for (int i = 0; i < j; ++i) {
              const int k_block = i_block;
              const size_t kG = k + nvl * k_block;
                const size_t iG = i + nvl * i_block;
                const size_t index_this = index + i + j*(size_t)(k-(j+1));
                COMET_ASSERT(index_this>=0 &&
                             index_this<metrics->num_metrics_local);
                metrics->data_coords_values_[index_this] =
                //  iG + metrics->num_vector * (jG + metrics->num_vector * (kG));
                  CoordsInfo::set(iG, jG, kG, *metrics, *env);
              }
            }
            index += j * (size_t)(nvl - (j+1));
          }
        } // if
      } // if
      ++section_block_num;
    } // section_step

    // Set index part 2: (triang prisms) i_block!=j_block==k_block part.

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
              const size_t jG = j + nvl * j_block;
              // don't use collapse because of overflow for large sizes
              //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
              #pragma omp parallel for schedule(dynamic,1000)
              for (int k = j+1; k < nvl; ++k) {
                for (int i = 0; i < nvl; ++i) {
                  const int k_block = j_block;
                  const size_t kG = k + nvl * k_block;
                  const size_t iG = i + nvl * i_block;
                  const size_t index_this = index + i + nvl*(size_t)(k-(j+1));
                  COMET_ASSERT(index_this>=0 &&
                               index_this<metrics->num_metrics_local);
                  metrics->data_coords_values_[index_this] =
                  //  iG + metrics->num_vector * (
                  //  jG + metrics->num_vector * (kG));
                    CoordsInfo::set(iG, jG, kG, *metrics, *env);
                } // for i
              } // for k
              index += nvl * (size_t)(nvl - (j+1));
            } // for j
          } // if
        } // if
        ++section_block_num;
      } // for j_i_offset
    } // section_step

    // Set index part 3: (block sections) i_block!=j_block!=k_block part.

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

                    const size_t jG = j + nvl * j_block;
                    const size_t kG = k + nvl * k_block;
                    const size_t iG = i + nvl * i_block;
                    COMET_ASSERT(iG >= 0 && iG < metrics->num_vector);
                    COMET_ASSERT(jG >= 0 && jG < metrics->num_vector);
                    COMET_ASSERT(kG >= 0 && kG < metrics->num_vector);
                    const size_t index_this = index + I + K * (size_t)nvl;
                    COMET_ASSERT(index_this>=0 &&
                             index_this<metrics->num_metrics_local);

                    COMET_ASSERT(Metrics_index_3_part3_permuted(*metrics,
                      I, J, K, i_block, j_block, k_block, *env) == index_this);

                    COMET_ASSERT(Metrics_index_3_part3(*metrics,
                      i, j, k, i_block, j_block, k_block, *env) == index_this);
                    metrics->data_coords_values_[index_this] =
                    //  iG + metrics->num_vector * (
                    //  jG + metrics->num_vector * (kG));
                      CoordsInfo::set(iG, jG, kG, *metrics, *env);
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

    COMET_INSIST(index == metrics->num_metrics_local &&
                 "Error in sizes calculation.");

  /*==================================================*/
  } else if (env->num_way() == NumWay::_2 && ! env->all2all()) {
  /*==================================================*/

    // Need store only strict upper triangular part of matrix.
    size_t index = 0;
    for (int j = 0; j < nvl; ++j) {
      const size_t jG = j + nvl * i_block;
      for (int i = 0; i < j; ++i) {
        const size_t iG = i + nvl * i_block;
        //metrics->data_coords_values_[index++] = iG + metrics->num_vector * jG;
        metrics->data_coords_values_[index++] =
          CoordsInfo::set(iG, jG, *metrics, *env);
      }
    }
    COMET_INSIST(index == metrics->num_metrics_local &&
                 "Error in sizes calculation.");

  /*==================================================*/
  } else { //  if (env->num_way() == NumWay::_3 && ! env->all2all())
  /*==================================================*/

    // Need store only strict interior of tetrahedron.
    size_t index = 0;
    for (int j = 0; j < nvl; ++j) {
      const int j_block = i_block;
      const size_t jG = j + nvl * j_block;
      for (int k = j+1; k < nvl; ++k) {
        const int k_block = i_block;
        const size_t kG = k + nvl * k_block;
        for (int i = 0; i < j; ++i) {
          const size_t iG = i + nvl * i_block;
          COMET_ASSERT(index < metrics->num_metrics_local);
          metrics->data_coords_values_[index++] =
          //  iG + metrics->num_vector * (jG + metrics->num_vector * (kG));
            CoordsInfo::set(iG, jG, kG, *metrics, *env);
        }
      }
    }
    COMET_INSIST(index == metrics->num_metrics_local &&
                 "Error in sizes calculation.");

  }
  
  //--------------------
  // Checks.
  //--------------------

  size_t num_elts = 0;
  COMET_MPI_SAFE_CALL(MPI_Allreduce(&metrics->num_metrics_local, &num_elts, 1,
    MPI_UNSIGNED_LONG_LONG, MPI_SUM, env->comm_repl_vector()));

  if (env->num_stage() == 1 && env->num_phase() == 1 && env->all2all())
    COMET_INSIST(utils::nchoosek(metrics->num_vector, env->num_way()) ==
                 num_elts);

  COMET_INSIST(env->all2all() ||
    utils::nchoosek(metrics->num_vector_local, env->num_way()) ==
                    metrics->num_metrics_local);
}

//=============================================================================
// Metrics pseudo-destructor.

void GMMetrics_destroy(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);
  COMET_INSIST(metrics->data || ! env->is_proc_active());

  if (! env->is_proc_active())
    return;

  *metrics = GMMetrics_null();
}

//=============================================================================
// Accessors: indexing: global coord from (contig) index: generic.

size_t Metrics_coords_getG(GMMetrics& metrics, size_t index, int ijk,
  CEnv& env) {
  COMET_ASSERT(index+1 >= 1 && index < metrics.num_metric_items_local_allocated);
  COMET_ASSERT(ijk >= 0 && ijk < env.num_way());

  const size_t result = CoordsInfo::getG(metrics.coords_value(index), ijk,
    metrics, env);

  return result;
}

//-----------------------------------------------------------------------------

void gm_metrics_pad_adjust(GMMetrics* metrics, MirroredBuf* metrics_buf,
                           CEnv* env, int weight) {
  COMET_INSIST(metrics && metrics_buf && env);

  if (!(env->is_metric_type_bitwise() && env->is_using_linalg()))
    return;

  if (env->is_using_tc())
    return;

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // TODO: should more of this be owned by decomp_mgr

  const int cbpe = env->counted_bits_per_elt();

  const bool is_cbpe_2 = 2 == cbpe;

  const int pad_adjustment = (is_cbpe_2 ? 4 : 1) * weight *
    metrics->dm->num_pad_field_local;

  // don't use collapse because of overflow for large sizes
  //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
  #pragma omp parallel for schedule(dynamic,1000)
  for (size_t j = 0; j < metrics_buf->dim1; ++j) {
    for (size_t i = 0; i < metrics_buf->dim0; ++i) {

#ifdef COMET_ASSERTIONS_ON
      const GMTally2x2 vold = metrics_buf->elt_const<GMTally2x2>(i, j);
#endif

      MFT::subtract(metrics_buf->elt<GMTally2x2>(i, j).data[0], pad_adjustment, 0);

#ifdef COMET_ASSERTIONS_ON
      const GMTally2x2 vnew = metrics_buf->elt_const<GMTally2x2>(i, j);
      COMET_ASSERT(GMTally2x2_get(vold, 0, 0) ==
               GMTally2x2_get(vnew, 0, 0) + pad_adjustment);
#endif

    } // for j
  }   // for i
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
