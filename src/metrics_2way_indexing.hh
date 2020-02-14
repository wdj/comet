//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, indexing.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_2way_indexing_hh_
#define _comet_metrics_2way_indexing_hh_

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Helper functions for 2-way case---*/

static int gm_bdiag_computed_max_allphase(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  /*---Max number of blocks of any block row computed on all phases---*/
  const int num_block = env->num_block_vector();
  return 1 + num_block / 2;
}

//-----------------------------------------------------------------------------

static int gm_bdiag_computed_min(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  // First block diag computed for this phase (min across vec procs)
  const int max_rectangle_width = gm_bdiag_computed_max_allphase(env);
  return (max_rectangle_width*env->phase_num()) / env->num_phase();
}

//-----------------------------------------------------------------------------

static int gm_bdiag_computed_max(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  // 1 + last block diag computed for this phase (max across vec procs)
  const int max_rectangle_width = gm_bdiag_computed_max_allphase(env);
  return (max_rectangle_width*(env->phase_num()+1)) / env->num_phase();
}

//-----------------------------------------------------------------------------

static int gm_block_computed_this_row_min(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  return gm_bdiag_computed_min(env);
}

//-----------------------------------------------------------------------------

static int gm_block_computed_this_row_max(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  const int num_block = env->num_block_vector();
  const int i_block = env->proc_num_vector();

  const bool is_row_short_by_1 = num_block % 2 == 0 && 2*i_block >= num_block;
  const bool is_last_phase = env->phase_num() == env->num_phase() - 1;

  const int diag_max = gm_bdiag_computed_max(env);

  // 1 + last block diag computed for this phase, all repl procs (this vec proc)
  const int n = is_last_phase && is_row_short_by_1 ? diag_max - 1 : diag_max;
  COMET_ASSERT(n >= 0);
  COMET_ASSERT(n <= num_block);
  return n;
}

//-----------------------------------------------------------------------------

static int gm_blocks_computed_this_row(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  // num block diags computed for this phase, all repl procs (this vec proc)
  const int n = gm_block_computed_this_row_max(env) -
                gm_block_computed_this_row_min(env);
  COMET_ASSERT(n >= 0);
  COMET_ASSERT(n <= env->num_block_vector());
  return n;
}

//=============================================================================
/*---Accessors: indexing: (contig) index from coord, 2-way---*/

static size_t gm_triang_(int i) {
  return (i * (size_t)(i-1)) >> 1;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0);
  COMET_ASSERT(i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0);
  COMET_ASSERT(j < metrics->num_vector_local);
  COMET_ASSERT(i < j);
  COMET_ASSERT(env->proc_num_repl() == 0);

  size_t index = gm_triang_(j) + i;
  COMET_ASSERT(i + metrics->num_vector_local *
               (size_t)env->proc_num_vector() ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  COMET_ASSERT(j + metrics->num_vector_local *
               (size_t)env->proc_num_vector() ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_maindiag_block_(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   CEnv* env) {
  COMET_ASSERT(j_block == env->proc_num_vector());

  return gm_triang_(j) + i;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_offdiag_block_(GMMetrics* metrics,
                                                  int i,
                                                  int j,
                                                  int j_block,
                                                  CEnv* env) {
  COMET_ASSERT(j_block != env->proc_num_vector());

  const int num_block = env->num_block_vector();

  const int num_proc_r = env->num_proc_repl();

  const int block_min = metrics->block_min;

  /* clang-format off */
  return metrics->index_offset_0_ +
      i + metrics->num_vector_local * (size_t)(
      j + metrics->num_vector_local * (
      ((j_block - block_min + num_block) % num_block) / num_proc_r ));
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_2(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0);
  COMET_ASSERT(i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0);
  COMET_ASSERT(j < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0);
  COMET_ASSERT(j_block < env->num_block_vector());
  COMET_ASSERT(i < j || j_block != env->proc_num_vector());
//  COMET_ASSERT(env->proc_num_repl() == 0 ||
//           j_block != env->proc_num_vector() // DEFUNCT
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  const int i_block = env->proc_num_vector();

  size_t index = j_block == i_block
           ? GMMetrics_helper2way_maindiag_block_(metrics, i, j, j_block, env)
           : GMMetrics_helper2way_offdiag_block_(metrics, i, j, j_block, env);

  COMET_ASSERT(index >= 0 && index < metrics->num_elts_local);
  return index;
}

//=============================================================================
//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: 2-way---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);


  const size_t i64 = metrics->coords_global_from_index[index] %
                     metrics->num_vector;
  const int i = (int)i64;
  COMET_ASSERT((size_t)i == i64);

  return i;
}

//-----------------------------------------------------------------------------

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  const size_t j64 = metrics->coords_global_from_index[index] /
                     metrics->num_vector;
  const int j = (int)j64;
  COMET_ASSERT((size_t)j == j64);

  return j;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_2way_indexing_hh_

//-----------------------------------------------------------------------------
