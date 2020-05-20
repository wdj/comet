//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, indexing.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_2way_indexing_i_hh_
#define _comet_metrics_2way_indexing_i_hh_

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Helper functions for 2-way case.

static int gm_bdiag_computed_max_allphase(CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());

  // Max number of blocks of any block row computed on all phases.
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
// Accessors: indexing: (contig) index from coord, 2-way.

//-----------------------------------------------------------------------------
/// \brief Helper: number of elements in triangle.

static size_t triangle_index(int i) {
  return (i * (size_t)(i-1)) >> 1;
}

//-----------------------------------------------------------------------------

static bool is_part1(int i_block, int j_block) {
  return i_block == j_block;
}

//-----------------------------------------------------------------------------

static bool is_part2(int i_block, int j_block) {
  return !is_part1(i_block, j_block);
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 2-way case, main diagonal block.

static size_t Metrics_index_2_part1(GMMetrics& metrics,
  int i, int j, int j_block, CEnv& env) {
  COMET_ASSERT(is_part1(env.proc_num_vector(), j_block));

  return triangle_index(j) + i;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 2-way case, off-diagonal block.

static size_t Metrics_index_2_part2(GMMetrics& metrics,
  int i, int j, int j_block, CEnv& env) {
  COMET_ASSERT(env.all2all());
  COMET_ASSERT(j_block != env.proc_num_vector());

  const int num_block = env.num_block_vector();
  const int num_proc_r = env.num_proc_repl();
  const int block_min_part2_ = metrics.block_min_part2_;

  /* clang-format off */
  return metrics.index_offset_part2_ +
      i + metrics.num_vector_local * (size_t)(
      j + metrics.num_vector_local * (
      ((j_block - block_min_part2_ + num_block) % num_block) / num_proc_r ));
  /* clang-format on */
}

//-----------------------------------------------------------------------------
/// \brief Convert element coordinates to index, 2-way case.

static size_t Metrics_index_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_2);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || env.proc_num_vector() == j_block);
  COMET_ASSERT(i >= 0 && i < metrics.num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics.num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env.num_block_vector());
  COMET_ASSERT(i < j || j_block != env.proc_num_vector());
  // WARNING: these conditions on j_block are not exhaustive.

  const int i_block = env.proc_num_vector();

  const int64_t index = is_part1(i_block, j_block)
           ? Metrics_index_2_part1(metrics, i, j, j_block, env)
           : Metrics_index_2_part2(metrics, i, j, j_block, env);

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  COMET_ASSERT(metrics.coords_global_from_index[index] % metrics.num_vector ==
           i + i_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT(metrics.coords_global_from_index[index] / metrics.num_vector ==
           j + j_block * (size_t)metrics.num_vector_local);

  return (size_t)index;
}

//=============================================================================
// Accessors: indexing: global coord from (contig) index: 2-way.

static int GMMetrics_iG_from_index_2(GMMetrics* metrics, size_t index,
  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);


  const size_t i = metrics->coords_global_from_index[index] %
                   metrics->num_vector;

  return safe_cast<int>(i);
}

//-----------------------------------------------------------------------------

static int GMMetrics_jG_from_index_2(GMMetrics* metrics, size_t index,
  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  const size_t j = metrics->coords_global_from_index[index] /
                   metrics->num_vector;

  return safe_cast<int>(j);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_2way_indexing_i_hh_

//-----------------------------------------------------------------------------
