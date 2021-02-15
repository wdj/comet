//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, indexing.
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

#ifndef _COMET_METRICS_2WAY_INDEXING_I_HH_
#define _COMET_METRICS_2WAY_INDEXING_I_HH_

//=============================================================================

namespace comet {

// Assumed ordering of block diagonals is RSP  - proc_repl axis varies fastest,
// then serial step axis, then phase axis.

//-----------------------------------------------------------------------------
/// \brief Num nonzero/stored block diags of full metrics matrix.

//static int gm_bdiag_computed_max_allphase(const CEnv* env) {
//  COMET_ASSERT(env);
//  COMET_ASSERT(env->num_way() == NumWay::_2);
//  COMET_ASSERT(env->all2all());
//
//  // Max number of blocks of any block row computed on all phases.
//  const int num_block = env->num_block_vector();
//  return 1 + num_block / 2;
//}

static int metrics_num_bdiag_(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  // Assuming here the full matrix has a circulant sparsity pattern that stores
  // all unique blocks, under assumption of symmetry of the full matrix.

  // = max (across block rows) of num blocks to be computed for all phases.
  return 1 + env.num_block_vector() / 2;
}


//-----------------------------------------------------------------------------
/// \brief Min limit of block diags to be computed this phase.

//static int gm_bdiag_computed_min(const CEnv* env) {
//  COMET_ASSERT(env);
//  COMET_ASSERT(env->num_way() == NumWay::_2);
//  COMET_ASSERT(env->all2all());
//
//  // First block diag computed for this phase (min across vec procs)
//  const int max_rectangle_width = metrics_num_bdiag_(*env);
//  return (max_rectangle_width*env->phase_num()) / env->num_phase();
//}

static int metrics_bdiag_thisphase_min(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  // First block diag computed for this phase (min across vec procs).
  const int num_bdiag = metrics_num_bdiag_(env);
  return (num_bdiag * env.phase_num()) / env.num_phase();
}

//-----------------------------------------------------------------------------
/// \brief Max limit of block diags to be computed this phase.

//static int gm_bdiag_computed_max(const CEnv* env) {
//  COMET_ASSERT(env);
//  COMET_ASSERT(env->num_way() == NumWay::_2);
//  COMET_ASSERT(env->all2all());
//
//  // (1 + ) last block diag computed for this phase (max across vec procs).
//  const int max_rectangle_width = metrics_num_bdiag_(*env);
//  return (max_rectangle_width*(env->phase_num()+1)) / env->num_phase();
//}

static int metrics_bdiag_thisphase_max(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  // (1 + ) last block diag computed for this phase (max across vec procs).
  const int num_bdiag = metrics_num_bdiag_(env);
  return (num_bdiag * (env.phase_num()+1)) / env.num_phase();
}

//-----------------------------------------------------------------------------
/// \brief Min limit of block diags (blocks) to compute this phase & block row.

static int gm_block_computed_this_row_min(const CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NumWay::_2);
  COMET_ASSERT(env->all2all());

  //return gm_bdiag_computed_min(env);
  return metrics_bdiag_thisphase_min(*env);
}


static int metrics_bdiag_thisphase_thisbrow_min(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  return metrics_bdiag_thisphase_min(env);
}

//-----------------------------------------------------------------------------
/// \brief Max limit of block diags (blocks) to compute this phase & block row.

static int gm_block_computed_this_row_max(const CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(env->num_way() == NumWay::_2);
  COMET_ASSERT(env->all2all());

  const int num_block = env->num_block_vector();
  const int i_block = env->proc_num_vector();

  const bool is_row_short_by_1 = num_block % 2 == 0 && 2*i_block >= num_block;
  const bool is_last_phase = env->phase_num() == env->num_phase() - 1;

  //const int diag_max = gm_bdiag_computed_max(env);
  const int diag_max = metrics_bdiag_thisphase_max(*env);

  // 1 + last block diag computed for this phase, all repl procs (this vec proc)
  const int n = is_last_phase && is_row_short_by_1 ? diag_max - 1 : diag_max;
  COMET_ASSERT(n >= 0);
  COMET_ASSERT(n <= num_block);
  return n;
}


static int metrics_bdiag_thisphase_thisbrow_max(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  const int num_block = env.num_block_vector();
  const int i_block = env.proc_num_vector();

  // The block in the last block diag may be in fact be zero, in consideration
  // of the set of blocks needed to cover all unique under symmetry. This would
  // potentially be encountered in the last phase, so account for this.

  const bool is_row_short_by_1 = num_block % 2 == 0 && 2 * i_block >= num_block;
  const bool is_last_phase = env.phase_num() == env.num_phase() - 1;

  //const int diag_max = gm_bdiag_computed_max(env);
  const int diag_max = metrics_bdiag_thisphase_max(env);

  // 1 + last block diag computed for this phase, all repl procs (this vec proc)
  const int result = is_last_phase && is_row_short_by_1 ?
     diag_max - 1 : diag_max;
  COMET_ASSERT(result >= 0 && result <= num_block);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Number of block diags (blocks) to compute this phase & block row.

//static int gm_blocks_computed_this_row(const CEnv* env) {
//  COMET_ASSERT(env);
//  COMET_ASSERT(env->num_way() == NumWay::_2);
//  COMET_ASSERT(env->all2all());
//
//  // num block diags computed for this phase, all repl procs (this vec proc)
//  const int n = gm_block_computed_this_row_max(env) -
//                gm_block_computed_this_row_min(env);
//  COMET_ASSERT(n >= 0);
//  COMET_ASSERT(n <= env->num_block_vector());
//  return n;
//}

static int metrics_num_bdiag_thisphase_thisbrow(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all());

  const int result = metrics_bdiag_thisphase_thisbrow_max(env) -
                     metrics_bdiag_thisphase_thisbrow_min(env);
  COMET_ASSERT(result >= 0 && result <= env.num_block_vector());
  return result;
}

//-----------------------------------------------------------------------------
/// \brief Number of steps needed fror 2-way methods.

static int metrics_2way_num_steps(const CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.all2all()); // assume this for now.

  const int num_bdiag_thisphase = metrics_bdiag_thisphase_max(env) -
                                  metrics_bdiag_thisphase_min(env);

  const int result = utils::ceil(num_bdiag_thisphase, env.num_proc_repl());

  return result;
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
  const int block_min_part2 = metrics.block_min_part2_;

//FIXRING
  const int block_part_2 =
    ((j_block - block_min_part2 + num_block) % num_block) / num_proc_r;

  /* clang-format off */
  return metrics.index_offset_part2_ +
      i + metrics.num_vector_local * (size_t)(
      j + metrics.num_vector_local * (
      block_part_2));
  /* clang-format on */
}

//-----------------------------------------------------------------------------
/// \brief Convert element coordinates to index, 2-way case.

static size_t Metrics_index_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
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

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_metrics_local);

  COMET_ASSERT(CoordsInfo::getiG(metrics.coords_value(index), metrics, env) ==
           i + i_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT(CoordsInfo::getjG(metrics.coords_value(index), metrics, env) ==
           j + j_block * (size_t)metrics.num_vector_local);

  return (size_t)index;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_METRICS_2WAY_INDEXING_I_HH_

//-----------------------------------------------------------------------------
