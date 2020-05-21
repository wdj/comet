//-----------------------------------------------------------------------------
/*!
 * \file   metrics_3way_indexing.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 3-way, indexing.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_3way_indexing_i_hh_
#define _comet_metrics_3way_indexing_i_hh_

#include "cstdint"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Helper functions for 3-way case.

//-----------------------------------------------------------------------------
// NOTE: the following does not specialize based on part1/2/3.

static int gm_num_section_steps(const CEnv* const env, int part_num) {
  COMET_ASSERT(env && part_num >= 1 && part_num <= 3);
  // Number of section steps to be executed for a given block.

  const bool collapse = ! env->all2all() || env->num_proc_repl() == 1;
  return collapse || part_num == 3 ? 1 : 6;
}

//-----------------------------------------------------------------------------

static int gm_num_sections(const CEnv* const env, int part_num) {
  COMET_ASSERT(env && part_num >= 1 && part_num <= 3);
  // Number of sections the block is divided into.

  return part_num == 3 ? 6 : gm_num_section_steps(env, part_num);
}

//-----------------------------------------------------------------------------

static int gm_num_section_blocks(const CEnv* const env) {
  COMET_ASSERT(env);
  // Total section steps across all blocks, phases.

  const int npv = env->num_proc_vector();

  const bool collapse = ! env->all2all() || env->num_proc_repl() == 1;
  return collapse ? npv*npv - 2*npv + 2 : (npv+1) * (npv+2);
}

//-----------------------------------------------------------------------------

static int gm_section_block_phase_min(const CEnv* const env) {
  COMET_ASSERT(env);

  return (gm_num_section_blocks(env)*env->phase_num()) / env->num_phase();
}

//-----------------------------------------------------------------------------

static int gm_section_block_phase_max(const CEnv* const env) {
  COMET_ASSERT(env);

  return (gm_num_section_blocks(env)*(env->phase_num()+1)) / env->num_phase();
}


//-----------------------------------------------------------------------------

static bool gm_is_section_block_in_phase(const CEnv* const env,
                                        int section_block) {
  COMET_ASSERT(env);
  COMET_ASSERT(section_block >= 0 && section_block < gm_num_section_blocks(env));

  return section_block >= gm_section_block_phase_min(env) &&
         section_block < gm_section_block_phase_max(env);
}

//-----------------------------------------------------------------------------

static bool is_part1(int i_block, int j_block, int k_block) {
  return i_block == j_block && j_block == k_block;
}

//-----------------------------------------------------------------------------

static bool is_part2(int i_block, int j_block, int k_block) {
  return j_block == k_block && !is_part1(i_block, j_block, k_block);
}

//-----------------------------------------------------------------------------

static bool is_part3(int i_block, int j_block, int k_block) {
  return i_block != j_block && j_block != k_block && i_block != k_block;
}

//-----------------------------------------------------------------------------

static int gm_section_axis_part3(int i_block, int j_block, int k_block) {
  // NOTE: this could possibly be implemented somewhat more efficiently.
  /* clang-format off */
  return i_block < j_block && i_block < k_block ? 0 : // i axis
         j_block < i_block && j_block < k_block ? 1 : // j axis
                                                  2;  // k axis
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static int gm_section_num_part3(int i_block, int j_block, int k_block) {
  // NOTE: this could possibly be implemented somewhat more efficiently.
  /* clang-format off */
  return i_block < k_block && k_block < j_block ?    0 :
         i_block < j_block && j_block < k_block ?    1 :
         j_block < i_block && i_block < k_block ?    2 :
         j_block < k_block && k_block < i_block ?    3 :
         k_block < j_block && j_block < i_block ?    4 :
       /*k_block < i_block && i_block < j_block ? */ 5;
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static int gm_J_lo(int section_num, int nvl, int part_num, CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(section_num >= 0 && section_num < 6);
  COMET_ASSERT(nvl >= 0);
  COMET_ASSERT(part_num >= 1 && part_num <= 3);
  const int num_sections = gm_num_sections(env, part_num);
  COMET_ASSERT(section_num >= 0 && section_num <= num_sections);
  COMET_ASSERT(env->num_stage() > 0);
  COMET_ASSERT(env->stage_num() >= 0 && env->stage_num() < env->num_stage());

  const int result = ((env->stage_num() + env->num_stage() * section_num)*nvl) /
                     (num_sections * env->num_stage());

  return result;
}

//-----------------------------------------------------------------------------

static int gm_J_hi(int section_num, int nvl, int part_num, CEnv* env) {
  COMET_ASSERT(env);
  COMET_ASSERT(section_num >= 0 && section_num < 6);
  COMET_ASSERT(nvl >= 0);
  COMET_ASSERT(part_num >= 1 && part_num <= 3);
  const int num_sections = gm_num_sections(env, part_num);
  COMET_ASSERT(section_num >= 0 && section_num <= num_sections);
  COMET_ASSERT(env->num_stage() > 0);
  COMET_ASSERT(env->stage_num() >= 0 && env->stage_num() < env->num_stage());

  const int result = ((env->stage_num() + 1 + env->num_stage() * section_num)*nvl) /
                     (num_sections * env->num_stage());

  return result;
}

//=============================================================================
// GMSectionInfo.

typedef struct {
  bool is_part1;
  bool is_part2;
  bool is_part3;
  bool sax0; // kij
  bool sax1; // ijk
  bool sax2; // jki
  int part_num;
  int section_axis;
  int section_num;
  int i_lb;
  int j_lb;
  int k_lb;
  int i_ub;
  int j_ub;
  int k_ub;
  int J_lb;
  int J_ub;
  int num_vector_local;
  bool no_perm() const {return !is_part3;}
  template<typename T> T perm0(T v0, T v1, T v2) const {
    return no_perm() ?   v0 :
           sax1      ?   v0 :
           sax0      ?   v2 :
       /*  sax2      ?*/ v1;
  }
  template<typename T> T perm1(T v0, T v1, T v2) const {
    return no_perm() ?   v1 :
           sax1      ?   v1 :
           sax0      ?   v0 :
       /*  sax2      ?*/ v2;
  }
  template<typename T> T perm2(T v0, T v1, T v2) const {
    return no_perm() ?   v2 :
           sax1      ?   v2 :
           sax0      ?   v1 :
       /*  sax2      ?*/ v0;
  }
  template<typename T> T unperm0(T v0, T v1, T v2) const {
    return no_perm() ?   v0 :
           sax0      ?   v1 :
           sax1      ?   v0 :
       /*  sax2      ?*/ v2;
  }
  template<typename T> T unperm1(T v0, T v1, T v2) const {
    return no_perm() ?   v1 :
           sax0      ?   v2 :
           sax1      ?   v1 :
       /*  sax2      ?*/ v0;
  }
  template<typename T> T unperm2(T v0, T v1, T v2) const {
    return no_perm() ?   v2 :
           sax0      ?   v0 :
           sax1      ?   v2 :
       /*  sax2      ?*/ v1;
  }
} GMSectionInfo;

//-----------------------------------------------------------------------------

static void GMSectionInfo_create(
  GMSectionInfo* si,
  int i_block,
  int j_block,
  int k_block,
  int section_step,
  int num_vector_local,
  CEnv* env) {
  COMET_INSIST(si && env);
  COMET_INSIST(i_block >= 0 && i_block < env->num_block_vector());
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env->num_block_vector());
  COMET_INSIST(num_vector_local >= 0);

  si->num_vector_local = num_vector_local;

  si->is_part1 = is_part1(i_block, j_block, k_block);
  si->is_part2 = is_part2(i_block, j_block, k_block);
  si->is_part3 = is_part3(i_block, j_block, k_block);

  si->part_num = si->is_part1 ? 1 :
                 si->is_part2 ? 2 : 3;

  const int num_section_steps = gm_num_section_steps(env, si->part_num);
  COMET_INSIST(section_step>=0);
  COMET_INSIST(section_step<num_section_steps);

  si->section_axis =
    ! si->is_part3 ? 1 : // j axis
    gm_section_axis_part3(i_block, j_block, k_block);

  si->section_num = ! si->is_part3 ? section_step :
                    gm_section_num_part3(i_block, j_block, k_block);

  // Define bounding box containing region to be computed.

  si->J_lb = gm_J_lo(si->section_num, num_vector_local, si->part_num, env);
  si->J_ub = gm_J_hi(si->section_num, num_vector_local, si->part_num, env);

  si->i_lb = si->section_axis == 0 ? si->J_lb : 0;
  si->j_lb = si->section_axis == 1 ? si->J_lb : 0;
  si->k_lb = si->section_axis == 2 ? si->J_lb : 0;

  si->i_ub = si->section_axis == 0 ? si->J_ub : num_vector_local;
  si->j_ub = si->section_axis == 1 ? si->J_ub : num_vector_local;
  si->k_ub = si->section_axis == 2 ? si->J_ub : num_vector_local;

  si->sax0 = si->section_axis == 0;
  si->sax1 = si->section_axis == 1;
  si->sax2 = si->section_axis == 2;
}

/*
- I_max = is_part1 ? J : nvl; // XXX can work same way if permuted or not
- K_min = is_part3 ? 0 : J + 1; // XXX can work same way if permuted or not
- put in functions for permuted (I, K) and nonpermuted (i, k)
- store I_ub, etc.
- should we always permute axes (I think not - perhaps only if parallel all2all)
- should we always slice into 6 sections (?)
- do lb/ub values apply for part1/2 - is there a permutation issue - ? OK
- should this be cognizant of all2all value

- * deploy section_num usage for part1/2
*/

//-----------------------------------------------------------------------------

static void GMSectionInfo_destroy(
  GMSectionInfo* si,
  CEnv* env) {
}

//-----------------------------------------------------------------------------

static int GMSectionInfo_k_min(
  GMSectionInfo* si,
  int j,
  CEnv* env) {
  COMET_INSIST(si && env);
  COMET_ASSERT(j >= 0 && j < si->num_vector_local);

  return si->is_part3 ? si->k_lb : j + 1;
}

//-----------------------------------------------------------------------------

static int GMSectionInfo_i_max(
  GMSectionInfo* si,
  int j,
  CEnv* env) {
  COMET_INSIST(si && env);
  COMET_ASSERT(j >= 0 && j < si->num_vector_local);

  return si->is_part1 ? j : si->i_ub;
}

//=============================================================================
// Accessors: indexing: global coord from (contig) index: 3-way.

static int GMMetrics_iG_from_index_3(GMMetrics* metrics, size_t index,
  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index >= 0 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  const GMMetrics::Coords_t i = metrics->coords_values(index) %
                                metrics->num_vector;

  return safe_cast<int>(i);
}

//-----------------------------------------------------------------------------

static int GMMetrics_jG_from_index_3(GMMetrics* metrics, size_t index,
  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index >= 0 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  const GMMetrics::Coords_t j = (metrics->coords_values(index) /
                                 metrics->num_vector) % metrics->num_vector;

  return safe_cast<int>(j);
}

//-----------------------------------------------------------------------------

static int GMMetrics_kG_from_index_3(GMMetrics* metrics, size_t index,
  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index >= 0 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  const GMMetrics::Coords_t k = (metrics->coords_values(index) /
                                 metrics->num_vector) / metrics->num_vector;

  return safe_cast<int>(k);
}

//=============================================================================
// Accessors: indexing: (contig) index from coord, 3-way.

//-----------------------------------------------------------------------------
/// \brief Helper: number of elts in part of a trapezoid, cut orthog to j axis.

static size_t trapezoid_size(int j, int nvl) {
  return ( j *(size_t) (j-1) *(size_t) (3*nvl-2*j-2) ) / 6;
}

//-----------------------------------------------------------------------------
/// \brief Helper: number of elts in part of a triangle, cut orthog to j axis.

static size_t triangle_size(int j, int nvl) {
  return triangle_index(nvl) - triangle_index(nvl-j);
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, main diag block.

static size_t Metrics_index_3_part1(GMMetrics& metrics,
  int i, int j, int k, int i_block, int j_block, int k_block, CEnv& env) {
  const int nvl = metrics.num_vector_local;

  const int num_section_steps = gm_num_section_steps(&env, 1);
  const int section_num = (j * num_section_steps) / nvl;
  // TODO: make this test work for !all2all case.
  COMET_ASSERT(metrics.is_section_num_valid_part1_[section_num] || !env.all2all());

  const int64_t elts_offset = metrics.index_offset_section_part1_[section_num];

  /* clang-format off */
  const int64_t index = elts_offset +
                        i +
                        (k-j-1)*(size_t)j +
                        trapezoid_size(j, nvl);
  /* clang-format on */

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Faster version of true mod, needed for special situation.

static int mod_fast(int i, int n) {
  return (i + n) % n;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, diag plane block.

static size_t Metrics_index_3_part2(GMMetrics& metrics,
  int i, int j, int k, int i_block, int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.all2all());
  const int nvl = metrics.num_vector_local;

  const int num_section_steps = gm_num_section_steps(&env, 2);
  const int section_num = (j * num_section_steps) / nvl;
  COMET_ASSERT(metrics.is_section_num_valid_part2_[section_num]);

  const int64_t elts_offset = metrics.index_offset_section_part2_[section_num];

  const int num_block = env.num_block_vector();
  const int j_i_offset = mod_fast(j_block - i_block, num_block);
  const int block_num_part2 = j_i_offset - 1
      - metrics.phase_block_start_part2_[section_num];

  // Packing offset for multiple section blocks for this proc_r and section_num
  const int num_proc_r = env.num_proc_repl();
  const int blocks_offset = block_num_part2 / num_proc_r;

  const size_t section_size = metrics.section_size_part2_[section_num];

  // Ordering: outer loop is section num, inner loop is block num.

  /* clang-format off */
  const int64_t index = elts_offset +
                        i + nvl*(
                        (k-j-1) +
                        triangle_size(j, nvl) + section_size*(
                        blocks_offset
                        ));
 /* clang-format on */

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, off-diag block.

static size_t Metrics_index_3_part3(GMMetrics& metrics,
  int i, int j, int k, int i_block, int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.all2all());
  const int nvl = metrics.num_vector_local;

  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int64_t elts_offset = metrics.index_offset_part3_;

  const int num_block = env.num_block_vector();
  const int j_i_offset = mod_fast(j_block - i_block, num_block);
  const int k_i_offset = mod_fast(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset))
    - metrics.phase_block_start_part3_;

  // Packing offset for multiple blocks for this proc_r
  const int num_proc_r = env.num_proc_repl();
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int section_axis = gm_section_axis_part3(i_block, j_block, k_block);
  const int J_lo = metrics.J_lo_part3_[section_num];
  const int J_wi = metrics.J_wi_part3_[section_num];

  const int I = section_axis == 1 ? i :
                section_axis == 0 ? k :
                                    j;

  const int J = section_axis == 1 ? j :
                section_axis == 0 ? i :
                                    k;

  const int K = section_axis == 1 ? k :
                section_axis == 0 ? j :
                                    i;

  /* clang-format off */
  const int64_t index = elts_offset +
                        I + nvl * (
                        K + nvl * (
                        J - J_lo + J_wi * (
                        (int64_t)blocks_offset
                        )));
  /* clang-format on */

#if 0
  const int64_t index = elts_offset +
                        i - ( section_axis == 0 ? J_lo : 0 ) +
                            ( section_axis == 0 ? J_wi : nvl ) * (
                        k - ( section_axis == 2 ? J_lo : 0 ) +
                            ( section_axis == 2 ? J_wi : nvl ) * (
                        j - ( section_axis == 1 ? J_lo : 0 ) +
                            ( section_axis == 1 ? J_wi : nvl ) * (
                        (int64_t)blocks_offset
        )));
#endif

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coordinates to index, 3-way case, not permuted/cached.

static size_t Metrics_index_3(GMMetrics& metrics, int i, int j, int k,
  int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  COMET_ASSERT(i >= 0 && j >= 0 && k >= 0);
  COMET_ASSERT(j_block >= 0 && j_block < env.num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env.num_block_vector());
  COMET_ASSERT(! (env.proc_num_vector() == j_block &&
                  env.proc_num_vector() != k_block));
  COMET_ASSERT(! (env.proc_num_vector() == k_block &&
                  env.proc_num_vector() != j_block));
  // WARNING: these conditions are not exhaustive.

  const int i_block = env.proc_num_vector();

  const int64_t index = is_part1(i_block, j_block, k_block) ?
    Metrics_index_3_part1(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
                        is_part2(i_block, j_block, k_block) ?
    Metrics_index_3_part2(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
    Metrics_index_3_part3(metrics, i, j, k,
                                i_block, j_block, k_block, env);

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  COMET_ASSERT((size_t)GMMetrics_iG_from_index_3(&metrics, index, &env) ==
           i + i_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT((size_t)GMMetrics_jG_from_index_3(&metrics, index, &env) ==
           j + j_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT((size_t)GMMetrics_kG_from_index_3(&metrics, index, &env) ==
           k + k_block * (size_t)metrics.num_vector_local);

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, main diag block, permuted.

static size_t Metrics_index_3_part1_permuted(GMMetrics& metrics,
    int I, int J, int K, int i_block, int j_block, int k_block, CEnv& env) {

  return Metrics_index_3_part1(metrics, I, J, K,
                               i_block, j_block, k_block, env);
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, diag plane block, permuted.

static size_t Metrics_index_3_part2_permuted(GMMetrics& metrics,
    int I, int J, int K, int i_block, int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.all2all());

  return Metrics_index_3_part2(metrics, I, J, K,
                               i_block, j_block, k_block, env);
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coords to index, 3-way case, off-diag block, permuted.

static size_t Metrics_index_3_part3_permuted(GMMetrics& metrics,
    int I, int J, int K, int i_block, int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.all2all());
  const int nvl = metrics.num_vector_local;

  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int64_t elts_offset = metrics.index_offset_part3_;

  const int num_block = env.num_block_vector();
  const int j_i_offset = mod_fast(j_block - i_block, num_block);
  const int k_i_offset = mod_fast(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset))
    - metrics.phase_block_start_part3_;

  const int num_proc_r = env.num_proc_repl();
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int J_lo = metrics.J_lo_part3_[section_num];
  const int J_wi = metrics.J_wi_part3_[section_num];

  /* clang-format off */

  const int64_t index = elts_offset +
                        I + nvl * (
                        K + nvl * (
                        J - J_lo + J_wi * (
                        (int64_t)blocks_offset
                        )));

  /* clang-format on */

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Convert elt coordinates to index, 3-way case, permuted.

static size_t Metrics_index_3_permuted_(GMMetrics& metrics,
    int I, int J, int K, int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  COMET_ASSERT(I >= 0 && J >= 0 && K >= 0);
  COMET_ASSERT(j_block >= 0 && j_block < env.num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env.num_block_vector());
  COMET_ASSERT(! (env.proc_num_vector() == j_block &&
              env.proc_num_vector() != k_block));
  COMET_ASSERT(! (env.proc_num_vector() == k_block &&
              env.proc_num_vector() != j_block));
  // WARNING: these conditions are not exhaustive.

  const int i_block = env.proc_num_vector();

  int64_t index = is_part1(i_block, j_block, k_block) ?
    Metrics_index_3_part1_permuted(metrics, I, J, K,
                                   i_block, j_block, k_block, env) :
                 is_part2(i_block, j_block, k_block) ?
    Metrics_index_3_part2_permuted(metrics, I, J, K,
                                   i_block, j_block, k_block, env) :
    Metrics_index_3_part3_permuted(metrics, I, J, K,
                                   i_block, j_block, k_block, env);

  COMET_ASSERT(index >= 0 && index < (int64_t)metrics.num_elts_local);

#ifdef COMET_ASSERTIONS_ON

  int section_step = 0;

  if (is_part3(i_block, j_block, k_block)) {
    section_step = 0;
  } else if (is_part1(i_block, j_block, k_block)) {
    const int num_section_steps = gm_num_section_steps(&env, 1);
    const int section_num = (J * num_section_steps) / metrics.num_vector_local;
    section_step = section_num;
  } else {
    const int num_section_steps = gm_num_section_steps(&env, 2);
    const int section_num = (J * num_section_steps) / metrics.num_vector_local;
    section_step = section_num;
  }

  GMSectionInfo si;
  GMSectionInfo_create(&si, i_block, j_block, k_block, section_step,
                       metrics.num_vector_local, &env);

  const int i = si.unperm0(I, J, K);
  const int j = si.unperm1(I, J, K);
  const int k = si.unperm2(I, J, K);

  GMSectionInfo_destroy(&si, &env);

  COMET_ASSERT((size_t)GMMetrics_iG_from_index_3(&metrics, index, &env) ==
           i + i_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT((size_t)GMMetrics_jG_from_index_3(&metrics, index, &env) ==
           j + j_block * (size_t)metrics.num_vector_local);

  COMET_ASSERT((size_t)GMMetrics_kG_from_index_3(&metrics, index, &env) ==
           k + k_block * (size_t)metrics.num_vector_local);
#endif

  return (size_t)index;
}

//-----------------------------------------------------------------------------
/// \brief Helper struct to speed up index calculations.

class MetricsIndexCache  {
  bool is_initialized_;
  int I_;
  int K_;
  size_t index_;
public:
  size_t index(GMMetrics& metrics, int I, int J, int K,
               int j_block, int k_block, CEnv& env) {

    const bool is_same_K_as_prev = (K == K_);

    if (is_initialized_ && is_same_K_as_prev) {
      // Fast calculation.
      index_ += I - I_;
      // Update the cache.
      I_ = I;
    } else {
      // Slow calculation.
      index_ = Metrics_index_3_permuted_(metrics, I, J, K, j_block, k_block,
        env);
      // Update the cache.
      I_ = I;
      K_ = K;
      is_initialized_ = true;
    }
    return index_;
  }
};

//-----------------------------------------------------------------------------
/// \brief Convert elt coordinates to index, 3-way case, permuted/cached.

static size_t Metrics_index_3( GMMetrics& metrics, int I, int J, int K,
    int j_block, int k_block, MetricsIndexCache& index_cache, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  COMET_ASSERT(I >= 0 && J >= 0 && K >= 0);
  COMET_ASSERT(j_block >= 0 && j_block < env.num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env.num_block_vector());
  COMET_ASSERT(! (env.proc_num_vector() == j_block &&
                  env.proc_num_vector() != k_block));
  COMET_ASSERT(! (env.proc_num_vector() == k_block &&
                  env.proc_num_vector() != j_block));
  // WARNING: these conditions are not exhaustive.

  return index_cache.index(metrics, I, J, K, j_block, k_block, env);;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_3way_indexing_i_hh_

//-----------------------------------------------------------------------------
