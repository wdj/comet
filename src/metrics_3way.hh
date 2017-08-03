/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics_3way.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 3-way.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_metrics_3way_hh_
#define _gm_metrics_3way_hh_

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Helper functions for 3-way case---*/

static bool gm_is_part1(int i_block, int j_block, int k_block) {
  return i_block == j_block && j_block == k_block;
}

/*---------------------------------------------------------------------------*/

static bool gm_is_part3(int i_block, int j_block, int k_block) {
  return i_block != j_block && j_block != k_block && i_block != k_block;
}

/*---------------------------------------------------------------------------*/

static int gm_section_axis_part3(int i_block, int j_block, int k_block) {
  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/
  /* clang-format off */
  return i_block < j_block && i_block < k_block ? 0 : /*---i axis---*/
         j_block < i_block && j_block < k_block ? 1 : /*---j axis---*/
                                                  2;  /*---k axis---*/
  /* clang-format on */
}

/*---------------------------------------------------------------------------*/

static int gm_section_num_part3(int i_block, int j_block, int k_block) {
  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/
  /* clang-format off */
  return i_block < k_block && k_block < j_block ?    0 :
         i_block < j_block && j_block < k_block ?    1 :
         j_block < i_block && i_block < k_block ?    2 :
         j_block < k_block && k_block < i_block ?    3 :
         k_block < j_block && j_block < i_block ?    4 :
       /*k_block < i_block && i_block < j_block ? */ 5;
  /* clang-format on */
}

/*---------------------------------------------------------------------------*/

static int gm_J_lo(int section_num, int nvl, int part_num, GMEnv* env) {
  GMAssert(env);
  GMAssert(section_num >= 0 && section_num < 6);
  GMAssert(nvl >= 0);
  GMAssert(part_num >= 1 && part_num <= 3);
  const int num_sections = GMEnv_num_sections(env, part_num);
  GMAssert(section_num >= 0 && section_num <= num_sections);
  GMAssert(env->num_stage > 0);
  GMAssert(env->stage_num >= 0 && env->stage_num < env->num_stage);

  const int result = ((env->stage_num + env->num_stage * section_num)*nvl) /
                     (num_sections * env->num_stage);

  return result;
}

/*---------------------------------------------------------------------------*/

static int gm_J_hi(int section_num, int nvl, int part_num, GMEnv* env) {
  GMAssert(env);
  GMAssert(section_num >= 0 && section_num < 6);
  GMAssert(nvl >= 0);
  GMAssert(part_num >= 1 && part_num <= 3);
  const int num_sections = GMEnv_num_sections(env, part_num);
  GMAssert(section_num >= 0 && section_num <= num_sections);
  GMAssert(env->num_stage > 0);
  GMAssert(env->stage_num >= 0 && env->stage_num < env->num_stage);

  const int result = ((env->stage_num + 1 + env->num_stage * section_num)*nvl) /
                     (num_sections * env->num_stage);

  return result;
}

/*===========================================================================*/
/*---GMSectionInfo---*/

typedef struct {
  bool is_part1;
  bool is_part2;
  bool is_part3;
  bool sax0;
  bool sax1;
  bool sax2;
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
} GMSectionInfo;

/*---------------------------------------------------------------------------*/

static void GMSectionInfo_create(
  GMSectionInfo* si,
  int i_block,
  int j_block,
  int k_block,
  int section_step,
  int num_vector_local,
  GMEnv* env) {
  GMAssertAlways(si && env);
  GMAssertAlways(i_block >= 0 && i_block < GMEnv_num_block_vector(env));
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssertAlways(num_vector_local >= 0);

  si->num_vector_local = num_vector_local;

  si->is_part1 = gm_is_part1(i_block, j_block, k_block);
  si->is_part3 = gm_is_part3(i_block, j_block, k_block);
  si->is_part2 = (!si->is_part1) && (!si->is_part3);

  si->part_num = si->is_part1 ? 1 :
                 si->is_part2 ? 2 : 3;

  const int num_section_steps = GMEnv_num_section_steps(env, si->part_num);
  GMAssertAlways(section_step>=0);
  GMAssertAlways(section_step<num_section_steps);

  si->section_axis =
    ! si->is_part3 ? 1 /*---j axis---*/ :
    gm_section_axis_part3(i_block, j_block, k_block);

  si->section_num = !si->is_part3 ? section_step :
                    gm_section_num_part3(i_block, j_block, k_block);

  /*---Define bounding box containing region to be computed---*/

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

#if 0

- I_max = is_part1 ? J : nvl; // XXX can work same way if permuted or not
- K_min = is_part3 ? 0 : J + 1; // XXX can work same way if permuted or not
- put in functions for permuted (I, K) and nonpermuted (i, k)
- store I_ub, etc.
- should we always permute axes (I think not - perhaps only if parallel all2all)
- should we always slice into 6 sections (?)
- do lb/ub values apply for part1/2 - is there a permutation issue - ? OK
- should this be cognizant of all2all value

- * deploy section_num usage for part1/2

#endif

/*---------------------------------------------------------------------------*/

static void GMSectionInfo_destroy(
  GMSectionInfo* si,
  GMEnv* env) {
}

/*---------------------------------------------------------------------------*/

static int GMSectionInfo_k_min(
  GMSectionInfo* si,
  int j,
  GMEnv* env) {
  GMAssertAlways(si && env);
  GMAssert(j >= 0 && j < si->num_vector_local);

  return si->is_part3 ? si->k_lb : j + 1;
}

/*---------------------------------------------------------------------------*/

static int GMSectionInfo_i_max(
  GMSectionInfo* si,
  int j,
  GMEnv* env) {
  GMAssertAlways(si && env);
  GMAssert(j >= 0 && j < si->num_vector_local);

  return si->is_part1 ? j : si->i_ub;
}

/*===========================================================================*/
/*---Accessors: indexing: (contig) index from coord, 3-way---*/

/*---------------------------------------------------------------------------*/
/*---elements in a part of a trapezoid, cut orthog to j axis---*/

static size_t gm_trap_size(int j, int nvl) {
  return ( j *(size_t) (j-1) *(size_t) (3*nvl-2*j-2) ) / 6;
}

/*---------------------------------------------------------------------------*/
/*---elements in a part of a triang, cut orthog to j axis---*/

static size_t gm_triang_size(int j, int nvl) {
  return gm_triang_(nvl) - gm_triang_(nvl-j);
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);

  const int nvl = metrics->num_vector_local;

  /* clang-format off */
  const GMInt64 index = i +
                        (k-j-1)*(size_t)j +
                        gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  GMAssert(i + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           (metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector);
  GMAssert(k + metrics->num_vector_local * (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] /
               (metrics->num_vector * (size_t)metrics->num_vector));

  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part1_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int num_section_steps = GMEnv_num_section_steps(env, 1);
  const int section_num = (j * num_section_steps) / nvl;
  GMAssert(metrics->section_num_valid_part1_[section_num]);

  const GMInt64 elts_offset = metrics->index_offset_section_part1_[section_num];

  /* clang-format off */
  const GMInt64 index = elts_offset +
                        i +
                        (k-j-1)*(size_t)j +
                        gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  return index;
  //return GMMetrics_index_from_coord_3(metrics, i, j, k, env);
}

/*---------------------------------------------------------------------------*/
/*---Faster version of true mod needed for special situation---*/

static int gm_mod1_(int i, int n) {
  return (i + n) % n;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part2_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int num_section_steps = GMEnv_num_section_steps(env, 2);
  const int section_num = (j * num_section_steps) / nvl;
  GMAssert(metrics->section_num_valid_part2_[section_num]);

  const GMInt64 elts_offset = metrics->index_offset_section_part2_[section_num];

  const int num_block = GMEnv_num_block_vector(env);
  const int block_num_part2 = gm_mod1_(j_block - i_block, num_block) - 1;
  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part2 / num_proc_r;

  const size_t section_size = metrics->section_size_part2[section_num];

  /* clang-format off */
  const GMInt64 index = elts_offset +
                        i + nvl*(
                        (k-j-1) +
                        gm_triang_size(j, nvl) + section_size*(
                        blocks_offset
                        ));
 /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part3_(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int i_block,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int section_axis = gm_section_axis_part3(i_block, j_block, k_block);
  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int num_block = GMEnv_num_block_vector(env);
  const int j_i_offset = gm_mod1_(j_block - i_block, num_block);
  const int k_i_offset = gm_mod1_(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset));
  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int J_lo = metrics->J_lo_part3_[section_num];
  const int J_wi = metrics->J_wi_part3_[section_num];

  const GMInt64 elts_offset = metrics->index_offset_01_;

  /* clang-format off */
  const GMInt64 index = elts_offset +
                        i - ( section_axis == 0 ? J_lo : 0 ) +
                            ( section_axis == 0 ? J_wi : nvl ) * (
                        k - ( section_axis == 2 ? J_lo : 0 ) +
                            ( section_axis == 2 ? J_wi : nvl ) * (
                        j - ( section_axis == 1 ? J_lo : 0 ) +
                            ( section_axis == 1 ? J_wi : nvl ) * (
                        (GMInt64)blocks_offset
        )));
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_3(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int k,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && j >= 0 && k >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == j_block &&
             GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == k_block &&
             GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  const GMInt64 index = j_block == i_block && k_block == i_block ?
    GMMetrics_helper3way_part1_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper3way_part2_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
    GMMetrics_helper3way_part3_(metrics, i, j, k,
                                i_block, j_block, k_block, env);

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  GMAssert(metrics->coords_global_from_index[index] %
             (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env)) ==
           i + i_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env))) %
               (metrics->num_vector_local * GMEnv_num_block_vector(env)) ==
           j + j_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)GMEnv_num_block_vector(env))) /
               (metrics->num_vector_local * GMEnv_num_block_vector(env)) ==
           k + k_block * (size_t)metrics->num_vector_local);

  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part1_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {

  return GMMetrics_helper3way_part1_(metrics, I, J, K,
                                     i_block, j_block, k_block, env);
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part2_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {

  return GMMetrics_helper3way_part2_(metrics, I, J, K,
                                     i_block, j_block, k_block, env);
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper3way_part3_permuted_(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int i_block,
    int j_block,
    int k_block,
    GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const int num_block = GMEnv_num_block_vector(env);
  const int j_i_offset = gm_mod1_(j_block - i_block, num_block);
  const int k_i_offset = gm_mod1_(k_block - i_block, num_block);
  const int block_num_part3 =
    ((num_block-2) * (k_i_offset - 1)) +
    (j_i_offset - 1 - (j_i_offset > k_i_offset));
  const int num_proc_r = GMEnv_num_proc_repl(env);
  const int blocks_offset = block_num_part3 / num_proc_r;

  const int J_lo = metrics->J_lo_part3_[section_num];
  const int J_wi = metrics->J_wi_part3_[section_num];

  const GMInt64 elts_offset = metrics->index_offset_01_;

  /* clang-format off */

  const GMInt64 index = elts_offset +
                        I + nvl * (
                        K + nvl * (
                        J - J_lo + J_wi * (
                        (GMInt64)blocks_offset
                        )));

  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < (GMInt64)metrics->num_elts_local);

  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_3_permuted(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && J >= 0 && K >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == j_block &&
             GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == k_block &&
             GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  size_t index = j_block == i_block && k_block == i_block ?
    GMMetrics_helper3way_part1_permuted_(metrics, I, J, K,
                                i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper3way_part2_permuted_(metrics, I, J, K,
                                i_block, j_block, k_block, env) :
    GMMetrics_helper3way_part3_permuted_(metrics, I, J, K,
                                         i_block, j_block, k_block, env);

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

#ifdef GM_ASSERTIONS_ON

  int section_step = 0;

  if (gm_is_part3(i_block, j_block, k_block)) {
    section_step = 0;
  } else if (gm_is_part1(i_block, j_block, k_block)) {
    const int num_section_steps = GMEnv_num_section_steps(env, 1);
    const int section_num = (J * num_section_steps) / metrics->num_vector_local;
    section_step = section_num;
  } else {
    const int num_section_steps = GMEnv_num_section_steps(env, 2);
    const int section_num = (J * num_section_steps) / metrics->num_vector_local;
    section_step = section_num;
  }

  GMSectionInfo si_value_;
  GMSectionInfo* si = &si_value_;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const bool is_part3 = si->is_part3;
  const bool sax0 = si->section_axis == 0;
  const bool sax1 = si->section_axis == 1;
  //const bool sax2 = si->section_axis == 2;

  /* clang-format off */
  const int i = !is_part3 ?   I :
                     sax0 ?   J :
                     sax1 ?   I :
                  /* sax2 ?*/ K;
  const int j = !is_part3 ?   J :
                     sax0 ?   K :
                     sax1 ?   J :
                  /* sax2 ?*/ I;
  const int k = !is_part3 ?   K :
                     sax0 ?   I :
                     sax1 ?   K :
                  /* sax2 ?*/ J;
  /* clang-format on */

  GMSectionInfo_destroy(si, env);

  GMAssert(metrics->coords_global_from_index[index] % metrics->num_vector ==
           i + i_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector ==
           j + j_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] / metrics->num_vector) /
               metrics->num_vector ==
           k + k_block * (size_t)metrics->num_vector_local);
#endif

  return index;
}

/*---------------------------------------------------------------------------*/

typedef struct {
  bool is_initialized;
  int I;
  int K;
  size_t index;
} GMIndexCache;

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && J >= 0 && K >= 0);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == j_block &&
             GMEnv_proc_num_vector_i(env) != k_block));
  GMAssert(!(GMEnv_proc_num_vector_i(env) == k_block &&
             GMEnv_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  if (index_cache->is_initialized && K == index_cache->K) {
      const size_t index = index_cache->index + (I-index_cache->I);
      index_cache->index = index;
      index_cache->I = I;
      return index;
  }

  const size_t index = GMMetrics_index_from_coord_all2all_3_permuted(
    metrics, I, J, K, j_block, k_block, env);

  index_cache->I = I;
  index_cache->K = K;
  index_cache->index = index;
  index_cache->is_initialized = true;

  return index;
}
/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: value from (contig) index: basic---*/

static GMFloat3 GMMetrics_float3_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  GMAssert(metrics->data_S);

  return ((GMFloat3*)(metrics->data_S))[index];
}

/*---------------------------------------------------------------------------*/

static GMFloat3 GMMetrics_float3_C_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  GMAssert(metrics->data_C);

  return ((GMFloat3*)(metrics->data_C))[index];
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_from_index(GMMetrics* metrics,
                                                    size_t index,
                                                    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);

  return ((GMTally4x2*)(metrics->data))[index];
}

/*===========================================================================*/
/*---Accessors: value from (contig) index: derived---*/

static GMFloat GMMetrics_ccc_value_3(GMMetrics* metrics,
                                    const GMTally1 rijk,
                                    const GMTally1 si,
                                    const GMTally1 sj,
                                    const GMTally1 sk,
                                    const GMFloat recip_ci,
                                    const GMFloat recip_cj,
                                    const GMFloat recip_ck,
                                    const GMFloat recip_sumcijk,
                                    GMEnv* env) {
  GMAssert(metrics && env);

  const GMFloat one = 1;

  const GMFloat fi = (one / 2) * recip_ci * si;
  const GMFloat fj = (one / 2) * recip_cj * sj;
  const GMFloat fk = (one / 2) * recip_ck * sk;

  const GMFloat fijk = recip_sumcijk * rijk;

  /*---Do the following to make floating point arithmetic order-independent---*/

  GMFloat fmin = 0;
  GMFloat fmid = 0;
  GMFloat fmax = 0;

  if (fi > fj) {
    if (fi > fk) {
      fmax = fi;
      if (fj > fk) {
        fmid = fj;
        fmin = fk;
      } else {
        fmid = fk;
        fmin = fj;
      }
    } else {
      fmid = fi;
      fmax = fk;
      fmin = fj;
    }
  } else {
    if (fj > fk) {
      fmax = fj;
      if (fi > fk) {
        fmid = fi;
        fmin = fk;
      } else {
        fmid = fk;
        fmin = fi;
      }
    } else {
      fmid = fj;
      fmax = fk;
      fmin = fi;
    }
  }

  GMAssert(fmin <= fmid);
  GMAssert(fmid <= fmax);

  const GMFloat ccc_multiplier = one; // TBD
  const GMFloat ccc_param = GMEnv_ccc_param(env);

  /* clang-format off */
  const GMFloat result = ccc_multiplier * fijk * (one - ccc_param * fmin) *
                                                 (one - ccc_param * fmid) *
                                                 (one - ccc_param * fmax);
  /* clang-format on */

  return result;
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_ccc_get_from_index_3(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              int i2,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);
  GMAssert(i2 >= 0 && i2 < 2);

  const GMFloat one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally4x2 t42 = GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(t42, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_S_get_from_index(metrics, index, env);

  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  GMTally1 ci, cj, ck;
  if (env->sparse) {
    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);
  } else {
    ci = metrics->num_field_active;
    cj = metrics->num_field_active;
    ck = metrics->num_field_active;
  }

  /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;
  const GMTally1 sk = i2 == 0 ? (2 * ck - sk1) : sk1;

  const GMFloat recip_ci = env->sparse ? one / ci : recip_m;
  const GMFloat recip_cj = env->sparse ? one / cj : recip_m;
  const GMFloat recip_ck = env->sparse ? one / ck : recip_m;

  const GMFloat recip_sumcijk = env->sparse ?
    one / (GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
           GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
           GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
           GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1)) :
    (one / 8) * recip_m;

  return GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk, recip_ci,
                               recip_cj, recip_ck, recip_sumcijk, env);
}

/*===========================================================================*/
/*---Accessors: value from (local) coord: set: 3-way---*/

static void GMMetrics_float_set_3(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  int k,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_S_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMFloat3 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_S);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_C_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMFloat3 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_C);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally4x2_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMTally4x2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_3(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int j_block,
                                          int k_block,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_S_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_C_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally4x2_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMTally4x2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_3_permuted(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat value,
    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(
    metrics, I, J, K, j_block, k_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_S_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_C_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally4x2_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMTally4x2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K,
                                                               j_block, k_block,
                                                               env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env && index_cache);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_S_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat3 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env && index_cache);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_C_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat3 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env && index_cache);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally4x2_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMTally4x2 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env && index_cache);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
      metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*===========================================================================*/
/*---Accessors: value from (local) coord: get: 3-way---*/

static GMFloat GMMetrics_float_get_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(!GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(i < j && j < k);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int k,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0 && i < metrics->num_vector_local);
  GMAssert(j >= 0 && j < metrics->num_vector_local);
  GMAssert(k >= 0 && k < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3_permuted(GMMetrics* metrics,
                                                   int I,
                                                   int J,
                                                   int K,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMIndexCache* index_cache,
    GMEnv* env) {
  GMAssert(metrics && env && index_cache);
  GMAssert(GMEnv_all2all(env));
  GMAssert(I >= 0 && I < metrics->num_vector_local);
  GMAssert(J >= 0 && J < metrics->num_vector_local);
  GMAssert(K >= 0 && K < metrics->num_vector_local);
  GMAssert(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssert(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
      metrics, I, J, K, j_block, k_block, index_cache, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: 3-way---*/

static int GMMetrics_coord0_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t i64 = metrics->coords_global_from_index[index] %
                     metrics->num_vector;
  const int i = (int)i64;
  GMAssert((size_t)i == i64);

  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t j64 =
      (metrics->coords_global_from_index[index] / metrics->num_vector) %
      metrics->num_vector;
  const int j = (int)j64;
  GMAssert((size_t)j == j64);

  return j;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord2_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index >= 0 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_3);

  const size_t k64 = metrics->coords_global_from_index[index] /
                     (metrics->num_vector * (size_t)metrics->num_vector);
  const int k = (int)k64;
  GMAssert((size_t)k == k64);

  return k;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_gm_metrics_3way_hh_---*/

/*---------------------------------------------------------------------------*/
