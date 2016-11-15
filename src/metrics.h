/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _metrics_h_
#define _metrics_h_

#include <stddef.h>

#include "env.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Lightweight functions---*/

static _Bool gm_is_part1(int i_block, int j_block, int k_block) {
  return i_block == j_block && j_block == k_block;
}

/*---------------------------------------------------------------------------*/

static _Bool gm_is_part3(int i_block, int j_block, int k_block) {
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

/*===========================================================================*/
/*---GMSectionInfo---*/

typedef struct {
  _Bool is_part1;
  _Bool is_part2;
  _Bool is_part3;
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
  GMAssertAlways(si != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(i_block >= 0 && i_block < Env_num_block_vector(env));
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(num_vector_local >= 0);

  si->num_vector_local = num_vector_local;

  si->is_part1 = gm_is_part1(i_block, j_block, k_block);
  si->is_part3 = gm_is_part3(i_block, j_block, k_block);
  si->is_part2 = (!si->is_part1) && (!si->is_part3);

  const int part_num = si->is_part1 ? 1 :
                       si->is_part2 ? 2 : 3;

  const int num_section_steps = GMEnv_num_section_steps(env, part_num);
  GMAssertAlways(section_step>=0);
  GMAssertAlways(section_step<num_section_steps);

  si->section_axis =
    ! si->is_part3 ? 1 /*---j axis---*/ :
    gm_section_axis_part3(i_block, j_block, k_block);

  si->section_num = !si->is_part3 ? section_step :
                    gm_section_num_part3(i_block, j_block, k_block);

  /*---Define bounding box containing region to be computed---*/

  const int nvl6 = si->num_vector_local / 6;

  si->i_lb = si->section_axis == 0 && (si->is_part3 || num_section_steps == 6) ?
             si->section_num * nvl6 : 0;

  si->j_lb = si->section_axis == 1 && (si->is_part3 || num_section_steps == 6) ?
             si->section_num * nvl6 : 0;

  si->k_lb = si->section_axis == 2 && (si->is_part3 || num_section_steps == 6) ?
             si->section_num * nvl6 : 0;

  si->i_ub = si->section_axis == 0 && (si->is_part3 || num_section_steps == 6) ?
             (si->section_num + 1) * nvl6 : num_vector_local;

  si->j_ub = si->section_axis == 1 && (si->is_part3 || num_section_steps == 6) ?
             (si->section_num + 1) * nvl6 : num_vector_local;

  si->k_ub = si->section_axis == 2 && (si->is_part3 || num_section_steps == 6) ?
             (si->section_num + 1) * nvl6 : num_vector_local;

  si->J_lb = si->is_part3 || num_section_steps == 6 ?
             si->section_num * nvl6 : 0;

  si->J_ub = si->is_part3 || num_section_steps == 6 ?
             (si->section_num + 1) * nvl6 : num_vector_local;
}

#if 0

- I_max = is_part1 ? J : numvecl; // XXX can work same way if permuted or not
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
  GMAssertAlways(si != NULL);
  GMAssert(j >= 0 && j < si->num_vector_local);
  GMAssertAlways(env != NULL);

  return si->is_part3 ? si->k_lb : j + 1;
}

/*---------------------------------------------------------------------------*/

static int GMSectionInfo_i_max(
  GMSectionInfo* si,
  int j,
  GMEnv* env) {
  GMAssertAlways(si != NULL);
  GMAssert(j >= 0 && j < si->num_vector_local);
  GMAssertAlways(env != NULL);

  return si->is_part1 ? j : si->i_ub;
}

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  /*---Logical sizes---*/
  int num_field;
  int num_field_local;
  int num_vector;
  int num_vector_local;
  int nvl6;
  int pad1;
  size_t num_elts_local;
  /*---Helper values---*/
  size_t index_offset_0_;
  size_t index_offset_01_;
  //int block_num_offset_0_;
  //int block_num_offset_01_;
  size_t index_offset_section_pt1_[6];
  size_t index_offset_section_pt2_[6];
  _Bool section_num_valid_pt1_[6];
  _Bool section_num_valid_pt2_[6];
  size_t section_size_pt2[6];
  GMFloat m;
  GMFloat recip_m;
  /*---map of (contig) index to linearized Cartesian coords---*/
  size_t* coords_global_from_index;
  /*---Other---*/
  int data_type_id;
  int data_type_num_values;
  void* __restrict__ data;
  void* __restrict__ data_M;
} GMMetrics;

/*===========================================================================*/
/*---Null object---*/

GMMetrics GMMetrics_null(void);

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void GMMetrics_create(GMMetrics* metrics,
                      int data_type_id,
                      int num_field,
                      int num_vector_local,
                      GMEnv* env);

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void GMMetrics_destroy(GMMetrics* metrics, GMEnv* env);

/*===========================================================================*/
/*---Metrics checksum---*/

GMChecksum GMMetrics_checksum(GMMetrics* metrics, GMEnv* env);

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: indexing: (contig) index from coord, 2-way---*/

static size_t gm_triang_(int i) {
  return (i * (size_t)(i-1)) >> 1;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(Env_proc_num_repl(env) == 0);

  size_t index = gm_triang_(j) + i;
  GMAssert(i + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper2way_maindiag_block_(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  //return GMMetrics_index_from_coord_2(metrics, i, j, env)
  return gm_triang_(j) + i;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper2way_offdiag_block_(GMMetrics* metrics,
                                                  int i,
                                                  int j,
                                                  int j_block,
                                                  GMEnv* env) {
  const int num_block = Env_num_block_vector(env);

  const int num_proc_r = Env_num_proc_repl(env);

  const int i_block = Env_proc_num_vector_i(env);

  /* clang-format off */
  return metrics->index_offset_0_ +
      i + metrics->num_vector_local * (
      j + metrics->num_vector_local * (
      ((j_block - i_block - 1 + 2*num_block) % num_block) / num_proc_r ));
  /* clang-format on */
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_index_from_coord_all2all_2(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(i < j || j_block != Env_proc_num_vector_i(env));
//  GMAssert(Env_proc_num_repl(env) == 0 ||
//           j_block != Env_proc_num_vector(env));
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  const int i_block = Env_proc_num_vector_i(env);

  size_t index = j_block == i_block
           ? GMMetrics_helper2way_maindiag_block_(metrics, i, j, j_block, env)
           : GMMetrics_helper2way_offdiag_block_(metrics, i, j, j_block, env);

  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);

  const int nvl = metrics->num_vector_local;

  /* clang-format off */
  size_t index = i +
                 (k-j-1)*(size_t)j +
                 gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

  GMAssert(i + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           (metrics->coords_global_from_index[index] / metrics->num_vector) %
               metrics->num_vector);
  GMAssert(k + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] /
               (metrics->num_vector * metrics->num_vector));

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
  GMAssert(metrics->section_num_valid_pt1_[section_num]);

  const size_t elts_offset = metrics->index_offset_section_pt1_[section_num];

  /* clang-format off */
  const size_t index = elts_offset +
                       i +
                       (k-j-1)*(size_t)j +
                       gm_trap_size(j, nvl);
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

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
  GMAssert(metrics->section_num_valid_pt2_[section_num]);

  const size_t elts_offset = metrics->index_offset_section_pt2_[section_num];

  const int num_block = Env_num_block_vector(env);
  const int block_num_pt2 = gm_mod1_(j_block - i_block, num_block) - 1;
  const int num_proc_r = Env_num_proc_repl(env);
  const int blocks_offset = block_num_pt2 / num_proc_r;

  const size_t section_size = metrics->section_size_pt2[section_num];

  /* clang-format off */
  size_t index = elts_offset +
                 i + nvl*(
                 (k-j-1) +
                 gm_triang_size(j, nvl) + section_size*(
                 blocks_offset
                 ));
 /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

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
  const int nvl6 = metrics->nvl6;

  const int section_axis = gm_section_axis_part3(i_block, j_block, k_block);
  const int section_num = gm_section_num_part3(i_block, j_block, k_block);

  const size_t elts_offset = metrics->index_offset_01_;

  const int num_block = Env_num_block_vector(env);
  const int j_i_block_delta = gm_mod1_(j_block - i_block, num_block);
  const int k_i_block_delta = gm_mod1_(k_block - i_block, num_block);
  const int block_num_pt3 =
    ((num_block-2) * (k_i_block_delta - 1)) +
    (j_i_block_delta - 1 - (j_i_block_delta > k_i_block_delta));
  const int num_proc_r = Env_num_proc_repl(env);
  const int blocks_offset = block_num_pt3 / num_proc_r;

  /* clang-format off */
  size_t index = elts_offset +
        i - ( section_axis == 0 ? section_num * nvl6 : 0 ) +
            ( section_axis == 0 ? nvl6 : nvl ) * (
        k - ( section_axis == 2 ? section_num * nvl6 : 0 ) +
            ( section_axis == 2 ? nvl6 : nvl ) * (
        j - ( section_axis == 1 ? section_num * nvl6 : 0 ) +
            ( section_axis == 1 ? nvl6 : nvl ) * (
        blocks_offset
        )));
  /* clang-format on */

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  GMAssert(!(Env_proc_num_vector_i(env) == j_block &&
             Env_proc_num_vector_i(env) != k_block));
  GMAssert(!(Env_proc_num_vector_i(env) == k_block &&
             Env_proc_num_vector_i(env) != j_block));
  /*---WARNING: these conditions are not exhaustive---*/

  const int i_block = Env_proc_num_vector_i(env);

  size_t index = j_block == i_block && k_block == i_block ?
    GMMetrics_helper3way_part1_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper3way_part2_(metrics, i, j, k,
                                i_block, j_block, k_block, env) :
    GMMetrics_helper3way_part3_(metrics, i, j, k,
                                i_block, j_block, k_block, env);

  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);

  GMAssert(metrics->coords_global_from_index[index] %
             (metrics->num_vector_local * (size_t)Env_num_block_vector(env)) ==
           i + i_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)Env_num_block_vector(env))) %
               (metrics->num_vector_local * Env_num_block_vector(env)) ==
           j + j_block * (size_t)metrics->num_vector_local);

  GMAssert((metrics->coords_global_from_index[index] /
            (metrics->num_vector_local * (size_t)Env_num_block_vector(env))) /
               (metrics->num_vector_local * Env_num_block_vector(env)) ==
           k + k_block * (size_t)metrics->num_vector_local);

  return index;
}

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: value from (contig) index: basic---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              int index,
                                              GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(metrics->data))[index];
}

/*---------------------------------------------------------------------------*/

static GMFloat2 GMMetrics_float2_M_get_from_index(GMMetrics* metrics,
                                                  int index,
                                                  GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  return ((GMFloat2*)(metrics->data_M))[index];
}

/*---------------------------------------------------------------------------*/

static GMFloat3 GMMetrics_float3_M_get_from_index(GMMetrics* metrics,
                                                  int index,
                                                  GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);

  return ((GMFloat3*)(metrics->data_M))[index];
}

/*---------------------------------------------------------------------------*/

static GMTally2x2 GMMetrics_tally2x2_get_from_index(GMMetrics* metrics,
                                                    int index,
                                                    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  return ((GMTally2x2*)(metrics->data))[index];
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_from_index(GMMetrics* metrics,
                                                    int index,
                                                    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);

  return ((GMTally4x2*)(metrics->data))[index];
}

/*===========================================================================*/
/*---Accessors: value from (contig) index: derived---*/


static GMFloat GMMetrics_czekanowski_get_from_index(GMMetrics* metrics,
                                                    int index,
                                                    GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);

  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_ccc_get_from_index_2(GMMetrics* metrics,
                                              int index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMTally2x2 tally2x2 =
      GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(tally2x2, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_M_get_from_index(metrics, index, env);

  GMTally1 si_1;
  GMTally1 sj_1;
  GMFloat2_decode(&si_1, &sj_1, si1_sj1);

  /*---Get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * metrics->num_field - si_1) : si_1;
  const GMTally1 sj = i1 == 0 ? (2 * metrics->num_field - sj_1) : sj_1;

  /*---Do the following to make floating point arithmetic order-independent---*/
  const GMTally1 smin = si < sj ? si : sj;
  const GMTally1 smax = si < sj ? sj : si;

  //---TODO: optimize
  const GMFloat one = 1;
  const GMFloat m = metrics->num_field;
  const GMFloat recip_m = metrics->recip_m;
  const GMFloat front_multiplier = 9 * one / 2;

  /*---Arrange so as to guarantee each factor nonnegative---*/
  /* clang-format off */
  const GMFloat result = (front_multiplier / 4) * recip_m * rij *
                         (3 * m - smin) * (one/3) * recip_m *
                         (3 * m - smax) * (one/3) * recip_m;
  /* clang-format on */

  return result;
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_ccc_get_from_index_3(GMMetrics* metrics,
                                              int index,
                                              int i0,
                                              int i1,
                                              int i2,
                                              GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index >= 0);
  GMAssert((size_t)index < metrics->num_elts_local);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);
  GMAssert(i2 >= 0 && i2 < 2);

  const GMTally4x2 tally4x2 =
      GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(tally4x2, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_M_get_from_index(metrics, index, env);

  GMTally1 si_1;
  GMTally1 sj_1;
  GMTally1 sk_1;
  GMFloat3_decode(&si_1, &sj_1, &sk_1, si1_sj1_sk1);

  /*---Get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * metrics->num_field - si_1) : si_1;
  const GMTally1 sj = i1 == 0 ? (2 * metrics->num_field - sj_1) : sj_1;
  const GMTally1 sk = i2 == 0 ? (2 * metrics->num_field - sk_1) : sk_1;

  /*---Do the following to make floating point arithmetic order-independent---*/

  GMTally1 smin = 0;
  GMTally1 smid = 0;
  GMTally1 smax = 0;

  if (si <= sj && si <= sk) {
    smin = si;
    if (sj < sk) {
      smid = sj;
      smax = sk;
    } else /*---sk <= sj---*/ {
      smid = sk;
      smax = sj;
    }
  } else if (sj <= si && sj <= sk) {
    smin = sj;
    if (si < sk) {
      smid = si;
      smax = sk;
    } else /*---sk <= si---*/ {
      smid = sk;
      smax = si;
    }
  } else /*---sk <= si && sk <= sj ...---*/ {
    smin = sk;
    if (si < sj) {
      smid = si;
      smax = sj;
    } else /*---sj <= si---*/ {
      smid = sj;
      smax = si;
    }
  }
  GMAssert(smin <= smid);
  GMAssert(smid <= smax);

  //---TODO: optimize
  const GMFloat one = 1;
  const GMFloat m = metrics->num_field;
  const GMFloat recip_m = metrics->recip_m;
  const GMFloat front_multiplier_TBD = 2 * one / 2;

  /*---Arrange so as to guarantee each factor nonnegative---*/
  /* clang-format off */
  const GMFloat result = (front_multiplier_TBD / 8) * recip_m * rijk *
                         (3 * m - smin) * (one/3) * recip_m *
                         (3 * m - smid) * (one/3) * recip_m *
                         (3 * m - smax) * (one/3) * recip_m;
  /* clang-format on */

  return result;
}

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: value from (local) coord: set: 2-way---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float2_M_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMFloat2 value,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat2*)(metrics->data_M))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally2x2_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMTally2x2 value,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_2(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int j_block,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(i < j || j_block != Env_proc_num_vector_i(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float2_M_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMFloat2 value,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(i < j || j_block != Env_proc_num_vector_i(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat2*)(metrics->data_M))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally2x2_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMTally2x2 value,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(i < j || j_block != Env_proc_num_vector_i(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

/*===========================================================================*/
/*---Accessors: value from (local) coord: set: 3-way---*/

static void GMMetrics_float_set_3(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  int k,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_M_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMFloat3 value,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat3*)(metrics->data_M))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_tally4x2_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMTally4x2 value,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float_set_all2all_3(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int j_block,
                                          int k_block,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

/*---------------------------------------------------------------------------*/

static void GMMetrics_float3_M_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_M))[index] = value;
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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: value from (local) coord: get: 2-way---*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(!Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(i < j || j_block != Env_proc_num_vector_i(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*===========================================================================*/
/*---Accessors: value from (local) coord: get: 3-way---*/

static GMFloat GMMetrics_float_get_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(env);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMFloat GMMetrics_float_get_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

/*---------------------------------------------------------------------------*/

static GMTally4x2 GMMetrics_tally4x2_get_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(k >= 0);
  GMAssert(k < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(j < k);
  GMAssert(env);
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);

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
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(Env_all2all(env));
  GMAssert(i >= 0);
  GMAssert(j >= 0);
  GMAssert(k >= 0);
  GMAssert(j_block >= 0);
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block >= 0);
  GMAssert(k_block < Env_num_block_vector(env));
  GMAssert(Env_data_type_metrics(env) == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: 2-way---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_2);

  const int j = metrics->coords_global_from_index[index] / metrics->num_vector;
  return j;
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: 3-way---*/

static int GMMetrics_coord0_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);

  const int i = metrics->coords_global_from_index[index] % metrics->num_vector;
  return i;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord1_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);

  const int j =
      (metrics->coords_global_from_index[index] / metrics->num_vector) %
      metrics->num_vector;
  return j;
}

/*---------------------------------------------------------------------------*/

static int GMMetrics_coord2_global_from_index_3(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);

  const int k = metrics->coords_global_from_index[index] /
                (metrics->num_vector * metrics->num_vector);
  return k;
}

/*===========================================================================*/
/*---Accessors: indexing: global coord from (contig) index: generic---*/

static int GMMetrics_coord_global_from_index(GMMetrics* metrics,
                                             size_t index,
                                             int coord_num,
                                             GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(index >= 0);
  GMAssert(index < metrics->num_elts_local);
  GMAssert(coord_num >= 0);
  GMAssert(coord_num < Env_num_way(env));

  int result = 0;

  GMAssert(Env_num_way(env) <= GM_NUM_NUM_WAY + 1
               ? "this num_way currently not supported"
               : 0);

  switch (Env_num_way(env) + 4 * coord_num) {
    case 2 + 4 * 0: /* 2-way, coord 0 */
      result = GMMetrics_coord0_global_from_index_2(metrics, index, env);
      break;
    case 2 + 4 * 1: /* 2-way, coord 1 */
      result = GMMetrics_coord1_global_from_index_2(metrics, index, env);
      break;
    case 3 + 4 * 0: /* 3-way, coord 0 */
      result = GMMetrics_coord0_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 1: /* 3-way, coord 1 */
      result = GMMetrics_coord1_global_from_index_3(metrics, index, env);
      break;
    case 3 + 4 * 2: /* 3-way, coord 2 */
      result = GMMetrics_coord2_global_from_index_3(metrics, index, env);
      break;
    default:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/

  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_metrics_h_---*/

/*---------------------------------------------------------------------------*/
