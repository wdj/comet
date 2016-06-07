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
/*---Struct declaration---*/

typedef struct {
  /*---Logical sizes---*/
  int num_field;
  int num_field_local;
  int num_vector;
  int num_vector_local;
  int pad1;
  size_t num_elts_local;
  /*---Helper values---*/
  size_t index_offset_0_;
  size_t index_offset_01_;
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
/*---Helper functions for 3-way all2all case---*/

static int gm_metrics_3way_section_axis(GMMetrics* metrics,
                                        int i_block,
                                        int j_block,
                                        int k_block,
                                        GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i_block >= 0);
  GMAssert(j_block >= 0);
  GMAssert(k_block >= 0);
  GMAssert(i_block < Env_num_block_vector(env));
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block < Env_num_block_vector(env));
  // GMAssert( i_block != j_block );
  // GMAssert( i_block != k_block );
  // GMAssert( k_block != j_block );

  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/

  /* clang-format off */
  return i_block < j_block && i_block < k_block ? 0 : /*---i axis---*/
         j_block < i_block && j_block < k_block ? 1 : /*---j axis---*/
                                                  2;  /*---k axis---*/
  /* clang-format on */
}

/*---------------------------------------------------------------------------*/

static int gm_metrics_3way_section_num(GMMetrics* metrics,
                                       int i_block,
                                       int j_block,
                                       int k_block,
                                       GMEnv* env) {
  GMAssert(metrics != NULL);
  GMAssert(env != NULL);
  GMAssert(Env_num_way(env) == GM_NUM_WAY_3);
  GMAssert(i_block >= 0);
  GMAssert(j_block >= 0);
  GMAssert(k_block >= 0);
  GMAssert(i_block < Env_num_block_vector(env));
  GMAssert(j_block < Env_num_block_vector(env));
  GMAssert(k_block < Env_num_block_vector(env));
  // GMAssert( i_block != j_block );
  // GMAssert( i_block != k_block );
  // GMAssert( k_block != j_block );

  /*---NOTE: this could possibly be implemented somewhat more efficiently---*/

  /* clang-format off */
  return i_block < k_block && k_block < j_block ?    0 :
         i_block < j_block && j_block < k_block ?    1 :
         j_block < i_block && i_block < k_block ?    2 :
         j_block < k_block && k_block < i_block ?    3 :
         k_block < j_block && j_block < i_block ?    4 :
      /* k_block < i_block && i_block < j_block ? */ 5;
  /* clang-format on */
}

/*===========================================================================*/
/*===========================================================================*/
/*---Accessors: indexing: (contig) index from coord, 2-way---*/

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
  GMAssert(Env_proc_num_vector_j(env) == 0);

  size_t index = ((j * (size_t)(j - 1)) >> 1) + i;
  GMAssert(i + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local * (size_t)Env_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper1_maindiag_block_(GMMetrics* metrics,
                                                int i,
                                                int j,
                                                int j_block,
                                                GMEnv* env) {
  //return GMMetrics_index_from_coord_2(metrics, i, j, env)
  return ((j * (size_t)(j - 1)) >> 1) + i;
}

/*---------------------------------------------------------------------------*/

static size_t GMMetrics_helper1_offdiag_block_(GMMetrics* metrics,
                                               int i,
                                               int j,
                                               int j_block,
                                               GMEnv* env) {
  const int i_block = Env_proc_num_vector_i(env);

  /* clang-format off */
  return metrics->index_offset_0_ +
      i + metrics->num_vector_local * (
      j + metrics->num_vector_local * (
      ((j_block - i_block - 1 + 2*Env_num_block_vector(env)) %
                                  Env_num_block_vector(env)
      ) / Env_num_proc_vector_j(env) ));
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
//  GMAssert(Env_proc_num_vector_j(env) == 0 ||
//           j_block != Env_proc_num_vector(env));
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  const int i_block = Env_proc_num_vector_i(env);

  size_t index = j_block == i_block
               ? GMMetrics_helper1_maindiag_block_(metrics, i, j, j_block, env)
               : GMMetrics_helper1_offdiag_block_(metrics, i, j, j_block, env);

  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
}

/*===========================================================================*/
/*---Accessors: indexing: (contig) index from coord, 3-way---*/

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

  size_t index = (k * (size_t)(k - 1) * (size_t)(k - 2)) / 6 +
                 (j * (size_t)(j - 1)) / 2 + i;
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

static size_t GMMetrics_helper1_part1(GMMetrics* metrics,
                                      int i,
                                      int j,
                                      int k,
                                      int i_block,
                                      int j_block,
                                      int k_block,
                                      GMEnv* env) {
  GMAssert(Env_proc_num_vector_j(env) == 0);
  return GMMetrics_index_from_coord_3(metrics, i, j, k, env);
}

/*---------------------------------------------------------------------------*/

//CHANGE
static size_t GMMetrics_helper1_part2(GMMetrics* metrics,
                                      int i,
                                      int j,
                                      int k,
                                      int i_block,
                                      int j_block,
                                      int k_block,
                                      GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  /* clang-format off */
  return metrics->index_offset_0_ +
         i + nvl * (
         ((k * (size_t)(k - 1)) >> 1) + j + 
             ((nvl * (size_t)(nvl - 1)) >> 1) * (
         j_block - ( j_block > i_block ) ));
  /* clang-format on */
}

/*---------------------------------------------------------------------------*/

//CHANGE
static size_t GMMetrics_helper1_part3(GMMetrics* metrics,
                                      int i,
                                      int j,
                                      int k,
                                      int i_block,
                                      int j_block,
                                      int k_block,
                                      GMEnv* env) {
  const int nvl = metrics->num_vector_local;

  const int section_axis =
      gm_metrics_3way_section_axis(metrics, i_block, j_block, k_block, env);
  const int section_num =
      gm_metrics_3way_section_num(metrics, i_block, j_block, k_block, env);

  /* clang-format off */
  return metrics->index_offset_01_ +
         i - ( section_axis == 0 ? section_num * nvl / 6 : 0 ) +
             ( section_axis == 0 ? nvl / 6 : nvl ) * (
         j - ( section_axis == 1 ? section_num * nvl / 6 : 0 ) +
             ( section_axis == 1 ? nvl / 6 : nvl ) * (
         k - ( section_axis == 2 ? section_num * nvl / 6 : 0 ) +
             ( section_axis == 2 ? nvl / 6 : nvl ) * (
         j_block - ( j_block > i_block ) - ( j_block > k_block ) +
                                       (Env_num_block_vector(env)-2) * (
         k_block - ( k_block > i_block ) ))));
  /* clang-format on */
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
    GMMetrics_helper1_part1(metrics, i, j, k, i_block, j_block, k_block, env) :
                 j_block == k_block ?
    GMMetrics_helper1_part2(metrics, i, j, k, i_block, j_block, k_block, env) :
    GMMetrics_helper1_part3(metrics, i, j, k, i_block, j_block, k_block, env);

  GMAssert(index >= 0 && index < metrics->num_elts_local);

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
