//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_metrics_2way_hh_
#define _gm_metrics_2way_hh_

//=============================================================================
/*---Helper functions for 2-way case---*/

static int gm_max_computed_blocks_per_row(GMEnv* env) {
  GMAssert(env);
  /*---Max number of blocks of any block row computed on all phases---*/
  const int num_block = GMEnv_num_block_vector(env);
  return 1 + num_block / 2;
}

//-----------------------------------------------------------------------------

static int gm_diag_computed_min(GMEnv* env) {
  GMAssert(env);
  const int max_rectangle_width = gm_max_computed_blocks_per_row(env);
  return (max_rectangle_width*env->phase_num) / env->num_phase;
}

//-----------------------------------------------------------------------------

static int gm_diag_computed_max(GMEnv* env) {
  GMAssert(env);
  const int max_rectangle_width = gm_max_computed_blocks_per_row(env);
  return (max_rectangle_width*(env->phase_num+1)) / env->num_phase;
}

//-----------------------------------------------------------------------------

static int gm_diag_computed_this_row_max(GMEnv* env) {
  GMAssert(env);

  const int num_block = GMEnv_num_block_vector(env);
  const int i_block = GMEnv_proc_num_vector_i(env);

  const bool is_row_short_by_1 = num_block % 2 == 0 && 2*i_block >= num_block;
  const bool is_last_phase = env->phase_num == env->num_phase - 1;

  const int diag_max = gm_diag_computed_max(env);

  const int n = is_last_phase && is_row_short_by_1 ? diag_max - 1 : diag_max;
  GMAssert(n >= 0);
  GMAssert(n <= num_block);
  return n;
}

//-----------------------------------------------------------------------------

static int gm_computed_blocks_this_row(GMEnv* env) {
  GMAssert(env);
  const int n = gm_diag_computed_this_row_max(env) - gm_diag_computed_min(env);
  GMAssert(n >= 0);
  GMAssert(n <= GMEnv_num_block_vector(env));
  return n;
}

//=============================================================================
/*---Accessors: indexing: (contig) index from coord, 2-way---*/

static size_t gm_triang_(int i) {
  return (i * (size_t)(i-1)) >> 1;
}

//-----------------------------------------------------------------------------
/*---Helper: is this (section_)block_num to be processed by this proc_r---*/

static bool gm_proc_r_active(int section_block_num, const GMEnv* const env) {
  GMAssert(env);
  GMAssert(section_block_num >= 0);
  return section_block_num % GMEnv_num_proc_repl(env)
         == GMEnv_proc_num_repl(env);
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_2(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_proc_num_repl(env) == 0);

  size_t index = gm_triang_(j) + i;
  GMAssert(i + metrics->num_vector_local *
               (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] % metrics->num_vector);
  GMAssert(j + metrics->num_vector_local *
               (size_t)GMEnv_proc_num_vector_i(env) ==
           metrics->coords_global_from_index[index] / metrics->num_vector);
  return index;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_maindiag_block_(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  GMAssert(j_block == GMEnv_proc_num_vector_i(env));

  return gm_triang_(j) + i;
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_helper2way_offdiag_block_(GMMetrics* metrics,
                                                  int i,
                                                  int j,
                                                  int j_block,
                                                  GMEnv* env) {
  GMAssert(j_block != GMEnv_proc_num_vector_i(env));

  const int num_block = GMEnv_num_block_vector(env);

  const int num_proc_r = GMEnv_num_proc_repl(env);

  const int block_min = metrics->block_min;

  /* clang-format off */
  return metrics->index_offset_0_ +
      i + metrics->num_vector_local * (
      j + metrics->num_vector_local * (
      ((j_block - block_min + num_block) % num_block) / num_proc_r ));
  /* clang-format on */
}

//-----------------------------------------------------------------------------

static size_t GMMetrics_index_from_coord_all2all_2(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int j_block,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
//  GMAssert(GMEnv_proc_num_repl(env) == 0 ||
//           j_block != GMEnv_proc_num_vector(env));
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  const int i_block = GMEnv_proc_num_vector_i(env);

  size_t index = j_block == i_block
           ? GMMetrics_helper2way_maindiag_block_(metrics, i, j, j_block, env)
           : GMMetrics_helper2way_offdiag_block_(metrics, i, j, j_block, env);

  GMAssert(index >= 0 && index < metrics->num_elts_local);
  return index;
}

//=============================================================================
//=============================================================================
/*---Accessors: value from (contig) index: basic---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              size_t index,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(metrics->data))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);

  return ((GMFloat2*)(metrics->data_S))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_C_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);

  return ((GMFloat2*)(metrics->data_C))[index];
}

//-----------------------------------------------------------------------------

static GMTally2x2 GMMetrics_tally2x2_get_from_index(GMMetrics* metrics,
                                                    size_t index,
                                                    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  return ((GMTally2x2*)(metrics->data))[index];
}

//=============================================================================
/*---Accessors: value from (contig) index: derived---*/


static GMFloat GMMetrics_czek_get_from_index(GMMetrics* metrics,
                                             size_t index,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);

  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_value_2(GMMetrics* metrics,
                                    const GMTally1 rij,
                                    const GMTally1 si,
                                    const GMTally1 sj,
                                    const GMFloat recip_ci,
                                    const GMFloat recip_cj,
                                    const GMFloat recip_sumcij,
                                    GMEnv* env) {
  GMAssert(metrics && env);

  const GMFloat f_one = 1;

  const GMFloat fi = (f_one / 2) * recip_ci * si;
  const GMFloat fj = (f_one / 2) * recip_cj * sj;

  const GMFloat fij = recip_sumcij * rij;

  /*---Do the following to make floating point arithmetic order-independent---*/
  const GMFloat fmin = fi < fj ? fi : fj;
  const GMFloat fmax = fi < fj ? fj : fi;

  const GMFloat ccc_multiplier = GMEnv_ccc_multiplier(env);
  const GMFloat ccc_param = GMEnv_ccc_param(env);

  /* clang-format off */
  const GMFloat result = ccc_multiplier * fij * (f_one - ccc_param * fmin) *
                                                (f_one - ccc_param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------

#ifndef XXX

static GMFloat GMMetrics_ccc_get_from_index_2(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMFloat f_one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  if (env->sparse) {

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
    if (ci == 0 || cj == 0 || cij == 0) {
      return 0;
    }
    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;

    const GMFloat f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const GMFloat f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const GMFloat f_cij = (GMFloat) cij;
    const GMFloat recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;

    const GMFloat recip_ci = f_cj * f_cij * recip_cicjcij;
    const GMFloat recip_cj = f_ci * f_cij * recip_cicjcij;

    const GMFloat recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

    return GMMetrics_ccc_value_2(metrics, rij, si, sj,
                                 recip_ci, recip_cj, recip_sumcij, env);
  } /*---if sparse---*/

  GMAssert(metrics->num_field_active > 0);

  /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * metrics->num_field_active - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * metrics->num_field_active - sj1) : sj1;

  const GMFloat recip_sumcij = (f_one / 4) * recip_m;

  return GMMetrics_ccc_value_2(metrics, rij, si, sj,
                               recip_m, recip_m, recip_sumcij, env);
}

#else
//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_get_from_index_2(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMFloat one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMTally1 ci, cj;
  if (env->sparse) {
    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
    if (ci == 0 || cj == 0 || cij == 0) {
      return 0;
    }
  } else {
    ci = metrics->num_field_active;
    cj = metrics->num_field_active;
    GMAssert(metrics->num_field_active > 0);
  }

  /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;

  const GMFloat recip_ci = env->sparse ? one / ci : recip_m;
  const GMFloat recip_cj = env->sparse ? one / cj : recip_m;

  const GMFloat recip_sumcij = env->sparse ?
    one / (GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
           GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1)) :
    (one / 4) * recip_m;

  return GMMetrics_ccc_value_2(metrics, rij, si, sj,
                               recip_ci, recip_cj, recip_sumcij, env);
}

#endif
//=============================================================================
//=============================================================================
/*---Accessors: value from (local) coord: set: 2-way---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMFloat2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat2*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMFloat2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat2*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMTally2x2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_2(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int j_block,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMFloat2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat2*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMFloat2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat2*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMTally2x2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

//=============================================================================
/*---Accessors: value from (local) coord: get: 2-way---*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//=============================================================================
/*---Accessors: indexing: global coord from (contig) index: 2-way---*/

static int GMMetrics_coord0_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);


  const size_t i64 = metrics->coords_global_from_index[index] %
                     metrics->num_vector;
  const int i = (int)i64;
  GMAssert((size_t)i == i64);

  return i;
}

//-----------------------------------------------------------------------------

static int GMMetrics_coord1_global_from_index_2(GMMetrics* metrics,
                                                size_t index,
                                                GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);

  const size_t j64 = metrics->coords_global_from_index[index] /
                     metrics->num_vector;
  const int j = (int)j64;
  GMAssert((size_t)j == j64);

  return j;
}

//=============================================================================

#endif /*---_gm_metrics_2way_hh_---*/

//-----------------------------------------------------------------------------
