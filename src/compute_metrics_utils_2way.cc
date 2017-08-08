//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_utils_2way.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, for 2-way case.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "env.hh"
#include "mirrored_buf.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "linalg.hh"
#include "compute_metrics_utils_2way.hh"

#ifdef __cplusplus
extern "C" {
#endif

//=============================================================================
/*---Start calculation of numerators, 2-way Czekanowski---*/

void gm_compute_czek_numerators_2way_start_(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_left_buf,
    GMMirroredBuf* vectors_right_buf,
    GMMirroredBuf* metrics_buf,
    int j_block,
    bool do_compute_triang_only,
    GMEnv* env) {
  GMAssertAlways(vectors_left && vectors_right && metrics && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU && GMEnv_all2all(env)) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for REF case not supported"
                      : 0);

    /*---Perform pseudo GEMM---*/

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int fl = 0; fl < vectors_left->num_field_local; ++fl) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, fl, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, fl, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Perform pseudo GEMM---*/

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int fl = 0; fl < vectors_left->num_field_local; ++fl) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, fl, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, fl, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_2(metrics, i, j, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else /* if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_linalg_set_matrix_zero_start(metrics_buf, env);

    /*---Perform pseudo GEMM---*/

    gm_linalg_gemm_start(vectors_left->num_vector_local,
                         vectors_left->num_vector_local,
                         vectors_left->num_field_local,
                         vectors_left_buf->d, vectors_left->num_field_local,
                         vectors_right_buf->d, vectors_left->num_field_local,
                         metrics_buf->d, vectors_left->num_vector_local, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Start calculation of numerators, 2-way CCC---*/

void gm_compute_ccc_numerators_2way_start_(GMVectors* vectors_left,
                                           GMVectors* vectors_right,
                                           GMMetrics* metrics,
                                           GMMirroredBuf* vectors_left_buf,
                                           GMMirroredBuf* vectors_right_buf,
                                           GMMirroredBuf* metrics_buf,
                                           int j_block,
                                           bool do_compute_triang_only,
                                           GMEnv* env) {
  GMAssertAlways(vectors_left && vectors_right && metrics && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_REF) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for REF case not supported"
                      : 0);

    /*---Perform pseudo GEMM---*/

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        for (int fl = 0; fl < vectors_left->num_field_local; ++fl) {
          const GMBits2 vi = GMVectors_bits2_get(vectors_left, fl, i, env);
          const GMBits2 vj = GMVectors_bits2_get(vectors_right, fl, j, env);
          const bool unknown_i = env->sparse ? vi == GM_2BIT_UNKNOWN
                                             : false;
          const bool unknown_j = env->sparse ? vj == GM_2BIT_UNKNOWN
                                             : false;

          if ((!unknown_i) && (!unknown_j)) {

            /* clang-format off */
            const int r00 = ( ( !(vi & 1) ) && ( !(vj & 1) ) ) +
                            ( ( !(vi & 1) ) && ( !(vj & 2) ) ) +
                            ( ( !(vi & 2) ) && ( !(vj & 1) ) ) +
                            ( ( !(vi & 2) ) && ( !(vj & 2) ) );
            const int r01 = ( ( !(vi & 1) ) && (  (vj & 1) ) ) +
                            ( ( !(vi & 1) ) && (  (vj & 2) ) ) +
                            ( ( !(vi & 2) ) && (  (vj & 1) ) ) +
                            ( ( !(vi & 2) ) && (  (vj & 2) ) );
            const int r10 = ( (  (vi & 1) ) && ( !(vj & 1) ) ) +
                            ( (  (vi & 1) ) && ( !(vj & 2) ) ) +
                            ( (  (vi & 2) ) && ( !(vj & 1) ) ) +
                            ( (  (vi & 2) ) && ( !(vj & 2) ) );
            const int r11 = ( (  (vi & 1) ) && (  (vj & 1) ) ) +
                            ( (  (vi & 1) ) && (  (vj & 2) ) ) +
                            ( (  (vi & 2) ) && (  (vj & 1) ) ) +
                            ( (  (vi & 2) ) && (  (vj & 2) ) );
            /* clang-format on */

            /*---NOTE: "since the sum of all 4 of these relative
                 co-occurences is 4, we really only need to compute 3 of them.
                 Then the last one is just 4 minus the rest." */

#if DOUG_WAY
//TODO: work on this as a possibly faster way.
            const int vi1 = (vi & 3) != 0;
            const int vi0 = ((~vi) & 3) != 0;
            const int vj1 = (vj & 3) != 0;
            const int vj0 = ((~vj) & 3) != 0;

            const int a11 = vi1 & vj1;

            const int r11 = a11 +
#endif
          /*---Accumulate---*/

            sum.data[0] += GMTally1_encode(r00, r01);
            sum.data[1] += GMTally1_encode(r10, r11);
          } /*---if !unknown---*/
        } /*---for fl---*/
        if (GMEnv_all2all(env)) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_CPU) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Perform pseudo GEMM---*/

    /*---Precompute masks for final (incomplete) packedval_field -
         can be 1 to 64 inclusive---*/

    /* clang-format off */

    const int nfl = vectors_left->num_field_local;
    const int num_field_active_local =
      GMEnv_proc_num_field(env) == GMEnv_num_proc_field(env) - 1
      ? nfl - (vectors_left->num_field - vectors_left->num_field_active) : nfl;
    const int num_packedval_field_active_local =
      (num_field_active_local + 64 - 1) / 64;

    const int num_seminibbles_edge = 1 + (num_field_active_local-1) % 64;

    const GMUInt64 nobits = 0;
    const GMUInt64 allbits = ~nobits;

    const GMUInt64 edgemask0 = num_seminibbles_edge >= 32 ?
                               allbits :
                               allbits >> (64 - 2*num_seminibbles_edge);

    const GMUInt64 edgemask1 = num_seminibbles_edge <= 32 ?
                               nobits :
                               num_seminibbles_edge == 64 ?
                               allbits :
                               allbits >> (128 - 2*num_seminibbles_edge);

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        const int pvfl_edge = num_packedval_field_active_local - 1;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {
          /*---Get masks for active seminibbles in each word---*/

          const GMUInt64 activebits0 = pvfl < pvfl_edge ? allbits :
                                       pvfl == pvfl_edge ? edgemask0 : nobits;
          const GMUInt64 activebits1 = pvfl < pvfl_edge ? allbits :
                                       pvfl == pvfl_edge ? edgemask1 : nobits;

          /*---Extract input values to process---*/

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const GMUInt64 vi0 = vi.data[0];
          const GMUInt64 vi1 = vi.data[1];
          const GMUInt64 vj0 = vj.data[0];
          const GMUInt64 vj1 = vj.data[1];

          /*---Compute masks---*/

          const GMUInt64 oddbits = 0x5555555555555555;

          const GMUInt64 vi0mask = activebits0 &
                       (env->sparse ? (vi0 | ~(vi0 >> 1)) & oddbits : oddbits);

          const GMUInt64 vi1mask = activebits1 &
                       (env->sparse ? (vi1 | ~(vi1 >> 1)) & oddbits : oddbits);

          const GMUInt64 vj0mask = activebits0 &
                       (env->sparse ? (vj0 | ~(vj0 >> 1)) & oddbits : oddbits);

          const GMUInt64 vj1mask = activebits1 &
                       (env->sparse ? (vj1 | ~(vj1 >> 1)) & oddbits : oddbits);

          const GMUInt64 v0mask = vi0mask & vj0mask;
          const GMUInt64 v1mask = vi1mask & vj1mask;

          /*---Get even/odd bits for each seminibble, masked to active---*/

          const GMUInt64 vi0_0 =  vi0       & v0mask;
          const GMUInt64 vi0_1 = (vi0 >> 1) & v0mask;
          const GMUInt64 vi1_0 =  vi1       & v1mask;
          const GMUInt64 vi1_1 = (vi1 >> 1) & v1mask;
          const GMUInt64 vj0_0 =  vj0       & v0mask;
          const GMUInt64 vj0_1 = (vj0 >> 1) & v0mask;
          const GMUInt64 vj1_0 =  vj1       & v1mask;
          const GMUInt64 vj1_1 = (vj1 >> 1) & v1mask;

          /*---Get complements of even/odd bits for each seminibble; mask---*/

          const GMUInt64 nvi0_0 = ~ vi0       & v0mask;
          const GMUInt64 nvi0_1 = ~(vi0 >> 1) & v0mask;
          const GMUInt64 nvi1_0 = ~ vi1       & v1mask;
          const GMUInt64 nvi1_1 = ~(vi1 >> 1) & v1mask;
          const GMUInt64 nvj0_0 = ~ vj0       & v0mask;
          const GMUInt64 nvj0_1 = ~(vj0 >> 1) & v0mask;
          const GMUInt64 nvj1_0 = ~ vj1       & v1mask;
          const GMUInt64 nvj1_1 = ~(vj1 >> 1) & v1mask;

          const int r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 & nvj0_0) |
                                      ( (nvi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 & nvj1_0) |
                                      ( (nvi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 & nvj1_0) |
                                      ( (nvi1_1 & nvj1_1) << 1 ));
          const int r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 &  vj0_0) |
                                      ( (nvi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 &  vj1_0) |
                                      ( (nvi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 &  vj1_0) |
                                      ( (nvi1_1 &  vj1_1) << 1 ));
          const int r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 & nvj0_0) |
                                      ( ( vi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 & nvj1_0) |
                                      ( ( vi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 & nvj1_0) |
                                      ( ( vi1_1 & nvj1_1) << 1 ));
          const int r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 &  vj0_0) |
                                      ( ( vi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 &  vj1_0) |
                                      ( ( vi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 &  vj1_0) |
                                      ( ( vi1_1 &  vj1_1) << 1 ));

          /*---Accumulate---*/

          sum.data[0] += GMTally1_encode(r00, r01);
          sum.data[1] += GMTally1_encode(r10, r11);
        } /*---for pvfl---*/
        if (GMEnv_all2all(env)) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_linalg_set_matrix_zero_start(metrics_buf, env);

    /*---Perform pseudo GEMM---*/

    gm_linalg_gemm_start(
        vectors_left->num_vector_local, vectors_left->num_vector_local,
        vectors_left->num_packedval_field_local,
        vectors_left_buf->d, vectors_left->num_packedval_field_local,
        vectors_right_buf->d, vectors_left->num_packedval_field_local,
        metrics_buf->d, vectors_left->num_vector_local, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Start calculation of numerators, 2-way generic---*/

// NOTE: unlike the 3-way case, this function does not retrieve the
// metrics_buf from the GPU.

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* metrics,
                                      GMMirroredBuf* vectors_left_buf,
                                      GMMirroredBuf* vectors_right_buf,
                                      GMMirroredBuf* metrics_buf,
                                      int j_block,
                                      bool do_compute_triang_only,
                                      GMEnv* env) {
  GMAssertAlways(vectors_left && vectors_right && metrics && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      gm_compute_czek_numerators_2way_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      gm_compute_ccc_numerators_2way_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);
    } break;
    default:
      GMInsist(env, false ? "Unimplemented." : 0);
  } /*---case---*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way Czek---*/

void gm_compute_czek_2way_combine_(
    GMMetrics* metrics,
    GMMirroredBuf* metrics_buf,
    const GMVectorSums* vector_sums_left,
    const GMVectorSums* vector_sums_right,
    int j_block, bool do_compute_triang_only, GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf);
  GMAssertAlways(vector_sums_left && vector_sums_right && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  // NOTE: here and elsewhere, vector_sums_left and vector_sums_right
  // may sometimes be aliases.  Since they are read-only, there should
  // be no danger of the compiler generating incorrect code.
  // (cf. https://en.wikipedia.org/wiki/Pointer_aliasing)
  // However by accounting for this one might be able to in principle
  // remove a load instruction to improve performance.

  const int nvl = metrics->num_vector_local;
  const GMVectorSums* vs_l = vector_sums_left;
  const GMVectorSums* vs_r = vector_sums_right;

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU && GMEnv_all2all(env)) {
    /*----------------------------------------*/

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vj = GMVectorSums_sum(vs_r, j, env);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
            GMMetrics_float_get_all2all_2(metrics, i, j, j_block, env);
        const GMFloat vi = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vi < vj ?  vi + vj : vj + vi;
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block,
                                      2 * numer / denom, env);
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/
        /*---TODO: here and elsewhere check for unlikely case denom is/nearly
         * zero---*/

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vj = GMVectorSums_sum(vs_r, j, env);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer = GMMetrics_float_get_2(metrics, i, j, env);
        const GMFloat vi = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vi < vj ?  vi + vj : vj + vi;
        GMMetrics_float_set_2(metrics, i, j, 2 * numer / denom, env);
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (GMEnv_all2all(env)) {
    /*----------------------------------------*/

    if (do_compute_triang_only) {
      for (int j = 0; j < nvl; ++j) {
        const GMFloat vj = GMVectorSums_sum(vs_r, j, env);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMFloat numer =
            GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
          const GMFloat vi = GMVectorSums_sum(vs_l, i, env);
          const GMFloat denom = vi < vj ? vi + vj : vj + vi;
          GMMetrics_float_set_all2all_2(metrics, i, j, j_block,
                                        2 * numer / denom, env);
        } /*---for i---*/
        metrics->num_elts_local_computed += i_max;
      }   /*---for j---*/
    } else {
#pragma omp parallel for collapse(2)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMFloat vj = GMVectorSums_sum(vs_r, j, env);
          const GMFloat numer =
            GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
          const GMFloat vi = GMVectorSums_sum(vs_l, i, env);
          const GMFloat denom = vi < vj ? vi + vj : vj + vi;
          GMMetrics_float_set_all2all_2(metrics, i, j, j_block,
                                        2 * numer / denom, env);
        } /*---for i---*/
      }   /*---for j---*/
      metrics->num_elts_local_computed += nvl * (size_t)nvl;
    }

    /*----------------------------------------*/
  } else {
    /*----------------------------------------*/

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vj = GMVectorSums_sum(vs_r, j, env);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
          GMMirroredBuf_elt<GMFloat>(metrics_buf, i, j);
        const GMFloat vi = GMVectorSums_sum(vs_l, i, env);
        const GMFloat denom = vi < vj ? vi + vj : vj + vi;
        GMMetrics_float_set_2(metrics, i, j, 2 * numer / denom, env);
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way CCC---*/

void gm_compute_ccc_2way_combine_(GMMetrics* metrics,
                                  GMMirroredBuf* metrics_buf,
                                  const GMVectorSums* vector_sums_left,
                                  const GMVectorSums* vector_sums_right,
                                  int j_block, bool do_compute_triang_only,
                                  GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf);
  GMAssertAlways(vector_sums_left && vector_sums_right && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  const int nvl = metrics->num_vector_local;
  const GMVectorSums* vs_l = vector_sums_left;
  const GMVectorSums* vs_r = vector_sums_right;

  /*---Copy from metrics_buffer for GPU case---*/

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*--------------------*/
    if (GMEnv_all2all(env)) {
      /*--------------------*/

      if (do_compute_triang_only) {
        for (int j = 0; j < nvl; ++j) {
          const int i_max = j;
          for (int i = 0; i < i_max; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
#ifdef GM_ASSERTIONS_ON
            if (!env->sparse) {
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              GMAssert((GMUInt64)r00 + (GMUInt64)r01 + (GMUInt64)r10 +
                           (GMUInt64)r11 ==
                       (GMUInt64)(4 * metrics->num_field_active));
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              GMAssert((GMUInt64)r10 + (GMUInt64)r11 == (GMUInt64)(2 * si1));
              GMAssert((GMUInt64)r01 + (GMUInt64)r11 == (GMUInt64)(2 * sj1));
            }
#endif
          } /*---for i---*/
        }   /*---for j---*/
      } else {
        #pragma omp parallel for collapse(2)
        for (int j = 0; j < nvl; ++j) {
          for (int i = 0; i < nvl; ++i) {
            const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
            GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
#ifdef GM_ASSERTIONS_ON
            if (!env->sparse) {
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              GMAssert((GMUInt64)r00 + (GMUInt64)r01 + (GMUInt64)r10 +
                           (GMUInt64)r11 ==
                       (GMUInt64)(4 * metrics->num_field_active));
              const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
              const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
              GMAssert((GMUInt64)r10 + (GMUInt64)r11 == (GMUInt64)(2 * si1));
              GMAssert((GMUInt64)r01 + (GMUInt64)r11 == (GMUInt64)(2 * sj1));
            }
#endif
          } /*---for i---*/
        }   /*---for j---*/
     }

      /*--------------------*/
    } else /*---(!GMEnv_all2all(env))---*/ {
      /*--------------------*/
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally2x2 value =
              GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j);
          GMMetrics_tally2x2_set_2(metrics, i, j, value, env);
#ifdef GM_ASSERTIONS_ON
          if (!env->sparse) {
            const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
            const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
            const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
            const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
            GMAssert((GMUInt64)r00 + (GMUInt64)r01 + (GMUInt64)r10 +
                         (GMUInt64)r11 ==
                     (GMUInt64)(4 * metrics->num_field_active));
            const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
            const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
            GMAssert((GMUInt64)r10 + (GMUInt64)r11 == (GMUInt64)(2 * si1));
            GMAssert((GMUInt64)r01 + (GMUInt64)r11 == (GMUInt64)(2 * sj1));
          }
#endif
        } /*---for i---*/
      }   /*---for j---*/
      /*--------------------*/
    } /*---if---*/
    /*--------------------*/
  }

  /*---Compute multipliers---*/

  /*--------------------*/
  if (GMEnv_all2all(env)) {
    /*--------------------*/

    if (do_compute_triang_only) {
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
        metrics->num_elts_local_computed += i_max;
      }   /*---for j---*/
    } else {
      #pragma omp parallel for collapse(2)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          GMMetrics_float2_S_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
          if (env->sparse) {
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            GMMetrics_float2_C_set_all2all_2(metrics, i, j, j_block, ci_cj, env);
          } /*---if sparse---*/
        }   /*---for i---*/
      }   /*---for j---*/
      metrics->num_elts_local_computed += nvl * (size_t)nvl;
   }

    /*--------------------*/
  } else /*---(!GMEnv_all2all(env))---*/ {
    /*--------------------*/
    for (int j = 0; j < nvl; ++j) {
      const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_r, j, env);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_l, i, env);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
        GMMetrics_float2_S_set_2(metrics, i, j, si1_sj1, env);
        if (env->sparse) {
          const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_r, j, env);
          const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_l, i, env);
          const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
          GMMetrics_float2_C_set_2(metrics, i, j, ci_cj, env);
        } /*---if sparse---*/
      } /*---for i---*/
      metrics->num_elts_local_computed += i_max;
    }   /*---for j---*/
    /*--------------------*/
  } /*---if---*/
  /*--------------------*/
}

//=============================================================================
/*---Combine nums and denoms on CPU to get final result, 2-way generic---*/

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredBuf* metrics_buf,
                             const GMVectorSums* vector_sums_left,
                             const GMVectorSums* vector_sums_right,
                             int j_block,
                             bool do_compute_triang_only,
                             GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf);
  GMAssertAlways(vector_sums_left && vector_sums_right && env);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_2);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      gm_compute_czek_2way_combine_(metrics, metrics_buf,
                                    vector_sums_left, vector_sums_right,
                                    j_block, do_compute_triang_only, env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      gm_compute_ccc_2way_combine_(metrics, metrics_buf,
                                   vector_sums_left, vector_sums_right,
                                   j_block, do_compute_triang_only, env);
    } break;
    default:
      GMInsist(env, false ? "Unimplemented." : 0);
  } /*---case---*/
}

//=============================================================================

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

//-----------------------------------------------------------------------------
