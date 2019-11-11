//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_nums.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate metrics numerators, 2-way, for a single block.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdint"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way_block_nums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Start calculation of numerators, 2-way Czekanowski---*/

void gm_compute_2way_proc_nums_czek_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(vectors_left && vectors_right && metrics && env);
  GMInsist(j_block >= 0 && j_block < env->num_block_vector());
  GMInsist(env->num_way() == NUM_WAY::_2);

  /*----------------------------------------*/
  if (env->compute_method() != ComputeMethod::GPU && env->all2all()) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    /*---Perform pseudo GEMM---*/

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (env->compute_method() != ComputeMethod::GPU) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    /*---Perform pseudo GEMM---*/

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_2(metrics, i, j, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else /* if (env->compute_method() == ComputeMethod::GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_linalg_set_matrix_zero_start(metrics_buf, env);

    /*---Perform pseudo GEMM---*/

    gm_linalg_gemm_start(
      vectors_left->num_vector_local,
      vectors_left->num_vector_local,
      vectors_left->num_field_local,
      vectors_left_buf->active, vectors_left->num_field_local,
      vectors_right_buf->active, vectors_left->num_field_local,
      metrics_buf->active, vectors_left->num_vector_local,
      vectors_left->dm, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Start calculation of numerators, 2-way CCC---*/

void gm_compute_2way_proc_nums_ccc_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(vectors_left && vectors_right && metrics && env);
  GMInsist(j_block >= 0 && j_block < env->num_block_vector());
  GMInsist(env->num_way() == NUM_WAY::_2);

  /*----------------------------------------*/
  if (env->compute_method() == ComputeMethod::REF) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    /*---Perform pseudo GEMM---*/

    const int nfal = vectors_left->dm->num_field_active_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        for (int f = 0; f < nfal; ++f) {
          const GMBits2 vi = GMVectors_bits2_get(vectors_left, f, i, env);
          const GMBits2 vj = GMVectors_bits2_get(vectors_right, f, j, env);
          const bool unknown_i = env->sparse() ? vi == GM_2BIT_UNKNOWN
                                             : false;
          const bool unknown_j = env->sparse() ? vj == GM_2BIT_UNKNOWN
                                             : false;

          if ( ! unknown_i && ! unknown_j ) {

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

            // NOTE: Since the sum of all 4 of these relative
            // cooccurences is 4, we really only need to compute 3 of them.
            //  Then the last one is just 4 minus the rest (non-sparse case)

#if DOUG_WAY
//TODO: work on this as a possibly faster way.
            const int vi1 = (vi & 3) != 0;
            const int vi0 = ((~vi) & 3) != 0;
            const int vj1 = (vj & 3) != 0;
            const int vj0 = ((~vj) & 3) != 0;

            const int a11 = vi1 & vj1;

            const int r11 = a11 +
#endif

            // Accumulate

            sum.data[0] += GMTally1_encode(r00, r01);
            sum.data[1] += GMTally1_encode(r10, r11);

          } /*---if ! unknown---*/
        } /*---for f---*/
        if (env->all2all()) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /*----------------------------------------*/
  } else if (env->compute_method() == ComputeMethod::CPU &&
             !env->is_using_linalg()) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
       "num_proc_field>1 for CPU compute_method for this case not supported");

    /*---Perform pseudo GEMM---*/

    /* clang-format off */

    const int cbpe = 2;

    const int pad_adjustment = cbpe * cbpe * metrics->dm->num_pad_field_local;
    const GMFloat float_pad_adjustment = GMTally1_encode(pad_adjustment, 0);

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          /*---Extract input values to process---*/

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          /*---Compute masks---*/

          const uint64_t oddbits = 0x5555555555555555;

          const uint64_t vi0mask =
                       (env->sparse() ? (vi0 | ~(vi0 >> 1)) & oddbits : oddbits);

          const uint64_t vi1mask =
                       (env->sparse() ? (vi1 | ~(vi1 >> 1)) & oddbits : oddbits);

          const uint64_t vj0mask =
                       (env->sparse() ? (vj0 | ~(vj0 >> 1)) & oddbits : oddbits);

          const uint64_t vj1mask =
                       (env->sparse() ? (vj1 | ~(vj1 >> 1)) & oddbits : oddbits);

          const uint64_t v0mask = vi0mask & vj0mask;
          const uint64_t v1mask = vi1mask & vj1mask;

          /*---Get even/odd bits for each seminibble, masked to active---*/

          const uint64_t vi0_0 =  vi0       & v0mask;
          const uint64_t vi0_1 = (vi0 >> 1) & v0mask;
          const uint64_t vi1_0 =  vi1       & v1mask;
          const uint64_t vi1_1 = (vi1 >> 1) & v1mask;
          const uint64_t vj0_0 =  vj0       & v0mask;
          const uint64_t vj0_1 = (vj0 >> 1) & v0mask;
          const uint64_t vj1_0 =  vj1       & v1mask;
          const uint64_t vj1_1 = (vj1 >> 1) & v1mask;

          /*---Get complements of even/odd bits for each seminibble; mask---*/

          const uint64_t nvi0_0 = ~ vi0       & v0mask;
          const uint64_t nvi0_1 = ~(vi0 >> 1) & v0mask;
          const uint64_t nvi1_0 = ~ vi1       & v1mask;
          const uint64_t nvi1_1 = ~(vi1 >> 1) & v1mask;
          const uint64_t nvj0_0 = ~ vj0       & v0mask;
          const uint64_t nvj0_1 = ~(vj0 >> 1) & v0mask;
          const uint64_t nvj1_0 = ~ vj1       & v1mask;
          const uint64_t nvj1_1 = ~(vj1 >> 1) & v1mask;

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

          // Accumulate

          sum.data[0] += GMTally1_encode(r00, r01);
          sum.data[1] += GMTally1_encode(r10, r11);

        } /*---for pvfl---*/

        // Adjust for pad

#ifdef GM_ASSERTIONS_ON
        GMTally2x2 sum_old = sum;
#endif
        sum.data[0] -= float_pad_adjustment;
#ifdef GM_ASSERTIONS_ON
        GMAssert(GMTally2x2_get(sum_old, 0, 0) ==
                 GMTally2x2_get(sum, 0, 0) + pad_adjustment);
#endif

        if (env->all2all()) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (env->is_using_linalg()) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_linalg_set_matrix_zero_start(metrics_buf, env);

    /*---Perform pseudo GEMM---*/

    gm_linalg_gemm_start(
      vectors_left->num_vector_local,
      vectors_left->num_vector_local,
      vectors_left->num_packedval_field_local,
      vectors_left_buf->active, vectors_left->num_packedval_field_local,
      vectors_right_buf->active, vectors_left->num_packedval_field_local,
      metrics_buf->active, vectors_left->num_vector_local,
      vectors_left->dm, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Start calculation of numerators, 2-way DUO---*/

void gm_compute_2way_proc_nums_duo_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(vectors_left && vectors_right && metrics && env);
  GMInsist(j_block >= 0 && j_block < env->num_block_vector());
  GMInsist(env->num_way() == NUM_WAY::_2);

  /*----------------------------------------*/
  if (env->compute_method() == ComputeMethod::REF) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    /*---Perform pseudo GEMM---*/

    const int nfal = vectors_left->dm->num_field_active_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        for (int f = 0; f < nfal; ++f) {
          const GMBits2 vi = GMVectors_bits2_get(vectors_left, f, i, env);
          const GMBits2 vj = GMVectors_bits2_get(vectors_right, f, j, env);
          const bool unknown_i = env->sparse() ? vi == GM_2BIT_UNKNOWN
                                             : false;
          const bool unknown_j = env->sparse() ? vj == GM_2BIT_UNKNOWN
                                             : false;

          if ( ! unknown_i && ! unknown_j ) {

            /* clang-format off */
            const int r00 = ( ( !(vi & 1) ) && ( !(vj & 1) ) );
            const int r01 = ( ( !(vi & 1) ) && (  (vj & 1) ) );
            const int r10 = ( (  (vi & 1) ) && ( !(vj & 1) ) );
            const int r11 = ( (  (vi & 1) ) && (  (vj & 1) ) );
            /* clang-format on */

            // NOTE: Since the sum of all 4 of these relative
            // cooccurences is 1, we really only need to compute 3 of them.
            //  Then the last one is just 1 minus the rest (non-sparse case)

            // Accumulate

            sum.data[0] += GMTally1_encode(r00, r01);
            sum.data[1] += GMTally1_encode(r10, r11);

          } /*---if ! unknown---*/
        } /*---for f---*/
        if (env->all2all()) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /*----------------------------------------*/
  } else if (env->compute_method() == ComputeMethod::CPU &&
             !env->is_using_linalg()) {
    /*----------------------------------------*/

    GMInsistInterface(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    /*---Perform pseudo GEMM---*/

    /* clang-format off */

    const int cbpe = 1;

    const int pad_adjustment = cbpe * cbpe * metrics->dm->num_pad_field_local;
    const GMFloat float_pad_adjustment = GMTally1_encode(pad_adjustment, 0);

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          /*---Extract input values to process---*/

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          /*---Compute masks---*/

          const uint64_t oddbits = 0x5555555555555555;

          const uint64_t vi0mask =
                       (env->sparse() ? (vi0 | ~(vi0 >> 1)) & oddbits : oddbits);

          const uint64_t vi1mask =
                       (env->sparse() ? (vi1 | ~(vi1 >> 1)) & oddbits : oddbits);

          const uint64_t vj0mask =
                       (env->sparse() ? (vj0 | ~(vj0 >> 1)) & oddbits : oddbits);

          const uint64_t vj1mask =
                       (env->sparse() ? (vj1 | ~(vj1 >> 1)) & oddbits : oddbits);

          const uint64_t v0mask = vi0mask & vj0mask;
          const uint64_t v1mask = vi1mask & vj1mask;

          /*---Get even/odd bits for each seminibble, masked to active---*/

          const uint64_t vi0_0 =  vi0       & v0mask;
          const uint64_t vi1_0 =  vi1       & v1mask;
          const uint64_t vj0_0 =  vj0       & v0mask;
          const uint64_t vj1_0 =  vj1       & v1mask;

          /*---Get complements of even/odd bits for each seminibble; mask---*/

          const uint64_t nvi0_0 = ~ vi0       & v0mask;
          const uint64_t nvi1_0 = ~ vi1       & v1mask;
          const uint64_t nvj0_0 = ~ vj0       & v0mask;
          const uint64_t nvj1_0 = ~ vj1       & v1mask;

          const int r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi1_0 & nvj1_0) << 1 ));
          const int r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi1_0 &  vj1_0) << 1 ));
          const int r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi1_0 & nvj1_0) << 1 ));
          const int r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi1_0 &  vj1_0) << 1 ));
//printf("cpu %i\n", (int)r01);

          // Accumulate

          sum.data[0] += GMTally1_encode(r00, r01);
          sum.data[1] += GMTally1_encode(r10, r11);

        } /*---for pvfl---*/

        // Adjust for pad

#ifdef GM_ASSERTIONS_ON
        GMTally2x2 sum_old = sum;
#endif
        sum.data[0] -= float_pad_adjustment;
#ifdef GM_ASSERTIONS_ON
        GMAssert(GMTally2x2_get(sum_old, 0, 0) ==
                 GMTally2x2_get(sum, 0, 0) + pad_adjustment);
#endif

        if (env->all2all()) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (env->is_using_linalg()) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_linalg_set_matrix_zero_start(metrics_buf, env);

    /*---Perform pseudo GEMM---*/

    gm_linalg_gemm_start(
      vectors_left->num_vector_local,
      vectors_left->num_vector_local,
      vectors_left->num_packedval_field_local,
      vectors_left_buf->active, vectors_left->num_packedval_field_local,
      vectors_right_buf->active, vectors_left->num_packedval_field_local,
      metrics_buf->active, vectors_left->num_vector_local,
      vectors_left->dm, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================
/*---Start calculation of numerators, 2-way generic---*/

// NOTE: unlike the 3-way case, this function does not retrieve the
// metrics_buf from the GPU.

void gm_compute_2way_proc_nums_start(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  GMMirroredBuf* vectors_left_buf,
  GMMirroredBuf* vectors_right_buf,
  GMMirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  GMEnv* env) {

  GMInsist(vectors_left && vectors_right && metrics && env);
  GMInsist(j_block >= 0 && j_block < env->num_block_vector());
  GMInsist(env->num_way() == NUM_WAY::_2);

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      gm_compute_2way_proc_nums_czek_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);
    } break;
    case MetricType::CCC: {
      gm_compute_2way_proc_nums_ccc_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);
    } break;
    case MetricType::DUO: {
      gm_compute_2way_proc_nums_duo_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);
    } break;
    default:
      GMInsistInterface(env, false && "Selected metric_type unimplemented.");
  } /*---case---*/
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
