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
// Start calculation of numerators, 2-way Czekanowski.

void gm_compute_2way_proc_nums_czek_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NUM_WAY::_2);
  COMET_INSIST(!env->is_using_linalg());

  // ----------------------------------
  if (env->all2all()) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } // for k
        // Update metrics array.
        Metrics_elt_2<GMFloat>(*metrics, i, j, j_block, *env) = metric;
      } // for i
    }   // for j

    // ----------------------------------
  } else {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    // Perform pseudo GEMM.

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (int f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } // for k
        // Update metrics array.
        Metrics_elt_2<GMFloat>(*metrics, i, j, env->proc_num_vector(), *env) = metric;
      } // for i
    }   // for j

    // ----------------------------------
  } // if
  // ----------------------------------
}

//=============================================================================
// Start calculation of numerators, 2-way CCC.

void gm_compute_2way_proc_nums_ccc_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NUM_WAY::_2);

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // ----------------------------------
  if (env->compute_method() == ComputeMethod::REF) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

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

            MFT::add(sum.data[0], r00, r01);
            MFT::add(sum.data[1], r10, r11);

          } // if ! unknown
        } // for f

        // Update metrics array.

        const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
        Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block_eff, *env) = sum;
        //}
      } // for j
    }   // for i

    // ----------------------------------
  } else { // ComputeMethod::CPU && !env->is_using_linalg()
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
       "num_proc_field>1 for CPU compute_method for this case not supported");

    // Perform pseudo GEMM.

    /* clang-format off */

    const int cbpe = env->counted_bits_per_elt();

    const int pad_adjustment = cbpe * cbpe * metrics->dm->num_pad_field_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          // Extract input values to process.

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          // Compute masks.

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

          // Get even/odd bits for each seminibble, masked to active.

          const uint64_t vi0_0 =  vi0       & v0mask;
          const uint64_t vi0_1 = (vi0 >> 1) & v0mask;
          const uint64_t vi1_0 =  vi1       & v1mask;
          const uint64_t vi1_1 = (vi1 >> 1) & v1mask;
          const uint64_t vj0_0 =  vj0       & v0mask;
          const uint64_t vj0_1 = (vj0 >> 1) & v0mask;
          const uint64_t vj1_0 =  vj1       & v1mask;
          const uint64_t vj1_1 = (vj1 >> 1) & v1mask;

          // Get complements of even/odd bits for each seminibble; mask.

          const uint64_t nvi0_0 = ~ vi0       & v0mask;
          const uint64_t nvi0_1 = ~(vi0 >> 1) & v0mask;
          const uint64_t nvi1_0 = ~ vi1       & v1mask;
          const uint64_t nvi1_1 = ~(vi1 >> 1) & v1mask;
          const uint64_t nvj0_0 = ~ vj0       & v0mask;
          const uint64_t nvj0_1 = ~(vj0 >> 1) & v0mask;
          const uint64_t nvj1_0 = ~ vj1       & v1mask;
          const uint64_t nvj1_1 = ~(vj1 >> 1) & v1mask;

          const int r00 = utils::popc64((nvi0_0 & nvj0_0) |
                                      ( (nvi0_0 & nvj0_1) << 1 )) +
                          utils::popc64((nvi0_1 & nvj0_0) |
                                      ( (nvi0_1 & nvj0_1) << 1 )) +
                          utils::popc64((nvi1_0 & nvj1_0) |
                                      ( (nvi1_0 & nvj1_1) << 1 )) +
                          utils::popc64((nvi1_1 & nvj1_0) |
                                      ( (nvi1_1 & nvj1_1) << 1 ));
          const int r01 = utils::popc64((nvi0_0 &  vj0_0) |
                                      ( (nvi0_0 &  vj0_1) << 1 )) +
                          utils::popc64((nvi0_1 &  vj0_0) |
                                      ( (nvi0_1 &  vj0_1) << 1 )) +
                          utils::popc64((nvi1_0 &  vj1_0) |
                                      ( (nvi1_0 &  vj1_1) << 1 )) +
                          utils::popc64((nvi1_1 &  vj1_0) |
                                      ( (nvi1_1 &  vj1_1) << 1 ));
          const int r10 = utils::popc64(( vi0_0 & nvj0_0) |
                                      ( ( vi0_0 & nvj0_1) << 1 )) +
                          utils::popc64(( vi0_1 & nvj0_0) |
                                      ( ( vi0_1 & nvj0_1) << 1 )) +
                          utils::popc64(( vi1_0 & nvj1_0) |
                                      ( ( vi1_0 & nvj1_1) << 1 )) +
                          utils::popc64(( vi1_1 & nvj1_0) |
                                      ( ( vi1_1 & nvj1_1) << 1 ));
          const int r11 = utils::popc64(( vi0_0 &  vj0_0) |
                                      ( ( vi0_0 &  vj0_1) << 1 )) +
                          utils::popc64(( vi0_1 &  vj0_0) |
                                      ( ( vi0_1 &  vj0_1) << 1 )) +
                          utils::popc64(( vi1_0 &  vj1_0) |
                                      ( ( vi1_0 &  vj1_1) << 1 )) +
                          utils::popc64(( vi1_1 &  vj1_0) |
                                      ( ( vi1_1 &  vj1_1) << 1 ));

          // Accumulate

          MFT::add(sum.data[0], r00, r01);
          MFT::add(sum.data[1], r10, r11);

        } // for pvfl

        // Adjust for pad

#ifdef COMET_ASSERTIONS_ON
        GMTally2x2 sum_old = sum;
#endif
        MFT::subtract(sum.data[0], pad_adjustment, 0);
#ifdef COMET_ASSERTIONS_ON
        COMET_ASSERT(GMTally2x2_get(sum_old, 0, 0) ==
                 GMTally2x2_get(sum, 0, 0) + pad_adjustment);
#endif

        // Update metrics array.

        const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
        Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block_eff, *env) = sum;
        //}
      } // for j
    }   // for i

    /* clang-format on */

    // ----------------------------------
  } // if
  // ----------------------------------
}

//=============================================================================
// Start calculation of numerators, 2-way DUO.

void gm_compute_2way_proc_nums_duo_start_(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NUM_WAY::_2);

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // ----------------------------------
  if (env->compute_method() == ComputeMethod::REF) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

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

            MFT::add(sum.data[0], r00, r01);
            MFT::add(sum.data[1], r10, r11);

          } // if ! unknown
        } // for f

        // Update metrics array.

        const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
        Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block_eff, *env) = sum;
      } // for j
    }   // for i

    // ----------------------------------
  } else if (!env->is_using_xor()) { // ComputeMethod::CPU && !env->is_using_linalg()
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    // Perform pseudo GEMM.

    /* clang-format off */

    const int cbpe = env->counted_bits_per_elt();

    const int pad_adjustment = cbpe * cbpe * metrics->dm->num_pad_field_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          // Extract input values to process.

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          // Compute masks.

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

          // Get even/odd bits for each seminibble, masked to active.

          const uint64_t vi0_0 =  vi0       & v0mask;
          const uint64_t vi1_0 =  vi1       & v1mask;
          const uint64_t vj0_0 =  vj0       & v0mask;
          const uint64_t vj1_0 =  vj1       & v1mask;

          // Get complements of even/odd bits for each seminibble; mask.

          const uint64_t nvi0_0 = ~ vi0       & v0mask;
          const uint64_t nvi1_0 = ~ vi1       & v1mask;
          const uint64_t nvj0_0 = ~ vj0       & v0mask;
          const uint64_t nvj1_0 = ~ vj1       & v1mask;

          const int r00 = utils::popc64((nvi0_0 & nvj0_0) |
                                      ( (nvi1_0 & nvj1_0) << 1 ));
          const int r01 = utils::popc64((nvi0_0 &  vj0_0) |
                                      ( (nvi1_0 &  vj1_0) << 1 ));
          const int r10 = utils::popc64(( vi0_0 & nvj0_0) |
                                      ( ( vi1_0 & nvj1_0) << 1 ));
          const int r11 = utils::popc64(( vi0_0 &  vj0_0) |
                                      ( ( vi1_0 &  vj1_0) << 1 ));

          // Accumulate

          MFT::add(sum.data[0], r00, r01);
          MFT::add(sum.data[1], r10, r11);

        } // for pvfl

        // Adjust for pad

#ifdef COMET_ASSERTIONS_ON
        GMTally2x2 sum_old = sum;
#endif
        MFT::subtract(sum.data[0], pad_adjustment, 0);
#ifdef COMET_ASSERTIONS_ON
        COMET_ASSERT(GMTally2x2_get(sum_old, 0, 0) ==
                 GMTally2x2_get(sum, 0, 0) + pad_adjustment);
#endif

        // Update metrics array.

        const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
        Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block_eff, *env) = sum;
      } // for j
    }   // for i

    /* clang-format on */

    // ----------------------------------
  } else { // ComputeMethod::CPU && !env->is_using_linalg() && env->is_using_xor()
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    // Perform pseudo GEMM.

    /* clang-format off */

    const int cbpe = env->counted_bits_per_elt();

    const int pad_adjustment = cbpe * cbpe * metrics->dm->num_pad_field_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        const int npvfl = vectors_left->num_packedval_field_local;
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          // Extract input values to process.

          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, pvfl, i,
                                                       env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, pvfl, j,
                                                       env);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          // Compute masks to sample the single needed bit from each seminibble,
          // and to ignore undefined entries.

          const uint64_t oddbits = 0x5555555555555555;

          const uint64_t vi0mask =
                       env->sparse() ? (vi0 | ~(vi0 >> 1)) & oddbits : oddbits;

          const uint64_t vi1mask =
                       env->sparse() ? (vi1 | ~(vi1 >> 1)) & oddbits : oddbits;

          const uint64_t vj0mask =
                       env->sparse() ? (vj0 | ~(vj0 >> 1)) & oddbits : oddbits;

          const uint64_t vj1mask =
                       env->sparse() ? (vj1 | ~(vj1 >> 1)) & oddbits : oddbits;

          // Extract elts that are a "1" bit (=01).

          const uint64_t pvi0 =  vi0  & vi0mask;
          const uint64_t pvi1 =  vi1  & vi1mask;
          const uint64_t pvj0 =  vj0  & vj0mask;
          const uint64_t pvj1 =  vj1  & vj1mask;

          // Extract elts that are an "0" bit (=00).

          const uint64_t nvi0 = ~vi0  & vi0mask;
          const uint64_t nvi1 = ~vi1  & vi1mask;
          const uint64_t nvj0 = ~vj0  & vj0mask;
          const uint64_t nvj1 = ~vj1  & vj1mask;

          // Combine lower, upper words - each only uses odd bits - make packed.

          const uint64_t pvi = pvi0 | (pvi1 << 1);
          const uint64_t pvj = pvj0 | (pvj1 << 1);
          const uint64_t nvi = nvi0 | (nvi1 << 1);
          const uint64_t nvj = nvj0 | (nvj1 << 1);

          // Compute values using xor.
          // NOTE: these will need adjustment later in the execution.

          const int r00 = utils::popc64(nvi ^ nvj);
          const int r01 = utils::popc64(nvi ^ pvj);
          const int r10 = utils::popc64(pvi ^ nvj);
          const int r11 = utils::popc64(pvi ^ pvj);

          // Accumulate

          MFT::add(sum.data[0], r00, r01);
          MFT::add(sum.data[1], r10, r11);

        } // for pvfl

        // Adjust for pad
        // NOTE this works differently from other cases because of xor

        MFT::subtract(sum.data[0], 0, pad_adjustment);
        MFT::subtract(sum.data[1], pad_adjustment, 0);

        // Update metrics array.

        const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
        Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block_eff, *env) = sum;
      } // for j
    }   // for i

    /* clang-format on */

    // ----------------------------------
  } // if
  // ----------------------------------
}

//=============================================================================
// Start calculation of numerators, 2-way generic.

// NOTE: unlike the 3-way case, this function does not retrieve the
// metrics_buf from the GPU.

void gm_compute_2way_proc_nums_start(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  VectorSums* vector_sums_left,
  VectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NUM_WAY::_2);

  if (env->is_using_linalg()) {

    // Perform pseudo GEMM.

    gm_linalg_gemm_start(
      vectors_left->num_vector_local,
      vectors_left->num_vector_local,
      vectors_left->num_packedval_field_local,
      vectors_left_buf,
      vectors_right_buf,
      metrics_buf,
      vector_sums_left->sums(), vector_sums_right->sums(),
      vector_sums_left->counts(), vector_sums_right->counts(),
      vectors_left->dm, env);

    return;

  } // if

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
      COMET_INSIST_INTERFACE(env, false && "Selected metric_type unimplemented.");
  } // case
}

//=============================================================================
// Finish calculation of numerators, 2-way generic.

void gm_compute_2way_proc_nums_wait(
  GMVectors* vectors_left,
  GMVectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  VectorSums* vector_sums_left,
  VectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NUM_WAY::_2);

  if (env->is_using_linalg()) {

    gm_linalg_gemm_wait(
      vectors_left->num_vector_local,
      vectors_left->num_vector_local,
      vectors_left->num_packedval_field_local,
      vectors_left_buf,
      vectors_right_buf,
      metrics_buf,
      vector_sums_left->sums(), vector_sums_right->sums(),
      vector_sums_left->counts(), vector_sums_right->counts(),
      vectors_left->dm, env);

  } // if
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
