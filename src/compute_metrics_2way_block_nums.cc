//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_nums.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate metrics numerators, 2-way, for a single block.
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

#include "cstdint"

#include "env.hh"
#include "tc.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Start calculation of numerators, 2-way Czekanowski.

static void compute_nums_nonlinalg_czek_start_(
  Vectors* vectors_left,
  Vectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);
  COMET_INSIST(!env->is_using_linalg());

  typedef GMFloat Float_t;

  // ----------------------------------
  if (env->all2all()) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        Float_t metric = 0;
        for (int f = 0; f < (int)vectors_left->dm()->num_field_local; ++f) {
          const Float_t value1 = vectors_left->elt_float_const<Float_t>(f, i);
          const Float_t value2 = vectors_right->elt_float_const<Float_t>(f, j);
          metric += value1 < value2 ? value1 : value2;
        } // for k
        // Update metrics array.
        Metrics_elt_2<Float_t>(*metrics, i, j, j_block, *env) = metric;
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
        Float_t metric = 0;
        for (int f = 0; f < (int)vectors_left->dm()->num_field_local; ++f) {
          const Float_t value1 = vectors_left->elt_float_const<Float_t>(f, i);
          const Float_t value2 = vectors_right->elt_float_const<Float_t>(f, j);
          metric += value1 < value2 ? value1 : value2;
        } // for k
        // Update metrics array.
        Metrics_elt_2<Float_t>(*metrics, i, j, env->proc_num_vector(), *env) =
          metric;
      } // for i
    }   // for j

    // ----------------------------------
  } // if (env->all2all())
  // ----------------------------------
}

//=============================================================================
// Start calculation of numerators, 2-way CCC.

static void compute_nums_nonlinalg_ccc_start_(
  Vectors* vectors_left,
  Vectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // ----------------------------------
  if (env->compute_method() == ComputeMethod::REF) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

    const int nfal = vectors_left->dm()->num_field_active_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        for (int f = 0; f < nfal; ++f) {
          const GMBits2 vi = vectors_left->bits2_get(f, i, *env);
          const GMBits2 vj = vectors_right->bits2_get(f, j, *env);
          //const GMBits2 vi = vectors_left->bits2_get(f, i, *env);
          //const GMBits2 vj = vectors_right->bits2_get(f, j, *env);
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

//printf("%i %i %f %f %f %f\n", i, j, (double)r00, (double)r01, (double)r10, (double)r11);
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
        const int npfl = vectors_left->num_packedfield_local();
        for (int pfl = 0; pfl < npfl; ++pfl) {

          // Extract input values to process.

          const GMBits2x64 vi = vectors_left->elt_bits2x64_const(pfl, i);
          const GMBits2x64 vj = vectors_right->elt_bits2x64_const(pfl, j);
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

        } // for pfl

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

static void compute_nums_nonlinalg_duo_start_(
  Vectors* vectors_left,
  Vectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // ----------------------------------
  if (env->compute_method() == ComputeMethod::REF) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    // Perform pseudo GEMM.

    const int nfal = vectors_left->dm()->num_field_active_local;

    for (int j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (int i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        for (int f = 0; f < nfal; ++f) {
          const GMBits2 vi = vectors_left->bits2_get(f, i, *env);
          const GMBits2 vj = vectors_right->bits2_get(f, j, *env);
          //const GMBits2 vi = Vectors_bits2_get(vectors_left, f, i, env);
          //const GMBits2 vj = Vectors_bits2_get(vectors_right, f, j, env);
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
        const int npfl = vectors_left->num_packedfield_local();
        for (int pfl = 0; pfl < npfl; ++pfl) {

          // Extract input values to process.

          const GMBits2x64 vi = vectors_left->elt_bits2x64_const(pfl, i);
          const GMBits2x64 vj = vectors_right->elt_bits2x64_const(pfl, j);
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

        } // for pfl

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
        const int npfl = vectors_left->num_packedfield_local();
        for (int pfl = 0; pfl < npfl; ++pfl) {

          // Extract input values to process.

          const GMBits2x64 vi = vectors_left->elt_bits2x64_const(pfl, i);
          const GMBits2x64 vj = vectors_right->elt_bits2x64_const(pfl, j);
          const uint64_t vi0 = vi.data[0];
          const uint64_t vi1 = vi.data[1];
          const uint64_t vj0 = vj.data[0];
          const uint64_t vj1 = vj.data[1];

          // Compute masks to sample the single needed bit from each seminibble,
          // and to ignore undefined vector entries.

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

        } // for pfl

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

void ComputeMetrics2WayBlock::compute_nums_start(
  Vectors* vectors_left,
  Vectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  VectorSums* vector_sums_left,
  VectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  MagmaWrapper& magma_wrapper,
  GemmShapes& gemm_shapes,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  if (env->is_using_linalg()) {

    //GemmShape2Way gemm_shape(do_compute_triang_only,
    //  metrics->dm()->num_vector_local, metrics->dm()->num_vector_active,
    //  env->proc_num_vector(), j_block);
//printf("%i  %i %i\n", System::proc_num(), (int)metrics->dm()->num_vector_active_local, do_compute_triang_only);
    ////GemmShapes gemm_shapes(&gemm_shape, NULL, NULL);
    //GemmShapes gemm_shapes(&gemm_shape);

    LinAlg::gemm_start(
      vectors_left->num_vector_local(),
      vectors_left->num_vector_local(),
      vectors_left->num_packedfield_local(),
      vectors_left_buf,
      vectors_right_buf,
      metrics_buf,
      vector_sums_left->sums(), vector_sums_right->sums(),
      vector_sums_left->counts(), vector_sums_right->counts(),
      *(vectors_left->dm()), magma_wrapper, gemm_shapes, *env);

  } else if(env->metric_type() == MetricType::CZEK) {

    compute_nums_nonlinalg_czek_start_(
        vectors_left, vectors_right, metrics, vectors_left_buf,
        vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
        env);

  } else if(env->metric_type() == MetricType::CCC) {

    compute_nums_nonlinalg_ccc_start_(
        vectors_left, vectors_right, metrics, vectors_left_buf,
        vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
        env);

  } else if(env->metric_type() == MetricType::DUO) {

      compute_nums_nonlinalg_duo_start_(
          vectors_left, vectors_right, metrics, vectors_left_buf,
          vectors_right_buf, metrics_buf, j_block, do_compute_triang_only,
          env);

  } else {

      COMET_INSIST_INTERFACE(env, false &&
        "Selected metric_type unimplemented.");

  } // if
}

//=============================================================================
// Finish calculation of numerators, 2-way generic.

void ComputeMetrics2WayBlock::compute_nums_wait(
  Vectors* vectors_left,
  Vectors* vectors_right,
  GMMetrics* metrics,
  MirroredBuf* vectors_left_buf,
  MirroredBuf* vectors_right_buf,
  MirroredBuf* metrics_buf,
  VectorSums* vector_sums_left,
  VectorSums* vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  GemmShapes& gemm_shapes,
  CEnv* env) {

  COMET_INSIST(vectors_left && vectors_right && metrics && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  if (env->is_using_linalg()) {

    //GemmShape2Way gemm_shape(do_compute_triang_only,
    //  metrics->dm()->num_vector_local, metrics->dm()->num_vector_active,
    //  env->proc_num_vector(), j_block);
    ////GemmShapes gemm_shapes(&gemm_shape, NULL, NULL);
    //GemmShapes gemm_shapes(&gemm_shape);

    LinAlg::gemm_wait(
      vectors_left->num_vector_local(),
      vectors_left->num_vector_local(),
      vectors_left->num_packedfield_local(),
      vectors_left_buf,
      vectors_right_buf,
      metrics_buf,
      vector_sums_left->sums(), vector_sums_right->sums(),
      vector_sums_left->counts(), vector_sums_right->counts(),
      *(vectors_left->dm()), gemm_shapes, *env);

  } // if
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
