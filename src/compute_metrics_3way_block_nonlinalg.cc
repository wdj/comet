//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block_nongpu.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block, non-GPU case.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdint"
#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "comm_xfer_utils.hh"
#include "compute_metrics_3way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Compute 3-way numerators Czek cases that don't use linalg package.

void ComputeMetrics3WayBlock::compute_czek_(VData vdata_i, VData vdata_j,
  VData vdata_k, GMMetrics& numerators,
  int j_block, int k_block, int section_step) {

  GMMetrics* metrics = &numerators;
  CEnv* env = &env_;
  GMVectors* vectors_i = vdata_i.vectors;
  GMVectors* vectors_j = vdata_j.vectors;
  GMVectors* vectors_k = vdata_k.vectors;
  MirroredBuf* vectors_i_buf = vdata_i.buf;
  MirroredBuf* vectors_j_buf = vdata_j.buf;
  MirroredBuf* vectors_k_buf = vdata_k.buf;
  VectorSums* vector_sums_i = vdata_i.sums;
  VectorSums* vector_sums_j = vdata_j.sums;
  VectorSums* vector_sums_k = vdata_k.sums;

  COMET_INSIST(metrics && env);
  COMET_INSIST(vectors_i && vectors_j && vectors_k);
  COMET_INSIST(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env->num_block_vector());
  COMET_INSIST(! (env->proc_num_vector() == j_block &&
              env->proc_num_vector() != k_block));
  COMET_INSIST(! (env->proc_num_vector() == k_block &&
              env->proc_num_vector() != j_block));
  COMET_INSIST(!env->is_compute_method_gpu());
  COMET_INSIST(env->num_way() == NumWay::_3);
  COMET_INSIST(vector_sums_i && vector_sums_j && vector_sums_k);

  // Initializations.

  const int nvl = metrics->num_vector_local;
  const int nfl = vectors_i->num_field_local;

  const int i_block = env->proc_num_vector();

  GMSectionInfo si_value, *si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const VectorSums* const vs_i = vector_sums_i;
  const VectorSums* const vs_j = vector_sums_j;
  const VectorSums* const vs_k = vector_sums_k;

  // ----------------------------------
  if (!env->is_compute_method_gpu() && !env->all2all()) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() && "num_proc_field>1 "
      "for REF/CPU compute_method for this case not supported");

    COMET_INSIST(gm_num_section_steps(env, 1) == 1 &&
             "not all2all case always has one section step.");

    // No off-proc all2all: compute tetrahedron of values.

    for (int j = 0; j < nvl; ++j) {
      for (int k = j+1; k < nvl; ++k) {
        for (int i = 0; i < j; ++i) {
          // Make arithmetic order-independent.
          GMFloat smin, smid, smax;
          const GMFloat si = vs_i->sum(i);
          const GMFloat sj = vs_i->sum(j);
          const GMFloat sk = vs_i->sum(k);
          utils::sort_3(smin, smid, smax, si, sj, sk);
          const GMFloat denom = smin + smid + smax;
          GMFloat numer = 0;
          for (int f = 0; f < nfl; ++f) {
            const GMFloat val1 = GMVectors_float_get(vectors_i, f, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_i, f, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_i, f, k, env);
            GMFloat min12 = val1 < val2 ? val1 : val2;
            numer += min12;
            numer += val1 < val3 ? val1 : val3;
            numer += val2 < val3 ? val2 : val3;
            numer -= min12 < val3 ? min12 : val3;
          } // for f
          const GMFloat value = ((GMFloat)3) * numer / (((GMFloat)2) * denom);

          Metrics_elt_3<GMFloat>(*metrics, i, j, k,
            env->proc_num_vector(), env->proc_num_vector(), *env) = value;
        }
        metrics->num_metrics_local_computed += j;
      }
    }

    // ----------------------------------
  } else if (!env->is_compute_method_gpu()) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() && "num_proc_field>1 "
      "for REF/CPU compute_method for this case not supported");

    // Compute tetrahedron, triang prism or block section.

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {

      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      MetricsIndexCache index_cache = {};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          const int i = si->unperm0(I, J, K);
          const int j = si->unperm1(I, J, K);
          const int k = si->unperm2(I, J, K);

          // Make arithmetic order-independent.
          GMFloat smin, smid, smax;
          const GMFloat si = vs_i->sum(i);
          const GMFloat sj = vs_j->sum(j);
          const GMFloat sk = vs_k->sum(k);
          utils::sort_3(smin, smid, smax, si, sj, sk);
          const GMFloat denom = smin + smid + smax;
          GMFloat numer = 0;
          for (int f = 0; f < nfl; ++f) {
            const GMFloat val1 = GMVectors_float_get(vectors_i, f, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_j, f, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_k, f, k, env);
            const GMFloat min_ij = val1 < val2 ? val1 : val2;
            const GMFloat min_ik = val1 < val3 ? val1 : val3;
            const GMFloat min_jk = val2 < val3 ? val2 : val3;
            const GMFloat min_ijk = min_ij < val3 ? min_ij : val3;
            numer += min_ij + min_ik + min_jk - min_ijk;
          } // for f

          const GMFloat value = ((GMFloat)3) * numer / (((GMFloat)2) * denom);

          Metrics_elt_3<GMFloat>(*metrics, I, J, K,
            j_block, k_block, index_cache, *env) = value;
        } //---I
        metrics->num_metrics_local_computed += I_max - I_min;
      } //---K
    } //---J

    // ----------------------------------
  } else /* if (env->is_compute_method_gpu()) */ {
    // ----------------------------------

    COMET_INSIST(false && "Invalid compute_method.");

  } // if GPU

  GMSectionInfo_destroy(si, env);
}

//-----------------------------------------------------------------------------
/// \brief Compute 3-way numerators CCC/DUO cases that don't use linalg package.

void ComputeMetrics3WayBlock::compute_ccc_duo_(VData vdata_i, VData vdata_j,
  VData vdata_k, GMMetrics& numerators,
  int j_block, int k_block, int section_step) {


  GMMetrics* metrics = &numerators;
  CEnv* env = &env_;
  GMVectors* vectors_i = vdata_i.vectors;
  GMVectors* vectors_j = vdata_j.vectors;
  GMVectors* vectors_k = vdata_k.vectors;
  MirroredBuf* vectors_i_buf = vdata_i.buf;
  MirroredBuf* vectors_j_buf = vdata_j.buf;
  MirroredBuf* vectors_k_buf = vdata_k.buf;
  VectorSums* vector_sums_i = vdata_i.sums;
  VectorSums* vector_sums_j = vdata_j.sums;
  VectorSums* vector_sums_k = vdata_k.sums;




  COMET_INSIST(metrics && env);
  COMET_INSIST(vectors_i && vectors_j && vectors_k);
  COMET_INSIST(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env->num_block_vector());
  COMET_INSIST(! (env->proc_num_vector() == j_block &&
              env->proc_num_vector() != k_block));
  COMET_INSIST(! (env->proc_num_vector() == k_block &&
              env->proc_num_vector() != j_block));
  COMET_INSIST(!env->is_using_linalg());
  COMET_INSIST(env->num_way() == NumWay::_3);
  COMET_INSIST(vector_sums_i && vector_sums_j && vector_sums_k);

  typedef MetricFormatTraits<MetricFormat::PACKED_DOUBLE> MFT;

  // Initializations.

  const int nvl = metrics->num_vector_local;

  const int i_block = env->proc_num_vector();

  GMSectionInfo si_value, *si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const VectorSums* const vs_i = vector_sums_i;
  const VectorSums* const vs_j = vector_sums_j;
  const VectorSums* const vs_k = vector_sums_k;

  // ----------------------------------
  if (env->compute_method() == ComputeMethod::REF) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for REF compute_method not supported");

    const int nfal = vectors_i->dm->num_field_active_local;

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {

      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      MetricsIndexCache index_cache = {};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          const int i = si->unperm0(I, J, K);
          const int j = si->unperm1(I, J, K);
          const int k = si->unperm2(I, J, K);

          GMTally4x2 sum = GMTally4x2_null();
          for (int f = 0; f < nfal; ++f) {
            const GMBits2 vi = GMVectors_bits2_get(vectors_i, f, i, env);
            const GMBits2 vj = GMVectors_bits2_get(vectors_j, f, j, env);
            const GMBits2 vk = GMVectors_bits2_get(vectors_k, f, k, env);

            const bool unknown_i = env->sparse() ? vi == GM_2BIT_UNKNOWN
                                               : false;
            const bool unknown_j = env->sparse() ? vj == GM_2BIT_UNKNOWN
                                               : false;
            const bool unknown_k = env->sparse() ? vk == GM_2BIT_UNKNOWN
                                               : false;

            if ( ! unknown_i && ! unknown_j  && ! unknown_k ) {

              if (env->metric_type() == MetricType::CCC) {

                /* clang-format off */
                const int r000 =
                  ((!(vi & 1)) && (!(vj & 1)) && (!(vk & 1))) +
                  ((!(vi & 1)) && (!(vj & 1)) && (!(vk & 2))) +
                  ((!(vi & 1)) && (!(vj & 2)) && (!(vk & 1))) +
                  ((!(vi & 1)) && (!(vj & 2)) && (!(vk & 2))) +
                  ((!(vi & 2)) && (!(vj & 1)) && (!(vk & 1))) +
                  ((!(vi & 2)) && (!(vj & 1)) && (!(vk & 2))) +
                  ((!(vi & 2)) && (!(vj & 2)) && (!(vk & 1))) +
                  ((!(vi & 2)) && (!(vj & 2)) && (!(vk & 2)));
                const int r001 =
                  ((!(vi & 1)) && (!(vj & 1)) && ( (vk & 1))) +
                  ((!(vi & 1)) && (!(vj & 1)) && ( (vk & 2))) +
                  ((!(vi & 1)) && (!(vj & 2)) && ( (vk & 1))) +
                  ((!(vi & 1)) && (!(vj & 2)) && ( (vk & 2))) +
                  ((!(vi & 2)) && (!(vj & 1)) && ( (vk & 1))) +
                  ((!(vi & 2)) && (!(vj & 1)) && ( (vk & 2))) +
                  ((!(vi & 2)) && (!(vj & 2)) && ( (vk & 1))) +
                  ((!(vi & 2)) && (!(vj & 2)) && ( (vk & 2)));
                const int r010 =
                  ((!(vi & 1)) && ( (vj & 1)) && (!(vk & 1))) +
                  ((!(vi & 1)) && ( (vj & 1)) && (!(vk & 2))) +
                  ((!(vi & 1)) && ( (vj & 2)) && (!(vk & 1))) +
                  ((!(vi & 1)) && ( (vj & 2)) && (!(vk & 2))) +
                  ((!(vi & 2)) && ( (vj & 1)) && (!(vk & 1))) +
                  ((!(vi & 2)) && ( (vj & 1)) && (!(vk & 2))) +
                  ((!(vi & 2)) && ( (vj & 2)) && (!(vk & 1))) +
                  ((!(vi & 2)) && ( (vj & 2)) && (!(vk & 2)));
                const int r011 =
                  ((!(vi & 1)) && ( (vj & 1)) && ( (vk & 1))) +
                  ((!(vi & 1)) && ( (vj & 1)) && ( (vk & 2))) +
                  ((!(vi & 1)) && ( (vj & 2)) && ( (vk & 1))) +
                  ((!(vi & 1)) && ( (vj & 2)) && ( (vk & 2))) +
                  ((!(vi & 2)) && ( (vj & 1)) && ( (vk & 1))) +
                  ((!(vi & 2)) && ( (vj & 1)) && ( (vk & 2))) +
                  ((!(vi & 2)) && ( (vj & 2)) && ( (vk & 1))) +
                  ((!(vi & 2)) && ( (vj & 2)) && ( (vk & 2)));
                const int r100 =
                  (( (vi & 1)) && (!(vj & 1)) && (!(vk & 1))) +
                  (( (vi & 1)) && (!(vj & 1)) && (!(vk & 2))) +
                  (( (vi & 1)) && (!(vj & 2)) && (!(vk & 1))) +
                  (( (vi & 1)) && (!(vj & 2)) && (!(vk & 2))) +
                  (( (vi & 2)) && (!(vj & 1)) && (!(vk & 1))) +
                  (( (vi & 2)) && (!(vj & 1)) && (!(vk & 2))) +
                  (( (vi & 2)) && (!(vj & 2)) && (!(vk & 1))) +
                  (( (vi & 2)) && (!(vj & 2)) && (!(vk & 2)));
                const int r101 =
                  (( (vi & 1)) && (!(vj & 1)) && ( (vk & 1))) +
                  (( (vi & 1)) && (!(vj & 1)) && ( (vk & 2))) +
                  (( (vi & 1)) && (!(vj & 2)) && ( (vk & 1))) +
                  (( (vi & 1)) && (!(vj & 2)) && ( (vk & 2))) +
                  (( (vi & 2)) && (!(vj & 1)) && ( (vk & 1))) +
                  (( (vi & 2)) && (!(vj & 1)) && ( (vk & 2))) +
                  (( (vi & 2)) && (!(vj & 2)) && ( (vk & 1))) +
                  (( (vi & 2)) && (!(vj & 2)) && ( (vk & 2)));
                const int r110 =
                  (( (vi & 1)) && ( (vj & 1)) && (!(vk & 1))) +
                  (( (vi & 1)) && ( (vj & 1)) && (!(vk & 2))) +
                  (( (vi & 1)) && ( (vj & 2)) && (!(vk & 1))) +
                  (( (vi & 1)) && ( (vj & 2)) && (!(vk & 2))) +
                  (( (vi & 2)) && ( (vj & 1)) && (!(vk & 1))) +
                  (( (vi & 2)) && ( (vj & 1)) && (!(vk & 2))) +
                  (( (vi & 2)) && ( (vj & 2)) && (!(vk & 1))) +
                  (( (vi & 2)) && ( (vj & 2)) && (!(vk & 2)));
                const int r111 =
                  (( (vi & 1)) && ( (vj & 1)) && ( (vk & 1))) +
                  (( (vi & 1)) && ( (vj & 1)) && ( (vk & 2))) +
                  (( (vi & 1)) && ( (vj & 2)) && ( (vk & 1))) +
                  (( (vi & 1)) && ( (vj & 2)) && ( (vk & 2))) +
                  (( (vi & 2)) && ( (vj & 1)) && ( (vk & 1))) +
                  (( (vi & 2)) && ( (vj & 1)) && ( (vk & 2))) +
                  (( (vi & 2)) && ( (vj & 2)) && ( (vk & 1))) +
                  (( (vi & 2)) && ( (vj & 2)) && ( (vk & 2)));
                /* clang-format on */

                // Accumulate

                MFT::add(sum.data[0], r000, r001);
                MFT::add(sum.data[1], r010, r011);
                MFT::add(sum.data[2], r100, r101);
                MFT::add(sum.data[3], r110, r111);

              } else { // (env->metric_type() == MetricType::DUO)

                /* clang-format off */
                const int r000 = ((!(vi & 1)) && (!(vj & 1)) && (!(vk & 1)));
                const int r001 = ((!(vi & 1)) && (!(vj & 1)) && ( (vk & 1)));
                const int r010 = ((!(vi & 1)) && ( (vj & 1)) && (!(vk & 1)));
                const int r011 = ((!(vi & 1)) && ( (vj & 1)) && ( (vk & 1)));
                const int r100 = (( (vi & 1)) && (!(vj & 1)) && (!(vk & 1)));
                const int r101 = (( (vi & 1)) && (!(vj & 1)) && ( (vk & 1)));
                const int r110 = (( (vi & 1)) && ( (vj & 1)) && (!(vk & 1)));
                const int r111 = (( (vi & 1)) && ( (vj & 1)) && ( (vk & 1)));
                /* clang-format on */

                // Accumulate

                MFT::add(sum.data[0], r000, r001);
                MFT::add(sum.data[1], r010, r011);
                MFT::add(sum.data[2], r100, r101);
                MFT::add(sum.data[3], r110, r111);

              } // if (env->metric_type() == MetricType::CCC)

            } // if ! unknown
          } // for f

          // Get denom

          const auto si1 = (GMTally1)vs_i->sum(i);
          const auto sj1 = (GMTally1)vs_j->sum(j);
          const auto sk1 = (GMTally1)vs_k->sum(k);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);

          const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
          const int k_block_eff = env->all2all() ? k_block : env->proc_num_vector();

          Metrics_elt_3<GMTally4x2>(*metrics, I, J, K,
            j_block_eff, k_block_eff, index_cache, *env) = sum;
          Metrics_elt_3<GMFloat3, MetricsArray::S>(*metrics, I, J, K,
            j_block_eff, k_block_eff, index_cache, *env) = si1_sj1_sk1;
          if (env->sparse()) {
            const auto ci1 = (GMTally1)vs_i->count(i);
            const auto cj1 = (GMTally1)vs_j->count(j);
            const auto ck1 = (GMTally1)vs_k->count(k);
            const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
            Metrics_elt_3<GMFloat3, MetricsArray::C>(*metrics, I, J, K,
              j_block_eff, k_block_eff, index_cache, *env) = ci1_cj1_ck1;
          } // if sparse

        } // for I
        metrics->num_metrics_local_computed += I_max - I_min;
      } // for K
    } // for J

    // ----------------------------------
  } else if (env->compute_method() == ComputeMethod::CPU &&
             !env->is_using_linalg()) {
    // ----------------------------------

    COMET_INSIST_INTERFACE(env, ! env->do_reduce() &&
      "num_proc_field>1 for CPU compute_method for this case not supported");

    /* clang-format off */

    const int cbpe = env->counted_bits_per_elt();

    const int pad_adjustment = cbpe * cbpe * cbpe * metrics->dm->num_pad_field_local;

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {

      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      MetricsIndexCache index_cache = {};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          const int i = si->unperm0(I, J, K);
          const int j = si->unperm1(I, J, K);
          const int k = si->unperm2(I, J, K);

          GMTally4x2 sum = GMTally4x2_null();
          const int npvfl = vectors_i->num_packedval_field_local;
          for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

            // Extract input values to process.

            const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_i, pvfl, i,
                                                         env);
            const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_j, pvfl, j,
                                                          env);
            const GMBits2x64 vk = GMVectors_bits2x64_get(vectors_k, pvfl, k,
                                                          env);

            const uint64_t vi0 = vi.data[0];
            const uint64_t vi1 = vi.data[1];
            const uint64_t vj0 = vj.data[0];
            const uint64_t vj1 = vj.data[1];
            const uint64_t vk0 = vk.data[0];
            const uint64_t vk1 = vk.data[1];

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

            const uint64_t vk0mask =
                       (env->sparse() ? (vk0 | ~(vk0 >> 1)) & oddbits : oddbits);

            const uint64_t vk1mask =
                       (env->sparse() ? (vk1 | ~(vk1 >> 1)) & oddbits : oddbits);

            const uint64_t v0mask = vi0mask & vj0mask & vk0mask;
            const uint64_t v1mask = vi1mask & vj1mask & vk1mask;

            // Get even/odd bits for each seminibble, masked to active.

            const uint64_t vi0_0 =  vi0       & v0mask;
            const uint64_t vi0_1 = (vi0 >> 1) & v0mask;
            const uint64_t vi1_0 =  vi1       & v1mask;
            const uint64_t vi1_1 = (vi1 >> 1) & v1mask;
            const uint64_t vj0_0 =  vj0       & v0mask;
            const uint64_t vj0_1 = (vj0 >> 1) & v0mask;
            const uint64_t vj1_0 =  vj1       & v1mask;
            const uint64_t vj1_1 = (vj1 >> 1) & v1mask;
            const uint64_t vk0_0 =  vk0       & v0mask;
            const uint64_t vk0_1 = (vk0 >> 1) & v0mask;
            const uint64_t vk1_0 =  vk1       & v1mask;
            const uint64_t vk1_1 = (vk1 >> 1) & v1mask;

            // Get complements of even/odd bits for each seminibble; mask.

            const uint64_t nvi0_0 = ~ vi0       & v0mask;
            const uint64_t nvi0_1 = ~(vi0 >> 1) & v0mask;
            const uint64_t nvi1_0 = ~ vi1       & v1mask;
            const uint64_t nvi1_1 = ~(vi1 >> 1) & v1mask;
            const uint64_t nvj0_0 = ~ vj0       & v0mask;
            const uint64_t nvj0_1 = ~(vj0 >> 1) & v0mask;
            const uint64_t nvj1_0 = ~ vj1       & v1mask;
            const uint64_t nvj1_1 = ~(vj1 >> 1) & v1mask;
            const uint64_t nvk0_0 = ~ vk0       & v0mask;
            const uint64_t nvk0_1 = ~(vk0 >> 1) & v0mask;
            const uint64_t nvk1_0 = ~ vk1       & v1mask;
            const uint64_t nvk1_1 = ~(vk1 >> 1) & v1mask;

            if (env->metric_type() == MetricType::CCC) {

              const int r000 = utils::popc64((nvi0_0 & nvj0_0 & nvk0_0) |
                                           ( (nvi0_0 & nvj0_0 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_0 & nvj0_1 & nvk0_0) |
                                           ( (nvi0_0 & nvj0_1 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_1 & nvj0_0 & nvk0_0) |
                                           ( (nvi0_1 & nvj0_0 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_1 & nvj0_1 & nvk0_0) |
                                           ( (nvi0_1 & nvj0_1 & nvk0_1) << 1 ))+
                               utils::popc64((nvi1_0 & nvj1_0 & nvk1_0) |
                                           ( (nvi1_0 & nvj1_0 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_0 & nvj1_1 & nvk1_0) |
                                           ( (nvi1_0 & nvj1_1 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_1 & nvj1_0 & nvk1_0) |
                                           ( (nvi1_1 & nvj1_0 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_1 & nvj1_1 & nvk1_0) |
                                           ( (nvi1_1 & nvj1_1 & nvk1_1) << 1 ));
              const int r001 = utils::popc64((nvi0_0 & nvj0_0 &  vk0_0) |
                                           ( (nvi0_0 & nvj0_0 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_0 & nvj0_1 &  vk0_0) |
                                           ( (nvi0_0 & nvj0_1 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_1 & nvj0_0 &  vk0_0) |
                                           ( (nvi0_1 & nvj0_0 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_1 & nvj0_1 &  vk0_0) |
                                           ( (nvi0_1 & nvj0_1 &  vk0_1) << 1 ))+
                               utils::popc64((nvi1_0 & nvj1_0 &  vk1_0) |
                                           ( (nvi1_0 & nvj1_0 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_0 & nvj1_1 &  vk1_0) |
                                           ( (nvi1_0 & nvj1_1 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_1 & nvj1_0 &  vk1_0) |
                                           ( (nvi1_1 & nvj1_0 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_1 & nvj1_1 &  vk1_0) |
                                           ( (nvi1_1 & nvj1_1 &  vk1_1) << 1 ));
              const int r010 = utils::popc64((nvi0_0 &  vj0_0 & nvk0_0) |
                                           ( (nvi0_0 &  vj0_0 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_0 &  vj0_1 & nvk0_0) |
                                           ( (nvi0_0 &  vj0_1 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_1 &  vj0_0 & nvk0_0) |
                                           ( (nvi0_1 &  vj0_0 & nvk0_1) << 1 ))+
                               utils::popc64((nvi0_1 &  vj0_1 & nvk0_0) |
                                           ( (nvi0_1 &  vj0_1 & nvk0_1) << 1 ))+
                               utils::popc64((nvi1_0 &  vj1_0 & nvk1_0) |
                                           ( (nvi1_0 &  vj1_0 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_0 &  vj1_1 & nvk1_0) |
                                           ( (nvi1_0 &  vj1_1 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_1 &  vj1_0 & nvk1_0) |
                                           ( (nvi1_1 &  vj1_0 & nvk1_1) << 1 ))+
                               utils::popc64((nvi1_1 &  vj1_1 & nvk1_0) |
                                           ( (nvi1_1 &  vj1_1 & nvk1_1) << 1 ));
              const int r011 = utils::popc64((nvi0_0 &  vj0_0 &  vk0_0) |
                                           ( (nvi0_0 &  vj0_0 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_0 &  vj0_1 &  vk0_0) |
                                           ( (nvi0_0 &  vj0_1 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_1 &  vj0_0 &  vk0_0) |
                                           ( (nvi0_1 &  vj0_0 &  vk0_1) << 1 ))+
                               utils::popc64((nvi0_1 &  vj0_1 &  vk0_0) |
                                           ( (nvi0_1 &  vj0_1 &  vk0_1) << 1 ))+
                               utils::popc64((nvi1_0 &  vj1_0 &  vk1_0) |
                                           ( (nvi1_0 &  vj1_0 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_0 &  vj1_1 &  vk1_0) |
                                           ( (nvi1_0 &  vj1_1 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_1 &  vj1_0 &  vk1_0) |
                                           ( (nvi1_1 &  vj1_0 &  vk1_1) << 1 ))+
                               utils::popc64((nvi1_1 &  vj1_1 &  vk1_0) |
                                           ( (nvi1_1 &  vj1_1 &  vk1_1) << 1 ));
              const int r100 = utils::popc64(( vi0_0 & nvj0_0 & nvk0_0) |
                                           ( ( vi0_0 & nvj0_0 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_0 & nvj0_1 & nvk0_0) |
                                           ( ( vi0_0 & nvj0_1 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_1 & nvj0_0 & nvk0_0) |
                                           ( ( vi0_1 & nvj0_0 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_1 & nvj0_1 & nvk0_0) |
                                           ( ( vi0_1 & nvj0_1 & nvk0_1) << 1 ))+
                               utils::popc64(( vi1_0 & nvj1_0 & nvk1_0) |
                                           ( ( vi1_0 & nvj1_0 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_0 & nvj1_1 & nvk1_0) |
                                           ( ( vi1_0 & nvj1_1 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_1 & nvj1_0 & nvk1_0) |
                                           ( ( vi1_1 & nvj1_0 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_1 & nvj1_1 & nvk1_0) |
                                           ( ( vi1_1 & nvj1_1 & nvk1_1) << 1 ));
              const int r101 = utils::popc64(( vi0_0 & nvj0_0 &  vk0_0) |
                                           ( ( vi0_0 & nvj0_0 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_0 & nvj0_1 &  vk0_0) |
                                           ( ( vi0_0 & nvj0_1 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_1 & nvj0_0 &  vk0_0) |
                                           ( ( vi0_1 & nvj0_0 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_1 & nvj0_1 &  vk0_0) |
                                           ( ( vi0_1 & nvj0_1 &  vk0_1) << 1 ))+
                               utils::popc64(( vi1_0 & nvj1_0 &  vk1_0) |
                                           ( ( vi1_0 & nvj1_0 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_0 & nvj1_1 &  vk1_0) |
                                           ( ( vi1_0 & nvj1_1 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_1 & nvj1_0 &  vk1_0) |
                                           ( ( vi1_1 & nvj1_0 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_1 & nvj1_1 &  vk1_0) |
                                           ( ( vi1_1 & nvj1_1 &  vk1_1) << 1 ));
              const int r110 = utils::popc64(( vi0_0 &  vj0_0 & nvk0_0) |
                                           ( ( vi0_0 &  vj0_0 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_0 &  vj0_1 & nvk0_0) |
                                           ( ( vi0_0 &  vj0_1 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_1 &  vj0_0 & nvk0_0) |
                                           ( ( vi0_1 &  vj0_0 & nvk0_1) << 1 ))+
                               utils::popc64(( vi0_1 &  vj0_1 & nvk0_0) |
                                           ( ( vi0_1 &  vj0_1 & nvk0_1) << 1 ))+
                               utils::popc64(( vi1_0 &  vj1_0 & nvk1_0) |
                                           ( ( vi1_0 &  vj1_0 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_0 &  vj1_1 & nvk1_0) |
                                           ( ( vi1_0 &  vj1_1 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_1 &  vj1_0 & nvk1_0) |
                                           ( ( vi1_1 &  vj1_0 & nvk1_1) << 1 ))+
                               utils::popc64(( vi1_1 &  vj1_1 & nvk1_0) |
                                           ( ( vi1_1 &  vj1_1 & nvk1_1) << 1 ));
              const int r111 = utils::popc64(( vi0_0 &  vj0_0 &  vk0_0) |
                                           ( ( vi0_0 &  vj0_0 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_0 &  vj0_1 &  vk0_0) |
                                           ( ( vi0_0 &  vj0_1 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_1 &  vj0_0 &  vk0_0) |
                                           ( ( vi0_1 &  vj0_0 &  vk0_1) << 1 ))+
                               utils::popc64(( vi0_1 &  vj0_1 &  vk0_0) |
                                           ( ( vi0_1 &  vj0_1 &  vk0_1) << 1 ))+
                               utils::popc64(( vi1_0 &  vj1_0 &  vk1_0) |
                                           ( ( vi1_0 &  vj1_0 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_0 &  vj1_1 &  vk1_0) |
                                           ( ( vi1_0 &  vj1_1 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_1 &  vj1_0 &  vk1_0) |
                                           ( ( vi1_1 &  vj1_0 &  vk1_1) << 1 ))+
                               utils::popc64(( vi1_1 &  vj1_1 &  vk1_0) |
                                           ( ( vi1_1 &  vj1_1 &  vk1_1) << 1 ));

              // Accumulate

              MFT::add(sum.data[0], r000, r001);
              MFT::add(sum.data[1], r010, r011);
              MFT::add(sum.data[2], r100, r101);
              MFT::add(sum.data[3], r110, r111);

            } else { // (env->metric_type() == MetricType::DUO)

              const int r000 = utils::popc64((nvi0_0 & nvj0_0 & nvk0_0) |
                                           ( (nvi1_0 & nvj1_0 & nvk1_0) << 1 ));
              const int r001 = utils::popc64((nvi0_0 & nvj0_0 &  vk0_0) |
                                           ( (nvi1_0 & nvj1_0 &  vk1_0) << 1 ));
              const int r010 = utils::popc64((nvi0_0 &  vj0_0 & nvk0_0) |
                                           ( (nvi1_0 &  vj1_0 & nvk1_0) << 1 ));
              const int r011 = utils::popc64((nvi0_0 &  vj0_0 &  vk0_0) |
                                           ( (nvi1_0 &  vj1_0 &  vk1_0) << 1 ));
              const int r100 = utils::popc64(( vi0_0 & nvj0_0 & nvk0_0) |
                                           ( ( vi1_0 & nvj1_0 & nvk1_0) << 1 ));
              const int r101 = utils::popc64(( vi0_0 & nvj0_0 &  vk0_0) |
                                           ( ( vi1_0 & nvj1_0 &  vk1_0) << 1 ));
              const int r110 = utils::popc64(( vi0_0 &  vj0_0 & nvk0_0) |
                                           ( ( vi1_0 &  vj1_0 & nvk1_0) << 1 ));
              const int r111 = utils::popc64(( vi0_0 &  vj0_0 &  vk0_0) |
                                           ( ( vi1_0 &  vj1_0 &  vk1_0) << 1 ));

              // Accumulate

              MFT::add(sum.data[0], r000, r001);
              MFT::add(sum.data[1], r010, r011);
              MFT::add(sum.data[2], r100, r101);
              MFT::add(sum.data[3], r110, r111);

            } // if (env->metric_type() == MetricType::CCC)

          } // for pvfl

          // Adjust for pad

#ifdef COMET_ASSERTIONS_ON
          GMTally4x2 sum_old = sum;
#endif
          MFT::subtract(sum.data[0], pad_adjustment, 0);
#ifdef COMET_ASSERTIONS_ON
          COMET_ASSERT(GMTally4x2_get(sum_old, 0, 0, 0) ==
                   GMTally4x2_get(sum, 0, 0, 0) + pad_adjustment);
#endif

          // Get denom

          const auto si1 = (GMTally1)vs_i->sum(i);
          const auto sj1 = (GMTally1)vs_j->sum(j);
          const auto sk1 = (GMTally1)vs_k->sum(k);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);

          const int j_block_eff = env->all2all() ? j_block : env->proc_num_vector();
          const int k_block_eff = env->all2all() ? k_block : env->proc_num_vector();

          Metrics_elt_3<GMTally4x2>(*metrics, I, J, K,
            j_block_eff, k_block_eff, index_cache, *env) = sum;
          Metrics_elt_3<GMFloat3, MetricsArray::S>(*metrics, I, J, K,
            j_block_eff, k_block_eff, index_cache, *env) = si1_sj1_sk1;
          if (env->sparse()) {
            const auto ci1 = (GMTally1)vs_i->count(i);
            const auto cj1 = (GMTally1)vs_j->count(j);
            const auto ck1 = (GMTally1)vs_k->count(k);
            const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
            Metrics_elt_3<GMFloat3, MetricsArray::C>(*metrics, I, J, K,
              j_block_eff, k_block_eff, index_cache, *env) = ci1_cj1_ck1;
          } // if sparse

        } //---I
        metrics->num_metrics_local_computed += I_max - I_min;
      } //---K
    } //---J
    /* clang-format on */

    // ----------------------------------
  } else /* if (env->is_using_linalg()) */ {
    // ----------------------------------
    COMET_INSIST(false && "Invalid compute_method");
    // ----------------------------------
  } // if
  // ----------------------------------

  GMSectionInfo_destroy(si, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
