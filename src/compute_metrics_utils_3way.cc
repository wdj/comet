//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_utils_3way.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, 3-way.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "string.h"

#include "env.hh"
#include "mirrored_buf.hh"
#include "linalg.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_metrics_utils.hh"
#include "compute_metrics_utils_3way.hh"

//=============================================================================
/*---Start calculation of numerators, 3-way Czekanowski non-gpu---*/

void gm_compute_czek_numerators_3way_nongpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_block,
    int k_block,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env) {
  GMAssertAlways(this_ && metrics && env);
  GMAssertAlways(vectors_i && vectors_j && vectors_k);
  GMAssertAlways(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == j_block &&
                   GMEnv_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == k_block &&
                   GMEnv_proc_num_vector_i(env) != j_block));
  GMAssertAlways(GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssertAlways(vector_sums_i && vector_sums_j && vector_sums_k);

  /*---Initializations---*/

  const int nvl = metrics->num_vector_local;
  const int nfl = vectors_i->num_field_local;

  const int i_block = GMEnv_proc_num_vector_i(env);

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const GMVectorSums* const vs_i = vector_sums_i;
  const GMVectorSums* const vs_j = vector_sums_j;
  const GMVectorSums* const vs_k = vector_sums_k;

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU &&
      !GMEnv_all2all(env)) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for REF/CPU cases not supported"
                      : 0);

    GMAssertAlways(gm_num_section_steps(env, 1) == 1);

    /*---No off-proc all2all: compute tetrahedron of values---*/

    //const int section_num = 0;
    //const int J_lo = gm_J_lo(section_num, nvl, 1, env);
    //const int J_hi = gm_J_hi(section_num, nvl, 1, env);
    //const int j_min = J_lo;
    //const int j_max = J_hi;
    //for (int j = j_min; j < j_max; ++j) {
    for (int j = 0; j < nvl; ++j) {
      for (int k = j+1; k < nvl; ++k) {
        for (int i = 0; i < j; ++i) {
          /*---Make arithmetic order-independent---*/
          GMFloat smin, smid, smax;
          const GMFloat si = GMVectorSums_sum(vs_i, i, env);
          const GMFloat sj = GMVectorSums_sum(vs_i, j, env);
          const GMFloat sk = GMVectorSums_sum(vs_i, k, env);
          GMFloat_sort_3(&smin, &smid, &smax, &si, &sj, &sk);
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
          } /*---for f---*/
          const GMFloat value = ((GMFloat)3) * numer / (((GMFloat)2) * denom);
          GMMetrics_float_set_3(metrics, i, j, k, value, env);
        }
        metrics->num_elts_local_computed += j;
      }
    }

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for REF/CPU cases not supported"
                      : 0);

    /*---Compute tetrahedron, triang prism or block section---*/

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {


      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      GMIndexCache index_cache = {0};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          /* clang-format off */
          const int i = !si->is_part3 ?   I :
                             si->sax0 ?   J :
                             si->sax1 ?   I :
                          /* si->sax2 ?*/ K;
          const int j = !si->is_part3 ?   J :
                             si->sax0 ?   K :
                             si->sax1 ?   J :
                          /* si->sax2 ?*/ I;
          const int k = !si->is_part3 ?   K :
                             si->sax0 ?   I :
                             si->sax1 ?   K :
                          /* si->sax2 ?*/ J;
          /* clang-format on */

          /*---Make arithmetic order-independent---*/
          GMFloat smin, smid, smax;
          const GMFloat si = GMVectorSums_sum(vs_i, i, env);
          const GMFloat sj = GMVectorSums_sum(vs_j, j, env);
          const GMFloat sk = GMVectorSums_sum(vs_k, k, env);
          GMFloat_sort_3(&smin, &smid, &smax, &si, &sj, &sk);
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
          } /*---for f---*/

          const GMFloat value = ((GMFloat)3) * numer / (((GMFloat)2) * denom);

          GMMetrics_float_set_all2all_3_permuted_cache(metrics, I, J, K,
                                   j_block, k_block, value, &index_cache, env);
        } //---I
        metrics->num_elts_local_computed += I_max - I_min;
      } //---K
    } //---J

    /*----------------------------------------*/
  } else /* if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    GMAssertAlways(false
                 ? "logic error - code branch should never be executed"
                 : 0);

  } /*---if GPU---*/

  GMSectionInfo_destroy(si, env);;
}

//=============================================================================
/*---Start calculation of numerators, 3-way CCC non-gpu---*/

void gm_compute_ccc_numerators_3way_nongpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_block,
    int k_block,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env) {
  GMAssertAlways(this_ && metrics && env);
  GMAssertAlways(vectors_i && vectors_j && vectors_k);
  GMAssertAlways(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == j_block &&
                   GMEnv_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == k_block &&
                   GMEnv_proc_num_vector_i(env) != j_block));
  GMAssertAlways(GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssertAlways(vector_sums_i && vector_sums_j && vector_sums_k);

  /*---Initializations---*/

  const int nvl = metrics->num_vector_local;
  const int nfl = vectors_i->num_field_local;

  const int i_block = GMEnv_proc_num_vector_i(env);

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const GMVectorSums* const vs_i = vector_sums_i;
  const GMVectorSums* const vs_j = vector_sums_j;
  const GMVectorSums* const vs_k = vector_sums_k;

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_REF) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {

      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      GMIndexCache index_cache = {0};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          /* clang-format off */
          const int i = !si->is_part3 ?   I :
                             si->sax0 ?   J :
                             si->sax1 ?   I :
                          /* si->sax2 ?*/ K;
          const int j = !si->is_part3 ?   J :
                             si->sax0 ?   K :
                             si->sax1 ?   J :
                          /* si->sax2 ?*/ I;
          const int k = !si->is_part3 ?   K :
                             si->sax0 ?   I :
                             si->sax1 ?   K :
                          /* si->sax2 ?*/ J;
          /* clang-format on */

          GMTally4x2 sum = GMTally4x2_null();
          for (int f = 0; f < nfl; ++f) {
            const GMBits2 vi = GMVectors_bits2_get(vectors_i, f, i, env);
            const GMBits2 vj = GMVectors_bits2_get(vectors_j, f, j, env);
            const GMBits2 vk = GMVectors_bits2_get(vectors_k, f, k, env);

            const bool unknown_i = env->sparse ? vi == GM_2BIT_UNKNOWN
                                               : false;
            const bool unknown_j = env->sparse ? vj == GM_2BIT_UNKNOWN
                                               : false;
            const bool unknown_k = env->sparse ? vk == GM_2BIT_UNKNOWN
                                               : false;

            if ((!unknown_i) && (!unknown_j) && (!unknown_k)) {

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
              sum.data[0] += GMTally1_encode(r000, r001);
              sum.data[1] += GMTally1_encode(r010, r011);
              sum.data[2] += GMTally1_encode(r100, r101);
              sum.data[3] += GMTally1_encode(r110, r111);
            } /*---if !unknown---*/
          } /*---for f---*/
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_i, i, env);
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_j, j, env);
          const GMTally1 sk1 = (GMTally1)GMVectorSums_sum(vs_k, k, env);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);
          if (GMEnv_all2all(env)) {
            GMMetrics_tally4x2_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, sum, &index_cache, env);
            GMMetrics_float3_S_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, si1_sj1_sk1, &index_cache, env);
            if (env->sparse) {
              const GMTally1 ci1 = (GMTally1)GMVectorSums_count(vs_i, i, env);
              const GMTally1 cj1 = (GMTally1)GMVectorSums_count(vs_j, j, env);
              const GMTally1 ck1 = (GMTally1)GMVectorSums_count(vs_k, k, env);
              const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
              GMMetrics_float3_C_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, ci1_cj1_ck1, &index_cache, env);
            } /*---if sparse---*/
          } else /*---! all2all---*/ {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
            GMMetrics_float3_S_set_3(metrics, i, j, k, si1_sj1_sk1, env);
            if (env->sparse) {
              const GMTally1 ci1 = (GMTally1)GMVectorSums_count(vs_i, i, env);
              const GMTally1 cj1 = (GMTally1)GMVectorSums_count(vs_j, j, env);
              const GMTally1 ck1 = (GMTally1)GMVectorSums_count(vs_k, k, env);
              const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
              GMMetrics_float3_C_set_3(metrics, i, j, k, ci1_cj1_ck1, env);
            } /*---if sparse---*/
          } /*---if all2all---*/
        } /*---for I---*/
        metrics->num_elts_local_computed += I_max - I_min;
      } /*---for K---*/
    } /*---for J---*/

    /*----------------------------------------*/
  } else if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_CPU) {
    /*----------------------------------------*/

    GMInsist(env, !env->do_reduce
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Precompute masks for final (incomplete) packedval_field -
         can be 1 to 64 inclusive---*/

    /* clang-format off */

    const int nfl = vectors_i->num_field_local;
    const int num_field_active_local =
      GMEnv_proc_num_field(env) == GMEnv_num_proc_field(env) - 1
      ? nfl - (vectors_i->num_field - vectors_i->num_field_active) : nfl;
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

    const int J_lo = si->J_lb;
    const int J_hi = si->J_ub;
    for (int J = J_lo; J < J_hi; ++J) {

      const int I_min = 0;
      const int I_max = si->is_part1 ? J : nvl;
      const int K_min = si->is_part3 ? 0 : J + 1;
      const int K_max = nvl;

      GMIndexCache index_cache = {0};
      for (int K=K_min; K<K_max; ++K) {
        for (int I=I_min; I<I_max; ++I) {

          /* clang-format off */
          const int i = !si->is_part3 ?   I :
                             si->sax0 ?   J :
                             si->sax1 ?   I :
                          /* si->sax2 ?*/ K;
          const int j = !si->is_part3 ?   J :
                             si->sax0 ?   K :
                             si->sax1 ?   J :
                          /* si->sax2 ?*/ I;
          const int k = !si->is_part3 ?   K :
                             si->sax0 ?   I :
                             si->sax1 ?   K :
                          /* si->sax2 ?*/ J;
          /* clang-format on */

          GMTally4x2 sum = GMTally4x2_null();
          const int npvfl = vectors_i->num_packedval_field_local;
          const int pvfl_edge = num_packedval_field_active_local - 1;
          for (int pvfl = 0; pvfl < npvfl; ++pvfl) {
            /*---Get masks for active seminibbles in each word---*/

            const GMUInt64 activebits0 = pvfl < pvfl_edge ? allbits :
                                         pvfl == pvfl_edge ? edgemask0 : nobits;
            const GMUInt64 activebits1 = pvfl < pvfl_edge ? allbits :
                                         pvfl == pvfl_edge ? edgemask1 : nobits;

            /*---Extract input values to process---*/

            const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_i, pvfl, i,
                                                         env);
            const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_j, pvfl, j,
                                                          env);
            const GMBits2x64 vk = GMVectors_bits2x64_get(vectors_k, pvfl, k,
                                                          env);

            const GMUInt64 vi0 = vi.data[0];
            const GMUInt64 vi1 = vi.data[1];
            const GMUInt64 vj0 = vj.data[0];
            const GMUInt64 vj1 = vj.data[1];
            const GMUInt64 vk0 = vk.data[0];
            const GMUInt64 vk1 = vk.data[1];

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

            const GMUInt64 vk0mask = activebits0 &
                       (env->sparse ? (vk0 | ~(vk0 >> 1)) & oddbits : oddbits);

            const GMUInt64 vk1mask = activebits1 &
                       (env->sparse ? (vk1 | ~(vk1 >> 1)) & oddbits : oddbits);

            const GMUInt64 v0mask = vi0mask & vj0mask & vk0mask;
            const GMUInt64 v1mask = vi1mask & vj1mask & vk1mask;

            /*---Get even/odd bits for each seminibble, masked to active---*/

            const GMUInt64 vi0_0 =  vi0       & v0mask;
            const GMUInt64 vi0_1 = (vi0 >> 1) & v0mask;
            const GMUInt64 vi1_0 =  vi1       & v1mask;
            const GMUInt64 vi1_1 = (vi1 >> 1) & v1mask;
            const GMUInt64 vj0_0 =  vj0       & v0mask;
            const GMUInt64 vj0_1 = (vj0 >> 1) & v0mask;
            const GMUInt64 vj1_0 =  vj1       & v1mask;
            const GMUInt64 vj1_1 = (vj1 >> 1) & v1mask;
            const GMUInt64 vk0_0 =  vk0       & v0mask;
            const GMUInt64 vk0_1 = (vk0 >> 1) & v0mask;
            const GMUInt64 vk1_0 =  vk1       & v1mask;
            const GMUInt64 vk1_1 = (vk1 >> 1) & v1mask;

            /*---Get complements of even/odd bits for each seminibble; mask---*/

            const GMUInt64 nvi0_0 = ~ vi0       & v0mask;
            const GMUInt64 nvi0_1 = ~(vi0 >> 1) & v0mask;
            const GMUInt64 nvi1_0 = ~ vi1       & v1mask;
            const GMUInt64 nvi1_1 = ~(vi1 >> 1) & v1mask;
            const GMUInt64 nvj0_0 = ~ vj0       & v0mask;
            const GMUInt64 nvj0_1 = ~(vj0 >> 1) & v0mask;
            const GMUInt64 nvj1_0 = ~ vj1       & v1mask;
            const GMUInt64 nvj1_1 = ~(vj1 >> 1) & v1mask;
            const GMUInt64 nvk0_0 = ~ vk0       & v0mask;
            const GMUInt64 nvk0_1 = ~(vk0 >> 1) & v0mask;
            const GMUInt64 nvk1_0 = ~ vk1       & v1mask;
            const GMUInt64 nvk1_1 = ~(vk1 >> 1) & v1mask;

            const int r000 = gm_popcount64((nvi0_0 & nvj0_0 & nvk0_0) |
                                         ( (nvi0_0 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 & nvj0_1 & nvk0_0) |
                                         ( (nvi0_0 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_0 & nvk0_0) |
                                         ( (nvi0_1 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_1 & nvk0_0) |
                                         ( (nvi0_1 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_0 & nvk1_0) |
                                         ( (nvi1_0 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_1 & nvk1_0) |
                                         ( (nvi1_0 & nvj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_0 & nvk1_0) |
                                         ( (nvi1_1 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_1 & nvk1_0) |
                                         ( (nvi1_1 & nvj1_1 & nvk1_1) << 1 ));
            const int r001 = gm_popcount64((nvi0_0 & nvj0_0 &  vk0_0) |
                                         ( (nvi0_0 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 & nvj0_1 &  vk0_0) |
                                         ( (nvi0_0 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_0 &  vk0_0) |
                                         ( (nvi0_1 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_1 &  vk0_0) |
                                         ( (nvi0_1 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_0 &  vk1_0) |
                                         ( (nvi1_0 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_1 &  vk1_0) |
                                         ( (nvi1_0 & nvj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_0 &  vk1_0) |
                                         ( (nvi1_1 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_1 &  vk1_0) |
                                         ( (nvi1_1 & nvj1_1 &  vk1_1) << 1 ));
            const int r010 = gm_popcount64((nvi0_0 &  vj0_0 & nvk0_0) |
                                         ( (nvi0_0 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 &  vj0_1 & nvk0_0) |
                                         ( (nvi0_0 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_0 & nvk0_0) |
                                         ( (nvi0_1 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_1 & nvk0_0) |
                                         ( (nvi0_1 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_0 & nvk1_0) |
                                         ( (nvi1_0 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_1 & nvk1_0) |
                                         ( (nvi1_0 &  vj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_0 & nvk1_0) |
                                         ( (nvi1_1 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_1 & nvk1_0) |
                                         ( (nvi1_1 &  vj1_1 & nvk1_1) << 1 ));
            const int r011 = gm_popcount64((nvi0_0 &  vj0_0 &  vk0_0) |
                                         ( (nvi0_0 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 &  vj0_1 &  vk0_0) |
                                         ( (nvi0_0 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_0 &  vk0_0) |
                                         ( (nvi0_1 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_1 &  vk0_0) |
                                         ( (nvi0_1 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_0 &  vk1_0) |
                                         ( (nvi1_0 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_1 &  vk1_0) |
                                         ( (nvi1_0 &  vj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_0 &  vk1_0) |
                                         ( (nvi1_1 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_1 &  vk1_0) |
                                         ( (nvi1_1 &  vj1_1 &  vk1_1) << 1 ));
            const int r100 = gm_popcount64(( vi0_0 & nvj0_0 & nvk0_0) |
                                         ( ( vi0_0 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 & nvj0_1 & nvk0_0) |
                                         ( ( vi0_0 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_0 & nvk0_0) |
                                         ( ( vi0_1 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_1 & nvk0_0) |
                                         ( ( vi0_1 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_0 & nvk1_0) |
                                         ( ( vi1_0 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_1 & nvk1_0) |
                                         ( ( vi1_0 & nvj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_0 & nvk1_0) |
                                         ( ( vi1_1 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_1 & nvk1_0) |
                                         ( ( vi1_1 & nvj1_1 & nvk1_1) << 1 ));
            const int r101 = gm_popcount64(( vi0_0 & nvj0_0 &  vk0_0) |
                                         ( ( vi0_0 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 & nvj0_1 &  vk0_0) |
                                         ( ( vi0_0 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_0 &  vk0_0) |
                                         ( ( vi0_1 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_1 &  vk0_0) |
                                         ( ( vi0_1 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_0 &  vk1_0) |
                                         ( ( vi1_0 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_1 &  vk1_0) |
                                         ( ( vi1_0 & nvj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_0 &  vk1_0) |
                                         ( ( vi1_1 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_1 &  vk1_0) |
                                         ( ( vi1_1 & nvj1_1 &  vk1_1) << 1 ));
            const int r110 = gm_popcount64(( vi0_0 &  vj0_0 & nvk0_0) |
                                         ( ( vi0_0 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 &  vj0_1 & nvk0_0) |
                                         ( ( vi0_0 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_0 & nvk0_0) |
                                         ( ( vi0_1 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_1 & nvk0_0) |
                                         ( ( vi0_1 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_0 & nvk1_0) |
                                         ( ( vi1_0 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_1 & nvk1_0) |
                                         ( ( vi1_0 &  vj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_0 & nvk1_0) |
                                         ( ( vi1_1 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_1 & nvk1_0) |
                                         ( ( vi1_1 &  vj1_1 & nvk1_1) << 1 ));
            const int r111 = gm_popcount64(( vi0_0 &  vj0_0 &  vk0_0) |
                                         ( ( vi0_0 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 &  vj0_1 &  vk0_0) |
                                         ( ( vi0_0 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_0 &  vk0_0) |
                                         ( ( vi0_1 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_1 &  vk0_0) |
                                         ( ( vi0_1 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_0 &  vk1_0) |
                                         ( ( vi1_0 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_1 &  vk1_0) |
                                         ( ( vi1_0 &  vj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_0 &  vk1_0) |
                                         ( ( vi1_1 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_1 &  vk1_0) |
                                         ( ( vi1_1 &  vj1_1 &  vk1_1) << 1 ));

            /*---Accumulate---*/

            sum.data[0] += GMTally1_encode(r000, r001);
            sum.data[1] += GMTally1_encode(r010, r011);
            sum.data[2] += GMTally1_encode(r100, r101);
            sum.data[3] += GMTally1_encode(r110, r111);
          } /*---for pvfl---*/
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_i, i, env);
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_j, j, env);
          const GMTally1 sk1 = (GMTally1)GMVectorSums_sum(vs_k, k, env);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);
          if (GMEnv_all2all(env)) {
            GMMetrics_tally4x2_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, sum, &index_cache, env);
            GMMetrics_float3_S_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, si1_sj1_sk1, &index_cache, env);
            if (env->sparse) {
              const GMTally1 ci1 = (GMTally1)GMVectorSums_count(vs_i, i, env);
              const GMTally1 cj1 = (GMTally1)GMVectorSums_count(vs_j, j, env);
              const GMTally1 ck1 = (GMTally1)GMVectorSums_count(vs_k, k, env);
              const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
              GMMetrics_float3_C_set_all2all_3_permuted_cache(metrics,
                I, J, K, j_block, k_block, ci1_cj1_ck1, &index_cache, env);
            } /*---if sparse---*/
          } else /*---! all2all---*/ {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
            GMMetrics_float3_S_set_3(metrics, i, j, k, si1_sj1_sk1, env);
            if (env->sparse) {
              const GMTally1 ci1 = (GMTally1)GMVectorSums_count(vs_i, i, env);
              const GMTally1 cj1 = (GMTally1)GMVectorSums_count(vs_j, j, env);
              const GMTally1 ck1 = (GMTally1)GMVectorSums_count(vs_k, k, env);
              const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
              GMMetrics_float3_C_set_3(metrics, i, j, k, ci1_cj1_ck1, env);
            } /*---if sparse---*/
          } /*---if all2all---*/
        } //---I
        metrics->num_elts_local_computed += I_max - I_min;
      } //---K
    } //---J
    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/
    GMAssertAlways(false ? "logic error; should never be executed" : 0);
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMSectionInfo_destroy(si, env);;
}

//=============================================================================

void GMComputeNumerators3Way_create(GMComputeNumerators3Way* this_,
                                    int nvl, int npvfl, GMEnv* env) {
  GMAssertAlways(this_ && env);
  GMAssertAlways(nvl >= 0 && npvfl >= 0);

  this_->matM_ij_buf = GMMirroredBuf_null();
  this_->matM_jk_buf = GMMirroredBuf_null();
  this_->matM_kik_buf = GMMirroredBuf_null();

  for (int i=0; i<2; ++i) {
    this_->tmp_buf[i] = GMMirroredBuf_null();
    this_->matX_buf[i] = GMMirroredBuf_null();
    this_->matB_buf[i] = GMMirroredBuf_null();
  }

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce) {
        GMMirroredBuf_create(&(this_->tmp_buf[i]), nvl, nvl, env);
      }
      GMMirroredBuf_create(&(this_->matX_buf[i]), npvfl, nvl, env);
      GMMirroredBuf_create(&(this_->matB_buf[i]), nvl, nvl, env);
    }
    if (env->need_2way) {
      GMMirroredBuf_create(&(this_->matM_ij_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_jk_buf), nvl, nvl, env);
      GMMirroredBuf_create(&(this_->matM_kik_buf), nvl, nvl, env);
    }
  }
}

//=============================================================================

void GMComputeNumerators3Way_destroy(GMComputeNumerators3Way* this_,
                                     GMEnv* env) {
  GMAssertAlways(this_ && env);

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    for (int i=0; i<2; ++i) {
      if (env->do_reduce) {
        GMMirroredBuf_destroy(&this_->tmp_buf[i], env);
      }
      GMMirroredBuf_destroy(&this_->matX_buf[i], env);
      GMMirroredBuf_destroy(&this_->matB_buf[i], env);
    }
    if (env->need_2way) {
      GMMirroredBuf_destroy(&this_->matM_ij_buf, env);
      GMMirroredBuf_destroy(&this_->matM_jk_buf, env);
      GMMirroredBuf_destroy(&this_->matM_kik_buf, env);
    }
  }
}

//=============================================================================
/*---Start calculation of numerators, 3-way generic---*/

void GMComputeNumerators3Way_start(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredBuf* vectors_i_buf,
    GMMirroredBuf* vectors_j_buf,
    GMMirroredBuf* vectors_k_buf,
    int j_block,
    int k_block,
    const GMVectorSums* vector_sums_i,
    const GMVectorSums* vector_sums_j,
    const GMVectorSums* vector_sums_k,
    int section_step,
    GMEnv* env) {
  GMAssertAlways(this_ && metrics && env);
  GMAssertAlways(vectors_i && vectors_j && vectors_k);
  GMAssertAlways(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMAssertAlways(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == j_block &&
                   GMEnv_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(GMEnv_proc_num_vector_i(env) == k_block &&
                   GMEnv_proc_num_vector_i(env) != j_block));
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssertAlways(vector_sums_i && vector_sums_j && vector_sums_k);

  /*----------------------------------------*/
  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/
    gm_compute_numerators_3way_gpu_start_(this_,
                                          vectors_i, vectors_j, vectors_k,
                                          metrics, vectors_i_buf, vectors_j_buf,
                                          vectors_k_buf, j_block, k_block,
                                          vector_sums_i, vector_sums_j,
                                          vector_sums_k,
                                          section_step, env);
    /*----------------------------------------*/
  } else /*---(GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU)---*/ {
    /*----------------------------------------*/
    switch (GMEnv_metric_type(env)) {
      case GM_METRIC_TYPE_CZEK: {
        gm_compute_czek_numerators_3way_nongpu_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      case GM_METRIC_TYPE_CCC: {
        gm_compute_ccc_numerators_3way_nongpu_start_(this_,
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      default:
        GMInsist(env, false ? "Unimplemented." : 0);
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

//=============================================================================

//-----------------------------------------------------------------------------
