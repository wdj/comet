/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils_3way.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, 3-way.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_utils_linalg.hh"
#include "compute_metrics_utils.hh"
#include "compute_metrics_utils_3way.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Start calculation of numerators, 3-way Czekanowski non-gpu---*/

void gm_compute_czekanowski_numerators_3way_nongpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
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

    GMInsist(env, GMEnv_num_proc_field(env) == 1
                      ? "num_proc_field>1 for REF/CPU cases not supported"
                      : 0);

    GMAssertAlways(GMEnv_num_section_steps(env, 1) == 1);

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

    GMInsist(env, GMEnv_num_proc_field(env) == 1
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

    GMAssertAlways(GM_BOOL_FALSE
                 ? "logic error - code branch should never be executed"
                 : 0);

  } /*---if GPU---*/

  GMSectionInfo_destroy(si, env);;
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way CCC non-gpu---*/

void gm_compute_ccc_numerators_3way_nongpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
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

    GMInsist(env, GMEnv_num_proc_field(env) == 1
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

            const _Bool unknown_i = env->sparse ? vi == GM_2BIT_UNKNOWN
                                                : GM_BOOL_FALSE;
            const _Bool unknown_j = env->sparse ? vj == GM_2BIT_UNKNOWN
                                                : GM_BOOL_FALSE;
            const _Bool unknown_k = env->sparse ? vk == GM_2BIT_UNKNOWN
                                                : GM_BOOL_FALSE;

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

    GMInsist(env, GMEnv_num_proc_field(env) == 1
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
    GMAssertAlways(GM_BOOL_FALSE ? "logic error; should never be executed" : 0);
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMSectionInfo_destroy(si, env);;
}

/*===========================================================================*/

void gm_compute_numerators_3way_gpu_form_matX_(
  const GMVectors* const vectors_i,
  const GMMirroredPointer* const vectors_I_buf,
  const GMMirroredPointer* const vectors_J_buf,
  GMMirroredPointer* const matX_buf,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  GMEnv* const env) {

  /*--------------------*/
  /*---Populate leading columns of matX---*/
  /*--------------------*/

  const int nfl = vectors_i->num_field_local;
  const int npvfl = vectors_i->num_packedval_field_local;

  /*----------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
    /*----------*/
#pragma omp parallel for collapse(2)
    for (int I = I_min; I < I_max; ++I) {
      /*---Operate on columns x_i and x_j elementwise---*/
      for (int f = 0; f < npvfl; ++f) {
        const GMFloat a = ((GMFloat*)(vectors_I_buf->h))[f + npvfl*I];
        const GMFloat b = ((GMFloat*)(vectors_J_buf->h))[f + npvfl*J];
        ((GMFloat*)(matX_buf->h))[f + npvfl * I] = a < b ? a : b;
      }  //---for f---//
    }    //---for I---//
    /*----------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC) {
    /*----------*/
    for (int I = I_min; I < I_max; ++I) {

      const int num_field_active_local =
        GMEnv_proc_num_field(env) == GMEnv_num_proc_field(env) - 1
        ? nfl - (vectors_i->num_field - vectors_i->num_field_active) : nfl;
      const int num_packedval_field_active_local =
        (num_field_active_local + 64 - 1) / 64;

      const int num_seminibbles_edge = 1 + (num_field_active_local-1) % 64;

      const GMUInt64 nobits = 0;
      const GMUInt64 allbits = ~nobits;
      const GMUInt64 oddbits = 0x5555555555555555;

      const GMUInt64 edgemask0 = num_seminibbles_edge >= 32 ?
                                 allbits :
                                 allbits >> (64 - 2*num_seminibbles_edge);
      const GMUInt64 edgemask1 = num_seminibbles_edge <= 32 ?
                                 nobits :
                                 num_seminibbles_edge == 64 ?
                                 allbits :
                                 allbits >> (128 - 2*num_seminibbles_edge);
      /*---Operate on columns v_i and v_j elementwise---*/
      const int pvfl_edge = num_packedval_field_active_local - 1;
      for (int pvfl = 0; pvfl < npvfl; ++pvfl) {
        const int indI = pvfl + npvfl * I;
        const int indJ = pvfl + npvfl * J;

        const GMUInt64 activebits0 = pvfl < pvfl_edge ? allbits :
                                     pvfl == pvfl_edge ? edgemask0 : nobits;
        const GMUInt64 activebits1 = pvfl < pvfl_edge ? allbits :
                                     pvfl == pvfl_edge ? edgemask1 : nobits;

        const GMUInt64 activebits[2] = {activebits0, activebits1};

        const _Bool sparse = env->sparse;

        for (int word = 0; word<2; ++word) {

          const GMUInt64 vI =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_I_buf->h))[indI].data[word]);
          const GMUInt64 vJ =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_J_buf->h))[indJ].data[word]);

          const GMUInt64 vIx = step_2way==0 ? vI | ~activebits[word] : vI;

          /*---Create word whose odd bits sample the lo or hi bit of interest
               of the seminibble.  Also create the complement thereof---*/

          const GMUInt64  vI_0 =   vIx        & oddbits;
          const GMUInt64  vI_1 =  (vIx >> 1)  & oddbits;
          const GMUInt64 nvI_0 = ~ vIx        & oddbits;
          const GMUInt64 nvI_1 = ~(vIx >> 1)  & oddbits;

          const GMUInt64  vJ_0 =   vJ            & oddbits;
          const GMUInt64  vJ_1 =  (vJ  >> 1)     & oddbits;

          /*---Create a mask whose odd bits denote whether each respective
               seminibble of vector I matches the case we are handling
               in this 2-way step (and the complement thereof)---*/

          const GMUInt64  vI_match =
            step_2way==0 ?  nvI_0 & nvI_1  & oddbits : // 00
            step_2way==1 && sparse ?
                           ( vI_0 & nvI_1) & oddbits : // 01
            step_2way==1 ? ( vI_0 ^  vI_1) & oddbits : // 01, 10
          /*step_2way==2*/   vI_0 &  vI_1  & oddbits;  // 11

          const GMUInt64 nvI_match =
            step_2way==0 ? ( vI_0 |  vI_1) & oddbits :
            step_2way==1 && sparse ?
                           (nvI_0 |  vI_1) & oddbits :
            step_2way==1 ? ( vI_0 ^ nvI_1) & oddbits :
          /*step_2way==2*/ (nvI_0 | nvI_1) & oddbits;

          /*---Construct the lo and hi bit of the result seminibble, based
               on truth table:
            - lo bit is 1 (00 or 01) if vJ is 01, 10 or 11 and vI is the
              case being handled for this 2-way step.
            - hi bit is 1 (10 or 11) if vJ is 11 or if vI is a case not
              being handled for this 2-way step.
          */

          const GMUInt64 r_0 =  vI_match & (sparse ? vJ_0 : vJ_0 | vJ_1);
          const GMUInt64 r_1 = nvI_match | (sparse ? vJ_1 : vJ_0 & vJ_1);

          const GMUInt64 r = r_0 | (r_1 << 1);
          ((GMBits2x64*)(matX_buf->h))[indI].data[word] = *(GMBits1_2x64*)&r;

        } /*---word---*/
      }  //---for f---//
    }    //---for I---//
    /*----------*/
  } /*---GMEnv_metric_type(env)---*/
  /*----------*/
}

/*===========================================================================*/

void gm_compute_numerators_3way_gpu_form_metrics_(
  const GMMirroredPointer* const matM_IJ_buf,
  const GMMirroredPointer* const matM_JK_buf,
  const GMMirroredPointer* const matM_KIK_buf,
  const GMMirroredPointer* const matB_buf,
  GMMetrics* metrics,
  const int nvl,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  const int K_min,
  const int K_max,
  const int j_block,
  const int k_block,
  const GMSectionInfo* const si,
  const GMVectorSums* vector_sums_i,
  const GMVectorSums* vector_sums_j,
  const GMVectorSums* vector_sums_k,
  GMEnv* const env) {

  GMAssertAlways(vector_sums_i && vector_sums_j && vector_sums_k);

  const _Bool is_part3 = si->is_part3;

  const GMVectorSums* const vs_i = vector_sums_i;
  const GMVectorSums* const vs_j = vector_sums_j;
  const GMVectorSums* const vs_k = vector_sums_k;

  /* clang-format off */
  const GMVectorSums* const vs_I =
          !si->is_part3 ?   vs_i :
               si->sax0 ?   vs_k :
               si->sax1 ?   vs_i :
            /* si->sax2 ?*/ vs_j;
  const GMVectorSums* const vs_J =
          !si->is_part3 ?   vs_j :
               si->sax0 ?   vs_i :
               si->sax1 ?   vs_j :
            /* si->sax2 ?*/ vs_k;
  const GMVectorSums* const vs_K =
          !si->is_part3 ?   vs_k :
               si->sax0 ?   vs_j :
               si->sax1 ?   vs_k :
            /* si->sax2 ?*/ vs_i;
  /* clang-format on */

  const size_t nvl64 = (size_t)nvl;
  const size_t I_max64 = (size_t)I_max;

  /*--------------------*/
  /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
  /*--------------------*/

  /*----------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI &&
      !GMEnv_all2all(env)) {
    /*----------*/

#pragma omp parallel for collapse(2)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
      const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + nvl * J];
        const GMFloat min_JK = ((GMFloat*)(matM_JK_buf->h))[J + nvl64*K];
        const GMFloat min_KIK = ((GMFloat*)(matM_KIK_buf->h))[K + nvl64*I];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max64*K];
        const GMFloat numer = min_IJ + min_JK + min_KIK - min_IJK;
        const int i = I;
        const int j = J;
        const int k = K;
        /*---Make arithmetic order-independent---*/
        GMFloat smin, smid, smax;
        const GMFloat si = GMVectorSums_sum(vs_i, i, env);
        const GMFloat sj = GMVectorSums_sum(vs_i, j, env);
        const GMFloat sk = GMVectorSums_sum(vs_i, k, env);
        GMFloat_sort_3(&smin, &smid, &smax, &si, &sj, &sk);
        const GMFloat denom = smin + smid + smax;
        const GMFloat value = ((GMFloat)1.5) * numer / denom;
        GMMetrics_float_set_3(metrics, i, j, k, value, env);
      } /*---for K---*/
    }   /*---for I---*/
    metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                        (K_max - K_min);
    /*----------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI &&
      GMEnv_all2all(env)) {
    /*----------*/

    GMIndexCache index_cache = {0};
#pragma omp parallel for collapse(2) firstprivate(index_cache)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + nvl64 * J];
        const GMFloat min_JK = ((GMFloat*)(matM_JK_buf->h))[J + nvl64 * K];
        const GMFloat min_KIK =
            is_part3 ? ((GMFloat*)(matM_KIK_buf->h))[K + nvl64 * I]
                     : ((GMFloat*)(matM_KIK_buf->h))[I + nvl64 * K];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max64 * K];
        const GMFloat numer = min_IJ + min_JK + min_KIK - min_IJK;
        /*---Make arithmetic order-independent---*/
        GMFloat smin, smid, smax;
        const GMFloat sI = GMVectorSums_sum(vs_I, I, env);
        const GMFloat sJ = GMVectorSums_sum(vs_J, J, env);
        const GMFloat sK = GMVectorSums_sum(vs_K, K, env);
        GMFloat_sort_3(&smin, &smid, &smax, &sI, &sJ, &sK);
        const GMFloat denom = smin + smid + smax;
        const GMFloat value = ((GMFloat)1.5) * numer / denom;
        GMMetrics_float_set_all2all_3_permuted_cache(metrics, I, J, K,
            j_block, k_block, value, &index_cache, env);
      } /*---for K---*/
    }   /*---for I---*/
    metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                        (K_max - K_min);
    /*----------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             !GMEnv_all2all(env)) {
    /*----------*/
#pragma omp parallel for collapse(2)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        /*---This is the notall2all case -- has no axis permutation---*/

        const int i = I;
        const int j = J;
        const int k = K;

        GMTally4x2 numer = step_2way==0 ? GMTally4x2_null() :
                        GMMetrics_tally4x2_get_3(metrics, i, j, k, env);

        GMTally1 r000, r001;
        GMTally1_decode(&r000, &r001, numer.data[0]);
        GMTally1 r010, r011;
        GMTally1_decode(&r010, &r011, numer.data[1]);
        GMTally1 r100, r101;
        GMTally1_decode(&r100, &r101, numer.data[2]);
        GMTally1 r110, r111;
        GMTally1_decode(&r110, &r111, numer.data[3]);

        const GMTally2x2 mB = ((GMTally2x2*)(matB_buf->h))[I + I_max * K];
        GMTally1 mB00, mB01;
        GMTally1_decode(&mB00, &mB01, mB.data[0]);
        GMTally1 mB10, mB11;
        GMTally1_decode(&mB10, &mB11, mB.data[1]);

        if (step_2way==0) {
          r000 += 2 * mB00;
          r001 += 2 * mB01;
          r010 += 2 * mB10;
          r011 += 2 * mB11;
        } else if (step_2way==1) {
          r000 += mB00;
          r001 += mB01;
          r010 += mB10;
          r011 += mB11;
          r100 += mB00;
          r101 += mB01;
          r110 += mB10;
          r111 += mB11;
        } else /*---step_2way==2---*/ {
          r100 += 2 * mB00;
          r101 += 2 * mB01;
          r110 += 2 * mB10;
          r111 += 2 * mB11;
        }
        /*---NOTE: pay attention to order here---*/
        numer.data[0] = GMTally1_encode(r000, r001);
        numer.data[1] = GMTally1_encode(r010, r011);
        numer.data[2] = GMTally1_encode(r100, r101);
        numer.data[3] = GMTally1_encode(r110, r111);
//        if (step_2way==0) {
          GMMetrics_tally4x2_set_3(metrics, i, j, k, numer, env);
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_i, i, env);
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_i, j, env);
          const GMTally1 sk1 = (GMTally1)GMVectorSums_sum(vs_i, k, env);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);
          GMMetrics_float3_S_set_3(metrics, i, j, k, si1_sj1_sk1, env);
          if (env->sparse) {
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_i, i, env);
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_i, j, env);
            const GMTally1 ck = (GMTally1)GMVectorSums_count(vs_i, k, env);
            const GMFloat3 ci_cj_ck = GMFloat3_encode(ci, cj, ck);
            GMMetrics_float3_C_set_3(metrics, i, j, k, ci_cj_ck, env);
          } /*---if sparse---*/
//        }
      } /*---for K---*/
    }   /*---for I---*/
    if (step_2way == 2) {
      metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                          (K_max - K_min);
    }
    /*----------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
             GMEnv_all2all(env)) {
    /*----------*/
    GMIndexCache index_cache = {0};
#pragma omp parallel for collapse(2) firstprivate(index_cache)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
/*---For the permuted case,
 1) pay attention to KIK access
 2) swap 01 and 10 if needed.
---*/

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

        GMTally4x2 numer = step_2way==0 ? GMTally4x2_null() :
          GMMetrics_tally4x2_get_all2all_3_permuted_cache(metrics, I, J, K,
                                           j_block, k_block, &index_cache, env);

        GMTally1 r000_permuted, r001_permuted;
        GMTally1_decode(&r000_permuted, &r001_permuted, numer.data[0]);
        GMTally1 r010_permuted, r011_permuted;
        GMTally1_decode(&r010_permuted, &r011_permuted, numer.data[1]);
        GMTally1 r100_permuted, r101_permuted;
        GMTally1_decode(&r100_permuted, &r101_permuted, numer.data[2]);
        GMTally1 r110_permuted, r111_permuted;
        GMTally1_decode(&r110_permuted, &r111_permuted, numer.data[3]);

        const GMTally2x2 mB = ((GMTally2x2*)(matB_buf->h))[I + I_max * K];
        GMTally1 mB00, mB01;
        GMTally1_decode(&mB00, &mB01, mB.data[0]);
        GMTally1 mB10, mB11;
        GMTally1_decode(&mB10, &mB11, mB.data[1]);

        /* clang-format off */
        int r000 = r000_permuted;

        int r100 = !si->is_part3 ?   r100_permuted :
                        si->sax0 ?   r001_permuted :
                        si->sax1 ?   r100_permuted :
                     /* si->sax2 ?*/ r010_permuted;
        int r010 = !si->is_part3 ?   r010_permuted :
                        si->sax0 ?   r100_permuted :
                        si->sax1 ?   r010_permuted :
                     /* si->sax2 ?*/ r001_permuted;
        int r001 = !si->is_part3 ?   r001_permuted :
                        si->sax0 ?   r010_permuted :
                        si->sax1 ?   r001_permuted :
                     /* si->sax2 ?*/ r100_permuted;

        int r011 = !si->is_part3 ?   r011_permuted :
                        si->sax0 ?   r110_permuted :
                        si->sax1 ?   r011_permuted :
                     /* si->sax2 ?*/ r101_permuted;
        int r101 = !si->is_part3 ?   r101_permuted :
                        si->sax0 ?   r011_permuted :
                        si->sax1 ?   r101_permuted :
                     /* si->sax2 ?*/ r110_permuted;
        int r110 = !si->is_part3 ?   r110_permuted :
                        si->sax0 ?   r101_permuted :
                        si->sax1 ?   r110_permuted :
                     /* si->sax2 ?*/ r011_permuted;

        int r111 = r111_permuted;
        /* clang-format on */

        if (step_2way==0) {
          r000 += 2 * mB00;
          r001 += 2 * mB01;
          r010 += 2 * mB10;
          r011 += 2 * mB11;
        } else if (step_2way==1) {
          r000 += mB00;
          r001 += mB01;
          r010 += mB10;
          r011 += mB11;
          r100 += mB00;
          r101 += mB01;
          r110 += mB10;
          r111 += mB11;
        } else /*---step_2way==2---*/ {
          r100 += 2 * mB00;
          r101 += 2 * mB01;
          r110 += 2 * mB10;
          r111 += 2 * mB11;
        }

        /* clang-format off */
        r000_permuted = r000;

        r100_permuted = !si->is_part3 ?   r100 :
                             si->sax0 ?   r010 :
                             si->sax1 ?   r100 :
                          /* si->sax2 ?*/ r001;
        r010_permuted = !si->is_part3 ?   r010 :
                             si->sax0 ?   r001 :
                             si->sax1 ?   r010 :
                          /* si->sax2 ?*/ r100;
        r001_permuted = !si->is_part3 ?   r001 :
                             si->sax0 ?   r100 :
                             si->sax1 ?   r001 :
                          /* si->sax2 ?*/ r010;

        r011_permuted = !si->is_part3 ?   r011 :
                             si->sax0 ?   r101 :
                             si->sax1 ?   r011 :
                          /* si->sax2 ?*/ r110;
        r101_permuted = !si->is_part3 ?   r101 :
                             si->sax0 ?   r110 :
                             si->sax1 ?   r101 :
                          /* si->sax2 ?*/ r011;
        r110_permuted = !si->is_part3 ?   r110 :
                             si->sax0 ?   r011 :
                             si->sax1 ?   r110 :
                          /* si->sax2 ?*/ r101;

        r111_permuted = r111;
        /* clang-format on */

        /*---NOTE: pay attention to order here---*/

        numer.data[0] = GMTally1_encode(r000_permuted, r001_permuted);
        numer.data[1] = GMTally1_encode(r010_permuted, r011_permuted);
        numer.data[2] = GMTally1_encode(r100_permuted, r101_permuted);
        numer.data[3] = GMTally1_encode(r110_permuted, r111_permuted);
        GMMetrics_tally4x2_set_all2all_3_permuted_cache(metrics, I, J, K,
                                   j_block, k_block, numer, &index_cache, env);
//        if (step_2way==0) {
          const GMTally1 si1 = (GMTally1)GMVectorSums_sum(vs_i, i, env);
          const GMTally1 sj1 = (GMTally1)GMVectorSums_sum(vs_j, j, env); 
          const GMTally1 sk1 = (GMTally1)GMVectorSums_sum(vs_k, k, env); 
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);
          GMMetrics_float3_S_set_all2all_3_permuted_cache(metrics, I, J, K,
              j_block, k_block, si1_sj1_sk1, &index_cache, env);
          if (env->sparse) {
            const GMTally1 ci = (GMTally1)GMVectorSums_count(vs_i, i, env);
            const GMTally1 cj = (GMTally1)GMVectorSums_count(vs_j, j, env); 
            const GMTally1 ck = (GMTally1)GMVectorSums_count(vs_k, k, env); 
            const GMFloat3 ci_cj_ck = GMFloat3_encode(ci, cj, ck);
            GMMetrics_float3_C_set_all2all_3_permuted_cache(metrics, I, J, K,
                j_block, k_block, ci_cj_ck, &index_cache, env);
          } /*---if sparse---*/
//        }
      } /*---for K---*/
    }   /*---for I---*/
    if (step_2way == 2) {
      metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                          (K_max - K_min);
    }
    /*----------*/
  } else {
    /*----------*/
    GMAssertAlways(GM_BOOL_FALSE);
    /*----------*/
  } /*---GMEnv_metric_type(env)---*/
  /*----------*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way gpu---*/

void gm_compute_numerators_3way_gpu_start_(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
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
  GMAssertAlways(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMAssertAlways(vector_sums_i && vector_sums_j && vector_sums_k);

  /*---Initializations---*/

  int mpi_code = 0;
  mpi_code *= 1; /*---Avoid unused variable warning---*/

  const int nvl = metrics->num_vector_local;
  const int npvfl = vectors_i->num_packedval_field_local;

  const int i_block = GMEnv_proc_num_vector_i(env);

  const _Bool need_2way = GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI;
  const _Bool need_reduce = GMEnv_num_proc_field(env) > 1;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const _Bool need_mat_ij = need_2way;
  const _Bool need_mat_jk = need_2way && !si->is_part1;
  const _Bool need_mat_kik = need_2way && si->is_part3;

  /*----------------------------------------*/
  /*---First get the required 2-way ij, jk, ik metrics---*/
  /*----------------------------------------*/

  /*--------------------*/
  /*---Compute i_block - j_block PROD---*/
  /*--------------------*/

  GMMirroredPointer* mat_buf_tmp[2] = {&this_->mat_buf_tmp[0],
                                       &this_->mat_buf_tmp[1]};

  GMMirroredPointer* const matM_ij_buf = need_mat_ij ? &this_->matM_ij_buf :
                                                      NULL;

  if (need_mat_ij) {
    GMMirroredPointer* matM_ij_buf_local =
       need_reduce ? mat_buf_tmp[0] : matM_ij_buf;

    gm_linalg_set_matrix_zero_start(matM_ij_buf_local, nvl, nvl, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_i_buf->d, npvfl,
                         vectors_j_buf->d, npvfl,
                         matM_ij_buf_local->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_ij_buf_local, env);
    gm_get_metrics_wait(metrics, matM_ij_buf_local, env);

    if (need_reduce) {
      gm_reduce_metrics(metrics, matM_ij_buf, matM_ij_buf_local, env);
    }
  }

  /*--------------------*/
  /*---Compute j_block - k_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  GMMirroredPointer* const matM_jk_buf =
      !si->is_part1 ? &this_->matM_jk_buf : matM_ij_buf;

  if (need_mat_jk) {
    GMMirroredPointer* matM_jk_buf_local =
        need_reduce ? mat_buf_tmp[0] : matM_jk_buf;

    gm_linalg_set_matrix_zero_start(matM_jk_buf_local, nvl, nvl, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_j_buf->d, npvfl,
                         vectors_k_buf->d, npvfl,
                         matM_jk_buf_local->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_jk_buf_local, env);
    gm_get_metrics_wait(metrics, matM_jk_buf_local, env);

    if (need_reduce) {
      gm_reduce_metrics(metrics, matM_jk_buf, matM_jk_buf_local, env);
    }
  }

  /*--------------------*/
  /*---Compute k_block - i_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  /*---NOTE: for Part 3, this is indexed directly as (k,i).
       Otherwise, it is indexed through an alias as (i,k)---*/

  GMMirroredPointer* const matM_kik_buf = si->is_part3
    ? &this_->matM_kik_buf : matM_ij_buf;

  if (need_mat_kik) {
    GMMirroredPointer* matM_kik_buf_local =
        need_reduce ? mat_buf_tmp[0] : matM_kik_buf;

    gm_linalg_set_matrix_zero_start(matM_kik_buf_local, nvl, nvl, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_k_buf->d, npvfl,
                         vectors_i_buf->d, npvfl,
                         matM_kik_buf_local->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_kik_buf_local, env);
    gm_get_metrics_wait(metrics, matM_kik_buf_local, env);

    if (need_reduce) {
      gm_reduce_metrics(metrics, matM_kik_buf, matM_kik_buf_local, env);
    }
  } /*---is_part3---*/

  /*----------------------------------------*/
  /*---Now compute ijk piece, via an outer loop over j values---*/
  /*----------------------------------------*/

  /*---Allocate magma CPU/GPU memory for matrices X and B---*/
  // X = elementwise OP of one vector with each of the rest of the vectors.
  // For the jth iteration, the ith column of X is the elementwise OP
  //   of vectors i and j.
  // B = X^T PROD V = three way PROD.

  GMMirroredPointer* matX_buf[2] = {&this_->matX_buf[0], &this_->matX_buf[1]};
  GMMirroredPointer* matB_buf[2] = {&this_->matB_buf[0], &this_->matB_buf[1]};

  /*---Set up pointers to permute the access of axes for Part 3---*/
  /*---We use capitals I, J, K here to denote the PERMUTED axes---*/

  const _Bool is_ijk = !si->is_part3 ? GM_BOOL_TRUE : si->sax1;
  const _Bool is_kij = !si->is_part3 ? GM_BOOL_FALSE : si->sax0;
  const _Bool is_jki = !si->is_part3 ? GM_BOOL_FALSE : si->sax2;

  /* clang-format off */
  GMMirroredPointer* const vectors_I_buf = is_ijk ? vectors_i_buf :
                                           is_kij ? vectors_k_buf :
                                           is_jki ? vectors_j_buf : 0;
 
  GMMirroredPointer* const vectors_J_buf = is_ijk ? vectors_j_buf :
                                           is_kij ? vectors_i_buf :
                                           is_jki ? vectors_k_buf : 0;
 
  GMMirroredPointer* const vectors_K_buf = is_ijk ? vectors_k_buf :
                                           is_kij ? vectors_j_buf :
                                           is_jki ? vectors_i_buf : 0;
  
  //TODO - use is_ijk etc. 
  GMMirroredPointer* const matM_IJ_buf  = !si->is_part3 ? matM_ij_buf  :
                                               si->sax0 ? matM_kik_buf :
                                               si->sax1 ? matM_ij_buf  :
                                               si->sax2 ? matM_jk_buf  : 0;
  
  GMMirroredPointer* const matM_JK_buf  = !si->is_part3 ? matM_jk_buf  :
                                               si->sax0 ? matM_ij_buf  :
                                               si->sax1 ? matM_jk_buf  :
                                               si->sax2 ? matM_kik_buf : 0;
  
  GMMirroredPointer* const matM_KIK_buf = !si->is_part3 ? matM_kik_buf :
                                               si->sax0 ? matM_jk_buf  :
                                               si->sax1 ? matM_kik_buf :
                                               si->sax2 ? matM_ij_buf  : 0;
  /* clang-format on */

  /*--------------------*/
  /*---Collapsed loops over J and over 2-way steps---*/
  /*--------------------*/

  const int J_min = si->J_lb;
  const int J_max = si->J_ub;
  const int J_count = J_max - J_min;

  const int num_step_2way = GMEnv_metric_type(env)==GM_METRIC_TYPE_CCC ? 3 : 1;
  const int num_step = J_count * num_step_2way;
  const int extra_step = 1;

  MPI_Request mpi_requests[2];
//  _Bool matX_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matX_d_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_d_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_ptr_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};

  typedef struct {
    int step_num;
    int step_2way;
    int J;
    int I_min;
    int I_max;
    int K_min;
    int K_max;
    _Bool empty;
    _Bool is_compute_step;
    _Bool do_compute;
    int index_01;
    GMMirroredPointer* matB_buf_ptr;
  } LoopVars;

  LoopVars vars = {0};
  LoopVars vars_prev = {0};
  LoopVars vars_prevprev = {0};
  LoopVars vars_next = {0};
  vars.do_compute = GM_BOOL_FALSE;
  vars_prev.do_compute = GM_BOOL_FALSE;
  vars_prevprev.do_compute = GM_BOOL_FALSE;
  vars_next.do_compute = GM_BOOL_FALSE;

  for (int step_num = 0-extra_step; step_num < num_step+extra_step*2;
       ++step_num) {

    vars_prevprev = vars_prev;
    vars_prev = vars;
    vars = vars_next;

    vars_next.step_num = step_num + 1;
    vars_next.step_2way = gm_mod_i(vars_next.step_num, num_step_2way);
    vars_next.J = J_min + gm_floor_i(vars_next.step_num, num_step_2way);
    vars_next.I_min = 0;
    vars_next.I_max = si->is_part1 ? vars_next.J : nvl;
    vars_next.K_min = si->is_part3 ? 0 : vars_next.J + 1;
    vars_next.K_max = nvl;
    vars_next.empty = vars_next.I_min >= vars_next.I_max ||
                      vars_next.K_min >= vars_next.K_max;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.do_compute = vars_next.is_compute_step && ! vars_next.empty;
    vars_next.index_01 = gm_mod_i(vars_next.step_num, 2);
    vars_next.matB_buf_ptr = need_reduce ?  mat_buf_tmp[vars_next.index_01] :
                                            matB_buf[vars_next.index_01];

//    vars_next.matX_h_updating = &matX_h_updating[vars_next.index_01];
//    vars_next.matX_d_updating = &matX_d_updating[vars_next.index_01];
//    vars_next.matB_h_updating = &matB_h_updating[vars_next.index_01];
//    vars_next.matB_d_updating = &matB_d_updating[vars_next.index_01];
//    vars_next.matB_ptr_h_updating = &matB_ptr_h_updating[vars_next.index_01];

    //==========

    if (vars_next.do_compute) {
      /*---Populate leading columns of matX---*/
//      GMAssertAlways(!*vars_next.matX_h_updating);
//      *vars_next.matX_h_updating = GM_BOOL_TRUE;
      gm_compute_numerators_3way_gpu_form_matX_(vectors_i,
          vectors_I_buf, vectors_J_buf, matX_buf[vars_next.index_01],
          vars_next.J, vars_next.step_2way,
          vars_next.I_min, vars_next.I_max, env);
      //matX_buf[vars_next.index_01]->h is now overwritten.
//      *vars_next.matX_h_updating = GM_BOOL_FALSE;
    }

    //==========

    if (vars.do_compute) {
      /*---Send matrix matX to GPU - WAIT---*/
      gm_linalg_set_matrix_wait(env);
      //matX_buf[vars.index_01]->h is now no longer needed.
      //matX_buf[vars.index_01]->d is now available.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Perform pseudo GEMM matB = matX^T PROD V - WAIT---*/
      gm_compute_wait(env);
      //matX_buf[vars_prev.index_01]->d is now no longer needed.
      //matB_buf_ptr[vars_prev.index_01]->d is now available.
    }

    //==========

    if (vars_next.do_compute) {
      /*---Send matrix matX to GPU - START---*/
      gm_linalg_set_matrix_start(matX_buf[vars_next.index_01], npvfl,
                                 vars_next.I_max, env);
      //matX_buf[vars_next.index_01]->d (= matX_buf[vars_prev.index_01]->d)
      //is now being overwritten.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - START---*/
      gm_linalg_get_matrix_start(vars_prev.matB_buf_ptr, vars_prev.I_max,
                                 nvl, env);
      //gm_get_metrics_start(metrics, vars_prev.matB_buf_ptr, env); //FIX
      //matB_buf_ptr[vars_prev.index_01]->h is now being overwritten.
    }

    //==========

    if (vars.do_compute) {
      /*---Initialize result matrix to zero (apparently magma requires)---*/
      gm_linalg_set_matrix_zero_start(vars.matB_buf_ptr, nvl,
                                      vars.I_max, env);
      /*---Perform pseudo GEMM matB = matX^T PROD V - START---*/
      gm_linalg_gemm_start(vars.I_max, nvl, npvfl,
                           matX_buf[vars.index_01]->d, npvfl,
                           vectors_K_buf->d, npvfl,
                           vars.matB_buf_ptr->d, vars.I_max, env);
      //matB_buf_ptr[vars.index_01]->d is now being overwritten.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - WAIT---*/
      gm_linalg_get_matrix_wait(env);
      //gm_get_metrics_wait(metrics, vars_prev.matB_buf_ptr, env); //FIX
      //matB_buf_ptr[vars_prev.index_01]->h is now available.
      //matB_buf_ptr[vars_prev.index_01]->d is now no longer needed.
    }

    //==========

    if (vars_prevprev.do_compute) {
      if (need_reduce) {
        gm_reduce_metrics_wait(&(mpi_requests[vars_prevprev.index_01]), env); 
        //matB_buf[vars.prevprev.index_01]->h is now available.
      }
    }

    //==========

    if (vars_prev.do_compute) {
      if (need_reduce) {
        mpi_requests[vars_prev.index_01] = gm_reduce_metrics_start(metrics,
            matB_buf[vars_prev.index_01], vars_prev.matB_buf_ptr, env);
        //matB_buf[vars.prev.index_01]->h is now being overwritten.
      }
    }

    //==========

    //---NOTE: matB_buf[vars_prevprev.index_01]->d is being updated now
    //---but matB_buf[vars_prevprev.index_01]->h is ok.

    if (vars_prevprev.do_compute) {
      /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
      gm_compute_numerators_3way_gpu_form_metrics_(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
          matB_buf[vars_prevprev.index_01],
          metrics, nvl,
          vars_prevprev.J,
          vars_prevprev.step_2way,
          vars_prevprev.I_min,
          vars_prevprev.I_max,
          vars_prevprev.K_min,
          vars_prevprev.K_max,
          j_block, k_block, si,
          vector_sums_i, vector_sums_j, vector_sums_k,
          env);
      //matB_buf[vars.prevprev.index_01]->h (=matB_buf[vars.index_01]->h)
      //now no longer needed.
    }

    //==========

  } /*---for step_num---*/

  /*--------------------*/
  /*---Free memory---*/
  /*--------------------*/

  GMSectionInfo_destroy(si, env);
}

/*===========================================================================*/

void GMComputeNumerators3Way_create(
    GMComputeNumerators3Way* this_,
    int nvl,
    int npvfl,
    GMEnv* env) {
  GMAssertAlways(this_ && env);

  const _Bool need_2way = GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI;
  const _Bool need_reduce = GMEnv_num_proc_field(env) > 1;

  const size_t metrics_buf_size = nvl * (size_t)nvl;
  const size_t vectors_buf_size = nvl * (size_t)npvfl;

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    for (int i=0; i<2; ++i) {
      this_->mat_buf_tmp[i] = need_reduce
                        ? gm_linalg_malloc(metrics_buf_size, env)
                        : GMMirroredPointer_null();
      this_->matX_buf[i] = gm_linalg_malloc(vectors_buf_size, env);
      this_->matB_buf[i] = gm_linalg_malloc(metrics_buf_size, env);
    }

    this_->matM_ij_buf = need_2way ? gm_linalg_malloc(metrics_buf_size, env)
                                   : GMMirroredPointer_null();
    this_->matM_jk_buf = need_2way ? gm_linalg_malloc(metrics_buf_size, env)
                                   : GMMirroredPointer_null();
    this_->matM_kik_buf = need_2way ? gm_linalg_malloc(metrics_buf_size, env)
                                    : GMMirroredPointer_null();
  }
}

/*===========================================================================*/

void GMComputeNumerators3Way_destroy(GMComputeNumerators3Way* this_,
                                     GMEnv* env) {
  GMAssertAlways(this_ && env);

  const _Bool need_2way = GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI;
  const _Bool need_reduce = GMEnv_num_proc_field(env) > 1;

  if (GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {

    for (int i=0; i<2; ++i) {
      if (need_reduce) {
        gm_linalg_free(&this_->mat_buf_tmp[i], env);
      }
      gm_linalg_free(&this_->matX_buf[i], env);
      gm_linalg_free(&this_->matB_buf[i], env);
    }

    if (need_2way) {
      gm_linalg_free(&this_->matM_ij_buf, env);
      gm_linalg_free(&this_->matM_jk_buf, env);
      gm_linalg_free(&this_->matM_kik_buf, env);
    }
  }
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way generic---*/

void GMComputeNumerators3Way_start(
    GMComputeNumerators3Way* this_,
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
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
      case GM_METRIC_TYPE_CZEKANOWSKI: {
        gm_compute_czekanowski_numerators_3way_nongpu_start_(this_,
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
        GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
