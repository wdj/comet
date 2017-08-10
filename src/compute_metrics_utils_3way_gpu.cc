//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_utils_3way_gpu.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities, 3-way, gpu.
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
#include "comm_xfer_utils.hh"
#include "compute_metrics_utils_3way.hh"

//=============================================================================

void gm_compute_numerators_3way_gpu_form_matX_(
  const GMVectors* vectors_i,
  const GMMirroredBuf* vectors_I_buf,
  const GMMirroredBuf* vectors_J_buf,
  GMMirroredBuf* const matX_buf,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  GMEnv* const env) {

  /*--------------------*/
  /*---Populate leading columns of matX---*/
  /*--------------------*/

  const int npvfl = vectors_i->num_packedval_field_local;

  /*----------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK) {
    /*----------*/
#pragma omp parallel for collapse(2)
    for (int I = I_min; I < I_max; ++I) {
      /*---Operate on columns x_i and x_j elementwise---*/
      for (int f = 0; f < npvfl; ++f) {
        const GMFloat a = GMMirroredBuf_elt_const<GMFloat>(vectors_I_buf, f, I);
        const GMFloat b = GMMirroredBuf_elt_const<GMFloat>(vectors_J_buf, f, J);
        GMMirroredBuf_elt<GMFloat>(matX_buf, f, I) = a < b ? a : b;
      }  //---for f---//
    }    //---for I---//
    /*----------*/
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC) {
    /*----------*/
    for (int I = I_min; I < I_max; ++I) {

      const GMUInt64 oddbits = 0x5555555555555555;

      /*---Operate on columns v_i and v_j elementwise---*/
      for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

        const bool sparse = env->sparse;

        for (int word = 0; word<2; ++word) {
          const GMUInt64 vI = GMMirroredBuf_elt_const<GMBits2x64>(
                                           vectors_I_buf, pvfl, I).data[word];
          const GMUInt64 vJ = GMMirroredBuf_elt_const<GMBits2x64>(
                                           vectors_J_buf, pvfl, J).data[word];

          // Create word whose odd bits sample the lo or hi bit of interest
          // of the seminibble.  Also create the complement thereof.

          const GMUInt64  vI_0 =   vI        & oddbits;
          const GMUInt64  vI_1 =  (vI >> 1)  & oddbits;
          const GMUInt64 nvI_0 = ~ vI        & oddbits;
          const GMUInt64 nvI_1 = ~(vI >> 1)  & oddbits;

          const GMUInt64  vJ_0 =   vJ        & oddbits;
          const GMUInt64  vJ_1 =  (vJ  >> 1) & oddbits;

          // Create a mask whose odd bits denote whether each respective
          // seminibble of vector I matches the case we are handling
          // in this 2-way step (and the complement thereof).

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

          // Construct the lo and hi bit of the result seminibble, based
          // on truth table:
          //  - lo bit is 1 (00 or 01) if vJ is 01, 10 or 11 and vI is the
          //    case being handled for this 2-way step.
          //  - hi bit is 1 (10 or 11) if vJ is 11 or if vI is a case not
          //    being handled for this 2-way step.

          const GMUInt64 r_0 =  vI_match & (sparse ? vJ_0 : vJ_0 | vJ_1);
          const GMUInt64 r_1 = nvI_match | (sparse ? vJ_1 : vJ_0 & vJ_1);

          const GMUInt64 r = r_0 | (r_1 << 1);

          GMMirroredBuf_elt<GMBits2x64>(matX_buf, pvfl, I).data[word] = r;
        } /*---word---*/
      }  //---for f---//
    }    //---for I---//
    /*----------*/
  } /*---GMEnv_metric_type(env)---*/
  /*----------*/
}

//=============================================================================

void gm_compute_numerators_3way_gpu_form_metrics_(
  GMMirroredBuf* const matM_IJ_buf,
  GMMirroredBuf* const matM_JK_buf,
  GMMirroredBuf* const matM_KIK_buf,
  GMMirroredBuf* const matB_buf,
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

  GMInsist(vector_sums_i && vector_sums_j && vector_sums_k);

  const bool is_part3 = si->is_part3;

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

  //const size_t nvl64 = (size_t)nvl;
  //const size_t I_max64 = (size_t)I_max;

  /*--------------------*/
  /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
  /*--------------------*/

  /*----------*/
  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK &&
      !GMEnv_all2all(env)) {
    /*----------*/

#pragma omp parallel for collapse(2)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = GMMirroredBuf_elt_const<GMFloat>(matM_IJ_buf, I, J);;
        const GMFloat min_JK = GMMirroredBuf_elt_const<GMFloat>(matM_JK_buf, J, K);;
        const GMFloat min_KIK = GMMirroredBuf_elt_const<GMFloat>(matM_KIK_buf, K, I);;
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = GMMirroredBuf_elt_const<GMFloat>(matB_buf, I, K);;
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
  } else if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEK &&
      GMEnv_all2all(env)) {
    /*----------*/

    GMIndexCache index_cache = {0};
#pragma omp parallel for collapse(2) firstprivate(index_cache)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = GMMirroredBuf_elt_const<GMFloat>(matM_IJ_buf, I, J);;
        const GMFloat min_JK = GMMirroredBuf_elt_const<GMFloat>(matM_JK_buf, J, K);;
        const GMFloat min_KIK = is_part3 ?
          GMMirroredBuf_elt_const<GMFloat>(matM_KIK_buf, K, I) :
          GMMirroredBuf_elt_const<GMFloat>(matM_KIK_buf, I, K);
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = GMMirroredBuf_elt_const<GMFloat>(matB_buf, I, K);;
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
//TODO: try to merge this case with the one below it
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

        const GMTally2x2 mB = GMMirroredBuf_elt_const<GMTally2x2>(matB_buf, I, K);;
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

        numer.data[0] = GMTally1_encode(r000, r001);
        numer.data[1] = GMTally1_encode(r010, r011);
        numer.data[2] = GMTally1_encode(r100, r101);
        numer.data[3] = GMTally1_encode(r110, r111);
        GMMetrics_tally4x2_set_3(metrics, i, j, k, numer, env);
        if (step_2way==2) {
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
        }
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

        const GMTally2x2 mB = GMMirroredBuf_elt_const<GMTally2x2>(matB_buf, I, K);;
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

        numer.data[0] = GMTally1_encode(r000_permuted, r001_permuted);
        numer.data[1] = GMTally1_encode(r010_permuted, r011_permuted);
        numer.data[2] = GMTally1_encode(r100_permuted, r101_permuted);
        numer.data[3] = GMTally1_encode(r110_permuted, r111_permuted);
        GMMetrics_tally4x2_set_all2all_3_permuted_cache(metrics, I, J, K,
                                   j_block, k_block, numer, &index_cache, env);
        if (step_2way==2) {
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
        }
      } /*---for K---*/
    }   /*---for I---*/
    if (step_2way == 2) {
      metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                          (K_max - K_min);
    }
    /*----------*/
  } else {
    /*----------*/
    GMInsist(false);
    /*----------*/
  } /*---GMEnv_metric_type(env)---*/
  /*----------*/
}

//=============================================================================

static void lock(bool& lock_val) {
  GMInsist(! lock_val);
  lock_val = true;
};

static void unlock(bool& lock_val) {
  GMInsist(lock_val);
  lock_val = false;
};

//=============================================================================
/*---Start calculation of numerators, 3-way gpu---*/

void gm_compute_numerators_3way_gpu_start_(
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
  GMInsist(this_ && metrics && env);
  GMInsist(vectors_i && vectors_j && vectors_k);
  GMInsist(vectors_i_buf && vectors_j_buf && vectors_k_buf);
  GMInsist(j_block >= 0 && j_block < GMEnv_num_block_vector(env));
  GMInsist(k_block >= 0 && k_block < GMEnv_num_block_vector(env));
  GMInsist(!(GMEnv_proc_num_vector_i(env) == j_block &&
                   GMEnv_proc_num_vector_i(env) != k_block));
  GMInsist(!(GMEnv_proc_num_vector_i(env) == k_block &&
                   GMEnv_proc_num_vector_i(env) != j_block));
  GMInsist(GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU);
  GMInsist(GMEnv_num_way(env) == GM_NUM_WAY_3);
  GMInsist(vector_sums_i && vector_sums_j && vector_sums_k);

  /*---Initializations---*/

  const int nvl = metrics->num_vector_local;
  const int npvfl = vectors_i->num_packedval_field_local;

  const int i_block = GMEnv_proc_num_vector_i(env);

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  const bool need_mat_ij = env->need_2way;
  const bool need_mat_jk = env->need_2way && !si->is_part1;
  const bool need_mat_kik = env->need_2way && si->is_part3;

  /*----------------------------------------*/
  /*---First get the required 2-way ij, jk, ik metrics---*/
  /*----------------------------------------*/

  /*--------------------*/
  /*---Compute i_block - j_block PROD---*/
  /*--------------------*/

  GMMirroredBuf* tmp_buf[2] = {&this_->tmp_buf[0],
                                   &this_->tmp_buf[1]};

  GMMirroredBuf* const matM_ij_buf = need_mat_ij ? &this_->matM_ij_buf :
                                                      NULL;

  if (need_mat_ij) {
    GMMirroredBuf* matM_ij_buf_ptr =
       env->do_reduce ? tmp_buf[0] : matM_ij_buf;

    gm_linalg_set_matrix_zero_start(matM_ij_buf_ptr, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_i_buf->d, npvfl,
                         vectors_j_buf->d, npvfl,
                         matM_ij_buf_ptr->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_ij_buf_ptr, env);
    gm_get_metrics_wait(metrics, matM_ij_buf_ptr, env);

    if (env->do_reduce) {
      gm_reduce_metrics(metrics, matM_ij_buf, matM_ij_buf_ptr, env);
    }
  }

  /*--------------------*/
  /*---Compute j_block - k_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  GMMirroredBuf* const matM_jk_buf =
      !si->is_part1 ? &this_->matM_jk_buf : matM_ij_buf;

  if (need_mat_jk) {
    GMMirroredBuf* matM_jk_buf_ptr =
        env->do_reduce ? tmp_buf[0] : matM_jk_buf;

    gm_linalg_set_matrix_zero_start(matM_jk_buf_ptr, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_j_buf->d, npvfl,
                         vectors_k_buf->d, npvfl,
                         matM_jk_buf_ptr->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_jk_buf_ptr, env);
    gm_get_metrics_wait(metrics, matM_jk_buf_ptr, env);

    if (env->do_reduce) {
      gm_reduce_metrics(metrics, matM_jk_buf, matM_jk_buf_ptr, env);
    }
  }

  /*--------------------*/
  /*---Compute k_block - i_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  /*---NOTE: for Part 3, this is indexed directly as (k,i).
       Otherwise, it is indexed through an alias as (i,k)---*/

  GMMirroredBuf* const matM_kik_buf = si->is_part3
    ? &this_->matM_kik_buf : matM_ij_buf;

  if (need_mat_kik) {
    GMMirroredBuf* matM_kik_buf_ptr =
        env->do_reduce ? tmp_buf[0] : matM_kik_buf;

    gm_linalg_set_matrix_zero_start(matM_kik_buf_ptr, env);

    gm_linalg_gemm_start(nvl, nvl, npvfl,
                         vectors_k_buf->d, npvfl,
                         vectors_i_buf->d, npvfl,
                         matM_kik_buf_ptr->d, nvl, env);
    gm_compute_wait(env);

    gm_get_metrics_start(metrics, matM_kik_buf_ptr, env);
    gm_get_metrics_wait(metrics, matM_kik_buf_ptr, env);

    if (env->do_reduce) {
      gm_reduce_metrics(metrics, matM_kik_buf, matM_kik_buf_ptr, env);
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

  GMMirroredBuf* matX_buf[2] = {&this_->matX_buf[0], &this_->matX_buf[1]};
  GMMirroredBuf* matB_buf[2] = {&this_->matB_buf[0], &this_->matB_buf[1]};

  /*---Set up pointers to permute the access of axes for Part 3---*/
  /*---We use capitals I, J, K here to denote the PERMUTED axes---*/

  const bool is_ijk = !si->is_part3 ? true : si->sax1;
  const bool is_kij = !si->is_part3 ? false : si->sax0;
  const bool is_jki = !si->is_part3 ? false : si->sax2;

  /* clang-format off */
  GMMirroredBuf* const vectors_I_buf = is_ijk ? vectors_i_buf :
                                           is_kij ? vectors_k_buf :
                                           is_jki ? vectors_j_buf : 0;
 
  GMMirroredBuf* const vectors_J_buf = is_ijk ? vectors_j_buf :
                                           is_kij ? vectors_i_buf :
                                           is_jki ? vectors_k_buf : 0;
 
  GMMirroredBuf* const vectors_K_buf = is_ijk ? vectors_k_buf :
                                           is_kij ? vectors_j_buf :
                                           is_jki ? vectors_i_buf : 0;
  
  //TODO - use is_ijk etc. 
  GMMirroredBuf* const matM_IJ_buf  = !si->is_part3 ? matM_ij_buf  :
                                               si->sax0 ? matM_kik_buf :
                                               si->sax1 ? matM_ij_buf  :
                                               si->sax2 ? matM_jk_buf  : 0;
  
  GMMirroredBuf* const matM_JK_buf  = !si->is_part3 ? matM_jk_buf  :
                                               si->sax0 ? matM_ij_buf  :
                                               si->sax1 ? matM_jk_buf  :
                                               si->sax2 ? matM_kik_buf : 0;
  
  GMMirroredBuf* const matM_KIK_buf = !si->is_part3 ? matM_kik_buf :
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

  MPI_Request mpi_requests[2] = {0, 0};

  typedef struct {
    int step_num;
    int step_2way;
    int J;
    int I_min;
    int I_max;
    int K_min;
    int K_max;
    bool empty;
    bool is_compute_step;
    bool do_compute;
    int index_01;
    GMMirroredBuf matB_buf;
    GMMirroredBuf tmp_buf;
  } LoopVars;

  LoopVars vars = {0};
  LoopVars vars_prev = {0};
  LoopVars vars_prevprev = {0};
  LoopVars vars_next = {0};
  vars.do_compute = false;
  vars_prev.do_compute = false;
  vars_prevprev.do_compute = false;
  vars_next.do_compute = false;

  // Use locks to verify no race condition on a buffer.
  // Lock buffer when in use for read or write, unlock when done.
  bool lock_tmp_buf_h[2] = {false, false};
  bool lock_tmp_buf_d[2] = {false, false};
  bool lock_matX_buf_h[2] = {false, false};
  bool lock_matX_buf_d[2] = {false, false};
  bool lock_matB_buf_h[2] = {false, false};
  bool lock_matB_buf_d[2] = {false, false};

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
    if (vars_next.I_max <= nvl) {
      // Create buffer aliases with required shape.
      if (env->do_reduce) {
        GMMirroredBuf_create(&vars_next.tmp_buf,
                                 tmp_buf[vars_next.index_01],
                                 vars_next.I_max, env);
      }
      GMMirroredBuf_create(&(vars_next.matB_buf),
                               matB_buf[vars_next.index_01],
                               vars_next.I_max, env);
    }

    GMMirroredBuf* matB_buf_ptr_prev = env->do_reduce ?
                                           &vars_prev.tmp_buf :
                                           &vars_prev.matB_buf;
    GMMirroredBuf* matB_buf_ptr = env->do_reduce ? &vars.tmp_buf :
                                                       &vars.matB_buf;

    bool& lock_matB_buf_ptr_h_prevprev = env->do_reduce ?
                                   lock_tmp_buf_h[vars_prevprev.index_01] :
                                   lock_matB_buf_h[vars_prevprev.index_01];
    bool& lock_matB_buf_ptr_h_prev = env->do_reduce ?
                                   lock_tmp_buf_h[vars_prev.index_01] :
                                   lock_matB_buf_h[vars_prev.index_01];
    bool& lock_matB_buf_ptr_d_prev = env->do_reduce ?
                                   lock_tmp_buf_d[vars_prev.index_01] :
                                   lock_matB_buf_d[vars_prev.index_01];
    bool& lock_matB_buf_ptr_d = env->do_reduce ?
                                   lock_tmp_buf_d[vars.index_01] :
                                   lock_matB_buf_d[vars.index_01];

    //==========

    if (vars_next.do_compute) {
      /*---Populate leading columns of matX---*/
      lock(lock_matX_buf_h[vars_next.index_01]);
      gm_compute_numerators_3way_gpu_form_matX_(vectors_i,
          vectors_I_buf, vectors_J_buf, matX_buf[vars_next.index_01],
          vars_next.J, vars_next.step_2way,
          vars_next.I_min, vars_next.I_max, env);
      unlock(lock_matX_buf_h[vars_next.index_01]);
    }

    //==========

    if (vars.do_compute) {
      /*---Send matrix matX to GPU - WAIT---*/
      gm_linalg_set_matrix_wait(env);
      unlock(lock_matX_buf_h[vars.index_01]);
      unlock(lock_matX_buf_d[vars.index_01]);
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Perform pseudo GEMM matB = matX^T PROD V - WAIT---*/
      gm_compute_wait(env);
      unlock(lock_matB_buf_ptr_d_prev);
      unlock(lock_matX_buf_d[vars_prev.index_01]);
    }

    //==========

    if (vars_next.do_compute) {
      /*---Send matrix matX to GPU - START---*/
      lock(lock_matX_buf_h[vars_next.index_01]);
      lock(lock_matX_buf_d[vars_next.index_01]);
      gm_linalg_set_matrix_start(matX_buf[vars_next.index_01], env);
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - START---*/
      lock(lock_matB_buf_ptr_d_prev);
      lock(lock_matB_buf_ptr_h_prev);
      gm_linalg_get_matrix_start(matB_buf_ptr_prev, env);
    }

    //==========

    if (vars.do_compute) {
      /*---Initialize result matrix to zero (apparently magma requires)---*/
      lock(lock_matB_buf_ptr_d);
      lock(lock_matX_buf_d[vars.index_01]);
      gm_linalg_set_matrix_zero_start(matB_buf_ptr, env);
      /*---Perform pseudo GEMM matB = matX^T PROD V - START---*/
      gm_linalg_gemm_start(vars.I_max, nvl, npvfl,
                           matX_buf[vars.index_01]->d, npvfl,
                           vectors_K_buf->d, npvfl,
                           matB_buf_ptr->d, vars.I_max, env);
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - WAIT---*/
      gm_linalg_get_matrix_wait(env);
      if (vars_prev.step_2way == 0) {
        gm_metrics_gpu_adjust(metrics, matB_buf_ptr_prev, env);
      }
      unlock(lock_matB_buf_ptr_d_prev);
      unlock(lock_matB_buf_ptr_h_prev);
    }

    //==========

    if (vars_prevprev.do_compute && env->do_reduce) {
      /*---Reduce along field procs - WAIT---*/
      gm_reduce_metrics_wait(&(mpi_requests[vars_prevprev.index_01]), env); 
      unlock(lock_matB_buf_ptr_h_prevprev);
      unlock(lock_matB_buf_h[vars_prevprev.index_01]);
    }

    //==========

    if (vars_prev.do_compute && env->do_reduce) {
      /*---Reduce along field procs - START---*/
      lock(lock_matB_buf_ptr_h_prev);
      lock(lock_matB_buf_h[vars_prev.index_01]);
      mpi_requests[vars_prev.index_01] = gm_reduce_metrics_start(metrics,
          &vars_prev.matB_buf, matB_buf_ptr_prev, env);
    }

    //==========

    //---NOTE: matB_buf[vars_prevprev.index_01]->d is locked now
    //---but matB_buf[vars_prevprev.index_01]->h is usable.

    if (vars_prevprev.do_compute) {
      /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
      lock(lock_matB_buf_ptr_h_prevprev);
      gm_compute_numerators_3way_gpu_form_metrics_(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
          &vars_prevprev.matB_buf,
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
      unlock(lock_matB_buf_ptr_h_prevprev);
    }

    //==========

  } /*---for step_num---*/

  /*--------------------*/
  /*---Free memory---*/
  /*--------------------*/

  GMSectionInfo_destroy(si, env);
}

//=============================================================================

//-----------------------------------------------------------------------------
