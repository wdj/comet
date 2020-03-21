//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block_gpu.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block, GPU case.
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
/// \brief Create matrix X from a vectors object and a single vector column.

static void compute_metrics_3way_block_linalg_form_matXitem_(
  const MirroredBuf* const vectors_I_buf,
  const MirroredBuf* const vectors_J_buf,
  MirroredBuf* const matXitem_buf,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  const int npvfl,
  CEnv& env) {

  matXitem_buf->lock_h();

  if (env.metric_type() == MetricType::CZEK) {

    // Populate leading columns of matX.

    // Don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
    #pragma omp parallel for schedule(dynamic,1000)
    for (int I = I_min; I < I_max; ++I) {
      // Operate on columns x_i and x_j elementwise.
      for (int f = 0; f < npvfl; ++f) {
        const GMFloat a = vectors_I_buf->elt_const<GMFloat>(f, I);
        const GMFloat b = vectors_J_buf->elt_const<GMFloat>(f, J);
        matXitem_buf->elt<GMFloat>(f, I) = a < b ? a : b;
      }  //---for f---//
    }    //---for I---//

  } else if (env.is_metric_type_bitwise()) {

    // Extract column J of vectors_J, later use to form matX.

    if (env.form_matX_on_accel()) {

      for (int word = 0; word<2; ++word) {
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

            matXitem_buf->elt<GMBits2x64>(pvfl, 0).data[word] =
              vectors_J_buf->elt_const<GMBits2x64>(pvfl, J).data[word];

        }
      }

    } else { // if (!env.form_matX_on_accel())

      COMET_INSIST(env.metric_type() == MetricType::CCC &&
                   "Case currently unimplemented.");

      // Populate leading columns of matX.

      for (int I = I_min; I < I_max; ++I) {

        // Mask for odd bits (starting at lowest-order bit: bit 0, bit 2, ...)

        const uint64_t oddbits = 0x5555555555555555;

        // Operate on columns v_i and v_j elementwise.
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

          const bool sparse = env.sparse();

          for (int word = 0; word<2; ++word) {
            const uint64_t vI = vectors_I_buf->elt_const<GMBits2x64>(
                                             pvfl, I).data[word];
            const uint64_t vJ = vectors_J_buf->elt_const<GMBits2x64>(
                                             pvfl, J).data[word];

            // Create word whose odd bits sample the lo (denoted here "..._0")
            // or hi ("..._1") bit of the seminibble.  Also create the
            // complement thereof (denoted "n...").

            const uint64_t  vI_0 =   vI        & oddbits;
            const uint64_t  vI_1 =  (vI >> 1)  & oddbits;
            const uint64_t nvI_0 = ~ vI        & oddbits;
            const uint64_t nvI_1 = ~(vI >> 1)  & oddbits;

            const uint64_t  vJ_0 =   vJ        & oddbits;
            const uint64_t  vJ_1 =  (vJ  >> 1) & oddbits;

            // Create a mask whose odd bits denote whether each respective
            // seminibble of vector I matches the case we are handling
            // in this 2-way step (and the complement thereof).
            // step 0: select (only) entries equal to 00
            // step 1: select (only) entries equal to 01 or 10 (nonsparse case)
            // step 1: select (only) entries equal to 01 (sparse) (ignore 10)
            // step 2: select (only) entries equal to 11
            // Note here that 10 is in some situations used as a special marker
            // meaning, ignore this seminiblle for the calculations.

            const uint64_t  vI_mask =
              step_2way==0 ?  nvI_0 & nvI_1  & oddbits : // 00
              step_2way==1 && sparse ?
                             ( vI_0 & nvI_1) & oddbits : // 01
              step_2way==1 ? ( vI_0 ^  vI_1) & oddbits : // 01, 10
            /*step_2way==2*/   vI_0 &  vI_1  & oddbits;  // 11

            const uint64_t nvI_mask =
              step_2way==0 ? ( vI_0 |  vI_1) & oddbits :
              step_2way==1 && sparse ?
                           (nvI_0 |  vI_1) & oddbits :
              step_2way==1 ? ( vI_0 ^ nvI_1) & oddbits :
            /*step_2way==2*/ (nvI_0 | nvI_1) & oddbits;

            // Construct the lo and hi bit of the result seminibble of matrix X.
            // This is best understood by looking at truth table (see paper).
            // case vI_mask = 1: vJ = 00 => X = 00
            // case vI_mask = 1: vJ = 01 => X = 01
            // case vI_mask = 1: vJ = 10 => X = 01 (nonsparse case)
            // case vI_mask = 1: vJ = 10 => X = 10 (sparse case)
            // case vI_mask = 1: vJ = 11 => X = 11
            // case vI_mask = 0: X = 10
            // Thus for the nonsparse case:
            //  - lo bit is 1 (11 or 01) if vJ is 01, 10 or 11 and vI is the
            //    case being handled for this 2-way step.
            //  - hi bit is 1 (10 or 11) if vJ is 11 or if vI is a case not
            //    being handled for this 2-way step.

            const uint64_t r_0 =  vI_mask & (sparse ? vJ_0 : vJ_0 | vJ_1);
            const uint64_t r_1 = nvI_mask | (sparse ? vJ_1 : vJ_0 & vJ_1);

            // Combine even and odd bits

            const uint64_t r = r_0 | (r_1 << 1);

            // Store result

            matXitem_buf->elt<GMBits2x64>(pvfl, I).data[word] = r;
          } // word
        }  // f
      }    // I

    } // if (env.form_matX_on_accel())

  } else {

    COMET_INSIST(false);

  } // env.metric_type()

  matXitem_buf->unlock_h();
}

//-----------------------------------------------------------------------------
/// \brief Set metrics numerators and denominators using GEMM results.

static void compute_metrics_3way_block_linalg_form_metrics_(
  MirroredBuf* const matM_IJ_buf,
  MirroredBuf* const matM_JK_buf,
  MirroredBuf* const matM_KIK_buf,
  MirroredBuf* const matB_buf,
  GMMetrics* metrics,
  int nvl, int J, int step_2way,
  int I_min, int I_max, int K_min, int K_max,
  int I_max_dim,
  int j_block, int k_block,
  GMSectionInfo* const si,
  VectorSums* vs_i, VectorSums* vs_j, VectorSums* vs_k,
  CEnv& env) {

  COMET_INSIST(vs_i && vs_j && vs_k);

  COMET_INSIST( ! (env.is_bitwise_3way_2step() && !env.form_matX_on_accel()) &&
               "Case currently unimplemented.");

  matB_buf->lock_h();

  const VectorSums* const vs_I = si->perm0(vs_i, vs_j, vs_k);
  const VectorSums* const vs_J = si->perm1(vs_i, vs_j, vs_k);
  const VectorSums* const vs_K = si->perm2(vs_i, vs_j, vs_k);

  //const size_t nvl64 = (size_t)nvl;
  //const size_t I_max64 = (size_t)I_max;

  //--------------------
  // Compute numerators using ijk piece and (if needed) 2-way pieces.
  //--------------------

  if (env.metric_type() == MetricType::CZEK && !env.all2all()) {

    // don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
    #pragma omp parallel for schedule(dynamic,1000)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = matM_IJ_buf->elt_const<GMFloat>(I, J);;
        const GMFloat min_JK = matM_JK_buf->elt_const<GMFloat>(J, K);;
        const GMFloat min_KIK = matM_KIK_buf->elt_const<GMFloat>(K, I);;
        // sum of mins vectors i, j, and k is matB(k,i).
        const GMFloat min_IJK = matB_buf->elt_const<GMFloat>(I, K);;
        const GMFloat numer = min_IJ + min_JK + min_KIK - min_IJK;
        const int i = I;
        const int j = J;
        const int k = K;
        // Make arithmetic order-independent.
        GMFloat smin, smid, smax;
        const auto si = vs_i->sum(i);
        const auto sj = vs_i->sum(j);
        const auto sk = vs_i->sum(k);
        utils::sort_3(smin, smid, smax, si, sj, sk);
        const GMFloat denom = smin + smid + smax;
        const GMFloat value = ((GMFloat)1.5) * numer / denom;
        GMMetrics_float_set_3(metrics, i, j, k, value, &env);
      } // K
    }   // I
    metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                        (K_max - K_min);

  } else if (env.metric_type() == MetricType::CZEK && env.all2all()) {

    GMIndexCache index_cache = {};
    // don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) firstprivate(index_cache) schedule(dynamic,1000)
    #pragma omp parallel for firstprivate(index_cache) schedule(dynamic,1000)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = matM_IJ_buf->elt_const<GMFloat>(I, J);;
        const GMFloat min_JK = matM_JK_buf->elt_const<GMFloat>(J, K);;
        const GMFloat min_KIK = si->is_part3 ?
          matM_KIK_buf->elt_const<GMFloat>(K, I) :
          matM_KIK_buf->elt_const<GMFloat>(I, K);
        // sum of mins vectors i, j, and k is matB(k,i).
        const GMFloat min_IJK = matB_buf->elt_const<GMFloat>(I, K);;
        const GMFloat numer = min_IJ + min_JK + min_KIK - min_IJK;
        // Make arithmetic order-independent.
        GMFloat smin, smid, smax;
        const auto sI = vs_I->sum(I);
        const auto sJ = vs_J->sum(J);
        const auto sK = vs_K->sum(K);
        utils::sort_3(smin, smid, smax, sI, sJ, sK);
        const GMFloat denom = smin + smid + smax;
        const GMFloat value = ((GMFloat)1.5) * numer / denom;
        GMMetrics_float_set_all2all_3_permuted_cache(metrics, I, J, K,
            j_block, k_block, value, &index_cache, &env);
      } // K
    }   // I
    metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                        (K_max - K_min);

  } else if (env.is_metric_type_bitwise()) {

    COMET_INSIST((env.metric_type() == MetricType::CCC ||
                  env.is_bitwise_3way_2step()) &&
                 "Case currently unimplemented.");

    const int nvle = I_max_dim;
    const int nvleD2 = nvle / 2;
    const int is_halved = env.is_bitwise_3way_2step() &&
                          env.is_vectors_halved();

    GMIndexCache index_cache = {};

    // don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) firstprivate(index_cache) schedule(dynamic,1000)
    #pragma omp parallel for firstprivate(index_cache) schedule(dynamic,1000)
    for (int K = K_min; K < K_max; ++K) {

      for (int I = I_min; I < I_max; ++I) {
        /*---For the permuted case,
         1) pay attention to KIK access
         2) swap 01 and 10 if needed.
        ---*/

        const bool is_I_in_range = !is_halved ? true :
          I >= nvleD2*step_2way && I < nvleD2*(step_2way+1);

        const bool init_numer = !is_halved ? 0 == step_2way : true;

        const int I0 = !is_halved ? I : I % nvleD2;
        const int I1 = !is_halved ? I : I % nvleD2 + nvleD2;
//        const int is_I0_in_range = I0 >= I_min && I0 < I_max;
//        const int is_I1_in_range = I1 >= I_min && I1 < I_max;

        const int i = si->unperm0(I, J, K);
        const int j = si->unperm1(I, J, K);
        const int k = si->unperm2(I, J, K);

        // Numerator.

if (is_I_in_range) {

        // NOTE: numer is accessed through a permuted index, but
        // the 001 etc. entries of numer are unpermuted.

        GMTally4x2 numer = init_numer ? GMTally4x2_null() :
          env.all2all() ?
          GMMetrics_tally4x2_get_all2all_3_permuted_cache(metrics, I, J, K,
                                        j_block, k_block, &index_cache, &env) :
          GMMetrics_tally4x2_get_3(metrics, i, j, k, &env);

        GMTally1 r000, r001, r010, r011, r100, r101, r110, r111;
        PackedDouble_decode(r000, r001, numer.data[0]);
        PackedDouble_decode(r010, r011, numer.data[1]);
        PackedDouble_decode(r100, r101, numer.data[2]);
        PackedDouble_decode(r110, r111, numer.data[3]);

        // NOTE: matB is generated by a GEMM from permuted vectors objects,
        // thus has permuted 001 etc. entries.

        const GMTally2x2 matB_permuted =
          matB_buf->elt_const<GMTally2x2>(I0, K);
// !is_I0_in_range ?  GMTally2x2_null() :
        GMTally1 matB00_permuted, matB01_permuted;
        GMTally1 matB10_permuted, matB11_permuted;
        PackedDouble_decode(matB00_permuted, matB01_permuted,
                            matB_permuted.data[0]);
        PackedDouble_decode(matB10_permuted, matB11_permuted,
                            matB_permuted.data[1]);

        int r000_permuted = r000;
        int r100_permuted = si->perm0(r100, r010, r001);
        int r010_permuted = si->perm1(r100, r010, r001);
        int r001_permuted = si->perm2(r100, r010, r001);
        int r011_permuted = si->perm0(r011, r101, r110);
        int r101_permuted = si->perm1(r011, r101, r110);
        int r110_permuted = si->perm2(r011, r101, r110);
        int r111_permuted = r111;

        // Add contribution from this sub-step.

        if (is_halved) {

          if (is_I_in_range) {

//printf("1--- %i %i %i    %i %i %i %i\n", I, J, K, matB00_permuted, matB01_permuted, matB10_permuted, matB11_permuted);
            r000_permuted += matB00_permuted;
            r001_permuted += matB01_permuted;
            r010_permuted += matB10_permuted;
            r011_permuted += matB11_permuted;

//printf("Hey1 %i %i %i   %i   %i   %i %i   %i \n", I, J, K, nvl, init_numer, I_min, I_max, step_2way);
            const GMTally2x2 matB_permuted =
              matB_buf->elt_const<GMTally2x2>(I1, K);
//printf("Hey2\n");
// !is_I1_in_range ?  GMTally2x2_null() :
            GMTally1 matB00_permuted, matB01_permuted;
            GMTally1 matB10_permuted, matB11_permuted;
            PackedDouble_decode(matB00_permuted, matB01_permuted,
                                matB_permuted.data[0]);
            PackedDouble_decode(matB10_permuted, matB11_permuted,
                                matB_permuted.data[1]);
//printf("2--- %i %i %i    %i %i %i %i\n", I, J, K, matB00_permuted, matB01_permuted, matB10_permuted, matB11_permuted);

            r100_permuted += matB00_permuted;
            r101_permuted += matB01_permuted;
            r110_permuted += matB10_permuted;
            r111_permuted += matB11_permuted;

          }

        } else if (env.is_bitwise_3way_2step()) {

          if (step_2way==0) {
//printf("1--- %i %i %i    %i %i %i %i\n", I, J, K, matB00_permuted, matB01_permuted, matB10_permuted, matB11_permuted);
            r000_permuted += matB00_permuted;
            r001_permuted += matB01_permuted;
            r010_permuted += matB10_permuted;
            r011_permuted += matB11_permuted;
          } else /*---step_2way==1---*/ {
//printf("2--- %i %i %i    %i %i %i %i\n", I, J, K, matB00_permuted, matB01_permuted, matB10_permuted, matB11_permuted);
            r100_permuted += matB00_permuted;
            r101_permuted += matB01_permuted;
            r110_permuted += matB10_permuted;
            r111_permuted += matB11_permuted;
          }

        } else { // if (!env.is_bitwise_3way_2step())

          if (step_2way==0) {
            r000_permuted += 2 * matB00_permuted;
            r001_permuted += 2 * matB01_permuted;
            r010_permuted += 2 * matB10_permuted;
            r011_permuted += 2 * matB11_permuted;
          } else if (step_2way==1) {
            r000_permuted += matB00_permuted;
            r001_permuted += matB01_permuted;
            r010_permuted += matB10_permuted;
            r011_permuted += matB11_permuted;
            r100_permuted += matB00_permuted;
            r101_permuted += matB01_permuted;
            r110_permuted += matB10_permuted;
            r111_permuted += matB11_permuted;
          } else /*---step_2way==2---*/ {
            r100_permuted += 2 * matB00_permuted;
            r101_permuted += 2 * matB01_permuted;
            r110_permuted += 2 * matB10_permuted;
            r111_permuted += 2 * matB11_permuted;
          }

        } // if (env.is_bitwise_3way_2step())

        r000 = r000_permuted;
        r100 = si->unperm0(r100_permuted, r010_permuted, r001_permuted);
        r010 = si->unperm1(r100_permuted, r010_permuted, r001_permuted);
        r001 = si->unperm2(r100_permuted, r010_permuted, r001_permuted);
        r011 = si->unperm0(r011_permuted, r101_permuted, r110_permuted);
        r101 = si->unperm1(r011_permuted, r101_permuted, r110_permuted);
        r110 = si->unperm2(r011_permuted, r101_permuted, r110_permuted);
        r111 = r111_permuted;
//printf("        %i %i %i    %i %i %i %i %i %i %i %i\n", I, J, K, r000, r001, r010, r011, r100, r101, r110, r111);

        numer.data[0] = GMTally1_encode(r000, r001);
        numer.data[1] = GMTally1_encode(r010, r011);
        numer.data[2] = GMTally1_encode(r100, r101);
        numer.data[3] = GMTally1_encode(r110, r111);
        if (env.all2all()) {
          GMMetrics_tally4x2_set_all2all_3_permuted_cache(metrics, I, J, K,
                                  j_block, k_block, numer, &index_cache, &env);
        } else {
          GMMetrics_tally4x2_set_3(metrics, i, j, k, numer, &env);
        }

} // if (is_I_in_range)

        // Denominator.

        if (step_2way == env.num_step_2way_for_3way() - 1) {
          const auto si1 = (GMTally1)vs_i->sum(i);
          const auto sj1 = (GMTally1)vs_j->sum(j);
          const auto sk1 = (GMTally1)vs_k->sum(k);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);
          if (env.all2all()) {
            GMMetrics_float3_S_set_all2all_3_permuted_cache(metrics, I, J, K,
              j_block, k_block, si1_sj1_sk1, &index_cache, &env);
          } else {
            GMMetrics_float3_S_set_3(metrics, i, j, k, si1_sj1_sk1, &env);
          }
          if (env.sparse()) {
            const auto ci = (GMTally1)vs_i->count(i);
            const auto cj = (GMTally1)vs_j->count(j);
            const auto ck = (GMTally1)vs_k->count(k);
            const GMFloat3 ci_cj_ck = GMFloat3_encode(ci, cj, ck);
            if (env.all2all()) {
              GMMetrics_float3_C_set_all2all_3_permuted_cache(metrics, I, J, K,
                j_block, k_block, ci_cj_ck, &index_cache, &env);
            } else {
              GMMetrics_float3_C_set_3(metrics, i, j, k, ci_cj_ck, &env);
            }
          } /*---if sparse---*/

        }
      } // K
    }   // I
    if (step_2way == env.num_step_2way_for_3way() - 1) {
      metrics->num_elts_local_computed += (I_max - I_min) * (size_t)
                                          (K_max - K_min);
    }

  } else { // if (env.metric_type() ...

    COMET_INSIST(false);

  } // if (env.metric_type() ...

  matB_buf->unlock_h();
}

//-----------------------------------------------------------------------------
/// \brief Compute 3-way metrics for cases that use the linalg package.

void ComputeMetrics3WayBlock::compute_linalg_(
  VData vdata_i, VData vdata_j, VData vdata_k, GMMetrics& metrics,
  int j_block, int k_block, int section_step) {

  COMET_INSIST(j_block >= 0 && j_block < env_.num_block_vector());
  COMET_INSIST(k_block >= 0 && k_block < env_.num_block_vector());
  COMET_INSIST(! (env_.proc_num_vector() == j_block &&
              env_.proc_num_vector() != k_block));
  COMET_INSIST(! (env_.proc_num_vector() == k_block &&
              env_.proc_num_vector() != j_block));
  COMET_INSIST(env_.is_using_linalg());
  COMET_INSIST(env_.num_way() == NUM_WAY::_3);

  // Initializations.

  const int nvl = metrics.num_vector_local;
  const int npvfl = vdata_i.vectors->num_packedval_field_local;
  const int i_block = env_.proc_num_vector();

  GMSectionInfo si_value, *si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step, nvl, &env_);

  const bool need_mat_ij = env_.does_3way_need_2way();
  const bool need_mat_jk = env_.does_3way_need_2way() && ! si->is_part1;
  const bool need_mat_kik = env_.does_3way_need_2way() && si->is_part3;

  GMDecompMgr* const dm = vdata_i.vectors->dm;

  //----------------------------------------
  // First get the required 2-way ij, jk, ik metrics.
  //----------------------------------------

  //--------------------
  // Compute i_block - j_block PROD.
  //--------------------

//CHANGE: remove this?
  MirroredBuf* tmp_buf[NUM_BUF] = {&tmp_buf_[0], &tmp_buf_[1]};

  MirroredBuf* const matM_ij_buf = need_mat_ij ? &matM_ij_buf_ : NULL;

  if (need_mat_ij) {
    MirroredBuf* matM_ij_buf_ptr = env_.do_reduce()
      ? tmp_buf[0] : matM_ij_buf;

    gm_linalg_gemm(nvl, nvl, npvfl,
                   vdata_i.buf, vdata_j.buf, matM_ij_buf_ptr,
                   vdata_i.sums->sums(), vdata_j.sums->sums(),
                   vdata_i.sums->counts(), vdata_j.sums->counts(),
                   dm, &env_);

    matM_ij_buf_ptr->from_accel();

    gm_reduce_metrics(&metrics, matM_ij_buf, matM_ij_buf_ptr, &env_);
  }

  //--------------------
  // Compute j_block - k_block PROD.
  //--------------------

  // Need to compute only if not identical to already computed values.

  MirroredBuf* const matM_jk_buf = ! si->is_part1
    ? &matM_jk_buf_ : matM_ij_buf;

  if (need_mat_jk) {
    MirroredBuf* matM_jk_buf_ptr = env_.do_reduce()
      ? tmp_buf[0] : matM_jk_buf;

    gm_linalg_gemm(nvl, nvl, npvfl,
                   vdata_j.buf, vdata_k.buf, matM_jk_buf_ptr,
                   vdata_j.sums->sums(), vdata_k.sums->sums(),
                   vdata_j.sums->counts(), vdata_k.sums->counts(),
                   dm, &env_);

    matM_jk_buf_ptr->from_accel();

    gm_reduce_metrics(&metrics, matM_jk_buf, matM_jk_buf_ptr, &env_);
  }

  //--------------------
  // Compute k_block - i_block PROD.
  //--------------------

  // Need to compute only if not identical to already computed values.

  // NOTE: for Part 3, this is indexed directly as (k,i).
  //  Otherwise, it is indexed through an alias as (i,k).

  MirroredBuf* const matM_kik_buf = si->is_part3
    ? &matM_kik_buf_ : matM_ij_buf;

  if (need_mat_kik) {
    MirroredBuf* matM_kik_buf_ptr =
        env_.do_reduce() ? tmp_buf[0] : matM_kik_buf;

    gm_linalg_gemm(nvl, nvl, npvfl,
                   vdata_k.buf, vdata_i.buf, matM_kik_buf_ptr,
                   vdata_k.sums->sums(), vdata_i.sums->sums(),
                   vdata_k.sums->counts(), vdata_i.sums->counts(),
                   dm, &env_);

    matM_kik_buf_ptr->from_accel();

    gm_reduce_metrics(&metrics, matM_kik_buf, matM_kik_buf_ptr, &env_);
  } /*---is_part3---*/

  //----------------------------------------
  // Now compute ijk piece, via an outer loop over j values.
  //----------------------------------------

  // Allocate magma CPU/GPU memory for matrices X and B.
  // X = elementwise OP of one vector with each of the rest of the vectors.
  // For the jth iteration, the ith column of X is the elementwise OP
  //   of vectors i and j.
  // B = X^T PROD V = three way PROD.

//CHANGE: these as class members, init in ctor?
  MirroredBuf* matXitem_buf[NUM_BUF] = {&matXitem_buf_[0], &matXitem_buf_[1]};
  MirroredBuf* matB_buf[NUM_BUF] = {&matB_buf_[0], &matB_buf_[1]};

  // Set up pointers to permute the access of axes for Part 3.
  // Use capitals I, J, K here to denote the PERMUTED axes.

  VectorSums* const vsums_I =
                        si->perm0(vdata_i.sums, vdata_j.sums, vdata_k.sums);
  VectorSums* const vsums_J =
                        si->perm1(vdata_i.sums, vdata_j.sums, vdata_k.sums);
  VectorSums* const vsums_K =
                        si->perm2(vdata_i.sums, vdata_j.sums, vdata_k.sums);

  MirroredBuf* const vectors_I_buf =
                        si->perm0(vdata_i.buf, vdata_j.buf, vdata_k.buf);
  MirroredBuf* const vectors_J_buf =
                        si->perm1(vdata_i.buf, vdata_j.buf, vdata_k.buf);
  MirroredBuf* const vectors_K_buf =
                        si->perm2(vdata_i.buf, vdata_j.buf, vdata_k.buf);

  MirroredBuf* const matM_IJ_buf =
                        si->perm0(matM_ij_buf, matM_jk_buf, matM_kik_buf);
  MirroredBuf* const matM_JK_buf =
                        si->perm1(matM_ij_buf, matM_jk_buf, matM_kik_buf);
  MirroredBuf* const matM_KIK_buf =
                        si->perm2(matM_ij_buf, matM_jk_buf, matM_kik_buf);

  if (env_.form_matX_on_accel()) {
    vectors_I_buf->to_accel();
  }

  //--------------------
  // Collapsed loops over J and over 2-way steps.
  //--------------------

  const int J_min = si->J_lb;
  const int J_max = si->J_ub;
  const int J_count = J_max - J_min;

  const int num_step_2way = env_.num_step_2way_for_3way();
  const int num_step = J_count * num_step_2way;
  const int extra_step = 1;

  MPI_Request mpi_requests[NUM_BUF] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL};

  struct LoopVars {
    int step_num;
    int step_2way;
    int J;
    int I_min;
    int I_max;
    int K_min;
    int K_max;
    int I_max_dim;
    bool empty;
    bool is_compute_step;
    bool do_compute;
    bool do_reduce;
    int index_01;
    MirroredBuf matB_buf;
    MirroredBuf tmp_buf;
    MirroredBuf* matB_buf_ptr() {return do_reduce ? &tmp_buf : &matB_buf;}
    LoopVars(CEnv& env)
      : do_compute(false)
      , do_reduce(env.do_reduce())
      , matB_buf(env)
      , tmp_buf(env) {}
    //LoopVars& operator=(LoopVars&) = default;
    //static void copy(LoopVars& w, const LoopVars* v) {
//      memcpy(&w, &v, sizeof(LoopVars));
//    }
  };

  const int num_buf = 4;

  std::vector<LoopVars> vars_buf;
  for (int i=0; i<num_buf; ++i)
    vars_buf.push_back(LoopVars(env_));

  //LoopVars vars(env_);
  //LoopVars vars_prev(env_);
  //LoopVars vars_prevprev(env_);
  //LoopVars vars_next(env_);

  const int first_step = 0 - extra_step;

  //========================================
  for (int step_num = first_step; step_num < num_step+extra_step*2;
       ++step_num) {
  //========================================

    // Set per-step variables.

    LoopVars& vars_prevprev = vars_buf[(step_num - first_step + 0) % num_buf];
    LoopVars& vars_prev = vars_buf[(step_num - first_step + 1) % num_buf];
    LoopVars& vars = vars_buf[(step_num - first_step + 2) % num_buf];
    LoopVars& vars_next = vars_buf[(step_num - first_step + 3) % num_buf];

    vars_next.step_num = step_num + 1;
    vars_next.step_2way = utils::mod_i(vars_next.step_num, num_step_2way);
    vars_next.J = J_min + utils::trunc(vars_next.step_num, num_step_2way);
    vars_next.I_min = 0;
    vars_next.I_max = si->is_part1 ? vars_next.J : nvl;
    vars_next.I_max_dim = tc_gemm_size_required(vars_next.I_max, env_);
    vars_next.K_min = si->is_part3 ? 0 : vars_next.J + 1;
    vars_next.K_max = nvl;
    vars_next.empty = vars_next.I_min >= vars_next.I_max ||
                      vars_next.K_min >= vars_next.K_max;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.do_compute = vars_next.is_compute_step && ! vars_next.empty;
    vars_next.index_01 = utils::mod_i(vars_next.step_num, (int)NUM_BUF);
    if (vars_next.I_max <= nvl) {
      COMET_INSIST(vars_next.I_max_dim <= nvl &&
               "Block size rounding-up error.");
      // Create buffer aliases with required shape.
      if (env_.do_reduce()) {
        vars_next.tmp_buf.allocate(*tmp_buf[vars_next.index_01],
                                   vars_next.I_max_dim);
      }
      vars_next.matB_buf.allocate(*matB_buf[vars_next.index_01],
                                  vars_next.I_max_dim);
    }

    //TODO: fix locks to work properly for CPU case.

    //========== Send matrix matXitem to GPU - WAIT

    if (vars.do_compute) {
      matXitem_buf[vars.index_01]->to_accel_wait();
    }

    //========== Perform pseudo GEMM matB = matX^T PROD V - WAIT

    if (vars_prev.do_compute) {
      gm_linalg_gemm_wait(vars_prev.I_max, nvl, npvfl,
          matXitem_buf[vars_prev.index_01], vectors_I_buf, vectors_K_buf,
          vars_prev.matB_buf_ptr(),
          vsums_I->sums(), vsums_J->sums(), vsums_K->sums(),
          vsums_I->counts(), vsums_J->counts(), vsums_K->counts(), vars_prev.J,
          vars_prev.step_2way, dm, &env_);
    }

    //========== Calculate matXitem

    if (vars_next.do_compute) {
      compute_metrics_3way_block_linalg_form_matXitem_(
          vectors_I_buf, vectors_J_buf, matXitem_buf[vars_next.index_01],
          vars_next.J, vars_next.step_2way,
          vars_next.I_min, vars_next.I_max, npvfl, env_);
    }

    //========== Send matrix matXitem to GPU - START

    if (vars_next.do_compute) {
      matXitem_buf[vars_next.index_01]->to_accel_start();
    }

    //========== Copy result matrix matB from GPU - START

    if (vars_prev.do_compute) {
      vars_prev.matB_buf_ptr()->from_accel_start();
    }

    //========== Perform pseudo GEMM matB = matX^T PROD V - START

    if (vars.do_compute) {
      gm_linalg_gemm_start(vars.I_max, nvl, npvfl,
          matXitem_buf[vars.index_01], vectors_I_buf, vectors_K_buf,
          vars.matB_buf_ptr(),
          vsums_I->sums(), vsums_J->sums(), vsums_K->sums(),
          vsums_I->counts(), vsums_J->counts(), vsums_K->counts(), vars_prev.J,
          vars.step_2way, dm, &env_);
    }

    //========== Copy result matrix matB from GPU - WAIT

    if (vars_prev.do_compute) {
      vars_prev.matB_buf_ptr()->from_accel_wait();
      if (vars_prev.step_2way == 0) {
        gm_metrics_pad_adjust(&metrics, vars_prev.matB_buf_ptr(), &env_,
          //CHECK
          env_.is_bitwise_3way_2step() && env_.metric_type() == MetricType::CCC ? 2 : 1);
      }
    }

    //========== Reduce along field procs - WAIT

    if (vars_prevprev.do_compute && env_.do_reduce()) {
      gm_reduce_metrics_wait(&(mpi_requests[vars_prevprev.index_01]),
          &vars_prevprev.matB_buf, vars_prevprev.matB_buf_ptr(), &env_);
    }

    //========== Reduce along field procs - START

    if (vars_prev.do_compute && env_.do_reduce()) {
      mpi_requests[vars_prev.index_01] = gm_reduce_metrics_start(&metrics,
          &vars_prev.matB_buf, vars_prev.matB_buf_ptr(), &env_);
    }

    //========== Compute numerators using ijk piece and (if needed) 2-way pieces

    //---NOTE: matB_buf[vars_prevprev.index_01]->d is locked now
    //---but matB_buf[vars_prevprev.index_01]->h is usable.

    if (vars_prevprev.do_compute) {
      compute_metrics_3way_block_linalg_form_metrics_(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
          &vars_prevprev.matB_buf,
          &metrics, nvl, vars_prevprev.J, vars_prevprev.step_2way,
          vars_prevprev.I_min, vars_prevprev.I_max,
          vars_prevprev.K_min, vars_prevprev.K_max,
          vars_prevprev.I_max_dim,
          j_block, k_block, si,
          vdata_i.sums, vdata_j.sums, vdata_k.sums,
          env_);
    }

  //========================================
  } // step_num
  //========================================

  // Terminations

  GMSectionInfo_destroy(si, &env_);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
