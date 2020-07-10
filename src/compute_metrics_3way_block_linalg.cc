//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way_block_gpu.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Calculate numerators, 3-way, for a single block, GPU case.
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
#include "string.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "compressed_buf.hh"
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

    if (env.form_matX_tc()) {

      for (int word = 0; word<2; ++word) {
        for (int pvfl = 0; pvfl < npvfl; ++pvfl) {

            matXitem_buf->elt<GMBits2x64>(pvfl, 0).data[word] =
              vectors_J_buf->elt_const<GMBits2x64>(pvfl, J).data[word];

        }
      }

    } else { // if (!env.form_matX_tc())

//if (env.metric_type() == MetricType::DUO)
//printf("%i\n", env.tc_eff());

      COMET_INSIST(env.metric_type() != MetricType::DUO &&
                   "Case currently unimplemented.");
      // TODO: implement DUO here.

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
            // step 0: select (only) vector entries equal to 00
            // step 1: select (only) vector entries equal to 01 or 10 (nonsparse case)
            // step 1: select (only) vector entries equal to 01 (sparse) (ignore 10)
            // step 2: select (only) vector entries equal to 11
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

    } // if (env.form_matX_tc())

  } else {

    COMET_INSIST(false);

  } // env.metric_type()

  matXitem_buf->unlock_h();
}

//-----------------------------------------------------------------------------
/// \brief Set metrics numerators and denominators using GEMM results.

template<int METRIC_FORMAT>
static void compute_metrics_3way_block_linalg_form_metrics_mf_(
  MirroredBuf* const matM_IJ_buf,
  MirroredBuf* const matM_JK_buf,
  MirroredBuf* const matM_KIK_buf,
  //MirroredBuf* const matB_buf,
  CompressedBuf* const matB_buf,
  GMMetrics* metrics,
  int nvl, int J, int step_2way,
  int I_min, int I_max, int K_min, int K_max,
  int I_max_dim,
  int j_block, int k_block,
  GMSectionInfo* const si,
  VectorSums* vs_i, VectorSums* vs_j, VectorSums* vs_k,
  CEnv& env) {

  COMET_INSIST(vs_i && vs_j && vs_k);

  COMET_INSIST( ! (env.is_bitwise_3way_2step() && !env.form_matX_tc()) &&
               "Case currently unimplemented.");

  matB_buf->lock_h();

  const VectorSums* const vs_I = si->perm0(vs_i, vs_j, vs_k);
  const VectorSums* const vs_J = si->perm1(vs_i, vs_j, vs_k);
  const VectorSums* const vs_K = si->perm2(vs_i, vs_j, vs_k);

  //--------------------
  // Compute numerators using ijk piece and (if needed) 2-way pieces.
  //--------------------
//if(env.proc_num_vector()==0) printf("Hey %i %i\n", J, step_2way);

  if (env.metric_type() == MetricType::CZEK && !env.all2all()) {

    // don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
    #pragma omp parallel for schedule(dynamic,1000)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = matM_IJ_buf->elt_const<GMFloat>(I, J);
        const GMFloat min_JK = matM_JK_buf->elt_const<GMFloat>(J, K);
        const GMFloat min_KIK = matM_KIK_buf->elt_const<GMFloat>(K, I);
        // sum of mins vectors i, j, and k is matB(k,i).
        const GMFloat min_IJK = matB_buf->elt_const<GMFloat>(I, K);
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
        Metrics_elt_3<GMFloat>(*metrics, i, j, k,
          env.proc_num_vector(), env.proc_num_vector(), env) = value;
      } // K
    }   // I
    metrics->num_metric_items_local_computed_inc((I_max - I_min) * (size_t)
                                                 (K_max - K_min));

  } else if (env.metric_type() == MetricType::CZEK && env.all2all()) {

    MetricsIndexCache index_cache = {};
    // don't use collapse because of overflow for large sizes
    //#pragma omp parallel for collapse(2) firstprivate(index_cache) schedule(dynamic,1000)
    #pragma omp parallel for firstprivate(index_cache) schedule(dynamic,1000)
    for (int K = K_min; K < K_max; ++K) {
      for (int I = I_min; I < I_max; ++I) {
        const GMFloat min_IJ = matM_IJ_buf->elt_const<GMFloat>(I, J);
        const GMFloat min_JK = matM_JK_buf->elt_const<GMFloat>(J, K);
        const GMFloat min_KIK = si->is_part3 ?
          matM_KIK_buf->elt_const<GMFloat>(K, I) :
          matM_KIK_buf->elt_const<GMFloat>(I, K);
        // sum of mins vectors i, j, and k is matB(k,i).
        const GMFloat min_IJK = matB_buf->elt_const<GMFloat>(I, K);
        const GMFloat numer = min_IJ + min_JK + min_KIK - min_IJK;
        // Make arithmetic order-independent.
        GMFloat smin, smid, smax;
        const auto sI = vs_I->sum(I);
        const auto sJ = vs_J->sum(J);
        const auto sK = vs_K->sum(K);
        utils::sort_3(smin, smid, smax, sI, sJ, sK);
        const GMFloat denom = smin + smid + smax;
        const GMFloat value = ((GMFloat)1.5) * numer / denom;
        Metrics_elt_3<GMFloat>(*metrics, I, J, K,
          j_block, k_block, index_cache, env) = value;
      } // K
    }   // I
    metrics->num_metric_items_local_computed_inc((I_max - I_min) * (size_t)
                                                 (K_max - K_min));

  } else if (env.is_metric_type_bitwise()) {

    COMET_INSIST((env.metric_type() == MetricType::CCC ||
                  env.is_bitwise_3way_2step()) &&
                 "Case of DUO + 3-step algorithm currently unimplemented.");

    enum {MF = METRIC_FORMAT};
    typedef MetricFormatTraits<MF> MFT;
    typedef typename MFT::TypeIn MFTypeIn;

    const int nvle = I_max_dim;
    const int nvleD2 = nvle / 2;
    const int is_halved = env.is_vectors_halved();
    const int num_halves = is_halved ? 2 : 1;

    matB_buf->elt_read_start();

    if (env.is_shrink()) {

      // NOTE: this may be a slight overestimate of the amount of mem needed.

#if 0
printf("%i %i %i %i\n", (int)metrics->num_metric_items_local_computed, (int)matB_buf->num_entries(), (int)metrics->num_metric_items_local_allocated, env.is_threshold_tc());
fflush(stdout);
#endif

      COMET_INSIST(metrics->num_metric_items_local_computed +
                   matB_buf->num_entries() <=
                   metrics->num_metric_items_local_allocated &&
                   "Insufficient metrics memory; "
                   "please decrease metrics_shrink.");

      const int i_block = env.proc_num_vector();

      for (size_t ind_entry = 0; ind_entry < matB_buf->num_entries();
           ++ind_entry) {

        const MFTypeIn metric_item = matB_buf->elt_const<MFTypeIn>(ind_entry);

        const size_t index = metrics->num_metric_items_local_computed;

        Metrics_elt<MFTypeIn>(*metrics, index, env) = metric_item;

        const size_t K = matB_buf->ind1_recent();
        const size_t I_mapped = matB_buf->ind0_recent();

        const int half_num = I_mapped / nvleD2;

        //const size_t I = I_mapped % nvleD2;
        const size_t I = I_mapped % nvleD2 + step_2way * nvleD2;

        const bool is_in_range = I >= (size_t)I_min && I < (size_t)I_max &&
                                 K >= (size_t)K_min && K < (size_t)K_max;

//if(I==1 && J==4 && K==6 && env.proc_num_vector()==0)
//printf("%i %i %i   %i %i %i   %i %i   %f\n", (int)I, (int)J, (int)K,
//       nvleD2, (int)I_mapped, (int)step_2way,
//       (int)ind_entry, is_in_range,
//       (double)metric_item);

        if (!is_in_range)
          continue;

        const size_t i = si->unperm0(I, (size_t)J, K);
        const size_t j = si->unperm1(I, (size_t)J, K);
        const size_t k = si->unperm2(I, (size_t)J, K);

        // TODO: accessor functions
        const size_t iG = i + nvl * i_block;
        const size_t jG = j + nvl * j_block;
        const size_t kG = k + nvl * k_block;

        const int IE = half_num;
        const int JE = matB_buf->iE_recent();
        const int KE = matB_buf->jE_recent();

//if(I==1 && J==4 && K==6 && env.proc_num_vector()==0)
//printf("%i %i\n", (int)JE, (int)KE);

        const int iE = si->unperm0(IE, JE, KE);
        const int jE = si->unperm1(IE, JE, KE);
        const int kE = si->unperm2(IE, JE, KE);

//if(iG==0 && jG==1 && kG==5)
//printf("%i %i %i  %i %i %i  %i  %i  %.20e\n", iE, jE, kE, IE, JE, KE, env.proc_num(), si->part_num, (double)metric_item);

        // TODO: accessor function
        metrics->data_coords_values_[index] =
          CoordsInfo::set(iG, jG, kG, iE, jE, kE, *metrics, env);

        metrics->num_metric_items_local_computed_inc(1);

//printf("%i  %i %i %i   %i %i %i  %f  %zu   %i %i %i %i\n", (int)ind_entry, (int)iG, (int)jG, (int)kG, iE, jE, kE, (double)metric_item, (size_t)CoordsInfo::set(iG, jG, kG, iE, jE, kE, *metrics, env), nvleD2, (int)I_mapped, (int)J, (int)step_2way);

      } // for ind_entry

    } else { // if (!env.is_shrink())

      MetricsIndexCache index_cache = {};

      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) firstprivate(index_cache) schedule(dynamic,1000)
      #pragma omp parallel for firstprivate(index_cache) schedule(dynamic,1000) if (!matB_buf->do_compress())
// if(!env.is_threshold_tc())
      for (int K = K_min; K < K_max; ++K) {

        for (int half_num = 0; half_num < num_halves; ++half_num) {

        // This "I" is the true I, the coordinate in the plane.

        for (int I = I_min; I < I_max; ++I) {

          const int i = si->unperm0(I, J, K);
          const int j = si->unperm1(I, J, K);
          const int k = si->unperm2(I, J, K);

          // For halved case:
          // first 2-way step: process lower half of I values;
          // second: upper half.

          const bool is_I_in_range = !is_halved ? true :
            I >= nvleD2 * step_2way && I < nvleD2 * (step_2way+1);

          // Numerator.

          if (is_I_in_range) {

            // Do we initialize relevant metric value on this 2-way step.

            const bool init_numer = is_halved ? 0 == half_num : 0 == step_2way;

            // NOTE: numer is accessed through a permuted index, but
            // the 000, 001 etc. indices of numer are unpermuted.

            const int j_block_eff = env.all2all() ? j_block : env.proc_num_vector();
            const int k_block_eff = env.all2all() ? k_block : env.proc_num_vector();

            auto numer = init_numer ? Tally4x2<MF>::null() :
              Metrics_elt_const_3<Tally4x2<MF>>(*metrics, I, J, K,
                j_block_eff, k_block_eff, index_cache, env);

            MFTypeIn r000, r001, r010, r011, r100, r101, r110, r111;
            MFT::decode(r000, r001, numer.data[0]);
            MFT::decode(r010, r011, numer.data[1]);
            MFT::decode(r100, r101, numer.data[2]);
            MFT::decode(r110, r111, numer.data[3]);

            MFTypeIn r000_perm = r000;
            MFTypeIn r100_perm = si->perm0(r100, r010, r001);
            MFTypeIn r010_perm = si->perm1(r100, r010, r001);
            MFTypeIn r001_perm = si->perm2(r100, r010, r001);
            MFTypeIn r011_perm = si->perm0(r011, r101, r110);
            MFTypeIn r101_perm = si->perm1(r011, r101, r110);
            MFTypeIn r110_perm = si->perm2(r011, r101, r110);
            MFTypeIn r111_perm = r111;

#if 0
    const size_t i = si->unperm0((size_t)I, (size_t)J, (size_t)K);
    const size_t j = si->unperm1((size_t)I, (size_t)J, (size_t)K);
    const size_t k = si->unperm2((size_t)I, (size_t)J, (size_t)K);

    const size_t iG = i + nvl * i_block;
    const size_t jG = j + nvl * j_block;
    const size_t kG = k + nvl * k_block;
#endif

            // Add contribution from this 2-way step.

            if (is_halved) {

              // For halved case, get the the 2 I values in matB to process.
              // For non-halved, just use I_mapped_lo (= I).

              const int I_mapped_lo = I % nvleD2;
              const int I_mapped_hi = I % nvleD2 + nvleD2;

              // Add in lo part.

              if (0 == half_num) {
                // NOTE: matB is generated by GEMM on permuted vectors objects,
                // thus has permuted 001 etc. table entries.

                const auto matB_perm =
                  matB_buf->elt_const<Tally2x2<MF>>(I_mapped_lo, K);
                MFTypeIn matB00_perm, matB01_perm, matB10_perm, matB11_perm;
                MFT::decode(matB00_perm, matB01_perm, matB_perm.data[0]);
                MFT::decode(matB10_perm, matB11_perm, matB_perm.data[1]);

                r000_perm += matB00_perm;
                r001_perm += matB01_perm;
                r010_perm += matB10_perm;
                r011_perm += matB11_perm;
//if (K==K_min) printf("lo %i %i %i %i\n", I, I_mapped_lo, J, step_2way);
//if(I==0 && J==1 && K==2 && env.proc_num_vector()==0)
//printf("lo %i %i %i   %f %f %f %f\n", (int)I, (int)J, (int)K,
//       (double)matB00_perm, (double)matB01_perm, (double)matB10_perm, (double)matB11_perm);
              }

              // Add in hi part.

              if (1 == half_num) {
                const Tally2x2<MF> matB_perm =
                  matB_buf->elt_const<Tally2x2<MF>>(I_mapped_hi, K);

                MFTypeIn matB00_perm, matB01_perm, matB10_perm, matB11_perm;
                MFT::decode(matB00_perm, matB01_perm, matB_perm.data[0]);
                MFT::decode(matB10_perm, matB11_perm, matB_perm.data[1]);

                r100_perm += matB00_perm;
                r101_perm += matB01_perm;
                r110_perm += matB10_perm;
                r111_perm += matB11_perm;
//if (K==K_min) printf("hi %i %i %i %i\n", I, I_mapped_hi, J, step_2way);
//if(I==0 && J==1 && K==2 && env.proc_num_vector()==0)
//printf("hi %i %i %i   %f %f %f %f\n", (int)I, (int)J, (int)K,
//       (double)matB00_perm, (double)matB01_perm, (double)matB10_perm, (double)matB11_perm);
              }

            } else if (env.is_bitwise_3way_2step()) { // && ! is_halved

              const auto matB_perm =
                matB_buf->elt_const<Tally2x2<MF>>(I, K);
              MFTypeIn matB00_perm, matB01_perm, matB10_perm, matB11_perm;
              MFT::decode(matB00_perm, matB01_perm, matB_perm.data[0]);
              MFT::decode(matB10_perm, matB11_perm, matB_perm.data[1]);

              if (0 == step_2way) {
                r000_perm += matB00_perm;
                r001_perm += matB01_perm;
                r010_perm += matB10_perm;
                r011_perm += matB11_perm;
              } else { // step_2way==1
                r100_perm += matB00_perm;
                r101_perm += matB01_perm;
                r110_perm += matB10_perm;
                r111_perm += matB11_perm;
              }

            } else { // if (! env.is_bitwise_3way_2step() && ! is_halved)

              const auto matB_perm =
                matB_buf->elt_const<Tally2x2<MF>>(I, K);
              MFTypeIn matB00_perm, matB01_perm, matB10_perm, matB11_perm;
              MFT::decode(matB00_perm, matB01_perm, matB_perm.data[0]);
              MFT::decode(matB10_perm, matB11_perm, matB_perm.data[1]);

              if (0 == step_2way) {
                r000_perm += 2 * matB00_perm;
                r001_perm += 2 * matB01_perm;
                r010_perm += 2 * matB10_perm;
                r011_perm += 2 * matB11_perm;
              } else if (1 == step_2way) {
                r000_perm += matB00_perm;
                r001_perm += matB01_perm;
                r010_perm += matB10_perm;
                r011_perm += matB11_perm;
                r100_perm += matB00_perm;
                r101_perm += matB01_perm;
                r110_perm += matB10_perm;
                r111_perm += matB11_perm;
              } else { // 2 == step_2way
                r100_perm += 2 * matB00_perm;
                r101_perm += 2 * matB01_perm;
                r110_perm += 2 * matB10_perm;
                r111_perm += 2 * matB11_perm;
              }

            } // if (is_halved) ...

            r000 = r000_perm;
            r100 = si->unperm0(r100_perm, r010_perm, r001_perm);
            r010 = si->unperm1(r100_perm, r010_perm, r001_perm);
            r001 = si->unperm2(r100_perm, r010_perm, r001_perm);
            r011 = si->unperm0(r011_perm, r101_perm, r110_perm);
            r101 = si->unperm1(r011_perm, r101_perm, r110_perm);
            r110 = si->unperm2(r011_perm, r101_perm, r110_perm);
            r111 = r111_perm;

            MFT::encode(numer.data[0], r000, r001);
            MFT::encode(numer.data[1], r010, r011);
            MFT::encode(numer.data[2], r100, r101);
            MFT::encode(numer.data[3], r110, r111);

            Metrics_elt_3<Tally4x2<MF>>(*metrics, I, J, K,
              j_block_eff, k_block_eff, index_cache, env) = numer;
//printf("%i %i %i  %f %f %f %f %f %f %f %f\n", I, J, K, (double)r000, (double)r001, (double)r010, (double)r011, (double)r100, (double)r101, (double)r110, (double)r111);

#if 0
int index_ = Metrics_index_3(*metrics, (int)I, J, (int)K, j_block, k_block, index_cache, env);
if (
(index_ == 204 && env.proc_num()==0 && env.compute_method()==2)
||
(index_ == 1120 && env.proc_num()==1 && env.compute_method()==0)
)
printf("0  %i   %i %i %i %i %i %i %i %i \n", si->part_num,
(int)r000, (int)r001, (int)r010, (int)r011,
(int)r100, (int)r101, (int)r110, (int)r111);
fflush(stdout);
#endif

          } // if (is_I_in_range)

          // Denominator.

          if ((!env.is_threshold_tc()) &&
              step_2way == env.num_step_2way_for_3way() - 1) {

            const auto si1 = (GMTally1)vs_i->sum(i);
            const auto sj1 = (GMTally1)vs_j->sum(j);
            const auto sk1 = (GMTally1)vs_k->sum(k);
            const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si1, sj1, sk1);

            const int j_block_eff = env.all2all() ? j_block : env.proc_num_vector();
            const int k_block_eff = env.all2all() ? k_block : env.proc_num_vector();
            
            Metrics_elt_3<GMFloat3, MetricsArray::S>(*metrics, I, J, K,
              j_block_eff, k_block_eff, index_cache, env) = si1_sj1_sk1;
            if (env.sparse()) {
              const auto ci1 = (GMTally1)vs_i->count(i);
              const auto cj1 = (GMTally1)vs_j->count(j);
              const auto ck1 = (GMTally1)vs_k->count(k);
              const GMFloat3 ci1_cj1_ck1 = GMFloat3_encode(ci1, cj1, ck1);
              Metrics_elt_3<GMFloat3, MetricsArray::C>(*metrics, I, J, K,
                j_block_eff, k_block_eff, index_cache, env) = ci1_cj1_ck1;
            } // if sparse

          } // if ((!env.is_threshold_tc()) && ...

        } // I
        } // for half_num
      } // K
      if (step_2way == env.num_step_2way_for_3way() - 1) {
        metrics->num_metric_items_local_computed_inc((I_max - I_min) * (size_t)
                                                     (K_max - K_min));
      }

    } // if (env.is_shrink())

  } else { // if (env.metric_type() ...

    COMET_INSIST(false);

  } // if (env.metric_type() ...

  matB_buf->unlock_h();
}

//-----------------------------------------------------------------------------
/// \brief Set metrics numerators and denominators using GEMM results.

static void compute_metrics_3way_block_linalg_form_metrics_(
  MirroredBuf* const matM_IJ_buf,
  MirroredBuf* const matM_JK_buf,
  MirroredBuf* const matM_KIK_buf,
  //MirroredBuf* const matB_buf,
  CompressedBuf* const matB_buf,
  GMMetrics* metrics,
  int nvl, int J, int step_2way,
  int I_min, int I_max, int K_min, int K_max,
  int I_max_dim,
  int j_block, int k_block,
  GMSectionInfo* const si,
  VectorSums* vs_i, VectorSums* vs_j, VectorSums* vs_k,
  CEnv& env) {

  if (env.is_threshold_tc()) {

    compute_metrics_3way_block_linalg_form_metrics_mf_<
      MetricFormat::SINGLE>(
      matM_IJ_buf, matM_JK_buf, matM_KIK_buf, matB_buf, metrics,
      nvl, J, step_2way, I_min, I_max, K_min, K_max, I_max_dim,
      j_block, k_block, si, vs_i, vs_j, vs_k, env);

  } else {

    compute_metrics_3way_block_linalg_form_metrics_mf_<
      MetricFormat::PACKED_DOUBLE>(
      matM_IJ_buf, matM_JK_buf, matM_KIK_buf, matB_buf, metrics,
      nvl, J, step_2way, I_min, I_max, K_min, K_max, I_max_dim,
      j_block, k_block, si, vs_i, vs_j, vs_k, env);

  }
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
  COMET_INSIST(env_.num_way() == NumWay::_3);

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

  MirroredBuf* const matM_ij_buf = need_mat_ij ? &matM_ij_buf_ : NULL;

  if (need_mat_ij) {
    MirroredBuf* matM_ij_buf_ptr = env_.do_reduce()
      ? tmp_buf_[0] : matM_ij_buf;

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
      ? tmp_buf_[0] : matM_jk_buf;

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
        env_.do_reduce() ? tmp_buf_[0] : matM_kik_buf;

    gm_linalg_gemm(nvl, nvl, npvfl,
                   vdata_k.buf, vdata_i.buf, matM_kik_buf_ptr,
                   vdata_k.sums->sums(), vdata_i.sums->sums(),
                   vdata_k.sums->counts(), vdata_i.sums->counts(),
                   dm, &env_);

    matM_kik_buf_ptr->from_accel();

    gm_reduce_metrics(&metrics, matM_kik_buf, matM_kik_buf_ptr, &env_);
  } // is_part3

  //----------------------------------------
  // Now compute ijk piece, via an outer loop over j values.
  //----------------------------------------

  // Allocate magma CPU/GPU memory for matrices X and B.
  // X = elementwise OP of one vector with each of the rest of the vectors.
  // For the jth iteration, the ith column of X is the elementwise OP
  //   of vectors i and j.
  // B = X^T PROD V = three way PROD.


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

  if (env_.form_matX_tc()) {
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

//  std::vector<LoopVars> vars_buf;
//  for (int i=0; i<num_buf; ++i)
//    vars_buf.push_back(LoopVars(env_));

  LoopVars vars_buf0(env_);
  LoopVars vars_buf1(env_);
  LoopVars vars_buf2(env_);
  LoopVars vars_buf3(env_);

  LoopVars* const vars_buf[num_buf] = {&vars_buf0, &vars_buf1, &vars_buf2, &vars_buf3};

  CompressedBuf matB_buf_compressed(*matB_buf_[0], env_);

  const int first_step = 0 - extra_step;

  //========================================
  for (int step_num = first_step; step_num < num_step+extra_step*2;
       ++step_num) {
  //========================================

    // Set per-step variables.

    LoopVars& vars_prevprev = *vars_buf[(step_num - first_step + 0) % num_buf];
    LoopVars& vars_prev = *vars_buf[(step_num - first_step + 1) % num_buf];
    LoopVars& vars = *vars_buf[(step_num - first_step + 2) % num_buf];
    LoopVars& vars_next = *vars_buf[(step_num - first_step + 3) % num_buf];

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
        vars_next.tmp_buf.allocate(*tmp_buf_[vars_next.index_01],
                                   vars_next.I_max_dim);
      }
      vars_next.matB_buf.allocate(*matB_buf_[vars_next.index_01],
                                  vars_next.I_max_dim);
    }

    //TODO: fix locks to work properly for CPU case.

    //========== Send matrix matXitem to GPU - WAIT

    if (vars.do_compute) {
      matXitem_buf_[vars.index_01]->to_accel_wait();
    }

    //========== Perform pseudo GEMM matB = matX^T PROD V - WAIT

    if (vars_prev.do_compute) {
      gm_linalg_gemm_wait(vars_prev.I_max, nvl, npvfl,
          matXitem_buf_[vars_prev.index_01], vectors_I_buf, vectors_K_buf,
//c
          vars_prev.matB_buf_ptr(),
          vsums_I->sums(), vsums_J->sums(), vsums_K->sums(),
          vsums_I->counts(), vsums_J->counts(), vsums_K->counts(),
          vars_prev.J, vars_prev.step_2way, dm, &env_);
      matB_buf_compressed.attach(*vars_prev.matB_buf_ptr());
      matB_buf_compressed.compress();
    }

    //========== Calculate matXitem

    if (vars_next.do_compute) {
      compute_metrics_3way_block_linalg_form_matXitem_(
          vectors_I_buf, vectors_J_buf, matXitem_buf_[vars_next.index_01],
          vars_next.J, vars_next.step_2way,
          vars_next.I_min, vars_next.I_max, npvfl, env_);
    }

    //========== Send matrix matXitem to GPU - START

    if (vars_next.do_compute) {
      matXitem_buf_[vars_next.index_01]->to_accel_start();
    }

    //========== Copy result matrix matB from GPU - START

    if (vars_prev.do_compute) {
//c
      //vars_prev.matB_buf_ptr()->from_accel_start();
      matB_buf_compressed.from_accel_start();
    }

    //========== Perform pseudo GEMM matB = matX^T PROD V - START

    if (vars.do_compute) {
      gm_linalg_gemm_start(vars.I_max, nvl, npvfl,
          matXitem_buf_[vars.index_01], vectors_I_buf, vectors_K_buf,
          vars.matB_buf_ptr(),
          vsums_I->sums(), vsums_J->sums(), vsums_K->sums(),
          vsums_I->counts(), vsums_J->counts(), vsums_K->counts(),
          vars.J, vars.step_2way, dm, &env_);
    }

    //========== Copy result matrix matB from GPU - WAIT

    if (vars_prev.do_compute) {
//c
      //vars_prev.matB_buf_ptr()->from_accel_wait();
      matB_buf_compressed.from_accel_wait();
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
      matB_buf_compressed.attach(vars_prevprev.matB_buf);
    }

    //========== Reduce along field procs - START

    if (vars_prev.do_compute && env_.do_reduce()) {
      mpi_requests[vars_prev.index_01] = gm_reduce_metrics_start(&metrics,
          &vars_prev.matB_buf, vars_prev.matB_buf_ptr(), &env_);
    }

    //========== Compute numerators using ijk piece and (if needed) 2-way pieces

    //---NOTE: matB_buf[vars_prevprev.index_01]->d is locked now
    //---but matB_buf[vars_prevprev.index_01]->h is usable.

    LoopVars& vars_tail = env_.do_reduce() ? vars_prevprev : vars_prev;

    if (vars_tail.do_compute) {
      compute_metrics_3way_block_linalg_form_metrics_(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
//c
          &matB_buf_compressed,
          //&vars_tail.matB_buf,
          &metrics, nvl, vars_tail.J, vars_tail.step_2way,
          vars_tail.I_min, vars_tail.I_max,
          vars_tail.K_min, vars_tail.K_max,
          vars_tail.I_max_dim,
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
