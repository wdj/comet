//-----------------------------------------------------------------------------
/*!
 * \file   tc_out.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, tc package: copying data from the accelerator.
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

#ifndef _COMET_TC_OUT_I_HH_
#define _COMET_TC_OUT_I_HH_

#include "formulas.hh"
#include "histograms.hh"

//#include "env.hh"
//#include "tc.hh"
//#include "tc_int.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Swizzle individual elements in buf.

template<typename GemmOut_t, int METRIC_FORMAT>
__host__ __device__ static void tc_repair_metrics_kernel_elt_(
  int nvl, int nvll, int nvllD2, void* vo, int thread_r, int thread_c) { 

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type MFType;

  // Considered as an array of floats, array is 2*nvl rows X 2*nvl cols.
  // Each thread manipulates a block of 4 rows and 2 cols.
  // Thus the dimensions of the metrics array in blocks is nvllD2 X nvl.
  // Each block viewed as an array of doubles is 2 X 2.

  // Two col numbers being processed of this (float) array.

  // ISSUE: does the compiler need to / understand that the pointers are aliased

  const size_t ivo_offset0 = 4*thread_r + thread_c * (size_t)(4*nvll);
  const size_t ivo_offset1 = 4*thread_r + thread_c * (size_t)(4*nvll) + 2*nvll;

  // Read the 8 values.

  GemmOut_t* const ivo = (GemmOut_t*)vo;

  const GemmOut_t i00 = ivo[ivo_offset0+0];
  const GemmOut_t i01 = ivo[ivo_offset0+1];
  const GemmOut_t i02 = ivo[ivo_offset0+2];
  const GemmOut_t i03 = ivo[ivo_offset0+3];

  const GemmOut_t i10 = ivo[ivo_offset1+0];
  const GemmOut_t i11 = ivo[ivo_offset1+1];
  const GemmOut_t i12 = ivo[ivo_offset1+2];
  const GemmOut_t i13 = ivo[ivo_offset1+3];

//printf("%i %i %i %i %i %i %i %i\n",
//(int)i00, (int)i01, (int)i02, (int)i03, (int)i10, (int)i11, (int)i12, (int)i13);

  // Apply the permutation:

  // [ i00  i10 ]  ->  [ i00  i02 ]
  // [ i01  i11 ]  ->  [ i01  i03 ]
  // [ i02  i12 ]  ->  [ i10  i12 ]
  // [ i03  i13 ]  ->  [ i11  i13 ]

  const GemmOut_t i00p = i00;
  const GemmOut_t i01p = i01;

  const GemmOut_t i02p = i10;
  const GemmOut_t i03p = i11;

  const GemmOut_t i10p = i02;
  const GemmOut_t i11p = i03;

  const GemmOut_t i12p = i12;
  const GemmOut_t i13p = i13;

  // Pack two 26-bit integers into mantissa of double.

  MFType o00 = {}, o01 = {}, o10 = {}, o11 = {};
  MFT::encode(o00, i00p, i02p);
  MFT::encode(o01, i01p, i03p);
  MFT::encode(o10, i10p, i12p);
  MFT::encode(o11, i11p, i13p);

  // Overwrite block with the new values.
  // All is isolated to a single thread, should be thread safe.

  const size_t o_offset0 = 2 * thread_r + thread_c * (size_t)(2*nvll);
  const size_t o_offset1 = 2 * thread_r + thread_c * (size_t)(2*nvll) + nvll;

  MFType* const ovo = (MFType*)vo;

  ovo[o_offset0+0] = o00;
  ovo[o_offset0+1] = o01;

  ovo[o_offset1+0] = o10;
  ovo[o_offset1+1] = o11;
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support tc_repair_metrics_.
///
///        This function has two purposes:
///        1. Convert the 2X2 table from each pair of compared vectors
///        from 4 32-bit (int32 or float32) values to the required
///        16-byte double complex packed format.
///        2. Permute the table elements to the required places.
///
///        The reason for the permutation is as follows.
///        For the output matrix of this function, each single 2X2 matrix
///        is arranged contiguously in memory as a double complex value.
///        However, the input matrices to the GEMM do not give a result
///        matrix that is consistent with this ordering.
///        Thus there needs to be a copy to rearrange.  Furthermore,
///        we want to make this an in-place rearrangement to save
///        space, and additionally we want to assign work to threads
///        with no race conditions and with coalesced memory accesses.
///
///        The method can be explained as follows.
///        1. The input "left" and "right" matrices to the modified GEMM
///        can be thought of each as a group of column vectors.
///        2. Each column (of 2-bit vector entries) is converted into two
///        columns, with elements being the counts of 0 bits and 1 bits of the
///        original vectors.  Each pair of vectors is kept together
///        side-by-side in these new left and right matrices L and R.
///        3. The columns of L are permuted, to give L' = L P
///        Example:
///          R  = [ G, G, H, H, I, I, J, J, K, K, L, L ]
///          L  = [ A, A, B, B, C, C, D, D, E, E, F, F ]
///          L' = [ A, A, D, D, B, B, E, E, C, C, F, F ]
///        (note L is used in 2 different senses here)
///        4. The GEMM is computed, M = (L')^T R = P^T L^T R.  Because of
///        the permutation of L, the rows of M are permuted.
///        Here, for brevity we drop the transpose, writing A^T G as AG, etc.
///          M = [ AG, AG, AH, AH, . . . ]
///              [ AG, AG, AH, AH, . . . ]
///              [ DG, DG, DH, DH, . . . ]
///              [ DG, DG, DH, DH, . . . ]
///              [ BG, BG, BH, BH, . . . ]
///              [ BG, BG, BH, BH, . . . ]
///              [ EG, EG, EH, EH, . . . ]
///              [ EG, EG, EH, EH, . . . ]
///              [ CG, CG, CH, CH, . . . ]
///              [ CG, CG, CH, CH, . . . ]
///              [ FG, FG, FH, FH, . . . ]
///              [ FG, FG, FH, FH, . . . ]
///        Here we are considering M to be stored in column-major order.
///        5. Next we consider this as composed of size 4X2 blocks,
///        assign a CUDA thread to each block and do an in-block
///        permutation. Note each thread loads 2 16-byte (double) words,
///        with stride between threads of 16 bytes.
///        (need to check on efficiency of this w.r.t. coalescing etc.)
///          [ AG, AG ] -> [ AG, DG ]
///          [ AG, AG ] -> [ AG, DG ]
///          [ DG, DG ] -> [ AG, DG ]
///          [ DG, DG ] -> [ AG, DG ]
///        As can be seen, all four entries AG of the table are now
///        contiguous in memory.

template<typename GemmOut_t, int METRIC_FORMAT>
__global__ static void tc_repair_metrics_kernel_(
  int nvl, int nvll, int nvllD2, void* vo) { 

  // Row and column threads of metrics array.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvllD2 || thread_c >= nvl)
    return;

  tc_repair_metrics_kernel_elt_<GemmOut_t, METRIC_FORMAT>(
    nvl, nvll, nvllD2, vo,
    thread_r, thread_c);
}

//-----------------------------------------------------------------------------
/// \brief Swizzle/cast values from cublas call into double complex format.
///
///        The cublas gemm poduces a matrix of scalars of 32 bit size
///        (int32 or float).  However the required format of the metrics
///        is a matrix of double complex values, with each double
///        containing two packed 26-bit integers.
///        This code does an in-place transformation from one to the other.

template<int TC_METHOD, int METRIC_FORMAT>
void tc_repair_metrics_( int nvll, int nvl, void* vo, CEnv& env) {
  COMET_INSIST(vo);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);

  // always true, because of tc_gemm_vaxis_divisibility_required()
  COMET_INSIST(nvll % 2 == 0 && "Failed divisibility condition for tc gemm.");
  const int nvllD2 = nvll / 2;

  typedef typename TCTraits<TC_METHOD>::GemmOut_t GemmOut_t;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

      const int threadblocksize = 256;
      COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                   "Current HIP limitation.");
      const int vll2_threadblocks = utils::ceil(nvllD2, threadblocksize);

      COMET_LAUNCH_KERNEL((tc_repair_metrics_kernel_<GemmOut_t, METRIC_FORMAT>),
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        nvl, nvll, nvllD2, vo);

      COMET_INSIST(System::accel_last_call_succeeded());

  } else { // (!env.is_compute_method_gpu())

    for (int thread_c=0; thread_c<nvl; ++thread_c) {
      for (int thread_r=0; thread_r<nvllD2; ++thread_r) {
        tc_repair_metrics_kernel_elt_<GemmOut_t, METRIC_FORMAT>(
          nvl, nvll, nvllD2, vo, thread_r, thread_c);
      }
    }

  } // if compute_method
}

//=============================================================================

class GemmShapeNone {
public:

  __host__ __device__ bool is_inside(int I, int J, int K = 0) const {
    return true;
  }
};

//=============================================================================

class HistogramsHelperNone {
public:

  typedef Histograms::Elt_t Elt_t;

  __host__ __device__ void add(double value, int indT_I, int indT_J, int indT_K = 0) {}

  __host__ __device__ void finalize(Elt_t* histograms_ptr, int num_buckets) {}

};

//=============================================================================

class HistogramsHelper2Way {

  double entries[4];

public:

  typedef Histograms::Elt_t Elt_t;

  __host__ __device__ void add(double value, int indT_I, int indT_J) {
    entries[indT_J + 2 * indT_I] = value;
  }

  __host__ __device__ void finalize(Elt_t* histograms_ptr, int num_buckets) {
    enum {LL=0, LH=1, HL=2, HH=3};

    Histograms::add(histograms_ptr, num_buckets, entries[LL], HistogramID::LL);
    Histograms::add(histograms_ptr, num_buckets, entries[LH], HistogramID::LH);
    Histograms::add(histograms_ptr, num_buckets, entries[HL], HistogramID::LH);
    Histograms::add(histograms_ptr, num_buckets, entries[HH], HistogramID::HH);
    Histograms::add(histograms_ptr, num_buckets, entries[LL] + entries[HH],
                                                            HistogramID::LLHH);
  }

};

//=============================================================================

class HistogramsHelper3Way {
public:

  typedef Histograms::Elt_t Elt_t;

  __host__ __device__ void add(double value, int indT_I, int indT_J, int indT_K) {}

  __host__ __device__ void finalize(Elt_t* histograms_ptr, int num_buckets) {}

};

//=============================================================================
/// \brief Do thresholding of metrics: individual elements, 2-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT, class HistogramsHelper,
         class GemmShape = GemmShapeNone>
__host__ __device__ void tc_threshold_2way_kernel_elt_(
  int thread_r, int thread_c,
  int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* counts_I, GMFloat* counts_J,
  double param, double multiplier, double threshold_eff, bool is_using_xor,
  Histograms::Elt_t* histograms_ptr, int num_buckets, GemmShape gemm_shape = {}) {
  COMET_ASSERT(vo && sums_I && sums_J && counts_I && counts_J);
  COMET_ASSERT(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_ASSERT(thread_r >= 0 && thread_r < nvll);
  COMET_ASSERT(thread_c >= 0 && thread_c < nvl);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  //----------
  // Initializations.
  //----------

  enum {CBPE = COUNTED_BITS_PER_ELT};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type MFType;
  typedef typename MFT::TypeIn MFTypeIn;

  // Metrics array, represented as 64-bit elements, each with 2 entries.

  MFType* const dvo = (MFType*)vo;

  //----------
  // Compute indexing.
  //----------

  // The incoming matrix is of dimension nvll X nvl.
  // Each matrix element is 4 table entries, stored as 16 bytes.
  // The matrix is stored column major.

  // NOMENCLATURE:
  // I, J - coordinate of element in the 2D plane.
  // indT_I, indT_J - indices into the 2x2 table at each coord

  // TODO: (maybe) use MirroredBuf class for indexing.

  const int I = thread_r;
  COMET_ASSERT(I >= 0 && I < nvll);

  const int J = thread_c;
  COMET_ASSERT(J >= 0 && J < nvl);

  //----------
  // Get some values needed for the CCC/DUO formula.
  //----------

  const GMTally1 cI = (GMTally1)counts_I[I];
  const GMTally1 cJ = (GMTally1)counts_J[J];

  const double recip_cI = 1e0 / cI;
  const double recip_cJ = 1e0 / cJ;

  const GMTally1 sI1 = (GMTally1)sums_I[I];
  const GMTally1 sJ1 = (GMTally1)sums_J[J];
  const GMTally1 sI0 = CBPE * cI - sI1;
  const GMTally1 sJ0 = CBPE * cJ - sJ1;

  //----------
  // Calculate cij, sum of all 4 metric entries at this I,J.
  //----------

  const int mftype_per_table = Tally2x2<METRIC_FORMAT>::NUM; // = 2
  enum {indT_J_0 = 0, indT_J_1 = 1};

  // Loop over entries in table.

  GMTally1 cij = 0;
  for (int indT_I = 0; indT_I < 2; ++indT_I) {
    // NOTE: indT_J runs across the two packed values (unrolled below).

    const size_t index = indT_I + mftype_per_table * (
                         thread_r + nvll * (
                         thread_c));

    // Access entry of metrics array.
    const MFType dvo_this = dvo[index];

    const GMTally1 sI = indT_I == 0 ? sI0 : sI1;

    // Access the two packed table entries from metrics array.
    MFTypeIn values[2];
    MFT::decode(values[indT_J_0], values[indT_J_1], dvo_this);

    // Accumulate result.  Use special formula if xor-gemm.
    cij += is_using_xor ?
      ((sI + sJ0 - (GMTally1)values[indT_J_0])/2) +
      ((sI + sJ1 - (GMTally1)values[indT_J_1])/2) :
                   (GMTally1)values[indT_J_0] +
                   (GMTally1)values[indT_J_1];
  }

  const double recip_sumcij = 1e0 / cij;

  HistogramsHelper helper;

  //----------
  // Loop to compute and threshold for all table entries.
  //----------

  for (int indT_I = 0; indT_I < 2; ++indT_I) {

    const GMTally1 sI = indT_I == 0 ? sI0 : sI1;

    const size_t index = indT_I + mftype_per_table * (
                         thread_r + nvll * (
                         thread_c));

    // Access entry of metrics array.
    MFType& dvo_this = dvo[index];

    // Access the two packed table entries from metrics array.
    MFTypeIn values[2];
    MFT::decode(values[indT_J_0], values[indT_J_1], dvo_this);

    for (int indT_J = 0; indT_J < 2; ++indT_J) {

      const GMTally1 sJ = indT_J == 0 ? sJ0 : sJ1;

      // Get numerator value.  Use special formula if xor-gemm.
      const auto rij = is_using_xor ?
        ((sI + sJ - (GMTally1)values[indT_J])/2) :
                    (GMTally1)values[indT_J];

      // Apply CCC/DUO formula.
      const bool is_zero_denom = 0 == cI || 0 == cJ || 0 == cij;
      const double metric_entry = is_zero_denom ? 0e0 :
        ccc_duo_value<CBPE, double>(
          rij, sI, sJ, recip_cI, recip_cJ,
          recip_sumcij, multiplier, param);

      // Save entry to store into histogram.
      helper.add(metric_entry, indT_I, indT_J);

      // Check against threshold.
      const bool pass_threshold = CEnv::pass_threshold(
        (double)(MFTypeIn)metric_entry, threshold_eff);

      // Apply threshold.
      values[indT_J] = pass_threshold ? (MFTypeIn)metric_entry :
                                        (MFTypeIn)0;
    } // indT_J

    // Store result back into metrics array.
    MFT::encode(dvo_this, values[indT_J_0], values[indT_J_1]);
  } // indT_I

  if (gemm_shape.is_inside(I, J))
    helper.finalize(histograms_ptr, num_buckets);
}

//=============================================================================
/// \brief Do thresholding of metrics: individual elements, 3-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT, class HistogramsHelper,
         class GemmShape = GemmShapeNone>
__host__ __device__ void tc_threshold_3way_kernel_elt_(
  int thread_r, int thread_c,
  int nvll, int nvllX2, int nvllD2,  int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K,
  uint32_t* matX_counts, int J, int step_2way, double param, double multiplier,
  double threshold_eff, bool is_using_xor,
  Histograms::Elt_t* histograms_ptr, int num_buckets, GemmShape gemm_shape = {}) {
  COMET_ASSERT(vo);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_ASSERT(nvll*2 == nvllX2 && nvll/2 == nvllD2);
  COMET_ASSERT(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_ASSERT(thread_r >= 0 && thread_r < nvllD2);
  COMET_ASSERT(thread_c >= 0 && thread_c < nvl);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  //----------
  // Initializations.
  //----------

  enum {CBPE = COUNTED_BITS_PER_ELT};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type MFType;
  typedef typename MFT::TypeIn MFTypeIn;

  // Metrics array, represented as 64-bit elements, each with 2 entries.

  MFType* const dvo = (MFType*)vo;

  //----------
  // Compute indexing.
  //----------

  // The incoming matrix is of dimension nvll X nvl.
  // Each matrix element is 4 table entries, stored as 16 bytes.
  // The matrix is stored column major.
  // Recall the matrix is "halved": the lower and upper half pertain
  // to different table entries (indT_I axis) for the same (I,J,K) coordinate

  // NOMENCLATURE:
  // I, J, K - (permuted) coordinate of element in the 3D block
  // indT_I, indT_J, indT_K - indices into the 2x2x2 table at each coord

  // TODO: (maybe) use MirroredBuf class for indexing.

  const int I = thread_r + nvllD2 * step_2way;
  COMET_ASSERT(I >= 0 && I < nvll);
  COMET_ASSERT(I/nvllD2 == step_2way);

  // NOTE: "J" is already defined as an argument to this function - denotes
  // the plane being processed.

  const int K = thread_c;
  COMET_ASSERT(K >= 0 && K < nvl);

  //----------
  // Get some values needed for the CCC/DUO formula.
  //----------

  const GMTally1 cI = (GMTally1)counts_I[I];
  const GMTally1 cJ = (GMTally1)counts_J[J];
  const GMTally1 cK = (GMTally1)counts_K[K];

  const double recip_cI = 1e0 / cI;
  const double recip_cJ = 1e0 / cJ;
  const double recip_cK = 1e0 / cK;

  const GMTally1 sI1 = (GMTally1)sums_I[I];
  const GMTally1 sJ1 = (GMTally1)sums_J[J];
  const GMTally1 sK1 = (GMTally1)sums_K[K];
  const GMTally1 sI0 = CBPE * cI - sI1;
  const GMTally1 sJ0 = CBPE * cJ - sJ1;
  const GMTally1 sK0 = CBPE * cK - sK1;

  const uint32_t ignore_me = 0;

  //----------
  // Calculate cijk, sum of all 8 metric entries at this I,J,K.
  //----------

  // NOTE for 3-way we process two 2X2 tables.
  const int mftype_per_table = Tally2x2<METRIC_FORMAT>::NUM; // = 2
  const size_t halves_per_whole = 2;
  enum {indT_K_0 = 0, indT_K_1 = 1};

  // Loop over entries in table.

  GMTally1 cijk = 0;
  for (int indT_I = 0; indT_I < 2; ++indT_I) {
    for (int indT_J = 0; indT_J < 2; ++indT_J) {
      // NOTE: indT_K runs across the two packed values (unrolled below).

      const size_t index = indT_J + mftype_per_table * (
                           thread_r + nvllD2 * (
                           indT_I + halves_per_whole * (
                           thread_c)));

      // Access entry of metrics array.
      const MFType dvo_this = dvo[index];

      // Need special vector element count to make xor-gemm work right.
      // NOTE: indexing here matches the row part of metrics index just above.
      const GMTally1 sIJ = is_using_xor ?
        matX_counts[indT_J + mftype_per_table * (
                    thread_r + nvllD2 *
                    indT_I)] : ignore_me;

      // Access the two packed table entries from metrics array.
      MFTypeIn values[2];
      MFT::decode(values[indT_K_0], values[indT_K_1], dvo_this);

      // Accumulate result.  Use special formula if xor-gemm.
      cijk += is_using_xor ?
        (sIJ + sK0 - (GMTally1)values[indT_K_0])/2 +
        (sIJ + sK1 - (GMTally1)values[indT_K_1])/2 :
                     (GMTally1)values[indT_K_0] +
                     (GMTally1)values[indT_K_1];
    }
  }

  const double recip_sumcijk = 1e0 / cijk;

  HistogramsHelper helper;

  //----------
  // Loop to compute and threshold for all table entries.
  //----------

  for (int indT_I = 0; indT_I < 2; ++indT_I) {

    const GMTally1 sI = indT_I == 0 ? sI0 : sI1;

    for (int indT_J = 0; indT_J < 2; ++indT_J) {

      const GMTally1 sJ = indT_J == 0 ? sJ0 : sJ1;

      const size_t index = indT_J + mftype_per_table * (
                           thread_r + nvllD2 * (
                           indT_I + halves_per_whole * (
                           thread_c)));

      // Access entry of metrics array.
      MFType& dvo_this = dvo[index];

      // Need special vector element count to make xor-gemm work right.
      const GMTally1 sIJ = is_using_xor ?
        matX_counts[indT_J + mftype_per_table * (
                    thread_r + nvllD2 *
                    indT_I)] : ignore_me;

      // Access the two packed table entries from metrics array.
      MFTypeIn values[2];
      MFT::decode(values[indT_K_0], values[indT_K_1], dvo_this);

      for (int indT_K = 0; indT_K < 2; ++indT_K) {

        const GMTally1 sK = indT_K == 0 ? sK0 : sK1;

        // Get numerator value.  Use special formula if xor-gemm.
        const auto rijk = is_using_xor ?
          ((sIJ + sK - (GMTally1)values[indT_K])/2) :
                       (GMTally1)values[indT_K];

        // Apply CCC/DUO formula.
        const bool is_zero_denom = 0 == cI || 0 == cJ || 0 == cK || 0 == cijk;
        const double metric_entry = is_zero_denom ? 0e0 :
          ccc_duo_value<CBPE, double>(
            rijk, sI, sJ, sK, recip_cI, recip_cJ, recip_cK,
            recip_sumcijk, multiplier, param);


        // Save entry to store into histogram.
        helper.add(metric_entry, indT_I, indT_J, indT_K);

        // Check against threshold.
        const bool pass_threshold = CEnv::pass_threshold(
          (double)(MFTypeIn)metric_entry, threshold_eff);

        // Apply threshold.
        values[indT_K] = pass_threshold ? (MFTypeIn)metric_entry:
                                          (MFTypeIn)0;
      } // indT_K

      // Store result back into metrics array.
      MFT::encode(dvo_this, values[indT_K_0], values[indT_K_1]);
    } // indT_J
  } // indT_I

  if (gemm_shape.is_inside(I, J, K))
    helper.finalize(histograms_ptr, num_buckets);
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics: GPU kernel, 2-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT, class HistogramsHelper,
         class GemmShape = GemmShapeNone>
__global__ void tc_threshold_2way_kernel_(
  int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* counts_I, GMFloat* counts_J,
  double param, double multiplier, double threshold_eff, bool is_using_xor,
  Histograms::Elt_t* histograms_ptr = 0, int num_buckets = 0,
  GemmShape gemm_shape = {}) {
  COMET_ASSERT(vo && sums_I && sums_J && counts_I && counts_J);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  // Row and column threads for metrics matrix.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvll || thread_c >= nvl)
    return;

  tc_threshold_2way_kernel_elt_<COUNTED_BITS_PER_ELT, METRIC_FORMAT,
                                HistogramsHelper>(
    thread_r, thread_c,
    nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
    param, multiplier, threshold_eff, is_using_xor,
    histograms_ptr, num_buckets, gemm_shape);
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics: GPU kernel, 3-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT, class HistogramsHelper,
         class GemmShape = GemmShapeNone>
__global__ void tc_threshold_3way_kernel_(
  int nvll, int nvllX2, int nvllD2,  int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K,
  uint32_t* matX_counts, int J, int step_2way, double param, double multiplier,
  double threshold_eff, bool is_using_xor,
  Histograms::Elt_t* histograms_ptr = 0, int num_buckets = 0,
  GemmShape gemm_shape = {}) {
  COMET_ASSERT(vo);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  // Row and column threads for metrics matrix.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvllD2 || thread_c >= nvl)
    return;

  tc_threshold_3way_kernel_elt_<COUNTED_BITS_PER_ELT, METRIC_FORMAT,
                                HistogramsHelper>(
    thread_r, thread_c,
    nvll, nvllX2, nvllD2, nvl, vo,
    sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
    J, step_2way, param, multiplier, threshold_eff, is_using_xor,
    histograms_ptr, num_buckets, gemm_shape);
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics if requested: templated on CBPE.

template<int CBPE, int MF>
void tc_threshold_per_CBPE_(int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K,
  uint32_t* matX_counts, int J, int step_2way,
  Histograms& histograms, CEnv& env) {

  const int nvllX2 = nvll * 2;
  const int nvllD2 = nvll / 2;

  // ? template specialization on num_way
  // 2way: i -> I, j -> K (???)

  const double param = env.ccc_param();
  const double multiplier = env.ccc_duo_multiplier();
  const double threshold_eff = env.threshold_eff();
  const bool is_using_xor = env.is_using_xor();
  COMET_INSIST_INTERFACE(&env, param <= 1 &&
                         "CCC/DUO param value not allowed.");
  COMET_INSIST_INTERFACE(&env, multiplier >= 0 &&
                         "CCC/DUO multiplier value not allowed.");

  const int num_threads_r = env.num_way() == NumWay::_2 ? nvll : nvllD2;
  const int num_threads_c = nvl;

  const bool is_computing_histograms = histograms.is_computing_histograms();

  if (env.is_compute_method_gpu()) {

    // Kernel call.

    const int threadblocksize = 256;
    COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                 "Current HIP limitation.");
    const int vll_threadblocks = utils::ceil(num_threads_r, threadblocksize);

    if (NumWay::_2 == env.num_way()) {

      if (is_computing_histograms)

        COMET_LAUNCH_KERNEL((tc_threshold_2way_kernel_<CBPE, MF,
                                                       HistogramsHelper2Way>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
          param, multiplier, threshold_eff, is_using_xor,
          histograms.get_ptr(), histograms.num_buckets());

      else

        COMET_LAUNCH_KERNEL((tc_threshold_2way_kernel_<CBPE, MF,
                                                       HistogramsHelperNone>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
          param, multiplier, threshold_eff, is_using_xor);


    } else { // if (NumWay::_3 == env.num_way())

      if (is_computing_histograms)

        COMET_LAUNCH_KERNEL((tc_threshold_3way_kernel_<CBPE, MF,
                                                       HistogramsHelper3Way>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvllX2, nvllD2, nvl, vo,
          sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
          J, step_2way, param, multiplier, threshold_eff, is_using_xor,
          histograms.get_ptr(), histograms.num_buckets());

      else

        COMET_LAUNCH_KERNEL((tc_threshold_3way_kernel_<CBPE, MF,
                                                       HistogramsHelperNone>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvllX2, nvllD2, nvl, vo,
          sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
          J, step_2way, param, multiplier, threshold_eff, is_using_xor);

    } // if env.num_way()

    COMET_INSIST(System::accel_last_call_succeeded());

  } else { // (!env.is_compute_method_gpu())

    if (NumWay::_2 == env.num_way()) {

      if (is_computing_histograms)

        for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
          for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
            tc_threshold_2way_kernel_elt_<CBPE, MF, HistogramsHelper2Way>(
              thread_r, thread_c,
              nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
              param, multiplier, threshold_eff, is_using_xor,
              histograms.get_ptr(), histograms.num_buckets());
          } // for
        } // for

      else

        for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
          for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
            tc_threshold_2way_kernel_elt_<CBPE, MF, HistogramsHelperNone>(
              thread_r, thread_c,
              nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
              param, multiplier, threshold_eff, is_using_xor,
              histograms.get_ptr(), histograms.num_buckets());
          } // for
        } // for

    } else { // if (NumWay::_3 == env.num_way())

      if (is_computing_histograms)

        for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
          for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
            tc_threshold_3way_kernel_elt_<CBPE, MF, HistogramsHelper3Way>(
              thread_r, thread_c,
              nvll, nvllX2, nvllD2, nvl, vo,
              sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
              J, step_2way, param, multiplier, threshold_eff, is_using_xor,
              histograms.get_ptr(), histograms.num_buckets());
          } // for
        } // for

      else

        for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
          for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
            tc_threshold_3way_kernel_elt_<CBPE, MF, HistogramsHelperNone>(
              thread_r, thread_c,
              nvll, nvllX2, nvllD2, nvl, vo,
              sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
              J, step_2way, param, multiplier, threshold_eff, is_using_xor,
              histograms.get_ptr(), histograms.num_buckets());
          } // for
        } // for

    } // if env.num_way()

  } // if compute_method
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics if requested.

template<int METRIC_FORMAT>
void tc_threshold_(int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K,
  uint32_t* matX_counts, int J, int step_2way,
  Histograms& histograms, CEnv& env) {
  COMET_INSIST(vo);
  COMET_INSIST(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);

  COMET_INSIST(env.sparse() && "Case not supported.");
  COMET_INSIST((env.is_vectors_halved() ==
                (NumWay::_3 == env.num_way())) &&
    "Case not supported.");
  COMET_INSIST(env.is_bitwise_3way_2step() && "Case not supported.");
  COMET_INSIST_INTERFACE(&env, env.num_proc_field() == 1 &&
    "Thresholding on accelerator currently requires num_proc_field = 1.");
  COMET_INSIST(METRIC_FORMAT == MetricFormat::SINGLE);

  enum {MF = METRIC_FORMAT};

  const int cbpe = env.counted_bits_per_elt();

  if (CBPE::DUO == cbpe) {

    enum {CBPE = CBPE::DUO};

    tc_threshold_per_CBPE_<CBPE, MF>(nvll, nvl, vo, sums_I, sums_J, sums_K,
      counts_I, counts_J, counts_K, matX_counts, J, step_2way,
      histograms, env);

  } else { // if (CBPE::CCC == cbpe)

    enum {CBPE = CBPE::CCC};

    tc_threshold_per_CBPE_<CBPE, MF>(nvll, nvl, vo, sums_I, sums_J, sums_K,
      counts_I, counts_J, counts_K, matX_counts, J, step_2way,
      histograms, env);

  } // if CBPE
}

//-----------------------------------------------------------------------------
/// \brief Postprocess metrics values previously computed by GEMMs.

template<int TC_METHOD, int METRIC_FORMAT>
void tc_out_( int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K,
  uint32_t* matX_counts, int J, int step_2way,
  Histograms& histograms, CEnv& env) {
  COMET_INSIST(vo);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);

  // Perform (1) swizzle and (2) reformatting to packed double format.

  tc_repair_metrics_<TC_METHOD, METRIC_FORMAT>(nvll, nvl, vo, env);

  // Apply thresholding of smaller values to zero, if requested.

  if (env.is_threshold_tc())
    tc_threshold_<METRIC_FORMAT>(nvll, nvl, vo,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, matX_counts,
      J, step_2way, histograms, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_OUT_I_HH_

//-----------------------------------------------------------------------------
