//-----------------------------------------------------------------------------
/*!
 * \file   tc_copyout.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, tc package: copying data from the accelerator.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_tc_copyout_i_hh_
#define _comet_tc_copyout_i_hh_

#include "formulas.hh"

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
//printf("%f %f %f %f %f %f %f %f\n", (float)i00, (float)i01, (float)i02, (float)i03, (float)i10, (float)i11, (float)i12, (float)i13);
//printf("%i %i %i %i %i %i %i %i\n", (int)i00, (int)i01, (int)i02, (int)i03, (int)i10, (int)i11, (int)i12, (int)i13);

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

//  // Use "shifter" to move a value to the upper half of the mantissa.
//
//  const double shifter = (((uint32_t)1) << GM_TALLY1_MAX_VALUE_BITS);
//
//  // Pack two 26-bit integers into mantissa of double.
//
//  const double o00 = (double)i00p + (double)i02p * shifter;
//  const double o01 = (double)i01p + (double)i03p * shifter;
//
//  const double o10 = (double)i10p + (double)i12p * shifter;
//  const double o11 = (double)i11p + (double)i13p * shifter;

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
///        2. Each column (of 2-bit entries) is converted into two columns,
///        with entries being the counts of 0 bits and 1 bits of the
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

  // always true, because of tc_gemm_divisibility_required()
  COMET_INSIST(nvll % 2 == 0 && "Failed divisibility condition for tc gemm.");
  const int nvllD2 = nvll / 2;

  typedef typename TCSelector<TC_METHOD>::GemmOut_t GemmOut_t;

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

      System::accel_last_call_succeeded();

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
/// \brief Do thresholding of metrics: individual elements, 2-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT>
__host__ __device__ void tc_threshold_2way_kernel_elt_(
  int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* counts_I, GMFloat* counts_J,
  double param, double multiplier,
  double threshold_eff, int thread_r, int thread_c) {
  COMET_ASSERT(vo && sums_I && sums_J && counts_I && counts_J);
  COMET_ASSERT(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_ASSERT(thread_r >= 0 && thread_r < nvll);
  COMET_ASSERT(thread_c >= 0 && thread_c < nvl);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type MFType;
  typedef typename MFT::TypeIn MFTypeIn;

  MFType* const dvo = (MFType*)vo;

  // Indexing.

  // The incoming matrix is of dimension nvll X nvl.
  // Each matrix element is 4 table entries, stored as 16 bytes.
  // The matrix is stored column major.
  // I, J - coordinate of element in the 2D plane.
  // indT_I, indT_J - indices into the 2x2 table at each coord

  // TODO: use MirroredBuf class for indexing.

  // Calculate cijk.

  const int mftype_per_tallytable = Tally2x2<METRIC_FORMAT>::NUM;

  GMTally1 cij = 0;
  for (int indT_I = 0; indT_I < 2; ++indT_I) {

    const size_t index_this = indT_I + mftype_per_tallytable * (
                              thread_r + nvll * (
                              thread_c));

    const MFType dvo_this = dvo[index_this];

    MFTypeIn values_this[2];
    MFT::decode(values_this[0], values_this[1], dvo_this);

    cij += (GMTally1)values_this[0] + (GMTally1)values_this[1];
  }

  const double recip_sumcij = 1e0 / cij;

  const int I = thread_r;
  COMET_ASSERT(I >= 0 && I < nvll);

  const int J = thread_c;
  COMET_ASSERT(J >= 0 && J < nvl);

  // (more) values to be plugged into CCC/DUO formula.

  const GMTally1 cI = (GMTally1)counts_I[I];
  const GMTally1 cJ = (GMTally1)counts_J[J];

  const double recip_cI = 1e0 / cI;
  const double recip_cJ = 1e0 / cJ;

  const GMTally1 sI1 = (GMTally1)sums_I[I];
  const GMTally1 sJ1 = (GMTally1)sums_J[J];

  // Loops to update all 8 table values.

  for (int indT_I = 0; indT_I < 2; ++indT_I) {

    const GMTally1 sI = indT_I == 0 ? CBPE * cI - sI1 : sI1;

    const size_t index_this = indT_I + mftype_per_tallytable * (
                              thread_r + nvll * (
                              thread_c));

    MFType& dvo_this = dvo[index_this];

    MFTypeIn values_this[2];
    MFT::decode(values_this[0], values_this[1], dvo_this);

    for (int indT_J = 0; indT_J < 2; ++indT_J) {

      const GMTally1 sJ = indT_J == 0 ? CBPE * cJ - sJ1 : sJ1;

      const MFTypeIn rij = values_this[indT_J];

      const double metric_value =
        0 == cI || 0 == cJ || 0 == cij ? 0e0 :
        ccc_duo_value<CBPE, double>(
          (GMTally1)rij, sI, sJ, recip_cI, recip_cJ,
          recip_sumcij, multiplier, param);

      const bool pass_threshold = CEnv::pass_threshold(metric_value,
                                                       threshold_eff);

      values_this[indT_J] = pass_threshold ? (MFTypeIn)metric_value :
                                             (MFTypeIn)0;
    } // indT_J

    MFT::encode(dvo_this, values_this[0], values_this[1]);
  } // indT_I
}

//=============================================================================
/// \brief Do thresholding of metrics: individual elements, 3-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT>
__host__ __device__ void tc_threshold_3way_kernel_elt_(
  int nvll, int nvllX2, int nvllD2,  int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int step_2way, double param, double multiplier,
  double threshold_eff, int thread_r, int thread_c) {
  COMET_ASSERT(vo);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_ASSERT(nvll*2 == nvllX2);
  COMET_ASSERT(nvll/2 == nvllD2);
  COMET_ASSERT(nvll >= 0 && nvl >= 0 && nvll <= nvl);
  COMET_ASSERT(thread_r >= 0 && thread_r < nvllD2);
  COMET_ASSERT(thread_c >= 0 && thread_c < nvl);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type MFType;
  typedef typename MFT::TypeIn MFTypeIn;

  MFType* const dvo = (MFType*)vo;

  // Indexing.

  // The incoming matrix is of dimension nvll X nvl.
  // Each matrix element is 4 table entries, stored as 16 bytes.
  // The matrix is stored column major.
  // Recall the matrix is "halved": the lower and upper half pertain
  // to different table entries (indT_I axis) for the same (I,J,K) coordinate
  // I, J, K - (permuted) coordinate of element in the 3D block
  // indT_I, indT_J, indT_K - indices into the 2x2x2 table at each coord

  // Calculate cijk.

  const int mftype_per_tallytable = Tally2x2<METRIC_FORMAT>::NUM;
  const size_t halves_per_whole = 2;

  GMTally1 cijk = 0;
  for (int indT_I = 0; indT_I < 2; ++indT_I) {
    for (int indT_J = 0; indT_J < 2; ++indT_J) {

      const MFType dvo_this = dvo[indT_J + mftype_per_tallytable * (
                                  thread_r + nvllD2 * (
                                  indT_I + halves_per_whole * (
                                  thread_c)))];

      MFTypeIn values_this[2];
      MFT::decode(values_this[0], values_this[1], dvo_this);

      cijk += (GMTally1)values_this[0] + (GMTally1)values_this[1];
    }
  }

  const double recip_sumcijk = 1e0 / cijk;

  const int I = thread_r + nvllD2 * step_2way;
  COMET_ASSERT(I >= 0 && I < nvll);
  COMET_ASSERT(I/nvllD2 == step_2way);

  const int K = thread_c;
  COMET_ASSERT(K >= 0 && K < nvl);

  // (more) values to be plugged into CCC/DUO formula.

  const GMTally1 cI = (GMTally1)counts_I[I];
  const GMTally1 cJ = (GMTally1)counts_J[J];
  const GMTally1 cK = (GMTally1)counts_K[K];

  const double recip_cI = 1e0 / cI;
  const double recip_cJ = 1e0 / cJ;
  const double recip_cK = 1e0 / cK;

  const GMTally1 sI1 = (GMTally1)sums_I[I];
  const GMTally1 sJ1 = (GMTally1)sums_J[J];
  const GMTally1 sK1 = (GMTally1)sums_K[K];

  // Loops to update all 8 table values.

  for (int indT_I = 0; indT_I < 2; ++indT_I) {

    const GMTally1 sI = indT_I == 0 ? CBPE * cI - sI1 : sI1;

    for (int indT_J = 0; indT_J < 2; ++indT_J) {

      const GMTally1 sJ = indT_J == 0 ? CBPE * cJ - sJ1 : sJ1;

      MFType& dvo_this = dvo[indT_J + mftype_per_tallytable * (
                             thread_r + nvllD2 * (
                             indT_I + halves_per_whole * (
                             thread_c)))];

      MFTypeIn values_this[2];
      MFT::decode(values_this[0], values_this[1], dvo_this);

      for (int indT_K = 0; indT_K < 2; ++indT_K) {

        const GMTally1 sK = indT_K == 0 ? CBPE * cK - sK1 : sK1;

        const MFTypeIn rijk = values_this[indT_K];

        const double metric_value =
          0 == cI || 0 == cJ || 0 == cK || 0 == cijk ? 0e0 :
          ccc_duo_value<CBPE, double>(
            (GMTally1)rijk, sI, sJ, sK, recip_cI, recip_cJ, recip_cK,
            recip_sumcijk, multiplier, param);

        const bool pass_threshold = CEnv::pass_threshold(metric_value,
                                                         threshold_eff);

        values_this[indT_K] = pass_threshold ? (MFTypeIn)metric_value :
                                               (MFTypeIn)0;
      } // indT_K

      MFT::encode(dvo_this, values_this[0], values_this[1]);
    } // indT_J
  } // indT_I
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics: GPU kernel, 2-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT>
__global__ void tc_threshold_2way_kernel_(
  int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* counts_I, GMFloat* counts_J,
  double param, double multiplier, double threshold_eff) {
  COMET_ASSERT(vo && sums_I && sums_J && counts_I && counts_J);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  // Row and column threads for metrics matrix.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvll || thread_c >= nvl)
    return;

  tc_threshold_2way_kernel_elt_<COUNTED_BITS_PER_ELT, METRIC_FORMAT>(
    nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
    param, multiplier, threshold_eff, thread_r, thread_c);
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics: GPU kernel, 3-way case.

template<int COUNTED_BITS_PER_ELT, int METRIC_FORMAT>
__global__ void tc_threshold_3way_kernel_(
  int nvll, int nvllX2, int nvllD2,  int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int step_2way, double param, double multiplier, double threshold_eff) {
  COMET_ASSERT(vo);
  COMET_ASSERT(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_ASSERT(METRIC_FORMAT == MetricFormat::SINGLE);

  // Row and column threads for metrics matrix.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvllD2 || thread_c >= nvl)
    return;

  tc_threshold_3way_kernel_elt_<COUNTED_BITS_PER_ELT, METRIC_FORMAT>(
    nvll, nvllX2, nvllD2, nvl, vo,
    sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
    step_2way, param, multiplier, threshold_eff, thread_r, thread_c);
}

//-----------------------------------------------------------------------------
/// \brief Do thresholding of metrics if requested.

template<int TC_METHOD, int METRIC_FORMAT>
void tc_threshold_(int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int step_2way, CEnv& env) {
  COMET_INSIST(vo);
  COMET_INSIST(sums_I && sums_J && sums_K && counts_I && counts_J && counts_K);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);

  COMET_INSIST(env.sparse() && "Case not supported.");
  COMET_INSIST((env.is_vectors_halved() ==
                (NUM_WAY::_3 == env.num_way())) &&
    "Case not supported.");
  COMET_INSIST(env.is_bitwise_3way_2step() && "Case not supported.");
  COMET_INSIST_INTERFACE(&env, env.num_proc_field() == 1 &&
    "Thresholding on accelerator currently requires num_proc_field = 1.");
  COMET_INSIST(METRIC_FORMAT == MetricFormat::SINGLE);

  enum {MF = METRIC_FORMAT};

  const int nvllX2 = nvll * 2;
  const int nvllD2 = nvll / 2;

  const int cbpe = env.counted_bits_per_elt();

  // ? template specialization on num_way
  // 2way: i -> I, j -> K (???)

  const double param = env.ccc_param();
  const double multiplier = env.ccc_duo_multiplier();
  const double threshold_eff = env.threshold_eff();

  const int num_threads_r = env.num_way() == NUM_WAY::_2 ? nvll : nvllD2;
  const int num_threads_c = nvl;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

      const int threadblocksize = 256;
      COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                   "Current HIP limitation.");
      const int vll_threadblocks = utils::ceil(num_threads_r, threadblocksize);

      if (CBPE::DUO == cbpe && NUM_WAY::_2 == env.num_way()) {

        enum {CBPE = CBPE::DUO};

        COMET_LAUNCH_KERNEL((tc_threshold_2way_kernel_<CBPE, MF>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
          param, multiplier, threshold_eff);

      } else if (NUM_WAY::_2 == env.num_way()) { // && CBPE::CCC == cbpe

        enum {CBPE = CBPE::CCC};

        COMET_LAUNCH_KERNEL((tc_threshold_2way_kernel_<CBPE, MF>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
          param, multiplier, threshold_eff);

      } else if (CBPE::DUO == cbpe) { // && (3 == env.num_way())

        enum {CBPE = CBPE::DUO};

        COMET_LAUNCH_KERNEL((tc_threshold_3way_kernel_<CBPE, MF>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvllX2, nvllD2, nvl, vo,
          sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
          step_2way, param, multiplier, threshold_eff);

      } else { // if (CBPE::CCC == cbpe && NUM_WAY::_3 == env.num_way())

        enum {CBPE = CBPE::CCC};

        COMET_LAUNCH_KERNEL((tc_threshold_3way_kernel_<CBPE, MF>),
          dim3(vll_threadblocks, num_threads_c, 1),
          dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
          nvll, nvllX2, nvllD2, nvl, vo,
          sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
          step_2way, param, multiplier, threshold_eff);

      } // if cbpe, env.num_way()

      System::accel_last_call_succeeded();

  } else { // (!env.is_compute_method_gpu())

    if (CBPE::DUO == cbpe && NUM_WAY::_2 == env.num_way()) {

      enum {CBPE = CBPE::DUO};

      for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
        for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
          tc_threshold_2way_kernel_elt_<CBPE, MF>(
            nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
            param, multiplier, threshold_eff, thread_r, thread_c);
        } // for
      } // for

    } else if (NUM_WAY::_2 == env.num_way()) { // && CBPE::CCC == cbpe

      enum {CBPE = CBPE::CCC};

      for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
        for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
          tc_threshold_2way_kernel_elt_<CBPE, MF>(
            nvll, nvl, vo, sums_I, sums_J, counts_I, counts_J,
            param, multiplier, threshold_eff, thread_r, thread_c);
        } // for
      } // for

    } else if (CBPE::DUO == cbpe) { // && (NUM_WAY::_3 == env.num_way())

      enum {CBPE = CBPE::DUO};

      for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
        for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
          tc_threshold_3way_kernel_elt_<CBPE, MF>(
            nvll, nvllX2, nvllD2, nvl, vo,
            sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
            step_2way, param, multiplier, threshold_eff, thread_r, thread_c);
        } // for
      } // for

    } else { // if (CBPE::CCC == cbpe && NUM_WAY::_3 == env.num_way())

      enum {CBPE = CBPE::CCC};

      for (int thread_c=0; thread_c<num_threads_c; ++thread_c) {
        for (int thread_r=0; thread_r<num_threads_r; ++thread_r) {
          tc_threshold_3way_kernel_elt_<CBPE, MF>(
            nvll, nvllX2, nvllD2, nvl, vo,
            sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J,
            step_2way, param, multiplier, threshold_eff, thread_r, thread_c);
        } // for
      } // for

    } // if cbpe, env.num_way()

  } // if compute_method
}

//-----------------------------------------------------------------------------
/// \brief Postprocess metrics values previously computed by GEMMs.

template<int TC_METHOD, int METRIC_FORMAT>
void tc_out_( int nvll, int nvl, void* vo,
  GMFloat* sums_I, GMFloat* sums_J, GMFloat* sums_K,
  GMFloat* counts_I, GMFloat* counts_J, GMFloat* counts_K, int J,
  int step_2way,
  CEnv& env) {
  COMET_INSIST(vo);
  COMET_INSIST(nvll >= 0 && nvl >= 0 && nvll <= nvl);

  // Perform (1) swizzle and (2) reformatting to packed double format.

  tc_repair_metrics_<TC_METHOD, METRIC_FORMAT>(nvll, nvl, vo, env);

  // Apply thresholding of smaller values to zero, if requested.

  if (env.threshold_tc())
    tc_threshold_<TC_METHOD, METRIC_FORMAT>(nvll, nvl, vo,
      sums_I, sums_J, sums_K, counts_I, counts_J, counts_K, J, step_2way, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_tc_copyout_i_hh_

//-----------------------------------------------------------------------------