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

//#include "env.hh"
//#include "tc.hh"
//#include "tc_int.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Swizzle individual elements in buf.

template<typename GemmOut_t>
__host__ __device__ static void tc_repair_metrics_kernel_elt_(
  int nvl, int nvll, int nvllD2, void* vo,
  int thread_r, int thread_c) { 

  // Considered as an array of floats, array is 2*nvl rows X 2*nvl cols.
  // Each thread manipulates a block of 4 rows and 2 cols.
  // Thus the dimensions of the metrics array in blocks is nvllD2 X nvl.
  // Each block viewed as an array of doubles is 2 X 2.

  // Two col numbers being processed of this (float) array.

  // ISSUE: does the compiler need to / understand that the pointers are aliased

  const size_t fcr_offset0 = 4*thread_r + thread_c * (size_t)(4*nvll);
  const size_t fcr_offset1 = 4*thread_r + thread_c * (size_t)(4*nvll) + 2*nvll;

  // Read the 8 values.

  GemmOut_t* const fvo = (GemmOut_t*)vo;

  const GemmOut_t f00 = fvo[fcr_offset0+0];
  const GemmOut_t f01 = fvo[fcr_offset0+1];
  const GemmOut_t f02 = fvo[fcr_offset0+2];
  const GemmOut_t f03 = fvo[fcr_offset0+3];

  const GemmOut_t f10 = fvo[fcr_offset1+0];
  const GemmOut_t f11 = fvo[fcr_offset1+1];
  const GemmOut_t f12 = fvo[fcr_offset1+2];
  const GemmOut_t f13 = fvo[fcr_offset1+3];
//printf("%f %f %f %f %f %f %f %f\n", (float)f00, (float)f01, (float)f02, (float)f03, (float)f10, (float)f11, (float)f12, (float)f13);

  // Apply the permutation:

  // [ f00  f10 ]  ->  [ f00  f02 ]
  // [ f01  f11 ]  ->  [ f01  f03 ]
  // [ f02  f12 ]  ->  [ f10  f12 ]
  // [ f03  f13 ]  ->  [ f11  f13 ]

  const GemmOut_t f00p = f00;
  const GemmOut_t f01p = f01;

  const GemmOut_t f02p = f10;
  const GemmOut_t f03p = f11;

  const GemmOut_t f10p = f02;
  const GemmOut_t f11p = f03;

  const GemmOut_t f12p = f12;
  const GemmOut_t f13p = f13;

  // Use "shifter" to move a value to the upper half of the mantissa.

  const double shifter = (((uint32_t)1) << GM_TALLY1_MAX_VALUE_BITS);

  // Pack two 26-bit integers into mantissa of double.

  const double d00 = (double)f00p + (double)f02p * shifter;
  const double d01 = (double)f01p + (double)f03p * shifter;

  const double d10 = (double)f10p + (double)f12p * shifter;
  const double d11 = (double)f11p + (double)f13p * shifter;

  // Overwrite block with the new values.
  // All is isolated to a single thread, should be thread safe.

  const size_t dc_offset0 = thread_c * (size_t)(2*nvll);
  const size_t dc_offset1 = thread_c * (size_t)(2*nvll) + nvll;

  const size_t dcr_offset0 = dc_offset0 + 2*thread_r;
  const size_t dcr_offset1 = dc_offset1 + 2*thread_r;

  double* const dvo = (double*)vo;

  dvo[dcr_offset0+0] = d00;
  dvo[dcr_offset0+1] = d01;

  dvo[dcr_offset1+0] = d10;
  dvo[dcr_offset1+1] = d11;
}

//-----------------------------------------------------------------------------
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

template<typename GemmOut_t>
__global__ static void tc_repair_metrics_kernel_(
  int nvl, int nvll, int nvllD2, void* vo) { 

  // Row and column threads of metrics array.
  const int thread_r = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_c = blockIdx_y_();

  if (thread_r >= nvllD2 || thread_c >= nvl) {
    return;
  }

  tc_repair_metrics_kernel_elt_<GemmOut_t>(
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

template<int TC_METHOD>
void tc_repair_metrics_(
  int nvll,
  int nvl,
  void* vo,
  TCBufs& tc_bufs,
  CEnv& env) {

  COMET_INSIST(vo);
  COMET_INSIST(nvll >= 0);
  COMET_INSIST(nvl >= 0);
  COMET_INSIST(nvll <= nvl);

  // always true, because of tc_gemm_divisibility_required()
  COMET_INSIST(nvll % 2 == 0 && "Failed divisibility condition for tc gemm.");
  const int nvllD2 = nvll / 2;

  typedef typename TCSelector<TC_METHOD>::GemmOut_t GemmOut_t;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

#ifdef COMET_USE_ACCEL

    const int threadblocksize = 256;
    COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                 "Current HIP limitation.");
    const int vll2_threadblocks = utils::ceil(nvllD2, threadblocksize);

#  ifdef COMET_USE_HIP
    hipLaunchKernelGGL(
#  endif
    tc_repair_metrics_kernel_<GemmOut_t>
#  ifdef COMET_USE_CUDA
        <<<
#  else
        ,
#  endif
        dim3(vll2_threadblocks, nvl, 1),
        dim3(threadblocksize, 1, 1),
        0,
        env.stream_compute()
#  ifdef COMET_USE_CUDA
        >>> (
#  else
        ,
#  endif
        nvl, nvll, nvllD2, vo);

    System::accel_last_call_succeeded();

#endif // COMET_USE_ACCEL

  } else { // (!env.is_compute_method_gpu())

//    for (int thread_c=0; thread_c<nvllD2; ++thread_c) {
//      for (int thread_r=0; thread_r<nvl; ++thread_r) {
    for (int thread_c=0; thread_c<nvl; ++thread_c) {
      for (int thread_r=0; thread_r<nvllD2; ++thread_r) {

        tc_repair_metrics_kernel_elt_<GemmOut_t>(
          nvl, nvll, nvllD2, vo, thread_r, thread_c);

      }
    }

  } // if compute_method
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_tc_copyout_i_hh_

//-----------------------------------------------------------------------------
