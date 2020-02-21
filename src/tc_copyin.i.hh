//-----------------------------------------------------------------------------
/*!
 * \file   tc_copyin.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, tc package: copying data to the accelerator.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_tc_copyin_i_hh_
#define _comet_tc_copyin_i_hh_

//#include "env.hh"
//#include "tc.hh"
//#include "tc_int.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Helper function: is a nonnegative integer a power of 2.

__host__ __device__ static bool is_po2(int x) {
  return x && (!(x&(x-1))); 
}

//-----------------------------------------------------------------------------
/// \brief Formula for element to write to buf.
///
/// Description of is_bitwise_3way_2step option:
///
/// For this method two passes are made (instead of three), each of which
/// calculates exactly 4 of the required 8 metrics values.
/// The main idea is this: the left matrix entries are calculated as
/// the combined number of paths through corresponding elements of the I
/// and J matrices.  For the first pass (step_2way == 0), paths 0-0 and
/// 0-1 (corresponding to I-J element values) are counted; for the second
/// pass, paths 1-0 and 1-1.
/// Below is a tabular representation of the values for the non-sparse CCC case.
/// The "10" cases are omtted here because they are the same as the "01" cases.
/// The 3-way CCC sparse case is easliy adapted from this.
///
/// step_2way ------>   0    0    1    1
/// i01       ------>   0    1    0    1
/// output    ------>  0-0  0-1  1-0  1-1
/// --------------------------------------
///  m (=I)   c (=J)
///  |        |
///  v        v
///  00       00        4    0    0    0
///  00       01        2    2    0    0
///  00       11        2    2    0    0
/// --------------------------------------
///  01       00        2    0    2    0
///  01       01        1    1    1    1
///  01       11        0    2    0    2
/// --------------------------------------
///  11       00        0    0    4    0
///  11       01        0    0    2    2
///  11       11        0    0    0    4
///
/// The 3-way non-sparse DUO method is similar, but we only
/// look at 1 bit of each seminibble, not 2:
///
/// step_2way ------>   0    0    1    1
/// i01       ------>   0    1    0    1
/// output    ------>  0-0  0-1  1-0  1-1
/// --------------------------------------
///  m (=I)   c (=J)
///  |        |
///  v        v
///  *0       *0        1    0    0    0
///  *0       *1        0    1    0    0
/// --------------------------------------
///  *1       *0        0    0    1    0
///  *1       *1        0    0    0    1

template<typename GemmIn_t>
__host__ __device__ static GemmIn_t tc_buf_write_kernel_value_(
  const int snm,
  const int snc,
  const int i01,
  const int num_way,
  const bool is_sparse,
  const bool is_right,
  const bool is_duo,
  const bool form_matX_on_accel,
  const int step_2way,
  const bool is_bitwise_3way_2step) {

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const bool is_left = ! is_right;
  const bool skip_10 = is_sparse || (num_way == 3 && is_left);

  // Possible counts, represented in target type.

  const GemmIn_t zero = TCBufTypes<GemmIn_t>::zero();
  const GemmIn_t one  = TCBufTypes<GemmIn_t>::one();
  const GemmIn_t two  = TCBufTypes<GemmIn_t>::two();
  const GemmIn_t four = TCBufTypes<GemmIn_t>::four();

  // Possible seminibble bit patterns.

  const int _00 = 0;
  const int _01 = 1;
  const int _10 = 2;
  const int _11 = 3;

  // Unimplemented cases:
  COMET_ASSERT( ! (is_duo && !is_sparse) );
  COMET_ASSERT( ! (is_duo && !is_bitwise_3way_2step) );
  COMET_ASSERT( ! (is_bitwise_3way_2step && !form_matX_on_accel) );

  const GemmIn_t out =  3 == num_way && is_left && is_duo ? (
                          snm == _10                             ? zero :
                          snc == _10                             ? zero :
                          (snm&1) == step_2way && (snc&1) == i01 ? one :
                                                                   zero
                        ) : //====================
                        is_duo ? (
                           snm == _10         ? zero :
                          (snm & 1) == i01    ? one :
                       /* (snm & 1) == 1-i01 */ zero
                        ) : //====================
                        3 == num_way && is_left &&
                        is_bitwise_3way_2step /* && is_ccc */ ? (
                          snm == _10 && is_sparse              ? zero :
                          snc == _10 && is_sparse              ? zero :
                          snm == _11*(1-step_2way)             ? zero :
                          snc == _11*(1-i01)                   ? zero :
                          is_po2(snm) && is_po2(snc)           ? one :
                          is_po2(snm) && snc == _11*i01        ? two :
                          is_po2(snc) && snm == _11*step_2way  ? two :
                        /* snm*(3-snm) + (snc)*(3-snc) == 0 */ four
                        ) : //====================
                        3 == num_way && is_left &&
                        form_matX_on_accel /* && is_ccc */ ? (
                          snm == _10 && is_sparse      ? zero :
                          snm == _00 && step_2way != 0 ? zero :
                          snm == _01 && step_2way != 1 ? zero :
                          snm == _10 && step_2way != 1 ? zero :
                          snm == _11 && step_2way != 2 ? zero :
                          snc == _11*i01               ? two :
                          snc == _11*(1-i01)           ? zero :
                          snc == _01                   ? one :
                                  is_sparse            ? zero :
                       /* snc == _10 */                  one
                        ) : //====================
                        /* is_ccc ... */ (
                          snm == _11*i01      ? two :
                          snm == _11*(1-i01)  ? zero :
                                  !skip_10    ? one :
                          snm == _01          ? one :
                       /* snm == _10 */         zero
                        );
  return out;
}

//-----------------------------------------------------------------------------
/// \brief Write individual elements to buf.
///

template<typename GemmIn_t>
__host__ __device__ static void tc_buf_write_kernel_elt_(
  GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_on_accel,
  int step_2way,
  bool is_bitwise_3way_2step,
  int nvlea,
  int nvle,
  int nvleD2,
  int nvleX2,
  int nfl,
  int nflD2,
  int nflD2_thisstep,
  int flD2_min,
  int vlX2,
  int flD2_thisstep) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int i01 = vlX2 % 2; // count either 0 bits or 1 bits.
  COMET_ASSERT(i01 == 0 || i01 == 1);
  const int vl = vlX2 / 2;

  const int flD2 = flD2_min + flD2_thisstep;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  const uint32_t* const vim_col = vim + vl * (size_t)vi_dim0;

  // Pick up two consecutive field values:
  // first field seminibble0, second field seminibble1
  // Set to zero if outside of active range.

  const int nibblem = vl<nvlea ? (vim_col[flD2/8] >> (4*(flD2%8))) & 15 : 0;
  const int snm0 = nibblem & 3;
  const int snm1 = (nibblem>>2) & 3;

  const int nibblec = vl<nvlea ? (vic[flD2/8] >> (4*(flD2%8))) & 15 : 0;
  const int snc0 = nibblec & 3;
  const int snc1 = (nibblec>>2) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const GemmIn_t out0 = tc_buf_write_kernel_value_<GemmIn_t>(snm0, snc0, i01,
    num_way, is_sparse, is_right, is_duo, form_matX_on_accel, step_2way,
    is_bitwise_3way_2step);

  const GemmIn_t out1 = tc_buf_write_kernel_value_<GemmIn_t>(snm1, snc1, i01,
    num_way, is_sparse, is_right, is_duo, form_matX_on_accel, step_2way,
    is_bitwise_3way_2step);

  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vl_index = is_right ? vl : vl < nvleD2 ? 2*vl : 2*vl - nvle + 1;
  const int vlX2_index = i01 + 2*vl_index;

  const int flD2_index = flD2_thisstep;

  const int fl_index_0 = 0 + 2 * flD2_index;
  const int fl_index_1 = 1 + 2 * flD2_index;

  const int vlX2_dim = nvleX2;

  vo[vlX2_index + vlX2_dim * (size_t)fl_index_0] = out0;
  vo[vlX2_index + vlX2_dim * (size_t)fl_index_1] = out1;
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support tc_buf_write_.

template<typename GemmIn_t>
__global__ static void tc_buf_write_kernel_(
  GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_on_accel,
  int step_2way,
  bool is_bitwise_3way_2step,
  int nvlea,
  int nvle,
  int nvleD2,
  int nvleX2,
  int nfl,
  int nflD2,
  int nflD2_thisstep,
  int flD2_min) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int vlX2 = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int flD2_thisstep = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  if (vlX2 >= nvleX2 || flD2_thisstep >= nflD2_thisstep) {
    return;
  }

  tc_buf_write_kernel_elt_<GemmIn_t>(vo, vim, vic, vi_dim0,
    num_way, is_sparse, is_right, is_duo, form_matX_on_accel, step_2way,
    is_bitwise_3way_2step,
    nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
    vlX2, flD2_thisstep);
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

template<int TC_METHOD>
void tc_buf_write_(
  bool is_right, int I_max, int I_max_dim, int nvl,
  int npvfl, int npvfl_thisstep, int pvfl_min,
  const uint32_t* vi1, const uint32_t* vi2, TCBufs& tc_bufs, int step_2way,
  CEnv& env) {

  COMET_INSIST(vi1 && vi2);
  COMET_INSIST(I_max_dim >= 0 && I_max_dim <= nvl);
  COMET_INSIST(I_max >= 0 && I_max <= I_max_dim);
  COMET_INSIST(nvl >= 0 && npvfl >= 0);
  COMET_INSIST(tc_bufs.tc_buf_left && tc_bufs.tc_buf_right);
  COMET_INSIST(npvfl_thisstep >= 0 && npvfl_thisstep <= npvfl);
  COMET_INSIST(pvfl_min >= 0 && pvfl_min + npvfl_thisstep <= npvfl);

  // num_vector-related dimensions.

  const int nvle = is_right ? nvl : I_max_dim; // effective nvl dimension
  const int nvleD2 = nvle / 2;
  const int nvleX2 = nvle * 2;
  const int nvlea = is_right ? nvl : I_max; // num active nvle; others zeroed
  // NOTE: ignoring here the issue from decomp_mgr that
  // num_vector_active_local may be strictly less than num_vector_local;
  // doesn't matter: just compute on fake values that will later be ignored.

  COMET_INSIST(nvle % 2 == 0 && nvl % 2 == 0 &&
               "tc method here requires num_vector_local multiple of 2.");

  // num_field-related dimensions.

  const int nfl = npvfl * 64;
  const int nflD2 = nfl / 2;
  const int nfl_thisstep = npvfl_thisstep * 64;
  const int nflD2_thisstep = nfl_thisstep / 2;
  const int fl_min = pvfl_min * 64;
  const int flD2_min = fl_min / 2;
  // Remember: end padding is set to zero; will correct zero counts later.

  // Arrays.

  typedef typename TCSelector<TC_METHOD>::GemmIn_t GemmIn_t;
  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)
  GemmIn_t* const tc_buf = is_right ? (GemmIn_t*)tc_bufs.tc_buf_right :
                                      (GemmIn_t*)tc_bufs.tc_buf_left;
  COMET_INSIST(nvleX2 * (size_t)(2*nflD2_thisstep) *
           sizeof(typename TCSelector<TC_METHOD>::GemmIn_t)
           <= tc_bufs.tc_buf_size &&
           "Subscriptrange error on tc buf.");

  const bool is_duo = env.metric_type() == MetricType::DUO;
  const bool form_matX_on_accel = env.form_matX_on_accel();
  const bool is_bitwise_3way_2step = env.is_bitwise_3way_2step();

  const uint32_t* unused_col = form_matX_on_accel ? NULL : vi1; // dummy
  const uint32_t* vim = form_matX_on_accel ? vi2 : vi1; // matrix
  const uint32_t* vic = form_matX_on_accel ? vi1 : unused_col; // column

  if (env.is_compute_method_gpu()) {

    // Kernel call.

#   ifdef COMET_USE_ACCEL


#ifdef COMET_USE_CUDA
# define COMET_LAUNCH_KERNEL(name, \
    numthreadblocks, threadblocksize, sharedmem, stream, ...) \
    name <<< numthreadblocks, threadblocksize, sharedmem, stream >>> \
      (__VA_ARGS__)
#endif
#ifdef COMET_USE_HIP
# define COMET_LAUNCH_KERNEL(name, \
    numthreadblocks, threadblocksize, sharedmem, stream, ...) \
    hipLaunchKernelGGL(name, \
      numthreadblocks, threadblocksize, sharedmem, stream, __VA_ARGS__)
#endif

      const int threadblocksize = 256;
      const int blockdim_y = 32768;
      const int num_threadblocks_0 = utils::ceil(nvleX2, threadblocksize);
      const int num_threadblocks_1 = utils::min(nflD2_thisstep, blockdim_y);
      const int num_threadblocks_2 = utils::ceil(nflD2_thisstep, blockdim_y);

      COMET_LAUNCH_KERNEL((tc_buf_write_kernel_<GemmIn_t>),
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
        is_duo, form_matX_on_accel, step_2way, is_bitwise_3way_2step,
        nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min);

#if 0

#     ifdef COMET_USE_HIP
        hipLaunchKernelGGL(
#     endif
        tc_buf_write_kernel_<GemmIn_t>
#     ifdef COMET_USE_CUDA
        <<<
#     else
        ,
#     endif
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute()
#     ifdef COMET_USE_CUDA
        >>> (
#     else
        ,
#     endif
        tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
        is_duo, form_matX_on_accel, step_2way, is_bitwise_3way_2step,
        nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min);

#endif

      System::accel_last_call_succeeded();

#   endif // COMET_USE_ACCEL

  } else { // (!env.is_compute_method_gpu())

    for (int flD2_thisstep=0; flD2_thisstep<nflD2_thisstep; ++flD2_thisstep) {
      for (int vlX2=0; vlX2<nvleX2; ++vlX2) {

        tc_buf_write_kernel_elt_<GemmIn_t>(
          tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
          is_duo, form_matX_on_accel, step_2way, is_bitwise_3way_2step,
          nvlea, nvle, nvleD2, nvleX2, nfl, nflD2, nflD2_thisstep, flD2_min,
          vlX2, flD2_thisstep);

      }
//printf("========================== %i\n", flD2_thisstep);
    }

  } // if (env.is_compute_method_gpu())
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_tc_copyin_i_hh_

//-----------------------------------------------------------------------------
