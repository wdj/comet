//-----------------------------------------------------------------------------
/*!
 * \file   tc_copyin.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  CUDA code, tc package: copying data to the accelerator.
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

#ifndef _COMET_TC_COPYIN_I_HH_
#define _COMET_TC_COPYIN_I_HH_

//#include "env.hh"
//#include "tc.hh"
//#include "tc_int.hh"

//#include <inttypes.h>

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Helper function: is a nonnegative integer a power of 2.

template<typename T>
__host__ __device__ static bool is_po2(T x) {
  return x && (!(x&(x-1))); 
}

//-----------------------------------------------------------------------------
/// \brief Formula for element to write to buf.
///
/// Description of is_bitwise_3way_2step option, for 3-way case:
///
/// For this method two passes are made (instead of three), each of which
/// calculates exactly 4 of the required 8 metrics values.
/// The main idea is this: the left matrix elements are calculated as
/// the combined number of paths through corresponding elements of the I
/// and J matrices.  For the first pass (step_2way (= kE) == 0), paths 0-0 and
/// 0-1 (corresponding to I-J element values) are counted; for the second
/// pass, paths 1-0 and 1-1.
/// Below is a tabular representation of the values for the non-sparse CCC case.
/// The "10" cases are omtted here because they are the same as the "01" cases.
/// The 3-way CCC sparse case is easliy adapted from this.
///
/// kE        ------>   0    0    1    1
/// jE        ------>   0    1    0    1
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
/// kE        ------>   0    0    1    1
/// jE        ------>   0    1    0    1
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
  const int jE,
  const int kE,
  const int step_2way,
  const int num_way,
  const bool is_sparse,
  const bool is_right,
  const bool is_duo,
  const bool form_matX_tc,
  const bool is_bitwise_3way_2step) {

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const bool is_left = ! is_right;
  const bool num_way_3 = 3 == num_way;
  const bool num_way_3_left = num_way_3 && is_left;

  // Possible counts, represented in target type.

  const GemmIn_t zero = TCBufTraits<GemmIn_t>::zero();
  const GemmIn_t one  = TCBufTraits<GemmIn_t>::one();
  const GemmIn_t two  = TCBufTraits<GemmIn_t>::two();
  const GemmIn_t four = TCBufTraits<GemmIn_t>::four();

  // Possible seminibble bit patterns.

  const int _00 = 0;
  const int _01 = 1;
  const int _10 = 2;
  const int _11 = 3;
  const int _UNDEF = _10;

  // Test for unimplemented cases:
  COMET_ASSERT( ! (is_duo && !is_sparse) &&
                "DUO only implemented for sparse case");
  COMET_ASSERT( ! (is_duo && !is_bitwise_3way_2step) &&
                "DUO 3way only implemented for 2step case");
  COMET_ASSERT( ! (is_bitwise_3way_2step && !form_matX_tc) &&
               "3way 2step requires form_matX_tc");

  const GemmIn_t out =
    is_duo && num_way_3_left ? (

             // Form X: combine the column and matrix
             snm == _UNDEF || snc == _UNDEF ? zero :
             (snm&1) == kE && (snc&1) == jE ? one :
                                              zero

    ) : is_duo /* && is_sparse */ ? (

             // DUO: pick up low order bit (unless undefined)
             snm == _UNDEF    ? zero :
             (snm&1) == jE    ? one :
          /* (snm&1) == 1-jE */ zero

    ) : num_way_3_left && is_bitwise_3way_2step
        /* && is_ccc && form_matX_tc */ ? (

             // Form X: combine the column and matrix
             snm == _UNDEF && is_sparse          ? zero :
             snc == _UNDEF && is_sparse          ? zero :
             snm == _11*(1-kE)                   ? zero :
             snc == _11*(1-jE)                   ? zero :
             is_po2(snm) && is_po2(snc)          ? one :
             is_po2(snm) && snc == _11*jE        ? two :
             is_po2(snc) && snm == _11*kE        ? two :
          /* snm*(3-snm) + (snc)*(3-snc) == 0 */   four

    ) : num_way_3_left && form_matX_tc
        /* && is_ccc && !is_bitwise_3way_2step */ ? (

             // Encode based on 2way step num
             snm == _UNDEF && is_sparse   ? zero :
             snm == _00 && step_2way != 0 ? zero :
             snm == _01 && step_2way != 1 ? zero :
             snm == _10 && step_2way != 1 ? zero :
             snm == _11 && step_2way != 2 ? zero :
             snc == _11*jE                ? two :
             snc == _11*(1-jE)            ? zero :
             snc == _01                   ? one :
             is_sparse                    ? zero :
          /* snc == _10 */                  one

    ) : /* is_ccc ... */ (

             snm == _11*jE                   ? two :
             snm == _11*(1-jE)               ? zero :
             snm == _01                      ? one :
             snm == _UNDEF && is_sparse      ? zero :
          /* snm == _UNDEF && num_way_3_left ? zero :*/
          /* snm == _10 && ... */              one

    );

  return out;
}

//-----------------------------------------------------------------------------
/// \brief Access tc_buf value to write that is associated with specific field.
///

template<typename GemmIn_t>
__host__ __device__ static GemmIn_t tc_buf_write_kernel_get_field_(
  GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_tc,
  int step_2way,
  bool is_bitwise_3way_2step,
  bool is_vectors_halved,

  int nvle,
  int nvleD2,
  int nvleX2_thread,
  int nvlea,

  int nfl,
  int nfal,

  int vlX2_thread,
  int fl) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int jE = vlX2_thread % 2; // count either 0 or 1 bits.
  const int vl_thread = vlX2_thread / 2;

  const int vl = is_vectors_halved && !is_right ?
                 vl_thread % nvleD2 + step_2way * nvleD2 :
                 vl_thread;

  const bool is_vector_inactive = vl >= nvlea;

  const int kE = is_vectors_halved && !is_right ?
                 vl_thread / nvleD2 :
                 step_2way;

  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

//  const int vl_index = is_right           ?   vl_thread :
//                       vl_thread < nvleD2 ? 2*vl_thread :
//                                            2*vl_thread - nvle + 1;
//  const int vlX2_index = jE + 2*vl_index;

//  const int vlX2_dim = nvleX2_thread;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  const uint32_t* const vim_col = vim + vl * (size_t)vi_dim0;

  // Sizes.

  enum {BITS_PER_BYTE = 8};
  enum {SNPW = 16}; // seminibbles per 32-bit word
  enum {BPSN = 2}; // bits per seminibble
  // assert(SNPW * BPSN == sizeof(*vim) * BITS_PER_BYTE);

  // Is this fl in active range.

  const bool is_field_inactive = fl >= nfal;

  // Pick up field value.
  // Set to zero if outside of active range.

  const uint32_t sn_m = is_vector_inactive ? 0 :
    ( (vim_col[fl/SNPW]) >> (BPSN*(fl%SNPW)) ) & 3;

  const uint32_t sn_c = is_vector_inactive ? 0 : /* is_right ? 0 : */
    ( (vic    [fl/SNPW]) >> (BPSN*(fl%SNPW)) ) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const GemmIn_t out = is_field_inactive ? 0 :
    tc_buf_write_kernel_value_<GemmIn_t>(sn_m, sn_c, jE, kE, step_2way,
      num_way, is_sparse, is_right, is_duo, form_matX_tc,
      is_bitwise_3way_2step);

  return out;
}

//-----------------------------------------------------------------------------
/// \brief Write individual elements to buf.
///

template<int TC_METHOD, int NGIPT, int NFPGI>
__host__ __device__ static void tc_buf_write_kernel_elt_(
  typename TCTraits<TC_METHOD>::GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_tc,
  int step_2way,
  bool is_bitwise_3way_2step,
  bool is_vectors_halved,

  int nvle,
  int nvleD2,
  int nvleX2_thread,
  int nvlea,

  int nfl,
  int nflT,
  int nflT_thread,
  int flT_min,
  int nfal,

  int vlX2_thread,
  int flT_thread) {

  typedef typename TCTraits<TC_METHOD>::GemmIn_t GemmIn_t;

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int jE = vlX2_thread % 2; // count either 0 or 1 bits.
  const int vl_thread = vlX2_thread / 2;

//  const int vl = is_vectors_halved && !is_right ?
//                 vl_thread % nvleD2 + step_2way * nvleD2 :
//                 vl_thread;

//  const bool is_vector_inactive = vl >= nvlea;

//  const int kE = is_vectors_halved && !is_right ?
//                 vl_thread / nvleD2 :
//                 step_2way;

  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vl_index = is_right           ?   vl_thread :
                       vl_thread < nvleD2 ? 2*vl_thread :
                                            2*vl_thread - nvle + 1;
  const int vlX2_index = jE + 2*vl_index;

  const int vlX2_dim = nvleX2_thread;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

//  const uint32_t* const vim_col = vim + vl * (size_t)vi_dim0;

  // "full" field-T local based on fl thread for this tc_step.
  // index into col, measured in GemmIn_t[NGIPT] sized elts.

  const int flT = flT_min + flT_thread;

  // Sizes.

  enum {BITS_PER_BYTE = 8};
//  enum {SNPW = 16}; // seminibbles per 32-bit word
//  enum {BPSN = 2}; // bits per seminibble
  // assert(SNPW * BPSN == sizeof(*vim) * BITS_PER_BYTE);

  // assert(NGIPT >= 1 && NFPGI >= 1);
  // assert(NFPGI * NGIPT <= SNPW);
  // assert NGIPT is power of 2
  // assert NFPGI is power of 2

  // Loop over GemmIn_t values per this thread.
  // ISSUE: would this be faster with unroll pragma or template recursion

  const size_t nflG_dim = NGIPT * nflT_thread;

  for (int igipt = 0; igipt < NGIPT; ++igipt) {

    const size_t flG_index = igipt + NGIPT * flT_thread;

    const int nb = 4;

    // see https://rocblas.readthedocs.io/en/latest/functions.html?highlight=rocblas_gemm_ex#blas-extensions

    GemmIn_t& vo_value =
      TC_METHOD == TC::B1 ?
        vo[flG_index + nflG_dim * vlX2_index] :
      (BuildHas::HIP && TC_METHOD == TC::INT8) ?
        vo[flG_index % nb + nb * (vlX2_index + vlX2_dim * (flG_index / nb))] :
        vo[vlX2_index + vlX2_dim * flG_index];

    for (int ifpgi = 0; ifpgi < NFPGI; ++ifpgi) {

      // Get active field_local number.

      const int fl = ifpgi + NFPGI * (igipt + NGIPT * flT);

      // Is this fl in active range.

//      const bool is_field_inactive = fl >= nfal;

      // Pick up field value.
      // Set to zero if outside of active range.


//      const uint32_t sn_m = is_vector_inactive ? 0 :
//        ( (vim_col[fl/SNPW]) >> (BPSN*(fl%SNPW)) ) & 3;

//      const uint32_t sn_c = is_vector_inactive ? 0 : /* is_right ? 0 : */
//        ( (vic    [fl/SNPW]) >> (BPSN*(fl%SNPW)) ) & 3;

      // Count number of 0 (or 1) bits in respective seminibble.
      // Determine whether to skip (1,0) null indicator value.
      // NOTE: does not work for all cases.

//      const GemmIn_t out = is_field_inactive ? 0 :
//        tc_buf_write_kernel_value_<GemmIn_t>(sn_m, sn_c, jE, kE, step_2way,
//          num_way, is_sparse, is_right, is_duo, form_matX_tc,
//          is_bitwise_3way_2step);

      const GemmIn_t out =
        tc_buf_write_kernel_get_field_<GemmIn_t>(vo, vim, vic, vi_dim0,
          num_way, is_sparse, is_right, is_duo, form_matX_tc, step_2way,
          is_bitwise_3way_2step, is_vectors_halved,
          nvle, nvleD2, nvleX2_thread, nvlea,
          nfl, nfal, 
          vlX2_thread, fl);

      // Store.

      if (ifpgi) {
        const int shift = ( ifpgi * BITS_PER_BYTE * sizeof(GemmIn_t) ) / NFPGI;
        vo_value = vo_value | ( out << shift);
      } else {
        vo_value = out;
      }

    } // ifpgi
  } // igipt
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support tc_buf_write_.

template<int TC_METHOD, int NGIPT, int NFPGI>
__global__ static void tc_buf_write_kernel_(
  typename TCTraits<TC_METHOD>::GemmIn_t* vo,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
  bool is_right,
  bool is_duo,
  bool form_matX_tc,
  int step_2way,
  bool is_bitwise_3way_2step,
  bool is_vectors_halved,

  int nvle,
  int nvleD2,
  int nvleX2_thread,
  int nvlea,

  int nfl,
  int nflT,
  int nflT_thread,
  int flT_min,
  int nfal) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int vlX2_thread = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int flT_thread = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  if (vlX2_thread >= nvleX2_thread || flT_thread >= nflT_thread)
    return;

  tc_buf_write_kernel_elt_<TC_METHOD, NGIPT, NFPGI>(vo, vim, vic, vi_dim0,
    num_way, is_sparse, is_right, is_duo, form_matX_tc, step_2way,
    is_bitwise_3way_2step, is_vectors_halved,
    nvle, nvleD2, nvleX2_thread, nvlea, nfl, nflT, nflT_thread, flT_min, nfal,
    vlX2_thread, flT_thread);
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

template<int TC_METHOD>
void tc_buf_write_(
  bool is_right, int I_max, int I_max_dim, int nvl,
  int npvfl, int npvfl_thisstep, int pvfl_min, int nfal,
  const uint32_t* vi1, const uint32_t* vi2, TCBufs& tc_bufs,
  int step_2way, CEnv& env) {

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

  enum {NUM_FIELD_PER_PACKEDVAL_FIELD = 64};
  enum {NUM_FL_PER_PVFL = NUM_FIELD_PER_PACKEDVAL_FIELD};

  const int nfl = npvfl * NUM_FL_PER_PVFL;
  const int nfl_thisstep = npvfl_thisstep * NUM_FL_PER_PVFL;
  const int fl_min = pvfl_min * NUM_FL_PER_PVFL;

  enum {NUM_GEMMIN_T_PER_THREAD = 2}; // tuning param
  enum {NGIPT = NUM_GEMMIN_T_PER_THREAD};

  enum {NUM_FIELD_PER_GEMMIN_T = TCTraits<TC_METHOD>::NFPGI};
  enum {NFPGI = NUM_FIELD_PER_GEMMIN_T};

  enum {NUM_FIELD_PER_THREAD = NFPGI * NGIPT};
  enum {NFPT = NUM_FIELD_PER_THREAD};

  //const int nGI = nfl / NFPGI;
  const int nGI_thisstep = nfl_thisstep / NFPGI;
  //const int GI_min = fl_min / NFPGI;

  const int nflT = nfl / NFPT;
  const int nflT_thisstep = nfl_thisstep / NFPT;
  const int flT_min = fl_min / NFPT;

  // Remember: end padding is set to zero; will correct zero counts later.

  // Arrays.

  typedef typename TCTraits<TC_METHOD>::GemmIn_t GemmIn_t;
  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)
  GemmIn_t* const tc_buf = is_right ? (GemmIn_t*)tc_bufs.tc_buf_right :
                                      (GemmIn_t*)tc_bufs.tc_buf_left;
  COMET_INSIST(nvleX2 * (size_t)(nGI_thisstep) *
           sizeof(typename TCTraits<TC_METHOD>::GemmIn_t)
           <= tc_bufs.tc_buf_size &&
           "Subscriptrange error on tc buf.");

  const bool is_duo = env.metric_type() == MetricType::DUO;
  const bool form_matX_tc = env.form_matX_tc();
  const bool is_bitwise_3way_2step = env.is_bitwise_3way_2step();

  const uint32_t* unused_col = form_matX_tc ? NULL : vi1; // dummy
  const uint32_t* vim = form_matX_tc ? vi2 : vi1; // matrix
  const uint32_t* vic = form_matX_tc ? vi1 : unused_col; // column

  const int nvleX2_thread = nvleX2;
  const int nflT_thread = nflT_thisstep;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

      const int threadblocksize = 256;
      COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                   "Current HIP limitation.");
      const int blockdim_y = 32768;
      const int num_threadblocks_0 = utils::ceil(nvleX2_thread, threadblocksize);
      const int num_threadblocks_1 = utils::min(nflT_thread, blockdim_y);
      const int num_threadblocks_2 = utils::ceil(nflT_thread, blockdim_y);

      COMET_LAUNCH_KERNEL((tc_buf_write_kernel_<TC_METHOD, NGIPT, NFPGI>),
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
        is_duo, form_matX_tc, step_2way, is_bitwise_3way_2step,
        env.is_vectors_halved(),
        nvle, nvleD2, nvleX2_thread, nvlea,
        nfl, nflT, nflT_thread, flT_min, nfal);

      System::accel_last_call_succeeded();

  } else { // (!env.is_compute_method_gpu())

    for (int flT_thread=0; flT_thread<nflT_thread; ++flT_thread) {
      for (int vlX2_thread=0; vlX2_thread<nvleX2_thread; ++vlX2_thread) {

        tc_buf_write_kernel_elt_<TC_METHOD, NGIPT, NFPGI>(
          tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(), is_right,
          is_duo, form_matX_tc, step_2way, is_bitwise_3way_2step,
          env.is_vectors_halved(),
          nvle, nvleD2, nvleX2_thread, nvlea,
          nfl, nflT, nflT_thread, flT_min, nfal,
          vlX2_thread, flT_thread);

      }
    }

  } // if (env.is_compute_method_gpu())
}

//=============================================================================

template<typename GemmIn_t>
__global__ static void tc_compute_matX_counts_kernel_(
  uint32_t* matX_counts,
  const uint32_t* vim,
  const uint32_t* vic,
  int vi_dim0,
  int step_2way,
  int nvle,
  int nvleD2,
  int nvleX2_thread,
  int nvlea,
  int nfl,
  int nfl_thread,
  int nfal) {

#if xxx




  const int fl_thread0 = threadIdx_x_();
  const int fl_thread = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int vlX2_thread = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  if (vlX2_thread >= nvleX2_thread || fl_thread >= nfl_thread)
    return;

  const int num_way = NumWay::_3;
  const bool is_sparse = true;
  const bool is_right = false;
  const bool is_duo = true;
  const bool form_matX_tc = true;
  const bool is_bitwise_3way_2step = true;
  const bool is_vectors_halved = true;

  __shared__ uint32_t sdata[];

  int tid = fl_thread0;
  sdata[tid] = 0;
  size_t idx = fl_thread;

  while (idx < N) { // grid stride loop to load data
    sdata[tid] += gdata[idx];
    idx += gridDim.x * blockDim.x;
  }

  for (unsigned int s=blockDim.x/2; s>0; s>>=1) {
    __syncthreads();
    if (tid < s) // parallel sweep reduction
      sdata[tid] += sdata[tid + s];
    }

  if (tid == 0) atomicAdd(out, sdata[0]);



#endif
}

//-----------------------------------------------------------------------------

void tc_compute_matX_counts(
  int I_max,
  int I_max_dim,
  int nvl,
  int npvfl,
  int nfal,
  const uint32_t* vi1,
  const uint32_t* vi2,
  TCBufs& tc_bufs,
  int step_2way,
  CEnv& env) {

  COMET_INSIST(vi1 && vi2);
  COMET_INSIST(I_max >= 0 && I_max <= I_max_dim);

  if (!(env.is_using_xor() && env.num_way() == NumWay::_3))
    return;

  COMET_INSIST(env.is_compute_method_gpu());
  COMET_INSIST(env.num_way() == NumWay::_3);
  COMET_INSIST(env.sparse());
  COMET_INSIST(env.metric_type() == MetricType::DUO);
  COMET_INSIST(env.form_matX_tc());
  COMET_INSIST(env.is_bitwise_3way_2step());
  COMET_INSIST(env.is_vectors_halved());

#ifdef COMET_USE_ACCEL

  // num_vector-related dimensions.

  const int nvle = I_max_dim; // effective nvl dimension
  const int nvleD2 = nvle / 2;
  const int nvleX2 = nvle * 2;
  const int nvlea = I_max; // num active nvle; others zeroed

  enum {NUM_FL_PER_PVFL = 64};
  const int nfl = npvfl * NUM_FL_PER_PVFL;
  COMET_INSIST(nfal <= nfl);

  const uint32_t* vim = vi2; // matrix
  const uint32_t* vic = vi1; // column

  enum {NUM_FL_PER_UINT32 = 16};

  const int nvleX2_thread = nvleX2;
  const int nfl_thread = nfl / NUM_FL_PER_UINT32;
  const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)

  typedef TCTraits<TC::B1>::GemmIn_t GemmIn_t;

  uint32_t* const matX_counts = tc_bufs.matX_counts;

  const int threadblocksize = 256;
  COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
               "Current HIP limitation.");
  const int blockdim_y = 32768;
  const int num_threadblocks_0 = utils::ceil(nfl_thread, threadblocksize);
  const int num_threadblocks_1 = utils::min(nvleX2_thread, blockdim_y);
  const int num_threadblocks_2 = utils::ceil(nvleX2_thread, blockdim_y);

  COMET_LAUNCH_KERNEL((tc_compute_matX_counts_kernel_<GemmIn_t>),
    dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
    dim3(threadblocksize, 1, 1),
    threadblocksize * sizeof(uint32_t),
    env.stream_compute(),
    matX_counts, vim, vic, vi_dim0, step_2way,
    nvle, nvleD2, nvleX2_thread, nvlea,
    nfl, nfl_thread, nfal);

  System::accel_last_call_succeeded();

#endif //COMET_USE_ACCEL
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_COPYIN_I_HH_

//-----------------------------------------------------------------------------
