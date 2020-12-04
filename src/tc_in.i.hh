//-----------------------------------------------------------------------------
/*!
 * \file   tc_in.i.hh
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

#ifndef _COMET_TC_IN_I_HH_
#define _COMET_TC_IN_I_HH_

//#include "env.hh"
//#include "tc.hh"
//#include "tc_int.hh"

//#include <inttypes.h>

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Word type to be used for manipulating tc_in input (fields) values.

//typedef uint32_t TCWord_t;
typedef uint64_t TCWord_t;

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

template<typename GemmIn_t, bool IS_LEFT>
__host__ __device__ static GemmIn_t tc_buf_write_kernel_value_(
  const int snm,
  const int snc,
  const int jE,
  const int kE,
  const int step_2way,
  const int num_way,
  const bool is_sparse,
  const bool is_duo,
  const bool form_matX_tc,
  const bool is_bitwise_3way_2step) {

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const bool num_way_3 = 3 == num_way;
  const bool num_way_3_left = num_way_3 && IS_LEFT;

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

template<typename GemmIn_t, bool IS_LEFT>
__host__ __device__ static GemmIn_t tc_buf_write_kernel_get_field_(
  const TCWord_t* vim,
  const TCWord_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
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

  int jE,
  int vl_thread,
  int fl) {

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  //const int jE = vlX2_thread % 2; // count either 0 or 1 bits.
  //const int vl_thread = vlX2_thread / 2;

  const int vl = is_vectors_halved && IS_LEFT ?
                 vl_thread % nvleD2 + step_2way * nvleD2 :
                 vl_thread;

  const bool is_vector_inactive = vl >= nvlea;

  const int kE = is_vectors_halved && IS_LEFT ?
                 vl_thread / nvleD2 :
                 step_2way;

  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  const TCWord_t* const vim_col = vim + vl * (size_t)vi_dim0;

  // Sizes.

  enum {BITS_PER_BYTE = 8};

  enum {BPSN = 2}; // bits per seminibble
  enum {SNPW = sizeof(TCWord_t) * BITS_PER_BYTE / BPSN}; // seminibbles / word

  // Is this fl in active range.

  const bool is_field_inactive = fl >= nfal;

  // Pick up field value.
  // Set to zero if outside of active range.

  const TCWord_t sn_m_word = is_vector_inactive ? 0 :
    vim_col[fl/SNPW];
  const TCWord_t sn_m = ( sn_m_word >> (BPSN*(fl%SNPW)) ) & 3;

  const bool num_way_3 = 3 == num_way;
  const bool num_way_3_left = num_way_3 && IS_LEFT;

  const TCWord_t sn_c_word = is_vector_inactive ? 0 : !num_way_3_left ? 0 :
    vic[fl/SNPW];
  const TCWord_t sn_c = ( sn_c_word >> (BPSN*(fl%SNPW)) ) & 3;

  // Count number of 0 (or 1) bits in respective seminibble.
  // Determine whether to skip (1,0) null indicator value.
  // NOTE: does not work for all cases.

  const GemmIn_t out = is_field_inactive ? 0 :
    tc_buf_write_kernel_value_<GemmIn_t, IS_LEFT>(sn_m, sn_c, jE, kE, step_2way,
      num_way, is_sparse, is_duo, form_matX_tc,
      is_bitwise_3way_2step);

  return out;
}

//-----------------------------------------------------------------------------
/// \brief Write individual elements to buf.
///

template<int TC_METHOD, bool IS_LEFT>
__host__ __device__ static void tc_buf_write_kernel_elt_(
  typename TCTraits<TC_METHOD>::GemmIn_t* vo,
  const TCWord_t* vim,
  const TCWord_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
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
  enum {NGIPT = TCTraits<TC_METHOD>::NGIPT};
  enum {NFPGI = TCTraits<TC_METHOD>::NFPGI};

  // Two fields (seminibbles) map to two halves of (2*sizeof(GemmIn_t))-bit word

  const int jE = vlX2_thread % 2; // count either 0 or 1 bits.
  const int vl_thread = vlX2_thread / 2;

  // Right case: straight copy of cols to cols in sequence.
  // Left case: interleave to make later swizzling of metrics array work:
  // [ A A B B C C D D E E F F ] -> [ A A D D B B E E C C F F]

  const int vl_index = (!IS_LEFT)         ?   vl_thread :
                       vl_thread < nvleD2 ? 2*vl_thread :
                                            2*vl_thread - nvle + 1;
  const int vlX2_index = jE + 2*vl_index;

  const int vlX2_dim = nvleX2_thread;

  // Output array interpreted as having GemmIn_t scalars has nfl rows.

  // "full" field-T local based on fl thread for this tc_step.
  // index into col, measured in GemmIn_t[NGIPT] sized elts.

  const int flT = flT_min + flT_thread;

  // Sizes.

  enum {BITS_PER_BYTE = 8};

  // assert(NGIPT >= 1 && NFPGI >= 1);
  // assert(NFPGI * NGIPT <= SNPW);
  // assert NGIPT is power of 2
  // assert NFPGI is power of 2

  // Loop over GemmIn_t values per this thread.
  // ISSUE: would this be faster with unroll pragma or template recursion

  const size_t nflG_dim = NGIPT * nflT_thread;

  //----------
  // Special-case code for a B1 case.

  if (TC_METHOD == TC::B1 && is_duo && NGIPT == 1 && sizeof(GemmIn_t) == 4 &&
      sizeof(TCWord_t) == 8 && TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR) {

    const int vl = is_vectors_halved && IS_LEFT ?
                   vl_thread % nvleD2 + step_2way * nvleD2 :
                   vl_thread;

    const bool is_vector_inactive = vl >= nvlea;

    const int kE = is_vectors_halved && IS_LEFT ?
                   vl_thread / nvleD2 :
                   step_2way;

    const TCWord_t* const vim_col = vim + vl * (size_t)vi_dim0;

    // Pick up fields value.
    // Set to zero if outside of active range.

    const TCWord_t m_word = is_vector_inactive ? 0 : vim_col[flT];

    const bool num_way_3 = 3 == num_way;
    const bool num_way_3_left = num_way_3 && IS_LEFT;

    const TCWord_t c_word = is_vector_inactive ? 0 : !num_way_3_left ? 0 :
      vic[flT];

    // Calculate target result on odd bits.

    const TCWord_t oddbits_mask = 0x5555555555555555;

    enum {BPSN = 2}; // bits per seminibble (field)
    enum {SNPW = sizeof(TCWord_t) * BITS_PER_BYTE / BPSN}; // seminibbles (fields) / word

    // Upper bound on fields computed by this thread.
    const int fl_max = (flT+1) * SNPW;

    // How many inactive fields in this word.
    const int fl_inactive = utils::min(32, utils::max(0, fl_max - nfal));

    const TCWord_t allbits = ~(TCWord_t)0;

    // NOTE: 64-bit word fills up the rightmost bits first - see vectors.hh.
    const TCWord_t field_active_mask = fl_inactive >= 32 ?
      ((TCWord_t)0) : allbits >> (2 * fl_inactive);

    // Use odd bits to denote whether each field is undef or not.
    const TCWord_t m_notundef_mask = (m_word | ~(m_word >> 1));
    const TCWord_t c_notundef_mask = (c_word | ~(c_word >> 1));

    const TCWord_t vo_value_64 = num_way_3_left ?
      (kE ? m_word : ~m_word) & (jE ? c_word : ~c_word) &
      oddbits_mask & field_active_mask & m_notundef_mask & c_notundef_mask :
      (jE ? m_word : ~m_word) &
      oddbits_mask & field_active_mask & m_notundef_mask;

    // Store.

    const size_t flG_index = 0 + NGIPT * flT_thread;

    GemmIn_t& vo_value = vo[flG_index + nflG_dim * vlX2_index];

    // Put odd bits of 64-bit word into 32-bit word.
    const GemmIn_t* const vo_value_32 = (GemmIn_t*)&vo_value_64;
    vo_value = vo_value_32[0] | (vo_value_32[1] << 1);

    return;
  }

  //----------
  // Special-case code for a HIP case.

  if (BuildHas::HIP && TC_METHOD == TC::INT8) {

    const int vl = is_vectors_halved && IS_LEFT ?
                   vl_thread % nvleD2 + step_2way * nvleD2 :
                   vl_thread;

    const bool is_vector_inactive = vl >= nvlea;

    const int kE = is_vectors_halved && IS_LEFT ?
                   vl_thread / nvleD2 :
                   step_2way;

    const TCWord_t* const vim_col = vim + vl * (size_t)vi_dim0;

    enum {BPSN = 2}; // bits per seminibble
    enum {SNPW = sizeof(TCWord_t) * BITS_PER_BYTE / BPSN}; // seminibbles / word

    for (int igipt = 0; igipt < NGIPT; igipt+=4) {

      const size_t flG_index = igipt + NGIPT * flT_thread;

      const int nb = 4;

      GemmIn_t* vo_value =
        &(vo[flG_index % nb + nb * (vlX2_index + vlX2_dim * (flG_index / nb))]);

      const int fl = igipt + NGIPT * flT;

      const bool num_way_3 = 3 == num_way;
      const bool num_way_3_left = num_way_3 && IS_LEFT;

      const TCWord_t sn_m_word = is_vector_inactive ? 0 : vim_col[fl/SNPW];
      const TCWord_t sn_c_word = is_vector_inactive ? 0 : !num_way_3_left ? 0 :
        vic[fl/SNPW];

// TODO: check that this is properly unrolled.

#     pragma unroll
      for (int i = 0; i < (int)(sizeof(TCWord_t)/sizeof(GemmIn_t)); ++i) {

        vo_value[i] = fl+i >= nfal ? 0 :
          tc_buf_write_kernel_value_<GemmIn_t, IS_LEFT>(
            ( sn_m_word >> (BPSN*((fl+i)%SNPW)) ) & 3,
            ( sn_c_word >> (BPSN*((fl+i)%SNPW)) ) & 3,
            jE, kE, step_2way,
            num_way, is_sparse, is_duo, form_matX_tc,
            is_bitwise_3way_2step);

      } // i
    } // igipt

    return;

  } // if (BuildHas::HIP && TC_METHOD == TC::INT8)

  //----------

  for (int igipt = 0; igipt < NGIPT; ++igipt) {

    const size_t flG_index = igipt + NGIPT * flT_thread;

    const int nb = 4;

    // see https://rocblas.readthedocs.io/en/latest/functions.html?highlight=rocblas_gemm_ex#blas-extensions

    GemmIn_t& vo_value =
      TCTraits<TC_METHOD>::IS_B_FIELD_MAJOR ?
        vo[flG_index + nflG_dim * vlX2_index] :
      (BuildHas::HIP && TC_METHOD == TC::INT8) ?
        vo[flG_index % nb + nb * (vlX2_index + vlX2_dim * (flG_index / nb))] :
        vo[vlX2_index + vlX2_dim * flG_index];

    GemmIn_t vo_value_tmp = vo_value;

    for (int ifpgi = 0; ifpgi < NFPGI; ++ifpgi) {

      // Get active field_local number.

      const int fl = ifpgi + NFPGI * (igipt + NGIPT * flT);

      // Count number of 0 (or 1) bits in respective seminibble.
      // Determine whether to skip (1,0) null indicator value.
      // NOTE: does not work for all cases.

      const GemmIn_t out =
        tc_buf_write_kernel_get_field_<GemmIn_t, IS_LEFT>(vim, vic, vi_dim0,
          num_way, is_sparse, is_duo, form_matX_tc, step_2way,
          is_bitwise_3way_2step, is_vectors_halved,
          nvle, nvleD2, nvleX2_thread, nvlea,
          nfl, nfal, 
          jE, vl_thread, fl);

      // Store.

      if (ifpgi) {
        const int shift = ifpgi * (BITS_PER_BYTE * sizeof(GemmIn_t) / NFPGI);
        vo_value_tmp = vo_value_tmp | ( out << shift);
      } else {
        vo_value_tmp = out;
      }

    } // ifpgi

    vo_value = vo_value_tmp;

  } // igipt
}

//-----------------------------------------------------------------------------
/// \brief GPU kernel to support tc_buf_write_.

template<int TC_METHOD, bool IS_LEFT>
__global__ static void tc_buf_write_kernel_(
  typename TCTraits<TC_METHOD>::GemmIn_t* vo,
  const TCWord_t* vim,
  const TCWord_t* vic,
  int vi_dim0,
  int num_way,
  bool is_sparse,
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

  const int thread_dim0 = threadIdx_x_() + blockIdx_x_() * blockDim_x_();
  const int thread_dim1 = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  enum {IS_THREAD_MAPPING_FIELD_MAJOR =
        TCTraits<TC_METHOD>::IS_THREAD_MAPPING_FIELD_MAJOR};

  const int flT_thread = IS_THREAD_MAPPING_FIELD_MAJOR ?
    thread_dim0 : thread_dim1;
  const int vlX2_thread = IS_THREAD_MAPPING_FIELD_MAJOR ?
    thread_dim1 : thread_dim0;

  if (vlX2_thread >= nvleX2_thread || flT_thread >= nflT_thread)
    return;

  tc_buf_write_kernel_elt_<TC_METHOD, IS_LEFT>(vo, vim, vic, vi_dim0,
    num_way, is_sparse, is_duo, form_matX_tc, step_2way,
    is_bitwise_3way_2step, is_vectors_halved,
    nvle, nvleD2, nvleX2_thread, nvlea, nfl, nflT, nflT_thread, flT_min, nfal,
    vlX2_thread, flT_thread);
}

//-----------------------------------------------------------------------------
/// \brief Convert bitwise matrix to required format for GEMM.

template<int TC_METHOD, bool IS_LEFT>
void tc_buf_write_(
  int I_max, int I_max_dim, int nvl,
  int npvfl, int npvfl_thisstep, int pvfl_min, int nfal,
  const TCWord_t* vi1, const TCWord_t* vi2, TCBufs& tc_bufs,
  int step_2way, CEnv& env) {

  COMET_INSIST(vi1 && vi2);
  COMET_INSIST(I_max_dim >= 0 && I_max_dim <= nvl);
  COMET_INSIST(I_max >= 0 && I_max <= I_max_dim);
  COMET_INSIST(nvl >= 0 && npvfl >= 0);
  COMET_INSIST(tc_bufs.tc_buf_left && tc_bufs.tc_buf_right);
  COMET_INSIST(npvfl_thisstep >= 0 && npvfl_thisstep <= npvfl);
  COMET_INSIST(pvfl_min >= 0 && pvfl_min + npvfl_thisstep <= npvfl);

  // num_vector-related dimensions.

  const int nvle = IS_LEFT ? I_max_dim : nvl; // effective nvl dimension
  const int nvleD2 = nvle / 2;
  const int nvleX2 = nvle * 2;
  const int nvlea = IS_LEFT ? I_max : nvl; // num active nvle; others zeroed
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

  enum {NGIPT = TCTraits<TC_METHOD>::NGIPT};
  enum {NFPGI = TCTraits<TC_METHOD>::NFPGI};

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

  const int vi_dim0 = npvfl * (2*sizeof(double)) / sizeof(TCWord_t);
  //const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)

  GemmIn_t* const tc_buf = IS_LEFT ? (GemmIn_t*)tc_bufs.tc_buf_left :
                                     (GemmIn_t*)tc_bufs.tc_buf_right;
  COMET_INSIST(nvleX2 * (size_t)(nGI_thisstep) * sizeof(GemmIn_t)
           <= tc_bufs.tc_buf_size &&
           "Subscriptrange error on tc buf.");

  const bool is_duo = env.metric_type() == MetricType::DUO;
  const bool form_matX_tc = env.form_matX_tc();
  const bool is_bitwise_3way_2step = env.is_bitwise_3way_2step();

  const TCWord_t* unused_col = form_matX_tc ? NULL : vi1; // dummy
  const TCWord_t* vim = form_matX_tc ? vi2 : vi1; // matrix
  const TCWord_t* vic = form_matX_tc ? vi1 : unused_col; // column

  const int nvleX2_thread = nvleX2;
  const int nflT_thread = nflT_thisstep;

  enum {IS_THREAD_MAPPING_FIELD_MAJOR =
        TCTraits<TC_METHOD>::IS_THREAD_MAPPING_FIELD_MAJOR};

  const int thread_dim0 = IS_THREAD_MAPPING_FIELD_MAJOR ?
    nflT_thread : nvleX2_thread;
  const int thread_dim1 = IS_THREAD_MAPPING_FIELD_MAJOR ?
    nvleX2_thread : nflT_thread;

  if (env.is_compute_method_gpu()) {

    // Kernel call.

      const int threadblocksize = 256;
      COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
                   "Current HIP limitation.");
      const int blockdim_y = 32768;
      const int num_threadblocks_0 = utils::ceil(thread_dim0, threadblocksize);
      const int num_threadblocks_1 = utils::min(thread_dim1, blockdim_y);
      const int num_threadblocks_2 = utils::ceil(thread_dim1, blockdim_y);

      COMET_LAUNCH_KERNEL((tc_buf_write_kernel_<TC_METHOD, IS_LEFT>),
        dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
        dim3(threadblocksize, 1, 1), 0, env.stream_compute(),
        tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(),
        is_duo, form_matX_tc, step_2way, is_bitwise_3way_2step,
        env.is_vectors_halved(),
        nvle, nvleD2, nvleX2_thread, nvlea,
        nfl, nflT, nflT_thread, flT_min, nfal);

      System::accel_last_call_succeeded();

  } else { // (!env.is_compute_method_gpu())

    for (int i1 = 0; i1 < thread_dim1; ++i1) {
      for (int i0 = 0; i0 < thread_dim0; ++i0) {

        const int flT_thread = IS_THREAD_MAPPING_FIELD_MAJOR ? i0 : i1;
        const int vlX2_thread = IS_THREAD_MAPPING_FIELD_MAJOR ? i1 : i0;

        tc_buf_write_kernel_elt_<TC_METHOD, IS_LEFT>(
          tc_buf, vim, vic, vi_dim0, env.num_way(), env.sparse(),
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
/// \brief GPU kernel to support tc_compute_matX_counts.

template<typename GemmIn_t>
__global__ static void tc_compute_matX_counts_kernel_(
  uint32_t* matX_counts,
  const TCWord_t* vim,
  const TCWord_t* vic,
  int vi_dim0,
  int step_2way,
  int nvle,
  int nvleD2,
  int nvleX2_thread,
  int nvlea,
  int nfl,
  int nfl_thread,
  int nfal) {

#ifdef COMET_USE_ACCEL
  const int vlX2_thread = blockIdx_y_() + gridDim_y_() * blockIdx_z_();

  if (vlX2_thread >= nvleX2_thread)
    return;

  const int jE = vlX2_thread % 2; // count either 0 or 1 bits.
  const int vl_thread = vlX2_thread / 2;

  const int fl_ind0 = threadIdx_x_();
  const int fl_dim0 = blockDim_x_();
  const int fl_ind1 = blockIdx_x_();
  const int fl_dim1 = gridDim_x_();

  const int fl_thread = fl_ind0 + fl_dim0 * fl_ind1;

  // NOTE: here nfl and nfal can be larger than the number of fl threads.

  const int num_way = NumWay::_3;
  const bool is_sparse = true;
  const bool is_duo = true;
  const bool form_matX_tc = true;
  const bool is_bitwise_3way_2step = true;
  const bool is_vectors_halved = true;
  enum {IS_LEFT = true};

  // Adapted from https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

  extern __shared__ uint32_t sdata[];

  sdata[fl_ind0] = 0;
  int fl = fl_thread;

  // Reduce across threadblocks along fl direction.

  while (fl < nfal) {
    const GemmIn_t elt = 
      tc_buf_write_kernel_get_field_<GemmIn_t, IS_LEFT>(vim, vic, vi_dim0,
        num_way, is_sparse, is_duo, form_matX_tc, step_2way,
        is_bitwise_3way_2step, is_vectors_halved,
        nvle, nvleD2, nvleX2_thread, nvlea,
        nfl, nfal, 
        jE, vl_thread, fl);

    sdata[fl_ind0] += elt;
    fl += fl_dim1 * fl_dim0;
  }

  // Reduce within threadblock.

  for (unsigned int s = fl_dim0/2; s > 0; s >>= 1) {
    __syncthreads();
    if (fl_ind0 < s) // parallel sweep reduction
      sdata[fl_ind0] += sdata[fl_ind0 + s];
  }

  // First thread of threadblock adds in its contribution.

  if (fl_ind0 == 0) {
    atomicAdd(&(matX_counts[vlX2_thread]), sdata[0]);
  }
#endif
}

//-----------------------------------------------------------------------------
/// \brief Compute vector (column) sums for matX to support 3-way xor gemm.

void tc_compute_matX_counts(
  int I_max,
  int I_max_dim,
  int nvl,
  int npvfl,
  int nfal,
  const TCWord_t* vi1,
  const TCWord_t* vi2,
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

  const TCWord_t* vim = vi2; // matrix
  const TCWord_t* vic = vi1; // column

  enum {NUM_FL_PER_UINT32 = 16};

  const int nvleX2_thread = nvleX2;
  const int nfl_thread = nfl / NUM_FL_PER_UINT32;
  const int vi_dim0 = npvfl * (2*sizeof(double)) / sizeof(TCWord_t);
  //const int vi_dim0 = npvfl * 4; // 4 = sizeof(doublecomplex) / sizeof(int32)

  typedef TCTraits<TC::B1>::GemmIn_t GemmIn_t;

  uint32_t* const matX_counts = tc_bufs.matX_counts;

# if defined COMET_USE_CUDA
    cudaMemsetAsync(matX_counts, 0, tc_bufs.matX_counts_size,
      env.stream_compute());
# elif defined COMET_USE_HIP
    hipMemsetAsync(matX_counts, 0, tc_bufs.matX_counts_size,
      env.stream_compute());
# endif

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
    threadblocksize * sizeof(int32_t),
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

#endif // _COMET_TC_IN_I_HH_

//-----------------------------------------------------------------------------
