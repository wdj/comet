//-----------------------------------------------------------------------------
/*!
 * \file   histograms.hh
 * \author Wayne Joubert
 * \date   Fri Jun 11 08:15:19 EDT 2021
 * \brief  Manage histograms, declarations.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2021, UT-Battelle, LLC

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

#ifndef _COMET_HISTOGRAMS_HH_
#define _COMET_HISTOGRAMS_HH_

#include "env.hh"
#include "mirrored_buf.hh"

//-----------------------------------------------------------------------------

namespace comet {

#if 0
TODO:

- V - 3way histograms helper
- write output
  - readable ascii file - csv?
  - fopen, fclose
- reduce function
  - NEED to retrieve from GPU - buf_.from_accel() ??
  - ? already have code for reduction - check
- properly pass histograms through call chain
  - metrics -> decomp_mgr -> tc_bufs
- implement for non threshold_tc cases
- CHECK HIP case - atomic add - ??
- command line argument --histogram_file
- env.histogram_file_ , env.is_computing_histogram()
- update documentation, code_levelization.txt
- write a unit test

#endif


//-----------------------------------------------------------------------------

struct HistogramID {

  enum {LL = 0,
        LH = 1,
        HH = 2,
        LLHH = 3};

  enum {LLL = 0,
        LLH = 1,
        LHH = 2,
        HHH = 3,
        LLLHHH = 4};
};

class Histograms {

// 2-way: LL, LH, HH, LL+HH (4)
// 3-way: LLL, LLH, LHH, HHH, LLL+HHH (5)

public:

  typedef double Elt_t;

  enum {RECIP_BUCKET_WIDTH = 1000};

  Histograms(CEnv& env);

  void reduce();

  void output();

  int num_histograms() const {return env_.num_way() + 2;}

  Elt_t* get_ptr() {return (Elt_t*)(buf_.active);}

  int num_buckets() const {return num_buckets_;}

  // Access an entry of the histogram
  static __host__ __device__ Elt_t& elt(Elt_t* ptr, int num_buckets,
    int bucket_num, int histogram_num) {
    return ptr[bucket_num + num_buckets * histogram_num];
  }

  // Add value to correct bucket of the histogram.
  template<typename T>
  static __host__ __device__ void add(Elt_t* ptr, int num_buckets,
    T value, int histogram_id) {

    const int bucket_num = value * RECIP_BUCKET_WIDTH;

    const int bucket_num_clamped = bucket_num <= 0 ? 0 :
                                   bucket_num > num_buckets-1 ? num_buckets-1 :
                                   bucket_num;

    Elt_t& elt_this = elt(ptr, num_buckets, bucket_num_clamped, histogram_id);

#   if defined COMET_USE_CUDA && defined __CUDA_ARCH__

      atomicAdd(&elt_this, 1e0);

#   elif defined COMET_USE_HIP && defined __HIPCC__

//FIX - use atomic CAS ?? does hip have atomic update double ?
      atomicAdd(&elt_this, 1e0);

#   else

      elt_this++;

#   endif

  }
 

private:

  CEnv& env_;

  Elt_t range_;

  int num_buckets_;

  // buf_ will represent a 2D array, dimensions are
  // num_buckets X num_histograms, and num_buckets is the stride-1 axis.
  MirroredBuf buf_;


}; // Histograms

//-----------------------------------------------------------------------------

} // namespace comet

//=============================================================================

#endif // _COMET_HISTOGRAMS_HH_

//-----------------------------------------------------------------------------
