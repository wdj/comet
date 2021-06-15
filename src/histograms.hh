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

#include "string"

#include "env.hh"
#include "utils.hh"
#include "mirrored_buf.hh"

//-----------------------------------------------------------------------------

namespace comet {

#if 0
TODO:

- V - 3-way histograms helper
- write output - DONE
  - readable ascii file - csv?
  - fopen, fclose
- reduce function - DONE
  - NEED to retrieve from GPU - buf_.from_accel() ??
  - ? already have code for reduction - check
- implement for non threshold_tc cases
- properly pass histograms through call chain - DONE
- CHECK HIP case - atomic add - DONE
- command line argument --histogram_file - DONE
- histogram_file_ , is_computing_histogram() - DONE
- update documentation, code_levelization.txt - DONE
- write a unit test
- function to return total number of histogram entries (ex. LLHH, LLLHHH) -
  in driver.cc try check against known metric count.
  only if is_computing_histograms.

#endif


//-----------------------------------------------------------------------------
/// \brief Helper class to specify individual histograms.

struct HistogramID {

  // 2-way
  enum {LL = 0,
        LH = 1,
        HH = 2,
        LLHH = 3};

  // 3-way
  enum {LLL = 0,
        LLH = 1,
        LHH = 2,
        HHH = 3,
        LLLHHH = 4};
};

//-----------------------------------------------------------------------------
/// \brief Class to manage histograms for metrics values.

class Histograms {

// 2-way: LL, LH, HH, LL+HH (4)
// 3-way: LLL, LLH, LHH, HHH, LLL+HHH (5)

public:

  typedef double Elt_t;

  enum {RECIP_BUCKET_WIDTH = 1000};

  Histograms(char* histograms_file, CEnv& env);

  void finalize();

  void output();

  bool is_computing_histograms() const {return is_computing_histograms_;}

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

    utils::atomic_add(&elt_this, 1e0);

  }

  Elt_t& elt(int bucket_num, int histogram_num) {
    return Histograms::elt((Elt_t*)buf_.h, num_buckets_, bucket_num, histogram_num);
  }

  Elt_t& elt_finalized(int bucket_num, int histogram_num) {
    return Histograms::elt((Elt_t*)buf_finalized_.h, num_buckets_, bucket_num,
                           histogram_num);
  }

  Elt_t bucket_min(int bucket_num) const {
    return (Elt_t)(bucket_num) / RECIP_BUCKET_WIDTH;}

private:

  CEnv& env_;

  const bool is_computing_histograms_;

  const std::string histograms_file_str_;

  const Elt_t range_;

  const int num_buckets_;

  // buf_ will represent a 2D array, dimensions are
  // num_buckets X num_histograms, and num_buckets is the stride-1 axis.
  MirroredBuf buf_;
  MirroredBuf buf_finalized_;

  bool is_finalized_;

}; // Histograms

//-----------------------------------------------------------------------------

} // namespace comet

//=============================================================================

#endif // _COMET_HISTOGRAMS_HH_

//-----------------------------------------------------------------------------
