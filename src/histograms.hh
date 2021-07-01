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

#include "types.hh"
#include "env.hh"
#include "utils.hh"
#include "formulas.hh"
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
- properly pass histograms through call chain - DONE
- CHECK HIP case - atomic add - DONE
- command line argument --histogram_file - DONE
- histogram_file_ , is_computing_histogram() - DONE
- update documentation, code_levelization.txt - DONE
- function to return total number of histogram entries (ex. LLHH, LLLHHH) -
  in driver.cc try check against known metric count.
  only if is_computing_histograms. - DONE

- write a unit test

- implement for non threshold_tc cases

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

  Histograms(const char* histograms_file, CEnv& env);

  void finalize();

  void output();

  bool is_computing_histograms() const {return is_computing_histograms_;}

  bool is_computing_on_accel() const {return env_.is_threshold_tc();}

  int num_histograms() const {return num_histograms_;}

  Elt_t* get_ptr() {return (Elt_t*)(buf_.active);}

  int num_buckets() const {return num_buckets_;}

  //----------

  // access an entry of the histogram (static)
  static __host__ __device__ Elt_t& elt(Elt_t* ptr, int num_buckets,
    int bucket_num, int histogram_num) {
    return ptr[bucket_num + num_buckets * histogram_num];
  }

  //----------

  // access an entry of the histogram (static) (const)
  static __host__ __device__ Elt_t elt_const(Elt_t* ptr, int num_buckets,
    int bucket_num, int histogram_num) {
    return ptr[bucket_num + num_buckets * histogram_num];
  }

  //----------

  // Add value to correct bucket of the histogram (static).
  template<typename T>
  static __host__ __device__ void add(Elt_t* ptr, int num_buckets,
    T value, int histogram_id) {

    const int bucket_num = value * RECIP_BUCKET_WIDTH;

    const int bucket_num_clamped = bucket_num <= 0 ? 0 :
                                   bucket_num > num_buckets-1 ? num_buckets-1 :
                                   bucket_num;

    Elt_t& elt_this = elt(ptr, num_buckets, bucket_num_clamped, histogram_id);

    // TODO: make this OMP 3.1 thread aware.
    utils::atomic_add(&elt_this, 1e0);

  }

  //----------

  // Add value to correct bucket of the histogram (non-static, host only).
  template<typename T>
  void add(T value, int histogram_id) {
    COMET_ASSERT(!is_computing_on_accel());
    Histograms::add((Elt_t*)buf_.h, num_buckets_, value, histogram_id);
  }

  //----------

  // Add 2X2 table entries to histogram buckets.
  template<int MF>
  void add(const Tally2x2<MF> ttable, GMTally1 si1, GMTally1 sj1,
           GMTally1 ci, GMTally1 cj, int nfa) {
    if (!is_computing_histograms_)
      return;

    const int cbpe = env_.counted_bits_per_elt();

    double metric[2][2] = {};

    for (int iE=0; iE<2; ++iE) {
      for (int jE=0; jE<2; ++jE) {

        metric[iE][jE] = cbpe == CBPE::DUO ?
          ccc_duo_value<CBPE::DUO, MF>(ttable, iE, jE, si1, sj1,
                                       ci, cj, nfa, env_) :
          ccc_duo_value<CBPE::CCC, MF>(ttable, iE, jE, si1, sj1,
                                       ci, cj, nfa, env_);
      }
    }

    add(metric[0][0], HistogramID::LL);
    add(metric[0][1], HistogramID::LH);
    add(metric[1][0], HistogramID::LH);
    add(metric[1][1], HistogramID::HH);
    add(metric[0][0] + metric[1][1],
                      HistogramID::LLHH);
  }

  //----------

  // Add 4X2 table entries to histogram buckets.
  template<int MF>
  void add(const Tally4x2<MF> ttable, GMTally1 si1, GMTally1 sj1, GMTally1 sk1,
           GMTally1 ci, GMTally1 cj, GMTally1 ck, int nfa) {
    if (!is_computing_histograms_)
      return;

    const int cbpe = env_.counted_bits_per_elt();

    double metric[2][2][2] = {};

    for (int iE=0; iE<2; ++iE) {
      for (int jE=0; jE<2; ++jE) {
        for (int kE=0; kE<2; ++kE) {
          metric[iE][jE][kE] = cbpe == CBPE::DUO ?
            ccc_duo_value<CBPE::DUO, MF>(ttable, iE, jE, kE, si1, sj1, sk1,
                                         ci, cj, ck,  nfa, env_) :
            ccc_duo_value<CBPE::CCC, MF>(ttable, iE, jE, kE, si1, sj1, sk1,
                                         ci, cj, ck,  nfa, env_);
        }
      }
    }

    add(metric[0][0][0], HistogramID::LLL);
    add(metric[0][0][1], HistogramID::LLH);
    add(metric[0][1][0], HistogramID::LLH);
    add(metric[1][0][0], HistogramID::LLH);
    add(metric[0][1][1], HistogramID::LHH);
    add(metric[1][0][1], HistogramID::LHH);
    add(metric[1][1][0], HistogramID::LHH);
    add(metric[1][1][1], HistogramID::HHH);
    add(metric[0][0][0] + metric[1][1][1],
                         HistogramID::LLLHHH);
  }

  //----------

  // accessor for element of histogram
  Elt_t& elt(int bucket_num, int histogram_num) {
    return Histograms::elt((Elt_t*)buf_.h, num_buckets_, bucket_num, histogram_num);
  }

  //----------

  // accessor for element of histogram (const)
  Elt_t elt_const(int bucket_num, int histogram_num) const {
    return Histograms::elt_const((Elt_t*)buf_.h, num_buckets_, bucket_num, histogram_num);
  }

  //----------

  // accessor for element of finalized histogram
  Elt_t& elt_finalized(int bucket_num, int histogram_num) {
    return Histograms::elt((Elt_t*)buf_finalized_.h, num_buckets_, bucket_num,
                           histogram_num);
  }

  //----------

  // accessor for element of finalized histogram (const)
  Elt_t elt_finalized_const(int bucket_num, int histogram_num) const {
    COMET_ASSERT(is_finalized_);
    return Histograms::elt_const((Elt_t*)buf_finalized_.h, num_buckets_, bucket_num,
                           histogram_num);
  }
  //----------

  // Return the lowest real-number value assigned to specified bucket.
  Elt_t bucket_min(int bucket_num) const {
    return (Elt_t)(bucket_num) / RECIP_BUCKET_WIDTH;}

  //----------

  void check(size_t num_vector);

private:

  CEnv& env_;

  const bool is_computing_histograms_;

  const std::string histograms_file_str_;

  const Elt_t range_;

  const int num_buckets_;

  const int num_histograms_;

  const int num_elts_;

  // buf_ will represent a 2D array, dimensions are
  // num_buckets X num_histograms, and num_buckets is the stride-1 axis.
  MirroredBuf buf_;
  MirroredBuf buf_finalized_;

  bool is_finalized_;

  // Compute sum of all (finalized) histogram elements.
  Elt_t sum_() const;
#if 0
  Elt_t sum_() {
    COMET_INSIST(is_finalized_);
    Elt_t result = 0;
    for (int col = 0; col < num_histograms_; ++col) {
      const int multiplier = env_.num_way() + 1 == col ? 0 : 1;
      for (int row = 0; row < num_buckets_ ; ++row) {
        result += multiplier * elt_finalized_const(row, col);
      }
    }
    return result;
  }
#endif

}; // Histograms

//-----------------------------------------------------------------------------

} // namespace comet

//=============================================================================

#endif // _COMET_HISTOGRAMS_HH_

//-----------------------------------------------------------------------------
