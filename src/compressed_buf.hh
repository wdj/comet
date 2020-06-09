//-----------------------------------------------------------------------------
/*!
 * \file   compressed_buf.hh
 * \author Wayne Joubert
 * \date   Tue May  5 20:11:12 EDT 2020
 * \brief  Mirrored buffer allowing for compression.
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

#ifndef _COMET_COMPRESSED_BUF_HH_
#define _COMET_COMPRESSED_BUF_HH_

#include <type_traits>

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Helper function for access to buffer elements.

// Forward declaration.
class CompressedBuf;

template<typename T>
struct CompressedBufAccessorUnspecialized_ {
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
  //static T elt_const(size_t ind_entry, const CompressedBuf* buf);
};

template<typename T>
struct CompressedBufAccessor_ {
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
  static T elt_const(size_t ind_entry, const CompressedBuf* buf);
};

template<>
struct CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>> {
  typedef Tally2x2<MetricFormat::SINGLE> T;
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
  static T elt_const(size_t ind_entry, const CompressedBuf* buf);
};

template<>
struct CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>::TypeIn> {
  typedef Tally2x2<MetricFormat::SINGLE>::TypeIn T;
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
  static T elt_const(size_t ind_entry, const CompressedBuf* buf);
};

//-----------------------------------------------------------------------------
/// \brief Mirrored cpu/accel buffer that allows for compression.

class CompressedBuf {

  enum {NUM_VALUES_PER_METRIC = 4};

  enum {METRIC_FORMAT = MetricFormat::SINGLE};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef MFT::TypeIn MFTTypeIn;
  typedef Tally2x2<MetricFormat::SINGLE> TTable_t;

  typedef void Workspace_t;

  typedef size_t Lengths_t;

  // Helper struct to ensure correct order of operations.

  struct State {
    enum {IDLE = 0,
          COMPRESSED = 1,
          TRANSFER_STARTED = 2};
  };

  // Helper struct to keep track of sequential read state of buffer.

  struct Reader {

    bool is_reading_started;
    size_t ind_run_recent;
    size_t ind_in_run_recent;
    size_t ind_entry_recent;
    size_t lengths_running_sum;
    size_t ind_typein_recent;
    size_t ind_recent;
    size_t ind0_recent;
    size_t ind1_recent;
    int iE_recent;
    int jE_recent;

    void init() {
      is_reading_started = false;;
      ind_run_recent = 0; // where to look next
      ind_in_run_recent = 0; // where to look next
      ind_entry_recent = 0;
      lengths_running_sum = 0;
      ind_typein_recent = 0;
      ind_recent = 0;
      ind0_recent = 0;
      ind1_recent = 0;
      iE_recent = 0;
      jE_recent = 0;
    }
  };

public:

  CompressedBuf(MirroredBuf& buf, CEnv& env);

  void compress();

  void from_accel_start();
  void from_accel_wait();

  static bool can_compress(CEnv& env) {return
    env.threshold_tc() &&
    env.num_way() == NumWay::_3 && // TODO: implement for 2-way
    env.is_compute_method_gpu() &&
    !env.do_reduce();}
//  static bool can_compress(CEnv& env) {return false;}

  void attach(MirroredBuf& buf);

  void lock_h();
  void unlock_h();

  // For accessing entries.

  template<typename T>
  T elt_const(size_t ind0, size_t ind1) const {
    return CompressedBufAccessor_<T>::elt_const(ind0, ind1, this);
  }

  template<typename T>
  T elt_const(size_t ind_entry) const {
    return CompressedBufAccessor_<T>::elt_const(ind_entry, this);
  }

  void elt_read_start() {reader_.init();}
  size_t ind0_recent() const {return reader_.ind0_recent;}
  size_t ind1_recent() const {return reader_.ind1_recent;}
  int iE_recent() const {return reader_.iE_recent;}
  int jE_recent() const {return reader_.jE_recent;}

private:

  CEnv& env_;

  MirroredBuf* buf_;

  MirroredBuf num_nonzeros_buf_;
  MirroredBuf num_runs_buf_;
  MirroredBuf keys_buf_;
  MirroredBuf lengths_buf_;

  MirroredBuf keys_alias_buf_;
  MirroredBuf lengths_alias_buf_;

  MirroredBuf reduce_workspace_buf_;
  MirroredBuf rle_workspace_buf_;

  double compression_factor_required_() const {return 1.;}

  const size_t buf_length_max_;
  size_t num_nonzeros_approx_;
  bool do_compress_;
  size_t num_runs_;

  int state_;

  size_t num_entries_;

  //----------

  mutable Reader reader_;

  void compute_num_nonzeros_();

  bool can_compress_() const {return can_compress(env_);}

  /// \brief Compute current len of (uncompressed) buffer, in MFTTypeIn vals.

  size_t buf_length_() const {
    // NOTE: is dynamic since buf_ can change size if elsewhere realloc/realias.
    const size_t length = buf_->dim0 * buf_->dim1 * NUM_VALUES_PER_METRIC;
    COMET_INSIST(length <= buf_length_max_);
    return length;
 }

  template<typename> friend struct CompressedBufAccessorUnspecialized_;
  template<typename> friend struct CompressedBufAccessor_;

  // Disallowed methods.
  CompressedBuf(const CompressedBuf&);
  void operator=(const CompressedBuf&);

public:

  size_t num_entries() const {return num_entries_;}

  /// \brief Internal class to support CUB reduction calls.

  struct ReductionOp {
    enum {INITIAL_VALUE = -1};

    // Native input values are nonneg; let counts be neg to differentiate.
    template <typename T> __device__ __forceinline__
    bool is_value_a_count_value(const T &a) const {return a < -.5;}
    template <typename T> __device__ __forceinline__
    bool is_zero(const T &a) const {return 0 == a;}

    template <typename T> __device__ __forceinline__
    T operator()(const T &a, const T &b) const {
        return (is_value_a_count_value(a) ? a : is_zero(a) ? 0 : -1) +
               (is_value_a_count_value(b) ? b : is_zero(b) ? 0 : -1);
    }
  };

}; // CompressedBuf

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
// Implementation include files.

#include "compressed_buf.i.hh"

#endif // _COMET_COMPRESSED_BUF_HH_

//-----------------------------------------------------------------------------
