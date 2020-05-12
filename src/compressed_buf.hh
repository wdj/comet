//-----------------------------------------------------------------------------
/*!
 * \file   compressed_buf.hh
 * \author Wayne Joubert
 * \date   Tue May  5 20:11:12 EDT 2020
 * \brief  Mirrored buffer allowing for compression.
 * \note   Copyright (C) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compressed_buf_hh_
#define _comet_compressed_buf_hh_

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class CompressedBuf {

  enum {NUM_VALUES_PER_METRIC = 4};

  enum {METRIC_FORMAT = MetricFormat::SINGLE};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef MFT::TypeIn MFTTypeIn;

  typedef void Workspace_t;

  struct State {
    enum {IDLE = 0,
          COMPRESSED = 1,
          TRANSFER_STARTED = 2};
  };

public:

  CompressedBuf(MirroredBuf& buf, CEnv& env);

  void compress();

  void from_accel_start();
  void from_accel_wait();

  template<typename T>
  T elt_const(size_t i0, size_t i1) const {
    if (do_compress_) {
      COMET_ASSERT(State::IDLE == state_);
      COMET_ASSERT(i0 < buf_->dim0 && i1 < buf_->dim1);
      const size_t i = i0 + buf_->dim0 * i1;
      COMET_ASSERT(is_reading_started_ ? i > i_prev_ : 0 == i);

printf("%i %i %i\n", (int)i, (int)i0, (int)i1); // FIX

      is_reading_started_ = true;
      i_prev_ = i;


#if 0

temporary 2x2 table; initialize to zero
static assertion to check T is correct type
loop over 4 table elements in proper order
  map from i and table index to TypeIn index
  march down rle from current location to get at or past sought entry, or to end
  manage run lengths vs. start pointers (scan-sums)
  is NOT an error if we skip something in the rle
  if found nonzero, insert

#endif


      // TODO
      // FIX
      return buf_->elt_const<T>(i0, i1);
      // FIX


    } else {
      return buf_->elt_const<T>(i0, i1);
    }
  }

  static bool can_compress(CEnv& env) {return env.threshold_tc() &&
                                              env.is_compute_method_gpu() &&
                                              !env.do_reduce();}
//FIX
//  static bool can_compress(CEnv& env) {return false;}

  void attach(MirroredBuf& buf);

  void lock_h();
  void unlock_h();

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

  double compress_threshold_() const {return .2;}

  const size_t length_max_;
  size_t num_nonzeros_approx_;
  bool do_compress_;
  size_t num_runs_;

  int state_;

  size_t mutable i_prev_;
  bool mutable is_reading_started_;
  size_t mutable ind_rle_;

  //bool is_open_;
  //size_t read_ptr_;

  void compute_num_nonzeros_();

  bool can_compress_() const {return can_compress(env_);}

  /// \brief Compute current len of (uncompressed) buffer, in MFTTypeIn vals.

  size_t length_() const {
    // NOTE: is dynamic since buf_ can change size if elsewhere realloc/realias.
    const size_t length = buf_->dim0 * buf_->dim1 * NUM_VALUES_PER_METRIC;
    COMET_INSIST(length <= length_max_);
    return length;
 }

  // Disallowed methods.
  CompressedBuf(const CompressedBuf&);
  void operator=(const CompressedBuf&);

public:

  // Internal class to support CUB reduction calls.

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

#endif // _comet_compressed_buf_hh_

//-----------------------------------------------------------------------------
