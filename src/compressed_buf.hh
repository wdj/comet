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

#include <type_traits>

#include "env.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Helper function for cartesian access to buffer elements.

// Forward declaration.
class CompressedBuf;

template<typename T>
struct CompressedBufAccessor_ {
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
};

template<>
struct CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>> {
  typedef Tally2x2<MetricFormat::SINGLE> T;
  static T elt_const(size_t ind0, size_t ind1, const CompressedBuf* buf);
};

//-----------------------------------------------------------------------------
/// \brief Mirrored cpu/accel buffer that allows for compression.

class CompressedBuf {

  enum {NUM_VALUES_PER_METRIC = 4};

  enum {METRIC_FORMAT = MetricFormat::SINGLE};

  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef MFT::TypeIn MFTTypeIn;

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
    size_t ind_prev;
    size_t ind_runs;
    size_t lengths_sum;

    void init() {
      is_reading_started = false;;
      ind_prev = 0;
      ind_runs = 0;
      lengths_sum = 0;
    }
  };

public:

  CompressedBuf(MirroredBuf& buf, CEnv& env);

  void compress();

  void from_accel_start();
  void from_accel_wait();

  template<typename T>
  T elt_const(size_t ind0, size_t ind1) const {
    return CompressedBufAccessor_<T>::elt_const(ind0, ind1, this);
  }

  static bool can_compress(CEnv& env) {return
    env.threshold_tc() &&
    env.num_way() == NUM_WAY::_3 && // TODO: implement for 2-way
    env.is_compute_method_gpu() &&
    !env.do_reduce();}
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

  double compression_factor_required_() const {return 1.;}

  const size_t length_max_;
  size_t num_nonzeros_approx_;
  bool do_compress_;
  size_t num_runs_;

  int state_;

  mutable Reader reader_;

  void compute_num_nonzeros_();

  bool can_compress_() const {return can_compress(env_);}

  /// \brief Compute current len of (uncompressed) buffer, in MFTTypeIn vals.

  size_t length_() const {
    // NOTE: is dynamic since buf_ can change size if elsewhere realloc/realias.
    const size_t length = buf_->dim0 * buf_->dim1 * NUM_VALUES_PER_METRIC;
    COMET_INSIST(length <= length_max_);
    return length;
 }

  template<typename> friend class CompressedBufAccessor_;

  // Disallowed methods.
  CompressedBuf(const CompressedBuf&);
  void operator=(const CompressedBuf&);

public:

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

//-----------------------------------------------------------------------------
// Inlinable implementations.

template<typename T> T CompressedBufAccessor_<T>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  return cbuf->buf_->elt_const<T>(ind0, ind1);
}

inline Tally2x2<MetricFormat::SINGLE>
CompressedBufAccessor_<Tally2x2<MetricFormat::SINGLE>>::
elt_const(size_t ind0, size_t ind1, const CompressedBuf* cbuf) {
  typedef Tally2x2<MetricFormat::SINGLE> T;
  typedef CompressedBuf CBuf;
  typedef CBuf::MFTTypeIn MFTTypeIn;
  const MirroredBuf* buf_ = cbuf->buf_;
  CompressedBuf::Reader& reader_ = cbuf->reader_;

  if (cbuf->do_compress_) {

    COMET_ASSERT(CBuf::State::IDLE == cbuf->state_);
    COMET_ASSERT(ind0 < buf_->dim0 && ind1 < buf_->dim1);
    const size_t ind = ind0 + buf_->dim0 * ind1;
    COMET_ASSERT(ind > reader_.ind_prev || ! reader_.is_reading_started);

    reader_.is_reading_started = true;
    reader_.ind_prev = ind;

    T result = T::null();

#ifdef COMET_ASSERTIONS_ON
    // number of TypeIn (float) entries of (uncompressed) buf.
    const size_t dim_typein = 4 * buf_->dim0 * buf_->dim1;
#endif

    const size_t dim_runs = cbuf->num_runs_;
    size_t& ind_runs = reader_.ind_runs;

    // Loop to look for, pick up 4 table entries.

    // NOTE: order of next 2 nested loops must match mem layout of Tally2x2,
    // since that is the order of elts submitted to the rle.

    for (int i0=0; i0<2; ++i0) {
      for (int i1=0; i1<2; ++i1) {

        // index for typein value in (uncompressed) array.
        const size_t ind_typein = i1 + 2 * ( i0 + 2 * ind );

        COMET_ASSERT(ind_typein >= reader_.lengths_sum);
        COMET_ASSERT(ind_typein < dim_typein);

        // Loop over rle, starting at current location, to seek element.

        for( ; ind_runs < dim_runs; ) {

            const size_t length_run = cbuf->lengths_alias_buf_.
              elt_const<CompressedBuf::Lengths_t>(ind_runs, 0);
            const size_t ind_typein_min = reader_.lengths_sum;
            const size_t ind_typein_max = reader_.lengths_sum + length_run;

            if (ind_typein >= ind_typein_min &&
                ind_typein < ind_typein_max) {
              const MFTTypeIn key_run =
                cbuf->keys_alias_buf_.elt_const<MFTTypeIn>(ind_runs, 0);
              T::set(result, i0, i1, key_run);
              break;
            }

            ++ind_runs;
            reader_.lengths_sum += length_run;

        } // for

      } // for i1
    } // for i0

    //COMET_ASSERT(buf_->elt_const<T>(ind0, ind1)==result);
    //return buf_->elt_const<T>(ind0, ind1);
    return result;

  } else {
    return buf_->elt_const<T>(ind0, ind1);
  }
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compressed_buf_hh_

//-----------------------------------------------------------------------------
