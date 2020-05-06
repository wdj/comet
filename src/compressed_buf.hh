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

public:

  CompressedBuf(MirroredBuf& buf, CEnv& env);



  struct ReductionOp {
    template <typename T> __device__ __forceinline__
    // Native input values are nonneg; let counts be neg to differentiate.
    bool is_a_count_value(const T &a) const {return a < -.5;}

    template <typename T> __device__ __forceinline__
    T operator()(const T &a, const T &b) const {
        return (is_a_count_value(a) ? a : -1) + 
               (is_a_count_value(b) ? b : -1);
    }
  };

private:

  CEnv& env_;
  MirroredBuf& buf_;

  MirroredBuf num_nonzeros_buf_;
  MirroredBuf num_runs_buf_;
  MirroredBuf keys_buf_;
  MirroredBuf lengths_buf_;

  MirroredBuf keys_alias_buf_;
  MirroredBuf lengths_alias_buf_;

  MirroredBuf reduce_workspace_buf_;
  MirroredBuf rle_workspace_buf_;

  size_t length_max_;

  bool is_open_;
  size_t read_ptr_;

  void compute_nonzeros_();

  bool try_compress() const {return env_.threshold_tc() &&
                                    env_.is_compute_method_gpu();}





  // Disallowed methods.
  CompressedBuf(const CompressedBuf&);
  void operator=(const CompressedBuf&);
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compressed_buf_hh_

//-----------------------------------------------------------------------------
