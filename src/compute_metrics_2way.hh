//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way.hh
 * \author Wayne Joubert
 * \date   Thu Jan  7 10:21:09 EST 2016
 * \brief  Calculate metrics, 2-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_2way_hh_
#define _comet_compute_metrics_2way_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class ComputeMetrics2Way {

  enum {NUM_BUF = 2};

public:

  ComputeMetrics2Way(GMDecompMgr& dm, GMEnv& env);
  ~ComputeMetrics2Way();

  void compute(GMMetrics& metrics, GMVectors& vectors);

private:

  Env& env_;

  GMVectors vectors_01_[NUM_BUF];
  MirroredBuf metrics_buf_01_[NUM_BUF];
  MirroredBuf vectors_buf_;
  MirroredBuf metrics_tmp_buf_;

  void compute_notall2all_(GMMetrics& metrics, GMVectors& vectors);
  void compute_all2all_(GMMetrics& metrics, GMVectors& vectors);

  void lock(bool& lock_val) {
    COMET_INSIST(! lock_val);
    lock_val = true;
  };

  void unlock(bool& lock_val) {
    COMET_INSIST(lock_val);
    lock_val = false;
  };

  // Disallowed methods.

  ComputeMetrics2Way(const ComputeMetrics2Way&);
  void operator=(const ComputeMetrics2Way&);

};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_2way_hh_

//-----------------------------------------------------------------------------
