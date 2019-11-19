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

  ComputeMetrics2Way(GMDecompMgr& dm, GMEnv& env)
    : env_(env)
    , vector_sums_onproc{0}
    , vector_sums_offproc{0}
    , vectors_01{0}
    , metrics_buf_01{0}
    , vectors_buf{0}
    , metrics_tmp_buf{0}
  {
    create(&dm, &env);
  }

  ~ComputeMetrics2Way() {
    destroy(&env_);
  }



  void create(
    GMDecompMgr* dm,
    GMEnv* env);
  
  void destroy(
    GMEnv* env);

  void compute(
    GMMetrics* metrics,
    GMVectors* vectors,
    GMEnv* env) {
    if (!env->all2all()) {
      compute_notall2all(metrics, vectors, env);
    } else {
      compute_all2all(metrics, vectors, env);
    }
  }

  void compute_notall2all(
    GMMetrics* metrics,
    GMVectors* vectors,
    GMEnv* env);

  void compute_all2all(
    GMMetrics* metrics,
    GMVectors* vectors,
    GMEnv* env);


private:

  Env& env_;

  GMVectorSums vector_sums_onproc;
  GMVectorSums vector_sums_offproc;
  GMVectors vectors_01[NUM_BUF];
  GMMirroredBuf metrics_buf_01[NUM_BUF];
  GMMirroredBuf vectors_buf;
  GMMirroredBuf metrics_tmp_buf;


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


//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_2way_hh_

//-----------------------------------------------------------------------------
