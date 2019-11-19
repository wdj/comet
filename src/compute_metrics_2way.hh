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

class GMComputeMetrics2Way {

  enum {NUM_BUF = 2};


public:

//  GMComputeMetrics2Way(GMDecompMgr& dm, GMEnv& env);
//  ~GMComputeMetrics2Way();


  GMVectorSums vector_sums_onproc;
  GMVectorSums vector_sums_offproc;
  GMVectors vectors_01[NUM_BUF];
  GMMirroredBuf metrics_buf_01[NUM_BUF];
  GMMirroredBuf vectors_buf;
  GMMirroredBuf metrics_tmp_buf;


//  GMComputeMetrics2Way(GMDecompMgr& dm, GMEnv& env) {
//    create(&dm, &env);
//  }



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

  void lock(bool& lock_val) {
    COMET_INSIST(! lock_val);
    lock_val = true;
  };

  void unlock(bool& lock_val) {
    COMET_INSIST(lock_val);
    lock_val = false;
  };

  // Disallowed methods.

//  GMComputeMetrics2Way(const GMComputeMetrics2Way&);
//  void operator=(const GMComputeMetrics2Way&);

};

//=============================================================================


//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_2way_hh_

//-----------------------------------------------------------------------------
