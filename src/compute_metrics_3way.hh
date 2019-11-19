//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_3way.hh
 * \author Wayne Joubert
 * \date   Thu Jan 21 19:07:47 EST 2016
 * \brief  Calculate metrics, 3-way.
 * \note   Copyright (C) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_compute_metrics_3way_hh_
#define _comet_compute_metrics_3way_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"
#include "metrics.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

class ComputeMetrics3Way {

public:

  ComputeMetrics3Way(GMDecompMgr& dm, GMEnv& env);
  ~ComputeMetrics3Way();

  void compute(GMMetrics& metrics, GMVectors& vectors);

private:

  GMEnv& env_;

  void compute_notall2all_(GMMetrics& metrics, GMVectors& vectors);
  void compute_all2all_(GMMetrics& metrics, GMVectors& vectors);

  // Disallowed methods.

  ComputeMetrics3Way(const ComputeMetrics3Way&);
  void operator=(const ComputeMetrics3Way&);

}; // ComputeMetrics3Way

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_compute_metrics_3way_hh_

//-----------------------------------------------------------------------------
