//-----------------------------------------------------------------------------
/*!
 * \file   test_problems.hh
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_test_problems_hh_
#define _comet_test_problems_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"

#include "driver.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

void set_vectors_synthetic(GMVectors* vectors, int problem_type, int verbosity,
                           GMEnv* env);

static int problem_type_default() {return GM_PROBLEM_TYPE_ANALYTIC;}
//static int problem_type_default() {return GM_PROBLEM_TYPE_RANDOM;}

void check_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env);

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_test_problems_hh_

//-----------------------------------------------------------------------------
