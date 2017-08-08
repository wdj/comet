//-----------------------------------------------------------------------------
/*!
 * \file   input_output.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  I/O functions used by driver, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_input_output_hh_
#define _gm_input_output_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "driver.hh"

//=============================================================================

void set_vectors_from_file(GMVectors* vectors, DriverOptions* do_, GMEnv* env);

void output_metrics_file(GMMetrics* metrics, DriverOptions* do_,
                         FILE* file, double threshold, GMEnv* env);

void output_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env);

//=============================================================================

#endif /*---_gm_input_output_hh_---*/

//-----------------------------------------------------------------------------
