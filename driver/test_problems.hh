/*---------------------------------------------------------------------------*/
/*!
 * \file   test_problems.hh
 * \author Wayne Joubert
 * \date   Mon Aug  7 17:02:51 EDT 2017
 * \brief  Generator for synthetic test problems, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#ifndef _gm_test_problems_hh_
#define _gm_test_problems_hh_

#include "env.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "driver.hh"

/*===========================================================================*/

void set_vectors_random(GMVectors* vectors, DriverOptions* do_, GMEnv* env);

void set_vectors_analytic(GMVectors* vectors, DriverOptions* do_, GMEnv* env);

void check_metrics(GMMetrics* metrics, DriverOptions* do_, GMEnv* env);

/*===========================================================================*/

#endif /*---_gm_test_problems_hh_---*/

/*---------------------------------------------------------------------------*/
