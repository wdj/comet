/*---------------------------------------------------------------------------*/
/*!
 * \file   checksums.hh
 * \author Wayne Joubert
 * \date   Mon Aug  7 14:47:01 EDT 2017
 * \brief  Checksums for metrics, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================
=============================================================================*/

#ifndef _gm_checksums_hh_
#define _gm_checksums_hh_

#include <float.h>

#include <env.hh>
#include <metrics.hh>

/*===========================================================================*/
/*---Checksums---*/

/*---Multiprecision integers---*/

enum { GM_MULTIPREC_INT_SIZE = 16 };

typedef struct {
  size_t data[GM_MULTIPREC_INT_SIZE];
} GMMultiprecInt;

/*---Struct with checksum info---*/

enum { GM_CHECKSUM_SIZE = 3 };

typedef struct {
  size_t data[GM_CHECKSUM_SIZE];
  bool is_overflowed;
  double value_max;
  GMMultiprecInt sum;
  double sum_d;
  bool is_started;
  bool computing_checksum;
} GMChecksum;

/*---------------------------------------------------------------------------*/

static GMChecksum GMChecksum_null() {
  GMChecksum result;
  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result.data[i] = 0;
  }
  result.is_overflowed = false;
  result.value_max = -DBL_MAX;
  GMMultiprecInt sum = {0};
  result.sum = sum;
  result.sum_d = 0;
  result.computing_checksum = true;
  return result;
}

/*---------------------------------------------------------------------------*/

static bool gm_are_checksums_equal(GMChecksum c1, GMChecksum c2) {
  if ((!c1.computing_checksum) || (!c2.computing_checksum)) {
    return true;
  }
  bool result = true;
  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result = result && c1.data[i] == c2.data[i];
  }
  result = result && c1.is_overflowed == c2.is_overflowed;
  result = result && c1.value_max == c2.value_max;
  return result;
}


/*===========================================================================*/
/*---Metrics checksum---*/

void GMMetrics_checksum(GMMetrics* metrics, GMChecksum* cs, GMEnv* env);

/*===========================================================================*/

#endif /*---_gm_checksums_hh_---*/

/*---------------------------------------------------------------------------*/
