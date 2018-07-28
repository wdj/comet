//-----------------------------------------------------------------------------
/*!
 * \file   checksums.hh
 * \author Wayne Joubert
 * \date   Mon Aug  7 14:47:01 EDT 2017
 * \brief  Checksums for metrics, header.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_checksums_hh_
#define _gm_checksums_hh_

#include "env.hh"
#include "metrics.hh"

//=============================================================================
// Multiprecision integers

enum { GM_MULTIPREC_INT_SIZE = 16 };

typedef struct {
  size_t data[GM_MULTIPREC_INT_SIZE];
} GMMultiprecInt;

//-----------------------------------------------------------------------------

// Checksum struct

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

//-----------------------------------------------------------------------------
// Set to null

GMChecksum GMChecksum_null();

//-----------------------------------------------------------------------------
// Check whether two checksums equal

bool GMChecksum_equal(GMChecksum* cksum1, GMChecksum* cksum2);

//-----------------------------------------------------------------------------
// Compute checksum of metrics object

void GMChecksum_metrics(GMChecksum* cksum, GMChecksum* cksum_local,
                        GMMetrics* metrics, GMEnv* env);

//-----------------------------------------------------------------------------
// Print ckecksum to stdout.

void GMChecksum_print(GMChecksum* cksum, GMEnv* env);

//=============================================================================

#endif // _gm_checksums_hh_

//-----------------------------------------------------------------------------
