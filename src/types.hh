//-----------------------------------------------------------------------------
/*!
 * \file   types.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Fundamental scalar types for algorithms; associated functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_types_hh_
#define _comet_types_hh_

#include "cstdint"
#include "cfloat"
#include <limits>
//#include <cstdio> // FIX

#include "assertions.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Type ids

enum {
  GM_DATA_TYPE_FLOAT = 1,
  GM_DATA_TYPE_BITS1 = 2, // Not implemented
  GM_DATA_TYPE_UINT64 = 3,  //---(design of this entry is not complete)
  GM_DATA_TYPE_BITS2 = 4,
  GM_DATA_TYPE_TALLY2X2 = 5,
  GM_DATA_TYPE_TALLY4X2 = 6
};

//=============================================================================
// Types: general

// TODO: use C++ types that have guaranteed sizes, e.g., int32_t.

// Integer types

//typedef unsigned char GMUInt8;
//typedef unsigned short int GMUInt16;
//typedef unsigned int GMUInt32;

//typedef signed long long int GMInt64;
//typedef unsigned long long int GMUInt64;

#ifdef COMET_USE_INT128
typedef unsigned __int128 GMUInt128;
#endif

// Floating point of explicit (double) precision

typedef float GMFp32;
typedef double GMFp64;

static void gm_check_type_sizes() {
  COMET_STATIC_ASSERT(sizeof(GMFp32) == 32/8);
  COMET_STATIC_ASSERT(sizeof(GMFp64) == 64/8);
  COMET_STATIC_ASSERT(sizeof(int) == 4);
  COMET_STATIC_ASSERT(sizeof(size_t) == 8);
//  COMET_STATIC_ASSERT(sizeof(GMUInt8) == 1);
//  COMET_STATIC_ASSERT(sizeof(GMUInt16) == 2);
//  COMET_STATIC_ASSERT(sizeof(GMUInt32) == 4);
//  COMET_STATIC_ASSERT(sizeof(GMInt64) == 8);
//  COMET_STATIC_ASSERT(sizeof(GMUInt64) == 8);
#ifdef COMET_USE_INT128
  COMET_STATIC_ASSERT(sizeof(GMUInt128) == 16);
#endif
}

//=============================================================================
// Types (mainly) for Czekanowski metric

//---TODO: revise nomenclature to be different from "GMFloat2" ...

#ifdef COMET_FP_PRECISION_SINGLE
  typedef float GMFloat;
  #define COMET_MPI_FLOAT MPI_FLOAT
  enum { FP_PRECISION_DOUBLE = false };
#ifdef COMET_FP_PRECISION_DOUBLE
#error Cannot set both COMET_FP_PRECISION_SINGLE and COMET_FP_PRECISION_DOUBLE.
#endif
#else
#ifdef COMET_FP_PRECISION_DOUBLE
  typedef double GMFloat;
  #define COMET_MPI_FLOAT MPI_DOUBLE
  enum { FP_PRECISION_DOUBLE = true };
#else
#error Must set COMET_FP_PRECISION_SINGLE or COMET_FP_PRECISION_DOUBLE.
#endif
#endif

//=============================================================================
// Types for CCC and DUO metrics

// For Vectors: single 2-bit value (seminibble):
// use unsigned int as a container for a single item
typedef unsigned int GMBits2;

// For Vectors: packed: 2 long integers, used to store 64 seminibbles
typedef unsigned long long int GMBits1_2x64;
typedef struct { GMBits1_2x64 data[2]; } GMBits2x64;

// For Vectors: largest allowed size of a data value
enum { GM_BITS2_MAX_VALUE_BITS = 2 };

// For metrics: single integer to store a tally result
typedef unsigned int GMTally1;

// For Metrics: 2 (4) doubles to store 4 (8) packed tally results:
// use 25 bits of each 52-bit mantissa to store a result
typedef struct { GMFp64 data[2]; } GMTally2x2;
typedef struct { GMFp64 data[4]; } GMTally4x2;

// For Metrics: largest allowed size of a data value
enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

// For Metrics: for packing of multipliers
typedef GMFp64 GMFloat2;
typedef struct { GMFp64 data[2]; } GMFloat3;

// Marker value for a missing or unknown 2-bit entry for sparse case

enum { GM_2BIT_UNKNOWN = 2 * 1 + 1 * 0 };

//=============================================================================
// Types for CCC and DUO metrics: functions

//-----------------------------------------------------------------------------
// Return null value; also use static asserts to check sizes

static GMBits2x64 GMBits2x64_null() {
  COMET_STATIC_ASSERT(2 * GM_TALLY1_MAX_VALUE_BITS <= DBL_MANT_DIG);
  COMET_STATIC_ASSERT(sizeof(GMBits2) * 8 >= GM_BITS2_MAX_VALUE_BITS);
  COMET_STATIC_ASSERT(sizeof(GMBits1_2x64) == 8);
  COMET_STATIC_ASSERT(sizeof(GMBits2x64) == 2 * sizeof(GMBits1_2x64));
  COMET_STATIC_ASSERT(sizeof(GMBits2x64) == 16);
  COMET_STATIC_ASSERT(sizeof(GMTally2x2) == sizeof(GMBits2x64)); // for Magma

  GMBits2x64 value;
  value.data[0] = 0;
  value.data[1] = 0;
  return value;
}

//-----------------------------------------------------------------------------
// Return null value; also use static asserts to check sizes

static GMTally2x2 GMTally2x2_null() {
  COMET_STATIC_ASSERT(sizeof(GMTally1) * 8 >= GM_TALLY1_MAX_VALUE_BITS);
  COMET_STATIC_ASSERT(sizeof(GMTally2x2) == 16);
  COMET_STATIC_ASSERT(sizeof(GMTally2x2) == sizeof(GMBits2x64)); // for Magma

  GMTally2x2 value;
  value.data[0] = 0;
  value.data[1] = 0;
  return value;
}

//-----

static GMTally4x2 GMTally4x2_null() {
  COMET_STATIC_ASSERT(sizeof(GMTally4x2) == 32);

  GMTally4x2 value;
  value.data[0] = 0;
  value.data[1] = 0;
  value.data[2] = 0;
  value.data[3] = 0;
  return value;
}

//-----------------------------------------------------------------------------
// Encode/decode between float and pair of tally values

__host__ __device__
static void GMTally1_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMFp64 v) {
  COMET_ASSERT(val0);
  COMET_ASSERT(val1);
  const uint64_t tally2 = (uint64_t)v;
  COMET_ASSERT(v == (GMFp64)tally2);
  const GMTally1 v0 =
      tally2 & ((((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) - 1);
  const GMTally1 v1 = tally2 >> GM_TALLY1_MAX_VALUE_BITS;
  *val0 = v0;
  *val1 = v1;
  COMET_ASSERT(v ==
           (GMFp64)(v0 + (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) * v1));
  //COMET_ASSERT(v0 >= 0);
  //COMET_ASSERT(v1 >= 0);
  COMET_ASSERT(v0 < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
  COMET_ASSERT(v1 < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
}

//----------

__host__ __device__
static void GMTally1_decode(GMTally1& val0,
                            GMTally1& val1,
                            GMFp64 v) {
  const uint64_t tally2 = (uint64_t)v;
  COMET_ASSERT(v == (GMFp64)tally2);
  const GMTally1 v0 =
      tally2 & ((((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) - 1);
  const GMTally1 v1 = tally2 >> GM_TALLY1_MAX_VALUE_BITS;
  val0 = v0;
  val1 = v1;
  COMET_ASSERT(v ==
           (GMFp64)(v0 + (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) * v1));
  //COMET_ASSERT(v0 >= 0);
  //COMET_ASSERT(v1 >= 0);
  COMET_ASSERT(v0 < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
  COMET_ASSERT(v1 < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
}

//----------

__host__ __device__
static GMFp64 GMTally1_encode(GMTally1 val0, GMTally1 val1) {
  const uint64_t tally2 =
      val0 + (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) * val1;
  const GMFp64 result = (GMFp64)tally2;
  COMET_ASSERT(val0 == (((uint64_t)result) &
                    ((((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS) - 1)));
  COMET_ASSERT(val1 == ((uint64_t)result) >> GM_TALLY1_MAX_VALUE_BITS);
  return result;
}

//-----------------------------------------------------------------------------
// Encode for multipliers/sums

static GMFloat2 GMFloat2_encode(GMTally1 val0, GMTally1 val1) {
  return GMTally1_encode(val0, val1);
}

//----------

static GMFloat3 GMFloat3_encode(GMTally1 val0, GMTally1 val1, GMTally1 val2) {
  GMFloat3 result; // here we should set = null to be super cautious
  const GMTally1 dummy = 0;
  result.data[0] = GMTally1_encode(val0, val1);
  result.data[1] = GMTally1_encode(val2, dummy);
  return result;
}

//-----------------------------------------------------------------------------
// Decode for multipliers/sums

static void GMFloat2_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMFloat2 v) {
  GMTally1_decode(val0, val1, v);
}

//----------

static void GMFloat3_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMTally1* __restrict__ val2,
                            GMFloat3 v) {
  GMTally1_decode(val0, val1, v.data[0]);
  GMTally1 dummy;
  GMTally1_decode(val2, &dummy, v.data[1]);
}

//-----------------------------------------------------------------------------
// Get an entry: 2x2

static GMTally1 GMTally2x2_get(GMTally2x2 tally2x2, int i0, int i1) {
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  const uint64_t tally2 = tally2x2.data[i0];

  const GMTally1 result =
      i1 == 0 ? tally2 % (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
              : tally2 / (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
  //COMET_ASSERT(result >= 0);
  COMET_ASSERT(result < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
  return result;
}

//-----------------------------------------------------------------------------
// Get an entry: 4x2

static GMTally1 GMTally4x2_get(GMTally4x2 tally4x2, int i0, int i1, int i2) {
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  const uint64_t tally2 = tally4x2.data[i1 + 2 * i0];

  const GMTally1 result =
      i2 == 0 ? tally2 % (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
              : tally2 / (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
  //COMET_ASSERT(result >= 0);
  COMET_ASSERT(result < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
  return result;
}

//=============================================================================

template<typename T> int mantissa_digits() {
  COMET_STATIC_ASSERT(std::numeric_limits<T>::radix == 2);
  return std::numeric_limits<T>::digits;
}

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_types_hh_

//-----------------------------------------------------------------------------
