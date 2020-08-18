//-----------------------------------------------------------------------------
/*!
 * \file   types.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Fundamental scalar types for algorithms; associated functions.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef _COMET_TYPES_HH_
#define _COMET_TYPES_HH_

#include "cstdint"
#include "cfloat"
#include <limits>

#include "assertions.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Type ids

enum {
  GM_DATA_TYPE_FLOAT = 1,
  GM_DATA_TYPE_BITS1 = 2, // Not implemented
  GM_DATA_TYPE_UINT64 = 3,  //---(design of this selection is not complete)
  GM_DATA_TYPE_BITS2 = 4,
  GM_DATA_TYPE_TALLY2X2 = 5,
  GM_DATA_TYPE_TALLY4X2 = 6
};

//=============================================================================
/// \brief Basic types.

// TODO: consistently use C++ types that have guaranteed sizes, e.g., int32_t.

struct BasicTypes {

  typedef float FP32;
  typedef double FP64;

# ifdef COMET_USE_INT128
    typedef unsigned __int128 BigUInt;
# else
    typedef uint64_t BigUInt;
# endif

private:

  enum {BITS_PER_BYTE_ = 8};

  static void check_type_sizes_() {
    COMET_STATIC_ASSERT(sizeof(FP32) == 32/BITS_PER_BYTE_);
    COMET_STATIC_ASSERT(sizeof(FP64) == 64/BITS_PER_BYTE_);
    COMET_STATIC_ASSERT(sizeof(int) == 4);
    COMET_STATIC_ASSERT(sizeof(size_t) == 8);
#   ifdef COMET_USE_INT128
      COMET_STATIC_ASSERT(sizeof(BigUInt) == 128/BITS_PER_BYTE_);
#   endif
  }

  // Disallowed methods.

  BasicTypes(const BasicTypes&);
  void operator=(const BasicTypes&);
}; // BasicTypes

//=============================================================================
// Type for storing coordinates of a metrics item.

typedef size_t MetricItemCoords_t;

//=============================================================================
// Types (mainly) for Czekanowski metric

//---TODO: revise nomenclature to be different from "GMFloat2" ...

#ifdef COMET_FP_PRECISION_DOUBLE
  typedef double GMFloat;
  #define COMET_MPI_FLOAT MPI_DOUBLE
#else
  typedef float GMFloat;
  #define COMET_MPI_FLOAT MPI_FLOAT
#endif

//-----------------------------------------------------------------------------
/// \brief Helper class for metric element storage format.

struct MetricFormat {
  enum {PACKED_DOUBLE = 0,
        SINGLE = 1};
};

//=============================================================================
// Types for CCC and DUO metrics

// For Metrics: largest allowed size of a data value
enum { GM_TALLY1_MAX_VALUE_BITS = 26 };

// For Vectors: single 2-bit value (seminibble):
// use unsigned int as a container for a single item
typedef unsigned int GMBits2;

// For Vectors: packed: 2 long integers, used to store 64 seminibbles
typedef unsigned long long int GMBits1_2x64;
typedef struct { GMBits1_2x64 data[2]; } GMBits2x64;

// For Vectors: largest allowed size of a data value
enum { GM_BITS2_MAX_VALUE_BITS = 2 };

// For Metrics: single integer to store a tally result
typedef unsigned int GMTally1;

// For Metrics: double used to store two metric numerator values.
typedef BasicTypes::FP64 PackedDouble;

// For Metrics: two floats used to store two metric numerator values.

typedef struct { BasicTypes::FP32 data[2]; } Single2;

// For Metrics: 2 (4) doubles to store 4 (8) packed tally results:
// use 25 bits of each 52-bit mantissa to store a result
typedef struct { PackedDouble data[2]; } GMTally2x2;
typedef struct { PackedDouble data[4]; } GMTally4x2;

// For Metrics: for packing of multipliers
typedef PackedDouble GMFloat2;
typedef struct { PackedDouble data[2]; } GMFloat3;

// Marker value for a missing or unknown 2-bit vector entry for sparse case

enum { GM_2BIT_UNKNOWN = 2 * 1 + 1 * 0 };

//=============================================================================
// Templatized types for CCC and DUO metrics

template<int METRIC_FORMAT> struct MetricFormatTraits;

//-----------------------------------------------------------------------------
/// \brief Metric format traits for packed double case.

template<> struct MetricFormatTraits<MetricFormat::PACKED_DOUBLE> {
  typedef PackedDouble Type;
  typedef GMTally1 TypeIn;

  __host__ __device__
  static void decode(TypeIn& __restrict__ val0,
                     TypeIn& __restrict__ val1,
                     const Type& v) {
    const uint64_t tally2 = (uint64_t)v;
    COMET_ASSERT(v == (Type)tally2);
    const uint64_t shifter = (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
    const GMTally1 v0 = tally2 & (shifter - 1);
    const GMTally1 v1 = tally2 >> GM_TALLY1_MAX_VALUE_BITS;
    val0 = v0;
    val1 = v1;
    COMET_ASSERT(v == (Type)(v0 + v1 * shifter));
    COMET_ASSERT(v0 < shifter);
    COMET_ASSERT(v1 < shifter);
  }

  __host__ __device__
  static void encode(Type& v,
                     const TypeIn& __restrict__ val0,
                     const TypeIn& __restrict__ val1) {
    const uint64_t shifter = (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
    const uint64_t tally2 = val0 + shifter * val1;
    v = (Type)tally2;
    COMET_ASSERT(val0 == (((uint64_t)v) & (shifter - 1)));
    COMET_ASSERT(val1 == ((uint64_t)v) >> GM_TALLY1_MAX_VALUE_BITS);
  }

  __host__ __device__
  static void add(Type& v,
                  const TypeIn& __restrict__ val0,
                  const TypeIn& __restrict__ val1) {
    const uint64_t shifter = (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
    const uint64_t tally2 = val0 + shifter * val1;
#ifdef COMET_ASSERTIONS_ON
    const Type vold = v;
#endif
    v += (Type)tally2;
    COMET_ASSERT(val0 == (((uint64_t)(v-vold)) & (shifter - 1)));
    COMET_ASSERT(val1 == ((uint64_t)(v-vold)) >> GM_TALLY1_MAX_VALUE_BITS);
  }

  __host__ __device__
  static void subtract(Type& v,
                       const TypeIn& __restrict__ val0,
                       const TypeIn& __restrict__ val1) {
    const uint64_t shifter = (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
    const uint64_t tally2 = val0 + shifter * val1;
#ifdef COMET_ASSERTIONS_ON
    const Type vold = v;
#endif
    v -= (Type)tally2;
    COMET_ASSERT(val0 == (((uint64_t)(vold-v)) & (shifter - 1)));
    COMET_ASSERT(val1 == ((uint64_t)(vold-v)) >> GM_TALLY1_MAX_VALUE_BITS);
  }

  __host__ __device__
  static Type null() {
    return 0;
  }
};

//-----------------------------------------------------------------------------
/// \brief Metric format traits for FP32 case.

template<> struct MetricFormatTraits<MetricFormat::SINGLE> {
  typedef Single2 Type;
  typedef float TypeIn;

  __host__ __device__
  static void decode(TypeIn& __restrict__ val0,
                     TypeIn& __restrict__ val1,
                     const Type& v) {
    val0 = v.data[0];
    val1 = v.data[1];
  }

  __host__ __device__
  static void encode(Type& v,
                     const TypeIn& __restrict__ val0,
                     const TypeIn& __restrict__ val1) {
    v.data[0] = val0;
    v.data[1] = val1;
  }

  __host__ __device__
  static void add(Type& v,
                  const TypeIn& __restrict__ val0,
                  const TypeIn& __restrict__ val1) {
    v.data[0] += val0;
    v.data[1] += val1;
  }

  __host__ __device__
  static void subtract(Type& v,
                       const TypeIn& __restrict__ val0,
                       const TypeIn& __restrict__ val1) {
    v.data[0] -= val0;
    v.data[1] -= val1;
  }

  __host__ __device__
  static Type null() {
    return {0, 0};
  }
};

//-----------------------------------------------------------------------------
/// \brief Tally table struct to support 2-way bitwise methods.

template<int METRIC_FORMAT> struct Tally2x2 {
  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type Type;
  typedef typename MFT::TypeIn TypeIn;
  enum {NUM = 2};
  Type data[NUM];
  typedef Tally2x2<METRIC_FORMAT> This_t;

  __host__ __device__ 
  static This_t null() {
    This_t result;
    result.data[0] = MFT::null();
    result.data[1] = MFT::null();
    return result;
  }

  __host__ __device__ 
  static TypeIn get(const This_t& value, int iE, int jE) {
    const Type data = value.data[iE];
    TypeIn results[2];
    MFT::decode(results[0], results[1], data);
    return results[jE];
  }

  __host__ __device__ 
  static void set(This_t& value, int iE, int jE, TypeIn v) {
    Type& data = value.data[iE];
    TypeIn results[2];
    MFT::decode(results[0], results[1], data);
    results[jE] = v;
    MFT::encode(data, results[0], results[1]);
  }
};

//-----------------------------------------------------------------------------
/// \brief Tally table struct to support 3-way bitwise methods.

template<int METRIC_FORMAT> struct Tally4x2 {
  typedef MetricFormatTraits<METRIC_FORMAT> MFT;
  typedef typename MFT::Type Type;
  typedef typename MFT::TypeIn TypeIn;
  enum {NUM = 4};
  Type data[NUM];
  typedef Tally4x2<METRIC_FORMAT> This_t;

  __host__ __device__ 
  static This_t null() {
    This_t result;
    result.data[0] = MFT::null();
    result.data[1] = MFT::null();
    result.data[2] = MFT::null();
    result.data[3] = MFT::null();
    return result;
  }

  __host__ __device__ 
  static TypeIn get(const This_t& value, int iE, int jE, int kE) {
    const Type data = value.data[jE + 2*iE];
    TypeIn results[2];
    MFT::decode(results[0], results[1], data);
    return results[kE];
  }
};

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
// Encode for multipliers/sums

static GMFloat2 GMFloat2_encode(GMTally1 val0, GMTally1 val1) {
  PackedDouble result = 0;
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::encode(result, val0, val1);
  return result;
}

//----------

static GMFloat3 GMFloat3_encode(GMTally1 val0, GMTally1 val1, GMTally1 val2) {
  GMFloat3 result; // here we should set = null to be super cautious
  const GMTally1 dummy = 0;
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::encode(result.data[0], val0, val1);
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::encode(result.data[1], val2, dummy);
  return result;
}

//-----------------------------------------------------------------------------
// Decode for multipliers/sums

static void GMFloat2_decode(GMTally1& __restrict__ val0,
                            GMTally1& __restrict__ val1,
                            const GMFloat2 v) {
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::decode(val0, val1, v);
}

//----------

static void GMFloat3_decode(GMTally1* __restrict__ val0,
                            GMTally1* __restrict__ val1,
                            GMTally1* __restrict__ val2,
                            GMFloat3 v) {
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::decode(*val0, *val1, v.data[0]);
  GMTally1 dummy;
  MetricFormatTraits<MetricFormat::PACKED_DOUBLE>::decode(*val2, dummy, v.data[1]);
}

//-----------------------------------------------------------------------------
// Get a table entry: 2x2

static GMTally1 GMTally2x2_get(GMTally2x2 tally2x2, int iE, int jE) {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);

  const uint64_t tally2 = tally2x2.data[iE];

  const GMTally1 result =
      jE == 0 ? tally2 % (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
              : tally2 / (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS);
  //COMET_ASSERT(result >= 0);
  COMET_ASSERT(result < (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS));
  return result;
}

//-----------------------------------------------------------------------------
// Get a table entry: 4x2

static GMTally1 GMTally4x2_get(GMTally4x2 tally4x2, int iE, int jE, int kE) {
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);

  const uint64_t tally2 = tally4x2.data[jE + 2 * iE];

  const GMTally1 result =
      kE == 0 ? tally2 % (((uint64_t)1) << GM_TALLY1_MAX_VALUE_BITS)
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

//=============================================================================
/// \brief Safely cast arithmetic value to strictly lower type size.

template<typename TO, typename TI>
static TO safe_cast(const TI v) {
#if ! defined(NDEBUG)
  static_assert(sizeof(TO) < sizeof(TI), "");
#endif
  // ISSUE: this doesn't compile under CUDA:
  // COMET_STATIC_ASSERT(sizeof(TO) < sizeof(TI)):
  COMET_ASSERT((TI)(TO)v == v);
  return (TO)v;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TYPES_HH_

//-----------------------------------------------------------------------------
