//-----------------------------------------------------------------------------
/*!
 * \file   utils.hh
 * \author Wayne Joubert
 * \date   Sat Nov 16 10:04:31 EST 2019
 * \brief  Miscellaneous utilities
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

#ifndef _COMET_UTILS_HH_
#define _COMET_UTILS_HH_

#include "cstdint"
#include "cstddef"
#include "cstring"
#include "assert.h"
#include "float.h"
#include "algorithm"
#include "vector"
#include "limits"

#include "env.hh"
#include "assertions.hh"
#include "types.hh"

//=============================================================================

namespace comet {

// Forward declaration.
class CEnv;

namespace System {

//-----------------------------------------------------------------------------
// System helper functions.

int num_proc();
int proc_num();
static bool is_proc_num_0() {
  return !proc_num();
}
int compute_capability();
int pci_bus_id();
int pci_domain_id();
double time();
bool accel_last_call_succeeded();

//-----------------------------------------------------------------------------
/*!
 * \brief Are we in an openmp parallel region.
 *
 */
static bool is_in_parallel_region() {
# if COMET_USE_OPENMP
    return omp_in_parallel();
# else
    return false;
# endif
};

//-----------------------------------------------------------------------------

} // namespace System

//=============================================================================

namespace utils {

//-----------------------------------------------------------------------------
/*!
 * \brief Minimum of two scalars, native implementation for speed.
 *
 */
template<typename T>
__host__ __device__
T min(const T& i, const T& j) {

  const T r = i < j ? i : j;

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Maximum of two scalars, native implementation for speed.
 *
 */
template<typename T>
__host__ __device__
T max(const T& i, const T& j) {

  const T r = i > j ? i : j;

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Truncate the quotient of two nonnegative integers to an integer.
 *
 */
template<typename T>
T trunc(const T& i, const T& n) {
  COMET_ASSERT(n > 0);
  COMET_ASSERT(i+1 >= 0+1);

  const T r = i / n;
  COMET_ASSERT(i >= r*n);
  COMET_ASSERT(i < r*n + n);

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Integer floor of quotient of integers.
 *
 */
template<typename T>
T floor(const T& i, const T& n) {
  COMET_STATIC_ASSERT(std::is_signed<T>::value);
  COMET_ASSERT(n > 0);

  const T r = i >= 0 ? i / n : (i + 1 - n) / n;
  COMET_ASSERT(i >= r*n);
  COMET_ASSERT(i < r*n + n);

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Integer ceiling of quotient of integers.
 *
 */
template<typename T>
T ceil(const T& i, const T& n) {
  COMET_ASSERT(n > 0);

  const T r = i > 0 ? (i + n - 1) / n : i / n;
  // WARNING: may fail if unsigned type.
  COMET_ASSERT(i + n > r*n && i <= r*n);

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Mathematical modulus of integer value.
 *
 */
template<typename T>
T mod_i(const T& i, const T& n) {
  COMET_STATIC_ASSERT((std::is_same<T,int>::value));
  COMET_ASSERT(n > 0);

  const T r = i - n * floor(i, n);
  COMET_ASSERT(r >= 0 && r < n);
  COMET_ASSERT((r-i) % n == 0);

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Faster version of true modulus, needed for a special situation.
 *
 */
static int mod_fast(const int& i, const int& n) {
  return (i + n) % n;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Upper bound of random numbers from generator.
 *
 */
static size_t randomize_max() {

  const size_t im = 714025;
  return im;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Random number genrator.
 *
 */
template<typename T>
static size_t randomize(const T& i) {

  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;
  return static_cast<size_t>((i * ia + ic) % im);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Ceiling of log2 of an integer.
 *
 */
template<class T>
static int log2(const T& n) {
  if (n <= 1)
    return 0;

  size_t n_ = n - 1;

  int r = 0; 
  for ( ; r <= static_cast<int>(sizeof(T)*BITS_PER_BYTE); ++r) {
    if (n_ == 0) {
      break;
    }
    n_ >>= 1;
  }

  COMET_ASSERT(r >= 1);

  COMET_ASSERT(n <= (static_cast<T>(1) << r));
  COMET_ASSERT(!(n <= (static_cast<T>(1) << (r-1))));

  return r;
}

//-----------------------------------------------------------------------------
/*!
 * \brief N choose K function.
 *
 */
template<class T>
static T nchoosek(const T& n, const int& k) {
  COMET_ASSERT(n+1 >= 0+1);
  COMET_ASSERT(k >= 0 && static_cast<T>(k) <= n);

  T numer = 1;
  int denom = 1;

  for (int i = 0; i < k; ++i) {
    // Subtract 2 because of possible sign bit
    // and also edge case of 2^i - 1 limit.
    COMET_ASSERT(log2(numer) + log2(n - i) <=
                 static_cast<int>(sizeof(T)*BITS_PER_BYTE) - 2);
    numer *= (n - i);
    denom *= (i + 1);
  }

  return numer / denom;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Convert string representing unsigned integer to arithmetic type.
 *
 */
template<class T>
T strtoarith(char* const str) {

  T result = 0;

  for (int i = 0; ; ++i) {

    if (!str[i])
      break;

    if (' ' == str[i] || '	' == str[i])
      continue;

    COMET_INSIST(str[i] >= '0' && str[i] <= '9');

    const int digit = str[i] - '0';
 
    const T result_new = 10 * result + digit;
    COMET_INSIST((result_new-digit) / 10 == result);
    result = result_new;
  }

  return result;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Count of 1-bits in 8-bit word, with templated version.
 *
 */
__host__ __device__
static int popc8(const uint8_t& x) {
  // Adapted from https://stackoverflow.com/questions/30688465/how-to-check-the-number-of-set-bits-in-an-8-bit-unsigned-char
  //x = x - ((x >> 1) & 0x55);
  //x = (x & 0x33) + ((x >> 2) & 0x33);
  //return (((x + (x >> 4)) & 0x0F) * 0x01);

  const auto ix = static_cast<int>(x);

  return (!!(ix&  1)) +
         (!!(ix&  2)) +
         (!!(ix&  4)) +
         (!!(ix&  8)) +
         (!!(ix& 16)) +
         (!!(ix& 32)) +
         (!!(ix& 64)) +
         (!!(ix&128));
}

template<typename In_t = uint8_t>
__host__ __device__
static int popc(In_t x) {return popc8(x);}

//-----------------------------------------------------------------------------
/*!
 * \brief Count of 1-bits in 32-bit word, with templated version.
 *
 */
__host__ __device__
static int popc32(int32_t x) {

  // Adapted from Hacker's Delight, 2nd ed.
   x = x - ((x >> 1) & 0x55555555);
   x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
   x = (x + (x >> 4)) & 0x0F0F0F0F;
   x = x + (x >> 8);
   x = x + (x >> 16);

   return x & 0x0000003F;
}

template<>
__host__ __device__
int popc<int32_t>(int32_t x) {return popc32(x);}

//-----------------------------------------------------------------------------
/*!
 * \brief Count of 1-bits in unsigned 32-bit word, with templated version.
 *
 */
__host__ __device__
static int popcu32(uint32_t x) {

  // Adapted from Hacker's Delight, 2nd ed.
   x = x - ((x >> 1) & 0x55555555);
   x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
   x = (x + (x >> 4)) & 0x0F0F0F0F;
   x = x + (x >> 8);
   x = x + (x >> 16);

   return x & 0x0000003F;
}

template<>
__host__ __device__
int popc<uint32_t>(uint32_t x) {return popcu32(x);}

//-----------------------------------------------------------------------------
/*!
 * \brief Count of 1-bits in 64-bit word, with templated version.
 *
 */
__host__ __device__
static int popc64(uint64_t x) {

  // Adapted from https://en.wikipedia.org/wiki/Hamming_weight
  const uint64_t m1 = 0x5555555555555555;
  const uint64_t m2 = 0x3333333333333333;
  const uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
  const uint64_t h01 = 0x0101010101010101;
  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;

  return (x * h01) >> 56;
}

template<>
__host__ __device__
int popc<uint64_t>(uint64_t x) {return popc64(x);}

//-----------------------------------------------------------------------------
/*!
 * \brief Fast sort of 3 values.
 *
 */
template<typename T>
static __host__ __device__ void sort_3(T& min_, T& mid_, T& max_,
                                       const T& a, const T& b, const T& c) {
  if (a > b) {
    if (a > c) {
      max_ = a;
      if (b > c) {
        mid_ = b;
        min_ = c;
      } else {
        mid_ = c;
        min_ = b;
      }
    } else {
      mid_ = a;
      max_ = c;
      min_ = b;
    }
  } else {
    if (b > c) {
      max_ = b;
      if (a > c) {
        mid_ = a;
        min_ = c;
      } else {
        mid_ = c;
        min_ = a;
      }
    } else {
      mid_ = b;
      max_ = c;
      min_ = a;
    }
  }

  COMET_ASSERT(min_ <= mid_);
  COMET_ASSERT(mid_ <= max_);

#if defined(COMET_ASSERTIONS_ON) && !defined(COMET_DEVICE_COMPILE)
  std::vector<T> v{a, b, c};
  std::sort(v.begin(), v.end());
  COMET_ASSERT(min_ == v[0]);
  COMET_ASSERT(mid_ == v[1]);
  COMET_ASSERT(max_ == v[2]);
#endif
}

//-----------------------------------------------------------------------------
/*!
 * \brief Atomic add with doubles.
 *
 */
__host__ __device__ static double atomic_add(double* address, double value) {

# if defined COMET_USE_CUDA && defined __CUDA_ARCH__

    return atomicAdd(address, value);

//# elif defined COMET_USE_HIP && defined __HIPCC__
# elif defined COMET_USE_HIP && defined __HIP_DEVICE_COMPILE__

    return atomicAdd(address, value);

#if 0
    // NOTE: this needs to be checked before using.

    double old = *address, assumed;
    do {
      assumed = old; old =
        __longlong_as_double(
          atomicCAS(reinterpret_cast<unsigned long long int*>(address),
                    __double_as_longlong(assumed),
      __double_as_longlong(value + assumed)));
    } while (assumed != old);
    return old;
#endif

# else

    COMET_ASSERT(!System::is_in_parallel_region() && "Unimplemented.");

    *address += value;
    return value;

# endif
}

//-----------------------------------------------------------------------------
/*!
 * \brief Fill an array with NaNs (for debug build only).
 *
 */
template<typename Float_t>
static void fill_nan(Float_t* a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);

# ifdef COMET_ASSERTIONS_ON
    const Float_t value = std::numeric_limits<Float_t>::quiet_NaN();
    for (size_t i=0; i<n; ++i)
      a[i] = value;
# endif
}

//-----------------------------------------------------------------------------

void* malloc(size_t n, CEnv& env);
void free(void* p, size_t n, CEnv& env);

size_t array_cksum(const unsigned char* const a, size_t n);

//=============================================================================

} // namespace utils
} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_UTILS_HH_

//-----------------------------------------------------------------------------
