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

#include "env.hh"
#include "assertions.hh"
#include "types.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Helper class for system functions.

struct System {
  static int num_proc();
  static int proc_num();
  static bool is_proc_num_0() {return !proc_num();}
  static int compute_capability();
  static int pci_bus_id();
  static int pci_domain_id();
  static double time();
  static bool accel_last_call_succeeded();

  static bool is_in_parallel_region();
#if 0
  static bool is_in_parallel_region() {
#   if COMET_USE_OPENMP
      return omp_in_parallel();
#   else
      return false;
#   endif
  };
#endif

private:

#if defined COMET_USE_CUDA
  typedef cudaDeviceProp accelDeviceProp_t;
#elif defined COMET_USE_HIP
  typedef hipDeviceProp_t accelDeviceProp_t;
#else
  typedef int accelDeviceProp_t;
#endif

  static accelDeviceProp_t& get_device_prop();
};

namespace utils {

//-----------------------------------------------------------------------------
/// \brief Minimum of two scalars, native implementation for speed.

template<typename T>
__host__ __device__
T min(const T& i, const T& j) {

  const T r = i < j ? i : j;

  return r;
}

//-----------------------------------------------------------------------------
/// \brief Maximum of two scalars, native implementation for speed.

template<typename T>
__host__ __device__
T max(const T& i, const T& j) {

  const T r = i > j ? i : j;

  return r;
}

//-----------------------------------------------------------------------------
/// \brief Truncate to integer the quotient of integers.

template<typename T>
T trunc(const T& i, const T& n) {
  COMET_ASSERT(n > 0);
  COMET_ASSERT(i+1 >= 1);

  const T r = i / n;
  COMET_ASSERT(i >= r*n);
  COMET_ASSERT(i < r*n + n);

  return r;
}

//-----------------------------------------------------------------------------
/// \brief Integer floor of quotient of integers.

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
/// \brief Integer ceiling of quotient of integers.

template<typename T>
T ceil(const T& i, const T& n) {
  COMET_ASSERT(n > 0);

  const T r = i > 0 ? (i + n - 1) / n : i / n;
  // WARNING: may fail if unsigned type.
  COMET_ASSERT(i + n > r*n && i <= r*n);

  return r;
}

//-----------------------------------------------------------------------------
/// \brief Mathematical modulus if int value.

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
/// \brief Faster version of true mod, needed for special situation.

static int mod_fast(int i, int n) {
  return (i + n) % n;
}

//-----------------------------------------------------------------------------
/// \brief Upper bound of random numbers from generator.

static size_t randomize_max() {

  const size_t im = 714025;
  return im;
}

//-----------------------------------------------------------------------------
/// \brief Random number genrator.

template<typename T>
static size_t randomize(T i) {

  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;
  return static_cast<size_t>((i * ia + ic) % im);
}

//-----------------------------------------------------------------------------
/// \brief Ceiling of log2 of an integer.

template<class T>
static int log2(T n) {
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
/// \brief N choose K function.

template<class T>
static T nchoosek(T n, int k) {
  COMET_ASSERT(n+1 >= 0+1);
  COMET_ASSERT(k >= 0 && static_cast<T>(k) <= n);

  T numer = 1;
  int denom = 1;

  for (int i = 0; i < k; ++i) {
    // Subtract 2 because of possible sign bit
    // and also edge case of 2^i - 1 limit.
    COMET_ASSERT(log2(numer) + log2(n - i) <=
                 static_cast<int>(sizeof(T)*BITS_PER_BYTE - 2));
    numer *= (n - i);
    denom *= (i + 1);
  }

  return numer / denom;
}

//-----------------------------------------------------------------------------
/// \brief Convert string representing unsigned integer to arithmetic type.

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
/// \brief Population count of 1-bits in 8-bit word.

__host__ __device__
static int popc8(uint8_t x) {
  // Adapted from https://stackoverflow.com/questions/30688465/how-to-check-the-number-of-set-bits-in-an-8-bit-unsigned-char
  //x = x - ((x >> 1) & 0x55);
  //x = (x & 0x33) + ((x >> 2) & 0x33);
  //return (((x + (x >> 4)) & 0x0F) * 0x01);

  return (!!(((int)x)&  1)) +
         (!!(((int)x)&  2)) +
         (!!(((int)x)&  4)) +
         (!!(((int)x)&  8)) +
         (!!(((int)x)& 16)) +
         (!!(((int)x)& 32)) +
         (!!(((int)x)& 64)) +
         (!!(((int)x)&128));
}

template<typename In_t = uint8_t>
__host__ __device__
static int popc(In_t x) {return popc8(x);}

//-----------------------------------------------------------------------------
/// \brief Population count of 1-bits in 32-bit word.

__host__ __device__
static int popcs32(int32_t x) {

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
int popc<int32_t>(int32_t x) {return popcs32(x);}

//-----------------------------------------------------------------------------
/// \brief Population count of 1-bits in 32-bit word.

__host__ __device__
static int popc32(uint32_t x) {

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
int popc<uint32_t>(uint32_t x) {return popc32(x);}

//-----------------------------------------------------------------------------
/// \brief Population count of 1-bits in 64-bit word.

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
/// \brief Fast sort of 3 values.

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
/// \brief Atomic add with doubles.

__host__ __device__ static double atomic_add(double* address, double value) {

# if defined COMET_USE_CUDA && defined __CUDA_ARCH__

    return atomicAdd(address, value);

//# elif defined COMET_USE_HIP && defined __HIPCC__
# elif defined COMET_USE_HIP && defined __HIP_DEVICE_COMPILE__

    // TODO: fix this to work under HIP.

    double old = *address, assumed;
    do {
      assumed = old; old =
        __longlong_as_double(
          atomicCAS((unsigned long long int*)address,
                    __double_as_longlong(assumed),
      __double_as_longlong(value + assumed)));
    } while (assumed != old);
    return old;

# else

    COMET_ASSERT(!System::is_in_parallel_region() && "Unimplemented.");

    *address += value;
    return value;

# endif
}


//=============================================================================

} // namespace utils
} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_UTILS_HH_

//-----------------------------------------------------------------------------
