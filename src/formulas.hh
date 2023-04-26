//-----------------------------------------------------------------------------
/*!
 * \file   formulas.hh
 * \author Wayne Joubert
 * \date   Tue Mar  3 09:59:06 EST 2020
 * \brief  Formulas to support metrics calculations.
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

#ifndef _COMET_FORMULAS_HH_
#define _COMET_FORMULAS_HH_

#include "env.hh"
#include "assertions.hh"
#include "types.hh"
#include "utils.hh"

//=============================================================================

namespace comet {

namespace formulas {

//-----------------------------------------------------------------------------
/*!
 * \brief Formula for 2-way CCC or DUO result, lower-level version.
 *
 */
template<int COUNTED_BITS_PER_ELT, typename Float_t = double>
static __host__ __device__ Float_t ccc_duo_value(
  const GMTally1 rij,
  const GMTally1 si,
  const GMTally1 sj,
  const Float_t recip_ci,
  const Float_t recip_cj,
  const Float_t recip_sumcij,
  const Float_t multiplier,
  const Float_t param) {

  const Float_t f_one = 1;

  // C = COUNTED_BITS_PER_ELT = is_CCC ? 2 : 1

  const Float_t fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
                // = si / (C * ci)
  const Float_t fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;
                // = sj / (C * cj)

  const Float_t fij = recip_sumcij * rij; // = rij / sumcij

  // Do the following to make floating point arithmetic order-independent.

  const Float_t fmin = fi < fj ? fi : fj;
  const Float_t fmax = fi < fj ? fj : fi;

  // Main CCC/DUO formula.

  /* clang-format off */
  const Float_t result = multiplier * fij * (f_one - param * fmin) *
                                            (f_one - param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Formula for 3-way CCC or DUO result, lower-level version.
 *
 */
template<int COUNTED_BITS_PER_ELT, typename Float_t = double>
static __host__ __device__ Float_t ccc_duo_value(
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const Float_t recip_ci,
  const Float_t recip_cj,
  const Float_t recip_ck,
  const Float_t recip_sumcijk,
  const Float_t multiplier,
  const Float_t param) {

  const Float_t f_one = 1;

  // C = COUNTED_BITS_PER_ELT = is_CCC ? 2 : 1

  const Float_t fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
                // = si / (C * ci)
  const Float_t fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;
                // = sj / (C * cj)
  const Float_t fk = (f_one / COUNTED_BITS_PER_ELT) * recip_ck * sk;
                // = sk / (C * ck)

  const Float_t fijk = recip_sumcijk * rijk; // = rijk / sumcijk

  // Do the following to make floating point arithmetic order-independent.

  Float_t fmin = 0, fmid = 0, fmax = 0;
  utils::sort_3(fmin, fmid, fmax, fi, fj, fk);

  // Main CCC/DUO formula.

  /* clang-format off */
  const Float_t result = multiplier * fijk * (f_one - param * fmin) *
                                             (f_one - param * fmid) *
                                             (f_one - param * fmax);
  /* clang-format on */

  return result;
}

//=============================================================================
/*!
 * \brief Formula for 2-way CCC or DUO result, higher-level version.
 *
 */
template<int CBPE, int MF, typename Float_t = double>
static Float_t ccc_duo_value(
  const Tally2x2<MF>& ttable,
  int iE,
  int jE,
  GMTally1 si1,
  GMTally1 sj1,
  GMTally1 ci,
  GMTally1 cj,
  int num_field_active,
  CEnv& env) {
  COMET_ASSERT(0 == iE || 1 == iE);
  COMET_ASSERT(0 == jE || 1 == jE);
  COMET_ASSERT(si1+1 >= 0+1 && sj1+1 >= 0+1);
  COMET_ASSERT(ci+1 >= 0+1 && cj+1 >= 0+1);
  COMET_ASSERT(num_field_active > 0);

  if (env.is_threshold_tc()) {
    // If was computed on GPU, then no further work to do.
    COMET_ASSERT(MF == MetricFormat::SINGLE);
    return static_cast<Float_t>(ttable.get(iE, jE));
  }

  COMET_ASSERT(MF == MetricFormat::PACKED_DOUBLE);

  if (0 == ci || 0 == cj)
    return static_cast<Float_t>(0);

  const Float_t f_one = 1;
  const Float_t recip_m = f_one / num_field_active;

  const GMTally1 rij = ttable.get(iE, jE);

  //==========
  if (env.sparse()) {
  //==========

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 si = 0 == iE ? si0 : si1;
    const GMTally1 sj = 0 == jE ? sj0 : sj1;

    // This complex code is to reduce the number of divides: Part 1.
    const Float_t f_ci = static_cast<Float_t>(ci);
    const Float_t f_cj = static_cast<Float_t>(cj);

    const Float_t f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const Float_t f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    // If GPU computes XOR-GEMM (e.g., NVIDIA Turing), then compensate.

    //----------
    if (env.is_using_xor()) {
    //----------

      // Make adjustment for xor gemm.
      // See notes for 8x8 system solve on how this to backs out the result.
      const GMTally1 cij = (si1 + sj1 - ttable.get(1, 1)) / 2 +
                           (si1 + sj0 - ttable.get(1, 0)) / 2 +
                           (si0 + sj1 - ttable.get(0, 1)) / 2 +
                           (si0 + sj0 - ttable.get(0, 0)) / 2;
      if (0 == cij)
        return static_cast<Float_t>(0);

      // This complex code is to reduce the number of divides: Part 2.
      const Float_t f_cij = static_cast<Float_t>(cij);
      const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

      const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
      const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

      const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

      // Make adjustment for xor gemm.
      // See notes for 8x8 system solve to back out this result.
      const GMTally1 rij_true = (si + sj - rij) / 2;

      // Apply formula.

      return ccc_duo_value<CBPE, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij,
        env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

    //----------
    } else { // if (!env.is_using_xor())
    //----------

      GMTally1 cij = ttable.get(0, 0) + ttable.get(0, 1) +
                     ttable.get(1, 0) + ttable.get(1, 1);
      if (0 == cij)
        return static_cast<Float_t>(0);

      // This complex code is to reduce the number of divides: Part 2.
      const Float_t f_cij = static_cast<Float_t>(cij);
      const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

      const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
      const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

      const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

      const GMTally1 rij_true = rij;

      // Apply formula.

      return ccc_duo_value<CBPE, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij,
        env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

    //----------
    }  // if (env.is_using_xor())
    //----------

  //==========
  } else { // !env.sparse()
  //==========

    COMET_ASSERT((!(env.is_using_xor() && !env.sparse())) &&
      "Selected set of options unimplemented."); // should never occur
    // TODO: implement this, if needed.

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * num_field_active - si1;
    const GMTally1 sj0 = CBPE * num_field_active - sj1;
    const GMTally1 si = 0 == iE ? si0 : si1;
    const GMTally1 sj = 0 == jE ? sj0 : sj1;

    const Float_t recip_sumcij = (f_one / (CBPE*CBPE)) * recip_m;

    // Apply formula.

    return ccc_duo_value<CBPE, Float_t>(
      rij, si, sj, recip_m, recip_m, recip_sumcij,
      env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

  //==========
  } // if (env.sparse())
  //==========
}

//=============================================================================
/*!
 * \brief Formula for 3-way CCC or DUO result, higher-level short version.
 *
 */
template<int CBPE, typename Float_t = double>
static Float_t ccc_duo_value(
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const double recip_ci,
  const double recip_cj,
  const double recip_ck,
  const double recip_sumcijk,
  CEnv& env) {

  return ccc_duo_value<CBPE, Float_t>(
    rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk,
    env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());
}

//-----------------------------------------------------------------------------
/*!
 * \brief Formula for 3-way CCC or DUO result, higher-level version.
 *
 */
template<int CBPE, int MF, typename Float_t = double>
static Float_t ccc_duo_value(
  const Tally4x2<MF>& ttable,
  int iE,
  int jE,
  int kE,
  GMTally1 si1,
  GMTally1 sj1,
  GMTally1 sk1,
  GMTally1 ci,
  GMTally1 cj,
  GMTally1 ck,
  int num_field_active,
  CEnv& env) {
  COMET_ASSERT(0 == iE || 1 == iE);
  COMET_ASSERT(0 == jE || 1 == jE);
  COMET_ASSERT(0 == kE || 1 == kE);
  COMET_ASSERT(si1+1 >= 0+1 && sj1+1 >= 0+1 && sk1+1 >= 0+1);
  COMET_ASSERT(ci+1 >= 0+1 && cj+1 >= 0+1 && ck+1 >= 0+1);
  COMET_ASSERT(num_field_active > 0);

  if (env.is_threshold_tc()) {
    // If was computed on GPU, then no further work to do.
    COMET_ASSERT(MF == MetricFormat::SINGLE);
    return static_cast<Float_t>(ttable.get(iE, jE, kE));
  }

  COMET_ASSERT(MF == MetricFormat::PACKED_DOUBLE);

  if (0 == ci || 0 == cj || 0 == ck)
    return static_cast<Float_t>(0);

  const Float_t f_one = 1;
  const Float_t recip_m = f_one / num_field_active;

  const GMTally1 rijk = ttable.get(iE, jE, kE);

  //==========
  if (env.sparse()) {
  //==========

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 sk0 = CBPE * ck - sk1;
    const GMTally1 si = 0 == iE ? si0 : si1;
    const GMTally1 sj = 0 == jE ? sj0 : sj1;
    const GMTally1 sk = 0 == kE ? sk0 : sk1;

    // TODO: it may be possible to decrease the number of divides
    // here - see GMMetrics_ccc_get_from_index_2.
    const Float_t recip_ci = f_one / ci;
    const Float_t recip_cj = f_one / cj;
    const Float_t recip_ck = f_one / ck;

    // NOTE: we are assuming here that xor gemm values have already been
    // adjusted to true gemm values.  Changing this in the future would
    // require computing and storing matX sum values somewhere and then using
    // them here.

    // NOTE: this is in fact currently never called if using xor gemm.
    COMET_ASSERT(!env.is_using_xor());

    const GMTally1 cijk = ttable.get(0, 0, 0) + ttable.get(0, 0, 1) +
                          ttable.get(0, 1, 0) + ttable.get(0, 1, 1) +
                          ttable.get(1, 0, 0) + ttable.get(1, 0, 1) +
                          ttable.get(1, 1, 0) + ttable.get(1, 1, 1);
    if (0 == cijk)
      return static_cast<Float_t>(0);

    const Float_t recip_sumcijk = f_one / cijk;

    // Apply formula.

    return ccc_duo_value<CBPE>(
      rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk, env);

  //==========
  } else { // !env.sparse()
  //==========

    COMET_ASSERT((!(env.is_using_xor() && !env.sparse())) &&
      "Selected set of options unimplemented."); // should never occur
    // TODO: implement this, if needed.

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si = 0 == iE ? (CBPE * num_field_active - si1) : si1;
    const GMTally1 sj = 0 == jE ? (CBPE * num_field_active - sj1) : sj1;
    const GMTally1 sk = 0 == kE ? (CBPE * num_field_active - sk1) : sk1;

    const Float_t recip_sumcijk = (f_one / (CBPE*CBPE*CBPE)) * recip_m;

    // Apply formula.

    return ccc_duo_value<CBPE>(
      rijk, si, sj, sk, recip_m, recip_m, recip_m, recip_sumcijk, env);

  //==========
  } // if (env.sparse())
  //==========
}

//=============================================================================

} // namespace formulas

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_FORMATS_HH_

//-----------------------------------------------------------------------------
