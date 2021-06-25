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

//-----------------------------------------------------------------------------
/// \brief Formula for a single 2-way CCC or DUO result.

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

  const Float_t fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const Float_t fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;

  const Float_t fij = recip_sumcij * rij;

  // Do the following to make floating point arithmetic order-independent.
  const Float_t fmin = fi < fj ? fi : fj;
  const Float_t fmax = fi < fj ? fj : fi;

  /* clang-format off */
  const Float_t result = multiplier * fij * (f_one - param * fmin) *
                                            (f_one - param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Formula for a single 3-way CCC or DUO result.

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

  const Float_t fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const Float_t fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;
  const Float_t fk = (f_one / COUNTED_BITS_PER_ELT) * recip_ck * sk;

  const Float_t fijk = recip_sumcijk * rijk;

  // Do the following to make floating point arithmetic order-independent.
  Float_t fmin = 0, fmid = 0, fmax = 0;
  utils::sort_3(fmin, fmid, fmax, fi, fj, fk);

  /* clang-format off */
  const Float_t result = multiplier * fijk * (f_one - param * fmin) *
                                             (f_one - param * fmid) *
                                             (f_one - param * fmax);
  /* clang-format on */

  return result;
}

//=============================================================================
/// \brief Formula for a single 2-way CCC or DUO result, higher level.

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
  COMET_ASSERT(iE == 0 || iE == 1);
  COMET_ASSERT(jE == 0 || jE == 1);
  COMET_ASSERT(si1+1 >= 0+1);
  COMET_ASSERT(sj1+1 >= 0+1);
  COMET_ASSERT(ci+1 >= 0+1);
  COMET_ASSERT(cj+1 >= 0+1);
  COMET_ASSERT(num_field_active > 0);

  if (env.is_threshold_tc()) {
    COMET_ASSERT(MF == MetricFormat::SINGLE);
    return (Float_t)ttable.get(iE, jE);
  }

  COMET_ASSERT(MF == MetricFormat::PACKED_DOUBLE);

  if (0 == ci || 0 == cj)
    return (Float_t)0;

  const Float_t f_one = 1;
  const Float_t recip_m = f_one / num_field_active;

  const GMTally1 rij = ttable.get(iE, jE);

  if (env.sparse()) {

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 si = iE == 0 ? si0 : si1;
    const GMTally1 sj = jE == 0 ? sj0 : sj1;

    // The following complex code is to reduce the number of divides.
    const Float_t f_ci = (Float_t) ci;
    const Float_t f_cj = (Float_t) cj;

    const Float_t f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const Float_t f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    if (env.is_using_xor()) {

      // Make adjustment for xor gemm.
      // See notes for 8x8 system solve to back out this result.
      const GMTally1 cij = (si1 + sj1 - ttable.get(1, 1)) / 2 +
                           (si1 + sj0 - ttable.get(1, 0)) / 2 +
                           (si0 + sj1 - ttable.get(0, 1)) / 2 +
                           (si0 + sj0 - ttable.get(0, 0)) / 2;
      if (0 == cij)
        return (Float_t)0;

      // The following complex code is to reduce the number of divides.
      const Float_t f_cij = (Float_t) cij;
      const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

      const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
      const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

      const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

      // Make adjustment for xor gemm.
      // See notes for 8x8 system solve to back out this result.
      const GMTally1 rij_true = (si + sj - rij) / 2;
      //const GMTally1 rij_true = iE &&  jE ? (si1 + sj1 - rij) / 2 :
      //                          iE && !jE ? (si1 + sj0 - rij) / 2 :
      //                         !iE &&  jE ? (si0 + sj1 - rij) / 2 :
      //                                      (si0 + sj0 - rij) / 2;

      return ccc_duo_value<CBPE, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij,
        env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

    } else { // if (!env.is_using_xor())

      GMTally1 cij = ttable.get(0, 0) +
                     ttable.get(0, 1) +
                     ttable.get(1, 0) +
                     ttable.get(1, 1);
      if (0 == cij)
        return (Float_t)0;

      // The following complex code is to reduce the number of divides.
      const Float_t f_cij = (Float_t) cij;
      const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

      const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
      const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

      const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

      const GMTally1 rij_true = rij;

      return ccc_duo_value<CBPE, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij,
        env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

    }  // if (env.is_using_xor())

  } else { // !env.sparse

    COMET_ASSERT(!(env.is_using_xor() && !env.sparse())); // should never occur

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * num_field_active - si1;
    const GMTally1 sj0 = CBPE * num_field_active - sj1;
    const GMTally1 si = iE == 0 ? si0 : si1;
    const GMTally1 sj = jE == 0 ? sj0 : sj1;

    const Float_t recip_sumcij = (f_one / (CBPE*CBPE)) * recip_m;

    return ccc_duo_value<CBPE, Float_t>(
      rij, si, sj, recip_m, recip_m, recip_sumcij,
      env_ccc_duo_multiplier<CBPE>(env), env.ccc_param());

  } // if (env.sparse())
}

//-----------------------------------------------------------------------------




//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_FORMATS_HH_

//-----------------------------------------------------------------------------
