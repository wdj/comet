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

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_FORMATS_HH_

//-----------------------------------------------------------------------------
