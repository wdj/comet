//-----------------------------------------------------------------------------
/*!
 * \file   formulas.hh
 * \author Wayne Joubert
 * \date   Tue Mar  3 09:59:06 EST 2020
 * \brief  Formulas to support metrics calculations.
 * \note   Copyright (C) 2020 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_formulas_hh_
#define _comet_formulas_hh_

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

//-----------------------------------------------------------------------------

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_formats_hh_

//-----------------------------------------------------------------------------
