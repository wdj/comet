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

template<typename FLOAT, int COUNTED_BITS_PER_ELT>
static __host__ __device__ FLOAT ccc_duo_value(
  const GMTally1 rij,
  const GMTally1 si,
  const GMTally1 sj,
  const FLOAT recip_ci,
  const FLOAT recip_cj,
  const FLOAT recip_sumcij,
  const FLOAT multiplier,
  const FLOAT param) {

  const FLOAT f_one = 1;

  const FLOAT fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const FLOAT fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;

  const FLOAT fij = recip_sumcij * rij;

  // Do the following to make floating point arithmetic order-independent.
  const FLOAT fmin = fi < fj ? fi : fj;
  const FLOAT fmax = fi < fj ? fj : fi;

  /* clang-format off */
  const FLOAT result = multiplier * fij * (f_one - param * fmin) *
                                          (f_one - param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Formula for a single 3-way CCC or DUO result.

template<typename FLOAT, int COUNTED_BITS_PER_ELT>
static __host__ __device__ FLOAT ccc_duo_value(
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const FLOAT recip_ci,
  const FLOAT recip_cj,
  const FLOAT recip_ck,
  const FLOAT recip_sumcijk,
  const FLOAT multiplier,
  const FLOAT param) {

  const FLOAT f_one = 1;

  const FLOAT fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const FLOAT fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;
  const FLOAT fk = (f_one / COUNTED_BITS_PER_ELT) * recip_ck * sk;

  const FLOAT fijk = recip_sumcijk * rijk;

  // Do the following to make floating point arithmetic order-independent.
  FLOAT fmin = 0, fmid = 0, fmax = 0;
  utils::sort_3(fmin, fmid, fmax, fi, fj, fk);

  /* clang-format off */
  const FLOAT result = multiplier * fijk * (f_one - param * fmin) *
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
