//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_accessors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, accessor functions.
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

#ifndef _COMET_METRICS_2WAY_ACCESSORS_I_HH_
#define _COMET_METRICS_2WAY_ACCESSORS_I_HH_

#include "metrics_2way_indexing.i.hh"
#include "formulas.hh"

//=============================================================================

namespace comet {

//=============================================================================
// Accessors: value from (contig) index: derived

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 2-way CCC value.

static void GMMetrics_ccc_check_size_nofp_2(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (!BuildHas::INT128)
    return;

  if (env->metric_type() != MetricType::CCC ||
      env->num_way() != NumWay::_2 || ! env->are_ccc_params_default())
    return;

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 2+lm + 2+lm + 2+lm; 

  const int shift = mantissa_digits<GMFloat>() - 3; // Note num/denom <= 4.5 < 1<<3

  // To make this less restrictive, see 3-way counterpart code.
  COMET_INSIST_INTERFACE(env, lnum + shift < 128 &&
                         "Number of fields too large.");
}

//-----------------------------------------------------------------------------
/// \brief Formula for CCC 2-way metric using 128 bit integer arithmetic.

static double GMMetrics_ccc_value_nofp_2(GMMetrics* metrics,
                                         const GMTally1 rij,
                                         const GMTally1 si,
                                         const GMTally1 sj,
                                         const GMTally1 ci,
                                         const GMTally1 cj,
                                         const GMTally1 cij,
                                         CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(BuildHas::INT128);

  typedef double Float_t;
  typedef BasicTypes::BigUInt UInt128;

  const UInt128 num = rij * (UInt128)(3 * ci - 1 * si) *
                            (UInt128)(3 * cj - 1 * sj);

  const UInt128 denom = 2 * cij * (UInt128)ci * (UInt128)cj;

  const int shift = mantissa_digits<Float_t>() - 3; // Note num/denom <= 4.5 < 1<<3

  // This should be an integer with no more than
  // mantissa_digits<Float_t>() binary digits
  const UInt128 int_ratio = (num << shift) / denom;

  // Convert to floting point and then adjust exponent.
  const Float_t result = ( (Float_t) int_ratio ) /
                         ( (Float_t) ( ((size_t)1) << shift ) );

  return result;
}


//-----------------------------------------------------------------------------
/// \brief Accessor for 2-way CCC metric computed with 128 bit int arithmetic.

static double GMMetrics_ccc_get_from_index_nofp_2(GMMetrics* metrics,
                                                  NML_t index,
                                                  int iE,
                                                  int jE,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_metrics_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NumWay::_2);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(env->are_ccc_params_default());
  COMET_ASSERT(BuildHas::INT128);

  typedef double Float_t;

  const auto ttable = Metrics_elt_const<GMTally2x2>(*metrics, index, *env);
  const GMTally1 rij = GMTally2x2_get(ttable, iE, jE);

  const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArrayId::S>(*metrics,
    index, *env);
  GMTally1 si1, sj1;
  GMFloat2_decode(si1, sj1, si1_sj1);

  GMTally1 ci = 0, cj = 0, cij = 0;

  if (env->sparse()) {
    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArrayId::C>(*metrics,
      index, *env);
    GMFloat2_decode(ci, cj, ci_cj);

    cij = GMTally2x2_get(ttable, 0, 0) + GMTally2x2_get(ttable, 0, 1) +
          GMTally2x2_get(ttable, 1, 0) + GMTally2x2_get(ttable, 1, 1);

    if (0 == ci || 0 == cj || 0 == cij) {
      return (Float_t)0;
    }
  } else {
    const int m = metrics->num_field_active;

    ci = m;
    cj = m;

    cij = 4 * m;
  }

  const GMTally1 si = iE == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = jE == 0 ? (2 * cj - sj1) : sj1;

  return GMMetrics_ccc_value_nofp_2(metrics, rij, si, sj, ci, cj, cij, env);
}

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 2-way CCC or DUO result, implementation.

template<int COUNTED_BITS_PER_ELT, typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_2_impl( GMMetrics& metrics,
  NML_t index, int iE, int jE, CEnv& env) {
  COMET_ASSERT(index < metrics.num_metrics_local); // && index >= 0
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.is_metric_type_bitwise());
  COMET_ASSERT(env.counted_bits_per_elt() == COUNTED_BITS_PER_ELT);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);

  if (env.is_threshold_tc()) {
    typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
    const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
    TTable_t::TypeIn result = TTable_t::get(ttable, iE, jE);
    return static_cast<FloatResult_t>(result);
  }

  typedef double Float_t; // Perform all calcs in double until return.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArrayId::S>(metrics,
    index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(si1, sj1, si1_sj1);

  enum {MF = MetricFormat::PACKED_DOUBLE};

  typedef Tally2x2<MF> TTable_t;
  const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);

  Float_t result_floatcalc = 0;

  if (env.sparse()) {

    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArrayId::C>(metrics,
      index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(ci, cj, ci_cj);

    result_floatcalc = formulas::ccc_duo_value<CBPE, MF>(ttable,
      iE, jE, si1, sj1, ci, cj, metrics.num_field_active, env);

  } else { // !env.sparse

    GMTally1 ci = CBPE * CBPE * metrics.num_field_active;
    GMTally1 cj = CBPE * CBPE * metrics.num_field_active;

    result_floatcalc = formulas::ccc_duo_value<CBPE, MF>(ttable,
      iE, jE, si1, sj1, ci, cj, metrics.num_field_active, env);

  } // if (env.sparse())

#if 0

  const Float_t f_one = 1;
  const Float_t recip_m = metrics.recip_m;

  const auto ttable = Metrics_elt_const<GMTally2x2>(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(ttable, iE, jE);

  const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArrayId::S>(metrics,
    index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(si1, sj1, si1_sj1);

  Float_t result_floatcalc = 0;

  if (env.sparse()) {

    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArrayId::C>(metrics,
      index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(ci, cj, ci_cj);
    if (0 == ci || 0 == cj)
      return (FloatResult_t)0;

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
      const GMTally1 cij = (si1 + sj1 - GMTally2x2_get(ttable, 1, 1)) / 2 +
                           (si1 + sj0 - GMTally2x2_get(ttable, 1, 0)) / 2 +
                           (si0 + sj1 - GMTally2x2_get(ttable, 0, 1)) / 2 +
                           (si0 + sj0 - GMTally2x2_get(ttable, 0, 0)) / 2;
      if (0 == cij)
        return (FloatResult_t)0;

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

      result_floatcalc = formulas::ccc_duo_value<COUNTED_BITS_PER_ELT, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij, 
        env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(env), env.ccc_param());

    } else { // if (!env.is_using_xor())

      GMTally1 cij = GMTally2x2_get(ttable, 0, 0) +
                     GMTally2x2_get(ttable, 0, 1) +
                     GMTally2x2_get(ttable, 1, 0) +
                     GMTally2x2_get(ttable, 1, 1);
      if (0 == cij)
        return (FloatResult_t)0;

      // The following complex code is to reduce the number of divides.
      const Float_t f_cij = (Float_t) cij;
      const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

      const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
      const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

      const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

      const GMTally1 rij_true = rij;

      result_floatcalc = formulas::ccc_duo_value<COUNTED_BITS_PER_ELT, Float_t>(
        rij_true, si, sj, recip_ci, recip_cj, recip_sumcij, 
        env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(env), env.ccc_param());

    }  // if (env.is_using_xor())

  } else { // !env.sparse

    COMET_ASSERT(!(env.is_using_xor() && !env.sparse())); // should never occur

    COMET_ASSERT(metrics.num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * metrics.num_field_active - si1;
    const GMTally1 sj0 = CBPE * metrics.num_field_active - sj1;
    const GMTally1 si = iE == 0 ? si0 : si1;
    const GMTally1 sj = jE == 0 ? sj0 : sj1;

    const Float_t recip_sumcij = (f_one / (CBPE*CBPE)) * recip_m;

    result_floatcalc = formulas::ccc_duo_value<COUNTED_BITS_PER_ELT, Float_t>(
      rij, si, sj, recip_m, recip_m, recip_sumcij, 
      env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(env), env.ccc_param());

  } // if (env.sparse())

#endif

  if (BuildHas::INT128 && env.metric_type() == MetricType::CCC &&
      env.are_ccc_params_default()) {
    const Float_t result_intcalc = GMMetrics_ccc_get_from_index_nofp_2(&metrics,
                                         index, iE, jE, &env);

    // TODO: CHECK floating point type here
    const double eps = 1. /
      ( ((size_t)1) << (mantissa_digits<FloatResult_t>() - 6) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      fprintf(stderr,
              "Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
              (double)result_floatcalc, (double)result_intcalc);
      COMET_INSIST(diff < eps);
    }
  }

  return static_cast<FloatResult_t>(result_floatcalc);
}

//-----------------------------------------------------------------------------
/// \brief Helper class for ccc_duo_get_from_index_2 template specialization.

// See https://stackoverflow.com/questions/12683165/partial-specialization-of-templates-with-integer-parameters

template<int CBPE, typename FloatResult_t>
struct MetricsIndexCCCDUO2Helper {
  static FloatResult_t impl(
  GMMetrics& metrics, NML_t index, int iE, int jE, CEnv& env) {
    return Metrics_ccc_duo_get_2_impl<
      CBPE, FloatResult_t>(metrics, index, iE, jE, env);;
  }
};

//-----------------------------------------------------------------------------

template<typename FloatResult_t>
struct MetricsIndexCCCDUO2Helper<CBPE::NONE, FloatResult_t> {
  static FloatResult_t impl(
  GMMetrics& metrics, NML_t index, int iE, int jE, CEnv& env) {
    return env.counted_bits_per_elt() == CBPE::CCC ?
      MetricsIndexCCCDUO2Helper<CBPE::CCC, FloatResult_t>::impl(
        metrics, index, iE, jE, env) :
      MetricsIndexCCCDUO2Helper<CBPE::DUO, FloatResult_t>::impl(
        metrics, index, iE, jE, env);
  }
};

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 2-way CCC or DUO result.

template<int CBPE = CBPE::NONE, typename FloatResult_t = GMFloat>
static FloatResult_t Metrics_ccc_duo_get_2(
  GMMetrics& metrics, NML_t index, int iE, int jE, CEnv& env) {
    return MetricsIndexCCCDUO2Helper<CBPE, FloatResult_t>::impl(
      metrics, index, iE, jE, env);
}

//-----------------------------------------------------------------------------

template<int CBPE = CBPE::NONE, typename FloatResult_t = GMFloat>
static FloatResult_t Metrics_ccc_duo_get_2(
  GMMetrics& metrics, NML_t index, int entry_num, CEnv& env) {

  if (env.is_shrink()) {
    COMET_ASSERT(0 == entry_num);
    typedef MetricFormatTraits<MetricFormat::SINGLE>::TypeIn TypeIn;
    return (FloatResult_t)Metrics_elt_const<TypeIn>(metrics, index, env);
  }

  COMET_ASSERT(0 <= entry_num && entry_num < env.pow2_num_way());
  const int iE = entry_num / 2;
  const int jE = entry_num % 2;;

  return MetricsIndexCCCDUO2Helper<CBPE, FloatResult_t>::impl(  
    metrics, index, iE, jE, env);
}

//-----------------------------------------------------------------------------

template<typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_2(
  GMMetrics& metrics, NML_t index, int entry_num, CEnv& env) {

  return Metrics_ccc_duo_get_2<CBPE::NONE, FloatResult_t>(metrics, index,
    entry_num, env);
}

//-----------------------------------------------------------------------------

template<typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_2(
  GMMetrics& metrics, NML_t index, int iE, int jE, CEnv& env) {

  return Metrics_ccc_duo_get_2<CBPE::NONE, FloatResult_t>(metrics, index,
    iE, jE, env);
}

//-----------------------------------------------------------------------------

//template<typename FloatResult_t>
//static FloatResult_t Metrics_ccc_duo_get_2(
//  GMMetrics& metrics, NML_t index, int iE, int jE, CEnv& env) {
//    return MetricsIndexCCCDUO2Helper<CBPE::NONE, FloatResult_t>::impl(
//      metrics, index, iE, jE, env);
//}

#if 0
//-----------------------------------------------------------------------------
/// \brief Templatized Check if any table value may exceed threshold.
/// NOTE: this function is currently unused.

template<int COUNTED_BITS_PER_ELT>
static bool Metrics_ccc_duo_threshold_detector_2(GMMetrics& metrics,
  const size_t index, CEnv& env) {
  COMET_ASSERT(index < metrics.num_metrics_local); // && index >= 0
  COMET_ASSERT(env.num_way() == NumWay::_2);

  // if no active threshold, then always pass threshold criterion.
  if (!env.thresholds_.is_active()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    return true;
  }

  // if is_shrink, assume a threshold pass my exist, don't take time to check.
  if (env.is_shrink()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    return true;
  }

  // If threshold_tc, then look for any prethresholded values not set to zero.
  if (env.is_threshold_tc()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
    const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
    for (int iE = 0; iE < 2; ++iE) {
      for (int jE = 0; jE < 2; ++jE) {
        if (TTable_t::get(ttable, iE, jE) != (TTable_t::TypeIn)0)
          return true;
      }
    }
    return false;
  }

  COMET_INSIST(env.is_using_threshold_detector());

  // this is here because xor stuff not implemented below.
  COMET_ASSERT(!env.is_using_xor()); // should never occur.

  typedef double Float_t; // Perform all calcs in double.

  enum {CBPE = COUNTED_BITS_PER_ELT};

// FIXTHRESHOLD - here and below

  const Float_t threshold_eff = env.threshold_eff();

  if (env.sparse()) {

    const Float_t f_one = 1;

    const auto ttable = Metrics_elt_const<GMTally2x2>(metrics, index, env);
    const GMTally1 rij00 = GMTally2x2_get(ttable, 0, 0);
    const GMTally1 rij01 = GMTally2x2_get(ttable, 0, 1);
    const GMTally1 rij10 = GMTally2x2_get(ttable, 1, 0);
    const GMTally1 rij11 = GMTally2x2_get(ttable, 1, 1);

    const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArrayId::S>(metrics,
      index, env);
    GMTally1 si1 = 0, sj1 = 0;
    GMFloat2_decode(si1, sj1, si1_sj1);

    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArrayId::C>(metrics,
      index, env);
    GMTally1 ci = 0, cj = 0;
    GMFloat2_decode(ci, cj, ci_cj);

    GMTally1 cij = rij00 + rij01 + rij10 + rij11;
    if (ci == 0 || cj == 0 || cij == 0) {
      return env.pass_threshold(0);
    }

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits

    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;

    const Float_t f_ci = (Float_t) ci;
    const Float_t f_cj = (Float_t) cj;
    const Float_t f_cij = (Float_t) cij;

    const Float_t f_cicj = f_ci * f_cj;

    const Float_t recip_2 = f_one / f_cicj;

    const Float_t recip_ci = f_cj * recip_2;
    const Float_t recip_cj = f_ci * recip_2;

    const Float_t fi0_2 = recip_ci * si0;
    const Float_t fi1_2 = recip_ci * si1;
    const Float_t fj0_2 = recip_cj * sj0;
    const Float_t fj1_2 = recip_cj * sj1;

    // Do simple algebra on threshold inequality to get all constnts on RHS

    const Float_t ccc_duo_multiplier = (Float_t)env_ccc_duo_multiplier<CBPE>(env);
    const Float_t ccc_param = (Float_t)env.ccc_param();

    const Float_t threshold_multiplier =
             ((f_one*CBPE*CBPE) / ccc_duo_multiplier) * f_cij;

    // Let a few extra values pass the thresold check just to make sure
    const Float_t roundoff_fuzz = 1.e-5;

    const Float_t f = threshold_eff * threshold_multiplier * (f_one-roundoff_fuzz);

    const Float_t v00 = rij00 * ((f_one*CBPE) - ccc_param * fi0_2) *
                                ((f_one*CBPE) - ccc_param * fj0_2);
    const Float_t v01 = rij01 * ((f_one*CBPE) - ccc_param * fi0_2) *
                                ((f_one*CBPE) - ccc_param * fj1_2);
    const Float_t v10 = rij10 * ((f_one*CBPE) - ccc_param * fi1_2) *
                                ((f_one*CBPE) - ccc_param * fj0_2);
    const Float_t v11 = rij11 * ((f_one*CBPE) - ccc_param * fi1_2) *
                                ((f_one*CBPE) - ccc_param * fj1_2);

    // Verify the algebra.
    COMET_ASSERT(fabs(v00 - threshold_multiplier *
      Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index, 0, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v01 - threshold_multiplier *
      Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index, 0, 1, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v10 - threshold_multiplier *
      Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index, 1, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v11 - threshold_multiplier *
      Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index, 1, 1, env) )
             < roundoff_fuzz * threshold_multiplier);

    return v00 > f || v01 > f || v10 > f || v11 > f;

  } // if (env.sparse())

  // Non-sparse case (less well-optimized).

  const Float_t v00 = Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index,
                                                           0, 0, env);
  const Float_t v01 = Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index,
                                                           0, 1, env);
  const Float_t v10 = Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index,
                                                           1, 0, env);
  const Float_t v11 = Metrics_ccc_duo_get_2<CBPE, Float_t>(metrics, index,
                                                           1, 1, env);

  return env.pass_threshold(v00) || env.pass_threshold(v01) ||
         env.pass_threshold(v10) || env.pass_threshold(v11);
}
#endif

//=============================================================================
// Accessors: value from (local) coord: set: 2-way.

template<typename T, int MA = MetricsArrayId::M>
static T& Metrics_elt_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || env.proc_num_vector() == j_block);
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_2(metrics, i, j, j_block, env);

  return Metrics_elt<T, MA>(metrics, index, env);
}

//=============================================================================
// Accessors: value from (local) coord: get: 2-way.

template<typename T>
static T Metrics_elt_const_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || env.proc_num_vector() == j_block);
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_2(metrics, i, j, j_block, env);

  return Metrics_elt_const<T>(metrics, index, env);
}

//=============================================================================
// Accessors: value from index: get: 2-way.

template<typename FloatResult_t>
FloatResult_t GMMetrics_get_2(GMMetrics& metrics,
  NML_t index, int iE, int jE, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(index >= 0 && index < metrics.num_metrics_local);
  COMET_ASSERT(iE >= 0 && iE < env.ijkE_max());
  COMET_ASSERT(jE >= 0 && jE < env.ijkE_max());

  COMET_ASSERT(!env.is_shrink()); //FIX
  const auto result =
    env.metric_type() == MetricType::CZEK ?
      Metrics_elt_const<FloatResult_t>(metrics, index, env) :
      Metrics_ccc_duo_get_2<FloatResult_t>(metrics, index, iE, jE, env);

  return result;
}

//=============================================================================
// Accessors: value from (global) coord: get: 2-way.

template<typename FloatResult_t>
FloatResult_t GMMetrics_get_2(GMMetrics& metrics,
  NV_t iG, NV_t jG, int iE, int jE, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(iE >= 0 && iE < env.ijkE_max());
  COMET_ASSERT(jE >= 0 && jE < env.ijkE_max());
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_2(metrics, iG, jG, env);

  const auto result = GMMetrics_get_2<FloatResult_t>(metrics, index, iE, jE, env);

  return result;
}

//=============================================================================
// Accessors: determine whether pass threshold: 2-way.

template<typename Float_t = GMFloat>
bool Metrics_is_pass_threshold_noshrink(GMMetrics& metrics,
  NML_t index, int iE, int jE, CEnv& env) {
  COMET_ASSERT(index+1 >= 0+1 && index < metrics.num_metrics_local);
  COMET_ASSERT(env.num_way() == NumWay::_2);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);

  bool result = true;

  if (!env.is_metric_type_bitwise()) {
    const auto metric_value = Metrics_elt_const<Float_t>(metrics, index, env);
    result = env.thresholds().is_pass(metric_value);
    return result;
  }

  // NOTE: If is_shrink, then it may not be possible on the CPU to check
  // whether passed threshold (specif., if thresholds.is_multi()).
  // So we will just mark as passed threshold, without the check.

  if (!env.is_shrink()) {
    if (env.is_threshold_tc()) {

      // NOTE: using MF::SINGLE if and only if is_threshold_tc()
      typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
      const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);

      // Check for zero since non-pass entries have already
      // been thresholded to zero in TC package.
      if (Thresholds::is_zero(TTable_t::get(ttable, iE, jE)))
        result = false;

    } else { // ! env.is_threshold_tc()

      // Convert fromMetricFormat::PACKED_DOUBLE to MetricFormat::SINGLE.
      // NOTE: for is_double / is_bitwise case, here we cast to float
      // before doing thresholding test. This is consistent throughout
      // (CHECK) so should work ok.

      typedef MetricFormatTraits<MetricFormat::SINGLE>::TypeIn TypeIn;
      typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
      TTable_t ttable = TTable_t::null();

      for (int iE_ = 0; iE_ < 2; ++iE_) {
        for (int jE_ = 0; jE_ < 2; ++jE_) {
          const auto metric
            = GMMetrics_get_2<Float_t>(metrics, index, iE_, jE_, env);
          const auto metric_single = static_cast<TypeIn>(metric);
          TTable_t::set(ttable, iE_, jE_, metric_single);
        } // for jE_
      } // for iE_

      if (!env.thresholds().is_pass(ttable, iE, jE))
        result = false;

    } // if (env.is_threshold_tc())
  } // if (!env.is_shrink())

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_METRICS_2WAY_ACCESSORS_I_HH_

//-----------------------------------------------------------------------------
