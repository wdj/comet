//-----------------------------------------------------------------------------
/*!
 * \file   metrics_3way_accessors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 3-way, accessor functions.
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

#ifndef _COMET_METRICS_3WAY_ACCESSORS_I_HH_
#define _COMET_METRICS_3WAY_ACCESSORS_I_HH_

#include "metrics_3way_indexing.i.hh"
#include "formulas.hh"

//=============================================================================

namespace comet {

//=============================================================================
// Accessors: value from (contig) index: derived

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 3-way CCC value.

static void GMMetrics_ccc_check_size_nofp_3(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (!BuildHas::INT128)
    return;

  if (env->metric_type() != MetricType::CCC || 
      env->num_way() != NumWay::_3 || ! env->are_ccc_params_default())
    return;

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 3+lm + 2+lm + 2+lm + 2+lm;

  COMET_INSIST_INTERFACE(env, lnum < 128 && "Number of fields too large.");
}

//-----------------------------------------------------------------------------
/// \brief Formula for CCC 3-way metric using 128 bit integer arithmetic.

static double GMMetrics_ccc_value_nofp_3(GMMetrics* metrics,
                                         const GMTally1 rijk,
                                         const GMTally1 si,
                                         const GMTally1 sj,
                                         const GMTally1 sk,
                                         const GMTally1 ci,
                                         const GMTally1 cj,
                                         const GMTally1 ck,
                                         const GMTally1 cijk,
                                         CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(BuildHas::INT128);

  typedef double Float_t;
  typedef BasicTypes::BigUInt UInt128;

  const UInt128 num = rijk * (UInt128)(3 * ci - 1 * si) *
                             (UInt128)(3 * cj - 1 * sj) *
                             (UInt128)(3 * ck - 1 * sk);

  const UInt128 denom = 6 * cijk * (UInt128)ci * (UInt128)cj
                                 * (UInt128)ck;

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 3+lm + 2+lm + 2+lm + 2+lm;

  const int shift = mantissa_digits<Float_t>() - 3; // Note num/denom <= 4.5 < 1<<3
                                                // always >= 0, < 128

  // Guarantee not to shift bits off the top.
  const int shift_limit = 128 - lnum; // always >= 0, < 128

  const int shift_left = utils::min(shift, shift_limit); // >= 0, < 128

  const int shift_right = shift - shift_left; // >= 0, < 128

  const Float_t result = ( (Float_t) ((num << shift_left) /
                                      (denom >> shift_right)) ) /
                         ( (Float_t)( ((size_t)1) << shift ) );

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Accessor for 3-way CCC metric computed with 128 bit int arithmetic.

static double GMMetrics_ccc_get_from_index_nofp_3(GMMetrics* metrics,
                                                  NML_t index,
                                                  int iE,
                                                  int jE,
                                                  int kE,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_metrics_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NumWay::_3);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);
  COMET_ASSERT(env->are_ccc_params_default());
  COMET_ASSERT(BuildHas::INT128);

  typedef double Float_t;

  const auto ttable = Metrics_elt_const<GMTally4x2>(*metrics, index, *env);
  const GMTally1 rijk = GMTally4x2_get(ttable, iE, jE, kE);

  const auto si1_sj1_sk1 = Metrics_elt_const<GMFloat3, MetricsArrayId::S>(*metrics,
    index, *env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  GMTally1 ci = 0, cj = 0, ck = 0, cijk = 0;

  if (env->sparse()) {
    const auto ci_cj_ck = Metrics_elt_const<GMFloat3, MetricsArrayId::C>(*metrics,
      index, *env);
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    cijk = GMTally4x2_get(ttable, 0, 0, 0) + GMTally4x2_get(ttable, 0, 0, 1) +
           GMTally4x2_get(ttable, 0, 1, 0) + GMTally4x2_get(ttable, 0, 1, 1) +
           GMTally4x2_get(ttable, 1, 0, 0) + GMTally4x2_get(ttable, 1, 0, 1) +
           GMTally4x2_get(ttable, 1, 1, 0) + GMTally4x2_get(ttable, 1, 1, 1);

    if (0 == ci || 0 == cj || 0 == ck || 0 == cijk) {
      return (Float_t)0;
    }
  } else {
    const int m = metrics->num_field_active;

    ci = m;
    cj = m;
    ck = m;

    cijk = 8 * m;
  }

  const GMTally1 si = iE == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = jE == 0 ? (2 * cj - sj1) : sj1;
  const GMTally1 sk = kE == 0 ? (2 * ck - sk1) : sk1;

  return GMMetrics_ccc_value_nofp_3(metrics, rijk, si, sj, sk, ci, cj, ck,
                                    cijk, env);
}

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 3-way CCC or DUO result, implementation.

template<int COUNTED_BITS_PER_ELT, typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_3_impl(GMMetrics& metrics,
  NML_t index, int iE, int jE, int kE, CEnv& env) {
  COMET_ASSERT(index < metrics.num_metric_items_local_allocated); // && index >= 0
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(env.is_metric_type_bitwise());
  COMET_ASSERT(env.counted_bits_per_elt() == COUNTED_BITS_PER_ELT);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);

  COMET_ASSERT(!env.is_shrink()); //FIX


//TODO: replace most of this with high-level ccc_duo_value call, if appropriate.


  if (env.is_threshold_tc()) {
    typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
    const auto table = Metrics_elt_const<TTable_t>(metrics, index, env);
    TTable_t::TypeIn result = TTable_t::get(table, iE, jE, kE);
    return (FloatResult_t)result;
  }





  typedef double Float_t; // Perform all calcs in double until return.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const Float_t f_one = 1;
  const Float_t recip_m = metrics.recip_m;

  const auto ttable = Metrics_elt_const<GMTally4x2>(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(ttable, iE, jE, kE);

  const auto si1_sj1_sk1 = Metrics_elt_const<GMFloat3, MetricsArrayId::S>(metrics,
    index, env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  Float_t result_floatcalc = 0;

  if (env.sparse()) {

    const auto ci_cj_ck = Metrics_elt_const<GMFloat3, MetricsArrayId::C>(metrics,
      index, env);
    GMTally1 ci, cj, ck;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);
    if (0 == ci || 0 == cj || 0 == ck)
      return (FloatResult_t)0;

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 sk0 = CBPE * ck - sk1;
    const GMTally1 si = iE == 0 ? si0 : si1;
    const GMTally1 sj = jE == 0 ? sj0 : sj1;
    const GMTally1 sk = kE == 0 ? sk0 : sk1;

    // TODO: it may be possible to decrease the number of divides
    // here - see GMMetrics_ccc_get_from_index_2.
    const Float_t recip_ci = f_one / ci;
    const Float_t recip_cj = f_one / cj;
    const Float_t recip_ck = f_one / ck;

    // NOTE: we are assuming here that xor gemm values have already been
    // adjusted to true gemm values.  Changing this in the future would
    // require computing and storing matX sum values somewhere and then using
    // them here.

    // NOTE: this is never called under conditions in which the xor fixup
    // must be done here - the CPU version does the xor fixup as part of
    // the core metrics calc process.

    const GMTally1 cijk =
        GMTally4x2_get(ttable, 0, 0, 0) + GMTally4x2_get(ttable, 0, 0, 1) +
        GMTally4x2_get(ttable, 0, 1, 0) + GMTally4x2_get(ttable, 0, 1, 1) +
        GMTally4x2_get(ttable, 1, 0, 0) + GMTally4x2_get(ttable, 1, 0, 1) +
        GMTally4x2_get(ttable, 1, 1, 0) + GMTally4x2_get(ttable, 1, 1, 1);
    if (0 == cijk)
      return (FloatResult_t)0;

    const Float_t recip_sumcijk = f_one / cijk;

    result_floatcalc = formulas::ccc_duo_value<CBPE>(
      rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk, env);

  } else { // !env.sparse

    COMET_ASSERT(!(env.is_using_xor() && !env.sparse())); // should never occur

    COMET_ASSERT(metrics.num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si =
      iE == 0 ? (CBPE * metrics.num_field_active - si1) : si1;
    const GMTally1 sj =
      jE == 0 ? (CBPE * metrics.num_field_active - sj1) : sj1;
    const GMTally1 sk =
      kE == 0 ? (CBPE * metrics.num_field_active - sk1) : sk1;

    const Float_t recip_sumcijk = (f_one / (CBPE*CBPE*CBPE)) * recip_m;

    result_floatcalc = formulas::ccc_duo_value<CBPE>(
      rijk, si, sj, sk, recip_m, recip_m, recip_m, recip_sumcijk, env);

  } // if (env.sparse())






  if (BuildHas::INT128 && env.metric_type() == MetricType::CCC &&
      env.are_ccc_params_default()) {
    const Float_t result_intcalc = GMMetrics_ccc_get_from_index_nofp_3(&metrics,
                                         index, iE, jE, kE, &env);

    // TODO: CHECK floating point type here
    const double eps = 1. / ( ((size_t)1) << (mantissa_digits<FloatResult_t>() - 6) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      fprintf(stderr,
              "Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      COMET_INSIST(diff < eps);
    }
  }

  return (FloatResult_t)result_floatcalc;
}

//-----------------------------------------------------------------------------
/// \brief Helper class for ccc_duo_get_from_index_3 template specialization.

// See https://stackoverflow.com/questions/12683165/partial-specialization-of-templates-with-integer-parameters

template<int CBPE, typename FloatResult_t>
struct MetricsIndexCCCDUO3Helper {
  static FloatResult_t impl(
  GMMetrics& metrics, NML_t index, int iE, int jE, int kE, CEnv& env) {
    return Metrics_ccc_duo_get_3_impl<
      CBPE, FloatResult_t>(metrics, index, iE, jE, kE, env);;
  }
};

template<typename FloatResult_t>
struct MetricsIndexCCCDUO3Helper<CBPE::NONE, FloatResult_t> {
  static FloatResult_t impl(
  GMMetrics& metrics, NML_t index, int iE, int jE, int kE, CEnv& env) {
    return env.counted_bits_per_elt() == CBPE::CCC ?
      MetricsIndexCCCDUO3Helper<CBPE::CCC, FloatResult_t>::impl( 
        metrics, index, iE, jE, kE, env) :
      MetricsIndexCCCDUO3Helper<CBPE::DUO, FloatResult_t>::impl(
        metrics, index, iE, jE, kE, env);
  }
};

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 3-way CCC or DUO result.

#if 1
template<int CBPE = CBPE::NONE, typename FloatResult_t = GMFloat>
static FloatResult_t Metrics_ccc_duo_get_3(
  GMMetrics& metrics, NML_t index, int iE, int jE, int kE, CEnv& env) {
    return MetricsIndexCCCDUO3Helper<CBPE, FloatResult_t>::impl( 
      metrics, index, iE, jE, kE, env);
}
#endif

//-----------------------------------------------------------------------------

template<int CBPE = CBPE::NONE, typename FloatResult_t = GMFloat>
static FloatResult_t Metrics_ccc_duo_get_3(
  GMMetrics& metrics, NML_t index, int entry_num, CEnv& env) {

  if (env.is_shrink()) {
    COMET_ASSERT(0 == entry_num);
    typedef MetricFormatTraits<MetricFormat::SINGLE>::TypeIn TypeIn;
    return (FloatResult_t)Metrics_elt_const<TypeIn>(metrics, index, env);
  }

  COMET_ASSERT(0 <= entry_num && entry_num < (1 << env.num_way()));
  const int iE = entry_num / 4;
  const int jE = (entry_num / 2) % 2;
  const int kE = entry_num % 2;

  return MetricsIndexCCCDUO3Helper<CBPE, FloatResult_t>::impl( 
    metrics, index, iE, jE, kE, env);
}

//-----------------------------------------------------------------------------

template<typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_3(
  GMMetrics& metrics, NML_t index, int entry_num, CEnv& env) {

  return Metrics_ccc_duo_get_3<CBPE::NONE, FloatResult_t>(metrics, index,
    entry_num, env);
}

//-----------------------------------------------------------------------------

template<typename FloatResult_t>
static FloatResult_t Metrics_ccc_duo_get_3(
  GMMetrics& metrics, NML_t index, int iE, int jE, int kE, CEnv& env) {

  return Metrics_ccc_duo_get_3<CBPE::NONE, FloatResult_t>(metrics, index,
    iE, jE, kE, env);
}

#if 0
//-----------------------------------------------------------------------------
/// \brief Templatized Check if any table value may exceed threshold.
/// NOTE: this function is currently unused.

template<int COUNTED_BITS_PER_ELT>
static bool Metrics_ccc_duo_threshold_detector_3(
  GMMetrics& metrics, const size_t index, CEnv& env) {
  COMET_ASSERT(index < metrics.num_metrics_local); // && index >= 0
  COMET_ASSERT(env.num_way() == NumWay::_3);

// FIXTHRESHOLD - here and below

  // if no active threshold, then always pass threshold criterion.
  if (!env.thresholds_.is_active()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    return true;
  }

  // if is_shrink, assume a threshold pass may exist, don't take time to check.
  if (env.is_shrink()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    return true;
  }

  // If threshold_tc, then look for any prethresholded values not set to zero.
  if (env.is_threshold_tc()) {
    COMET_INSIST(!env.is_using_threshold_detector());
    typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
    const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);
    for (int iE = 0; iE < 2; ++iE) {
      for (int jE = 0; jE < 2; ++jE) {
        for (int kE = 0; kE < 2; ++kE) {
          if (TTable_t::get(ttable, iE, jE, kE) != (TTable_t::TypeIn)0)
            return true;
        }
      }
    }
    return false;
  }

  COMET_INSIST(env.is_using_threshold_detector());

  // this is here because xor stuff not implemented below.
  COMET_ASSERT(!env.is_using_xor()); // should never occur.

  typedef double Float_t; // Perform all calcs in double.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  if (env.sparse()) {

    const Float_t f_one = 1;

    const auto ttable = Metrics_elt_const<GMTally4x2>(metrics, index, env);
    const GMTally1 rijk000 = GMTally4x2_get(ttable, 0, 0, 0);
    const GMTally1 rijk001 = GMTally4x2_get(ttable, 0, 0, 1);
    const GMTally1 rijk010 = GMTally4x2_get(ttable, 0, 1, 0);
    const GMTally1 rijk011 = GMTally4x2_get(ttable, 0, 1, 1);
    const GMTally1 rijk100 = GMTally4x2_get(ttable, 1, 0, 0);
    const GMTally1 rijk101 = GMTally4x2_get(ttable, 1, 0, 1);
    const GMTally1 rijk110 = GMTally4x2_get(ttable, 1, 1, 0);
    const GMTally1 rijk111 = GMTally4x2_get(ttable, 1, 1, 1);

    const auto si1_sj1_sk1 = Metrics_elt_const<GMFloat3, MetricsArrayId::S>(metrics,
      index, env);
    GMTally1 si1 = 0, sj1 = 0, sk1 = 0;
    GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

    const auto ci_cj_ck = Metrics_elt_const<GMFloat3, MetricsArrayId::C>(metrics,
      index, env);
    GMTally1 ci = 0, cj = 0, ck = 0;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    GMTally1 cijk = rijk000 + rijk001 + rijk010 + rijk011 +
                    rijk100 + rijk101 + rijk110 + rijk111;
    if (ci == 0 || cj == 0 || ck == 0 || cijk == 0) {
      return env.pass_threshold(0);
    }

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits

    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 sk0 = CBPE * ck - sk1;

    // TODO: optimize this further

    const Float_t recip_ci = f_one / ci;
    const Float_t recip_cj = f_one / cj;
    const Float_t recip_ck = f_one / ck;

    const Float_t recip_sumcijk = f_one / cijk;

    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk000,
      si0, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 0, 0, 0, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk001,
      si0, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 0, 0, 1, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk010,
      si0, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 0, 1, 0, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk011,
      si0, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 0, 1, 1, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk100,
      si1, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 1, 0, 0, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk101,
      si1, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 1, 0, 1, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk110,
      si1, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 1, 1, 0, env)));
    COMET_ASSERT((formulas::ccc_duo_value<CBPE>(rijk111,
      si1, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ==
      (Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index, 1, 1, 1, env)));

    return env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk000,
           si0, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk001,
           si0, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk010,
           si0, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk011,
           si0, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk100,
           si1, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk101,
           si1, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk110,
           si1, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env)) ||
           env.pass_threshold(formulas::ccc_duo_value<CBPE>(rijk111,
           si1, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env));

  } // if (env.sparse())

  // Non-sparse case (less well-optimized).

  const Float_t v000 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            0, 0, 0, env);
  const Float_t v001 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            0, 0, 1, env);
  const Float_t v010 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            0, 1, 0, env);
  const Float_t v011 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            0, 1, 1, env);
  const Float_t v100 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            1, 0, 0, env);
  const Float_t v101 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            1, 0, 1, env);
  const Float_t v110 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            1, 1, 0, env);
  const Float_t v111 = Metrics_ccc_duo_get_3<CBPE, Float_t>(metrics, index,
                                                            1, 1, 1, env);

  return env.pass_threshold(v000) || env.pass_threshold(v001) ||
         env.pass_threshold(v010) || env.pass_threshold(v011) ||
         env.pass_threshold(v100) || env.pass_threshold(v101) ||
         env.pass_threshold(v110) || env.pass_threshold(v111);
}
#endif

//=============================================================================
// Accessors: value from (local) coord: set: 3-way.

//-----------------------------------------------------------------------------

template<typename T, int MA = MetricsArrayId::M>
static T& Metrics_elt_3(GMMetrics& metrics, int i, int j, int k,
  int j_block, int k_block, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_3(metrics, i, j, k, j_block, k_block, env);

  return Metrics_elt<T, MA>(metrics, index, env);
}

//-----------------------------------------------------------------------------

template<typename T, int MA = MetricsArrayId::M>
static T& Metrics_elt_3(GMMetrics& metrics, int i, int j, int k,
  int j_block, int k_block, MetricsIndexCache& index_cache, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_3(metrics, i, j, k, j_block, k_block,
    index_cache, env);

  return Metrics_elt<T, MA>(metrics, index, env);
}

//=============================================================================
// Accessors: value from (local) coord: get: 3-way.

//-----------------------------------------------------------------------------

template<typename T>
static T Metrics_elt_const_3(GMMetrics& metrics, int i, int j, int k,
  int j_block, int k_block, MetricsIndexCache& index_cache, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || (env.proc_num_vector() == j_block &&
                                 env.proc_num_vector() == k_block));
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_3(metrics, i, j, k,
    j_block, k_block, index_cache, env);

  return Metrics_elt_const<T>(metrics, index, env);
}

//=============================================================================
// Accessors: value from index: get: 3-way.

template<typename FloatResult_t>
FloatResult_t GMMetrics_get_3(GMMetrics& metrics,
  NML_t index, int iE, int jE, int kE, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(index >= 0 && index < metrics.num_metrics_local);
  COMET_ASSERT(iE >= 0 && iE < env.ijkE_max());
  COMET_ASSERT(jE >= 0 && jE < env.ijkE_max());
  COMET_ASSERT(kE >= 0 && kE < env.ijkE_max());

  COMET_ASSERT(!env.is_shrink()); //FIX
  const auto result =
    env.metric_type() == MetricType::CZEK ?
      Metrics_elt_const<FloatResult_t>(metrics, index, env) :
      Metrics_ccc_duo_get_3<FloatResult_t>(metrics, index, iE, jE, kE, env);

  return result;
}

//=============================================================================
// Accessors: value from (global) coord: get: 3-way.

template<typename FloatResult_t>
FloatResult_t GMMetrics_get_3(GMMetrics& metrics,
  NV_t iG, NV_t jG, NV_t kG, int iE, int jE, int kE, CEnv& env) {
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(iE >= 0 && iE < env.ijkE_max());
  COMET_ASSERT(jE >= 0 && jE < env.ijkE_max());
  COMET_ASSERT(kE >= 0 && kE < env.ijkE_max());
  COMET_ASSERT(!env.is_shrink());
  // WARNING: these conditions are not exhaustive.

  const NML_t index = Metrics_index_3(metrics, iG, jG, kG, env);

  const auto result = GMMetrics_get_3<FloatResult_t>(metrics, index, iE, jE, kE, env);

  return result;
}

//=============================================================================
// Accessors: determine whether pass threshold: 3-way.

template<typename Float_t = GMFloat>
bool Metrics_is_pass_threshold_noshrink(GMMetrics& metrics,
  NML_t index, int iE, int jE, int kE, CEnv& env) {
  COMET_ASSERT(index < metrics.num_metrics_local); // && index >= 0
  COMET_ASSERT(env.num_way() == NumWay::_3);
  COMET_ASSERT(iE >= 0 && iE < 2);
  COMET_ASSERT(jE >= 0 && jE < 2);
  COMET_ASSERT(kE >= 0 && kE < 2);

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
      typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
      const auto ttable = Metrics_elt_const<TTable_t>(metrics, index, env);

      // Check for zero since non-pass entries have already
      // been thresholded to zero in TC package.
      if (Thresholds::is_zero(TTable_t::get(ttable, iE, jE, kE)))
        result = false;

    } else { // ! env.is_threshold_tc()

      // Convert fromMetricFormat::PACKED_DOUBLE to MetricFormat::SINGLE.
      // NOTE: for is_double / is_bitwise case, here we cast to float
      // before doing thresholding test. This is consistent throughout
      // (CHECK) so should work ok.

      typedef MetricFormatTraits<MetricFormat::SINGLE>::TypeIn TypeIn;
      typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
      TTable_t ttable = TTable_t::null();

      for (int iE_ = 0; iE_ < 2; ++iE_) {
        for (int jE_ = 0; jE_ < 2; ++jE_) {
          for (int kE_ = 0; kE_ < 2; ++kE_) {
            const auto metric
              = GMMetrics_get_3<Float_t>(metrics, index, iE_, jE_, kE_, env);
            const auto metric_single = static_cast<TypeIn>(metric);
            TTable_t::set(ttable, iE_, jE_, kE_, metric_single);
          }
        }
      }

      if (!env.thresholds().is_pass(ttable, iE, jE, kE))
        result = false;

    } // if (env.is_threshold_tc())
  } // if (!envs_shrink())

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_METRICS_3WAY_ACCESSORS_I_HH_

//-----------------------------------------------------------------------------
