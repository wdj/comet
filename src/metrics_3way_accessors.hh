//-----------------------------------------------------------------------------
/*!
 * \file   metrics_3way_accessors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 3-way, accessor functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_3way_accessors_hh_
#define _comet_metrics_3way_accessors_hh_

#include "metrics_3way_indexing.hh"
#include "formulas.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Accessors: value from (contig) index: basic---*/

static GMFloat3 GMMetrics_float3_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  COMET_ASSERT(metrics->data_S);

  return ((GMFloat3*)(metrics->data_S))[index];
}

//-----------------------------------------------------------------------------

static GMFloat3 GMMetrics_float3_C_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  COMET_ASSERT(metrics->data_C);

  return ((GMFloat3*)(metrics->data_C))[index];
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_from_index(GMMetrics* metrics,
                                                    size_t index,
                                                    CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);

  return ((GMTally4x2*)(metrics->data))[index];
}

//=============================================================================
/*---Accessors: value from (contig) index: derived---*/

#if 0
//-----------------------------------------------------------------------------
/// \brief Formula for a single 3-way CCC or DUO result.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_value_3_impl_(
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const GMFloat recip_ci,
  const GMFloat recip_cj,
  const GMFloat recip_ck,
  const GMFloat recip_sumcijk,
  const GMFloat multiplier,
  const GMFloat param) {

  const GMFloat f_one = 1;

  const GMFloat fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const GMFloat fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;
  const GMFloat fk = (f_one / COUNTED_BITS_PER_ELT) * recip_ck * sk;

  const GMFloat fijk = recip_sumcijk * rijk;

  /*---Do the following to make floating point arithmetic order-independent---*/
  GMFloat fmin = 0, fmid = 0, fmax = 0;
  utils::sort_3(fmin, fmid, fmax, fi, fj, fk);

  /* clang-format off */
  const GMFloat result = multiplier * fijk * (f_one - param * fmin) *
                                             (f_one - param * fmid) *
                                             (f_one - param * fmax);
  /* clang-format on */

  return result;
}
#endif

//-----------------------------------------------------------------------------
/// \brief Templatized access to the CCC or DUO formula.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_value_3(
  GMMetrics* metrics,
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const GMFloat recip_ci,
  const GMFloat recip_cj,
  const GMFloat recip_ck,
  const GMFloat recip_sumcijk,
  CEnv* env) {
  COMET_ASSERT(metrics && env);

//  return GMMetrics_ccc_duo_value_3_impl_<COUNTED_BITS_PER_ELT>(
//    rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk,
//    env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env),
//    env->ccc_param());

  return ccc_duo_value_3<GMFloat, COUNTED_BITS_PER_ELT>(
    rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk,
    env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env),
    env->ccc_param());
}

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 3-way CCC value.

static void GMMetrics_ccc_check_size_nofp_3(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (env->metric_type() != MetricType::CCC || 
      env->num_way() != NUM_WAY::_3 || ! env->are_ccc_params_default()) {
    return;
  }

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 3+lm + 2+lm + 2+lm + 2+lm;

  COMET_INSIST_INTERFACE(env, lnum < 128 && "Number of fields too large.");
}

//-----------------------------------------------------------------------------

#ifdef COMET_USE_INT128

//-----------------------------------------------------------------------------
/// \brief Formula for CCC 3-way metric using 128 bit integer arithmetic.

static GMFloat GMMetrics_ccc_value_nofp_3(GMMetrics* metrics,
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

  const GMUInt128 num = rijk * (GMUInt128)(3 * ci - 1 * si) *
                               (GMUInt128)(3 * cj - 1 * sj) *
                               (GMUInt128)(3 * ck - 1 * sk);

  const GMUInt128 denom = 6 * cijk * (GMUInt128)ci * (GMUInt128)cj
                                   * (GMUInt128)ck;

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 3+lm + 2+lm + 2+lm + 2+lm;

  const int shift = mantissa_digits<GMFloat>() - 3; // Note num/denom <= 4.5 < 1<<3
                                                // always >= 0, < 128

  // Guarantee not to shift bits off the top.
  const int shift_limit = 128 - lnum; // always >= 0, < 128

  const int shift_left = utils::min(shift, shift_limit); // >= 0, < 128

  const int shift_right = shift - shift_left; // >= 0, < 128

  const GMFloat result = ( (GMFloat) ((num << shift_left) /
                                      (denom >> shift_right)) ) /
                         ( (GMFloat)( ((size_t)1) << shift ) );

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Accessor for 3-way CCC metric computed with 128 bit int arithmetic.

static GMFloat GMMetrics_ccc_get_from_index_nofp_3(GMMetrics* metrics,
                                                   size_t index,
                                                   int i0,
                                                   int i1,
                                                   int i2,
                                                   CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);
  COMET_ASSERT(env->are_ccc_params_default());

  const GMTally4x2 t42 = GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(t42, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  GMTally1 ci = 0, cj = 0, ck = 0, cijk = 0;

  if (env->sparse()) {
    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    cijk = GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
           GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
           GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
           GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1);

    if (0 == ci || 0 == cj || 0 == ck || 0 == cijk) {
      return (GMFloat)0;
    }
  } else {
    const int m = metrics->num_field_active;

    ci = m;
    cj = m;
    ck = m;

    cijk = 8 * m;
  }

  const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;
  const GMTally1 sk = i2 == 0 ? (2 * ck - sk1) : sk1;

  return GMMetrics_ccc_value_nofp_3(metrics, rijk, si, sj, sk, ci, cj, ck,
                                    cijk, env);
}

#endif

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 3-way CCC or DUO result.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_get_from_index_3(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  int i2,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const GMFloat f_one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally4x2 t42 = GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(t42, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  GMFloat result_floatcalc = 0;

  if (env->sparse()) {

    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj, ck;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    // TODO: it may be possible to decrease the number of divides
    // here - see GMMetrics_ccc_get_from_index_2.
    const GMFloat recip_ci = f_one / ci;
    const GMFloat recip_cj = f_one / cj;
    const GMFloat recip_ck = f_one / ck;

    GMTally1 cijk =
           GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
           GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
           GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
           GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1);
    if (0 == ci || 0 == cj || 0 == ck || 0 == cijk)
      return (GMFloat)0;

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si = i0 == 0 ? (CBPE * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (CBPE * cj - sj1) : sj1;
    const GMTally1 sk = i2 == 0 ? (CBPE * ck - sk1) : sk1;

    const GMFloat recip_sumcijk = f_one / cijk;
//      f_one / (GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
//               GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
//               GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
//               GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1));

    result_floatcalc = GMMetrics_ccc_duo_value_3<COUNTED_BITS_PER_ELT>(metrics,
      rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk, env);

  } else { // !env->sparse

    COMET_ASSERT(metrics->num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si =
      i0 == 0 ? (CBPE * metrics->num_field_active - si1) : si1;
    const GMTally1 sj =
      i1 == 0 ? (CBPE * metrics->num_field_active - sj1) : sj1;
    const GMTally1 sk =
      i2 == 0 ? (CBPE * metrics->num_field_active - sk1) : sk1;

    const GMFloat recip_sumcijk = (f_one / (CBPE*CBPE*CBPE)) * recip_m;

    result_floatcalc = GMMetrics_ccc_duo_value_3<CBPE>(metrics,
      rijk, si, sj, sk, recip_m, recip_m, recip_m, recip_sumcijk, env);

  } // if (env->sparse())

#ifdef COMET_USE_INT128
  if (env->metric_type() == MetricType::CCC && env->are_ccc_params_default()) {
    const GMFloat result_intcalc = GMMetrics_ccc_get_from_index_nofp_3(metrics,
                                         index, i0, i1, i2, env);

    const double eps = 1. / ( ((size_t)1) << (mantissa_digits<GMFloat>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      printf("Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      COMET_INSIST(diff < eps);
    }
  }
#endif
  return result_floatcalc;
}

//-----------------------------------------------------------------------------
/// \brief Accessor for 3-way CCC result.

static GMFloat GMMetrics_ccc_get_from_index_3(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  int i2,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  enum { COUNTED_BITS_PER_ELT = 2 };

  return GMMetrics_ccc_duo_get_from_index_3<COUNTED_BITS_PER_ELT>(metrics,
    index, i0, i1, i2, env);
}

//-----------------------------------------------------------------------------
/// \brief Accessor for 3-way DUO result.

static GMFloat GMMetrics_duo_get_from_index_3(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  int i2,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  enum { COUNTED_BITS_PER_ELT = 1 };

  return GMMetrics_ccc_duo_get_from_index_3<COUNTED_BITS_PER_ELT>(metrics,
    index, i0, i1, i2, env);
}

//-----------------------------------------------------------------------------
/// \brief Templatized Check if any table value may exceed threshold.

template<int COUNTED_BITS_PER_ELT>
static bool GMMetrics_ccc_duo_get_from_index_3_threshold(
  GMMetrics* metrics,
  const size_t index,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const double threshold_eff = env->threshold_eff();

  if (env->sparse()) {

    const GMFloat f_one = 1;

    const GMTally4x2 t42 = GMMetrics_tally4x2_get_from_index(metrics, index,
                                                             env);
    const GMTally1 rijk000 = GMTally4x2_get(t42, 0, 0, 0);
    const GMTally1 rijk001 = GMTally4x2_get(t42, 0, 0, 1);
    const GMTally1 rijk010 = GMTally4x2_get(t42, 0, 1, 0);
    const GMTally1 rijk011 = GMTally4x2_get(t42, 0, 1, 1);
    const GMTally1 rijk100 = GMTally4x2_get(t42, 1, 0, 0);
    const GMTally1 rijk101 = GMTally4x2_get(t42, 1, 0, 1);
    const GMTally1 rijk110 = GMTally4x2_get(t42, 1, 1, 0);
    const GMTally1 rijk111 = GMTally4x2_get(t42, 1, 1, 1);

    const GMFloat3 si1_sj1_sk1 =
        GMMetrics_float3_S_get_from_index(metrics, index, env);
    GMTally1 si1 = 0, sj1 = 0, sk1 = 0;
    GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMTally1 ci = 0, cj = 0, ck = 0;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    GMTally1 cijk = rijk000 + rijk001 + rijk010 + rijk011 +
                    rijk100 + rijk101 + rijk110 + rijk111;
    if (ci == 0 || cj == 0 || ck == 0 || cijk == 0) {
      return 0 > threshold_eff;
    }

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits

    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;
    const GMTally1 sk0 = CBPE * ck - sk1;

    // TODO: optimize this further

    const GMFloat recip_ci = f_one / ci;
    const GMFloat recip_cj = f_one / cj;
    const GMFloat recip_ck = f_one / ck;

    const GMFloat recip_sumcijk = f_one / cijk;

    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk000,
      si0, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 0, 0, 0, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk001,
      si0, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 0, 0, 1, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk010,
      si0, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 0, 1, 0, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk011,
      si0, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 0, 1, 1, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk100,
      si1, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 1, 0, 0, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk101,
      si1, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 1, 0, 1, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk110,
      si1, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 1, 1, 0, env));
    COMET_ASSERT(GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk111,
      si1, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
      GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index, 1, 1, 1, env));

    return GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk000, si0, sj0, sk0,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk001, si0, sj0, sk1,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk010, si0, sj1, sk0,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk011, si0, sj1, sk1,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk100, si1, sj0, sk0,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk101, si1, sj0, sk1,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk110, si1, sj1, sk0,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff ||
          GMMetrics_ccc_duo_value_3<CBPE>(metrics, rijk111, si1, sj1, sk1,
            recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold_eff;

  } // if (env->sparse())

  // Non-sparse case (less well-optimized).

  const GMFloat v000 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                0, 0, 0, env);
  const GMFloat v001 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                0, 0, 1, env);
  const GMFloat v010 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                0, 1, 0, env);
  const GMFloat v011 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                0, 1, 1, env);
  const GMFloat v100 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                1, 0, 0, env);
  const GMFloat v101 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                1, 0, 1, env);
  const GMFloat v110 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                1, 1, 0, env);
  const GMFloat v111 = GMMetrics_ccc_duo_get_from_index_3<CBPE>(metrics, index,
                                                                1, 1, 1, env);

  return v000 > threshold_eff || v001 > threshold_eff ||
         v010 > threshold_eff || v011 > threshold_eff ||
         v100 > threshold_eff || v101 > threshold_eff ||
         v110 > threshold_eff || v111 > threshold_eff;
}

//-----------------------------------------------------------------------------
/// \brief Check if any 3-way CCC table value may exceed threshold.

static bool GMMetrics_ccc_get_from_index_3_threshold(
  GMMetrics* metrics, const size_t index, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  enum { COUNTED_BITS_PER_ELT = 2 };

  return GMMetrics_ccc_duo_get_from_index_3_threshold<COUNTED_BITS_PER_ELT>(
    metrics, index, env);
}

//-----------------------------------------------------------------------------
/// \brief Check if any 3-way DUO table value may exceed threshold.

static bool GMMetrics_duo_get_from_index_3_threshold(
  GMMetrics* metrics, const size_t index, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  enum { COUNTED_BITS_PER_ELT = 1 };

  return GMMetrics_ccc_duo_get_from_index_3_threshold<COUNTED_BITS_PER_ELT>(
    metrics, index, env);
}

//=============================================================================
/*---Accessors: value from (local) coord: set: 3-way---*/

template<typename T>
static void GMMetrics_set_3(GMMetrics* metrics, void* p, int i, int j, int k,
  T value, CEnv* env) {    
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < k);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(!env->all2all());

  const size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((T*)p)[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_3(GMMetrics* metrics, int i, int j, int k,
  GMFloat value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  GMMetrics_set_3<GMFloat>(metrics, metrics->data, i, j, k, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_3(GMMetrics* metrics, int i, int j, int k,
  GMFloat3 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_3<GMFloat3>(metrics, metrics->data_S, i, j, k, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_3(GMMetrics* metrics, int i, int j, int k,
  GMFloat3 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_3<GMFloat3>(metrics, metrics->data_C, i, j, k, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_3(GMMetrics* metrics, int i, int j, int k,
  GMTally4x2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_3<GMTally4x2>(metrics, metrics->data, i, j, k, value, env);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T>
static void GMMetrics_set_all2all_3(GMMetrics* metrics, void* p,
  int i, int j, int k, int j_block, int k_block, T value, CEnv* env) {
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  // WARNING: these conditions are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k,
    j_block, k_block, env);
  ((T*)p)[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_3(GMMetrics* metrics,
  int i, int j, int k, int j_block, int k_block, GMFloat value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  GMMetrics_set_all2all_3<GMFloat>(metrics, metrics->data, i, j, k,
    j_block, k_block, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_all2all_3(GMMetrics* metrics,
  int i, int j, int k, int j_block, int k_block, GMFloat3 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3<GMFloat3>(metrics, metrics->data_S, i, j, k,
    j_block, k_block, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_all2all_3(GMMetrics* metrics,
  int i, int j, int k, int j_block, int k_block, GMFloat3 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3<GMFloat3>(metrics, metrics->data_C, i, j, k,
    j_block, k_block, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_all2all_3(GMMetrics* metrics,
  int i, int j, int k, int j_block, int k_block, GMTally4x2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3<GMTally4x2>(metrics, metrics->data, i, j, k,
    j_block, k_block, value, env);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T>
static void GMMetrics_set_all2all_3_permuted_cache(GMMetrics* metrics, void* p,
  int I, int J, int K, int j_block, int k_block, T value,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  // WARNING: these conditions are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((T*)p)[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_3_permuted_cache(GMMetrics* metrics,
  int I, int J, int K, int j_block, int k_block, GMFloat value,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  GMMetrics_set_all2all_3_permuted_cache<GMFloat>(metrics, metrics->data,
    I, J, K, j_block, k_block, value, index_cache, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_all2all_3_permuted_cache(GMMetrics* metrics,
  int I, int J, int K, int j_block, int k_block, GMFloat3 value,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3_permuted_cache<GMFloat3>(metrics, metrics->data_S,
    I, J, K, j_block, k_block, value, index_cache, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_all2all_3_permuted_cache(GMMetrics* metrics,
  int I, int J, int K, int j_block, int k_block, GMFloat3 value,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3_permuted_cache<GMFloat3>(metrics, metrics->data_C,
    I, J, K, j_block, k_block, value, index_cache, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_all2all_3_permuted_cache(GMMetrics* metrics,
  int I, int J, int K, int j_block, int k_block, GMTally4x2 value,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  GMMetrics_set_all2all_3_permuted_cache<GMTally4x2>(metrics, metrics->data,
    I, J, K, j_block, k_block, value, index_cache, env);
}

//=============================================================================
/*---Accessors: value from (local) coord: get: 3-way---*/

static GMFloat GMMetrics_float_get_3(GMMetrics* metrics,
  int i, int j, int k, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  const size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_3(GMMetrics* metrics,
  int i, int j, int k, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  const size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3( GMMetrics* metrics,
  int i, int j, int k, int j_block, int k_block, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  // WARNING: these conditions are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k,
    j_block, k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3_permuted_cache(
  GMMetrics* metrics, int I, int J, int K, int j_block, int k_block,
  GMIndexCache* index_cache, CEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  // WARNING: these conditions are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_3way_accessors_hh_

//-----------------------------------------------------------------------------
