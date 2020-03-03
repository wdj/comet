//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_accessors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, accessor functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_metrics_2way_accessors_hh_
#define _comet_metrics_2way_accessors_hh_

#include "cstdio"

#include "metrics_2way_indexing.hh"
#include "formulas.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Accessors: value from (contig) index: basic---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              size_t index,
                                              CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(metrics->data))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);
  COMET_ASSERT(metrics->data_S);

  return ((GMFloat2*)(metrics->data_S))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_C_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);
  COMET_ASSERT(metrics->data_C);

  return ((GMFloat2*)(metrics->data_C))[index];
}

//-----------------------------------------------------------------------------

static GMTally2x2 GMMetrics_tally2x2_get_from_index(GMMetrics* metrics,
                                                    size_t index,
                                                    CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  return ((GMTally2x2*)(metrics->data))[index];
}

//=============================================================================
/*---Accessors: value from (contig) index: derived---*/

static GMFloat GMMetrics_czek_get_from_index(GMMetrics* metrics,
                                             size_t index,
                                             CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);

  return GMMetrics_float_get_from_index(metrics, index, env);
}

#if 0
//-----------------------------------------------------------------------------
/// \brief Formula for a single 2-way CCC or DUO result.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_value_2_impl_(
  const GMTally1 rij,
  const GMTally1 si,
  const GMTally1 sj,
  const GMFloat recip_ci,
  const GMFloat recip_cj,
  const GMFloat recip_sumcij,
  const GMFloat multiplier,
  const GMFloat param) {

  const GMFloat f_one = 1;

  const GMFloat fi = (f_one / COUNTED_BITS_PER_ELT) * recip_ci * si;
  const GMFloat fj = (f_one / COUNTED_BITS_PER_ELT) * recip_cj * sj;

  const GMFloat fij = recip_sumcij * rij;

  /*---Do the following to make floating point arithmetic order-independent---*/
  const GMFloat fmin = fi < fj ? fi : fj;
  const GMFloat fmax = fi < fj ? fj : fi;

  /* clang-format off */
  const GMFloat result = multiplier * fij * (f_one - param * fmin) *
                                            (f_one - param * fmax);
  /* clang-format on */

  return result;
}
#endif

//-----------------------------------------------------------------------------
/// \brief Templatized access to the CCC or DUO formula.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_value_2(
  GMMetrics* metrics,
  const GMTally1 rij,
  const GMTally1 si,
  const GMTally1 sj,
  const GMFloat recip_ci,
  const GMFloat recip_cj,
  const GMFloat recip_sumcij,
  CEnv* env) {
  COMET_ASSERT(metrics && env);

//  return GMMetrics_ccc_duo_value_2_impl_<COUNTED_BITS_PER_ELT>(
//    rij, si, sj, recip_ci, recip_cj, recip_sumcij,
//    env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env),
//    env->ccc_param());

  return ccc_duo_value_2<GMFloat, COUNTED_BITS_PER_ELT>(
    rij, si, sj, recip_ci, recip_cj, recip_sumcij,
    env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env),
    env->ccc_param());
}

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 2-way CCC value.

static void GMMetrics_ccc_check_size_nofp_2(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (env->metric_type() != MetricType::CCC ||
      env->num_way() != NUM_WAY::_2 || ! env->are_ccc_params_default()) {
    return;
  }

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

#ifdef COMET_USE_INT128

//-----------------------------------------------------------------------------
/// \brief Formula for CCC 2-way metric using 128 bit integer arithmetic.

static GMFloat GMMetrics_ccc_value_nofp_2(GMMetrics* metrics,
                                          const GMTally1 rij,
                                          const GMTally1 si,
                                          const GMTally1 sj,
                                          const GMTally1 ci,
                                          const GMTally1 cj,
                                          const GMTally1 cij,
                                          CEnv* env) {
  COMET_ASSERT(metrics && env);

  const GMUInt128 num = rij * (GMUInt128)(3 * ci - 1 * si) *
                              (GMUInt128)(3 * cj - 1 * sj);

  const GMUInt128 denom = 2 * cij * (GMUInt128)ci * (GMUInt128)cj;

  const int shift = mantissa_digits<GMFloat>() - 3; // Note num/denom <= 4.5 < 1<<3

  // This should be an integer with no more than
  // mantissa_digits<GMFloat>() binary digits
  const GMUInt128 int_ratio = (num << shift) / denom;

  // Convert to floting point and then adjust exponent.
  const GMFloat result = ( (GMFloat) int_ratio ) /
                         ( (GMFloat) ( ((size_t)1) << shift ) );

  return result;
}


//-----------------------------------------------------------------------------
/// \brief Accessor for 2-way CCC metric computed with 128 bit int arithmetic.

static GMFloat GMMetrics_ccc_get_from_index_nofp_2(GMMetrics* metrics,
                                                   size_t index,
                                                   int i0,
                                                   int i1,
                                                   CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(env->are_ccc_params_default());

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMTally1 ci = 0, cj = 0, cij = 0;

  if (env->sparse()) {
    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMFloat2_decode(&ci, &cj, ci_cj);

    cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
          GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);

    if (0 == ci || 0 == cj || 0 == cij) {
      return (GMFloat)0;
    }
  } else {
    const int m = metrics->num_field_active;

    ci = m;
    cj = m;

    cij = 4 * m;
  }

  const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;

  return GMMetrics_ccc_value_nofp_2(metrics, rij, si, sj, ci, cj, cij, env);
}

#endif

//-----------------------------------------------------------------------------
/// \brief Templatized accessor for 2-way CCC or DUO result.

template<int COUNTED_BITS_PER_ELT>
static GMFloat GMMetrics_ccc_duo_get_from_index_2(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const GMFloat f_one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMFloat result_floatcalc = 0;

  if (env->sparse()) {

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
    if (0 == ci || 0 == cj || 0 == cij)
      return (GMFloat)0;

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si = i0 == 0 ? (CBPE * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (CBPE * cj - sj1) : sj1;

    // The following complex code is to reduce the number of divides.
    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;

    const GMFloat f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const GMFloat f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const GMFloat f_cij = (GMFloat) cij;
    const GMFloat recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

    const GMFloat recip_ci = f_cj * f_cij * recip_cicjcij;
    const GMFloat recip_cj = f_ci * f_cij * recip_cicjcij;

    const GMFloat recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

    result_floatcalc = GMMetrics_ccc_duo_value_2<COUNTED_BITS_PER_ELT>(
      metrics, rij, si, sj, recip_ci, recip_cj, recip_sumcij, env);

  } else { // !env->sparse

    COMET_ASSERT(metrics->num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si =
      i0 == 0 ? (CBPE * metrics->num_field_active - si1) : si1;
    const GMTally1 sj =
      i1 == 0 ? (CBPE * metrics->num_field_active - sj1) : sj1;

    const GMFloat recip_sumcij = (f_one / (CBPE*CBPE)) * recip_m;

    result_floatcalc = GMMetrics_ccc_duo_value_2<COUNTED_BITS_PER_ELT>(
      metrics, rij, si, sj, recip_m, recip_m, recip_sumcij, env);

  } // if (env->sparse())

#ifdef COMET_USE_INT128
  if (env->metric_type() == MetricType::CCC && env->are_ccc_params_default()) {
    const GMFloat result_intcalc = GMMetrics_ccc_get_from_index_nofp_2(metrics,
                                         index, i0, i1, env);

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
/// \brief Accessor for 2-way CCC result.

static GMFloat GMMetrics_ccc_get_from_index_2(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  enum { COUNTED_BITS_PER_ELT = 2 };

  return GMMetrics_ccc_duo_get_from_index_2<COUNTED_BITS_PER_ELT>(
    metrics, index, i0, i1, env);
}

//-----------------------------------------------------------------------------
/// \brief Accessor for 2-way DUO result.

static GMFloat GMMetrics_duo_get_from_index_2(
  GMMetrics* metrics,
  size_t index,
  int i0,
  int i1,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  enum { COUNTED_BITS_PER_ELT = 1 };

  return GMMetrics_ccc_duo_get_from_index_2<COUNTED_BITS_PER_ELT>(
    metrics, index, i0, i1, env);
}

//-----------------------------------------------------------------------------
/// \brief Templatized Check if any table value may exceed threshold.

template<int COUNTED_BITS_PER_ELT>
static bool GMMetrics_ccc_duo_get_from_index_2_threshold(
  GMMetrics* metrics,
  const size_t index,
  GMFloat threshold,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  enum {CBPE = COUNTED_BITS_PER_ELT};

  if (env->sparse()) {

    const GMFloat f_one = 1;

    const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index,
                                                             env);
    const GMTally1 rij00 = GMTally2x2_get(t22, 0, 0);
    const GMTally1 rij01 = GMTally2x2_get(t22, 0, 1);
    const GMTally1 rij10 = GMTally2x2_get(t22, 1, 0);
    const GMTally1 rij11 = GMTally2x2_get(t22, 1, 1);

    const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
    GMTally1 si1 = 0, sj1 = 0;
    GMFloat2_decode(&si1, &sj1, si1_sj1);

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci = 0, cj = 0;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = rij00 + rij01 + rij10 + rij11;
    if (ci == 0 || cj == 0 || cij == 0) {
      return 0 > threshold;
    }

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits

    const GMTally1 si0 = CBPE * ci - si1;
    const GMTally1 sj0 = CBPE * cj - sj1;

    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;
    const GMFloat f_cij = (GMFloat) cij;

    const GMFloat f_cicj = f_ci * f_cj;

    const GMFloat recip_2 = f_one / f_cicj;

    const GMFloat recip_ci = f_cj * recip_2;
    const GMFloat recip_cj = f_ci * recip_2;

    const GMFloat fi0_2 = recip_ci * si0;
    const GMFloat fi1_2 = recip_ci * si1;
    const GMFloat fj0_2 = recip_cj * sj0;
    const GMFloat fj1_2 = recip_cj * sj1;

    // Do simple algebra on threshold inequality to get all constnts on RHS

    const GMFloat ccc_duo_multiplier = env_ccc_duo_multiplier<CBPE>(*env);
    const GMFloat ccc_param = env->ccc_param();

    const GMFloat threshold_multiplier =
             ((f_one*CBPE*CBPE) / ccc_duo_multiplier) * f_cij;

    // Let a few extra values pass the thresold check just to make sure
    const GMFloat roundoff_fuzz = 1.e-5;

    const GMFloat f = threshold * threshold_multiplier * (f_one-roundoff_fuzz);

    const GMFloat v00 = rij00 * ((f_one*CBPE) - ccc_param * fi0_2) *
                                ((f_one*CBPE) - ccc_param * fj0_2);
    const GMFloat v01 = rij01 * ((f_one*CBPE) - ccc_param * fi0_2) *
                                ((f_one*CBPE) - ccc_param * fj1_2);
    const GMFloat v10 = rij10 * ((f_one*CBPE) - ccc_param * fi1_2) *
                                ((f_one*CBPE) - ccc_param * fj0_2);
    const GMFloat v11 = rij11 * ((f_one*CBPE) - ccc_param * fi1_2) *
                                ((f_one*CBPE) - ccc_param * fj1_2);

    // Verify the algebra.
    COMET_ASSERT(fabs(v00 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index, 0, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v01 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index, 0, 1, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v10 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index, 1, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v11 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index, 1, 1, env) )
             < roundoff_fuzz * threshold_multiplier);

    return v00 > f || v01 > f || v10 > f || v11 > f;

  } // if (env->sparse())

  // Non-sparse case (less well-optimized).

  const GMFloat v00 = GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index,
                                                               0, 0, env);
  const GMFloat v01 = GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index,
                                                               0, 1, env);
  const GMFloat v10 = GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index,
                                                               1, 0, env);
  const GMFloat v11 = GMMetrics_ccc_duo_get_from_index_2<CBPE>(metrics, index,
                                                               1, 1, env);

  return v00 > threshold || v01 > threshold ||
         v10 > threshold || v11 > threshold;
}

//-----------------------------------------------------------------------------
/// \brief Check if any 2-way CCC table value may exceed threshold.

static bool GMMetrics_ccc_get_from_index_2_threshold(
  GMMetrics* metrics,
  const size_t index,
  GMFloat threshold,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  enum { COUNTED_BITS_PER_ELT = 2 };

  return GMMetrics_ccc_duo_get_from_index_2_threshold<COUNTED_BITS_PER_ELT>(
    metrics, index, threshold, env);
}

//-----------------------------------------------------------------------------
/// \brief Check if any 2-way DUO table value may exceed threshold.

static bool GMMetrics_duo_get_from_index_2_threshold(
  GMMetrics* metrics,
  const size_t index,
  GMFloat threshold,
  CEnv* env) {

  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  enum { COUNTED_BITS_PER_ELT = 1 };

  return GMMetrics_ccc_duo_get_from_index_2_threshold<COUNTED_BITS_PER_ELT>(
    metrics, index, threshold, env);
}

//=============================================================================
/*---Accessors: value from (local) coord: set: 2-way---*/

template<typename T>
static void GMMetrics_set_2(GMMetrics* metrics, void* p, int i, int j,
  T value, CEnv* env) {
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(!env->all2all());

  const size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((T*)p)[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_2(GMMetrics* metrics, int i, int j,
  GMFloat value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  GMMetrics_set_2<GMFloat>(metrics, metrics->data, i, j, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_2(GMMetrics* metrics, int i, int j,
  GMFloat2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_2<GMFloat2>(metrics, metrics->data_S, i, j, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_2(GMMetrics* metrics, int i, int j,
  GMFloat2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_2<GMFloat2>(metrics, metrics->data_C, i, j, value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_2(GMMetrics* metrics, int i, int j,
  GMTally2x2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_2<GMTally2x2>(metrics, metrics->data, i, j, value, env);
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

template<typename T>
static void GMMetrics_set_all2all_2(GMMetrics* metrics, void* p, int i, int j,
  int j_block, T value, CEnv* env) {
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(i < j || j_block != env->proc_num_vector());
  // WARNING: these conditions on j_block are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_2(metrics, i, j,
                                                            j_block, env);
  ((T*)p)[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_2(GMMetrics* metrics, int i, int j,
  int j_block, GMFloat value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  GMMetrics_set_all2all_2<GMFloat>(metrics, metrics->data, i, j, j_block,
    value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_all2all_2(GMMetrics* metrics, int i, int j,
  int j_block, GMFloat2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_all2all_2<GMFloat2>(metrics, metrics->data_S, i, j, j_block,
    value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_all2all_2(GMMetrics* metrics, int i, int j,
  int j_block, GMFloat2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_all2all_2<GMFloat2>(metrics, metrics->data_C, i, j, j_block,
    value, env);
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_all2all_2(GMMetrics* metrics, int i, int j,
  int j_block, GMTally2x2 value, CEnv* env) {
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY2X2);

  GMMetrics_set_all2all_2<GMTally2x2>(metrics, metrics->data, i, j, j_block,
    value, env);
}

//=============================================================================
/*---Accessors: value from (local) coord: get: 2-way---*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
  int i, int j, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  const size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
  int i, int j, int j_block, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(i < j || j_block != env->proc_num_vector());
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);
  /// WARNING: these conditions on j_block are not exhaustive.

  const size_t index = GMMetrics_index_from_coord_all2all_2(metrics, i, j,
    j_block, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_2way_accessors_hh_

//-----------------------------------------------------------------------------
