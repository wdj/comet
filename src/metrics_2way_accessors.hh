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

#include "metrics_2way_indexing.hh"
#include "formulas.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Accessors: value from (contig) index: basic

#if 0
template<typename T>
static T Metrics_get(GMMetrics& metrics, size_t index, CEnv& env) {
  COMET_ASSERT(index+1 >= 1 && index < metrics.num_elts_local);

  return ((T*)(metrics.data))[index];
}

//-----------------------------------------------------------------------------

// TODO: consolidate these accessors.

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
#endif

//=============================================================================
// Accessors: value from (contig) index: derived

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 2-way CCC value.

static void GMMetrics_ccc_check_size_nofp_2(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (env->metric_type() != MetricType::CCC ||
      env->num_way() != NUM_WAY::_2 || ! env->are_ccc_params_default())
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

  //const GMTally2x2 ttable = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const auto ttable = Metrics_elt_const<GMTally2x2>(*metrics, index, *env);
  const GMTally1 rij = GMTally2x2_get(ttable, i0, i1);

  //const GMFloat2 si1_sj1 =
  //    GMMetrics_float2_S_get_from_index(metrics, index, env);
  const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArray::S>(*metrics,
    index, *env);
  GMTally1 si1, sj1;
  GMFloat2_decode(si1, sj1, si1_sj1);

  GMTally1 ci = 0, cj = 0, cij = 0;

  if (env->sparse()) {
    //const GMFloat2 ci_cj =
    //  GMMetrics_float2_C_get_from_index(metrics, index, env);
    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArray::C>(*metrics,
      index, *env);
    GMFloat2_decode(ci, cj, ci_cj);

    cij = GMTally2x2_get(ttable, 0, 0) + GMTally2x2_get(ttable, 0, 1) +
          GMTally2x2_get(ttable, 1, 0) + GMTally2x2_get(ttable, 1, 1);

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

template<int COUNTED_BITS_PER_ELT, typename FloatResult_t = GMFloat>
static FloatResult_t GMMetrics_ccc_duo_get_from_index_2(
  GMMetrics* metrics, size_t index, int i0, int i1, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->is_metric_type_bitwise());
  COMET_ASSERT(env->counted_bits_per_elt() == COUNTED_BITS_PER_ELT);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);

  if (env->threshold_tc()) {
    typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
    //const auto ttable = Metrics_get<TTable_t>(*metrics, index, *env);
    const auto ttable = Metrics_elt_const<TTable_t>(*metrics, index, *env);
    TTable_t::TypeIn result = TTable_t::get(ttable, i0, i1);
    return (FloatResult_t)result;
  }

  typedef double Float_t; // Perform all calcs in double until return.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const Float_t f_one = 1;
  const Float_t recip_m = metrics->recip_m;

  //const GMTally2x2 ttable = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const auto ttable = Metrics_elt_const<GMTally2x2>(*metrics, index, *env);
  const GMTally1 rij = GMTally2x2_get(ttable, i0, i1);

  //const GMFloat2 si1_sj1 =
  //    GMMetrics_float2_S_get_from_index(metrics, index, env);
  const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArray::S>(*metrics,
    index, *env);
  GMTally1 si1, sj1;
  GMFloat2_decode(si1, sj1, si1_sj1);

  Float_t result_floatcalc = 0;

  if (env->sparse()) {

    //const GMFloat2 ci_cj =
    //  GMMetrics_float2_C_get_from_index(metrics, index, env);
    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArray::C>(*metrics,
      index, *env);
    GMTally1 ci, cj;
    GMFloat2_decode(ci, cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(ttable, 0, 0) + GMTally2x2_get(ttable, 0, 1) +
                   GMTally2x2_get(ttable, 1, 0) + GMTally2x2_get(ttable, 1, 1);
    if (0 == ci || 0 == cj || 0 == cij)
      return (FloatResult_t)0;

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si = i0 == 0 ? (CBPE * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (CBPE * cj - sj1) : sj1;

    // The following complex code is to reduce the number of divides.
    const Float_t f_ci = (Float_t) ci;
    const Float_t f_cj = (Float_t) cj;

    const Float_t f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const Float_t f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const Float_t f_cij = (Float_t) cij;
    const Float_t recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

    const Float_t recip_ci = f_cj * f_cij * recip_cicjcij;
    const Float_t recip_cj = f_ci * f_cij * recip_cicjcij;

    const Float_t recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

    result_floatcalc = ccc_duo_value<COUNTED_BITS_PER_ELT, Float_t>(
      rij, si, sj, recip_ci, recip_cj, recip_sumcij, 
      env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env), env->ccc_param());

  } else { // !env->sparse

    COMET_ASSERT(metrics->num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si =
      i0 == 0 ? (CBPE * metrics->num_field_active - si1) : si1;
    const GMTally1 sj =
      i1 == 0 ? (CBPE * metrics->num_field_active - sj1) : sj1;

    const Float_t recip_sumcij = (f_one / (CBPE*CBPE)) * recip_m;

    result_floatcalc = ccc_duo_value<COUNTED_BITS_PER_ELT, Float_t>(
      rij, si, sj, recip_m, recip_m, recip_sumcij, 
      env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(*env), env->ccc_param());

  } // if (env->sparse())

#ifdef COMET_USE_INT128
  if (env->metric_type() == MetricType::CCC && env->are_ccc_params_default()) {
    const Float_t result_intcalc = GMMetrics_ccc_get_from_index_nofp_2(metrics,
                                         index, i0, i1, env);

    // TODO: CHECK floating point type here
    const double eps = 1. / ( ((size_t)1) << (mantissa_digits<FloatResult_t>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      fprintf(stderr, "Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      COMET_INSIST(diff < eps);
    }
  }
#endif
  return (FloatResult_t)result_floatcalc;
}

//-----------------------------------------------------------------------------
/// \brief Templatized Check if any table value may exceed threshold.

template<int COUNTED_BITS_PER_ELT>
static bool GMMetrics_ccc_duo_get_from_index_2_threshold(
  GMMetrics* metrics, const size_t index, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);

  if (env->threshold_tc()) {
    typedef Tally2x2<MetricFormat::SINGLE> TTable_t;
    //const auto ttable = Metrics_get<TTable_t>(*metrics, index, *env);
    const auto ttable = Metrics_elt_const<TTable_t>(*metrics, index, *env);
    for (int i0 = 0; i0 < 2; ++i0) {
      for (int i1 = 0; i1 < 2; ++i1) {
        if (TTable_t::get(ttable, i0, i1) != (TTable_t::TypeIn)0)
          return true;
      }
    }
    return false;
  }

  typedef double Float_t; // Perform all calcs in double.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const Float_t threshold_eff = env->threshold_eff();

  if (env->sparse()) {

    const Float_t f_one = 1;

    //const GMTally2x2 ttable = GMMetrics_tally2x2_get_from_index(metrics, index,
    //                                                         env);
    const auto ttable = Metrics_elt_const<GMTally2x2>(*metrics, index, *env);
    const GMTally1 rij00 = GMTally2x2_get(ttable, 0, 0);
    const GMTally1 rij01 = GMTally2x2_get(ttable, 0, 1);
    const GMTally1 rij10 = GMTally2x2_get(ttable, 1, 0);
    const GMTally1 rij11 = GMTally2x2_get(ttable, 1, 1);

    //const GMFloat2 si1_sj1 =
    //  GMMetrics_float2_S_get_from_index(metrics, index, env);
    const auto si1_sj1 = Metrics_elt_const<GMFloat2, MetricsArray::S>(*metrics,
      index, *env);
    GMTally1 si1 = 0, sj1 = 0;
    GMFloat2_decode(si1, sj1, si1_sj1);

    //const GMFloat2 ci_cj =
    //  GMMetrics_float2_C_get_from_index(metrics, index, env);
    const auto ci_cj = Metrics_elt_const<GMFloat2, MetricsArray::C>(*metrics,
      index, *env);
    GMTally1 ci = 0, cj = 0;
    GMFloat2_decode(ci, cj, ci_cj);

    GMTally1 cij = rij00 + rij01 + rij10 + rij11;
    if (ci == 0 || cj == 0 || cij == 0) {
      return env->pass_threshold(0);
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

    const Float_t ccc_duo_multiplier = (Float_t)env_ccc_duo_multiplier<CBPE>(*env);
    const Float_t ccc_param = (Float_t)env->ccc_param();

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
      GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index, 0, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v01 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index, 0, 1, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v10 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index, 1, 0, env) )
             < roundoff_fuzz * threshold_multiplier);
    COMET_ASSERT(fabs(v11 - threshold_multiplier *
      GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index, 1, 1, env) )
             < roundoff_fuzz * threshold_multiplier);

    return v00 > f || v01 > f || v10 > f || v11 > f;

  } // if (env->sparse())

  // Non-sparse case (less well-optimized).

  const Float_t v00 = GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index,
                                                               0, 0, env);
  const Float_t v01 = GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index,
                                                               0, 1, env);
  const Float_t v10 = GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index,
                                                               1, 0, env);
  const Float_t v11 = GMMetrics_ccc_duo_get_from_index_2<CBPE, Float_t>(metrics, index,
                                                               1, 1, env);

  return env->pass_threshold(v00) || env->pass_threshold(v01) ||
         env->pass_threshold(v10) || env->pass_threshold(v11);
}

//=============================================================================
// Accessors: value from (local) coord: set: 2-way.









#if 0
template<typename T>
static void GMMetrics_set_2(GMMetrics* metrics, void* p, int i, int j,
  T value, CEnv* env) {
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(!env->all2all());

  //const size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  const size_t index = Metrics_index_2(*metrics, i, j, env->proc_num_vector(), *env);
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
#endif

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

#if 0
template<typename T>
static void GMMetrics_set_all2all_2(GMMetrics* metrics, void* p, int i, int j,
  int j_block, T value, CEnv* env) {
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(i < j || j_block != env->proc_num_vector());
  // WARNING: these conditions on j_block are not exhaustive.

  const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
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
#endif

//-----------------------------------------------------------------------------

template<typename T, int MA = MetricsArray::_>
static T& Metrics_elt_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_2);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || env.proc_num_vector() == j_block);
  // WARNING: these conditions are not exhaustive.

  const size_t index = Metrics_index_2(metrics, i, j, j_block, env);

  return Metrics_elt<T, MA>(metrics, index, env);
}

//=============================================================================
// Accessors: value from (local) coord: get: 2-way.

#if 0
static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
  int i, int j, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_2);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  const size_t index = Metrics_index_2(*metrics, i, j, env->proc_num_vector(), *env);
  return Metrics_elt_const<GMFloat>(*metrics, index, *env);
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

  const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
  return Metrics_elt_const<GMFloat>(*metrics, index, *env);
}

#endif

//-----------------------------------------------------------------------------

template<typename T>
static T Metrics_elt_const_2(GMMetrics& metrics, int i, int j, int j_block,
  CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_2);
  COMET_ASSERT(env.proc_num_repl() == 0 || env.all2all());
  COMET_ASSERT(env.all2all() || env.proc_num_vector() == j_block);
  // WARNING: these conditions are not exhaustive.

  const size_t index = Metrics_index_2(metrics, i, j, j_block, env);

  return Metrics_elt_const<T>(metrics, index, env);
}

//=============================================================================
// Accessors: value from index: get: 2-way.

static GMFloat GMMetrics_get_2(GMMetrics& metrics,
  size_t index, int i0, int i1, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_2);
  COMET_ASSERT(index >= 0 && index < metrics.num_elts_local);
  COMET_ASSERT(i0 >= 0 && i0 < env.i012_max());
  COMET_ASSERT(i1 >= 0 && i1 < env.i012_max());

  const GMFloat result =
    env.metric_type() == MetricType::CZEK ?
      Metrics_elt_const<GMFloat>(metrics, index, env) :
    env.metric_type() == MetricType::CCC ?
      (GMFloat)GMMetrics_ccc_duo_get_from_index_2<CBPE::CCC>(
        &metrics, index, i0, i1, &env) :
  //env.metric_type() == MetricType::DUO ?
      (GMFloat)GMMetrics_ccc_duo_get_from_index_2<CBPE::DUO>(
        &metrics, index, i0, i1, &env);

  return result;
}

//=============================================================================
// Accessors: value from (global) coord: get: 2-way.

static GMFloat GMMetrics_get_2(GMMetrics& metrics,
  size_t ig, size_t jg, int i0, int i1, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_2);
  COMET_ASSERT(i0 >= 0 && i0 < env.i012_max());
  COMET_ASSERT(i1 >= 0 && i1 < env.i012_max());

  // WARNING: these conditions are not exhaustive.

  const size_t i = GMDecompMgr_get_vector_local_from_vector_active(
    metrics.dm, ig, &env);
  const size_t j = GMDecompMgr_get_vector_local_from_vector_active(
    metrics.dm, jg, &env);
  COMET_ASSERT(i >= 0 && i < metrics.dm->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics.dm->num_vector_local);

  const int i_proc = GMDecompMgr_get_proc_vector_from_vector_active(
    metrics.dm, ig, &env);
  const int j_proc = GMDecompMgr_get_proc_vector_from_vector_active(
    metrics.dm, jg, &env);
  no_unused_variable_warning(i_proc);
  COMET_ASSERT(env.proc_num_vector() == i_proc);
  COMET_ASSERT(j_proc >= 0 && j_proc < env.num_proc_vector());

  const int j_proc_eff = env.all2all() ? j_proc : env.proc_num_vector();

  const size_t index = Metrics_index_2(metrics, i, j, j_proc_eff, env);

  const GMFloat result = GMMetrics_get_2(metrics, index, i0, i1, env);

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_2way_accessors_hh_

//-----------------------------------------------------------------------------
