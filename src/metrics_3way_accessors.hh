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
// Accessors: value from (contig) index: basic

// TODO: consolidate these accessors.

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
// Accessors: value from (contig) index: derived

//-----------------------------------------------------------------------------
/// \brief Templatized access to the CCC or DUO formula.

template<int COUNTED_BITS_PER_ELT>
static double Metrics_ccc_duo_value(
  GMMetrics& metrics,
  const GMTally1 rijk,
  const GMTally1 si,
  const GMTally1 sj,
  const GMTally1 sk,
  const double recip_ci,
  const double recip_cj,
  const double recip_ck,
  const double recip_sumcijk,
  CEnv& env) {

  return ccc_duo_value<COUNTED_BITS_PER_ELT, double>(
    rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk,
    env_ccc_duo_multiplier<COUNTED_BITS_PER_ELT>(env),
    env.ccc_param());
}

//-----------------------------------------------------------------------------
/// \brief Check to ensure 128 bits is enough to store 3-way CCC value.

static void GMMetrics_ccc_check_size_nofp_3(GMMetrics* metrics, CEnv* env) {
  COMET_INSIST(metrics && env);

  if (env->metric_type() != MetricType::CCC || 
      env->num_way() != NUM_WAY::_3 || ! env->are_ccc_params_default())
    return;

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

  const GMTally4x2 ttable = GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(ttable, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  GMTally1 ci = 0, cj = 0, ck = 0, cijk = 0;

  if (env->sparse()) {
    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    cijk = GMTally4x2_get(ttable, 0, 0, 0) + GMTally4x2_get(ttable, 0, 0, 1) +
           GMTally4x2_get(ttable, 0, 1, 0) + GMTally4x2_get(ttable, 0, 1, 1) +
           GMTally4x2_get(ttable, 1, 0, 0) + GMTally4x2_get(ttable, 1, 0, 1) +
           GMTally4x2_get(ttable, 1, 1, 0) + GMTally4x2_get(ttable, 1, 1, 1);

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

template<int COUNTED_BITS_PER_ELT, typename FloatResult_t = GMFloat>
static FloatResult_t GMMetrics_ccc_duo_get_from_index_3(
  GMMetrics* metrics, size_t index, int i0, int i1, int i2, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->is_metric_type_bitwise());
  COMET_ASSERT(env->counted_bits_per_elt() == COUNTED_BITS_PER_ELT);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

  if (env->threshold_tc()) {
    typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
    const auto table = Metrics_get<TTable_t>(*metrics, index, *env);
    TTable_t::TypeIn result = TTable_t::get(table, i0, i1, i2);
    return (FloatResult_t)result;
  }

  typedef double Float_t; // Perform all calcs in double until return.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  const Float_t f_one = 1;
  const Float_t recip_m = metrics->recip_m;

  const GMTally4x2 ttable = GMMetrics_tally4x2_get_from_index(metrics, index, env);
  const GMTally1 rijk = GMTally4x2_get(ttable, i0, i1, i2);

  const GMFloat3 si1_sj1_sk1 =
      GMMetrics_float3_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1, sk1;
  GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

  Float_t result_floatcalc = 0;

  if (env->sparse()) {

    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj, ck;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    // TODO: it may be possible to decrease the number of divides
    // here - see GMMetrics_ccc_get_from_index_2.
    const Float_t recip_ci = f_one / ci;
    const Float_t recip_cj = f_one / cj;
    const Float_t recip_ck = f_one / ck;

    GMTally1 cijk =
           GMTally4x2_get(ttable, 0, 0, 0) + GMTally4x2_get(ttable, 0, 0, 1) +
           GMTally4x2_get(ttable, 0, 1, 0) + GMTally4x2_get(ttable, 0, 1, 1) +
           GMTally4x2_get(ttable, 1, 0, 0) + GMTally4x2_get(ttable, 1, 0, 1) +
           GMTally4x2_get(ttable, 1, 1, 0) + GMTally4x2_get(ttable, 1, 1, 1);
    if (0 == ci || 0 == cj || 0 == ck || 0 == cijk)
      return (FloatResult_t)0;

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si = i0 == 0 ? (CBPE * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (CBPE * cj - sj1) : sj1;
    const GMTally1 sk = i2 == 0 ? (CBPE * ck - sk1) : sk1;

    const Float_t recip_sumcijk = f_one / cijk;
//      f_one / (GMTally4x2_get(ttable, 0, 0, 0) + GMTally4x2_get(ttable, 0, 0, 1) +
//               GMTally4x2_get(ttable, 0, 1, 0) + GMTally4x2_get(ttable, 0, 1, 1) +
//               GMTally4x2_get(ttable, 1, 0, 0) + GMTally4x2_get(ttable, 1, 0, 1) +
//               GMTally4x2_get(ttable, 1, 1, 0) + GMTally4x2_get(ttable, 1, 1, 1));

    result_floatcalc = Metrics_ccc_duo_value<CBPE>(*metrics,
      rijk, si, sj, sk, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env);

  } else { // !env->sparse

    COMET_ASSERT(metrics->num_field_active > 0);

    // Get number of 1 bits OR get number of 0 bits from number of 1 bits.
    const GMTally1 si =
      i0 == 0 ? (CBPE * metrics->num_field_active - si1) : si1;
    const GMTally1 sj =
      i1 == 0 ? (CBPE * metrics->num_field_active - sj1) : sj1;
    const GMTally1 sk =
      i2 == 0 ? (CBPE * metrics->num_field_active - sk1) : sk1;

    const Float_t recip_sumcijk = (f_one / (CBPE*CBPE*CBPE)) * recip_m;

    result_floatcalc = Metrics_ccc_duo_value<CBPE>(*metrics,
      rijk, si, sj, sk, recip_m, recip_m, recip_m, recip_sumcijk, *env);

  } // if (env->sparse())

#ifdef COMET_USE_INT128
  if (env->metric_type() == MetricType::CCC && env->are_ccc_params_default()) {
    const Float_t result_intcalc = GMMetrics_ccc_get_from_index_nofp_3(metrics,
                                         index, i0, i1, i2, env);

    // TODO: CHECK floating point type here
    const double eps = 1. / ( ((size_t)1) << (mantissa_digits<FloatResult_t>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      printf("Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
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
static bool GMMetrics_ccc_duo_get_from_index_3_threshold(
  GMMetrics* metrics, const size_t index, CEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index < metrics->num_elts_local); // && index >= 0
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  if (env->threshold_tc()) {
    typedef Tally4x2<MetricFormat::SINGLE> TTable_t;
    const auto ttable = Metrics_get<TTable_t>(*metrics, index, *env);
    for (int i0 = 0; i0 < 2; ++i0) {
      for (int i1 = 0; i1 < 2; ++i1) {
        for (int i2 = 0; i2 < 2; ++i2) {
          if (TTable_t::get(ttable, i0, i1, i2) != (TTable_t::TypeIn)0)
            return true;
        }
      }
    }
    return false;
  }

  typedef double Float_t; // Perform all calcs in double.

  enum {CBPE = COUNTED_BITS_PER_ELT};

  if (env->sparse()) {

    const Float_t f_one = 1;

    const GMTally4x2 ttable = GMMetrics_tally4x2_get_from_index(metrics, index,
                                                             env);
    const GMTally1 rijk000 = GMTally4x2_get(ttable, 0, 0, 0);
    const GMTally1 rijk001 = GMTally4x2_get(ttable, 0, 0, 1);
    const GMTally1 rijk010 = GMTally4x2_get(ttable, 0, 1, 0);
    const GMTally1 rijk011 = GMTally4x2_get(ttable, 0, 1, 1);
    const GMTally1 rijk100 = GMTally4x2_get(ttable, 1, 0, 0);
    const GMTally1 rijk101 = GMTally4x2_get(ttable, 1, 0, 1);
    const GMTally1 rijk110 = GMTally4x2_get(ttable, 1, 1, 0);
    const GMTally1 rijk111 = GMTally4x2_get(ttable, 1, 1, 1);

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
      return env->pass_threshold(0);
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

    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk000,
      si0, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 0, 0, 0, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk001,
      si0, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 0, 0, 1, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk010,
      si0, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 0, 1, 0, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk011,
      si0, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 0, 1, 1, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk100,
      si1, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 1, 0, 0, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk101,
      si1, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 1, 0, 1, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk110,
      si1, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 1, 1, 0, env)));
    COMET_ASSERT((Metrics_ccc_duo_value<CBPE>(*metrics, rijk111,
      si1, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ==
      (GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index, 1, 1, 1, env)));

    return env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk000,
           si0, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk001,
           si0, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk010,
           si0, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk011,
           si0, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk100,
           si1, sj0, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk101,
           si1, sj0, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk110,
           si1, sj1, sk0, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env)) ||
           env->pass_threshold(Metrics_ccc_duo_value<CBPE>(*metrics, rijk111,
           si1, sj1, sk1, recip_ci, recip_cj, recip_ck, recip_sumcijk, *env));

  } // if (env->sparse())

  // Non-sparse case (less well-optimized).

  const Float_t v000 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                0, 0, 0, env);
  const Float_t v001 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                0, 0, 1, env);
  const Float_t v010 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                0, 1, 0, env);
  const Float_t v011 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                0, 1, 1, env);
  const Float_t v100 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                1, 0, 0, env);
  const Float_t v101 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                1, 0, 1, env);
  const Float_t v110 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                1, 1, 0, env);
  const Float_t v111 = GMMetrics_ccc_duo_get_from_index_3<CBPE, Float_t>(metrics, index,
                                                                1, 1, 1, env);

  return env->pass_threshold(v000) || env->pass_threshold(v001) ||
         env->pass_threshold(v010) || env->pass_threshold(v011) ||
         env->pass_threshold(v100) || env->pass_threshold(v101) ||
         env->pass_threshold(v110) || env->pass_threshold(v111);
}

//=============================================================================
// Accessors: value from (local) coord: set: 3-way.

template<typename T>
static void GMMetrics_set_3(GMMetrics* metrics, void* p, int i, int j, int k,
  T value, CEnv* env) {    
  COMET_ASSERT(metrics && p && env);
  COMET_ASSERT(i >= 0 && i < j);
  COMET_ASSERT(j >= 0 && j < k);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(!env->all2all());

  //const size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  const size_t index = Metrics_index_3(*metrics, i, j, k, *env);
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

template<int MF>
static void Metrics_set(
  GMMetrics& metrics, int i, int j, int k, Tally4x2<MF> value, CEnv& env) {

  const size_t index = Metrics_index_3(metrics, i, j, k, env);
  ((Tally4x2<MF>*)(metrics.data))[index] = value;
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

  //const size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
  //  metrics, I, J, K, j_block, k_block, index_cache, env);
  const size_t index = Metrics_index_3(
    *metrics, I, J, K, j_block, k_block, *index_cache, *env);
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

//-----------------------------------------------------------------------------

template<int MF>
static void Metrics_set(
  GMMetrics& metrics, int I, int J, int K, int j_block, int k_block,
  Tally4x2<MF> value, GMIndexCache& index_cache, CEnv& env) {

  const size_t index = Metrics_index_3(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((Tally4x2<MF>*)(metrics.data))[index] = value;
}

//=============================================================================
// Accessors: value from (local) coord: get: 3-way.

template<int MF>
static Tally4x2<MF> Metrics_get(GMMetrics& metrics,
  int i, int j, int k, CEnv& env) {

  const size_t index = Metrics_index_3(metrics, i, j, k, env);
  return Metrics_get<Tally4x2<MF>>(metrics, index, env);
}

//-----------------------------------------------------------------------------

template<int MF>
static Tally4x2<MF> Metrics_get(
  GMMetrics& metrics, int I, int J, int K, int j_block, int k_block,
  GMIndexCache& index_cache, CEnv& env) {

  const size_t index = Metrics_index_3(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  return Metrics_get<Tally4x2<MF>>(metrics, index, env);
}

//=============================================================================
// Accessors: value from index: get: 3-way.

static GMFloat GMMetrics_get_3(GMMetrics& metrics,
  size_t index, int i0, int i1, int i2, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_3);
  COMET_ASSERT(index >= 0 && index < metrics.num_elts_local);
  COMET_ASSERT(i0 >= 0 && i0 < env.i012_max());
  COMET_ASSERT(i1 >= 0 && i1 < env.i012_max());
  COMET_ASSERT(i2 >= 0 && i2 < env.i012_max());

  const GMFloat result =
    env.metric_type() == MetricType::CZEK ?
      Metrics_get<GMFloat>(metrics, index, env) :
    env.metric_type() == MetricType::CCC ?
      (GMFloat)GMMetrics_ccc_duo_get_from_index_3<CBPE::CCC>(
        &metrics, index, i0, i1, i2, &env) :
  //env.metric_type() == MetricType::DUO ?
      (GMFloat)GMMetrics_ccc_duo_get_from_index_3<CBPE::DUO>(
        &metrics, index, i0, i1, i2, &env);

  return result;
}

//=============================================================================
// Accessors: value from (global) coord: get: 3-way.

static GMFloat GMMetrics_get_3(GMMetrics& metrics,
  size_t ig, size_t jg, size_t kg, int i0, int i1, int i2, CEnv& env) {
  COMET_ASSERT(env.num_way() == NUM_WAY::_3);
  COMET_ASSERT(i0 >= 0 && i0 < env.i012_max());
  COMET_ASSERT(i1 >= 0 && i1 < env.i012_max());
  COMET_ASSERT(i2 >= 0 && i2 < env.i012_max());
  // WARNING: these conditions are not exhaustive.

  const size_t i = GMDecompMgr_get_vector_local_from_vector_active(
    metrics.dm, ig, &env);
  const size_t j = GMDecompMgr_get_vector_local_from_vector_active(
    metrics.dm, jg, &env);
  const size_t k = GMDecompMgr_get_vector_local_from_vector_active(
    metrics.dm, kg, &env);
  COMET_ASSERT(i >= 0 && i < metrics.dm->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics.dm->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics.dm->num_vector_local);

  const int i_proc = GMDecompMgr_get_proc_vector_from_vector_active(
    metrics.dm, ig, &env);
  const int j_proc = GMDecompMgr_get_proc_vector_from_vector_active(
    metrics.dm, jg, &env);
  const int k_proc = GMDecompMgr_get_proc_vector_from_vector_active(
    metrics.dm, kg, &env);
  no_unused_variable_warning(i_proc);
  COMET_ASSERT(env.proc_num_vector() == i_proc);
  COMET_ASSERT(j_proc >= 0 && j_proc < env.num_proc_vector());
  COMET_ASSERT(k_proc >= 0 && k_proc < env.num_proc_vector());

  const size_t index = env.all2all() ?
    GMMetrics_index_from_coord_all2all_3(&metrics,
      i, j, k, j_proc, k_proc, &env) :
    Metrics_index_3(metrics, i, j, k, env);
    //GMMetrics_index_from_coord_3(&metrics, i, j, k, &env);

  const GMFloat result = GMMetrics_get_3(metrics, index, i0, i1, i2, env);

  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_3way_accessors_hh_

//-----------------------------------------------------------------------------
