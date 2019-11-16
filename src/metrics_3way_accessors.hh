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

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Accessors: value from (contig) index: basic---*/

static GMFloat3 GMMetrics_float3_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
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
                                                  GMEnv* env) {
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
                                                    GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);

  return ((GMTally4x2*)(metrics->data))[index];
}

//=============================================================================
/*---Accessors: value from (contig) index: derived---*/

static GMFloat GMMetrics_ccc_value_3(GMMetrics* metrics,
                                    const GMTally1 rijk,
                                    const GMTally1 si,
                                    const GMTally1 sj,
                                    const GMTally1 sk,
                                    const GMFloat recip_ci,
                                    const GMFloat recip_cj,
                                    const GMFloat recip_ck,
                                    const GMFloat recip_sumcijk,
                                    GMEnv* env) {
  COMET_ASSERT(metrics && env);

  const GMFloat one = 1;

  const GMFloat fi = (one / 2) * recip_ci * si;
  const GMFloat fj = (one / 2) * recip_cj * sj;
  const GMFloat fk = (one / 2) * recip_ck * sk;

  const GMFloat fijk = recip_sumcijk * rijk;

  /*---Do the following to make floating point arithmetic order-independent---*/

  GMFloat fmin = 0;
  GMFloat fmid = 0;
  GMFloat fmax = 0;

  if (fi > fj) {
    if (fi > fk) {
      fmax = fi;
      if (fj > fk) {
        fmid = fj;
        fmin = fk;
      } else {
        fmid = fk;
        fmin = fj;
      }
    } else {
      fmid = fi;
      fmax = fk;
      fmin = fj;
    }
  } else {
    if (fj > fk) {
      fmax = fj;
      if (fi > fk) {
        fmid = fi;
        fmin = fk;
      } else {
        fmid = fk;
        fmin = fi;
      }
    } else {
      fmid = fj;
      fmax = fk;
      fmin = fi;
    }
  }

  COMET_ASSERT(fmin <= fmid);
  COMET_ASSERT(fmid <= fmax);

  const GMFloat ccc_multiplier = env->ccc_multiplier();
  const GMFloat ccc_param = env->ccc_param();

  /* clang-format off */
  const GMFloat result = ccc_multiplier * fijk * (one - ccc_param * fmin) *
                                                 (one - ccc_param * fmid) *
                                                 (one - ccc_param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------

static void GMMetrics_ccc_check_size_nofp_3(GMMetrics* metrics, GMEnv* env) {
  COMET_INSIST(metrics && env);

#ifdef USE_INT128
  if (env->metric_type() != MetricType::CCC || 
      env->num_way() != NUM_WAY::_3 || ! env->are_ccc_params_default()) {
    return;
  }

  const size_t m = metrics->num_field_active;
  const int lm = utils::log2(m);

  // Bound on log2(numerator)
  const int lnum = 3+lm + 2+lm + 2+lm + 2+lm;

  COMET_INSIST_INTERFACE(env, lnum < 128 && "Number of fields too large.");
#endif
}

//-----------------------------------------------------------------------------

#ifdef USE_INT128

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_value_nofp_3(GMMetrics* metrics,
                                          const GMTally1 rijk,
                                          const GMTally1 si,
                                          const GMTally1 sj,
                                          const GMTally1 sk,
                                          const GMTally1 ci,
                                          const GMTally1 cj,
                                          const GMTally1 ck,
                                          const GMTally1 cijk,
                                          GMEnv* env) {
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

static GMFloat GMMetrics_ccc_get_from_index_nofp_3(GMMetrics* metrics,
                                                   size_t index,
                                                   int i0,
                                                   int i1,
                                                   int i2,
                                                   GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
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

  GMTally1 ci, cj, ck, cijk;

  if (env->sparse()) {
    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    cijk =
           GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
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

static GMFloat GMMetrics_ccc_get_from_index_3(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              int i2,
                                              GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(i0 >= 0 && i0 < 2);
  COMET_ASSERT(i1 >= 0 && i1 < 2);
  COMET_ASSERT(i2 >= 0 && i2 < 2);

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

    GMTally1 cijk =
           GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
           GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
           GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
           GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1);

    if (0 == ci || 0 == cj || 0 == ck || 0 == cijk) {
      return (GMFloat)0;
    }

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;
    const GMTally1 sk = i2 == 0 ? (2 * ck - sk1) : sk1;

    // TODO: it may be possible to decrease the number of divides
    // here - see GMMetrics_ccc_get_from_index_2.
    const GMFloat recip_ci = f_one / ci;
    const GMFloat recip_cj = f_one / cj;
    const GMFloat recip_ck = f_one / ck;

    const GMFloat recip_sumcijk =
      f_one / (GMTally4x2_get(t42, 0, 0, 0) + GMTally4x2_get(t42, 0, 0, 1) +
               GMTally4x2_get(t42, 0, 1, 0) + GMTally4x2_get(t42, 0, 1, 1) +
               GMTally4x2_get(t42, 1, 0, 0) + GMTally4x2_get(t42, 1, 0, 1) +
               GMTally4x2_get(t42, 1, 1, 0) + GMTally4x2_get(t42, 1, 1, 1));

    result_floatcalc = GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk,
                             recip_ci, recip_cj, recip_ck, recip_sumcijk, env);
  } else { /*---if sparse---*/

    COMET_ASSERT(metrics->num_field_active > 0);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (2 * metrics->num_field_active - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (2 * metrics->num_field_active - sj1) : sj1;
    const GMTally1 sk = i2 == 0 ? (2 * metrics->num_field_active - sk1) : sk1;

    const GMFloat recip_sumcijk = (f_one / 8) * recip_m;

    result_floatcalc = GMMetrics_ccc_value_3(metrics, rijk, si, sj, sk, recip_m,
                               recip_m, recip_m, recip_sumcijk, env);

  } /*---if sparse---*/

#ifdef USE_INT128
  if (env->are_ccc_params_default()) {
    const GMFloat result_intcalc = GMMetrics_ccc_get_from_index_nofp_3(metrics,
                                         index, i0, i1, i2, env);

    const double eps = 1. / ( ((size_t)1) << (mantissa_digits<GMFloat>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      printf("Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      COMET_INSIST(diff < eps);
    }

    //COMET_INSIST(diff < eps);

    //return result_intcalc;
    return result_floatcalc;
  } else {
    return result_floatcalc;
  }
#else
  return result_floatcalc;
#endif
}

//-----------------------------------------------------------------------------

static bool GMMetrics_ccc_get_from_index_3_threshold(GMMetrics* metrics,
                                                     const size_t index,
                                                     GMFloat threshold,
                                                     GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(index+1 >= 1 && index < metrics->num_elts_local);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);

  if (env->sparse()) {

    const GMFloat f_one = 1;

    const GMTally4x2 t42 = GMMetrics_tally4x2_get_from_index(metrics, index, env);
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
    GMTally1 si1, sj1, sk1;
    GMFloat3_decode(&si1, &sj1, &sk1, si1_sj1_sk1);

    const GMFloat3 ci_cj_ck =
      GMMetrics_float3_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj, ck;
    GMFloat3_decode(&ci, &cj, &ck, ci_cj_ck);

    GMTally1 cijk = rijk000 + rijk001 + rijk010 + rijk011 +
                    rijk100 + rijk101 + rijk110 + rijk111;
    if (ci == 0 || cj == 0 || ck == 0 || cijk == 0) {
      return 0 > threshold;
    }

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/

    const GMTally1 si0 = 2 * ci - si1;
    const GMTally1 sj0 = 2 * cj - sj1;
    const GMTally1 sk0 = 2 * ck - sk1;

    // TODO: optimize this further

    const GMFloat recip_ci = f_one / ci;
    const GMFloat recip_cj = f_one / cj;
    const GMFloat recip_ck = f_one / ck;

    const GMFloat recip_sumcijk = f_one / cijk;

    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk000, si0, sj0, sk0,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 0, 0, 0, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk001, si0, sj0, sk1,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 0, 0, 1, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk010, si0, sj1, sk0,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 0, 1, 0, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk011, si0, sj1, sk1,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 0, 1, 1, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk100, si1, sj0, sk0,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 1, 0, 0, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk101, si1, sj0, sk1,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 1, 0, 1, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk110, si1, sj1, sk0,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 1, 1, 0, env));
    COMET_ASSERT(GMMetrics_ccc_value_3(metrics, rijk111, si1, sj1, sk1,
                       recip_ci, recip_cj, recip_ck, recip_sumcijk, env) ==
             GMMetrics_ccc_get_from_index_3(metrics, index, 1, 1, 1, env));

    return GMMetrics_ccc_value_3(metrics, rijk000, si0, sj0, sk0,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk001, si0, sj0, sk1,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk010, si0, sj1, sk0,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk011, si0, sj1, sk1,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk100, si1, sj0, sk0,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk101, si1, sj0, sk1,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk110, si1, sj1, sk0,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold ||
           GMMetrics_ccc_value_3(metrics, rijk111, si1, sj1, sk1,
               recip_ci, recip_cj, recip_ck, recip_sumcijk, env) > threshold;

  } /*---if sparse---*/

  GMFloat v000, v001, v010, v011, v100, v101, v110, v111;

  v000 = GMMetrics_ccc_get_from_index_3(metrics, index, 0, 0, 0, env);
  v001 = GMMetrics_ccc_get_from_index_3(metrics, index, 0, 0, 1, env);
  v010 = GMMetrics_ccc_get_from_index_3(metrics, index, 0, 1, 0, env);
  v011 = GMMetrics_ccc_get_from_index_3(metrics, index, 0, 1, 1, env);
  v100 = GMMetrics_ccc_get_from_index_3(metrics, index, 1, 0, 0, env);
  v101 = GMMetrics_ccc_get_from_index_3(metrics, index, 1, 0, 1, env);
  v110 = GMMetrics_ccc_get_from_index_3(metrics, index, 1, 1, 0, env);
  v111 = GMMetrics_ccc_get_from_index_3(metrics, index, 1, 1, 1, env);
  return v000 > threshold || v001 > threshold ||
         v010 > threshold || v011 > threshold ||
         v100 > threshold || v101 > threshold ||
         v110 > threshold || v111 > threshold;
}

//=============================================================================
/*---Accessors: value from (local) coord: set: 3-way---*/

static void GMMetrics_float_set_3(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  int k,
                                  GMFloat value,
                                  GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMFloat3 value,
                                     GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_S);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMFloat3 value,
                                     GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_C);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMTally4x2 value,
                                     GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_3(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int k,
                                          int j_block,
                                          int k_block,
                                          GMFloat value,
                                          GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_all2all_3(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int k,
                                             int j_block,
                                             int k_block,
                                             GMTally4x2 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_3_permuted(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat value,
    GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(
    metrics, I, J, K, j_block, k_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMFloat3 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_all2all_3_permuted(GMMetrics* metrics,
                                             int I,
                                             int J,
                                             int K,
                                             int j_block,
                                             int k_block,
                                             GMTally4x2 value,
                                             GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K,
                                                               j_block, k_block,
                                                               env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_S_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat3 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_S);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat3*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float3_C_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMFloat3 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(metrics->data_C);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
    metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMFloat3*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally4x2_set_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMTally4x2 value,
    GMIndexCache* index_cache,
    GMEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
      metrics, I, J, K, j_block, k_block, index_cache, env);
  ((GMTally4x2*)(metrics->data))[index] = value;
}

//=============================================================================
/*---Accessors: value from (local) coord: get: 3-way---*/

static GMFloat GMMetrics_float_get_3(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     int k,
                                     GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_3(GMMetrics* metrics,
                                           int i,
                                           int j,
                                           int k,
                                           GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(! env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(i < j && j < k);
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);

  size_t index = GMMetrics_index_from_coord_3(metrics, i, j, k, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3(GMMetrics* metrics,
                                                   int i,
                                                   int j,
                                                   int k,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(i >= 0 && i < metrics->num_vector_local);
  COMET_ASSERT(j >= 0 && j < metrics->num_vector_local);
  COMET_ASSERT(k >= 0 && k < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3(metrics, i, j, k, j_block,
                                                      k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3_permuted(GMMetrics* metrics,
                                                   int I,
                                                   int J,
                                                   int K,
                                                   int j_block,
                                                   int k_block,
                                                   GMEnv* env) {
  COMET_ASSERT(metrics && env);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted(metrics, I, J, K, j_block,
                                                      k_block, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMTally4x2 GMMetrics_tally4x2_get_all2all_3_permuted_cache(
    GMMetrics* metrics,
    int I,
    int J,
    int K,
    int j_block,
    int k_block,
    GMIndexCache* index_cache,
    GMEnv* env) {
  COMET_ASSERT(metrics && env && index_cache);
  COMET_ASSERT(env->all2all());
  COMET_ASSERT(I >= 0 && I < metrics->num_vector_local);
  COMET_ASSERT(J >= 0 && J < metrics->num_vector_local);
  COMET_ASSERT(K >= 0 && K < metrics->num_vector_local);
  COMET_ASSERT(j_block >= 0 && j_block < env->num_block_vector());
  COMET_ASSERT(k_block >= 0 && k_block < env->num_block_vector());
  COMET_ASSERT(env->num_way() == NUM_WAY::_3);
  COMET_ASSERT(env->data_type_metrics() == GM_DATA_TYPE_TALLY4X2);
  /*---WARNING: these conditions are not exhaustive---*/

  size_t index = GMMetrics_index_from_coord_all2all_3_permuted_cache(
      metrics, I, J, K, j_block, k_block, index_cache, env);
  return GMMetrics_tally4x2_get_from_index(metrics, index, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_metrics_3way_accessors_hh_

//-----------------------------------------------------------------------------
