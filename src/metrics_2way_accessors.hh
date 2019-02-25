//-----------------------------------------------------------------------------
/*!
 * \file   metrics_2way_accessors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class, header, 2-way, accessor functions.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_metrics_2way_accessors_hh_
#define _gm_metrics_2way_accessors_hh_

#include "metrics_2way_indexing.hh"

//=============================================================================
/*---Accessors: value from (contig) index: basic---*/

static GMFloat GMMetrics_float_get_from_index(GMMetrics* metrics,
                                              size_t index,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(metrics->data))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_S_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);

  return ((GMFloat2*)(metrics->data_S))[index];
}

//-----------------------------------------------------------------------------

static GMFloat2 GMMetrics_float2_C_get_from_index(GMMetrics* metrics,
                                                  size_t index,
                                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);

  return ((GMFloat2*)(metrics->data_C))[index];
}

//-----------------------------------------------------------------------------

static GMTally2x2 GMMetrics_tally2x2_get_from_index(GMMetrics* metrics,
                                                    size_t index,
                                                    GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  return ((GMTally2x2*)(metrics->data))[index];
}

//=============================================================================
/*---Accessors: value from (contig) index: derived---*/


static GMFloat GMMetrics_czek_get_from_index(GMMetrics* metrics,
                                             size_t index,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);

  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_value_2(GMMetrics* metrics,
                                     const GMTally1 rij,
                                     const GMTally1 si,
                                     const GMTally1 sj,
                                     const GMFloat recip_ci,
                                     const GMFloat recip_cj,
                                     const GMFloat recip_sumcij,
                                     GMEnv* env) {
  GMAssert(metrics && env);

  const GMFloat f_one = 1;

  const GMFloat fi = (f_one / 2) * recip_ci * si;
  const GMFloat fj = (f_one / 2) * recip_cj * sj;

  const GMFloat fij = recip_sumcij * rij;

  /*---Do the following to make floating point arithmetic order-independent---*/
  const GMFloat fmin = fi < fj ? fi : fj;
  const GMFloat fmax = fi < fj ? fj : fi;

  const GMFloat ccc_multiplier = GMEnv_ccc_multiplier(env);
  const GMFloat ccc_param = GMEnv_ccc_param(env);

  /* clang-format off */
  const GMFloat result = ccc_multiplier * fij * (f_one - ccc_param * fmin) *
                                                (f_one - ccc_param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_duo_value_2(GMMetrics* metrics,
                                     const GMTally1 rij,
                                     const GMTally1 si,
                                     const GMTally1 sj,
                                     const GMFloat recip_ci,
                                     const GMFloat recip_cj,
                                     const GMFloat recip_sumcij,
                                     GMEnv* env) {
  GMAssert(metrics && env);

  const GMFloat f_one = 1;

  const GMFloat fi = recip_ci * si;
  const GMFloat fj = recip_cj * sj;
//printf("fi  %f\n", (double)fi);
//printf("fj  %f\n", (double)fj);

  const GMFloat fij = recip_sumcij * rij;
//printf("fij  %f\n", (double)fij);

  /*---Do the following to make floating point arithmetic order-independent---*/
  const GMFloat fmin = fi < fj ? fi : fj;
  const GMFloat fmax = fi < fj ? fj : fi;

  const GMFloat ccc_multiplier = GMEnv_ccc_multiplier(env);
  const GMFloat ccc_param = GMEnv_ccc_param(env);
//printf("ccc_multiplier  %f\n", (double)ccc_multiplier);
//printf("ccc_param  %f\n", (double)ccc_param);

  /* clang-format off */
  const GMFloat result = ccc_multiplier * fij * (f_one - ccc_param * fmin) *
                                                (f_one - ccc_param * fmax);
  /* clang-format on */

  return result;
}

//-----------------------------------------------------------------------------

static void GMMetrics_ccc_check_size_nofp_2(GMMetrics* metrics, GMEnv* env) {
  GMInsist(metrics && env);

#ifdef HAVE_INT128
  if (GMEnv_num_way(env) != GM_NUM_WAY_2 || ! env->are_ccc_params_default) {
    return;
  }

  const size_t m = metrics->num_field_active;
  const int lm = gm_log2(m);

  // Bound on log2(numerator)
  const int lnum = 2+lm + 2+lm + 2+lm; 

  const int shift = gm_mant_dig<GMFloat>() - 3; // Note num/denom <= 4.5 < 1<<3

  // To make this less restrictive, see 3-way counterpart code.
  GMInsistInterface(env, lnum + shift < 128 && "Number of fields too large.");
#endif
}

//-----------------------------------------------------------------------------

#ifdef HAVE_INT128

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_value_nofp_2(GMMetrics* metrics,
                                          const GMTally1 rij,
                                          const GMTally1 si,
                                          const GMTally1 sj,
                                          const GMTally1 ci,
                                          const GMTally1 cj,
                                          const GMTally1 cij,
                                          GMEnv* env) {
  GMAssert(metrics && env);

  const GMUInt128 num = rij * (GMUInt128)(3 * ci - 1 * si) *
                              (GMUInt128)(3 * cj - 1 * sj);

  const GMUInt128 denom = 2 * cij * (GMUInt128)ci * (GMUInt128)cj;

  const int shift = gm_mant_dig<GMFloat>() - 3; // Note num/denom <= 4.5 < 1<<3

  // This should be an integer with no more than
  // gm_mant_dig<GMFloat>() binary digits
  const GMUInt128 int_ratio = (num << shift) / denom;

  // Convert to floting point and then adjust exponent.
  const GMFloat result = ( (GMFloat) int_ratio ) /
                         ( (GMFloat) ( ((size_t)1) << shift ) );

  return result;
}


//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_get_from_index_nofp_2(GMMetrics* metrics,
                                                   size_t index,
                                                   int i0,
                                                   int i1,
                                                   GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);
  GMAssert(env->are_ccc_params_default);

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMTally1 ci, cj, cij;

  if (env->sparse) {
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

static GMFloat GMMetrics_ccc_get_from_index_2(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMFloat f_one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMFloat result_floatcalc = 0;

  if (env->sparse) {

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
    if (0 == ci || 0 == cj || 0 == cij) {
      return (GMFloat)0;
    }

    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;

    const GMFloat f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const GMFloat f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const GMFloat f_cij = (GMFloat) cij;
    const GMFloat recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;

    const GMFloat recip_ci = f_cj * f_cij * recip_cicjcij;
    const GMFloat recip_cj = f_ci * f_cij * recip_cicjcij;

    const GMFloat recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

    result_floatcalc = GMMetrics_ccc_value_2(metrics, rij, si, sj,
                                 recip_ci, recip_cj, recip_sumcij, env);

  } else { /*---if sparse---*/

    GMAssert(metrics->num_field_active > 0);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (2 * metrics->num_field_active - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (2 * metrics->num_field_active - sj1) : sj1;

    const GMFloat recip_sumcij = (f_one / 4) * recip_m;

    result_floatcalc = GMMetrics_ccc_value_2(metrics, rij, si, sj,
                                   recip_m, recip_m, recip_sumcij, env);
    //GMFloat v00, v01, v10, v11;
    //GMMetrics_ccc_get_from_index_2_all(metrics, index, v00, v01, v10, v11, env);
    //const GMFloat value = i0 ? (i1 ? v11 : v10) : (i1 ? v01 : v00);
    //GMAssert(result == value);

  } /*---if sparse---*/

#ifdef HAVE_INT128
  if (env->are_ccc_params_default) {
    const GMFloat result_intcalc = GMMetrics_ccc_get_from_index_nofp_2(metrics,
                                         index, i0, i1, env);

    const double eps = 1. / ( ((size_t)1) << (gm_mant_dig<GMFloat>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      printf("Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      GMInsist(diff < eps);
    }

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

static GMFloat GMMetrics_duo_get_from_index_2(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMFloat f_one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMFloat result_floatcalc = 0;

  if (env->sparse) {

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
//FIX - something maybe not right here, divison by 4 - ??
    if (0 == ci || 0 == cj || 0 == cij) {
      return (GMFloat)0;
    }

    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;

    const GMFloat f_cicj_min = f_ci < f_cj ? f_ci : f_cj;
    const GMFloat f_cicj_max = f_ci > f_cj ? f_ci : f_cj;

    const GMFloat f_cij = (GMFloat) cij;
    const GMFloat recip_cicjcij = f_one / (f_cicj_min * f_cicj_max * f_cij);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (ci - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (cj - sj1) : sj1;
//printf("%i %i   %i %i   %i %i\n", i0, i1, ci, cj, si, sj);

    const GMFloat recip_ci = f_cj * f_cij * recip_cicjcij;
    const GMFloat recip_cj = f_ci * f_cij * recip_cicjcij;

    const GMFloat recip_sumcij = f_cicj_min * f_cicj_max * recip_cicjcij;

    result_floatcalc = GMMetrics_duo_value_2(metrics, rij, si, sj,
                                 recip_ci, recip_cj, recip_sumcij, env);

  } else { /*---if sparse---*/

    GMAssert(metrics->num_field_active > 0);

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
    const GMTally1 si = i0 == 0 ? (metrics->num_field_active - si1) : si1;
    const GMTally1 sj = i1 == 0 ? (metrics->num_field_active - sj1) : sj1;

    const GMFloat recip_sumcij = recip_m;

    result_floatcalc = GMMetrics_duo_value_2(metrics, rij, si, sj,
                                   recip_m, recip_m, recip_sumcij, env);
    //GMFloat v00, v01, v10, v11;
    //GMMetrics_ccc_get_from_index_2_all(metrics, index, v00, v01, v10, v11, env);
    //const GMFloat value = i0 ? (i1 ? v11 : v10) : (i1 ? v01 : v00);
    //GMAssert(result == value);

  } /*---if sparse---*/

#if 0
#ifdef HAVE_INT128
  if (env->are_ccc_params_default) {
    const GMFloat result_intcalc = GMMetrics_ccc_get_from_index_nofp_2(metrics,
                                         index, i0, i1, env);

    const double eps = 1. / ( ((size_t)1) << (gm_mant_dig<GMFloat>() - 5) );

    const double diff = fabs(result_intcalc - result_floatcalc);

    if (!(diff < eps)) {
      printf("Error: mismatch result_floatcalc %.16e result_intcalc %.16e\n",
             (double)result_floatcalc, (double)result_intcalc);
      GMInsist(diff < eps);
    }

    //return result_intcalc;
    return result_floatcalc;
  } else {
    return result_floatcalc;
  }
#else
  return result_floatcalc;
#endif
#endif
  return result_floatcalc;
}

//-----------------------------------------------------------------------------

//FIX: make DUO version of this.
static bool GMMetrics_ccc_get_from_index_2_threshold(GMMetrics* metrics,
                                                     const size_t index,
                                                     GMFloat threshold,
                                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);

  if (env->sparse) {

    const GMFloat f_one = 1;

    const GMFloat ccc_multiplier = GMEnv_ccc_multiplier(env);
    const GMFloat ccc_param = GMEnv_ccc_param(env);

    const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
    const GMTally1 rij00 = GMTally2x2_get(t22, 0, 0);
    const GMTally1 rij01 = GMTally2x2_get(t22, 0, 1);
    const GMTally1 rij10 = GMTally2x2_get(t22, 1, 0);
    const GMTally1 rij11 = GMTally2x2_get(t22, 1, 1);

    const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
    GMTally1 si1, sj1;
    GMFloat2_decode(&si1, &sj1, si1_sj1);

    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMTally1 ci, cj;
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = rij00 + rij01 + rij10 + rij11;
    if (ci == 0 || cj == 0 || cij == 0) {
      return 0 > threshold;
    }

    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;
    const GMFloat f_cij = (GMFloat) cij;

    //const GMFloat f_cicj = f_ci < f_cj ? f_ci * f_cj : f_cj * f_ci;
    const GMFloat f_cicj = f_ci * f_cj;

    const GMFloat recip_2 = f_one / f_cicj;

    const GMFloat recip_ci = f_cj * recip_2;
    const GMFloat recip_cj = f_ci * recip_2;

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/

    const GMTally1 si0 = 2 * ci - si1;
    const GMTally1 sj0 = 2 * cj - sj1;

    const GMFloat fi0_2 = recip_ci * si0;
    const GMFloat fi1_2 = recip_ci * si1;
    const GMFloat fj0_2 = recip_cj * sj0;
    const GMFloat fj1_2 = recip_cj * sj1;

    /*---Make floating point arithmetic order-independent---*/

    const GMFloat threshold_multiplier = ((f_one*4) / ccc_multiplier) * f_cij;


    const GMFloat f = threshold * threshold_multiplier * (f_one - 1.e-5);;

    const GMFloat v00 = rij00 * ((f_one*2) - ccc_param * fi0_2) *
                                ((f_one*2) - ccc_param * fj0_2);
    const GMFloat v01 = rij01 * ((f_one*2) - ccc_param * fi0_2) *
                                ((f_one*2) - ccc_param * fj1_2);
    const GMFloat v10 = rij10 * ((f_one*2) - ccc_param * fi1_2) *
                                ((f_one*2) - ccc_param * fj0_2);
    const GMFloat v11 = rij11 * ((f_one*2) - ccc_param * fi1_2) *
                                ((f_one*2) - ccc_param * fj1_2);

    GMAssert(fabs(v00 - threshold_multiplier * GMMetrics_ccc_get_from_index_2(metrics, index, 0, 0, env) )
             < 1.e-5 * threshold_multiplier);
    GMAssert(fabs(v01 - threshold_multiplier * GMMetrics_ccc_get_from_index_2(metrics, index, 0, 1, env) )
             < 1.e-5 * threshold_multiplier);
    GMAssert(fabs(v10 - threshold_multiplier * GMMetrics_ccc_get_from_index_2(metrics, index, 1, 0, env) )
             < 1.e-5 * threshold_multiplier);
    GMAssert(fabs(v11 - threshold_multiplier * GMMetrics_ccc_get_from_index_2(metrics, index, 1, 1, env) )
             < 1.e-5 * threshold_multiplier);

    return v00 > f || v01 > f || v10 > f || v11 > f;

#if 0
    const GMFloat recip_m = metrics->recip_m;
    const GMFloat f_ci = (GMFloat) ci;
    const GMFloat f_cj = (GMFloat) cj;

    const GMFloat f_cicj = f_ci < f_cj ? f_ci * f_cj : f_cj * f_ci;

    const GMFloat f_cij = (GMFloat) cij;
    const GMFloat recip_3 = f_one / (f_cicj * f_cij);

    const GMFloat recip_sumcij = f_cicj * recip_3;

    const GMFloat recip_ci = f_cj * f_cij * recip_3;
    const GMFloat recip_cj = f_ci * f_cij * recip_3;

    /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/

    const GMTally1 si0 = 2 * ci - si1;
    const GMTally1 sj0 = 2 * cj - sj1;

    const GMFloat fi0 = (f_one / 2) * recip_ci * si0;
    const GMFloat fi1 = (f_one / 2) * recip_ci * si1;
    const GMFloat fj0 = (f_one / 2) * recip_cj * sj0;
    const GMFloat fj1 = (f_one / 2) * recip_cj * sj1;

    const GMFloat fij00 = recip_sumcij * rij00;
    const GMFloat fij01 = recip_sumcij * rij01;
    const GMFloat fij10 = recip_sumcij * rij10;
    const GMFloat fij11 = recip_sumcij * rij11;

    /*---Make floating point arithmetic order-independent---*/

    GMFloat v00, v01, v10, v11;

    if (fi0 < fj0) {
      v00 = ccc_multiplier * fij00 * (f_one - ccc_param * fi0) *
                                     (f_one - ccc_param * fj0);
    } else {
      v00 = ccc_multiplier * fij00 * (f_one - ccc_param * fj0) *
                                     (f_one - ccc_param * fi0);
    }

    if (fi0 < fj1) {
      v01 = ccc_multiplier * fij01 * (f_one - ccc_param * fi0) *
                                     (f_one - ccc_param * fj1);
    } else {
      v01 = ccc_multiplier * fij01 * (f_one - ccc_param * fj1) *
                                     (f_one - ccc_param * fi0);
    }

    if (fi1 < fj0) {
      v10 = ccc_multiplier * fij10 * (f_one - ccc_param * fi1) *
                                     (f_one - ccc_param * fj0);
    } else {
      v10 = ccc_multiplier * fij10 * (f_one - ccc_param * fj0) *
                                     (f_one - ccc_param * fi1);
    }

    if (fi1 < fj1) {
      v11 = ccc_multiplier * fij11 * (f_one - ccc_param * fi1) *
                                     (f_one - ccc_param * fj1);
    } else {
      v11 = ccc_multiplier * fij11 * (f_one - ccc_param * fj1) *
                                     (f_one - ccc_param * fi1);
    }

    return v00 > threshold || v01 > threshold ||
           v10 > threshold || v11 > threshold;;
#endif


  } /*---if sparse---*/

  GMFloat v00, v01, v10, v11;

  v00 = GMMetrics_ccc_get_from_index_2(metrics, index, 0, 0, env);
  v01 = GMMetrics_ccc_get_from_index_2(metrics, index, 0, 1, env);
  v10 = GMMetrics_ccc_get_from_index_2(metrics, index, 1, 0, env);
  v11 = GMMetrics_ccc_get_from_index_2(metrics, index, 1, 1, env);
  return v00 > threshold || v01 > threshold ||
         v10 > threshold || v11 > threshold;

}

#if 0
//-----------------------------------------------------------------------------

static GMFloat GMMetrics_ccc_get_from_index_2(GMMetrics* metrics,
                                              size_t index,
                                              int i0,
                                              int i1,
                                              GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(index+1 >= 1 && index < metrics->num_elts_local);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(i0 >= 0 && i0 < 2);
  GMAssert(i1 >= 0 && i1 < 2);

  const GMFloat one = 1;
  const GMFloat recip_m = metrics->recip_m;

  const GMTally2x2 t22 = GMMetrics_tally2x2_get_from_index(metrics, index, env);
  const GMTally1 rij = GMTally2x2_get(t22, i0, i1);

  const GMFloat2 si1_sj1 =
      GMMetrics_float2_S_get_from_index(metrics, index, env);
  GMTally1 si1, sj1;
  GMFloat2_decode(&si1, &sj1, si1_sj1);

  GMTally1 ci, cj;
  if (env->sparse) {
    const GMFloat2 ci_cj =
      GMMetrics_float2_C_get_from_index(metrics, index, env);
    GMFloat2_decode(&ci, &cj, ci_cj);

    GMTally1 cij = GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
                   GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1);
    if (ci == 0 || cj == 0 || cij == 0) {
      return 0;
    }
  } else {
    ci = metrics->num_field_active;
    cj = metrics->num_field_active;
    GMAssert(metrics->num_field_active > 0);
  }

  /*---Get number of 1 bits OR get number of 0 bits from number of 1 bits---*/
  const GMTally1 si = i0 == 0 ? (2 * ci - si1) : si1;
  const GMTally1 sj = i1 == 0 ? (2 * cj - sj1) : sj1;

  const GMFloat recip_ci = env->sparse ? one / ci : recip_m;
  const GMFloat recip_cj = env->sparse ? one / cj : recip_m;

  const GMFloat recip_sumcij = env->sparse ?
    one / (GMTally2x2_get(t22, 0, 0) + GMTally2x2_get(t22, 0, 1) +
           GMTally2x2_get(t22, 1, 0) + GMTally2x2_get(t22, 1, 1)) :
    (one / 4) * recip_m;

  return GMMetrics_ccc_value_2(metrics, rij, si, sj,
                               recip_ci, recip_cj, recip_sumcij, env);
}
#endif

//=============================================================================
//=============================================================================
/*---Accessors: value from (local) coord: set: 2-way---*/

static void GMMetrics_float_set_2(GMMetrics* metrics,
                                  int i,
                                  int j,
                                  GMFloat value,
                                  GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMFloat2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat2*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMFloat2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMFloat2*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMTally2x2 value,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float_set_all2all_2(GMMetrics* metrics,
                                          int i,
                                          int j,
                                          int j_block,
                                          GMFloat value,
                                          GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat*)(metrics->data))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_S_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMFloat2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_S);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat2*)(metrics->data_S))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_float2_C_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMFloat2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  GMAssert(metrics->data_C);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMFloat2*)(metrics->data_C))[index] = value;
}

//-----------------------------------------------------------------------------

static void GMMetrics_tally2x2_set_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMTally2x2 value,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_TALLY2X2);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  ((GMTally2x2*)(metrics->data))[index] = value;
}

//=============================================================================
/*---Accessors: value from (local) coord: get: 2-way---*/

static GMFloat GMMetrics_float_get_2(GMMetrics* metrics,
                                     int i,
                                     int j,
                                     GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(! GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(i < j);
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);

  size_t index = GMMetrics_index_from_coord_2(metrics, i, j, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//-----------------------------------------------------------------------------

static GMFloat GMMetrics_float_get_all2all_2(GMMetrics* metrics,
                                             int i,
                                             int j,
                                             int j_block,
                                             GMEnv* env) {
  GMAssert(metrics && env);
  GMAssert(GMEnv_num_way(env) == GM_NUM_WAY_2);
  GMAssert(GMEnv_all2all(env));
  GMAssert(i >= 0);
  GMAssert(i < metrics->num_vector_local);
  GMAssert(j >= 0);
  GMAssert(j < metrics->num_vector_local);
  GMAssert(j_block >= 0);
  GMAssert(j_block < GMEnv_num_block_vector(env));
  GMAssert(i < j || j_block != GMEnv_proc_num_vector_i(env));
  GMAssert(GMEnv_data_type_metrics(env) == GM_DATA_TYPE_FLOAT);
  /*---WARNING: these conditions on j_block are not exhaustive---*/

  size_t index =
      GMMetrics_index_from_coord_all2all_2(metrics, i, j, j_block, env);
  return GMMetrics_float_get_from_index(metrics, index, env);
}

//=============================================================================

#endif // _gm_metrics_2way_accessors_hh_

//-----------------------------------------------------------------------------
