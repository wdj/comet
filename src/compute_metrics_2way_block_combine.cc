//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_combine.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Combine numerators and denominators, 2-way, for a single block.
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

#include "cstdint"
#define __STDC_FORMAT_MACROS 1
#include "inttypes.h"

#include "env.hh"
#include "linalg.hh"
#include "mirrored_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way_block_combine.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Combine nums and denoms on CPU to get final result, 2-way Czek.

void gm_compute_2way_proc_combine_czek_(
  GMMetrics* metrics,
  MirroredBuf* metrics_buf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
   bool do_compute_triang_only,
   CEnv* env) {

  COMET_INSIST(metrics && metrics_buf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  // NOTE: here and elsewhere, vector_sums_left and vector_sums_right
  // may sometimes be aliases.  Since they are read-only, there should
  // be no danger of the compiler generating incorrect code.
  // (cf. https://en.wikipedia.org/wiki/Pointer_aliasing)
  // However by accounting for this one might be able to in principle
  // remove a load instruction to improve performance.

  const int nvl = metrics->num_vector_local;
  const VectorSums* const vs_l = vector_sums_left;
  const VectorSums* const vs_r = vector_sums_right;

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  // ----------------------------------
  if (!env->is_compute_method_gpu() && env->all2all()) {
    // ----------------------------------

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = vs_r->sum(j);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer = Metrics_elt_const_2<GMFloat>(*metrics,
          i, j, j_block, *env);
        const GMFloat vs_i = vs_l->sum(i);
        const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        Metrics_elt_2<GMFloat>(*metrics, i, j, j_block, *env) = value;
      } // for i
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j
        /*---TODO: here and elsewhere check for unlikely case denom is/nearly
         * zero---*/

    // ----------------------------------
  } else if (!env->is_compute_method_gpu()) {
    // ----------------------------------

    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = vs_r->sum(j);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer = Metrics_elt_const_2<GMFloat>(*metrics,
          i, j, env->proc_num_vector(), *env);
        const GMFloat vs_i = vs_l->sum(i);
        const GMFloat denom = vs_i < vs_j ?  vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        Metrics_elt_2<GMFloat>(*metrics, i, j, env->proc_num_vector(), *env) = value;
      } // for i
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j

    // ----------------------------------
  } else if (env->all2all()) {
    // ----------------------------------

    if (do_compute_triang_only) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMFloat vs_j = vs_r->sum(j);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMFloat numer =
            metrics_buf->elt_const<GMFloat>(i, j);
          const GMFloat vs_i = vs_l->sum(i);
          const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
          const GMFloat multiplier = (GMFloat)2;
          const GMFloat value = (multiplier * numer) / denom;
          Metrics_elt_2<GMFloat>(*metrics, i, j, j_block, *env) = value;
        } // for i
      }   // for j
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_metric_items_local_computed_inc(i_max);
      }   // for j
    } else {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMFloat vs_j = vs_r->sum(j);
          const GMFloat numer =
            metrics_buf->elt_const<GMFloat>(i, j);
          const GMFloat vs_i = vs_l->sum(i);
          const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
          const GMFloat multiplier = (GMFloat)2;
          const GMFloat value = (multiplier * numer) / denom;
          Metrics_elt_2<GMFloat>(*metrics, i, j, j_block, *env) = value;
        } // for i
      }   // for j
      metrics->num_metric_items_local_computed_inc(nvl * (size_t)nvl);
    }

    // ----------------------------------
  } else {
    // ----------------------------------

    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = vs_r->sum(j);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
          metrics_buf->elt_const<GMFloat>(i, j);
        const GMFloat vs_i = vs_l->sum(i);
        const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        Metrics_elt_2<GMFloat>(*metrics, i, j, env->proc_num_vector(), *env) = value;
      } // for i
    }   // for j
    for (int j = 0; j < nvl; ++j) {
      const int i_max = j;
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j

    // ----------------------------------
  } // if
  // ----------------------------------
}

//=============================================================================
// Combine nums and denoms on CPU to get final result, 2-way CCC.

void gm_compute_2way_proc_combine_ccc_(
  GMMetrics* metrics,
  MirroredBuf* metrics_buf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {
  COMET_INSIST(metrics && metrics_buf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  const int nvl = metrics->num_vector_local;
  const VectorSums* const vs_l = vector_sums_left;
  const VectorSums* const vs_r = vector_sums_right;

  // Copy from metrics_buffer for GPU case; perform checks.

  if (env->is_using_linalg()) {
    // --------------
    if (env->all2all()) {
      // --------------

      if (do_compute_triang_only) {
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          const int i_max = j;
          for (int i = 0; i < i_max; ++i) {
            const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
            Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block, *env) = value;
#if 1
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse()) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(4 * metrics->num_field_active);
              if (error1) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " m %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       coords, env->proc_num());
                COMET_INSIST((! error1) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 si1 = vs_l->sum(i);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(2 * si1);
              if (error2) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(2 * si1),
                       coords, env->proc_num());
                COMET_INSIST((! error2) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 sj1 = vs_r->sum(j);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(2 * sj1);
              if (error3) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(2 * sj1),
                       coords, env->proc_num());
                COMET_INSIST((! error3) && "Violation of algorithm computational invariant.");
              }
            }
#endif
#ifdef COMET_ASSERTIONS_ON
            if (! env->sparse()) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              COMET_ASSERT((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                           (uint64_t)r11 ==
                       (uint64_t)(4 * metrics->num_field_active));
              // 2-sum checks.
              const GMTally1 si1 = vs_l->sum(i);
              const GMTally1 sj1 = vs_r->sum(j);
              COMET_ASSERT((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
              COMET_ASSERT((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
            }
#endif
          } // for i
        }   // for j
      } else {
        // don't use collapse because of overflow for large sizes
        //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          for (int i = 0; i < nvl; ++i) {
            const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
            Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block, *env) = value;
#if 1
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse()) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(4 * metrics->num_field_active);
              if (error1) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " m %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       coords, env->proc_num());
                COMET_INSIST((! error1) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 si1 = vs_l->sum(i);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(2 * si1);
              if (error2) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(2 * si1),
                       coords, env->proc_num());
                COMET_INSIST((! error2) && "Violation of algorithm computational invariant.");
              }
              // 2-sum check.
              const GMTally1 sj1 = vs_r->sum(j);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(2 * sj1);
              if (error3) {
                const size_t index = Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error: r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(2 * sj1),
                       coords, env->proc_num());
                COMET_INSIST((! error3) && "Violation of algorithm computational invariant.");
              }
            }
#endif
#ifdef COMET_ASSERTIONS_ON
            if (! env->sparse()) {
              // 4-sum check.
              const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
              const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
              const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
              const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
              COMET_ASSERT((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                           (uint64_t)r11 ==
                       (uint64_t)(4 * metrics->num_field_active));
              // 2-sum checks.
              const GMTally1 si1 = vs_l->sum(i);
              const GMTally1 sj1 = vs_r->sum(j);
              COMET_ASSERT((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
              COMET_ASSERT((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
            }
#endif
          } // for i
        }   // for j
     }

      // --------------
    } else { // (! env->all2all())
      // --------------
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
          Metrics_elt_2<GMTally2x2>(*metrics, i, j, env->proc_num_vector(), *env) = value;
#ifdef COMET_ASSERTIONS_ON
          if (! env->sparse()) {
            // 4-sum check.
            const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
            const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
            const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
            const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
            COMET_ASSERT((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                         (uint64_t)r11 ==
                     (uint64_t)(4 * metrics->num_field_active));
            // 2-sum checks.
            const GMTally1 si1 = vs_l->sum(i);
            const GMTally1 sj1 = vs_r->sum(j);
            COMET_ASSERT((uint64_t)r10 + (uint64_t)r11 == (uint64_t)(2 * si1));
            COMET_ASSERT((uint64_t)r01 + (uint64_t)r11 == (uint64_t)(2 * sj1));
          }
#endif
        } // for i
      }   // for j
      // --------------
    } // if
    // --------------
  }

  // Compute multipliers.

  // --------------
  if (env->all2all()) {
    // --------------

    if (do_compute_triang_only) {
      if (!env->is_threshold_tc()) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = vs_r->sum(j);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, j_block, *env) = si1_sj1;
          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse
        }   // for i
      }   // for j
      } // if (!env->is_threshold_tc())
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_metric_items_local_computed_inc(i_max);
      }   // for j
    } else {
      if (!env->is_threshold_tc()) {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMTally1 sj1 = vs_r->sum(j);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, j_block, *env) = si1_sj1;
          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse
        }   // for i
      }   // for j
      } // if (!env->is_threshold_tc())
      metrics->num_metric_items_local_computed_inc(nvl * (size_t)nvl);
   }

    // --------------
  } else { // (! env->all2all())
    // --------------
    if (!env->is_threshold_tc()) {
    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMTally1 sj1 = vs_r->sum(j);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMTally1 si1 = vs_l->sum(i);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
        Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, env->proc_num_vector(), *env) = si1_sj1;
        if (env->sparse()) {
          const GMTally1 ci = vs_l->count(i);
          const GMTally1 cj = vs_r->count(j);
          const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
          Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, env->proc_num_vector(), *env) = ci_cj;
        } // if sparse
      } // for i
    }   // for j
    } // if (!env->is_threshold_tc())
    for (int j = 0; j < nvl; ++j) {
      const int i_max = do_compute_triang_only ? j : nvl;
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j
    // --------------
  } // if
  // --------------
}

//=============================================================================
// Combine nums and denoms on CPU to get final result, 2-way DUO.

void gm_compute_2way_proc_combine_duo_(
  GMMetrics* metrics,
  MirroredBuf* metrics_buf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && metrics_buf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  const int nvl = metrics->num_vector_local;
  const VectorSums* const vs_l = vector_sums_left;
  const VectorSums* const vs_r = vector_sums_right;

  // Copy from metrics_buffer for GPU case.

  if (env->is_using_linalg()) {

    // NOTE: this technically uses the wrong type for copying if env->is_threshold_tc(),
    // but still works ok.

    // --------------
    if (env->all2all()) {
      // --------------

      if (do_compute_triang_only) {
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          const int i_max = j;
          for (int i = 0; i < i_max; ++i) {
            const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
            Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block, *env) = value;
          } // for i
        }   // for j
      } else {
        // don't use collapse because of overflow for large sizes
        //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
        #pragma omp parallel for schedule(dynamic,1000)
        for (int j = 0; j < nvl; ++j) {
          for (int i = 0; i < nvl; ++i) {
            const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
            Metrics_elt_2<GMTally2x2>(*metrics, i, j, j_block, *env) = value;
          } // for i
        }   // for j
     }

      // --------------
    } else { // (! env->all2all())
      // --------------
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally2x2 value =
              metrics_buf->elt_const<GMTally2x2>(i, j);
          Metrics_elt_2<GMTally2x2>(*metrics, i, j, env->proc_num_vector(), *env) = value;
        } // for i
      }   // for j
      // --------------
    } // if
    // --------------
  }

  // Compute multipliers.

  // --------------
  if (env->all2all()) {
    // --------------

    if (do_compute_triang_only) {
      if (!env->is_threshold_tc()) {
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = vs_r->sum(j);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, j_block, *env) = si1_sj1;
          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse
        }   // for i
      }   // for j
      } // if (!env->is_threshold_tc())
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_metric_items_local_computed_inc(i_max);
      }   // for j
    } else {
      if (!env->is_threshold_tc()) {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMTally1 sj1 = vs_r->sum(j);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, j_block, *env) = si1_sj1;
          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse
        }   // for i
      }   // for j
      } // if (!env->is_threshold_tc())
      metrics->num_metric_items_local_computed_inc(nvl * (size_t)nvl);
   }

    // --------------
  } else { // (! env->all2all())
    // --------------
    if (!env->is_threshold_tc()) {
    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMTally1 sj1 = vs_r->sum(j);
      const int i_max = do_compute_triang_only ? j : nvl;
      for (int i = 0; i < i_max; ++i) {
        const GMTally1 si1 = vs_l->sum(i);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
        Metrics_elt_2<GMFloat2, MetricsArray::S>(*metrics, i, j, env->proc_num_vector(), *env) = si1_sj1;
        if (env->sparse()) {
          const GMTally1 ci = vs_l->count(i);
          const GMTally1 cj = vs_r->count(j);
          const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
          Metrics_elt_2<GMFloat2, MetricsArray::C>(*metrics, i, j, env->proc_num_vector(), *env) = ci_cj;
        } // if sparse
      } // for i
    }   // for j
    } // if (!env->is_threshold_tc())
    for (int j = 0; j < nvl; ++j) {
      const int i_max = do_compute_triang_only ? j : nvl;
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j
    // --------------
  } // if
  // --------------
}

//=============================================================================
// Combine nums and denoms on CPU to get final result, 2-way generic.

void gm_compute_2way_proc_combine(
  GMMetrics* metrics,
  MirroredBuf* metrics_buf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && metrics_buf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      gm_compute_2way_proc_combine_czek_(metrics, metrics_buf,
                                         vector_sums_left, vector_sums_right,
                                         j_block, do_compute_triang_only, env);
    } break;
    case MetricType::CCC: {
      gm_compute_2way_proc_combine_ccc_(metrics, metrics_buf,
                                        vector_sums_left, vector_sums_right,
                                        j_block, do_compute_triang_only, env);
    } break;
    case MetricType::DUO: {
      gm_compute_2way_proc_combine_duo_(metrics, metrics_buf,
                                        vector_sums_left, vector_sums_right,
                                        j_block, do_compute_triang_only, env);
    } break;
    default:
      COMET_INSIST_INTERFACE(env, false && "Selected metric_type unimplemented.");
  } // case
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
