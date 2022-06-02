//-----------------------------------------------------------------------------
/*!
 * \file   compute_metrics_2way_block_finalize.cc
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
#include "compressed_buf.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "vector_sums.hh"
#include "compute_metrics_2way_block.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Combine nums and denoms on CPU to get final result, 2-way Czek.

template<int METRIC_FORMAT>
static void finalize_czek_(
  GMMetrics* metrics,
  CompressedBuf* matB_cbuf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && matB_cbuf);
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
       For GPU case, directly access the matB_cbuf holding the numerators.
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
  } else if (!env->is_compute_method_gpu()) { // && !env->all2all()
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
        Metrics_elt_2<GMFloat>(*metrics, i, j, env->proc_num_vector(), *env) =
          value;
      } // for i
      metrics->num_metric_items_local_computed_inc(i_max);
    }   // for j

    // ----------------------------------
  } else if (env->all2all()) { // && env->is_compute_method_gpu()
    // ----------------------------------

    if (do_compute_triang_only) {
      COMET_INSIST(!matB_cbuf->do_compress());
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        const GMFloat vs_j = vs_r->sum(j);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMFloat numer =
            matB_cbuf->elt_const<GMFloat>(i, j);
          const GMFloat vs_i = vs_l->sum(i);
          const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
          const GMFloat multiplier = (GMFloat)2;
          const GMFloat value = (multiplier * numer) / denom;

//if (i==0 && j==4)
//printf("1 %.20e %.20e %.20e\n",
//(double)numer,
//(double)denom,
//(double)value);

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
      COMET_INSIST(!matB_cbuf->do_compress());
      #pragma omp parallel for schedule(dynamic,1000)
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMFloat vs_j = vs_r->sum(j);
          const GMFloat numer =
            matB_cbuf->elt_const<GMFloat>(i, j);
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
  } else { // !env->all2all()) && env->is_compute_method_gpu()
    // ----------------------------------

    COMET_INSIST(!matB_cbuf->do_compress());
    #pragma omp parallel for schedule(dynamic,1000)
    for (int j = 0; j < nvl; ++j) {
      const GMFloat vs_j = vs_r->sum(j);
      const int i_max = j;
      for (int i = 0; i < i_max; ++i) {
        const GMFloat numer =
          matB_cbuf->elt_const<GMFloat>(i, j);
        const GMFloat vs_i = vs_l->sum(i);
        const GMFloat denom = vs_i < vs_j ? vs_i + vs_j : vs_j + vs_i;
        const GMFloat multiplier = (GMFloat)2;
        const GMFloat value = (multiplier * numer) / denom;
        Metrics_elt_2<GMFloat>(*metrics, i, j, env->proc_num_vector(), *env) =
          value;
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

template<int METRIC_FORMAT>
static void finalize_ccc_duo_(
  GMMetrics* metrics,
  CompressedBuf* matB_cbuf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && matB_cbuf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  const int nvl = metrics->num_vector_local;
  const VectorSums* const vs_l = vector_sums_left;
  const VectorSums* const vs_r = vector_sums_right;

  enum {MF = METRIC_FORMAT};
  typedef MetricFormatTraits<MF> MFT;
  typedef typename MFT::TypeIn MFTypeIn;

  matB_cbuf->elt_read_start();

  // Copy from matB_cbuf for linalg case; perform checks.

  if (env->is_shrink()) { // && env->is_using_linalg() -- this always true here

    COMET_ASSERT(env->is_threshold_tc());

    // NOTE: this may be slight overestimate of amt of mem that will be needed.

    COMET_INSIST(metrics->num_metric_items_local_computed +
      matB_cbuf->num_entries() <=
      metrics->num_metric_items_local_allocated && 
      "Insufficient metrics memory; please decrease metrics_shrink.");

    const int i_block = env->proc_num_vector();

    // Loop over all table entries stored in compressed buffer.

    for (size_t ind_entry = 0; ind_entry < matB_cbuf->num_entries();
         ++ind_entry) {

      // Read current item (i.e., entry).
      const MFTypeIn metric_item = matB_cbuf->elt_const<MFTypeIn>(ind_entry);

      // If this buf did not do_compress, may actually have zeros.
      // Make guarantee that if is_shrink, all zeros
      // (or to be exact, failed-threshold) values are removed from metrics.

// FIXTHRESHOLD: this needs to be essentially a check for zero.
// equivalent (CHECK) to pass_threshold here if single thresh case, but may not otherwise.
// this should be EXACTLY a check for zero
      //if (!env->pass_threshold(metric_item))
      //if (!env->thresholds().is_pass(metric_item))
      if (Thresholds::is_zero(metric_item))
        continue;

      // Location to store it (item number in metrics array).
      const size_t index = metrics->num_metric_items_local_computed;
      COMET_ASSERT(index < metrics->num_metric_items_local_allocated);

      // Get row, col nums of item just read.

      const size_t j = matB_cbuf->ind1_recent();
      const size_t i = matB_cbuf->ind0_recent();

      // It was computed by GEMM; check is it an entry we need.

      const bool is_in_range = do_compute_triang_only ?
        j >= 0 && j < (size_t)nvl && i >= 0 && i < j :
        j >= 0 && j < (size_t)nvl && i >= 0 && i < (size_t)nvl;

      if (!is_in_range)
        continue;

      // TODO: accessor functions
      const size_t iG = i + nvl * i_block;
      const size_t jG = j + nvl * j_block;

      const int iE = matB_cbuf->iE_recent();
      const int jE = matB_cbuf->jE_recent();

      // Store metric item.

      Metrics_elt<MFTypeIn>(*metrics, index, *env) = metric_item;

      // Store the coords information for this metric item.
      // TODO: accessor function
      metrics->coords_[index] =
        CoordsInfo::set(iG, jG, iE, jE, *metrics, *env);

      metrics->num_metric_items_local_computed_inc(1);

    } // for ind_entry

  } else if (env->is_using_linalg()) { // && ! env->is_shrink()

    // --------------
    if (env->all2all() && do_compute_triang_only) {
    // --------------

      // here and below don't use collapse because of overflow for large sizes
#     pragma omp parallel for schedule(dynamic,1000) if (!matB_cbuf->do_compress())
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const auto value = matB_cbuf->elt_const<Tally2x2<MF>>(i, j);
//printf("finalize   %i %i %f %f %f %f\n", i, j,
//(double)value.get(0,0), (double)value.get(0,1), (double)value.get(1,0), (double)value.get(1,1));
          Metrics_elt_2<Tally2x2<MF>>(*metrics, i, j, j_block, *env) = value;
#         ifdef COMET_ASSERTIONS_ON
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse() && ! env->is_threshold_tc()) {
              const int cbpe = env->counted_bits_per_elt();
              //-----4-sum check.
              const auto r00 = Tally2x2<MF>::get(value, 0, 0);
              const auto r01 = Tally2x2<MF>::get(value, 0, 1);
              const auto r10 = Tally2x2<MF>::get(value, 1, 0);
              const auto r11 = Tally2x2<MF>::get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(cbpe * cbpe * metrics->num_field_active);
              if (error1) {
                const size_t index =
                  Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64
                        " m %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error1 &&
                  "Violation of algorithm computational invariant.");
              }
              //-----2-sum check.
              const GMTally1 si1 = vs_l->sum(i);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(cbpe * si1);
              if (error2) {
                const size_t index =
                  Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(cbpe * si1),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error2 &&
                  "Violation of algorithm computational invariant.");
              }
              //-----2-sum check.
              const GMTally1 sj1 = vs_r->sum(j);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(cbpe * sj1);
              if (error3) {
                const size_t index =
                  Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(cbpe * sj1),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error3 &&
                  "Violation of algorithm computational invariant.");
              }
            } // if (! env->sparse() && ! env->is_threshold_tc())
#         endif // COMET_ASSERTIONS_ON
        } // for i
      }   // for j

    // --------------
    } else if (env->all2all()) { // && ! do_compute_triang_only
    // --------------

#     pragma omp parallel for schedule(dynamic,1000) if (!matB_cbuf->do_compress())
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const Tally2x2<MF> value =
            matB_cbuf->elt_const<Tally2x2<MF>>(i, j);
          Metrics_elt_2<Tally2x2<MF>>(*metrics, i, j, j_block, *env) = value;
#         ifdef COMET_ASSERTIONS_ON
            // ISSUE: this check may increase runtime nontrivially
            if (! env->sparse() && ! env->is_threshold_tc()) {
              const int cbpe = env->counted_bits_per_elt();
              //-----4-sum check.
              const auto r00 = Tally2x2<MF>::get(value, 0, 0);
              const auto r01 = Tally2x2<MF>::get(value, 0, 1);
              const auto r10 = Tally2x2<MF>::get(value, 1, 0);
              const auto r11 = Tally2x2<MF>::get(value, 1, 1);
              const bool error1 = (uint64_t)r00 + (uint64_t)r01 +
                                  (uint64_t)r10 + (uint64_t)r11 !=
                       (uint64_t)(cbpe * cbpe * metrics->num_field_active);
              if (error1) {
                const size_t index =
                   Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " m %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)(metrics->num_field_active),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error1 &&
                  "Violation of algorithm computational invariant.");
              }
              //-----2-sum check.
              const GMTally1 si1 = vs_l->sum(i);
              const bool error2 = (uint64_t)r10 + (uint64_t)r11 !=
                                  (uint64_t)(cbpe * si1);
              if (error2) {
                const size_t index =
                  Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " si1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)si1,
                       (uint64_t)r10 + (uint64_t)r11, (uint64_t)(cbpe * si1),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error2 &&
                  "Violation of algorithm computational invariant.");
              }
              //-----2-sum check.
              const GMTally1 sj1 = vs_r->sum(j);
              const bool error3 = (uint64_t)r01 + (uint64_t)r11 !=
                                  (uint64_t)(cbpe * sj1);
              if (error3) {
                const size_t index =
                  Metrics_index_2(*metrics, i, j, j_block, *env);
                const MetricItemCoords_t coords = metrics->coords_value(index);
                fprintf(stderr, "Error:"
                        " r00 %" PRIu64 " r01 %" PRIu64
                        " r10 %" PRIu64 " r11 %" PRIu64 " sj1 %" PRIu64
                        " actual %" PRIu64 " expected %" PRIu64
                        " coords %zu %zu"
                        " rank %i\n",
                       (uint64_t)r00, (uint64_t)r01, (uint64_t)r10,
                       (uint64_t)r11, (uint64_t)sj1,
                       (uint64_t)r01 + (uint64_t)r11, (uint64_t)(cbpe * sj1),
                       CoordsInfo::getiG(coords, *metrics, *env),
                       CoordsInfo::getjG(coords, *metrics, *env),
                       env->proc_num());
                COMET_INSIST(! error3 &&
                  "Violation of algorithm computational invariant.");
              }
            } // if (! env->sparse() && ! env->is_threshold_tc())
#         endif // COMET_ASSERTIONS_ON
        } // for i
      }   // for j

    // --------------
    } else { // if (! env->all2all())
    // --------------

      const int j_block = env->proc_num_vector();
#     pragma omp parallel for schedule(dynamic,1000) if (!matB_cbuf->do_compress())
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const Tally2x2<MF> value =
              matB_cbuf->elt_const<Tally2x2<MF>>(i, j);
          Metrics_elt_2<Tally2x2<MF>>(*metrics, i, j, j_block, *env) = value;
#         ifdef COMET_ASSERTIONS_ON
            if (! env->sparse() && ! env->is_threshold_tc()) {
              const int cbpe = env->counted_bits_per_elt();
              //-----4-sum check.
              const auto r00 = Tally2x2<MF>::get(value, 0, 0);
              const auto r01 = Tally2x2<MF>::get(value, 0, 1);
              const auto r10 = Tally2x2<MF>::get(value, 1, 0);
              const auto r11 = Tally2x2<MF>::get(value, 1, 1);
              COMET_ASSERT((uint64_t)r00 + (uint64_t)r01 + (uint64_t)r10 +
                           (uint64_t)r11 ==
                       (uint64_t)(cbpe * cbpe * metrics->num_field_active));
              //-----2-sum checks.
              const GMTally1 si1 = vs_l->sum(i);
              const GMTally1 sj1 = vs_r->sum(j);
              COMET_ASSERT((uint64_t)r10 + (uint64_t)r11 ==
                           (uint64_t)(cbpe * si1));
              COMET_ASSERT((uint64_t)r01 + (uint64_t)r11 ==
                           (uint64_t)(cbpe * sj1));
            }
#         endif
        } // for i
      }   // for j

    // --------------
    } // if (env->all2all() && do_compute_triang_only)
    // --------------

  } // if (env->is_shrink()) // if (env->is_using_linalg())

  // Compute multipliers; update counts.

  enum {S = MetricsArray::S};
  enum {C = MetricsArray::C};

  // --------------
  if (env->all2all() && do_compute_triang_only) {
  // --------------

    if (!env->is_threshold_tc()) {
      COMET_INSIST(!matB_cbuf->do_compress());
#     pragma omp parallel for schedule(dynamic,1000) if (!metrics->is_computing_histograms())
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = vs_r->sum(j);
        const int i_max = j;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, S>(*metrics, i, j, j_block, *env) = si1_sj1;

          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse

          if (metrics->is_computing_histograms() &&
              CoordsInfo::is_active(
                metrics->coords_value(Metrics_index_2(*metrics, i, j, j_block, *env)),
                *metrics, *env)) {
            const int cbpe = env->counted_bits_per_elt();
            const int nfa = metrics->num_field_active;
            const GMTally1 ci = env->sparse() ? vs_l->count(i) :
                                                cbpe * cbpe * nfa;
            const GMTally1 cj = env->sparse() ? vs_r->count(j) :
                                                cbpe * cbpe * nfa;
            Tally2x2<MF> ttable = Metrics_elt_2<Tally2x2<MF>>(*metrics,
              i, j, j_block, *env);
            metrics->dm->histograms()->add(ttable, si1, sj1, ci, cj, nfa);
          } // if is_computing_histograms
        }   // for i
      }   // for j
    } // if (!env->is_threshold_tc())

    // Update counts.
    if (!env->is_shrink()) {
      for (int j = 0; j < nvl; ++j) {
        const int i_max = j;
        metrics->num_metric_items_local_computed_inc(i_max);
      }   // for j
    }

  // --------------
  } else if (env->all2all()) { // && ! do_compute_triang_only
  // --------------

    if (!env->is_threshold_tc()) {
      COMET_INSIST(!matB_cbuf->do_compress());
#     pragma omp parallel for schedule(dynamic,1000) if (!metrics->is_computing_histograms())
      for (int j = 0; j < nvl; ++j) {
        for (int i = 0; i < nvl; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMTally1 sj1 = vs_r->sum(j);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, S>(*metrics, i, j, j_block, *env) = si1_sj1;

          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse

          if (metrics->is_computing_histograms() &&
              CoordsInfo::is_active(
                metrics->coords_value(Metrics_index_2(*metrics, i, j, j_block, *env)),
                *metrics, *env)) {
//printf("B %i %i %i %i %i %i %i\n",
//env->proc_num_vector(),
//i, j, j_block, (int)metrics->num_vector_active, (int)nvl, (int)metrics->dm->num_vector_active_local) ;
            const int cbpe = env->counted_bits_per_elt();
            const int nfa = metrics->num_field_active;
            const GMTally1 ci = env->sparse() ? vs_l->count(i) :
                                                cbpe * cbpe * nfa;
            const GMTally1 cj = env->sparse() ? vs_r->count(j) :
                                                cbpe * cbpe * nfa;
            Tally2x2<MF> ttable = Metrics_elt_2<Tally2x2<MF>>(*metrics,
              i, j, j_block, *env);
            metrics->dm->histograms()->add(ttable, si1, sj1, ci, cj, nfa);
          } // if is_computing_histograms
        }   // for i
      }   // for j
    } // if (!env->is_threshold_tc())

    // Update counts.
    if (!env->is_shrink()) {
      metrics->num_metric_items_local_computed_inc(nvl * (size_t)nvl);
    }

  // --------------
  } else { // ! env->all2all()
  // --------------

    if (!env->is_threshold_tc()) {
      const int j_block = env->proc_num_vector();
      COMET_INSIST(!matB_cbuf->do_compress());
#     pragma omp parallel for schedule(dynamic,1000) if (!metrics->is_computing_histograms())
      for (int j = 0; j < nvl; ++j) {
        const GMTally1 sj1 = vs_r->sum(j);
        const int i_max = do_compute_triang_only ? j : nvl;
        for (int i = 0; i < i_max; ++i) {
          const GMTally1 si1 = vs_l->sum(i);
          const GMFloat2 si1_sj1 = GMFloat2_encode(si1, sj1);
          Metrics_elt_2<GMFloat2, S>(*metrics, i, j, j_block, *env) = si1_sj1;

          if (env->sparse()) {
            const GMTally1 ci = vs_l->count(i);
            const GMTally1 cj = vs_r->count(j);
            const GMFloat2 ci_cj = GMFloat2_encode(ci, cj);
            Metrics_elt_2<GMFloat2, C>(*metrics, i, j, j_block, *env) = ci_cj;
          } // if sparse

          if (metrics->is_computing_histograms() &&
              (size_t)i < metrics->dm->num_vector_active_local &&
              (size_t)j < metrics->dm->num_vector_active_local) {
            const int cbpe = env->counted_bits_per_elt();
            const int nfa = metrics->num_field_active;
            const GMTally1 ci = env->sparse() ? vs_l->count(i) :
                                                cbpe * cbpe * nfa;
            const GMTally1 cj = env->sparse() ? vs_r->count(j) :
                                                cbpe * cbpe * nfa;
            Tally2x2<MF> ttable = Metrics_elt_2<Tally2x2<MF>>(*metrics,
              i, j, j_block, *env);
            metrics->dm->histograms()->add(ttable, si1, sj1, ci, cj, nfa);
          } // if is_computing_histograms
        } // for i
      }   // for j
    } // if (!env->is_threshold_tc())

    // Update counts.
    if (!env->is_shrink()) {
      for (int j = 0; j < nvl; ++j) {
        const int i_max = do_compute_triang_only ? j : nvl;
        metrics->num_metric_items_local_computed_inc(i_max);
      }   // for j
    }

  // --------------
  } // if (env->all2all() && do_compute_triang_only)
  // --------------
}

//=============================================================================
// Combine nums and denoms on CPU to get final result, 2-way generic, templated.

template<int METRIC_FORMAT>
static void finalize_(
  GMMetrics* metrics,
  CompressedBuf* matB_cbuf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && matB_cbuf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  // Select action based on metric_type.

  switch (env->metric_type()) {
    case MetricType::CZEK: {

      finalize_czek_<METRIC_FORMAT>(metrics, matB_cbuf,
         vector_sums_left, vector_sums_right,
         j_block, do_compute_triang_only, env);

    } break;
    case MetricType::CCC: {

      finalize_ccc_duo_<METRIC_FORMAT>(metrics, matB_cbuf,
        vector_sums_left, vector_sums_right,
        j_block, do_compute_triang_only, env);

    } break;
    case MetricType::DUO: {

      finalize_ccc_duo_<METRIC_FORMAT>(metrics, matB_cbuf,
        vector_sums_left, vector_sums_right,
        j_block, do_compute_triang_only, env);

    } break;
    default:
      COMET_INSIST_INTERFACE(env, false &&
        "Selected metric_type unimplemented.");
  } // case
}

//=============================================================================
// Combine nums and denoms on CPU to get final result, 2-way generic.

void ComputeMetrics2WayBlock::finalize(
  GMMetrics* metrics,
  CompressedBuf* matB_cbuf,
  const VectorSums* const vector_sums_left,
  const VectorSums* const vector_sums_right,
  int j_block,
  bool do_compute_triang_only,
  CEnv* env) {

  COMET_INSIST(metrics && matB_cbuf);
  COMET_INSIST(vector_sums_left && vector_sums_right && env);
  COMET_INSIST(j_block >= 0 && j_block < env->num_block_vector());
  COMET_INSIST(env->num_way() == NumWay::_2);

  // Select action based on MetricFormat type.

  if (env->is_threshold_tc())

    finalize_<MetricFormat::SINGLE>(metrics, matB_cbuf,
      vector_sums_left, vector_sums_right,
      j_block, do_compute_triang_only, env);

  else // ! env->is_threshold_tc()

    finalize_<MetricFormat::PACKED_DOUBLE>(metrics, matB_cbuf,
      vector_sums_left, vector_sums_right,
      j_block, do_compute_triang_only, env);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
