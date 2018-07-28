//-----------------------------------------------------------------------------
/*!
 * \file   checksums.cc
 * \author Wayne Joubert
 * \date   Mon Aug  7 14:47:01 EDT 2017
 * \brief  Checksums for metrics.
 * \note   Copyright (C) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "mpi.h"

#include "env.hh"
#include "metrics.hh"
#include "checksums.hh"

//=============================================================================
// Set to null

GMChecksum GMChecksum_null() {
  GMChecksum result;
  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result.data[i] = 0;
  }
  result.is_overflowed = false;
  result.value_max = -DBL_MAX;
  GMMultiprecInt sum = {0};
  result.sum = sum;
  result.sum_d = 0;
  result.computing_checksum = true;
  return result;
}

//-----------------------------------------------------------------------------
// Check whether two checksums equal

bool GMChecksum_equal(GMChecksum* cksum1, GMChecksum* cksum2) {
  if ( ! cksum1->computing_checksum || ! cksum2->computing_checksum ) {
    return true;
  }
  bool result = true;
  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    result = result && cksum1->data[i] == cksum2->data[i];
  }
  result = result && cksum1->is_overflowed == cksum2->is_overflowed;
  result = result && cksum1->value_max == cksum2->value_max;
  return result;
}

//-----------------------------------------------------------------------------
// Helper function - perform one bubble sort step

static void makegreater(size_t& i, size_t& j, int& ind_i, int& ind_j) {
  if (i < j) {
    const size_t tmp = i;
    i = j;
    j = tmp;
    const int tmp2 = ind_i;
    ind_i = ind_j;
    ind_j = tmp2;
  }
}

//-----------------------------------------------------------------------------
// Helper function - left shift that works for any shift amount

static size_t lshift(size_t a, int j) {
  if (j >= 64 || j <= -64) {
    return 0;
  }
  return j > 0 ? a << j : a >> (-j);
}

//-----------------------------------------------------------------------------
// Compute checksum of metrics object

void GMChecksum_metrics(GMChecksum* cksum, GMChecksum* cksum_local,
                        GMMetrics* metrics, GMEnv* env) {
  GMInsist(cksum && cksum_local && metrics && env);
  GMInsist(metrics->data || ! GMEnv_is_proc_active(env));

  //--------------------
  // Initializations
  //--------------------

  if (! GMEnv_is_proc_active(env)) {
    return;
  }

  enum { NUM_WAY_MAX = GM_NUM_NUM_WAY + 1 };
  GMInsist(GMEnv_num_way(env) <= NUM_WAY_MAX && "This num_way not supported.");

  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  //--------------------
  // Check for NaNs if appropriate
  //--------------------

  switch (metrics->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat_check((GMFloat*)(metrics->data), metrics->num_elts_local);
    } break;
  }

  //--------------------
  // Calculate the global largest value
  //--------------------

  double value_max_this = cksum->value_max;
  #pragma omp parallel for reduction(max:value_max_this)
  for (UI64 index = 0; index < metrics->num_elts_local; ++index) {
    bool is_active = true;
    for (int i = 0; i < GMEnv_num_way(env); ++i) {
      const UI64 coord = GMMetrics_coord_global_from_index(metrics, index,
                                                           i, env);
      is_active = is_active && coord < metrics->num_vector_active;
    }
    // Loop over data values at this index
    for (int i_value = 0; i_value < metrics->data_type_num_values; ++i_value) {
      // Pick up value of this metrics elt
      double value = 0;
      switch (metrics->data_type_id) {
        // --------------
        case GM_DATA_TYPE_FLOAT: {
          value = GMMetrics_czek_get_from_index(metrics, index, env);
        } break;
        // --------------
        case GM_DATA_TYPE_TALLY2X2: {
          const int i0 = i_value / 2;
          const int i1 = i_value % 2;
          value = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
        } break;
        // --------------
        case GM_DATA_TYPE_TALLY4X2: {
          const int i0 = i_value / 4;
          const int i1 = (i_value / 2) % 2;
          const int i2 = i_value % 2;
          value =
              GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2, env);
        } break;
        // --------------
        default:
          GMInsist(false && "Invalid data type.");
      } // switch
      if (is_active) {
        value_max_this = value > value_max_this ? value : value_max_this;
      }
    }
  } // for index

  int mpi_code = 0;
  mpi_code = MPI_Allreduce(&value_max_this, &cksum->value_max, 1,
                           MPI_DOUBLE, MPI_MAX, GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  cksum_local->value_max = cksum->value_max;

  // The largest we expect any value to be if using "special" inputs
  const int log2_value_max_allowed = 4;
  const double value_max_allowed = 1 << log2_value_max_allowed;

  cksum->is_overflowed = cksum->is_overflowed && cksum->value_max > value_max_allowed;
  cksum_local->is_overflowed = cksum->is_overflowed;

  const double scaling = value_max_allowed;

  //const double scaling = cksum->is_overflowed ? cksum->value_max : value_max_allowed;

  //--------------------
  // Calculate checksum
  //--------------------

  const int w = 30;
  GMInsist(64 - 2 * w >= 4);
  const UI64 one64 = 1;

  const UI64 lomask = (one64 << w) - 1;
  const UI64 lohimask = (one64 << (2 * w)) - 1;

  GMMultiprecInt sum_this = {0};
  double sum_d_this = 0;

#pragma omp parallel
{
  GMMultiprecInt sum_this_private = {0};
  double sum_d_this_private = 0;
  #pragma omp for collapse(2)
  for (UI64 index = 0; index < metrics->num_elts_local; ++index) {
    // Loop over data values at this index
    for (int i_value = 0; i_value < metrics->data_type_num_values; ++i_value) {

      // Obtain global coords of metrics elt
      UI64 coords[NUM_WAY_MAX];
      int ind_coords[NUM_WAY_MAX];
      for (int i = 0; i < NUM_WAY_MAX; ++i) {
        coords[i] = 0;
        ind_coords[i] = i;
      }
      bool is_active = true;
      for (int i = 0; i < GMEnv_num_way(env); ++i) {
        const UI64 coord = GMMetrics_coord_global_from_index(metrics, index,
                                                             i, env);
        is_active = is_active && coord < metrics->num_vector_active;
        coords[i] = coord;
      }
      // Reflect coords by symmetry to get uniform result -
      //   sort into descending order
      makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);
      makegreater(coords[0], coords[1], ind_coords[0], ind_coords[1]);
      makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);

      // Pick up value of this metrics elt
      double value = 0;
      switch (metrics->data_type_id) {
        // --------------
        case GM_DATA_TYPE_FLOAT: {
          value = GMMetrics_czek_get_from_index(metrics, index, env);
        } break;
        // --------------
        case GM_DATA_TYPE_TALLY2X2: {
          const int i0_unpermuted = i_value / 2;
          const int i1_unpermuted = i_value % 2;
          const int i0 = ind_coords[0] == 0 ? i0_unpermuted : i1_unpermuted;
          const int i1 = ind_coords[0] == 0 ? i1_unpermuted : i0_unpermuted;
          value = GMMetrics_ccc_get_from_index_2(metrics, index, i0, i1, env);
        } break;
        // --------------
        case GM_DATA_TYPE_TALLY4X2: {
          const int i0_unpermuted = i_value / 4;
          const int i1_unpermuted = (i_value / 2) % 2;
          const int i2_unpermuted = i_value % 2;
          const int i0 = ind_coords[0] == 0 ? i0_unpermuted :
                         ind_coords[1] == 0 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i1 = ind_coords[0] == 1 ? i0_unpermuted :
                         ind_coords[1] == 1 ? i1_unpermuted :
                                              i2_unpermuted;
          const int i2 = ind_coords[0] == 2 ? i0_unpermuted :
                         ind_coords[1] == 2 ? i1_unpermuted :
                                              i2_unpermuted;
          value =
              GMMetrics_ccc_get_from_index_3(metrics, index, i0, i1, i2, env);
        } break;
        // --------------
        default:
          GMInsist(false && "Invalid data type.");
      } // switch
      // Convert to integer.  Store only 2*w bits max
      UI64 ivalue = (value / scaling) * (one64 << (2 * w));
      // Construct global id for metrics data value
      UI64 uid = coords[0];
      for (int i = 1; i < GMEnv_num_way(env); ++i) {
        uid = uid * metrics->num_vector_active + coords[i];
      }
      uid = uid * metrics->data_type_num_values + i_value;
      // Randomize
      const UI64 rand1 = gm_randomize(uid + 956158765);
      const UI64 rand2 = gm_randomize(uid + 842467637);
      UI64 rand_value = rand1 + gm_randomize_max() * rand2;
      rand_value &= lohimask;
      // Multiply the two values
      const UI64 a = rand_value;
      const UI64 alo = a & lomask;
      const UI64 ahi = a >> w;
      const UI64 b = ivalue;
      const UI64 blo = b & lomask;
      const UI64 bhi = b >> w;
      const UI64 cx = alo * bhi + ahi * blo;
      UI64 clo = alo * blo + ((cx & lomask) << w);
      UI64 chi = ahi * bhi + (cx >> w);
      // (move the carry bits)
      chi += clo >> (2 * w);
      clo &= lohimask;
      const double value_d =
          ivalue * (double)rand_value / ((double)(one64 << (2 * w)));
      if (is_active) {
        sum_d_this_private += value_d; // Reduction
        // Split the product into one-char chunks, accumulate to sums
        for (int i = 0; i < 8; ++i) {
          const UI64 value0 = (clo << (64 - 8 - 8 * i)) >> (64 - 8);
          const UI64 value1 = (chi << (64 - 8 - 8 * i)) >> (64 - 8);
          sum_this_private.data[0 + i] += value0; // Reduction
          sum_this_private.data[8 + i] += value1; // Reduction
        }
      }
    } // for i_value
  }   // for index

  #pragma omp critical
  {
      sum_d_this += sum_d_this_private; // Reduction
      for (int i = 0; i < 8; ++i) {
        sum_this.data[0 + i] += sum_this_private.data[0 + i];// Reduction
        sum_this.data[8 + i] += sum_this_private.data[8 + i];// Reduction
      }
  }
} // omp parallel

  // Global sum

  GMMultiprecInt sum = {0};
  mpi_code = MPI_Allreduce(sum_this.data, sum.data, GM_MULTIPREC_INT_SIZE,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  for (int i = 0; i < GM_MULTIPREC_INT_SIZE; ++i) {
    cksum->sum.data[i] += sum.data[i];
    cksum_local->sum.data[i] += sum_this.data[i];
  }

  // Combine results

  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    cksum->data[i] = 0;
    cksum_local->data[i] = 0;
    for (int j = 0; j < 8; ++j) {
      cksum->data[i] +=
          lshift(cksum->sum.data[0 + j], 8 * j - 2 * w * i) & lohimask;
      cksum->data[i] +=
          lshift(cksum->sum.data[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
      cksum_local->data[i] +=
          lshift(cksum_local->sum.data[0 + j], 8 * j - 2 * w * i) & lohimask;
      cksum_local->data[i] +=
          lshift(cksum_local->sum.data[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
    }
  }
  /// (move the carry bits
  cksum->data[1] += cksum->data[0] >> (2 * w);
  cksum->data[0] &= lohimask;
  cksum->data[2] += cksum->data[1] >> (2 * w);
  cksum->data[1] &= lohimask;
  cksum_local->data[1] += cksum_local->data[0] >> (2 * w);
  cksum_local->data[0] &= lohimask;
  cksum_local->data[2] += cksum_local->data[1] >> (2 * w);
  cksum_local->data[1] &= lohimask;

  //--------------------
  // Check against floating point result
  //--------------------

  double sum_d;
  mpi_code = MPI_Allreduce(&sum_d_this, &sum_d, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  cksum->sum_d += sum_d;
  cksum_local->sum_d += sum_d_this;

  double result_d = cksum->data[0] / ((double)(one64 << (2 * w))) +
                    cksum->data[1] +
                    cksum->data[2] * ((double)(one64 << (2 * w)));
  GMInsist(fabs(cksum->sum_d - result_d) <= cksum->sum_d * 1.e-10);
}

//=============================================================================
// Print ckecksum to stdout.

void GMChecksum_print(GMChecksum* cksum, GMEnv* env) {
  GMInsist(cksum && env);

  for (int i = 0; i < GM_CHECKSUM_SIZE; ++i) {
    printf("%s%li", i == 0 ? "" : "-",
           cksum->data[GM_CHECKSUM_SIZE - 1 - i]);
  }
  if (cksum->is_overflowed) {
    printf("-OVFL");
    printf("-%e", cksum->value_max);
  }
}


//-----------------------------------------------------------------------------
