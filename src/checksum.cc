//-----------------------------------------------------------------------------
/*!
 * \file   checksum.cc
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
#include "checksum.hh"

//-----------------------------------------------------------------------------

namespace CoMet {

//=============================================================================
/// \brief Checksum default constructor.

Checksum::Checksum()
  : is_overflowed(false)
  , value_max(-DBL_MAX)
  , is_started(false)
  , computing_checksum(true)
  {}

//-----------------------------------------------------------------------------
/// \brief Check whether two checksums are equal.

bool Checksum::is_equal(const Checksum& cksum2) const {
  bool result = true;

  // Don't perform this test if not computing both checksums.
  if ( ! this->computing_checksum || ! cksum2.computing_checksum ) {
    for (int i = 0; i < CKSUM_SIZE; ++i) {
      result = result && this->data[i] == cksum2.data[i];
    }
    // ISSUE: for now, check overflow like this because
    // overflow is not ncessarily an error (e.g., setting of
    // ccc_multiplier).
    // TODO: fix this better.
    result = result && this->is_overflowed == cksum2.is_overflowed;
    result = result && this->value_max == cksum2.value_max;
  }

  return result;
}

//-----------------------------------------------------------------------------
/// \brief Print ckecksum to stdout.

void Checksum::print(GMEnv& env) {

  for (int i = 0; i < CKSUM_SIZE; ++i) {
    printf("%s%li", i == 0 ? "" : "-",
           this->data[CKSUM_SIZE - 1 - i]);
  }
  if (this->is_overflowed) {
    printf("-OVFL");
    printf("-%e", this->value_max);
  }
}

//-----------------------------------------------------------------------------
/// \brief Checksum helper function: perform one bubble sort step.

inline static void makegreater(size_t& i, size_t& j, int& ind_i, int& ind_j) {
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
/// \brief Checksum helper function: left shift that works for any shift amount.

inline static size_t lshift(size_t a, int j) {
  if (j >= 64 || j <= -64) {
    return 0;
  }
  return j > 0 ? a << j : a >> (-j);
}

//-----------------------------------------------------------------------------
/// \brief Checksum helper: return largest value in metrics object.
///
///        NOTE: this computes the max value on proc, not across procs,
///        so each proc's result in general can have a different value.

double Checksum::metrics_max_value(GMMetrics& metrics, GMEnv& env) {

  double result = -DBL_MAX;

  // TODO: make this unnecessary.
  if (! GMEnv_is_proc_active(&env)) {
    return result;
  }

  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  // Loop over metrics indices to find max.
  #pragma omp parallel for reduction(max:result)
  for (UI64 index = 0; index < metrics.num_elts_local; ++index) {
    // Determine whether this cell is active.
    bool is_active = true;
    for (int i = 0; i < GMEnv_num_way(&env); ++i) {
      const UI64 coord = GMMetrics_coord_global_from_index(&metrics, index,
                                                           i, &env);
      is_active = is_active && coord < metrics.num_vector_active;
    }
    double value_max = -DBL_MAX;
    if (is_active) {
      // Loop over data values at this index
      for (int i_value = 0; i_value < metrics.data_type_num_values; ++i_value) {
        // Pick up value of this metrics elt
        double value = 0;
        switch (metrics.data_type_id) {
          // --------------
          case GM_DATA_TYPE_FLOAT: {
            value = GMMetrics_czek_get_from_index(&metrics, index, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY2X2: {
            const int i0 = i_value / 2;
            const int i1 = i_value % 2;
            value =
              GMMetrics_ccc_get_from_index_2(&metrics, index, i0, i1, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY4X2: {
            const int i0 = i_value / 4;
            const int i1 = (i_value / 2) % 2;
            const int i2 = i_value % 2;
            value =
              GMMetrics_ccc_get_from_index_3(&metrics, index, i0, i1, i2, &env);
          } break;
          // --------------
        default:
          GMInsist(false && "Invalid data type.");
        } // switch
        // value_max is the largest of the values at this index.
        value_max = value > value_max ? value : value_max;
      } // for i_value
    } // if is_active
    result = value_max > result ? value_max : result;
  } // for index

  return result;
} // Checksum::metrics_max_value

//-----------------------------------------------------------------------------
/// \brief compute (global and local) checksum of metrics object.
///
///        The global checksum (across all procs) is generally of most
///        interest.  The local checksum is to debug difficult cases
///        when checksum errors need to be narrowed down to the
///        specific proc.
///
///        Note that cksum and cksum_local are input/output
///        variables; they are added to for multiple stages or phases
///        of the calculation.

void Checksum::compute(Checksum& cksum, Checksum& cksum_local,
                       GMMetrics& metrics, GMEnv& env){
  // TODO: make this check unnecessary.
  GMInsist(metrics.data || ! GMEnv_is_proc_active(&env));

  // TODO: make this unnecessary.
  if (! GMEnv_is_proc_active(&env)) {
    return;
  }

  // TODO: is there a standard C++ type for this.
  typedef size_t UI64;
  GMStaticAssert(sizeof(UI64) == 8);

  enum { NUM_WAY_MAX = GM_NUM_NUM_WAY + 1 };
  GMInsist(GMEnv_num_way(&env) <= NUM_WAY_MAX && "This num_way not supported.");

  // Check for NaNs if appropriate

  // TODO: put this in metrics class - a heavyweight validity check function
  switch (metrics.data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat_check((GMFloat*)(metrics.data), metrics.num_elts_local);
    } break;
  }

  // Get max metrics value.

  // The only impact this has on the output is to determine whether
  // an "overflow" in size of the metric value beyond what the checksum
  // functionality can compute has occurred.
  // It might also be useful for debugging.
  // Note there is no effort to strictly differentiate the
  // local (per proc) vs. global value.

  double value_max_tmp = Checksum::metrics_max_value(metrics, env);
  value_max_tmp = value_max_tmp > cksum.value_max ?
                  value_max_tmp : cksum.value_max;

  int mpi_code = MPI_Allreduce(&value_max_tmp, &cksum.value_max, 1,
                        MPI_DOUBLE, MPI_MAX, GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS);
  cksum_local.value_max = cksum.value_max;

  // Check whether values are within a range for which we can compute
  // the checksum with this code.

  // The largest we expect any value to be if using "special" inputs.
  //
  // This size constraint may be violated if there is an error in the
  // calculation.  It may also be violated by setting of ccc_multiplier large.
  // TODO: examine whether the bound here could be made tighter.
  // TODO: consider polling the metrics object for what the max value
  // upper bound should be - e.g., ccc_multiplier * (1 + roundoff_fuzz)
  const int log2_value_max_allowed = 4;
  const double value_max_allowed = 1 << log2_value_max_allowed;

  cksum.is_overflowed = cksum_local.is_overflowed =
    cksum.is_overflowed && cksum.value_max > value_max_allowed;

  // Scaling factor for values - so that after scaling, value is <= 1.
  const double scaling = value_max_allowed;

  //--------------------
  // Calculate checksum
  //--------------------

  const int w = 30; // 2*w is the integer size
  GMInsist(64 - 2 * w >= 4); // fits into uint64, with some headroom
  const UI64 one64 = 1; // the constant "1"

  const UI64 lomask = (one64 << w) - 1; // masks for lo and hi parts of integer
  const UI64 lohimask = (one64 << (2 * w)) - 1;

  GMMultiprecInt sum_local; // = 0 // checksum valune on this proc
  double sum_d_local = 0; // floating point representation of the same, as check

  #pragma omp parallel
  {
    GMMultiprecInt sum_local_private; // = 0
    double sum_d_local_private = 0;
    // Loop over metrics indices to get checksum contribution.
    #pragma omp for collapse(2)
    for (UI64 index = 0; index < metrics.num_elts_local; ++index) {
      // Loop over data values at this index
      for (int i_value = 0; i_value < metrics.data_type_num_values; ++i_value) {
        // Obtain global coords of metrics elt
        UI64 coords[NUM_WAY_MAX];
        int ind_coords[NUM_WAY_MAX]; // permutation index
        for (int i = 0; i < NUM_WAY_MAX; ++i) {
          coords[i] = 0;
          ind_coords[i] = i;
        }
        bool is_active = true;
        for (int i = 0; i < GMEnv_num_way(&env); ++i) {
          const UI64 coord = GMMetrics_coord_global_from_index(&metrics, index,
                                                               i, &env);
          // Ignore padding vectors.
          is_active = is_active && coord < metrics.num_vector_active;
          coords[i] = coord;
        }
        // Reflect coords by symmetry to get uniform result -
        //   sort into descending order
        //
        // The idea here is that, because the tensor has reflective
        // symmetries, a different equivlent reflected tensor value may be
        // computed based on the parallel decomposition.
        // This permutation puts the indices into a uniform order
        // so that this is not viewed as a difference in the results.
        // Note also below we will permute i0 / i1 / i2 as needed.
        makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);
        makegreater(coords[0], coords[1], ind_coords[0], ind_coords[1]);
        makegreater(coords[1], coords[2], ind_coords[1], ind_coords[2]);

        // Pick up value of this metrics elt
        double value = 0;
        switch (metrics.data_type_id) {
          // --------------
          case GM_DATA_TYPE_FLOAT: {
            value = GMMetrics_czek_get_from_index(&metrics, index, &env);
          } break;
          // --------------
          case GM_DATA_TYPE_TALLY2X2: {
            const int i0_unpermuted = i_value / 2;
            const int i1_unpermuted = i_value % 2;
            const int i0 = ind_coords[0] == 0 ? i0_unpermuted : i1_unpermuted;
            const int i1 = ind_coords[0] == 0 ? i1_unpermuted : i0_unpermuted;
            value =
              GMMetrics_ccc_get_from_index_2(&metrics, index, i0, i1, &env);
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
              GMMetrics_ccc_get_from_index_3(&metrics, index, i0, i1, i2, &env);
          } break;
          // --------------
          default:
            GMInsist(false && "Invalid data type.");
        } // switch
        // Convert to uint64.  Store only 2*w+1 bits, at most -
        // if (value / scaling) <= 1, which it should be if
        // floating point arithmetic works as expected,
        // must have ivalue <= (1<<(2*w)).
        // Note: it would have been better to set this to be 2*w bits max
        // instead - would need to subtract 1 from the second multiplicand
        // and make sure floating point arith works as expected.
        // HOWEVER, see note below.
        UI64 ivalue = (UI64)( (value / scaling) * (one64 << (2 * w)) );
        // Construct an id that is a single number representing the coord
        // and value number.
        UI64 uid = coords[0];
        for (int i = 1; i < GMEnv_num_way(&env); ++i) {
          uid = uid * metrics.num_vector_active + coords[i];
        }
        uid = uid * metrics.data_type_num_values + i_value;
        // Randomize this id
        const UI64 rand1 = gm_randomize(uid + 956158765);
        const UI64 rand2 = gm_randomize(uid + 842467637);
        UI64 rand_value = rand1 + gm_randomize_max() * rand2;
        // Truncate to 2*w bits.
        rand_value &= lohimask;
        // Multiply the two values.
        const UI64 a = rand_value;
        const UI64 alo = a & lomask;
        const UI64 ahi = a >> w;
        const UI64 b = ivalue;
        const UI64 blo = b & lomask;
        const UI64 bhi = b >> w;
        const UI64 cx = alo * bhi + ahi * blo;
        // Note: since a < (1<<(2*w)) and b <= (1<<(2*w)),
        // it is guaranteed that c < (1<<(4*w)),
        // so the result is in the right range of bits.
        UI64 clo = alo * blo + ((cx & lomask) << w);
        UI64 chi = ahi * bhi + (cx >> w);
        // (move the carry bits)
        chi += clo >> (2 * w);
        clo &= lohimask;
        // The checksum is the product of the metric value and the
        // global coord / value number id, this summed across all.
        const double value_d =
            ivalue * (double)rand_value / ((double)(one64 << (2 * w)));
        // Accumulate
        if (is_active) {
          sum_d_local_private += value_d; // (private) reduction
          // Split the product into one-char chunks, accumulate to sums
          for (int i = 0; i < 8; ++i) {
            const UI64 value0 = (clo << (64 - 8 - 8 * i)) >> (64 - 8);
            const UI64 value1 = (chi << (64 - 8 - 8 * i)) >> (64 - 8);
            sum_local_private.data[0 + i] += value0; // (private) reduction
            sum_local_private.data[8 + i] += value1; // (private) reduction
          }
        }
      } // for i_value
    } // for index
    // omp for collapse

    #pragma omp critical
    {
        sum_d_local += sum_d_local_private; // (critical) reduction
        for (int i = 0; i < 8; ++i) {
          sum_local.data[0 + i] += sum_local_private.data[0 + i]; // (critical)
          sum_local.data[8 + i] += sum_local_private.data[8 + i]; // reduction
        }
    }
  } // omp parallel

  // Global sum of multiprecision int

  GMMultiprecInt sum; // = 0
  mpi_code = MPI_Allreduce(sum_local.data, sum.data,
                           MultiprecInt::SIZE,
                           MPI_UNSIGNED_LONG_LONG, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS);
  // Add to multiprec int data we have so far.
  for (int i = 0; i < MultiprecInt::SIZE; ++i) {
    cksum.sum.data[i] += sum.data[i];
    cksum_local.sum.data[i] += sum_local.data[i];
  }

  // Condense results into smaller number of 64 bit ints.

  for (int i = 0; i < CKSUM_SIZE; ++i) {
    cksum.data[i] = 0;
    cksum_local.data[i] = 0;
    for (int j = 0; j < 8; ++j) {
      cksum.data[i] +=
        lshift(cksum.sum.data[0 + j], 8 * j - 2 * w * i) & lohimask;
      cksum.data[i] +=
        lshift(cksum.sum.data[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
      cksum_local.data[i] +=
        lshift(cksum_local.sum.data[0 + j], 8 * j - 2 * w * i) & lohimask;
      cksum_local.data[i] +=
        lshift(cksum_local.sum.data[8 + j], 8 * j - 2 * w * (i - 1)) & lohimask;
    }
  }
  // Adjustments: move the carry bits
  cksum.data[1] += cksum.data[0] >> (2 * w);
  cksum.data[0] &= lohimask;
  cksum.data[2] += cksum.data[1] >> (2 * w);
  cksum.data[1] &= lohimask;
  cksum_local.data[1] += cksum_local.data[0] >> (2 * w);
  cksum_local.data[0] &= lohimask;
  cksum_local.data[2] += cksum_local.data[1] >> (2 * w);
  cksum_local.data[1] &= lohimask;

  // Validate by checking against floating point result

  double sum_d = 0;
  mpi_code = MPI_Allreduce(&sum_d_local, &sum_d, 1, MPI_DOUBLE, MPI_SUM,
                           GMEnv_mpi_comm_repl_vector(&env));
  GMInsist(mpi_code == MPI_SUCCESS);
  cksum.sum_d += sum_d;
  cksum_local.sum_d += sum_d_local;

  double result_d = cksum.data[0] / ((double)(one64 << (2 * w))) +
                    cksum.data[1] +
                    cksum.data[2] * ((double)(one64 << (2 * w)));
  GMInsist(fabs(cksum.sum_d - result_d) <= cksum.sum_d * 1.e-10);
} // Checksum::compute

//=============================================================================

} // namespace CoMet






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
