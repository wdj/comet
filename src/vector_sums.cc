//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.cc
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Per-vector computed quantities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdlib.h"

#include "env.hh"
#include "vectors.hh"
#include "vector_sums.hh"

//=============================================================================
/*---Null object---*/

GMVectorSums GMVectorSums_null(void) {
  GMVectorSums v;
  v.sums = NULL;
  v.counts = NULL;
  v.sums_tmp_ = NULL;
  v.counts_tmp_ = NULL;
  v.size_ = 0;
  //v.num_field_ = 0;
  return v;
}

//=============================================================================
/*---Pseudo-constructor---*/

void GMVectorSums_create(GMVectorSums* this_,
                         int num_vector_local,
                         GMEnv* env) {
  GMInsist(this_ && env);
  GMInsist(num_vector_local >= 0);

  this_->size_ = num_vector_local;
  //this_->num_field_ = vectors->num_field;
  const int num_proc = GMEnv_num_proc_field(env);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      this_->sums = GMFloat_malloc(num_vector_local, env);
      this_->sums_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      this_->sums = GMFloat_malloc(num_vector_local, env);
      this_->sums_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      if (env->sparse) {
        this_->counts = GMFloat_malloc(num_vector_local, env);
        this_->counts_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      }
    } break;
    default:
      GMInsistInterface(env, false && "Unimplemented.");
  } /*---case---*/
}

//=============================================================================
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* this_, GMEnv* env) {
  GMInsist(this_ && env);

  GMFloat_free((GMFloat*)this_->sums, this_->size_, env);
  if (this_->sums_tmp_) {
    GMFloat_free((GMFloat*)this_->sums_tmp_, this_->size_, env);
  }
  if (this_->counts) {
    GMFloat_free((GMFloat*)this_->counts, this_->size_, env);
  }
  if (this_->counts_tmp_) {
    GMFloat_free((GMFloat*)this_->counts_tmp_, this_->size_, env);
  }

  *this_ = GMVectorSums_null();
}

//=============================================================================
/*---Compute the sum of elements of each vector on CPU, for denom---*/

void GMVectorSums_compute_float_(GMVectorSums* this_,
                                 GMVectors* vectors,
                                 GMEnv* env) {
  GMInsist(this_ && vectors && env);

  GMFloat* const __restrict__ sums = this_->sums;
  GMFloat* const __restrict__ sums_tmp = this_->sums_tmp_;

  const int num_proc = GMEnv_num_proc_field(env);
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;

  /*---Sum up all values in each vector---*/

  #pragma omp parallel for schedule(dynamic,1000)
  for (int i = 0; i < vectors->num_vector_local; ++i) {
    GMFloat sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (int f = 0; f < vectors->num_field_local; ++f) {
      const GMFloat value = GMVectors_float_get(vectors, f, i, env);
      sum += value;
    }
    sums_local[i] = sum;
  }

  /*---Do reduction across field procs if needed---*/

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = MPI_Allreduce(sums_local, sums, vectors->num_vector_local,
                             GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
    GMInsist(mpi_code == MPI_SUCCESS);
  }

  env->ops_local += 2 * vectors->num_vector_local *
                    (double)vectors->num_field_local;
}

//-----------------------------------------------------------------------------

void GMVectorSums_compute_bits2_(GMVectorSums* this_,
                                 GMVectors* vectors,
                                 GMEnv* env) {
  GMInsist(this_ && vectors && env);

  GMFloat* const __restrict__ sums = this_->sums;
  GMFloat* const __restrict__ sums_tmp = this_->sums_tmp_;
  GMFloat* const __restrict__ counts = this_->counts;
  GMFloat* const __restrict__ counts_tmp = this_->counts_tmp_;

  const int num_proc = GMEnv_num_proc_field(env);
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;
  GMFloat* const counts_local = num_proc == 1 ? counts : counts_tmp;

  // Count number of 1-bits in each vector

  //----------
  if (env->compute_method_ == GM_COMPUTE_METHOD_REF) {
    //----------
    #pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      if (env->sparse) {
        GMFloat count = 0;
        for (int f = 0; f < vectors->num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 value = GMVectors_bits2_get(vectors, f, i, env);
          if (value != GM_2BIT_UNKNOWN){
            sum += ((value & 1) != 0) + ((value & 2) != 0);
            count++;
          }
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        GMAssert(count >= 0 && count <= vectors->num_field);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { // ! sparse
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors->num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 value = GMVectors_bits2_get(vectors, f, i, env);
          sum += ((value & 1) != 0) + ((value & 2) != 0);
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->num_field);
        sums_local[i] = sum;
      } // if sparse
    } // for i
    //----------
  } else {
    //----------
    for (int i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      if (env->sparse) {
        const GMUInt64 oddbits = 0x5555555555555555;
        GMFloat count = 0;
        for (int f = 0; f < vectors->num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          const GMUInt64 data0 = value.data[0];
          const GMUInt64 data1 = value.data[1];
          const GMUInt64 v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const GMUInt64 v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const GMUInt64 v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const GMUInt64 v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          sum += (GMFloat)gm_popcount64(data0 & v10_mask0);
          sum += (GMFloat)gm_popcount64(data1 & v10_mask1);
          count += (GMFloat)gm_popcount64(v10_oddmask0 | (v10_oddmask1 << 1));
        }
        // Adjust for end pad
        const int num_fields_pad =
          vectors->dm->num_packedfield_local *
          vectors->dm->num_field_per_packedfield -
          vectors->dm->num_field_active_local;
        count -= num_fields_pad;
        // Finish
        GMAssert(sum >= 0 && sum <= 2 * vectors->dm->num_field_active_local);
        GMAssert(count >= 0 && count <= vectors->dm->num_field_active_local);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { // ! sparse
        for (int f = 0; f < vectors->num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          sum += (GMFloat)gm_popcount64(value.data[0]);
          sum += (GMFloat)gm_popcount64(value.data[1]);
          // NOTE: for this case pad entries are all zero so no effect on sum
        }
        GMAssert(sum >= 0 && sum <= 2 * vectors->dm->num_field_active_local);
        sums_local[i] = sum;
      } // if sparse
    } // for i
    //----------
  } // if
  //----------

  // Do reduction across field procs if needed

  if (num_proc > 1) {
    int mpi_code = 0;
    mpi_code = MPI_Allreduce(sums_local, sums, vectors->num_vector_local,
                      GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
    GMInsist(mpi_code == MPI_SUCCESS);
    if (env->sparse) {
      mpi_code = MPI_Allreduce(counts_local, counts, vectors->num_vector_local,
                        GM_MPI_FLOAT, MPI_SUM, GMEnv_mpi_comm_field(env));
      GMInsist(mpi_code == MPI_SUCCESS);
    } /*---if sparse---*/
  }
}

//-----------------------------------------------------------------------------

void GMVectorSums_compute(GMVectorSums* this_, GMVectors* vectors, GMEnv* env) {
  GMInsist(this_ && vectors && env);

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
      GMVectorSums_compute_float_(this_, vectors, env);
    } break;
    case GM_METRIC_TYPE_CCC: {
      GMVectorSums_compute_bits2_(this_, vectors, env);
    } break;
    default:
      GMInsistInterface(env, false && "Unimplemented.");
  } /*---case---*/
}

//-----------------------------------------------------------------------------

GMFloat GMVectorSums_sum(const GMVectorSums* this_, int i,  GMEnv* env) {
  GMAssert(this_ && env);
  GMAssert(i >= 0 && (size_t)i < this_->size_);

  return this_->sums[i];
}

//-----------------------------------------------------------------------------

GMFloat GMVectorSums_count(const GMVectorSums* this_, int i,  GMEnv* env) {
  GMAssert(this_ && env);
  GMAssert(i >= 0 && (size_t)i < this_->size_);
  GMAssert(env->sparse);

  //return this_->counts ? this_->counts[i] : this_->num_field_;
  return this_->counts[i];
}

//-----------------------------------------------------------------------------
