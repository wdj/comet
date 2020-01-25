//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.cc
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Compute the denominators needed by the methods.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdint"
#include "cstdlib"
#include <iostream> //FIX

#include "env.hh"
#include "vectors.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

VectorSums::VectorSums(size_t num_vector_local, Env& env)
  : env_(env)
  , num_vector_local_(num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
}

//-----------------------------------------------------------------------------

VectorSums::VectorSums(const GMVectors& vectors, Env& env)
  : env_(env)
  , num_vector_local_(vectors.num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
  compute(vectors);
}

//-----------------------------------------------------------------------------

void VectorSums::allocate_() {

  if (!num_vector_local_)
    return;

  sums_.allocate(num_vector_local_, 1, sizeof(Float_t));
  if (env_.do_reduce())
    sums_tmp_.allocate(num_vector_local_, 1, sizeof(Float_t));

  if (env_.sparse() && env_.is_metric_type_bitwise()) {
    counts_.allocate(num_vector_local_, 1, sizeof(Float_t));
    if (env_.do_reduce())
      counts_tmp_.allocate(num_vector_local_, 1, sizeof(Float_t));
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  if (env_.is_metric_type_bitwise()) {
    compute_bits2(vectors);
  } else {
    compute_float(vectors);
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute_float(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  auto sums_h_local = (Float_t* const __restrict__)(env_.do_reduce() ?
       sums_tmp_.h : sums_.h);

  // Sum up all values in each vector.

# pragma omp parallel for schedule(dynamic,1000)
  for (int i = 0; i < num_vector_local_; ++i) {
    GMFloat sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (int f = 0; f < vectors.num_field_local; ++f) {
      const Float_t value = GMVectors_float_get(&vectors, f, i, &env_);
      sum += value;
    }
    sums_h_local[i] = sum;
  }

  // Do reduction across field procs if needed.

  if (env_.do_reduce())
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_h_local, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));

  env_.ops_local_inc(2 * num_vector_local_ * (double)vectors.num_field_local);
}

//-----------------------------------------------------------------------------

void VectorSums::compute_bits2(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  auto sums_h_local = (Float_t* const __restrict__)(env_.do_reduce() ?
       sums_tmp_.h : sums_.h);
  auto counts_h_local = (Float_t* const __restrict__)(env_.do_reduce() ?
       counts_tmp_.h : counts_.h);

  // Count number of 1-bits in each vector

  const bool do_count_2 = env_.metric_type() == MetricType::CCC;
# ifdef COMET_ASSERTIONS_ON
    const int count_2_1 = do_count_2 ? 2 : 1;
# endif

  //----------
  if (env_.compute_method() == ComputeMethod::REF) {
  //----------

#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (env_.sparse()) {
        Float_t count = 0;
        for (int f = 0; f < (int)vectors.dm->num_field_active_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          if (GM_2BIT_UNKNOWN != v){
            sum += do_count_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                              : ((v & 1) != 0);
            count++;
          }
        }
        COMET_ASSERT(sum >= 0 && sum <= count_2_1 * count);
        COMET_ASSERT(count >= 0 &&
                     count <= vectors.dm->num_field_active_local);
        sums_h_local[i] = sum;
        counts_h_local[i] = count;
      } else { // ! sparse
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors.num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          sum += do_count_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                            : ((v & 1) != 0);
        }
        COMET_ASSERT(sum >= 0 &&
                     sum <= count_2_1 * vectors.dm->num_field_active_local);
        sums_h_local[i] = sum;
      } // if sparse
    } // for i

    //----------
  } else { // REF
    //----------

    //TODO: should decomp_mgr own more of this
    const int num_fields_pad = vectors.dm->num_packedfield_local *
                               vectors.dm->num_field_per_packedfield -
                               vectors.dm->num_field_active_local;
#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (env_.sparse()) {
        const uint64_t oddbits = 0x5555555555555555;
        Float_t count = 0;
        for (int f = 0; f < vectors.num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          const uint64_t data0 = v.data[0];
          const uint64_t data1 = v.data[1];
          const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          sum += do_count_2 ? utils::popc64(data0 & v10_mask0)
                            : utils::popc64(data0 & oddbits & v10_mask0);
          sum += do_count_2 ? utils::popc64(data1 & v10_mask1)
                            : utils::popc64(data1 & oddbits & v10_mask1);
          // NOTE: the code below interlaces half the bits of each of the two
          // 64-bit words being processed here.
          // In fact, "count" counts the VECTOR ELEMENTS that are defined, not
          // the number of BITS for all the defined elements.
          count += utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));
        }
        // Adjust for end pad
        count -= num_fields_pad;
        // Finish
        COMET_ASSERT(sum >= 0 && sum <= count_2_1 * count);
        COMET_ASSERT(count >= 0 &&
                     count <= vectors.dm->num_field_active_local);
        sums_h_local[i] = sum;
        counts_h_local[i] = count;
      } else { // ! sparse
        const uint64_t oddbits = 0x5555555555555555;
        for (int f = 0; f < vectors.num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          sum += do_count_2 ? utils::popc64(v.data[0])
                            : utils::popc64(v.data[0] & oddbits);
          sum += do_count_2 ? utils::popc64(v.data[1])
                            : utils::popc64(v.data[1] & oddbits);
          // NOTE: for this case pad entries are all zero so no effect on sum
        }
        COMET_ASSERT(sum >= 0 &&
                     sum <= count_2_1 * vectors.dm->num_field_active_local);
        sums_h_local[i] = sum;
      } // if sparse
    } // for i

    //----------
  } // if
  //----------

  // Do reduction across field procs if needed

  if (env_.do_reduce()) {
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_h_local, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
    if (env_.sparse())
      COMET_MPI_SAFE_CALL(MPI_Allreduce(counts_h_local, counts_.h,
        num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
  }

}

//-----------------------------------------------------------------------------



#if 0
//-----------------------------------------------------------------------------
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
  COMET_INSIST(this_ && env);
  COMET_INSIST(num_vector_local >= 0);

  this_->size_ = num_vector_local;
  //this_->num_field_ = vectors->num_field;
  const int num_proc = env->num_proc_field();

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      this_->sums = GMFloat_malloc(num_vector_local, env);
      this_->sums_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
    } break;
    case MetricType::CCC: {
      this_->sums = GMFloat_malloc(num_vector_local, env);
      this_->sums_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      if (env->sparse()) {
        this_->counts = GMFloat_malloc(num_vector_local, env);
        this_->counts_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      }
    } break;
    case MetricType::DUO: {
      this_->sums = GMFloat_malloc(num_vector_local, env);
      this_->sums_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      if (env->sparse()) {
        this_->counts = GMFloat_malloc(num_vector_local, env);
        this_->counts_tmp_ = num_proc == 1
                              ? NULL
                              : GMFloat_malloc(num_vector_local, env);
      }
    } break;
    default:
      COMET_INSIST_INTERFACE(env, false && "Unimplemented metric_type.");
  } /*---case---*/
}

//=============================================================================
/*---Pseudo-destructor---*/

void GMVectorSums_destroy(GMVectorSums* this_, GMEnv* env) {
  COMET_INSIST(this_ && env);

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
  COMET_INSIST(this_ && vectors && env);

  GMFloat* const __restrict__ sums = this_->sums;
  GMFloat* const __restrict__ sums_tmp = this_->sums_tmp_;

  const int num_proc = env->num_proc_field();
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;

  /*---Sum up all values in each vector---*/

# pragma omp parallel for schedule(dynamic,1000)
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
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local, sums,
      vectors->num_vector_local, COMET_MPI_FLOAT, MPI_SUM, env->comm_field()));
  }

  env->ops_local_inc(2 * vectors->num_vector_local *
                     (double)vectors->num_field_local);
}

//-----------------------------------------------------------------------------

void GMVectorSums_compute_bits2_(GMVectorSums* this_,
                                 GMVectors* vectors,
                                 GMEnv* env) {
  COMET_INSIST(this_ && vectors && env);

  GMFloat* const __restrict__ sums = this_->sums;
  GMFloat* const __restrict__ sums_tmp = this_->sums_tmp_;
  GMFloat* const __restrict__ counts = this_->counts;
  GMFloat* const __restrict__ counts_tmp = this_->counts_tmp_;

  const int num_proc = env->num_proc_field();
  GMFloat* const sums_local = num_proc == 1 ? sums : sums_tmp;
  GMFloat* const counts_local = num_proc == 1 ? counts : counts_tmp;

  // Count number of 1-bits in each vector

  const bool count_2 = env->metric_type() == MetricType::CCC;

  //----------
  if (env->compute_method() == ComputeMethod::REF) {
    //----------
#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      if (env->sparse()) {
        GMFloat count = 0;
        for (int f = 0; f < vectors->num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 value = GMVectors_bits2_get(vectors, f, i, env);
          if (value != GM_2BIT_UNKNOWN){
            sum += count_2 ? ((value & 1) != 0) + ((value & 2) != 0)
                           : ((value & 1) != 0);
            count++;
          }
        }
        COMET_ASSERT(sum >= 0 && sum <= 2 * vectors->num_field);
        COMET_ASSERT(count >= 0 && count <= vectors->num_field);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { // ! sparse
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors->num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 value = GMVectors_bits2_get(vectors, f, i, env);
          sum += count_2 ? ((value & 1) != 0) + ((value & 2) != 0)
                         : ((value & 1) != 0);
        }
        COMET_ASSERT(sum >= 0 && sum <= (count_2 ? 2 : 1) * vectors->num_field);
        sums_local[i] = sum;
      } // if sparse
    } // for i
    //----------
  } else { // REF
    //----------
    //TODO: should decomp_mgr own more of this
    const int num_fields_pad =
      vectors->dm->num_packedfield_local *
      vectors->dm->num_field_per_packedfield -
      vectors->dm->num_field_active_local;
    for (int i = 0; i < vectors->num_vector_local; ++i) {
      GMFloat sum = 0;
      if (env->sparse()) {
        const uint64_t oddbits = 0x5555555555555555;
        GMFloat count = 0;
        for (int f = 0; f < vectors->num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          const uint64_t data0 = value.data[0];
          const uint64_t data1 = value.data[1];
          const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          sum += count_2 ? (GMFloat)utils::popc64(data0 & v10_mask0)
                         : (GMFloat)utils::popc64(data0 & oddbits & v10_mask0);
          sum += count_2 ? (GMFloat)utils::popc64(data1 & v10_mask1)
                         : (GMFloat)utils::popc64(data1 & oddbits & v10_mask1);
          // NOTE: the code below interlaces half the bits of each of the two
          // 64-bit words being processed here.
          // In fact, "count" counts the VECTOR ELEMENTS that are defined, not
          // the number of BITS for all the defined elements.
          count += (GMFloat)utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));
        }
        // Adjust for end pad
        count -= num_fields_pad;
        // Finish
        COMET_ASSERT(sum >= 0 && sum <= (count_2 ? 2 : 1) * vectors->dm->num_field_active_local);
        COMET_ASSERT(count >= 0 && count <= vectors->dm->num_field_active_local);
        sums_local[i] = sum;
        counts_local[i] = count;
      } else { // ! sparse
        const uint64_t oddbits = 0x5555555555555555;
        for (int f = 0; f < vectors->num_packedval_field_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 value = GMVectors_bits2x64_get(vectors, f, i, env);
          sum += count_2 ? (GMFloat)utils::popc64(value.data[0])
                         : (GMFloat)utils::popc64(value.data[0] & oddbits);
          sum += count_2 ? (GMFloat)utils::popc64(value.data[1])
                         : (GMFloat)utils::popc64(value.data[1] & oddbits);
          // NOTE: for this case pad entries are all zero so no effect on sum
        }
        COMET_ASSERT(sum >= 0 && sum <= (count_2 ? 2 : 1) * vectors->dm->num_field_active_local);
        sums_local[i] = sum;
      } // if sparse
    } // for i
    //----------
  } // if
  //----------

  // Do reduction across field procs if needed

  if (num_proc > 1) {
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local, sums,
      vectors->num_vector_local, COMET_MPI_FLOAT, MPI_SUM, env->comm_field()));
    if (env->sparse()) {
      COMET_MPI_SAFE_CALL(MPI_Allreduce(counts_local, counts,
        vectors->num_vector_local, COMET_MPI_FLOAT, MPI_SUM, env->comm_field()));
    } /*---if sparse---*/
  }
}

//-----------------------------------------------------------------------------

void GMVectorSums_compute(GMVectorSums* this_, GMVectors* vectors, GMEnv* env) {
  COMET_INSIST(this_ && vectors && env);

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      GMVectorSums_compute_float_(this_, vectors, env);
    } break;
    case MetricType::CCC: {
      GMVectorSums_compute_bits2_(this_, vectors, env);
    } break;
    case MetricType::DUO: {
      GMVectorSums_compute_bits2_(this_, vectors, env);
    } break;
    default:
      COMET_INSIST_INTERFACE(env, false && "Unimplemented metric_type.");
  } /*---case---*/
}

//-----------------------------------------------------------------------------

GMFloat GMVectorSums_sum(const GMVectorSums* this_, int i,  GMEnv* env) {
  COMET_ASSERT(this_ && env);
  COMET_ASSERT(i >= 0 && (size_t)i < this_->size_);

  return this_->sums[i];
}

//-----------------------------------------------------------------------------

GMFloat GMVectorSums_count(const GMVectorSums* this_, int i,  GMEnv* env) {
  COMET_ASSERT(this_ && env);
  COMET_ASSERT(i >= 0 && (size_t)i < this_->size_);
  COMET_ASSERT(env->sparse());

  //return this_->counts ? this_->counts[i] : this_->num_field_;
  return this_->counts[i];
}
#endif

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
