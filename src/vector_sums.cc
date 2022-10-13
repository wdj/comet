//-----------------------------------------------------------------------------
/*!
 * \file   vector_sums.cc
 * \author Wayne Joubert
 * \date   Sat Oct 24 10:48:32 EDT 2015
 * \brief  Compute the denominators needed by the methods.
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
#include "cstdlib"

#include "env.hh"
#include "vectors.hh"
#include "vector_sums.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------

VectorSums::VectorSums(size_t num_vector_local, CEnv& env)
  : env_(env)
  , num_vector_local_(num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
}

//-----------------------------------------------------------------------------

VectorSums::VectorSums(const GMVectors& vectors, CEnv& env)
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

  if (need_counts_()) {
    counts_.allocate(num_vector_local_, 1, sizeof(Float_t));
    if (env_.do_reduce())
      counts_tmp_.allocate(num_vector_local_, 1, sizeof(Float_t));
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  // Compute the vector sums.

  if (env_.is_metric_type_bitwise()) {
    compute_bits2_(vectors);
  } else {
    compute_float_(vectors);
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute_accel(const GMVectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  if (!env_.is_compute_method_gpu())
    return;

  COMET_INSIST(!env_.do_reduce() && "Not implemented.");

  if (env_.is_metric_type_bitwise()) {
    compute_bits2_accel_(vectors, accel_stream);
  } else {
    compute_float_accel_(vectors, accel_stream);
  }
}

//-----------------------------------------------------------------------------

void VectorSums::compute_float_(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;

  // Sum up all values in each vector.

# pragma omp parallel for schedule(dynamic,1000)
  for (int i = 0; i < num_vector_local_; ++i) {
    GMFloat sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (int f = 0; f < vectors.num_field_local; ++f) {
      const Float_t value = GMVectors_float_get(&vectors, f, i, &env_);
      sum += value;
    }
    elt_ref_(sums_local, i) = sum;
  }

  // Do reduction across field procs if needed.

  if (env_.do_reduce())
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));

  env_.ops_local_inc(2 * num_vector_local_ * (double)vectors.num_field_local);
}

//-----------------------------------------------------------------------------

void VectorSums::compute_bits2_(const GMVectors& vectors) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;
  MirroredBuf& counts_local = env_.do_reduce() ? counts_tmp_ : counts_;

  // Count number of 1-bits in each vector

  const int cbpe = env_.counted_bits_per_elt();
  const bool is_cbpe_2 = cbpe == 2;

  //----------
  if (env_.compute_method() == ComputeMethod::REF) {
  //----------

#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (need_counts_()) {
        Float_t count = 0;
        for (int f = 0; f < (int)vectors.dm()->num_field_active_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          if (GM_2BIT_UNKNOWN != v){
            sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                             : ((v & 1) != 0);
            count++;
          } // for f
        }
        COMET_ASSERT(count >= 0 && count <= vectors.dm()->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_(sums_local, i) = sum;
        elt_ref_(counts_local, i) = count;
      } else { // ! need_counts_()
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors.num_field_local; ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = GMVectors_bits2_get(&vectors, f, i, &env_);
          sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                           : ((v & 1) != 0);
        } // for f
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm()->num_field_active_local);
        elt_ref_(sums_local, i) = sum;
      } // if need_counts_()
    } // for i

    //----------
  } else { // env_.compute_method() != ComputeMethod::REF
    //----------

    const int num_pad_field_local = vectors.dm()->num_pad_field_local;

#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (need_counts_()) {
        const uint64_t oddbits = 0x5555555555555555;
        Float_t count = 0;
        for (int f = 0; f < vectors.num_packedfield_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          const uint64_t data0 = v.data[0];
          const uint64_t data1 = v.data[1];
          const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          sum += is_cbpe_2 ? utils::popc64(data0 & v10_mask0)
                           : utils::popc64(data0 & v10_mask0 & oddbits);
          sum += is_cbpe_2 ? utils::popc64(data1 & v10_mask1)
                           : utils::popc64(data1 & v10_mask1 & oddbits);
          // NOTE: the code below interlaces half the bits of each of the two
          // 64-bit words being processed here.
          // In fact, "count" counts the VECTOR ELEMENTS that are defined, not
          // the number of BITS for all the defined elements.
          count += utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));
        } // for f
        // Adjust for end pad
        count -= num_pad_field_local;
        // Finish
        COMET_ASSERT(count >= 0 && count <= vectors.dm()->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_(sums_local, i) = sum;
        elt_ref_(counts_local, i) = count;
//printf("1111 %f %f\n", sum, count);
      } else { // ! need_counts_()
        const uint64_t oddbits = 0x5555555555555555;
        for (int f = 0; f < vectors.num_packedfield_local; ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = GMVectors_bits2x64_get(&vectors, f, i, &env_);
          sum += is_cbpe_2 ? utils::popc64(v.data[0])
                            : utils::popc64(v.data[0] & oddbits);
          sum += is_cbpe_2 ? utils::popc64(v.data[1])
                            : utils::popc64(v.data[1] & oddbits);
          // NOTE: for this case pad entries of vec all zero so no effect on sum
        } // for f
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm()->num_field_active_local);
        elt_ref_(sums_local, i) = sum;
      } // if need_counts_()
    } // for i

    //----------
  } // if
  //----------

  // Do reduction across field procs if needed

  if (env_.do_reduce()) {
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
    if (need_counts_())
      COMET_MPI_SAFE_CALL(MPI_Allreduce(counts_local.h, counts_.h,
        num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
  }
}

//-----------------------------------------------------------------------------

template<typename Out_t>
__global__
static void VectorSums_compute_float_accel_kernel_(
  Out_t* sums,
  VectorSums::Float_t* vectors_data,
  int nvl_thread,
  int npfl_thread) {
#ifdef COMET_USE_ACCEL

  //const int nvl = nvl_thread;
  const int npfl = npfl_thread;

  const int vl_thread = blockIdx_y_() + gridDim_y_() * blockIdx_z_();
  if (vl_thread >= nvl_thread)
    return;

  const int pfl_ind0 = threadIdx_x_();
  const int pfl_dim0 = blockDim_x_();
  const int pfl_ind1 = blockIdx_x_();
  const int pfl_dim1 = gridDim_x_();

  const int pfl_thread = pfl_ind0 + pfl_dim0 * pfl_ind1;

  // Adapted from https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

  extern __shared__ Out_t sdata_vector_sums_1[];

  Out_t* sdata = &(sdata_vector_sums_1[0]);

  sdata[pfl_ind0] = (Out_t)0;

  const int vl = vl_thread;
  int pfl = pfl_thread;

  // Serially reduce across (threadblock+grid)-sized blocks along pfl direction.

  while (pfl < npfl_thread) {

    const VectorSums::Float_t v = vectors_data[GMVectors_index(pfl, vl, npfl)];
//printf("%f\n", (double)v);

    const Out_t sum = v;

    sdata[pfl_ind0] += sum;
    pfl += pfl_dim0 * pfl_dim1;

  } // while

  // Reduce within threadblock.

  for (unsigned int s = pfl_dim0/2; s > 0; s >>= 1) {
    __syncthreads();
    if (pfl_ind0 < s) { // parallel sweep reduction
      sdata[pfl_ind0] += sdata[pfl_ind0 + s];
    }
  }

  // First thread of each threadblock adds in its contribution.

  if (pfl_ind0 == 0)
    atomicAdd(&(sums[vl]), sdata[0]);
//if (pfl_ind0 == 0) printf("%f\n", (double)sums[vl]);

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------

template<typename Out_t>
__global__
static void VectorSums_compute_bits2_accel_kernel_(
  Out_t* sums,
  Out_t* counts,
  GMBits2x64* vectors_data,
  int nvl_thread,
  int npfl_thread,
  int cbpe,
  int num_pad_field_local,
  bool need_counts) {
#ifdef COMET_USE_ACCEL

  //const int nvl = nvl_thread;
  const int npfl = npfl_thread;

  const int vl_thread = blockIdx_y_() + gridDim_y_() * blockIdx_z_();
  if (vl_thread >= nvl_thread)
    return;

  const int pfl_ind0 = threadIdx_x_();
  const int pfl_dim0 = blockDim_x_();
  const int pfl_ind1 = blockIdx_x_();
  const int pfl_dim1 = gridDim_x_();

  const int pfl_thread = pfl_ind0 + pfl_dim0 * pfl_ind1;

  const bool is_cbpe_2 = 2 == cbpe;

  // Adapted from https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

  extern __shared__ Out_t sdata_vector_sums_2[];

  Out_t* sdata = &(sdata_vector_sums_2[0]);

  Out_t* sdata_sums = &(sdata[0]);
  Out_t* sdata_counts = sdata_sums + pfl_dim0;

  sdata_sums[pfl_ind0] = (Out_t)0;
  if (need_counts)
    sdata_counts[pfl_ind0] = (Out_t)0;

  const int vl = vl_thread;
  int pfl = pfl_thread;

  const uint64_t oddbits = 0x5555555555555555;

  // Serially reduce across (threadblock+grid)-sized blocks along pfl direction.

  while (pfl < npfl_thread) {

    const GMBits2x64 v = vectors_data[GMVectors_index(pfl, vl, npfl)];

    if (need_counts) {

      const uint64_t data0 = v.data[0];
      const uint64_t data1 = v.data[1];
      const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
      const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
      const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
      const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);

      const Out_t sum = is_cbpe_2 ?
          utils::popc64(data0 & v10_mask0) +
          utils::popc64(data1 & v10_mask1) :
          utils::popc64(data0 & v10_mask0 & oddbits) +
          utils::popc64(data1 & v10_mask1 & oddbits);

      const Out_t count = utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));

      sdata_sums[pfl_ind0] += sum;
      sdata_counts[pfl_ind0] += count;
      pfl += pfl_dim0 * pfl_dim1;

    } else { // if (!need_counts)

      const Out_t sum = is_cbpe_2 ?
          utils::popc64(v.data[0]) +
          utils::popc64(v.data[1]) :
          utils::popc64(v.data[0] & oddbits) +
          utils::popc64(v.data[1] & oddbits);
          // NOTE: for this case pad entries of vec all zero so no effect on sum

      sdata_sums[pfl_ind0] += sum;
      pfl += pfl_dim0 * pfl_dim1;

    } // if (need_counts)

  } // while

  // Reduce within threadblock.

  for (unsigned int s = pfl_dim0/2; s > 0; s >>= 1) {
    __syncthreads();
    if (pfl_ind0 < s) { // parallel sweep reduction
      sdata_sums[pfl_ind0] += sdata_sums[pfl_ind0 + s];
      if (need_counts)
        sdata_counts[pfl_ind0] += sdata_counts[pfl_ind0 + s];
    }
  }

  // First thread of each threadblock adds in its contribution.

  if (pfl_ind0 == 0) {
    atomicAdd(&(sums[vl]), sdata_sums[0]);
    if (need_counts)
      atomicAdd(&(counts[vl]), sdata_counts[0]);
  }

  // Adjust for pad if needed.

  if (pfl_thread == 0 && need_counts)
    atomicAdd(&(counts[vl]), -num_pad_field_local);

//if (pfl_thread < npfl_thread)
//if (vl <= 2) {
//printf("2222s %f  %i %i\n", (double)sums[vl], vl, pfl);
////printf("2222s %f  %i %i  %i\n", (double)sdata_sums[pfl_ind0], vl, pfl, pfl_ind0);
//if (need_counts)
//printf("2222s %f  %i %i\n", (double)counts[vl], vl, pfl);
////printf("2222c %f  %i %i  %i\n", (double)sdata_counts[pfl_ind0], vl, pfl, pfl_ind0);
//}

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------

void VectorSums::compute_float_accel_(const GMVectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  COMET_INSIST(env_.is_compute_method_gpu());

#ifdef COMET_USE_ACCEL

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;

  // Sum up all values in each vector.

  const int nvl_thread = num_vector_local_;
  const int npfl_thread = vectors.num_packedfield_local;

# if defined COMET_USE_CUDA
    cudaMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# elif defined COMET_USE_HIP
    hipMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# endif

  // ISSUE: this may need to be a power of 2.
  const int threadblocksize = 256;
  COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
               "Current HIP limitation.");
  const int blockdim_y = 32768;
  const int num_threadblocks_0 = utils::ceil(npfl_thread, threadblocksize);
  const int num_threadblocks_1 = utils::min(nvl_thread, blockdim_y);
  const int num_threadblocks_2 = utils::ceil(nvl_thread, blockdim_y);

  COMET_LAUNCH_KERNEL((VectorSums_compute_float_accel_kernel_<Float_t>),
    dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
    dim3(threadblocksize, 1, 1),
    threadblocksize * sizeof(Float_t),
    accel_stream,
    (Float_t*)sums_local.d,
    //(Float_t*)vectors.data,
    (Float_t*)vectors.buf()->d,
    nvl_thread,
    npfl_thread);

  System::accel_last_call_succeeded();

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------

void VectorSums::compute_bits2_accel_(const GMVectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(num_vector_local_ == vectors.num_vector_local);

  COMET_INSIST(env_.is_compute_method_gpu());

#ifdef COMET_USE_ACCEL

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;
  MirroredBuf& counts_local = env_.do_reduce() ? counts_tmp_ : counts_;

  // Count number of 1-bits in each vector

  const int cbpe = env_.counted_bits_per_elt();

  const int num_pad_field_local = vectors.dm()->num_pad_field_local;
  const bool need_counts = need_counts_();

  const int nvl_thread = num_vector_local_;
  const int npfl_thread = vectors.num_packedfield_local;

# if defined COMET_USE_CUDA
    cudaMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
    if (need_counts)
      cudaMemsetAsync(counts_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# elif defined COMET_USE_HIP
    hipMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
    if (need_counts)
      hipMemsetAsync(counts_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# endif

  // ISSUE: this may need to be a power of 2.
  const int threadblocksize = 256;
  COMET_INSIST((threadblocksize <= 256 || ! BuildHas::HIP) &&
               "Current HIP limitation.");
  const int blockdim_y = 32768;
  const int num_threadblocks_0 = utils::ceil(npfl_thread, threadblocksize);
  const int num_threadblocks_1 = utils::min(nvl_thread, blockdim_y);
  const int num_threadblocks_2 = utils::ceil(nvl_thread, blockdim_y);

  COMET_LAUNCH_KERNEL((VectorSums_compute_bits2_accel_kernel_<Float_t>),
    dim3(num_threadblocks_0, num_threadblocks_1, num_threadblocks_2),
    dim3(threadblocksize, 1, 1),
    threadblocksize * sizeof(Float_t) * (1+need_counts),
    accel_stream,
    (Float_t*)sums_local.d,
    (Float_t*)counts_local.d,
    //(GMBits2x64*)vectors.data,
    (GMBits2x64*)vectors.buf()->d,
    nvl_thread,
    npfl_thread,
    cbpe,
    num_pad_field_local,
    need_counts);

  System::accel_last_call_succeeded();

#endif // COMET_USE_ACCEL
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
