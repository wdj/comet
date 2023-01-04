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
/*!
 * \brief Constructor. Allocate only.
 *
 */
VectorSums::VectorSums(size_t num_vector_local, CEnv& env)
  : env_(env)
  , is_allocated_(false)
  , num_vector_local_(num_vector_local)
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
}

//-----------------------------------------------------------------------------
/*!
 * \brief Constructor. Allocate and compute.
 *
 */
VectorSums::VectorSums(const Vectors& vectors, CEnv& env)
  : env_(env)
  , num_vector_local_(vectors.num_vector_local())
  , sums_(env)
  , sums_tmp_(env)
  , counts_(env)
  , counts_tmp_(env) {

  allocate_();
  compute(vectors);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Constructor. Allocate and compute.
 *
 */
VectorSums::~VectorSums() {
  if (!num_vector_local_)
    return;

  COMET_INSIST(is_allocated_);

  is_allocated_ = false;
  // NOTE MirroredBuf members are deallocated automatically.
}

//-----------------------------------------------------------------------------
/*!
 * \brief Allocate required arrays.
 *
 */
void VectorSums::allocate_() {

  if (!num_vector_local_)
    return;

  const size_t sizeof_float = env_.is_double_prec() ?
                              sizeof(double) : sizeof(float);

  sums_.allocate(num_vector_local_, 1, sizeof_float);
  if (env_.do_reduce())
    sums_tmp_.allocate(num_vector_local_, 1, sizeof_float);

  if (need_counts_()) {
    counts_.allocate(num_vector_local_, 1, sizeof_float);
    if (env_.do_reduce())
      counts_tmp_.allocate(num_vector_local_, 1, sizeof_float);
  }

  is_allocated_ = true;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts (on CPU).
 *
 */
void VectorSums::compute(const Vectors& vectors) {
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

  if (!env_.is_double_prec()) {
    if (env_.is_metric_type_bitwise())
      compute_bits2_<float>(vectors);
    else
      compute_float_<float>(vectors);
  } else { // env_.is_double_prec()
    if (env_.is_metric_type_bitwise())
      compute_bits2_<double>(vectors);
    else
      compute_float_<double>(vectors);
  } // if (!env_.is_double_prec())
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts on GPU.
 *
 */
void VectorSums::compute_accel(const Vectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

  if (!env_.is_compute_method_gpu())
    return;

  COMET_INSIST(!env_.do_reduce() && "This case currently not implemented.");

  if (!env_.is_double_prec()) {
    if (env_.is_metric_type_bitwise())
      compute_bits2_accel_<float>(vectors, accel_stream);
    else
      compute_float_accel_<float>(vectors, accel_stream);
  } else { // env_.is_double_prec()
    if (env_.is_metric_type_bitwise())
      compute_bits2_accel_<double>(vectors, accel_stream);
    else
      compute_float_accel_<double>(vectors, accel_stream);
  } // if (!env_.is_double_prec())
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts (on CPU), float case.
 *
 */
template<typename Float_t>
void VectorSums::compute_float_(const Vectors& vectors) {
  COMET_INSIST(env_.is_double_prec() == (8 == sizeof(Float_t)));
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

  // Where to put the (possibly intermediate) results.

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;

  // Sum up all values in each vector.

# pragma omp parallel for schedule(dynamic,1000)
  for (int i = 0; i < num_vector_local_; ++i) {
    Float_t sum = 0;
    //#pragma omp parallel for reduction(+:sum)
    for (int f = 0; f < (int)vectors.dm()->num_field_local; ++f) {
      const Float_t value = vectors.elt_float_const<Float_t>(f, i);
      sum += value;
    }
    elt_ref_<Float_t>(sums_local, i) = sum;
  }

  // Do reduction across field procs if needed.

  if (env_.do_reduce())
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));

  // Increment operations count.

  env_.ops_local_inc(num_vector_local_ *
                     static_cast<double>(vectors.dm()->num_field_local) * 2);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts (on CPU), bits2 case.
 *
 */
template<typename Float_t>
void VectorSums::compute_bits2_(const Vectors& vectors) {
  COMET_INSIST(env_.is_double_prec() == (8 == sizeof(Float_t)));
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

  // Where to put the (possibly intermediate) results.

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;
  MirroredBuf& counts_local = env_.do_reduce() ? counts_tmp_ : counts_;

  // Count num of 1-bits in each vector, and num nonzero entries if needed.
  // NOTE is_cbpe is used because must handle CCC and DUO differently.

  const int cbpe = env_.counted_bits_per_elt();
  const bool is_cbpe_2 = 2 == cbpe;

  //----------
  if (env_.compute_method() == ComputeMethod::REF) {
  //----------

    // Use reference method, slower but simpler.

    // This should be enough parallelism.
#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (need_counts_()) {
        Float_t count = 0;
        //#pragma omp parallel for reduction(+:sum) reduction(+:count)
        for (int f = 0;
             f < static_cast<int>(vectors.dm()->num_field_active_local); ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = vectors.bits2_get(f, i, env_);
          if (GM_2BIT_UNKNOWN != v){
            sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                             : ((v & 1) != 0);
            count++;
          } // for f
        }
        // Finish
        COMET_ASSERT(count >= 0 && count <= vectors.dm()->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_<Float_t>(sums_local, i) = sum;
        elt_ref_<Float_t>(counts_local, i) = count;
      } else { // ! need_counts_()
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0;
             f < static_cast<int>(vectors.dm()->num_field_active_local); ++f) {
          // Slow way: sum each seminibble individually
          const GMBits2 v = vectors.bits2_get(f, i, env_);
          //const GMBits2 v = Vectors_bits2_get(&vectors, f, i, &env_);
          sum += is_cbpe_2 ? ((v & 1) != 0) + ((v & 2) != 0)
                           : ((v & 1) != 0);
        } // for f
        // Finish
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm()->num_field_active_local);
        elt_ref_<Float_t>(sums_local, i) = sum;
      } // (if need_counts_())
    } // for i

    //----------
  } else { // env_.compute_method() != ComputeMethod::REF
    //----------

    // Use faster method.

    const int num_pad_field_local = vectors.dm()->num_pad_field_local;

    // This should be enough parallelism.
#   pragma omp parallel for schedule(dynamic,1000)
    for (int i = 0; i < num_vector_local_; ++i) {
      Float_t sum = 0;
      if (need_counts_()) {
        const uint64_t oddbits = 0x5555555555555555;
        Float_t count = 0;
        //#pragma omp parallel for reduction(+:sum) reduction(+:count)
        for (int f = 0; f < vectors.num_packedfield_local(); ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = vectors.elt_bits2x64_const(f, i);
          const uint64_t data0 = v.data[0];
          const uint64_t data1 = v.data[1];
          const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
          const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
          const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
          const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);
          // NOTE: for this case, pad entries of vector are all zero,
          // so no effect on sum.
          sum += is_cbpe_2 ? utils::popc64(data0 & v10_mask0)
                           : utils::popc64(data0 & v10_mask0 & oddbits);
          sum += is_cbpe_2 ? utils::popc64(data1 & v10_mask1)
                           : utils::popc64(data1 & v10_mask1 & oddbits);
          // Count number of occurrences of GM_2BIT_UNKNOWN (= binary "10").
          // NOTE: the code below interlaces half the bits of each of the two
          // 64-bit words being processed here.
          // Notice here "count" counts the VECTOR ELEMENTS that are defined,
          // not the number of BITS for all the defined elements.
          count += utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));
        } // for f
        // Adjust for end pad - zero padding entries should not be counted.
        count -= num_pad_field_local;
        // Finish
        COMET_ASSERT(count >= 0 &&
                     count <= vectors.dm()->num_field_active_local);
        COMET_ASSERT(sum >= 0 && sum <= cbpe * count);
        elt_ref_<Float_t>(sums_local, i) = sum;
        elt_ref_<Float_t>(counts_local, i) = count;
      } else { // ! need_counts_()
        const uint64_t oddbits = 0x5555555555555555;
        //#pragma omp parallel for reduction(+:sum)
        for (int f = 0; f < vectors.num_packedfield_local(); ++f) {
          // Fast way: sum all 64 bits of each word immediately
          const GMBits2x64 v = vectors.elt_bits2x64_const(f, i);
          // NOTE: for this case, pad entries of vector are all zero,
          // so no effect on sum.
          sum += is_cbpe_2 ? utils::popc64(v.data[0])
                            : utils::popc64(v.data[0] & oddbits);
          sum += is_cbpe_2 ? utils::popc64(v.data[1])
                            : utils::popc64(v.data[1] & oddbits);
        } // for f
        // Finish
        COMET_ASSERT(sum >= 0 &&
                     sum <= cbpe * vectors.dm()->num_field_active_local);
        elt_ref_<Float_t>(sums_local, i) = sum;
      } // if (need_counts_())
    } // for i

    //----------
  } // if (env_.compute_method() == ComputeMethod::REF)
  //----------

  // Do reduction across field procs if needed.

  if (env_.do_reduce()) {
    COMET_MPI_SAFE_CALL(MPI_Allreduce(sums_local.h, sums_.h,
      num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
    if (need_counts_())
      COMET_MPI_SAFE_CALL(MPI_Allreduce(counts_local.h, counts_.h,
        num_vector_local_, COMET_MPI_FLOAT, MPI_SUM, env_.comm_field()));
  }
}

//-----------------------------------------------------------------------------
/*!
 * \brief File-scope GPU kernel to support float case.
 *
 */
template<typename Float_t>
__global__
static void VectorSums_compute_float_accel_kernel_(
  Float_t* sums,
  Float_t* vectors_data,
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

  // NOTE: each vector has a GPU parallel reduction applied to compute the resut.
  // Adapted from https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

  typedef int Dummy_t;
  extern __shared__ Dummy_t sdata_vector_sums[];
  Float_t* sdata = reinterpret_cast<Float_t*>(&(sdata_vector_sums[0]));

  sdata[pfl_ind0] = static_cast<Float_t>(0);

  const int vl = vl_thread;
  int pfl = pfl_thread;

  // Serially reduce across (threadblock+grid)-sized blocks along pfl direction.

  while (pfl < npfl_thread) {

    const Float_t v = vectors_data[Vectors::index(pfl, vl, npfl)];

    const Float_t sum = v;

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

  if (0 == pfl_ind0)
    atomicAdd(&(sums[vl]), sdata[0]);

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------
/*!
 * \brief File-scope GPU kernel to support bits2 case.
 *
 */
template<typename Float_t>
__global__
static void VectorSums_compute_bits2_accel_kernel_(
  Float_t* sums,
  Float_t* counts,
  GMBits2x64* vectors_data,
  int nvl_thread,
  int npfl_thread,
  int cbpe,
  int num_pad_field_local,
  bool need_counts) {

#ifdef COMET_USE_ACCEL

  // NOTE: no compelling need to support REF case here.

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

  // NOTE: each vector has a GPU parallel reduction applied to compute the resut.
  // Adapted from https://www.olcf.ornl.gov/wp-content/uploads/2019/12/05_Atomics_Reductions_Warp_Shuffle.pdf

  typedef int Dummy_t;
  extern __shared__ Dummy_t sdata_vector_sums[];
  Float_t* sdata = reinterpret_cast<Float_t*>(&(sdata_vector_sums[0]));

  Float_t* sdata_sums = &(sdata[0]);
  Float_t* sdata_counts = sdata_sums + pfl_dim0;

  sdata_sums[pfl_ind0] = static_cast<Float_t>(0);
  if (need_counts)
    sdata_counts[pfl_ind0] = static_cast<Float_t>(0);

  const int vl = vl_thread;
  int pfl = pfl_thread;

  const uint64_t oddbits = 0x5555555555555555;

  // Serially reduce across (threadblock+grid)-sized blocks along pfl direction.

  while (pfl < npfl_thread) {

    const GMBits2x64 v = vectors_data[Vectors::index(pfl, vl, npfl)];

    if (need_counts) {

      const uint64_t data0 = v.data[0];
      const uint64_t data1 = v.data[1];
      const uint64_t v10_oddmask0 = (data0 | ~(data0 >> 1)) & oddbits;
      const uint64_t v10_oddmask1 = (data1 | ~(data1 >> 1)) & oddbits;
      const uint64_t v10_mask0 = v10_oddmask0 | (v10_oddmask0 << 1);
      const uint64_t v10_mask1 = v10_oddmask1 | (v10_oddmask1 << 1);

      const Float_t sum = is_cbpe_2 ?
          utils::popc64(data0 & v10_mask0) +
          utils::popc64(data1 & v10_mask1) :
          utils::popc64(data0 & v10_mask0 & oddbits) +
          utils::popc64(data1 & v10_mask1 & oddbits);

      const Float_t count = utils::popc64(v10_oddmask0 | (v10_oddmask1 << 1));

      sdata_sums[pfl_ind0] += sum;
      sdata_counts[pfl_ind0] += count;
      pfl += pfl_dim0 * pfl_dim1;

    } else { // ! need_counts

      const Float_t sum = is_cbpe_2 ?
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

  if (0 == pfl_ind0) {
    atomicAdd(&(sums[vl]), sdata_sums[0]);
    if (need_counts)
      atomicAdd(&(counts[vl]), sdata_counts[0]);
  }

  // Adjust for pad if needed.

  if (0 == pfl_thread && need_counts)
    atomicAdd(&(counts[vl]), -num_pad_field_local);

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts on GPU, float case.
 *
 */
template<typename Float_t> 
void VectorSums::compute_float_accel_(const Vectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(env_.is_compute_method_gpu());
  COMET_INSIST(env_.is_double_prec() == (8 == sizeof(Float_t)));
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

#ifdef COMET_USE_ACCEL

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;

  // Sum up all values in each vector.

  const int nvl_thread = num_vector_local_;
  const int npfl_thread = vectors.num_packedfield_local();

  // Initialize GPU result arrays to (integer) zero.

# if defined COMET_USE_CUDA
    cudaMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# elif defined COMET_USE_HIP
    hipMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# endif

  // Set up kernel launch.

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
    reinterpret_cast<Float_t*>(sums_local.d),
    reinterpret_cast<Float_t*>(vectors.buf()->d),
    nvl_thread,
    npfl_thread);

  System::accel_last_call_succeeded();

#endif // COMET_USE_ACCEL
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute the sums/counts on GPU, bits2 case.
 *
 */
template<typename Float_t> 
void VectorSums::compute_bits2_accel_(const Vectors& vectors,
  AccelStream_t accel_stream) {
  COMET_INSIST(env_.is_compute_method_gpu());
  COMET_INSIST(env_.is_double_prec() == (8 == sizeof(Float_t)));
  COMET_INSIST(vectors.num_vector_local() == num_vector_local_);

#ifdef COMET_USE_ACCEL

  MirroredBuf& sums_local = env_.do_reduce() ? sums_tmp_ : sums_;
  MirroredBuf& counts_local = env_.do_reduce() ? counts_tmp_ : counts_;

  // Count number of 1-bits in each vector

  const int cbpe = env_.counted_bits_per_elt();

  const int num_pad_field_local = vectors.dm()->num_pad_field_local;
  const bool need_counts = need_counts_();

  const int nvl_thread = num_vector_local_;
  const int npfl_thread = vectors.num_packedfield_local();

  // Initialize GPU result arrays to (integer) zero.

# if defined COMET_USE_CUDA
    cudaMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
    if (need_counts)
      cudaMemsetAsync(counts_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# elif defined COMET_USE_HIP
    hipMemsetAsync(sums_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
    if (need_counts)
      hipMemsetAsync(counts_local.d, 0, nvl_thread*sizeof(Float_t), accel_stream);
# endif

  // Set up kernel launch.

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
    reinterpret_cast<Float_t*>(sums_local.d),
    reinterpret_cast<Float_t*>(counts_local.d),
    reinterpret_cast<GMBits2x64*>(vectors.buf()->d),
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
