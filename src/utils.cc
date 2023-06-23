//-----------------------------------------------------------------------------
/*!
 * \file   utils.cc
 * \author Wayne Joubert
 * \date   Sat Nov 16 10:04:31 EST 2019
 * \brief  Miscellaneous utilities
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

#include "cstdio"
#include "cstdlib"
#include "cstddef"
#include "cstring"
#include "math.h"
#include "limits"

#include "errno.h"
#include "sys/time.h"
#include "signal.h"

#include "mpi.h"

#if defined COMET_USE_CUDA
#include "cuda.h"
#endif

#include "env.hh"

//=============================================================================

namespace comet {

namespace System {

//-----------------------------------------------------------------------------
/*!
 * \brief System (wallclock) timer.
 *
 */
double time() {

  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = static_cast<double>(tv.tv_sec) +
                  static_cast<double>(tv.tv_usec) * 1.e-6;
  return result;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Number of processors (MPI ranks) available.
 *
 */
int num_proc() {
  int num_proc = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_size(MPI_COMM_WORLD, &num_proc));
  return num_proc;
}

//-----------------------------------------------------------------------------
/*!
 * \brief MPI rank in comm world.
 *
 */
int proc_num() {
  int proc_num = 0;
  COMET_MPI_SAFE_CALL(MPI_Comm_rank(MPI_COMM_WORLD, &proc_num));
  return proc_num;
}

//-----------------------------------------------------------------------------
// GPU device properties.

#if defined COMET_USE_CUDA
  typedef cudaDeviceProp accelDeviceProp_t;
#elif defined COMET_USE_HIP
  typedef hipDeviceProp_t accelDeviceProp_t;
#else
  typedef int accelDeviceProp_t;
#endif

//-----------------------------------------------------------------------------
/*!
 * \brief Get GPU device properties.
 *
 */
accelDeviceProp_t& get_device_prop_() {
  // NOTE: local static variable.
  static accelDeviceProp_t device_prop = {};
  static bool is_initialized = false;
  if (!is_initialized) {
#if defined COMET_USE_CUDA
    const cudaError_t error = cudaGetDeviceProperties(&device_prop, 0);
    no_unused_variable_warning(error);
#elif defined COMET_USE_HIP
    hipGetDeviceProperties(&device_prop, 0); // Assume only one GPU per rank.
#else
    device_prop = 0;
#endif
    //COMET_INSIST(accel_last_call_succeeded());
    is_initialized = true;
  }
  return device_prop;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accelerator compute capability.
 *
 */
int compute_capability() {
  const accelDeviceProp_t device_prop = get_device_prop_();
#if defined COMET_USE_CUDA
  //// Assume only one GPU per rank.
  //cudaDeviceProp device_prop;
  //cudaError_t error = cudaGetDeviceProperties(&device_prop, 0);
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  const int compute_capability = num_devices == 0 ? 0 :
    device_prop.major * 100 + device_prop.minor;
#elif defined COMET_USE_HIP
  //hipDeviceProp_t device_prop;
  //hipGetDeviceProperties(&device_prop, 0); // Assume only one GPU per rank.
  //const int compute_capability = device_prop.major * 100 + device_prop.minor;
  // This seems more stable than major/minor.
  int num_devices;
  hipGetDeviceCount(&num_devices);
  const int compute_capability = num_devices == 0 ? 0 : device_prop.gcnArch;
#else
  no_unused_variable_warning(device_prop);
  const int compute_capability = 0;
#endif
  COMET_INSIST(accel_last_call_succeeded());
  return compute_capability;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accelerator pci bus id.
 *
 */
int pci_bus_id() {
  const accelDeviceProp_t device_prop = get_device_prop_();
#if defined COMET_USE_CUDA
  //// Assume only one GPU per rank.
  //cudaDeviceProp device_prop;
  //cudaError_t error = cudaGetDeviceProperties(&device_prop, 0);
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  const int pci_bus_id = num_devices == 0 ? 0 : device_prop.pciBusID;
#elif defined COMET_USE_HIP
  //hipDeviceProp_t device_prop;
  //hipGetDeviceProperties(&device_prop, 0); // Assume only one GPU per rank.
  int num_devices;
  hipGetDeviceCount(&num_devices);
  const int pci_bus_id = num_devices == 0 ? 0 :device_prop.pciBusID;
#else
  no_unused_variable_warning(device_prop);
  const int pci_bus_id = 0;
#endif
  COMET_INSIST(accel_last_call_succeeded());
  return pci_bus_id;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accelerator pci domain id.
 *
 */
int pci_domain_id() {
  const accelDeviceProp_t device_prop = get_device_prop_();
#if defined COMET_USE_CUDA
  //// Assume only one GPU per rank.
  //cudaDeviceProp device_prop;
  //cudaError_t error = cudaGetDeviceProperties(&device_prop, 0);
  int num_devices;
  cudaGetDeviceCount(&num_devices);
  const int pci_domain_id = num_devices == 0 ? 0 : device_prop.pciDomainID;
#elif defined COMET_USE_HIP
  //hipDeviceProp_t device_prop;
  //hipGetDeviceProperties(&device_prop, 0); // Assume only one GPU per rank.
  int num_devices;
  hipGetDeviceCount(&num_devices);
  const int pci_domain_id = num_devices == 0 ? 0 : device_prop.pciDomainID;
#else
  no_unused_variable_warning(device_prop);
  const int pci_domain_id = 0;
#endif
  COMET_INSIST(accel_last_call_succeeded());
  return pci_domain_id;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Accelerator did most recent call succeed.
 *
 */
bool accel_last_call_succeeded() {

#if defined COMET_USE_CUDA
  // NOTE: this read of the last error is a destructive read.
  cudaError_t error = cudaGetLastError();
  const bool result = error == cudaSuccess;

  if (!result) {
    fprintf(stderr, "CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
#elif defined COMET_USE_HIP
  // NOTE: this read of the last error is (apparently) a destructive read.
  hipError_t error = hipGetLastError();
  const bool result = error == hipSuccess;

  if (!result) {
    fprintf(stderr, "HIP error detected: %s\n", hipGetErrorString(error));
  }

  return result;
#endif

  return true;
}

} // namespace System

//-----------------------------------------------------------------------------

namespace utils {

//-----------------------------------------------------------------------------
/*!
 * \brief Memory allocation with memory usage tracking.
 *
 */
void* malloc(size_t n, CEnv& env) {
  COMET_INSIST(n+1 >= 1);

  void* p = ::malloc(n);
  COMET_INSIST(p &&
    "Invalid pointer from malloc, possibly due to insufficient memory.");
  env.cpu_mem_local_inc(n);
  return p;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Memory allocation without memory usage tracking.
 *
 */
void* malloc_nometer(size_t n, CEnv& env) {
  COMET_INSIST(n+1 >= 1);

  void* p = ::malloc(n);
  COMET_INSIST(p &&
    "Invalid pointer from malloc, possibly due to insufficient memory.");
  return p;
}

//-----------------------------------------------------------------------------
/*!
 * \brief Memory deallocation with memory usage tracking.
 *
 */
void free(void* p, size_t n, CEnv& env) {
  COMET_INSIST(p);
  COMET_INSIST(n+1 >= 1);

  ::free(p);
  env.cpu_mem_local_dec(n);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Memory deallocation without memory usage tracking.
 *
 */
void free_nometer(void* p, size_t n, CEnv& env) {
  COMET_INSIST(p);
  COMET_INSIST(n+1 >= 1);

  ::free(p);
}

//-----------------------------------------------------------------------------
/*!
 * \brief Compute checksum of a simple array.
 *
 */
size_t array_cksum(const unsigned char* const a, size_t n) {
  COMET_INSIST(a);
  COMET_INSIST(n+1 >= 1);

  size_t result = 0;
  const size_t mask = (((size_t)1) << 32) - 1;

# pragma omp parallel for schedule(dynamic,1000) reduction(+:result)
  for (size_t i=0; i<n; ++i) {
    result += (a[i] * i) & mask;
  }

  return result;
}

//-----------------------------------------------------------------------------

} // namespace utils

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
