//-----------------------------------------------------------------------------
/*!
 * \file   tc_helpers.i.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  Declarations needed internally for the tc package.
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

#ifndef _COMET_TC_HELPERS_I_HH_
#define _COMET_TC_HELPERS_I_HH_

//=============================================================================

namespace comet {

//=============================================================================
// HELPERS

//-----------------------------------------------------------------------------
/// \brief Select types etc. based on the setting of the tc param.

struct TCTraitsBase {
  //enum {IS_THREAD_MAPPING_FIELD_MAJOR = BuildHas::HIP ? false : false}; // tuning param
  enum {IS_THREAD_MAPPING_FIELD_MAJOR = BuildHas::HIP ? false : false}; // tuning param
  //enum {NUM_GEMMIN_T_PER_THREAD = BuildHas::HIP ? 64 : 2}; // tuning param
  enum {NUM_GEMMIN_T_PER_THREAD = BuildHas::HIP ? 64 : 2}; // tuning param
  enum {NGIPT = NUM_GEMMIN_T_PER_THREAD};

  enum {NUM_FIELDS_PER_GEMMIN_T_DEFAULT = 1};
  enum {NFPGI = NUM_FIELDS_PER_GEMMIN_T_DEFAULT};

  enum {IS_B_FIELD_MAJOR_DEFAULT = false};
  enum {IS_B_FIELD_MAJOR = IS_B_FIELD_MAJOR_DEFAULT};
};

template<int TC_METHOD> struct TCTraits;

//----------

template<> struct TCTraits<TC::FP32> : public TCTraitsBase {
  //typedef GMFp32 GemmIn_t; // don't use this because harder to access bits.
  typedef uint32_t GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined COMET_USE_CUDA
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_32F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f32_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
};

//----------

template<> struct TCTraits<TC::FP16> : public TCTraitsBase {
  typedef uint16_t GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined COMET_USE_CUDA
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_16F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f16_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
};

//----------

template<> struct TCTraits<TC::INT8> : public TCTraitsBase {
  typedef int8_t GemmIn_t;
  typedef int32_t GemmOut_t;
#if defined COMET_USE_CUDA
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_8I;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32I;}
#elif defined COMET_USE_HIP
  static rocblas_datatype __host__ __device__ gemm_type_in() {
   return rocblas_datatype_i8_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
   return rocblas_datatype_i32_r;
  }
#endif
};

//----------

template<> struct TCTraits<TC::B1> : public TCTraitsBase {
  typedef int8_t GemmIn_t;
  typedef int32_t GemmOut_t;
#if defined COMET_USE_CUDA
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_8I;} // UNUSED
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32I;}
#elif defined COMET_USE_HIP
  static rocblas_datatype __host__ __device__ gemm_type_in() {
   return rocblas_datatype_u8_r; // UNUSED
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
   return rocblas_datatype_u32_r;
  }
#endif
  enum {NFPGI = 8};
  enum {IS_B_FIELD_MAJOR = true};
};

//-----------------------------------------------------------------------------
/// \brief Provide needed constants of GemmIn_t type.
///
///        Provide constants "0", "1". "2", "4" in the datatype
///        required to input to the reduced precision GEMM.

// Note: we could use __half here instead of uint16_t.  The intent here
// was to use a type based on standard C/C++.  No actual computations
// are done in this code based on the specifics of the type, so it doesn't
// matter.  Important thing is that sizeof(uint16_t) == sizeof(__half) == 2.

template<typename GemmIn_t> struct TCBufTraits;

//----------

template<> struct TCBufTraits<uint32_t> {
private:
  //static __host__ __device__ uint32_t mycast(GMFp32 v) {
  //  // ISSUE: the result here is, "by the book", undefined; the proper
  //  // way is to use char* to access the component bytes of the types
  //  // to copy from one type to another.
  //  const uint32_t* const p = (uint32_t*)&v;
  //  return *p;
  //}
  static __host__ __device__ uint32_t mycast(const GMFp32 v) {
    const char* const pv = (const char*)&v;
    uint32_t r = 0;
    char* const pr = (char*)&r;
    pr[0] = pv[0];
    pr[1] = pv[1];
    pr[2] = pv[2];
    pr[3] = pv[3];
    return r;
  }
public:
  static __host__ __device__ uint32_t zero() {return mycast(0);}
  static __host__ __device__ uint32_t one() {return mycast(1);}
  static __host__ __device__ uint32_t two() {return mycast(2);}
  static __host__ __device__ uint32_t four() {return mycast(4);}
};

//----------

template<> struct TCBufTraits<uint16_t> {
  static __host__ __device__ uint16_t zero() {return (uint16_t)0x0000;}
  static __host__ __device__ uint16_t one() {return (uint16_t)0x3c00;}
  static __host__ __device__ uint16_t two() {return (uint16_t)0x4000;}
  static __host__ __device__ uint16_t four() {return (uint16_t)0x4400;}
                                        // = *(uint16_t*)&__float2half(0.);
                                        // = *(uint16_t*)&__float2half(1.);
                                        // = *(uint16_t*)&__float2half(2.);
                                        // = *(uint16_t*)&__float2half(4.);
};

//----------

template<> struct TCBufTraits<int8_t> {
  static __host__ __device__ int8_t zero() {return (int8_t)0;}
  static __host__ __device__ int8_t one() {return (int8_t)1;}
  static __host__ __device__ int8_t two() {return (int8_t)2;}
  static __host__ __device__ int8_t four() {return (int8_t)4;}
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_HELPERS_I_HH_

//-----------------------------------------------------------------------------
