//-----------------------------------------------------------------------------
/*!
 * \file   tc_int.hh
 * \author Wayne Joubert
 * \date   Tue May 15 12:03:55 EDT 2018
 * \brief  Declarations needed internally for the tc package.
 * \note   Copyright (C) 2018 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _comet_tc_int_hh_
#define _comet_tc_int_hh_

//=============================================================================

namespace comet {

//=============================================================================
// HELPERS
//=============================================================================

//-----------------------------------------------------------------------------
/// \brief Provide needed constants of GemmIn_t type.
///
///        Provide constants "0", "1" and "2" in the datatype
///        required to input to the reduced precision GEMM.

// Note: we could use __half here instead of uint16_t.  The intent here
// was to use a type based on standard C/C++.  No actual computations
// are done in this code based on the specifics of the type, so it doesn't
// matter.  Important thing is that sizeof(uint16_t) == sizeof(__half) == 2.

template<typename GemmIn_t> struct TCBufTypes;

template<> struct TCBufTypes<uint16_t> {
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

template<> struct TCBufTypes<int8_t> {
  static __host__ __device__ int8_t zero() {return (int8_t)0;}
  static __host__ __device__ int8_t one() {return (int8_t)1;}
  static __host__ __device__ int8_t two() {return (int8_t)2;}
  static __host__ __device__ int8_t four() {return (int8_t)4;}
};

//----------

template<> struct TCBufTypes<GMFp32> {
  static __host__ __device__ GMFp32 zero() {return (GMFp32)0;}
  static __host__ __device__ GMFp32 one() {return (GMFp32)1;}
  static __host__ __device__ GMFp32 two() {return (GMFp32)2;}
  static __host__ __device__ GMFp32 four() {return (GMFp32)4;}
};

//-----------------------------------------------------------------------------
/// \brief Select types etc. based on the setting of the tc param.

template<int TC_METHOD> struct TCSelector;

template<> struct TCSelector<TC::INT8> {
  // types.
  typedef int8_t GemmIn_t;
  typedef int32_t GemmOut_t;
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_8I;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32I;}
#elif defined COMET_USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
   return rocblas_datatype_u8_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
   return rocblas_datatype_u32_r;
  }
#endif
  enum { COUNT = 1 };
};

//----------

template<> struct TCSelector<TC::FP16> {
  // types.
  typedef uint16_t GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_16F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f16_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
  enum { COUNT = 1 };
};

//----------

template<> struct TCSelector<TC::FP32> {
  // types.
  typedef GMFp32 GemmIn_t;
  typedef GMFp32 GemmOut_t;
#if defined COMET_USE_CUDA
  // type selector parameters.
  static cudaDataType __host__ __device__ gemm_type_in() {return CUDA_R_32F;}
  static cudaDataType __host__ __device__ gemm_type_out() {return CUDA_R_32F;}
#elif defined COMET_USE_HIP
  // type selector parameters.
  static rocblas_datatype __host__ __device__ gemm_type_in() {
    return rocblas_datatype_f32_r;
  }
  static rocblas_datatype __host__ __device__ gemm_type_out() {
    return rocblas_datatype_f32_r;
  }
#endif
  enum { COUNT = 1 };
};

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _comet_tc_int_hh_

//-----------------------------------------------------------------------------
