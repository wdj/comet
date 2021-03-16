//-----------------------------------------------------------------------------
/*!
 * \file   tc_solve_cutlass_general.i.hh
 * \author Paul Eller
 * \date   Tue Nov  3 08:26:29 EST 2020
 * \brief  CUDA code, gemm operation, 
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

#ifndef _COMET_TC_SOLVE_CUTLASS_GENERAL_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_GENERAL_I_HH_

// Defines cutlass::gemm::device::Gemm, the generic Gemm computation template class.
#include "cutlass/gemm/device/gemm.h"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/// \brief Cutlass GEMM driver routine

template <typename Gemm>
void CutlassGemmRun(int M, int N, int K, cutlass::uint1b_t const *A,
  int lda, cutlass::uint1b_t const *B, int ldb, int32_t beta, int32_t *C,
  int ldc, AccelStream_t accel_stream) {

  // Define a CUTLASS GEMM type
  Gemm gemm_operator;
  int32_t alpha = 1;

  // Construct the CUTLASS GEMM arguments object.
  typename Gemm::Arguments args({M, N, K},      // Gemm Problem dimensions
                                {A, lda},       // Tensor-ref for source matrix A
                                {B, ldb},       // Tensor-ref for source matrix B
                                {C, ldc},       // Tensor-ref for source matrix C
                                {C, ldc},       // Tensor-ref for destination matrix C
                                {alpha,beta});  // Scalars used in the Epilogue

  // Launch the CUTLASS GEMM kernel
  cutlass::Status status = gemm_operator(args, nullptr, accel_stream);
  COMET_INSIST(status == cutlass::Status::kSuccess);
  System::accel_last_call_succeeded();
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 128 x 256 blocks

struct TCTBlockType {
  enum {_128_256_512  = 0,
	_256_128_512  = 1,
	_128_128_512  = 2,
	_128_64_512   = 3,
	_64_128_512   = 4,
	_64_64_512    = 5,
        _128_256_1024 = 6
  };
};

template <int TBType> struct TBlockType;

// Turing Settings
template<> struct TBlockType<TCTBlockType::_128_256_512> {
  enum {t0 = 128, t1 = 256, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_256_128_512> {
  enum {t0 = 256, t1 = 128, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_128_128_512> {
  enum {t0 = 128, t1 = 128, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_128_64_512> {
  enum {t0 = 128, t1 = 64, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_64_128_512> {
  enum {t0 = 64, t1 = 128, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_64_64_512> {
  enum {t0 = 64, t1 = 64, t2 = 512};
};

// Ampere Settings
#if defined COMET_USE_AMPERE
template<> struct TBlockType<TCTBlockType::_128_256_1024> {
  enum {t0 = 128, t1 = 256, t2 = 1024};
};
#endif

struct TCWarpType {
  enum {_64_64_512  = 0,
	_64_32_512  = 1,
	_32_64_512  = 2,
	_32_32_512  = 3,
        _64_64_1024 = 4
  };
};

template <int WType> struct WarpType;

// Turing Settings
template<> struct WarpType<TCWarpType::_64_64_512> {
  enum {w0 = 64, w1 = 64, w2 = 512};
};

template<> struct WarpType<TCWarpType::_64_32_512> {
  enum {w0 = 64, w1 = 32, w2 = 512};
};

template<> struct WarpType<TCWarpType::_32_64_512> {
  enum {w0 = 32, w1 = 64, w2 = 512};
};

template<> struct WarpType<TCWarpType::_32_32_512> {
  enum {w0 = 32, w1 = 32, w2 = 512};
};

// Ampere Settings
#if defined COMET_USE_AMPERE
template<> struct WarpType<TCWarpType::_64_64_1024> {
  enum {w0 = 64, w1 = 64, w2 = 1024};
};
#endif

struct TCInstType {
  enum{_8_8_128 = 0,
       _16_8_256 = 1
  };
};

template <int IType> struct InstType;

// Turing Settings
template<> struct InstType<TCInstType::_8_8_128> {
  enum {i0 = 8, i1 = 8, i2 = 128};
};

// Ampere Settings
#if defined COMET_USE_AMPERE
template<> struct InstType<TCInstType::_16_8_256> {
  enum {i0 = 16, i1 = 8, i2 = 256};
};
#endif

//-----------------------------------------------------------------------------
/// \brief 1-bit Int Tensor Core GEMM

template <int TBType, int WType, int IType, int NStages>
void CutlassTCGemm1B(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc, 
  AccelStream_t accel_stream) {

  /*printf("In CutlassTCGemm1B TB=%d,%d,%d W=%d,%d,%d I=%d,%d,%d\n",
    TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2,
    WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2,
    InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2);*/

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2>;
  using WarpShape        = cutlass::gemm::GemmShape<WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2>;
  using InstructionShape = cutlass::gemm::GemmShape<InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

#if defined COMET_USE_AMPERE
  //printf("Using Ampere Sm80\n");
  using ArchType = cutlass::arch::Sm80;
#elif defined COMET_USE_TURING
  //printf("Using Turing Sm75\n");
  using ArchType = cutlass::arch::Sm75;
#else
  //printf("Using default Sm50\n");
  using ArchType = cutlass::arch::Sm50;
#endif

  using Gemm = cutlass::gemm::device::Gemm<prec, RowMajor,    // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           ArchType,          // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           NStages, 128, 128, false, cutlass::arch::OpXorPopc>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int WMMA Tensor Core GEMM

template <int TBType, int WType, int IType, int NStages>
void CutlassTCGemm1BWmma(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc,
  AccelStream_t accel_stream) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2>;
  using WarpShape        = cutlass::gemm::GemmShape<WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2>;
  using InstructionShape = cutlass::gemm::GemmShape<InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

#if defined COMET_USE_AMPERE
  using ArchType = cutlass::arch::Sm80;
#elif defined COMET_USE_TURING
  using ArchType = cutlass::arch::Sm75;
#else
  using ArchType = cutlass::arch::Sm50;
#endif

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassWmmaTensorOp, // Operator class
                                           ArchType,          // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           NStages, 128, 128, false, cutlass::arch::OpXorPopc>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_CUTLASS_GENERAL_I_HH_

//-----------------------------------------------------------------------------
