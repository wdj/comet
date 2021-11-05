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
/// \brief 1-bit Int GEMM with variable Thread Block, Warp, and Instruction
//         sizes

struct TCTBlockType {
  enum {_64_64_512    = 0,
	_64_128_512   = 1,
	_64_256_512   = 2,
        _128_64_512   = 3,
        _128_128_512  = 4,
	_128_256_512  = 5,
	_256_64_512   = 6,
	_256_128_512  = 7,
	_64_64_1024   = 8,
	_64_128_1024  = 9,
	_64_256_1024  = 10,
	_128_64_1024  = 11,
	_128_128_1024 = 12,
        _128_256_1024 = 13,
	_256_64_1024  = 14,
	_256_128_1024 = 15,
	_160_256_1024 = 16,
	_160_288_1024 = 17,
	_192_224_1024 = 18,
	_192_192_1024 = 19
	// Unlisted - Very slow
	/*_128_288_1024 = 17,
	// Unlisted - Runs out of memory
        _160_288_1024 = 17
	_192_256_1024 = 18,
	_256_256_1024 = 19*/
  };
};

template <int TBType> struct TBlockType;

template<> struct TBlockType<TCTBlockType::_192_192_1024> {
  enum {t0 = 192, t1 = 192, t2 = 1024};
};

// Turing Settings
template<> struct TBlockType<TCTBlockType::_64_64_512> {
  enum {t0 = 64, t1 = 64, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_64_128_512> {
  enum {t0 = 64, t1 = 128, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_128_64_512> {
  enum {t0 = 128, t1 = 64, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_128_128_512> {
  enum {t0 = 128, t1 = 128, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_128_256_512> {
  enum {t0 = 128, t1 = 256, t2 = 512};
};

template<> struct TBlockType<TCTBlockType::_256_128_512> {
  enum {t0 = 256, t1 = 128, t2 = 512};
};

// Ampere Settings
#if defined COMET_USE_AMPERE
template<> struct TBlockType<TCTBlockType::_64_64_1024> {
  enum {t0 = 64, t1 = 64, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_64_256_1024> {
  enum {t0 = 64, t1 = 256, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_128_128_1024> {
  enum {t0 = 128, t1 = 128, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_128_256_1024> {
  enum {t0 = 128, t1 = 256, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_256_64_1024> {
  enum {t0 = 256, t1 = 64, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_256_128_1024> {
  enum {t0 = 256, t1 = 128, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_160_256_1024> {
  enum {t0 = 160, t1 = 256, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_192_224_1024> {
  enum {t0 = 192, t1 = 224, t2 = 1024};
};

// Very slow
/*template<> struct TBlockType<TCTBlockType::_128_288_1024> {
  enum {t0 = 128, t1 = 288, t2 = 1024};
};

// Runs out of memory
template<> struct TBlockType<TCTBlockType::_160_288_1024> {
  enum {t0 = 160, t1 = 288, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_192_256_1024> {
  enum {t0 = 192, t1 = 256, t2 = 1024};
};

template<> struct TBlockType<TCTBlockType::_256_256_1024> {
  enum {t0 = 256, t1 = 256, t2 = 1024};
};*/
#endif

struct TCWarpType {
  enum {_32_32_512  = 0,
	_64_32_512  = 1,
	_64_64_512  = 2,
	_32_64_512  = 3,
	_32_32_1024 = 4,
	_32_64_1024 = 5,
	_64_32_1024 = 6,
        _64_64_1024 = 7,
	_64_128_512 = 8,
        _80_64_1024 = 9,
	_80_72_1024 = 10,
	_96_56_1024 = 11,
	_32_128_1024 = 12
	// Unlisted - Very slow/ran out of memory
	/*_80_72_1024 = 9
	_64_96_1024 = 9,
	_96_64_1024 = 10,
	_64_128_1024 = 11,
	_128_64_1024 = 12,
	_128_128_1024 = 13*/
  };
};

template <int WType> struct WarpType;

template<> struct WarpType<TCWarpType::_32_128_1024> {
  enum {w0 = 32, w1 = 128, w2 = 1024};
};

// Turing Settings
template<> struct WarpType<TCWarpType::_32_32_512> {
  enum {w0 = 32, w1 = 32, w2 = 512};
};

template<> struct WarpType<TCWarpType::_32_64_512> {
  enum {w0 = 32, w1 = 64, w2 = 512};
};

template<> struct WarpType<TCWarpType::_64_32_512> {
  enum {w0 = 64, w1 = 32, w2 = 512};
};

template<> struct WarpType<TCWarpType::_64_64_512> {
  enum {w0 = 64, w1 = 64, w2 = 512};
};

// Ampere Settings
#if defined COMET_USE_AMPERE
template<> struct WarpType<TCWarpType::_32_32_1024> {
  enum {w0 = 32, w1 = 32, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_32_64_1024> {
  enum {w0 = 32, w1 = 64, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_64_32_1024> {
  enum {w0 = 64, w1 = 32, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_64_64_1024> {
  enum {w0 = 64, w1 = 64, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_64_128_512> {
  enum {w0 = 64, w1 = 128, w2 = 512};
};

template<> struct WarpType<TCWarpType::_80_64_1024> {
  enum {w0 = 80, w1 = 64, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_80_72_1024> {
  enum {w0 = 80, w1 = 72, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_96_56_1024> {
  enum {w0 = 96, w1 = 56, w2 = 1024};
};

// Very slow settings
/*template<> struct WarpType<TCWarpType::_64_96_1024> {
  enum {w0 = 64, w1 = 96, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_96_64_1024> {
  enum {w0 = 96, w1 = 64, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_64_128_1024> {
  enum {w0 = 64, w1 = 128, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_128_64_1024> {
  enum {w0 = 128, w1 = 64, w2 = 1024};
};

template<> struct WarpType<TCWarpType::_128_128_1024> {
  enum {w0 = 128, w1 = 128, w2 = 1024};
};*/
#endif

// Options are defined in include/cutlass/arch/mma_sm80.h
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

struct TCOpType {
  enum{Xor  = 0,
       Mult = 1
  };
};

template <int OType> struct OpType;

template<> struct OpType<TCOpType::Xor> {
  typedef cutlass::arch::OpXorPopc Op;
};

template<> struct OpType<TCOpType::Mult> {
  typedef cutlass::arch::OpMultiplyAdd Op;
};

//-----------------------------------------------------------------------------
/// \brief 1-bit Int Tensor Core GEMM

template <int TBType, int WType, int IType, int OType, int NStages>
void CutlassTCGemm1B(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc, 
  AccelStream_t accel_stream) {

  printf("In CutlassTCGemm1B TB=%d,%d,%d W=%d,%d,%d I=%d,%d,%d\n",
    TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2,
    WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2,
    InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2);

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2>;
  using WarpShape        = cutlass::gemm::GemmShape<WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2>;
  using InstructionShape = cutlass::gemm::GemmShape<InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;
  using OperationType    = typename OpType<OType>::Op;

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
                                           NStages,
					   128, // AlignmentA
					   128, // AlignmentB
					   false, // SplitKSerial
					   OperationType>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int Tensor Core GEMM

template <int TBType, int WType, int IType, int OType, int NStages>
void CutlassTCGemm1BTest(int M, int N, int K, uint8_t const *A,
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
  using OperationType    = typename OpType<OType>::Op;

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
                                           //cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
					   cutlass::gemm::threadblock::GemmHorizontalThreadblockSwizzle,
                                           NStages,
                                           128, // AlignmentA
                                           128, // AlignmentB
                                           false, // SplitKSerial
                                           OperationType>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int WMMA Tensor Core GEMM

template <int TBType, int WType, int IType, int OType, int NStages>
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
  using OperationType    = typename OpType<OType>::Op;

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
                                           NStages, 128, 128, false, OperationType>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int Tensor Core GEMM

template <int TBType, int WType, int IType, int OType, int NStages>
void CutlassTCGemm1BSplitK(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc,
  AccelStream_t accel_stream) {

  //printf("In CutlassTCGemm1B TB=%d,%d,%d W=%d,%d,%d I=%d,%d,%d\n",
  //  TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2,
  //  WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2,
  //  InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2);

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<TBlockType<TBType>::t0,TBlockType<TBType>::t1,TBlockType<TBType>::t2>;
  using WarpShape        = cutlass::gemm::GemmShape<WarpType<WType>::w0,WarpType<WType>::w1,WarpType<WType>::w2>;
  using InstructionShape = cutlass::gemm::GemmShape<InstType<IType>::i0,InstType<IType>::i1,InstType<IType>::i2>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;
  using OperationType    = typename OpType<OType>::Op;

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
                                           ThreadBlockShape,
					   WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           //cutlass::gemm::threadblock::GemmHorizontalThreadblockSwizzle<>,
                                           NStages,
                                           128, // AlignmentA
                                           128, // AlignmentB
                                           true, // SplitKSerial
                                           OperationType>;

  CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                       ldb,beta,(accprec*)C,ldc,accel_stream);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_CUTLASS_GENERAL_I_HH_

//-----------------------------------------------------------------------------
