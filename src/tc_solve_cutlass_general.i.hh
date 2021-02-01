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
cudaError_t CutlassGemmRun(int M, int N, int K, cutlass::uint1b_t const *A,
  int lda, cutlass::uint1b_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  // Define a CUTLASS GEMM type
  Gemm gemm_operator;
  int32_t alpha = 1;

  // Construct the CUTLASS GEMM arguments object.
  typename Gemm::Arguments args({M, N, K},      // Gemm Problem dimensions
                                {A, lda},       // Tensor-ref for source matrix A
                                {B, ldb},       // Tensor-ref for source matrix B
                                {C, ldc},       // Tensor-ref for source matrix C
                                {C, ldc},       // Tensor-ref for destination matrix C
                                {alpha,beta}); // Scalars used in the Epilogue

  // Launch the CUTLASS GEMM kernel
  cutlass::Status status = gemm_operator(args);

  // Return a cudaError_t if the CUTLASS GEMM operator returned an error code.
  if (status != cutlass::Status::kSuccess) return cudaErrorUnknown;

  // Return success, if no errors were encountered.
  return cudaSuccess;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 256 x 128 blocks

cudaError_t CutlassTCGemm1B_256x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<256,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 128 x 256 blocks

cudaError_t CutlassTCGemm1B_128x256(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,256,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 128 x 128 blocks

cudaError_t CutlassTCGemm1B_128x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 128 x 64 blocks

cudaError_t CutlassTCGemm1B_128x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 64 x 128 blocks

cudaError_t CutlassTCGemm1B_64x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int GEMM with 64 x 64 blocks

cudaError_t CutlassTCGemm1B_64x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//-----------------------------------------------------------------------------
/// \brief 1-bit Int WMMA GEMM with 64 x 64 blocks

cudaError_t CutlassTCGemm1BWmma_64x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t beta, int32_t *C, int ldc) {

  using prec             = cutlass::uint1b_t;
  using accprec          = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprec, 128/cutlass::sizeof_bits<accprec>::value, accprec, accprec>;

  using Gemm = cutlass::gemm::device::Gemm<prec,RowMajor,     // Data-type/Layout of A
                                           prec, ColumnMajor, // Data-type/Layout of B
                                           accprec, RowMajor, // Data-type/Layout of C
                                           accprec,           // Data-type of accumulator
                                           cutlass::arch::OpClassWmmaTensorOp, // Operator class
                                           cutlass::arch::Sm75,                // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(prec const*)A,lda,(prec const*)B,
                                            ldb,beta,(accprec*)C,ldc);
  return result;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_TC_SOLVE_CUTLASS_GENERAL_I_HH_

//-----------------------------------------------------------------------------
