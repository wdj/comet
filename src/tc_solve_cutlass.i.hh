#ifndef _COMET_TC_SOLVE_CUTLASS_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_I_HH_

// Defines cutlass::gemm::device::Gemm, the generic Gemm computation template class.
#include "cutlass/gemm/device/gemm.h"

/********************************************************
 * Cutlass GEMM driver routine
 *******************************************************/
template <typename Gemm>
cudaError_t CutlassGemmRun(int M, int N, int K, cutlass::uint1b_t const *A,
  int lda, cutlass::uint1b_t const *B, int ldb, int32_t *C, int ldc) {

  // Define a CUTLASS GEMM type
  Gemm gemm_operator;
  int32_t alpha = 1, beta = 0;

  // Construct the CUTLASS GEMM arguments object.
  //
  // One of CUTLASS's design patterns is to define gemm argument objects that are constructible
  // in host code and passed to kernels by value. These may include pointers, strides, scalars,
  // and other arguments needed by Gemm and its components.
  //
  // The benefits of this pattern are (1.) a structured, composable strategy for passing host-constructible
  // arguments to kernels and (2.) minimized initialization overhead on kernel entry.
  typename Gemm::Arguments args({M, N, K},      // Gemm Problem dimensions
                                {A, lda},       // Tensor-ref for source matrix A
                                {B, ldb},       // Tensor-ref for source matrix B
                                {C, ldc},       // Tensor-ref for source matrix C
                                {C, ldc},       // Tensor-ref for destination matrix D (may be different memory than C)
                                {alpha,beta}); // Scalars used in the Epilogue

  // Launch the CUTLASS GEMM kernel
  cutlass::Status status = gemm_operator(args);

  // Return a cudaError_t if the CUTLASS GEMM operator returned an error code.
  if (status != cutlass::Status::kSuccess) return cudaErrorUnknown;

  // Return success, if no errors were encountered.
  return cudaSuccess;
}

/***********************************
 * 1-bit Int GEMM
 **********************************/
cudaError_t CutlassTCGemm1B_256x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<256,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1B_128x256(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,256,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1B_128x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1B_128x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<128,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<64,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1B_64x128(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,128,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,64,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1B_64x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassTensorOp, // Operator class
                                           cutlass::arch::Sm75,            // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

cudaError_t CutlassTCGemm1BWmma_64x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

  using precision        = cutlass::uint1b_t;
  using accprecision     = int32_t;
  using RowMajor         = cutlass::layout::RowMajor;
  using ColumnMajor      = cutlass::layout::ColumnMajor;
  using ThreadBlockShape = cutlass::gemm::GemmShape<64,64,512>;
  using WarpShape        = cutlass::gemm::GemmShape<32,32,512>;
  using InstructionShape = cutlass::gemm::GemmShape<8,8,128>;
  using EpilogueOutputOp = cutlass::epilogue::thread::LinearCombination<
    accprecision, 128/cutlass::sizeof_bits<accprecision>::value, accprecision, accprecision>;

  using Gemm = cutlass::gemm::device::Gemm<precision,RowMajor,     // Data-type/Layout of A
                                           precision, ColumnMajor, // Data-type/Layout of B
                                           accprecision, RowMajor, // Data-type/Layout of C
                                           accprecision,           // Data-type of accumulator
                                           cutlass::arch::OpClassWmmaTensorOp, // Operator class
                                           cutlass::arch::Sm75,                // Architecture type
                                           ThreadBlockShape, WarpShape,
                                           InstructionShape, EpilogueOutputOp,
                                           cutlass::gemm::threadblock::GemmIdentityThreadblockSwizzle<>,
                                           2, 128, 128, false, cutlass::arch::OpXorPopc>;

  cudaError_t result = CutlassGemmRun<Gemm>(M,N,K,(precision const*)A,lda,(precision const*)B,ldb,(accprecision*)C,ldc);
  return result;
}

#endif

