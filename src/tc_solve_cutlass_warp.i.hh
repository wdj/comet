#ifndef _COMET_TC_SOLVE_CUTLASS_WARP_I_HH_
#define _COMET_TC_SOLVE_CUTLASS_WARP_I_HH_

//-----------------------------------------------------------------------------
/// \brief 1-bit Int WMMA GEMM with 64 x 64 blocks

cudaError_t CutlassTCWarp_64x64(int M, int N, int K, uint8_t const *A,
  int lda, uint8_t const *B, int ldb, int32_t *C, int ldc) {

}
#endif

