/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Jakub Kurzak
       @author Stan Tomov
       @author Mark Gates
       @author Azzam Haidar

       See [zcds]gemm_fermi.cu for description of related files.
*/

///////////////////////////////////////////////////////////////////////////////////////////////////
// reset variables from previous includes of this file.
#undef TRANS_A
#undef TRANS_B
#undef CONJ_A
#undef CONJ_B

#undef BLK_M
#undef BLK_N

#undef THR_M
#undef THR_N

#undef batched_herk_kernel_name_
#undef batched_herk_kernel_name
#undef batched_gemm_kernel_name_
#undef batched_gemm_kernel_name
#undef gemm_kernel_name_
#undef gemm_kernel_name

#undef devfunc_name_
#undef devfunc_name

///////////////////////////////////////////////////////////////////////////////////////////////////

#if   (version == trans_nn)
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nn_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nn
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_nn
  #define BLK_M BLK_M_nn
  #define BLK_N BLK_N_nn

#elif (version == trans_nt)
  #define TRANS_B
  #define batched_herk_kernel_name_(p)  magma_tally3blas_ ## p ## _herk_kernel_fermi_nt_batched
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nt_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nt
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_nt
  #define BLK_M BLK_M_nt
  #define BLK_N BLK_N_nt

#elif (version == trans_nc)
  #define TRANS_B
  #define CONJ_B
  #define batched_herk_kernel_name_(p)  magma_tally3blas_ ## p ## _herk_kernel_fermi_nc_batched
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nc_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_nc
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_nc
  #define BLK_M BLK_M_nc
  #define BLK_N BLK_N_nc

#elif (version == trans_tn)
  #define TRANS_A
  #define batched_herk_kernel_name_(p)  magma_tally3blas_ ## p ## _herk_kernel_fermi_tn_batched
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tn_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tn
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_tn
  #define BLK_M BLK_M_tn
  #define BLK_N BLK_N_tn

#elif (version == trans_tt)
  #define TRANS_A
  #define TRANS_B
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tt_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tt
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_tt
  #define BLK_M BLK_M_tt
  #define BLK_N BLK_N_tt

#elif (version == trans_tc)
  #define TRANS_A
  #define TRANS_B
  #define CONJ_B
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tc_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_tc
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_tc
  #define BLK_M BLK_M_tc
  #define BLK_N BLK_N_tc

#elif (version == trans_cn)
  #define TRANS_A
  #define CONJ_A
  #define batched_herk_kernel_name_(p)  magma_tally3blas_ ## p ## _herk_kernel_fermi_cn_batched
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_cn_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_cn
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_cn
  #define BLK_M BLK_M_cn
  #define BLK_N BLK_N_cn

#elif (version == trans_ct)
  #define TRANS_A
  #define CONJ_A
  #define TRANS_B
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_ct_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_ct
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_ct
  #define BLK_M BLK_M_ct
  #define BLK_N BLK_N_ct

#elif (version == trans_cc)
  #define TRANS_A
  #define CONJ_A
  #define TRANS_B
  #define CONJ_B
  #define batched_gemm_kernel_name_(p)  p ## gemm_kernel_fermi_cc_batched
  #define gemm_kernel_name_(p)  p ## gemm_kernel_fermi_cc
  #define devfunc_name_(p) p ## gemm_devfunc_fermi_cc
  #define BLK_M BLK_M_cc
  #define BLK_N BLK_N_cc

#endif

// need a second macro in order to expand precision;
// see http://gcc.gnu.org/onlinedocs/cpp/Argument-Prescan.html
#define batched_herk_kernel_name(p) batched_herk_kernel_name_(p)
#define batched_gemm_kernel_name(p) batched_gemm_kernel_name_(p)
#define gemm_kernel_name(p) gemm_kernel_name_(p)
#define devfunc_name(p) devfunc_name_(p)

///////////////////////////////////////////////////////////////////////////////////////////////////

// size of work for a thread
#define THR_M ( BLK_M / DIM_X )
#define THR_N ( BLK_N / DIM_Y )

//---BEGIN ADDED GENOMICS_METRICS
//-----------------------------------------------------------------------------

//http://docs.nvidia.com/cuda/cuda-math-api/group__CUDA__MATH__INTRINSIC__INT.html
//__DEVICE_FUNCTIONS_DECL__ int __popc ( unsigned int  x )
//Count the number of bits that are set to 1 in a 32 bit integer.
//__DEVICE_FUNCTIONS_DECL__ int __popcll ( unsigned long long int x )
//Count the number of bits that are set to 1 in a 64 bit integer.

#if (version == trans_tn) && (precision == z)

#define gm_popcount64(x) __popcll(x)

extern "C" static __device__
FloatingPoint_t tally3_compute(
    FloatingPoint_t A,
    FloatingPoint_t B,
    FloatingPoint_t C)
{
    typedef unsigned long long int GMUInt64;
    typedef GMUInt64 GMBits1_2x64;
    typedef struct { GMBits1_2x64 data[2]; } GMBits2x64;
    enum { GM_TALLY1_MAX_VALUE_BITS = 25 };
    typedef struct { double data[2]; } GMTally2x2;

    GMTally2x2 sum = *(GMTally2x2*)&C;

    /*---Extract input values to process---*/

    const GMBits2x64 vi = *(GMBits2x64*)&A;
    const GMBits2x64 vj = *(GMBits2x64*)&B;

          const GMUInt64 vi0 = vi.data[0];
          const GMUInt64 vi1 = vi.data[1];
          const GMUInt64 vj0 = vj.data[0];
          const GMUInt64 vj1 = vj.data[1];

          /*---Get mask to ignore vi seminibbles with value of 1,0---*/
          /*---NOTE: check that this handles pad properly---*/

          const GMUInt64 oddbits = 0x5555555555555555;

          const GMUInt64 vi0_10mask = (vi0 | ~(vi0 >> 1)) & oddbits;
          const GMUInt64 vi1_10mask = (vi1 | ~(vi1 >> 1)) & oddbits;

          /*---Get even, odd bits for each semi-nibble, then mask---*/

          //--------------------
          // Nomenclature:
          //
          // ( ) v(i)(0)_(0)
          // (n) v(j)(1)_(1)
          //  ^    ^ ^   ^
          //  |    | |   |--- lower or upper bit of each seminibble
          //  |    | |--- lower or upper word
          //  |    |--- left or right vector
          //  |---test for value or for its negative/complement
          //--------------------

          const GMUInt64 vi0_0 =  vi0       & vi0_10mask;
          const GMUInt64 vi0_1 = (vi0 >> 1) & vi0_10mask;
          const GMUInt64 vi1_0 =  vi1       & vi1_10mask;
          const GMUInt64 vi1_1 = (vi1 >> 1) & vi1_10mask;
          const GMUInt64 vj0_0 =  vj0       & oddbits;
          const GMUInt64 vj0_1 = (vj0 >> 1) & oddbits;
          const GMUInt64 vj1_0 =  vj1       & oddbits;
          const GMUInt64 vj1_1 = (vj1 >> 1) & oddbits;

          /*---Get complements of the same bits, then mask---*/

          const GMUInt64 nvi0_0 = ~ vi0       & vi0_10mask;
          const GMUInt64 nvi0_1 = ~(vi0 >> 1) & vi0_10mask;
          const GMUInt64 nvi1_0 = ~ vi1       & vi1_10mask;
          const GMUInt64 nvi1_1 = ~(vi1 >> 1) & vi1_10mask;
          const GMUInt64 nvj0_0 = ~ vj0       & oddbits;
          const GMUInt64 nvj0_1 = ~(vj0 >> 1) & oddbits;
          const GMUInt64 nvj1_0 = ~ vj1       & oddbits;
          const GMUInt64 nvj1_1 = ~(vj1 >> 1) & oddbits;

          const int r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 & nvj0_0) |
                                      ( (nvi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 & nvj1_0) |
                                      ( (nvi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 & nvj1_0) |
                                      ( (nvi1_1 & nvj1_1) << 1 ));
          const int r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 &  vj0_0) |
                                      ( (nvi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 &  vj1_0) |
                                      ( (nvi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 &  vj1_0) |
                                      ( (nvi1_1 &  vj1_1) << 1 ));
          const int r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 & nvj0_0) |
                                      ( ( vi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 & nvj1_0) |
                                      ( ( vi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 & nvj1_0) |
                                      ( ( vi1_1 & nvj1_1) << 1 ));
          const int r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 &  vj0_0) |
                                      ( ( vi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 &  vj1_0) |
                                      ( ( vi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 &  vj1_0) |
                                      ( ( vi1_1 &  vj1_1) << 1 ));

          /*---Accumulate---*/

          sum.data[0] += r00 + (((GMUInt64)1)<<GM_TALLY1_MAX_VALUE_BITS) * r01;
          sum.data[1] += r10 + (((GMUInt64)1)<<GM_TALLY1_MAX_VALUE_BITS) * r11;

//printf("%i %i %i %i %e %e\n", r00, r01, r10, r11, sum.data[0], sum.data[1]);
    return *(FloatingPoint_t*)&sum;
}
#endif

//-----------------------------------------------------------------------------
//---END ADDED GENOMICS_METRICS
///////////////////////////////////////////////////////////////////////////////////////////////////
extern "C" static __device__
void devfunc_name(precision) (
    int M, int N, int K,
    const FloatingPoint_t* __restrict__ A, int LDA,
    const FloatingPoint_t* __restrict__ B, int LDB,
    FloatingPoint_t*       __restrict__ C, int LDC,
    FloatingPoint_t alpha, FloatingPoint_t beta,
    int offsetA, int offsetB )
{
#if (__CUDA_ARCH__ >= 200)
    int idx = threadIdx.x;  // thread's m dimension
    int idy = threadIdx.y;  // thread's n dimension

    int idt = DIM_X * idy + idx;    // thread's global number

    int idxA = idt % DIM_XA;    // idx within A
    int idyA = idt / DIM_XA;    // idy within A

    int idxB = idt % DIM_XB;    // idx within B
    int idyB = idt / DIM_XB;    // idy within B

    int blx = blockIdx.x;   // block's m dimension
    int bly = blockIdx.y;   // block's n dimension

    __shared__ FloatingPoint_t sA[BLK_K][BLK_M+1];      // +1 only required if A is transposed
    __shared__ FloatingPoint_t sB[BLK_N][BLK_K+1];      // +1 always required

    // Registers for the innermost loop
    FloatingPoint_t rC[THR_N][THR_M];
    FloatingPoint_t rA[THR_M];
    FloatingPoint_t rB[THR_N];

    // Registers for the dev->shmem copy
    #ifdef TRANS_A
        FloatingPoint_t ra[BLK_M/DIM_YA][BLK_K/DIM_XA];
    #else
        FloatingPoint_t ra[BLK_K/DIM_YA][BLK_M/DIM_XA];
    #endif
    #ifdef TRANS_B
        FloatingPoint_t rb[BLK_K/DIM_YB][BLK_N/DIM_XB];
    #else
        FloatingPoint_t rb[BLK_N/DIM_YB][BLK_K/DIM_XB];
    #endif

    #ifdef TEXTURE_1D
        #ifdef TRANS_A
            int coord_A = offsetA + blx*BLK_M*LDA + idyA*LDA + idxA;
        #else
            int coord_A = offsetA + blx*BLK_M     + idyA*LDA + idxA;
        #endif
        #ifdef TRANS_B
            int coord_B = offsetB + bly*BLK_N     + idyB*LDB + idxB;
        #else
            int coord_B = offsetB + bly*BLK_N*LDB + idyB*LDB + idxB;
        #endif
    #else
        // bound is the correction to offs_d in order to not get out of memory bound
        // so bound could be negative value since offs_d could be out of bound
        #ifdef TRANS_A
            const FloatingPoint_t *offs_dA = A + blx*BLK_M*LDA + idyA*LDA + idxA;
            ptrdiff_t boundA = (LDA*(M-1) + K) - ( blx*BLK_M*LDA + idyA*LDA + idxA ) -1;
        #else
            const FloatingPoint_t *offs_dA = A + blx*BLK_M     + idyA*LDA + idxA;
            ptrdiff_t boundA = (LDA*(K-1) + M) - ( blx*BLK_M  + idyA*LDA + idxA ) -1;
        #endif
        #ifdef TRANS_B
            const FloatingPoint_t *offs_dB = B + bly*BLK_N     + idyB*LDB + idxB;
            ptrdiff_t boundB = (LDB*(K-1) + N) - ( bly*BLK_N     + idyB*LDB + idxB ) -1;
        #else
            const FloatingPoint_t *offs_dB = B + bly*BLK_N*LDB + idyB*LDB + idxB;
            ptrdiff_t boundB = (LDB*(N-1) + K) - ( bly*BLK_N*LDB + idyB*LDB + idxB ) -1;
        #endif
    #endif

    int m, n, k, kk;


    // Zero C
    #pragma unroll
    for (n = 0; n < THR_N; n++)
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rC[n][m] = make_FloatingPoint(0.0, 0.0);

    // Load A dev->shmem
    #ifdef TRANS_A
        #pragma unroll
        for (n = 0; n < BLK_M; n += DIM_YA)
            #pragma unroll
            for (m = 0; m < BLK_K; m += DIM_XA)
                sA[m+idxA][n+idyA] = fetch(A, m, n, boundA);
    #else
        #pragma unroll
        for (n = 0; n < BLK_K; n += DIM_YA)
            #pragma unroll
            for (m = 0; m < BLK_M; m += DIM_XA)
                sA[n+idyA][m+idxA] = fetch(A, m, n, boundA);
    #endif

    // Load B dev->shmem
    #ifdef TRANS_B
        #pragma unroll
        for (n = 0; n < BLK_K; n += DIM_YB)
            #pragma unroll
            for (m = 0; m < BLK_N; m += DIM_XB)
                sB[m+idxB][n+idyB] = fetch(B, m, n, boundB);
    #else
        #pragma unroll
        for (n = 0; n < BLK_N; n += DIM_YB)
            #pragma unroll
            for (m = 0; m < BLK_K; m += DIM_XB)
                sB[n+idyB][m+idxB] = fetch(B, m, n, boundB);
    #endif

    __syncthreads();

    for (kk = 0; kk < K-BLK_K; kk += BLK_K)
    {
        #ifdef TEXTURE_1D
            #ifdef TRANS_A
                coord_A += BLK_K;
            #else
                coord_A += BLK_K*LDA;
            #endif
            #ifdef TRANS_B
                coord_B += BLK_K*LDB;
            #else
                coord_B += BLK_K;
            #endif
        #else
            #ifdef TRANS_A
                offs_dA += BLK_K;
                boundA  -= BLK_K;
            #else
                offs_dA += BLK_K*LDA;
                boundA  -= BLK_K*LDA;
            #endif
            #ifdef TRANS_B
                offs_dB += BLK_K*LDB;
                boundB  -= BLK_K*LDB;
            #else
                offs_dB += BLK_K;
                boundB  -= BLK_K;
            #endif
        #endif

        // Load A dev->regs
        #ifdef TRANS_A
            #pragma unroll
            for (n = 0; n < BLK_M/DIM_YA; n++)
                #pragma unroll
                for (m = 0; m < BLK_K/DIM_XA; m++)
                    ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);
        #else
            #pragma unroll
            for (n = 0; n < BLK_K/DIM_YA; n++)
                #pragma unroll
                for (m = 0; m < BLK_M/DIM_XA; m++)
                    ra[n][m] = fetch(A, m*DIM_XA, n*DIM_YA, boundA);
        #endif

        // Load B dev->regs
        #ifdef TRANS_B
            #pragma unroll
            for (n = 0; n < BLK_K/DIM_YB; n++)
                #pragma unroll
                for (m = 0; m < BLK_N/DIM_XB; m++)
                    rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);
        #else
            #pragma unroll
            for (n = 0; n < BLK_N/DIM_YB; n++)
                #pragma unroll
                for (m = 0; m < BLK_K/DIM_XB; m++)
                    rb[n][m] = fetch(B, m*DIM_XB, n*DIM_YB, boundB);
        #endif

        // Multiply
        #pragma unroll
        for (k = 0; k < BLK_K; k++)
        {
            // Load A shmem->regs
            #pragma unroll
            for (m = 0; m < THR_M; m++)
                rA[m] = sA[k][m*DIM_X+idx];

            // Load B shmem->regs
            #pragma unroll
            for (n = 0; n < THR_N; n++)
                rB[n] = sB[n*DIM_Y+idy][k];

            // Compute
            #pragma unroll
            for (n = 0; n < THR_N; n++) {
                #pragma unroll
                for (m = 0; m < THR_M; m++) {
                    #ifdef CONJ_A
                      #ifdef CONJ_B
                        fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
                      #else
                        fma(conj(rA[m]), rB[n], rC[n][m]);
                      #endif
                    #else
                      #ifdef CONJ_B
                        fma(rA[m], conj(rB[n]), rC[n][m]);
                      #else
                        //---BEGIN MODIFIED GENOMICS_METRICS
                        #if (version == trans_tn) && (precision == z)
                          rC[n][m] = tally3_compute(rA[m], rB[n], rC[n][m]);
                        #else
                          fma(rA[m], rB[n], rC[n][m]);
                        #endif
                        //---END MODIFIED GENOMICS_METRICS
                      #endif
                    #endif
                }
            }
        }

        __syncthreads();

        // Load A regs->shmem
        #ifdef TRANS_A
            #pragma unroll
            for (n = 0; n < BLK_M/DIM_YA; n++)
                #pragma unroll
                for (m = 0; m < BLK_K/DIM_XA; m++)
                    sA[m*DIM_XA+idxA][n*DIM_YA+idyA] = ra[n][m];
        #else
            #pragma unroll
            for (n = 0; n < BLK_K/DIM_YA; n++)
                #pragma unroll
                for (m = 0; m < BLK_M/DIM_XA; m++)
                    sA[n*DIM_YA+idyA][m*DIM_XA+idxA] = ra[n][m];
            #endif

        // Load B regs->shmem
        #ifdef TRANS_B
            #pragma unroll
            for (n = 0; n < BLK_K/DIM_YB; n++)
                #pragma unroll
                for (m = 0; m < BLK_N/DIM_XB; m++)
                    sB[m*DIM_XB+idxB][n*DIM_YB+idyB] = rb[n][m];
        #else
            #pragma unroll
            for (n = 0; n < BLK_N/DIM_YB; n++)
                #pragma unroll
                for (m = 0; m < BLK_K/DIM_XB; m++)
                    sB[n*DIM_YB+idyB][m*DIM_XB+idxB] = rb[n][m];
        #endif

        __syncthreads();
    }

    // Multiply last full (BLK_K) or partial block of
    // columns of op(A) and rows of op(B).
    // It's okay that m,n exceed matrix bounds as all work is in registers
    // or shared memory, and out-of-bounds rC[n][m] will not be saved later.
    kk = K - kk;
    #pragma unroll
    for (k = 0; k < kk; k++)
    {
        // Load A shmem->regs
        #pragma unroll
        for (m = 0; m < THR_M; m++)
            rA[m] = sA[k][m*DIM_X+idx];

        // Load B shmem->regs
        #pragma unroll
        for (n = 0; n < THR_N; n++)
            rB[n] = sB[n*DIM_Y+idy][k];

        // Compute
        #pragma unroll
        for (n = 0; n < THR_N; n++) {
            #pragma unroll
            for (m = 0; m < THR_M; m++) {
                #ifdef CONJ_A
                  #ifdef CONJ_B
                    fma(conj(rA[m]), conj(rB[n]), rC[n][m]);
                  #else
                    fma(conj(rA[m]), rB[n], rC[n][m]);
                  #endif
                #else
                  #ifdef CONJ_B
                    fma(rA[m], conj(rB[n]), rC[n][m]);
                  #else
                    //---BEGIN MODIFIED GENOMICS_METRICS
                    #if (version == trans_tn) && (precision == z)
                      rC[n][m] = tally3_compute(rA[m], rB[n], rC[n][m]);
                    #else
                      fma(rA[m], rB[n], rC[n][m]);
                    #endif
                    //---END MODIFIED GENOMICS_METRICS
                  #endif
                #endif
            }
        }
    }

    // Store C regs->dev
    #pragma unroll
    for (n = 0; n < THR_N; n++) {
        int coord_dCn = bly*BLK_N + n*DIM_Y + idy;
        #pragma unroll
        for (m = 0; m < THR_M; m++) {
            int coord_dCm = blx*BLK_M + m*DIM_X + idx;
            if (coord_dCm < M && coord_dCn < N) {
                int offsC = coord_dCn*LDC + coord_dCm;

                FloatingPoint_t &regC = rC[n][m];
                FloatingPoint_t &memC = C[offsC];

                memC = add(mul(alpha, regC), mul(beta, memC));
            }
        }
    }
#endif /* (__CUDA_ARCH__ >= 200) */
}
