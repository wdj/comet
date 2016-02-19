/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @generated from zhegst_gpu.cpp normal z -> s, Fri Jan 30 19:00:18 2015
*/

#include "common_magma_tally3.h"

/**
    Purpose
    -------
    SSYGST_GPU reduces a real symmetric-definite generalized
    eigenproblem to standard form.
    
    If ITYPE = 1, the problem is A*x = lambda*B*x,
    and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
    
    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
    B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
    
    B must have been previously factorized as U**H*U or L*L**H by SPOTRF.
    
    Arguments
    ---------
    @param[in]
    itype   INTEGER
            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
            = 2 or 3: compute U*A*U**H or L**H*A*L.
    
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of A is stored and B is factored as
                    U**H*U;
      -     = Magma_tally3Lower:  Lower triangle of A is stored and B is factored as
                    L*L**H.
    
    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.
    
    @param[in,out]
    dA      REAL array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally3Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally3Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the transformed matrix, stored in the
            same format as A.
    
    @param[in]
    ldda    INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    @param[in]
    dB      REAL array, dimension (LDB,N)
            The triangular factor from the Cholesky factorization of B,
            as returned by SPOTRF.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally3_ssyev_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_ssygst_gpu(
    magma_tally3_int_t itype, magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3Float_ptr dA, magma_tally3_int_t ldda,
    magma_tally3Float_ptr dB, magma_tally3_int_t lddb,
    magma_tally3_int_t *info)
{
#define A(i, j) (w + (j)*lda + (i))
#define B(i, j) (w + nb*lda + (j)*ldb + (i))

#define dA(i, j) (dA + (j)*ldda + (i))
#define dB(i, j) (dB + (j)*lddb + (i))

    const char* uplo_ = lapack_uplo_const_tally3( uplo );
    magma_tally3_int_t        nb;
    magma_tally3_int_t        k, kb, kb2;
    float    c_one      = MAGMA_tally3_S_ONE;
    float    c_neg_one  = MAGMA_tally3_S_NEG_ONE;
    float    c_half     = MAGMA_tally3_S_HALF;
    float    c_neg_half = MAGMA_tally3_S_NEG_HALF;
    float   *w;
    magma_tally3_int_t        lda;
    magma_tally3_int_t        ldb;
    float             d_one = 1.0;
    int upper = (uplo == Magma_tally3Upper);
    
    /* Test the input parameters. */
    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally3Lower) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }
    
    /* Quick return */
    if ( n == 0 )
        return *info;
    
    nb = magma_tally3_get_ssygst_nb(n);
    
    lda = nb;
    ldb = nb;
    
    if (MAGMA_tally3_SUCCESS != magma_tally3_smalloc_pinned( &w, 2*nb*nb )) {
        *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    magma_tally3_queue_t stream[3];
    magma_tally3_queue_create( &stream[0] );
    magma_tally3_queue_create( &stream[1] );
    magma_tally3_queue_create( &stream[2] );
    
    /* Use hybrid blocked code */
    if (itype == 1) {
        if (upper) {
            kb = min(n,nb);
            
            /* Compute inv(U')*A*inv(U) */
            magma_tally3_sgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_tally3_sgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for (k = 0; k < n; k += nb) {
                kb = min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the upper triangle of A(k:n,k:n) */
                
                magma_tally3_queue_sync( stream[2] );
                magma_tally3_queue_sync( stream[1] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0,0), &lda, B(0,0), &ldb, info);
                
                magma_tally3_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    // Start copying the new B block
                    magma_tally3_sgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_tally3_strsm(Magma_tally3Left, Magma_tally3Upper, Magma_tally3ConjTrans, Magma_tally3NonUnit,
                                kb, n-k-kb,
                                c_one, dB(k,k), lddb,
                                dA(k,k+kb), ldda);
                    
                    magma_tally3_queue_sync( stream[0] );
                    
                    magma_tally3_ssymm(Magma_tally3Left, Magma_tally3Upper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_tally3_ssyr2k(Magma_tally3Upper, Magma_tally3ConjTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k,k+kb), ldda,
                                 dB(k,k+kb), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_tally3_sgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_tally3_ssymm(Magma_tally3Left, Magma_tally3Upper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_tally3_strsm(Magma_tally3Right, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit,
                                kb, n-k-kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k,k+kb), ldda);
                }
            }
            
            magma_tally3_queue_sync( stream[0] );
        }
        else {
            kb = min(n,nb);
            
            /* Compute inv(L)*A*inv(L') */
            magma_tally3_sgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_tally3_sgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the lower triangle of A(k:n,k:n) */
                
                magma_tally3_queue_sync( stream[2] );
                magma_tally3_queue_sync( stream[1] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_tally3_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    // Start copying the new B block
                    magma_tally3_sgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_tally3_strsm(Magma_tally3Right, Magma_tally3Lower, Magma_tally3ConjTrans, Magma_tally3NonUnit,
                                n-k-kb, kb,
                                c_one, dB(k,k), lddb,
                                dA(k+kb,k), ldda);
                    
                    magma_tally3_queue_sync( stream[0] );
                    
                    magma_tally3_ssymm(Magma_tally3Right, Magma_tally3Lower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_tally3_ssyr2k(Magma_tally3Lower, Magma_tally3NoTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k+kb,k), ldda,
                                 dB(k+kb,k), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_tally3_sgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_tally3_ssymm(Magma_tally3Right, Magma_tally3Lower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_tally3_strsm(Magma_tally3Left, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3NonUnit,
                                n-k-kb, kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k+kb,k), ldda);
                }
            }
        }
        
        magma_tally3_queue_sync( stream[0] );
    }
    else {
        if (upper) {
            /* Compute U*A*U' */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_tally3_sgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    magma_tally3_strmm(Magma_tally3Left, Magma_tally3Upper, Magma_tally3NoTrans, Magma_tally3NonUnit,
                                k, kb,
                                c_one, dB(0,0), lddb,
                                dA(0,k), ldda);
                    
                    magma_tally3_ssymm(Magma_tally3Right, Magma_tally3Upper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_tally3_queue_sync( stream[1] );
                }
                
                magma_tally3_sgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if (k > 0) {
                    magma_tally3_ssyr2k(Magma_tally3Upper, Magma_tally3NoTrans,
                                 k, kb,
                                 c_one, dA(0,k), ldda,
                                 dB(0,k), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_tally3_ssymm(Magma_tally3Right, Magma_tally3Upper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_tally3_strmm(Magma_tally3Right, Magma_tally3Upper, Magma_tally3ConjTrans, Magma_tally3NonUnit,
                                k, kb,
                                c_one, dB(k,k), lddb,
                                dA(0,k), ldda);
                }
                
                magma_tally3_queue_sync( stream[2] );
                magma_tally3_queue_sync( stream[0] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_tally3_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_tally3_queue_sync( stream[1] );
        }
        else {
            /* Compute L'*A*L */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_tally3_sgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    magma_tally3_strmm(Magma_tally3Right, Magma_tally3Lower, Magma_tally3NoTrans, Magma_tally3NonUnit,
                                kb, k,
                                c_one, dB(0,0), lddb,
                                dA(k,0), ldda);
                    
                    magma_tally3_ssymm(Magma_tally3Left, Magma_tally3Lower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_tally3_queue_sync( stream[1] );
                }
                
                magma_tally3_sgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if (k > 0) {
                    magma_tally3_ssyr2k(Magma_tally3Lower, Magma_tally3ConjTrans,
                                 k, kb,
                                 c_one, dA(k,0), ldda,
                                 dB(k,0), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_tally3_ssymm(Magma_tally3Left, Magma_tally3Lower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_tally3_strmm(Magma_tally3Left, Magma_tally3Lower, Magma_tally3ConjTrans, Magma_tally3NonUnit,
                                kb, k,
                                c_one, dB(k,k), lddb,
                                dA(k,0), ldda);
                }
                
                magma_tally3_queue_sync( stream[2] );
                magma_tally3_queue_sync( stream[0] );
                
                lapackf77_ssygst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_tally3_ssetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_tally3_queue_sync( stream[1] );
        }
    }
    magma_tally3_queue_destroy( stream[0] );
    magma_tally3_queue_destroy( stream[1] );
    magma_tally3_queue_destroy( stream[2] );
    
    magma_tally3_free_pinned( w );
    
    return *info;
} /* magma_tally3_ssygst_gpu */

#undef A
#undef B
#undef dA
#undef dB
