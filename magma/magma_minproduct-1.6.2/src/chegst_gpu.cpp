/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @generated from zhegst_gpu.cpp normal z -> c, Fri Jan 30 19:00:18 2015
*/

#include "common_magma_minproduct.h"

/**
    Purpose
    -------
    CHEGST_GPU reduces a complex Hermitian-definite generalized
    eigenproblem to standard form.
    
    If ITYPE = 1, the problem is A*x = lambda*B*x,
    and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
    
    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
    B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
    
    B must have been previously factorized as U**H*U or L*L**H by CPOTRF.
    
    Arguments
    ---------
    @param[in]
    itype   INTEGER
            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
            = 2 or 3: compute U*A*U**H or L**H*A*L.
    
    @param[in]
    uplo    magma_minproduct_uplo_t
      -     = Magma_minproductUpper:  Upper triangle of A is stored and B is factored as
                    U**H*U;
      -     = Magma_minproductLower:  Lower triangle of A is stored and B is factored as
                    L*L**H.
    
    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.
    
    @param[in,out]
    dA      COMPLEX array, dimension (LDA,N)
            On entry, the Hermitian matrix A.  If UPLO = Magma_minproductUpper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_minproductLower, the
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
    dB      COMPLEX array, dimension (LDB,N)
            The triangular factor from the Cholesky factorization of B,
            as returned by CPOTRF.
    
    @param[in]
    lddb    INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_minproduct_cheev_comp
    ********************************************************************/
extern "C" magma_minproduct_int_t
magma_minproduct_chegst_gpu(
    magma_minproduct_int_t itype, magma_minproduct_uplo_t uplo, magma_minproduct_int_t n,
    magma_minproductFloatComplex_ptr dA, magma_minproduct_int_t ldda,
    magma_minproductFloatComplex_ptr dB, magma_minproduct_int_t lddb,
    magma_minproduct_int_t *info)
{
#define A(i, j) (w + (j)*lda + (i))
#define B(i, j) (w + nb*lda + (j)*ldb + (i))

#define dA(i, j) (dA + (j)*ldda + (i))
#define dB(i, j) (dB + (j)*lddb + (i))

    const char* uplo_ = lapack_uplo_const( uplo );
    magma_minproduct_int_t        nb;
    magma_minproduct_int_t        k, kb, kb2;
    magma_minproductFloatComplex    c_one      = MAGMA_minproduct_C_ONE;
    magma_minproductFloatComplex    c_neg_one  = MAGMA_minproduct_C_NEG_ONE;
    magma_minproductFloatComplex    c_half     = MAGMA_minproduct_C_HALF;
    magma_minproductFloatComplex    c_neg_half = MAGMA_minproduct_C_NEG_HALF;
    magma_minproductFloatComplex   *w;
    magma_minproduct_int_t        lda;
    magma_minproduct_int_t        ldb;
    float             d_one = 1.0;
    int upper = (uplo == Magma_minproductUpper);
    
    /* Test the input parameters. */
    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! upper && uplo != Magma_minproductLower) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (ldda < max(1,n)) {
        *info = -5;
    } else if (lddb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_minproduct_xerbla( __func__, -(*info) );
        return *info;
    }
    
    /* Quick return */
    if ( n == 0 )
        return *info;
    
    nb = magma_minproduct_get_chegst_nb(n);
    
    lda = nb;
    ldb = nb;
    
    if (MAGMA_minproduct_SUCCESS != magma_minproduct_cmalloc_pinned( &w, 2*nb*nb )) {
        *info = MAGMA_minproduct_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    magma_minproduct_queue_t stream[3];
    magma_minproduct_queue_create( &stream[0] );
    magma_minproduct_queue_create( &stream[1] );
    magma_minproduct_queue_create( &stream[2] );
    
    /* Use hybrid blocked code */
    if (itype == 1) {
        if (upper) {
            kb = min(n,nb);
            
            /* Compute inv(U')*A*inv(U) */
            magma_minproduct_cgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_minproduct_cgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for (k = 0; k < n; k += nb) {
                kb = min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the upper triangle of A(k:n,k:n) */
                
                magma_minproduct_queue_sync( stream[2] );
                magma_minproduct_queue_sync( stream[1] );
                
                lapackf77_chegst( &itype, uplo_, &kb, A(0,0), &lda, B(0,0), &ldb, info);
                
                magma_minproduct_csetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    // Start copying the new B block
                    magma_minproduct_cgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                kb, n-k-kb,
                                c_one, dB(k,k), lddb,
                                dA(k,k+kb), ldda);
                    
                    magma_minproduct_queue_sync( stream[0] );
                    
                    magma_minproduct_chemm(Magma_minproductLeft, Magma_minproductUpper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_minproduct_cher2k(Magma_minproductUpper, Magma_minproductConjTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k,k+kb), ldda,
                                 dB(k,k+kb), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_minproduct_cgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_minproduct_chemm(Magma_minproductLeft, Magma_minproductUpper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_minproduct_ctrsm(Magma_minproductRight, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit,
                                kb, n-k-kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k,k+kb), ldda);
                }
            }
            
            magma_minproduct_queue_sync( stream[0] );
        }
        else {
            kb = min(n,nb);
            
            /* Compute inv(L)*A*inv(L') */
            magma_minproduct_cgetmatrix_async( kb, kb,
                                    dB(0, 0), lddb,
                                    B(0, 0),  nb, stream[2] );
            magma_minproduct_cgetmatrix_async( kb, kb,
                                    dA(0, 0), ldda,
                                    A(0, 0),  nb, stream[1] );
            
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the lower triangle of A(k:n,k:n) */
                
                magma_minproduct_queue_sync( stream[2] );
                magma_minproduct_queue_sync( stream[1] );
                
                lapackf77_chegst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_minproduct_csetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    // Start copying the new B block
                    magma_minproduct_cgetmatrix_async( kb2, kb2,
                                            dB(k+kb, k+kb), lddb,
                                            B(0, 0),        nb, stream[2] );
                    
                    magma_minproduct_ctrsm(Magma_minproductRight, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                n-k-kb, kb,
                                c_one, dB(k,k), lddb,
                                dA(k+kb,k), ldda);
                    
                    magma_minproduct_queue_sync( stream[0] );
                    
                    magma_minproduct_chemm(Magma_minproductRight, Magma_minproductLower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_minproduct_cher2k(Magma_minproductLower, Magma_minproductNoTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k+kb,k), ldda,
                                 dB(k+kb,k), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_minproduct_cgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(0, 0),        lda, stream[1] );
                    
                    magma_minproduct_chemm(Magma_minproductRight, Magma_minproductLower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_minproduct_ctrsm(Magma_minproductLeft, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductNonUnit,
                                n-k-kb, kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k+kb,k), ldda);
                }
            }
        }
        
        magma_minproduct_queue_sync( stream[0] );
    }
    else {
        if (upper) {
            /* Compute U*A*U' */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_minproduct_cgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    magma_minproduct_ctrmm(Magma_minproductLeft, Magma_minproductUpper, Magma_minproductNoTrans, Magma_minproductNonUnit,
                                k, kb,
                                c_one, dB(0,0), lddb,
                                dA(0,k), ldda);
                    
                    magma_minproduct_chemm(Magma_minproductRight, Magma_minproductUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_minproduct_queue_sync( stream[1] );
                }
                
                magma_minproduct_cgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if (k > 0) {
                    magma_minproduct_cher2k(Magma_minproductUpper, Magma_minproductNoTrans,
                                 k, kb,
                                 c_one, dA(0,k), ldda,
                                 dB(0,k), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_minproduct_chemm(Magma_minproductRight, Magma_minproductUpper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_minproduct_ctrmm(Magma_minproductRight, Magma_minproductUpper, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                k, kb,
                                c_one, dB(k,k), lddb,
                                dA(0,k), ldda);
                }
                
                magma_minproduct_queue_sync( stream[2] );
                magma_minproduct_queue_sync( stream[0] );
                
                lapackf77_chegst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_minproduct_csetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_minproduct_queue_sync( stream[1] );
        }
        else {
            /* Compute L'*A*L */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_minproduct_cgetmatrix_async( kb, kb,
                                        dB(k, k), lddb,
                                        B(0, 0),  nb, stream[2] );
                
                /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    magma_minproduct_ctrmm(Magma_minproductRight, Magma_minproductLower, Magma_minproductNoTrans, Magma_minproductNonUnit,
                                kb, k,
                                c_one, dB(0,0), lddb,
                                dA(k,0), ldda);
                    
                    magma_minproduct_chemm(Magma_minproductLeft, Magma_minproductLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_minproduct_queue_sync( stream[1] );
                }
                
                magma_minproduct_cgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(0, 0),  lda, stream[0] );
                
                if (k > 0) {
                    magma_minproduct_cher2k(Magma_minproductLower, Magma_minproductConjTrans,
                                 k, kb,
                                 c_one, dA(k,0), ldda,
                                 dB(k,0), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_minproduct_chemm(Magma_minproductLeft, Magma_minproductLower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_minproduct_ctrmm(Magma_minproductLeft, Magma_minproductLower, Magma_minproductConjTrans, Magma_minproductNonUnit,
                                kb, k,
                                c_one, dB(k,k), lddb,
                                dA(k,0), ldda);
                }
                
                magma_minproduct_queue_sync( stream[2] );
                magma_minproduct_queue_sync( stream[0] );
                
                lapackf77_chegst( &itype, uplo_, &kb, A(0, 0), &lda, B(0, 0), &ldb, info);
                
                magma_minproduct_csetmatrix_async( kb, kb,
                                        A(0, 0),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_minproduct_queue_sync( stream[1] );
        }
    }
    magma_minproduct_queue_destroy( stream[0] );
    magma_minproduct_queue_destroy( stream[1] );
    magma_minproduct_queue_destroy( stream[2] );
    
    magma_minproduct_free_pinned( w );
    
    return *info;
} /* magma_minproduct_chegst_gpu */

#undef A
#undef B
#undef dA
#undef dB
