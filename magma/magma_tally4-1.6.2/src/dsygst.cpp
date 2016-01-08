/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Raffaele Solca
       @author Azzam Haidar

       @generated from zhegst.cpp normal z -> d, Fri Jan 30 19:00:18 2015
*/

#include "common_magma_tally4.h"

/**
    Purpose
    -------
    DSYGST reduces a real symmetric-definite generalized
    eigenproblem to standard form.
    
    If ITYPE = 1, the problem is A*x = lambda*B*x,
    and A is overwritten by inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H)
    
    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or
    B*A*x = lambda*x, and A is overwritten by U*A*U**H or L**H*A*L.
    
    B must have been previously factorized as U**H*U or L*L**H by DPOTRF.
    
    Arguments
    ---------
    @param[in]
    itype   INTEGER
            = 1: compute inv(U**H)*A*inv(U) or inv(L)*A*inv(L**H);
            = 2 or 3: compute U*A*U**H or L**H*A*L.
    
    @param[in]
    uplo    magma_tally4_uplo_t
      -     = Magma_tally4Upper:  Upper triangle of A is stored and B is factored as U**H*U;
      -     = Magma_tally4Lower:  Lower triangle of A is stored and B is factored as L*L**H.
    
    @param[in]
    n       INTEGER
            The order of the matrices A and B.  N >= 0.
    
    @param[in,out]
    A       DOUBLE_PRECISION array, dimension (LDA,N)
            On entry, the symmetric matrix A.  If UPLO = Magma_tally4Upper, the leading
            N-by-N upper triangular part of A contains the upper
            triangular part of the matrix A, and the strictly lower
            triangular part of A is not referenced.  If UPLO = Magma_tally4Lower, the
            leading N-by-N lower triangular part of A contains the lower
            triangular part of the matrix A, and the strictly upper
            triangular part of A is not referenced.
    \n
            On exit, if INFO = 0, the transformed matrix, stored in the
            same format as A.
    
    @param[in]
    lda     INTEGER
            The leading dimension of the array A.  LDA >= max(1,N).
    
    @param[in]
    B       DOUBLE_PRECISION array, dimension (LDB,N)
            The triangular factor from the Cholesky factorization of B,
            as returned by DPOTRF.
    
    @param[in]
    ldb     INTEGER
            The leading dimension of the array B.  LDB >= max(1,N).
    
    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value

    @ingroup magma_tally4_dsyev_comp
    ********************************************************************/
extern "C" magma_tally4_int_t
magma_tally4_dsygst(
    magma_tally4_int_t itype, magma_tally4_uplo_t uplo, magma_tally4_int_t n,
    double *A, magma_tally4_int_t lda,
    double *B, magma_tally4_int_t ldb,
    magma_tally4_int_t *info)
{
#define A(i, j) (A + (j)*lda + (i))
#define B(i, j) (B + (j)*ldb + (i))

#define dA(i, j) (dw + (j)*ldda + (i))
#define dB(i, j) (dw + n*ldda + (j)*lddb + (i))

    const char* uplo_ = lapack_uplo_const( uplo );
    magma_tally4_int_t        nb;
    magma_tally4_int_t        k, kb, kb2;
    double    c_one      = MAGMA_tally4_D_ONE;
    double    c_neg_one  = MAGMA_tally4_D_NEG_ONE;
    double    c_half     = MAGMA_tally4_D_HALF;
    double    c_neg_half = MAGMA_tally4_D_NEG_HALF;
    double   *dw;
    magma_tally4_int_t        ldda = n;
    magma_tally4_int_t        lddb = n;
    double             d_one = 1.0;
    int upper = (uplo == Magma_tally4Upper);
    
    /* Test the input parameters. */
    *info = 0;
    if (itype < 1 || itype > 3) {
        *info = -1;
    } else if (! upper && uplo != Magma_tally4Lower) {
        *info = -2;
    } else if (n < 0) {
        *info = -3;
    } else if (lda < max(1,n)) {
        *info = -5;
    } else if (ldb < max(1,n)) {
        *info = -7;
    }
    if (*info != 0) {
        magma_tally4_xerbla( __func__, -(*info) );
        return *info;
    }
    
    /* Quick return */
    if ( n == 0 )
        return *info;
    
    if (MAGMA_tally4_SUCCESS != magma_tally4_dmalloc( &dw, 2*n*n )) {
        *info = MAGMA_tally4_ERR_DEVICE_ALLOC;
        return *info;
    }
    
    nb = magma_tally4_get_dsygst_nb(n);
    
    magma_tally4_queue_t stream[2];
    magma_tally4_queue_create( &stream[0] );
    magma_tally4_queue_create( &stream[1] );
    
    magma_tally4_dsetmatrix( n, n, A(0, 0), lda, dA(0, 0), ldda );
    magma_tally4_dsetmatrix( n, n, B(0, 0), ldb, dB(0, 0), lddb );
    
    /* Use hybrid blocked code */
    
    if (itype == 1) {
        if (upper) {
            /* Compute inv(U')*A*inv(U) */
            
            for (k = 0; k < n; k += nb) {
                kb = min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the upper triangle of A(k:n,k:n) */
                
                lapackf77_dsygst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
                
                magma_tally4_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                kb, n-k-kb,
                                c_one, dB(k,k), lddb,
                                dA(k,k+kb), ldda);
                    
                    magma_tally4_queue_sync( stream[0] );
                    
                    magma_tally4_dsymm(Magma_tally4Left, Magma_tally4Upper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_tally4_dsyr2k(Magma_tally4Upper, Magma_tally4ConjTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k,k+kb), ldda,
                                 dB(k,k+kb), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_tally4_dgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(k+kb, k+kb),  lda, stream[1] );
                    
                    magma_tally4_dsymm(Magma_tally4Left, Magma_tally4Upper,
                                kb, n-k-kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k,k+kb), lddb,
                                c_one, dA(k, k+kb), ldda);
                    
                    magma_tally4_dtrsm(Magma_tally4Right, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit,
                                kb, n-k-kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k,k+kb), ldda);
                    
                    magma_tally4_queue_sync( stream[1] );
                }
            }
            
            magma_tally4_queue_sync( stream[0] );
        }
        else {
            /* Compute inv(L)*A*inv(L') */
            
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                kb2= min(n-k-nb,nb);
                
                /* Update the lower triangle of A(k:n,k:n) */
                
                lapackf77_dsygst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
                
                magma_tally4_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[0] );
                
                if (k+kb < n) {
                    magma_tally4_dtrsm(Magma_tally4Right, Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                n-k-kb, kb,
                                c_one, dB(k,k), lddb,
                                dA(k+kb,k), ldda);
                    
                    magma_tally4_queue_sync( stream[0] );
                    
                    magma_tally4_dsymm(Magma_tally4Right, Magma_tally4Lower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_tally4_dsyr2k(Magma_tally4Lower, Magma_tally4NoTrans,
                                 n-k-kb, kb,
                                 c_neg_one, dA(k+kb,k), ldda,
                                 dB(k+kb,k), lddb,
                                 d_one, dA(k+kb,k+kb), ldda);
                    
                    magma_tally4_dgetmatrix_async( kb2, kb2,
                                            dA(k+kb, k+kb), ldda,
                                            A(k+kb, k+kb),  lda, stream[1] );
                    
                    magma_tally4_dsymm(Magma_tally4Right, Magma_tally4Lower,
                                n-k-kb, kb,
                                c_neg_half, dA(k,k), ldda,
                                dB(k+kb,k), lddb,
                                c_one, dA(k+kb, k), ldda);
                    
                    magma_tally4_dtrsm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4NonUnit,
                                n-k-kb, kb,
                                c_one, dB(k+kb,k+kb), lddb,
                                dA(k+kb,k), ldda);
                }
                
                magma_tally4_queue_sync( stream[1] );
            }
        }
        
        magma_tally4_queue_sync( stream[0] );
    }
    else {
        if (upper) {
            /* Compute U*A*U' */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_tally4_dgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(k, k),  lda, stream[0] );
                
                /* Update the upper triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    magma_tally4_dtrmm(Magma_tally4Left, Magma_tally4Upper, Magma_tally4NoTrans, Magma_tally4NonUnit,
                                k, kb,
                                c_one, dB(0,0), lddb,
                                dA(0,k), ldda);
                    
                    magma_tally4_dsymm(Magma_tally4Right, Magma_tally4Upper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_tally4_queue_sync( stream[1] );
                    
                    magma_tally4_dsyr2k(Magma_tally4Upper, Magma_tally4NoTrans,
                                 k, kb,
                                 c_one, dA(0,k), ldda,
                                 dB(0,k), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_tally4_dsymm(Magma_tally4Right, Magma_tally4Upper,
                                k, kb,
                                c_half, dA(k,k), ldda,
                                dB(0,k), lddb,
                                c_one, dA(0, k), ldda);
                    
                    magma_tally4_dtrmm(Magma_tally4Right, Magma_tally4Upper, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                k, kb,
                                c_one, dB(k,k), lddb,
                                dA(0,k), ldda);
                }
                
                magma_tally4_queue_sync( stream[0] );
                
                lapackf77_dsygst( &itype, uplo_, &kb, A(k, k), &lda, B(k, k), &ldb, info);
                
                magma_tally4_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_tally4_queue_sync( stream[1] );
        }
        else {
            /* Compute L'*A*L */
            for (k = 0; k < n; k += nb) {
                kb= min(n-k,nb);
                
                magma_tally4_dgetmatrix_async( kb, kb,
                                        dA(k, k), ldda,
                                        A(k, k),  lda, stream[0] );
                
                /* Update the lower triangle of A(1:k+kb-1,1:k+kb-1) */
                if (k > 0) {
                    
                    magma_tally4_dtrmm(Magma_tally4Right, Magma_tally4Lower, Magma_tally4NoTrans, Magma_tally4NonUnit,
                                kb, k,
                                c_one, dB(0,0), lddb,
                                dA(k,0), ldda);
                    
                    magma_tally4_dsymm(Magma_tally4Left, Magma_tally4Lower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_tally4_queue_sync( stream[1] );
                    
                    magma_tally4_dsyr2k(Magma_tally4Lower, Magma_tally4ConjTrans,
                                 k, kb,
                                 c_one, dA(k,0), ldda,
                                 dB(k,0), lddb,
                                 d_one, dA(0,0), ldda);
                    
                    magma_tally4_dsymm(Magma_tally4Left, Magma_tally4Lower,
                                kb, k,
                                c_half, dA(k,k), ldda,
                                dB(k,0), lddb,
                                c_one, dA(k, 0), ldda);
                    
                    magma_tally4_dtrmm(Magma_tally4Left, Magma_tally4Lower, Magma_tally4ConjTrans, Magma_tally4NonUnit,
                                kb, k,
                                c_one, dB(k,k), lddb,
                                dA(k,0), ldda);
                }
                
                magma_tally4_queue_sync( stream[0] );
                
                lapackf77_dsygst( &itype, uplo_, &kb, A(k,k), &lda, B(k,k), &ldb, info);
                
                magma_tally4_dsetmatrix_async( kb, kb,
                                        A(k, k),  lda,
                                        dA(k, k), ldda, stream[1] );
            }
            
            magma_tally4_queue_sync( stream[1] );
        }
    }
    
    magma_tally4_dgetmatrix( n, n, dA(0, 0), ldda, A(0, 0), lda );
    
    magma_tally4_queue_destroy( stream[0] );
    magma_tally4_queue_destroy( stream[1] );
    
    magma_tally4_free( dw );
    
    return *info;
} /* magma_tally4_dsygst_gpu */

#undef A
#undef B
#undef dA
#undef dB
