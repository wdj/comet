/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @precisions normal z -> s d c
 
*/
#include "common_magma_minproduct.h"
#define PRECISION_z

#define  A(i, j) ( A[(j)*lda  + (i)])
#define  C(i, j) ( C[(j)*ldc  + (i)])
#define  D(i)    ( D[(i)*incD] )

// trailing submatrix update with inner-blocking 
int zherk_d(magma_minproduct_uplo_t uplo, magma_minproduct_int_t m, magma_minproduct_int_t n,
            magma_minproductDoubleComplex alpha, magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
            magma_minproductDoubleComplex beta,  magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
            magma_minproductDoubleComplex *D, magma_minproduct_int_t incD)
{
    magma_minproductDoubleComplex *Aik;
    magma_minproductDoubleComplex *Dkk;
    magma_minproductDoubleComplex *Akj;

    /* Check input arguments */
    if ((uplo != Magma_minproductLower) && (uplo != Magma_minproductUpper)) {
        return -1;
    }
    if (m < 0) {
        return -3;
    }
    if (n < 0) {
        return -4;
    }
    if ((lda < max(1, m)) && (m > 0)) {
        return -7;
    }
    if ((ldc < max(1, m)) && (m > 0)) {
        return -10;
    }
    if ( incD < 0 ) {
        return -12;
    }

    /* Quick return */
    if (m == 0 || n == 0 ||
        ((alpha == 0.0 || m == 0) && beta == 1.0) ) {
        return MAGMA_minproduct_SUCCESS;
    }

    if ( uplo == Magma_minproductLower )
    {
        for(int j=0; j<m; j++)
        {
            for(int i=j; i<m; i++)
            {
                magma_minproductDoubleComplex tmp = MAGMA_minproduct_Z_ZERO;
                Aik = A+i;
                Dkk = D;
                Akj = A+j;
                for(int k=0; k<n; k++, Aik+=lda, Dkk+=incD, Akj+=lda )
                {
                    tmp += (*Aik) * (*Dkk) * conj( *Akj );
                }
                C(i, j) = beta * C(i, j) + alpha * tmp;
            }
        }
    }
    else
    {
        for(int j=0; j<m; j++)
        {
            for(int i=0; i<=j; i++)
            {
                magma_minproductDoubleComplex tmp = MAGMA_minproduct_Z_ZERO;
                for(int k=0; k<n; k++)
                {
                    tmp += A(i, k) * D( k ) * conj( A(k, j) );
                }
                C(i, j) = beta * C(i, j) + alpha * tmp;
            }
        }
    }
    return MAGMA_minproduct_SUCCESS;
}

// trailing submatrix update with inner-blocking, using workshpace that
// stores D*L'
int zherk_d_workspace(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t k,
                      magma_minproductDoubleComplex alpha, magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
                      magma_minproductDoubleComplex beta,  magma_minproductDoubleComplex *C, magma_minproduct_int_t ldc,
                      magma_minproductDoubleComplex *work, magma_minproduct_int_t ldw)
{
    magma_minproductDoubleComplex c_one  =  MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex c_mone = -MAGMA_minproduct_Z_ONE;

    /* Check input arguments */
    if ((uplo != Magma_minproductLower) && (uplo != Magma_minproductUpper)) {
        return -1;
    }
    if (n < 0) {
        return -2;
    }
    if (k < 0) {
        return -3;
    }
    if ((lda < max(1,n)) && (n > 0)) {
        return -6;
    }
    if ((ldc < max(1,n)) && (n > 0)) {
        return -9;
    }

    /* Quick return */
    if (n == 0 || k == 0 ||
        ((alpha == 0.0 || k == 0) && beta == 1.0) ) {
        return MAGMA_minproduct_SUCCESS;
    }

    if ( uplo == Magma_minproductLower )
    {
         blasf77_zgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, 
                        &n, &n, &k,
                        &c_mone, A,    &lda,
                                 work, &ldw,
                        &c_one,  C,    &ldc );
    }
    else
    {
         blasf77_zgemm( Magma_minproductNoTransStr, Magma_minproductNoTransStr, 
                        &n, &n, &k,
                        &c_mone, work, &ldw,
                                 A,    &lda,
                        &c_one,  C,    &ldc );
    }
    return MAGMA_minproduct_SUCCESS;
}

// diagonal factorization with inner-block
int zhetrf_diag_nopiv(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, 
                      magma_minproductDoubleComplex *A, magma_minproduct_int_t lda)
{
    /* Quick return */
    if (n == 1)
        return 0;
    if (lda < n) 
        return -1;

    /**/
    magma_minproduct_int_t info = 0, ione = 1;
    magma_minproductDoubleComplex *Ak1k = NULL;
    magma_minproductDoubleComplex Akk;
    double done = 1.0;
    double alpha;

    if ( uplo == Magma_minproductLower )
    {
        /* Diagonal element */
        Akk  = *A;

        /* Pointer on first extra diagonal element */
        Ak1k = A + 1;

        for (magma_minproduct_int_t k=n-1; k>0; k--) {
            if ( fabs(Akk) < lapackf77_dlamch("Epsilon") ) {
                info = k;
                return info;
            }

            // scale off-diagonals
            alpha = done / MAGMA_minproduct_Z_REAL( Akk );
            blasf77_zdscal(&k, &alpha, Ak1k, &ione);

            // update remaining
            alpha = - MAGMA_minproduct_Z_REAL( Akk );
            blasf77_zher(Magma_minproductLowerStr, &k, 
                         &alpha, Ak1k, &ione, Ak1k + lda, &lda);

            /* Move to next diagonal element */
            Ak1k += lda;
            Akk = *Ak1k;
            Ak1k++;
        }
    } else {
        /* Diagonal element */
        Akk  = *A;

        /* Pointer on first extra diagonal element */
        Ak1k = A + lda;

        for (magma_minproduct_int_t k=n-1; k>0; k--) {
            if ( fabs(Akk) < lapackf77_dlamch("Epsilon") ) {
                info = k;
                return info;
            }

            // scale off-diagonals
            alpha = done / MAGMA_minproduct_Z_REAL( Akk );
            blasf77_zdscal(&k, &alpha, Ak1k, &lda);

            // update remaining
            alpha = - MAGMA_minproduct_Z_REAL( Akk );

            #if defined(PRECISION_z) | defined(PRECISION_c)
            lapackf77_zlacgv(&k, Ak1k, &lda);
            #endif
            blasf77_zher(Magma_minproductUpperStr, &k, 
                         &alpha, Ak1k, &lda, Ak1k + 1, &lda);
            #if defined(PRECISION_z) | defined(PRECISION_c)
            lapackf77_zlacgv(&k, Ak1k, &lda);
            #endif

            /* Move to next diagonal element */
            Ak1k ++;
            Akk = *Ak1k;
            Ak1k += lda;
        }
    }
    return info;
}


// main routine
extern "C" magma_minproduct_int_t
zhetrf_nopiv_cpu(magma_minproduct_uplo_t uplo, magma_minproduct_int_t n, magma_minproduct_int_t ib,
                 magma_minproductDoubleComplex *A, magma_minproduct_int_t lda,
                 magma_minproduct_int_t *info)
{
    magma_minproduct_int_t ione = 1;
    double alpha;
    double done = 1.0;
    magma_minproductDoubleComplex zone  =  MAGMA_minproduct_Z_ONE;
    magma_minproductDoubleComplex mzone = -MAGMA_minproduct_Z_ONE;

    /* Check input arguments */
    if (lda < n) {
        *info = -1;
        return *info;
    }

    *info = 0;
    /* Quick return */
    if (n == 1) {
        return *info;
    }

    if ( uplo == Magma_minproductLower ) {
        for(magma_minproduct_int_t i = 0; i < n; i += ib) {
            magma_minproduct_int_t sb = min(n-i, ib);

            /* Factorize the diagonal block */
            *info = zhetrf_diag_nopiv(uplo, sb, &A(i, i), lda);
            if (*info != 0) return *info;

            if ( i + sb < n ) {
                magma_minproduct_int_t height = n - i - sb;

                /* Solve the lower panel ( L21*D11 )*/
                blasf77_ztrsm(
                    Magma_minproductRightStr, Magma_minproductLowerStr, 
                    Magma_minproductConjTransStr, Magma_minproductUnitStr,
                    &height, &sb, 
                    &zone, &A(i, i),    &lda,
                           &A(i+sb, i), &lda);

                /* Scale the block to divide by D */
                for (magma_minproduct_int_t k=0; k<sb; k++) {
                    #define ZHERK_D_WORKSPACE
                    #ifdef ZHERK_D_WORKSPACE
                    for (magma_minproduct_int_t ii=i+sb; ii<n; ii++) A(i+k, ii) = MAGMA_minproduct_Z_CNJG(A(ii, i+k));
                    #endif
                    alpha = done / MAGMA_minproduct_Z_REAL(A(i+k, i+k));
                    blasf77_zdscal(&height, &alpha, &A(i+sb, i+k), &ione);
                    A(i+k, i+k) = MAGMA_minproduct_Z_MAKE(MAGMA_minproduct_Z_REAL(A(i+k, i+k)), 0.0);
                }

                /* Update the trailing submatrix A22 = A22 - A21 * D11 * A21' */
                #ifdef ZHERK_D_WORKSPACE
                zherk_d_workspace(Magma_minproductLower, height, sb,
                                  mzone, &A(i+sb, i), lda,    // A21
                                  zone,  &A(i+sb, i+sb), lda, // A22
                                         &A(i, i+sb), lda);   // workspace, I am writing on upper part :)
                #else
                zherk_d(Magma_minproductLower, height, sb,
                        mzone, &A(i+sb, i), lda,    // A21
                        zone,  &A(i+sb, i+sb), lda, // A22
                               &A(i, i), lda+1);    // D11
                #endif
            }
        }
    } else {
        for(magma_minproduct_int_t i = 0; i < n; i += ib) {
            magma_minproduct_int_t sb = min(n-i, ib);

            /* Factorize the diagonal block */
            *info = zhetrf_diag_nopiv(uplo, sb, &A(i, i), lda);
            if (*info != 0) return *info;

            if ( i + sb < n ) {
                magma_minproduct_int_t height = n - i - sb;

                /* Solve the lower panel ( L21*D11 )*/
                blasf77_ztrsm(
                    Magma_minproductLeftStr, Magma_minproductUpperStr, 
                    Magma_minproductConjTransStr, Magma_minproductUnitStr,
                    &sb, &height, 
                    &zone, &A(i, i),    &lda,
                           &A(i, i+sb), &lda);

                /* Scale the block to divide by D */
                for (magma_minproduct_int_t k=0; k<sb; k++) {
                    #define ZHERK_D_WORKSPACE
                    #ifdef ZHERK_D_WORKSPACE
                    for (magma_minproduct_int_t ii=i+sb; ii<n; ii++) A(ii, i+k) = MAGMA_minproduct_Z_CNJG(A(i+k, ii));
                    #endif
                    alpha = done / MAGMA_minproduct_Z_REAL(A(i+k, i+k));
                    blasf77_zdscal(&height, &alpha, &A(i+k, i+sb), &lda);
                    A(i+k, i+k) = MAGMA_minproduct_Z_MAKE(MAGMA_minproduct_Z_REAL(A(i+k, i+k)), 0.0);
                }

                /* Update the trailing submatrix A22 = A22 - A21 * D11 * A21' */
                #ifdef ZHERK_D_WORKSPACE
                zherk_d_workspace(Magma_minproductUpper, height, sb,
                                  mzone, &A(i, i+sb), lda,    // A21
                                  zone,  &A(i+sb, i+sb), lda, // A22
                                         &A(i+sb, i), lda);   // workspace, I am writing on upper part :)
                #else
                zherk_d(Magma_minproductUpper, height, sb,
                        mzone, &A(i, i+sb), lda,    // A21
                        zone,  &A(i+sb, i+sb), lda, // A22
                               &A(i, i), lda+1);    // D11
                #endif
            }
        }
    }

    return *info;
}

