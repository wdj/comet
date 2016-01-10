/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       November 2011

       @generated from zhetrf_nopiv_tally4_cpu.cpp normal z -> c, Fri Jan 30 19:00:17 2015
 
*/
#include "common_magma_tally4.h"
#define PRECISION_c

#define  A(i, j) ( A[(j)*lda  + (i)])
#define  C(i, j) ( C[(j)*ldc  + (i)])
#define  D(i)    ( D[(i)*incD] )

// trailing submatrix update with inner-blocking 
int cherk_d(magma_tally4_uplo_t uplo, magma_tally4_int_t m, magma_tally4_int_t n,
            magma_tally4FloatComplex alpha, magma_tally4FloatComplex *A, magma_tally4_int_t lda,
            magma_tally4FloatComplex beta,  magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
            magma_tally4FloatComplex *D, magma_tally4_int_t incD)
{
    magma_tally4FloatComplex *Aik;
    magma_tally4FloatComplex *Dkk;
    magma_tally4FloatComplex *Akj;

    /* Check input arguments */
    if ((uplo != Magma_tally4Lower) && (uplo != Magma_tally4Upper)) {
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
        return MAGMA_tally4_SUCCESS;
    }

    if ( uplo == Magma_tally4Lower )
    {
        for(int j=0; j<m; j++)
        {
            for(int i=j; i<m; i++)
            {
                magma_tally4FloatComplex tmp = MAGMA_tally4_C_ZERO;
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
                magma_tally4FloatComplex tmp = MAGMA_tally4_C_ZERO;
                for(int k=0; k<n; k++)
                {
                    tmp += A(i, k) * D( k ) * conj( A(k, j) );
                }
                C(i, j) = beta * C(i, j) + alpha * tmp;
            }
        }
    }
    return MAGMA_tally4_SUCCESS;
}

// trailing submatrix update with inner-blocking, using workshpace that
// stores D*L'
int cherk_d_workspace(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t k,
                      magma_tally4FloatComplex alpha, magma_tally4FloatComplex *A, magma_tally4_int_t lda,
                      magma_tally4FloatComplex beta,  magma_tally4FloatComplex *C, magma_tally4_int_t ldc,
                      magma_tally4FloatComplex *work, magma_tally4_int_t ldw)
{
    magma_tally4FloatComplex c_one  =  MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex c_mone = -MAGMA_tally4_C_ONE;

    /* Check input arguments */
    if ((uplo != Magma_tally4Lower) && (uplo != Magma_tally4Upper)) {
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
        return MAGMA_tally4_SUCCESS;
    }

    if ( uplo == Magma_tally4Lower )
    {
         blasf77_cgemm( Magma_tally4NoTransStr, Magma_tally4NoTransStr, 
                        &n, &n, &k,
                        &c_mone, A,    &lda,
                                 work, &ldw,
                        &c_one,  C,    &ldc );
    }
    else
    {
         blasf77_cgemm( Magma_tally4NoTransStr, Magma_tally4NoTransStr, 
                        &n, &n, &k,
                        &c_mone, work, &ldw,
                                 A,    &lda,
                        &c_one,  C,    &ldc );
    }
    return MAGMA_tally4_SUCCESS;
}

// diagonal factorization with inner-block
int chetrf_diag_nopiv(magma_tally4_uplo_t uplo, magma_tally4_int_t n, 
                      magma_tally4FloatComplex *A, magma_tally4_int_t lda)
{
    /* Quick return */
    if (n == 1)
        return 0;
    if (lda < n) 
        return -1;

    /**/
    magma_tally4_int_t info = 0, ione = 1;
    magma_tally4FloatComplex *Ak1k = NULL;
    magma_tally4FloatComplex Akk;
    float done = 1.0;
    float alpha;

    if ( uplo == Magma_tally4Lower )
    {
        /* Diagonal element */
        Akk  = *A;

        /* Pointer on first extra diagonal element */
        Ak1k = A + 1;

        for (magma_tally4_int_t k=n-1; k>0; k--) {
            if ( fabs(Akk) < lapackf77_slamch("Epsilon") ) {
                info = k;
                return info;
            }

            // scale off-diagonals
            alpha = done / MAGMA_tally4_C_REAL( Akk );
            blasf77_csscal(&k, &alpha, Ak1k, &ione);

            // update remaining
            alpha = - MAGMA_tally4_C_REAL( Akk );
            blasf77_cher(Magma_tally4LowerStr, &k, 
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

        for (magma_tally4_int_t k=n-1; k>0; k--) {
            if ( fabs(Akk) < lapackf77_slamch("Epsilon") ) {
                info = k;
                return info;
            }

            // scale off-diagonals
            alpha = done / MAGMA_tally4_C_REAL( Akk );
            blasf77_csscal(&k, &alpha, Ak1k, &lda);

            // update remaining
            alpha = - MAGMA_tally4_C_REAL( Akk );

            #if defined(PRECISION_z) | defined(PRECISION_c)
            lapackf77_clacgv(&k, Ak1k, &lda);
            #endif
            blasf77_cher(Magma_tally4UpperStr, &k, 
                         &alpha, Ak1k, &lda, Ak1k + 1, &lda);
            #if defined(PRECISION_z) | defined(PRECISION_c)
            lapackf77_clacgv(&k, Ak1k, &lda);
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
extern "C" magma_tally4_int_t
chetrf_nopiv_tally4_cpu(magma_tally4_uplo_t uplo, magma_tally4_int_t n, magma_tally4_int_t ib,
                 magma_tally4FloatComplex *A, magma_tally4_int_t lda,
                 magma_tally4_int_t *info)
{
    magma_tally4_int_t ione = 1;
    float alpha;
    float done = 1.0;
    magma_tally4FloatComplex zone  =  MAGMA_tally4_C_ONE;
    magma_tally4FloatComplex mzone = -MAGMA_tally4_C_ONE;

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

    if ( uplo == Magma_tally4Lower ) {
        for(magma_tally4_int_t i = 0; i < n; i += ib) {
            magma_tally4_int_t sb = min(n-i, ib);

            /* Factorize the diagonal block */
            *info = chetrf_diag_nopiv(uplo, sb, &A(i, i), lda);
            if (*info != 0) return *info;

            if ( i + sb < n ) {
                magma_tally4_int_t height = n - i - sb;

                /* Solve the lower panel ( L21*D11 )*/
                blasf77_ctrsm(
                    Magma_tally4RightStr, Magma_tally4LowerStr, 
                    Magma_tally4ConjTransStr, Magma_tally4UnitStr,
                    &height, &sb, 
                    &zone, &A(i, i),    &lda,
                           &A(i+sb, i), &lda);

                /* Scale the block to divide by D */
                for (magma_tally4_int_t k=0; k<sb; k++) {
                    #define CHERK_D_WORKSPACE
                    #ifdef CHERK_D_WORKSPACE
                    for (magma_tally4_int_t ii=i+sb; ii<n; ii++) A(i+k, ii) = MAGMA_tally4_C_CNJG(A(ii, i+k));
                    #endif
                    alpha = done / MAGMA_tally4_C_REAL(A(i+k, i+k));
                    blasf77_csscal(&height, &alpha, &A(i+sb, i+k), &ione);
                    A(i+k, i+k) = MAGMA_tally4_C_MAKE(MAGMA_tally4_C_REAL(A(i+k, i+k)), 0.0);
                }

                /* Update the trailing submatrix A22 = A22 - A21 * D11 * A21' */
                #ifdef CHERK_D_WORKSPACE
                cherk_d_workspace(Magma_tally4Lower, height, sb,
                                  mzone, &A(i+sb, i), lda,    // A21
                                  zone,  &A(i+sb, i+sb), lda, // A22
                                         &A(i, i+sb), lda);   // workspace, I am writing on upper part :)
                #else
                cherk_d(Magma_tally4Lower, height, sb,
                        mzone, &A(i+sb, i), lda,    // A21
                        zone,  &A(i+sb, i+sb), lda, // A22
                               &A(i, i), lda+1);    // D11
                #endif
            }
        }
    } else {
        for(magma_tally4_int_t i = 0; i < n; i += ib) {
            magma_tally4_int_t sb = min(n-i, ib);

            /* Factorize the diagonal block */
            *info = chetrf_diag_nopiv(uplo, sb, &A(i, i), lda);
            if (*info != 0) return *info;

            if ( i + sb < n ) {
                magma_tally4_int_t height = n - i - sb;

                /* Solve the lower panel ( L21*D11 )*/
                blasf77_ctrsm(
                    Magma_tally4LeftStr, Magma_tally4UpperStr, 
                    Magma_tally4ConjTransStr, Magma_tally4UnitStr,
                    &sb, &height, 
                    &zone, &A(i, i),    &lda,
                           &A(i, i+sb), &lda);

                /* Scale the block to divide by D */
                for (magma_tally4_int_t k=0; k<sb; k++) {
                    #define CHERK_D_WORKSPACE
                    #ifdef CHERK_D_WORKSPACE
                    for (magma_tally4_int_t ii=i+sb; ii<n; ii++) A(ii, i+k) = MAGMA_tally4_C_CNJG(A(i+k, ii));
                    #endif
                    alpha = done / MAGMA_tally4_C_REAL(A(i+k, i+k));
                    blasf77_csscal(&height, &alpha, &A(i+k, i+sb), &lda);
                    A(i+k, i+k) = MAGMA_tally4_C_MAKE(MAGMA_tally4_C_REAL(A(i+k, i+k)), 0.0);
                }

                /* Update the trailing submatrix A22 = A22 - A21 * D11 * A21' */
                #ifdef CHERK_D_WORKSPACE
                cherk_d_workspace(Magma_tally4Upper, height, sb,
                                  mzone, &A(i, i+sb), lda,    // A21
                                  zone,  &A(i+sb, i+sb), lda, // A22
                                         &A(i+sb, i), lda);   // workspace, I am writing on upper part :)
                #else
                cherk_d(Magma_tally4Upper, height, sb,
                        mzone, &A(i, i+sb), lda,    // A21
                        zone,  &A(i+sb, i+sb), lda, // A22
                               &A(i, i), lda+1);    // D11
                #endif
            }
        }
    }

    return *info;
}
