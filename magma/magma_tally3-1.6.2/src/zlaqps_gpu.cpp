/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"

#define PRECISION_z

/**
    @deprecated
    
    Purpose
    -------
    ZLAQPS computes a step of QR factorization with column pivoting
    of a complex M-by-N matrix A by using Blas-3.  It tries to factorize
    NB columns from A starting from the row OFFSET+1, and updates all
    of the matrix with Blas-3 xGEMM.

    In some cases, due to catastrophic cancellations, it cannot
    factorize NB columns.  Hence, the actual number of factorized
    columns is returned in KB.

    Block A(1:OFFSET,1:N) is accordingly pivoted, but not factorized.

    Arguments
    ---------
    @param[in]
    m       INTEGER
            The number of rows of the matrix A. M >= 0.

    @param[in]
    n       INTEGER
            The number of columns of the matrix A. N >= 0

    @param[in]
    offset  INTEGER
            The number of rows of A that have been factorized in
            previous steps.

    @param[in]
    nb      INTEGER
            The number of columns to factorize.

    @param[out]
    kb      INTEGER
            The number of columns actually factorized.

    @param[in,out]
    dA      COMPLEX_16 array, dimension (LDDA,N), on the GPU.
            On entry, the M-by-N matrix A.
            On exit, block A(OFFSET+1:M,1:KB) is the triangular
            factor obtained and block A(1:OFFSET,1:N) has been
            accordingly pivoted, but no factorized.
            The rest of the matrix, block A(OFFSET+1:M,KB+1:N) has
            been updated.

    @param[in]
    ldda    INTEGER
            The leading dimension of the array A. LDDA >= max(1,M).

    @param[in,out]
    jpvt    INTEGER array, dimension (N)
            JPVT(I) = K <==> Column K of the full matrix A has been
            permuted into position I in AP.

    @param[out]
    tau     COMPLEX_16 array, dimension (KB)
            The scalar factors of the elementary reflectors.

    @param[in,out]
    vn1     DOUBLE PRECISION array, dimension (N)
            The vector with the partial column norms.

    @param[in,out]
    vn2     DOUBLE PRECISION array, dimension (N)
            The vector with the exact column norms.

    @param[in,out]
    dauxv   COMPLEX_16 array, dimension (NB), on the GPU
            Auxiliary vector.

    @param[in,out]
    dF      COMPLEX_16 array, dimension (LDDF,NB), on the GPU
            Matrix F' = L*Y'*A.

    @param[in]
    lddf    INTEGER
            The leading dimension of the array F. LDDF >= max(1,N).

    @ingroup magma_tally3_zgeqp3_aux
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zlaqps_gpu(
    magma_tally3_int_t m, magma_tally3_int_t n, magma_tally3_int_t offset,
    magma_tally3_int_t nb, magma_tally3_int_t *kb,
    magma_tally3DoubleComplex_ptr dA,  magma_tally3_int_t ldda,
    magma_tally3_int_t *jpvt, magma_tally3DoubleComplex *tau,
    double *vn1, double *vn2,
    magma_tally3DoubleComplex_ptr dauxv,
    magma_tally3DoubleComplex_ptr dF,  magma_tally3_int_t lddf)
{
#define  dA(i, j) (dA  + (i) + (j)*(ldda))
#define  dF(i, j) (dF  + (i) + (j)*(lddf))

    magma_tally3DoubleComplex c_zero    = MAGMA_tally3_Z_MAKE( 0.,0.);
    magma_tally3DoubleComplex c_one     = MAGMA_tally3_Z_MAKE( 1.,0.);
    magma_tally3DoubleComplex c_neg_one = MAGMA_tally3_Z_MAKE(-1.,0.);
    magma_tally3_int_t ione = 1;
    
    magma_tally3_int_t i__1, i__2;
    //double d__1;
    magma_tally3DoubleComplex z__1;
    
    //magma_tally3_int_t j;
    magma_tally3_int_t k, rk;
    //magma_tally3DoubleComplex Akk;
    magma_tally3DoubleComplex_ptr dAks;
    magma_tally3DoubleComplex tauk = MAGMA_tally3_Z_ZERO;
    magma_tally3_int_t pvt;
    //double temp, temp2;
    double tol3z;
    magma_tally3_int_t itemp;

    double lsticc;
    magma_tally3Double_ptr dlsticcs;
    magma_tally3_dmalloc( &dlsticcs, 1+256*(n+255)/256 );

    //lastrk = min( m, n + offset );
    tol3z = magma_tally3_dsqrt( lapackf77_dlamch("Epsilon"));

    lsticc = 0;
    k = 0;
    magma_tally3_zmalloc( &dAks, nb );

    while( k < nb && lsticc == 0 ) {
        rk = offset + k;
        
        /* Determine ith pivot column and swap if necessary */
        // subtract 1 from Fortran/CUBLAS idamax; pvt, k are 0-based.
        pvt = k + magma_tally3_idamax( n-k, &vn1[k], ione ) - 1;
        
        if (pvt != k) {
            /*if (pvt >= nb) {
                // 1. Start copy from GPU
                magma_tally3_zgetmatrix_async( m - offset - nb, 1,
                                        dA(offset + nb, pvt), ldda,
                                        A (offset + nb, pvt), lda, stream );
            }*/

            /* F gets swapped so F must be sent at the end to GPU   */
            i__1 = k;
            /*if (pvt < nb) {
                // no need of transfer if pivot is within the panel
                blasf77_zswap( &m, A(0, pvt), &ione, A(0, k), &ione );
            }
            else {
                // 1. Finish copy from GPU
                magma_tally3_queue_sync( stream );

                // 2. Swap as usual on CPU
                blasf77_zswap(&m, A(0, pvt), &ione, A(0, k), &ione);

                // 3. Restore the GPU
                magma_tally3_zsetmatrix_async( m - offset - nb, 1,
                                        A (offset + nb, pvt), lda,
                                        dA(offset + nb, pvt), ldda, stream);
            }*/
            magma_tally3blas_zswap( m, dA(0, pvt), ione, dA(0, k), ione );

            //blasf77_zswap( &i__1, F(pvt,0), &ldf, F(k,0), &ldf );
            magma_tally3blas_zswap( i__1, dF(pvt, 0), lddf, dF(k, 0), lddf);
            itemp     = jpvt[pvt];
            jpvt[pvt] = jpvt[k];
            jpvt[k]   = itemp;
            //vn1[pvt] = vn1[k];
            //vn2[pvt] = vn2[k];
            #if defined(PRECISION_d) || defined(PRECISION_z)
                //magma_tally3_dswap( 1, &vn1[pvt], 1, &vn1[k], 1 );
                //magma_tally3_dswap( 1, &vn2[pvt], 1, &vn2[k], 1 );
                magma_tally3_dswap( 2, &vn1[pvt], n+offset, &vn1[k], n+offset );
            #else
                //magma_tally3_sswap( 1, &vn1[pvt], 1, &vn1[k], 1 );
                //magma_tally3_sswap( 1, &vn2[pvt], 1, &vn2[k], 1 );
                magma_tally3_sswap(2, &vn1[pvt], n+offset, &vn1[k], n+offset);
            #endif
        }

        /* Apply previous Householder reflectors to column K:
           A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)'.
           Optimization: multiply with beta=0; wait for vector and subtract */
        if (k > 0) {
            /*#if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 0; j < k; ++j) {
                *F(k,j) = MAGMA_tally3_Z_CNJG( *F(k,j) );
            }
            #endif*/

//#define RIGHT_UPDATE
#ifdef RIGHT_UPDATE
            i__1 = m - offset - nb;
            i__2 = k;
            magma_tally3_zgemv( Magma_tally3NoTrans, i__1, i__2,
                         c_neg_one, A(offset+nb, 0), lda,
                                    F(k,         0), ldf,
                         c_one,     A(offset+nb, k), ione );
#else
            i__1 = m - rk;
            i__2 = k;
            /*blasf77_zgemv( Magma_tally3NoTransStr, &i__1, &i__2,
                           &c_neg_one, A(rk, 0), &lda,
                                       F(k,  0), &ldf,
                           &c_one,     A(rk, k), &ione ); */
            magma_tally3_zgemv( Magma_tally3NoTrans, i__1, i__2,
                         c_neg_one, dA(rk, 0), ldda,
                                    dF(k,  0), lddf,
                         c_one,     dA(rk, k), ione );
#endif

            /*#if (defined(PRECISION_c) || defined(PRECISION_z))
            for (j = 0; j < k; ++j) {
                *F(k,j) = MAGMA_tally3_Z_CNJG( *F(k,j) );
            }
            #endif*/
        }
        
        /*  Generate elementary reflector H(k). */
        magma_tally3_zlarfg_gpu( m-rk, dA(rk, k), dA(rk + 1, k), &tau[k], &vn1[k], &dAks[k]);

        //Akk = *A(rk, k);
        //*A(rk, k) = c_one;
        //magma_tally3_zgetvector( 1, &dAks[k],  1, &Akk,     1 );

        /* needed to avoid the race condition */
        if (k == 0) magma_tally3_zsetvector(  1,    &c_one,        1, dA(rk, k), 1 );
        else        magma_tally3_zcopymatrix( 1, 1, dA(offset, 0), 1, dA(rk, k), 1 );

        /* Compute Kth column of F:
           Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K) on the GPU */
        if (k < n-1 || k > 0) magma_tally3_zgetvector( 1, &tau[k], 1, &tauk, 1 );
        if (k < n-1) {
            i__1 = m - rk;
            i__2 = n - k - 1;

            /* Send the vector to the GPU */
            //magma_tally3_zsetmatrix( i__1, 1, A(rk, k), lda, dA(rk,k), ldda );

            /* Multiply on GPU */
            // was CALL ZGEMV( 'Conjugate transpose', M-RK+1, N-K,
            //                 TAU( K ), A( RK,  K+1 ), LDA,
            //                           A( RK,  K   ), 1,
            //                 CZERO,    F( K+1, K   ), 1 )
            //magma_tally3_zgetvector( 1, &tau[k], 1, &tauk, 1 );
            magma_tally3_zgemv( Magma_tally3ConjTrans, m-rk, n-k-1,
                         tauk,   dA( rk,  k+1 ), ldda,
                                 dA( rk,  k   ), 1,
                         c_zero, dF( k+1, k   ), 1 );
            //magma_tally3_zscal( m-rk, tau[k], F( k+1, k), 1 );
            //magma_tally3_int_t i__3 = nb-k-1;
            //magma_tally3_int_t i__4 = i__2 - i__3;
            //magma_tally3_int_t i__5 = nb-k;
            //magma_tally3_zgemv( Magma_tally3ConjTrans, i__1 - i__5, i__2 - i__3,
            //             tau[k], dA(rk +i__5, k+1+i__3), ldda,
            //                     dA(rk +i__5, k       ), ione,
            //             c_zero, dF(k+1+i__3, k       ), ione );
            
            //magma_tally3_zgetmatrix_async( i__2-i__3, 1,
            //                        dF(k + 1 +i__3, k), i__2,
            //                        F (k + 1 +i__3, k), i__2, stream );
            
            //blasf77_zgemv( Magma_tally3ConjTransStr, &i__1, &i__3,
            //               &tau[k], A(rk,  k+1), &lda,
            //                        A(rk,  k  ), &ione,
            //               &c_zero, F(k+1, k  ), &ione );
            
            //magma_tally3_queue_sync( stream );
            //blasf77_zgemv( Magma_tally3ConjTransStr, &i__5, &i__4,
            //               &tau[k], A(rk, k+1+i__3), &lda,
            //                        A(rk, k       ), &ione,
            //               &c_one,  F(k+1+i__3, k ), &ione );
        }
        
        /* Padding F(1:K,K) with zeros.
        for (j = 0; j <= k; ++j) {
            magma_tally3_zsetvector( 1, &c_zero, 1, F(j, k), 1 );
        }*/
        
        /* Incremental updating of F:
           F(1:N,K) := F(1:N,K)                        - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'*A(RK:M,K).
           F(1:N,K) := tau(K)*A(RK:M,K+1:N)'*A(RK:M,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)'*A(RK:M,K)
                    := tau(K)(A(RK:M,K+1:N)' - F(1:N,1:K-1)*A(RK:M,1:K-1)') A(RK:M,K)
           so, F is (updated A)*V */
        //if (k > 0 && k < n-1) {
        if (k > 0) {
            //magma_tally3_zgetvector( 1, &tau[k], 1, &tauk, 1 );
            z__1 = MAGMA_tally3_Z_NEGATE( tauk );
#ifdef RIGHT_UPDATE
            i__1 = m - offset - nb;
            i__2 = k;
            magma_tally3_zgemv( Magma_tally3ConjTrans, i__1, i__2,
                         z__1,   dA(offset+nb, 0), lda,
                                 dA(offset+nb, k), ione,
                         c_zero, dauxv, ione );
            
            i__1 = k;
            magma_tally3_zgemv( Magma_tally3NoTrans, n-k-1, i__1,
                         c_one, F(k+1,0), ldf,
                                dauxv,     ione,
                         c_one, F(k+1,k), ione );
#else
            i__1 = m - rk;
            i__2 = k;
            //blasf77_zgemv( Magma_tally3ConjTransStr, &i__1, &i__2,
            //               &z__1,   A(rk, 0), &lda,
            //                        A(rk, k), &ione,
            //               &c_zero, auxv, &ione );

            magma_tally3_zgemv( Magma_tally3ConjTrans, i__1, i__2,
                         z__1,   dA(rk, 0), ldda,
                                 dA(rk, k), ione,
                         c_zero, dauxv, ione );
            
            //i__1 = k;
            //blasf77_zgemv( Magma_tally3NoTransStr, &n, &i__1,
            //               &c_one, F(0,0), &ldf,
            //                       auxv,   &ione,
            //               &c_one, F(0,k), &ione );
            /*magma_tally3_zgemv( Magma_tally3NoTrans, n, i__1,
                           c_one, F(0,0), ldf,
                                  auxv,   ione,
                           c_one, F(0,k), ione ); */
            /* I think we only need stricly lower-triangular part :) */
            magma_tally3_zgemv( Magma_tally3NoTrans, n-k-1, i__2,
                         c_one, dF(k+1,0), lddf,
                                dauxv,     ione,
                         c_one, dF(k+1,k), ione );
#endif
        }
        
        /* Optimization: On the last iteration start sending F back to the GPU */
        
        /* Update the current row of A:
           A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)'.               */
        if (k < n-1) {
            i__1 = n - k - 1;
            i__2 = k + 1;
            //blasf77_zgemm( Magma_tally3NoTransStr, Magma_tally3ConjTransStr, &ione, &i__1, &i__2,
            //               &c_neg_one, A(rk, 0  ), &lda,
            //                           F(k+1,0  ), &ldf,
            //               &c_one,     A(rk, k+1), &lda );
#ifdef RIGHT_UPDATE
            /* right-looking update of rows,                     */
            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, nb-k, i__1, ione,
                         c_neg_one, dA(rk,  k  ), ldda,
                                    dF(k+1, k  ), lddf,
                         c_one,     dA(rk,  k+1), ldda );
#else
            /* left-looking update of rows,                     *
             * since F=A'v with original A, so no right-looking */
            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, ione, i__1, i__2,
                         c_neg_one, dA(rk, 0  ), ldda,
                                    dF(k+1,0  ), lddf,
                         c_one,     dA(rk, k+1), ldda );
#endif
        }
        
        /* Update partial column norms. */
        if (rk < min(m, n+offset)-1 ) {
            magma_tally3blas_dznrm2_row_check_adjust(n-k-1, tol3z, &vn1[k+1], &vn2[k+1], dA(rk,k+1), ldda, dlsticcs);

            magma_tally3_device_sync();
            #if defined(PRECISION_d) || defined(PRECISION_z)
            magma_tally3_dgetvector( 1, &dlsticcs[0], 1, &lsticc, 1 );
            #else
            magma_tally3_sgetvector( 1, &dlsticcs[0], 1, &lsticc, 1 );
            #endif
        }


        /*if (rk < lastrk) {
            for (j = k + 1; j < n; ++j) {
                if (vn1[j] != 0.) {
                    // NOTE: The following 4 lines follow from the analysis in
                    //   Lapack Working Note 176.
                    temp = MAGMA_tally3_Z_ABS( *A(rk,j) ) / vn1[j];
                    temp = max( 0., ((1. + temp) * (1. - temp)) );

                    d__1 = vn1[j] / vn2[j];
                    temp2 = temp * (d__1 * d__1);

                    if (temp2 <= tol3z) {
                        vn2[j] = (double) lsticc;
                        lsticc = j;
                    } else {
                        vn1[j] *= magma_tally3_dsqrt(temp);
                    }
                }
            }
        }*/
        
        //*A(rk, k) = Akk;
        //magma_tally3_zsetvector( 1, &Akk, 1, A(rk, k), 1 );
        //magma_tally3_zswap( 1, &dAks[k], 1, A(rk, k), 1 );
        
        ++k;
    }
    magma_tally3_zcopymatrix( 1, k, dAks, 1, dA(offset, 0), ldda+1 );

    // leave k as the last column done
    --k;
    *kb = k + 1;
    rk = offset + *kb - 1;

    /* Apply the block reflector to the rest of the matrix:
       A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) - A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)'  */
    if (*kb < min(n, m - offset)) {
        i__1 = m - rk - 1;
        i__2 = n - *kb;
        
        /* Send F to the GPU
        magma_tally3_zsetmatrix( i__2, *kb,
                          F (*kb, 0), ldf,
                          dF(*kb, 0), i__2 ); */

        magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans, i__1, i__2, *kb,
                     c_neg_one, dA(rk+1, 0  ), ldda,
                                dF(*kb,  0  ), lddf,
                     c_one,     dA(rk+1, *kb), ldda );
    }
    /* Recomputation of difficult columns. */
    if ( lsticc > 0 ) {
        // printf( " -- recompute dnorms --\n" );
        magma_tally3blas_dznrm2_check( m-rk-1, n-*kb, dA(rk+1,*kb), ldda,
                                &vn1[*kb], dlsticcs );
        magma_tally3_dcopymatrix( n-*kb, 1, &vn1[*kb], *kb, &vn2[*kb], *kb );
    /*while( lsticc > 0 ) {
        itemp = (magma_tally3_int_t)(vn2[lsticc] >= 0. ? floor(vn2[lsticc] + .5) : -floor(.5 - vn2[lsticc]));
        i__1 = m - rk - 1;
        if (lsticc <= nb)
            vn1[lsticc] = magma_tally3_cblas_dznrm2( i__1, A(rk+1,lsticc), ione );
        else {
            // Where is the data, CPU or GPU ?
            double r1, r2;
            
            r1 = magma_tally3_cblas_dznrm2( nb-k, A(rk+1,lsticc), ione );
            r2 = magma_tally3_dznrm2(m-offset-nb, dA(offset + nb + 1, lsticc), ione);
            
            vn1[lsticc] = magma_tally3_dsqrt(r1*r1+r2*r2);
        }
        
        // NOTE: The computation of VN1( LSTICC ) relies on the fact that
        //   SNRM2 does not fail on vectors with norm below the value of SQRT(DLAMCH('S'))
        vn2[lsticc] = vn1[lsticc];
        lsticc = itemp; */
    }
    magma_tally3_free(dAks);
    magma_tally3_free(dlsticcs);

    return MAGMA_tally3_SUCCESS;
} /* magma_tally3_zlaqps */
