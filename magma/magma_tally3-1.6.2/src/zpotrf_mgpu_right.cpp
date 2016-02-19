/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015
       
       @author Azzam Haidar
       @author Ichi Yamazaki

       @precisions normal z -> s d c

*/
#include "common_magma_tally3.h"
#include "trace.h"


/**
    Purpose
    -------
    ZPOTRF computes the Cholesky factorization of a complex Hermitian
    positive definite matrix dA.

    The factorization has the form
       dA = U**H * U,   if UPLO = Magma_tally3Upper, or
       dA = L  * L**H,  if UPLO = Magma_tally3Lower,
    where U is an upper triangular matrix and L is lower triangular.

    This is the block version of the algorithm, calling Level 3 BLAS.

    Arguments
    ---------
    @param[in]
    uplo    magma_tally3_uplo_t
      -     = Magma_tally3Upper:  Upper triangle of dA is stored;
      -     = Magma_tally3Lower:  Lower triangle of dA is stored.

    @param[in]
    n       INTEGER
            The order of the matrix dA.  N >= 0.

    @param[in,out]
    d_lA    COMPLEX_16 array of pointers on the GPU, dimension (ngpu)
            On entry, the Hermitian matrix dA distributed over GPUs
            (dl_A[d] points to the local matrix on the d-th GPU).
            It is distributed in 1D block column or row cyclic (with the
            block size of nb) if UPLO = Magma_tally3Upper or Magma_tally3Lower, respectively.
            If UPLO = Magma_tally3Upper, the leading N-by-N upper triangular
            part of dA contains the upper triangular part of the matrix dA,
            and the strictly lower triangular part of dA is not referenced.
            If UPLO = Magma_tally3Lower, the leading N-by-N lower triangular part
            of dA contains the lower triangular part of the matrix dA, and
            the strictly upper triangular part of dA is not referenced.
    \n
            On exit, if INFO = 0, the factor U or L from the Cholesky
            factorization dA = U**H * U or dA = L * L**H.

    @param[in]
    ldda     INTEGER
            The leading dimension of the array dA.  LDDA >= max(1,N).
            To benefit from coalescent memory accesses LDDA must be
            divisible by 16.

    @param[out]
    info    INTEGER
      -     = 0:  successful exit
      -     < 0:  if INFO = -i, the i-th argument had an illegal value
      -     > 0:  if INFO = i, the leading minor of order i is not
                  positive definite, and the factorization could not be
                  completed.

    @ingroup magma_tally3_zposv_comp
    ********************************************************************/
extern "C" magma_tally3_int_t
magma_tally3_zpotrf_mgpu_right(
    magma_tally3_int_t ngpu,
    magma_tally3_uplo_t uplo, magma_tally3_int_t n,
    magma_tally3DoubleComplex_ptr d_lA[], magma_tally3_int_t ldda,
    magma_tally3_int_t *info )
{
    #define dlA(id, i, j)  (d_lA[(id)] + (j) * ldda + (i))
    #define dlP(id, i, j)  (d_lP[(id)] + (j) * ldda + (i))

    #define panel(j)  (panel + (j))
    #define tmppanel(j)  (tmppanel + (j))
    #define tmpprevpanel(j)  (tmpprevpanel + (j))
    #define STREAM_ID(i) (nqueue > 1 ? 1+((i)/nb)%(nqueue-1) : 0)

    magma_tally3DoubleComplex z_one = MAGMA_tally3_Z_MAKE(  1.0, 0.0 );
    magma_tally3DoubleComplex mz_one = MAGMA_tally3_Z_MAKE( -1.0, 0.0 );
    double             one =  1.0;
    double             m_one = -1.0;
    const char* uplo_ = lapack_uplo_const_tally3( uplo );

    magma_tally3_int_t j, nb, d, id, j_local, blkid, crosspoint, prevtrsmrows=0, nqueue = 5;
    magma_tally3DoubleComplex *panel, *tmppanel0, *tmppanel1, *tmppanel, *tmpprevpanel;
    magma_tally3DoubleComplex *d_lP[Magma_tally3MaxGPUs], *dlpanel, *dlpanels[Magma_tally3MaxGPUs];
    magma_tally3_int_t rows, trsmrows, igpu, n_local[Magma_tally3MaxGPUs], ldpanel;
    magma_tally3_queue_t queues[Magma_tally3MaxGPUs][10];

    *info = 0;
    if ( uplo != Magma_tally3Upper && uplo != Magma_tally3Lower ) {
        *info = -1;
    } else if (n < 0) {
        *info = -2;
    } else if (ldda < max(1,n)) {
        *info = -4;
    }
    if (*info != 0) {
        magma_tally3_xerbla( __func__, -(*info) );
        return *info;
    }

    magma_tally3_device_t orig_dev;
    magma_tally3_getdevice( &orig_dev );
    magma_tally3_queue_t orig_stream;
    magma_tally3blasGetKernelStream( &orig_stream );

    nb = magma_tally3_get_zpotrf_nb(n);

    ldpanel = ldda;
    magma_tally3_setdevice(0);
    if (MAGMA_tally3_SUCCESS != magma_tally3_zmalloc_pinned( &panel, 2 * nb * ldpanel )) {
        *info = MAGMA_tally3_ERR_HOST_ALLOC;
        return *info;
    }

    tmppanel0 = panel;
    tmppanel1 = tmppanel0 + nb * ldpanel;

    if ((nb <= 1) || (nb >= n)) {
        // Use unblocked code.
        magma_tally3_zgetmatrix( n, n, dlA(0, 0, 0), ldda, panel, ldpanel);
        lapackf77_zpotrf( uplo_, &n, panel, &ldpanel, info);
        magma_tally3_zsetmatrix( n, n, panel, ldpanel, dlA(0, 0, 0), ldda );
    } else {
        for( d = 0; d < ngpu; d++ ) {
            // local-n and local-ld
            n_local[d] = ((n / nb) / ngpu) * nb;
            if (d < (n / nb) % ngpu)
                n_local[d] += nb;
            else if (d == (n / nb) % ngpu)
                n_local[d] += n % nb;

            magma_tally3_setdevice(d);
            magma_tally3_device_sync();
            if (MAGMA_tally3_SUCCESS != magma_tally3_zmalloc( &d_lP[d], nb * ldda )) {
                for( j = 0; j < d; j++ ) {
                    magma_tally3_setdevice(j);
                    magma_tally3_free( d_lP[d] );
                }
                *info = MAGMA_tally3_ERR_DEVICE_ALLOC;
                return *info;
            }
            for( j=0; j < nqueue; j++ ) {
                magma_tally3_queue_create( &queues[d][j] );
            }
        }

        //#define ENABLE_TIMER
        #if defined (ENABLE_TIMER)
        real_Double_t therk[4], tmtc, tcchol, tctrsm, tctm, tmnp, tcnp;
        real_Double_t ttot_herk[4] = {0,0,0,0}, ttot_mtc = 0, ttot_cchol = 0, ttot_ctrsm = 0, ttot_ctm = 0, ttot_mnp = 0, ttot_cnp = 0;
        printf("\n\n %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s %10s\n",
                "j", "nb", "row", "mtc", "CPU_np", "panel", "ctrsm", "CH+TRSM", "CPU", "dsyrk[0]", "dsyrk[1]", "dsyrk[2]", "dsyrk[3]", "ctm P", "gpu_np");
        printf("     ====================================================================================================\n");
        #endif

        // Use blocked code.
        if (uplo == Magma_tally3Upper) {
            printf( " === not supported, yet ===\n" );
        } else {
            blkid = -1;
            if (ngpu == 4)
                crosspoint = n;
            else if (ngpu == 3)
                crosspoint = n;
            else if (ngpu == 2)
                crosspoint = 20160;
            else
                crosspoint = 0;
            crosspoint = 0; //n; //n -- > gpu always does next panel, 0 --> cpu always does next panel
            crosspoint = n;

            #if defined (ENABLE_TIMER)
            real_Double_t tget = magma_tally3_wtime(), tset = 0.0, ttot = 0.0;
            #endif
            if ( n > nb ) {
                // send first panel to cpu
                magma_tally3_setdevice(0);
                tmppanel = tmppanel0;
                magma_tally3_zgetmatrix_async(n, nb,
                        dlA(0, 0, 0), ldda,
                        tmppanel(0),  ldpanel,
                        queues[0][0] );
            }
            #if defined (ENABLE_TIMER)
            for( d=0; d < ngpu; d++ ) {
                magma_tally3_setdevice(d);
                magma_tally3_device_sync();
            }
            tget = magma_tally3_wtime()-tget;
            #endif

            // Compute the Cholesky factorization A = L*L'
            for (j = 0; (j + nb) < n; j += nb) {
                #if defined (ENABLE_TIMER)
                therk[0] = therk[1] = therk[2] = therk[3] = tmtc = tcchol = tctrsm = tctm = tmnp = tcnp = 0.0;
                #endif

                blkid += 1;
                tmppanel = (blkid % 2 == 0) ? tmppanel0 : tmppanel1;
                // Set the gpu number that holds the current panel
                id = (j / nb) % ngpu;
                magma_tally3_setdevice(id);

                // Set the local index where the current panel is
                j_local = j / (nb * ngpu) * nb;
                
                rows = n - j;
                // Wait for the panel on cpu
                magma_tally3_queue_sync( queues[id][0] );
                if (j > 0 && prevtrsmrows > crosspoint) {
                    #if defined (ENABLE_TIMER)
                    tcnp = magma_tally3_wtime();
                    #endif

                    tmpprevpanel = ((blkid - 1) % 2) == 0 ? tmppanel0 : tmppanel1;

                    blasf77_zgemm( Magma_tally3NoTransStr, Magma_tally3ConjTransStr,
                            &rows, &nb, &nb,
                            &mz_one, tmpprevpanel(j), &ldpanel,
                                     tmpprevpanel(j), &ldpanel,
                            &z_one,      tmppanel(j), &ldpanel );

                    #if defined (ENABLE_TIMER)
                    tcnp = magma_tally3_wtime() - tcnp;
                    ttot_cnp += tcnp;
                    #endif
                }

                #if defined (ENABLE_TIMER)
                tcchol = magma_tally3_wtime();
                #endif
                lapackf77_zpotrf(Magma_tally3LowerStr, &nb, tmppanel(j), &ldpanel, info);
                if (*info != 0) {
                    *info = *info + j;
                    break;
                }

                #if defined (ENABLE_TIMER)
                tcchol = magma_tally3_wtime() - tcchol;
                ttot_cchol += tcchol;
                tctrsm = magma_tally3_wtime();
                #endif

                trsmrows = rows - nb;

                if (trsmrows > 0) {
                    blasf77_ztrsm(Magma_tally3RightStr, Magma_tally3LowerStr, Magma_tally3ConjTransStr, Magma_tally3NonUnitStr,
                                  &trsmrows, &nb,
                                  &z_one, tmppanel(j), &ldpanel,
                                          tmppanel(j + nb), &ldpanel);
                }

                #if defined (ENABLE_TIMER)
                tctrsm = magma_tally3_wtime() - tctrsm;
                ttot_ctrsm += tctrsm;
                tctm = magma_tally3_wtime();
                #endif

                d = (id + 1) % ngpu;
                // send current panel to gpus
                for (igpu = 0; igpu < ngpu; igpu++, d = (d + 1) % ngpu ) {
                    magma_tally3_int_t myrows = 0;
                    magma_tally3_int_t row_offset = 0;
                    if ( d == id ) {
                        dlpanel = dlA(d, j, j_local);
                        myrows = rows;
                        row_offset = 0;
                    } else {
                        dlpanel = dlP(d, 0, 0);
                        myrows = trsmrows;
                        row_offset = nb;
                    }

                    if (myrows > 0) {
                        magma_tally3_setdevice(d);
                        magma_tally3_zsetmatrix_async(myrows, nb,
                                tmppanel(j + row_offset),    ldpanel,
                                dlpanel, ldda, queues[d][0] );
                    }
                }
                /* make sure panel is on GPUs */
                d = (id + 1) % ngpu;
                for (igpu = 0; igpu < ngpu; igpu++, d = (d + 1) % ngpu ) {
                    magma_tally3_setdevice(d);
                    magma_tally3_queue_sync( queues[d][0] );
                }

                #if defined (ENABLE_TIMER)
                tctm = magma_tally3_wtime() - tctm;
                ttot_ctm += tctm;
                #endif

                if ( (j + nb) < n) {
                    magma_tally3_int_t offset = 0;
                    magma_tally3_int_t row_offset = 0;
                    if (j + nb + nb < n) {
                        d = (id + 1) % ngpu;
                        magma_tally3_setdevice(d);
                        magma_tally3_int_t j_local2 = (j + nb) / (nb * ngpu) * nb;
                        if (trsmrows <= crosspoint) {
                            #if defined (ENABLE_TIMER)
                            tmnp = magma_tally3_wtime();
                            #endif

                            // do gemm on look ahead panel
                            if ( d == id ) {
                                dlpanel = dlA(d, j + nb, j_local);
                            } else {
                                dlpanel = dlP(d, 0, 0);
                            }

                            magma_tally3blasSetKernelStream( queues[d][STREAM_ID(j_local2)] );
                            #define ZHERK_ON_DIAG
                            #ifdef  ZHERK_ON_DIAG
                            magma_tally3_zherk( Magma_tally3Lower, Magma_tally3NoTrans,
                                    nb, nb,
                                    m_one, dlpanel, ldda,
                                     one,  dlA(d, j + nb, j_local2), ldda);
                            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans,
                                    trsmrows-nb, nb, nb,
                                    mz_one, dlpanel+nb, ldda,
                                            dlpanel,    ldda,
                                     z_one, dlA(d, j + nb +nb, j_local2), ldda);
                            #else
                            magma_tally3_zgemm( Magma_tally3NoTrans, Magma_tally3ConjTrans,
                                    trsmrows, nb, nb,
                                    mz_one, dlpanel, ldda,
                                            dlpanel, ldda,
                                     z_one, dlA(d, j + nb, j_local2), ldda);
                            #endif

                            #if defined (ENABLE_TIMER)
                            magma_tally3_device_sync();
                            tmnp = magma_tally3_wtime() - tmnp;
                            ttot_mnp += tmnp;
                            #endif
                        }
                        // send next panel to cpu
                        magma_tally3_queue_sync( queues[d][STREAM_ID(j_local2)] ); // make sure lookahead is done
                        tmppanel = ((blkid+1) % 2 == 0) ? tmppanel0 : tmppanel1;
                        magma_tally3_zgetmatrix_async(rows-nb, nb,
                                dlA(d, j+nb, j_local2), ldda,
                                tmppanel(j+nb),  ldpanel,
                                queues[d][0] );
                        tmppanel = (blkid % 2 == 0) ? tmppanel0 : tmppanel1;

                        offset = j + nb + nb;
                        row_offset = nb;
                    } else {
                        offset = j + nb;
                        row_offset = 0;
                    }

                    if (n - offset > 0) {
                        // syrk on multiple gpu
                        for (d = 0; d < ngpu; d++ ) {
                            if ( d == id ) {
                                dlpanels[d] = dlA(d, j + nb + row_offset, j_local);
                            } else {
                                dlpanels[d] = dlP(d, row_offset, 0);
                            }
                        }

                        #if defined (ENABLE_TIMER)
                        for( d=0; d < ngpu; d++ ) therk[d] = magma_tally3_wtime();
                        #endif

                        //magma_tally3blasSetKernelStream( queues[d] );
                        //magma_tally3_zherk(Magma_tally3Lower, Magma_tally3NoTrans, n - offset, nb,
                        //        m_one, dlpanel, ldda,
                        //        one, &d_lA[d][offset + offset*ldda], ldda );
                        #ifdef  ZHERK_ON_DIAG
                        magma_tally3_zherk_mgpu
                        #else
                        magma_tally3_zherk_mgpu2
                        #endif
                                        (ngpu, Magma_tally3Lower, Magma_tally3NoTrans,
                                         nb, n - offset, nb,
                                         m_one, dlpanels, ldda, 0,
                                         one,   d_lA,     ldda, offset,
                                         nqueue, queues );
                        #if defined (ENABLE_TIMER)
                        for( d=0; d < ngpu; d++ ) {
                            magma_tally3_setdevice(d);
                            magma_tally3_device_sync();
                            therk[d] = magma_tally3_wtime() - therk[d];
                            ttot_herk[d] += therk[d];
                        }
                        #endif
                    }

                    prevtrsmrows = trsmrows;

                    #if defined (ENABLE_TIMER)
                    ttot += (tcnp+tcchol+tctrsm+therk[0]+therk[1]+therk[2]+tctm+tmnp);
                    printf("%10d %10d %10d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf(%d) %10.3lf\n",
                            j, nb, rows, tmtc,
                            tcnp,     // gemm
                            tcchol,   // potrf
                            tctrsm,   // trsm
                            (tcchol + tctrsm),
                            (tmtc+tcnp+tcchol+tctrsm),
                            therk[0], therk[1], therk[2], therk[3], // syrk
                            tctm, // copy panel to GPU
                            tmnp, // lookahead on GPU
                            (id + 1) % ngpu,
                            (tcnp+tcchol+tctrsm+therk[0]+therk[1]+therk[2]+tctm+tmnp));
                    fflush(0);
                    #endif
                }
            }
            for( d = 0; d < ngpu; d++ ) {
                magma_tally3_setdevice(d);
                for( id=0; id < nqueue; id++ ) {
                    magma_tally3_queue_sync( queues[d][id] );
                }
            }
            #if defined (ENABLE_TIMER)
            printf("\n%10d %10d %10d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf(-) %10.3lf\n",
                    n, n, 0, ttot_mtc,
                    ttot_cnp,     // gemm
                    ttot_cchol,   // potrf
                    ttot_ctrsm,   // trsm
                    (ttot_cchol + ttot_ctrsm),
                    (ttot_mtc+ttot_cnp+ttot_cchol+ttot_ctrsm),
                    ttot_herk[0], ttot_herk[1], ttot_herk[2], ttot_herk[3], // syrk
                    ttot_ctm, // copy panel to GPU
                    ttot_mnp, // lookahead on GPU
                    (ttot_cnp+ttot_cchol+ttot_ctrsm+ttot_herk[0]+ttot_herk[1]+ttot_herk[2]+ttot_ctm+ttot_mnp));
            printf("%10d %10d %10d %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf %10.3lf(-) %10.3lf (ratio)\n",
                    n, n, 0, ttot_mtc/ttot,
                    ttot_cnp/ttot,     // gemm
                    ttot_cchol/ttot,   // potrf
                    ttot_ctrsm/ttot,   // trsm
                    (ttot_cchol + ttot_ctrsm)/ttot,
                    (ttot_mtc+ttot_cnp+ttot_cchol+ttot_ctrsm)/ttot,
                    ttot_herk[0]/ttot, ttot_herk[1]/ttot, ttot_herk[2]/ttot, ttot_herk[3]/ttot, // syrk
                    ttot_ctm/ttot, // copy panel to GPU
                    ttot_mnp/ttot, // lookahead on GPU
                    (ttot_cnp+ttot_cchol+ttot_ctrsm+ttot_herk[0]+ttot_herk[1]+ttot_herk[2]+ttot_ctm+ttot_mnp)/ttot);
            #endif

            // cholesky for the last block
            if (j < n && *info == 0) {
                rows = n - j;
                id = (j / nb) % ngpu;

                // Set the local index where the current panel is
                j_local = j / (nb * ngpu) * nb;
                
                magma_tally3_setdevice(id);
                #if defined (ENABLE_TIMER)
                tset = magma_tally3_wtime();
                #endif
                magma_tally3_zgetmatrix(rows, rows, dlA(id, j, j_local), ldda, panel(j), ldpanel);
                lapackf77_zpotrf(Magma_tally3LowerStr, &rows, panel(j), &ldpanel, info);
                magma_tally3_zsetmatrix(rows, rows, panel(j), ldpanel, dlA(id, j, j_local), ldda);
                #if defined (ENABLE_TIMER)
                tset = magma_tally3_wtime() - tset;
                #endif
            }
            #if defined (ENABLE_TIMER)
            printf( " matrix_get,set: %10.3lf %10.3lf -> %10.3lf\n",tget,tset,ttot+tget+tset );
            #endif
        } // end of else not upper

        // clean up
        for( d = 0; d < ngpu; d++ ) {
            magma_tally3_setdevice(d);
            for( j=0; j < nqueue; j++ ) {
                magma_tally3_queue_destroy( queues[d][j] );
            }
            magma_tally3_free( d_lP[d] );
        }
    } // end of not lapack

    // free workspace
    magma_tally3_free_pinned( panel );
    magma_tally3_setdevice( orig_dev );
    magma_tally3blasSetKernelStream( orig_stream );

    return *info;
} /* magma_tally3_zpotrf_mgpu_right */
#undef dlA
#undef dlP
#undef panel
#undef tmppanel
#undef tmpprevpanel
#undef STREAM_ID

