!
!   -- MAGMA_tally4 (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!
!   @precisions normal z -> c d s
!

#define PRECISION_z

module magma_tally4_zfortran

  use magma_tally4_param, only: sizeof_complex_16

  implicit none

  !---- Fortran interfaces to MAGMA_tally4 subroutines ----
  interface

     subroutine magma_tally4f_zgetptr( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       double precision:: d(*)
       double precision:: e(*)
       complex*16    :: tauq(*)
       complex*16    :: taup(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgetptr


     subroutine magma_tally4f_zgebrd( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       double precision:: d(*)
       double precision:: e(*)
       complex*16    :: tauq(*)
       complex*16    :: taup(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgebrd

     subroutine magma_tally4f_zgehrd2(n, ilo, ihi,A, lda, tau, work, lwork, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgehrd2

     subroutine magma_tally4f_zgehrd(n, ilo, ihi,A, lda, tau, work, lwork, d_T, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       complex*16    :: d_T(*)
       integer       :: info
     end subroutine magma_tally4f_zgehrd

     subroutine magma_tally4f_zgelqf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgelqf

     subroutine magma_tally4f_zgeqlf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgeqlf

     subroutine magma_tally4f_zgeqrf( m, n, A, lda, tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgeqrf

     subroutine magma_tally4f_zgesv(  n, nrhs, A, lda, ipiv, B, ldb, info)
       integer       :: n
       integer       :: nrhs
       complex*16    :: A
       integer       :: lda
       integer       :: ipiv(*)
       complex*16    :: B
       integer       :: ldb
       integer       :: info
     end subroutine magma_tally4f_zgesv

     subroutine magma_tally4f_zgetrf( m, n, A, lda, ipiv, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magma_tally4f_zgetrf

     subroutine magma_tally4f_zposv(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character     :: uplo
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_tally4f_zposv
     
     subroutine magma_tally4f_zpotrf( uplo, n, A, lda, info)
       character          :: uplo
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       integer       :: info
     end subroutine magma_tally4f_zpotrf

     subroutine magma_tally4f_zhetrd( uplo, n, A, lda, d, e, tau, work, lwork, info)
       character          :: uplo
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       double precision:: d(*)
       double precision:: e(*)
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zhetrd

     subroutine magma_tally4f_zunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: c(*)
       integer       :: ldc
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zunmqr

     subroutine magma_tally4f_zunmtr( side, uplo, trans, m, n, a, lda,tau,c,    ldc,work, lwork,info)
       character          :: side
       character          :: uplo
       character          :: trans
       integer       :: m
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: c(*)
       integer       :: ldc
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zunmtr
#if defined(PRECISION_z) || defined(PRECISION_c)

     subroutine magma_tally4f_zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: w(*)
       complex*16    :: vl(*)
       integer       :: ldvl
       complex*16    :: vr(*)
       integer       :: ldvr
       complex*16    :: work(*)
       integer       :: lwork
       double precision:: rwork(*)
       integer       :: info
     end subroutine magma_tally4f_zgeev

     subroutine magma_tally4f_zgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       double precision:: s(*)
       complex*16    :: u(*)
       integer       :: ldu
       complex*16    :: vt(*)
       integer       :: ldvt
       complex*16    :: work(*)
       integer       :: lwork
       double precision:: rwork(*)
       integer       :: info
     end subroutine magma_tally4f_zgesvd

     subroutine magma_tally4f_zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       double precision:: w(*)
       complex*16    :: work(*)
       integer       :: lwork
       double precision:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magma_tally4f_zheevd

     subroutine magma_tally4f_zhegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: b(*)
       integer       :: ldb
       double precision:: w(*)
       complex*16    :: work(*)
       integer       :: lwork
       double precision:: rwork(*)
       integer       :: lrwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magma_tally4f_zhegvd

#else
     subroutine magma_tally4f_zgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
       character          :: jobvl
       character          :: jobvr
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: wr(*)
       complex*16    :: wi(*)
       complex*16    :: vl(*)
       integer       :: ldvl
       complex*16    :: vr(*)
       integer       :: ldvr
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgeev

     subroutine magma_tally4f_zgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
       character          :: jobu
       character          :: jobvt
       integer       :: m
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       double precision:: s(*)
       complex*16    :: u(*)
       integer       :: ldu
       complex*16    :: vt(*)
       integer       :: ldvt
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgesvd

     subroutine magma_tally4f_zheevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
       character          :: jobz
       character          :: uplo
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       double precision:: w(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magma_tally4f_zheevd

     subroutine magma_tally4f_zhegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
       integer       :: itype
       character     :: jobz
       character     :: uplo
       integer       :: n
       complex*16    :: a(*)
       integer       :: lda
       complex*16    :: b(*)
       integer       :: ldb
       double precision:: w(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: iwork(*)
       integer       :: liwork
       integer       :: info
     end subroutine magma_tally4f_zhegvd
#endif

     subroutine magma_tally4f_zgels_gpu(  trans, m, n, nrhs, dA, ldda, dB, lddb, hwork, lwork, info)
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_tally4f_zgels_gpu

     subroutine magma_tally4f_zgeqrf_gpu( m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       magma_tally4_devptr_t:: dT
       integer       :: info
     end subroutine magma_tally4f_zgeqrf_gpu

     subroutine magma_tally4f_zgeqrf2_gpu(m, n, dA, ldda, tau, info)
       integer       :: m
       integer       :: n
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       integer       :: info
     end subroutine magma_tally4f_zgeqrf2_gpu

     subroutine magma_tally4f_zgeqrf3_gpu(m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       magma_tally4_devptr_t:: dT
       integer       :: info
     end subroutine magma_tally4f_zgeqrf3_gpu

     subroutine magma_tally4f_zgeqrs_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau
       magma_tally4_devptr_t:: dT
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magma_tally4f_zgeqrs_gpu

     subroutine magma_tally4f_zgeqrs3_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau
       magma_tally4_devptr_t:: dT
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magma_tally4f_zgeqrs3_gpu

     subroutine magma_tally4f_zgessm_gpu( storev, m, n, k, ib, ipiv, dL1, lddl1, dL,  lddl, dA,  ldda, info)
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       integer       :: ib
       integer       :: ipiv(*)
       magma_tally4_devptr_t:: dL1
       integer       :: lddl1
       magma_tally4_devptr_t:: dL
       integer       :: lddl
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magma_tally4f_zgessm_gpu

     subroutine magma_tally4f_zgesv_gpu(  n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_tally4f_zgesv_gpu

     subroutine magma_tally4f_zgetrf_gpu( m, n, dA, ldda, ipiv, info)
       integer       :: m
       integer       :: n
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magma_tally4f_zgetrf_gpu

     subroutine magma_tally4f_zgetrs_gpu( trans, n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       character          :: trans
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_tally4f_zgetrs_gpu

     subroutine magma_tally4f_zlabrd_gpu( m, n, nb, a, lda, da, ldda, d, e, tauq, taup, x, ldx, dx, lddx, y, ldy, dy, lddy)
       integer       :: m
       integer       :: n
       integer       :: nb
       complex*16    :: a(*)
       integer       :: lda
       magma_tally4_devptr_t:: da
       integer       :: ldda
       double precision:: d(*)
       double precision:: e(*)
       complex*16    :: tauq(*)
       complex*16    :: taup(*)
       complex*16    :: x(*)
       integer       :: ldx
       magma_tally4_devptr_t:: dx
       integer       :: lddx
       complex*16    :: y(*)
       integer       :: ldy
       magma_tally4_devptr_t:: dy
       integer       :: lddy
     end subroutine magma_tally4f_zlabrd_gpu

     subroutine magma_tally4f_zlarfb_gpu( side, trans, direct, storev, m, n, k, dv, ldv, dt, ldt, dc, ldc, dowrk, ldwork)
       character          :: side
       character          :: trans
       character          :: direct
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       magma_tally4_devptr_t:: dv
       integer       :: ldv
       magma_tally4_devptr_t:: dt
       integer       :: ldt
       magma_tally4_devptr_t:: dc
       integer       :: ldc
       magma_tally4_devptr_t:: dowrk
       integer       :: ldwork
     end subroutine magma_tally4f_zlarfb_gpu

     subroutine magma_tally4f_zposv_gpu(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_tally4f_zposv_gpu

     subroutine magma_tally4f_zpotrf_gpu( uplo, n, dA, ldda, info)
       character          :: uplo
       integer       :: n
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magma_tally4f_zpotrf_gpu

     subroutine magma_tally4f_zpotrs_gpu( uplo,  n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_tally4_devptr_t:: dA
       integer       :: ldda
       magma_tally4_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_tally4f_zpotrs_gpu

     subroutine magma_tally4f_zssssm_gpu( storev, m1, n1, m2, n2, k, ib, dA1, ldda1, dA2, ldda2, dL1, lddl1, dL2, lddl2, IPIV, info)
       character          :: storev
       integer       :: m1
       integer       :: n1
       integer       :: m2
       integer       :: n2
       integer       :: k
       integer       :: ib
       magma_tally4_devptr_t:: dA1
       integer       :: ldda1
       magma_tally4_devptr_t:: dA2
       integer       :: ldda2
       magma_tally4_devptr_t:: dL1
       integer       :: lddl1
       magma_tally4_devptr_t:: dL2
       integer       :: lddl2
       integer       :: IPIV(*)
       integer       :: info
     end subroutine magma_tally4f_zssssm_gpu

     subroutine magma_tally4f_zungqr_gpu( m, n, k, da, ldda, tau, dwork, nb, info)
       integer       :: m
       integer       :: n
       integer       :: k
       magma_tally4_devptr_t:: da
       integer       :: ldda
       complex*16    :: tau(*)
       magma_tally4_devptr_t:: dwork
       integer       :: nb
       integer       :: info
     end subroutine magma_tally4f_zungqr_gpu

     subroutine magma_tally4f_zunmqr_gpu( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, td, nb, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       magma_tally4_devptr_t:: a
       integer       :: lda
       complex*16    :: tau(*)
       magma_tally4_devptr_t:: c
       integer       :: ldc
       magma_tally4_devptr_t:: work
       integer       :: lwork
       magma_tally4_devptr_t:: td
       integer       :: nb
       integer       :: info
     end subroutine magma_tally4f_zunmqr_gpu

  end interface

contains
  
  subroutine magma_tally4f_zoff1d( ptrNew, ptrOld, inc, i)
    magma_tally4_devptr_t :: ptrNew
    magma_tally4_devptr_t :: ptrOld
    integer        :: inc, i
    
    ptrNew = ptrOld + (i-1) * inc * sizeof_complex_16
    
  end subroutine magma_tally4f_zoff1d
  
  subroutine magma_tally4f_zoff2d( ptrNew, ptrOld, lda, i, j)
    magma_tally4_devptr_t :: ptrNew
    magma_tally4_devptr_t :: ptrOld
    integer        :: lda, i, j
    
    ptrNew = ptrOld + ((j-1) * lda + (i-1)) * sizeof_complex_16
    
  end subroutine magma_tally4f_zoff2d
  
end module magma_tally4_zfortran
