!
!   -- MAGMA_minproduct (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!
!   @precisions normal z -> c d s
!

#define PRECISION_z

module magma_minproduct_zfortran

  use magma_minproduct_param, only: sizeof_complex_16

  implicit none

  !---- Fortran interfaces to MAGMA_minproduct subroutines ----
  interface

     subroutine magma_minproductf_zgetptr( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
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
     end subroutine magma_minproductf_zgetptr


     subroutine magma_minproductf_zgebrd( m, n, A, lda, d, e,tauq, taup, work, lwork, info)
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
     end subroutine magma_minproductf_zgebrd

     subroutine magma_minproductf_zgehrd2(n, ilo, ihi,A, lda, tau, work, lwork, info)
       integer       :: n
       integer       :: ilo
       integer       :: ihi
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_minproductf_zgehrd2

     subroutine magma_minproductf_zgehrd(n, ilo, ihi,A, lda, tau, work, lwork, d_T, info)
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
     end subroutine magma_minproductf_zgehrd

     subroutine magma_minproductf_zgelqf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_minproductf_zgelqf

     subroutine magma_minproductf_zgeqlf( m, n, A,    lda,   tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_minproductf_zgeqlf

     subroutine magma_minproductf_zgeqrf( m, n, A, lda, tau, work, lwork, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       complex*16    :: tau(*)
       complex*16    :: work(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_minproductf_zgeqrf

     subroutine magma_minproductf_zgesv(  n, nrhs, A, lda, ipiv, B, ldb, info)
       integer       :: n
       integer       :: nrhs
       complex*16    :: A
       integer       :: lda
       integer       :: ipiv(*)
       complex*16    :: B
       integer       :: ldb
       integer       :: info
     end subroutine magma_minproductf_zgesv

     subroutine magma_minproductf_zgetrf( m, n, A, lda, ipiv, info)
       integer       :: m
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magma_minproductf_zgetrf

     subroutine magma_minproductf_zposv(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character     :: uplo
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_minproductf_zposv
     
     subroutine magma_minproductf_zpotrf( uplo, n, A, lda, info)
       character          :: uplo
       integer       :: n
       complex*16    :: A(*)
       integer       :: lda
       integer       :: info
     end subroutine magma_minproductf_zpotrf

     subroutine magma_minproductf_zhetrd( uplo, n, A, lda, d, e, tau, work, lwork, info)
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
     end subroutine magma_minproductf_zhetrd

     subroutine magma_minproductf_zunmqr( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, info)
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
     end subroutine magma_minproductf_zunmqr

     subroutine magma_minproductf_zunmtr( side, uplo, trans, m, n, a, lda,tau,c,    ldc,work, lwork,info)
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
     end subroutine magma_minproductf_zunmtr
#if defined(PRECISION_z) || defined(PRECISION_c)

     subroutine magma_minproductf_zgeev( jobvl, jobvr, n, a, lda, w, vl, ldvl, vr, ldvr, work, lwork, rwork, info)
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
     end subroutine magma_minproductf_zgeev

     subroutine magma_minproductf_zgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info)
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
     end subroutine magma_minproductf_zgesvd

     subroutine magma_minproductf_zheevd( jobz, uplo, n, a, lda, w, work, lwork, rwork, lrwork, iwork, liwork, info)
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
     end subroutine magma_minproductf_zheevd

     subroutine magma_minproductf_zhegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, lrwork, iwork, liwork, info)
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
     end subroutine magma_minproductf_zhegvd

#else
     subroutine magma_minproductf_zgeev( jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info)
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
     end subroutine magma_minproductf_zgeev

     subroutine magma_minproductf_zgesvd( jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info)
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
     end subroutine magma_minproductf_zgesvd

     subroutine magma_minproductf_zheevd( jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
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
     end subroutine magma_minproductf_zheevd

     subroutine magma_minproductf_zhegvd( itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, iwork, liwork, info)
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
     end subroutine magma_minproductf_zhegvd
#endif

     subroutine magma_minproductf_zgels_gpu(  trans, m, n, nrhs, dA, ldda, dB, lddb, hwork, lwork, info)
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lwork
       integer       :: info
     end subroutine magma_minproductf_zgels_gpu

     subroutine magma_minproductf_zgeqrf_gpu( m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       magma_minproduct_devptr_t:: dT
       integer       :: info
     end subroutine magma_minproductf_zgeqrf_gpu

     subroutine magma_minproductf_zgeqrf2_gpu(m, n, dA, ldda, tau, info)
       integer       :: m
       integer       :: n
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       integer       :: info
     end subroutine magma_minproductf_zgeqrf2_gpu

     subroutine magma_minproductf_zgeqrf3_gpu(m, n, dA, ldda, tau, dT, info)
       integer       :: m
       integer       :: n
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau(*)
       magma_minproduct_devptr_t:: dT
       integer       :: info
     end subroutine magma_minproductf_zgeqrf3_gpu

     subroutine magma_minproductf_zgeqrs_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau
       magma_minproduct_devptr_t:: dT
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magma_minproductf_zgeqrs_gpu

     subroutine magma_minproductf_zgeqrs3_gpu( m, n, nrhs, dA, ldda, tau, dT, dB, lddb, hwork, lhwork, info)
       integer       :: m
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       complex*16    :: tau
       magma_minproduct_devptr_t:: dT
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       complex*16    :: hwork(*)
       integer       :: lhwork
       integer       :: info
     end subroutine magma_minproductf_zgeqrs3_gpu

     subroutine magma_minproductf_zgessm_gpu( storev, m, n, k, ib, ipiv, dL1, lddl1, dL,  lddl, dA,  ldda, info)
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       integer       :: ib
       integer       :: ipiv(*)
       magma_minproduct_devptr_t:: dL1
       integer       :: lddl1
       magma_minproduct_devptr_t:: dL
       integer       :: lddl
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magma_minproductf_zgessm_gpu

     subroutine magma_minproductf_zgesv_gpu(  n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_minproductf_zgesv_gpu

     subroutine magma_minproductf_zgetrf_gpu( m, n, dA, ldda, ipiv, info)
       integer       :: m
       integer       :: n
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       integer       :: info
     end subroutine magma_minproductf_zgetrf_gpu

     subroutine magma_minproductf_zgetrs_gpu( trans, n, nrhs, dA, ldda, ipiv, dB, lddb, info)
       character          :: trans
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       integer       :: ipiv(*)
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_minproductf_zgetrs_gpu

     subroutine magma_minproductf_zlabrd_gpu( m, n, nb, a, lda, da, ldda, d, e, tauq, taup, x, ldx, dx, lddx, y, ldy, dy, lddy)
       integer       :: m
       integer       :: n
       integer       :: nb
       complex*16    :: a(*)
       integer       :: lda
       magma_minproduct_devptr_t:: da
       integer       :: ldda
       double precision:: d(*)
       double precision:: e(*)
       complex*16    :: tauq(*)
       complex*16    :: taup(*)
       complex*16    :: x(*)
       integer       :: ldx
       magma_minproduct_devptr_t:: dx
       integer       :: lddx
       complex*16    :: y(*)
       integer       :: ldy
       magma_minproduct_devptr_t:: dy
       integer       :: lddy
     end subroutine magma_minproductf_zlabrd_gpu

     subroutine magma_minproductf_zlarfb_gpu( side, trans, direct, storev, m, n, k, dv, ldv, dt, ldt, dc, ldc, dowrk, ldwork)
       character          :: side
       character          :: trans
       character          :: direct
       character          :: storev
       integer       :: m
       integer       :: n
       integer       :: k
       magma_minproduct_devptr_t:: dv
       integer       :: ldv
       magma_minproduct_devptr_t:: dt
       integer       :: ldt
       magma_minproduct_devptr_t:: dc
       integer       :: ldc
       magma_minproduct_devptr_t:: dowrk
       integer       :: ldwork
     end subroutine magma_minproductf_zlarfb_gpu

     subroutine magma_minproductf_zposv_gpu(  uplo, n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_minproductf_zposv_gpu

     subroutine magma_minproductf_zpotrf_gpu( uplo, n, dA, ldda, info)
       character          :: uplo
       integer       :: n
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       integer       :: info
     end subroutine magma_minproductf_zpotrf_gpu

     subroutine magma_minproductf_zpotrs_gpu( uplo,  n, nrhs, dA, ldda, dB, lddb, info)
       character          :: uplo
       integer       :: n
       integer       :: nrhs
       magma_minproduct_devptr_t:: dA
       integer       :: ldda
       magma_minproduct_devptr_t:: dB
       integer       :: lddb
       integer       :: info
     end subroutine magma_minproductf_zpotrs_gpu

     subroutine magma_minproductf_zssssm_gpu( storev, m1, n1, m2, n2, k, ib, dA1, ldda1, dA2, ldda2, dL1, lddl1, dL2, lddl2, IPIV, info)
       character          :: storev
       integer       :: m1
       integer       :: n1
       integer       :: m2
       integer       :: n2
       integer       :: k
       integer       :: ib
       magma_minproduct_devptr_t:: dA1
       integer       :: ldda1
       magma_minproduct_devptr_t:: dA2
       integer       :: ldda2
       magma_minproduct_devptr_t:: dL1
       integer       :: lddl1
       magma_minproduct_devptr_t:: dL2
       integer       :: lddl2
       integer       :: IPIV(*)
       integer       :: info
     end subroutine magma_minproductf_zssssm_gpu

     subroutine magma_minproductf_zungqr_gpu( m, n, k, da, ldda, tau, dwork, nb, info)
       integer       :: m
       integer       :: n
       integer       :: k
       magma_minproduct_devptr_t:: da
       integer       :: ldda
       complex*16    :: tau(*)
       magma_minproduct_devptr_t:: dwork
       integer       :: nb
       integer       :: info
     end subroutine magma_minproductf_zungqr_gpu

     subroutine magma_minproductf_zunmqr_gpu( side, trans, m, n, k, a, lda, tau, c, ldc, work, lwork, td, nb, info)
       character          :: side
       character          :: trans
       integer       :: m
       integer       :: n
       integer       :: k
       magma_minproduct_devptr_t:: a
       integer       :: lda
       complex*16    :: tau(*)
       magma_minproduct_devptr_t:: c
       integer       :: ldc
       magma_minproduct_devptr_t:: work
       integer       :: lwork
       magma_minproduct_devptr_t:: td
       integer       :: nb
       integer       :: info
     end subroutine magma_minproductf_zunmqr_gpu

  end interface

contains
  
  subroutine magma_minproductf_zoff1d( ptrNew, ptrOld, inc, i)
    magma_minproduct_devptr_t :: ptrNew
    magma_minproduct_devptr_t :: ptrOld
    integer        :: inc, i
    
    ptrNew = ptrOld + (i-1) * inc * sizeof_complex_16
    
  end subroutine magma_minproductf_zoff1d
  
  subroutine magma_minproductf_zoff2d( ptrNew, ptrOld, lda, i, j)
    magma_minproduct_devptr_t :: ptrNew
    magma_minproduct_devptr_t :: ptrOld
    integer        :: lda, i, j
    
    ptrNew = ptrOld + ((j-1) * lda + (i-1)) * sizeof_complex_16
    
  end subroutine magma_minproductf_zoff2d
  
end module magma_minproduct_zfortran
