!
!   -- MAGMA_minproduct (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!

module magma_minproduct

  use magma_minproduct_param
  use magma_minproduct_zfortran
  use magma_minproduct_dfortran
  use magma_minproduct_cfortran
  use magma_minproduct_sfortran

  interface

  subroutine magma_minproductf_init( )
  end subroutine
  
  subroutine magma_minproductf_finalize(  )
  end subroutine
  
  end interface
  
  ! parameter constants from magma_minproduct_types.h
  integer, parameter :: &
        Magma_minproductFalse         = 0,    &
        Magma_minproductTrue          = 1,    &
        Magma_minproductRowMajor      = 101,  &
        Magma_minproductColMajor      = 102,  &
        Magma_minproductNoTrans       = 111,  &
        Magma_minproductTrans         = 112,  &
        Magma_minproductConjTrans     = 113,  &
        Magma_minproductUpper         = 121,  &
        Magma_minproductLower         = 122,  &
        Magma_minproductUpperLower    = 123,  &
        Magma_minproductFull          = 123,  &
        Magma_minproductNonUnit       = 131,  &
        Magma_minproductUnit          = 132,  &
        Magma_minproductLeft          = 141,  &
        Magma_minproductRight         = 142,  &
        Magma_minproductBothSides     = 143,  &
        Magma_minproductOneNorm       = 171,  &
        Magma_minproductRealOneNorm   = 172,  &
        Magma_minproductTwoNorm       = 173,  &
        Magma_minproductFrobeniusNorm = 174,  &
        Magma_minproductInfNorm       = 175,  &
        Magma_minproductRealInfNorm   = 176,  &
        Magma_minproductMaxNorm       = 177,  &
        Magma_minproductRealMaxNorm   = 178,  &
        Magma_minproductDistUniform   = 201,  &
        Magma_minproductDistSymmetric = 202,  &
        Magma_minproductDistNormal    = 203,  &
        Magma_minproductHermGeev      = 241,  &
        Magma_minproductHermPoev      = 242,  &
        Magma_minproductNonsymPosv    = 243,  &
        Magma_minproductSymPosv       = 244,  &
        Magma_minproductNoPacking     = 291,  &
        Magma_minproductPackSubdiag   = 292,  &
        Magma_minproductPackSupdiag   = 293,  &
        Magma_minproductPackColumn    = 294,  &
        Magma_minproductPackRow       = 295,  &
        Magma_minproductPackLowerBand = 296,  &
        Magma_minproductPackUpeprBand = 297,  &
        Magma_minproductPackAll       = 298,  &
        Magma_minproductNoVec         = 301,  &
        Magma_minproductVec           = 302,  &
        Magma_minproductIVec          = 303,  &
        Magma_minproductAllVec        = 304,  &
        Magma_minproductSomeVec       = 305,  &
        Magma_minproductOverwriteVec  = 306,  &
        Magma_minproductBacktransVec  = 307,  &
        Magma_minproductRangeAll      = 311,  &
        Magma_minproductRangeV        = 312,  &
        Magma_minproductRangeI        = 313,  &
        Magma_minproductQ             = 322,  &
        Magma_minproductP             = 323,  &
        Magma_minproductForward       = 391,  &
        Magma_minproductBackward      = 392,  &
        Magma_minproductColumnwise    = 401,  &
        Magma_minproductRowwise       = 402

end module magma_minproduct
