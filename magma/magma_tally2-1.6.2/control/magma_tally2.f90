!
!   -- MAGMA_tally2 (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!

module magma_tally2

  use magma_tally2_param
  use magma_tally2_zfortran
  use magma_tally2_dfortran
  use magma_tally2_cfortran
  use magma_tally2_sfortran

  interface

  subroutine magma_tally2f_init( )
  end subroutine
  
  subroutine magma_tally2f_finalize(  )
  end subroutine
  
  end interface
  
  ! parameter constants from magma_tally2_types.h
  integer, parameter :: &
        Magma_tally2False         = 0,    &
        Magma_tally2True          = 1,    &
        Magma_tally2RowMajor      = 101,  &
        Magma_tally2ColMajor      = 102,  &
        Magma_tally2NoTrans       = 111,  &
        Magma_tally2Trans         = 112,  &
        Magma_tally2ConjTrans     = 113,  &
        Magma_tally2Upper         = 121,  &
        Magma_tally2Lower         = 122,  &
        Magma_tally2UpperLower    = 123,  &
        Magma_tally2Full          = 123,  &
        Magma_tally2NonUnit       = 131,  &
        Magma_tally2Unit          = 132,  &
        Magma_tally2Left          = 141,  &
        Magma_tally2Right         = 142,  &
        Magma_tally2BothSides     = 143,  &
        Magma_tally2OneNorm       = 171,  &
        Magma_tally2RealOneNorm   = 172,  &
        Magma_tally2TwoNorm       = 173,  &
        Magma_tally2FrobeniusNorm = 174,  &
        Magma_tally2InfNorm       = 175,  &
        Magma_tally2RealInfNorm   = 176,  &
        Magma_tally2MaxNorm       = 177,  &
        Magma_tally2RealMaxNorm   = 178,  &
        Magma_tally2DistUniform   = 201,  &
        Magma_tally2DistSymmetric = 202,  &
        Magma_tally2DistNormal    = 203,  &
        Magma_tally2HermGeev      = 241,  &
        Magma_tally2HermPoev      = 242,  &
        Magma_tally2NonsymPosv    = 243,  &
        Magma_tally2SymPosv       = 244,  &
        Magma_tally2NoPacking     = 291,  &
        Magma_tally2PackSubdiag   = 292,  &
        Magma_tally2PackSupdiag   = 293,  &
        Magma_tally2PackColumn    = 294,  &
        Magma_tally2PackRow       = 295,  &
        Magma_tally2PackLowerBand = 296,  &
        Magma_tally2PackUpeprBand = 297,  &
        Magma_tally2PackAll       = 298,  &
        Magma_tally2NoVec         = 301,  &
        Magma_tally2Vec           = 302,  &
        Magma_tally2IVec          = 303,  &
        Magma_tally2AllVec        = 304,  &
        Magma_tally2SomeVec       = 305,  &
        Magma_tally2OverwriteVec  = 306,  &
        Magma_tally2BacktransVec  = 307,  &
        Magma_tally2RangeAll      = 311,  &
        Magma_tally2RangeV        = 312,  &
        Magma_tally2RangeI        = 313,  &
        Magma_tally2Q             = 322,  &
        Magma_tally2P             = 323,  &
        Magma_tally2Forward       = 391,  &
        Magma_tally2Backward      = 392,  &
        Magma_tally2Columnwise    = 401,  &
        Magma_tally2Rowwise       = 402

end module magma_tally2
