!
!   -- MAGMA_tally3 (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!

module magma_tally3

  use magma_tally3_param
  use magma_tally3_zfortran
  use magma_tally3_dfortran
  use magma_tally3_cfortran
  use magma_tally3_sfortran

  interface

  subroutine magma_tally3f_init( )
  end subroutine
  
  subroutine magma_tally3f_finalize(  )
  end subroutine
  
  end interface
  
  ! parameter constants from magma_tally3_types.h
  integer, parameter :: &
        Magma_tally3False         = 0,    &
        Magma_tally3True          = 1,    &
        Magma_tally3RowMajor      = 101,  &
        Magma_tally3ColMajor      = 102,  &
        Magma_tally3NoTrans       = 111,  &
        Magma_tally3Trans         = 112,  &
        Magma_tally3ConjTrans     = 113,  &
        Magma_tally3Upper         = 121,  &
        Magma_tally3Lower         = 122,  &
        Magma_tally3UpperLower    = 123,  &
        Magma_tally3Full          = 123,  &
        Magma_tally3NonUnit       = 131,  &
        Magma_tally3Unit          = 132,  &
        Magma_tally3Left          = 141,  &
        Magma_tally3Right         = 142,  &
        Magma_tally3BothSides     = 143,  &
        Magma_tally3OneNorm       = 171,  &
        Magma_tally3RealOneNorm   = 172,  &
        Magma_tally3TwoNorm       = 173,  &
        Magma_tally3FrobeniusNorm = 174,  &
        Magma_tally3InfNorm       = 175,  &
        Magma_tally3RealInfNorm   = 176,  &
        Magma_tally3MaxNorm       = 177,  &
        Magma_tally3RealMaxNorm   = 178,  &
        Magma_tally3DistUniform   = 201,  &
        Magma_tally3DistSymmetric = 202,  &
        Magma_tally3DistNormal    = 203,  &
        Magma_tally3HermGeev      = 241,  &
        Magma_tally3HermPoev      = 242,  &
        Magma_tally3NonsymPosv    = 243,  &
        Magma_tally3SymPosv       = 244,  &
        Magma_tally3NoPacking     = 291,  &
        Magma_tally3PackSubdiag   = 292,  &
        Magma_tally3PackSupdiag   = 293,  &
        Magma_tally3PackColumn    = 294,  &
        Magma_tally3PackRow       = 295,  &
        Magma_tally3PackLowerBand = 296,  &
        Magma_tally3PackUpeprBand = 297,  &
        Magma_tally3PackAll       = 298,  &
        Magma_tally3NoVec         = 301,  &
        Magma_tally3Vec           = 302,  &
        Magma_tally3IVec          = 303,  &
        Magma_tally3AllVec        = 304,  &
        Magma_tally3SomeVec       = 305,  &
        Magma_tally3OverwriteVec  = 306,  &
        Magma_tally3BacktransVec  = 307,  &
        Magma_tally3RangeAll      = 311,  &
        Magma_tally3RangeV        = 312,  &
        Magma_tally3RangeI        = 313,  &
        Magma_tally3Q             = 322,  &
        Magma_tally3P             = 323,  &
        Magma_tally3Forward       = 391,  &
        Magma_tally3Backward      = 392,  &
        Magma_tally3Columnwise    = 401,  &
        Magma_tally3Rowwise       = 402

end module magma_tally3
