!
!   -- MAGMA_tally4 (version 1.6.1) --
!      Univ. of Tennessee, Knoxville
!      Univ. of California, Berkeley
!      Univ. of Colorado, Denver
!      @date January 2015
!

module magma_tally4

  use magma_tally4_param
  use magma_tally4_zfortran
  use magma_tally4_dfortran
  use magma_tally4_cfortran
  use magma_tally4_sfortran

  interface

  subroutine magma_tally4f_init( )
  end subroutine
  
  subroutine magma_tally4f_finalize(  )
  end subroutine
  
  end interface
  
  ! parameter constants from magma_tally4_types.h
  integer, parameter :: &
        Magma_tally4False         = 0,    &
        Magma_tally4True          = 1,    &
        Magma_tally4RowMajor      = 101,  &
        Magma_tally4ColMajor      = 102,  &
        Magma_tally4NoTrans       = 111,  &
        Magma_tally4Trans         = 112,  &
        Magma_tally4ConjTrans     = 113,  &
        Magma_tally4Upper         = 121,  &
        Magma_tally4Lower         = 122,  &
        Magma_tally4UpperLower    = 123,  &
        Magma_tally4Full          = 123,  &
        Magma_tally4NonUnit       = 131,  &
        Magma_tally4Unit          = 132,  &
        Magma_tally4Left          = 141,  &
        Magma_tally4Right         = 142,  &
        Magma_tally4BothSides     = 143,  &
        Magma_tally4OneNorm       = 171,  &
        Magma_tally4RealOneNorm   = 172,  &
        Magma_tally4TwoNorm       = 173,  &
        Magma_tally4FrobeniusNorm = 174,  &
        Magma_tally4InfNorm       = 175,  &
        Magma_tally4RealInfNorm   = 176,  &
        Magma_tally4MaxNorm       = 177,  &
        Magma_tally4RealMaxNorm   = 178,  &
        Magma_tally4DistUniform   = 201,  &
        Magma_tally4DistSymmetric = 202,  &
        Magma_tally4DistNormal    = 203,  &
        Magma_tally4HermGeev      = 241,  &
        Magma_tally4HermPoev      = 242,  &
        Magma_tally4NonsymPosv    = 243,  &
        Magma_tally4SymPosv       = 244,  &
        Magma_tally4NoPacking     = 291,  &
        Magma_tally4PackSubdiag   = 292,  &
        Magma_tally4PackSupdiag   = 293,  &
        Magma_tally4PackColumn    = 294,  &
        Magma_tally4PackRow       = 295,  &
        Magma_tally4PackLowerBand = 296,  &
        Magma_tally4PackUpeprBand = 297,  &
        Magma_tally4PackAll       = 298,  &
        Magma_tally4NoVec         = 301,  &
        Magma_tally4Vec           = 302,  &
        Magma_tally4IVec          = 303,  &
        Magma_tally4AllVec        = 304,  &
        Magma_tally4SomeVec       = 305,  &
        Magma_tally4OverwriteVec  = 306,  &
        Magma_tally4BacktransVec  = 307,  &
        Magma_tally4RangeAll      = 311,  &
        Magma_tally4RangeV        = 312,  &
        Magma_tally4RangeI        = 313,  &
        Magma_tally4Q             = 322,  &
        Magma_tally4P             = 323,  &
        Magma_tally4Forward       = 391,  &
        Magma_tally4Backward      = 392,  &
        Magma_tally4Columnwise    = 401,  &
        Magma_tally4Rowwise       = 402

end module magma_tally4
