/*
    -- MAGMA_tally4 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/
#include <stdio.h>

#include "testings.h"
#include "magma_tally4.h"

int gStatus;

void check_( bool flag, const char* msg, int line )
{
    if ( ! flag ) {
        gStatus += 1;
        printf( "line %d: %s failed\n", line, msg );
    }
}

#define check( flag ) check_( flag, #flag, __LINE__ )


/* ////////////////////////////////////////////////////////////////////////////
   -- Testing lapack_xxxxx_const and related
*/
int main( int argc, char** argv )
{
    gStatus = 0;
    int s;

    // ------------------------------------------------------------
    s = gStatus;
    check( lapack_bool_const_tally4(   Magma_tally4False         )[0] == 'N' );
    check( lapack_bool_const_tally4(   Magma_tally4True          )[0] == 'Y' );

    check( lapack_order_const_tally4(  Magma_tally4RowMajor      )[0] == 'R' );
    check( lapack_order_const_tally4(  Magma_tally4ColMajor      )[0] == 'C' );

    check( lapack_trans_const_tally4(  Magma_tally4NoTrans       )[0] == 'N' );
    check( lapack_trans_const_tally4(  Magma_tally4Trans         )[0] == 'T' );
    check( lapack_trans_const_tally4(  Magma_tally4ConjTrans     )[0] == 'C' );

    check( lapack_uplo_const_tally4(   Magma_tally4Upper         )[0] == 'U' );
    check( lapack_uplo_const_tally4(   Magma_tally4Lower         )[0] == 'L' );
    check( lapack_uplo_const_tally4(   Magma_tally4Full          )[0] == 'G' );

    check( lapack_diag_const_tally4(   Magma_tally4NonUnit       )[0] == 'N' );
    check( lapack_diag_const_tally4(   Magma_tally4Unit          )[0] == 'U' );

    check( lapack_side_const_tally4(   Magma_tally4Left          )[0] == 'L' );
    check( lapack_side_const_tally4(   Magma_tally4Right         )[0] == 'R' );
    check( lapack_side_const_tally4(   Magma_tally4BothSides     )[0] == 'B' );

    check( lapack_norm_const_tally4(   Magma_tally4OneNorm       )[0] == '1' );
    check( lapack_norm_const_tally4(   Magma_tally4TwoNorm       )[0] == '2' );
    check( lapack_norm_const_tally4(   Magma_tally4FrobeniusNorm )[0] == 'F' );
    check( lapack_norm_const_tally4(   Magma_tally4InfNorm       )[0] == 'I' );
    check( lapack_norm_const_tally4(   Magma_tally4MaxNorm       )[0] == 'M' );

    check( lapack_dist_const_tally4(   Magma_tally4DistUniform   )[0] == 'U' );
    check( lapack_dist_const_tally4(   Magma_tally4DistSymmetric )[0] == 'S' );
    check( lapack_dist_const_tally4(   Magma_tally4DistNormal    )[0] == 'N' );

    check( lapack_sym_const_tally4(    Magma_tally4HermGeev      )[0] == 'H' );
    check( lapack_sym_const_tally4(    Magma_tally4HermPoev      )[0] == 'P' );
    check( lapack_sym_const_tally4(    Magma_tally4NonsymPosv    )[0] == 'N' );
    check( lapack_sym_const_tally4(    Magma_tally4SymPosv       )[0] == 'S' );

    check( lapack_pack_const_tally4(   Magma_tally4NoPacking     )[0] == 'N' );
    check( lapack_pack_const_tally4(   Magma_tally4PackSubdiag   )[0] == 'U' );
    check( lapack_pack_const_tally4(   Magma_tally4PackSupdiag   )[0] == 'L' );
    check( lapack_pack_const_tally4(   Magma_tally4PackColumn    )[0] == 'C' );
    check( lapack_pack_const_tally4(   Magma_tally4PackRow       )[0] == 'R' );
    check( lapack_pack_const_tally4(   Magma_tally4PackLowerBand )[0] == 'B' );
    check( lapack_pack_const_tally4(   Magma_tally4PackUpeprBand )[0] == 'Q' );
    check( lapack_pack_const_tally4(   Magma_tally4PackAll       )[0] == 'Z' );

    check( lapack_vec_const_tally4(    Magma_tally4NoVec         )[0] == 'N' );
    check( lapack_vec_const_tally4(    Magma_tally4Vec           )[0] == 'V' );
    check( lapack_vec_const_tally4(    Magma_tally4IVec          )[0] == 'I' );
    check( lapack_vec_const_tally4(    Magma_tally4AllVec        )[0] == 'A' );
    check( lapack_vec_const_tally4(    Magma_tally4SomeVec       )[0] == 'S' );
    check( lapack_vec_const_tally4(    Magma_tally4OverwriteVec  )[0] == 'O' );

    check( lapack_range_const_tally4(  Magma_tally4RangeAll      )[0] == 'A' );
    check( lapack_range_const_tally4(  Magma_tally4RangeV        )[0] == 'V' );
    check( lapack_range_const_tally4(  Magma_tally4RangeI        )[0] == 'I' );

    check( lapack_vect_const_tally4(   Magma_tally4Q             )[0] == 'Q' );
    check( lapack_vect_const_tally4(   Magma_tally4P             )[0] == 'P' );

    check( lapack_direct_const_tally4( Magma_tally4Forward       )[0] == 'F' );
    check( lapack_direct_const_tally4( Magma_tally4Backward      )[0] == 'B' );

    check( lapack_storev_const_tally4( Magma_tally4Columnwise    )[0] == 'C' );
    check( lapack_storev_const_tally4( Magma_tally4Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally4  -> lapack_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_bool_const_tally4(   Magma_tally4False         ) == 'N' );
    check( lapacke_bool_const_tally4(   Magma_tally4True          ) == 'Y' );

    check( lapacke_order_const_tally4(  Magma_tally4RowMajor      ) == 'R' );
    check( lapacke_order_const_tally4(  Magma_tally4ColMajor      ) == 'C' );

    check( lapacke_trans_const_tally4(  Magma_tally4NoTrans       ) == 'N' );
    check( lapacke_trans_const_tally4(  Magma_tally4Trans         ) == 'T' );
    check( lapacke_trans_const_tally4(  Magma_tally4ConjTrans     ) == 'C' );

    check( lapacke_uplo_const_tally4(   Magma_tally4Upper         ) == 'U' );
    check( lapacke_uplo_const_tally4(   Magma_tally4Lower         ) == 'L' );
    check( lapacke_uplo_const_tally4(   Magma_tally4Full          ) == 'G' );

    check( lapacke_diag_const_tally4(   Magma_tally4NonUnit       ) == 'N' );
    check( lapacke_diag_const_tally4(   Magma_tally4Unit          ) == 'U' );

    check( lapacke_side_const_tally4(   Magma_tally4Left          ) == 'L' );
    check( lapacke_side_const_tally4(   Magma_tally4Right         ) == 'R' );
    check( lapacke_side_const_tally4(   Magma_tally4BothSides     ) == 'B' );

    check( lapacke_norm_const_tally4(   Magma_tally4OneNorm       ) == '1' );
    check( lapacke_norm_const_tally4(   Magma_tally4TwoNorm       ) == '2' );
    check( lapacke_norm_const_tally4(   Magma_tally4FrobeniusNorm ) == 'F' );
    check( lapacke_norm_const_tally4(   Magma_tally4InfNorm       ) == 'I' );
    check( lapacke_norm_const_tally4(   Magma_tally4MaxNorm       ) == 'M' );

    check( lapacke_dist_const_tally4(   Magma_tally4DistUniform   ) == 'U' );
    check( lapacke_dist_const_tally4(   Magma_tally4DistSymmetric ) == 'S' );
    check( lapacke_dist_const_tally4(   Magma_tally4DistNormal    ) == 'N' );

    check( lapacke_sym_const_tally4(    Magma_tally4HermGeev      ) == 'H' );
    check( lapacke_sym_const_tally4(    Magma_tally4HermPoev      ) == 'P' );
    check( lapacke_sym_const_tally4(    Magma_tally4NonsymPosv    ) == 'N' );
    check( lapacke_sym_const_tally4(    Magma_tally4SymPosv       ) == 'S' );

    check( lapacke_pack_const_tally4(   Magma_tally4NoPacking     ) == 'N' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackSubdiag   ) == 'U' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackSupdiag   ) == 'L' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackColumn    ) == 'C' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackRow       ) == 'R' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackLowerBand ) == 'B' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackUpeprBand ) == 'Q' );
    check( lapacke_pack_const_tally4(   Magma_tally4PackAll       ) == 'Z' );

    check( lapacke_vec_const_tally4(    Magma_tally4NoVec         ) == 'N' );
    check( lapacke_vec_const_tally4(    Magma_tally4Vec           ) == 'V' );
    check( lapacke_vec_const_tally4(    Magma_tally4IVec          ) == 'I' );
    check( lapacke_vec_const_tally4(    Magma_tally4AllVec        ) == 'A' );
    check( lapacke_vec_const_tally4(    Magma_tally4SomeVec       ) == 'S' );
    check( lapacke_vec_const_tally4(    Magma_tally4OverwriteVec  ) == 'O' );

    check( lapacke_range_const_tally4(  Magma_tally4RangeAll      ) == 'A' );
    check( lapacke_range_const_tally4(  Magma_tally4RangeV        ) == 'V' );
    check( lapacke_range_const_tally4(  Magma_tally4RangeI        ) == 'I' );

    check( lapacke_vect_const_tally4(   Magma_tally4Q             ) == 'Q' );
    check( lapacke_vect_const_tally4(   Magma_tally4P             ) == 'P' );

    check( lapacke_direct_const_tally4( Magma_tally4Forward       ) == 'F' );
    check( lapacke_direct_const_tally4( Magma_tally4Backward      ) == 'B' );

    check( lapacke_storev_const_tally4( Magma_tally4Columnwise    ) == 'C' );
    check( lapacke_storev_const_tally4( Magma_tally4Rowwise       ) == 'R' );
    printf( "MAGMA_tally4  -> lapacke_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapack_const_tally4( Magma_tally4False         )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4True          )[0] == 'Y' );

    check( lapack_const_tally4( Magma_tally4RowMajor      )[0] == 'R' );
    check( lapack_const_tally4( Magma_tally4ColMajor      )[0] == 'C' );

    check( lapack_const_tally4( Magma_tally4NoTrans       )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4Trans         )[0] == 'T' );
    check( lapack_const_tally4( Magma_tally4ConjTrans     )[0] == 'C' );

    check( lapack_const_tally4( Magma_tally4Upper         )[0] == 'U' );
    check( lapack_const_tally4( Magma_tally4Lower         )[0] == 'L' );
    check( lapack_const_tally4( Magma_tally4Full          )[0] == 'G' );

    check( lapack_const_tally4( Magma_tally4NonUnit       )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4Unit          )[0] == 'U' );

    check( lapack_const_tally4( Magma_tally4Left          )[0] == 'L' );
    check( lapack_const_tally4( Magma_tally4Right         )[0] == 'R' );
    check( lapack_const_tally4( Magma_tally4BothSides     )[0] == 'B' );

    check( lapack_const_tally4( Magma_tally4OneNorm       )[0] == '1' );
    check( lapack_const_tally4( Magma_tally4TwoNorm       )[0] == '2' );
    check( lapack_const_tally4( Magma_tally4FrobeniusNorm )[0] == 'F' );
    check( lapack_const_tally4( Magma_tally4InfNorm       )[0] == 'I' );
    check( lapack_const_tally4( Magma_tally4MaxNorm       )[0] == 'M' );

    check( lapack_const_tally4( Magma_tally4DistUniform   )[0] == 'U' );
    check( lapack_const_tally4( Magma_tally4DistSymmetric )[0] == 'S' );
    check( lapack_const_tally4( Magma_tally4DistNormal    )[0] == 'N' );

    check( lapack_const_tally4( Magma_tally4HermGeev      )[0] == 'H' );
    check( lapack_const_tally4( Magma_tally4HermPoev      )[0] == 'P' );
    check( lapack_const_tally4( Magma_tally4NonsymPosv    )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4SymPosv       )[0] == 'S' );

    check( lapack_const_tally4( Magma_tally4NoPacking     )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4PackSubdiag   )[0] == 'U' );
    check( lapack_const_tally4( Magma_tally4PackSupdiag   )[0] == 'L' );
    check( lapack_const_tally4( Magma_tally4PackColumn    )[0] == 'C' );
    check( lapack_const_tally4( Magma_tally4PackRow       )[0] == 'R' );
    check( lapack_const_tally4( Magma_tally4PackLowerBand )[0] == 'B' );
    check( lapack_const_tally4( Magma_tally4PackUpeprBand )[0] == 'Q' );
    check( lapack_const_tally4( Magma_tally4PackAll       )[0] == 'Z' );

    check( lapack_const_tally4( Magma_tally4NoVec         )[0] == 'N' );
    check( lapack_const_tally4( Magma_tally4Vec           )[0] == 'V' );
    check( lapack_const_tally4( Magma_tally4IVec          )[0] == 'I' );
    check( lapack_const_tally4( Magma_tally4AllVec        )[0] == 'A' );
    check( lapack_const_tally4( Magma_tally4SomeVec       )[0] == 'S' );
    check( lapack_const_tally4( Magma_tally4OverwriteVec  )[0] == 'O' );

    check( lapack_const_tally4( Magma_tally4RangeAll      )[0] == 'A' );
    check( lapack_const_tally4( Magma_tally4RangeV        )[0] == 'V' );
    check( lapack_const_tally4( Magma_tally4RangeI        )[0] == 'I' );

    check( lapack_const_tally4( Magma_tally4Q             )[0] == 'Q' );
    check( lapack_const_tally4( Magma_tally4P             )[0] == 'P' );

    check( lapack_const_tally4( Magma_tally4Forward       )[0] == 'F' );
    check( lapack_const_tally4( Magma_tally4Backward      )[0] == 'B' );

    check( lapack_const_tally4( Magma_tally4Columnwise    )[0] == 'C' );
    check( lapack_const_tally4( Magma_tally4Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally4  -> lapack_const_tally4          %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_const_tally4( Magma_tally4False         ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4True          ) == 'Y' );

    check( lapacke_const_tally4( Magma_tally4RowMajor      ) == 'R' );
    check( lapacke_const_tally4( Magma_tally4ColMajor      ) == 'C' );

    check( lapacke_const_tally4( Magma_tally4NoTrans       ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4Trans         ) == 'T' );
    check( lapacke_const_tally4( Magma_tally4ConjTrans     ) == 'C' );

    check( lapacke_const_tally4( Magma_tally4Upper         ) == 'U' );
    check( lapacke_const_tally4( Magma_tally4Lower         ) == 'L' );
    check( lapacke_const_tally4( Magma_tally4Full          ) == 'G' );

    check( lapacke_const_tally4( Magma_tally4NonUnit       ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4Unit          ) == 'U' );

    check( lapacke_const_tally4( Magma_tally4Left          ) == 'L' );
    check( lapacke_const_tally4( Magma_tally4Right         ) == 'R' );
    check( lapacke_const_tally4( Magma_tally4BothSides     ) == 'B' );

    check( lapacke_const_tally4( Magma_tally4OneNorm       ) == '1' );
    check( lapacke_const_tally4( Magma_tally4TwoNorm       ) == '2' );
    check( lapacke_const_tally4( Magma_tally4FrobeniusNorm ) == 'F' );
    check( lapacke_const_tally4( Magma_tally4InfNorm       ) == 'I' );
    check( lapacke_const_tally4( Magma_tally4MaxNorm       ) == 'M' );

    check( lapacke_const_tally4( Magma_tally4DistUniform   ) == 'U' );
    check( lapacke_const_tally4( Magma_tally4DistSymmetric ) == 'S' );
    check( lapacke_const_tally4( Magma_tally4DistNormal    ) == 'N' );

    check( lapacke_const_tally4( Magma_tally4HermGeev      ) == 'H' );
    check( lapacke_const_tally4( Magma_tally4HermPoev      ) == 'P' );
    check( lapacke_const_tally4( Magma_tally4NonsymPosv    ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4SymPosv       ) == 'S' );

    check( lapacke_const_tally4( Magma_tally4NoPacking     ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4PackSubdiag   ) == 'U' );
    check( lapacke_const_tally4( Magma_tally4PackSupdiag   ) == 'L' );
    check( lapacke_const_tally4( Magma_tally4PackColumn    ) == 'C' );
    check( lapacke_const_tally4( Magma_tally4PackRow       ) == 'R' );
    check( lapacke_const_tally4( Magma_tally4PackLowerBand ) == 'B' );
    check( lapacke_const_tally4( Magma_tally4PackUpeprBand ) == 'Q' );
    check( lapacke_const_tally4( Magma_tally4PackAll       ) == 'Z' );

    check( lapacke_const_tally4( Magma_tally4NoVec         ) == 'N' );
    check( lapacke_const_tally4( Magma_tally4Vec           ) == 'V' );
    check( lapacke_const_tally4( Magma_tally4IVec          ) == 'I' );
    check( lapacke_const_tally4( Magma_tally4AllVec        ) == 'A' );
    check( lapacke_const_tally4( Magma_tally4SomeVec       ) == 'S' );
    check( lapacke_const_tally4( Magma_tally4OverwriteVec  ) == 'O' );

    check( lapacke_const_tally4( Magma_tally4RangeAll      ) == 'A' );
    check( lapacke_const_tally4( Magma_tally4RangeV        ) == 'V' );
    check( lapacke_const_tally4( Magma_tally4RangeI        ) == 'I' );

    check( lapacke_const_tally4( Magma_tally4Q             ) == 'Q' );
    check( lapacke_const_tally4( Magma_tally4P             ) == 'P' );

    check( lapacke_const_tally4( Magma_tally4Forward       ) == 'F' );
    check( lapacke_const_tally4( Magma_tally4Backward      ) == 'B' );

    check( lapacke_const_tally4( Magma_tally4Columnwise    ) == 'C' );
    check( lapacke_const_tally4( Magma_tally4Rowwise       ) == 'R' );
    printf( "MAGMA_tally4  -> lapacke_const_tally4         %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( magma_tally4_bool_const('N') == Magma_tally4False );
    check( magma_tally4_bool_const('n') == Magma_tally4False );
    check( magma_tally4_bool_const('Y') == Magma_tally4True  );
    check( magma_tally4_bool_const('y') == Magma_tally4True  );

    check( magma_tally4_order_const( 'R' ) == Magma_tally4RowMajor  );
    check( magma_tally4_order_const( 'r' ) == Magma_tally4RowMajor  );
    check( magma_tally4_order_const( 'C' ) == Magma_tally4ColMajor  );
    check( magma_tally4_order_const( 'c' ) == Magma_tally4ColMajor  );

    check( magma_tally4_trans_const( 'N' ) == Magma_tally4NoTrans   );
    check( magma_tally4_trans_const( 'n' ) == Magma_tally4NoTrans   );
    check( magma_tally4_trans_const( 'T' ) == Magma_tally4Trans     );
    check( magma_tally4_trans_const( 't' ) == Magma_tally4Trans     );
    check( magma_tally4_trans_const( 'C' ) == Magma_tally4ConjTrans );
    check( magma_tally4_trans_const( 'c' ) == Magma_tally4ConjTrans );

    check( magma_tally4_uplo_const( 'U' ) == Magma_tally4Upper      );
    check( magma_tally4_uplo_const( 'u' ) == Magma_tally4Upper      );
    check( magma_tally4_uplo_const( 'L' ) == Magma_tally4Lower      );
    check( magma_tally4_uplo_const( 'l' ) == Magma_tally4Lower      );
    check( magma_tally4_uplo_const( 'A' ) == Magma_tally4Full       );  // anything else
    check( magma_tally4_uplo_const( 'a' ) == Magma_tally4Full       );
    check( magma_tally4_uplo_const( 'G' ) == Magma_tally4Full       );
    check( magma_tally4_uplo_const( 'g' ) == Magma_tally4Full       );
    check( magma_tally4_uplo_const( 'F' ) == Magma_tally4Full       );
    check( magma_tally4_uplo_const( 'f' ) == Magma_tally4Full       );

    check( magma_tally4_diag_const( 'N' ) == Magma_tally4NonUnit    );
    check( magma_tally4_diag_const( 'n' ) == Magma_tally4NonUnit    );
    check( magma_tally4_diag_const( 'U' ) == Magma_tally4Unit       );
    check( magma_tally4_diag_const( 'u' ) == Magma_tally4Unit       );

    check( magma_tally4_side_const( 'L' ) == Magma_tally4Left       );
    check( magma_tally4_side_const( 'l' ) == Magma_tally4Left       );
    check( magma_tally4_side_const( 'R' ) == Magma_tally4Right      );
    check( magma_tally4_side_const( 'r' ) == Magma_tally4Right      );

    check( magma_tally4_norm_const( 'O' ) == Magma_tally4OneNorm       );
    check( magma_tally4_norm_const( 'o' ) == Magma_tally4OneNorm       );
    check( magma_tally4_norm_const( '1' ) == Magma_tally4OneNorm       );
    check( magma_tally4_norm_const( '2' ) == Magma_tally4TwoNorm       );
    check( magma_tally4_norm_const( 'F' ) == Magma_tally4FrobeniusNorm );
    check( magma_tally4_norm_const( 'f' ) == Magma_tally4FrobeniusNorm );
    check( magma_tally4_norm_const( 'E' ) == Magma_tally4FrobeniusNorm );
    check( magma_tally4_norm_const( 'e' ) == Magma_tally4FrobeniusNorm );
    check( magma_tally4_norm_const( 'I' ) == Magma_tally4InfNorm       );
    check( magma_tally4_norm_const( 'i' ) == Magma_tally4InfNorm       );
    check( magma_tally4_norm_const( 'M' ) == Magma_tally4MaxNorm       );
    check( magma_tally4_norm_const( 'm' ) == Magma_tally4MaxNorm       );

    check( magma_tally4_dist_const( 'U' ) == Magma_tally4DistUniform   );
    check( magma_tally4_dist_const( 'u' ) == Magma_tally4DistUniform   );
    check( magma_tally4_dist_const( 'S' ) == Magma_tally4DistSymmetric );
    check( magma_tally4_dist_const( 's' ) == Magma_tally4DistSymmetric );
    check( magma_tally4_dist_const( 'N' ) == Magma_tally4DistNormal    );
    check( magma_tally4_dist_const( 'n' ) == Magma_tally4DistNormal    );

    //check( magma_tally4_xxxx_const( 'H' ) == Magma_tally4HermGeev      );
    //check( magma_tally4_xxxx_const( 'P' ) == Magma_tally4HermPoev      );
    //check( magma_tally4_xxxx_const( 'N' ) == Magma_tally4NonsymPosv    );
    //check( magma_tally4_xxxx_const( 'S' ) == Magma_tally4SymPosv       );

    check( magma_tally4_pack_const( 'N' ) == Magma_tally4NoPacking     );
    check( magma_tally4_pack_const( 'n' ) == Magma_tally4NoPacking     );
    check( magma_tally4_pack_const( 'U' ) == Magma_tally4PackSubdiag   );
    check( magma_tally4_pack_const( 'u' ) == Magma_tally4PackSubdiag   );
    check( magma_tally4_pack_const( 'L' ) == Magma_tally4PackSupdiag   );
    check( magma_tally4_pack_const( 'l' ) == Magma_tally4PackSupdiag   );
    check( magma_tally4_pack_const( 'C' ) == Magma_tally4PackColumn    );
    check( magma_tally4_pack_const( 'c' ) == Magma_tally4PackColumn    );
    check( magma_tally4_pack_const( 'R' ) == Magma_tally4PackRow       );
    check( magma_tally4_pack_const( 'r' ) == Magma_tally4PackRow       );
    check( magma_tally4_pack_const( 'B' ) == Magma_tally4PackLowerBand );
    check( magma_tally4_pack_const( 'b' ) == Magma_tally4PackLowerBand );
    check( magma_tally4_pack_const( 'Q' ) == Magma_tally4PackUpeprBand );
    check( magma_tally4_pack_const( 'q' ) == Magma_tally4PackUpeprBand );
    check( magma_tally4_pack_const( 'Z' ) == Magma_tally4PackAll       );
    check( magma_tally4_pack_const( 'z' ) == Magma_tally4PackAll       );

    check( magma_tally4_vec_const( 'N' )  == Magma_tally4NoVec         );
    check( magma_tally4_vec_const( 'n' )  == Magma_tally4NoVec         );
    check( magma_tally4_vec_const( 'V' )  == Magma_tally4Vec           );
    check( magma_tally4_vec_const( 'v' )  == Magma_tally4Vec           );
    check( magma_tally4_vec_const( 'I' )  == Magma_tally4IVec          );
    check( magma_tally4_vec_const( 'i' )  == Magma_tally4IVec          );
    check( magma_tally4_vec_const( 'A' )  == Magma_tally4AllVec        );
    check( magma_tally4_vec_const( 'a' )  == Magma_tally4AllVec        );
    check( magma_tally4_vec_const( 'S' )  == Magma_tally4SomeVec       );
    check( magma_tally4_vec_const( 's' )  == Magma_tally4SomeVec       );
    check( magma_tally4_vec_const( 'O' )  == Magma_tally4OverwriteVec  );
    check( magma_tally4_vec_const( 'o' )  == Magma_tally4OverwriteVec  );

    check( magma_tally4_range_const( 'A' )  == Magma_tally4RangeAll    );
    check( magma_tally4_range_const( 'a' )  == Magma_tally4RangeAll    );
    check( magma_tally4_range_const( 'V' )  == Magma_tally4RangeV      );
    check( magma_tally4_range_const( 'v' )  == Magma_tally4RangeV      );
    check( magma_tally4_range_const( 'I' )  == Magma_tally4RangeI      );
    check( magma_tally4_range_const( 'i' )  == Magma_tally4RangeI      );

    check( magma_tally4_vect_const( 'Q' )   == Magma_tally4Q           );
    check( magma_tally4_vect_const( 'q' )   == Magma_tally4Q           );
    check( magma_tally4_vect_const( 'P' )   == Magma_tally4P           );
    check( magma_tally4_vect_const( 'p' )   == Magma_tally4P           );

    check( magma_tally4_direct_const( 'F' ) == Magma_tally4Forward     );
    check( magma_tally4_direct_const( 'f' ) == Magma_tally4Forward     );
    check( magma_tally4_direct_const( 'B' ) == Magma_tally4Backward    );
    check( magma_tally4_direct_const( 'b' ) == Magma_tally4Backward    );

    check( magma_tally4_storev_const( 'C' ) == Magma_tally4Columnwise  );
    check( magma_tally4_storev_const( 'c' ) == Magma_tally4Columnwise  );
    check( magma_tally4_storev_const( 'R' ) == Magma_tally4Rowwise     );
    check( magma_tally4_storev_const( 'r' ) == Magma_tally4Rowwise     );
    printf( "LAPACK -> magma_tally4_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    #ifdef HAVE_clAmdBlas
    s = gStatus;
    check( amdblas_order_const( Magma_tally4RowMajor      ) == clAmdBlasRowMajor    );
    check( amdblas_order_const( Magma_tally4ColMajor      ) == clAmdBlasColumnMajor );

    check( amdblas_trans_const( Magma_tally4NoTrans       ) == clAmdBlasNoTrans     );
    check( amdblas_trans_const( Magma_tally4Trans         ) == clAmdBlasTrans       );
    check( amdblas_trans_const( Magma_tally4ConjTrans     ) == clAmdBlasConjTrans   );

    check( amdblas_uplo_const(  Magma_tally4Upper         ) == clAmdBlasUpper       );
    check( amdblas_uplo_const(  Magma_tally4Lower         ) == clAmdBlasLower       );

    check( amdblas_diag_const(  Magma_tally4NonUnit       ) == clAmdBlasNonUnit     );
    check( amdblas_diag_const(  Magma_tally4Unit          ) == clAmdBlasUnit        );

    check( amdblas_side_const(  Magma_tally4Left          ) == clAmdBlasLeft        );
    check( amdblas_side_const(  Magma_tally4Right         ) == clAmdBlasRight       );
    printf( "MAGMA_tally4  -> amdblas_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef CUBLAS_V2_H_
    s = gStatus;
    check( cublas_trans_const_tally4( Magma_tally4NoTrans       ) == CUBLAS_OP_N            );
    check( cublas_trans_const_tally4( Magma_tally4Trans         ) == CUBLAS_OP_T            );
    check( cublas_trans_const_tally4( Magma_tally4ConjTrans     ) == CUBLAS_OP_C            );

    check( cublas_uplo_const_tally4(  Magma_tally4Upper         ) == CUBLAS_FILL_MODE_UPPER );
    check( cublas_uplo_const_tally4(  Magma_tally4Lower         ) == CUBLAS_FILL_MODE_LOWER );

    check( cublas_diag_const_tally4(  Magma_tally4NonUnit       ) == CUBLAS_DIAG_NON_UNIT   );
    check( cublas_diag_const_tally4(  Magma_tally4Unit          ) == CUBLAS_DIAG_UNIT       );

    check( cublas_side_const_tally4(  Magma_tally4Left          ) == CUBLAS_SIDE_LEFT       );
    check( cublas_side_const_tally4(  Magma_tally4Right         ) == CUBLAS_SIDE_RIGHT      );
    printf( "MAGMA_tally4  -> cublas_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef HAVE_CBLAS
    s = gStatus;
    check( cblas_order_const( Magma_tally4RowMajor      ) == CblasRowMajor  );
    check( cblas_order_const( Magma_tally4ColMajor      ) == CblasColMajor  );

    check( cblas_trans_const( Magma_tally4NoTrans       ) == CblasNoTrans   );
    check( cblas_trans_const( Magma_tally4Trans         ) == CblasTrans     );
    check( cblas_trans_const( Magma_tally4ConjTrans     ) == CblasConjTrans );

    check( cblas_uplo_const(  Magma_tally4Upper         ) == CblasUpper     );
    check( cblas_uplo_const(  Magma_tally4Lower         ) == CblasLower     );

    check( cblas_diag_const(  Magma_tally4NonUnit       ) == CblasNonUnit   );
    check( cblas_diag_const(  Magma_tally4Unit          ) == CblasUnit      );

    check( cblas_side_const(  Magma_tally4Left          ) == CblasLeft      );
    check( cblas_side_const(  Magma_tally4Right         ) == CblasRight     );
    printf( "MAGMA_tally4  -> cblas_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));
    #endif

    return gStatus;
}
