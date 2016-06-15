/*
    -- MAGMA_tally3 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/
#include <stdio.h>

#include "testings.h"
#include "magma_tally3.h"

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
    check( lapack_bool_const_tally3(   Magma_tally3False         )[0] == 'N' );
    check( lapack_bool_const_tally3(   Magma_tally3True          )[0] == 'Y' );

    check( lapack_order_const_tally3(  Magma_tally3RowMajor      )[0] == 'R' );
    check( lapack_order_const_tally3(  Magma_tally3ColMajor      )[0] == 'C' );

    check( lapack_trans_const_tally3(  Magma_tally3NoTrans       )[0] == 'N' );
    check( lapack_trans_const_tally3(  Magma_tally3Trans         )[0] == 'T' );
    check( lapack_trans_const_tally3(  Magma_tally3ConjTrans     )[0] == 'C' );

    check( lapack_uplo_const_tally3(   Magma_tally3Upper         )[0] == 'U' );
    check( lapack_uplo_const_tally3(   Magma_tally3Lower         )[0] == 'L' );
    check( lapack_uplo_const_tally3(   Magma_tally3Full          )[0] == 'G' );

    check( lapack_diag_const_tally3(   Magma_tally3NonUnit       )[0] == 'N' );
    check( lapack_diag_const_tally3(   Magma_tally3Unit          )[0] == 'U' );

    check( lapack_side_const_tally3(   Magma_tally3Left          )[0] == 'L' );
    check( lapack_side_const_tally3(   Magma_tally3Right         )[0] == 'R' );
    check( lapack_side_const_tally3(   Magma_tally3BothSides     )[0] == 'B' );

    check( lapack_norm_const_tally3(   Magma_tally3OneNorm       )[0] == '1' );
    check( lapack_norm_const_tally3(   Magma_tally3TwoNorm       )[0] == '2' );
    check( lapack_norm_const_tally3(   Magma_tally3FrobeniusNorm )[0] == 'F' );
    check( lapack_norm_const_tally3(   Magma_tally3InfNorm       )[0] == 'I' );
    check( lapack_norm_const_tally3(   Magma_tally3MaxNorm       )[0] == 'M' );

    check( lapack_dist_const_tally3(   Magma_tally3DistUniform   )[0] == 'U' );
    check( lapack_dist_const_tally3(   Magma_tally3DistSymmetric )[0] == 'S' );
    check( lapack_dist_const_tally3(   Magma_tally3DistNormal    )[0] == 'N' );

    check( lapack_sym_const_tally3(    Magma_tally3HermGeev      )[0] == 'H' );
    check( lapack_sym_const_tally3(    Magma_tally3HermPoev      )[0] == 'P' );
    check( lapack_sym_const_tally3(    Magma_tally3NonsymPosv    )[0] == 'N' );
    check( lapack_sym_const_tally3(    Magma_tally3SymPosv       )[0] == 'S' );

    check( lapack_pack_const_tally3(   Magma_tally3NoPacking     )[0] == 'N' );
    check( lapack_pack_const_tally3(   Magma_tally3PackSubdiag   )[0] == 'U' );
    check( lapack_pack_const_tally3(   Magma_tally3PackSupdiag   )[0] == 'L' );
    check( lapack_pack_const_tally3(   Magma_tally3PackColumn    )[0] == 'C' );
    check( lapack_pack_const_tally3(   Magma_tally3PackRow       )[0] == 'R' );
    check( lapack_pack_const_tally3(   Magma_tally3PackLowerBand )[0] == 'B' );
    check( lapack_pack_const_tally3(   Magma_tally3PackUpeprBand )[0] == 'Q' );
    check( lapack_pack_const_tally3(   Magma_tally3PackAll       )[0] == 'Z' );

    check( lapack_vec_const_tally3(    Magma_tally3NoVec         )[0] == 'N' );
    check( lapack_vec_const_tally3(    Magma_tally3Vec           )[0] == 'V' );
    check( lapack_vec_const_tally3(    Magma_tally3IVec          )[0] == 'I' );
    check( lapack_vec_const_tally3(    Magma_tally3AllVec        )[0] == 'A' );
    check( lapack_vec_const_tally3(    Magma_tally3SomeVec       )[0] == 'S' );
    check( lapack_vec_const_tally3(    Magma_tally3OverwriteVec  )[0] == 'O' );

    check( lapack_range_const_tally3(  Magma_tally3RangeAll      )[0] == 'A' );
    check( lapack_range_const_tally3(  Magma_tally3RangeV        )[0] == 'V' );
    check( lapack_range_const_tally3(  Magma_tally3RangeI        )[0] == 'I' );

    check( lapack_vect_const_tally3(   Magma_tally3Q             )[0] == 'Q' );
    check( lapack_vect_const_tally3(   Magma_tally3P             )[0] == 'P' );

    check( lapack_direct_const_tally3( Magma_tally3Forward       )[0] == 'F' );
    check( lapack_direct_const_tally3( Magma_tally3Backward      )[0] == 'B' );

    check( lapack_storev_const_tally3( Magma_tally3Columnwise    )[0] == 'C' );
    check( lapack_storev_const_tally3( Magma_tally3Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally3  -> lapack_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_bool_const_tally3(   Magma_tally3False         ) == 'N' );
    check( lapacke_bool_const_tally3(   Magma_tally3True          ) == 'Y' );

    check( lapacke_order_const_tally3(  Magma_tally3RowMajor      ) == 'R' );
    check( lapacke_order_const_tally3(  Magma_tally3ColMajor      ) == 'C' );

    check( lapacke_trans_const_tally3(  Magma_tally3NoTrans       ) == 'N' );
    check( lapacke_trans_const_tally3(  Magma_tally3Trans         ) == 'T' );
    check( lapacke_trans_const_tally3(  Magma_tally3ConjTrans     ) == 'C' );

    check( lapacke_uplo_const_tally3(   Magma_tally3Upper         ) == 'U' );
    check( lapacke_uplo_const_tally3(   Magma_tally3Lower         ) == 'L' );
    check( lapacke_uplo_const_tally3(   Magma_tally3Full          ) == 'G' );

    check( lapacke_diag_const_tally3(   Magma_tally3NonUnit       ) == 'N' );
    check( lapacke_diag_const_tally3(   Magma_tally3Unit          ) == 'U' );

    check( lapacke_side_const_tally3(   Magma_tally3Left          ) == 'L' );
    check( lapacke_side_const_tally3(   Magma_tally3Right         ) == 'R' );
    check( lapacke_side_const_tally3(   Magma_tally3BothSides     ) == 'B' );

    check( lapacke_norm_const_tally3(   Magma_tally3OneNorm       ) == '1' );
    check( lapacke_norm_const_tally3(   Magma_tally3TwoNorm       ) == '2' );
    check( lapacke_norm_const_tally3(   Magma_tally3FrobeniusNorm ) == 'F' );
    check( lapacke_norm_const_tally3(   Magma_tally3InfNorm       ) == 'I' );
    check( lapacke_norm_const_tally3(   Magma_tally3MaxNorm       ) == 'M' );

    check( lapacke_dist_const_tally3(   Magma_tally3DistUniform   ) == 'U' );
    check( lapacke_dist_const_tally3(   Magma_tally3DistSymmetric ) == 'S' );
    check( lapacke_dist_const_tally3(   Magma_tally3DistNormal    ) == 'N' );

    check( lapacke_sym_const_tally3(    Magma_tally3HermGeev      ) == 'H' );
    check( lapacke_sym_const_tally3(    Magma_tally3HermPoev      ) == 'P' );
    check( lapacke_sym_const_tally3(    Magma_tally3NonsymPosv    ) == 'N' );
    check( lapacke_sym_const_tally3(    Magma_tally3SymPosv       ) == 'S' );

    check( lapacke_pack_const_tally3(   Magma_tally3NoPacking     ) == 'N' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackSubdiag   ) == 'U' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackSupdiag   ) == 'L' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackColumn    ) == 'C' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackRow       ) == 'R' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackLowerBand ) == 'B' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackUpeprBand ) == 'Q' );
    check( lapacke_pack_const_tally3(   Magma_tally3PackAll       ) == 'Z' );

    check( lapacke_vec_const_tally3(    Magma_tally3NoVec         ) == 'N' );
    check( lapacke_vec_const_tally3(    Magma_tally3Vec           ) == 'V' );
    check( lapacke_vec_const_tally3(    Magma_tally3IVec          ) == 'I' );
    check( lapacke_vec_const_tally3(    Magma_tally3AllVec        ) == 'A' );
    check( lapacke_vec_const_tally3(    Magma_tally3SomeVec       ) == 'S' );
    check( lapacke_vec_const_tally3(    Magma_tally3OverwriteVec  ) == 'O' );

    check( lapacke_range_const_tally3(  Magma_tally3RangeAll      ) == 'A' );
    check( lapacke_range_const_tally3(  Magma_tally3RangeV        ) == 'V' );
    check( lapacke_range_const_tally3(  Magma_tally3RangeI        ) == 'I' );

    check( lapacke_vect_const_tally3(   Magma_tally3Q             ) == 'Q' );
    check( lapacke_vect_const_tally3(   Magma_tally3P             ) == 'P' );

    check( lapacke_direct_const_tally3( Magma_tally3Forward       ) == 'F' );
    check( lapacke_direct_const_tally3( Magma_tally3Backward      ) == 'B' );

    check( lapacke_storev_const_tally3( Magma_tally3Columnwise    ) == 'C' );
    check( lapacke_storev_const_tally3( Magma_tally3Rowwise       ) == 'R' );
    printf( "MAGMA_tally3  -> lapacke_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapack_const_tally3( Magma_tally3False         )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3True          )[0] == 'Y' );

    check( lapack_const_tally3( Magma_tally3RowMajor      )[0] == 'R' );
    check( lapack_const_tally3( Magma_tally3ColMajor      )[0] == 'C' );

    check( lapack_const_tally3( Magma_tally3NoTrans       )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3Trans         )[0] == 'T' );
    check( lapack_const_tally3( Magma_tally3ConjTrans     )[0] == 'C' );

    check( lapack_const_tally3( Magma_tally3Upper         )[0] == 'U' );
    check( lapack_const_tally3( Magma_tally3Lower         )[0] == 'L' );
    check( lapack_const_tally3( Magma_tally3Full          )[0] == 'G' );

    check( lapack_const_tally3( Magma_tally3NonUnit       )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3Unit          )[0] == 'U' );

    check( lapack_const_tally3( Magma_tally3Left          )[0] == 'L' );
    check( lapack_const_tally3( Magma_tally3Right         )[0] == 'R' );
    check( lapack_const_tally3( Magma_tally3BothSides     )[0] == 'B' );

    check( lapack_const_tally3( Magma_tally3OneNorm       )[0] == '1' );
    check( lapack_const_tally3( Magma_tally3TwoNorm       )[0] == '2' );
    check( lapack_const_tally3( Magma_tally3FrobeniusNorm )[0] == 'F' );
    check( lapack_const_tally3( Magma_tally3InfNorm       )[0] == 'I' );
    check( lapack_const_tally3( Magma_tally3MaxNorm       )[0] == 'M' );

    check( lapack_const_tally3( Magma_tally3DistUniform   )[0] == 'U' );
    check( lapack_const_tally3( Magma_tally3DistSymmetric )[0] == 'S' );
    check( lapack_const_tally3( Magma_tally3DistNormal    )[0] == 'N' );

    check( lapack_const_tally3( Magma_tally3HermGeev      )[0] == 'H' );
    check( lapack_const_tally3( Magma_tally3HermPoev      )[0] == 'P' );
    check( lapack_const_tally3( Magma_tally3NonsymPosv    )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3SymPosv       )[0] == 'S' );

    check( lapack_const_tally3( Magma_tally3NoPacking     )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3PackSubdiag   )[0] == 'U' );
    check( lapack_const_tally3( Magma_tally3PackSupdiag   )[0] == 'L' );
    check( lapack_const_tally3( Magma_tally3PackColumn    )[0] == 'C' );
    check( lapack_const_tally3( Magma_tally3PackRow       )[0] == 'R' );
    check( lapack_const_tally3( Magma_tally3PackLowerBand )[0] == 'B' );
    check( lapack_const_tally3( Magma_tally3PackUpeprBand )[0] == 'Q' );
    check( lapack_const_tally3( Magma_tally3PackAll       )[0] == 'Z' );

    check( lapack_const_tally3( Magma_tally3NoVec         )[0] == 'N' );
    check( lapack_const_tally3( Magma_tally3Vec           )[0] == 'V' );
    check( lapack_const_tally3( Magma_tally3IVec          )[0] == 'I' );
    check( lapack_const_tally3( Magma_tally3AllVec        )[0] == 'A' );
    check( lapack_const_tally3( Magma_tally3SomeVec       )[0] == 'S' );
    check( lapack_const_tally3( Magma_tally3OverwriteVec  )[0] == 'O' );

    check( lapack_const_tally3( Magma_tally3RangeAll      )[0] == 'A' );
    check( lapack_const_tally3( Magma_tally3RangeV        )[0] == 'V' );
    check( lapack_const_tally3( Magma_tally3RangeI        )[0] == 'I' );

    check( lapack_const_tally3( Magma_tally3Q             )[0] == 'Q' );
    check( lapack_const_tally3( Magma_tally3P             )[0] == 'P' );

    check( lapack_const_tally3( Magma_tally3Forward       )[0] == 'F' );
    check( lapack_const_tally3( Magma_tally3Backward      )[0] == 'B' );

    check( lapack_const_tally3( Magma_tally3Columnwise    )[0] == 'C' );
    check( lapack_const_tally3( Magma_tally3Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally3  -> lapack_const_tally3          %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_const_tally3( Magma_tally3False         ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3True          ) == 'Y' );

    check( lapacke_const_tally3( Magma_tally3RowMajor      ) == 'R' );
    check( lapacke_const_tally3( Magma_tally3ColMajor      ) == 'C' );

    check( lapacke_const_tally3( Magma_tally3NoTrans       ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3Trans         ) == 'T' );
    check( lapacke_const_tally3( Magma_tally3ConjTrans     ) == 'C' );

    check( lapacke_const_tally3( Magma_tally3Upper         ) == 'U' );
    check( lapacke_const_tally3( Magma_tally3Lower         ) == 'L' );
    check( lapacke_const_tally3( Magma_tally3Full          ) == 'G' );

    check( lapacke_const_tally3( Magma_tally3NonUnit       ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3Unit          ) == 'U' );

    check( lapacke_const_tally3( Magma_tally3Left          ) == 'L' );
    check( lapacke_const_tally3( Magma_tally3Right         ) == 'R' );
    check( lapacke_const_tally3( Magma_tally3BothSides     ) == 'B' );

    check( lapacke_const_tally3( Magma_tally3OneNorm       ) == '1' );
    check( lapacke_const_tally3( Magma_tally3TwoNorm       ) == '2' );
    check( lapacke_const_tally3( Magma_tally3FrobeniusNorm ) == 'F' );
    check( lapacke_const_tally3( Magma_tally3InfNorm       ) == 'I' );
    check( lapacke_const_tally3( Magma_tally3MaxNorm       ) == 'M' );

    check( lapacke_const_tally3( Magma_tally3DistUniform   ) == 'U' );
    check( lapacke_const_tally3( Magma_tally3DistSymmetric ) == 'S' );
    check( lapacke_const_tally3( Magma_tally3DistNormal    ) == 'N' );

    check( lapacke_const_tally3( Magma_tally3HermGeev      ) == 'H' );
    check( lapacke_const_tally3( Magma_tally3HermPoev      ) == 'P' );
    check( lapacke_const_tally3( Magma_tally3NonsymPosv    ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3SymPosv       ) == 'S' );

    check( lapacke_const_tally3( Magma_tally3NoPacking     ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3PackSubdiag   ) == 'U' );
    check( lapacke_const_tally3( Magma_tally3PackSupdiag   ) == 'L' );
    check( lapacke_const_tally3( Magma_tally3PackColumn    ) == 'C' );
    check( lapacke_const_tally3( Magma_tally3PackRow       ) == 'R' );
    check( lapacke_const_tally3( Magma_tally3PackLowerBand ) == 'B' );
    check( lapacke_const_tally3( Magma_tally3PackUpeprBand ) == 'Q' );
    check( lapacke_const_tally3( Magma_tally3PackAll       ) == 'Z' );

    check( lapacke_const_tally3( Magma_tally3NoVec         ) == 'N' );
    check( lapacke_const_tally3( Magma_tally3Vec           ) == 'V' );
    check( lapacke_const_tally3( Magma_tally3IVec          ) == 'I' );
    check( lapacke_const_tally3( Magma_tally3AllVec        ) == 'A' );
    check( lapacke_const_tally3( Magma_tally3SomeVec       ) == 'S' );
    check( lapacke_const_tally3( Magma_tally3OverwriteVec  ) == 'O' );

    check( lapacke_const_tally3( Magma_tally3RangeAll      ) == 'A' );
    check( lapacke_const_tally3( Magma_tally3RangeV        ) == 'V' );
    check( lapacke_const_tally3( Magma_tally3RangeI        ) == 'I' );

    check( lapacke_const_tally3( Magma_tally3Q             ) == 'Q' );
    check( lapacke_const_tally3( Magma_tally3P             ) == 'P' );

    check( lapacke_const_tally3( Magma_tally3Forward       ) == 'F' );
    check( lapacke_const_tally3( Magma_tally3Backward      ) == 'B' );

    check( lapacke_const_tally3( Magma_tally3Columnwise    ) == 'C' );
    check( lapacke_const_tally3( Magma_tally3Rowwise       ) == 'R' );
    printf( "MAGMA_tally3  -> lapacke_const_tally3         %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( magma_tally3_bool_const('N') == Magma_tally3False );
    check( magma_tally3_bool_const('n') == Magma_tally3False );
    check( magma_tally3_bool_const('Y') == Magma_tally3True  );
    check( magma_tally3_bool_const('y') == Magma_tally3True  );

    check( magma_tally3_order_const( 'R' ) == Magma_tally3RowMajor  );
    check( magma_tally3_order_const( 'r' ) == Magma_tally3RowMajor  );
    check( magma_tally3_order_const( 'C' ) == Magma_tally3ColMajor  );
    check( magma_tally3_order_const( 'c' ) == Magma_tally3ColMajor  );

    check( magma_tally3_trans_const( 'N' ) == Magma_tally3NoTrans   );
    check( magma_tally3_trans_const( 'n' ) == Magma_tally3NoTrans   );
    check( magma_tally3_trans_const( 'T' ) == Magma_tally3Trans     );
    check( magma_tally3_trans_const( 't' ) == Magma_tally3Trans     );
    check( magma_tally3_trans_const( 'C' ) == Magma_tally3ConjTrans );
    check( magma_tally3_trans_const( 'c' ) == Magma_tally3ConjTrans );

    check( magma_tally3_uplo_const( 'U' ) == Magma_tally3Upper      );
    check( magma_tally3_uplo_const( 'u' ) == Magma_tally3Upper      );
    check( magma_tally3_uplo_const( 'L' ) == Magma_tally3Lower      );
    check( magma_tally3_uplo_const( 'l' ) == Magma_tally3Lower      );
    check( magma_tally3_uplo_const( 'A' ) == Magma_tally3Full       );  // anything else
    check( magma_tally3_uplo_const( 'a' ) == Magma_tally3Full       );
    check( magma_tally3_uplo_const( 'G' ) == Magma_tally3Full       );
    check( magma_tally3_uplo_const( 'g' ) == Magma_tally3Full       );
    check( magma_tally3_uplo_const( 'F' ) == Magma_tally3Full       );
    check( magma_tally3_uplo_const( 'f' ) == Magma_tally3Full       );

    check( magma_tally3_diag_const( 'N' ) == Magma_tally3NonUnit    );
    check( magma_tally3_diag_const( 'n' ) == Magma_tally3NonUnit    );
    check( magma_tally3_diag_const( 'U' ) == Magma_tally3Unit       );
    check( magma_tally3_diag_const( 'u' ) == Magma_tally3Unit       );

    check( magma_tally3_side_const( 'L' ) == Magma_tally3Left       );
    check( magma_tally3_side_const( 'l' ) == Magma_tally3Left       );
    check( magma_tally3_side_const( 'R' ) == Magma_tally3Right      );
    check( magma_tally3_side_const( 'r' ) == Magma_tally3Right      );

    check( magma_tally3_norm_const( 'O' ) == Magma_tally3OneNorm       );
    check( magma_tally3_norm_const( 'o' ) == Magma_tally3OneNorm       );
    check( magma_tally3_norm_const( '1' ) == Magma_tally3OneNorm       );
    check( magma_tally3_norm_const( '2' ) == Magma_tally3TwoNorm       );
    check( magma_tally3_norm_const( 'F' ) == Magma_tally3FrobeniusNorm );
    check( magma_tally3_norm_const( 'f' ) == Magma_tally3FrobeniusNorm );
    check( magma_tally3_norm_const( 'E' ) == Magma_tally3FrobeniusNorm );
    check( magma_tally3_norm_const( 'e' ) == Magma_tally3FrobeniusNorm );
    check( magma_tally3_norm_const( 'I' ) == Magma_tally3InfNorm       );
    check( magma_tally3_norm_const( 'i' ) == Magma_tally3InfNorm       );
    check( magma_tally3_norm_const( 'M' ) == Magma_tally3MaxNorm       );
    check( magma_tally3_norm_const( 'm' ) == Magma_tally3MaxNorm       );

    check( magma_tally3_dist_const( 'U' ) == Magma_tally3DistUniform   );
    check( magma_tally3_dist_const( 'u' ) == Magma_tally3DistUniform   );
    check( magma_tally3_dist_const( 'S' ) == Magma_tally3DistSymmetric );
    check( magma_tally3_dist_const( 's' ) == Magma_tally3DistSymmetric );
    check( magma_tally3_dist_const( 'N' ) == Magma_tally3DistNormal    );
    check( magma_tally3_dist_const( 'n' ) == Magma_tally3DistNormal    );

    //check( magma_tally3_xxxx_const( 'H' ) == Magma_tally3HermGeev      );
    //check( magma_tally3_xxxx_const( 'P' ) == Magma_tally3HermPoev      );
    //check( magma_tally3_xxxx_const( 'N' ) == Magma_tally3NonsymPosv    );
    //check( magma_tally3_xxxx_const( 'S' ) == Magma_tally3SymPosv       );

    check( magma_tally3_pack_const( 'N' ) == Magma_tally3NoPacking     );
    check( magma_tally3_pack_const( 'n' ) == Magma_tally3NoPacking     );
    check( magma_tally3_pack_const( 'U' ) == Magma_tally3PackSubdiag   );
    check( magma_tally3_pack_const( 'u' ) == Magma_tally3PackSubdiag   );
    check( magma_tally3_pack_const( 'L' ) == Magma_tally3PackSupdiag   );
    check( magma_tally3_pack_const( 'l' ) == Magma_tally3PackSupdiag   );
    check( magma_tally3_pack_const( 'C' ) == Magma_tally3PackColumn    );
    check( magma_tally3_pack_const( 'c' ) == Magma_tally3PackColumn    );
    check( magma_tally3_pack_const( 'R' ) == Magma_tally3PackRow       );
    check( magma_tally3_pack_const( 'r' ) == Magma_tally3PackRow       );
    check( magma_tally3_pack_const( 'B' ) == Magma_tally3PackLowerBand );
    check( magma_tally3_pack_const( 'b' ) == Magma_tally3PackLowerBand );
    check( magma_tally3_pack_const( 'Q' ) == Magma_tally3PackUpeprBand );
    check( magma_tally3_pack_const( 'q' ) == Magma_tally3PackUpeprBand );
    check( magma_tally3_pack_const( 'Z' ) == Magma_tally3PackAll       );
    check( magma_tally3_pack_const( 'z' ) == Magma_tally3PackAll       );

    check( magma_tally3_vec_const( 'N' )  == Magma_tally3NoVec         );
    check( magma_tally3_vec_const( 'n' )  == Magma_tally3NoVec         );
    check( magma_tally3_vec_const( 'V' )  == Magma_tally3Vec           );
    check( magma_tally3_vec_const( 'v' )  == Magma_tally3Vec           );
    check( magma_tally3_vec_const( 'I' )  == Magma_tally3IVec          );
    check( magma_tally3_vec_const( 'i' )  == Magma_tally3IVec          );
    check( magma_tally3_vec_const( 'A' )  == Magma_tally3AllVec        );
    check( magma_tally3_vec_const( 'a' )  == Magma_tally3AllVec        );
    check( magma_tally3_vec_const( 'S' )  == Magma_tally3SomeVec       );
    check( magma_tally3_vec_const( 's' )  == Magma_tally3SomeVec       );
    check( magma_tally3_vec_const( 'O' )  == Magma_tally3OverwriteVec  );
    check( magma_tally3_vec_const( 'o' )  == Magma_tally3OverwriteVec  );

    check( magma_tally3_range_const( 'A' )  == Magma_tally3RangeAll    );
    check( magma_tally3_range_const( 'a' )  == Magma_tally3RangeAll    );
    check( magma_tally3_range_const( 'V' )  == Magma_tally3RangeV      );
    check( magma_tally3_range_const( 'v' )  == Magma_tally3RangeV      );
    check( magma_tally3_range_const( 'I' )  == Magma_tally3RangeI      );
    check( magma_tally3_range_const( 'i' )  == Magma_tally3RangeI      );

    check( magma_tally3_vect_const( 'Q' )   == Magma_tally3Q           );
    check( magma_tally3_vect_const( 'q' )   == Magma_tally3Q           );
    check( magma_tally3_vect_const( 'P' )   == Magma_tally3P           );
    check( magma_tally3_vect_const( 'p' )   == Magma_tally3P           );

    check( magma_tally3_direct_const( 'F' ) == Magma_tally3Forward     );
    check( magma_tally3_direct_const( 'f' ) == Magma_tally3Forward     );
    check( magma_tally3_direct_const( 'B' ) == Magma_tally3Backward    );
    check( magma_tally3_direct_const( 'b' ) == Magma_tally3Backward    );

    check( magma_tally3_storev_const( 'C' ) == Magma_tally3Columnwise  );
    check( magma_tally3_storev_const( 'c' ) == Magma_tally3Columnwise  );
    check( magma_tally3_storev_const( 'R' ) == Magma_tally3Rowwise     );
    check( magma_tally3_storev_const( 'r' ) == Magma_tally3Rowwise     );
    printf( "LAPACK -> magma_tally3_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    #ifdef HAVE_clAmdBlas
    s = gStatus;
    check( amdblas_order_const( Magma_tally3RowMajor      ) == clAmdBlasRowMajor    );
    check( amdblas_order_const( Magma_tally3ColMajor      ) == clAmdBlasColumnMajor );

    check( amdblas_trans_const( Magma_tally3NoTrans       ) == clAmdBlasNoTrans     );
    check( amdblas_trans_const( Magma_tally3Trans         ) == clAmdBlasTrans       );
    check( amdblas_trans_const( Magma_tally3ConjTrans     ) == clAmdBlasConjTrans   );

    check( amdblas_uplo_const(  Magma_tally3Upper         ) == clAmdBlasUpper       );
    check( amdblas_uplo_const(  Magma_tally3Lower         ) == clAmdBlasLower       );

    check( amdblas_diag_const(  Magma_tally3NonUnit       ) == clAmdBlasNonUnit     );
    check( amdblas_diag_const(  Magma_tally3Unit          ) == clAmdBlasUnit        );

    check( amdblas_side_const(  Magma_tally3Left          ) == clAmdBlasLeft        );
    check( amdblas_side_const(  Magma_tally3Right         ) == clAmdBlasRight       );
    printf( "MAGMA_tally3  -> amdblas_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef CUBLAS_V2_H_
    s = gStatus;
    check( cublas_trans_const_tally3( Magma_tally3NoTrans       ) == CUBLAS_OP_N            );
    check( cublas_trans_const_tally3( Magma_tally3Trans         ) == CUBLAS_OP_T            );
    check( cublas_trans_const_tally3( Magma_tally3ConjTrans     ) == CUBLAS_OP_C            );

    check( cublas_uplo_const_tally3(  Magma_tally3Upper         ) == CUBLAS_FILL_MODE_UPPER );
    check( cublas_uplo_const_tally3(  Magma_tally3Lower         ) == CUBLAS_FILL_MODE_LOWER );

    check( cublas_diag_const_tally3(  Magma_tally3NonUnit       ) == CUBLAS_DIAG_NON_UNIT   );
    check( cublas_diag_const_tally3(  Magma_tally3Unit          ) == CUBLAS_DIAG_UNIT       );

    check( cublas_side_const_tally3(  Magma_tally3Left          ) == CUBLAS_SIDE_LEFT       );
    check( cublas_side_const_tally3(  Magma_tally3Right         ) == CUBLAS_SIDE_RIGHT      );
    printf( "MAGMA_tally3  -> cublas_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef HAVE_CBLAS
    s = gStatus;
    check( cblas_order_const( Magma_tally3RowMajor      ) == CblasRowMajor  );
    check( cblas_order_const( Magma_tally3ColMajor      ) == CblasColMajor  );

    check( cblas_trans_const( Magma_tally3NoTrans       ) == CblasNoTrans   );
    check( cblas_trans_const( Magma_tally3Trans         ) == CblasTrans     );
    check( cblas_trans_const( Magma_tally3ConjTrans     ) == CblasConjTrans );

    check( cblas_uplo_const(  Magma_tally3Upper         ) == CblasUpper     );
    check( cblas_uplo_const(  Magma_tally3Lower         ) == CblasLower     );

    check( cblas_diag_const(  Magma_tally3NonUnit       ) == CblasNonUnit   );
    check( cblas_diag_const(  Magma_tally3Unit          ) == CblasUnit      );

    check( cblas_side_const(  Magma_tally3Left          ) == CblasLeft      );
    check( cblas_side_const(  Magma_tally3Right         ) == CblasRight     );
    printf( "MAGMA_tally3  -> cblas_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));
    #endif

    return gStatus;
}
