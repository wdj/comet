/*
    -- MAGMA_tally2 (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/
#include <stdio.h>

#include "testings.h"
#include "magma_tally2.h"

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
    check( lapack_bool_const_tally2(   Magma_tally2False         )[0] == 'N' );
    check( lapack_bool_const_tally2(   Magma_tally2True          )[0] == 'Y' );

    check( lapack_order_const_tally2(  Magma_tally2RowMajor      )[0] == 'R' );
    check( lapack_order_const_tally2(  Magma_tally2ColMajor      )[0] == 'C' );

    check( lapack_trans_const_tally2(  Magma_tally2NoTrans       )[0] == 'N' );
    check( lapack_trans_const_tally2(  Magma_tally2Trans         )[0] == 'T' );
    check( lapack_trans_const_tally2(  Magma_tally2ConjTrans     )[0] == 'C' );

    check( lapack_uplo_const_tally2(   Magma_tally2Upper         )[0] == 'U' );
    check( lapack_uplo_const_tally2(   Magma_tally2Lower         )[0] == 'L' );
    check( lapack_uplo_const_tally2(   Magma_tally2Full          )[0] == 'G' );

    check( lapack_diag_const_tally2(   Magma_tally2NonUnit       )[0] == 'N' );
    check( lapack_diag_const_tally2(   Magma_tally2Unit          )[0] == 'U' );

    check( lapack_side_const_tally2(   Magma_tally2Left          )[0] == 'L' );
    check( lapack_side_const_tally2(   Magma_tally2Right         )[0] == 'R' );
    check( lapack_side_const_tally2(   Magma_tally2BothSides     )[0] == 'B' );

    check( lapack_norm_const_tally2(   Magma_tally2OneNorm       )[0] == '1' );
    check( lapack_norm_const_tally2(   Magma_tally2TwoNorm       )[0] == '2' );
    check( lapack_norm_const_tally2(   Magma_tally2FrobeniusNorm )[0] == 'F' );
    check( lapack_norm_const_tally2(   Magma_tally2InfNorm       )[0] == 'I' );
    check( lapack_norm_const_tally2(   Magma_tally2MaxNorm       )[0] == 'M' );

    check( lapack_dist_const_tally2(   Magma_tally2DistUniform   )[0] == 'U' );
    check( lapack_dist_const_tally2(   Magma_tally2DistSymmetric )[0] == 'S' );
    check( lapack_dist_const_tally2(   Magma_tally2DistNormal    )[0] == 'N' );

    check( lapack_sym_const_tally2(    Magma_tally2HermGeev      )[0] == 'H' );
    check( lapack_sym_const_tally2(    Magma_tally2HermPoev      )[0] == 'P' );
    check( lapack_sym_const_tally2(    Magma_tally2NonsymPosv    )[0] == 'N' );
    check( lapack_sym_const_tally2(    Magma_tally2SymPosv       )[0] == 'S' );

    check( lapack_pack_const_tally2(   Magma_tally2NoPacking     )[0] == 'N' );
    check( lapack_pack_const_tally2(   Magma_tally2PackSubdiag   )[0] == 'U' );
    check( lapack_pack_const_tally2(   Magma_tally2PackSupdiag   )[0] == 'L' );
    check( lapack_pack_const_tally2(   Magma_tally2PackColumn    )[0] == 'C' );
    check( lapack_pack_const_tally2(   Magma_tally2PackRow       )[0] == 'R' );
    check( lapack_pack_const_tally2(   Magma_tally2PackLowerBand )[0] == 'B' );
    check( lapack_pack_const_tally2(   Magma_tally2PackUpeprBand )[0] == 'Q' );
    check( lapack_pack_const_tally2(   Magma_tally2PackAll       )[0] == 'Z' );

    check( lapack_vec_const_tally2(    Magma_tally2NoVec         )[0] == 'N' );
    check( lapack_vec_const_tally2(    Magma_tally2Vec           )[0] == 'V' );
    check( lapack_vec_const_tally2(    Magma_tally2IVec          )[0] == 'I' );
    check( lapack_vec_const_tally2(    Magma_tally2AllVec        )[0] == 'A' );
    check( lapack_vec_const_tally2(    Magma_tally2SomeVec       )[0] == 'S' );
    check( lapack_vec_const_tally2(    Magma_tally2OverwriteVec  )[0] == 'O' );

    check( lapack_range_const_tally2(  Magma_tally2RangeAll      )[0] == 'A' );
    check( lapack_range_const_tally2(  Magma_tally2RangeV        )[0] == 'V' );
    check( lapack_range_const_tally2(  Magma_tally2RangeI        )[0] == 'I' );

    check( lapack_vect_const_tally2(   Magma_tally2Q             )[0] == 'Q' );
    check( lapack_vect_const_tally2(   Magma_tally2P             )[0] == 'P' );

    check( lapack_direct_const_tally2( Magma_tally2Forward       )[0] == 'F' );
    check( lapack_direct_const_tally2( Magma_tally2Backward      )[0] == 'B' );

    check( lapack_storev_const_tally2( Magma_tally2Columnwise    )[0] == 'C' );
    check( lapack_storev_const_tally2( Magma_tally2Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally2  -> lapack_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_bool_const_tally2(   Magma_tally2False         ) == 'N' );
    check( lapacke_bool_const_tally2(   Magma_tally2True          ) == 'Y' );

    check( lapacke_order_const_tally2(  Magma_tally2RowMajor      ) == 'R' );
    check( lapacke_order_const_tally2(  Magma_tally2ColMajor      ) == 'C' );

    check( lapacke_trans_const_tally2(  Magma_tally2NoTrans       ) == 'N' );
    check( lapacke_trans_const_tally2(  Magma_tally2Trans         ) == 'T' );
    check( lapacke_trans_const_tally2(  Magma_tally2ConjTrans     ) == 'C' );

    check( lapacke_uplo_const_tally2(   Magma_tally2Upper         ) == 'U' );
    check( lapacke_uplo_const_tally2(   Magma_tally2Lower         ) == 'L' );
    check( lapacke_uplo_const_tally2(   Magma_tally2Full          ) == 'G' );

    check( lapacke_diag_const_tally2(   Magma_tally2NonUnit       ) == 'N' );
    check( lapacke_diag_const_tally2(   Magma_tally2Unit          ) == 'U' );

    check( lapacke_side_const_tally2(   Magma_tally2Left          ) == 'L' );
    check( lapacke_side_const_tally2(   Magma_tally2Right         ) == 'R' );
    check( lapacke_side_const_tally2(   Magma_tally2BothSides     ) == 'B' );

    check( lapacke_norm_const_tally2(   Magma_tally2OneNorm       ) == '1' );
    check( lapacke_norm_const_tally2(   Magma_tally2TwoNorm       ) == '2' );
    check( lapacke_norm_const_tally2(   Magma_tally2FrobeniusNorm ) == 'F' );
    check( lapacke_norm_const_tally2(   Magma_tally2InfNorm       ) == 'I' );
    check( lapacke_norm_const_tally2(   Magma_tally2MaxNorm       ) == 'M' );

    check( lapacke_dist_const_tally2(   Magma_tally2DistUniform   ) == 'U' );
    check( lapacke_dist_const_tally2(   Magma_tally2DistSymmetric ) == 'S' );
    check( lapacke_dist_const_tally2(   Magma_tally2DistNormal    ) == 'N' );

    check( lapacke_sym_const_tally2(    Magma_tally2HermGeev      ) == 'H' );
    check( lapacke_sym_const_tally2(    Magma_tally2HermPoev      ) == 'P' );
    check( lapacke_sym_const_tally2(    Magma_tally2NonsymPosv    ) == 'N' );
    check( lapacke_sym_const_tally2(    Magma_tally2SymPosv       ) == 'S' );

    check( lapacke_pack_const_tally2(   Magma_tally2NoPacking     ) == 'N' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackSubdiag   ) == 'U' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackSupdiag   ) == 'L' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackColumn    ) == 'C' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackRow       ) == 'R' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackLowerBand ) == 'B' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackUpeprBand ) == 'Q' );
    check( lapacke_pack_const_tally2(   Magma_tally2PackAll       ) == 'Z' );

    check( lapacke_vec_const_tally2(    Magma_tally2NoVec         ) == 'N' );
    check( lapacke_vec_const_tally2(    Magma_tally2Vec           ) == 'V' );
    check( lapacke_vec_const_tally2(    Magma_tally2IVec          ) == 'I' );
    check( lapacke_vec_const_tally2(    Magma_tally2AllVec        ) == 'A' );
    check( lapacke_vec_const_tally2(    Magma_tally2SomeVec       ) == 'S' );
    check( lapacke_vec_const_tally2(    Magma_tally2OverwriteVec  ) == 'O' );

    check( lapacke_range_const_tally2(  Magma_tally2RangeAll      ) == 'A' );
    check( lapacke_range_const_tally2(  Magma_tally2RangeV        ) == 'V' );
    check( lapacke_range_const_tally2(  Magma_tally2RangeI        ) == 'I' );

    check( lapacke_vect_const_tally2(   Magma_tally2Q             ) == 'Q' );
    check( lapacke_vect_const_tally2(   Magma_tally2P             ) == 'P' );

    check( lapacke_direct_const_tally2( Magma_tally2Forward       ) == 'F' );
    check( lapacke_direct_const_tally2( Magma_tally2Backward      ) == 'B' );

    check( lapacke_storev_const_tally2( Magma_tally2Columnwise    ) == 'C' );
    check( lapacke_storev_const_tally2( Magma_tally2Rowwise       ) == 'R' );
    printf( "MAGMA_tally2  -> lapacke_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapack_const_tally2( Magma_tally2False         )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2True          )[0] == 'Y' );

    check( lapack_const_tally2( Magma_tally2RowMajor      )[0] == 'R' );
    check( lapack_const_tally2( Magma_tally2ColMajor      )[0] == 'C' );

    check( lapack_const_tally2( Magma_tally2NoTrans       )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2Trans         )[0] == 'T' );
    check( lapack_const_tally2( Magma_tally2ConjTrans     )[0] == 'C' );

    check( lapack_const_tally2( Magma_tally2Upper         )[0] == 'U' );
    check( lapack_const_tally2( Magma_tally2Lower         )[0] == 'L' );
    check( lapack_const_tally2( Magma_tally2Full          )[0] == 'G' );

    check( lapack_const_tally2( Magma_tally2NonUnit       )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2Unit          )[0] == 'U' );

    check( lapack_const_tally2( Magma_tally2Left          )[0] == 'L' );
    check( lapack_const_tally2( Magma_tally2Right         )[0] == 'R' );
    check( lapack_const_tally2( Magma_tally2BothSides     )[0] == 'B' );

    check( lapack_const_tally2( Magma_tally2OneNorm       )[0] == '1' );
    check( lapack_const_tally2( Magma_tally2TwoNorm       )[0] == '2' );
    check( lapack_const_tally2( Magma_tally2FrobeniusNorm )[0] == 'F' );
    check( lapack_const_tally2( Magma_tally2InfNorm       )[0] == 'I' );
    check( lapack_const_tally2( Magma_tally2MaxNorm       )[0] == 'M' );

    check( lapack_const_tally2( Magma_tally2DistUniform   )[0] == 'U' );
    check( lapack_const_tally2( Magma_tally2DistSymmetric )[0] == 'S' );
    check( lapack_const_tally2( Magma_tally2DistNormal    )[0] == 'N' );

    check( lapack_const_tally2( Magma_tally2HermGeev      )[0] == 'H' );
    check( lapack_const_tally2( Magma_tally2HermPoev      )[0] == 'P' );
    check( lapack_const_tally2( Magma_tally2NonsymPosv    )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2SymPosv       )[0] == 'S' );

    check( lapack_const_tally2( Magma_tally2NoPacking     )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2PackSubdiag   )[0] == 'U' );
    check( lapack_const_tally2( Magma_tally2PackSupdiag   )[0] == 'L' );
    check( lapack_const_tally2( Magma_tally2PackColumn    )[0] == 'C' );
    check( lapack_const_tally2( Magma_tally2PackRow       )[0] == 'R' );
    check( lapack_const_tally2( Magma_tally2PackLowerBand )[0] == 'B' );
    check( lapack_const_tally2( Magma_tally2PackUpeprBand )[0] == 'Q' );
    check( lapack_const_tally2( Magma_tally2PackAll       )[0] == 'Z' );

    check( lapack_const_tally2( Magma_tally2NoVec         )[0] == 'N' );
    check( lapack_const_tally2( Magma_tally2Vec           )[0] == 'V' );
    check( lapack_const_tally2( Magma_tally2IVec          )[0] == 'I' );
    check( lapack_const_tally2( Magma_tally2AllVec        )[0] == 'A' );
    check( lapack_const_tally2( Magma_tally2SomeVec       )[0] == 'S' );
    check( lapack_const_tally2( Magma_tally2OverwriteVec  )[0] == 'O' );

    check( lapack_const_tally2( Magma_tally2RangeAll      )[0] == 'A' );
    check( lapack_const_tally2( Magma_tally2RangeV        )[0] == 'V' );
    check( lapack_const_tally2( Magma_tally2RangeI        )[0] == 'I' );

    check( lapack_const_tally2( Magma_tally2Q             )[0] == 'Q' );
    check( lapack_const_tally2( Magma_tally2P             )[0] == 'P' );

    check( lapack_const_tally2( Magma_tally2Forward       )[0] == 'F' );
    check( lapack_const_tally2( Magma_tally2Backward      )[0] == 'B' );

    check( lapack_const_tally2( Magma_tally2Columnwise    )[0] == 'C' );
    check( lapack_const_tally2( Magma_tally2Rowwise       )[0] == 'R' );
    printf( "MAGMA_tally2  -> lapack_const_tally2          %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_const_tally2( Magma_tally2False         ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2True          ) == 'Y' );

    check( lapacke_const_tally2( Magma_tally2RowMajor      ) == 'R' );
    check( lapacke_const_tally2( Magma_tally2ColMajor      ) == 'C' );

    check( lapacke_const_tally2( Magma_tally2NoTrans       ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2Trans         ) == 'T' );
    check( lapacke_const_tally2( Magma_tally2ConjTrans     ) == 'C' );

    check( lapacke_const_tally2( Magma_tally2Upper         ) == 'U' );
    check( lapacke_const_tally2( Magma_tally2Lower         ) == 'L' );
    check( lapacke_const_tally2( Magma_tally2Full          ) == 'G' );

    check( lapacke_const_tally2( Magma_tally2NonUnit       ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2Unit          ) == 'U' );

    check( lapacke_const_tally2( Magma_tally2Left          ) == 'L' );
    check( lapacke_const_tally2( Magma_tally2Right         ) == 'R' );
    check( lapacke_const_tally2( Magma_tally2BothSides     ) == 'B' );

    check( lapacke_const_tally2( Magma_tally2OneNorm       ) == '1' );
    check( lapacke_const_tally2( Magma_tally2TwoNorm       ) == '2' );
    check( lapacke_const_tally2( Magma_tally2FrobeniusNorm ) == 'F' );
    check( lapacke_const_tally2( Magma_tally2InfNorm       ) == 'I' );
    check( lapacke_const_tally2( Magma_tally2MaxNorm       ) == 'M' );

    check( lapacke_const_tally2( Magma_tally2DistUniform   ) == 'U' );
    check( lapacke_const_tally2( Magma_tally2DistSymmetric ) == 'S' );
    check( lapacke_const_tally2( Magma_tally2DistNormal    ) == 'N' );

    check( lapacke_const_tally2( Magma_tally2HermGeev      ) == 'H' );
    check( lapacke_const_tally2( Magma_tally2HermPoev      ) == 'P' );
    check( lapacke_const_tally2( Magma_tally2NonsymPosv    ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2SymPosv       ) == 'S' );

    check( lapacke_const_tally2( Magma_tally2NoPacking     ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2PackSubdiag   ) == 'U' );
    check( lapacke_const_tally2( Magma_tally2PackSupdiag   ) == 'L' );
    check( lapacke_const_tally2( Magma_tally2PackColumn    ) == 'C' );
    check( lapacke_const_tally2( Magma_tally2PackRow       ) == 'R' );
    check( lapacke_const_tally2( Magma_tally2PackLowerBand ) == 'B' );
    check( lapacke_const_tally2( Magma_tally2PackUpeprBand ) == 'Q' );
    check( lapacke_const_tally2( Magma_tally2PackAll       ) == 'Z' );

    check( lapacke_const_tally2( Magma_tally2NoVec         ) == 'N' );
    check( lapacke_const_tally2( Magma_tally2Vec           ) == 'V' );
    check( lapacke_const_tally2( Magma_tally2IVec          ) == 'I' );
    check( lapacke_const_tally2( Magma_tally2AllVec        ) == 'A' );
    check( lapacke_const_tally2( Magma_tally2SomeVec       ) == 'S' );
    check( lapacke_const_tally2( Magma_tally2OverwriteVec  ) == 'O' );

    check( lapacke_const_tally2( Magma_tally2RangeAll      ) == 'A' );
    check( lapacke_const_tally2( Magma_tally2RangeV        ) == 'V' );
    check( lapacke_const_tally2( Magma_tally2RangeI        ) == 'I' );

    check( lapacke_const_tally2( Magma_tally2Q             ) == 'Q' );
    check( lapacke_const_tally2( Magma_tally2P             ) == 'P' );

    check( lapacke_const_tally2( Magma_tally2Forward       ) == 'F' );
    check( lapacke_const_tally2( Magma_tally2Backward      ) == 'B' );

    check( lapacke_const_tally2( Magma_tally2Columnwise    ) == 'C' );
    check( lapacke_const_tally2( Magma_tally2Rowwise       ) == 'R' );
    printf( "MAGMA_tally2  -> lapacke_const_tally2         %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( magma_tally2_bool_const('N') == Magma_tally2False );
    check( magma_tally2_bool_const('n') == Magma_tally2False );
    check( magma_tally2_bool_const('Y') == Magma_tally2True  );
    check( magma_tally2_bool_const('y') == Magma_tally2True  );

    check( magma_tally2_order_const( 'R' ) == Magma_tally2RowMajor  );
    check( magma_tally2_order_const( 'r' ) == Magma_tally2RowMajor  );
    check( magma_tally2_order_const( 'C' ) == Magma_tally2ColMajor  );
    check( magma_tally2_order_const( 'c' ) == Magma_tally2ColMajor  );

    check( magma_tally2_trans_const( 'N' ) == Magma_tally2NoTrans   );
    check( magma_tally2_trans_const( 'n' ) == Magma_tally2NoTrans   );
    check( magma_tally2_trans_const( 'T' ) == Magma_tally2Trans     );
    check( magma_tally2_trans_const( 't' ) == Magma_tally2Trans     );
    check( magma_tally2_trans_const( 'C' ) == Magma_tally2ConjTrans );
    check( magma_tally2_trans_const( 'c' ) == Magma_tally2ConjTrans );

    check( magma_tally2_uplo_const( 'U' ) == Magma_tally2Upper      );
    check( magma_tally2_uplo_const( 'u' ) == Magma_tally2Upper      );
    check( magma_tally2_uplo_const( 'L' ) == Magma_tally2Lower      );
    check( magma_tally2_uplo_const( 'l' ) == Magma_tally2Lower      );
    check( magma_tally2_uplo_const( 'A' ) == Magma_tally2Full       );  // anything else
    check( magma_tally2_uplo_const( 'a' ) == Magma_tally2Full       );
    check( magma_tally2_uplo_const( 'G' ) == Magma_tally2Full       );
    check( magma_tally2_uplo_const( 'g' ) == Magma_tally2Full       );
    check( magma_tally2_uplo_const( 'F' ) == Magma_tally2Full       );
    check( magma_tally2_uplo_const( 'f' ) == Magma_tally2Full       );

    check( magma_tally2_diag_const( 'N' ) == Magma_tally2NonUnit    );
    check( magma_tally2_diag_const( 'n' ) == Magma_tally2NonUnit    );
    check( magma_tally2_diag_const( 'U' ) == Magma_tally2Unit       );
    check( magma_tally2_diag_const( 'u' ) == Magma_tally2Unit       );

    check( magma_tally2_side_const( 'L' ) == Magma_tally2Left       );
    check( magma_tally2_side_const( 'l' ) == Magma_tally2Left       );
    check( magma_tally2_side_const( 'R' ) == Magma_tally2Right      );
    check( magma_tally2_side_const( 'r' ) == Magma_tally2Right      );

    check( magma_tally2_norm_const( 'O' ) == Magma_tally2OneNorm       );
    check( magma_tally2_norm_const( 'o' ) == Magma_tally2OneNorm       );
    check( magma_tally2_norm_const( '1' ) == Magma_tally2OneNorm       );
    check( magma_tally2_norm_const( '2' ) == Magma_tally2TwoNorm       );
    check( magma_tally2_norm_const( 'F' ) == Magma_tally2FrobeniusNorm );
    check( magma_tally2_norm_const( 'f' ) == Magma_tally2FrobeniusNorm );
    check( magma_tally2_norm_const( 'E' ) == Magma_tally2FrobeniusNorm );
    check( magma_tally2_norm_const( 'e' ) == Magma_tally2FrobeniusNorm );
    check( magma_tally2_norm_const( 'I' ) == Magma_tally2InfNorm       );
    check( magma_tally2_norm_const( 'i' ) == Magma_tally2InfNorm       );
    check( magma_tally2_norm_const( 'M' ) == Magma_tally2MaxNorm       );
    check( magma_tally2_norm_const( 'm' ) == Magma_tally2MaxNorm       );

    check( magma_tally2_dist_const( 'U' ) == Magma_tally2DistUniform   );
    check( magma_tally2_dist_const( 'u' ) == Magma_tally2DistUniform   );
    check( magma_tally2_dist_const( 'S' ) == Magma_tally2DistSymmetric );
    check( magma_tally2_dist_const( 's' ) == Magma_tally2DistSymmetric );
    check( magma_tally2_dist_const( 'N' ) == Magma_tally2DistNormal    );
    check( magma_tally2_dist_const( 'n' ) == Magma_tally2DistNormal    );

    //check( magma_tally2_xxxx_const( 'H' ) == Magma_tally2HermGeev      );
    //check( magma_tally2_xxxx_const( 'P' ) == Magma_tally2HermPoev      );
    //check( magma_tally2_xxxx_const( 'N' ) == Magma_tally2NonsymPosv    );
    //check( magma_tally2_xxxx_const( 'S' ) == Magma_tally2SymPosv       );

    check( magma_tally2_pack_const( 'N' ) == Magma_tally2NoPacking     );
    check( magma_tally2_pack_const( 'n' ) == Magma_tally2NoPacking     );
    check( magma_tally2_pack_const( 'U' ) == Magma_tally2PackSubdiag   );
    check( magma_tally2_pack_const( 'u' ) == Magma_tally2PackSubdiag   );
    check( magma_tally2_pack_const( 'L' ) == Magma_tally2PackSupdiag   );
    check( magma_tally2_pack_const( 'l' ) == Magma_tally2PackSupdiag   );
    check( magma_tally2_pack_const( 'C' ) == Magma_tally2PackColumn    );
    check( magma_tally2_pack_const( 'c' ) == Magma_tally2PackColumn    );
    check( magma_tally2_pack_const( 'R' ) == Magma_tally2PackRow       );
    check( magma_tally2_pack_const( 'r' ) == Magma_tally2PackRow       );
    check( magma_tally2_pack_const( 'B' ) == Magma_tally2PackLowerBand );
    check( magma_tally2_pack_const( 'b' ) == Magma_tally2PackLowerBand );
    check( magma_tally2_pack_const( 'Q' ) == Magma_tally2PackUpeprBand );
    check( magma_tally2_pack_const( 'q' ) == Magma_tally2PackUpeprBand );
    check( magma_tally2_pack_const( 'Z' ) == Magma_tally2PackAll       );
    check( magma_tally2_pack_const( 'z' ) == Magma_tally2PackAll       );

    check( magma_tally2_vec_const( 'N' )  == Magma_tally2NoVec         );
    check( magma_tally2_vec_const( 'n' )  == Magma_tally2NoVec         );
    check( magma_tally2_vec_const( 'V' )  == Magma_tally2Vec           );
    check( magma_tally2_vec_const( 'v' )  == Magma_tally2Vec           );
    check( magma_tally2_vec_const( 'I' )  == Magma_tally2IVec          );
    check( magma_tally2_vec_const( 'i' )  == Magma_tally2IVec          );
    check( magma_tally2_vec_const( 'A' )  == Magma_tally2AllVec        );
    check( magma_tally2_vec_const( 'a' )  == Magma_tally2AllVec        );
    check( magma_tally2_vec_const( 'S' )  == Magma_tally2SomeVec       );
    check( magma_tally2_vec_const( 's' )  == Magma_tally2SomeVec       );
    check( magma_tally2_vec_const( 'O' )  == Magma_tally2OverwriteVec  );
    check( magma_tally2_vec_const( 'o' )  == Magma_tally2OverwriteVec  );

    check( magma_tally2_range_const( 'A' )  == Magma_tally2RangeAll    );
    check( magma_tally2_range_const( 'a' )  == Magma_tally2RangeAll    );
    check( magma_tally2_range_const( 'V' )  == Magma_tally2RangeV      );
    check( magma_tally2_range_const( 'v' )  == Magma_tally2RangeV      );
    check( magma_tally2_range_const( 'I' )  == Magma_tally2RangeI      );
    check( magma_tally2_range_const( 'i' )  == Magma_tally2RangeI      );

    check( magma_tally2_vect_const( 'Q' )   == Magma_tally2Q           );
    check( magma_tally2_vect_const( 'q' )   == Magma_tally2Q           );
    check( magma_tally2_vect_const( 'P' )   == Magma_tally2P           );
    check( magma_tally2_vect_const( 'p' )   == Magma_tally2P           );

    check( magma_tally2_direct_const( 'F' ) == Magma_tally2Forward     );
    check( magma_tally2_direct_const( 'f' ) == Magma_tally2Forward     );
    check( magma_tally2_direct_const( 'B' ) == Magma_tally2Backward    );
    check( magma_tally2_direct_const( 'b' ) == Magma_tally2Backward    );

    check( magma_tally2_storev_const( 'C' ) == Magma_tally2Columnwise  );
    check( magma_tally2_storev_const( 'c' ) == Magma_tally2Columnwise  );
    check( magma_tally2_storev_const( 'R' ) == Magma_tally2Rowwise     );
    check( magma_tally2_storev_const( 'r' ) == Magma_tally2Rowwise     );
    printf( "LAPACK -> magma_tally2_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    #ifdef HAVE_clAmdBlas
    s = gStatus;
    check( amdblas_order_const( Magma_tally2RowMajor      ) == clAmdBlasRowMajor    );
    check( amdblas_order_const( Magma_tally2ColMajor      ) == clAmdBlasColumnMajor );

    check( amdblas_trans_const( Magma_tally2NoTrans       ) == clAmdBlasNoTrans     );
    check( amdblas_trans_const( Magma_tally2Trans         ) == clAmdBlasTrans       );
    check( amdblas_trans_const( Magma_tally2ConjTrans     ) == clAmdBlasConjTrans   );

    check( amdblas_uplo_const(  Magma_tally2Upper         ) == clAmdBlasUpper       );
    check( amdblas_uplo_const(  Magma_tally2Lower         ) == clAmdBlasLower       );

    check( amdblas_diag_const(  Magma_tally2NonUnit       ) == clAmdBlasNonUnit     );
    check( amdblas_diag_const(  Magma_tally2Unit          ) == clAmdBlasUnit        );

    check( amdblas_side_const(  Magma_tally2Left          ) == clAmdBlasLeft        );
    check( amdblas_side_const(  Magma_tally2Right         ) == clAmdBlasRight       );
    printf( "MAGMA_tally2  -> amdblas_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef CUBLAS_V2_H_
    s = gStatus;
    check( cublas_trans_const_tally2( Magma_tally2NoTrans       ) == CUBLAS_OP_N            );
    check( cublas_trans_const_tally2( Magma_tally2Trans         ) == CUBLAS_OP_T            );
    check( cublas_trans_const_tally2( Magma_tally2ConjTrans     ) == CUBLAS_OP_C            );

    check( cublas_uplo_const_tally2(  Magma_tally2Upper         ) == CUBLAS_FILL_MODE_UPPER );
    check( cublas_uplo_const_tally2(  Magma_tally2Lower         ) == CUBLAS_FILL_MODE_LOWER );

    check( cublas_diag_const_tally2(  Magma_tally2NonUnit       ) == CUBLAS_DIAG_NON_UNIT   );
    check( cublas_diag_const_tally2(  Magma_tally2Unit          ) == CUBLAS_DIAG_UNIT       );

    check( cublas_side_const_tally2(  Magma_tally2Left          ) == CUBLAS_SIDE_LEFT       );
    check( cublas_side_const_tally2(  Magma_tally2Right         ) == CUBLAS_SIDE_RIGHT      );
    printf( "MAGMA_tally2  -> cublas_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef HAVE_CBLAS
    s = gStatus;
    check( cblas_order_const( Magma_tally2RowMajor      ) == CblasRowMajor  );
    check( cblas_order_const( Magma_tally2ColMajor      ) == CblasColMajor  );

    check( cblas_trans_const( Magma_tally2NoTrans       ) == CblasNoTrans   );
    check( cblas_trans_const( Magma_tally2Trans         ) == CblasTrans     );
    check( cblas_trans_const( Magma_tally2ConjTrans     ) == CblasConjTrans );

    check( cblas_uplo_const(  Magma_tally2Upper         ) == CblasUpper     );
    check( cblas_uplo_const(  Magma_tally2Lower         ) == CblasLower     );

    check( cblas_diag_const(  Magma_tally2NonUnit       ) == CblasNonUnit   );
    check( cblas_diag_const(  Magma_tally2Unit          ) == CblasUnit      );

    check( cblas_side_const(  Magma_tally2Left          ) == CblasLeft      );
    check( cblas_side_const(  Magma_tally2Right         ) == CblasRight     );
    printf( "MAGMA_tally2  -> cblas_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));
    #endif

    return gStatus;
}
