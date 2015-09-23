/*
    -- MAGMA_minproduct (version 1.6.1) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date January 2015

       @author Mark Gates
*/
#include <stdio.h>

#include "testings.h"
#include "magma_minproduct.h"

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
    check( lapack_bool_const(   Magma_minproductFalse         )[0] == 'N' );
    check( lapack_bool_const(   Magma_minproductTrue          )[0] == 'Y' );

    check( lapack_order_const(  Magma_minproductRowMajor      )[0] == 'R' );
    check( lapack_order_const(  Magma_minproductColMajor      )[0] == 'C' );

    check( lapack_trans_const(  Magma_minproductNoTrans       )[0] == 'N' );
    check( lapack_trans_const(  Magma_minproductTrans         )[0] == 'T' );
    check( lapack_trans_const(  Magma_minproductConjTrans     )[0] == 'C' );

    check( lapack_uplo_const(   Magma_minproductUpper         )[0] == 'U' );
    check( lapack_uplo_const(   Magma_minproductLower         )[0] == 'L' );
    check( lapack_uplo_const(   Magma_minproductFull          )[0] == 'G' );

    check( lapack_diag_const(   Magma_minproductNonUnit       )[0] == 'N' );
    check( lapack_diag_const(   Magma_minproductUnit          )[0] == 'U' );

    check( lapack_side_const(   Magma_minproductLeft          )[0] == 'L' );
    check( lapack_side_const(   Magma_minproductRight         )[0] == 'R' );
    check( lapack_side_const(   Magma_minproductBothSides     )[0] == 'B' );

    check( lapack_norm_const(   Magma_minproductOneNorm       )[0] == '1' );
    check( lapack_norm_const(   Magma_minproductTwoNorm       )[0] == '2' );
    check( lapack_norm_const(   Magma_minproductFrobeniusNorm )[0] == 'F' );
    check( lapack_norm_const(   Magma_minproductInfNorm       )[0] == 'I' );
    check( lapack_norm_const(   Magma_minproductMaxNorm       )[0] == 'M' );

    check( lapack_dist_const(   Magma_minproductDistUniform   )[0] == 'U' );
    check( lapack_dist_const(   Magma_minproductDistSymmetric )[0] == 'S' );
    check( lapack_dist_const(   Magma_minproductDistNormal    )[0] == 'N' );

    check( lapack_sym_const(    Magma_minproductHermGeev      )[0] == 'H' );
    check( lapack_sym_const(    Magma_minproductHermPoev      )[0] == 'P' );
    check( lapack_sym_const(    Magma_minproductNonsymPosv    )[0] == 'N' );
    check( lapack_sym_const(    Magma_minproductSymPosv       )[0] == 'S' );

    check( lapack_pack_const(   Magma_minproductNoPacking     )[0] == 'N' );
    check( lapack_pack_const(   Magma_minproductPackSubdiag   )[0] == 'U' );
    check( lapack_pack_const(   Magma_minproductPackSupdiag   )[0] == 'L' );
    check( lapack_pack_const(   Magma_minproductPackColumn    )[0] == 'C' );
    check( lapack_pack_const(   Magma_minproductPackRow       )[0] == 'R' );
    check( lapack_pack_const(   Magma_minproductPackLowerBand )[0] == 'B' );
    check( lapack_pack_const(   Magma_minproductPackUpeprBand )[0] == 'Q' );
    check( lapack_pack_const(   Magma_minproductPackAll       )[0] == 'Z' );

    check( lapack_vec_const(    Magma_minproductNoVec         )[0] == 'N' );
    check( lapack_vec_const(    Magma_minproductVec           )[0] == 'V' );
    check( lapack_vec_const(    Magma_minproductIVec          )[0] == 'I' );
    check( lapack_vec_const(    Magma_minproductAllVec        )[0] == 'A' );
    check( lapack_vec_const(    Magma_minproductSomeVec       )[0] == 'S' );
    check( lapack_vec_const(    Magma_minproductOverwriteVec  )[0] == 'O' );

    check( lapack_range_const(  Magma_minproductRangeAll      )[0] == 'A' );
    check( lapack_range_const(  Magma_minproductRangeV        )[0] == 'V' );
    check( lapack_range_const(  Magma_minproductRangeI        )[0] == 'I' );

    check( lapack_vect_const(   Magma_minproductQ             )[0] == 'Q' );
    check( lapack_vect_const(   Magma_minproductP             )[0] == 'P' );

    check( lapack_direct_const( Magma_minproductForward       )[0] == 'F' );
    check( lapack_direct_const( Magma_minproductBackward      )[0] == 'B' );

    check( lapack_storev_const( Magma_minproductColumnwise    )[0] == 'C' );
    check( lapack_storev_const( Magma_minproductRowwise       )[0] == 'R' );
    printf( "MAGMA_minproduct  -> lapack_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_bool_const(   Magma_minproductFalse         ) == 'N' );
    check( lapacke_bool_const(   Magma_minproductTrue          ) == 'Y' );

    check( lapacke_order_const(  Magma_minproductRowMajor      ) == 'R' );
    check( lapacke_order_const(  Magma_minproductColMajor      ) == 'C' );

    check( lapacke_trans_const(  Magma_minproductNoTrans       ) == 'N' );
    check( lapacke_trans_const(  Magma_minproductTrans         ) == 'T' );
    check( lapacke_trans_const(  Magma_minproductConjTrans     ) == 'C' );

    check( lapacke_uplo_const(   Magma_minproductUpper         ) == 'U' );
    check( lapacke_uplo_const(   Magma_minproductLower         ) == 'L' );
    check( lapacke_uplo_const(   Magma_minproductFull          ) == 'G' );

    check( lapacke_diag_const(   Magma_minproductNonUnit       ) == 'N' );
    check( lapacke_diag_const(   Magma_minproductUnit          ) == 'U' );

    check( lapacke_side_const(   Magma_minproductLeft          ) == 'L' );
    check( lapacke_side_const(   Magma_minproductRight         ) == 'R' );
    check( lapacke_side_const(   Magma_minproductBothSides     ) == 'B' );

    check( lapacke_norm_const(   Magma_minproductOneNorm       ) == '1' );
    check( lapacke_norm_const(   Magma_minproductTwoNorm       ) == '2' );
    check( lapacke_norm_const(   Magma_minproductFrobeniusNorm ) == 'F' );
    check( lapacke_norm_const(   Magma_minproductInfNorm       ) == 'I' );
    check( lapacke_norm_const(   Magma_minproductMaxNorm       ) == 'M' );

    check( lapacke_dist_const(   Magma_minproductDistUniform   ) == 'U' );
    check( lapacke_dist_const(   Magma_minproductDistSymmetric ) == 'S' );
    check( lapacke_dist_const(   Magma_minproductDistNormal    ) == 'N' );

    check( lapacke_sym_const(    Magma_minproductHermGeev      ) == 'H' );
    check( lapacke_sym_const(    Magma_minproductHermPoev      ) == 'P' );
    check( lapacke_sym_const(    Magma_minproductNonsymPosv    ) == 'N' );
    check( lapacke_sym_const(    Magma_minproductSymPosv       ) == 'S' );

    check( lapacke_pack_const(   Magma_minproductNoPacking     ) == 'N' );
    check( lapacke_pack_const(   Magma_minproductPackSubdiag   ) == 'U' );
    check( lapacke_pack_const(   Magma_minproductPackSupdiag   ) == 'L' );
    check( lapacke_pack_const(   Magma_minproductPackColumn    ) == 'C' );
    check( lapacke_pack_const(   Magma_minproductPackRow       ) == 'R' );
    check( lapacke_pack_const(   Magma_minproductPackLowerBand ) == 'B' );
    check( lapacke_pack_const(   Magma_minproductPackUpeprBand ) == 'Q' );
    check( lapacke_pack_const(   Magma_minproductPackAll       ) == 'Z' );

    check( lapacke_vec_const(    Magma_minproductNoVec         ) == 'N' );
    check( lapacke_vec_const(    Magma_minproductVec           ) == 'V' );
    check( lapacke_vec_const(    Magma_minproductIVec          ) == 'I' );
    check( lapacke_vec_const(    Magma_minproductAllVec        ) == 'A' );
    check( lapacke_vec_const(    Magma_minproductSomeVec       ) == 'S' );
    check( lapacke_vec_const(    Magma_minproductOverwriteVec  ) == 'O' );

    check( lapacke_range_const(  Magma_minproductRangeAll      ) == 'A' );
    check( lapacke_range_const(  Magma_minproductRangeV        ) == 'V' );
    check( lapacke_range_const(  Magma_minproductRangeI        ) == 'I' );

    check( lapacke_vect_const(   Magma_minproductQ             ) == 'Q' );
    check( lapacke_vect_const(   Magma_minproductP             ) == 'P' );

    check( lapacke_direct_const( Magma_minproductForward       ) == 'F' );
    check( lapacke_direct_const( Magma_minproductBackward      ) == 'B' );

    check( lapacke_storev_const( Magma_minproductColumnwise    ) == 'C' );
    check( lapacke_storev_const( Magma_minproductRowwise       ) == 'R' );
    printf( "MAGMA_minproduct  -> lapacke_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapack_const( Magma_minproductFalse         )[0] == 'N' );
    check( lapack_const( Magma_minproductTrue          )[0] == 'Y' );

    check( lapack_const( Magma_minproductRowMajor      )[0] == 'R' );
    check( lapack_const( Magma_minproductColMajor      )[0] == 'C' );

    check( lapack_const( Magma_minproductNoTrans       )[0] == 'N' );
    check( lapack_const( Magma_minproductTrans         )[0] == 'T' );
    check( lapack_const( Magma_minproductConjTrans     )[0] == 'C' );

    check( lapack_const( Magma_minproductUpper         )[0] == 'U' );
    check( lapack_const( Magma_minproductLower         )[0] == 'L' );
    check( lapack_const( Magma_minproductFull          )[0] == 'G' );

    check( lapack_const( Magma_minproductNonUnit       )[0] == 'N' );
    check( lapack_const( Magma_minproductUnit          )[0] == 'U' );

    check( lapack_const( Magma_minproductLeft          )[0] == 'L' );
    check( lapack_const( Magma_minproductRight         )[0] == 'R' );
    check( lapack_const( Magma_minproductBothSides     )[0] == 'B' );

    check( lapack_const( Magma_minproductOneNorm       )[0] == '1' );
    check( lapack_const( Magma_minproductTwoNorm       )[0] == '2' );
    check( lapack_const( Magma_minproductFrobeniusNorm )[0] == 'F' );
    check( lapack_const( Magma_minproductInfNorm       )[0] == 'I' );
    check( lapack_const( Magma_minproductMaxNorm       )[0] == 'M' );

    check( lapack_const( Magma_minproductDistUniform   )[0] == 'U' );
    check( lapack_const( Magma_minproductDistSymmetric )[0] == 'S' );
    check( lapack_const( Magma_minproductDistNormal    )[0] == 'N' );

    check( lapack_const( Magma_minproductHermGeev      )[0] == 'H' );
    check( lapack_const( Magma_minproductHermPoev      )[0] == 'P' );
    check( lapack_const( Magma_minproductNonsymPosv    )[0] == 'N' );
    check( lapack_const( Magma_minproductSymPosv       )[0] == 'S' );

    check( lapack_const( Magma_minproductNoPacking     )[0] == 'N' );
    check( lapack_const( Magma_minproductPackSubdiag   )[0] == 'U' );
    check( lapack_const( Magma_minproductPackSupdiag   )[0] == 'L' );
    check( lapack_const( Magma_minproductPackColumn    )[0] == 'C' );
    check( lapack_const( Magma_minproductPackRow       )[0] == 'R' );
    check( lapack_const( Magma_minproductPackLowerBand )[0] == 'B' );
    check( lapack_const( Magma_minproductPackUpeprBand )[0] == 'Q' );
    check( lapack_const( Magma_minproductPackAll       )[0] == 'Z' );

    check( lapack_const( Magma_minproductNoVec         )[0] == 'N' );
    check( lapack_const( Magma_minproductVec           )[0] == 'V' );
    check( lapack_const( Magma_minproductIVec          )[0] == 'I' );
    check( lapack_const( Magma_minproductAllVec        )[0] == 'A' );
    check( lapack_const( Magma_minproductSomeVec       )[0] == 'S' );
    check( lapack_const( Magma_minproductOverwriteVec  )[0] == 'O' );

    check( lapack_const( Magma_minproductRangeAll      )[0] == 'A' );
    check( lapack_const( Magma_minproductRangeV        )[0] == 'V' );
    check( lapack_const( Magma_minproductRangeI        )[0] == 'I' );

    check( lapack_const( Magma_minproductQ             )[0] == 'Q' );
    check( lapack_const( Magma_minproductP             )[0] == 'P' );

    check( lapack_const( Magma_minproductForward       )[0] == 'F' );
    check( lapack_const( Magma_minproductBackward      )[0] == 'B' );

    check( lapack_const( Magma_minproductColumnwise    )[0] == 'C' );
    check( lapack_const( Magma_minproductRowwise       )[0] == 'R' );
    printf( "MAGMA_minproduct  -> lapack_const          %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( lapacke_const( Magma_minproductFalse         ) == 'N' );
    check( lapacke_const( Magma_minproductTrue          ) == 'Y' );

    check( lapacke_const( Magma_minproductRowMajor      ) == 'R' );
    check( lapacke_const( Magma_minproductColMajor      ) == 'C' );

    check( lapacke_const( Magma_minproductNoTrans       ) == 'N' );
    check( lapacke_const( Magma_minproductTrans         ) == 'T' );
    check( lapacke_const( Magma_minproductConjTrans     ) == 'C' );

    check( lapacke_const( Magma_minproductUpper         ) == 'U' );
    check( lapacke_const( Magma_minproductLower         ) == 'L' );
    check( lapacke_const( Magma_minproductFull          ) == 'G' );

    check( lapacke_const( Magma_minproductNonUnit       ) == 'N' );
    check( lapacke_const( Magma_minproductUnit          ) == 'U' );

    check( lapacke_const( Magma_minproductLeft          ) == 'L' );
    check( lapacke_const( Magma_minproductRight         ) == 'R' );
    check( lapacke_const( Magma_minproductBothSides     ) == 'B' );

    check( lapacke_const( Magma_minproductOneNorm       ) == '1' );
    check( lapacke_const( Magma_minproductTwoNorm       ) == '2' );
    check( lapacke_const( Magma_minproductFrobeniusNorm ) == 'F' );
    check( lapacke_const( Magma_minproductInfNorm       ) == 'I' );
    check( lapacke_const( Magma_minproductMaxNorm       ) == 'M' );

    check( lapacke_const( Magma_minproductDistUniform   ) == 'U' );
    check( lapacke_const( Magma_minproductDistSymmetric ) == 'S' );
    check( lapacke_const( Magma_minproductDistNormal    ) == 'N' );

    check( lapacke_const( Magma_minproductHermGeev      ) == 'H' );
    check( lapacke_const( Magma_minproductHermPoev      ) == 'P' );
    check( lapacke_const( Magma_minproductNonsymPosv    ) == 'N' );
    check( lapacke_const( Magma_minproductSymPosv       ) == 'S' );

    check( lapacke_const( Magma_minproductNoPacking     ) == 'N' );
    check( lapacke_const( Magma_minproductPackSubdiag   ) == 'U' );
    check( lapacke_const( Magma_minproductPackSupdiag   ) == 'L' );
    check( lapacke_const( Magma_minproductPackColumn    ) == 'C' );
    check( lapacke_const( Magma_minproductPackRow       ) == 'R' );
    check( lapacke_const( Magma_minproductPackLowerBand ) == 'B' );
    check( lapacke_const( Magma_minproductPackUpeprBand ) == 'Q' );
    check( lapacke_const( Magma_minproductPackAll       ) == 'Z' );

    check( lapacke_const( Magma_minproductNoVec         ) == 'N' );
    check( lapacke_const( Magma_minproductVec           ) == 'V' );
    check( lapacke_const( Magma_minproductIVec          ) == 'I' );
    check( lapacke_const( Magma_minproductAllVec        ) == 'A' );
    check( lapacke_const( Magma_minproductSomeVec       ) == 'S' );
    check( lapacke_const( Magma_minproductOverwriteVec  ) == 'O' );

    check( lapacke_const( Magma_minproductRangeAll      ) == 'A' );
    check( lapacke_const( Magma_minproductRangeV        ) == 'V' );
    check( lapacke_const( Magma_minproductRangeI        ) == 'I' );

    check( lapacke_const( Magma_minproductQ             ) == 'Q' );
    check( lapacke_const( Magma_minproductP             ) == 'P' );

    check( lapacke_const( Magma_minproductForward       ) == 'F' );
    check( lapacke_const( Magma_minproductBackward      ) == 'B' );

    check( lapacke_const( Magma_minproductColumnwise    ) == 'C' );
    check( lapacke_const( Magma_minproductRowwise       ) == 'R' );
    printf( "MAGMA_minproduct  -> lapacke_const         %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    s = gStatus;
    check( magma_minproduct_bool_const('N') == Magma_minproductFalse );
    check( magma_minproduct_bool_const('n') == Magma_minproductFalse );
    check( magma_minproduct_bool_const('Y') == Magma_minproductTrue  );
    check( magma_minproduct_bool_const('y') == Magma_minproductTrue  );

    check( magma_minproduct_order_const( 'R' ) == Magma_minproductRowMajor  );
    check( magma_minproduct_order_const( 'r' ) == Magma_minproductRowMajor  );
    check( magma_minproduct_order_const( 'C' ) == Magma_minproductColMajor  );
    check( magma_minproduct_order_const( 'c' ) == Magma_minproductColMajor  );

    check( magma_minproduct_trans_const( 'N' ) == Magma_minproductNoTrans   );
    check( magma_minproduct_trans_const( 'n' ) == Magma_minproductNoTrans   );
    check( magma_minproduct_trans_const( 'T' ) == Magma_minproductTrans     );
    check( magma_minproduct_trans_const( 't' ) == Magma_minproductTrans     );
    check( magma_minproduct_trans_const( 'C' ) == Magma_minproductConjTrans );
    check( magma_minproduct_trans_const( 'c' ) == Magma_minproductConjTrans );

    check( magma_minproduct_uplo_const( 'U' ) == Magma_minproductUpper      );
    check( magma_minproduct_uplo_const( 'u' ) == Magma_minproductUpper      );
    check( magma_minproduct_uplo_const( 'L' ) == Magma_minproductLower      );
    check( magma_minproduct_uplo_const( 'l' ) == Magma_minproductLower      );
    check( magma_minproduct_uplo_const( 'A' ) == Magma_minproductFull       );  // anything else
    check( magma_minproduct_uplo_const( 'a' ) == Magma_minproductFull       );
    check( magma_minproduct_uplo_const( 'G' ) == Magma_minproductFull       );
    check( magma_minproduct_uplo_const( 'g' ) == Magma_minproductFull       );
    check( magma_minproduct_uplo_const( 'F' ) == Magma_minproductFull       );
    check( magma_minproduct_uplo_const( 'f' ) == Magma_minproductFull       );

    check( magma_minproduct_diag_const( 'N' ) == Magma_minproductNonUnit    );
    check( magma_minproduct_diag_const( 'n' ) == Magma_minproductNonUnit    );
    check( magma_minproduct_diag_const( 'U' ) == Magma_minproductUnit       );
    check( magma_minproduct_diag_const( 'u' ) == Magma_minproductUnit       );

    check( magma_minproduct_side_const( 'L' ) == Magma_minproductLeft       );
    check( magma_minproduct_side_const( 'l' ) == Magma_minproductLeft       );
    check( magma_minproduct_side_const( 'R' ) == Magma_minproductRight      );
    check( magma_minproduct_side_const( 'r' ) == Magma_minproductRight      );

    check( magma_minproduct_norm_const( 'O' ) == Magma_minproductOneNorm       );
    check( magma_minproduct_norm_const( 'o' ) == Magma_minproductOneNorm       );
    check( magma_minproduct_norm_const( '1' ) == Magma_minproductOneNorm       );
    check( magma_minproduct_norm_const( '2' ) == Magma_minproductTwoNorm       );
    check( magma_minproduct_norm_const( 'F' ) == Magma_minproductFrobeniusNorm );
    check( magma_minproduct_norm_const( 'f' ) == Magma_minproductFrobeniusNorm );
    check( magma_minproduct_norm_const( 'E' ) == Magma_minproductFrobeniusNorm );
    check( magma_minproduct_norm_const( 'e' ) == Magma_minproductFrobeniusNorm );
    check( magma_minproduct_norm_const( 'I' ) == Magma_minproductInfNorm       );
    check( magma_minproduct_norm_const( 'i' ) == Magma_minproductInfNorm       );
    check( magma_minproduct_norm_const( 'M' ) == Magma_minproductMaxNorm       );
    check( magma_minproduct_norm_const( 'm' ) == Magma_minproductMaxNorm       );

    check( magma_minproduct_dist_const( 'U' ) == Magma_minproductDistUniform   );
    check( magma_minproduct_dist_const( 'u' ) == Magma_minproductDistUniform   );
    check( magma_minproduct_dist_const( 'S' ) == Magma_minproductDistSymmetric );
    check( magma_minproduct_dist_const( 's' ) == Magma_minproductDistSymmetric );
    check( magma_minproduct_dist_const( 'N' ) == Magma_minproductDistNormal    );
    check( magma_minproduct_dist_const( 'n' ) == Magma_minproductDistNormal    );

    //check( magma_minproduct_xxxx_const( 'H' ) == Magma_minproductHermGeev      );
    //check( magma_minproduct_xxxx_const( 'P' ) == Magma_minproductHermPoev      );
    //check( magma_minproduct_xxxx_const( 'N' ) == Magma_minproductNonsymPosv    );
    //check( magma_minproduct_xxxx_const( 'S' ) == Magma_minproductSymPosv       );

    check( magma_minproduct_pack_const( 'N' ) == Magma_minproductNoPacking     );
    check( magma_minproduct_pack_const( 'n' ) == Magma_minproductNoPacking     );
    check( magma_minproduct_pack_const( 'U' ) == Magma_minproductPackSubdiag   );
    check( magma_minproduct_pack_const( 'u' ) == Magma_minproductPackSubdiag   );
    check( magma_minproduct_pack_const( 'L' ) == Magma_minproductPackSupdiag   );
    check( magma_minproduct_pack_const( 'l' ) == Magma_minproductPackSupdiag   );
    check( magma_minproduct_pack_const( 'C' ) == Magma_minproductPackColumn    );
    check( magma_minproduct_pack_const( 'c' ) == Magma_minproductPackColumn    );
    check( magma_minproduct_pack_const( 'R' ) == Magma_minproductPackRow       );
    check( magma_minproduct_pack_const( 'r' ) == Magma_minproductPackRow       );
    check( magma_minproduct_pack_const( 'B' ) == Magma_minproductPackLowerBand );
    check( magma_minproduct_pack_const( 'b' ) == Magma_minproductPackLowerBand );
    check( magma_minproduct_pack_const( 'Q' ) == Magma_minproductPackUpeprBand );
    check( magma_minproduct_pack_const( 'q' ) == Magma_minproductPackUpeprBand );
    check( magma_minproduct_pack_const( 'Z' ) == Magma_minproductPackAll       );
    check( magma_minproduct_pack_const( 'z' ) == Magma_minproductPackAll       );

    check( magma_minproduct_vec_const( 'N' )  == Magma_minproductNoVec         );
    check( magma_minproduct_vec_const( 'n' )  == Magma_minproductNoVec         );
    check( magma_minproduct_vec_const( 'V' )  == Magma_minproductVec           );
    check( magma_minproduct_vec_const( 'v' )  == Magma_minproductVec           );
    check( magma_minproduct_vec_const( 'I' )  == Magma_minproductIVec          );
    check( magma_minproduct_vec_const( 'i' )  == Magma_minproductIVec          );
    check( magma_minproduct_vec_const( 'A' )  == Magma_minproductAllVec        );
    check( magma_minproduct_vec_const( 'a' )  == Magma_minproductAllVec        );
    check( magma_minproduct_vec_const( 'S' )  == Magma_minproductSomeVec       );
    check( magma_minproduct_vec_const( 's' )  == Magma_minproductSomeVec       );
    check( magma_minproduct_vec_const( 'O' )  == Magma_minproductOverwriteVec  );
    check( magma_minproduct_vec_const( 'o' )  == Magma_minproductOverwriteVec  );

    check( magma_minproduct_range_const( 'A' )  == Magma_minproductRangeAll    );
    check( magma_minproduct_range_const( 'a' )  == Magma_minproductRangeAll    );
    check( magma_minproduct_range_const( 'V' )  == Magma_minproductRangeV      );
    check( magma_minproduct_range_const( 'v' )  == Magma_minproductRangeV      );
    check( magma_minproduct_range_const( 'I' )  == Magma_minproductRangeI      );
    check( magma_minproduct_range_const( 'i' )  == Magma_minproductRangeI      );

    check( magma_minproduct_vect_const( 'Q' )   == Magma_minproductQ           );
    check( magma_minproduct_vect_const( 'q' )   == Magma_minproductQ           );
    check( magma_minproduct_vect_const( 'P' )   == Magma_minproductP           );
    check( magma_minproduct_vect_const( 'p' )   == Magma_minproductP           );

    check( magma_minproduct_direct_const( 'F' ) == Magma_minproductForward     );
    check( magma_minproduct_direct_const( 'f' ) == Magma_minproductForward     );
    check( magma_minproduct_direct_const( 'B' ) == Magma_minproductBackward    );
    check( magma_minproduct_direct_const( 'b' ) == Magma_minproductBackward    );

    check( magma_minproduct_storev_const( 'C' ) == Magma_minproductColumnwise  );
    check( magma_minproduct_storev_const( 'c' ) == Magma_minproductColumnwise  );
    check( magma_minproduct_storev_const( 'R' ) == Magma_minproductRowwise     );
    check( magma_minproduct_storev_const( 'r' ) == Magma_minproductRowwise     );
    printf( "LAPACK -> magma_minproduct_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));


    // ------------------------------------------------------------
    #ifdef HAVE_clAmdBlas
    s = gStatus;
    check( amdblas_order_const( Magma_minproductRowMajor      ) == clAmdBlasRowMajor    );
    check( amdblas_order_const( Magma_minproductColMajor      ) == clAmdBlasColumnMajor );

    check( amdblas_trans_const( Magma_minproductNoTrans       ) == clAmdBlasNoTrans     );
    check( amdblas_trans_const( Magma_minproductTrans         ) == clAmdBlasTrans       );
    check( amdblas_trans_const( Magma_minproductConjTrans     ) == clAmdBlasConjTrans   );

    check( amdblas_uplo_const(  Magma_minproductUpper         ) == clAmdBlasUpper       );
    check( amdblas_uplo_const(  Magma_minproductLower         ) == clAmdBlasLower       );

    check( amdblas_diag_const(  Magma_minproductNonUnit       ) == clAmdBlasNonUnit     );
    check( amdblas_diag_const(  Magma_minproductUnit          ) == clAmdBlasUnit        );

    check( amdblas_side_const(  Magma_minproductLeft          ) == clAmdBlasLeft        );
    check( amdblas_side_const(  Magma_minproductRight         ) == clAmdBlasRight       );
    printf( "MAGMA_minproduct  -> amdblas_xxxxx_const   %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef CUBLAS_V2_H_
    s = gStatus;
    check( cublas_trans_const( Magma_minproductNoTrans       ) == CUBLAS_OP_N            );
    check( cublas_trans_const( Magma_minproductTrans         ) == CUBLAS_OP_T            );
    check( cublas_trans_const( Magma_minproductConjTrans     ) == CUBLAS_OP_C            );

    check( cublas_uplo_const(  Magma_minproductUpper         ) == CUBLAS_FILL_MODE_UPPER );
    check( cublas_uplo_const(  Magma_minproductLower         ) == CUBLAS_FILL_MODE_LOWER );

    check( cublas_diag_const(  Magma_minproductNonUnit       ) == CUBLAS_DIAG_NON_UNIT   );
    check( cublas_diag_const(  Magma_minproductUnit          ) == CUBLAS_DIAG_UNIT       );

    check( cublas_side_const(  Magma_minproductLeft          ) == CUBLAS_SIDE_LEFT       );
    check( cublas_side_const(  Magma_minproductRight         ) == CUBLAS_SIDE_RIGHT      );
    printf( "MAGMA_minproduct  -> cublas_xxxxx_const    %s\n", (s == gStatus ? "ok" : "failed"));
    #endif


    // ------------------------------------------------------------
    #ifdef HAVE_CBLAS
    s = gStatus;
    check( cblas_order_const( Magma_minproductRowMajor      ) == CblasRowMajor  );
    check( cblas_order_const( Magma_minproductColMajor      ) == CblasColMajor  );

    check( cblas_trans_const( Magma_minproductNoTrans       ) == CblasNoTrans   );
    check( cblas_trans_const( Magma_minproductTrans         ) == CblasTrans     );
    check( cblas_trans_const( Magma_minproductConjTrans     ) == CblasConjTrans );

    check( cblas_uplo_const(  Magma_minproductUpper         ) == CblasUpper     );
    check( cblas_uplo_const(  Magma_minproductLower         ) == CblasLower     );

    check( cblas_diag_const(  Magma_minproductNonUnit       ) == CblasNonUnit   );
    check( cblas_diag_const(  Magma_minproductUnit          ) == CblasUnit      );

    check( cblas_side_const(  Magma_minproductLeft          ) == CblasLeft      );
    check( cblas_side_const(  Magma_minproductRight         ) == CblasRight     );
    printf( "MAGMA_minproduct  -> cblas_xxxxx_const     %s\n", (s == gStatus ? "ok" : "failed"));
    #endif

    return gStatus;
}
