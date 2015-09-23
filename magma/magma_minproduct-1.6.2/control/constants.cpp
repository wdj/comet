#include <assert.h>
#include <stdio.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#endif

#include "magma_minproduct_types.h"

// ----------------------------------------
// Convert LAPACK character constants to MAGMA_minproduct constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
// These functions and cases are in the same order as the constants are
// declared in magma_minproduct_types.h

extern "C"
magma_minproduct_bool_t   magma_minproduct_bool_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_minproductFalse;
        case 'Y': case 'y': return Magma_minproductTrue;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductFalse;
    }
}

extern "C"
magma_minproduct_order_t  magma_minproduct_order_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'R': case 'r': return Magma_minproductRowMajor;
        case 'C': case 'c': return Magma_minproductColMajor;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductRowMajor;
    }
}

extern "C"
magma_minproduct_trans_t  magma_minproduct_trans_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_minproductNoTrans;
        case 'T': case 't': return Magma_minproductTrans;
        case 'C': case 'c': return Magma_minproductConjTrans;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductNoTrans;
    }
}

extern "C"
magma_minproduct_uplo_t   magma_minproduct_uplo_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_minproductUpper;
        case 'L': case 'l': return Magma_minproductLower;
        default:            return Magma_minproductFull;        // see laset
    }
}

extern "C"
magma_minproduct_diag_t   magma_minproduct_diag_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_minproductNonUnit;
        case 'U': case 'u': return Magma_minproductUnit;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductNonUnit;
    }
}

extern "C"
magma_minproduct_side_t   magma_minproduct_side_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'L': case 'l': return Magma_minproductLeft;
        case 'R': case 'r': return Magma_minproductRight;
        case 'B': case 'b': return Magma_minproductBothSides;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductLeft;
    }
}

extern "C"
magma_minproduct_norm_t   magma_minproduct_norm_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'O': case 'o': case '1': return Magma_minproductOneNorm;
        case '2':           return Magma_minproductTwoNorm;
        case 'F': case 'f': case 'E': case 'e': return Magma_minproductFrobeniusNorm;
        case 'I': case 'i': return Magma_minproductInfNorm;
        case 'M': case 'm': return Magma_minproductMaxNorm;
        // Magma_minproductRealOneNorm
        // Magma_minproductRealInfNorm
        // Magma_minproductRealMaxNorm
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductOneNorm;
    }
}

extern "C"
magma_minproduct_dist_t   magma_minproduct_dist_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_minproductDistUniform;
        case 'S': case 's': return Magma_minproductDistSymmetric;
        case 'N': case 'n': return Magma_minproductDistNormal;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductDistUniform;
    }
}

extern "C"
magma_minproduct_sym_t    magma_minproduct_sym_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'H': case 'h': return Magma_minproductHermGeev;
        case 'P': case 'p': return Magma_minproductHermPoev;
        case 'N': case 'n': return Magma_minproductNonsymPosv;
        case 'S': case 's': return Magma_minproductSymPosv;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductHermGeev;
    }
}

extern "C"
magma_minproduct_pack_t   magma_minproduct_pack_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_minproductNoPacking;
        case 'U': case 'u': return Magma_minproductPackSubdiag;
        case 'L': case 'l': return Magma_minproductPackSupdiag;
        case 'C': case 'c': return Magma_minproductPackColumn;
        case 'R': case 'r': return Magma_minproductPackRow;
        case 'B': case 'b': return Magma_minproductPackLowerBand;
        case 'Q': case 'q': return Magma_minproductPackUpeprBand;
        case 'Z': case 'z': return Magma_minproductPackAll;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductNoPacking;
    }
}

extern "C"
magma_minproduct_vec_t    magma_minproduct_vec_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_minproductNoVec;
        case 'V': case 'v': return Magma_minproductVec;
        case 'I': case 'i': return Magma_minproductIVec;
        case 'A': case 'a': return Magma_minproductAllVec;
        case 'S': case 's': return Magma_minproductSomeVec;
        case 'O': case 'o': return Magma_minproductOverwriteVec;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductNoVec;
    }
}

extern "C"
magma_minproduct_range_t  magma_minproduct_range_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'A': case 'a': return Magma_minproductRangeAll;
        case 'V': case 'v': return Magma_minproductRangeV;
        case 'I': case 'i': return Magma_minproductRangeI;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductRangeAll;
    }
}

extern "C"
magma_minproduct_vect_t magma_minproduct_vect_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'Q': case 'q': return Magma_minproductQ;
        case 'P': case 'p': return Magma_minproductP;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductQ;
    }
}

extern "C"
magma_minproduct_direct_t magma_minproduct_direct_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'F': case 'f': return Magma_minproductForward;
        case 'B': case 'b': return Magma_minproductBackward;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductForward;
    }
}

extern "C"
magma_minproduct_storev_t magma_minproduct_storev_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'C': case 'c': return Magma_minproductColumnwise;
        case 'R': case 'r': return Magma_minproductRowwise;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_minproductColumnwise;
    }
}


// ----------------------------------------
// Convert MAGMA_minproduct constants to LAPACK constants.

const char *magma_minproduct2lapack_constants[] =
{
    "No",                                    //  0: Magma_minproductFalse
    "Yes",                                   //  1: Magma_minproductTrue (zlatrs)
    "", "", "", "", "", "", "", "", "",      //  2-10
    "", "", "", "", "", "", "", "", "", "",  // 11-20
    "", "", "", "", "", "", "", "", "", "",  // 21-30
    "", "", "", "", "", "", "", "", "", "",  // 31-40
    "", "", "", "", "", "", "", "", "", "",  // 41-50
    "", "", "", "", "", "", "", "", "", "",  // 51-60
    "", "", "", "", "", "", "", "", "", "",  // 61-70
    "", "", "", "", "", "", "", "", "", "",  // 71-80
    "", "", "", "", "", "", "", "", "", "",  // 81-90
    "", "", "", "", "", "", "", "", "", "",  // 91-100
    "Row",                                   // 101: Magma_minproductRowMajor
    "Column",                                // 102: Magma_minproductColMajor
    "", "", "", "", "", "", "", "",          // 103-110
    "No transpose",                          // 111: Magma_minproductNoTrans
    "Transpose",                             // 112: Magma_minproductTrans
    "Conjugate transpose",                   // 113: Magma_minproductConjTrans
    "", "", "", "", "", "", "",              // 114-120
    "Upper",                                 // 121: Magma_minproductUpper
    "Lower",                                 // 122: Magma_minproductLower
    "GFull",                                 // 123: Magma_minproductFull; see lascl for "G"
    "", "", "", "", "", "", "",              // 124-130
    "Non-unit",                              // 131: Magma_minproductNonUnit
    "Unit",                                  // 132: Magma_minproductUnit
    "", "", "", "", "", "", "", "",          // 133-140
    "Left",                                  // 141: Magma_minproductLeft
    "Right",                                 // 142: Magma_minproductRight
    "Both",                                  // 143: Magma_minproductBothSides (dtrevc)
    "", "", "", "", "", "", "",              // 144-150
    "", "", "", "", "", "", "", "", "", "",  // 151-160
    "", "", "", "", "", "", "", "", "", "",  // 161-170
    "1 norm",                                // 171: Magma_minproductOneNorm
    "",                                      // 172: Magma_minproductRealOneNorm
    "2 norm",                                // 173: Magma_minproductTwoNorm
    "Frobenius norm",                        // 174: Magma_minproductFrobeniusNorm
    "Infinity norm",                         // 175: Magma_minproductInfNorm
    "",                                      // 176: Magma_minproductRealInfNorm
    "Maximum norm",                          // 177: Magma_minproductMaxNorm
    "",                                      // 178: Magma_minproductRealMaxNorm
    "", "",                                  // 179-180
    "", "", "", "", "", "", "", "", "", "",  // 181-190
    "", "", "", "", "", "", "", "", "", "",  // 191-200
    "Uniform",                               // 201: Magma_minproductDistUniform
    "Symmetric",                             // 202: Magma_minproductDistSymmetric
    "Normal",                                // 203: Magma_minproductDistNormal
    "", "", "", "", "", "", "",              // 204-210
    "", "", "", "", "", "", "", "", "", "",  // 211-220
    "", "", "", "", "", "", "", "", "", "",  // 221-230
    "", "", "", "", "", "", "", "", "", "",  // 231-240
    "Hermitian",                             // 241 Magma_minproductHermGeev
    "Positive ev Hermitian",                 // 242 Magma_minproductHermPoev
    "NonSymmetric pos sv",                   // 243 Magma_minproductNonsymPosv
    "Symmetric pos sv",                      // 244 Magma_minproductSymPosv
    "", "", "", "", "", "",                  // 245-250
    "", "", "", "", "", "", "", "", "", "",  // 251-260
    "", "", "", "", "", "", "", "", "", "",  // 261-270
    "", "", "", "", "", "", "", "", "", "",  // 271-280
    "", "", "", "", "", "", "", "", "", "",  // 281-290
    "No Packing",                            // 291 Magma_minproductNoPacking
    "U zero out subdiag",                    // 292 Magma_minproductPackSubdiag
    "L zero out superdiag",                  // 293 Magma_minproductPackSupdiag
    "C",                                     // 294 Magma_minproductPackColumn
    "R",                                     // 295 Magma_minproductPackRow
    "B",                                     // 296 Magma_minproductPackLowerBand
    "Q",                                     // 297 Magma_minproductPackUpeprBand
    "Z",                                     // 298 Magma_minproductPackAll
    "", "",                                  // 299-300
    "No vectors",                            // 301 Magma_minproductNoVec
    "Vectors needed",                        // 302 Magma_minproductVec
    "I",                                     // 303 Magma_minproductIVec
    "All",                                   // 304 Magma_minproductAllVec
    "Some",                                  // 305 Magma_minproductSomeVec
    "Overwrite",                             // 306 Magma_minproductOverwriteVec
    "", "", "", "",                          // 307-310
    "All",                                   // 311 Magma_minproductRangeAll
    "V",                                     // 312 Magma_minproductRangeV
    "I",                                     // 313 Magma_minproductRangeI
    "", "", "", "", "", "", "",              // 314-320
    "",                                      // 321
    "Q",                                     // 322
    "P",                                     // 323
    "", "", "", "", "", "", "",              // 324-330
    "", "", "", "", "", "", "", "", "", "",  // 331-340
    "", "", "", "", "", "", "", "", "", "",  // 341-350
    "", "", "", "", "", "", "", "", "", "",  // 351-360
    "", "", "", "", "", "", "", "", "", "",  // 361-370
    "", "", "", "", "", "", "", "", "", "",  // 371-380
    "", "", "", "", "", "", "", "", "", "",  // 381-390
    "Forward",                               // 391: Magma_minproductForward
    "Backward",                              // 392: Magma_minproductBackward
    "", "", "", "", "", "", "", "",          // 393-400
    "Columnwise",                            // 401: Magma_minproductColumnwise
    "Rowwise",                               // 402: Magma_minproductRowwise
    "", "", "", "", "", "", "", ""           // 403-410
    // Remember to add a comma!
};

extern "C"
const char* lapack_const( int magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproduct2lapack_Min );
    assert( magma_minproduct_const <= Magma_minproduct2lapack_Max );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_bool_const( magma_minproduct_bool_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductFalse );
    assert( magma_minproduct_const <= Magma_minproductTrue  );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_order_const( magma_minproduct_order_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductRowMajor );
    assert( magma_minproduct_const <= Magma_minproductColMajor );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_trans_const( magma_minproduct_trans_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoTrans   );
    assert( magma_minproduct_const <= Magma_minproductConjTrans );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_uplo_const ( magma_minproduct_uplo_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductUpper );
    assert( magma_minproduct_const <= Magma_minproductFull  );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_diag_const ( magma_minproduct_diag_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNonUnit );
    assert( magma_minproduct_const <= Magma_minproductUnit    );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_side_const ( magma_minproduct_side_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductLeft  );
    assert( magma_minproduct_const <= Magma_minproductBothSides );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_norm_const  ( magma_minproduct_norm_t   magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductOneNorm     );
    assert( magma_minproduct_const <= Magma_minproductRealMaxNorm );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_dist_const  ( magma_minproduct_dist_t   magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductDistUniform );
    assert( magma_minproduct_const <= Magma_minproductDistNormal );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_sym_const   ( magma_minproduct_sym_t    magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductHermGeev );
    assert( magma_minproduct_const <= Magma_minproductSymPosv  );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_pack_const  ( magma_minproduct_pack_t   magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoPacking );
    assert( magma_minproduct_const <= Magma_minproductPackAll   );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_vec_const   ( magma_minproduct_vec_t    magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoVec );
    assert( magma_minproduct_const <= Magma_minproductOverwriteVec );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_range_const ( magma_minproduct_range_t  magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductRangeAll );
    assert( magma_minproduct_const <= Magma_minproductRangeI   );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_vect_const( magma_minproduct_vect_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductQ );
    assert( magma_minproduct_const <= Magma_minproductP );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_direct_const( magma_minproduct_direct_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductForward );
    assert( magma_minproduct_const <= Magma_minproductBackward );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}

extern "C"
const char* lapack_storev_const( magma_minproduct_storev_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductColumnwise );
    assert( magma_minproduct_const <= Magma_minproductRowwise    );
    return magma_minproduct2lapack_constants[ magma_minproduct_const ];
}


// ----------------------------------------
// Convert magma_minproduct constants to clAmdBlas constants.

#ifdef HAVE_clAmdBlas
const int magma_minproduct2amdblas_constants[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 100
    clAmdBlasRowMajor,      // 101: Magma_minproductRowMajor
    clAmdBlasColumnMajor,   // 102: Magma_minproductColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNoTrans,       // 111: Magma_minproductNoTrans
    clAmdBlasTrans,         // 112: Magma_minproductTrans
    clAmdBlasConjTrans,     // 113: Magma_minproductConjTrans
    0, 0, 0, 0, 0, 0, 0,
    clAmdBlasUpper,         // 121: Magma_minproductUpper
    clAmdBlasLower,         // 122: Magma_minproductLower
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNonUnit,       // 131: Magma_minproductNonUnit
    clAmdBlasUnit,          // 132: Magma_minproductUnit
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasLeft,          // 141: Magma_minproductLeft
    clAmdBlasRight,         // 142: Magma_minproductRight
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
clAmdBlasOrder       amdblas_order_const( magma_minproduct_order_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductRowMajor );
    assert( magma_minproduct_const <= Magma_minproductColMajor );
    return (clAmdBlasOrder)     magma_minproduct2amdblas_constants[ magma_minproduct_const ];
}

extern "C"
clAmdBlasTranspose   amdblas_trans_const( magma_minproduct_trans_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoTrans   );
    assert( magma_minproduct_const <= Magma_minproductConjTrans );
    return (clAmdBlasTranspose) magma_minproduct2amdblas_constants[ magma_minproduct_const ];
}

extern "C"
clAmdBlasUplo        amdblas_uplo_const ( magma_minproduct_uplo_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductUpper );
    assert( magma_minproduct_const <= Magma_minproductLower );
    return (clAmdBlasUplo)      magma_minproduct2amdblas_constants[ magma_minproduct_const ];
}

extern "C"
clAmdBlasDiag        amdblas_diag_const ( magma_minproduct_diag_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNonUnit );
    assert( magma_minproduct_const <= Magma_minproductUnit    );
    return (clAmdBlasDiag)      magma_minproduct2amdblas_constants[ magma_minproduct_const ];
}

extern "C"
clAmdBlasSide        amdblas_side_const ( magma_minproduct_side_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductLeft  );
    assert( magma_minproduct_const <= Magma_minproductRight );
    return (clAmdBlasSide)      magma_minproduct2amdblas_constants[ magma_minproduct_const ];
}
#endif  // HAVE_clAmdBlas


// ----------------------------------------
// Convert magma_minproduct constants to Nvidia CUBLAS constants.

#ifdef HAVE_CUBLAS
const int magma_minproduct2cublas_constants[] =
{
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0,                      // 100
    0,                      // 101: Magma_minproductRowMajor
    0,                      // 102: Magma_minproductColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: Magma_minproductNoTrans
    CUBLAS_OP_T,            // 112: Magma_minproductTrans
    CUBLAS_OP_C,            // 113: Magma_minproductConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: Magma_minproductUpper
    CUBLAS_FILL_MODE_LOWER, // 122: Magma_minproductLower
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: Magma_minproductNonUnit
    CUBLAS_DIAG_UNIT,       // 132: Magma_minproductUnit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: Magma_minproductLeft
    CUBLAS_SIDE_RIGHT,      // 142: Magma_minproductRight
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
cublasOperation_t    cublas_trans_const ( magma_minproduct_trans_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoTrans   );
    assert( magma_minproduct_const <= Magma_minproductConjTrans );
    return (cublasOperation_t)  magma_minproduct2cublas_constants[ magma_minproduct_const ];
}

extern "C"
cublasFillMode_t     cublas_uplo_const  ( magma_minproduct_uplo_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductUpper );
    assert( magma_minproduct_const <= Magma_minproductLower );
    return (cublasFillMode_t)   magma_minproduct2cublas_constants[ magma_minproduct_const ];
}

extern "C"
cublasDiagType_t     cublas_diag_const  ( magma_minproduct_diag_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNonUnit );
    assert( magma_minproduct_const <= Magma_minproductUnit    );
    return (cublasDiagType_t)   magma_minproduct2cublas_constants[ magma_minproduct_const ];
}

extern "C"
cublasSideMode_t     cublas_side_const  ( magma_minproduct_side_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductLeft  );
    assert( magma_minproduct_const <= Magma_minproductRight );
    return (cublasSideMode_t)   magma_minproduct2cublas_constants[ magma_minproduct_const ];
}
#endif  // HAVE_CUBLAS


// ----------------------------------------
// Convert magma_minproduct constants to CBLAS constants.
// We assume that magma_minproduct constants are consistent with cblas constants,
// so verify that with asserts.

#ifdef HAVE_CBLAS
extern "C"
enum CBLAS_ORDER     cblas_order_const  ( magma_minproduct_order_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductRowMajor );
    assert( magma_minproduct_const <= Magma_minproductColMajor );
    assert( (int)Magma_minproductRowMajor == CblasRowMajor );
    return (enum CBLAS_ORDER)     magma_minproduct_const;
}

extern "C"
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_minproduct_trans_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNoTrans   );
    assert( magma_minproduct_const <= Magma_minproductConjTrans );
    assert( (int)Magma_minproductNoTrans == CblasNoTrans );
    return (enum CBLAS_TRANSPOSE) magma_minproduct_const;
}

extern "C"
enum CBLAS_UPLO      cblas_uplo_const   ( magma_minproduct_uplo_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductUpper );
    assert( magma_minproduct_const <= Magma_minproductLower );
    assert( (int)Magma_minproductUpper == CblasUpper );
    return (enum CBLAS_UPLO)      magma_minproduct_const;
}

extern "C"
enum CBLAS_DIAG      cblas_diag_const   ( magma_minproduct_diag_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductNonUnit );
    assert( magma_minproduct_const <= Magma_minproductUnit    );
    assert( (int)Magma_minproductUnit == CblasUnit );
    return (enum CBLAS_DIAG)      magma_minproduct_const;
}

extern "C"
enum CBLAS_SIDE      cblas_side_const   ( magma_minproduct_side_t magma_minproduct_const )
{
    assert( magma_minproduct_const >= Magma_minproductLeft  );
    assert( magma_minproduct_const <= Magma_minproductRight );
    assert( (int)Magma_minproductLeft == CblasLeft );
    return (enum CBLAS_SIDE)      magma_minproduct_const;
}
#endif  // HAVE_CBLAS
