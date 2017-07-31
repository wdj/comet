#include <assert.h>
#include <stdio.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#endif

#include "magma_tally2_types.h"

// ----------------------------------------
// Convert LAPACK character constants to MAGMA_tally2 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
// These functions and cases are in the same order as the constants are
// declared in magma_tally2_types.h

extern "C"
magma_tally2_bool_t   magma_tally2_bool_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally2False;
        case 'Y': case 'y': return Magma_tally2True;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2False;
    }
}

extern "C"
magma_tally2_order_t  magma_tally2_order_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'R': case 'r': return Magma_tally2RowMajor;
        case 'C': case 'c': return Magma_tally2ColMajor;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2RowMajor;
    }
}

extern "C"
magma_tally2_trans_t  magma_tally2_trans_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally2NoTrans;
        case 'T': case 't': return Magma_tally2Trans;
        case 'C': case 'c': return Magma_tally2ConjTrans;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2NoTrans;
    }
}

extern "C"
magma_tally2_uplo_t   magma_tally2_uplo_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally2Upper;
        case 'L': case 'l': return Magma_tally2Lower;
        default:            return Magma_tally2Full;        // see laset
    }
}

extern "C"
magma_tally2_diag_t   magma_tally2_diag_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally2NonUnit;
        case 'U': case 'u': return Magma_tally2Unit;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2NonUnit;
    }
}

extern "C"
magma_tally2_side_t   magma_tally2_side_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'L': case 'l': return Magma_tally2Left;
        case 'R': case 'r': return Magma_tally2Right;
        case 'B': case 'b': return Magma_tally2BothSides;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2Left;
    }
}

extern "C"
magma_tally2_norm_t   magma_tally2_norm_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'O': case 'o': case '1': return Magma_tally2OneNorm;
        case '2':           return Magma_tally2TwoNorm;
        case 'F': case 'f': case 'E': case 'e': return Magma_tally2FrobeniusNorm;
        case 'I': case 'i': return Magma_tally2InfNorm;
        case 'M': case 'm': return Magma_tally2MaxNorm;
        // Magma_tally2RealOneNorm
        // Magma_tally2RealInfNorm
        // Magma_tally2RealMaxNorm
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2OneNorm;
    }
}

extern "C"
magma_tally2_dist_t   magma_tally2_dist_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally2DistUniform;
        case 'S': case 's': return Magma_tally2DistSymmetric;
        case 'N': case 'n': return Magma_tally2DistNormal;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2DistUniform;
    }
}

extern "C"
magma_tally2_sym_t    magma_tally2_sym_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'H': case 'h': return Magma_tally2HermGeev;
        case 'P': case 'p': return Magma_tally2HermPoev;
        case 'N': case 'n': return Magma_tally2NonsymPosv;
        case 'S': case 's': return Magma_tally2SymPosv;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2HermGeev;
    }
}

extern "C"
magma_tally2_pack_t   magma_tally2_pack_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally2NoPacking;
        case 'U': case 'u': return Magma_tally2PackSubdiag;
        case 'L': case 'l': return Magma_tally2PackSupdiag;
        case 'C': case 'c': return Magma_tally2PackColumn;
        case 'R': case 'r': return Magma_tally2PackRow;
        case 'B': case 'b': return Magma_tally2PackLowerBand;
        case 'Q': case 'q': return Magma_tally2PackUpeprBand;
        case 'Z': case 'z': return Magma_tally2PackAll;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2NoPacking;
    }
}

extern "C"
magma_tally2_vec_t    magma_tally2_vec_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally2NoVec;
        case 'V': case 'v': return Magma_tally2Vec;
        case 'I': case 'i': return Magma_tally2IVec;
        case 'A': case 'a': return Magma_tally2AllVec;
        case 'S': case 's': return Magma_tally2SomeVec;
        case 'O': case 'o': return Magma_tally2OverwriteVec;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2NoVec;
    }
}

extern "C"
magma_tally2_range_t  magma_tally2_range_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'A': case 'a': return Magma_tally2RangeAll;
        case 'V': case 'v': return Magma_tally2RangeV;
        case 'I': case 'i': return Magma_tally2RangeI;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2RangeAll;
    }
}

extern "C"
magma_tally2_vect_t magma_tally2_vect_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'Q': case 'q': return Magma_tally2Q;
        case 'P': case 'p': return Magma_tally2P;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2Q;
    }
}

extern "C"
magma_tally2_direct_t magma_tally2_direct_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'F': case 'f': return Magma_tally2Forward;
        case 'B': case 'b': return Magma_tally2Backward;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2Forward;
    }
}

extern "C"
magma_tally2_storev_t magma_tally2_storev_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'C': case 'c': return Magma_tally2Columnwise;
        case 'R': case 'r': return Magma_tally2Rowwise;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally2Columnwise;
    }
}


// ----------------------------------------
// Convert MAGMA_tally2 constants to LAPACK constants.

const char *magma_tally22lapack_const_tally2ants[] =
{
    "No",                                    //  0: Magma_tally2False
    "Yes",                                   //  1: Magma_tally2True (zlatrs)
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
    "Row",                                   // 101: Magma_tally2RowMajor
    "Column",                                // 102: Magma_tally2ColMajor
    "", "", "", "", "", "", "", "",          // 103-110
    "No transpose",                          // 111: Magma_tally2NoTrans
    "Transpose",                             // 112: Magma_tally2Trans
    "Conjugate transpose",                   // 113: Magma_tally2ConjTrans
    "", "", "", "", "", "", "",              // 114-120
    "Upper",                                 // 121: Magma_tally2Upper
    "Lower",                                 // 122: Magma_tally2Lower
    "GFull",                                 // 123: Magma_tally2Full; see lascl for "G"
    "", "", "", "", "", "", "",              // 124-130
    "Non-unit",                              // 131: Magma_tally2NonUnit
    "Unit",                                  // 132: Magma_tally2Unit
    "", "", "", "", "", "", "", "",          // 133-140
    "Left",                                  // 141: Magma_tally2Left
    "Right",                                 // 142: Magma_tally2Right
    "Both",                                  // 143: Magma_tally2BothSides (dtrevc)
    "", "", "", "", "", "", "",              // 144-150
    "", "", "", "", "", "", "", "", "", "",  // 151-160
    "", "", "", "", "", "", "", "", "", "",  // 161-170
    "1 norm",                                // 171: Magma_tally2OneNorm
    "",                                      // 172: Magma_tally2RealOneNorm
    "2 norm",                                // 173: Magma_tally2TwoNorm
    "Frobenius norm",                        // 174: Magma_tally2FrobeniusNorm
    "Infinity norm",                         // 175: Magma_tally2InfNorm
    "",                                      // 176: Magma_tally2RealInfNorm
    "Maximum norm",                          // 177: Magma_tally2MaxNorm
    "",                                      // 178: Magma_tally2RealMaxNorm
    "", "",                                  // 179-180
    "", "", "", "", "", "", "", "", "", "",  // 181-190
    "", "", "", "", "", "", "", "", "", "",  // 191-200
    "Uniform",                               // 201: Magma_tally2DistUniform
    "Symmetric",                             // 202: Magma_tally2DistSymmetric
    "Normal",                                // 203: Magma_tally2DistNormal
    "", "", "", "", "", "", "",              // 204-210
    "", "", "", "", "", "", "", "", "", "",  // 211-220
    "", "", "", "", "", "", "", "", "", "",  // 221-230
    "", "", "", "", "", "", "", "", "", "",  // 231-240
    "Hermitian",                             // 241 Magma_tally2HermGeev
    "Positive ev Hermitian",                 // 242 Magma_tally2HermPoev
    "NonSymmetric pos sv",                   // 243 Magma_tally2NonsymPosv
    "Symmetric pos sv",                      // 244 Magma_tally2SymPosv
    "", "", "", "", "", "",                  // 245-250
    "", "", "", "", "", "", "", "", "", "",  // 251-260
    "", "", "", "", "", "", "", "", "", "",  // 261-270
    "", "", "", "", "", "", "", "", "", "",  // 271-280
    "", "", "", "", "", "", "", "", "", "",  // 281-290
    "No Packing",                            // 291 Magma_tally2NoPacking
    "U zero out subdiag",                    // 292 Magma_tally2PackSubdiag
    "L zero out superdiag",                  // 293 Magma_tally2PackSupdiag
    "C",                                     // 294 Magma_tally2PackColumn
    "R",                                     // 295 Magma_tally2PackRow
    "B",                                     // 296 Magma_tally2PackLowerBand
    "Q",                                     // 297 Magma_tally2PackUpeprBand
    "Z",                                     // 298 Magma_tally2PackAll
    "", "",                                  // 299-300
    "No vectors",                            // 301 Magma_tally2NoVec
    "Vectors needed",                        // 302 Magma_tally2Vec
    "I",                                     // 303 Magma_tally2IVec
    "All",                                   // 304 Magma_tally2AllVec
    "Some",                                  // 305 Magma_tally2SomeVec
    "Overwrite",                             // 306 Magma_tally2OverwriteVec
    "", "", "", "",                          // 307-310
    "All",                                   // 311 Magma_tally2RangeAll
    "V",                                     // 312 Magma_tally2RangeV
    "I",                                     // 313 Magma_tally2RangeI
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
    "Forward",                               // 391: Magma_tally2Forward
    "Backward",                              // 392: Magma_tally2Backward
    "", "", "", "", "", "", "", "",          // 393-400
    "Columnwise",                            // 401: Magma_tally2Columnwise
    "Rowwise",                               // 402: Magma_tally2Rowwise
    "", "", "", "", "", "", "", ""           // 403-410
    // Remember to add a comma!
};

extern "C"
const char* lapack_const_tally2( int magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally22lapack_Min );
    assert( magma_tally2_const <= Magma_tally22lapack_Max );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_bool_const_tally2( magma_tally2_bool_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2False );
    assert( magma_tally2_const <= Magma_tally2True  );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_order_const_tally2( magma_tally2_order_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2RowMajor );
    assert( magma_tally2_const <= Magma_tally2ColMajor );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_trans_const_tally2( magma_tally2_trans_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoTrans   );
    assert( magma_tally2_const <= Magma_tally2ConjTrans );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_uplo_const_tally2 ( magma_tally2_uplo_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Upper );
    assert( magma_tally2_const <= Magma_tally2Full  );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_diag_const_tally2 ( magma_tally2_diag_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NonUnit );
    assert( magma_tally2_const <= Magma_tally2Unit    );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_side_const_tally2 ( magma_tally2_side_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Left  );
    assert( magma_tally2_const <= Magma_tally2BothSides );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_norm_const_tally2  ( magma_tally2_norm_t   magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2OneNorm     );
    assert( magma_tally2_const <= Magma_tally2RealMaxNorm );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_dist_const_tally2  ( magma_tally2_dist_t   magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2DistUniform );
    assert( magma_tally2_const <= Magma_tally2DistNormal );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_sym_const_tally2   ( magma_tally2_sym_t    magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2HermGeev );
    assert( magma_tally2_const <= Magma_tally2SymPosv  );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_pack_const_tally2  ( magma_tally2_pack_t   magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoPacking );
    assert( magma_tally2_const <= Magma_tally2PackAll   );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_vec_const_tally2   ( magma_tally2_vec_t    magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoVec );
    assert( magma_tally2_const <= Magma_tally2OverwriteVec );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_range_const_tally2 ( magma_tally2_range_t  magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2RangeAll );
    assert( magma_tally2_const <= Magma_tally2RangeI   );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_vect_const_tally2( magma_tally2_vect_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Q );
    assert( magma_tally2_const <= Magma_tally2P );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_direct_const_tally2( magma_tally2_direct_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Forward );
    assert( magma_tally2_const <= Magma_tally2Backward );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}

extern "C"
const char* lapack_storev_const_tally2( magma_tally2_storev_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Columnwise );
    assert( magma_tally2_const <= Magma_tally2Rowwise    );
    return magma_tally22lapack_const_tally2ants[ magma_tally2_const ];
}


// ----------------------------------------
// Convert magma_tally2 constants to clAmdBlas constants.

#ifdef HAVE_clAmdBlas
const int magma_tally22amdblas_constants[] =
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
    clAmdBlasRowMajor,      // 101: Magma_tally2RowMajor
    clAmdBlasColumnMajor,   // 102: Magma_tally2ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNoTrans,       // 111: Magma_tally2NoTrans
    clAmdBlasTrans,         // 112: Magma_tally2Trans
    clAmdBlasConjTrans,     // 113: Magma_tally2ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    clAmdBlasUpper,         // 121: Magma_tally2Upper
    clAmdBlasLower,         // 122: Magma_tally2Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNonUnit,       // 131: Magma_tally2NonUnit
    clAmdBlasUnit,          // 132: Magma_tally2Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasLeft,          // 141: Magma_tally2Left
    clAmdBlasRight,         // 142: Magma_tally2Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
clAmdBlasOrder       amdblas_order_const( magma_tally2_order_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2RowMajor );
    assert( magma_tally2_const <= Magma_tally2ColMajor );
    return (clAmdBlasOrder)     magma_tally22amdblas_constants[ magma_tally2_const ];
}

extern "C"
clAmdBlasTranspose   amdblas_trans_const( magma_tally2_trans_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoTrans   );
    assert( magma_tally2_const <= Magma_tally2ConjTrans );
    return (clAmdBlasTranspose) magma_tally22amdblas_constants[ magma_tally2_const ];
}

extern "C"
clAmdBlasUplo        amdblas_uplo_const ( magma_tally2_uplo_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Upper );
    assert( magma_tally2_const <= Magma_tally2Lower );
    return (clAmdBlasUplo)      magma_tally22amdblas_constants[ magma_tally2_const ];
}

extern "C"
clAmdBlasDiag        amdblas_diag_const ( magma_tally2_diag_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NonUnit );
    assert( magma_tally2_const <= Magma_tally2Unit    );
    return (clAmdBlasDiag)      magma_tally22amdblas_constants[ magma_tally2_const ];
}

extern "C"
clAmdBlasSide        amdblas_side_const ( magma_tally2_side_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Left  );
    assert( magma_tally2_const <= Magma_tally2Right );
    return (clAmdBlasSide)      magma_tally22amdblas_constants[ magma_tally2_const ];
}
#endif  // HAVE_clAmdBlas


// ----------------------------------------
// Convert magma_tally2 constants to Nvidia CUBLAS constants.

#ifdef HAVE_CUBLAS
const int magma_tally22cublas_constants[] =
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
    0,                      // 101: Magma_tally2RowMajor
    0,                      // 102: Magma_tally2ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: Magma_tally2NoTrans
    CUBLAS_OP_T,            // 112: Magma_tally2Trans
    CUBLAS_OP_C,            // 113: Magma_tally2ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: Magma_tally2Upper
    CUBLAS_FILL_MODE_LOWER, // 122: Magma_tally2Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: Magma_tally2NonUnit
    CUBLAS_DIAG_UNIT,       // 132: Magma_tally2Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: Magma_tally2Left
    CUBLAS_SIDE_RIGHT,      // 142: Magma_tally2Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
cublasOperation_t    cublas_trans_const_tally2 ( magma_tally2_trans_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoTrans   );
    assert( magma_tally2_const <= Magma_tally2ConjTrans );
    return (cublasOperation_t)  magma_tally22cublas_constants[ magma_tally2_const ];
}

extern "C"
cublasFillMode_t     cublas_uplo_const_tally2  ( magma_tally2_uplo_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Upper );
    assert( magma_tally2_const <= Magma_tally2Lower );
    return (cublasFillMode_t)   magma_tally22cublas_constants[ magma_tally2_const ];
}

extern "C"
cublasDiagType_t     cublas_diag_const_tally2  ( magma_tally2_diag_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NonUnit );
    assert( magma_tally2_const <= Magma_tally2Unit    );
    return (cublasDiagType_t)   magma_tally22cublas_constants[ magma_tally2_const ];
}

extern "C"
cublasSideMode_t     cublas_side_const_tally2  ( magma_tally2_side_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Left  );
    assert( magma_tally2_const <= Magma_tally2Right );
    return (cublasSideMode_t)   magma_tally22cublas_constants[ magma_tally2_const ];
}
#endif  // HAVE_CUBLAS


// ----------------------------------------
// Convert magma_tally2 constants to CBLAS constants.
// We assume that magma_tally2 constants are consistent with cblas constants,
// so verify that with asserts.

#ifdef HAVE_CBLAS
extern "C"
enum CBLAS_ORDER     cblas_order_const  ( magma_tally2_order_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2RowMajor );
    assert( magma_tally2_const <= Magma_tally2ColMajor );
    assert( (int)Magma_tally2RowMajor == CblasRowMajor );
    return (enum CBLAS_ORDER)     magma_tally2_const;
}

extern "C"
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally2_trans_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NoTrans   );
    assert( magma_tally2_const <= Magma_tally2ConjTrans );
    assert( (int)Magma_tally2NoTrans == CblasNoTrans );
    return (enum CBLAS_TRANSPOSE) magma_tally2_const;
}

extern "C"
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally2_uplo_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Upper );
    assert( magma_tally2_const <= Magma_tally2Lower );
    assert( (int)Magma_tally2Upper == CblasUpper );
    return (enum CBLAS_UPLO)      magma_tally2_const;
}

extern "C"
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally2_diag_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2NonUnit );
    assert( magma_tally2_const <= Magma_tally2Unit    );
    assert( (int)Magma_tally2Unit == CblasUnit );
    return (enum CBLAS_DIAG)      magma_tally2_const;
}

extern "C"
enum CBLAS_SIDE      cblas_side_const   ( magma_tally2_side_t magma_tally2_const )
{
    assert( magma_tally2_const >= Magma_tally2Left  );
    assert( magma_tally2_const <= Magma_tally2Right );
    assert( (int)Magma_tally2Left == CblasLeft );
    return (enum CBLAS_SIDE)      magma_tally2_const;
}
#endif  // HAVE_CBLAS
