#include <assert.h>
#include <stdio.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#endif

#include "magma_tally3_types.h"

// ----------------------------------------
// Convert LAPACK character constants to MAGMA_tally3 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
// These functions and cases are in the same order as the constants are
// declared in magma_tally3_types.h

extern "C"
magma_tally3_bool_t   magma_tally3_bool_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally3False;
        case 'Y': case 'y': return Magma_tally3True;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3False;
    }
}

extern "C"
magma_tally3_order_t  magma_tally3_order_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'R': case 'r': return Magma_tally3RowMajor;
        case 'C': case 'c': return Magma_tally3ColMajor;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3RowMajor;
    }
}

extern "C"
magma_tally3_trans_t  magma_tally3_trans_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally3NoTrans;
        case 'T': case 't': return Magma_tally3Trans;
        case 'C': case 'c': return Magma_tally3ConjTrans;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3NoTrans;
    }
}

extern "C"
magma_tally3_uplo_t   magma_tally3_uplo_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally3Upper;
        case 'L': case 'l': return Magma_tally3Lower;
        default:            return Magma_tally3Full;        // see laset
    }
}

extern "C"
magma_tally3_diag_t   magma_tally3_diag_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally3NonUnit;
        case 'U': case 'u': return Magma_tally3Unit;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3NonUnit;
    }
}

extern "C"
magma_tally3_side_t   magma_tally3_side_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'L': case 'l': return Magma_tally3Left;
        case 'R': case 'r': return Magma_tally3Right;
        case 'B': case 'b': return Magma_tally3BothSides;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3Left;
    }
}

extern "C"
magma_tally3_norm_t   magma_tally3_norm_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'O': case 'o': case '1': return Magma_tally3OneNorm;
        case '2':           return Magma_tally3TwoNorm;
        case 'F': case 'f': case 'E': case 'e': return Magma_tally3FrobeniusNorm;
        case 'I': case 'i': return Magma_tally3InfNorm;
        case 'M': case 'm': return Magma_tally3MaxNorm;
        // Magma_tally3RealOneNorm
        // Magma_tally3RealInfNorm
        // Magma_tally3RealMaxNorm
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3OneNorm;
    }
}

extern "C"
magma_tally3_dist_t   magma_tally3_dist_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally3DistUniform;
        case 'S': case 's': return Magma_tally3DistSymmetric;
        case 'N': case 'n': return Magma_tally3DistNormal;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3DistUniform;
    }
}

extern "C"
magma_tally3_sym_t    magma_tally3_sym_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'H': case 'h': return Magma_tally3HermGeev;
        case 'P': case 'p': return Magma_tally3HermPoev;
        case 'N': case 'n': return Magma_tally3NonsymPosv;
        case 'S': case 's': return Magma_tally3SymPosv;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3HermGeev;
    }
}

extern "C"
magma_tally3_pack_t   magma_tally3_pack_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally3NoPacking;
        case 'U': case 'u': return Magma_tally3PackSubdiag;
        case 'L': case 'l': return Magma_tally3PackSupdiag;
        case 'C': case 'c': return Magma_tally3PackColumn;
        case 'R': case 'r': return Magma_tally3PackRow;
        case 'B': case 'b': return Magma_tally3PackLowerBand;
        case 'Q': case 'q': return Magma_tally3PackUpeprBand;
        case 'Z': case 'z': return Magma_tally3PackAll;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3NoPacking;
    }
}

extern "C"
magma_tally3_vec_t    magma_tally3_vec_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally3NoVec;
        case 'V': case 'v': return Magma_tally3Vec;
        case 'I': case 'i': return Magma_tally3IVec;
        case 'A': case 'a': return Magma_tally3AllVec;
        case 'S': case 's': return Magma_tally3SomeVec;
        case 'O': case 'o': return Magma_tally3OverwriteVec;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3NoVec;
    }
}

extern "C"
magma_tally3_range_t  magma_tally3_range_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'A': case 'a': return Magma_tally3RangeAll;
        case 'V': case 'v': return Magma_tally3RangeV;
        case 'I': case 'i': return Magma_tally3RangeI;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3RangeAll;
    }
}

extern "C"
magma_tally3_vect_t magma_tally3_vect_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'Q': case 'q': return Magma_tally3Q;
        case 'P': case 'p': return Magma_tally3P;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3Q;
    }
}

extern "C"
magma_tally3_direct_t magma_tally3_direct_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'F': case 'f': return Magma_tally3Forward;
        case 'B': case 'b': return Magma_tally3Backward;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3Forward;
    }
}

extern "C"
magma_tally3_storev_t magma_tally3_storev_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'C': case 'c': return Magma_tally3Columnwise;
        case 'R': case 'r': return Magma_tally3Rowwise;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally3Columnwise;
    }
}


// ----------------------------------------
// Convert MAGMA_tally3 constants to LAPACK constants.

const char *magma_tally32lapack_const_tally3ants[] =
{
    "No",                                    //  0: Magma_tally3False
    "Yes",                                   //  1: Magma_tally3True (zlatrs)
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
    "Row",                                   // 101: Magma_tally3RowMajor
    "Column",                                // 102: Magma_tally3ColMajor
    "", "", "", "", "", "", "", "",          // 103-110
    "No transpose",                          // 111: Magma_tally3NoTrans
    "Transpose",                             // 112: Magma_tally3Trans
    "Conjugate transpose",                   // 113: Magma_tally3ConjTrans
    "", "", "", "", "", "", "",              // 114-120
    "Upper",                                 // 121: Magma_tally3Upper
    "Lower",                                 // 122: Magma_tally3Lower
    "GFull",                                 // 123: Magma_tally3Full; see lascl for "G"
    "", "", "", "", "", "", "",              // 124-130
    "Non-unit",                              // 131: Magma_tally3NonUnit
    "Unit",                                  // 132: Magma_tally3Unit
    "", "", "", "", "", "", "", "",          // 133-140
    "Left",                                  // 141: Magma_tally3Left
    "Right",                                 // 142: Magma_tally3Right
    "Both",                                  // 143: Magma_tally3BothSides (dtrevc)
    "", "", "", "", "", "", "",              // 144-150
    "", "", "", "", "", "", "", "", "", "",  // 151-160
    "", "", "", "", "", "", "", "", "", "",  // 161-170
    "1 norm",                                // 171: Magma_tally3OneNorm
    "",                                      // 172: Magma_tally3RealOneNorm
    "2 norm",                                // 173: Magma_tally3TwoNorm
    "Frobenius norm",                        // 174: Magma_tally3FrobeniusNorm
    "Infinity norm",                         // 175: Magma_tally3InfNorm
    "",                                      // 176: Magma_tally3RealInfNorm
    "Maximum norm",                          // 177: Magma_tally3MaxNorm
    "",                                      // 178: Magma_tally3RealMaxNorm
    "", "",                                  // 179-180
    "", "", "", "", "", "", "", "", "", "",  // 181-190
    "", "", "", "", "", "", "", "", "", "",  // 191-200
    "Uniform",                               // 201: Magma_tally3DistUniform
    "Symmetric",                             // 202: Magma_tally3DistSymmetric
    "Normal",                                // 203: Magma_tally3DistNormal
    "", "", "", "", "", "", "",              // 204-210
    "", "", "", "", "", "", "", "", "", "",  // 211-220
    "", "", "", "", "", "", "", "", "", "",  // 221-230
    "", "", "", "", "", "", "", "", "", "",  // 231-240
    "Hermitian",                             // 241 Magma_tally3HermGeev
    "Positive ev Hermitian",                 // 242 Magma_tally3HermPoev
    "NonSymmetric pos sv",                   // 243 Magma_tally3NonsymPosv
    "Symmetric pos sv",                      // 244 Magma_tally3SymPosv
    "", "", "", "", "", "",                  // 245-250
    "", "", "", "", "", "", "", "", "", "",  // 251-260
    "", "", "", "", "", "", "", "", "", "",  // 261-270
    "", "", "", "", "", "", "", "", "", "",  // 271-280
    "", "", "", "", "", "", "", "", "", "",  // 281-290
    "No Packing",                            // 291 Magma_tally3NoPacking
    "U zero out subdiag",                    // 292 Magma_tally3PackSubdiag
    "L zero out superdiag",                  // 293 Magma_tally3PackSupdiag
    "C",                                     // 294 Magma_tally3PackColumn
    "R",                                     // 295 Magma_tally3PackRow
    "B",                                     // 296 Magma_tally3PackLowerBand
    "Q",                                     // 297 Magma_tally3PackUpeprBand
    "Z",                                     // 298 Magma_tally3PackAll
    "", "",                                  // 299-300
    "No vectors",                            // 301 Magma_tally3NoVec
    "Vectors needed",                        // 302 Magma_tally3Vec
    "I",                                     // 303 Magma_tally3IVec
    "All",                                   // 304 Magma_tally3AllVec
    "Some",                                  // 305 Magma_tally3SomeVec
    "Overwrite",                             // 306 Magma_tally3OverwriteVec
    "", "", "", "",                          // 307-310
    "All",                                   // 311 Magma_tally3RangeAll
    "V",                                     // 312 Magma_tally3RangeV
    "I",                                     // 313 Magma_tally3RangeI
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
    "Forward",                               // 391: Magma_tally3Forward
    "Backward",                              // 392: Magma_tally3Backward
    "", "", "", "", "", "", "", "",          // 393-400
    "Columnwise",                            // 401: Magma_tally3Columnwise
    "Rowwise",                               // 402: Magma_tally3Rowwise
    "", "", "", "", "", "", "", ""           // 403-410
    // Remember to add a comma!
};

extern "C"
const char* lapack_const_tally3( int magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally32lapack_Min );
    assert( magma_tally3_const <= Magma_tally32lapack_Max );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_bool_const_tally3( magma_tally3_bool_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3False );
    assert( magma_tally3_const <= Magma_tally3True  );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_order_const_tally3( magma_tally3_order_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3RowMajor );
    assert( magma_tally3_const <= Magma_tally3ColMajor );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_trans_const_tally3( magma_tally3_trans_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoTrans   );
    assert( magma_tally3_const <= Magma_tally3ConjTrans );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_uplo_const_tally3 ( magma_tally3_uplo_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Upper );
    assert( magma_tally3_const <= Magma_tally3Full  );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_diag_const_tally3 ( magma_tally3_diag_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NonUnit );
    assert( magma_tally3_const <= Magma_tally3Unit    );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_side_const_tally3 ( magma_tally3_side_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Left  );
    assert( magma_tally3_const <= Magma_tally3BothSides );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_norm_const_tally3  ( magma_tally3_norm_t   magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3OneNorm     );
    assert( magma_tally3_const <= Magma_tally3RealMaxNorm );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_dist_const_tally3  ( magma_tally3_dist_t   magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3DistUniform );
    assert( magma_tally3_const <= Magma_tally3DistNormal );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_sym_const_tally3   ( magma_tally3_sym_t    magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3HermGeev );
    assert( magma_tally3_const <= Magma_tally3SymPosv  );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_pack_const_tally3  ( magma_tally3_pack_t   magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoPacking );
    assert( magma_tally3_const <= Magma_tally3PackAll   );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_vec_const_tally3   ( magma_tally3_vec_t    magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoVec );
    assert( magma_tally3_const <= Magma_tally3OverwriteVec );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_range_const_tally3 ( magma_tally3_range_t  magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3RangeAll );
    assert( magma_tally3_const <= Magma_tally3RangeI   );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_vect_const_tally3( magma_tally3_vect_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Q );
    assert( magma_tally3_const <= Magma_tally3P );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_direct_const_tally3( magma_tally3_direct_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Forward );
    assert( magma_tally3_const <= Magma_tally3Backward );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}

extern "C"
const char* lapack_storev_const_tally3( magma_tally3_storev_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Columnwise );
    assert( magma_tally3_const <= Magma_tally3Rowwise    );
    return magma_tally32lapack_const_tally3ants[ magma_tally3_const ];
}


// ----------------------------------------
// Convert magma_tally3 constants to clAmdBlas constants.

#ifdef HAVE_clAmdBlas
const int magma_tally32amdblas_constants[] =
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
    clAmdBlasRowMajor,      // 101: Magma_tally3RowMajor
    clAmdBlasColumnMajor,   // 102: Magma_tally3ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNoTrans,       // 111: Magma_tally3NoTrans
    clAmdBlasTrans,         // 112: Magma_tally3Trans
    clAmdBlasConjTrans,     // 113: Magma_tally3ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    clAmdBlasUpper,         // 121: Magma_tally3Upper
    clAmdBlasLower,         // 122: Magma_tally3Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNonUnit,       // 131: Magma_tally3NonUnit
    clAmdBlasUnit,          // 132: Magma_tally3Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasLeft,          // 141: Magma_tally3Left
    clAmdBlasRight,         // 142: Magma_tally3Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
clAmdBlasOrder       amdblas_order_const( magma_tally3_order_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3RowMajor );
    assert( magma_tally3_const <= Magma_tally3ColMajor );
    return (clAmdBlasOrder)     magma_tally32amdblas_constants[ magma_tally3_const ];
}

extern "C"
clAmdBlasTranspose   amdblas_trans_const( magma_tally3_trans_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoTrans   );
    assert( magma_tally3_const <= Magma_tally3ConjTrans );
    return (clAmdBlasTranspose) magma_tally32amdblas_constants[ magma_tally3_const ];
}

extern "C"
clAmdBlasUplo        amdblas_uplo_const ( magma_tally3_uplo_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Upper );
    assert( magma_tally3_const <= Magma_tally3Lower );
    return (clAmdBlasUplo)      magma_tally32amdblas_constants[ magma_tally3_const ];
}

extern "C"
clAmdBlasDiag        amdblas_diag_const ( magma_tally3_diag_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NonUnit );
    assert( magma_tally3_const <= Magma_tally3Unit    );
    return (clAmdBlasDiag)      magma_tally32amdblas_constants[ magma_tally3_const ];
}

extern "C"
clAmdBlasSide        amdblas_side_const ( magma_tally3_side_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Left  );
    assert( magma_tally3_const <= Magma_tally3Right );
    return (clAmdBlasSide)      magma_tally32amdblas_constants[ magma_tally3_const ];
}
#endif  // HAVE_clAmdBlas


// ----------------------------------------
// Convert magma_tally3 constants to Nvidia CUBLAS constants.

#ifdef HAVE_CUBLAS
const int magma_tally32cublas_constants[] =
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
    0,                      // 101: Magma_tally3RowMajor
    0,                      // 102: Magma_tally3ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: Magma_tally3NoTrans
    CUBLAS_OP_T,            // 112: Magma_tally3Trans
    CUBLAS_OP_C,            // 113: Magma_tally3ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: Magma_tally3Upper
    CUBLAS_FILL_MODE_LOWER, // 122: Magma_tally3Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: Magma_tally3NonUnit
    CUBLAS_DIAG_UNIT,       // 132: Magma_tally3Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: Magma_tally3Left
    CUBLAS_SIDE_RIGHT,      // 142: Magma_tally3Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
cublasOperation_t    cublas_trans_const_tally3 ( magma_tally3_trans_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoTrans   );
    assert( magma_tally3_const <= Magma_tally3ConjTrans );
    return (cublasOperation_t)  magma_tally32cublas_constants[ magma_tally3_const ];
}

extern "C"
cublasFillMode_t     cublas_uplo_const_tally3  ( magma_tally3_uplo_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Upper );
    assert( magma_tally3_const <= Magma_tally3Lower );
    return (cublasFillMode_t)   magma_tally32cublas_constants[ magma_tally3_const ];
}

extern "C"
cublasDiagType_t     cublas_diag_const_tally3  ( magma_tally3_diag_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NonUnit );
    assert( magma_tally3_const <= Magma_tally3Unit    );
    return (cublasDiagType_t)   magma_tally32cublas_constants[ magma_tally3_const ];
}

extern "C"
cublasSideMode_t     cublas_side_const_tally3  ( magma_tally3_side_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Left  );
    assert( magma_tally3_const <= Magma_tally3Right );
    return (cublasSideMode_t)   magma_tally32cublas_constants[ magma_tally3_const ];
}
#endif  // HAVE_CUBLAS


// ----------------------------------------
// Convert magma_tally3 constants to CBLAS constants.
// We assume that magma_tally3 constants are consistent with cblas constants,
// so verify that with asserts.

#ifdef HAVE_CBLAS
extern "C"
enum CBLAS_ORDER     cblas_order_const  ( magma_tally3_order_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3RowMajor );
    assert( magma_tally3_const <= Magma_tally3ColMajor );
    assert( (int)Magma_tally3RowMajor == CblasRowMajor );
    return (enum CBLAS_ORDER)     magma_tally3_const;
}

extern "C"
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally3_trans_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NoTrans   );
    assert( magma_tally3_const <= Magma_tally3ConjTrans );
    assert( (int)Magma_tally3NoTrans == CblasNoTrans );
    return (enum CBLAS_TRANSPOSE) magma_tally3_const;
}

extern "C"
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally3_uplo_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Upper );
    assert( magma_tally3_const <= Magma_tally3Lower );
    assert( (int)Magma_tally3Upper == CblasUpper );
    return (enum CBLAS_UPLO)      magma_tally3_const;
}

extern "C"
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally3_diag_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3NonUnit );
    assert( magma_tally3_const <= Magma_tally3Unit    );
    assert( (int)Magma_tally3Unit == CblasUnit );
    return (enum CBLAS_DIAG)      magma_tally3_const;
}

extern "C"
enum CBLAS_SIDE      cblas_side_const   ( magma_tally3_side_t magma_tally3_const )
{
    assert( magma_tally3_const >= Magma_tally3Left  );
    assert( magma_tally3_const <= Magma_tally3Right );
    assert( (int)Magma_tally3Left == CblasLeft );
    return (enum CBLAS_SIDE)      magma_tally3_const;
}
#endif  // HAVE_CBLAS
