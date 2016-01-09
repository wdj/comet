#include <assert.h>
#include <stdio.h>

#ifdef HAVE_CUBLAS
#include <cublas_v2.h>
#endif

#include "magma_tally4_types.h"

// ----------------------------------------
// Convert LAPACK character constants to MAGMA_tally4 constants.
// This is a one-to-many mapping, requiring multiple translators
// (e.g., "N" can be NoTrans or NonUnit or NoVec).
// These functions and cases are in the same order as the constants are
// declared in magma_tally4_types.h

extern "C"
magma_tally4_bool_t   magma_tally4_bool_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally4False;
        case 'Y': case 'y': return Magma_tally4True;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4False;
    }
}

extern "C"
magma_tally4_order_t  magma_tally4_order_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'R': case 'r': return Magma_tally4RowMajor;
        case 'C': case 'c': return Magma_tally4ColMajor;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4RowMajor;
    }
}

extern "C"
magma_tally4_trans_t  magma_tally4_trans_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally4NoTrans;
        case 'T': case 't': return Magma_tally4Trans;
        case 'C': case 'c': return Magma_tally4ConjTrans;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4NoTrans;
    }
}

extern "C"
magma_tally4_uplo_t   magma_tally4_uplo_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally4Upper;
        case 'L': case 'l': return Magma_tally4Lower;
        default:            return Magma_tally4Full;        // see laset
    }
}

extern "C"
magma_tally4_diag_t   magma_tally4_diag_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally4NonUnit;
        case 'U': case 'u': return Magma_tally4Unit;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4NonUnit;
    }
}

extern "C"
magma_tally4_side_t   magma_tally4_side_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'L': case 'l': return Magma_tally4Left;
        case 'R': case 'r': return Magma_tally4Right;
        case 'B': case 'b': return Magma_tally4BothSides;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4Left;
    }
}

extern "C"
magma_tally4_norm_t   magma_tally4_norm_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'O': case 'o': case '1': return Magma_tally4OneNorm;
        case '2':           return Magma_tally4TwoNorm;
        case 'F': case 'f': case 'E': case 'e': return Magma_tally4FrobeniusNorm;
        case 'I': case 'i': return Magma_tally4InfNorm;
        case 'M': case 'm': return Magma_tally4MaxNorm;
        // Magma_tally4RealOneNorm
        // Magma_tally4RealInfNorm
        // Magma_tally4RealMaxNorm
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4OneNorm;
    }
}

extern "C"
magma_tally4_dist_t   magma_tally4_dist_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'U': case 'u': return Magma_tally4DistUniform;
        case 'S': case 's': return Magma_tally4DistSymmetric;
        case 'N': case 'n': return Magma_tally4DistNormal;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4DistUniform;
    }
}

extern "C"
magma_tally4_sym_t    magma_tally4_sym_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'H': case 'h': return Magma_tally4HermGeev;
        case 'P': case 'p': return Magma_tally4HermPoev;
        case 'N': case 'n': return Magma_tally4NonsymPosv;
        case 'S': case 's': return Magma_tally4SymPosv;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4HermGeev;
    }
}

extern "C"
magma_tally4_pack_t   magma_tally4_pack_const  ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally4NoPacking;
        case 'U': case 'u': return Magma_tally4PackSubdiag;
        case 'L': case 'l': return Magma_tally4PackSupdiag;
        case 'C': case 'c': return Magma_tally4PackColumn;
        case 'R': case 'r': return Magma_tally4PackRow;
        case 'B': case 'b': return Magma_tally4PackLowerBand;
        case 'Q': case 'q': return Magma_tally4PackUpeprBand;
        case 'Z': case 'z': return Magma_tally4PackAll;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4NoPacking;
    }
}

extern "C"
magma_tally4_vec_t    magma_tally4_vec_const   ( char lapack_char )
{
    switch( lapack_char ) {
        case 'N': case 'n': return Magma_tally4NoVec;
        case 'V': case 'v': return Magma_tally4Vec;
        case 'I': case 'i': return Magma_tally4IVec;
        case 'A': case 'a': return Magma_tally4AllVec;
        case 'S': case 's': return Magma_tally4SomeVec;
        case 'O': case 'o': return Magma_tally4OverwriteVec;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4NoVec;
    }
}

extern "C"
magma_tally4_range_t  magma_tally4_range_const ( char lapack_char )
{
    switch( lapack_char ) {
        case 'A': case 'a': return Magma_tally4RangeAll;
        case 'V': case 'v': return Magma_tally4RangeV;
        case 'I': case 'i': return Magma_tally4RangeI;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4RangeAll;
    }
}

extern "C"
magma_tally4_vect_t magma_tally4_vect_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'Q': case 'q': return Magma_tally4Q;
        case 'P': case 'p': return Magma_tally4P;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4Q;
    }
}

extern "C"
magma_tally4_direct_t magma_tally4_direct_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'F': case 'f': return Magma_tally4Forward;
        case 'B': case 'b': return Magma_tally4Backward;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4Forward;
    }
}

extern "C"
magma_tally4_storev_t magma_tally4_storev_const( char lapack_char )
{
    switch( lapack_char ) {
        case 'C': case 'c': return Magma_tally4Columnwise;
        case 'R': case 'r': return Magma_tally4Rowwise;
        default:
            fprintf( stderr, "Error in %s: unexpected value %c\n", __func__, lapack_char );
            return Magma_tally4Columnwise;
    }
}


// ----------------------------------------
// Convert MAGMA_tally4 constants to LAPACK constants.

const char *magma_tally42lapack_const_tally4ants[] =
{
    "No",                                    //  0: Magma_tally4False
    "Yes",                                   //  1: Magma_tally4True (zlatrs)
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
    "Row",                                   // 101: Magma_tally4RowMajor
    "Column",                                // 102: Magma_tally4ColMajor
    "", "", "", "", "", "", "", "",          // 103-110
    "No transpose",                          // 111: Magma_tally4NoTrans
    "Transpose",                             // 112: Magma_tally4Trans
    "Conjugate transpose",                   // 113: Magma_tally4ConjTrans
    "", "", "", "", "", "", "",              // 114-120
    "Upper",                                 // 121: Magma_tally4Upper
    "Lower",                                 // 122: Magma_tally4Lower
    "GFull",                                 // 123: Magma_tally4Full; see lascl for "G"
    "", "", "", "", "", "", "",              // 124-130
    "Non-unit",                              // 131: Magma_tally4NonUnit
    "Unit",                                  // 132: Magma_tally4Unit
    "", "", "", "", "", "", "", "",          // 133-140
    "Left",                                  // 141: Magma_tally4Left
    "Right",                                 // 142: Magma_tally4Right
    "Both",                                  // 143: Magma_tally4BothSides (dtrevc)
    "", "", "", "", "", "", "",              // 144-150
    "", "", "", "", "", "", "", "", "", "",  // 151-160
    "", "", "", "", "", "", "", "", "", "",  // 161-170
    "1 norm",                                // 171: Magma_tally4OneNorm
    "",                                      // 172: Magma_tally4RealOneNorm
    "2 norm",                                // 173: Magma_tally4TwoNorm
    "Frobenius norm",                        // 174: Magma_tally4FrobeniusNorm
    "Infinity norm",                         // 175: Magma_tally4InfNorm
    "",                                      // 176: Magma_tally4RealInfNorm
    "Maximum norm",                          // 177: Magma_tally4MaxNorm
    "",                                      // 178: Magma_tally4RealMaxNorm
    "", "",                                  // 179-180
    "", "", "", "", "", "", "", "", "", "",  // 181-190
    "", "", "", "", "", "", "", "", "", "",  // 191-200
    "Uniform",                               // 201: Magma_tally4DistUniform
    "Symmetric",                             // 202: Magma_tally4DistSymmetric
    "Normal",                                // 203: Magma_tally4DistNormal
    "", "", "", "", "", "", "",              // 204-210
    "", "", "", "", "", "", "", "", "", "",  // 211-220
    "", "", "", "", "", "", "", "", "", "",  // 221-230
    "", "", "", "", "", "", "", "", "", "",  // 231-240
    "Hermitian",                             // 241 Magma_tally4HermGeev
    "Positive ev Hermitian",                 // 242 Magma_tally4HermPoev
    "NonSymmetric pos sv",                   // 243 Magma_tally4NonsymPosv
    "Symmetric pos sv",                      // 244 Magma_tally4SymPosv
    "", "", "", "", "", "",                  // 245-250
    "", "", "", "", "", "", "", "", "", "",  // 251-260
    "", "", "", "", "", "", "", "", "", "",  // 261-270
    "", "", "", "", "", "", "", "", "", "",  // 271-280
    "", "", "", "", "", "", "", "", "", "",  // 281-290
    "No Packing",                            // 291 Magma_tally4NoPacking
    "U zero out subdiag",                    // 292 Magma_tally4PackSubdiag
    "L zero out superdiag",                  // 293 Magma_tally4PackSupdiag
    "C",                                     // 294 Magma_tally4PackColumn
    "R",                                     // 295 Magma_tally4PackRow
    "B",                                     // 296 Magma_tally4PackLowerBand
    "Q",                                     // 297 Magma_tally4PackUpeprBand
    "Z",                                     // 298 Magma_tally4PackAll
    "", "",                                  // 299-300
    "No vectors",                            // 301 Magma_tally4NoVec
    "Vectors needed",                        // 302 Magma_tally4Vec
    "I",                                     // 303 Magma_tally4IVec
    "All",                                   // 304 Magma_tally4AllVec
    "Some",                                  // 305 Magma_tally4SomeVec
    "Overwrite",                             // 306 Magma_tally4OverwriteVec
    "", "", "", "",                          // 307-310
    "All",                                   // 311 Magma_tally4RangeAll
    "V",                                     // 312 Magma_tally4RangeV
    "I",                                     // 313 Magma_tally4RangeI
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
    "Forward",                               // 391: Magma_tally4Forward
    "Backward",                              // 392: Magma_tally4Backward
    "", "", "", "", "", "", "", "",          // 393-400
    "Columnwise",                            // 401: Magma_tally4Columnwise
    "Rowwise",                               // 402: Magma_tally4Rowwise
    "", "", "", "", "", "", "", ""           // 403-410
    // Remember to add a comma!
};

extern "C"
const char* lapack_const_tally4( int magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally42lapack_Min );
    assert( magma_tally4_const <= Magma_tally42lapack_Max );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_bool_const_tally4( magma_tally4_bool_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4False );
    assert( magma_tally4_const <= Magma_tally4True  );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_order_const_tally4( magma_tally4_order_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4RowMajor );
    assert( magma_tally4_const <= Magma_tally4ColMajor );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_trans_const_tally4( magma_tally4_trans_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoTrans   );
    assert( magma_tally4_const <= Magma_tally4ConjTrans );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_uplo_const_tally4 ( magma_tally4_uplo_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Upper );
    assert( magma_tally4_const <= Magma_tally4Full  );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_diag_const_tally4 ( magma_tally4_diag_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NonUnit );
    assert( magma_tally4_const <= Magma_tally4Unit    );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_side_const_tally4 ( magma_tally4_side_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Left  );
    assert( magma_tally4_const <= Magma_tally4BothSides );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_norm_const_tally4  ( magma_tally4_norm_t   magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4OneNorm     );
    assert( magma_tally4_const <= Magma_tally4RealMaxNorm );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_dist_const_tally4  ( magma_tally4_dist_t   magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4DistUniform );
    assert( magma_tally4_const <= Magma_tally4DistNormal );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_sym_const_tally4   ( magma_tally4_sym_t    magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4HermGeev );
    assert( magma_tally4_const <= Magma_tally4SymPosv  );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_pack_const_tally4  ( magma_tally4_pack_t   magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoPacking );
    assert( magma_tally4_const <= Magma_tally4PackAll   );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_vec_const_tally4   ( magma_tally4_vec_t    magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoVec );
    assert( magma_tally4_const <= Magma_tally4OverwriteVec );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_range_const_tally4 ( magma_tally4_range_t  magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4RangeAll );
    assert( magma_tally4_const <= Magma_tally4RangeI   );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_vect_const_tally4( magma_tally4_vect_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Q );
    assert( magma_tally4_const <= Magma_tally4P );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_direct_const_tally4( magma_tally4_direct_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Forward );
    assert( magma_tally4_const <= Magma_tally4Backward );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}

extern "C"
const char* lapack_storev_const_tally4( magma_tally4_storev_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Columnwise );
    assert( magma_tally4_const <= Magma_tally4Rowwise    );
    return magma_tally42lapack_const_tally4ants[ magma_tally4_const ];
}


// ----------------------------------------
// Convert magma_tally4 constants to clAmdBlas constants.

#ifdef HAVE_clAmdBlas
const int magma_tally42amdblas_constants[] =
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
    clAmdBlasRowMajor,      // 101: Magma_tally4RowMajor
    clAmdBlasColumnMajor,   // 102: Magma_tally4ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNoTrans,       // 111: Magma_tally4NoTrans
    clAmdBlasTrans,         // 112: Magma_tally4Trans
    clAmdBlasConjTrans,     // 113: Magma_tally4ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    clAmdBlasUpper,         // 121: Magma_tally4Upper
    clAmdBlasLower,         // 122: Magma_tally4Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasNonUnit,       // 131: Magma_tally4NonUnit
    clAmdBlasUnit,          // 132: Magma_tally4Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    clAmdBlasLeft,          // 141: Magma_tally4Left
    clAmdBlasRight,         // 142: Magma_tally4Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
clAmdBlasOrder       amdblas_order_const( magma_tally4_order_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4RowMajor );
    assert( magma_tally4_const <= Magma_tally4ColMajor );
    return (clAmdBlasOrder)     magma_tally42amdblas_constants[ magma_tally4_const ];
}

extern "C"
clAmdBlasTranspose   amdblas_trans_const( magma_tally4_trans_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoTrans   );
    assert( magma_tally4_const <= Magma_tally4ConjTrans );
    return (clAmdBlasTranspose) magma_tally42amdblas_constants[ magma_tally4_const ];
}

extern "C"
clAmdBlasUplo        amdblas_uplo_const ( magma_tally4_uplo_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Upper );
    assert( magma_tally4_const <= Magma_tally4Lower );
    return (clAmdBlasUplo)      magma_tally42amdblas_constants[ magma_tally4_const ];
}

extern "C"
clAmdBlasDiag        amdblas_diag_const ( magma_tally4_diag_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NonUnit );
    assert( magma_tally4_const <= Magma_tally4Unit    );
    return (clAmdBlasDiag)      magma_tally42amdblas_constants[ magma_tally4_const ];
}

extern "C"
clAmdBlasSide        amdblas_side_const ( magma_tally4_side_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Left  );
    assert( magma_tally4_const <= Magma_tally4Right );
    return (clAmdBlasSide)      magma_tally42amdblas_constants[ magma_tally4_const ];
}
#endif  // HAVE_clAmdBlas


// ----------------------------------------
// Convert magma_tally4 constants to Nvidia CUBLAS constants.

#ifdef HAVE_CUBLAS
const int magma_tally42cublas_constants[] =
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
    0,                      // 101: Magma_tally4RowMajor
    0,                      // 102: Magma_tally4ColMajor
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_OP_N,            // 111: Magma_tally4NoTrans
    CUBLAS_OP_T,            // 112: Magma_tally4Trans
    CUBLAS_OP_C,            // 113: Magma_tally4ConjTrans
    0, 0, 0, 0, 0, 0, 0,
    CUBLAS_FILL_MODE_UPPER, // 121: Magma_tally4Upper
    CUBLAS_FILL_MODE_LOWER, // 122: Magma_tally4Lower
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_DIAG_NON_UNIT,   // 131: Magma_tally4NonUnit
    CUBLAS_DIAG_UNIT,       // 132: Magma_tally4Unit
    0, 0, 0, 0, 0, 0, 0, 0,
    CUBLAS_SIDE_LEFT,       // 141: Magma_tally4Left
    CUBLAS_SIDE_RIGHT,      // 142: Magma_tally4Right
    0, 0, 0, 0, 0, 0, 0, 0
};

extern "C"
cublasOperation_t    cublas_trans_const_tally4 ( magma_tally4_trans_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoTrans   );
    assert( magma_tally4_const <= Magma_tally4ConjTrans );
    return (cublasOperation_t)  magma_tally42cublas_constants[ magma_tally4_const ];
}

extern "C"
cublasFillMode_t     cublas_uplo_const_tally4  ( magma_tally4_uplo_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Upper );
    assert( magma_tally4_const <= Magma_tally4Lower );
    return (cublasFillMode_t)   magma_tally42cublas_constants[ magma_tally4_const ];
}

extern "C"
cublasDiagType_t     cublas_diag_const_tally4  ( magma_tally4_diag_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NonUnit );
    assert( magma_tally4_const <= Magma_tally4Unit    );
    return (cublasDiagType_t)   magma_tally42cublas_constants[ magma_tally4_const ];
}

extern "C"
cublasSideMode_t     cublas_side_const_tally4  ( magma_tally4_side_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Left  );
    assert( magma_tally4_const <= Magma_tally4Right );
    return (cublasSideMode_t)   magma_tally42cublas_constants[ magma_tally4_const ];
}
#endif  // HAVE_CUBLAS


// ----------------------------------------
// Convert magma_tally4 constants to CBLAS constants.
// We assume that magma_tally4 constants are consistent with cblas constants,
// so verify that with asserts.

#ifdef HAVE_CBLAS
extern "C"
enum CBLAS_ORDER     cblas_order_const  ( magma_tally4_order_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4RowMajor );
    assert( magma_tally4_const <= Magma_tally4ColMajor );
    assert( (int)Magma_tally4RowMajor == CblasRowMajor );
    return (enum CBLAS_ORDER)     magma_tally4_const;
}

extern "C"
enum CBLAS_TRANSPOSE cblas_trans_const  ( magma_tally4_trans_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NoTrans   );
    assert( magma_tally4_const <= Magma_tally4ConjTrans );
    assert( (int)Magma_tally4NoTrans == CblasNoTrans );
    return (enum CBLAS_TRANSPOSE) magma_tally4_const;
}

extern "C"
enum CBLAS_UPLO      cblas_uplo_const   ( magma_tally4_uplo_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Upper );
    assert( magma_tally4_const <= Magma_tally4Lower );
    assert( (int)Magma_tally4Upper == CblasUpper );
    return (enum CBLAS_UPLO)      magma_tally4_const;
}

extern "C"
enum CBLAS_DIAG      cblas_diag_const   ( magma_tally4_diag_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4NonUnit );
    assert( magma_tally4_const <= Magma_tally4Unit    );
    assert( (int)Magma_tally4Unit == CblasUnit );
    return (enum CBLAS_DIAG)      magma_tally4_const;
}

extern "C"
enum CBLAS_SIDE      cblas_side_const   ( magma_tally4_side_t magma_tally4_const )
{
    assert( magma_tally4_const >= Magma_tally4Left  );
    assert( magma_tally4_const <= Magma_tally4Right );
    assert( (int)Magma_tally4Left == CblasLeft );
    return (enum CBLAS_SIDE)      magma_tally4_const;
}
#endif  // HAVE_CBLAS
