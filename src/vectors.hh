//-----------------------------------------------------------------------------
/*!
 * \file   vectors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Vectors pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#ifndef _gm_vectors_hh_
#define _gm_vectors_hh_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================
/*---Struct declaration---*/

typedef struct {
  /*---Logical sizes---*/
  int num_field;
  int num_field_local;
  size_t num_field_active;
  int num_vector;
  int num_vector_local;
  /*---Stored sizes---*/
  int num_bits_per_val;
  int num_bits_per_packedval;
  int num_val_per_packedval;
  int num_packedval_field_local;
  size_t num_packedval_local;
  /*---Other---*/
  int data_type_id;
  int pad1;
  void* __restrict__ data;
  size_t data_size;
  bool has_buf;
  GMMirroredBuf buf;
  GMDecompMgr* dm;
} GMVectors;

//=============================================================================
/*---Null object---*/

GMVectors GMVectors_null(void);

//=============================================================================
/*---Vectors pseudo-constructor---*/

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      GMEnv* env);

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               GMEnv* env);

//=============================================================================
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, GMEnv* env);

//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(GMMirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env);

//=============================================================================
/*---Accessors: Float---*/

static GMFloat* GMVectors_float_ptr(GMVectors* const vectors,
                                   int field_local,
                                   int vector_local,
                                   GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(env);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(vectors->data)) + (field_local +
                        vectors->num_field_local*(size_t)vector_local);
}

//-----------------------------------------------------------------------------

static void GMVectors_float_set(GMVectors* vectors,
                                int field_local,
                                int vector_local,
                                GMFloat value,
                                GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT);

  *(GMVectors_float_ptr(vectors, field_local, vector_local, env)) = value;

  //((GMFloat*)(vectors->data))[field_local +
  //                      vectors->num_field_local*(size_t)vector_local] = value;
}

//-----------------------------------------------------------------------------

static GMFloat GMVectors_float_get_from_index(GMVectors* const vectors,
                                              size_t index,
                                              GMEnv* env) {
  GMAssert(vectors);
  GMAssert(index >= 0);
  GMAssert(index < vectors->num_vector_local*(size_t)vectors->num_field_local);
  GMAssert(env);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(vectors->data))[index];
}

//-----------------------------------------------------------------------------

static GMFloat GMVectors_float_get(GMVectors* const vectors,
                                   int field_local,
                                   int vector_local,
                                   GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(env);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_FLOAT);

  return GMVectors_float_get_from_index(vectors,
    field_local + vectors->num_field_local*(size_t)vector_local, env);
}

//=============================================================================
/*---Accessors: Bits2, Bits2x64---*/

static GMBits2 GMVectors_bits2_get(GMVectors* vectors,
                                   int field_local,
                                   int vector_local,
                                   GMEnv* env) {
  /*---This function gets a single 2-bit value---*/
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_BITS2);

  /*---The field address is expressible as a tuple:
       which GMBits2x64 value,
       which of the 2 (size2) data entries,
       which of the 32 (size1) 2-bit (size0) fields in the data entry
  ---*/

  const int size0 = 2;
  const int size1 = 32;
  const int size2 = 2;

  int field_index0 = field_local % size1;
  int field_index1 = (field_local / size1) % size2;
  size_t field_index2 = field_local / (size1 * size2);

  GMBits1_2x64* const __restrict__ address =
      &(((GMBits2x64*)(vectors->data))[field_index2 +
                                       vectors->num_packedval_field_local *
                                       (size_t)vector_local]
            .data[field_index1]);

  return (GMBits2)(((*address) >> (size0 * field_index0)) & ((GMBits1_2x64)3));
}

//-----------------------------------------------------------------------------

static void GMVectors_bits2_set(GMVectors* vectors,
                                int field_local,
                                int vector_local,
                                GMBits2 value,
                                GMEnv* env) {
  /*---This function sets a single 2-bit value---*/
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(value >= 0 && value < (1 << GM_BITS2_MAX_VALUE_BITS));
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_BITS2);

  /*---The field address is expressible as a tuple:
       which GMBits2x64 value,
       which of the 2 (size2) data entries,
       which of the 32 (size1) 2-bit (size0) fields in the data entry
  ---*/

  const int size0 = 2;
  const int size1 = 32;
  const int size2 = 2;

  int field_index0 = field_local % size1;
  int field_index1 = (field_local / size1) % size2;
  size_t field_index2 = field_local / (size1 * size2);

  GMBits1_2x64* const __restrict__ address =
      &(((GMBits2x64*)(vectors->data))[field_index2 +
                                       vectors->num_packedval_field_local *
                                       (size_t)vector_local]
            .data[field_index1]);

  *address &= ~(((GMBits1_2x64)3) << (size0 * field_index0));
  *address |= ((GMBits1_2x64)value) << (size0 * field_index0);

  GMAssert(value ==
           GMVectors_bits2_get(vectors, field_local, vector_local, env));
}

//-----------------------------------------------------------------------------

static GMBits2x64* GMVectors_bits2x64_ptr(GMVectors* vectors,
                                          int packedval_field_local,
                                          int vector_local,
                                          GMEnv* env) {
  /*---This function accesses an entire packed value containing 2-bit values---*/
  GMAssert(vectors);
  GMAssert(packedval_field_local >= 0);
  GMAssert(packedval_field_local < vectors->num_packedval_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_BITS2);

  const size_t index = packedval_field_local +
    vectors->num_packedval_field_local*(size_t)vector_local;

  return ((GMBits2x64*)(vectors->data)) + index;
}

//-----------------------------------------------------------------------------

static void GMVectors_bits2x64_set(GMVectors* vectors,
                                   int packedval_field_local,
                                   int vector_local,
                                   GMBits2x64 value,
                                   GMEnv* env) {
  /*---This function sets an entire packed value containing 2-bit values---*/
  GMAssert(vectors);
  GMAssert(packedval_field_local >= 0);
  GMAssert(packedval_field_local < vectors->num_packedval_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_BITS2);

  const size_t index = packedval_field_local +
    vectors->num_packedval_field_local*(size_t)vector_local;

  ((GMBits2x64*)(vectors->data))[index].data[0] = value.data[0];
  ((GMBits2x64*)(vectors->data))[index].data[1] = value.data[1];
}

//-----------------------------------------------------------------------------

static GMBits2x64 GMVectors_bits2x64_get(GMVectors* vectors,
                                         int packedval_field_local,
                                         int vector_local,
                                         GMEnv* env) {
  /*---This function gets an entire packed value containing 2-bit values---*/
  GMAssert(vectors);
  GMAssert(packedval_field_local >= 0);
  GMAssert(packedval_field_local < vectors->num_packedval_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(GMEnv_data_type_vectors(env) == GM_DATA_TYPE_BITS2);

  const size_t index = packedval_field_local +
    vectors->num_packedval_field_local*(size_t)vector_local;

  const GMBits2x64 value = ((GMBits2x64*)(vectors->data))[index];

  return value;
}

//=============================================================================
/*---Misc---*/

static size_t GMVectors_num_local_required(size_t num_vector_active,
                                           GMEnv* const env) {

  const int num_proc_vector = GMEnv_num_proc_vector_i(env);

  const size_t nvl_1 = gm_ceil_i8(num_vector_active, num_proc_vector);

  const bool need_divisible_by_6 = GMEnv_num_way(env) == GM_NUM_WAY_3 &&
                                   GMEnv_all2all(env) && num_proc_vector > 2;

  const size_t num_vector_local = need_divisible_by_6 ?
                                  gm_ceil_i8(nvl_1, 6)*6 : nvl_1;

  return num_vector_local;
}

//=============================================================================

#endif /*---_gm_vectors_hh_---*/

//-----------------------------------------------------------------------------
