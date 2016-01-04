/*---------------------------------------------------------------------------*/
/*!
 * \file   vectors.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Vectors pseudo-class, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _vectors_h_
#define _vectors_h_

#include "env.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  /*---Logical sizes---*/
  int num_field;
  int num_field_local;
  int num_vector;
  int num_vector_local;
  /*---Stored sizes---*/
  int num_bits_per_val;
  int num_bits_per_packedval;
  int num_packedval_field_local;
  size_t num_packedval_local;
  /*---Other---*/
  int data_type_id;
  void* __restrict__ data;
} GMVectors;

/*===========================================================================*/
/*---Null object---*/

GMVectors GMVectors_null(void);

/*===========================================================================*/
/*---Vectors pseudo-constructor---*/

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      int num_field,
                      int num_vector_local,
                      GMEnv* env);

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, GMEnv* env);

/*===========================================================================*/
/*---Accessors: Float---*/

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

  ((GMFloat*)(vectors->data))[field_local + vectors->num_field_local *
                              vector_local] = value;
}

/*---------------------------------------------------------------------------*/

static GMFloat GMVectors_float_get_from_index(GMVectors* const vectors,
                                              int index,
                                              GMEnv* env) {
  GMAssert(vectors);
  GMAssert(index >= 0);
  GMAssert(index < vectors->num_vector_local * vectors->num_field_local);
  GMAssert(env);

  return ((GMFloat*)(vectors->data))[index];
}

/*---------------------------------------------------------------------------*/

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

  return GMVectors_float_get_from_index(
      vectors, field_local + vectors->num_field_local * vector_local, env);
}

/*===========================================================================*/
/*---Accessors: Bits2, Bits2x64---*/

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

  GMBits1_2x64* const __restrict__ address = &( 
    ((GMBits2x64*)(vectors->data))[
      field_index2 + vectors->num_packedval_field_local *
      vector_local].data[field_index1] );

  *address &= ~ ( ((GMBits1_2x64)3) << (size0*field_index0) );
  *address |=     value             << (size0*field_index0);
}

/*---------------------------------------------------------------------------*/

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

  GMBits1_2x64* const __restrict__ address = &( 
    ((GMBits2x64*)(vectors->data))[
      field_index2 + vectors->num_packedval_field_local *
      vector_local].data[field_index1] );

  return (GMBits2)(( (*address) >> (size0*field_index0) ) & ((GMBits1_2x64)3));
}

/*---------------------------------------------------------------------------*/

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

  const int index = packedval_field_local + vectors->num_packedval_field_local *
                    vector_local;

  ((GMBits2x64*)(vectors->data))[index].data[0] = value.data[0];
  ((GMBits2x64*)(vectors->data))[index].data[1] = value.data[1];
}

/*===========================================================================*/
/*---Accessors: Bits1---*/
//---(design is not complete)

static void GMVectors_bits1x64_set(GMVectors* vectors,
                                   int packedval_field_local,
                                   int vector_local,
                                   GMBits1x64 value,
                                   GMEnv* env) {
  GMAssert(vectors);
  GMAssert(packedval_field_local >= 0);
  GMAssert(packedval_field_local < vectors->num_packedval_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);

  const int index = packedval_field_local + vectors->num_packedval_field_local *
                    vector_local;

  ((GMBits1x64*)(vectors->data))[index] = value;
}

/*---------------------------------------------------------------------------*/



static int GMVectors_bit_dataval_num(GMVectors* vectors,
                                     int field_local,
                                     int vector_local,
                                     GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);

  const int field_dataval_num = field_local / (8 * sizeof(GMBits));

  const int dataval_num =
      field_dataval_num + vectors->num_packedval_field_local * vector_local;

  return dataval_num;
}

/*---------------------------------------------------------------------------*/

static void GMVectors_bit_set(GMVectors* vectors,
                              int field_local,
                              int vector_local,
                              _Bool value,
                              GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);

  const int bit_num = field_local % (8 * sizeof(GMBits));

  const int dataval_num =
      GMVectors_bit_dataval_num(vectors, field_local, vector_local, env);

  GMAssert(sizeof(GMULInt) == sizeof(GMBits));
  GMAssert(sizeof(GMBits) == 8);

  const GMULInt mask = ((GMULInt)1) << bit_num;
  GMULInt* const p = &(((GMULInt*)(vectors->data))[dataval_num]);
  *p = ((*p) & (~mask)) | mask;
}

/*---------------------------------------------------------------------------*/

static _Bool GMVectors_bit_get(GMVectors* const vectors,
                               int field_local,
                               int vector_local,
                               GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);
  GMAssert(env);

  const int bit_num = field_local % (8 * sizeof(GMBits));

  const int dataval_num =
      GMVectors_bit_dataval_num(vectors, field_local, vector_local, env);

  GMAssert(sizeof(GMULInt) == sizeof(GMBits));
  GMAssert(sizeof(GMBits) == 8);

  _Bool result = 0;

  GMULInt* const p = &(((GMULInt*)(vectors->data))[dataval_num]);
  result = ((*p) >> bit_num) & ((GMULInt)1) ? GM_BOOL_TRUE : GM_BOOL_FALSE;

  return result;
}

/*---------------------------------------------------------------------------*/

static _Bool GMVectors_bit_get_from_index(GMVectors* const vectors,
                                          int index,
                                          GMEnv* env) {
  GMAssert(vectors);
  GMAssert(index >= 0);
  GMAssert(index < vectors->num_vector_local * vectors->num_field_local);
  GMAssert(env);

  return GMVectors_bit_get(vectors, index % vectors->num_field_local,
                           index / vectors->num_field_local, env);
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*===========================================================================*/

#endif /*---_vectors_h_---*/

/*---------------------------------------------------------------------------*/
