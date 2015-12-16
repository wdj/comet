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
  int num_field;
  int num_field_local;
  int num_vector;
  int num_vector_local;
  size_t num_dataval_local;
  int num_field_dataval;
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
/*---Accessors---*/

#if 0
static void GMVectors_bitd128_set(GMVectors* vectors,
                                  int field_local,
                                  int vector_local,
                                  int value,
                                  GMEnv* env) {
  GMAssert(vectors);
  GMAssert(field_local >= 0);
  GMAssert(field_local < vectors->num_field_local);
  GMAssert(vector_local >= 0);
  GMAssert(vector_local < vectors->num_vector_local);

  GMAssert(value >= 0);
  GMAssert(value < 4);

  // the field is mapped to which GMBits128 entry is used,
  // which of the two data value is used (2 choices),
  // and the location in that data value (32 choices).
  int field_index2 = field_local / 64;
  int field_index1 = (field_local % 64) / 2;
  int field_index0 = field_local / 128;

// MUST modify to account for the vector number
  ((GMBits128*)(vectors->data))[field_index2].data[field_index1]

    |= value << field_index0;



  ((GMFloat*)(vectors->data))[field_local + vectors->num_field_local *
                              vector_local] = value;
}

/*---------------------------------------------------------------------------*/
#endif

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
      field_dataval_num + vectors->num_field_dataval * vector_local;

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
