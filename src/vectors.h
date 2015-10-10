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

/*===========================================================================*/
/*---Struct declaration---*/

typedef struct {
  int num_field;
  int num_vector;
  int num_vector_local;
  int num_vector_local_max;
  int num_field_dataval;
  int data_type_id;
  void* __restrict__ data;
} Vectors;

/*===========================================================================*/
/*---Null object---*/

Vectors Vectors_null(void);

/*===========================================================================*/
/*---Vectors pseudo-constructor---*/

void Vectors_create(Vectors* vectors,
                    int data_type_id,
                    int num_field,
                    int num_vector_local,
                    Env* env);

/*===========================================================================*/
/*---Vectors pseudo-destructor---*/

void Vectors_destroy(Vectors* vectors, Env* env);

/*===========================================================================*/
/*---Accessors---*/

static void Vectors_float_set(Vectors* vectors,
                              int field,
                              int vector_local,
                              Float_t value,
                              Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);

  ((Float_t*)(vectors->data))[field + vectors->num_field * vector_local] =
      value;
}

/*---------------------------------------------------------------------------*/

static Float_t Vectors_float_get_from_index(Vectors* const vectors,
                                            int index,
                                            Env* env) {
  Assert(vectors);
  Assert(index >= 0);
  Assert(index < vectors->num_vector_local * vectors->num_field);
  Assert(env);

  return ((Float_t*)(vectors->data))[index];
}

/*---------------------------------------------------------------------------*/

static Float_t Vectors_float_get(Vectors* const vectors,
                                 int field,
                                 int vector_local,
                                 Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);
  Assert(env);

  return Vectors_float_get_from_index( vectors,
              field + vectors->num_field * vector_local, env );
}

/*---------------------------------------------------------------------------*/

static int Vectors_bit_dataval_num(Vectors* vectors,
                                   int field,
                                   int vector_local,
                                   Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);

  const int field_dataval_num = field / ( 8 * sizeof(Bits_t) );

  const int dataval_num = field_dataval_num + vectors->num_field_dataval *
                          vector_local;

  return dataval_num;
}

/*---------------------------------------------------------------------------*/

static void Vectors_bit_set(Vectors* vectors,
                            int field,
                            int vector_local,
                            Bool_t value,
                            Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);

  const int bit_num = field % ( 8 * sizeof(Bits_t) );

  const int dataval_num = Vectors_bit_dataval_num( vectors, field,
                                                   vector_local, env );

  Assert( sizeof(ULInt_t) == sizeof(Bits_t) );
  Assert( sizeof(Bits_t) == 8 );

  const ULInt_t mask = ( (ULInt_t) 1 ) << bit_num;
  ULInt_t* const p = &(((ULInt_t*)(vectors->data))[dataval_num]);
  *p = ( ( *p ) & ( ~ mask ) ) | mask;
}

/*---------------------------------------------------------------------------*/

static Bool_t Vectors_bit_get(Vectors* const vectors,
                              int field,
                              int vector_local,
                              Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);
  Assert(env);

  const int bit_num = field % ( 8 * sizeof(Bits_t) );

  const int dataval_num = Vectors_bit_dataval_num( vectors, field,
                                                   vector_local, env );

  Assert( sizeof(ULInt_t) == sizeof(Bits_t) );
  Assert( sizeof(Bits_t) == 8 );

  Bool_t result = 0;

  ULInt_t* const p = &(((ULInt_t*)(vectors->data))[dataval_num]);
  result = ( ( *p ) >> bit_num ) & ((ULInt_t)1) ? Bool_true : Bool_false;

  return result;
}

/*---------------------------------------------------------------------------*/

static Bool_t Vectors_bit_get_from_index(Vectors* const vectors,
                                         int index,
                                         Env* env) {
  Assert(vectors);
  Assert(index >= 0);
  Assert(index < vectors->num_vector_local * vectors->num_field);
  Assert(env);

  return Vectors_bit_get( vectors, index % vectors->num_field,
                                   index / vectors->num_field, env );
}

/*===========================================================================*/

#endif /*---_vectors_h_---*/

/*---------------------------------------------------------------------------*/
