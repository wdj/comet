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
  int data_type_id;
  void* __restrict__ data;
} Vectors;

/*===========================================================================*/
/*---Functions---*/

void Vectors_create(Vectors* vectors,
                    int data_type_id,
                    int num_field,
                    int num_vector_local,
                    Env* env);

void Vectors_destroy(Vectors* vectors, Env* env);

/*---------------------------------------------------------------------------*/

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
/*---Accessors---*/

__attribute__((unused)) static void Vectors_set(Vectors* vectors,
                                                int field,
                                                int vector_local,
                                                double value,
                                                Env* env) {
  /*---WARNING: this is inefficient when called in a high-tripcount loop---*/
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);

  switch (vectors->data_type_id) {
    case DATA_TYPE_ID_FLOAT:
      Vectors_float_set(vectors, field, vector_local, (Float_t)value, env);
      return;
    default:
      Insist(Bool_false ? "Invalid data_type_id." : 0);
  }
}

/*---------------------------------------------------------------------------*/

__attribute__((unused)) static Float_t Vectors_float_get(Vectors* const vectors,
                                                         int field,
                                                         int vector_local,
                                                         Env* env) {
  Assert(vectors);
  Assert(field >= 0);
  Assert(field < vectors->num_field);
  Assert(vector_local >= 0);
  Assert(vector_local < vectors->num_vector_local);
  Assert(env);

  return ((Float_t*)(vectors->data))[field + vectors->num_field * vector_local];
}

/*===========================================================================*/

#endif /*---_vectors_h_---*/

/*---------------------------------------------------------------------------*/
