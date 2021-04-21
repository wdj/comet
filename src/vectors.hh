//-----------------------------------------------------------------------------
/*!
 * \file   vectors.hh
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the set of vectors taken as input to the methods.
 */
//-----------------------------------------------------------------------------
/*-----------------------------------------------------------------------------

Copyright 2020, UT-Battelle, LLC

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

-----------------------------------------------------------------------------*/

#ifndef _COMET_VECTORS_HH_
#define _COMET_VECTORS_HH_

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
// Struct declaration.

struct GMVectors {
  // Logical sizes.
  int num_field;
  int num_field_local;
  size_t num_field_active;
  int num_vector;
  int num_vector_local;
  // Stored sizes.
  int num_packedfield_local;
  size_t num_packedfield_vector_local;
  // Other.
  int data_type_id;
  int pad1;
  void* __restrict__ data;
  size_t data_size;
  bool has_buf_;
  bool has_buf() const {return has_buf_;}
  MirroredBuf* buf_;
  MirroredBuf* buf() const {
    COMET_INSIST(has_buf_);
    return buf_;
  };
  GMDecompMgr* dm;
};

//=============================================================================
// Null object.

GMVectors GMVectors_null(void);

//=============================================================================
// Vectors pseudo-constructor.

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      CEnv* env);

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               CEnv* env);

//-----------------------------------------------------------------------------

void GMVectors_initialize(GMVectors* vectors, CEnv* env);

void GMVectors_initialize_pad(GMVectors* vectors, CEnv* env);

//=============================================================================
// Vectors pseudo-destructor.

void GMVectors_destroy(GMVectors* vectors, CEnv* env);

//=============================================================================

size_t GMVectors_cksum(GMVectors* vectors, CEnv* env);

//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(MirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       CEnv* env);

//=============================================================================
// Accessors: Float.

static GMFloat* GMVectors_float_ptr(GMVectors* const vectors,
                                   int field_local,
                                   int vector_local,
                                   CEnv* env) {
  COMET_ASSERT(vectors);
  COMET_ASSERT(field_local >= 0);
  COMET_ASSERT(field_local < vectors->num_field_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(vectors->data)) + (field_local +
                        vectors->num_field_local*(size_t)vector_local);
}

//-----------------------------------------------------------------------------

static void GMVectors_float_set(GMVectors* vectors,
                                int field_local,
                                int vector_local,
                                GMFloat value,
                                CEnv* env) {
  COMET_ASSERT(vectors);
  COMET_ASSERT(field_local >= 0);
  COMET_ASSERT(field_local < vectors->num_field_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_FLOAT);

  *(GMVectors_float_ptr(vectors, field_local, vector_local, env)) = value;
}

//-----------------------------------------------------------------------------

static GMFloat GMVectors_float_get_from_index(const GMVectors* const vectors,
                                              size_t index,
                                              CEnv* env) {
  COMET_ASSERT(vectors);
  //COMET_ASSERT(index >= 0);
  COMET_ASSERT(index < vectors->num_vector_local*(size_t)vectors->num_field_local);
  COMET_ASSERT(env);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_FLOAT);

  return ((GMFloat*)(vectors->data))[index];
}

//-----------------------------------------------------------------------------

static GMFloat GMVectors_float_get(const GMVectors* const vectors,
                                   int field_local,
                                   int vector_local,
                                   CEnv* env) {
  COMET_ASSERT(vectors);
  COMET_ASSERT(field_local >= 0);
  COMET_ASSERT(field_local < vectors->num_field_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_FLOAT);

  return GMVectors_float_get_from_index(vectors,
    field_local + vectors->num_field_local*(size_t)vector_local, env);
}

//=============================================================================
// Accessors: Bits2, Bits2x64.

static GMBits2 GMVectors_bits2_get(const GMVectors* vectors,
                                   int field_local,
                                   int vector_local,
                                   CEnv* env) {
  // This function gets a single 2-bit value.
  COMET_ASSERT(vectors);
  COMET_ASSERT(field_local >= 0);
  COMET_ASSERT(field_local < vectors->num_field_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_BITS2);

  /*---The field address is expressible as a tuple:
       which GMBits2x64 value,
       which of the 2 (size2) vector data entries,
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
                                       vectors->num_packedfield_local *
                                       (size_t)vector_local]
            .data[field_index1]);

  return (GMBits2)(((*address) >> (size0 * field_index0)) & ((GMBits1_2x64)3));
}

//-----------------------------------------------------------------------------

static void GMVectors_bits2_set(GMVectors* vectors,
                                int field_local,
                                int vector_local,
                                GMBits2 value,
                                CEnv* env) {
  // This function sets a single 2-bit value.
  COMET_ASSERT(vectors);
  COMET_ASSERT(field_local >= 0);
  COMET_ASSERT(field_local < vectors->num_field_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(value+1 >= 1 && value < (1 << GM_BITS2_MAX_VALUE_BITS));
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_BITS2);

  /*---The field address is expressible as a tuple:
       which GMBits2x64 value,
       which of the 2 (size2) vector data entries,
       which of the 32 (size1) 2-bit (size0) fields in the data entry
  ---*/

  const int size0 = 2;
  const int size1 = 32;
  const int size2 = 2;

  const int field_index0 = field_local % size1;
  const int field_index1 = (field_local / size1) % size2;
  const size_t field_index2 = field_local / (size1 * size2);

  GMBits1_2x64* const __restrict__ address =
      &(((GMBits2x64*)(vectors->data))[field_index2 +
                                       vectors->num_packedfield_local *
                                       (size_t)vector_local]
            .data[field_index1]);

  *address &= ~(((GMBits1_2x64)3) << (size0 * field_index0));
  *address |= ((GMBits1_2x64)value) << (size0 * field_index0);

  COMET_ASSERT(value ==
           GMVectors_bits2_get(vectors, field_local, vector_local, env));
}

//-----------------------------------------------------------------------------

static GMBits2x64* GMVectors_bits2x64_ptr(GMVectors* vectors,
                                          int packedfield_local,
                                          int vector_local,
                                          CEnv* env) {
  // This function accesses an entire packed value containing 2-bit values.
  COMET_ASSERT(vectors);
  COMET_ASSERT(packedfield_local >= 0);
  COMET_ASSERT(packedfield_local < vectors->num_packedfield_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_BITS2);

  const size_t index = packedfield_local +
    vectors->num_packedfield_local*(size_t)vector_local;

  return ((GMBits2x64*)(vectors->data)) + index;
}

//-----------------------------------------------------------------------------

__host__ __device__
static size_t GMVectors_index(int packedfield_local,
                              int vector_local,
                              int num_packedfield_local) {
  // This function accesses an entire packed value.

  const size_t index = packedfield_local +
    num_packedfield_local*(size_t)vector_local;

  return index;
}

//-----------------------------------------------------------------------------

static void GMVectors_bits2x64_set(GMVectors* vectors,
                                   int packedfield_local,
                                   int vector_local,
                                   GMBits2x64 value,
                                   CEnv* env) {
  // This function sets an entire packed value containing 2-bit values.
  COMET_ASSERT(vectors);
  COMET_ASSERT(packedfield_local >= 0);
  COMET_ASSERT(packedfield_local < vectors->num_packedfield_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_BITS2);

  const size_t index = packedfield_local +
    vectors->num_packedfield_local*(size_t)vector_local;

  ((GMBits2x64*)(vectors->data))[index].data[0] = value.data[0];
  ((GMBits2x64*)(vectors->data))[index].data[1] = value.data[1];
}

//-----------------------------------------------------------------------------

static GMBits2x64 GMVectors_bits2x64_get(const GMVectors* vectors,
                                         int packedfield_local,
                                         int vector_local,
                                         CEnv* env) {
  // This function gets an entire packed value containing 2-bit values.
  COMET_ASSERT(vectors);
  COMET_ASSERT(packedfield_local >= 0);
  COMET_ASSERT(packedfield_local < vectors->num_packedfield_local);
  COMET_ASSERT(vector_local >= 0);
  COMET_ASSERT(vector_local < vectors->num_vector_local);
  COMET_ASSERT(env->data_type_vectors() == GM_DATA_TYPE_BITS2);

  const size_t index = packedfield_local +
    vectors->num_packedfield_local*(size_t)vector_local;

  const GMBits2x64 value = ((GMBits2x64*)(vectors->data))[index];

  return value;
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------

#endif // _COMET_VECTORS_HH_

//-----------------------------------------------------------------------------
