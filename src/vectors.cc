//-----------------------------------------------------------------------------
/*!
 * \file   vectors.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Vectors pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "stdlib.h"
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "mirrored_buf.hh"
#include "linalg.hh"
#include "vectors.hh"

//=============================================================================
/*---Null object---*/

GMVectors GMVectors_null() {
  GMVectors result;
  memset((void*)&result, 0, sizeof(GMVectors));
  return result;
}

//=============================================================================
/*---Set unused (pad) vector entries to zero---*/

void GMVectors_initialize_pad_(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);

  /*---Ensure final pad words/bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  const size_t pfl_min = vectors->dm->num_field_active_local /
                         vectors->dm->num_field_per_packedfield;
  const size_t pfl_max = vectors->dm->num_packedfield_local;

  switch (vectors->data_type_id) {
    case GM_DATA_TYPE_FLOAT: {
      GMFloat zero = 0;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          GMVectors_float_set(vectors, pfl, vl, zero, env);
        }
      }
    } break;
    case GM_DATA_TYPE_BITS2: {
      GMBits2x64 zero = GMBits2x64_null();
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          /*---Doesn't hurt to set partly active words to zero here---*/
          GMVectors_bits2x64_set(vectors, pfl, vl, zero, env);
        }
      }
    } break;
    default:
      GMInsist(false ? "Invalid data type." : 0);
  } /*---switch---*/
}

//=============================================================================
/*---Vectors pseudo-constructor---*/

void GMVectors_create_imp_(GMVectors* vectors,
                           int data_type_id,
                           GMDecompMgr* dm,
                           GMEnv* env) {
  GMInsist(vectors && dm && env);

  vectors->data_type_id = data_type_id;
  vectors->dm = dm;
  vectors->num_field = dm->num_field;
  vectors->num_field_active = dm->num_field_active;
  vectors->num_field_local = dm->num_field_local;
  vectors->num_vector = dm->num_vector;
  vectors->num_vector_local = dm->num_vector_local;

  // Set element sizes

  vectors->num_bits_per_val = dm->num_bits_per_field;
  vectors->num_bits_per_packedval = dm->num_bits_per_packedfield;
  vectors->num_val_per_packedval = dm->num_field_per_packedfield;

  const int bits_per_byte = 8;

  // Allocation size for vector storage

  vectors->num_packedval_field_local = dm->num_packedfield_local;

  vectors->num_packedval_local =
      vectors->num_packedval_field_local * dm->num_vector_local;

  vectors->data_size = vectors->num_packedval_local *
                       (dm->num_bits_per_packedfield / bits_per_byte);

  // Set up vector storage, mirrored buffer

  vectors->buf = GMMirroredBuf_null();

  if (vectors->has_buf) {
    GMMirroredBuf_create(&(vectors->buf),vectors->num_packedval_field_local,
                         vectors->num_vector_local, env);
    vectors->data = vectors->buf.h; // alias vector storage to buf
  } else {
    vectors->data = gm_malloc(vectors->data_size, env);
  }

  // Set pad entries to zero

  GMVectors_initialize_pad_(vectors, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      GMEnv* env) {
  GMInsist(vectors && dm && env);

  *vectors = GMVectors_null();

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = false;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               GMEnv* env) {
  GMInsist(vectors && dm && env);

  *vectors = GMVectors_null();

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//=============================================================================
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(vectors->data || !GMEnv_is_proc_active(env));

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  if (vectors->has_buf) {
    GMMirroredBuf_destroy(&vectors->buf, env);
  } else {
    gm_free(vectors->data, vectors->data_size, env);
  }

  *vectors = GMVectors_null();
}


//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(GMMirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env) {
  GMInsist(vectors && vectors_buf && env);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  /*---Copy vectors into GPU buffers if needed---*/

  switch (GMEnv_metric_type(env)) {
    case GM_METRIC_TYPE_CZEK: {
#pragma omp parallel for collapse(2)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          GMMirroredBuf_elt<GMFloat>(vectors_buf, fl, i) =
            GMVectors_float_get(vectors, fl, i, env);
        }
      }
    } break;
    case GM_METRIC_TYPE_CCC: {
#pragma omp parallel for collapse(2)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedval_field_local; ++fl) {
          GMMirroredBuf_elt<GMBits2x64>(vectors_buf, fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    default:
      GMInsistInterface(env, false ? "Unimplemented." : 0);
  } /*---case---*/
}

//-----------------------------------------------------------------------------
