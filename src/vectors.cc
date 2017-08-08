//-----------------------------------------------------------------------------
/*!
 * \file   vectors.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Vectors pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "env.hh"
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

void GMVectors_initialize_pad(GMVectors* vectors,
                              GMEnv* env) {
  GMInsist(vectors && env);

  /*---Ensure final pad bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  const int nfl = vectors->num_field_local;
  const int fl_min = GMEnv_proc_num_field(env) == GMEnv_num_proc_field(env)-1 ?
    nfl - (vectors->num_field - vectors->num_field_active) : nfl;
  GMInsist(fl_min >= 0);

  switch (vectors->data_type_id) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      GMFloat zero = 0;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (int fl = fl_min; fl < nfl; ++fl) {
          GMVectors_float_set(vectors, fl, vl, zero, env);
        }
      }
      /*---NO-OP---*/
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
      const int pvfl_min = fl_min / vectors->num_val_per_packedval;
      GMBits2x64 zero = GMBits2x64_null();
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (int pvfl=pvfl_min; pvfl<vectors->num_packedval_field_local;
             ++pvfl) {
          /*---Doesn't hurt to set whole words to zero here---*/
          GMVectors_bits2x64_set(vectors, pvfl, vl, zero, env);
        }
      }
    } break;
    /*--------------------*/
    default:
      GMInsist(false ? "Invalid data type." : 0);
  } /*---switch---*/
}

//=============================================================================
/*---Vectors pseudo-constructor---*/

void GMVectors_create_imp_(GMVectors* vectors,
                           int data_type_id,
                           int num_field,
                           size_t num_field_active,
                           int num_vector_local,
                           GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(num_field >= 0);
  GMInsist(num_field_active >= 0);
  GMInsist(num_field_active <= (size_t)num_field);
  GMInsist(num_vector_local >= 0);

  GMInsistInterface(env,
           num_field % GMEnv_num_proc_field(env) == 0
               ? "num_proc_field must exactly divide the total number of fields"
               : 0);

  vectors->data_type_id = data_type_id;
  vectors->num_field = num_field;
  vectors->num_field_active = num_field_active;
  vectors->num_field_local = num_field / GMEnv_num_proc_field(env);
  vectors->num_vector_local = num_vector_local;

  /*---Compute global values---*/

  const int num_block = GMEnv_num_block_vector(env);

  const size_t num_vector_bound = vectors->num_vector_local * 
                            (size_t)num_block * (size_t)GMEnv_num_proc_repl(env);
  GMInsist(num_vector_bound == (size_t)(int)num_vector_bound
    ? "Vector count too large to store in 32-bit int; please modify code." : 0);

  int mpi_code = 0;
  mpi_code = MPI_Allreduce(&(vectors->num_vector_local), &(vectors->num_vector),
                           1, MPI_INT, MPI_SUM, GMEnv_mpi_comm_vector(env));
  GMInsist(mpi_code == MPI_SUCCESS);
  GMInsist((size_t)(vectors->num_vector) == num_vector_bound);
  vectors->num_vector /= GMEnv_num_proc_repl(env);

  /*---Set element sizes---*/

  const int bits_per_byte = 8;

  switch (vectors->data_type_id) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
      vectors->num_bits_per_val = bits_per_byte * sizeof(GMFloat);
      vectors->num_bits_per_packedval = bits_per_byte * sizeof(GMFloat);
    } break;
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
      vectors->num_bits_per_val = GM_BITS2_MAX_VALUE_BITS;
      vectors->num_bits_per_packedval = bits_per_byte * sizeof(GMBits2x64);
      /*---By design can only store this number of fields for this metric---*/
      GMInsistInterface(env,
               ((GMUInt64)(4 * num_field)) <
                       (((GMUInt64)1) << GM_TALLY1_MAX_VALUE_BITS)
                   ? "Number of fields requested is too large for this metric"
                   : 0);
    } break;
    /*--------------------*/
    default:
      GMInsist(false ? "Invalid data type." : 0);
  } /*---switch---*/

  vectors->num_val_per_packedval = vectors->num_bits_per_packedval /
                                   vectors->num_bits_per_val;

  /*---Calculate number of (packed) values to set aside storage for---*/

  vectors->num_packedval_field_local =
      gm_ceil_i8(vectors->num_field_local * (size_t)(vectors->num_bits_per_val),
                 vectors->num_bits_per_packedval);
  vectors->num_packedval_local =
      vectors->num_packedval_field_local * (size_t)num_vector_local;

  /*---Allocation for vector storage---*/

  GMInsist(vectors->num_bits_per_packedval % bits_per_byte == 0);

  vectors->data_size = vectors->num_packedval_local *
                       (vectors->num_bits_per_packedval / bits_per_byte);

  vectors->buf = GMMirroredBuf_null();

  if (vectors->has_buf) {
    GMMirroredBuf_create(&(vectors->buf),vectors->num_packedval_field_local,
                         num_vector_local, env);
    vectors->data = vectors->buf.h;
  } else {
    vectors->data = gm_malloc(vectors->data_size, env);
  }

  /*---Set pad entries to zero---*/

  GMVectors_initialize_pad(vectors, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      int num_field,
                      size_t num_field_active,
                      int num_vector_local,
                      GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(num_field >= 0);
  GMInsist(num_field_active >= 0);
  GMInsist(num_field_active <= (size_t)num_field);
  GMInsist(num_vector_local >= 0);

  *vectors = GMVectors_null();

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = false;

  GMVectors_create_imp_(vectors, data_type_id, num_field, num_field_active,
    num_vector_local, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               int num_field,
                               size_t num_field_active,
                               int num_vector_local,
                               GMEnv* env) {
  GMInsist(vectors && env);
  GMInsist(num_field >= 0);
  GMInsist(num_field_active >= 0);
  GMInsist(num_field_active <= (size_t)num_field);
  GMInsist(num_vector_local >= 0);

  *vectors = GMVectors_null();

  if (!GMEnv_is_proc_active(env)) {
    return;
  }

  vectors->has_buf = GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU;

  GMVectors_create_imp_(vectors, data_type_id, num_field, num_field_active,
    num_vector_local, env);
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

//=============================================================================

#ifdef __cplusplus
extern "C" {
#endif

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

//-----------------------------------------------------------------------------
