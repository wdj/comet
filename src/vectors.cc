//-----------------------------------------------------------------------------
/*!
 * \file   vectors.cc
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Class to manage the set of vectors taken as input to the methods.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//-----------------------------------------------------------------------------

#include "cstdlib"
#include "cstdint"
#include "cstdio"
#include "string.h"

#include "mpi.h"

#include "env.hh"
#include "decomp_mgr.hh"
#include "vectors.hh"

//=============================================================================

namespace comet {

//-----------------------------------------------------------------------------
/*---Null object---*/

GMVectors GMVectors_null() {
  GMVectors result;
  memset((void*)&result, 0, sizeof(GMVectors));
  return result;
}

//=============================================================================
/*---Set unused (pad) vector entries to zero---*/

void GMVectors_initialize_pad(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  /*---Ensure final pad words/bits of each vector are set to zero so that
       word-wise summations of bits aren't corrupted with bad trailing data---*/

  if (! env->is_proc_active()) {
    return;
  }

  const size_t nfal = vectors->dm->num_field_active_local;

  const size_t pfl_min = nfal / vectors->dm->num_field_per_packedfield;
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
      const GMBits2x64 zero = GMBits2x64_null();
      const uint64_t allbits = 0xffffffffffffffff;
      for (int vl = 0; vl < vectors->num_vector_local; ++vl) {
        for (size_t pfl = pfl_min; pfl < pfl_max; ++pfl) {
          const size_t fl = 64 * pfl;

          if (fl >= nfal) {

            GMVectors_bits2x64_set(vectors, pfl, vl, zero, env);

          } else if (fl + 32 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[0] &= allbits >> shift_dist;
            val.data[1] = 0;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } else if (fl + 64 >= nfal) {

            GMBits2x64 val = GMVectors_bits2x64_get(vectors, pfl, vl, env);
            const int shift_dist = 64 - 2*(nfal-fl-32);
            COMET_ASSERT(shift_dist >= 0 && shift_dist < 64);
            val.data[1] &= allbits >> shift_dist;
            GMVectors_bits2x64_set(vectors, pfl, vl, val, env);

          } // if
        } // for pfl
      }
    } break;
    default:
      COMET_INSIST(false && "Invalid vectors data_type_id.");
  } /*---switch---*/
}

//=============================================================================
/*---Vectors pseudo-constructor---*/

void GMVectors_create_imp_(GMVectors* vectors,
                           int data_type_id,
                           GMDecompMgr* dm,
                           CEnv* env) {
  COMET_INSIST(vectors && dm && env);

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
                       (vectors->num_bits_per_packedval / bits_per_byte);

  // Set up vector storage, mirrored buffer

  vectors->buf = new MirroredBuf(*env);

  if (vectors->has_buf) {
    vectors->buf->allocate(vectors->num_packedval_field_local,
                           vectors->num_vector_local);
    vectors->data = vectors->buf->h; // alias vector storage to buf
  } else {
    vectors->data = gm_malloc(vectors->data_size, env);
  }

  // Set pad entries to zero

  GMVectors_initialize_pad(vectors, env);

//  gm_tc_bufs_malloc(env, vectors->num_vector_local,
//                    vectors->num_packedval_field_local);
}

//-----------------------------------------------------------------------------

void GMVectors_create(GMVectors* vectors,
                      int data_type_id,
                      GMDecompMgr* dm,
                      CEnv* env) {
  COMET_INSIST(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! env->is_proc_active()) {
    return;
  }

  vectors->has_buf = false;

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//-----------------------------------------------------------------------------

void GMVectors_create_with_buf(GMVectors* vectors,
                               int data_type_id,
                               GMDecompMgr* dm,
                               CEnv* env) {
  COMET_INSIST(vectors && dm && env);

  *vectors = GMVectors_null();

  if (! env->is_proc_active()) {
    return;
  }

  vectors->has_buf = env->is_using_linalg();

  GMVectors_create_imp_(vectors, data_type_id, dm, env);
}

//=============================================================================
/*---Vectors pseudo-destructor---*/

void GMVectors_destroy(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);
  COMET_INSIST(vectors->data || ! env->is_proc_active());

  if (! env->is_proc_active()) {
    return;
  }

  if (!vectors->has_buf) {
    gm_free(vectors->data, vectors->data_size, env);
  }

  delete vectors->buf;

  *vectors = GMVectors_null();
}


//=============================================================================
// Copy vectors to mirrored buffer

void gm_vectors_to_buf(MirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       CEnv* env) {
  COMET_INSIST(vectors && vectors_buf && env);

  if (!env->is_using_linalg())
    return;

  /*---Copy vectors into GPU buffers if needed---*/

  switch (env->metric_type()) {
    case MetricType::CZEK: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_field_local; ++fl) {
          vectors_buf->elt<GMFloat>(fl, i) =
            GMVectors_float_get(vectors, fl, i, env);
        }
      }
    } break;
    case MetricType::CCC: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedval_field_local; ++fl) {
          vectors_buf->elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    case MetricType::DUO: {
      // don't use collapse because of overflow for large sizes
      //#pragma omp parallel for collapse(2) schedule(dynamic,1000)
      #pragma omp parallel for schedule(dynamic,1000)
      for (int i = 0; i < vectors->num_vector_local; ++i) {
        for (int fl = 0; fl < vectors->num_packedval_field_local; ++fl) {
          vectors_buf->elt<GMBits2x64>(fl, i) =
            GMVectors_bits2x64_get(vectors, fl, i, env);
        }
      }
    } break;
    default:
      COMET_INSIST_INTERFACE(env, false && "Unimplemented metric_type.");
  } /*---case---*/
}

//=============================================================================
// Print entries of vectors

void GMVectors_print(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active()) {
    return;
  }

  const int nval = vectors->dm->num_vector_active_local;
  const int nfal = vectors->dm->num_field_active_local;
  
  switch (env->data_type_vectors()) {
    /*--------------------*/
    case GM_DATA_TYPE_FLOAT: {
    /*--------------------*/
      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMFloat float_value = GMVectors_float_get(vectors, fl, vl, env);
            printf("vec_proc %i vec %i field_proc %i field %i value %e\n",
                   env->proc_num_vector(), vl,
                   env->proc_num_field(), fl, float_value);
        } /*---fl---*/
      }   /*---vl---*/
    } break;       
    /*--------------------*/
    case GM_DATA_TYPE_BITS2: {
    /*--------------------*/
      for (int vl = 0; vl < nval; ++vl) {
        for (int fl = 0; fl < nfal; ++fl) {
          const GMBits2 value = GMVectors_bits2_get(vectors, fl, vl, env);
            printf("vec_proc %i vec %i "
                   "field_proc %i field %i value %.1i%.1i\n",
                   env->proc_num_vector(), vl,
                   env->proc_num_field(), fl, value / 2, value % 2);
        } /*---fl---*/
      }   /*---vl---*/
    } break;
    /*--------------------*/
    default:
    /*--------------------*/
      COMET_INSIST(false && "Invalid data_type_vectors.");
  } /*---switch---*/
}

//=============================================================================
// hecksum of entries.

size_t GMVectors_cksum(GMVectors* vectors, CEnv* env) {
  COMET_INSIST(vectors && env);

  if (! env->is_proc_active()) {
    return 0;
  }

  return gm_array_cksum((unsigned char*)(vectors->data), vectors->data_size);
}

//=============================================================================

} // namespace comet

//-----------------------------------------------------------------------------
