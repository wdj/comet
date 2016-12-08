/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.cc
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "compute_utils_linalg.hh"
#include "compute_metrics_utils.hh"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Start/end MPI send/receive of vectors data---*/

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(proc_num >= 0 && proc_num < GMEnv_num_proc_vector_total(env));

  //const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /* clang-format off */
  const int mpi_type = GMEnv_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  mpi_code =
      MPI_Isend((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, GMEnv_mpi_comm_vector(env), &mpi_request);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  int mpi_tag,
                                  GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(proc_num >= 0 && proc_num < GMEnv_num_proc_vector_total(env));

  //const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /* clang-format off */
  const int mpi_type = GMEnv_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  mpi_code =
      MPI_Irecv((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, GMEnv_mpi_comm_vector(env), &mpi_request);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssertAlways(mpi_request != NULL);
  GMAssertAlways(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssertAlways(mpi_request != NULL);
  GMAssertAlways(env != NULL);

  MPI_Status mpi_status;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---MPI allreduce operations---*/

void gm_allreduce_metrics(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf_target,
                          GMMirroredPointer* metrics_buf_source,
                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(metrics_buf_target != NULL);
  GMAssertAlways(metrics_buf_source != NULL);
  GMAssertAlways(env != NULL);

  const int numvecl = metrics->num_vector_local;

  /* clang-format off */
  const int mpi_type = GMEnv_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(metrics_buf_source->h, metrics_buf_target->h,
                           numvecl * (size_t)numvecl, mpi_type, MPI_SUM,
                           GMEnv_mpi_comm_field(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Start/end transfer of vectors data to GPU---*/

void gm_set_vectors_start(GMVectors* vectors,
                          GMMirroredPointer* vectors_buf,
                          GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(vectors_buf != NULL);
  GMAssertAlways(env != NULL);

  gm_linalg_set_matrix_start(vectors_buf, vectors->num_packedval_field_local,
                             vectors->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_set_vectors_wait(GMEnv* env) {
  GMAssertAlways(env != NULL);

  gm_linalg_set_matrix_wait(env);
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(metrics_buf != NULL);
  GMAssertAlways(env != NULL);

  gm_linalg_get_matrix_start(metrics_buf, metrics->num_vector_local,
                             metrics->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_get_metrics_wait(GMMetrics* metrics,
                         GMMirroredPointer* metrics_buf,
                         GMEnv* env) {
  GMAssertAlways(env != NULL);

  gm_linalg_get_matrix_wait(env);

  if (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
      GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*---Adjust entries because of computation on pad values.
         EXPLANATION: the final word of each vector may have zero-pad bits
         to fill out the word.  The Magma call will tally these into the
         GMTally2x2 data[0] entry, because this is here the zero X zero
         semi-nibble pairs are tallied.  The code here fixes this by
         subtracting off this unwanted tally result.---*/
    /*---NOTE: this should work for both 2-way and 3-way---*/

    const int num_seminibbles_pad =
        64 - (1 + (metrics->num_field_local - 1) % 64);
    const GMFloat adjustment = 4 * num_seminibbles_pad;
    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      int i = 0;
      for (i = 0; i < metrics->num_vector_local; ++i) {
        ((GMTally2x2*)(metrics_buf->h))[i + metrics->num_vector_local * j]
            .data[0] -= adjustment;
      } /*---for j---*/
    }   /*---for i---*/
  }     /*---if---*/
}

/*===========================================================================*/
/*---CPU-GPU transfer buffer manipulation---*/

void gm_vectors_to_buf(GMMirroredPointer* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(vectors_buf != NULL);
  GMAssertAlways(env != NULL);

  if (GMEnv_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  int i = 0;
  int f = 0;

  switch (GMEnv_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      /*---Copy vectors into GPU buffers if needed---*/
      for (i = 0; i < vectors->num_vector_local; ++i) {
        const size_t nfl = vectors->num_field_local;
        for (f = 0; f < vectors->num_field_local; ++f) {
          ((GMFloat*)vectors_buf->h)[f + nfl * i] =
              GMVectors_float_get(vectors, f, i, env);
        }
      }
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      /*---Copy vectors into GPU buffers if needed---*/
      for (i = 0; i < vectors->num_vector_local; ++i) {
        const size_t npfl = vectors->num_packedval_field_local;
        for (f = 0; f < vectors->num_packedval_field_local; ++f) {
          ((GMBits2x64*)vectors_buf->h)[f + npfl * i] =
              GMVectors_bits2x64_get(vectors, f, i, env);
        }
      }
    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
