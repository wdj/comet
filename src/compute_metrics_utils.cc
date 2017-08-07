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
#include "mirrored_buf.hh"
#include "vector_sums.hh"
#include "vectors.hh"
#include "metrics.hh"
#include "linalg.hh"
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
  GMAssertAlways(vectors && env);
  GMAssertAlways(proc_num >= 0 && proc_num < GMEnv_num_proc_vector_total(env));

  MPI_Request mpi_request;

  const int mpi_type = gm_mpi_type(env);

  const int mpi_code =
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
  GMAssertAlways(vectors && env);
  GMAssertAlways(proc_num >= 0 && proc_num < GMEnv_num_proc_vector_total(env));

  MPI_Request mpi_request;

  const int mpi_type = gm_mpi_type(env);

  const int mpi_code =
      MPI_Irecv((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, GMEnv_mpi_comm_vector(env), &mpi_request);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

void gm_send_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssertAlways(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

void gm_recv_vectors_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssertAlways(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---MPI reduce operations---*/

void gm_reduce_metrics(GMMetrics* metrics,
                       GMMirroredBuf* metrics_buf_target,
                       GMMirroredBuf* metrics_buf_source,
                       GMEnv* env) {
  GMAssertAlways(metrics && env);
  GMAssertAlways(metrics_buf_target && metrics_buf_source);

  const int nvl = metrics->num_vector_local;

  const int mpi_type = gm_mpi_type(env);

  const int mpi_code = MPI_Allreduce(metrics_buf_source->h,
                                     metrics_buf_target->h,
                                     nvl * (size_t)nvl, mpi_type, MPI_SUM,
                                     GMEnv_mpi_comm_field(env));
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

MPI_Request gm_reduce_metrics_start(GMMetrics* metrics,
                                    GMMirroredBuf* metrics_buf_target,
                                    GMMirroredBuf* metrics_buf_source,
                                    GMEnv* env) {
  GMAssertAlways(metrics && env);
  GMAssertAlways(metrics_buf_target && metrics_buf_source);

  const int nvl = metrics->num_vector_local;

  const int mpi_type = gm_mpi_type(env);

  MPI_Request mpi_request;
  const int mpi_code = MPI_Iallreduce(metrics_buf_source->h,
                                      metrics_buf_target->h,
                                      nvl * (size_t)nvl, mpi_type, MPI_SUM,
                                      GMEnv_mpi_comm_field(env),
                                      &mpi_request);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

void gm_reduce_metrics_wait(MPI_Request* mpi_request, GMEnv* env) {
  GMAssertAlways(mpi_request && env);

  MPI_Status mpi_status;

  const int mpi_code = MPI_Wait(mpi_request, &mpi_status);
  GMAssertAlways(mpi_code == MPI_SUCCESS);
}

/*===========================================================================*/
/*---Start/end transfer of vectors data to GPU---*/

void gm_set_vectors_start(GMVectors* vectors, GMMirroredBuf* vectors_buf,
                          GMEnv* env) {
  GMAssertAlways(vectors && vectors_buf && env);

  gm_linalg_set_matrix_start(vectors_buf, env);
}

/*---------------------------------------------------------------------------*/

void gm_set_vectors_wait(GMEnv* env) {
  GMAssertAlways(env);

  gm_linalg_set_matrix_wait(env);
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_metrics_start(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                          GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf && env);

  gm_linalg_get_matrix_start(metrics_buf, env);
}

/*---------------------------------------------------------------------------*/

void gm_get_metrics_wait(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                         GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf && env);

  gm_linalg_get_matrix_wait(env);
}

/*---------------------------------------------------------------------------*/

int gm_num_seminibbles_pad(GMMetrics* metrics, GMEnv* env) {
  GMAssertAlways(metrics && env);

  const int num_bits_per_val = 2;
  const int num_bits_per_packedval = 128;
  const int num_val_per_packedval = 64;

  const int nfl = metrics->num_field_local;

  const int num_packedval_field_local
     = gm_ceil_i8(nfl * (size_t)num_bits_per_val, num_bits_per_packedval);

  const int num_field_calculated = num_packedval_field_local *
                                   num_val_per_packedval;

  //const int num_packedval_field_local = (nfl + 64 - 1) / 64;
  //const int num_field_calculated = num_packedval_field_local * 64;

  const bool final_proc = GMEnv_proc_num_field(env) ==
                          GMEnv_num_proc_field(env) - 1;

  const int num_field_active_local = final_proc
    ? nfl - (metrics->num_field - metrics->num_field_active) : nfl;

  GMAssertAlways(num_field_active_local >= 0);

  const int num_seminibbles_pad = num_field_calculated -
                                  num_field_active_local;

  return num_seminibbles_pad;
}

/*---------------------------------------------------------------------------*/

void gm_metrics_gpu_adjust(GMMetrics* metrics, GMMirroredBuf* metrics_buf,
                           GMEnv* env) {
  GMAssertAlways(metrics && metrics_buf && env);

  if (! (GMEnv_metric_type(env) == GM_METRIC_TYPE_CCC &&
      GMEnv_compute_method(env) == GM_COMPUTE_METHOD_GPU)) {
    return;
  }

  /*---Adjust entries because of computation on pad values.
       EXPLANATION: the final word of each vector may have zero-pad bits
       to fill out the word.  The Magma call will tally these into the
       GMTally2x2 data[0] entry, because this is here the zero X zero
       seminibble pairs are tallied.  The code here fixes this by
       subtracting off this unwanted tally result.---*/
  /*---NOTE: this should work for both 2-way and 3-way---*/

//  const int nfl = metrics->num_field_local;
//  const int num_packedval_field_local = (nfl + 64 - 1) / 64;
//  const int num_field_calculated = num_packedval_field_local * 64;
//  const int num_field_active_local =
//    GMEnv_proc_num_field(env) == GMEnv_num_proc_field(env)-1
//    ? nfl - (metrics->num_field - metrics->num_field_active) : nfl;
//  GMAssertAlways(num_field_active_local >= 0);
//  const int num_seminibbles_pad = num_field_calculated -
//                                  num_field_active_local;
  //const int num_seminibbles_pad =
  //    64 - (1 + (metrics->num_field_local - 1) % 64);

  const int num_seminibbles_pad = gm_num_seminibbles_pad(metrics, env);

  const GMFloat adjustment = 4 * num_seminibbles_pad;
#pragma omp parallel for collapse(2)
  for (int j = 0; j < metrics->num_vector_local; ++j) {
    for (int i = 0; i < metrics->num_vector_local; ++i) {
      GMMirroredBuf_elt<GMTally2x2>(metrics_buf, i, j)
      //((GMTally2x2*)(metrics_buf->h))[i + metrics->num_vector_local * j]
          .data[0] -= adjustment;


    } /*---for j---*/
  }   /*---for i---*/
}

/*===========================================================================*/
/*---CPU-GPU transfer buffer manipulation---*/

void gm_vectors_to_buf(GMMirroredBuf* vectors_buf,
                       GMVectors* vectors,
                       GMEnv* env) {
  GMAssertAlways(vectors && vectors_buf && env);

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
      GMInsist(env, false ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
