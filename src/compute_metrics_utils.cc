/*---------------------------------------------------------------------------*/
/*!
 * \file   compute_metrics_utils.c
 * \author Wayne Joubert, James Nance
 * \date   Fri Oct  9 14:06:44 EDT 2015
 * \brief  Functions for computing metrics, utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include "env.h"
#include "vector_sums.h"
#include "vectors.h"
#include "metrics.h"
#include "compute_utils_magma.h"
#include "compute_metrics_utils.h"

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Start/end MPI send/receive of vectors data---*/

MPI_Request gm_send_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(proc_num >= 0 && proc_num < Env_num_proc_vector_total(env));

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /* clang-format off */
  const int mpi_type = Env_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       Env_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  mpi_code =
      MPI_Isend((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, Env_mpi_comm_vector(env), &mpi_request);
  GMAssertAlways(mpi_code == MPI_SUCCESS);

  return mpi_request;
}

/*---------------------------------------------------------------------------*/

MPI_Request gm_recv_vectors_start(GMVectors* vectors,
                                  int proc_num,
                                  GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(proc_num >= 0 && proc_num < Env_num_proc_vector_total(env));

  const int mpi_tag = 0;
  MPI_Request mpi_request;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /* clang-format off */
  const int mpi_type = Env_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       Env_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  mpi_code =
      MPI_Irecv((void*)vectors->data, vectors->num_packedval_local, mpi_type,
                proc_num, mpi_tag, Env_mpi_comm_vector(env), &mpi_request);
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
  const int mpi_type = Env_metric_type(env) == GM_METRIC_TYPE_SORENSON ?
                         GM_MPI_FLOAT : /*---NOTE: not fully designed---*/
                       Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI ?
                         GM_MPI_FLOAT :
                       Env_metric_type(env) == GM_METRIC_TYPE_CCC ?
                         MPI_DOUBLE_COMPLEX :
                       0;
  /* clang-format on */

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Allreduce(metrics_buf_source->h, metrics_buf_target->h,
                           numvecl * (size_t)numvecl, mpi_type, MPI_SUM,
                           Env_mpi_comm_field(env));
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

  gm_set_matrix_start(vectors_buf, vectors->num_packedval_field_local,
                      vectors->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_set_vectors_wait(GMEnv* env) {
  GMAssertAlways(env != NULL);

  gm_set_matrix_wait(env);
}

/*===========================================================================*/
/*---Start/end transfer of metrics data from GPU---*/

void gm_get_metrics_start(GMMetrics* metrics,
                          GMMirroredPointer* metrics_buf,
                          GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(metrics_buf != NULL);
  GMAssertAlways(env != NULL);

  gm_get_matrix_start(metrics_buf, metrics->num_vector_local,
                      metrics->num_vector_local, env);
}

/*---------------------------------------------------------------------------*/

void gm_get_metrics_wait(GMMetrics* metrics,
                         GMMirroredPointer* metrics_buf,
                         GMEnv* env) {
  GMAssertAlways(env != NULL);

  gm_get_matrix_wait(env);

  if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
      Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
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

void gm_vectors_to_buf(GMVectors* vectors,
                       GMMirroredPointer* vectors_buf,
                       GMEnv* env) {
  GMAssertAlways(vectors != NULL);
  GMAssertAlways(vectors_buf != NULL);
  GMAssertAlways(env != NULL);

  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    return;
  }

  int i = 0;
  int f = 0;

  switch (Env_metric_type(env)) {
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
/*---Start calculation of numerators, 2-way Czekanowski---*/

void gm_compute_czekanowski_numerators_2way_start(
    GMVectors* vectors_left,
    GMVectors* vectors_right,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_left_buf,
    GMMirroredPointer* vectors_right_buf,
    GMMirroredPointer* metrics_buf,
    int j_block,
    _Bool do_compute_triang_only,
    GMEnv* env) {
  GMAssertAlways(vectors_left != NULL);
  GMAssertAlways(vectors_right != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  int i = 0;
  int j = 0;
  int f = 0;

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && Env_all2all(env)) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Perform pseudo matrix-matrix product---*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Perform pseudo matrix-matrix product---*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      for (i = 0; i < i_max; ++i) {
        GMFloat metric = 0;
        for (f = 0; f < vectors_left->num_field_local; ++f) {
          const GMFloat value1 = GMVectors_float_get(vectors_left, f, i, env);
          const GMFloat value2 = GMVectors_float_get(vectors_right, f, j, env);
          metric += value1 < value2 ? value1 : value2;
        } /*---for k---*/
        GMMetrics_float_set_2(metrics, i, j, metric, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_magma_set_matrix_zero_start(metrics_buf, metrics->num_vector_local,
                                   metrics->num_vector_local, env);

    /*---Perform pseudo matrix-matrix product---*/

    /* .63 / 1.56 */
    gm_magma_gemm_start(vectors_left->num_vector_local,
                        vectors_left->num_vector_local,
                        vectors_left->num_field_local,
                        vectors_left_buf->d, vectors_left->num_field_local,
                        vectors_right_buf->d, vectors_left->num_field_local,
                        metrics_buf->d, vectors_left->num_vector_local, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 2-way CCC---*/

void gm_compute_ccc_numerators_2way_start(GMVectors* vectors_left,
                                          GMVectors* vectors_right,
                                          GMMetrics* metrics,
                                          GMMirroredPointer* vectors_left_buf,
                                          GMMirroredPointer* vectors_right_buf,
                                          GMMirroredPointer* metrics_buf,
                                          int j_block,
                                          _Bool do_compute_triang_only,
                                          GMEnv* env) {
  GMAssertAlways(vectors_left != NULL);
  GMAssertAlways(vectors_right != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  /*----------------------------------------*/
  if (Env_compute_method(env) == GM_COMPUTE_METHOD_REF) {
    /*----------------------------------------*/

    /*---Perform pseudo matrix-matrix product---*/

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        int f = 0;
        for (f = 0; f < vectors_left->num_field_local; ++f) {
          const GMBits2 value_i = GMVectors_bits2_get(vectors_left, f, i, env);
          const GMBits2 value_j = GMVectors_bits2_get(vectors_right, f, j, env);

          /* clang-format off */
          const int r00 = ( ( !(value_i & 1) ) && ( !(value_j & 1) ) ) +
                          ( ( !(value_i & 1) ) && ( !(value_j & 2) ) ) +
                          ( ( !(value_i & 2) ) && ( !(value_j & 1) ) ) +
                          ( ( !(value_i & 2) ) && ( !(value_j & 2) ) );
          const int r01 = ( ( !(value_i & 1) ) && (  (value_j & 1) ) ) +
                          ( ( !(value_i & 1) ) && (  (value_j & 2) ) ) +
                          ( ( !(value_i & 2) ) && (  (value_j & 1) ) ) +
                          ( ( !(value_i & 2) ) && (  (value_j & 2) ) );
          const int r10 = ( (  (value_i & 1) ) && ( !(value_j & 1) ) ) +
                          ( (  (value_i & 1) ) && ( !(value_j & 2) ) ) +
                          ( (  (value_i & 2) ) && ( !(value_j & 1) ) ) +
                          ( (  (value_i & 2) ) && ( !(value_j & 2) ) );
          const int r11 = ( (  (value_i & 1) ) && (  (value_j & 1) ) ) +
                          ( (  (value_i & 1) ) && (  (value_j & 2) ) ) +
                          ( (  (value_i & 2) ) && (  (value_j & 1) ) ) +
                          ( (  (value_i & 2) ) && (  (value_j & 2) ) );
          /* clang-format on */

          /*---NOTE: "since the sum of all 4 of these relative
               co-occurences is 1, we really only need to compute 3 of them.
               Then the last one is just 1 minus the rest." */

#if DOUG_WAY
//TODO: work on this as a possibly faster way.
          const int vi1 = (value_i & 3) != 0;
          const int vi0 = ((~value_i) & 3) != 0;
          const int vj1 = (value_j & 3) != 0;
          const int vj0 = ((~value_j) & 3) != 0;

          const int a11 = vi1 & vj1;

          const int r11 = a11 +
#endif

          /*---Accumulate---*/

          sum.data[0] += GMTally1_encode(r00, r01);
          sum.data[1] += GMTally1_encode(r10, r11);
        } /*---for f---*/
        if (Env_all2all(env)) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) == GM_COMPUTE_METHOD_CPU) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Perform pseudo matrix-matrix product---*/

    /*---Precompute masks for final packedval_field -
         can be 1 to 64 inclusive---*/

    /* clang-format off */

    const int num_seminibbles = 1 + (vectors_left->num_field_local-1) % 64;

    const GMUInt64 nobits = 0;
    const GMUInt64 allbits = 0xffffffffffffffff;

    const GMUInt64 lastmask0 = num_seminibbles >= 32 ?
                               allbits :
                               allbits >> (64 - 2*num_seminibbles);

    const GMUInt64 lastmask1 = num_seminibbles <= 32 ?
                               nobits :
                               num_seminibbles == 64 ?
                               allbits :
                               allbits >> (128 - 2*num_seminibbles);

    int j = 0;
    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max =
          do_compute_triang_only ? j : metrics->num_vector_local;
      int i = 0;
      for (i = 0; i < i_max; ++i) {
        GMTally2x2 sum = GMTally2x2_null();
        int f = 0;
        const int f_last = vectors_left->num_packedval_field_local - 1;
        for (f = 0; f <= f_last ; ++f) {
          /*---Get masks for active seminibbles in each word---*/

          const GMUInt64 activebits0 = f < f_last ? allbits : lastmask0;
          const GMUInt64 activebits1 = f < f_last ? allbits : lastmask1;

          /*---Extract input values to process---*/
          const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_left, f, i, env);
          const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_right, f, j,
                                                       env);
          const GMUInt64 vi0 = vi.data[0];
          const GMUInt64 vi1 = vi.data[1];
          const GMUInt64 vj0 = vj.data[0];
          const GMUInt64 vj1 = vj.data[1];

          /*---Get even, odd bits for each semi-nibble, masked to active---*/

          const GMUInt64 oddbits = 0x5555555555555555;

          const GMUInt64 vi0_0 =  vi0       & oddbits & activebits0;
          const GMUInt64 vi0_1 = (vi0 >> 1) & oddbits & activebits0;
          const GMUInt64 vi1_0 =  vi1       & oddbits & activebits1;
          const GMUInt64 vi1_1 = (vi1 >> 1) & oddbits & activebits1;
          const GMUInt64 vj0_0 =  vj0       & oddbits & activebits0;
          const GMUInt64 vj0_1 = (vj0 >> 1) & oddbits & activebits0;
          const GMUInt64 vj1_0 =  vj1       & oddbits & activebits1;
          const GMUInt64 vj1_1 = (vj1 >> 1) & oddbits & activebits1;

          /*---Get complements of the same bits, set other bits zero---*/

          const GMUInt64 nvi0_0 = ~ vi0       & oddbits & activebits0;
          const GMUInt64 nvi0_1 = ~(vi0 >> 1) & oddbits & activebits0;
          const GMUInt64 nvi1_0 = ~ vi1       & oddbits & activebits1;
          const GMUInt64 nvi1_1 = ~(vi1 >> 1) & oddbits & activebits1;
          const GMUInt64 nvj0_0 = ~ vj0       & oddbits & activebits0;
          const GMUInt64 nvj0_1 = ~(vj0 >> 1) & oddbits & activebits0;
          const GMUInt64 nvj1_0 = ~ vj1       & oddbits & activebits1;
          const GMUInt64 nvj1_1 = ~(vj1 >> 1) & oddbits & activebits1;

          const int r00 = gm_popcount64((nvi0_0 & nvj0_0) |
                                      ( (nvi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 & nvj0_0) |
                                      ( (nvi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 & nvj1_0) |
                                      ( (nvi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 & nvj1_0) |
                                      ( (nvi1_1 & nvj1_1) << 1 ));
          const int r01 = gm_popcount64((nvi0_0 &  vj0_0) |
                                      ( (nvi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi0_1 &  vj0_0) |
                                      ( (nvi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64((nvi1_0 &  vj1_0) |
                                      ( (nvi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64((nvi1_1 &  vj1_0) |
                                      ( (nvi1_1 &  vj1_1) << 1 ));
          const int r10 = gm_popcount64(( vi0_0 & nvj0_0) |
                                      ( ( vi0_0 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 & nvj0_0) |
                                      ( ( vi0_1 & nvj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 & nvj1_0) |
                                      ( ( vi1_0 & nvj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 & nvj1_0) |
                                      ( ( vi1_1 & nvj1_1) << 1 ));
          const int r11 = gm_popcount64(( vi0_0 &  vj0_0) |
                                      ( ( vi0_0 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi0_1 &  vj0_0) |
                                      ( ( vi0_1 &  vj0_1) << 1 )) +
                          gm_popcount64(( vi1_0 &  vj1_0) |
                                      ( ( vi1_0 &  vj1_1) << 1 )) +
                          gm_popcount64(( vi1_1 &  vj1_0) |
                                      ( ( vi1_1 &  vj1_1) << 1 ));

          /*---Accumulate---*/

          sum.data[0] += GMTally1_encode(r00, r01);
          sum.data[1] += GMTally1_encode(r10, r11);
        } /*---for f---*/
        if (Env_all2all(env)) {
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, sum, env);
        } else {
          GMMetrics_tally2x2_set_2(metrics, i, j, sum, env);
        }
      } /*---for j---*/
    }   /*---for i---*/

    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_magma_set_matrix_zero_start(metrics_buf, metrics->num_vector_local,
                                   metrics->num_vector_local, env);

    /*---Perform pseudo matrix-matrix product---*/

    gm_magma_gemm_start(
        vectors_left->num_vector_local, vectors_left->num_vector_local,
        vectors_left->num_packedval_field_local,
        vectors_left_buf->d, vectors_left->num_packedval_field_local,
        vectors_right_buf->d, vectors_left->num_packedval_field_local,
        metrics_buf->d, vectors_left->num_vector_local, env);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 2-way generic---*/

void gm_compute_numerators_2way_start(GMVectors* vectors_left,
                                      GMVectors* vectors_right,
                                      GMMetrics* numerators,
                                      GMMirroredPointer* vectors_left_buf,
                                      GMMirroredPointer* vectors_right_buf,
                                      GMMirroredPointer* numerators_buf,
                                      int j_block,
                                      _Bool do_compute_triang_only,
                                      GMEnv* env) {
  GMAssertAlways(vectors_left != NULL);
  GMAssertAlways(vectors_right != NULL);
  GMAssertAlways(numerators != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/

      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);

    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_czekanowski_numerators_2way_start(
          vectors_left, vectors_right, numerators, vectors_left_buf,
          vectors_right_buf, numerators_buf, j_block, do_compute_triang_only,
          env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      gm_compute_ccc_numerators_2way_start(vectors_left, vectors_right,
                                           numerators, vectors_left_buf,
                                           vectors_right_buf, numerators_buf,
                                           j_block, do_compute_triang_only, env);
    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 2-way Czek---*/

void gm_compute_czekanowski_2way_combine(
    GMMetrics* metrics,
    GMMirroredPointer* metrics_buf,
    GMFloat* __restrict__ vector_sums_left,
    GMFloat* __restrict__ vector_sums_right,
    int j_block,
    _Bool do_compute_triang_only,
    GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_left != NULL);
  GMAssertAlways(vector_sums_right != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  /*---For CPU case, copy numerator out of metrics struct which is temporarily
       holding numerators.
       For GPU case, directly access the metrics_buf holding the numerators.
  ---*/

  const _Bool are_vector_sums_aliased = vector_sums_left == vector_sums_right;

  int i = 0;
  int j = 0;

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && Env_all2all(env)) {
    /*----------------------------------------*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            GMMetrics_float_get_all2all_2(metrics, i, j, j_block, env);
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator =
            are_vector_sums_aliased
                ? vector_sums_left[i] + vector_sums_left[j]
                : vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/
        /*---TODO: here and elsewhere check for unlikely case denom is/nearly
         * zero---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator = GMMetrics_float_get_2(metrics, i, j, env);
        /*---Don't use two different pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else if (Env_all2all(env)) {
    /*----------------------------------------*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = do_compute_triang_only ? j : metrics->num_vector_local;
      const size_t nvl = metrics->num_vector_local;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            ((GMFloat*)metrics_buf->h)[i + nvl * j];
        /*---Don't use two pointers pointing to the same thing---*/
        const GMFloat denominator =
            are_vector_sums_aliased
                ? vector_sums_left[i] + vector_sums_left[j]
                : vector_sums_left[i] + vector_sums_right[j];
        GMMetrics_float_set_all2all_2(metrics, i, j, j_block,
                                      2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } else {
    /*----------------------------------------*/

    for (j = 0; j < metrics->num_vector_local; ++j) {
      const int i_max = j;
      const size_t nvl = metrics->num_vector_local;
      for (i = 0; i < i_max; ++i) {
        const GMFloat numerator =
            ((GMFloat*)metrics_buf->h)[i + nvl * j];
        /*---Don't use two different pointers pointing to the same thing---*/
        const GMFloat denominator = vector_sums_left[i] + vector_sums_left[j];
        GMMetrics_float_set_2(metrics, i, j, 2 * numerator / denominator, env);
      } /*---for i---*/
    }   /*---for j---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 2-way CCC---*/

void gm_compute_ccc_2way_combine(GMMetrics* metrics,
                                 GMMirroredPointer* metrics_buf,
                                 GMFloat* __restrict__ vector_sums_left,
                                 GMFloat* __restrict__ vector_sums_right,
                                 int j_block,
                                 _Bool do_compute_triang_only,
                                 GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_left != NULL);
  GMAssertAlways(vector_sums_right != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  const int numvecl = metrics->num_vector_local;

  int i = 0;
  int j = 0;

  /*---Copy from metrics_buffer for GPU case---*/

  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*--------------------*/
    if (Env_all2all(env)) {
      /*--------------------*/
      for (j = 0; j < numvecl; ++j) {
        const int i_max = do_compute_triang_only ? j : numvecl;
        for (i = 0; i < i_max; ++i) {
          const GMTally2x2 value = ((
              GMTally2x2*)(metrics_buf->h))[i + metrics->num_vector_local * j];
          GMMetrics_tally2x2_set_all2all_2(metrics, i, j, j_block, value, env);
#ifdef GM_ASSERTIONS_ON
          const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
          const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
          const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
          const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
          GMAssert((GMUInt64)r00 + (GMUInt64)r01 + (GMUInt64)r10 +
                       (GMUInt64)r11 ==
                   (GMUInt64)(4 * metrics->num_field));
          const GMTally1 si_1 = (GMTally1)(vector_sums_left[i]);
          const GMTally1 sj_1 = (GMTally1)(vector_sums_right[j]);
          GMAssert((GMUInt64)r10 + (GMUInt64)r11 == (GMUInt64)(2 * si_1));
          GMAssert((GMUInt64)r01 + (GMUInt64)r11 == (GMUInt64)(2 * sj_1));
#endif
        } /*---for i---*/
      }   /*---for j---*/
      /*--------------------*/
    } else /*---(!Env_all2all(env))---*/ {
      /*--------------------*/
      for (j = 0; j < numvecl; ++j) {
        const int i_max = do_compute_triang_only ? j : numvecl;
        for (i = 0; i < i_max; ++i) {
          const GMTally2x2 value = ((
              GMTally2x2*)(metrics_buf->h))[i + metrics->num_vector_local * j];
          GMMetrics_tally2x2_set_2(metrics, i, j, value, env);
#ifdef GM_ASSERTIONS_ON
          const GMTally1 r00 = GMTally2x2_get(value, 0, 0);
          const GMTally1 r01 = GMTally2x2_get(value, 0, 1);
          const GMTally1 r10 = GMTally2x2_get(value, 1, 0);
          const GMTally1 r11 = GMTally2x2_get(value, 1, 1);
          GMAssert((GMUInt64)r00 + (GMUInt64)r01 + (GMUInt64)r10 +
                       (GMUInt64)r11 ==
                   (GMUInt64)(4 * metrics->num_field));
          const GMTally1 si_1 = (GMTally1)(vector_sums_left[i]);
          const GMTally1 sj_1 = (GMTally1)(vector_sums_right[j]);
          GMAssert((GMUInt64)r10 + (GMUInt64)r11 == (GMUInt64)(2 * si_1));
          GMAssert((GMUInt64)r01 + (GMUInt64)r11 == (GMUInt64)(2 * sj_1));
#endif
        } /*---for i---*/
      }   /*---for j---*/
      /*--------------------*/
    } /*---if---*/
    /*--------------------*/
  }

  /*---Compute multipliers---*/

  /*--------------------*/
  if (Env_all2all(env)) {
    /*--------------------*/
    for (j = 0; j < numvecl; ++j) {
      const GMTally1 sj_1 = (GMTally1)(vector_sums_right[j]);
      const int i_max = do_compute_triang_only ? j : numvecl;
      for (i = 0; i < i_max; ++i) {
        const GMTally1 si_1 = (GMTally1)(vector_sums_left[i]);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si_1, sj_1);
        GMMetrics_float2_M_set_all2all_2(metrics, i, j, j_block, si1_sj1, env);
      } /*---for i---*/
    }   /*---for j---*/
    /*--------------------*/
  } else /*---(!Env_all2all(env))---*/ {
    /*--------------------*/
    for (j = 0; j < numvecl; ++j) {
      const GMTally1 sj_1 = (GMTally1)(vector_sums_right[j]);
      const int i_max = do_compute_triang_only ? j : numvecl;
      for (i = 0; i < i_max; ++i) {
        const GMTally1 si_1 = (GMTally1)(vector_sums_left[i]);
        const GMFloat2 si1_sj1 = GMFloat2_encode(si_1, sj_1);
        GMMetrics_float2_M_set_2(metrics, i, j, si1_sj1, env);
      } /*---for i---*/
    }   /*---for j---*/
    /*--------------------*/
  } /*---if---*/
  /*--------------------*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 2-way generic---*/

void gm_compute_2way_combine(GMMetrics* metrics,
                             GMMirroredPointer* metrics_buf,
                             GMVectorSums* vector_sums_left,
                             GMVectorSums* vector_sums_right,
                             int j_block,
                             _Bool do_compute_triang_only,
                             GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_left != NULL);
  GMAssertAlways(vector_sums_right != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_2);

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_czekanowski_2way_combine(metrics, metrics_buf,
                                          (GMFloat*)vector_sums_left->data,
                                          (GMFloat*)vector_sums_right->data,
                                          j_block, do_compute_triang_only, env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      gm_compute_ccc_2way_combine(metrics, metrics_buf,
                                  (GMFloat*)vector_sums_left->data,
                                  (GMFloat*)vector_sums_right->data, j_block,
                                  do_compute_triang_only, env);
    } break;
    /*----------------------------------------*/
    default:
      /*----------------------------------------*/
      /*---Should never get here---*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  } /*---case---*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way Czekanowski non-gpu---*/

void gm_compute_czekanowski_numerators_3way_nongpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_block,
    int k_block,
    GMEnv* env) {
  GMAssertAlways(vectors_i != NULL);
  GMAssertAlways(vectors_j != NULL);
  GMAssertAlways(vectors_k != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors_i_buf != NULL);
  GMAssertAlways(vectors_j_buf != NULL);
  GMAssertAlways(vectors_k_buf != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == j_block &&
                   Env_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == k_block &&
                   Env_proc_num_vector_i(env) != j_block));
  GMAssertAlways(Env_compute_method(env) != GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  int i = 0;
  int j = 0;
  int k = 0;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && !Env_all2all(env)) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---No off-proc all2all: compute tetrahedron of values---*/

    for (k = 0; k < numvecl; ++k) {
      for (j = 0; j < k; ++j) {
        for (i = 0; i < j; ++i) {
          GMFloat sum = 0;
          int f = 0;
          for (f = 0; f < numfieldl; ++f) {
            const GMFloat val1 = GMVectors_float_get(vectors_i, f, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_j, f, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_k, f, k, env);
            GMFloat min12 = val1 < val2 ? val1 : val2;
            sum += min12;
            sum += val1 < val3 ? val1 : val3;
            sum += val2 < val3 ? val2 : val3;
            sum -= min12 < val3 ? min12 : val3;
          } /*---for f---*/
          GMMetrics_float_set_3(metrics, i, j, k, sum, env);
        }
      }
    }

    /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Compute tetrahedron, triang prism or block section---*/

    for (k = si->k_lb; k < si->k_ub; ++k) {
      const int j_max = si->is_part3 ? si->j_ub : k;
      for (j = si->j_lb; j < j_max; ++j) {
        const int i_max = si->is_part1 ? j : si->i_ub;
        for (i = si->i_lb; i < i_max; ++i) {
          GMFloat numerator = 0;
          int f = 0;
          for (f = 0; f < numfieldl; ++f) {
            const GMFloat val1 = GMVectors_float_get(vectors_i, f, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_j, f, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_k, f, k, env);
            const GMFloat min_ij = val1 < val2 ? val1 : val2;
            const GMFloat min_ik = val1 < val3 ? val1 : val3;
            const GMFloat min_jk = val2 < val3 ? val2 : val3;
            const GMFloat min_ijk = min_ij < val3 ? min_ij : val3;
            numerator += min_ij + min_ik + min_jk - min_ijk;
          } /*---for f---*/
          GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                        numerator, env);
        }
      }
    }

    /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    GMAssertAlways(GM_BOOL_FALSE
                 ? "logic error - code branch should never be executed"
                 : 0);

  } /*---if GPU---*/

  GMSectionInfo_destroy(si, env);;
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way CCC non-gpu---*/

void gm_compute_ccc_numerators_3way_nongpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_block,
    int k_block,
    GMEnv* env) {
  GMAssertAlways(vectors_i != NULL);
  GMAssertAlways(vectors_j != NULL);
  GMAssertAlways(vectors_k != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors_i_buf != NULL);
  GMAssertAlways(vectors_j_buf != NULL);
  GMAssertAlways(vectors_k_buf != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == j_block &&
                   Env_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == k_block &&
                   Env_proc_num_vector_i(env) != j_block));
  GMAssertAlways(Env_compute_method(env) != GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  if (Env_compute_method(env) == GM_COMPUTE_METHOD_REF) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---No off-proc all2all: compute tetrahedron of values---*/

    int k = 0;
    const int k_min = (!Env_all2all(env)) ? 0 : si->k_lb;
    const int k_max = (!Env_all2all(env)) ? numvecl : si->k_ub;
    for (k = k_min; k < k_max; ++k) {
      const int j_min = (!Env_all2all(env)) ? 0 : si->j_lb;
      const int j_max = (!Env_all2all(env)) ? k : si->is_part3 ? si->j_ub : k;
      int j = 0;
      for (j = j_min; j < j_max; ++j) {
        const int i_min = (!Env_all2all(env)) ? 0 : si->i_lb;
        const int i_max = (!Env_all2all(env)) ? j : si->is_part1 ? j : si->i_ub;
        int i = 0;
        for (i = i_min; i < i_max; ++i) {
          GMTally4x2 sum = GMTally4x2_null();
          int f = 0;
          for (f = 0; f < numfieldl; ++f) {
            const GMBits2 value_i = GMVectors_bits2_get(vectors_i, f, i, env);
            const GMBits2 value_j = GMVectors_bits2_get(vectors_j, f, j, env);
            const GMBits2 value_k = GMVectors_bits2_get(vectors_k, f, k, env);
            /* clang-format off */
            const int r000 =
              (( !(value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 1) )) +
              (( !(value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 2) )) +
              (( !(value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 1) )) +
              (( !(value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 2) )) +
              (( !(value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 1) )) +
              (( !(value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 2) )) +
              (( !(value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 1) )) +
              (( !(value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ));
            const int r001 =
              (( !(value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 1) )) +
              (( !(value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 2) )) +
              (( !(value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 1) )) +
              (( !(value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 2) )) +
              (( !(value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 1) )) +
              (( !(value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 2) )) +
              (( !(value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 1) )) +
              (( !(value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 2) ));
            const int r010 =
              (( !(value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 1) )) +
              (( !(value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 2) )) +
              (( !(value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 1) )) +
              (( !(value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 2) )) +
              (( !(value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 1) )) +
              (( !(value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 2) )) +
              (( !(value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 1) )) +
              (( !(value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 2) ));
            const int r011 =
              (( !(value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 1) )) +
              (( !(value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 2) )) +
              (( !(value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 1) )) +
              (( !(value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 2) )) +
              (( !(value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 1) )) +
              (( !(value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 2) )) +
              (( !(value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 1) )) +
              (( !(value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 2) ));
            const int r100 =
              ((  (value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 1) )) +
              ((  (value_i & 1) ) && ( !(value_j & 1) ) && ( !(value_k & 2) )) +
              ((  (value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 1) )) +
              ((  (value_i & 1) ) && ( !(value_j & 2) ) && ( !(value_k & 2) )) +
              ((  (value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 1) )) +
              ((  (value_i & 2) ) && ( !(value_j & 1) ) && ( !(value_k & 2) )) +
              ((  (value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 1) )) +
              ((  (value_i & 2) ) && ( !(value_j & 2) ) && ( !(value_k & 2) ));
            const int r101 =
              ((  (value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 1) )) +
              ((  (value_i & 1) ) && ( !(value_j & 1) ) && (  (value_k & 2) )) +
              ((  (value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 1) )) +
              ((  (value_i & 1) ) && ( !(value_j & 2) ) && (  (value_k & 2) )) +
              ((  (value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 1) )) +
              ((  (value_i & 2) ) && ( !(value_j & 1) ) && (  (value_k & 2) )) +
              ((  (value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 1) )) +
              ((  (value_i & 2) ) && ( !(value_j & 2) ) && (  (value_k & 2) ));
            const int r110 =
              ((  (value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 1) )) +
              ((  (value_i & 1) ) && (  (value_j & 1) ) && ( !(value_k & 2) )) +
              ((  (value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 1) )) +
              ((  (value_i & 1) ) && (  (value_j & 2) ) && ( !(value_k & 2) )) +
              ((  (value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 1) )) +
              ((  (value_i & 2) ) && (  (value_j & 1) ) && ( !(value_k & 2) )) +
              ((  (value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 1) )) +
              ((  (value_i & 2) ) && (  (value_j & 2) ) && ( !(value_k & 2) ));
            const int r111 =
              ((  (value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 1) )) +
              ((  (value_i & 1) ) && (  (value_j & 1) ) && (  (value_k & 2) )) +
              ((  (value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 1) )) +
              ((  (value_i & 1) ) && (  (value_j & 2) ) && (  (value_k & 2) )) +
              ((  (value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 1) )) +
              ((  (value_i & 2) ) && (  (value_j & 1) ) && (  (value_k & 2) )) +
              ((  (value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 1) )) +
              ((  (value_i & 2) ) && (  (value_j & 2) ) && (  (value_k & 2) ));
            /* clang-format on */

//              printf("%i %i %i %i %i %i %i %i\n",
//                (int)r000, (int)r001, (int)r010, (int)r011,
//                (int)r100, (int)r101, (int)r110, (int)r111);
            sum.data[0] += GMTally1_encode(r000, r001);
            sum.data[1] += GMTally1_encode(r010, r011);
            sum.data[2] += GMTally1_encode(r100, r101);
            sum.data[3] += GMTally1_encode(r110, r111);
          } /*---for f---*/
          if (Env_all2all(env)) {
            GMMetrics_tally4x2_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                             sum, env);
          } else {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
          }
        } /*---for i---*/
      }   /*---for j---*/
    }     /*---for k---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    /*---Precompute masks for final packedval_field -
         can be 1 to 64 inclusive---*/

    /* clang-format off */

    const int num_seminibbles = 1 + (vectors_i->num_field_local-1) % 64;

    const GMUInt64 nobits = 0;
    const GMUInt64 allbits = 0xffffffffffffffff;

    const GMUInt64 lastmask0 = num_seminibbles >= 32 ?
                               allbits :
                               allbits >> (64 - 2*num_seminibbles);

    const GMUInt64 lastmask1 = num_seminibbles <= 32 ?
                               nobits :
                               num_seminibbles == 64 ?
                               allbits :
                               allbits >> (128 - 2*num_seminibbles);
    int k = 0;
    const int k_min = (!Env_all2all(env)) ? 0 : si->k_lb;
    const int k_max = (!Env_all2all(env)) ? numvecl : si->k_ub;
    for (k = k_min; k < k_max; ++k) {
      const int j_min = (!Env_all2all(env)) ? 0 : si->j_lb;
      const int j_max = (!Env_all2all(env)) ? k : si->is_part3 ? si->j_ub : k;
      int j = 0;
      for (j = j_min; j < j_max; ++j) {
        const int i_min = (!Env_all2all(env)) ? 0 : si->i_lb;
        const int i_max = (!Env_all2all(env)) ? j : si->is_part1 ? j : si->i_ub;
        int i = 0;
        for (i = i_min; i < i_max; ++i) {
          GMTally4x2 sum = GMTally4x2_null();
          int f = 0;
          const int f_last = vectors_i->num_packedval_field_local - 1;
          for (f = 0; f <= f_last ; ++f) {
            /*---Get masks for active seminibbles in each word---*/

            const GMUInt64 activebits0 = f < f_last ? allbits : lastmask0;
            const GMUInt64 activebits1 = f < f_last ? allbits : lastmask1;

            /*---Extract input values to process---*/

            const GMBits2x64 vi = GMVectors_bits2x64_get(vectors_i, f, i, env);
            const GMBits2x64 vj = GMVectors_bits2x64_get(vectors_j, f, j, env);
            const GMBits2x64 vk = GMVectors_bits2x64_get(vectors_k, f, k, env);

            const GMUInt64 vi0 = vi.data[0];
            const GMUInt64 vi1 = vi.data[1];
            const GMUInt64 vj0 = vj.data[0];
            const GMUInt64 vj1 = vj.data[1];
            const GMUInt64 vk0 = vk.data[0];
            const GMUInt64 vk1 = vk.data[1];

            /*---Get even, odd bits for each semi-nibble, masked to active---*/

            const GMUInt64 oddbits = 0x5555555555555555;

            const GMUInt64 vi0_0 =  vi0       & oddbits & activebits0;
            const GMUInt64 vi0_1 = (vi0 >> 1) & oddbits & activebits0;
            const GMUInt64 vi1_0 =  vi1       & oddbits & activebits1;
            const GMUInt64 vi1_1 = (vi1 >> 1) & oddbits & activebits1;
            const GMUInt64 vj0_0 =  vj0       & oddbits & activebits0;
            const GMUInt64 vj0_1 = (vj0 >> 1) & oddbits & activebits0;
            const GMUInt64 vj1_0 =  vj1       & oddbits & activebits1;
            const GMUInt64 vj1_1 = (vj1 >> 1) & oddbits & activebits1;
            const GMUInt64 vk0_0 =  vk0       & oddbits & activebits0;
            const GMUInt64 vk0_1 = (vk0 >> 1) & oddbits & activebits0;
            const GMUInt64 vk1_0 =  vk1       & oddbits & activebits1;
            const GMUInt64 vk1_1 = (vk1 >> 1) & oddbits & activebits1;

            /*---Get complements of the same bits, set other bits zero---*/

            const GMUInt64 nvi0_0 = ~ vi0       & oddbits & activebits0;
            const GMUInt64 nvi0_1 = ~(vi0 >> 1) & oddbits & activebits0;
            const GMUInt64 nvi1_0 = ~ vi1       & oddbits & activebits1;
            const GMUInt64 nvi1_1 = ~(vi1 >> 1) & oddbits & activebits1;
            const GMUInt64 nvj0_0 = ~ vj0       & oddbits & activebits0;
            const GMUInt64 nvj0_1 = ~(vj0 >> 1) & oddbits & activebits0;
            const GMUInt64 nvj1_0 = ~ vj1       & oddbits & activebits1;
            const GMUInt64 nvj1_1 = ~(vj1 >> 1) & oddbits & activebits1;
            const GMUInt64 nvk0_0 = ~ vk0       & oddbits & activebits0;
            const GMUInt64 nvk0_1 = ~(vk0 >> 1) & oddbits & activebits0;
            const GMUInt64 nvk1_0 = ~ vk1       & oddbits & activebits1;
            const GMUInt64 nvk1_1 = ~(vk1 >> 1) & oddbits & activebits1;

            const int r000 = gm_popcount64((nvi0_0 & nvj0_0 & nvk0_0) |
                                         ( (nvi0_0 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 & nvj0_1 & nvk0_0) |
                                         ( (nvi0_0 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_0 & nvk0_0) |
                                         ( (nvi0_1 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_1 & nvk0_0) |
                                         ( (nvi0_1 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_0 & nvk1_0) |
                                         ( (nvi1_0 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_1 & nvk1_0) |
                                         ( (nvi1_0 & nvj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_0 & nvk1_0) |
                                         ( (nvi1_1 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_1 & nvk1_0) |
                                         ( (nvi1_1 & nvj1_1 & nvk1_1) << 1 ));
            const int r001 = gm_popcount64((nvi0_0 & nvj0_0 &  vk0_0) |
                                         ( (nvi0_0 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 & nvj0_1 &  vk0_0) |
                                         ( (nvi0_0 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_0 &  vk0_0) |
                                         ( (nvi0_1 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 & nvj0_1 &  vk0_0) |
                                         ( (nvi0_1 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_0 &  vk1_0) |
                                         ( (nvi1_0 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 & nvj1_1 &  vk1_0) |
                                         ( (nvi1_0 & nvj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_0 &  vk1_0) |
                                         ( (nvi1_1 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 & nvj1_1 &  vk1_0) |
                                         ( (nvi1_1 & nvj1_1 &  vk1_1) << 1 ));
            const int r010 = gm_popcount64((nvi0_0 &  vj0_0 & nvk0_0) |
                                         ( (nvi0_0 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 &  vj0_1 & nvk0_0) |
                                         ( (nvi0_0 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_0 & nvk0_0) |
                                         ( (nvi0_1 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_1 & nvk0_0) |
                                         ( (nvi0_1 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_0 & nvk1_0) |
                                         ( (nvi1_0 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_1 & nvk1_0) |
                                         ( (nvi1_0 &  vj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_0 & nvk1_0) |
                                         ( (nvi1_1 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_1 & nvk1_0) |
                                         ( (nvi1_1 &  vj1_1 & nvk1_1) << 1 ));
            const int r011 = gm_popcount64((nvi0_0 &  vj0_0 &  vk0_0) |
                                         ( (nvi0_0 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_0 &  vj0_1 &  vk0_0) |
                                         ( (nvi0_0 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_0 &  vk0_0) |
                                         ( (nvi0_1 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi0_1 &  vj0_1 &  vk0_0) |
                                         ( (nvi0_1 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_0 &  vk1_0) |
                                         ( (nvi1_0 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_0 &  vj1_1 &  vk1_0) |
                                         ( (nvi1_0 &  vj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_0 &  vk1_0) |
                                         ( (nvi1_1 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64((nvi1_1 &  vj1_1 &  vk1_0) |
                                         ( (nvi1_1 &  vj1_1 &  vk1_1) << 1 ));
            const int r100 = gm_popcount64(( vi0_0 & nvj0_0 & nvk0_0) |
                                         ( ( vi0_0 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 & nvj0_1 & nvk0_0) |
                                         ( ( vi0_0 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_0 & nvk0_0) |
                                         ( ( vi0_1 & nvj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_1 & nvk0_0) |
                                         ( ( vi0_1 & nvj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_0 & nvk1_0) |
                                         ( ( vi1_0 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_1 & nvk1_0) |
                                         ( ( vi1_0 & nvj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_0 & nvk1_0) |
                                         ( ( vi1_1 & nvj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_1 & nvk1_0) |
                                         ( ( vi1_1 & nvj1_1 & nvk1_1) << 1 ));
            const int r101 = gm_popcount64(( vi0_0 & nvj0_0 &  vk0_0) |
                                         ( ( vi0_0 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 & nvj0_1 &  vk0_0) |
                                         ( ( vi0_0 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_0 &  vk0_0) |
                                         ( ( vi0_1 & nvj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 & nvj0_1 &  vk0_0) |
                                         ( ( vi0_1 & nvj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_0 &  vk1_0) |
                                         ( ( vi1_0 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 & nvj1_1 &  vk1_0) |
                                         ( ( vi1_0 & nvj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_0 &  vk1_0) |
                                         ( ( vi1_1 & nvj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 & nvj1_1 &  vk1_0) |
                                         ( ( vi1_1 & nvj1_1 &  vk1_1) << 1 ));
            const int r110 = gm_popcount64(( vi0_0 &  vj0_0 & nvk0_0) |
                                         ( ( vi0_0 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 &  vj0_1 & nvk0_0) |
                                         ( ( vi0_0 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_0 & nvk0_0) |
                                         ( ( vi0_1 &  vj0_0 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_1 & nvk0_0) |
                                         ( ( vi0_1 &  vj0_1 & nvk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_0 & nvk1_0) |
                                         ( ( vi1_0 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_1 & nvk1_0) |
                                         ( ( vi1_0 &  vj1_1 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_0 & nvk1_0) |
                                         ( ( vi1_1 &  vj1_0 & nvk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_1 & nvk1_0) |
                                         ( ( vi1_1 &  vj1_1 & nvk1_1) << 1 ));
            const int r111 = gm_popcount64(( vi0_0 &  vj0_0 &  vk0_0) |
                                         ( ( vi0_0 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_0 &  vj0_1 &  vk0_0) |
                                         ( ( vi0_0 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_0 &  vk0_0) |
                                         ( ( vi0_1 &  vj0_0 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi0_1 &  vj0_1 &  vk0_0) |
                                         ( ( vi0_1 &  vj0_1 &  vk0_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_0 &  vk1_0) |
                                         ( ( vi1_0 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_0 &  vj1_1 &  vk1_0) |
                                         ( ( vi1_0 &  vj1_1 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_0 &  vk1_0) |
                                         ( ( vi1_1 &  vj1_0 &  vk1_1) << 1 )) +
                             gm_popcount64(( vi1_1 &  vj1_1 &  vk1_0) |
                                         ( ( vi1_1 &  vj1_1 &  vk1_1) << 1 ));

            /*---Accumulate---*/

//              printf("%i %i %i %i %i %i %i %i\n",
//              (int)r000, (int)r001, (int)r010, (int)r011,
//              (int)r100, (int)r101, (int)r110, (int)r111);
            sum.data[0] += GMTally1_encode(r000, r001);
            sum.data[1] += GMTally1_encode(r010, r011);
            sum.data[2] += GMTally1_encode(r100, r101);
            sum.data[3] += GMTally1_encode(r110, r111);
          } /*---for f---*/
          if (Env_all2all(env)) {
            GMMetrics_tally4x2_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                             sum, env);
          } else {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
          }
        }
      }
    }

    /* clang-format on */

    /*----------------------------------------*/
  } else /* if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) */ {
    /*----------------------------------------*/

    GMAssertAlways(GM_BOOL_FALSE
                 ? "logic error - code branch should never be executed"
                 : 0);

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMSectionInfo_destroy(si, env);;
}

/*===========================================================================*/

void gm_compute_numerators_3way_gpu_form_matV(
  const GMMirroredPointer* const vectors_I_buf,
  const GMMirroredPointer* const vectors_J_buf,
  GMMirroredPointer* const matV_buf,
  const int numfieldl,
  const int numpfieldl,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  GMEnv* const env) {

  /*--------------------*/
  /*---Populate leading columns of matV---*/
  /*--------------------*/

  /*----------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      /*---Operate on columns x_i and x_j elementwise---*/
      int f = 0;
      for (f = 0; f < numpfieldl; ++f) {
        const GMFloat a = ((GMFloat*)(vectors_I_buf->h))[f + numpfieldl*I];
        const GMFloat b = ((GMFloat*)(vectors_J_buf->h))[f + numpfieldl*J];
        ((GMFloat*)(matV_buf->h))[f + numpfieldl * I] = a < b ? a : b;
      }  //---for f---//
    }    //---for i---//
    /*----------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      const GMUInt64 nobits = 0;
      const GMUInt64 allbits = 0xffffffffffffffff;
      const GMUInt64 oddbits = 0x5555555555555555;
      const int num_seminibbles = 1 + (numfieldl-1) % 64;
      const GMUInt64 lastmask0 = num_seminibbles >= 32 ?
                                 allbits :
                                 allbits >> (64 - 2*num_seminibbles);
      const GMUInt64 lastmask1 = num_seminibbles <= 32 ?
                                 nobits :
                                 num_seminibbles == 64 ?
                                 allbits :
                                 allbits >> (128 - 2*num_seminibbles);
      /*---Operate on columns x_i and x_j elementwise---*/
      int f = 0;
      const int f_last = numpfieldl - 1;
      for (f = 0; f <= f_last; ++f) {
        const int indI = f + numpfieldl * I;
        const int indJ = f + numpfieldl * J;

        /*-----*/

        const GMUInt64 vI0 =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_I_buf->h))[indI].data[0]);
        const GMUInt64 vJ0 =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_J_buf->h))[indJ].data[0]);

        const GMUInt64 activebits0 = f < f_last ? allbits : lastmask0;
        const GMUInt64 vI0x = step_2way==0 ? vI0 | ~activebits0 : vI0;

        const GMUInt64  vI0_0 =   vI0x       & oddbits;
        const GMUInt64  vI0_1 =  (vI0x >> 1) & oddbits;
        const GMUInt64  vJ0_0 =   vJ0        & oddbits;
        const GMUInt64  vJ0_1 =  (vJ0  >> 1) & oddbits;
        const GMUInt64 nvI0_0 = ~ vI0x       & oddbits;
        const GMUInt64 nvI0_1 = ~(vI0x >> 1) & oddbits;
        /*
        const GMUInt64 nvJ0_0 = ~ vJ0        & oddbits;
        const GMUInt64 nvJ0_1 = ~(vJ0  >> 1) & oddbits;
        */

        const GMUInt64  vI0_match =
          step_2way==0 ?  nvI0_0 & nvI0_1  & oddbits : 
          step_2way==1 ? ( vI0_0 ^  vI0_1) & oddbits : 
                           vI0_0 &  vI0_1  & oddbits;
        const GMUInt64 nvI0_match =
          step_2way==0 ? ( vI0_0 |  vI0_1) & oddbits : 
          step_2way==1 ? ( vI0_0 ^ nvI0_1) & oddbits : 
                         (nvI0_0 | nvI0_1) & oddbits;

        const GMUInt64 r0_0 =  vI0_match & (vJ0_0 |  vJ0_1);
        const GMUInt64 r0_1 = nvI0_match | (vJ0_0 &  vJ0_1);
        const GMUInt64 r0 = r0_0 | (r0_1 << 1);
        ((GMBits2x64*)(matV_buf->h))[indI].data[0] = *(GMBits1_2x64*)&r0;

        /*-----*/

        const GMUInt64 vI1 =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_I_buf->h))[indI].data[1]);
        const GMUInt64 vJ1 =
            *(GMUInt64*)&(((GMBits2x64*)(vectors_J_buf->h))[indJ].data[1]);

        const GMUInt64 activebits1 = f < f_last ? allbits : lastmask1;
        const GMUInt64 vI1x = step_2way==0 ? vI1 | ~activebits1 : vI1;

        const GMUInt64  vI1_0 =   vI1x       & oddbits;
        const GMUInt64  vI1_1 =  (vI1x >> 1) & oddbits;
        const GMUInt64  vJ1_0 =   vJ1        & oddbits;
        const GMUInt64  vJ1_1 =  (vJ1  >> 1) & oddbits;
        const GMUInt64 nvI1_0 = ~ vI1x       & oddbits;
        const GMUInt64 nvI1_1 = ~(vI1x >> 1) & oddbits;
        /*
        const GMUInt64 nvJ1_0 = ~ vJ1        & oddbits;
        const GMUInt64 nvJ1_1 = ~(vJ1  >> 1) & oddbits;
        */

        const GMUInt64  vI1_match =
          step_2way==0 ?  nvI1_0 & nvI1_1  & oddbits : 
          step_2way==1 ? ( vI1_0 ^  vI1_1) & oddbits : 
                               vI1_0 &  vI1_1  & oddbits;
            const GMUInt64 nvI1_match =
          step_2way==0 ? ( vI1_0 |  vI1_1) & oddbits : 
          step_2way==1 ? ( vI1_0 ^ nvI1_1) & oddbits : 
                             (nvI1_0 | nvI1_1) & oddbits;

        const GMUInt64 r1_0 =  vI1_match & (vJ1_0 |  vJ1_1);
        const GMUInt64 r1_1 = nvI1_match | (vJ1_0 &  vJ1_1);
        const GMUInt64 r1 = r1_0 | (r1_1 << 1);
        ((GMBits2x64*)(matV_buf->h))[indI].data[1] = *(GMBits1_2x64*)&r1;
      }  //---for f---//
    }    //---for i---//
    /*----------*/
  } /*---Env_metric_type(env)---*/
  /*----------*/
}

/*===========================================================================*/

void gm_compute_numerators_3way_gpu_form_metrics(
  const GMMirroredPointer* const matM_IJ_buf,
  const GMMirroredPointer* const matM_JK_buf,
  const GMMirroredPointer* const matM_KIK_buf,
  const GMMirroredPointer* const matB_buf,
  GMMetrics* metrics,
  const int numvecl,
  const int J,
  const int step_2way,
  const int I_min,
  const int I_max,
  const int K_min,
  const int K_max,
  const int j_block,
  const int k_block,
  const _Bool sax0,
  const _Bool sax1,
  const _Bool is_part3,
  GMEnv* const env) {

  /*--------------------*/
  /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
  /*--------------------*/

  /*----------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI &&
      !Env_all2all(env)) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvecl * J];
      int K = 0;
      for (K = K_min; K < K_max; ++K) {
        const GMFloat min_JK =
            ((GMFloat*)(matM_JK_buf->h))[J + numvecl * K];
        const GMFloat min_KIK =
            ((GMFloat*)(matM_KIK_buf->h))[K + numvecl * I];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max * K];
        const GMFloat numerator = min_IJ + min_JK + min_KIK - min_IJK;
        const int i = I;
        const int j = J;
        const int k = K;
        GMMetrics_float_set_3(metrics, i, j, k, numerator, env);
      } /*---for K---*/
    }   /*---for I---*/
    /*----------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI &&
      Env_all2all(env)) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvecl*J];
      int K = 0;
      for (K = K_min; K < K_max; ++K) {
        const GMFloat min_JK =
            ((GMFloat*)(matM_JK_buf->h))[J + numvecl * K];
        const GMFloat min_KIK =
            is_part3 ? ((GMFloat*)(matM_KIK_buf->h))[K + numvecl * I]
                     : ((GMFloat*)(matM_KIK_buf->h))[I + numvecl * K];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max * K];
        const GMFloat numerator = min_IJ + min_JK + min_KIK - min_IJK;
        /* clang-format off */
        const int i = !is_part3 ?   I :
                           sax0 ?   J :
                           sax1 ?   I :
                        /* sax2 ?*/ K;
        const int j = !is_part3 ?   J :
                           sax0 ?   K :
                           sax1 ?   J :
                        /* sax2 ?*/ I;
        const int k = !is_part3 ?   K :
                           sax0 ?   I :
                           sax1 ?   K :
                        /* sax2 ?*/ J;
        /* clang-format on */
        GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                      numerator, env);
      } /*---for K---*/
    }   /*---for I---*/
    /*----------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             !Env_all2all(env)) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {

      int K = 0;
      for (K = K_min; K < K_max; ++K) {
        /*---This is the notall2all case -- has no axis permutation---*/

        const int i = I;
        const int j = J;
        const int k = K;

        GMTally4x2 numerator = step_2way==0 ? GMTally4x2_null() :
                        GMMetrics_tally4x2_get_3(metrics, i, j, k, env);

        GMTally1 r000, r001;
        GMTally1_decode(&r000, &r001, numerator.data[0]);
        GMTally1 r010, r011;
        GMTally1_decode(&r010, &r011, numerator.data[1]);
        GMTally1 r100, r101;
        GMTally1_decode(&r100, &r101, numerator.data[2]);
        GMTally1 r110, r111;
        GMTally1_decode(&r110, &r111, numerator.data[3]);

        const GMTally2x2 mB = ((GMTally2x2*)(matB_buf->h))[I + I_max * K];
        GMTally1 mB00, mB01;
        GMTally1_decode(&mB00, &mB01, mB.data[0]);
        GMTally1 mB10, mB11;
        GMTally1_decode(&mB10, &mB11, mB.data[1]);

        if (step_2way==0) {
          r000 += 2 * mB00;
          r001 += 2 * mB01;
          r010 += 2 * mB10;
          r011 += 2 * mB11;
        } else if (step_2way==1) {
          r000 += mB00;
          r001 += mB01;
          r010 += mB10;
          r011 += mB11;
          r100 += mB00;
          r101 += mB01;
          r110 += mB10;
          r111 += mB11;
        } else /*---step_2way==2---*/ {
          r100 += 2 * mB00;
          r101 += 2 * mB01;
          r110 += 2 * mB10;
          r111 += 2 * mB11;
        }
        /*---NOTE: pay attention to order here---*/
        numerator.data[0] = GMTally1_encode(r000, r001);
        numerator.data[1] = GMTally1_encode(r010, r011);
        numerator.data[2] = GMTally1_encode(r100, r101);
        numerator.data[3] = GMTally1_encode(r110, r111);
        GMMetrics_tally4x2_set_3(metrics, i, j, k, numerator, env);
      } /*---for K---*/
    }   /*---for I---*/
    /*----------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CCC &&
             Env_all2all(env)) {
    /*----------*/
    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      int K = 0;
      for (K = K_min; K < K_max; ++K) {
/*---For the permuted case,
 1) pay attention to KIK access
 2) swap 01 and 10 if needed.
---*/

        /* clang-format off */
        const int i = !is_part3 ?   I :
                           sax0 ?   J :
                           sax1 ?   I :
                        /* sax2 ?*/ K;
        const int j = !is_part3 ?   J :
                           sax0 ?   K :
                           sax1 ?   J :
                        /* sax2 ?*/ I;
        const int k = !is_part3 ?   K :
                           sax0 ?   I :
                           sax1 ?   K :
                        /* sax2 ?*/ J;
        /* clang-format on */

        GMTally4x2 numer = step_2way==0 ? GMTally4x2_null() :
          GMMetrics_tally4x2_get_all2all_3(metrics, i, j, k,
                                           j_block, k_block, env);
        GMTally1 r000_permuted, r001_permuted;
        GMTally1_decode(&r000_permuted, &r001_permuted, numer.data[0]);
        GMTally1 r010_permuted, r011_permuted;
        GMTally1_decode(&r010_permuted, &r011_permuted, numer.data[1]);
        GMTally1 r100_permuted, r101_permuted;
        GMTally1_decode(&r100_permuted, &r101_permuted, numer.data[2]);
        GMTally1 r110_permuted, r111_permuted;
        GMTally1_decode(&r110_permuted, &r111_permuted, numer.data[3]);

        const GMTally2x2 mB = ((GMTally2x2*)(matB_buf->h))[I + I_max * K];
        GMTally1 mB00, mB01;
        GMTally1_decode(&mB00, &mB01, mB.data[0]);
        GMTally1 mB10, mB11;
        GMTally1_decode(&mB10, &mB11, mB.data[1]);

        /* clang-format off */
        int r000 = r000_permuted;

        int r100 = !is_part3 ?   r100_permuted :
                        sax0 ?   r001_permuted :
                        sax1 ?   r100_permuted :
                     /* sax2 ?*/ r010_permuted;
        int r010 = !is_part3 ?   r010_permuted :
                        sax0 ?   r100_permuted :
                        sax1 ?   r010_permuted :
                     /* sax2 ?*/ r001_permuted;
        int r001 = !is_part3 ?   r001_permuted :
                        sax0 ?   r010_permuted :
                        sax1 ?   r001_permuted :
                     /* sax2 ?*/ r100_permuted;

        int r011 = !is_part3 ?   r011_permuted :
                        sax0 ?   r110_permuted :
                        sax1 ?   r011_permuted :
                     /* sax2 ?*/ r101_permuted;
        int r101 = !is_part3 ?   r101_permuted :
                        sax0 ?   r011_permuted :
                        sax1 ?   r101_permuted :
                     /* sax2 ?*/ r110_permuted;
        int r110 = !is_part3 ?   r110_permuted :
                        sax0 ?   r101_permuted :
                        sax1 ?   r110_permuted :
                     /* sax2 ?*/ r011_permuted;

        int r111 = r111_permuted;
        /* clang-format on */

        if (step_2way==0) {
          r000 += 2 * mB00;
          r001 += 2 * mB01;
          r010 += 2 * mB10;
          r011 += 2 * mB11;
        } else if (step_2way==1) {
          r000 += mB00;
          r001 += mB01;
          r010 += mB10;
          r011 += mB11;
          r100 += mB00;
          r101 += mB01;
          r110 += mB10;
          r111 += mB11;
        } else /*---step_2way==2---*/ {
          r100 += 2 * mB00;
          r101 += 2 * mB01;
          r110 += 2 * mB10;
          r111 += 2 * mB11;
        }

        /* clang-format off */
        r000_permuted = r000;

        r100_permuted = !is_part3 ?   r100 :
                             sax0 ?   r010 :
                             sax1 ?   r100 :
                          /* sax2 ?*/ r001;
        r010_permuted = !is_part3 ?   r010 :
                             sax0 ?   r001 :
                             sax1 ?   r010 :
                          /* sax2 ?*/ r100;
        r001_permuted = !is_part3 ?   r001 :
                             sax0 ?   r100 :
                             sax1 ?   r001 :
                          /* sax2 ?*/ r010;

        r011_permuted = !is_part3 ?   r011 :
                             sax0 ?   r101 :
                             sax1 ?   r011 :
                          /* sax2 ?*/ r110;
        r101_permuted = !is_part3 ?   r101 :
                             sax0 ?   r110 :
                             sax1 ?   r101 :
                          /* sax2 ?*/ r011;
        r110_permuted = !is_part3 ?   r110 :
                             sax0 ?   r011 :
                             sax1 ?   r110 :
                          /* sax2 ?*/ r101;

        r111_permuted = r111;
        /* clang-format on */

        /*---NOTE: pay attention to order here---*/

        numer.data[0] = GMTally1_encode(r000_permuted, r001_permuted);
        numer.data[1] = GMTally1_encode(r010_permuted, r011_permuted);
        numer.data[2] = GMTally1_encode(r100_permuted, r101_permuted);
        numer.data[3] = GMTally1_encode(r110_permuted, r111_permuted);
        GMMetrics_tally4x2_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                         numer, env);

      } /*---for K---*/
    }   /*---for I---*/
    /*----------*/
  } else {
    /*----------*/
    GMAssertAlways(GM_BOOL_FALSE);
    /*----------*/
  } /*---Env_metric_type(env)---*/
  /*----------*/
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way gpu---*/

void gm_compute_numerators_3way_gpu_start(GMVectors* vectors_i,
                                          GMVectors* vectors_j,
                                          GMVectors* vectors_k,
                                          GMMetrics* metrics,
                                          GMMirroredPointer* vectors_i_buf,
                                          GMMirroredPointer* vectors_j_buf,
                                          GMMirroredPointer* vectors_k_buf,
                                          int j_block,
                                          int k_block,
                                          GMEnv* env) {
  GMAssertAlways(vectors_i != NULL);
  GMAssertAlways(vectors_j != NULL);
  GMAssertAlways(vectors_k != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors_i_buf != NULL);
  GMAssertAlways(vectors_j_buf != NULL);
  GMAssertAlways(vectors_k_buf != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == j_block &&
                   Env_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == k_block &&
                   Env_proc_num_vector_i(env) != j_block));
  GMAssertAlways(Env_compute_method(env) == GM_COMPUTE_METHOD_GPU);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;
  const int numpfieldl = vectors_i->num_packedval_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  const _Bool need_2way = Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  /*---First get the required 2-way ij, jk, ik metrics---*/
  /*----------------------------------------*/

  /*--------------------*/
  /*---Compute i_block - j_block PROD---*/
  /*--------------------*/

  /*---Allocate magma CPU/GPU memory for M = X^T PROD X---*/

  GMMirroredPointer mat_buf_tmp_value0 =
    gm_malloc_magma(numvecl * (size_t)numvecl, env);
  GMMirroredPointer mat_buf_tmp_value1 =
    gm_malloc_magma(numvecl * (size_t)numvecl, env);
  GMMirroredPointer* mat_buf_tmp[2] =
    {&mat_buf_tmp_value0, &mat_buf_tmp_value1};

  GMMirroredPointer matM_ij_buf_value = need_2way
    ? gm_malloc_magma(numvecl * (size_t)numvecl, env)
    : GMMirroredPointer_null();  // M = X^T PROD X
  GMMirroredPointer* const matM_ij_buf = &matM_ij_buf_value;

  if (need_2way) {
    GMMirroredPointer* matM_ij_buf_local =
        Env_num_proc_field(env) == 1 ? matM_ij_buf : mat_buf_tmp[0];

    /*---Initialize result matrix to zero (apparently magma requires)---*/

    gm_magma_set_matrix_zero_start(matM_ij_buf_local, numvecl, numvecl, env);

    /*---Perform pseudo matrix-matrix min product for M = X^T PROD X---*/

    gm_magma_gemm_start(numvecl, numvecl, numpfieldl,
                        vectors_i_buf->d, numpfieldl,
                        vectors_j_buf->d, numpfieldl,
                        matM_ij_buf_local->d, numvecl, env);
    gm_compute_wait(env);

    /*---Copy matM_ij from GPU---*/

    gm_get_matrix_start(matM_ij_buf_local, numvecl, numvecl, env);
    gm_get_matrix_wait(env);

    // TODO - outline into gm_allreduce_metrics
    if (Env_num_proc_field(env) > 1) {
      int mpi_code = 0;
      mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
      mpi_code = MPI_Allreduce(matM_ij_buf_local->h, matM_ij_buf->h,
                               numvecl * (size_t)numvecl, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
    // TODO - outline into gm_allreduce_metrics
  }

  /*--------------------*/
  /*---Compute j_block - k_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  GMMirroredPointer matM_jk_buf_value = need_2way && !si->is_part1
    ? gm_malloc_magma(numvecl * (size_t)numvecl, env)
    : GMMirroredPointer_null();
  GMMirroredPointer* const matM_jk_buf =
      !si->is_part1 ? &matM_jk_buf_value : matM_ij_buf;

  if (need_2way && !si->is_part1) {
    GMMirroredPointer* matM_jk_buf_local =
        Env_num_proc_field(env) == 1 ? matM_jk_buf : mat_buf_tmp[0];

    gm_magma_set_matrix_zero_start(matM_jk_buf_local, numvecl, numvecl, env);

    gm_magma_gemm_start(numvecl, numvecl, numpfieldl,
                        vectors_j_buf->d, numpfieldl,
                        vectors_k_buf->d, numpfieldl,
                        matM_jk_buf_local->d, numvecl, env);
    gm_compute_wait(env);

    gm_get_matrix_start(matM_jk_buf_local, numvecl, numvecl, env);
    gm_get_matrix_wait(env);

    // TODO - outline into gm_allreduce_metrics
    if (Env_num_proc_field(env) > 1) {
      int mpi_code = 0;
      mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
      mpi_code = MPI_Allreduce(matM_jk_buf_local->h, matM_jk_buf->h,
                               numvecl * (size_t)numvecl, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
    // TODO - outline into gm_allreduce_metrics
  }

  /*--------------------*/
  /*---Compute k_block - i_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  /*---NOTE: for Part 3, this is indexed directly as (k,i).
       Otherwise, it is indexed through an alias as (i,k)---*/

  GMMirroredPointer matM_kik_buf_value = need_2way && si->is_part3
    ? gm_malloc_magma(numvecl * (size_t)numvecl, env)
    : GMMirroredPointer_null();

  GMMirroredPointer* const matM_kik_buf = si->is_part3
    ? &matM_kik_buf_value : matM_ij_buf;

  if (need_2way && si->is_part3) {
    GMMirroredPointer* matM_kik_buf_local =
        Env_num_proc_field(env) == 1 ? matM_kik_buf : mat_buf_tmp[0];

    gm_magma_set_matrix_zero_start(matM_kik_buf_local, numvecl, numvecl, env);

    gm_magma_gemm_start(numvecl, numvecl, numpfieldl,
                        vectors_k_buf->d, numpfieldl,
                        vectors_i_buf->d, numpfieldl,
                        matM_kik_buf_local->d, numvecl, env);
    gm_compute_wait(env);

    gm_get_matrix_start(matM_kik_buf_local, numvecl, numvecl, env);
    gm_get_matrix_wait(env);

    // TODO - outline into gm_allreduce_metrics
    if (Env_num_proc_field(env) > 1) {
      int mpi_code = 0;
      mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
      mpi_code = MPI_Allreduce(matM_kik_buf_local->h, matM_kik_buf->h,
                               numvecl * (size_t)numvecl, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
    // TODO - outline into gm_allreduce_metrics
  } /*---is_part3---*/

  /*----------------------------------------*/
  /*---Now compute ijk piece, via an outer loop over j values---*/
  /*----------------------------------------*/

  /*---Allocate magma CPU/GPU memory for matrices V and B---*/
  /*
     V = elementwise OP of one vector with each of the rest of the vectors.
     For the jth iteration, the ith column of V is the elementwise OP
       of vectors i and j.
     B = X^T PROD V = three way PROD.
  */
  GMMirroredPointer matV_buf_value0 =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);
  GMMirroredPointer matV_buf_value1 =
      gm_malloc_magma(numvecl * (size_t)numpfieldl, env);
  GMMirroredPointer* matV_buf[2] = {&matV_buf_value0, &matV_buf_value1};

  GMMirroredPointer matB_buf_value0 =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);
  GMMirroredPointer matB_buf_value1 =
      gm_malloc_magma(numvecl * (size_t)numvecl, env);
  GMMirroredPointer* matB_buf[2] = {&matB_buf_value0, &matB_buf_value1};

  /*---Set up pointers to permute the access of axes for Part 3---*/
  /*---We use capitals I, J, K here to denote the PERMUTED axes---*/

  const _Bool sax0 = si->section_axis == 0;
  const _Bool sax1 = si->section_axis == 1;
  const _Bool sax2 = si->section_axis == 2;

  const _Bool is_ijk = !si->is_part3 ? GM_BOOL_TRUE : sax1;
  const _Bool is_kij = !si->is_part3 ? GM_BOOL_FALSE : sax0;
  const _Bool is_jki = !si->is_part3 ? GM_BOOL_FALSE : sax2;

#if 0
  GMMirroredPointer* const vectors_I_buf = !is_part3 ? vectors_i_buf :
                                                sax0 ? vectors_k_buf :
                                                sax1 ? vectors_i_buf :
                                                sax2 ? vectors_j_buf : 0;
  
  GMMirroredPointer* const vectors_J_buf = !is_part3 ? vectors_j_buf :
                                                sax0 ? vectors_i_buf :
                                                sax1 ? vectors_j_buf :
                                                sax2 ? vectors_k_buf : 0;
  
  GMMirroredPointer* const vectors_K_buf = !is_part3 ? vectors_k_buf :
                                                sax0 ? vectors_j_buf :
                                                sax1 ? vectors_k_buf :
                                                sax2 ? vectors_i_buf : 0;
#endif

  /* clang-format off */
  GMMirroredPointer* const vectors_I_buf = is_ijk ? vectors_i_buf :
                                           is_kij ? vectors_k_buf :
                                           is_jki ? vectors_j_buf : 0;
 
  GMMirroredPointer* const vectors_J_buf = is_ijk ? vectors_j_buf :
                                           is_kij ? vectors_i_buf :
                                           is_jki ? vectors_k_buf : 0;
 
  GMMirroredPointer* const vectors_K_buf = is_ijk ? vectors_k_buf :
                                           is_kij ? vectors_j_buf :
                                           is_jki ? vectors_i_buf : 0;
  
  /*---NOTE: must pay attention that these permuted matrices
       are indexed the right way by the permuted indices---*/
 
//TODO - use is_ijk etc. 
  GMMirroredPointer* const matM_IJ_buf  = !si->is_part3 ? matM_ij_buf  :
                                                   sax0 ? matM_kik_buf :
                                                   sax1 ? matM_ij_buf  :
                                                   sax2 ? matM_jk_buf  : 0;
  
  GMMirroredPointer* const matM_JK_buf  = !si->is_part3 ? matM_jk_buf  :
                                                   sax0 ? matM_ij_buf  :
                                                   sax1 ? matM_jk_buf  :
                                                   sax2 ? matM_kik_buf : 0;
  
  GMMirroredPointer* const matM_KIK_buf = !si->is_part3 ? matM_kik_buf :
                                                   sax0 ? matM_jk_buf  :
                                                   sax1 ? matM_kik_buf :
                                                   sax2 ? matM_ij_buf  : 0;
  /* clang-format on */

  /*---Process all combinations with nested loops in j, i, k---*/

  const int J_min = si->is_part3 ? (si->section_num + 0) * numvecl / 6 : 0;
  const int J_max = si->is_part3 ? (si->section_num + 1) * numvecl / 6 : numvecl;

  //NEW const int J_min = (si->section_num + 0) * numvecl / 6;
  //NEW const int J_max = (si->section_num + 1) * numvecl / 6;

  const int J_count = J_max - J_min;

  /*--------------------*/
  /*---Collapsed loops over J and over 2-way steps---*/
  /*--------------------*/

  const int num_step_2way = Env_metric_type(env) == GM_METRIC_TYPE_CCC ?  3 : 1;
  const int num_step = J_count*num_step_2way;

  const int extra_step = 1;
//  const int fill_step = 0 - extra_step;
//  const int drain_step = num_step + extra_step;

  MPI_Request mpi_requests[2];

  int step_num = 0;

  for (step_num = 0-extra_step; step_num < num_step+extra_step*2; ++step_num) {

    //==========

    //const int step_2way = (step_num + 2*num_step_2way) % num_step_2way;
    const int J = J_min + (step_num + 2*num_step_2way) / num_step_2way - 2;
    const int I_min = 0;
    const int I_max = si->is_part1 ? J : numvecl;
    const int K_min = si->is_part3 ? 0 : J + 1;
    const int K_max = numvecl;
    const _Bool empty = I_min >= I_max || K_min >= K_max;
    const _Bool is_compute_step = step_num >= 0 && step_num < num_step;
    const _Bool do_compute = is_compute_step && !empty;
    const int index_01 = (step_num + 2) % 2;
    GMMirroredPointer* matB_buf_local =
      Env_num_proc_field(env) == 1 ? matB_buf[index_01] : mat_buf_tmp[index_01];

    //==========

    const int step_num_prev = step_num - 1;
    //const int step_2way_prev =
    //  (step_num_prev + 2*num_step_2way) % num_step_2way;
    const int J_prev = J_min +
      (step_num_prev + 2*num_step_2way) / num_step_2way - 2;
    const int I_min_prev = 0;
    const int I_max_prev = si->is_part1 ? J_prev : numvecl;
    const int K_min_prev = si->is_part3 ? 0 : J_prev + 1;
    const int K_max_prev = numvecl;
    const _Bool empty_prev = I_min_prev >= I_max_prev ||
                             K_min_prev >= K_max_prev;
    const _Bool is_compute_step_prev = step_num_prev >= 0 &&
                                       step_num_prev < num_step;
    const _Bool do_compute_prev = is_compute_step_prev && !empty_prev;
    const int index_01_prev = (step_num_prev + 2) % 2;
    GMMirroredPointer* matB_buf_local_prev = Env_num_proc_field(env) == 1 ?
      matB_buf[index_01_prev] : mat_buf_tmp[index_01_prev];

    //==========

    const int step_num_prevprev = step_num - 2;
    const int step_2way_prevprev =
      (step_num_prevprev + 3*num_step_2way) % num_step_2way;
    const int J_prevprev = J_min +
      (step_num_prevprev + 3*num_step_2way) / num_step_2way - 3;
    const int I_min_prevprev = 0;
    const int I_max_prevprev = si->is_part1 ? J_prevprev : numvecl;
    const int K_min_prevprev = si->is_part3 ? 0 : J_prevprev + 1;
    const int K_max_prevprev = numvecl;
    const _Bool empty_prevprev = I_min_prevprev >= I_max_prevprev ||
                                 K_min_prevprev >= K_max_prevprev;
    const _Bool is_compute_step_prevprev = step_num_prevprev >= 0 &&
                                           step_num_prevprev < num_step;
    const _Bool do_compute_prevprev = is_compute_step_prevprev && !empty_prevprev;
    const int index_01_prevprev = (step_num_prevprev + 4) % 2;
    //GMMirroredPointer* matB_buf_local_prevprev = Env_num_proc_field(env) == 1 ?
    //  matB_buf[index_01_prevprev] : mat_buf_tmp[index_01_prevprev];

    //==========

    const int step_num_next = step_num + 1;
    const int step_2way_next = step_num_next % num_step_2way;
    const int J_next = J_min + step_num_next / num_step_2way;
    const int I_min_next = 0;
    const int I_max_next = si->is_part1 ? J_next : numvecl;
    const int K_min_next = si->is_part3 ? 0 : J_next + 1;
    const int K_max_next = numvecl;
    const _Bool empty_next = I_min_next >= I_max_next ||
                             K_min_next >= K_max_next;
    const _Bool is_compute_step_next = step_num_next >= 0 &&
                                       step_num_next < num_step;
    const _Bool do_compute_next = is_compute_step_next && !empty_next;
    const int index_01_next = (step_num_next + 2) % 2;

    //==========

    if (do_compute_next) {
      /*---Populate leading columns of matV---*/
      gm_compute_numerators_3way_gpu_form_matV(
          vectors_I_buf, vectors_J_buf,
          matV_buf[index_01_next],
          numfieldl, numpfieldl,
          J_next, step_2way_next, I_min_next, I_max_next,
          env);
    } /*---if do_compute---*/

    //==========

    if (do_compute) {
      /*---Send matrix matV to GPU - WAIT---*/
      gm_set_matrix_wait(env);
    } /*---if do_compute---*/

    //==========

    if (do_compute_prev) {
      /*---Perform pseudo mat X matt matB = matV^T PROD X - WAIT---*/
      gm_compute_wait(env);
    } /*---if do_compute---*/

    //==========

    if (do_compute_next) {
      /*---Send matrix matV to GPU - START---*/
      gm_set_matrix_start(matV_buf[index_01_next], numpfieldl, I_max_next, env);
    } /*---if do_compute---*/

    //==========

    if (do_compute_prev) {
      /*---Copy result matrix matB from GPU - START---*/
      gm_get_matrix_start(matB_buf_local_prev, I_max_prev, numvecl, env);
    } /*---if do_compute---*/

    //==========

    if (do_compute) {
      /*---Initialize result matrix to zero (apparently magma requires)---*/
      gm_magma_set_matrix_zero_start(matB_buf_local, numvecl, I_max, env);
      /*---Perform pseudo mat X mat matB = matV^T PROD X - START---*/
      gm_magma_gemm_start(I_max, numvecl, numpfieldl,
                          matV_buf[index_01]->d, numpfieldl,
                          vectors_K_buf->d, numpfieldl,
                          matB_buf_local->d, I_max, env);
    } /*---if do_compute---*/

    //==========

    if (do_compute_prev) {
      /*---Copy result matrix matB from GPU - WAIT---*/
      gm_get_matrix_wait(env);
    } /*---if do_compute---*/

    //==========

    if (do_compute_prevprev) {
      if (Env_num_proc_field(env) > 1) {
        MPI_Status mpi_status;
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Wait(&(mpi_requests[index_01_prevprev]), &mpi_status);
        GMAssertAlways(mpi_code == MPI_SUCCESS);
      }
    }

    //==========

    if (do_compute_prev) {
      // TODO - outline into gm_allreduce_metrics
      if (Env_num_proc_field(env) > 1) {
        //TODO: fix this properly.
        int multiplier = Env_metric_type(env) == GM_METRIC_TYPE_CCC ? 2 : 1;
        int mpi_code = 0;
        mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
        mpi_code = MPI_Iallreduce(matB_buf_local_prev->h,
                                  matB_buf[index_01_prev]->h,
                                  numvecl * (size_t)numvecl * multiplier,
                                  GM_MPI_FLOAT, MPI_SUM,
                                  Env_mpi_comm_field(env),
                                  &(mpi_requests[index_01_prev]));
        GMAssertAlways(mpi_code == MPI_SUCCESS);
      }
      // TODO - outline into gm_allreduce_metrics
    } /*---if do_compute---*/

    //==========

    //---NOTE: matB_buf[index_01_prevprev]->d is being updated now but matB_buf[index_01_prevprev]->h is ok.

    if (do_compute_prevprev) {
      /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
      gm_compute_numerators_3way_gpu_form_metrics(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
          matB_buf[index_01_prevprev],
          metrics, numvecl,
          J_prevprev,
          step_2way_prevprev,
          I_min_prevprev,
          I_max_prevprev,
          K_min_prevprev,
          K_max_prevprev,
          j_block, k_block, sax0, sax1, si->is_part3,
          env);
    } /*---if do_compute---*/

    //==========

  } /*---for step_num---*/

  /*--------------------*/
  /*---Free memory---*/
  /*--------------------*/

  if (need_2way) {
    gm_free_magma(&matM_ij_buf_value, env);
    if (!si->is_part1) {
      gm_free_magma(&matM_jk_buf_value, env);
    }
    if (si->is_part3) {
      gm_free_magma(&matM_kik_buf_value, env);
    }
  }
  gm_free_magma(mat_buf_tmp[0], env);
  gm_free_magma(mat_buf_tmp[1], env);
  gm_free_magma(matV_buf[0], env);
  gm_free_magma(matV_buf[1], env);
  gm_free_magma(matB_buf[0], env);
  gm_free_magma(matB_buf[1], env);
  GMSectionInfo_destroy(si, env);
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way generic---*/

void gm_compute_numerators_3way_start(GMVectors* vectors_i,
                                      GMVectors* vectors_j,
                                      GMVectors* vectors_k,
                                      GMMetrics* metrics,
                                      GMMirroredPointer* vectors_i_buf,
                                      GMMirroredPointer* vectors_j_buf,
                                      GMMirroredPointer* vectors_k_buf,
                                      int j_block,
                                      int k_block,
                                      GMEnv* env) {
  GMAssertAlways(vectors_i != NULL);
  GMAssertAlways(vectors_j != NULL);
  GMAssertAlways(vectors_k != NULL);
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vectors_i_buf != NULL);
  GMAssertAlways(vectors_j_buf != NULL);
  GMAssertAlways(vectors_k_buf != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == j_block &&
                   Env_proc_num_vector_i(env) != k_block));
  GMAssertAlways(!(Env_proc_num_vector_i(env) == k_block &&
                   Env_proc_num_vector_i(env) != j_block));
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  /*----------------------------------------*/
  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    gm_compute_numerators_3way_gpu_start(vectors_i, vectors_j, vectors_k,
                                         metrics, vectors_i_buf, vectors_j_buf,
                                         vectors_k_buf, j_block, k_block, env);
    /*----------------------------------------*/
  } else /*---(Env_compute_method(env) != GM_COMPUTE_METHOD_GPU)---*/ {
    /*----------------------------------------*/
    switch (Env_metric_type(env)) {
      /*----------------------------------------*/
      case GM_METRIC_TYPE_SORENSON: {
        /*----------------------------------------*/
        GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
      } break;
      /*----------------------------------------*/
      case GM_METRIC_TYPE_CZEKANOWSKI: {
        /*----------------------------------------*/
        gm_compute_czekanowski_numerators_3way_nongpu_start(
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block, env);
      } break;
      /*----------------------------------------*/
      case GM_METRIC_TYPE_CCC: {
        /*----------------------------------------*/
        gm_compute_ccc_numerators_3way_nongpu_start(
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block, env);
      } break;
      /*----------------------------------------*/
      default:
        /*----------------------------------------*/
        /*---Should never get here---*/
        GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    } /*---case---*/
    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 3-way Czek---*/

void gm_compute_czekanowski_3way_combine(GMMetrics* metrics,
                                         GMFloat* __restrict__ vector_sums_i,
                                         GMFloat* __restrict__ vector_sums_j,
                                         GMFloat* __restrict__ vector_sums_k,
                                         int j_block,
                                         int k_block,
                                         GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(Env_proc_num_vector_i(env) != j_block || j_block == k_block);
  GMAssertAlways(Env_proc_num_vector_i(env) != k_block || j_block == k_block);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  const int i_block = Env_proc_num_vector_i(env);

  const int numvecl = metrics->num_vector_local;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  if (Env_all2all(env)) {
    /*----------------------------------------*/

//TODO: combine the following 3 cases using same logic as gm_compute_czekanowski_numerators_3way_nongpu_start
    /*----------------------------------------*/
    if (si->is_part1) {
      /*----------------------------------------*/

      int k = 0;
      for (k = 0; k < numvecl; ++k) {
        int j = 0;
        for (j = 0; j < k; ++j) {
          int i = 0;
          for (i = 0; i < j; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_block, k_block, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_i[j] + vector_sums_i[k];
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                          3 * numerator / (2 * denominator),
                                          env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/

    /*----------------------------------------*/
    } else if (si->is_part3) {
      /*----------------------------------------*/

      int k = 0;
      for (k = si->k_lb; k < si->k_ub; ++k) {
        int j = 0;
        for (j = si->j_lb; j < si->j_ub; ++j) {
          int i = 0;
          for (i = si->i_lb; i < si->i_ub; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_block, k_block, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_j[j] + vector_sums_k[k];
            const GMFloat value =
                ((GMFloat)3) * numerator / (((GMFloat)2) * denominator);
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                          value, env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/

    /*----------------------------------------*/
    } else /* if (j_block == k_block) */ {
    /*----------------------------------------*/

      int k = 0;
      for (k = 0; k < numvecl; ++k) {
        int j = 0;
        for (j = 0; j < k; ++j) {
          int i = 0;
          for (i = 0; i < numvecl; ++i) {
            const GMFloat numerator = GMMetrics_float_get_all2all_3(
                metrics, i, j, k, j_block, k_block, env);
            const GMFloat denominator =
                vector_sums_i[i] + vector_sums_j[j] + vector_sums_j[k];
            GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                          3 * numerator / (2 * denominator),
                                          env);
          } /*---for i---*/
        }   /*---for j---*/
      }     /*---for k---*/

    }       /*---if---*/

    /*----------------------------------------*/
  } else /*---! Env_all2all(env)---*/ {
    /*----------------------------------------*/

    int i = 0;
    for (i = 0; i < numvecl; ++i) {
      int j = 0;
      for (j = i + 1; j < numvecl; ++j) {
        int k = 0;
        for (k = j + 1; k < numvecl; ++k) {
          const GMFloat numerator =
              GMMetrics_float_get_3(metrics, i, j, k, env);
          /*---Don't use two different pointers pointing to the same thing---*/
          const GMFloat denominator =
              vector_sums_i[i] + vector_sums_i[j] + vector_sums_i[k];
          GMMetrics_float_set_3(metrics, i, j, k,
                                3 * numerator / (2 * denominator), env);
        } /*---for k---*/
      }   /*---for j---*/
    }     /*---for i---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMSectionInfo_destroy(si, env);
}
//TODO: here and throughout, check that pointer aliasing is ok.

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 3-way CCC---*/

void gm_compute_ccc_3way_combine(GMMetrics* metrics,
                                 GMFloat* __restrict__ vector_sums_i,
                                 GMFloat* __restrict__ vector_sums_j,
                                 GMFloat* __restrict__ vector_sums_k,
                                 int j_block,
                                 int k_block,
                                 GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(Env_proc_num_vector_i(env) != j_block || j_block == k_block);
  GMAssertAlways(Env_proc_num_vector_i(env) != k_block || j_block == k_block);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;

  const int i_block = Env_proc_num_vector_i(env);

  int i = 0;
  int j = 0;
  int k = 0;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  if (Env_all2all(env)) {
    /*----------------------------------------*/

    /*---Compute tetrahedron, triang prism or block section---*/
    /*---Store multipliers---*/

    for (k = si->k_lb; k < si->k_ub; ++k) {
      const GMTally1 sk_1 = (GMTally1)(vector_sums_k[k]);
      GMAssert((GMFloat)sk_1 == vector_sums_k[k]);
      const int j_max = si->is_part3 ? si->j_ub : k;
      for (j = si->j_lb; j < j_max; ++j) {
        const GMTally1 sj_1 = (GMTally1)(vector_sums_j[j]);
        GMAssert((GMFloat)sj_1 == vector_sums_j[j]);
        const int i_max = si->is_part1 ? j : si->i_ub;
        for (i = si->i_lb; i < i_max; ++i) {
          const GMTally1 si_1 = (GMTally1)(vector_sums_i[i]);
          GMAssert((GMFloat)si_1 == vector_sums_i[i]);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
          GMMetrics_float3_M_set_all2all_3(metrics, i, j, k, j_block, k_block,
              si1_sj1_sk1, env);
        } /*---for k---*/
      }   /*---for j---*/
    }     /*---for i---*/

    /*----------------------------------------*/
  } else /*---! Env_all2all(env)---*/ {
    /*----------------------------------------*/

    /*---No off-proc all2all: compute tetrahedron of values---*/
    /*---Store multipliers---*/

    for (k = 0; k < numvecl; ++k) {
      const GMTally1 sk_1 = (GMTally1)(vector_sums_i[k]);
      GMAssert((GMFloat)sk_1 == vector_sums_i[k]);
      for (j = 0; j < k; ++j) {
        const GMTally1 sj_1 = (GMTally1)(vector_sums_i[j]);
        GMAssert((GMFloat)sj_1 == vector_sums_i[j]);
        for (i = 0; i < j; ++i) {
          const GMTally1 si_1 = (GMTally1)(vector_sums_i[i]);
          GMAssert((GMFloat)si_1 == vector_sums_i[i]);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
          GMMetrics_float3_M_set_3(metrics, i, j, k, si1_sj1_sk1, env);
        } /*---for k---*/
      }   /*---for j---*/
    }     /*---for i---*/

    /*----------------------------------------*/
  } /*---if---*/
  /*----------------------------------------*/

  GMSectionInfo_destroy(si, env);
}

/*===========================================================================*/
/*---Combine nums and denoms on CPU to get final result, 3-way generic---*/

void gm_compute_3way_combine(GMMetrics* metrics,
                             GMVectorSums* vector_sums_i,
                             GMVectorSums* vector_sums_j,
                             GMVectorSums* vector_sums_k,
                             int j_block,
                             int k_block,
                             GMEnv* env) {
  GMAssertAlways(metrics != NULL);
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);
  GMAssertAlways(env != NULL);
  GMAssertAlways(j_block >= 0 && j_block < Env_num_block_vector(env));
  GMAssertAlways(k_block >= 0 && k_block < Env_num_block_vector(env));
  GMAssertAlways(Env_proc_num_vector_i(env) != j_block || j_block == k_block);
  GMAssertAlways(Env_proc_num_vector_i(env) != k_block || j_block == k_block);
  GMAssertAlways(Env_num_way(env) == GM_NUM_WAY_3);

  switch (Env_metric_type(env)) {
    /*----------------------------------------*/
    case GM_METRIC_TYPE_SORENSON: {
      /*----------------------------------------*/
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CZEKANOWSKI: {
      /*----------------------------------------*/
      gm_compute_czekanowski_3way_combine(
          metrics, (GMFloat*)vector_sums_i->data, (GMFloat*)vector_sums_j->data,
          (GMFloat*)vector_sums_k->data, j_block, k_block, env);
    } break;
    /*----------------------------------------*/
    case GM_METRIC_TYPE_CCC: {
      /*----------------------------------------*/
      gm_compute_ccc_3way_combine(
          metrics, (GMFloat*)vector_sums_i->data, (GMFloat*)vector_sums_j->data,
          (GMFloat*)vector_sums_k->data, j_block, k_block, env);
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
