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

void gm_vectors_to_buf(GMMirroredPointer* vectors_buf,
                       GMVectors* vectors,
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
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
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
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);

  /*---Initializations---*/

  const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  int i = 0;
  int j = 0;
  int k = 0;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  /*---TODO: address aliasing---*/
  GMFloat* __restrict__ vector_sums_i_ = (GMFloat*)(vector_sums_i->data);
  GMFloat* __restrict__ vector_sums_j_ = (GMFloat*)(vector_sums_j->data);
  GMFloat* __restrict__ vector_sums_k_ = (GMFloat*)(vector_sums_k->data);

  //---ISSUE: are there aliasing issues with vectors_i/j/k

  /*----------------------------------------*/
  if (Env_compute_method(env) != GM_COMPUTE_METHOD_GPU && !Env_all2all(env)) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    GMAssertAlways(GMEnv_num_section_steps(env, 1) == 1);

    /*---No off-proc all2all: compute tetrahedron of values---*/

    for (j = 0; j < numvecl; ++j) {
      for (k = j+1; k < numvecl; ++k) {
        for (i = 0; i < j; ++i) {

          const GMFloat denominator =
            vector_sums_i_[i] + vector_sums_i_[j] + vector_sums_i_[k];

          GMFloat numerator = 0;
          int f = 0;
          for (f = 0; f < numfieldl; ++f) {
            const GMFloat val1 = GMVectors_float_get(vectors_i, f, i, env);
            const GMFloat val2 = GMVectors_float_get(vectors_i, f, j, env);
            const GMFloat val3 = GMVectors_float_get(vectors_i, f, k, env);
            GMFloat min12 = val1 < val2 ? val1 : val2;
            numerator += min12;
            numerator += val1 < val3 ? val1 : val3;
            numerator += val2 < val3 ? val2 : val3;
            numerator -= min12 < val3 ? min12 : val3;
          } /*---for f---*/

          const GMFloat value =
              ((GMFloat)3) * numerator / (((GMFloat)2) * denominator);

          GMMetrics_float_set_3(metrics, i, j, k, value, env);
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

    for (j = si->j_lb; j < si->j_ub; ++j) {
      const int k_min = GMSectionInfo_k_min(si, j, env);
      for (k = k_min; k < si->k_ub; ++k) {
        const int i_max = GMSectionInfo_i_max(si, j, env);
        for (i = si->i_lb; i < i_max; ++i) {

          const GMFloat denominator =
            vector_sums_i_[i] + vector_sums_j_[j] + vector_sums_k_[k];

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

          const GMFloat value =
              ((GMFloat)3) * numerator / (((GMFloat)2) * denominator);

          GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                        value, env);
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
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
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
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);

  /*---Initializations---*/

  //const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;

  const int i_block = Env_proc_num_vector_i(env);

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  int i = 0;
  int j = 0;
  int k = 0;

  /*---TODO: address aliasing---*/
  GMFloat* __restrict__ vector_sums_i_ = (GMFloat*)(vector_sums_i->data);
  GMFloat* __restrict__ vector_sums_j_ = (GMFloat*)(vector_sums_j->data);
  GMFloat* __restrict__ vector_sums_k_ = (GMFloat*)(vector_sums_k->data);

  /*----------------------------------------*/
  if (Env_compute_method(env) == GM_COMPUTE_METHOD_REF) {
    /*----------------------------------------*/

    GMInsist(env, Env_num_proc_field(env) == 1
                      ? "num_proc_field>1 for CPU case not supported"
                      : 0);

    const int j_min = si->j_lb;
    const int j_max = si->j_ub;
    for (j = j_min; j < j_max; ++j) {
      const GMTally1 sj_1 = (si->is_part1 || !Env_all2all(env)) ?
                            (GMTally1)(vector_sums_i_[j]) :
                            (GMTally1)(vector_sums_j_[j]);
      const int k_min = GMSectionInfo_k_min(si, j, env);
      const int k_max = si->k_ub;
      for (k = k_min; k < k_max; ++k) {
        const GMTally1 sk_1 = (si->is_part1 || !Env_all2all(env)) ?
                              (GMTally1)(vector_sums_i_[k]) :
                              si->is_part2 ?
                              (GMTally1)(vector_sums_j_[k]) :
                              (GMTally1)(vector_sums_k_[k]);
        const int i_min = si->i_lb;
        const int i_max = GMSectionInfo_i_max(si, j, env);
        for (i = i_min; i < i_max; ++i) {
          const GMTally1 si_1 = (GMTally1)(vector_sums_i_[i]);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
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
            GMMetrics_float3_M_set_all2all_3(metrics, i, j, k, j_block, k_block,
                si1_sj1_sk1, env);
          } else {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
            GMMetrics_float3_M_set_3(metrics, i, j, k, si1_sj1_sk1, env);
          }
        } /*---for i---*/
      }   /*---for j---*/
    }     /*---for k---*/

    /*----------------------------------------*/
  } else if (Env_compute_method(env) == GM_COMPUTE_METHOD_CPU) {
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
    const int j_min = si->j_lb;
    const int j_max = si->j_ub;
    for (j = j_min; j < j_max; ++j) {
      const GMTally1 sj_1 = (si->is_part1 || !Env_all2all(env)) ?
                            (GMTally1)(vector_sums_i_[j]) :
                            (GMTally1)(vector_sums_j_[j]);
      const int k_min = GMSectionInfo_k_min(si, j, env);
      const int k_max = si->k_ub;
      for (k = k_min; k < k_max; ++k) {
        const GMTally1 sk_1 = (si->is_part1 || !Env_all2all(env)) ?
                              (GMTally1)(vector_sums_i_[k]) :
                              si->is_part2 ?
                              (GMTally1)(vector_sums_j_[k]) :
                              (GMTally1)(vector_sums_k_[k]);
        const int i_min = si->i_lb;
        const int i_max = GMSectionInfo_i_max(si, j, env);
        for (i = i_min; i < i_max; ++i) {
          const GMTally1 si_1 = (GMTally1)(vector_sums_i_[i]);
          const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
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
            GMMetrics_float3_M_set_all2all_3(metrics, i, j, k, j_block, k_block,
                si1_sj1_sk1, env);
          } else {
            GMMetrics_tally4x2_set_3(metrics, i, j, k, sum, env);
            GMMetrics_float3_M_set_3(metrics, i, j, k, si1_sj1_sk1, env);
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
  const GMSectionInfo* const si,
  GMVectorSums* vector_sums_i,
  GMVectorSums* vector_sums_j,
  GMVectorSums* vector_sums_k,
  GMEnv* const env) {

  const _Bool is_part3 = si->is_part3;

  /*---TODO: address aliasing---*/
  GMFloat* __restrict__ vector_sums_i_ = (GMFloat*)(vector_sums_i->data);
  GMFloat* __restrict__ vector_sums_j_ = (GMFloat*)(vector_sums_j->data);
  GMFloat* __restrict__ vector_sums_k_ = (GMFloat*)(vector_sums_k->data);

  const size_t numvecl64 = (size_t)numvecl;
  const size_t I_max64 = (size_t)I_max;

  /*--------------------*/
  /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
  /*--------------------*/

  /*----------*/
  if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI && !Env_all2all(env)) {
    /*----------*/

    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvecl * J];
      int K = 0;
      for (K = K_min; K < K_max; ++K) {
        const GMFloat min_JK = ((GMFloat*)(matM_JK_buf->h))[J + numvecl64*K];
        const GMFloat min_KIK = ((GMFloat*)(matM_KIK_buf->h))[K + numvecl64*I];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max64*K];
        const GMFloat numerator = min_IJ + min_JK + min_KIK - min_IJK;
        const int i = I;
        const int j = J;
        const int k = K;
        const GMFloat denominator =
            vector_sums_i_[i] + vector_sums_i_[j] + vector_sums_i_[k];
        const GMFloat value = ((GMFloat)1.5) * numerator / denominator;
        GMMetrics_float_set_3(metrics, i, j, k, value, env);
      } /*---for K---*/
    }   /*---for I---*/
    /*----------*/
  } else if (Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI &&
      Env_all2all(env)) {
    /*----------*/

    int I = 0;
    for (I = I_min; I < I_max; ++I) {
      const GMFloat min_IJ = ((GMFloat*)(matM_IJ_buf->h))[I + numvecl64 * J];
      int K = 0;
      for (K = K_min; K < K_max; ++K) {
        const GMFloat min_JK = ((GMFloat*)(matM_JK_buf->h))[J + numvecl64 * K];
        const GMFloat min_KIK =
            is_part3 ? ((GMFloat*)(matM_KIK_buf->h))[K + numvecl64 * I]
                     : ((GMFloat*)(matM_KIK_buf->h))[I + numvecl64 * K];
        // sum of mins vectors i, j, and k is matB(k,i)
        const GMFloat min_IJK = ((GMFloat*)(matB_buf->h))[I + I_max64 * K];
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

        const GMFloat denominator =
            vector_sums_i_[i] + vector_sums_j_[j] + vector_sums_k_[k];
        const GMFloat value = ((GMFloat)1.5) * numerator / denominator;
        GMMetrics_float_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                      value, env);
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

        const GMTally1 sj_1 = (GMTally1)(vector_sums_i_[j]);
        const GMTally1 sk_1 = (GMTally1)(vector_sums_i_[k]);
        const GMTally1 si_1 = (GMTally1)(vector_sums_i_[i]);
        const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
        GMMetrics_float3_M_set_3(metrics, i, j, k, si1_sj1_sk1, env);

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
        const GMTally1 sj_1 = si->is_part1 ?
                              (GMTally1)(vector_sums_i_[j]) :
                              (GMTally1)(vector_sums_j_[j]); 
        const GMTally1 sk_1 = si->is_part1 ?
                              (GMTally1)(vector_sums_i_[k]) :
                              si->is_part2 ?
                              (GMTally1)(vector_sums_j_[k]) :
                              (GMTally1)(vector_sums_k_[k]); 
        const GMTally1 si_1 = (GMTally1)(vector_sums_i_[i]);
        const GMFloat3 si1_sj1_sk1 = GMFloat3_encode(si_1, sj_1, sk_1);
        GMMetrics_float3_M_set_all2all_3(metrics, i, j, k, j_block, k_block,
                                         si1_sj1_sk1, env);

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

void gm_compute_numerators_3way_gpu_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_block,
    int k_block,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
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
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);

  /*---Initializations---*/

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  const int numvecl = metrics->num_vector_local;
  const int numfieldl = vectors_i->num_field_local;
  const int numpfieldl = vectors_i->num_packedval_field_local;
  const size_t metrics_buf_size = numvecl * (size_t)numvecl;
  const size_t vectors_buf_size = numvecl * (size_t)numpfieldl;

  const int i_block = Env_proc_num_vector_i(env);

  const _Bool need_2way = Env_metric_type(env) == GM_METRIC_TYPE_CZEKANOWSKI;

  GMSectionInfo si_value;
  GMSectionInfo* si = &si_value;
  GMSectionInfo_create(si, i_block, j_block, k_block, section_step,
                       metrics->num_vector_local, env);

  /*----------------------------------------*/
  /*---First get the required 2-way ij, jk, ik metrics---*/
  /*----------------------------------------*/

  /*--------------------*/
  /*---Compute i_block - j_block PROD---*/
  /*--------------------*/

  GMMirroredPointer mat_buf_tmp_value0 = gm_malloc_magma(metrics_buf_size, env);
  GMMirroredPointer mat_buf_tmp_value1 = gm_malloc_magma(metrics_buf_size, env);
  GMMirroredPointer* mat_buf_tmp[2] =
    {&mat_buf_tmp_value0, &mat_buf_tmp_value1};

  GMMirroredPointer matM_ij_buf_value = need_2way
    ? gm_malloc_magma(metrics_buf_size, env)
    : GMMirroredPointer_null();  // M = X^T PROD X
  GMMirroredPointer* const matM_ij_buf = need_2way ? &matM_ij_buf_value : NULL;

  if (need_2way) {
    GMMirroredPointer* matM_ij_buf_local =
        Env_num_proc_field(env) == 1 ? matM_ij_buf : mat_buf_tmp[0];

    gm_magma_set_matrix_zero_start(matM_ij_buf_local, numvecl, numvecl, env);

    gm_magma_gemm_start(numvecl, numvecl, numpfieldl,
                        vectors_i_buf->d, numpfieldl,
                        vectors_j_buf->d, numpfieldl,
                        matM_ij_buf_local->d, numvecl, env);
    gm_compute_wait(env);

    gm_get_matrix_start(matM_ij_buf_local, numvecl, numvecl, env);
    gm_get_matrix_wait(env);

    if (Env_num_proc_field(env) > 1) {
      GMAssertAlways(metrics_buf_size == (size_t)(int)metrics_buf_size);
      mpi_code = MPI_Allreduce(matM_ij_buf_local->h, matM_ij_buf->h,
                               metrics_buf_size, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
  }

  /*--------------------*/
  /*---Compute j_block - k_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  GMMirroredPointer matM_jk_buf_value = need_2way && !si->is_part1
    ? gm_malloc_magma(metrics_buf_size, env)
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

    if (Env_num_proc_field(env) > 1) {
      GMAssertAlways(metrics_buf_size == (size_t)(int)metrics_buf_size);
      mpi_code = MPI_Allreduce(matM_jk_buf_local->h, matM_jk_buf->h,
                               metrics_buf_size, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
  }

  /*--------------------*/
  /*---Compute k_block - i_block PROD---*/
  /*--------------------*/

  /*---Need to compute only if not identical to already computed values---*/

  /*---NOTE: for Part 3, this is indexed directly as (k,i).
       Otherwise, it is indexed through an alias as (i,k)---*/

  GMMirroredPointer matM_kik_buf_value = need_2way && si->is_part3
    ? gm_malloc_magma(metrics_buf_size, env)
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

    if (Env_num_proc_field(env) > 1) {
      GMAssertAlways(metrics_buf_size == (size_t)(int)metrics_buf_size);
      mpi_code = MPI_Allreduce(matM_kik_buf_local->h, matM_kik_buf->h,
                               metrics_buf_size, GM_MPI_FLOAT, MPI_SUM,
                               Env_mpi_comm_field(env));
      GMAssertAlways(mpi_code == MPI_SUCCESS);
    }
  } /*---is_part3---*/

  /*----------------------------------------*/
  /*---Now compute ijk piece, via an outer loop over j values---*/
  /*----------------------------------------*/

  /*---Allocate magma CPU/GPU memory for matrices V and B---*/
  // V = elementwise OP of one vector with each of the rest of the vectors.
  // For the jth iteration, the ith column of V is the elementwise OP
  //   of vectors i and j.
  // B = X^T PROD V = three way PROD.

  GMMirroredPointer matV_buf_value0 = gm_malloc_magma(vectors_buf_size, env);
  GMMirroredPointer matV_buf_value1 = gm_malloc_magma(vectors_buf_size, env);
  GMMirroredPointer* matV_buf[2] = {&matV_buf_value0, &matV_buf_value1};

  GMMirroredPointer matB_buf_value0 = gm_malloc_magma(metrics_buf_size, env);
  GMMirroredPointer matB_buf_value1 = gm_malloc_magma(metrics_buf_size, env);
  GMMirroredPointer* matB_buf[2] = {&matB_buf_value0, &matB_buf_value1};

  /*---Set up pointers to permute the access of axes for Part 3---*/
  /*---We use capitals I, J, K here to denote the PERMUTED axes---*/

  const _Bool sax0 = si->section_axis == 0;
  const _Bool sax1 = si->section_axis == 1;
  const _Bool sax2 = si->section_axis == 2;

  const _Bool is_ijk = !si->is_part3 ? GM_BOOL_TRUE : sax1;
  const _Bool is_kij = !si->is_part3 ? GM_BOOL_FALSE : sax0;
  const _Bool is_jki = !si->is_part3 ? GM_BOOL_FALSE : sax2;

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

  /*--------------------*/
  /*---Collapsed loops over J and over 2-way steps---*/
  /*--------------------*/

  const int J_min = si->J_lb;
  const int J_max = si->J_ub;
  const int J_count = J_max - J_min;

  const int num_step_2way = Env_metric_type(env) == GM_METRIC_TYPE_CCC ?  3 : 1;
  const int num_step = J_count * num_step_2way;
  const int extra_step = 1;
  int step_num = 0;

  MPI_Request mpi_requests[2];
//  _Bool matV_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matV_d_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_d_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};
//  _Bool matB_ptr_h_updating[2] = {GM_BOOL_FALSE, GM_BOOL_FALSE};

  typedef struct {
    int step_num;
    int step_2way;
    int J;
    int I_min;
    int I_max;
    int K_min;
    int K_max;
    _Bool empty;
    _Bool is_compute_step;
    _Bool do_compute;
    int index_01;
    GMMirroredPointer* matB_buf_ptr;
  } LoopVars;

  LoopVars vars;
  LoopVars vars_prev;
  LoopVars vars_prevprev;
  LoopVars vars_next;
  vars.do_compute = GM_BOOL_FALSE;
  vars_prev.do_compute = GM_BOOL_FALSE;
  vars_prevprev.do_compute = GM_BOOL_FALSE;
  vars_next.do_compute = GM_BOOL_FALSE;

  for (step_num = 0-extra_step; step_num < num_step+extra_step*2; ++step_num) {

    vars_prevprev = vars_prev;
    vars_prev = vars;
    vars = vars_next;

    vars_next.step_num = step_num + 1;
    vars_next.step_2way = gm_mod_i(vars_next.step_num, num_step_2way);
    vars_next.J = J_min + gm_floor_i(vars_next.step_num, num_step_2way);
    vars_next.I_min = 0;
    vars_next.I_max = si->is_part1 ? vars_next.J : numvecl;
    vars_next.K_min = si->is_part3 ? 0 : vars_next.J + 1;
    vars_next.K_max = numvecl;
    vars_next.empty = vars_next.I_min >= vars_next.I_max ||
                      vars_next.K_min >= vars_next.K_max;
    vars_next.is_compute_step = vars_next.step_num >= 0 &&
                                vars_next.step_num < num_step;
    vars_next.do_compute = vars_next.is_compute_step && ! vars_next.empty;
    vars_next.index_01 = gm_mod_i(vars_next.step_num, 2);
    vars_next.matB_buf_ptr = Env_num_proc_field(env) == 1 ?
      matB_buf[vars_next.index_01] : mat_buf_tmp[vars_next.index_01];

//    vars_next.matV_h_updating = &matV_h_updating[vars_next.index_01];
//    vars_next.matV_d_updating = &matV_d_updating[vars_next.index_01];
//    vars_next.matB_h_updating = &matB_h_updating[vars_next.index_01];
//    vars_next.matB_d_updating = &matB_d_updating[vars_next.index_01];
//    vars_next.matB_ptr_h_updating = &matB_ptr_h_updating[vars_next.index_01];

    //==========

    if (vars_next.do_compute) {
      /*---Populate leading columns of matV---*/
//      GMAssertAlways(!*vars_next.matV_h_updating);
//      *vars_next.matV_h_updating = GM_BOOL_TRUE;
      gm_compute_numerators_3way_gpu_form_matV(
          vectors_I_buf, vectors_J_buf, matV_buf[vars_next.index_01],
          numfieldl, numpfieldl, vars_next.J, vars_next.step_2way,
          vars_next.I_min, vars_next.I_max, env);
      //matV_buf[vars_next.index_01]->h is now overwritten.
//      *vars_next.matV_h_updating = GM_BOOL_FALSE;
    }

    //==========

    if (vars.do_compute) {
      /*---Send matrix matV to GPU - WAIT---*/
      gm_set_matrix_wait(env);
      //matV_buf[vars.index_01]->h is now no longer needed.
      //matV_buf[vars.index_01]->d is now available.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Perform pseudo mat X matt matB = matV^T PROD X - WAIT---*/
      gm_compute_wait(env);
      //matV_buf[vars_prev.index_01]->d is now no longer needed.
      //matB_buf_ptr[vars_prev.index_01]->d is now available.
    }

    //==========

    if (vars_next.do_compute) {
      /*---Send matrix matV to GPU - START---*/
      gm_set_matrix_start(matV_buf[vars_next.index_01], numpfieldl,
                          vars_next.I_max, env);
      //matV_buf[vars_next.index_01]->d (= matV_buf[vars_prev.index_01]->d)
      //is now being overwritten.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - START---*/
      gm_get_matrix_start(vars_prev.matB_buf_ptr, vars_prev.I_max,
                          numvecl, env);
      //matB_buf_ptr[vars_prev.index_01]->h is now being overwritten.
    }

    //==========

    if (vars.do_compute) {
      /*---Initialize result matrix to zero (apparently magma requires)---*/
      gm_magma_set_matrix_zero_start(vars.matB_buf_ptr, numvecl,
                                     vars.I_max, env);
      /*---Perform pseudo mat X mat matB = matV^T PROD X - START---*/
      gm_magma_gemm_start(vars.I_max, numvecl, numpfieldl,
                          matV_buf[vars.index_01]->d, numpfieldl,
                          vectors_K_buf->d, numpfieldl,
                          vars.matB_buf_ptr->d, vars.I_max, env);
      //matB_buf_ptr[vars.index_01]->d is now being overwritten.
    }

    //==========

    if (vars_prev.do_compute) {
      /*---Copy result matrix matB from GPU - WAIT---*/
      gm_get_matrix_wait(env);
      //matB_buf_ptr[vars_prev.index_01]->h is now available.
      //matB_buf_ptr[vars_prev.index_01]->d is now no longer needed.
    }

    //==========

    if (vars_prevprev.do_compute) {
      if (Env_num_proc_field(env) > 1) {
        MPI_Status mpi_status;
        mpi_code = MPI_Wait(&(mpi_requests[vars_prevprev.index_01]),
                            &mpi_status);
        GMAssertAlways(mpi_code == MPI_SUCCESS);
        //matB_buf[vars.prevprev.index_01]->h is now available.
      }
    }

    //==========

    if (vars_prev.do_compute) {
      if (Env_num_proc_field(env) > 1) {
        //TODO: fix this properly.
        const int factor = Env_metric_type(env) == GM_METRIC_TYPE_CCC ? 2 : 1;
        const size_t metrics_buf_size_this = metrics_buf_size * factor;
        GMAssertAlways(metrics_buf_size_this ==
                       (size_t)(int)metrics_buf_size_this);
        mpi_code = MPI_Iallreduce(vars_prev.matB_buf_ptr->h,
                                  matB_buf[vars_prev.index_01]->h,
                                  metrics_buf_size_this,
                                  GM_MPI_FLOAT, MPI_SUM,
                                  Env_mpi_comm_field(env),
                                  &(mpi_requests[vars_prev.index_01]));
        GMAssertAlways(mpi_code == MPI_SUCCESS);
        //matB_buf[vars.prev.index_01]->h is now being overwritten.
      }
    }

    //==========

    //---NOTE: matB_buf[vars_prevprev.index_01]->d is being updated now
    //---but matB_buf[vars_prevprev.index_01]->h is ok.

    if (vars_prevprev.do_compute) {
      /*---Compute numerators using ijk piece and (if needed) 2-way pieces---*/
      gm_compute_numerators_3way_gpu_form_metrics(
          matM_IJ_buf, matM_JK_buf, matM_KIK_buf,
          matB_buf[vars_prevprev.index_01],
          metrics, numvecl,
          vars_prevprev.J,
          vars_prevprev.step_2way,
          vars_prevprev.I_min,
          vars_prevprev.I_max,
          vars_prevprev.K_min,
          vars_prevprev.K_max,
          j_block, k_block, sax0, sax1, si,
          vector_sums_i, vector_sums_j, vector_sums_k,
          env);
      //matB_buf[vars.prevprev.index_01]->h (=matB_buf[vars.index_01]->h)
      //now no longer needed.
    }

    //==========

  } /*---for step_num---*/

  /*--------------------*/
  /*---Free memory---*/
  /*--------------------*/

  if (matM_ij_buf != NULL) {
    gm_free_magma(matM_ij_buf, env);
  }

  if (need_2way) {
    if (!si->is_part1) {
      gm_free_magma(&matM_jk_buf_value, env);
    }
    if (si->is_part3) {
      gm_free_magma(&matM_kik_buf_value, env);
    }
  }
  int i = 0;
  for (i=0; i<2; ++i) {
    gm_free_magma(mat_buf_tmp[i], env);
    gm_free_magma(matV_buf[i], env);
    gm_free_magma(matB_buf[i], env);
  }
  GMSectionInfo_destroy(si, env);
}

/*===========================================================================*/
/*---Start calculation of numerators, 3-way generic---*/

void gm_compute_numerators_3way_start(
    GMVectors* vectors_i,
    GMVectors* vectors_j,
    GMVectors* vectors_k,
    GMMetrics* metrics,
    GMMirroredPointer* vectors_i_buf,
    GMMirroredPointer* vectors_j_buf,
    GMMirroredPointer* vectors_k_buf,
    int j_block,
    int k_block,
    GMVectorSums* vector_sums_i,
    GMVectorSums* vector_sums_j,
    GMVectorSums* vector_sums_k,
    int section_step,
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
  GMAssertAlways(vector_sums_i != NULL);
  GMAssertAlways(vector_sums_j != NULL);
  GMAssertAlways(vector_sums_k != NULL);

  /*----------------------------------------*/
  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    /*----------------------------------------*/

    gm_compute_numerators_3way_gpu_start(vectors_i, vectors_j, vectors_k,
                                         metrics, vectors_i_buf, vectors_j_buf,
                                         vectors_k_buf, j_block, k_block,
                                         vector_sums_i, vector_sums_j,
                                         vector_sums_k,
                                         section_step, env);
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
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
      } break;
      /*----------------------------------------*/
      case GM_METRIC_TYPE_CCC: {
        /*----------------------------------------*/
        gm_compute_ccc_numerators_3way_nongpu_start(
            vectors_i, vectors_j, vectors_k, metrics, vectors_i_buf,
            vectors_j_buf, vectors_k_buf, j_block, k_block,
            vector_sums_i, vector_sums_j, vector_sums_k,
            section_step, env);
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

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
