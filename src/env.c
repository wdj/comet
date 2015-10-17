/*---------------------------------------------------------------------------*/
/*!
 * \file   env.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "mpi.h"
#include "cuda.h"

#include "env.h"

/*===========================================================================*/
/*---Null object---*/

GMEnv GMEnv_null() {
  GMEnv result;
  memset((void*)&result, 0, sizeof(GMEnv));
  return result;
}

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create(GMEnv* env) {
  env->mpi_comm = MPI_COMM_WORLD;

  int mpi_code = MPI_Comm_rank(env->mpi_comm, &(env->proc_num));
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  GMAssert(mpi_code == MPI_SUCCESS);
  mpi_code = MPI_Comm_size(env->mpi_comm, &(env->num_proc));
  GMAssert(mpi_code == MPI_SUCCESS);

  /*---Set default values---*/
  env->metric_type = GM_METRIC_TYPE_CZEKANOWSKI;
  env->num_way = 2;
  env->all2all = GM_BOOL_FALSE;
  env->compute_method = GM_COMPUTE_METHOD_GPU;

  /*---Prepare for lazy intialization of cuda streams---*/
  env->are_cuda_streams_initialized = GM_BOOL_FALSE;
  env->stream_compute = 0;
  env->stream_vectors = 0;
  env->stream_metrics = 0;
}

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create_from_args(GMEnv* env, int argc, char** argv) {
  /*---Initialize with standard constructor---*/
  GMEnv_create(env);

  /*---Modify based on user options---*/
  int i = 0;
  for (i = 0; i < argc; ++i) {
    if (strcmp(argv[i], "--metric_type") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for metric_type." : 0);
      if (strcmp(argv[i], "sorenson") == 0) {
        env->metric_type = GM_METRIC_TYPE_SORENSON;
      } else if (strcmp(argv[i], "czekanowski") == 0) {
        env->metric_type = GM_METRIC_TYPE_CZEKANOWSKI;
      } else if (strcmp(argv[i], "ccc") == 0) {
        env->metric_type = GM_METRIC_TYPE_CCC;
      } else {
        GMInsist(env, GM_BOOL_FALSE ? "Invalid setting for metric_type." : 0);
      }
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_way." : 0);
      env->num_way = atoi(argv[i]);
      GMInsist(env, env->num_way == 2 || env->num_way == 3
                        ? "Invalid setting for num_way."
                        : 0);

    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for all2all." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        env->all2all = GM_BOOL_TRUE;
      } else if (strcmp(argv[i], "no") == 0) {
        env->all2all = GM_BOOL_FALSE;
      } else {
        GMInsist(env, GM_BOOL_FALSE ? "Invalid setting for all2all." : 0);
      }

    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for compute_method." : 0);
      if (strcmp(argv[i], "CPU") == 0) {
        env->compute_method = GM_COMPUTE_METHOD_CPU;
      } else if (strcmp(argv[i], "GPU") == 0) {
        env->compute_method = GM_COMPUTE_METHOD_GPU;
      } else if (strcmp(argv[i], "REF") == 0) {
        env->compute_method = GM_COMPUTE_METHOD_REF;
      } else {
        GMInsist(env,
                 GM_BOOL_FALSE ? "Invalid setting for compute_method." : 0);
      }
    } /*---if/else---*/
  }   /*---for i---*/
}

/*===========================================================================*/
/*---Initialize cuda streams---*/

void GMEnv_initialize_streams(GMEnv* env) {
  GMAssert(env != NULL );

  if ( ! env->are_cuda_streams_initialized ) {
    env->stream_compute = 0;
    cudaStreamCreate( & env->stream_vectors );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamCreate( & env->stream_metrics );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    env->are_cuda_streams_initialized = GM_BOOL_TRUE;
  }
}

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* env) {
  GMAssert(env != NULL );

  /*---Make sure no communicator to destroy---*/
  GMAssert(env->mpi_comm == MPI_COMM_WORLD);

  if ( env->compute_method == GM_COMPUTE_METHOD_GPU &&
       env->are_cuda_streams_initialized ) {
    cudaStreamDestroy( env->stream_vectors );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamDestroy( env->stream_metrics );
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }

  *env = GMEnv_null();
}

/*===========================================================================*/
/*---Timer functions---*/

double GMEnv_get_time(GMEnv* env) {
  GMAssert(env);
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);
  return result;
}

/*---------------------------------------------------------------------------*/

double GMEnv_get_synced_time(GMEnv* env) {
  GMAssert(env);

  /*
  cudaThreadSynchronize();
  */

  if (env->compute_method == GM_COMPUTE_METHOD_GPU) {
    cudaDeviceSynchronize();
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }

  int mpi_code = MPI_Barrier(env->mpi_comm);
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  GMAssert(mpi_code == MPI_SUCCESS);
  return GMEnv_get_time(env);
}

/*===========================================================================*/
/*---Assertions---*/

void gm_insist(GMEnv* env,
               const char* condition_string,
               const char* file,
               int line) {
  if (env->proc_num == 0) {
    fprintf(stderr, "GM insist error: \"%s\". At file %s, line %i.\n",
            condition_string, file, line);
  }
  exit(EXIT_FAILURE);
}

/*===========================================================================*/
/*---Misc.---*/

int gm_data_type_from_metric_type(int metric_type, GMEnv* env) {
  switch (metric_type) {
    case GM_METRIC_TYPE_SORENSON:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    case GM_METRIC_TYPE_CZEKANOWSKI:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
  }
  GMAssert(GM_BOOL_FALSE ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

GMBool GMEnv_cuda_last_call_succeeded(GMEnv* env) {
  GMBool result = GM_BOOL_TRUE;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if (error != cudaSuccess) {
    result = GM_BOOL_FALSE;
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
}

/*---------------------------------------------------------------------------*/

GMFloat* GMFloat_malloc( size_t n ) {
  GMFloat* result = (GMFloat*)malloc( n * sizeof(GMFloat) );
  GMAssert(result);
  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
