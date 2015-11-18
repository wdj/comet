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

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Types---*/

GMMirroredPointer GMMirroredPointer_null(void) {
  GMMirroredPointer p;
  p.h = NULL;
  p.d = NULL;
  return p;
}

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
  GMAssert(env != NULL);

  /*---Initialize MPI info---*/

  env->mpi_comm_ = MPI_COMM_WORLD;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Comm_size(MPI_COMM_WORLD, &env->num_proc_world_);
  GMAssert(mpi_code == MPI_SUCCESS);

  Env_set_num_proc(env, env->num_proc_world_);

  /*---Set default values---*/
  env->metric_type = GM_METRIC_TYPE_CZEKANOWSKI;
  env->num_way = 2;
  env->all2all = GM_BOOL_FALSE;
  env->are_cuda_streams_initialized_ = GM_BOOL_FALSE;
  Env_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
}

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create_from_args(GMEnv* env, int argc, char** argv) {
  GMAssert(env != NULL);

  /*---First initialize with standard constructor---*/
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
        Env_set_compute_method(env, GM_COMPUTE_METHOD_CPU);
      } else if (strcmp(argv[i], "GPU") == 0) {
        Env_set_compute_method(env, GM_COMPUTE_METHOD_GPU);
      } else if (strcmp(argv[i], "REF") == 0) {
        Env_set_compute_method(env, GM_COMPUTE_METHOD_REF);
      } else {
        GMInsist(env,
                 GM_BOOL_FALSE ? "Invalid setting for compute_method." : 0);
      }
    } else if (strcmp(argv[i], "--num_proc") == 0) {
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_proc." : 0);
      const int num_proc = atoi(argv[i]);
      Env_set_num_proc(env, num_proc);
    } /*---if/else---*/
  }   /*---for i---*/
}

/*===========================================================================*/
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* env) {
  GMAssert(env != NULL);

  if (!env->are_cuda_streams_initialized_) {
    cudaStreamCreate(&env->stream_compute);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamCreate(&env->stream_togpu);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamCreate(&env->stream_fromgpu);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    env->are_cuda_streams_initialized_ = GM_BOOL_TRUE;
  }
}

/*---------------------------------------------------------------------------*/

void GMEnv_terminate_streams(GMEnv* env) {
  GMAssert(env != NULL);

  if (env->are_cuda_streams_initialized_) {
    cudaStreamDestroy(env->stream_compute);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamDestroy(env->stream_togpu);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    cudaStreamDestroy(env->stream_fromgpu);
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
    env->are_cuda_streams_initialized_ = GM_BOOL_FALSE;
  }
}

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* env) {
  GMAssert(env != NULL);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /*---Destroy any nontrivial communicator---*/
  if (env->mpi_comm_ != MPI_COMM_WORLD) {
    mpi_code = MPI_Comm_free(&(env->mpi_comm_));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  GMEnv_terminate_streams(env);

  *env = GMEnv_null();
}

/*===========================================================================*/
/*---Accessors---*/

void Env_set_num_proc(GMEnv* env, int num_proc) {
  GMAssert(env != NULL);
  GMAssert(num_proc >= 0);
  GMAssert(num_proc <= env->num_proc_world_);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /*---Destroy old communictor if necessary---*/

  if (env->mpi_comm_ != MPI_COMM_WORLD) {
    mpi_code = MPI_Comm_free(&(env->mpi_comm_));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  /*---Make new communicator---*/

  mpi_code = MPI_Comm_rank(MPI_COMM_WORLD, &env->proc_num_);
  GMAssert(mpi_code == MPI_SUCCESS);
  env->is_proc_active_ = env->proc_num_ < num_proc;
  mpi_code = MPI_Comm_split(MPI_COMM_WORLD, env->is_proc_active_,
                                         env->proc_num_, &env->mpi_comm_);
  GMAssert(mpi_code == MPI_SUCCESS);

  /*---Get info---*/

  mpi_code = MPI_Comm_size(env->mpi_comm_, &env->num_proc_);
  GMAssert(mpi_code == MPI_SUCCESS);
  mpi_code = MPI_Comm_rank(env->mpi_comm_, &env->proc_num_);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

void Env_set_compute_method(GMEnv* env, int compute_method) {
  GMAssert(env != NULL);
  GMAssert(compute_method >= 0);
  GMAssert(compute_method < GM_NUM_COMPUTE_METHOD);

  if (compute_method == GM_COMPUTE_METHOD_GPU) {
    GMEnv_initialize_streams(env);
  } else {
    GMEnv_terminate_streams(env);
  }

  env->compute_method_ = compute_method;
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
  GMAssert(env != NULL);

  /*
  cudaThreadSynchronize();
  */

  if (Env_compute_method(env) == GM_COMPUTE_METHOD_GPU) {
    cudaDeviceSynchronize();
    GMAssert(GMEnv_cuda_last_call_succeeded(env));
  }

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  mpi_code = MPI_Barrier(Env_mpi_comm(env));
  GMAssert(mpi_code == MPI_SUCCESS);
  return GMEnv_get_time(env);
}

/*===========================================================================*/
/*---Assertions---*/

void gm_insist(GMEnv* env,
               const char* condition_string,
               const char* file,
               int line) {
  if (Env_proc_num(env) == 0) {
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

_Bool GMEnv_cuda_last_call_succeeded(GMEnv* env) {
  _Bool result = GM_BOOL_TRUE;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if (error != cudaSuccess) {
    result = GM_BOOL_FALSE;
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
}

/*---------------------------------------------------------------------------*/

GMFloat* GMFloat_malloc(size_t n) {
  GMFloat* result = (GMFloat*)malloc(n * sizeof(GMFloat));
  GMAssert(result);
  return result;
}

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
