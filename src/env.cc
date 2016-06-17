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

#ifdef TESTING
#include "gtest/gtest.h"
#endif

#include "mpi.h"
#include "cuda.h"

#include "env.h"

/*---------------------------------------------------------------------------*/

static void gm_test_wrapper() {
#ifdef TESTING
  ASSERT_TRUE(0);
#endif
}

/*---------------------------------------------------------------------------*/

#ifdef __cplusplus
extern "C" {
#endif

/*===========================================================================*/
/*---Assertions---*/

void gm_assert(const char* condition_string, const char* file, int line) {
  fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Assertion error",
          condition_string, file, line);
  gm_test_wrapper();
  exit(EXIT_FAILURE);
}

/*---------------------------------------------------------------------------*/

void gm_insist(const void* env,
               const char* condition_string,
               const char* file,
               int line) {
  if (Env_proc_num((const GMEnv*)env) == 0) {
    fprintf(stderr, "%s: \"%s\". At file %s, line %i.\n", "Interface error",
            condition_string, file, line);
  }
  gm_test_wrapper();
  exit(EXIT_FAILURE);
}

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

  GMStaticAssert(sizeof(GMBits) == 8);
  GMStaticAssert(sizeof(GMULInt) == 8);

  *env = GMEnv_null();

  /*---Initialize MPI info---*/

  env->mpi_comm_ = MPI_COMM_WORLD;
  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/
  mpi_code = MPI_Comm_size(MPI_COMM_WORLD, &env->num_proc_world_);
  GMAssert(mpi_code == MPI_SUCCESS);

  Env_set_num_proc(env, env->num_proc_world_, 1, 1);

  /*---Set default values---*/
  env->metric_type_ = GM_METRIC_TYPE_CZEKANOWSKI;
  env->num_way_ = GM_NUM_WAY_2;
  env->all2all_ = GM_BOOL_FALSE;
  env->are_cuda_streams_initialized_ = GM_BOOL_FALSE;
  Env_set_compute_method(env, GM_COMPUTE_METHOD_GPU);

  env->time = 0;
  env->ops = 0;
}

/*===========================================================================*/
/*---Initialize environment---*/

void GMEnv_create_from_args(GMEnv* env, int argc, const char** argv) {
  GMAssert(env != NULL);

  /*---First initialize with standard constructor---*/
  GMEnv_create(env);

  /*---Modify based on user options---*/
  int i = 0;
  for (i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "--metric_type") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for metric_type." : 0);
      if (strcmp(argv[i], "sorenson") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_SORENSON;
      } else if (strcmp(argv[i], "czekanowski") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CZEKANOWSKI;
      } else if (strcmp(argv[i], "ccc") == 0) {
        env->metric_type_ = GM_METRIC_TYPE_CCC;
      } else {
        GMInsist(env, GM_BOOL_FALSE ? "Invalid setting for metric_type." : 0);
      }
      /*----------*/
    } else if (strcmp(argv[i], "--num_way") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_way." : 0);
      env->num_way_ = atoi(argv[i]);
      GMInsist(env, env->num_way_ == GM_NUM_WAY_2 ||
                            env->num_way_ == GM_NUM_WAY_3
                        ? "Invalid setting for num_way."
                        : 0);
      Env_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                       env->num_proc_field_);
      /*----------*/
    } else if (strcmp(argv[i], "--all2all") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for all2all." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        env->all2all_ = GM_BOOL_TRUE;
      } else if (strcmp(argv[i], "no") == 0) {
        env->all2all_ = GM_BOOL_FALSE;
      } else {
        GMInsist(env, GM_BOOL_FALSE ? "Invalid setting for all2all." : 0);
      }
      /*----------*/
    } else if (strcmp(argv[i], "--compute_method") == 0) {
      /*----------*/
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
      /*----------*/
    } else if (strcmp(argv[i], "--num_proc_vector") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_proc_vector." : 0);
      const int num_proc_vector_i = atoi(argv[i]);
      Env_set_num_proc(env, num_proc_vector_i, env->num_proc_repl_,
                       env->num_proc_field_);
      /*----------*/
    } else if (strcmp(argv[i], "--num_proc_field") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_proc_field." : 0);
      const int num_proc_field = atoi(argv[i]);
      Env_set_num_proc(env, env->num_proc_vector_i_, env->num_proc_repl_,
                       num_proc_field);
      /*----------*/
    } else if (strcmp(argv[i], "--num_proc_repl") == 0) {
      /*----------*/
      ++i;
      GMInsist(env, i < argc ? "Missing value for num_proc_repl." : 0);
      const int num_proc_repl = atoi(argv[i]);
      Env_set_num_proc(env, env->num_proc_vector_i_, num_proc_repl,
                       env->num_proc_field_);
      /*----------*/
    } /*---if/else---*/
  }   /*---for i---*/
}

/*===========================================================================*/
/*---Manage cuda streams---*/

void GMEnv_initialize_streams(GMEnv* env) {
  GMAssert(env != NULL);

  /*---NOTE: this is used for lazy initialization---*/

  if (!env->are_cuda_streams_initialized_) {
    if (env->compute_method_ == GM_COMPUTE_METHOD_GPU) {
      cudaStreamCreate(&env->stream_compute_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));

      cudaStreamCreate(&env->stream_togpu_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));

      cudaStreamCreate(&env->stream_fromgpu_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));
    }
  }
  env->are_cuda_streams_initialized_ = GM_BOOL_TRUE;
}

/*---------------------------------------------------------------------------*/

void GMEnv_terminate_streams(GMEnv* env) {
  GMAssert(env != NULL);

  if (env->are_cuda_streams_initialized_) {
    if (env->compute_method_ == GM_COMPUTE_METHOD_GPU) {
      cudaStreamDestroy(env->stream_compute_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));

      cudaStreamDestroy(env->stream_togpu_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));

      cudaStreamDestroy(env->stream_fromgpu_);
      GMAssert(GMEnv_cuda_last_call_succeeded(env));
    }
  }
  env->are_cuda_streams_initialized_ = GM_BOOL_FALSE;
}

/*===========================================================================*/
/*---Finalize environment---*/

void GMEnv_destroy(GMEnv* env) {
  GMAssert(env != NULL);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /*---Destroy any nontrivial communicators---*/
  if (env->mpi_comm_ != MPI_COMM_WORLD) {
    mpi_code = MPI_Comm_free(&(env->mpi_comm_));
    GMAssert(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_free(&(env->mpi_comm_vector_));
    GMAssert(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_free(&(env->mpi_comm_field_));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  GMEnv_terminate_streams(env);

  *env = GMEnv_null();
}

/*===========================================================================*/
/*---Accessors---*/

void Env_set_compute_method(GMEnv* env, int compute_method) {
  GMAssert(env != NULL);
  GMAssert(compute_method >= 0);
  GMAssert(compute_method < GM_NUM_COMPUTE_METHOD);

  env->compute_method_ = compute_method;
}

/*---------------------------------------------------------------------------*/

int Env_data_type_vectors(const GMEnv* env) {
  GMAssert(env != NULL);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_SORENSON:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    case GM_METRIC_TYPE_CZEKANOWSKI:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return GM_DATA_TYPE_BITS2;
  }
  GMAssert(GM_BOOL_FALSE ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

int Env_data_type_metrics(const GMEnv* env) {
  GMAssert(env != NULL);

  switch (env->metric_type_) {
    case GM_METRIC_TYPE_SORENSON:
      GMInsist(env, GM_BOOL_FALSE ? "Unimplemented." : 0);
    case GM_METRIC_TYPE_CZEKANOWSKI:
      return GM_DATA_TYPE_FLOAT;
    case GM_METRIC_TYPE_CCC:
      return env->num_way_ == GM_NUM_WAY_2 ? GM_DATA_TYPE_TALLY2X2
                                           : GM_DATA_TYPE_TALLY4X2;
  }
  GMAssert(GM_BOOL_FALSE ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

void Env_set_num_proc(GMEnv* env, int num_proc_vector_i, int num_proc_repl,
                      int num_proc_field) {
  GMAssert(env != NULL);
  GMAssertAlways(num_proc_vector_i > 0);
  GMAssertAlways(num_proc_repl > 0);
  GMAssertAlways(num_proc_field >= 0);

  int mpi_code = 0;
  mpi_code = mpi_code * 1; /*---Avoid unused variable warning---*/

  /*---Destroy old communicators if necessary---*/

  if (env->mpi_comm_ != MPI_COMM_WORLD) {
    mpi_code = MPI_Comm_free(&(env->mpi_comm_));
    GMAssert(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_free(&(env->mpi_comm_vector_));
    GMAssert(mpi_code == MPI_SUCCESS);
    mpi_code = MPI_Comm_free(&(env->mpi_comm_field_));
    GMAssert(mpi_code == MPI_SUCCESS);
  }

  /*---Set proc counts---*/

  env->num_proc_vector_i_ = num_proc_vector_i;
  env->num_proc_repl_ = num_proc_repl;
  env->num_proc_field_ = num_proc_field;

  env->num_proc_vector_ = env->num_proc_vector_i_ * env->num_proc_repl_;

  env->num_proc_ = env->num_proc_vector_ * num_proc_field;
  GMAssertAlways(env->num_proc_ <= env->num_proc_world_);

  /*---Set proc nums---*/

  mpi_code = MPI_Comm_rank(MPI_COMM_WORLD, &env->proc_num_);
  GMAssert(mpi_code == MPI_SUCCESS);

  int itmp = env->proc_num_;

  env->proc_num_repl_ = itmp % env->num_proc_repl_;
  itmp /= env->num_proc_repl_;

  env->proc_num_vector_i_ = itmp % env->num_proc_vector_i_;
  itmp /= env->num_proc_vector_i_;

  env->proc_num_field_ = itmp % env->num_proc_field_;
  env->proc_num_vector_ = env->proc_num_ % env->num_proc_vector_;

  /*---Make new communicators---*/

  env->is_proc_active_ = env->proc_num_ < env->num_proc_;
  mpi_code = MPI_Comm_split(MPI_COMM_WORLD, env->is_proc_active_,
                            env->proc_num_, &env->mpi_comm_);
  GMAssert(mpi_code == MPI_SUCCESS);
  mpi_code = MPI_Comm_split(MPI_COMM_WORLD,
      env->is_proc_active_ ? env->proc_num_field_ : env->num_proc_,
      env->proc_num_, &env->mpi_comm_vector_);
  GMAssert(mpi_code == MPI_SUCCESS);
  mpi_code = MPI_Comm_split(MPI_COMM_WORLD,
      env->is_proc_active_ ? env->proc_num_vector_ : env->num_proc_,
      env->proc_num_, &env->mpi_comm_field_);
  GMAssert(mpi_code == MPI_SUCCESS);
}

/*---------------------------------------------------------------------------*/

cudaStream_t Env_stream_compute(GMEnv* env) {
  GMAssert(env != NULL);
  // GMAssert(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_compute_;
}

/*---------------------------------------------------------------------------*/

cudaStream_t Env_stream_togpu(GMEnv* env) {
  GMAssert(env != NULL);
  // GMAssert(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_togpu_;
}

/*---------------------------------------------------------------------------*/

cudaStream_t Env_stream_fromgpu(GMEnv* env) {
  GMAssert(env != NULL);
  // GMAssert(env->are_cuda_streams_initialized_);

  GMEnv_initialize_streams(env);

  return env->stream_fromgpu_;
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
/*---Misc.---*/

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
  GMAssertAlways(result);
  //TODO: in debug case, fill with NANs
  return result;
}

//TODO: GMInt_malloc ...

/*===========================================================================*/

#ifdef __cplusplus
} /*---extern "C"---*/
#endif

/*---------------------------------------------------------------------------*/
