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
#include "cuda_runtime_api.h"

#include "env.h"

/*===========================================================================*/
/*---Null object---*/

Env Env_null() {
  Env result;
  memset((void*)&result, 0, sizeof(Env));
  return result;
}

/*===========================================================================*/
/*---Initialize environment---*/

void Env_create(Env* env) {
  env->mpi_comm = MPI_COMM_WORLD;

  int mpi_code = MPI_Comm_rank(env->mpi_comm, &(env->proc_num));
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  Assert(mpi_code == MPI_SUCCESS);
  mpi_code = MPI_Comm_size(env->mpi_comm, &(env->num_proc));
  Assert(mpi_code == MPI_SUCCESS);

  /*---Set default values---*/
  env->metric_type = METRIC_TYPE_CZEKANOWSKI;
  env->num_way = 2;
  env->all2all = Bool_false;
  env->compute_method = COMPUTE_METHOD_GPU;
}

/*===========================================================================*/
/*---Initialize environment---*/

void Env_create_from_args(Env* env, int argc, char** argv) {
  /*---Initialize with standard constructor---*/
  Env_create(env);

  /*---Modify based on user options---*/
  int i = 0;
  for (i = 0; i < argc; ++i) {
    if (strcmp(argv[i], "--metric_type") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for metric_type." : 0);
      if (strcmp(argv[i], "sorenson") == 0) {
        env->metric_type = METRIC_TYPE_SORENSON;
      } else if (strcmp(argv[i], "czekanowski") == 0) {
        env->metric_type = METRIC_TYPE_CZEKANOWSKI;
      } else if (strcmp(argv[i], "ccc") == 0) {
        env->metric_type = METRIC_TYPE_CCC;
      } else {
        Insist(env, Bool_false ? "Invalid setting for metric_type." : 0);
      }
    } else if (strcmp(argv[i], "--num_way") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for num_way." : 0);
      env->num_way = atoi(argv[i]);
      Insist(env, env->num_way == 2 || env->num_way == 3
                      ? "Invalid setting for num_way."
                      : 0);

    } else if (strcmp(argv[i], "--all2all") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for all2all." : 0);
      if (strcmp(argv[i], "yes") == 0) {
        env->all2all = Bool_true;
      } else if (strcmp(argv[i], "no") == 0) {
        env->all2all = Bool_false;
      } else {
        Insist(env, Bool_false ? "Invalid setting for all2all." : 0);
      }

    } else if (strcmp(argv[i], "--compute_method") == 0) {
      ++i;
      Insist(env, i < argc ? "Missing value for compute_method." : 0);
      if (strcmp(argv[i], "CPU") == 0) {
        env->compute_method = COMPUTE_METHOD_CPU;
      } else if (strcmp(argv[i], "GPU") == 0) {
        env->compute_method = COMPUTE_METHOD_GPU;
      } else if (strcmp(argv[i], "REF") == 0) {
        env->compute_method = COMPUTE_METHOD_REFERENCE;
      } else {
        Insist(env, Bool_false ? "Invalid setting for compute_method." : 0);
      }
    } /*---if/else---*/
  }   /*---for i---*/
}

/*===========================================================================*/
/*---Finalize environment---*/

void Env_destroy(Env* env) {
  /*---Make sure no communicator to destroy---*/
  Assert(env->mpi_comm == MPI_COMM_WORLD);
  *env = Env_null();
}

/*===========================================================================*/
/*---Timer functions---*/

double Env_get_time(Env* env) {
  Assert(env);
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double result = ((double)tv.tv_sec + (double)tv.tv_usec * 1.e-6);
  return result;
}

/*---------------------------------------------------------------------------*/

double Env_get_synced_time(Env* env) {
  Assert(env);

  /*
  cudaThreadSynchronize();
  */
  cudaDeviceSynchronize();
  Assert(Env_cuda_last_call_succeeded(env));

  int mpi_code = MPI_Barrier(env->mpi_comm);
  if (mpi_code) {
  } /*---Avoid unused variable warning---*/
  Assert(mpi_code == MPI_SUCCESS);
  return Env_get_time(env);
}

/*===========================================================================*/
/*---Assertions---*/

void insist_(Env* env,
             const char* condition_string,
             const char* file,
             int line) {
  if (env->proc_num == 0) {
    fprintf(stderr, "Insist error: \"%s\". At file %s, line %i.\n",
            condition_string, file, line);
  }
  exit(EXIT_FAILURE);
}

/*===========================================================================*/
/*---Misc.---*/

int data_type_id_from_metric_type(int metric_type, Env* env) {
  switch (metric_type) {
    case METRIC_TYPE_SORENSON:
      Insist(env, Bool_false ? "Unimplemented." : 0);
    case METRIC_TYPE_CZEKANOWSKI:
      return DATA_TYPE_ID_FLOAT;
    case METRIC_TYPE_CCC:
      Insist(env, Bool_false ? "Unimplemented." : 0);
  }
  Assert(Bool_false ? "Invalid metric type." : 0);
  return 0;
}

/*---------------------------------------------------------------------------*/

Bool_t Env_cuda_last_call_succeeded(Env* env) {
  Bool_t result = Bool_true;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if (error != cudaSuccess) {
    result = Bool_false;
    printf("CUDA error detected: %s\n", cudaGetErrorString(error));
  }

  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
