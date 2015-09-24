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
/*---Assertions---*/

void insist_(const char *condition_string, const char *file, int line) {
    fprintf( stderr, "Insist error: \"%s\". At file %s, line %i.\n",
                     condition_string, file, line );
    exit( EXIT_FAILURE );
}

/*===========================================================================*/
/*---Initialize environment---*/

void Env_create(Env* env) {

    int mpi_code = 0;
    env->mpi_comm = MPI_COMM_WORLD;
    mpi_code = MPI_Comm_rank(env->mpi_comm, &(env->proc_num));
    Assert( mpi_code == MPI_SUCCESS );
    mpi_code = MPI_Comm_size(env->mpi_comm, &(env->num_proc));
    Assert( mpi_code == MPI_SUCCESS );

    /*---Set default values---*/
    env->metric_type = METRIC_TYPE_CZEKANOWSKI;
    env->num_way = 2;
    env->global_all2all = Bool_false;
    env->compute_method = COMPUTE_METHOD_GPU;
}

/*===========================================================================*/
/*---Initialize environment---*/

void Env_create_from_args(Env* env, int argc, char** argv) {

    Env_create(env);

    int i = 0;
    for (i=0; i<argc; ++i) {

        if (strcmp(argv[i], "--metric_type")==0) {
            ++i;
            Insist(i<argc ? "Missing value for metric_type." : 0 );
            env->metric_type = atoi(argv[i]);
            Insist(env->metric_type==METRIC_TYPE_SORENSON ||
                   env->metric_type==METRIC_TYPE_CZEKANOWSKI ||
                   env->metric_type==METRIC_TYPE_CCC ?
                   "Invalid setting for metric_type." : 0);

        } else if (strcmp(argv[i], "--num_way")==0) {
            ++i;
            Insist(i<argc ? "Missing value for num_way." : 0 );
            env->num_way = atoi(argv[i]);
            Insist(env->num_way==2 ||
                   env->num_way==3 ?
                   "Invalid setting for num_way." : 0);

        } else if (strcmp(argv[i], "--global_all2all")==0) {
            ++i;
            Insist(i<argc ? "Missing value for global_all2all." : 0 );
            env->global_all2all = atoi(argv[i]);
            Insist(env->global_all2all==Bool_false ||
                   env->global_all2all==Bool_true ?
                   "Invalid setting for global_all2all." : 0);

        } else if (strcmp(argv[i], "--compute_method")==0) {
            ++i;
            Insist(i<argc ? "Missing value for compute_method." : 0 );
            env->compute_method = atoi(argv[i]);
            Insist(env->compute_method==COMPUTE_METHOD_CPU ||
                   env->compute_method==COMPUTE_METHOD_GPU ?
                   "Invalid setting for compute_method." : 0);

        } /*---if/else---*/
    } /*---for i---*/
}

/*===========================================================================*/
/*---Finalize environment---*/

void Env_destroy(Env* env) {
    /*---Nothing to do presently---*/
}

/*===========================================================================*/
/*---Timer functions---*/

double Env_get_time(Env* env) {
    Assert(env);
    struct timeval tv;
    double result;
    gettimeofday( &tv, NULL );
    result = ( (double) tv.tv_sec +
               (double) tv.tv_usec * 1.e-6 );
    return result;
}

/*---------------------------------------------------------------------------*/

double Env_get_synced_time(Env* env) {
    Assert(env);

    cudaThreadSynchronize();
    Assert( Env_cuda_last_call_succeeded(env) );

    int mpi_code = MPI_Barrier(env->mpi_comm);
    Assert( mpi_code == MPI_SUCCESS );
    return Env_get_time(env);
}

/*===========================================================================*/
/*---Misc.---*/

int data_type_id_from_metric_type(int metric_type) {
    if (metric_type==METRIC_TYPE_CZEKANOWSKI) {
        return DATA_TYPE_ID_FLOAT;
    }
    Insist(Bool_false ? "Invalid metric type." : 0);
    return 0;
}

/*---------------------------------------------------------------------------*/

Bool_t Env_cuda_last_call_succeeded(Env* env) {
  Bool_t result = Bool_true;

  /*---NOTE: this read of the last error is a destructive read---*/
  cudaError_t error = cudaGetLastError();

  if ( error != cudaSuccess )
  {
      result = Bool_false;
      printf( "CUDA error detected: %s\n", cudaGetErrorString( error ) );
  }

  return result;
}

/*===========================================================================*/

/*---------------------------------------------------------------------------*/
