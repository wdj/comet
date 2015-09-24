/*---------------------------------------------------------------------------*/
/*!
 * \file   env.h
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Environment settings and general utilities, header.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

/*=============================================================================


=============================================================================*/

#ifndef _env_h_
#define _env_h_

#include <assert.h>

#include "mpi.h"

/*===========================================================================*/
/*---Assertions---*/

#define Assert(v) assert(v)

#ifndef Insist
#define Insist(condition) \
    (void)((condition) || (insist_ (#condition, __FILE__, __LINE__),0))
#endif

void insist_(const char *condition_string, const char *file, int line);

/*===========================================================================*/
/*---Types---*/

/*---Boolean type---*/

typedef int Bool_t;

enum { Bool_true = 1,
       Bool_false = 0
     };

/*---Default floating point type---*/

typedef double Float_t;
enum { MPI_Float_t = MPI_DOUBLE };

/*---Type ids---*/

enum { DATA_TYPE_ID_FLOAT = 1,
       DATA_TYPE_ID_BIT = 2
     };

/*===========================================================================*/
/*---Environment struct declarations---*/

typedef struct {
    int    metric_type;
    int    num_way;
    Bool_t global_all2all;
    int    compute_method;
    int    mpi_comm;
    int    num_proc;
    int    proc_num;
} Env;

enum { METRIC_TYPE_SORENSON = 1,
       METRIC_TYPE_CZEKANOWSKI = 2,
       METRIC_TYPE_CCC = 3
     };

enum { COMPUTE_METHOD_CPU = 0,
       COMPUTE_METHOD_GPU = 1
     };

void Env_create(Env* env);

void Env_create_from_args(Env* env, int argc, char** argv);

void Env_destroy(Env* env);

int data_type_id_from_metric_type(int metric_type);

/*===========================================================================*/
/*---Timer functions---*/

double Env_get_time(Env* env);
double Env_get_synced_time(Env* env);

/*===========================================================================*/
/*---Misc.---*/

Bool_t Env_cuda_last_call_succeeded(Env* env);

/*===========================================================================*/

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
