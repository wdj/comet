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

#include <stddef.h>
#include <assert.h>

#include "mpi.h"

/*===========================================================================*/
/*---Types---*/

/*---Boolean type---*/

typedef int Bool_t;

enum { Bool_true = (1 == 1), Bool_false = (1 == 0) };

/*---Default floating point type---*/

#ifdef FP_PRECISION_SINGLE
  typedef float Float_t;
  enum { MPI_Float_t = MPI_FLOAT };
#endif

#ifdef FP_PRECISION_DOUBLE
  typedef double Float_t;
  enum { MPI_Float_t = MPI_DOUBLE };
#endif

/*---Type ids---*/

enum { DATA_TYPE_ID_FLOAT = 1, DATA_TYPE_ID_BIT = 2 };

/*===========================================================================*/
/*---Environment struct declarations---*/

typedef struct {
  int metric_type;
  int num_way;
  Bool_t all2all;
  int compute_method;
  int mpi_comm;
  int num_proc;
  int proc_num;
} Env;

enum {
  METRIC_TYPE_SORENSON = 0,
  METRIC_TYPE_CZEKANOWSKI = 1,
  METRIC_TYPE_CCC = 2,
  NUM_METRIC_TYPE = 3
};

enum { COMPUTE_METHOD_CPU = 0,
       COMPUTE_METHOD_GPU = 1,
       NUM_COMPUTE_METHOD = 2 };

/*===========================================================================*/
/*---Null object---*/

Env Env_null(void);

/*===========================================================================*/
/*---Initialize environment---*/

void Env_create(Env* env);

void Env_create_from_args(Env* env, int argc, char** argv);

/*===========================================================================*/
/*---Finalize environment---*/

void Env_destroy(Env* env);

/*===========================================================================*/
/*---Timer functions---*/

double Env_get_time(Env* env);
double Env_get_synced_time(Env* env);

/*===========================================================================*/
/*---Assertions---*/

#define Assert(v) assert(v)

#ifndef Insist
#define Insist(env, condition) \
  (void)((condition) || (insist_(env, #condition, __FILE__, __LINE__), 0))
#endif

void insist_(Env* env,
             const char* condition_string,
             const char* file,
             int line);

#ifndef NDEBUG
/*---Fail compilation and (hopefully) give a filename/line number---*/
#define Static_Assert( condition ) { int a[ ( condition ) ? 1 : -1 ]; (void)a; }
#else
#define Static_Assert( condition )
#endif

/*===========================================================================*/
/*---Misc utility functions---*/

static int imin( const int i,
                 const int j )
{
    return i < j ? i : j;
}

/*---------------------------------------------------------------------------*/

static int imax( const int i,
                 const int j )
{
    return i > j ? i : j;
}

/*---------------------------------------------------------------------------*/

static int ifloor( const int i,
                   const int j )
{
    Assert( j > 0 );

    return i >= 0 ? i/j : (i-j+1)/j;
}

/*---------------------------------------------------------------------------*/

static int iceil( const int i,
                  const int j )
{
    Assert( j > 0 );

    return - ifloor( -i, j );
}

/*---------------------------------------------------------------------------*/

static size_t randomize( size_t i ) {
  const size_t im = 714025;
  const size_t ia = 4096;
  const size_t ic = 150889;

  return ( i * ia + ic ) % im;
}

/*---------------------------------------------------------------------------*/

static size_t randomize_max() {
  const size_t im = 714025;

  return im;
}

/*---------------------------------------------------------------------------*/

static size_t nchoosek( int n, int k ) {
  Assert( n >= 1 );
  Assert( k >= 0 && k <= n );
  int i;
  size_t num = 1;
  size_t denom = 1;
  for ( i = 0; i < k; ++i ) {
    num *= ( n - i );
    denom *= ( i + 1 );
  }
  return num / denom;
}

/*===========================================================================*/
/*---Misc.---*/

int data_type_id_from_metric_type(int metric_type, Env* env);

Bool_t Env_cuda_last_call_succeeded(Env* env);

/*===========================================================================*/

#endif /*---_env_h_---*/

/*---------------------------------------------------------------------------*/
