/*---------------------------------------------------------------------------*/
/*!
 * \file   metrics.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Metrics pseudo-class.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdlib.h>

#include "mpi.h"

#include "env.h"
#include "metrics.h"

/*===========================================================================*/
/*---Metrics pseudo-constructor---*/

void Metrics_create(Metrics* metrics, int data_type_id, int num_vector_local,
                    Env* env) {
    Assert(metrics);
    Assert(num_vector_local>=0);
    Assert(env);

    if (data_type_id == DATA_TYPE_ID_FLOAT) {
/*
FIX
        metrics->data = malloc(num_field*num_vector_local*sizeof(Float_t));
*/
    } else {
        Insist(Bool_false ? "Invalid data type." : 0);
    }
    metrics->data_type_id = data_type_id;
    metrics->num_vector_local = num_vector_local;

    int mpi_code;
    mpi_code = MPI_Allreduce(&(metrics->num_vector_local),
        &(metrics->num_vector), 1, MPI_INT, MPI_SUM, env->mpi_comm);
    Assert( mpi_code == MPI_SUCCESS );
    mpi_code = MPI_Allreduce(&(metrics->num_vector_local),
        &(metrics->num_vector_local_max), 1, MPI_INT, MPI_MAX, env->mpi_comm);
    Assert( mpi_code == MPI_SUCCESS );
}

/*===========================================================================*/
/*---Metrics pseudo-destructor---*/

void Metrics_destroy(Metrics* metrics, Env* env) {
    Assert(metrics);
    Assert(metrics->data);
    Assert(env);

    free(metrics->data);
    metrics->data = 0;
}

/*---------------------------------------------------------------------------*/
