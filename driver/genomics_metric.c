/*---------------------------------------------------------------------------*/
/*!
 * \file   driver.c
 * \author Wayne Joubert
 * \date   Wed Sep 23 12:39:13 EDT 2015
 * \brief  Main driver code for genomics metric calculation.
 * \note   Copyright (C) 2015 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
/*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi.h"

#include "env.h"
#include "vectors.h"
#include "metrics.h"

/*===========================================================================*/

void load_vectors(Vectors* vectors, Env* env) {
    int vector_local;
    for (vector_local=0; vector_local<vectors->num_vector_local;
         ++vector_local) {
        int field;
        for (field=0; field<vectors->num_field; ++field) {
            /*---Make up a number between 0 and 2---*/
            double value = 2. *
                ( field + vectors->num_field * (double)(
                  vector_local + vectors->num_vector_local_max * (double)(
                  env->proc_num ) ) ) /
                ( vectors->num_field * (double)
                  vectors->num_vector_local_max * (double)
                  env->num_proc );
            Vectors_set(vectors, field, vector_local, value, env);
        }
    }
}

/*===========================================================================*/

#if 0
void store_metrics(Metrics* metrics, Env* env) {


}
#endif

/*===========================================================================*/

int main(int argc, char** argv) {

    /*---Initialize---*/

    MPI_Init(&argc, &argv);

    Env env;
    Env_create_from_args(&env, argc, argv);

    /*---Parse remaining unprocessed arguments---*/

    int num_field = -1;
    int num_vector_local = -1;

    int i;
    for (i=0; i<argc; ++i) {

        if (strcmp(argv[i], "--num_field")==0) {
            ++i;
            Insist(i<argc ? "Missing value for num_field." : 0 );
            num_field = atoi(argv[i]);
            Insist(num_field>=0 ? "Invalid setting for num_field." : 0);

        } else if (strcmp(argv[i], "--num_vector_local")==0) {
            ++i;
            Insist(i<argc ? "Missing value for num_vector_local." : 0 );
            num_vector_local = atoi(argv[i]);
            Insist(num_field>=0 ? "Invalid setting for num_vector_local." : 0);

        } /*---if/else---*/
    } /*---for i---*/

    Insist(num_field != -1 ? "Error: num_field not set." : 0);
    Insist(num_vector_local != -1 ? "Error: num_vector_local not set." : 0);

    /*---Initialize vectors---*/

    Vectors vectors;
    Vectors_create(&vectors, data_type_id_from_metric_type(env.metric_type),
                   num_field, num_vector_local, &env);

    load_vectors(&vectors, &env);

    /*---Set up metrics---*/

#if 0
    Metrics metrics;
#endif



    /*---Calculate metrics---*/

#if 0
    compute_metrics(&metrics, &vectors, &env);
#endif


    /*---Output results---*/






    /*---Finalize---*/

    Vectors_destroy(&vectors, &env);
    Env_destroy(&env);
    MPI_Finalize();
    return 0;
}

/*---------------------------------------------------------------------------*/
