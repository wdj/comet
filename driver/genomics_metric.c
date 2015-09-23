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

#include "mpi.h"

#include "env.h"

/*===========================================================================*/

int main(int argc, char** argv) {

  MPI_Init(&argc, &argv);

  Env env;
  Env_create_from_args(&env, argc, argv);











    Env_destroy(&env);
    MPI_Finalize();
}

/*---------------------------------------------------------------------------*/
