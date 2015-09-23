/*
    -- MAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_zutil_sparse.cpp normal z -> s, Sun May  3 11:23:01 2015

       @author Hartwig Anzt

       Utilities for testing MAGMA_minproduct-sparse.
*/
#include "common_magma_minproductsparse.h"


// --------------------
static const char *usage_sparse_short =
"Usage: %s [options] [-h|--help]  matrices\n\n";

static const char *usage_sparse =
"Options are:\n"
" --format x     Possibility to choose a format for the sparse matrix:\n"
"               0   CSR\n"
"               1   ELL\n"
"               2   SELL-P\n"
" --blocksize x Set a specific blocksize for SELL-P format.\n"
" --alignment x Set a specific alignment for SELL-P format.\n"
" --mscale      Possibility to scale the original matrix:\n"
"               0   no scaling\n"
"               1   symmetric scaling to unit diagonal\n"
" --solver      Possibility to choose a solver:\n"
"               0   CG\n"
"               1   merged CG\n"
"               2   preconditioned CG\n"
"               3   BiCGSTAB\n"
"               4   merged BiCGSTAB\n"
"               5   preconditioned BiCGSTAB\n"
"               6   GMRES\n"
"               7   preconditioned GMRES\n"
"               8   LOBPCG\n"
"               9   Jacobi\n"
"               10  Block-asynchronous Iteration\n"
"               21  Iterative Refinement\n"
" --restart     For GMRES: possibility to choose the restart.\n"
" --precond x   Possibility to choose a preconditioner:\n"
"               0   no preconditioner\n"
"               1   Jacobi\n"
"               2   ILU(0) / IC(0)\n"
"               -2   iterative ILU(0) / IC(0)\n"
"                   For Iterative Refinement also possible: \n"
"                   3   CG\n"
"                   4   BiCGSTAB\n"
"                   5   GMRES\n"
"                   6   Block-asynchronous Iteration\n"
"                   --ptol eps    Relative resiudal stopping criterion for preconditioner.\n"
"                   --psweeps k   Iteration count for iterative incomplete factorizations.\n"
"                   --piter k     Iteration count for iterative preconditioner.\n"
"                   --plevels k   Number of ILU levels.\n"
" --ev x        For eigensolvers, set number of eigenvalues/eigenvectors to compute.\n"
" --verbose x   Possibility to print intermediate residuals every x iteration.\n"
" --maxiter x   Set an upper limit for the iteration count.\n"
" --tol eps     Set a relative residual stopping criterion.\n";



/**
    Purpose
    -------

    Parses input options for a solver

    Arguments
    ---------

    @param[in]
    argc            int
                    command line input
                
    @param[in]
    argv            char**
                    command line input

    @param[in,out]
    opts            magma_minproduct_sopts *
                    magma_minproduct solver options

    @param[out]
    matrices        int
                    counter how many linear systems to process

    @param[in]
    queue           magma_minproduct_queue_t
                    Queue to execute in.

    @ingroup magma_minproductsparse_saux
    ********************************************************************/

extern "C"
magma_minproduct_int_t
magma_minproduct_sparse_opts(
    int argc,
    char** argv,
    magma_minproduct_sopts *opts,
    int *matrices,
    magma_minproduct_queue_t queue )
{

    
    // fill in default values
    opts->input_format = Magma_minproduct_CSR;
    opts->blocksize = 8;
    opts->alignment = 8;
    opts->output_format = Magma_minproduct_CSR;
    opts->input_location = Magma_minproduct_CPU;
    opts->output_location = Magma_minproduct_CPU;
    opts->scaling = Magma_minproduct_NOSCALE;
    opts->solver_par.epsilon = 10e-16;
    opts->solver_par.maxiter = 1000;
    opts->solver_par.verbose = 0;
    opts->solver_par.version = 0;
    opts->solver_par.restart = 30;
    opts->solver_par.num_eigenvalues = 0;
    opts->precond_par.solver = Magma_minproduct_NONE;
    opts->precond_par.epsilon = 0.01;
    opts->precond_par.maxiter = 100;
    opts->precond_par.restart = 10;
    opts->precond_par.levels = 0;
    opts->precond_par.sweeps = 5;
    opts->solver_par.solver = Magma_minproduct_CG;
    
    printf( usage_sparse_short, argv[0] );
    
    int ndevices;
    cudaGetDeviceCount( &ndevices );
    
    int info;
    

    for( int i = 1; i < argc; ++i ) {
        if ( strcmp("--format", argv[i]) == 0 && i+1 < argc ) {
            info = atoi( argv[++i] );
            switch( info ) {
                case 0: opts->output_format = Magma_minproduct_CSR; break;
                case 1: opts->output_format = Magma_minproduct_ELL; break;
                case 2: opts->output_format = Magma_minproduct_SELLP; break;
                //case 2: opts->output_format = Magma_minproduct_ELLRT; break;
            }
        } else if ( strcmp("--mscale", argv[i]) == 0 && i+1 < argc ) {
            info = atoi( argv[++i] );
            switch( info ) {
                case 0: opts->scaling = Magma_minproduct_NOSCALE; break;
                case 1: opts->scaling = Magma_minproduct_UNITDIAG; break;
                case 2: opts->scaling = Magma_minproduct_UNITROW; break;
            }

        } else if ( strcmp("--solver", argv[i]) == 0 && i+1 < argc ) {
            info = atoi( argv[++i] );
            switch( info ) {
                case 0: opts->solver_par.solver = Magma_minproduct_CG; break;
                case 1: opts->solver_par.solver = Magma_minproduct_CGMERGE; break;
                case 2: opts->solver_par.solver = Magma_minproduct_PCG; break;
                case 3: opts->solver_par.solver = Magma_minproduct_BICGSTAB; break;
                case 4: opts->solver_par.solver = Magma_minproduct_BICGSTABMERGE; break;
                case 5: opts->solver_par.solver = Magma_minproduct_PBICGSTAB; break;
                case 6: opts->solver_par.solver = Magma_minproduct_GMRES; break;
                case 7: opts->solver_par.solver = Magma_minproduct_PGMRES; break;
                case 8: opts->solver_par.solver = Magma_minproduct_LOBPCG;
                            opts->solver_par.num_eigenvalues = 16;break;
                case 9: opts->solver_par.solver = Magma_minproduct_JACOBI; break;
                case 10: opts->solver_par.solver = Magma_minproduct_BAITER; break;
                case 21: opts->solver_par.solver = Magma_minproduct_ITERREF; break;
            }
        } else if ( strcmp("--restart", argv[i]) == 0 && i+1 < argc ) {
            opts->solver_par.restart = atoi( argv[++i] );
        } else if ( strcmp("--precond", argv[i]) == 0 && i+1 < argc ) {
            info = atoi( argv[++i] );
            switch( info ) {
                case 0: opts->precond_par.solver = Magma_minproduct_NONE; break;
                case 1: opts->precond_par.solver = Magma_minproduct_JACOBI; break;
                case 2: opts->precond_par.solver = Magma_minproduct_ILU; break;
                case -2: opts->precond_par.solver = Magma_minproduct_AILU; break;
                case 3: opts->precond_par.solver = Magma_minproduct_CG; break;
                case 4: opts->precond_par.solver = Magma_minproduct_BICGSTAB; break;
                case 5: opts->precond_par.solver = Magma_minproduct_GMRES; break;
                case 6: opts->precond_par.solver = Magma_minproduct_BAITER; break;

            }
        } else if ( strcmp("--ptol", argv[i]) == 0 && i+1 < argc ) {
            sscanf( argv[++i], "%f", &opts->precond_par.epsilon );
        } else if ( strcmp("--piter", argv[i]) == 0 && i+1 < argc ) {
            opts->precond_par.maxiter = atoi( argv[++i] );
        } else if ( strcmp("--psweeps", argv[i]) == 0 && i+1 < argc ) {
            opts->precond_par.sweeps = atoi( argv[++i] );
        } else if ( strcmp("--plevels", argv[i]) == 0 && i+1 < argc ) {
            opts->precond_par.levels = atoi( argv[++i] );
        } else if ( strcmp("--blocksize", argv[i]) == 0 && i+1 < argc ) {
            opts->blocksize = atoi( argv[++i] );
        } else if ( strcmp("--alignment", argv[i]) == 0 && i+1 < argc ) {
            opts->alignment = atoi( argv[++i] );
        } else if ( strcmp("--verbose", argv[i]) == 0 && i+1 < argc ) {
            opts->solver_par.verbose = atoi( argv[++i] );
        }  else if ( strcmp("--maxiter", argv[i]) == 0 && i+1 < argc ) {
            opts->solver_par.maxiter = atoi( argv[++i] );
        } else if ( strcmp("--tol", argv[i]) == 0 && i+1 < argc ) {
            sscanf( argv[++i], "%f", &opts->solver_par.epsilon );
        } else if ( strcmp("--ev", argv[i]) == 0 && i+1 < argc ) {
            opts->solver_par.num_eigenvalues = atoi( argv[++i] );
        } else if ( strcmp("--version", argv[i]) == 0 && i+1 < argc ) {
            opts->solver_par.version = atoi( argv[++i] );
        }
        // ----- usage
        else if ( strcmp("-h",     argv[i]) == 0 ||
                  strcmp("--help", argv[i]) == 0 ) {
            fprintf( stderr, usage_sparse, argv[0] );
        } else {
            *matrices = i;
            break;
        }
    }
    // ensure to take a symmetric preconditioner for the PCG
    if ( opts->solver_par.solver == Magma_minproduct_PCG
        && opts->precond_par.solver == Magma_minproduct_ILU )
            opts->precond_par.solver = Magma_minproduct_ICC;
    if ( opts->solver_par.solver == Magma_minproduct_PCG
        && opts->precond_par.solver == Magma_minproduct_AILU )
            opts->precond_par.solver = Magma_minproduct_AICC;
            
    return MAGMA_minproduct_SUCCESS;
}

    

