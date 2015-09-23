/*
    -- micMAGMA_minproduct (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_minproduct_ziteriluutils.cpp normal z -> d, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_minproductsparse.h"

#define PRECISION_d

/**
    Purpose
    -------

    Computes the Frobenius norm of the difference between the CSR matrices A
    and B. They need to share the same sparsity pattern!


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix in CSR

    @param[in]
    B           magma_minproduct_d_matrix
                sparse matrix in CSR
                
    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dfrobenius(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix B,
    real_Double_t *res,
    magma_minproduct_queue_t queue ){

    real_Double_t tmp2;
    magma_minproduct_int_t i,j,k;
    *res = 0.0;
    
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_minproduct_index_t localcol = A.col[j];
            for( k=B.row[i]; k<B.row[i+1]; k++){
                if(B.col[k] == localcol){
                    tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(A.val[j] )
                                                    - MAGMA_minproduct_D_REAL(B.val[k]) );

                    (*res) = (*res) + tmp2* tmp2;
                }
            }
        }
    }

    (*res) =  sqrt((*res));

    return MAGMA_minproduct_SUCCESS;
}



/**
    Purpose
    -------

    Computes the nonlinear residual A - LU and returns the difference as
    well es the Frobenius norm of the difference


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    L           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    U           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[out]
    LU          magma_minproduct_d_matrix*
                output sparse matrix in A-LU in CSR

    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dnonlinres(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix L,
    magma_minproduct_d_matrix U,
    magma_minproduct_d_matrix *LU,
    real_Double_t *res,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    real_Double_t tmp2;
    magma_minproduct_int_t i,j,k;
        
    double one = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );

    magma_minproduct_d_matrix L_d={Magma_minproduct_CSR}, U_d={Magma_minproduct_CSR}, LU_d={Magma_minproduct_CSR}, A_t={Magma_minproduct_CSR};

    CHECK( magma_minproduct_dmtransfer( L, &L_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue  ));
    CHECK( magma_minproduct_dmtransfer( U, &U_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue  ));
    CHECK( magma_minproduct_dmtransfer( A, &A_t, Magma_minproduct_CPU, Magma_minproduct_CPU, queue  ));
    CHECK( magma_minproduct_d_spmm( one, L_d, U_d, &LU_d, queue ));

    CHECK( magma_minproduct_dmtransfer(LU_d, LU, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    magma_minproduct_dmfree( &L_d, queue  );
    magma_minproduct_dmfree( &U_d, queue  );
    magma_minproduct_dmfree( &LU_d, queue  );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_minproduct_index_t lcol = A.col[j];
            double newval = MAGMA_minproduct_D_MAKE(0.0, 0.0);
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){
                    newval = MAGMA_minproduct_D_MAKE(
                        MAGMA_minproduct_D_REAL( LU->val[k] )- MAGMA_minproduct_D_REAL( A.val[j] )
                                                , 0.0 );
                }
            }
            A_t.val[j] = newval;
        }
    }

    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(A_t.val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }

    magma_minproduct_dmfree( LU, queue  );
    magma_minproduct_dmfree( &A_t, queue  );

    (*res) =  sqrt((*res));
    
cleanup:
    if( info !=0 ){
        magma_minproduct_dmfree( LU, queue  );
    }
    magma_minproduct_dmfree( &A_t, queue  );
    magma_minproduct_dmfree( &L_d, queue  );
    magma_minproduct_dmfree( &U_d, queue  );
    magma_minproduct_dmfree( &LU_d, queue  );
    return info;
}

/**
    Purpose
    -------

    Computes the ILU residual A - LU and returns the difference as
    well es the Frobenius norm of the difference


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    L           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    U           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[out]
    LU          magma_minproduct_d_matrix*
                output sparse matrix in A-LU in CSR
                
    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[out]
    nonlinres   real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dilures(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix L,
    magma_minproduct_d_matrix U,
    magma_minproduct_d_matrix *LU,
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    double tmp;
    real_Double_t tmp2;
    magma_minproduct_int_t i,j,k;
    
    double one = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );

    magma_minproduct_d_matrix LL={Magma_minproduct_CSR}, L_d={Magma_minproduct_CSR}, U_d={Magma_minproduct_CSR}, LU_d={Magma_minproduct_CSR};

    if( L.row[1]==1 ){        // lower triangular with unit diagonal
        //printf("L lower triangular.\n");
        LL.diagorder_type = Magma_minproduct_UNITY;
        CHECK( magma_minproduct_dmconvert( L, &LL, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
    }
    else if( L.row[1]==0 ){ // strictly lower triangular
        //printf("L strictly lower triangular.\n");
        CHECK( magma_minproduct_dmtransfer( L, &LL, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
        magma_minproduct_free_cpu( LL.col );
        magma_minproduct_free_cpu( LL.val );
        LL.nnz = L.nnz+L.num_rows;
        CHECK( magma_minproduct_dmalloc_cpu( &LL.val, LL.nnz ));
        CHECK( magma_minproduct_index_malloc_cpu( &LL.col, LL.nnz ));
        magma_minproduct_int_t z=0;
        for( magma_minproduct_int_t i=0; i<L.num_rows; i++){
            LL.row[i] = z;
            for( magma_minproduct_int_t j=L.row[i]; j<L.row[i+1]; j++){
                LL.val[z] = L.val[j];
                LL.col[z] = L.col[j];
                z++;
            }
            // add unit diagonal
            LL.val[z] = MAGMA_minproduct_D_MAKE(1.0, 0.0);
            LL.col[z] = i;
            z++;
        }
        LL.row[LL.num_rows] = z;
    }
    else{
        printf("error: L neither lower nor strictly lower triangular!\n");
    }

    CHECK( magma_minproduct_dmtransfer( LL, &L_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue  ));
    CHECK( magma_minproduct_dmtransfer( U, &U_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue  ));
    magma_minproduct_dmfree( &LL, queue );
    CHECK( magma_minproduct_d_spmm( one, L_d, U_d, &LU_d, queue ));



    CHECK( magma_minproduct_dmtransfer(LU_d, LU, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    magma_minproduct_dmfree( &L_d, queue );
    magma_minproduct_dmfree( &U_d, queue );
    magma_minproduct_dmfree( &LU_d, queue );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_minproduct_index_t lcol = A.col[j];
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){

                    tmp = MAGMA_minproduct_D_MAKE(
                        MAGMA_minproduct_D_REAL( LU->val[k] )- MAGMA_minproduct_D_REAL( A.val[j] )
                                                , 0.0 );
                    LU->val[k] = tmp;

                    tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(tmp) );
                    (*nonlinres) = (*nonlinres) + tmp2*tmp2;
                }

            }
        }
    }

    for(i=0; i<LU->num_rows; i++){
        for(j=LU->row[i]; j<LU->row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(LU->val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }

    (*res) =  sqrt((*res));
    (*nonlinres) =  sqrt((*nonlinres));

cleanup:
    if( info !=0 ){
        magma_minproduct_dmfree( LU, queue  );
    }
    magma_minproduct_dmfree( &LL, queue );
    magma_minproduct_dmfree( &L_d, queue  );
    magma_minproduct_dmfree( &U_d, queue  );
    magma_minproduct_dmfree( &LU_d, queue  );
    return info;
}



/**
    Purpose
    -------

    Computes the IC residual A - CC^T and returns the difference as
    well es the Frobenius norm of the difference


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    C           magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    CT          magma_minproduct_d_matrix
                input sparse matrix in CSR

    @param[in]
    LU          magma_minproduct_d_matrix*
                output sparse matrix in A-LU in CSR

    @param[out]
    res         real_Double_t*
                IC residual
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dicres(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix C,
    magma_minproduct_d_matrix CT,
    magma_minproduct_d_matrix *LU,
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;
    
    double tmp;
    real_Double_t tmp2;
    magma_minproduct_int_t i,j,k;

    double one = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );
    
    magma_minproduct_d_matrix L_d={Magma_minproduct_CSR}, U_d={Magma_minproduct_CSR}, LU_d={Magma_minproduct_CSR};
    
    *res = 0.0;
    *nonlinres = 0.0;

    CHECK( magma_minproduct_dmtransfer( C, &L_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_dmtransfer( CT, &U_d, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_d_spmm( one, L_d, U_d, &LU_d, queue ));
    CHECK( magma_minproduct_dmtransfer(LU_d, LU, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));

    magma_minproduct_dmfree( &LU_d, queue );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_minproduct_index_t lcol = A.col[j];
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){

                    tmp = MAGMA_minproduct_D_MAKE(
                        MAGMA_minproduct_D_REAL( LU->val[k] )- MAGMA_minproduct_D_REAL( A.val[j] )
                                                , 0.0 );
                    LU->val[k] = tmp;

                    tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(tmp) );
                    (*nonlinres) = (*nonlinres) + tmp2*tmp2;
                }
            }
        }
    }

    for(i=0; i<LU->num_rows; i++){
        for(j=LU->row[i]; j<LU->row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_minproduct_D_REAL(LU->val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }


    (*res) =  sqrt((*res));
    (*nonlinres) =  sqrt((*nonlinres));

cleanup:
    if( info !=0 ){
        magma_minproduct_dmfree( LU, queue  );
    }
    magma_minproduct_dmfree( &L_d, queue  );
    magma_minproduct_dmfree( &U_d, queue  );
    magma_minproduct_dmfree( &LU_d, queue  );
    return info;
}



/**
    Purpose
    -------

    Computes an initial guess for the iterative ILU/IC


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix
                sparse matrix in CSR

    @param[out]
    L           magma_minproduct_d_matrix*
                sparse matrix in CSR

    @param[out]
    U           magma_minproduct_d_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.


    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dinitguess(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix *L,
    magma_minproduct_d_matrix *U,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    double one = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );
    
    magma_minproduct_d_matrix hAL={Magma_minproduct_CSR}, hAU={Magma_minproduct_CSR}, dAL={Magma_minproduct_CSR}, 
    dAU={Magma_minproduct_CSR}, dALU={Magma_minproduct_CSR}, hALU={Magma_minproduct_CSR}, hD={Magma_minproduct_CSR}, 
    dD={Magma_minproduct_CSR}, dL={Magma_minproduct_CSR}, hL={Magma_minproduct_CSR};
    magma_minproduct_int_t i,j;
    
    magma_minproduct_int_t offdiags = 0;
    magma_minproduct_index_t *diag_offset;
    double *diag_vals=NULL;

    // need only lower triangular
    hAL.diagorder_type = Magma_minproduct_VALUE;
    CHECK( magma_minproduct_dmconvert( A, &hAL, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
    //magma_minproduct_dmconvert( hAL, &hALCOO, Magma_minproduct_CSR, Magma_minproduct_CSRCOO );

    // need only upper triangular
    //magma_minproduct_dmconvert( A, &hAU, Magma_minproduct_CSR, Magma_minproduct_CSRU );
    CHECK( magma_minproduct_d_cucsrtranspose(  hAL, &hAU, queue ));
    //magma_minproduct_dmconvert( hAU, &hAUCOO, Magma_minproduct_CSR, Magma_minproduct_CSRCOO );
    CHECK( magma_minproduct_dmtransfer( hAL, &dAL, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    CHECK( magma_minproduct_d_spmm( one, dAL, dAU, &dALU, queue ));
    CHECK( magma_minproduct_dmtransfer( dALU, &hALU, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));

    magma_minproduct_dmfree( &dAU, queue);
    magma_minproduct_dmfree( &dALU, queue);


    CHECK( magma_minproduct_dmalloc_cpu( &diag_vals, offdiags+1 ));
    CHECK( magma_minproduct_index_malloc_cpu( &diag_offset, offdiags+1 ));
    diag_offset[0] = 0;
    diag_vals[0] = MAGMA_minproduct_D_MAKE( 1.0, 0.0 );
    CHECK( magma_minproduct_dmgenerator( hALU.num_rows, offdiags, diag_offset, diag_vals, &hD, queue ));
    magma_minproduct_dmfree( &hALU, queue );

    
    for(i=0; i<hALU.num_rows; i++){
        for(j=hALU.row[i]; j<hALU.row[i+1]; j++){
            if( hALU.col[j] == i ){
                //printf("%d %d  %d == %d -> %f   -->", i, j, hALU.col[j], i, hALU.val[j]);
                hD.val[i] = MAGMA_minproduct_D_MAKE(
                        1.0 / sqrt(fabs(MAGMA_minproduct_D_REAL(hALU.val[j])))  , 0.0 );
                //printf("insert %f at %d\n", hD.val[i], i);
            }
        }
    }


    CHECK( magma_minproduct_dmtransfer( hD, &dD, Magma_minproduct_CPU, Magma_minproduct_DEV, queue ));
    magma_minproduct_dmfree( &hD, queue);

    CHECK( magma_minproduct_d_spmm( one, dD, dAL, &dL, queue ));
    magma_minproduct_dmfree( &dAL, queue );
    magma_minproduct_dmfree( &dD, queue );



/*
    // check for diagonal = 1
    magma_minproduct_d_matrix dLt={Magma_minproduct_CSR}, dLL={Magma_minproduct_CSR}, LL={Magma_minproduct_CSR};
    CHECK( magma_minproduct_d_cucsrtranspose(  dL, &dLt ));
    CHECK( magma_minproduct_dcuspmm( dL, dLt, &dLL ));
    CHECK( magma_minproduct_dmtransfer( dLL, &LL, Magma_minproduct_DEV, Magma_minproduct_CPU ));
    //for(i=0; i < hALU.num_rows; i++) {
    for(i=0; i < 100; i++) {
        for(j=hALU.row[i]; j < hALU.row[i+1]; j++) {
            if( hALU.col[j] == i ){
                printf("%d %d -> %f   -->", i, i, LL.val[j]);
            }
        }
    }
*/
    CHECK( magma_minproduct_dmtransfer( dL, &hL, Magma_minproduct_DEV, Magma_minproduct_CPU, queue ));
    CHECK( magma_minproduct_dmconvert( hL, L, Magma_minproduct_CSR, Magma_minproduct_CSRCOO, queue ));



cleanup:
    if( info !=0 ){
        magma_minproduct_dmfree( L, queue  );
        magma_minproduct_dmfree( U, queue  );
    }
    magma_minproduct_dmfree( &dAU, queue);
    magma_minproduct_dmfree( &dALU, queue);
    magma_minproduct_dmfree( &dL, queue );
    magma_minproduct_dmfree( &hL, queue );
    magma_minproduct_dmfree( &dAL, queue );
    magma_minproduct_dmfree( &dD, queue );
    magma_minproduct_dmfree( &hD, queue);
    magma_minproduct_dmfree( &hALU, queue );
    return info;
}



/**
    Purpose
    -------

    Using the iterative approach of computing ILU factorizations with increasing
    fill-in, it takes the input matrix A, containing the approximate factors,
    ( L and U as well )
    computes a matrix with one higher level of fill-in, inserts the original
    approximation as initial guess, and provides the factors L and U also
    filled with the scaled initial guess.


    Arguments
    ---------

    @param[in]
    A           magma_minproduct_d_matrix*
                sparse matrix in CSR

    @param[out]
    L           magma_minproduct_d_matrix*
                sparse matrix in CSR

    @param[out]
    U           magma_minproduct_d_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.


    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dinitrecursiveLU(
    magma_minproduct_d_matrix A,
    magma_minproduct_d_matrix *B,
    magma_minproduct_queue_t queue ){

    magma_minproduct_int_t i,j,k;

    for(i=0; i<A.num_rows; i++){
        for(j=B->row[i]; j<B->row[i+1]; j++){
            B->val[j] = MAGMA_minproduct_D_MAKE(0.0, 0.0);
            magma_minproduct_index_t localcol = B->col[j];
            for( k=A.row[i]; k<A.row[i+1]; k++){
                if(A.col[k] == localcol){
                    B->val[j] = A.val[k];
                }
            }
        }
    }

    return MAGMA_minproduct_SUCCESS; 
}



/**
    Purpose
    -------

    Checks for a lower triangular matrix whether it is strictly lower triangular
    and in the negative case adds a unit diagonal. It does this in-place.


    Arguments
    ---------

    @param[in,out]
    L           magma_minproduct_d_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_minproduct_queue_t
                Queue to execute in.

    @ingroup magma_minproductsparse_daux
    ********************************************************************/

magma_minproduct_int_t
magma_minproduct_dmLdiagadd(
    magma_minproduct_d_matrix *L,
    magma_minproduct_queue_t queue )
{
    magma_minproduct_int_t info = 0;

    magma_minproduct_d_matrix LL={Magma_minproduct_CSR};

    if( L->row[1]==1 ){        // lower triangular with unit diagonal
        //printf("L lower triangular.\n");
        LL.diagorder_type = Magma_minproduct_UNITY;
        CHECK( magma_minproduct_dmconvert( *L, &LL, Magma_minproduct_CSR, Magma_minproduct_CSRL, queue ));
    }
    else if( L->row[1]==0 ){ // strictly lower triangular
        //printf("L strictly lower triangular.\n");
        CHECK( magma_minproduct_dmtransfer( *L, &LL, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));
        magma_minproduct_free_cpu( LL.col );
        magma_minproduct_free_cpu( LL.val );
        LL.nnz = L->nnz+L->num_rows;
        CHECK( magma_minproduct_dmalloc_cpu( &LL.val, LL.nnz ));
        CHECK( magma_minproduct_index_malloc_cpu( &LL.col, LL.nnz ));
        magma_minproduct_int_t z=0;
        for( magma_minproduct_int_t i=0; i<L->num_rows; i++){
            LL.row[i] = z;
            for( magma_minproduct_int_t j=L->row[i]; j<L->row[i+1]; j++){
                LL.val[z] = L->val[j];
                LL.col[z] = L->col[j];
                z++;
            }
            // add unit diagonal
            LL.val[z] = MAGMA_minproduct_D_MAKE(1.0, 0.0);
            LL.col[z] = i;
            z++;
        }
        LL.row[LL.num_rows] = z;
        LL.nnz = z;
    }
    else{
        printf("error: L neither lower nor strictly lower triangular!\n");
    }
    magma_minproduct_dmfree( L, queue );
    CHECK( magma_minproduct_dmtransfer(LL, L, Magma_minproduct_CPU, Magma_minproduct_CPU, queue ));

cleanup:
    if( info != 0 ){
        magma_minproduct_dmfree( L, queue );
    }
    magma_minproduct_dmfree( &LL, queue );
    return info;
}


