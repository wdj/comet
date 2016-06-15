/*
    -- micMAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_ziteriluutils.cpp normal z -> s, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/
#include "common_magma_tally3sparse.h"

#define PRECISION_s

/**
    Purpose
    -------

    Computes the Frobenius norm of the difference between the CSR matrices A
    and B. They need to share the same sparsity pattern!


    Arguments
    ---------

    @param[in]
    A           magma_tally3_s_matrix
                sparse matrix in CSR

    @param[in]
    B           magma_tally3_s_matrix
                sparse matrix in CSR
                
    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_sfrobenius(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix B,
    real_Double_t *res,
    magma_tally3_queue_t queue ){

    real_Double_t tmp2;
    magma_tally3_int_t i,j,k;
    *res = 0.0;
    
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_tally3_index_t localcol = A.col[j];
            for( k=B.row[i]; k<B.row[i+1]; k++){
                if(B.col[k] == localcol){
                    tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(A.val[j] )
                                                    - MAGMA_tally3_S_REAL(B.val[k]) );

                    (*res) = (*res) + tmp2* tmp2;
                }
            }
        }
    }

    (*res) =  sqrt((*res));

    return MAGMA_tally3_SUCCESS;
}



/**
    Purpose
    -------

    Computes the nonlinear residual A - LU and returns the difference as
    well es the Frobenius norm of the difference


    Arguments
    ---------

    @param[in]
    A           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    L           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    U           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[out]
    LU          magma_tally3_s_matrix*
                output sparse matrix in A-LU in CSR

    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_snonlinres(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U,
    magma_tally3_s_matrix *LU,
    real_Double_t *res,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    real_Double_t tmp2;
    magma_tally3_int_t i,j,k;
        
    float one = MAGMA_tally3_S_MAKE( 1.0, 0.0 );

    magma_tally3_s_matrix L_d={Magma_tally3_CSR}, U_d={Magma_tally3_CSR}, LU_d={Magma_tally3_CSR}, A_t={Magma_tally3_CSR};

    CHECK( magma_tally3_smtransfer( L, &L_d, Magma_tally3_CPU, Magma_tally3_DEV, queue  ));
    CHECK( magma_tally3_smtransfer( U, &U_d, Magma_tally3_CPU, Magma_tally3_DEV, queue  ));
    CHECK( magma_tally3_smtransfer( A, &A_t, Magma_tally3_CPU, Magma_tally3_CPU, queue  ));
    CHECK( magma_tally3_s_spmm( one, L_d, U_d, &LU_d, queue ));

    CHECK( magma_tally3_smtransfer(LU_d, LU, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
    magma_tally3_smfree( &L_d, queue  );
    magma_tally3_smfree( &U_d, queue  );
    magma_tally3_smfree( &LU_d, queue  );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_tally3_index_t lcol = A.col[j];
            float newval = MAGMA_tally3_S_MAKE(0.0, 0.0);
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){
                    newval = MAGMA_tally3_S_MAKE(
                        MAGMA_tally3_S_REAL( LU->val[k] )- MAGMA_tally3_S_REAL( A.val[j] )
                                                , 0.0 );
                }
            }
            A_t.val[j] = newval;
        }
    }

    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(A_t.val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }

    magma_tally3_smfree( LU, queue  );
    magma_tally3_smfree( &A_t, queue  );

    (*res) =  sqrt((*res));
    
cleanup:
    if( info !=0 ){
        magma_tally3_smfree( LU, queue  );
    }
    magma_tally3_smfree( &A_t, queue  );
    magma_tally3_smfree( &L_d, queue  );
    magma_tally3_smfree( &U_d, queue  );
    magma_tally3_smfree( &LU_d, queue  );
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
    A           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    L           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    U           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[out]
    LU          magma_tally3_s_matrix*
                output sparse matrix in A-LU in CSR
                
    @param[out]
    res         real_Double_t*
                Frobenius norm of difference
                
    @param[out]
    nonlinres   real_Double_t*
                Frobenius norm of difference
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_silures(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix L,
    magma_tally3_s_matrix U,
    magma_tally3_s_matrix *LU,
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    float tmp;
    real_Double_t tmp2;
    magma_tally3_int_t i,j,k;
    
    float one = MAGMA_tally3_S_MAKE( 1.0, 0.0 );

    magma_tally3_s_matrix LL={Magma_tally3_CSR}, L_d={Magma_tally3_CSR}, U_d={Magma_tally3_CSR}, LU_d={Magma_tally3_CSR};

    if( L.row[1]==1 ){        // lower triangular with unit diagonal
        //printf("L lower triangular.\n");
        LL.diagorder_type = Magma_tally3_UNITY;
        CHECK( magma_tally3_smconvert( L, &LL, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
    }
    else if( L.row[1]==0 ){ // strictly lower triangular
        //printf("L strictly lower triangular.\n");
        CHECK( magma_tally3_smtransfer( L, &LL, Magma_tally3_CPU, Magma_tally3_CPU, queue ));
        magma_tally3_free_cpu( LL.col );
        magma_tally3_free_cpu( LL.val );
        LL.nnz = L.nnz+L.num_rows;
        CHECK( magma_tally3_smalloc_cpu( &LL.val, LL.nnz ));
        CHECK( magma_tally3_index_malloc_cpu( &LL.col, LL.nnz ));
        magma_tally3_int_t z=0;
        for( magma_tally3_int_t i=0; i<L.num_rows; i++){
            LL.row[i] = z;
            for( magma_tally3_int_t j=L.row[i]; j<L.row[i+1]; j++){
                LL.val[z] = L.val[j];
                LL.col[z] = L.col[j];
                z++;
            }
            // add unit diagonal
            LL.val[z] = MAGMA_tally3_S_MAKE(1.0, 0.0);
            LL.col[z] = i;
            z++;
        }
        LL.row[LL.num_rows] = z;
    }
    else{
        printf("error: L neither lower nor strictly lower triangular!\n");
    }

    CHECK( magma_tally3_smtransfer( LL, &L_d, Magma_tally3_CPU, Magma_tally3_DEV, queue  ));
    CHECK( magma_tally3_smtransfer( U, &U_d, Magma_tally3_CPU, Magma_tally3_DEV, queue  ));
    magma_tally3_smfree( &LL, queue );
    CHECK( magma_tally3_s_spmm( one, L_d, U_d, &LU_d, queue ));



    CHECK( magma_tally3_smtransfer(LU_d, LU, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
    magma_tally3_smfree( &L_d, queue );
    magma_tally3_smfree( &U_d, queue );
    magma_tally3_smfree( &LU_d, queue );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_tally3_index_t lcol = A.col[j];
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){

                    tmp = MAGMA_tally3_S_MAKE(
                        MAGMA_tally3_S_REAL( LU->val[k] )- MAGMA_tally3_S_REAL( A.val[j] )
                                                , 0.0 );
                    LU->val[k] = tmp;

                    tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(tmp) );
                    (*nonlinres) = (*nonlinres) + tmp2*tmp2;
                }

            }
        }
    }

    for(i=0; i<LU->num_rows; i++){
        for(j=LU->row[i]; j<LU->row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(LU->val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }

    (*res) =  sqrt((*res));
    (*nonlinres) =  sqrt((*nonlinres));

cleanup:
    if( info !=0 ){
        magma_tally3_smfree( LU, queue  );
    }
    magma_tally3_smfree( &LL, queue );
    magma_tally3_smfree( &L_d, queue  );
    magma_tally3_smfree( &U_d, queue  );
    magma_tally3_smfree( &LU_d, queue  );
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
    A           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    C           magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    CT          magma_tally3_s_matrix
                input sparse matrix in CSR

    @param[in]
    LU          magma_tally3_s_matrix*
                output sparse matrix in A-LU in CSR

    @param[out]
    res         real_Double_t*
                IC residual
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_sicres(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix C,
    magma_tally3_s_matrix CT,
    magma_tally3_s_matrix *LU,
    real_Double_t *res,
    real_Double_t *nonlinres,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    float tmp;
    real_Double_t tmp2;
    magma_tally3_int_t i,j,k;

    float one = MAGMA_tally3_S_MAKE( 1.0, 0.0 );
    
    magma_tally3_s_matrix L_d={Magma_tally3_CSR}, U_d={Magma_tally3_CSR}, LU_d={Magma_tally3_CSR};
    
    *res = 0.0;
    *nonlinres = 0.0;

    CHECK( magma_tally3_smtransfer( C, &L_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_smtransfer( CT, &U_d, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_s_spmm( one, L_d, U_d, &LU_d, queue ));
    CHECK( magma_tally3_smtransfer(LU_d, LU, Magma_tally3_DEV, Magma_tally3_CPU, queue ));

    magma_tally3_smfree( &LU_d, queue );

    // compute Frobenius norm of A-LU
    for(i=0; i<A.num_rows; i++){
        for(j=A.row[i]; j<A.row[i+1]; j++){
            magma_tally3_index_t lcol = A.col[j];
            for(k=LU->row[i]; k<LU->row[i+1]; k++){
                if( LU->col[k] == lcol ){

                    tmp = MAGMA_tally3_S_MAKE(
                        MAGMA_tally3_S_REAL( LU->val[k] )- MAGMA_tally3_S_REAL( A.val[j] )
                                                , 0.0 );
                    LU->val[k] = tmp;

                    tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(tmp) );
                    (*nonlinres) = (*nonlinres) + tmp2*tmp2;
                }
            }
        }
    }

    for(i=0; i<LU->num_rows; i++){
        for(j=LU->row[i]; j<LU->row[i+1]; j++){
            tmp2 = (real_Double_t) fabs( MAGMA_tally3_S_REAL(LU->val[j]) );
            (*res) = (*res) + tmp2* tmp2;
        }
    }


    (*res) =  sqrt((*res));
    (*nonlinres) =  sqrt((*nonlinres));

cleanup:
    if( info !=0 ){
        magma_tally3_smfree( LU, queue  );
    }
    magma_tally3_smfree( &L_d, queue  );
    magma_tally3_smfree( &U_d, queue  );
    magma_tally3_smfree( &LU_d, queue  );
    return info;
}



/**
    Purpose
    -------

    Computes an initial guess for the iterative ILU/IC


    Arguments
    ---------

    @param[in]
    A           magma_tally3_s_matrix
                sparse matrix in CSR

    @param[out]
    L           magma_tally3_s_matrix*
                sparse matrix in CSR

    @param[out]
    U           magma_tally3_s_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.


    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_sinitguess(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix *L,
    magma_tally3_s_matrix *U,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    float one = MAGMA_tally3_S_MAKE( 1.0, 0.0 );
    
    magma_tally3_s_matrix hAL={Magma_tally3_CSR}, hAU={Magma_tally3_CSR}, dAL={Magma_tally3_CSR}, 
    dAU={Magma_tally3_CSR}, dALU={Magma_tally3_CSR}, hALU={Magma_tally3_CSR}, hD={Magma_tally3_CSR}, 
    dD={Magma_tally3_CSR}, dL={Magma_tally3_CSR}, hL={Magma_tally3_CSR};
    magma_tally3_int_t i,j;
    
    magma_tally3_int_t offdiags = 0;
    magma_tally3_index_t *diag_offset;
    float *diag_vals=NULL;

    // need only lower triangular
    hAL.diagorder_type = Magma_tally3_VALUE;
    CHECK( magma_tally3_smconvert( A, &hAL, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
    //magma_tally3_smconvert( hAL, &hALCOO, Magma_tally3_CSR, Magma_tally3_CSRCOO );

    // need only upper triangular
    //magma_tally3_smconvert( A, &hAU, Magma_tally3_CSR, Magma_tally3_CSRU );
    CHECK( magma_tally3_s_cucsrtranspose(  hAL, &hAU, queue ));
    //magma_tally3_smconvert( hAU, &hAUCOO, Magma_tally3_CSR, Magma_tally3_CSRCOO );
    CHECK( magma_tally3_smtransfer( hAL, &dAL, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    CHECK( magma_tally3_s_spmm( one, dAL, dAU, &dALU, queue ));
    CHECK( magma_tally3_smtransfer( dALU, &hALU, Magma_tally3_DEV, Magma_tally3_CPU, queue ));

    magma_tally3_smfree( &dAU, queue);
    magma_tally3_smfree( &dALU, queue);


    CHECK( magma_tally3_smalloc_cpu( &diag_vals, offdiags+1 ));
    CHECK( magma_tally3_index_malloc_cpu( &diag_offset, offdiags+1 ));
    diag_offset[0] = 0;
    diag_vals[0] = MAGMA_tally3_S_MAKE( 1.0, 0.0 );
    CHECK( magma_tally3_smgenerator( hALU.num_rows, offdiags, diag_offset, diag_vals, &hD, queue ));
    magma_tally3_smfree( &hALU, queue );

    
    for(i=0; i<hALU.num_rows; i++){
        for(j=hALU.row[i]; j<hALU.row[i+1]; j++){
            if( hALU.col[j] == i ){
                //printf("%d %d  %d == %d -> %f   -->", i, j, hALU.col[j], i, hALU.val[j]);
                hD.val[i] = MAGMA_tally3_S_MAKE(
                        1.0 / sqrt(fabs(MAGMA_tally3_S_REAL(hALU.val[j])))  , 0.0 );
                //printf("insert %f at %d\n", hD.val[i], i);
            }
        }
    }


    CHECK( magma_tally3_smtransfer( hD, &dD, Magma_tally3_CPU, Magma_tally3_DEV, queue ));
    magma_tally3_smfree( &hD, queue);

    CHECK( magma_tally3_s_spmm( one, dD, dAL, &dL, queue ));
    magma_tally3_smfree( &dAL, queue );
    magma_tally3_smfree( &dD, queue );



/*
    // check for diagonal = 1
    magma_tally3_s_matrix dLt={Magma_tally3_CSR}, dLL={Magma_tally3_CSR}, LL={Magma_tally3_CSR};
    CHECK( magma_tally3_s_cucsrtranspose(  dL, &dLt ));
    CHECK( magma_tally3_scuspmm( dL, dLt, &dLL ));
    CHECK( magma_tally3_smtransfer( dLL, &LL, Magma_tally3_DEV, Magma_tally3_CPU ));
    //for(i=0; i < hALU.num_rows; i++) {
    for(i=0; i < 100; i++) {
        for(j=hALU.row[i]; j < hALU.row[i+1]; j++) {
            if( hALU.col[j] == i ){
                printf("%d %d -> %f   -->", i, i, LL.val[j]);
            }
        }
    }
*/
    CHECK( magma_tally3_smtransfer( dL, &hL, Magma_tally3_DEV, Magma_tally3_CPU, queue ));
    CHECK( magma_tally3_smconvert( hL, L, Magma_tally3_CSR, Magma_tally3_CSRCOO, queue ));



cleanup:
    if( info !=0 ){
        magma_tally3_smfree( L, queue  );
        magma_tally3_smfree( U, queue  );
    }
    magma_tally3_smfree( &dAU, queue);
    magma_tally3_smfree( &dALU, queue);
    magma_tally3_smfree( &dL, queue );
    magma_tally3_smfree( &hL, queue );
    magma_tally3_smfree( &dAL, queue );
    magma_tally3_smfree( &dD, queue );
    magma_tally3_smfree( &hD, queue);
    magma_tally3_smfree( &hALU, queue );
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
    A           magma_tally3_s_matrix*
                sparse matrix in CSR

    @param[out]
    L           magma_tally3_s_matrix*
                sparse matrix in CSR

    @param[out]
    U           magma_tally3_s_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.


    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_sinitrecursiveLU(
    magma_tally3_s_matrix A,
    magma_tally3_s_matrix *B,
    magma_tally3_queue_t queue ){

    magma_tally3_int_t i,j,k;

    for(i=0; i<A.num_rows; i++){
        for(j=B->row[i]; j<B->row[i+1]; j++){
            B->val[j] = MAGMA_tally3_S_MAKE(0.0, 0.0);
            magma_tally3_index_t localcol = B->col[j];
            for( k=A.row[i]; k<A.row[i+1]; k++){
                if(A.col[k] == localcol){
                    B->val[j] = A.val[k];
                }
            }
        }
    }

    return MAGMA_tally3_SUCCESS; 
}



/**
    Purpose
    -------

    Checks for a lower triangular matrix whether it is strictly lower triangular
    and in the negative case adds a unit diagonal. It does this in-place.


    Arguments
    ---------

    @param[in,out]
    L           magma_tally3_s_matrix*
                sparse matrix in CSR
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_saux
    ********************************************************************/

magma_tally3_int_t
magma_tally3_smLdiagadd(
    magma_tally3_s_matrix *L,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;

    magma_tally3_s_matrix LL={Magma_tally3_CSR};

    if( L->row[1]==1 ){        // lower triangular with unit diagonal
        //printf("L lower triangular.\n");
        LL.diagorder_type = Magma_tally3_UNITY;
        CHECK( magma_tally3_smconvert( *L, &LL, Magma_tally3_CSR, Magma_tally3_CSRL, queue ));
    }
    else if( L->row[1]==0 ){ // strictly lower triangular
        //printf("L strictly lower triangular.\n");
        CHECK( magma_tally3_smtransfer( *L, &LL, Magma_tally3_CPU, Magma_tally3_CPU, queue ));
        magma_tally3_free_cpu( LL.col );
        magma_tally3_free_cpu( LL.val );
        LL.nnz = L->nnz+L->num_rows;
        CHECK( magma_tally3_smalloc_cpu( &LL.val, LL.nnz ));
        CHECK( magma_tally3_index_malloc_cpu( &LL.col, LL.nnz ));
        magma_tally3_int_t z=0;
        for( magma_tally3_int_t i=0; i<L->num_rows; i++){
            LL.row[i] = z;
            for( magma_tally3_int_t j=L->row[i]; j<L->row[i+1]; j++){
                LL.val[z] = L->val[j];
                LL.col[z] = L->col[j];
                z++;
            }
            // add unit diagonal
            LL.val[z] = MAGMA_tally3_S_MAKE(1.0, 0.0);
            LL.col[z] = i;
            z++;
        }
        LL.row[LL.num_rows] = z;
        LL.nnz = z;
    }
    else{
        printf("error: L neither lower nor strictly lower triangular!\n");
    }
    magma_tally3_smfree( L, queue );
    CHECK( magma_tally3_smtransfer(LL, L, Magma_tally3_CPU, Magma_tally3_CPU, queue ));

cleanup:
    if( info != 0 ){
        magma_tally3_smfree( L, queue );
    }
    magma_tally3_smfree( &LL, queue );
    return info;
}

