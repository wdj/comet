/*
    -- MAGMA_tally3 (version 1.6.2) --
       Univ. of Tennessee, Knoxville
       Univ. of California, Berkeley
       Univ. of Colorado, Denver
       @date May 2015

       @generated from magma_tally3_zmilustruct.cpp normal z -> c, Sun May  3 11:23:01 2015
       @author Hartwig Anzt
*/

//  in this file, many routines are taken from
//  the IO functions provided by MatrixMarket

#include "common_magma_tally3sparse.h"


/******************************************************************************
 * ILU mex function from MATLAB:
 *
 * [l, u] = ilu_mex(a, level, omega, storage);
 *****************************************************************************/

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define mwIndex magma_tally3_index_t

void magma_tally3_cshell_sort(
  const magma_tally3_int_t n,
  magma_tally3_int_t x[]);

void magma_tally3_csymbolic_ilu(
  const magma_tally3_int_t levinc,
  const magma_tally3_int_t n,
  magma_tally3_int_t *nzl,
  magma_tally3_int_t *nzu,
  const mwIndex *ia,
  const mwIndex *ja,
  mwIndex *ial,
  mwIndex *jal,
  mwIndex *iau,
  mwIndex *jau);


/******************************************************************************
 *
 * MEX function
 *
 *****************************************************************************/

void magma_tally3_cmexFunction(magma_tally3_int_t nlhs, magma_tally3_int_t n, magma_tally3FloatComplex omega,
                 magma_tally3_int_t levfill, magma_tally3_int_t storage,
                magma_tally3_index_t * ial, magma_tally3_index_t *jal, magma_tally3FloatComplex *al,
                magma_tally3_index_t * iau, magma_tally3_index_t *jau, magma_tally3FloatComplex *au,
                magma_tally3_int_t nrhs,
                magma_tally3_index_t * ia, magma_tally3_index_t *ja, magma_tally3FloatComplex *a ){

    /* matrix is stored in CSC format, 0-based */

    magma_tally3_int_t nzl, nzu;

    
    nzl = storage;
    nzu = storage;


    /* the following will fail and return to matlab if insufficient storage */
    magma_tally3_csymbolic_ilu(levfill, n, &nzl, &nzu, ia, ja, ial, jal, iau, jau);

}

/* shell sort
// stable, so it is fast if already sorted
// sorts x[0:n-1] in place, ascending order.
*/

void magma_tally3_cshell_sort(
  const magma_tally3_int_t n, magma_tally3_index_t *x)
{
    magma_tally3_int_t m, max, j, k, itemp;

    m = n/2;

    while (m > 0) {
        max = n - m;
        for (j=0; j<max; j++)
        {
            for (k=j; k>=0; k-=m)
            {
                if (x[k+m] >= x[k])
                    break;
                itemp = x[k+m];
                x[k+m] = x[k];
                x[k] = itemp;
            }
        }
        m = m/2;
    }
}

/*
// symbolic level ILU
// factors magma_tally3_int_to separate upper and lower parts
// sorts the entries in each row of A by index
// assumes no zero rows
*/

void magma_tally3_csymbolic_ilu(
  const magma_tally3_int_t levfill,                 /* level of fill */
  const magma_tally3_int_t n,                       /* order of matrix */
  magma_tally3_int_t *nzl,                          /* input-output */
  magma_tally3_int_t *nzu,                          /* input-output */
  const mwIndex *ia,
  const mwIndex *ja,    /* input */
  mwIndex *ial,
  mwIndex *jal,              /* output lower factor structure */
  mwIndex *iau,
  mwIndex *jau)              /* output upper factor structure */
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_int_t i;
    magma_tally3_index_t *lnklst=NULL;
    magma_tally3_index_t *curlev=NULL;
    magma_tally3_index_t *levels=NULL;
    magma_tally3_index_t *iwork=NULL;
    
    magma_tally3_int_t knzl = 0;
    magma_tally3_int_t knzu = 0;

    CHECK( magma_tally3_index_malloc_cpu( &lnklst, n ));
    CHECK( magma_tally3_index_malloc_cpu( &curlev, n ));
    CHECK( magma_tally3_index_malloc_cpu( &levels, *nzu ));
    CHECK( magma_tally3_index_malloc_cpu( &iwork, n ));
    
    for(magma_tally3_int_t t=0; t<n; t++){
        lnklst[t] = 0;
        curlev[t] = 0;
        iwork[t] = 0;
    }

    for(magma_tally3_int_t t=0; t<*nzu; t++){
        levels[t] = 0;
    }

    ial[0] = 0;
    iau[0] = 0;

    for (i=0; i<n; i++)
    {

     //   printf("check line %d\n", i);
        magma_tally3_int_t first, next, j;

        /* copy column indices of row into workspace and sort them */

        magma_tally3_int_t len = ia[i+1] - ia[i];
        next = 0;
        for (j=ia[i]; j<ia[i+1]; j++)
            iwork[next++] = ja[j];
        magma_tally3_cshell_sort(len, iwork);
     //   printf("check2 line %d\n", i);
        /* construct implied linked list for row */

        first = iwork[0];
        curlev[first] = 0;

        for (j=0; j<=len-2; j++)
        {
            lnklst[iwork[j]] = iwork[j+1];
            curlev[iwork[j]] = 0;
        }
       // printf("check3 line %d iwork[len-1]:%d\n", i, iwork[len-1]);
        lnklst[iwork[len-1]] = n;
        curlev[iwork[len-1]] = 0;

        /* merge with rows in U */
       // printf("check4 line %d lnklst[iwork[len-1]]:%d\n", i, lnklst[iwork[len-1]]);
        next = first;
       // printf("next:%d (!<) first:%d\n", next, i);
        while (next < i)
        {
          //  printf("check line %d while %d\n", i, next);
            magma_tally3_int_t oldlst = next;
            magma_tally3_int_t nxtlst = lnklst[next];
            magma_tally3_int_t row = next;
            magma_tally3_int_t ii;

            /* scan row */

            for (ii=iau[row]+1; ii<iau[row+1]; /*nop*/)
            {
                if (jau[ii] < nxtlst)
                {
                    /* new fill-in */
                    magma_tally3_int_t newlev = curlev[row] + levels[ii] + 1;
                    if (newlev <= levfill)
                    {
                        lnklst[oldlst]  = jau[ii];
                        lnklst[jau[ii]] = nxtlst;
                        oldlst = jau[ii];
                        curlev[jau[ii]] = newlev;
                    }
                    ii++;
                }
                else if (jau[ii] == nxtlst)
                {
            magma_tally3_int_t newlev;
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                    newlev = curlev[row] + levels[ii] + 1;
                    curlev[jau[ii]] = MIN(curlev[jau[ii]], newlev);
                    ii++;
                }
                else /* (jau[ii] > nxtlst) */
                {
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                }
            }
            next = lnklst[next];
        }
        
        /* gather the pattern magma_tally3_int_to L and U */
       // printf("check line5 %d\n", i);
        next = first;
        while (next < i)
        {
            if (knzl >= *nzl)
        {
            printf("ILU: STORAGE parameter value %d<%d too small.\n", *nzl, knzl);
                printf("Increase STORAGE parameter.");
                    info = -1;
                    goto cleanup;
        }
            jal[knzl++] = next;
            next = lnklst[next];
        }
        ial[i+1] = knzl;
      //  printf("check line6 %d\n", i);
        if (next != i)
        {
        printf("ILU structurally singular.\n");
        /*
            assert(knzu < *nzu);
            levels[knzu] = 2*n;
            jau[knzu++] = i;
        */
        }
      //  printf("check line7 %d\n", i);
           //                 printf("next:%d  n:%d \n", next, n);
        while (next < n)
        {
            if (knzu >= *nzu)
        {
            printf("ILU: STORAGE parameter value %d<%d too small.\n", *nzu, knzu);
                printf("Increase STORAGE parameter.");
                info = -1;
                goto cleanup;
        }
                   // printf("1 knzu:%d  next:%d \n", knzu, next );
            levels[knzu] = curlev[next];
                  //  printf("2 knzu:%d  next:%d \n", knzu, next );
            jau[knzu++] = next;
                  //  printf("3 knzu:%d  next:%d \n", knzu, next );
            next = lnklst[next];
                  //  printf("4 next:%d  n:%d \n", next, n);
        }
        iau[i+1] = knzu;
    }



    *nzl = knzl;
    *nzu = knzu;

   // printf("ende\n");

#if 0
    printf( "Actual nnz for ILU: %d\n", *nzl + *nzu );
#endif

cleanup:
    magma_tally3_free_cpu(lnklst);
    magma_tally3_free_cpu(curlev);
    magma_tally3_free_cpu(levels);
    magma_tally3_free_cpu(iwork);

}



/**
    Purpose
    -------

    This routine performs a symbolic ILU factorization.
    The algorithm is taken from an implementation written by Edmond Chow.

    Arguments
    ---------
    @param[in,out]
    A           magma_tally3_c_matrix*
                matrix in magma_tally3 sparse matrix format containing the original
                matrix on input, and L,U on output

    @param[in]
    levels      magma_tally3_magma_tally3_int_t_t
                fill in level

    @param[out]
    L           magma_tally3_c_matrix*
                output lower triangular matrix in magma_tally3 sparse matrix format
                empty on function call

    @param[out]
    U           magma_tally3_c_matrix*
                output upper triangular matrix in magma_tally3 sparse matrix format
                empty on function call
                
    @param[in]
    queue       magma_tally3_queue_t
                Queue to execute in.

    @ingroup magma_tally3sparse_caux
    ********************************************************************/

extern "C"
magma_tally3_int_t
magma_tally3_csymbilu(
    magma_tally3_c_matrix *A,
    magma_tally3_int_t levels,
    magma_tally3_c_matrix *L,
    magma_tally3_c_matrix *U,
    magma_tally3_queue_t queue )
{
    magma_tally3_int_t info = 0;
    
    magma_tally3_c_matrix A_copy={Magma_tally3_CSR}, B={Magma_tally3_CSR};
    magma_tally3_c_matrix hA={Magma_tally3_CSR}, CSRCOOA={Magma_tally3_CSR};
    
    if( A->memory_location == Magma_tally3_CPU && A->storage_type == Magma_tally3_CSR ){

        CHECK( magma_tally3_cmtransfer( *A, &A_copy, Magma_tally3_CPU, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_cmtransfer( *A, &B, Magma_tally3_CPU, Magma_tally3_CPU, queue ));

        // possibility to scale to unit diagonal
        //magma_tally3_cmscale( &B, Magma_tally3_UNITDIAG );

        CHECK( magma_tally3_cmconvert( B, L, Magma_tally3_CSR, Magma_tally3_CSR , queue));
        CHECK( magma_tally3_cmconvert( B, U, Magma_tally3_CSR, Magma_tally3_CSR, queue ));

        magma_tally3_int_t num_lnnz = (levels > 0 ) ? B.nnz/2*(2*levels+50) : B.nnz;
        magma_tally3_int_t num_unnz = (levels > 0 ) ? B.nnz/2*(2*levels+50) : B.nnz;

        magma_tally3_free_cpu( L->col );
        magma_tally3_free_cpu( U->col );
        CHECK( magma_tally3_index_malloc_cpu( &L->col, num_lnnz ));
        CHECK( magma_tally3_index_malloc_cpu( &U->col, num_unnz ));

        magma_tally3_csymbolic_ilu( levels, A->num_rows, &num_lnnz, &num_unnz, B.row, B.col,
                                            L->row, L->col, U->row, U->col );
        L->nnz = num_lnnz;
        U->nnz = num_unnz;
        magma_tally3_free_cpu( L->val );
        magma_tally3_free_cpu( U->val );
        CHECK( magma_tally3_cmalloc_cpu( &L->val, L->nnz ));
        CHECK( magma_tally3_cmalloc_cpu( &U->val, U->nnz ));
        for( magma_tally3_int_t i=0; i<L->nnz; i++ )
            L->val[i] = MAGMA_tally3_C_MAKE( 0.0, 0.0 );

        for( magma_tally3_int_t i=0; i<U->nnz; i++ )
            U->val[i] = MAGMA_tally3_C_MAKE( 0.0, 0.0 );
        // take the original values (scaled) as initial guess for L
        for(magma_tally3_int_t i=0; i<L->num_rows; i++){
            for(magma_tally3_int_t j=B.row[i]; j<B.row[i+1]; j++){
                magma_tally3_index_t lcol = B.col[j];
                for(magma_tally3_int_t k=L->row[i]; k<L->row[i+1]; k++){
                    if( L->col[k] == lcol ){
                        L->val[k] =  B.val[j];
                    }
                }
            }
        }

        // take the original values (scaled) as initial guess for U
        for(magma_tally3_int_t i=0; i<U->num_rows; i++){
            for(magma_tally3_int_t j=B.row[i]; j<B.row[i+1]; j++){
                magma_tally3_index_t lcol = B.col[j];
                for(magma_tally3_int_t k=U->row[i]; k<U->row[i+1]; k++){
                    if( U->col[k] == lcol ){
                        U->val[k] =  B.val[j];
                    }
                }
            }
        }
        magma_tally3_cmfree( &B, queue );
        // fill A with the new structure;
        magma_tally3_free_cpu( A->col );
        magma_tally3_free_cpu( A->val );
        CHECK( magma_tally3_index_malloc_cpu( &A->col, L->nnz+U->nnz ));
        CHECK( magma_tally3_cmalloc_cpu( &A->val, L->nnz+U->nnz ));
        A->nnz = L->nnz+U->nnz ;
        
        magma_tally3_int_t z = 0;
        for(magma_tally3_int_t i=0; i<A->num_rows; i++){
            A->row[i] = z;
            for(magma_tally3_int_t j=L->row[i]; j<L->row[i+1]; j++){
                A->col[z] = L->col[j];
                A->val[z] = L->val[j];
                z++;
            }
            for(magma_tally3_int_t j=U->row[i]; j<U->row[i+1]; j++){
                A->col[z] = U->col[j];
                A->val[z] = U->val[j];
                z++;
            }
        }
        A->row[A->num_rows] = z;
        // reset the values of A to the original entries
        for(magma_tally3_int_t i=0; i<A->num_rows; i++){
            for(magma_tally3_int_t j=A_copy.row[i]; j<A_copy.row[i+1]; j++){
                magma_tally3_index_t lcol = A_copy.col[j];
                for(magma_tally3_int_t k=A->row[i]; k<A->row[i+1]; k++){
                    if( A->col[k] == lcol ){
                        A->val[k] =  A_copy.val[j];
                    }
                }
            }
        }
    }
    else{
        magma_tally3_storage_t A_storage = A->storage_type;
        magma_tally3_location_t A_location = A->memory_location;
        CHECK( magma_tally3_cmtransfer( *A, &hA, A->memory_location, Magma_tally3_CPU, queue ));
        CHECK( magma_tally3_cmconvert( hA, &CSRCOOA, hA.storage_type, Magma_tally3_CSR, queue ));

        CHECK( magma_tally3_csymbilu( &CSRCOOA, levels, L, U, queue ));

        magma_tally3_cmfree( &hA, queue );
        magma_tally3_cmfree( A, queue );
        CHECK( magma_tally3_cmconvert( CSRCOOA, &hA, Magma_tally3_CSR, A_storage, queue ));
        CHECK( magma_tally3_cmtransfer( hA, A, Magma_tally3_CPU, A_location, queue ));
    }
    
cleanup:
    if( info != 0 ){
        magma_tally3_cmfree( L, queue );
        magma_tally3_cmfree( U, queue );
    }
    magma_tally3_cmfree( &A_copy, queue );
    magma_tally3_cmfree( &B, queue );
    magma_tally3_cmfree( &hA, queue );
    magma_tally3_cmfree( &CSRCOOA, queue );
    return info;
}
