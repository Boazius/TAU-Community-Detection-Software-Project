/*
 * symSparseMat.c
 * symmetric sparse matrix module summary:
 *
 * module that contains a struct of a symmteric csr sparse matrix, implemented by arrays
 * saves only the upper triangle of values, because it is symmetric.
 *
 * spmat_allocate_symmetric_array				- creates a sparse matrix
 * unAllocate						- frees memory taken by the sparse matrix
 * pop						- pops a node from the stack
 * add_row					- reads an input of nodes , saves the data in itself.
 */

#include <stdio.h>
#include <stdlib.h>
#include "symSparseMat.h"
#include "errorChecking.h"

/*allocate memory for the sparse matrix - three arrays and n*/
symSpArrMat* spmat_allocate_symmetric_array(int n, int nnz)
{
	symSpArrMat *res;

	res = malloc(sizeof(symSpArrMat));
	checkNull(res);

    res->values = calloc(nnz/2, sizeof(double));
    checkNull(res->values);
    res->colind = calloc(nnz/2, sizeof(int));
    checkNull(res->colind);
    res->rowptr = calloc(n + 1, sizeof(int));
    checkNull(res->rowptr);

    res->n = n;
    res->rowptr[n] = 0; /*just to emphasize it is a counter of nnz*/
    return res;
}

/*free resources of spmat based on arrays*/
void unAllocate(symSpArrMat *A) {
    free(A->values);
    free(A->colind);
    free(A->rowptr);
    free(A);
}

/* read input row array into a symmetric sparse matrix, only upper triangle area */
void add_row(symSpArrMat *A, int *inputRow, int i)
{
	int *colind,*rowptr;
	int n,currIndex;
	double *values;

    colind = A->colind;
    rowptr = A->rowptr;
    n = A->n;
    values = A->values;

    currIndex = rowptr[n];
    /*update the row array in spmat with current nnz*/
    rowptr[i] = rowptr[n];
    /*insert values to values, colind, rowptr and update nnz counter*/
    while((*inputRow) != -1)
    {
    	/*save only the upper triangular area of the matrix, because its symmetric */
    	if(*inputRow > i)
    	{
    		values[currIndex] = EDGE_VALUE;
    		colind[currIndex] = *inputRow;
    		currIndex++;
    		rowptr[n]++;
    	}
    	inputRow++;
    }
}
