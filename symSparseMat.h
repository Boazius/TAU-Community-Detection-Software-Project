/*
 * symSparseMat.h
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
#ifndef SYMSPARSEMAT_H_
#define SYMSPARSEMAT_H_

/* no magic numbers, value of edges between nodes, default: 1 */
#define EDGE_VALUE  1


/* Definition of symmetric sparse matrix represented by arrays*/
typedef struct symSpArrMat {
	/* Matrix size (n*n) */
	int		n;
	/*csr matrix arrays: */
    double *values;
    int *colind;
    int *rowptr;
} symSpArrMat;


/* Adds row i the symmetric matrix, from input vector */
void add_row(symSpArrMat *A, int *inputRow, int i);

/* Frees all resources used by A */
void	unAllocate(symSpArrMat *A);


/* Allocates a new arrays symmetric sparse matrix of size n with nnz non-zero elements*/
symSpArrMat* spmat_allocate_symmetric_array(int n, int nnz);

#endif /* SYMSPARSEMAT_H_ */
