/*
 * network.c
 * Network module summary:
 *
 * module that contains a struct for the sparse matrix A and the degrees array for B.
 *
 *
 * read_network_from_file	- reads an input file and saves the relevant data into data structures
 * free_network				- frees allocated memory of the network structure
 * copyIntArrToDouble		- static helper function to convert an int array into a double array.
 */

#ifndef NETWORK_H_
#define NETWORK_H_

#include "symSparseMat.h"
#include "errorChecking.h"

typedef struct _Network{

	/* Adjacency matrix */
	symSpArrMat* A;

	/* Degrees array */
	double* degreeArray;

	/* Sum of degrees M */
	double M;

	/* Num of vertices in the network */
	int n;

} Network;

/* Read from open input stream and build a network with the sparse matrix etc*/
Network* read_network_from_file(FILE* inputFile);

/* Free all Network memory */
void free_network(Network* network);


#endif /* NETWORK_H_ */
