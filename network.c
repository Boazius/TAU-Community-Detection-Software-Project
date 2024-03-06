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

#ifndef NETWORK_C_
#define NETWORK_C_

#include <stdio.h>
#include <stdlib.h>
#include "network.h"

/*cast int array to double array */
static void copyIntArrToDouble(int* intArr,double* doubleArr,int size)
{
	int i;
	for(i=0;i<size;i++)
	{
		doubleArr[i] = (double)intArr[i];
	}
}


/* allocates and constructs the network using the input file */
Network* read_network_from_file(FILE* inputFile)
{
	Network* network;
	symSpArrMat *A;
	int num_of_nodes,intM,i,currDegree;
	double *degreeArray;
	int *intDegreeArray, *inputVector;;

	/* Allocate the network */
	network = malloc(sizeof(Network));
	checkNull(network);

	/*read first int from file - this is number of nodes in graph n*/

	fread(&num_of_nodes, sizeof(int), 1, inputFile);
	network->n = num_of_nodes;

	/*seek and save k vector of degrees, and M size of all degrees.*/
	degreeArray = malloc(num_of_nodes*sizeof(double));
	checkNull(degreeArray);
	intDegreeArray = malloc(num_of_nodes*sizeof(int));
	checkNull(intDegreeArray);

	intM = 0;

	/*read all degrees of nodes in graph, put in k array degreeArray */
	for (i=0;i<num_of_nodes; i++)
	{
		/*read the i_th degree*/
		fread(&currDegree,sizeof(int),1,inputFile);
		/*jump forward k edges, to next degree*/
		fseek(inputFile,currDegree*sizeof(int),SEEK_CUR);

		/*save the current k into array and into M */
		intDegreeArray[i]=currDegree;
		intM += currDegree;
	}

	/*cast integers to doubles */
	copyIntArrToDouble(intDegreeArray,degreeArray,num_of_nodes);
	network->degreeArray = degreeArray;
	network->M = (double)intM;

	/*like rewind(inputFile): start from beginning of file, and skip first M value:*/
	fseek(inputFile,sizeof(int),SEEK_SET);
	/*create sparse matrix with this M and initial value of non zero elements. this is A matrix*/
	A = spmat_allocate_symmetric_array(num_of_nodes, intM);
	network->A = A;

	/*create vector to read input file edges - size is +1 because of "-1" at the end of input */
	inputVector = malloc((num_of_nodes+1)*sizeof(int));
	checkNull(inputVector);

	/*file reading loop: read k edges into vector, put row into sparse*/
	for(i=0;i<num_of_nodes;i++)
	{
		/*jump over k, to edges list*/
		fseek(inputFile,sizeof(int),SEEK_CUR);

		/*read k_i edges into inputVector - could be 0!*/
		fread(inputVector,sizeof(int),intDegreeArray[i],inputFile);

		/*put -1 as end of input*/
		inputVector[intDegreeArray[i]] = -1;

		/*insert array into symmetric sparse matrix*/
		add_row(A,inputVector,i);
	}
	free(intDegreeArray);
	free(inputVector);
	return network;
}



/*free memory of network and of its sparse matrix */
void free_network(Network* network){
	unAllocate(network->A);
	free(network->degreeArray);
	free(network);
}


#endif /* NETWORK_C_ */
