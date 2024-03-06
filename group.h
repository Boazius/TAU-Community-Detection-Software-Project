/*
 * group.h
 * Group module summary:
 *
 * module that has a few helpful functions that iterate an array using a given group (cluster)
 * a group array is built like so: first index is its size, the rest are a sorted list of its nodes.
 *
 *
 * writeRandomValues			- writes random values from -RAND_MAX/2 to RAND_MAX/2 into an array
 * initializeByGroup
 * initializeByGroupInt			- writes 0 into an array, in indicies denoted by group
 * normalizeByGroup				- divides a vector by its norma, only in group related indicies
 * dotProductByGroup			- multiplies vectors only in group related indicies
 * maxEpsByGroup				- calculates maximum difference between two vectors indicies in group
 * trivialDivision				- writes 1 into an array, in indicies denoted by group
 * createStartingGroup			- creates the first group, with all nodes from 0 to n
 * create_eigen_division		- using a given eigenVector, calculates an s division of a gven cluster
 * split_group					- gets two stacks, and an s division, and
 * 								- pushes a division of the given group into the relevant stacks
 * altGroupify			- creates an alternate representation of the group - of size n, with 1's and 0's
 */

#ifndef GROUP_H_
#define GROUP_H_

#include "stack.h"
#include "errorChecking.h"



/*writes random values to starting vector in power iteration. using indices of given group */
void writeRandomValues(double* vector, int* group);

/*writes 0 to array in all indicies denoted by group */
void initializeByGroup(double* arr, int* group);

/*writes 0 to (int) array in all indicies denoted by group */
void initializeByGroupInt(int* arr,int* group);

/*divides vector by its norma, only in indicies denoted by group */
void normalizeByGroup(double* vector, int* group);

/*calculates dot product of two vectors only by indicies denoted by group */
double dotProductByGroup(double* vector1,double* vector2, int* group);

/*return the maximum between the diffrences of the vectors cordinates, by indices denoted by group*/
double maxEpsByGroup(double *vec1, double *vec2, int* group);

/*helper function - puts double(1) in all indexes for s multiplication*/
void trivialDivision(double* s, int* group);

/*helper function - returns a group of size n with all nodes. index 0 = size of group*/
int* createStartingGroup(int n);

/* using given eigenvector, create s division of g group */
void create_eigen_division(double* eigenVector, double* s, int* group);

/* split g into two groups to push to stacks. frees group memory if needed */
void split_group(int* group, double* s, node** stackP, node** stackO, int* stackOCount);

/* converts the group to [0,0,..,1,0,1,..] (size n) etc instead of [1,4,6,8,..] (size group) etc
 * writes into pre Allocated altGroup of size n*/
void altGroupify(int* group, int* altGroup,int num_of_nodes);


#endif /* GROUP_H_ */
