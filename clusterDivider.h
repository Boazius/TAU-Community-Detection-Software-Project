/*
 * clusterDivider.h
 * Cluster Divider Module Summary:
 *
 * The central module that controls the flow of the program.
 * gets input file with a graph with nodes and edges.
 * divides the graph into community structures using its leading eigen pair,
 * also tries to move nodes between the resulting division in order to maximize modularity
 * outputs the resulting division into output file in binary.
 *
 *
 * divide_network_to_clusters	- gets network in input file and divides into clusters
 * write_communities_to_output 	- writes into output file, the clusters calculated.
 * getLeadingEigenPair			- calculates for given clusters its leading eigen pair
 * getMatrixShift 				- used for the power iteration, calculates shift amount for the given cluster
 * powerIterate					- calculates one iteration of the power iteration algorithm
 * getEigenValue 				- approximates eigen value from eigen vector
 * getDeltaQ					- calculates s^T * B * s. on a given cluster
 * maximizeModularity			- algorithm that moves nodes in a given division to maximize its modularity
 * matVecMult					- matrix vector multiplication for a given cluster
 * getMaxI, getMaxJ				- used for maximize modularity
 *
 */

#ifndef CLUSTERDIVIDER_H_
#define CLUSTERDIVIDER_H_

#include <stdio.h>
#include <stdlib.h>
#include "network.h"
#include "stack.h"
#include "group.h"
#include "errorChecking.h"


/* Algorithm 3 , makes the maximal modularity community structures from input, writes to output*/
void divide_network_to_clusters(FILE* inputFile, FILE* outputFile);

/* write O stack into output file */
void write_communities_to_output(node* O, FILE* outputFile,int stackOCount);



/*eigenvalue finding: ***************************************************************/

/*gets eigenpair of B_hat[g], using power iteration and matrix shift */
void getLeadingEigenPair(Network* network,int* group,int* altGroup,double* columnSums,
		double** eigenVector,double* eigenValue,double* absColumnSums,double* prevB, double* currB);

/* returns matrix shift amount for power iteration. also calculates f_i^g
 * writes to columnSums which is of size n.
 * writes to absColumnSums which is of size n */
double getMatrixShift(Network* network,int* group,int* altGroup,double* columnSums,double* absColumnSums);

/* from prevB (b_k) , calculate currB (b_k+1), using power iteration */
void powerIterate(double* prevB,double* currB,Network* network, int* group,
		int* altGroup,double* columnSums,double matrixShift);

/*approximates eigenval from eigenVector using formula: λ1=bk·B_hatShifted[g] bk / bk·bk */
double getEigenValue(double* eigenVector,Network* network, int* group, int* altGroup,
		double matrixShift,double* columnSums);

/*delta Q finding: ******************************************************************/

/* calculate deltaQ = s^T*B^[g]*s, uses precalculated f_i^g */
double getDeltaQ(Network* network, double* s, int* group,int* altGroup, double* columnSums);

/*modularity optimization: **********************************************************/

/* starting with given s, divide g into two groups (in s) such that it maximizes mod of network
 * writes to s and BVector of size n*/
void maximizeModularity(int* group, Network* network, int* altGroup, double* s,
		double* score, int* unMovedArray,double *BVector);

/*multiply B matrix by s vector, in rows and columns denoted by group. writes to BsVec vector */
void matVecMult(Network* network, int* group,int* altGroup, double* s, double* BsVec);

/* argmax helper function for improve array in modularity maximization */
int getMaxI(double* improve,int groupSize);

/* argmax helper function for score array in modularity maxization
 * writes to MaxJIdx*/
int getMaxJ(double* score,int* group, int* unMovedArray);



#endif /* CLUSTERDIVIDER_H_ */
