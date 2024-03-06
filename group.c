/*
 * group.c
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

#ifndef GROUP_C_
#define GROUP_C_

#include "group.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


/*writes random values to starting vector in power iteration. using indices of given group */
void writeRandomValues(double* vector, int* group)
{
	int size,i;
	int* currGroupPtr = group;
	size = *currGroupPtr;
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		vector[*currGroupPtr] = rand() - RAND_MAX/2 ;
	}
}

/*writes 0 to array in all indicies denoted by group */
void initializeByGroup(double* arr, int* group)
{
	int i,size;
	int* currGroupPtr = group;
	size = *currGroupPtr;
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		arr[*currGroupPtr] = 0;
	}
}

/*writes 0 to (int) array in all indicies denoted by group. C..*/
void initializeByGroupInt(int* arr,int* group)
{
	int i,size;
	int* currGroupPtr = group;
	size = *currGroupPtr;
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		arr[*currGroupPtr] = 0;
	}
}

/*divides vector by its norma, only in indicies denoted by group */
void normalizeByGroup(double* vector, int* group)
{
	double norm;
	int size,i;
	int* currGroupPtr = group;
	norm = sqrt(dotProductByGroup(vector, vector, group));
	size = *currGroupPtr;
    for(i = 0;i < size;i++)
    {
    	currGroupPtr++;
    	vector[*currGroupPtr] /= norm;
    }

}

/*calculates dot product of two vectors only by indicies denoted by group */
double dotProductByGroup(double* vector1,double* vector2, int* group)
{
	int size,currIndex,i;
	double res;

	int* currGroupPtr = group;
	size = *currGroupPtr;
    res = 0;

    for(i = 0;i < size;i++)
    {
    	currGroupPtr++;
    	currIndex = *currGroupPtr;
    	res += vector1[currIndex] * vector2[currIndex];
    }
    return res;
}

/*return the maximum between the diffrences of the vectors cordinates, by indices denoted by group*/
double maxEpsByGroup(double *vec1, double *vec2, int* group)
{
    double max,currDiff;
    int size,currIndex,i;


	int* currGroupPtr = group;
	size = *currGroupPtr;
	max = 0;

	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		currIndex = *currGroupPtr;
		currDiff = fabs ( vec1[currIndex] - vec2[currIndex]);
		if(max < currDiff)
			max = currDiff;
	}
    return max;
}

/*helper function - puts double(1) in all indexes for s multiplication*/
void trivialDivision(double* s, int* group)
{
	int size,i;
	int* currGroupPtr = group;
	size = *currGroupPtr;
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		s[*currGroupPtr] = 1.0;
	}
}

/*helper function - returns a group of size n with all nodes. index 0 = size of group*/
int* createStartingGroup(int n)
{
	int *group, *currGroupPtr;
	int i;

	group = malloc((n+1)*sizeof(int));
	checkNull(group);
	currGroupPtr = group;
	*currGroupPtr = n;
	for(i=0; i<n;i++)
	{
		currGroupPtr++;
		*currGroupPtr=i;
	}
	return group;
}

/* using given eigenvector, create s division of g group */
void create_eigen_division(double* eigenVector, double* s, int* group)
{
	int size,i,currIndex;
	int* currGroupPtr = group;
	size = *currGroupPtr;

	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		currIndex = *currGroupPtr;
		if(ISPOSITIVE(eigenVector[currIndex]))
			s[currIndex] = 1.0;
		else
			s[currIndex] = -1.0;
	}
}

/* split g into two groups to push to stacks. frees group memory if needed */
void split_group(int* group, double* s, node** stackP, node** stackO,int* stackOCount)
{
	int size, g1Size,i;
	int *g1, *g2, *currg1Ptr, *currg2Ptr;
	int* currGroupPtr = group;
	size = *currGroupPtr;
	g1Size = 0;
	/*check size of g1: */
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		if(s[*currGroupPtr]>0) /* i.e 1.0 or -1.0 */
			g1Size++;
	}
	/*if no need to split:	 */
	if(g1Size == size || g1Size == 0)
	{
		*stackO = push(*stackO,group);
		(*stackOCount)++;
	}
	else
	{
		/* allocate two groups */
		g1 = malloc(sizeof(int)*(g1Size+1));
		checkNull(g1);
		g1[0] = g1Size;
		g2 = malloc(sizeof(int)*( (size-g1Size)+1));
		checkNull(g2);
		g2[0] = size-g1Size;

		/*iterate group and put values into g1 or g2 */
		currg1Ptr = g1 +1;
		currg2Ptr = g2 +1;
		currGroupPtr = group;
		for(i=0;i<size;i++)
		{
			currGroupPtr++;
			if(s[*currGroupPtr] > 0) /* i.e 1.0 or -1.0 */
			{
				*currg1Ptr = *currGroupPtr;
				currg1Ptr++;
			}
			else
			{
				*currg2Ptr = *currGroupPtr;
				currg2Ptr++;
			}
		}
		if(g1[0] == 1)
		{
			*stackO = push(*stackO,g1);
			(*stackOCount)++;
		}
		else
			*stackP = push(*stackP,g1);
		if(g2[0] == 1)
		{
			*stackO = push(*stackO,g2);
			(*stackOCount)++;
		}
		else
			*stackP = push(*stackP,g2);
		free(group);
	}
}

/* converts the group to [0,0,..,1,0,1,..] (size n) etc instead of [1,4,6,8,..] (size group) etc
 * writes into pre Allocated altGroup of size n*/
void altGroupify(int* group, int* altGroup,int num_of_nodes)
{
	/*int currJ,j; */
	int size,i;

	int* currGroupPtr = group;
	size = *currGroupPtr;

	for(i=0;i<num_of_nodes;i++)
	{
		altGroup[i] = 0;
	}
	for(i=0;i<size;i++)
	{
		currGroupPtr++;
		altGroup[*currGroupPtr] = 1;
	}
}

#endif /* GROUP_C_ */
