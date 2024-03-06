/*
 * clusterDivider.c
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


#include "clusterDivider.h"
#include <math.h>
#include "network.h"
#include "stack.h"
#include <time.h>

/* Algorithm 3 , makes the maximal modularity community structures from input, writes to output*/
void divide_network_to_clusters(FILE* inputFile, FILE* outputFile)
{
	Network *network;
	node *stackP, *stackO;
	int stackOCount,n;
	int *group, *altGroup,*unMovedArray;
	double *eigenVector, *columnSums, *absColumnSums, *prevB, *currB, *s,*score,*BVector;
	double eigenValue,deltaQ;


	/* read input file into a new network */
	network = read_network_from_file(inputFile);

	/*create O stack and P stack. */
	stackP = NULL;
	stackP = NULL;
	initialize(&stackP);
	initialize(&stackO);
	stackOCount = 0;
	n = network->n;

	/*starting group vector - has all of the nodes, of size n+1 */
	group = createStartingGroup(network->n);
	stackP = push(stackP,group);

	/* eigenpair allocation */
	eigenValue = 0;

	/* f_i^g allocation - this is the sum of each column in B[g] */
	columnSums = malloc(n * sizeof(double));
	checkNull(columnSums);
	/*save sums of absolute values of each row in matrix , for matrix shift */
	absColumnSums = malloc(n * sizeof(double));
	checkNull(absColumnSums);
	/*allocate n sized vectors for power iteration */
	prevB = malloc(n * sizeof(double));
	checkNull(prevB);
	currB = malloc(n * sizeof(double));
	checkNull(currB);
	/* alternative group allocation: save the group as [0,0,..,1,0,1,..] (size n) */
	altGroup = malloc(n * sizeof(int));
	checkNull(altGroup);
	/*allocate s which is division of size of g. s_i = 1 or -1. double for calculations */
	s = malloc(n * sizeof(double));
	checkNull(s);
	/*allocate arrays score and unMovedArray and BVector for modularity maximization */
	unMovedArray = malloc(n * sizeof(int));
	checkNull(unMovedArray);
	score = malloc(n * sizeof(double));
	checkNull(score);
	BVector = malloc(n * sizeof(double));
	checkNull(BVector);

	/*Algorithm 3 loop: */
	while(!isEmpty(stackP))
	{
		/*pop from P into group and process it */
		stackP = pop(stackP,&group);
		altGroupify(group,altGroup,network->n);

		/*calculate leading eigenpair of B[g], also save f_i^g array of size g */
		getLeadingEigenPair(network,group,altGroup,columnSums,&eigenVector,&eigenValue,
				absColumnSums,prevB,currB);

		if(!ISPOSITIVE(eigenValue))
		{
			/* s = [1,1,1,1..] */
			trivialDivision(s,group);
		}
		else
		{
			/*calculate s such that s_i = 1 if u_i > 0 else s_i = -1 */
			create_eigen_division(eigenVector,s,group);
			/* calculate sT̂B[g]s, using columnSums */
			deltaQ = getDeltaQ(network,s,group,altGroup,columnSums);
			if(!ISPOSITIVE(deltaQ))
			{
				trivialDivision(s,group);
			}
		}
		/* modularity maximize g with calculated s */
		maximizeModularity(group, network, altGroup, s, score, unMovedArray,BVector);

		/*using u, create two new group g1 and g2 */
		split_group(group,s,&stackP,&stackO,&stackOCount);

	}

	/*pops out all groups in stackO, writes them, and frees them */
	write_communities_to_output(stackO,outputFile,stackOCount);
	free(altGroup);
	free(s);
	free(columnSums);
	free(absColumnSums);
	free(prevB);
	free(currB);
	free(score);
	free(unMovedArray);
	free(BVector);
	free_network(network);
}

/* write O stack into output file, frees memory of O */
void write_communities_to_output(node* stackO, FILE* outputFile,int stackOCount)
{
	int *currGroup;
	currGroup = NULL;
	/*write number of clusters to start of file */
	fwrite(&stackOCount,sizeof(int),1,outputFile);
	while(!isEmpty(stackO))
	{
		/*write group in stack O and free its memory*/
		stackO = pop(stackO,&currGroup);
		fwrite(currGroup,sizeof(int),*currGroup+1,outputFile);
		free(currGroup);
	}
}

/*eigenvalue finding: ***********************/

/*gets eigenpair of B_hat[g], using power iteration and matrix shift.
 * also gets columnSums for BHat[g]
 * writes to columnSums, eigenVector, and eigenValue */
void getLeadingEigenPair(Network* network,int* group,int* altGroup,double* columnSums,
		double** eigenVector,double* eigenValue,double* absColumnSums,double* prevB, double* currB)
{
	double *temp;
	double maxDiff,matrixShift;

	/*create b_0 for power iteration */
	writeRandomValues(prevB,group);

	/*get matrix shift amount, and f_i^g for future iteration */
	matrixShift = getMatrixShift(network,group,altGroup,columnSums,absColumnSums);

	/* power iterate B_hat[g] matrix shifted */
	do
	{
		/* from prevB (b_k) , calculate currB (b_k+1), using power iteration */
		powerIterate(prevB,currB,network, group, altGroup, columnSums,matrixShift);
		maxDiff = maxEpsByGroup(prevB,currB,group);
		/*switch currB and prevB for next iteration, prevB is now e.g malloced vector */
		temp = prevB;
		prevB = currB;
		currB = temp;

	} while(ISPOSITIVE(fabs(maxDiff)));

	/* save b_k as eigenVector */
	*eigenVector = prevB;
	*eigenValue = getEigenValue(*eigenVector,network,group,altGroup,matrixShift,columnSums)
				-matrixShift;
}

/* returns matrix shift amount for power iteration. also calculates f_i^g
 * writes to columnSums which is of size n.
 * writes to absColumnSums which is of size n */
double getMatrixShift(Network* network,int* group,int* altGroup,double* columnSums,double* absColumnSums)
{
	double *degreeArray;
	int *colind, *rowptr,*currGroupPtr, *currGroupPtr2,*currColindPtr, *currRowPtr;
	double matrixShift,currVal,currBVal,M;
	int currRow,currCol,i,j,groupSize;

	/*get pointers to network variables, makes code readable */
	colind = network->A->colind;
	rowptr = network->A->rowptr;
	degreeArray = network->degreeArray;
	M = network->M;

	matrixShift=0;
	/*initialize f_i^g and abs(f_i^g) at zeroes*/
	initializeByGroup(absColumnSums,group);
	initializeByGroup(columnSums,group);
	groupSize = *group;
	currGroupPtr = group;
	/*iterate over A*/
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currRowPtr = rowptr + currRow;
		j = *currRowPtr;
		currColindPtr = colind +  *currRowPtr;
		currRowPtr++;
		/*iterate every value in current row in A. it is symmetric so add both to cols and rows */
		for(j=rowptr[currRow];j<*currRowPtr;j++)
		{

			currCol = *currColindPtr;
			/* check if current column is in our group using altGroup */
			if(altGroup[currCol])
			{
				currVal = EDGE_VALUE;
				currBVal = degreeArray[currRow]*degreeArray[currCol]/M ;
				columnSums[currRow] += currVal;
				columnSums[currCol] += currVal;
				absColumnSums[currRow] += fabs(currVal - currBVal) - currBVal;
				absColumnSums[currCol] += fabs(currVal - currBVal) - currBVal;
			}
			currColindPtr++;
		}
	}

	currGroupPtr = group;

	/*iterate over B, without the diagonal */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currGroupPtr2 = group;
		for(j=i+1;j<groupSize;j++)
		{
			currGroupPtr2++;
			currCol = *currGroupPtr2;
			currBVal = degreeArray[currRow]*degreeArray[currCol]/M ;
			columnSums[currRow] -= currBVal;
			columnSums[currCol] -= currBVal;
			absColumnSums[currRow] += currBVal;
			absColumnSums[currCol] += currBVal;
		}
	}
	currGroupPtr = group;
	/*iterate over the diagonal of B */
	for(i=0;i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currBVal = (degreeArray[currRow]*degreeArray[currRow])/M ;
		columnSums[currRow] -= currBVal;
		absColumnSums[currRow] += fabs(currBVal + columnSums[currRow]);
		/* find matrix shift amount - max_j (Sum_i |B_hat[g]| ) */
		if(absColumnSums[currRow]>matrixShift)
			matrixShift = absColumnSums[currRow];
	}
	return matrixShift;
}

/* from prevB (b_k) , calculate currB (b_k+1), using power iteration */
void powerIterate(double* prevB,double* currB,Network* network, int* group,
		int* altGroup,double* columnSums,double matrixShift)
{
	double *degreeArray;
	int *colind, *rowptr, *currGroupPtr, *currGroupPtr2, *currRowPtr, *currColindPtr;
	int currRow,currCol,i,j,groupSize;
	double currVal,currBVal,M;

	/*get pointers to network variables, makes code readable */
	colind = network->A->colind;
	rowptr = network->A->rowptr;
	/*int n = network->A->n;*/
	degreeArray = network->degreeArray;
	M = network->M;
	initializeByGroup(currB,group);

	groupSize = *group;
	currGroupPtr = group;

	/*iterate over A */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currRowPtr = rowptr + currRow;
		j = *currRowPtr;
		currColindPtr = colind +  *currRowPtr;
		currRowPtr++;
		/*iterate every value in current row in A. it is symmetric so add both to cols and rows */
		for(j=rowptr[currRow];j<*currRowPtr;j++)
		{
			currCol = *currColindPtr;
			/* check if current column is in our group using altGroup */
			if(altGroup[currCol])
			{
				currVal = EDGE_VALUE;
				currB[currRow] += currVal * (prevB[currCol]);
				currB[currCol] += currVal * (prevB[currRow]);
			}
			currColindPtr++;
		}
	}

	currGroupPtr = group;
	/*iterate over B without diagonal*/
	for(i=0; i<groupSize ;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currGroupPtr2 = group;

		for(j=i+1;j<groupSize ;j++)
		{
			currGroupPtr2++;
			currCol = *currGroupPtr2;
			currBVal = degreeArray[currRow]*degreeArray[currCol]/M ;
			currB[currRow] -= currBVal * (prevB[currCol]);
			currB[currCol] -= currBVal * (prevB[currRow]);
		}
	}

	currGroupPtr = group;
	/*iterate over the diagonal of B */
	for(i=0;i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currBVal = -(degreeArray[currRow]*degreeArray[currRow])/M ;
		currB[currRow] += (currBVal - columnSums[currRow] + matrixShift) * (prevB[currRow]);
	}
	/* divide currB by its norma */
	normalizeByGroup(currB,group);
}

/*approximates eigenval from eigenVector using formula: λ1=bk·B_hatShifted[g] bk / bk·bk */
double getEigenValue(double* eigenVector,Network* network, int* group, int* altGroup,
		double matrixShift,double* columnSums)
{
	double *degreeArray;
	int *colind, *rowptr, *currGroupPtr, *currGroupPtr2, *currRowPtr, *currColindPtr;
	int currRow,currCol,i,j,groupSize;
	double currVal,currBVal,numerator,denominator,M;

	/*get pointers to network variables, makes code readable */
	colind = network->A->colind;
	rowptr = network->A->rowptr;
	degreeArray = network->degreeArray;
	M = network->M;

	/* division of bk* BHat[g] *bk and bk*bk */
	numerator = 0;

	groupSize = *group;
	currGroupPtr = group;
	/*iterate over A */
	for (i = 0; i < groupSize; i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currRowPtr = rowptr + currRow;
		j = *currRowPtr;
		currColindPtr = colind + *currRowPtr;
		currRowPtr++;
		/*iterate every value in current row in A. it is symmetric so add both to cols and rows */
		for (j = rowptr[currRow]; j < *currRowPtr; j++)
		{
			currCol = *currColindPtr;
			/* check if current column is in our group using altGroup */
			if(altGroup[currCol])
			{
				currVal = EDGE_VALUE*eigenVector[currRow]*eigenVector[currCol];
				numerator += 2 * currVal;
			}
		}
	}

	currGroupPtr = group;
	/*iterate over B */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currGroupPtr2 = group;
		for(j=i+1;j<groupSize ;j++)
		{
			currGroupPtr2++;
			currCol = *currGroupPtr2;
			currBVal = -degreeArray[currRow]*degreeArray[currCol]/M ;
			currVal = currBVal*eigenVector[currRow]*eigenVector[currCol];
			numerator += 2 * currVal;
		}
	}

	currGroupPtr = group;
	/*iterate over the diagonal of B */
	for(i=0;i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currBVal = -(degreeArray[currRow]*degreeArray[currRow])/M - columnSums[currRow] + matrixShift;
		currVal = currBVal*eigenVector[currRow]*eigenVector[currRow];
		numerator += currVal;
	}
	denominator = dotProductByGroup(eigenVector, eigenVector, group);
	/* dont divide by zero. norma of eigen vector must be positive */
	if(!ISPOSITIVE(denominator))
	{
		printf("Zero Division Failure!");
		exit(EXIT_FAILURE);
	}
	return numerator / denominator ;
}

/*delta Q finding: ************************/

/* calculate deltaQ = s^T*BHat^[g]*s, uses precalculated f_i^g */
double getDeltaQ(Network* network, double* s, int* group,int* altGroup, double* columnSums)
{
	double *degreeArray;
	int *colind, *rowptr, * currGroupPtr, * currGroupPtr2, * currColindPtr, * currRowPtr;
	double M,currVal,currBVal,deltaQ;
	int currRow,currCol,i,j, groupSize;

	/*get pointers to network variables, makes code readable */
	colind = network->A->colind;
	rowptr = network->A->rowptr;
	degreeArray = network->degreeArray;
	M = network->M;
	deltaQ = 0;

	groupSize = *group;
	currGroupPtr = group;
	/*iterate over A */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currRowPtr = rowptr + currRow;
		j = *currRowPtr;
		currColindPtr = colind + *currRowPtr;
		currRowPtr++;
		/*iterate every value in current row in A. it is symmetric so add both to cols and rows */
		for (j = rowptr[currRow]; j < *currRowPtr; j++)
		{
			currCol = *currColindPtr;
			/* check if current column is in our group using altGroup */
			if(altGroup[currCol])
			{
				currVal = EDGE_VALUE * s[currRow]*s[currCol];
				deltaQ += 2 * currVal;
			}
		}
	}

	currGroupPtr = group;
	/*iterate over B */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currGroupPtr2 = group;
		for(j=i+1;j<groupSize ;j++)
		{
			currGroupPtr2++;
			currCol = *currGroupPtr2;
			currBVal = -degreeArray[currRow]*degreeArray[currCol]/M ;
			currVal = currBVal*s[currRow]*s[currCol];
			deltaQ += 2 * currVal;
		}
	}

	currGroupPtr = group;
	/*iterate over the diagonal of B */
	for(i=0;i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currBVal = -(degreeArray[currRow]*degreeArray[currRow])/M - columnSums[currRow];
		currVal = currBVal*s[currRow]*s[currRow];
		deltaQ += currVal;
	}
	return deltaQ;
}


/*modularity "optimization": **************************/

/* starting with given s, divide g into two groups (in s) such that it maximizes mod of network
 * writes to s and BVector of size n*/
void maximizeModularity(int* group, Network* network, int* altGroup, double* s,
		double* score, int* unMovedArray, double* BVector)
{
	double deltaQ,M,currBVal,prevImprovePtr;
	int groupSize,maxJ,maxI,i,j,z,currIndex;
	int *indices, *currGroupPtr, *currIndicesPtr ;
	double *improve, *degreeArray, *currImprovePtr;

	degreeArray = network->degreeArray;
	M = network->M;
	groupSize = *group;

	/* array of moved indicies, sorted by time (indices[0] moved first) */
	indices = malloc(sizeof(int)*groupSize);
	/*array that shows the Q (modularity) when indices[i] was moved */
	improve = malloc(sizeof(double)*groupSize);

	do
	{
		currImprovePtr = improve;
		prevImprovePtr = 0;
		currIndicesPtr = indices;
		initializeByGroupInt(unMovedArray,group);
		for(z=0;z<groupSize;z++)
		{
			/*calculate scores smartly:	*/
			matVecMult(network,group,altGroup,s,BVector);
			currGroupPtr = group;
			for(i=0;i<groupSize;i++)
			{
				currGroupPtr++;
				currIndex = *currGroupPtr;
				currBVal = degreeArray[currIndex] * degreeArray[currIndex] / M;
				/*//-4((BS)_k * d_k + D_kk )*/
				score[currIndex] = -4* (BVector[currIndex] * s[currIndex] + currBVal);
			}

			/*move the maximal score node */
			maxJ = getMaxJ(score,group,unMovedArray);
			s[maxJ] = -s[maxJ];
			*currIndicesPtr = maxJ;
			*currImprovePtr = prevImprovePtr + score[maxJ];
			prevImprovePtr = *currImprovePtr;
			currImprovePtr++;
			currIndicesPtr++;
			unMovedArray[maxJ] = 1;
		}

		/*find the state that has the maximum modularity , and go back to it */
		maxI = getMaxI(improve,groupSize);
		currIndicesPtr = indices + groupSize - 1;
		for(i=groupSize-1;i>maxI;i--)
		{
			j = *currIndicesPtr;
			s[j] = -s[j];
			currIndicesPtr--;
		}
		if(maxI == groupSize-1)
			deltaQ = 0;
		else
			deltaQ = improve[maxI];

	} while(ISPOSITIVE(deltaQ));

	free(improve);
	free(indices);
}

/*multiply B matrix by s vector, in rows and columns denoted by group. writes to BsVec vector */
void matVecMult(Network* network, int* group,int* altGroup, double* s, double* BVector)
{
	double *degreeArray;
	int *colind, *rowptr;
	int currRow,currCol,i,j, groupSize;
	double currVal,currBVal,M;

	int* currGroupPtr, *currRowPtr, *currColindPtr;
	initializeByGroup(BVector,group);
	/*get pointers to network variables, makes code readable */
	colind = network->A->colind;
	rowptr = network->A->rowptr;
	degreeArray = network->degreeArray;
	M = network->M;

	groupSize = *group;
	currGroupPtr = group;

	/*iterate over A */
	for(i=0; i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currRowPtr = rowptr + currRow;
		j = *currRowPtr;
		currColindPtr = colind + *currRowPtr;
		currRowPtr++;
		/*iterate every value in current row in A. it is symmetric so add both to cols and rows */
		for(j=rowptr[currRow];j<*currRowPtr;j++)
		{
			currCol = *currColindPtr;
			/* check if current column is in our group using altGroup */
			if(altGroup[currCol])
			{
				currVal = EDGE_VALUE;
				BVector[currRow] += currVal * (s[currCol]);
				BVector[currCol] += currVal * (s[currRow]);
			}
			currColindPtr++;
		}
	}

	currGroupPtr = group;
	/*iterate over B without diagonal*/
	for(i=0; i<groupSize ;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		for(j=i+1;j<groupSize ;j++)
		{
			currCol = group[j+1];
			currBVal = -degreeArray[currRow]*degreeArray[currCol]/M ;
			BVector[currRow] += currBVal * (s[currCol]);
			BVector[currCol] += currBVal * (s[currRow]);
		}
	}

	currGroupPtr = group;

	/*iterate over the diagonal of B */
	for(i=0;i<groupSize;i++)
	{
		currGroupPtr++;
		currRow = *currGroupPtr;
		currBVal = -(degreeArray[currRow]*degreeArray[currRow])/M ;
		BVector[currRow] += currBVal * (s[currRow]);
	}
}


/* argmax helper function for improve array in modularity maximization */
int getMaxI(double* improve,int groupSize)
{
	double maxImprove,currImprove;
	int maxI,i;
	double* currImprovePtr;
	currImprovePtr = improve;

	maxImprove = *currImprovePtr;
	maxI = 0;
	for(i=1;i<groupSize;i++)
	{
		currImprovePtr++;
		currImprove = *currImprovePtr;
		if(currImprove> maxImprove)
		{
			maxImprove = currImprove;
			maxI = i;
		}
	}
	return maxI;

}

/* argmax helper function for score array in modularity maxization*/
int getMaxJ(double* score,int* group, int* unMovedArray)
{
	double maxScore,currScore;
	int flag,groupSize,maxJ,currIndex,i;
	int *currGroupPtr;

	flag = 1;

	currGroupPtr = group;
	groupSize = *currGroupPtr;
	currGroupPtr++;
	maxJ = *currGroupPtr;

	for(i=0;i<groupSize;i++)
	{
		currIndex = *currGroupPtr;
		currScore = score[currIndex];
		if(!unMovedArray[currIndex])
		{
			if(flag)
			{
				maxScore = currScore;
				maxJ = currIndex;
				flag = 0;
			}
			else
			{
				if(maxScore<currScore)
				{
					maxScore = currScore;
					maxJ = currIndex;
				}
			}
		}
		currGroupPtr++;
	}
	return maxJ;
}
