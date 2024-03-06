/*
 * main.c
 *
 *  contains the main function of the program. opens and closes the files, and sends them to the module.
 */
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "clusterDivider.h"
#define timecheck


/*recieves two arguments - 1st input filename and 2nd output filename */
int main(int argc, char* argv[])
{
	char *inputFileName, *outputFileName;
	FILE *inputFile, *outputFile;

	/*timing clock: */
	#ifdef timecheck
	clock_t start, end;
	start = clock();
	#endif



	/*initialize random number generator with random seed */
	srand(time(NULL));

	/*check if argc is 3 - i.e two arguments*/
	if (argc != 3)
	{
		printf("number of arguments is not 2\n");
		exit(EXIT_FAILURE);
	}

	/*try open files*/
	inputFileName = argv[1];
	inputFile = fopen(inputFileName, "rb");
	if (inputFile == NULL)
	{
		printf("cannot find or open input file - make sure it is not in use and has the correct name\n");
		return 2;
	}
	outputFileName = argv[2];
	outputFile = fopen(outputFileName, "wb");
	if (outputFile == NULL)
	{
		printf("cannot find or open input file - make sure it is not in use and has the correct name\n");
		return 2;
	}

	/*send file streams to clusterDivider */
	divide_network_to_clusters(inputFile, outputFile);

	#ifdef timecheck
	end = clock();
	printf("cluster took %f seconds\n", ((double)(end-start) / CLOCKS_PER_SEC));
	#endif

	/*close files */
	if(fclose(inputFile))
	{
		printf("failed to close input file");
		exit(EXIT_FAILURE);
	}
	if(fclose(outputFile))
	{
		printf("failed to close output file");
		exit(EXIT_FAILURE);
	}

	return EXIT_SUCCESS;

}
