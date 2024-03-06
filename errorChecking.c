/*
 * errorChecking.c
 * tiny module just to check if mallocs/callocs are successful
 *
 * checkNull - exits the program if a given memory is null
 */

#include "errorChecking.h"
#include <stdio.h>
#include <stdlib.h>

/*function to check if mallocs and callocs etc are successful
 * exits the program if null */
void checkNull(void* Memory)
{
	if (Memory== NULL)
	{
		printf("Error in memory allocation, Exiting Program");
		exit(EXIT_FAILURE);
	}
}
