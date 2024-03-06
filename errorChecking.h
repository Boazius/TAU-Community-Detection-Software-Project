/*
 * errorChecking.h
 * tiny module just to check if mallocs/callocs are successful
 *
 * checkNull - exits the program if a given memory is null
 */

#ifndef ERRORCHECKING_H_
#define ERRORCHECKING_H_

#define epsilon 0.00001
#define ISPOSITIVE(X) ((X) > 0.00001)

/*function to check if mallocs and callocs etc are successful
 * exits the program if null */
void checkNull(void* Memory);


#endif /* ERRORCHECKING_H_ */
